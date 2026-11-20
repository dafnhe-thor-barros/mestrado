"""auditoria_marco4_robusta.py -- Auditoria robusta do Marco 4 (sem alterar modelo_pinn.py):
  1. 500 casos Latin Hypercube (300 no dominio v2 + 200 fora/fronteira por distancia normalizada a caixa de treino)
  2. Metricas oficiais com sinal (pior subestimacao do pico de Tj, p05, p95, mediana, fracao <= 0.5 C)
  3. Curva OOD: erro vs. distancia normalizada a caixa de treino
  4. Concorrentes naturais: modelos concentrados RC de 2 nos e 3 nos (mesmas constantes fisicas, sem dados)
  5. Benchmark de tempo do solver EDP (N=80 e N=160) vs. PINN v2c.pt (fator de aceleracao real)
"""
import json
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats.qmc import LatinHypercube
import torch

DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(DIR, "..", ".."))
FIG_DIR = os.path.join(ROOT, "figuras")
sys.path.insert(0, DIR)

import modelo_termico as mt
import modelo_pinn as mp


def box_distance_v2(case):
    """Distancia euclidiana normalizada a caixa de treino correlacionada da v2 (0 se dentro do dominio)."""
    f, qa, qr, tj0, tr0, dtip = case
    djr = tj0 - tr0
    bounds = [
        (f, mp.RNG_["f"][0], mp.RNG_["f"][1]),
        (qa, mp.RNG_["q_abs"][0], mp.RNG_["q_abs"][1]),
        (qr, mp.RNG_["q_rest"][0], mp.RNG_["q_rest"][1]),
        (djr, mp.DJR[0], mp.DJR[1]),
        (tr0, mp.RNG_["tr0"][0], mp.RNG_["tr0"][1]),
        (dtip, mp.RNG_["dtip"][0], mp.RNG_["dtip"][1]),
    ]
    sq = 0.0
    for v, lo, hi in bounds:
        span = hi - lo
        if v < lo:
            sq += ((lo - v) / span) ** 2
        elif v > hi:
            sq += ((v - hi) / span) ** 2
    return float(np.sqrt(sq))


def sample_lhs_in_domain(n=300, seed=42):
    u = LatinHypercube(d=6, seed=seed).random(n)
    cases = []
    for row in u:
        f = mp.RNG_["f"][0] + row[0] * (mp.RNG_["f"][1] - mp.RNG_["f"][0])
        qa = mp.RNG_["q_abs"][0] + row[1] * (mp.RNG_["q_abs"][1] - mp.RNG_["q_abs"][0])
        qr = mp.RNG_["q_rest"][0] + row[2] * (mp.RNG_["q_rest"][1] - mp.RNG_["q_rest"][0])
        djr = mp.DJR[0] + row[3] * (mp.DJR[1] - mp.DJR[0])
        tr0 = mp.RNG_["tr0"][0] + row[4] * (mp.RNG_["tr0"][1] - mp.RNG_["tr0"][0])
        dtip = mp.RNG_["dtip"][0] + row[5] * (mp.RNG_["dtip"][1] - mp.RNG_["dtip"][0])
        cases.append((float(f), float(qa), float(qr), float(tr0 + djr), float(tr0), float(dtip)))
    return cases


def sample_lhs_ood(n=200, seed=99):
    """Amostragem LHS na caixa v1 independente + extensoes em (Tj0 - Tr0) para mapear OOD por distancia."""
    u = LatinHypercube(d=6, seed=seed).random(n)
    cases = []
    lo = mp.LO.numpy()
    hi = mp.HI.numpy()
    for row in u:
        p = lo + row * (hi - lo)
        # Forca alguns casos fora do intervalo DJR [5, 40] K (inclusive Tj0 < Tr0 ou Tj0 - Tr0 > 40)
        c = tuple(float(x) for x in p)
        if box_distance_v2(c) > 0.0:
            cases.append(c)
    # Completa ate n casos estritamente OOD
    rng = np.random.default_rng(seed + 1)
    while len(cases) < n:
        p = lo + rng.random(6) * (hi - lo)
        c = tuple(float(x) for x in p)
        if box_distance_v2(c) > 0.0:
            cases.append(c)
    return cases[:n]


def rc2_window(case, ph, t_eval):
    """Concorrente natural 1: modelo RC concentrado de 2 nos (Tj + no unico radiador/coletor), fisica pura."""
    f, qa, qr, tj0, tr0, dtip = case
    # Temperatura media inicial da aleta com perfil quadratico: T_bar = tr0 - (2/3)*dtip
    t_fin0 = tr0 - (2.0 / 3.0) * dtip
    c_fin = ph.a_pl * ph.ca
    c_tot = ph.cfl + c_fin
    t_rad0 = (ph.cfl * tr0 + c_fin * t_fin0) / c_tot
    n_mod = ph.p["n_mod"]

    def rhs(t, y):
        tj, trad = y
        pm = mt.ia_p(f, tj - 273.15)[0]
        q_jm = (tj - trad) / ph.R
        dtj = (pm - q_jm) / ph.cmod
        q_in = n_mod * (q_jm + ph.kap * pm) + qr + ph.a_pl * qa
        q_out = ph.a_pl * ph.es * (trad ** 4 - ph.ts ** 4)
        dtrad = (q_in - q_out) / c_tot
        return [dtj, dtrad]

    sol = solve_ivp(rhs, (0.0, mp.W_S), [tj0, t_rad0], t_eval=t_eval, method="Radau", rtol=1e-6, atol=1e-8)
    return sol.y[0]


def rc3_window(case, ph, t_eval):
    """Concorrente natural 2: modelo RC concentrado de 3 nos (Tj + raiz/coletor + aleta media), fisica pura."""
    f, qa, qr, tj0, tr0, dtip = case
    t_fin0 = tr0 - (2.0 / 3.0) * dtip
    c_fin = ph.a_pl * ph.ca
    # Condutancia efetiva raiz -> centroide da aleta (distancia L/3 para perfil parabolico de fluxo)
    g_fin = ph.ell * ph.kd / (ph.L / 3.0)
    n_mod = ph.p["n_mod"]

    def rhs(t, y):
        tj, tr, tf = y
        pm = mt.ia_p(f, tj - 273.15)[0]
        q_jm = (tj - tr) / ph.R
        dtj = (pm - q_jm) / ph.cmod
        q_rf = g_fin * (tr - tf)
        dtr = (n_mod * (q_jm + ph.kap * pm) + qr - q_rf) / ph.cfl
        q_rad = ph.a_pl * ph.es * (tf ** 4 - ph.ts ** 4)
        dtf = (q_rf + ph.a_pl * qa - q_rad) / c_fin
        return [dtj, dtr, dtf]

    sol = solve_ivp(rhs, (0.0, mp.W_S), [tj0, tr0, t_fin0], t_eval=t_eval, method="Radau", rtol=1e-6, atol=1e-8)
    return sol.y[0]


def benchmark_solver(case, n_runs=20):
    res = {}
    for label, N, dt in [("N80_dt2", 80, 2.0), ("N160_dt1", 160, 1.0)]:
        ts = []
        for _ in range(n_runs):
            t0 = time.perf_counter()
            mp.solver_window(case, mp.A_PL, N=N, dt=dt)
            ts.append((time.perf_counter() - t0) * 1e3)
        ts = np.array(ts)
        res[label] = {"mean_ms": float(ts.mean()), "p50_ms": float(np.percentile(ts, 50)), "p99_ms": float(np.percentile(ts, 99))}
    return res


def main():
    torch.set_num_threads(1)
    ph = mp.Phys()
    model = mp.load(os.path.join(DIR, "v2c.pt"))

    in_cases = sample_lhs_in_domain(300, seed=42)
    ood_cases = sample_lhs_ood(200, seed=99)

    t_eval = np.linspace(0.0, mp.W_S, 151)  # dt = 2.0 s
    taus = torch.tensor(t_eval / mp.W_S, dtype=torch.float32)

    def eval_set(cases, run_rc=False):
        out = []
        for c in cases:
            _, Y = mp.solver_window(c, mp.A_PL, N=80, dt=2.0)
            tjr = Y[:, 0]
            tjp, trp, _ = mp.predict_state(model, torch.tensor(c, dtype=torch.float32), taus)
            tjp = tjp.numpy().astype(float)
            dj = tjp - tjr
            rec = {
                "dist": box_distance_v2(c),
                "tj_max_abs_C": float(np.abs(dj).max()),
                "tj_rms_K": float(np.sqrt(np.mean(dj ** 2))),
                "tj_peak_err_C": float(tjp.max() - tjr.max()),
                "tj_final_err_C": float(dj[-1]),
                "tj_peak_ref_C": float(tjr.max() - 273.15),
            }
            if run_rc:
                tj_rc2 = rc2_window(c, ph, t_eval)
                tj_rc3 = rc3_window(c, ph, t_eval)
                rec["rc2_max_abs_C"] = float(np.abs(tj_rc2 - tjr).max())
                rec["rc2_peak_err_C"] = float(tj_rc2.max() - tjr.max())
                rec["rc3_max_abs_C"] = float(np.abs(tj_rc3 - tjr).max())
                rec["rc3_peak_err_C"] = float(tj_rc3.max() - tjr.max())
            out.append(rec)
        return out

    t0 = time.perf_counter()
    res_in = eval_set(in_cases, run_rc=True)
    res_ood = eval_set(ood_cases, run_rc=False)
    elapsed = time.perf_counter() - t0

    def summarize_pinn(rows):
        ma = np.array([r["tj_max_abs_C"] for r in rows])
        pk = np.array([r["tj_peak_err_C"] for r in rows])
        rm = np.array([r["tj_rms_K"] for r in rows])
        return {
            "n": len(rows),
            "tj_max_abs_mediana_C": float(np.median(ma)),
            "tj_max_abs_p95_C": float(np.percentile(ma, 95)),
            "tj_max_abs_pior_C": float(ma.max()),
            "tj_rms_mediana_K": float(np.median(rm)),
            "frac_le_0p5C": float(np.mean(ma <= 0.5)),
            "n_le_0p5C": int(np.sum(ma <= 0.5)),
            "pico_Tj_subestimacao_p05_C": float(np.percentile(pk, 5)),
            "pico_Tj_pior_subestimacao_C": float(pk.min()),
            "pico_Tj_vies_medio_C": float(pk.mean()),
        }

    def summarize_rc(rows, key_ma, key_pk):
        ma = np.array([r[key_ma] for r in rows])
        pk = np.array([r[key_pk] for r in rows])
        return {
            "tj_max_abs_mediana_C": float(np.median(ma)),
            "tj_max_abs_p95_C": float(np.percentile(ma, 95)),
            "tj_max_abs_pior_C": float(ma.max()),
            "frac_le_0p5C": float(np.mean(ma <= 0.5)),
            "pico_Tj_pior_subestimacao_C": float(pk.min()),
        }

    # Bins OOD por distancia normalizada
    bins = [(0.0, 0.0, "d=0 (in-domain)"), (1e-9, 0.15, "0 < d <= 0.15"), (0.15, 0.35, "0.15 < d <= 0.35"), (0.35, 10.0, "d > 0.35")]
    all_rows = res_in + res_ood
    ood_curve = []
    for lo_d, hi_d, lbl in bins:
        if lo_d == 0.0 and hi_d == 0.0:
            sub = [r for r in all_rows if r["dist"] == 0.0]
        else:
            sub = [r for r in all_rows if lo_d <= r["dist"] <= hi_d]
        if sub:
            ma = np.array([r["tj_max_abs_C"] for r in sub])
            pk = np.array([r["tj_peak_err_C"] for r in sub])
            ood_curve.append({
                "faixa": lbl, "n": len(sub),
                "mediana_C": float(np.median(ma)), "p95_C": float(np.percentile(ma, 95)),
                "pior_C": float(ma.max()), "pior_subestimacao_pico_C": float(pk.min()),
            })

    case_a = (1353.752, 135.2, 5551.99, 351.0, 321.0, 6.0)
    solver_bench = benchmark_solver(case_a, n_runs=15)
    pinn_lat = mp.measure_latency(model, n=1000, nt=16, threads=1)

    report = {
        "tempo_auditoria_s": round(elapsed, 2),
        "pinn_in_domain_LHS300": summarize_pinn(res_in),
        "pinn_ood_LHS200": summarize_pinn(res_ood),
        "rc2_in_domain_LHS300": summarize_rc(res_in, "rc2_max_abs_C", "rc2_peak_err_C"),
        "rc3_in_domain_LHS300": summarize_rc(res_in, "rc3_max_abs_C", "rc3_peak_err_C"),
        "curva_ood_por_distancia": ood_curve,
        "benchmark_solver_cpu_ms": solver_bench,
        "benchmark_pinn_cpu_ms": pinn_lat,
        "speedup_vs_N80_dt2": float(solver_bench["N80_dt2"]["mean_ms"] / pinn_lat["mean_ms"]),
        "speedup_vs_N160_dt1": float(solver_bench["N160_dt1"]["mean_ms"] / pinn_lat["mean_ms"]),
    }

    out_json = os.path.join(DIR, "saida_marco4_auditoria_500_rc.json")
    open(out_json, "w", encoding="utf-8", newline="\n").write(json.dumps(report, indent=1, ensure_ascii=False) + "\n")
    print(json.dumps(report, indent=1, ensure_ascii=False))

    # Figura de publicacao (3 paineis): (a) Histograma erro com sinal (n=300); (b) OOD vs distancia; (c) PINN vs RC2 vs RC3
    plt.rcParams.update({"font.family": "serif", "font.size": 9.5, "axes.grid": True, "grid.alpha": 0.3, "grid.linestyle": "--"})
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12.5, 3.8))

    pk_in = np.array([r["tj_peak_err_C"] for r in res_in])
    ax1.hist(pk_in, bins=25, color="#1f77b4", alpha=0.85, edgecolor="white")
    ax1.axvline(0.0, color="black", lw=1.0)
    ax1.axvline(pk_in.min(), color="#d62728", ls="--", lw=1.3, label=rf"Pior subest.: {pk_in.min():.2f}$^\circ$C")
    ax1.axvline(np.percentile(pk_in, 5), color="#ff7f0e", ls=":", lw=1.3, label=rf"p05: {np.percentile(pk_in,5):.2f}$^\circ$C")
    ax1.set_xlabel(r"Erro no pico de $T_j$ ($\hat{T}_j - T_{j,\mathrm{ref}}$) [$^\circ$C]")
    ax1.set_ylabel("Número de casos (LHS, n = 300)")
    ax1.set_title("(a) Distribuição do erro com sinal")
    ax1.legend(fontsize=7.8, loc="upper left")

    ds = np.array([r["dist"] for r in all_rows])
    ms = np.array([r["tj_max_abs_C"] for r in all_rows])
    ax2.scatter(ds[ds == 0], ms[ds == 0], s=10, alpha=0.5, color="#2ca02c", label="In-domain (d = 0)")
    ax2.scatter(ds[ds > 0], ms[ds > 0], s=12, alpha=0.6, color="#d62728", label="OOD (d > 0)")
    ax2.axhline(0.5, color="black", ls="--", lw=1.0, label=r"Meta $0{,}5^\circ$C")
    ax2.set_xlabel("Distância normalizada à caixa de treino ($d_{\\mathrm{box}}$)")
    ax2.set_ylabel(r"Erro máximo absoluto em $T_j$ ($^\circ$C)")
    ax2.set_title("(b) Degradação fora do domínio (n = 500)")
    ax2.legend(fontsize=7.8, loc="upper left")

    models = ["PINN v2c\n(física pura)", "RC 3 nós\n(física pura)", "RC 2 nós\n(física pura)"]
    meds = [report["pinn_in_domain_LHS300"]["tj_max_abs_mediana_C"],
            report["rc3_in_domain_LHS300"]["tj_max_abs_mediana_C"],
            report["rc2_in_domain_LHS300"]["tj_max_abs_mediana_C"]]
    p95s = [report["pinn_in_domain_LHS300"]["tj_max_abs_p95_C"],
            report["rc3_in_domain_LHS300"]["tj_max_abs_p95_C"],
            report["rc2_in_domain_LHS300"]["tj_max_abs_p95_C"]]
    x = np.arange(3)
    ax3.bar(x - 0.18, meds, width=0.35, color="#1f77b4", label="Mediana")
    ax3.bar(x + 0.18, p95s, width=0.35, color="#ff7f0e", label="Percentil 95%")
    ax3.axhline(0.5, color="#2ca02c", ls="--", lw=1.1, label=r"Meta $0{,}5^\circ$C")
    ax3.set_xticks(x)
    ax3.set_xticklabels(models, fontsize=8.5)
    ax3.set_ylabel(r"Erro máximo em $T_j$ ($^\circ$C)")
    ax3.set_title("(c) PINN vs. modelos RC reduzidos")
    ax3.legend(fontsize=7.8, loc="upper left")

    fig.tight_layout()
    out_fig = os.path.join(FIG_DIR, "fig_4_2_pinn_500casos_ood_rc.png")
    fig.savefig(out_fig, dpi=300)
    plt.close(fig)
    print("Saved figure:", out_fig)


if __name__ == "__main__":
    main()
