"""gerar_figuras_marcos1_a_4.py -- Gera figuras de publicacao (300 DPI) diretamente dos modulos auditados (Marcos 1 a 4)."""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(DIR, "..", ".."))
FIG_DIR = os.path.join(ROOT, "figuras")
os.makedirs(FIG_DIR, exist_ok=True)
sys.path.insert(0, DIR)

import modelo_orbital as mo
import modelo_ia as mia
import modelo_termico as mt
import modelo_pinn as mp

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "legend.fontsize": 8.5,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": "--",
    "figure.dpi": 300,
    "savefig.bbox": "tight",
})


def fig_1_orbita():
    """Marco 1: Angulo beta e duracao de eclipse ao longo de 2026 para SSO 700 km."""
    days = np.arange(0, 365, 5)
    t0 = mo.datetime(2026, 1, 1, 12, 0, tzinfo=mo.timezone.utc)
    ltans = [(18.0, "LTAN 18:00 (nominal dawn-dusk)", "#1f77b4", "-"),
             (16.5, "LTAN 16:30", "#ff7f0e", "--"),
             (12.0, "LTAN 12:00 (meio-dia)", "#d62728", ":")]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.0))
    for ltan, lbl, col, ls in ltans:
        betas, ecls = [], []
        for d in days:
            dt = t0 + mo.timedelta(days=int(d))
            ds = mo.design_summary(dt, ltan, 700e3, n=4000)
            betas.append(ds["beta_deg"])
            ecls.append(ds["ecl_elipsoide_cone_s"] / 60.0)
        ax1.plot(days, betas, label=lbl, color=col, ls=ls, lw=1.8)
        ax2.plot(days, ecls, label=lbl, color=col, ls=ls, lw=1.8)

    ax1.axhline(64.5, color="gray", ls="-.", lw=1.0, alpha=0.8, label=r"Limite geométrico $\beta^* \approx 64{,}5^\circ$")
    ax1.set_xlabel("Dia do ano (2026)")
    ax1.set_ylabel(r"Ângulo solar $\beta$ ($^\circ$)")
    ax1.set_title("(a) Evolução anual do ângulo solar $\\beta$")
    ax1.legend(loc="lower right", frameon=True)

    ax2.annotate("Solstício jun/2026\n18,42 min (1105,4 s)",
                 xy=(172, 1105.39 / 60.0), xytext=(200, 23.0),
                 arrowprops=dict(arrowstyle="->", color="#1f77b4", lw=1.1),
                 fontsize=8.5, color="#1f77b4", fontweight="bold")
    ax2.set_xlabel("Dia do ano (2026)")
    ax2.set_ylabel("Duração do eclipse por órbita (min)")
    ax2.set_title("(b) Eclipse orbital (elipsoide WGS84 + cone)")
    ax2.legend(loc="upper right", frameon=True)

    fig.tight_layout()
    out = os.path.join(FIG_DIR, "fig_3_1_orbita_eclipse_ltan.png")
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print("Saved:", out)


def fig_2_ia_dvfs():
    """Marco 2d: Potencia por modulo A100 SXM4 vs. frequencia e temperatura de juncao."""
    fs = np.linspace(300.0, 1410.0, 120)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.0))
    for tj, col, ls in [(40.0, "#1f77b4", "-"), (65.0, "#2ca02c", "--"), (80.0, "#ff7f0e", "-."), (95.0, "#d62728", ":")]:
        ptot = [mt.ia_p(float(f), tj)[0] for f in fs]
        ax1.plot(fs, ptot, label=rf"$T_j = {int(tj)}^\circ\mathrm{{C}}$", color=col, ls=ls, lw=1.8)
    ax1.scatter([1410.0], [400.0], color="black", zorder=5, s=35, label="Ponto nominal (1410 MHz, 400 W)")
    ax1.set_xlabel("Frequência de operação $f$ (MHz)")
    ax1.set_ylabel("Potência por módulo A100 (W)")
    ax1.set_title("(a) Potência total por módulo sob DVFS")
    ax1.legend(loc="upper left", frameon=True)

    # Decomposicao em T_j = 80 C
    ia = mt.IA
    fv = ia["fv"] * ia["f_nom"]
    pfix = np.full_like(fs, ia["fix"] * ia["tdp"])
    pdyn, pleak = [], []
    for f in fs:
        v = ia["v_min"] + (ia["v_nom"] - ia["v_min"]) * np.clip((f - fv) / (ia["f_nom"] - fv), 0.0, 1.0)
        pd = (1.0 - ia["fix"] - ia["leak"]) * ia["tdp"] * (v / ia["v_nom"]) ** 2 * (f / ia["f_nom"])
        il = ia["leak"] * ia["tdp"] / ia["v_nom"] * (v / ia["v_nom"]) ** (ia["kv"] - 1.0) * np.exp(ia["beta"] * (80.0 - ia["t_ref"]))
        pdyn.append(pd)
        pleak.append(v * il)
    ax2.stackplot(fs, pfix, pdyn, pleak,
                  labels=[r"$P_{\mathrm{fix}}$ (estática)", r"$P_{\mathrm{dyn}}(V,f)$ (dinâmica)", r"$P_{\mathrm{leak}}(V,T_j)$ (fuga a $80^\circ\mathrm{C}$)"],
                  colors=["#aec7e8", "#1f77b4", "#ff7f0e"], alpha=0.88)
    ax2.set_xlabel("Frequência de operação $f$ (MHz)")
    ax2.set_ylabel("Parcelas de potência por módulo (W)")
    ax2.set_title(r"(b) Decomposição de potência em $T_j = 80^\circ\mathrm{C}$")
    ax2.legend(loc="upper left", frameon=True)

    fig.tight_layout()
    out = os.path.join(FIG_DIR, "fig_3_2_ia_dvfs_fuga.png")
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print("Saved:", out)


def fig_3_planta_termica():
    """Marco 3: Evolucao orbital em regime periodico S0 (sem limite) vs S1 (despacho sustentavel)."""
    T = mt.T_ORB
    t_ecl = mt.T_ECL
    n_orb = 6
    p_sun = mt.analytic_orbit(t_ecl, T, 18000.0)["P_TI_sun_sust_W"]
    des = mt.design_area(p_sun)
    a, n = des["a_pl"], mt.PAR["n_mod"]
    f_sun, f_ecl = mt.ia_f_for_power(p_sun / n, 80.0)[0], mt.ia_f_for_power(18000.0 / n, 80.0)[0]
    hs, he = mt.bus_heat(50000.0, p_sun), mt.bus_heat(0.0, 18000.0)
    qs, qe = hs["heat_aux"] + hs["q_buck"] + hs["q_bd"], he["heat_aux"] + he["q_buck"] + he["q_bd"]
    _, r1 = mt._orbit_run(a, f_sun, f_ecl, qs, qe, n_orb, 25.0)
    h4, h4e = mt.bus_heat(50000.0, 40000.0), mt.bus_heat(0.0, 40000.0)
    _, r0 = mt._orbit_run(a, mt.IA["f_nom"], mt.IA["f_nom"],
                         h4["heat_aux"] + h4["q_buck"] + h4["q_bd"],
                         h4e["heat_aux"] + h4e["q_buck"] + h4e["q_bd"], n_orb, 25.0)

    def last_orb(r):
        m = r["t"] >= (n_orb - 1) * T - 1e-9
        return (r["t"][m] - (n_orb - 1) * T) / 60.0, r["tj"][m] - 273.15, r["t0"][m] - 273.15, r["tN"][m] - 273.15

    t0_m, tj0, tr0, tp0 = last_orb(r0)
    t1_m, tj1, tr1, tp1 = last_orb(r1)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.0), sharey=True)
    for ax, tm, tj, tr, tp, title in [
        (ax1, t0_m, tj0, tr0, tp0, r"(a) Cenário S0 — 1410 MHz contínuo ($T_{j,\max} = %.2f^\circ\mathrm{C}$)" % tj0.max()),
        (ax2, t1_m, tj1, tr1, tp1, r"(b) Cenário S1 — Despacho sustentável ($T_{j,\max} = %.2f^\circ\mathrm{C}$)" % tj1.max()),
    ]:
        ax.axvspan((T - t_ecl) / 60.0, T / 60.0, color="#d9d9d9", alpha=0.55, label="Eclipse (18,42 min)")
        ax.plot(tm, tj, label=r"$T_j$ (junção GPU)", color="#d62728", lw=2.0)
        ax.plot(tm, tr, label=r"$T_{\mathrm{raiz}}$ (coletor)", color="#1f77b4", lw=1.7, ls="--")
        ax.plot(tm, tp, label=r"$T_{\mathrm{ponta}}$ (aleta)", color="#2ca02c", lw=1.5, ls=":")
        ax.axhline(80.0, color="black", ls="-.", lw=1.2, label=r"Limite $T_j = 80^\circ\mathrm{C}$")
        ax.set_xlabel("Tempo na órbita em regime periódico (min)")
        ax.set_title(title)
        ax.legend(loc="lower left", frameon=True, fontsize=8.0)
    ax1.set_ylabel(r"Temperatura ($^\circ\mathrm{C}$)")

    fig.tight_layout()
    out = os.path.join(FIG_DIR, "fig_3_3_planta_termica_s0_s1.png")
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print("Saved:", out)


def fig_4_pinn_v2c():
    """Marco 4: (a) Ablacoes Base -> E1 -> E2 -> E3 -> V2a -> V2b -> V2c; (b) v2c vs Solver em A_sol_quente."""
    labels = ["Base\n(v1)", "E1\n(hard IC)", "E2\n(w_bc=25)", "E3\n(20k ep.)", "V2a\n(K-eq indep)", "V2b\n(K-eq corr)", "V2c\n(+L-BFGS)"]
    # Valores exatos verificados nos JSONs (pior caso e mediana: v1 geral / v2 dominio de treino)
    def load_j(fn):
        return json.load(open(os.path.join(DIR, fn), encoding="utf-8"))

    j_base = load_j("saida_marco4_pinn_eval.json")
    j_e1 = load_j("saida_marco4_e1_eval.json")
    j_e2 = load_j("saida_marco4_e2_eval.json")
    j_e3 = load_j("saida_marco4_e3_v2_eval_reduced.json")
    j_v2a = load_j("saida_marco4_v2a_eval_reduced.json")
    j_v2b = load_j("saida_marco4_v2b_eval_reduced.json")
    j_v2c = load_j("saida_marco4_v2c_eval_reduced.json")

    def med_all(j):
        return float(np.median([v["max_abs_C"] for v in j["casos"].values()]))

    pior = [j_base["pior_caso"]["max_abs_C"],
            j_e1["pior_caso"]["max_abs_C"],
            j_e2["pior_caso"]["max_abs_C"],
            j_e3["resumo"]["no_dominio_treino_v2"]["pior_max_abs_C"],
            j_v2a["resumo"]["no_dominio_treino_v2"]["pior_max_abs_C"],
            j_v2b["resumo"]["no_dominio_treino_v2"]["pior_max_abs_C"],
            j_v2c["resumo"]["no_dominio_treino_v2"]["pior_max_abs_C"]]
    med = [med_all(j_base), med_all(j_e1), med_all(j_e2),
           j_e3["resumo"]["no_dominio_treino_v2"]["mediana_max_abs_C"],
           j_v2a["resumo"]["no_dominio_treino_v2"]["mediana_max_abs_C"],
           j_v2b["resumo"]["no_dominio_treino_v2"]["mediana_max_abs_C"],
           j_v2c["resumo"]["no_dominio_treino_v2"]["mediana_max_abs_C"]]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.2))
    x = np.arange(len(labels))
    w = 0.36
    b1 = ax1.bar(x - w / 2, pior, width=w, label="Erro máximo absoluto (pior caso)", color="#1f77b4", alpha=0.9)
    b2 = ax1.bar(x + w / 2, med, width=w, label="Erro máximo absoluto (mediana)", color="#ff7f0e", alpha=0.9)
    ax1.axhline(0.5, color="#2ca02c", ls="--", lw=1.5, label=r"Meta de projeto ($\leq 0{,}5^\circ\mathrm{C}$)")
    ax1.set_yscale("log")
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, fontsize=8.0)
    ax1.set_ylabel(r"Erro absoluto vs. solver EDP ($^\circ\mathrm{C}$, escala log)")
    ax1.set_title("(a) Evolução das ablações da PINN pura")
    ax1.set_ylim(0.15, 16.0)
    ax1.annotate("0,47 °C\n(10/10 ≤ 0,5 °C)", xy=(6 - w / 2, pior[-1]), xytext=(4.3, 0.22),
                 arrowprops=dict(arrowstyle="->", color="#1f77b4", lw=1.1),
                 fontsize=8.2, fontweight="bold", color="#1f77b4")
    ax1.legend(loc="upper right", frameon=True, fontsize=8.0)

    # Painel (b): Trajetoria A_sol_quente: Solver EDP vs PINN v2c
    case_a = (1353.752, 135.2, 5551.99, 351.0, 321.0, 6.0)
    m = mp.load(os.path.join(DIR, "v2c.pt"))
    t, Y = mp.solver_window(case_a, mp.A_PL, N=160)
    _, tjp, Tp = mp.eval_grid(m, case_a, N=160, nt=len(t))

    ax2.plot(t, Y[:, 0] - 273.15, label=r"$T_j$ — Solver EDP (Marco 3)", color="#d62728", lw=2.0)
    ax2.plot(t, tjp - 273.15, label=r"$T_j$ — PINN pura (v2c.pt)", color="black", ls="--", lw=1.6)
    ax2.plot(t, Y[:, 1] - 273.15, label=r"$T_{\mathrm{raiz}}$ — Solver EDP", color="#1f77b4", lw=2.0)
    ax2.plot(t, Tp[:, 0] - 273.15, label=r"$T_{\mathrm{raiz}}$ — PINN pura (v2c.pt)", color="#17becf", ls="--", lw=1.6)
    ax2.axhline(80.0, color="gray", ls="-.", lw=1.0, label=r"Limite $T_j = 80^\circ\mathrm{C}$")
    ax2.set_xlabel("Tempo na janela de predição $W$ (s)")
    ax2.set_ylabel(r"Temperatura ($^\circ\mathrm{C}$)")
    ax2.set_title(r"(b) Caso crítico A (sol quente): erro máx. $0{,}43^\circ\mathrm{C}$")
    ax2.legend(loc="center right", frameon=True, fontsize=8.0)

    fig.tight_layout()
    out = os.path.join(FIG_DIR, "fig_4_1_pinn_ablacoes_e_validacao.png")
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print("Saved:", out)


if __name__ == "__main__":
    fig_1_orbita()
    fig_2_ia_dvfs()
    fig_3_planta_termica()
    fig_4_pinn_v2c()
