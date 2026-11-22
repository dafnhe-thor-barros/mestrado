"""teste_convergencia_orbitas.py -- convergencia por orbita (transiente) e recarga do BESS a partir de varias cargas iniciais.
ESCRITO SEM TER SIDO EXECUTADO. Nao edita nenhum arquivo congelado: so importa marco5_tems.run_scenario e passa cfg_over={"soc0": x}.
(marco5_tems.py L503: soc0 = cfg.get("soc0", c["soc_max"]); run_scenario L646-647: cfg.update(cfg_over).)

Por orbita k (1..N) mede: SoC no fim, SoC min/max, vazao (f/f_max), T_j max [C] e |dE|/E_TI (criterio do PUE periodico, 1 orbita).
Convergencia = primeira orbita (contagem 1-based) a partir da qual a variacao em relacao a orbita anterior fica abaixo da tolerancia ATE O FIM.
As tolerancias abaixo sao ESCOLHA DESTE SCRIPT (nao pre-registradas): 0,002 SoC; 0,001 vazao; 0,1 K; 1e-3 para |dE|/E_TI (valor, nao variacao).
Com soc0 < 0,70 e supervisor ativo (floor = 0,70), no eclipse nao ha descarga (bus_step L113-120) e o controlador manda f_min; o modelo entrega TI = 0
mas a planta termica recebe f_eff = f_min (marco5_tems.py L530): possivel artefato, NAO verificado. Interprete essas partidas como comportamento do supervisor.

Uso (dentro de simulacoes/v2_auditado, com o venv do torch ativo):
  python teste_convergencia_orbitas.py --autoteste
  python teste_convergencia_orbitas.py --norb 8 --soc0 0.2 0.95 --seeds 10 --jobs 2 --out conv_curto        (teste curto)
  python teste_convergencia_orbitas.py --norb 20 --soc0 0.0 0.2 0.7 0.95 --seeds 10 --jobs 4 --out conv_n20
  python teste_convergencia_orbitas.py --norb 50 --soc0 0.0 0.95 --seeds 10 --jobs 2 --out conv_n50
Pasta --out NOVA por --norb (o nome do arquivo inclui norb; nao reaproveite pastas de outras campanhas).
"""
import argparse
import json
import os
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)

import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

TOL = {"soc_end": 0.002, "tops": 0.001, "tj_max_C": 0.1, "dE_over_ti": 1e-3}


def por_orbita(ser, T, ecap, fmax, soc0, norb, t0k=273.15):
    """Metricas por orbita a partir da serie completa (dict de listas: t, h, tj [K], soc, f, p_ti)."""
    t, h = np.array(ser["t"], float), np.array(ser["h"], float)
    soc, f, p_ti, tj = np.array(ser["soc"], float), np.array(ser["f"], float), np.array(ser["p_ti"], float), np.array(ser["tj"], float) - t0k
    rows, prev = [], float(soc0)
    for k in range(norb):
        m = (t > k * T + 1e-9) & (t <= (k + 1) * T + 1e-9)
        if not m.any():
            break
        ti_wh = float(np.sum(p_ti[m] * h[m]) / 3600.0)
        end = float(soc[m][-1])
        de = (end - prev) * ecap
        rows.append({"orbita": k + 1, "soc_end": end, "soc_min": float(soc[m].min()), "soc_max": float(soc[m].max()),
                     "tops": float(np.sum(f[m] * h[m]) / np.sum(h[m]) / fmax), "tj_max_C": float(tj[m].max()),
                     "ti_Wh": ti_wh, "dE_Wh": float(de), "dE_over_ti": (abs(de) / ti_wh) if ti_wh > 0 else None})
        prev = end
    return rows


def conv_variacao(vals, tol):
    """Primeira orbita (1-based) cuja variacao |v[k]-v[k-1]| e todas as seguintes ficam < tol. None se nunca (ou se so ha 1 orbita)."""
    v = [x for x in vals]
    if len(v) < 2 or any(x is None for x in v):
        return None
    d = [abs(v[i] - v[i - 1]) for i in range(1, len(v))]
    for k in range(len(d)):
        if all(x < tol for x in d[k:]):
            return k + 2
    return None


def conv_valor(vals, tol):
    """Primeira orbita (1-based) a partir da qual |v| < tol ate o fim."""
    if not vals or any(x is None for x in vals):
        return None
    for k in range(len(vals)):
        if all(abs(x) < tol for x in vals[k:]):
            return k + 1
    return None


def resume(rows):
    return {"orbitas": len(rows),
            "conv_soc_end": conv_variacao([r["soc_end"] for r in rows], TOL["soc_end"]),
            "conv_tops": conv_variacao([r["tops"] for r in rows], TOL["tops"]),
            "conv_tj_max": conv_variacao([r["tj_max_C"] for r in rows], TOL["tj_max_C"]),
            "conv_dE_over_ti": conv_valor([r["dE_over_ti"] for r in rows], TOL["dE_over_ti"]),
            "soc_end_final": rows[-1]["soc_end"] if rows else None,
            "soc_min_global": min(r["soc_min"] for r in rows) if rows else None,
            "soc_max_global": max(r["soc_max"] for r in rows) if rows else None}


def roda(a):
    ctrl, sensor, seed, soc0, norb, ckpt, out = a
    from marco5_tems import run_scenario, T0K
    from modelo_termico import T_ORB
    from contrato_parametros import PROJETO_1P0
    dest = os.path.join(out, "conv_%s_%s_s%d_soc%03d_n%d.json" % (ctrl, sensor, seed, int(round(soc0 * 100)), norb))
    if os.path.exists(dest):
        return dest, "pulado", 0.0
    t0 = time.time()
    try:
        r = run_scenario("nominal", ctrl, seed, n_orb=norb, sensor=sensor, ckpt=ckpt, cfg_over={"soc0": float(soc0)}, series=True)
    except RuntimeError as e:
        res = {"ctrl": ctrl, "sensor": sensor, "seed": seed, "soc0": soc0, "norb": norb, "erro": str(e)[:300]}
        with open(dest, "w", encoding="utf-8", newline="\n") as fh:
            json.dump(res, fh)
        return dest, "ERRO_RUNTIME: " + str(e)[:120], time.time() - t0
    rows = por_orbita(r["series"], T_ORB, r["ecap_Wh"], PROJETO_1P0["f_max"], soc0, norb, T0K)
    res = {"ctrl": ctrl, "sensor": sensor, "seed": seed, "soc0": soc0, "norb": norb, "ecap_Wh": r["ecap_Wh"], "T_orb_s": T_ORB,
           "tolerancias_escolha_do_script": TOL, "resumo": resume(rows), "por_orbita": rows,
           "tj_max_C": r["tj_max_C"], "t_above_80_s": r["t_above_80_s"], "n_infeasible": r["n_infeasible"], "shed_Wh": r["shed_Wh"],
           "soc_min": r["soc_min"], "t_soc_below_min_s": r["t_soc_below_min_s"]}
    with open(dest + ".tmp", "w", encoding="utf-8", newline="\n") as fh:
        json.dump(res, fh)
    os.replace(dest + ".tmp", dest)
    return dest, "ok", time.time() - t0


def autoteste():
    """Teste puro (sem simulacao): serie sintetica com decaimento exponencial do SoC."""
    T, ecap, n = 100.0, 1000.0, 8
    t = np.arange(1.0, n * T + 1.0, 1.0)
    soc_end_k = 0.8 - 0.8 * (0.5 ** np.arange(1, n + 1))
    soc = np.interp(t, np.concatenate(([0.0], np.arange(1, n + 1) * T)), np.concatenate(([0.0], soc_end_k)))
    ser = {"t": list(t), "h": [1.0] * len(t), "tj": [350.0] * len(t), "soc": list(soc), "f": [705.0] * len(t), "p_ti": [1000.0] * len(t)}
    rows = por_orbita(ser, T, ecap, 1410.0, 0.0, n)
    assert len(rows) == n and abs(rows[0]["soc_end"] - 0.4) < 1e-6, rows[0]
    assert abs(rows[0]["tops"] - 0.5) < 1e-9 and abs(rows[0]["tj_max_C"] - 76.85) < 1e-6
    assert conv_variacao([0.4, 0.6, 0.7, 0.75, 0.7501, 0.7501], 0.002) == 5
    assert conv_variacao([0.4, 0.6, 0.7], 0.002) is None
    assert conv_valor([0.1, 0.01, 0.0005, 0.0004], 1e-3) == 3
    assert conv_valor([0.1, 0.0005, 0.01, 0.0004], 1e-3) == 4
    print("autoteste OK")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--norb", type=int, default=20)
    ap.add_argument("--soc0", nargs="+", type=float, default=[0.0, 0.2, 0.7, 0.95])
    ap.add_argument("--seeds", nargs="+", type=int, default=[10])
    ap.add_argument("--sensor", default="sem_calibracao")
    ap.add_argument("--ctrl", default="TEMS_RC3", help="TEMS_RC3 (sem ckpt) ou TEMS_PINN (exige --ckpt)")
    ap.add_argument("--ckpt", default=None)
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--out", default="conv_orbitas")
    ap.add_argument("--autoteste", action="store_true")
    a = ap.parse_args(argv)
    if a.autoteste:
        autoteste()
        return 0
    if any(not (0.0 <= s <= 0.95) for s in a.soc0):
        raise SystemExit("--soc0 deve estar em [0, 0,95] (soc_max = 0,95)")
    if a.norb < 4:
        raise SystemExit("--norb >= 4")
    os.makedirs(a.out, exist_ok=True)
    tarefas = [(a.ctrl, a.sensor, sd, s0, a.norb, a.ckpt, a.out) for sd in a.seeds for s0 in a.soc0]
    saidas = []
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        fs = [ex.submit(roda, t) for t in tarefas]
        for fu in as_completed(fs):
            dest, st, dt = fu.result()
            print("%-60s %6.0fs %s" % (os.path.basename(dest), dt, st), flush=True)
            saidas.append(dest)
    L = ["ctrl,sensor,seed,soc0,norb,conv_soc_end,conv_tops,conv_tj_max,conv_dE_over_ti,soc_end_final,soc_min_global,soc_max_global,n_infeasible,t_above_80_s,shed_Wh"]
    for dest in sorted(saidas):
        if not os.path.exists(dest):
            continue
        d = json.load(open(dest, encoding="utf-8"))
        if "erro" in d:
            L.append("%s,%s,%s,%s,%s,ERRO %s" % (d["ctrl"], d["sensor"], d["seed"], d["soc0"], d["norb"], d["erro"][:60].replace(",", ";")))
            continue
        r = d["resumo"]
        L.append("%s,%s,%d,%.2f,%d,%s,%s,%s,%s,%.4f,%.4f,%.4f,%d,%.1f,%.1f" % (
            d["ctrl"], d["sensor"], d["seed"], d["soc0"], d["norb"], r["conv_soc_end"], r["conv_tops"], r["conv_tj_max"], r["conv_dE_over_ti"],
            r["soc_end_final"], r["soc_min_global"], r["soc_max_global"], d["n_infeasible"], d["t_above_80_s"], d["shed_Wh"]))
    txt = "\n".join(L)
    with open(os.path.join(a.out, "resumo_convergencia.csv"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(txt)
    print("\n##### CONVERGENCIA POR ORBITA (cole este bloco na conversa; None = nao convergiu ate o fim) #####")
    print(txt)
    print("\n##### SoC no fim de cada orbita (primeira semente, cada soc0) #####")
    for dest in sorted(saidas):
        if os.path.exists(dest):
            d = json.load(open(dest, encoding="utf-8"))
            if "por_orbita" in d and d["seed"] == a.seeds[0]:
                print("soc0=%.2f: %s" % (d["soc0"], " ".join("%.3f" % r["soc_end"] for r in d["por_orbita"])))
    return 0


if __name__ == "__main__":
    sys.exit(main())
