"""analise_fase_janela.py -- le os JSONs gerados por teste_convergencia_orbitas.py (SO LEITURA: nao simula, nao importa marco5_tems).
ESCRITO SEM TER SIDO EXECUTADO.

Para cada orbita n mostra:
  - fase do INICIO e do FIM do eclipse em relacao a grade de decisao (dt_ctrl, padrao 60 s);
  - SoC no fim da orbita;
  - |dE|/E_TI em 1 orbita e na janela de 3 orbitas (mesmo criterio de dE_over_ti_last3, marco5_tems.summarize L559-571) e se passa em 1e-3.
No fim: fracao das orbitas finais n >= 4 em que a janela de 3 orbitas passaria (isto e, o que o criterio de PUE periodico diria se a
simulacao terminasse na orbita n), e a correlacao de Spearman entre a fase do fim do eclipse e o SoC no fim da orbita (orbitas >= nmin).
Premissa NAO verificada aqui: eclipse centrado em T/2 (docstring de modelo_termico.py L175), logo inicio_k = T/2 - T_ECL/2 + k*T.
T_ORB e T_ECL vem de modelo_termico.py L21 (5926.379071134441 e 1105.3908).
Uso:
  python analise_fase_janela.py conv_curto            (pasta com conv_*_n20.json; use --dt se dt_ctrl != 60)
  python analise_fase_janela.py --autoteste
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

T_ORB, T_ECL = 5926.379071134441, 1105.3908
LIM = 1e-3


def fases(n, dt, T=T_ORB, te=T_ECL):
    """Lista de (fase_inicio_eclipse, fase_fim_eclipse) em segundos, modulo dt, para as orbitas k = 0..n-1."""
    ini0 = 0.5 * T - 0.5 * te
    out = []
    for k in range(n):
        ini = ini0 + k * T
        out.append((ini % dt, (ini + te) % dt))
    return out


def _rank(x):
    x = np.asarray(x, float)
    order = np.argsort(x, kind="mergesort")
    r = np.empty(len(x))
    r[order] = np.arange(len(x), dtype=float)
    for v in np.unique(x):
        m = x == v
        if m.sum() > 1:
            r[m] = r[m].mean()
    return r


def spearman(a, b):
    ra, rb = _rank(a), _rank(b)
    if len(ra) < 3 or ra.std() == 0 or rb.std() == 0:
        return None
    return float(np.corrcoef(ra, rb)[0, 1])


def analisa(path, dt, nmin):
    d = json.load(open(path, encoding="utf-8"))
    if "por_orbita" not in d:
        print("%s: sem por_orbita (%s)" % (os.path.basename(path), d.get("erro", "?")))
        return
    rows, ecap, s0 = d["por_orbita"], d["ecap_Wh"], d["soc0"]
    n_tot = len(rows)
    ph = fases(n_tot, dt)
    soc = [r["soc_end"] for r in rows]
    ti = [r["ti_Wh"] for r in rows]
    print("\n== %s | soc0=%.2f norb=%d ecap=%.1f Wh ==" % (os.path.basename(path), s0, n_tot, ecap))
    print("orb  ini_ecl(s)  fim_ecl(s)  soc_end   rel1     rel3     ok1  ok3")
    ok3_list, fim_ph, soc_sel = [], [], []
    for n in range(1, n_tot + 1):
        prev = soc[n - 2] if n >= 2 else s0
        r1 = abs(soc[n - 1] - prev) * ecap / ti[n - 1] if ti[n - 1] > 0 else None
        r3 = None
        if n >= 4:
            tt = sum(ti[n - 3:n])
            r3 = abs(soc[n - 1] - soc[n - 4]) * ecap / tt if tt > 0 else None
            if r3 is not None:
                ok3_list.append((n, r3 < LIM))
        if n >= nmin:
            fim_ph.append(ph[n - 1][1])
            soc_sel.append(soc[n - 1])
        print("%3d  %9.2f  %9.2f  %.4f  %-8s %-8s %-4s %-4s" % (
            n, ph[n - 1][0], ph[n - 1][1], soc[n - 1],
            "-" if r1 is None else "%.2e" % r1, "-" if r3 is None else "%.2e" % r3,
            "-" if r1 is None else ("sim" if r1 < LIM else "nao"), "-" if r3 is None else ("sim" if r3 < LIM else "nao")))
    if ok3_list:
        ok = [n for n, v in ok3_list if v]
        print("janela de 3 orbitas (n >= 4): passa em %d de %d finais possiveis; nao passa em n = %s" % (
            len(ok), len(ok3_list), [n for n, v in ok3_list if not v]))
    rho = spearman(fim_ph, soc_sel)
    print("Spearman(fase do fim do eclipse, SoC no fim da orbita), orbitas >= %d: %s" % (nmin, "n/d" if rho is None else "%.3f" % rho))


def autoteste():
    f = fases(3, 60.0)
    assert abs(f[0][0] - 10.4941) < 1e-3, f[0]
    assert abs(f[1][0] - 56.8731) < 1e-3, f[1]
    assert abs(f[0][1] - ((2410.4941355672 + T_ECL) % 60.0)) < 1e-6
    assert spearman([1, 2, 3, 4], [4, 3, 2, 1]) == -1.0
    assert abs(spearman([1, 2, 3, 4], [1, 2, 3, 4]) - 1.0) < 1e-12
    print("autoteste OK")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("pasta", nargs="?", default=None)
    ap.add_argument("--dt", type=float, default=60.0)
    ap.add_argument("--nmin", type=int, default=1, help="primeira orbita usada na correlacao de Spearman")
    ap.add_argument("--autoteste", action="store_true")
    a = ap.parse_args(argv)
    if a.autoteste:
        autoteste()
        return 0
    if not a.pasta:
        raise SystemExit("informe a pasta com os conv_*.json")
    arqs = sorted(glob.glob(os.path.join(a.pasta, "conv_*.json")))
    if not arqs:
        raise SystemExit("nenhum conv_*.json em %s" % a.pasta)
    for p in arqs:
        analisa(p, a.dt, a.nmin)
    return 0


if __name__ == "__main__":
    sys.exit(main())
