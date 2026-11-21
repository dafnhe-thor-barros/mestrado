"""diag_rint_oraculo.py -- diagnostico de erro de modelo x limite fisico nos cenarios com R_int perturbado.
ESCRITO SEM TER SIDO EXECUTADO. Nao altera nenhum arquivo congelado; so chama marco5_tems.py (CLI ja usada pela campanha) e le os JSONs.

Pergunta: no piloto (R_int alto) o TEMS violou 80 C e o erro de predicao foi de 9 a 17 K (hipotese: o preditor usa o R_int nominal).
Roda, em cada cenario, TEMS_RC3 e TEMS_SOLVER (oraculo), com a margem da regra e com --mtotal 0 (limite superior de vazao).
Como ler (HIPOTESES a confirmar pela saida, nao conclusoes):
  * se no cenario nominal o TEMS_SOLVER tem pred_err pequeno (no piloto: 0,02 K) mas em r_int_* ele continua com pred_err de varios K,
    o solver tambem usa o R_int nominal: este diagnostico NAO separa erro de modelo de limite fisico, e seria preciso passar o R_int real ao
    preditor (mudanca em marco5_tems.py, congelado; nao feita);
  * se em r_int_* o SOLVER tem pred_err pequeno e T_j max <= 80 C enquanto o RC3 viola, o erro de modelo (R_int) explica a violacao;
  * se o SOLVER com --mtotal 0 tambem viola 80 C, ha limite fisico (nem o oraculo segura a temperatura sem sacrificar a vazao).
Uso:  python diag_rint_oraculo.py --jobs 4 --out diag_rint_oraculo
      python diag_rint_oraculo.py --cenarios r_int_alto --seeds 10 --norb 4 --jobs 2 --out diag_teste     (teste curto)
Rode com a maquina livre (nada da campanha em execucao).
"""
import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

AQUI = os.path.dirname(os.path.abspath(__file__))
ENV = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
DIVERG_TXT = "Newton nao convergiu"


def tarefas(a):
    L = []
    for cen in a.cenarios:
        for ct in ("TEMS_RC3", "TEMS_SOLVER"):
            for m in (None, 0):
                L.append((cen, ct, m))
    return L


def nome(a, cen, ct, m):
    return "%s_%s_%s_n%d_s%s.json" % (cen, ct, "regra" if m is None else "m%s" % m, a.norb, "-".join(str(s) for s in a.seeds))


def roda(a, cen, ct, m):
    dest = os.path.join(a.out, nome(a, cen, ct, m))
    if os.path.exists(dest):
        return cen, ct, m, "pulado", dest
    tmp = dest + ".tmp"
    cmd = [sys.executable, "marco5_tems.py", "--scenarios", cen, "--ctrls", ct, "--seeds"] + [str(s) for s in a.seeds] + \
          ["--sensor", a.sensor, "--norb", str(a.norb), "--ckpt", a.ckpt, "--rtune", a.rtune] + (["--mtotal", str(m)] if m is not None else []) + ["--out", tmp]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True, env=ENV, cwd=AQUI)
    err = (r.stderr or "")[-1500:]
    with open(dest + ".err", "w", encoding="utf-8", newline="\n") as fh:
        fh.write(err)
    if r.returncode == 0 and os.path.exists(tmp):
        os.replace(tmp, dest)
        return cen, ct, m, "ok %.0fs" % (time.time() - t0), dest
    if DIVERG_TXT in err:
        return cen, ct, m, "DIVERGIU (Newton nao convergiu)", None
    return cen, ct, m, "ERRO: " + err[-300:].replace("\n", " | "), None


def media(xs):
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None


def f(x, nd=3):
    return "-" if x is None else ("%.*f" % (nd, x))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="diag_rint_oraculo")
    ap.add_argument("--cenarios", nargs="+", default=["nominal", "r_int_alto", "r_int_muito_alto", "combinada"])
    ap.add_argument("--seeds", nargs="+", type=int, default=[10, 11, 12])
    ap.add_argument("--norb", type=int, default=6)
    ap.add_argument("--sensor", default="sem_calibracao")
    ap.add_argument("--ckpt", default="v2c.pt")
    ap.add_argument("--rtune", default="m5_final/sintonia_reativo.json")
    ap.add_argument("--jobs", type=int, default=2)
    a = ap.parse_args(argv)
    os.chdir(AQUI)
    os.makedirs(a.out, exist_ok=True)
    res = []
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        fs = [ex.submit(roda, a, *t) for t in tarefas(a)]
        for fu in as_completed(fs):
            c = fu.result()
            print("%-18s %-12s %-6s %s" % (c[0], c[1], "regra" if c[2] is None else "m%s" % c[2], c[3][:200]), flush=True)
            res.append(c)
    L = ["cenario,ctrl,margem,n,tj_max_pior_C,n_viol_80C,tops_last_medio,pred_err_rms_last_medio_K,n_infeasible_total,status"]
    for cen in a.cenarios:
        for ct in ("TEMS_RC3", "TEMS_SOLVER"):
            for m in (None, 0):
                p = os.path.join(a.out, nome(a, cen, ct, m))
                rot = "regra" if m is None else "m%s" % m
                if not os.path.exists(p):
                    st = [c[3] for c in res if c[0] == cen and c[1] == ct and c[2] == m]
                    L.append("%s,%s,%s,0,-,-,-,-,-,%s" % (cen, ct, rot, (st[0] if st else "nao rodou")[:60].replace(",", ";")))
                    continue
                rs = json.load(open(p, encoding="utf-8"))["resultados"]
                tj = [r.get("tj_max_C") for r in rs if r.get("tj_max_C") is not None]
                L.append("%s,%s,%s,%d,%s,%d,%s,%s,%d,ok" % (cen, ct, rot, len(rs), f(max(tj) if tj else None), sum(1 for x in tj if x > 80.0),
                                                           f(media([r.get("tops_frac_last") for r in rs]), 4), f(media([r.get("pred_err_rms_last_K") for r in rs])),
                                                           int(sum(r.get("n_infeasible") or 0 for r in rs))))
    txt = "\n".join(L)
    with open(os.path.join(a.out, "diag_rint_oraculo_resumo.csv"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(txt)
    print("\n##### DIAGNOSTICO R_int (cole este bloco na conversa) #####")
    print("seeds", a.seeds, "| norb", a.norb, "| sensor", a.sensor, "| maquina: informar na conversa")
    print(txt)
    return 0


if __name__ == "__main__":
    sys.exit(main())
