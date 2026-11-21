"""rodar_campanha_marco5.py -- orquestrador da campanha do Marco 5 (Colab ou local). Chama marco5_tems.py e sintonia_reativo.py como processos.
Uso:   python rodar_campanha_marco5.py --quick --jobs 2 --out m5_quick          (teste de fumaca: 2 sementes, 4 orbitas, 1 sensor)
       python rodar_campanha_marco5.py --jobs 4 --out m5_saidas                  (campanha completa)
       python rodar_campanha_marco5.py --blocos 0 1 --jobs 4 --out m5_saidas      (so alguns blocos; retomavel: parte pronta e pulada)
Saidas (UTF-8) em --out: sintonia_reativo.json, m5_campanha_{nominal,burst,quant,pert,pareto,ood}.json, m5_traces_figuras.json, m5_agregado.json.
NAO produz m5_veredito_P1_P5.json: o texto dos criterios P1-P5 nao esta neste fluxo (ver m5_agregado.json, que traz os numeros de entrada do veredito).
Estado: escrito sem ter sido executado (o sandbox do assistente nao tem os modulos do repositorio). Rode --dry-run e --quick antes da campanha completa.
"""
import argparse
import glob
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

AQUI = os.path.dirname(os.path.abspath(__file__))
PERT = ["r_int_alto", "r_int_muito_alto", "c_mod_menos30", "c_mod_mais30", "eps_baixo", "combinada"]
ENV = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
NOMES = {1: "nominal", 2: "burst", 3: "quant", 4: "pert", 5: "pareto", 6: "ood"}


def chunks(xs, n):
    return [xs[i:i + n] for i in range(0, len(xs), n)]


def monta(cfg):
    """Retorna {bloco: [(arquivo_parte, args_do_marco5_tems)]}."""
    S, S3 = cfg["seeds"], cfg["seeds"][:3]
    base = ["--norb", str(cfg["norb"]), "--ckpt", cfg["ckpt"], "--rtune", cfg["rtune"]]
    B = {k: [] for k in range(1, 8)}

    def add(b, tag, cen, ctrls, sensor, seeds, extra=()):
        for ch in chunks(seeds, cfg["chunk"]):
            B[b].append(("%s_s%d-%d.json" % (tag, ch[0], ch[-1]), ["--scenarios"] + cen + ["--ctrls"] + ctrls + ["--seeds"] + [str(s) for s in ch] + ["--sensor", sensor] + base + list(extra)))
    for sn in cfg["sensores"]:
        add(1, "nom_" + sn, ["nominal"], ["R_SINT", "TEMS_RC3", "TEMS_PINN"], sn, S)
        add(1, "nomdt300_" + sn, ["nominal"], ["TEMS_PINN"], sn, S, ["--dt_ctrl", "300"])    # sem R_SINT: nao ha sintonia para dt = 300 s
        add(1, "nomsolver_" + sn, ["nominal"], ["TEMS_SOLVER"], sn, S3)
        add(2, "burst_" + sn, ["burst"], ["S1", "R_SINT", "TEMS_RC3", "TEMS_PINN"], sn, cfg["seeds_burst"], ["--burst_per_orbit", str(cfg["burst_per_orbit"])])
        for q in cfg["quants"]:
            add(3, "q%s_%s" % (q, sn), ["nominal", "burst"], ["R_SINT", "TEMS_RC3", "TEMS_PINN"], sn, S, ["--quant", q])
        for m in cfg["mt"]:
            add(5, "pareto_m%s_%s" % (m, sn), ["nominal"], ["TEMS_RC3", "TEMS_PINN"], sn, S, ["--mtotal", str(m)])
    add(1, "ref_det", ["nominal"], ["S0", "S1"], cfg["sensores"][-1], [0])              # S0 e S1 nao usam sensores nem sementes: 1 execucao basta (confirmado nos JSONs do Colab)
    sp = cfg["sensores"][-1]
    add(4, "pert", cfg["pert"], ["S1", "R_SINT", "TEMS_RC3", "TEMS_PINN"], sp, S)
    add(6, "ood", ["ood"], ["S1", "R_SINT", "TEMS_RC3", "TEMS_PINN"], sp, S)
    for cen in ("nominal", cfg["pert"][0]):
        B[7].append(("traces_%s.json" % cen, ["--scenarios", cen, "--ctrls", "S0", "S1", "R_SINT", "TEMS_RC3", "TEMS_PINN", "--seeds", "0", "--sensor", sp] + base + ["--save_traces", "@PARTE"]))
    return B


def roda_um(parte, args, cwd):
    if os.path.exists(parte):
        return parte, "pulado", 0.0
    t0 = time.time()
    tmp = parte + ".tmp"
    a = list(args)
    if "@PARTE" in a:                                   # bloco 7: o arquivo principal vai para .res e as series para a parte
        a[a.index("@PARTE")] = tmp
        out = parte + ".res.tmp"
    else:
        out = tmp
    r = subprocess.run([sys.executable, "marco5_tems.py"] + a + ["--out", out], capture_output=True, text=True, env=ENV, cwd=cwd)
    if r.returncode == 0 and os.path.exists(tmp if "@PARTE" in args else out):
        os.replace(tmp, parte)
        return parte, "ok", time.time() - t0
    return parte, "ERRO: " + (r.stderr or r.stdout)[-700:], time.time() - t0


def junta(bloco_dir, destino):
    linhas, contratos, cks, rts = [], set(), set(), set()
    for f in sorted(glob.glob(os.path.join(bloco_dir, "*.json"))):
        d = json.load(open(f, encoding="utf-8"))
        contratos.add(d["contrato_sha256"])
        cks.add(d.get("ckpt_sha256"))
        rts.add(d.get("rtune_sha256"))
        for r in d["resultados"]:
            linhas.append(dict(r, contrato_sha256=d["contrato_sha256"], parte=os.path.basename(f)))
    with open(destino, "w", encoding="utf-8", newline="\n") as fh:
        json.dump({"contratos": sorted(contratos), "ckpt_sha256": sorted(x for x in cks if x), "rtune_sha256": sorted(x for x in rts if x), "n": len(linhas), "resultados": linhas}, fh, ensure_ascii=False)
    return linhas


def agrega(todas):
    g = {}
    for r in todas:
        k = (r["cenario"], r["sensor"], r.get("quant"), r.get("mtotal"), r.get("dt_ctrl"), r["ctrl"])
        g.setdefault(k, []).append(r)
    out = []
    for k, rs in g.items():
        n = len(rs)
        out.append({"cenario": k[0], "sensor": k[1], "quant": k[2], "mtotal": k[3], "dt_ctrl": k[4], "ctrl": k[5], "n": n,
                    "tj_max_last_C_pior": max(x["tj_max_last_C"] for x in rs), "runs_com_violacao": sum(1 for x in rs if x["viol_int_Ks"] > 0),
                    "t_above_80_s_total": sum(x["t_above_80_s"] for x in rs), "tops_last_media": sum(x["tops_frac_last"] for x in rs) / n,
                    "soc_min": min(x["soc_min"] for x in rs), "n_infeasible": sum(x["n_infeasible"] for x in rs), "n_fallback": sum(x["n_fallback"] for x in rs),
                    "frac_pue_periodico_last3": (sum(1 for x in rs if x.get("dE_over_ti_last3") is not None and x["dE_over_ti_last3"] < 1e-3) / n),
                    "pred_err_rms_last_K_media": (lambda v: sum(v) / len(v) if v else None)([x["pred_err_rms_last_K"] for x in rs if x.get("pred_err_rms_last_K") is not None]),
                    "decision_ms_mean_media": sum(x["decision_ms_mean"] for x in rs) / n})
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description="Campanha do Marco 5")
    ap.add_argument("--out", default="m5_saidas")
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--blocos", nargs="+", type=int, default=list(range(0, 8)))
    ap.add_argument("--ckpt", default="v2c.pt")
    ap.add_argument("--norb", type=int, default=None)
    ap.add_argument("--norb-sint", type=int, default=3)
    ap.add_argument("--chunk", type=int, default=5)
    ap.add_argument("--burst-per-orbit", type=int, default=4)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)
    os.chdir(AQUI)
    out = os.path.abspath(a.out)
    os.makedirs(out, exist_ok=True)
    rtune = os.path.join(out, "sintonia_reativo.json")
    cfg = {"seeds": [0, 1, 2, 3, 4], "seeds_burst": list(range(20)), "norb": a.norb or 6, "ckpt": a.ckpt, "rtune": rtune, "chunk": a.chunk, "burst_per_orbit": a.burst_per_orbit,
           "sensores": ["calibrado", "sem_calibracao"], "quants": ["trunc", "round"], "mt": [0, 0.5, 1, 2, 3], "pert": PERT}
    sint = ["--sensors", "calibrado", "sem_calibracao", "--quants", "none", "trunc", "round", "--seeds-cal", "100", "101", "102", "--norb", str(a.norb_sint)]
    if a.quick:
        cfg.update({"seeds": [0, 1], "seeds_burst": [0, 1], "norb": a.norb or 4, "sensores": ["sem_calibracao"], "quants": ["trunc"], "mt": [0, 2], "pert": PERT[:1], "chunk": 2})
        sint = ["--quick", "--sensors", "sem_calibracao", "--quants", "none", "trunc", "--seeds-cal", "100", "--norb", "2"]
    B = monta(cfg)
    if a.dry_run:
        print({NOMES.get(b, "traces"): len(v) for b, v in B.items()}, "| norb", cfg["norb"], "| seeds", cfg["seeds"], "| sensores", cfg["sensores"])
        return B
    from contrato_parametros import PROJETO_1P0, contract_hash, verify_checkpoint
    print("hash do contrato:", contract_hash(PROJETO_1P0), "| checkpoint:", verify_checkpoint(a.ckpt))
    if 0 in a.blocos and not os.path.exists(rtune):
        t0 = time.time()
        r = subprocess.run([sys.executable, "sintonia_reativo.py"] + sint + ["--jobs", str(a.jobs), "--out", rtune], capture_output=True, text=True, env=ENV)
        print(r.stdout[-1500:], r.stderr[-500:] if r.returncode else "", "| bloco 0: %.1f min" % ((time.time() - t0) / 60))
        if r.returncode:
            raise SystemExit("sintonia falhou")
    erros = []
    for b in sorted(x for x in a.blocos if x >= 1):
        d = os.path.join(out, "parts", str(b))
        os.makedirs(d, exist_ok=True)
        t0 = time.time()
        with ThreadPoolExecutor(max_workers=a.jobs) as ex:
            fs = [ex.submit(roda_um, os.path.join(d, p), args, AQUI) for p, args in B[b]]
            for f in as_completed(fs):
                p, st, dt = f.result()
                print("bloco %d %-40s %6.0fs %s" % (b, os.path.basename(p), dt, st[:300]))
                if st.startswith("ERRO"):
                    erros.append((b, p, st))
        print("bloco %d concluido em %.1f min" % (b, (time.time() - t0) / 60))
        if b in NOMES:
            junta(d, os.path.join(out, "m5_campanha_%s.json" % NOMES[b]))
        elif b == 7:
            tr = []
            for f in sorted(glob.glob(os.path.join(d, "traces_*.json"))):
                tr += json.load(open(f, encoding="utf-8"))["traces"]
            with open(os.path.join(out, "m5_traces_figuras.json"), "w", encoding="utf-8", newline="\n") as fh:
                json.dump({"traces": tr}, fh, ensure_ascii=False)
    todas = []
    for n in NOMES.values():
        f = os.path.join(out, "m5_campanha_%s.json" % n)
        if os.path.exists(f):
            todas += json.load(open(f, encoding="utf-8"))["resultados"]
    with open(os.path.join(out, "m5_agregado.json"), "w", encoding="utf-8", newline="\n") as fh:
        json.dump({"agregado": agrega(todas)}, fh, ensure_ascii=False)
    print("blocos com erro:", len(erros))
    for e in erros:
        print(e)
    return 0 if not erros else 1


if __name__ == "__main__":
    sys.exit(main())
