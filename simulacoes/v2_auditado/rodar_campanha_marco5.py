"""rodar_campanha_marco5.py -- orquestrador da campanha do Marco 5 (Colab ou local). Chama marco5_tems.py e sintonia_reativo.py como processos.
Uso:   python rodar_campanha_marco5.py --dry-run                         (so lista quantas execucoes cada bloco tem)
       python rodar_campanha_marco5.py --quick --jobs 2 --out m5_quick    (teste de fumaca: 2 sementes, 4 orbitas, 1 sensor)
       python rodar_campanha_marco5.py --blocos 0 --jobs 4 --out m5_saidas
       python rodar_campanha_marco5.py --blocos 1 2 --jobs 4 --out m5_saidas   (retomavel: partes prontas sao puladas)
       python rodar_campanha_marco5.py --latencia --out m5_saidas         (rodar SOZINHO, sem outros processos)
       python rodar_campanha_marco5.py --resumo --out m5_saidas           (imprime e grava m5_resumo_para_colar.txt)
Saidas (UTF-8) em --out: sintonia_reativo.json, m5_campanha_{nominal,burst,quant,pert,pareto,ood}.json, m5_traces_figuras.json, tempos.json, m5_resumo_para_colar.txt.
NAO produz o veredito P1-P5: o texto dos criterios nao foi fornecido. O resumo traz os numeros de entrada dele.
Estado: escrito sem ter sido executado. Rode --dry-run e --quick antes da campanha completa.
"""
import argparse
import glob
import hashlib
import json
import math
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

AQUI = os.path.dirname(os.path.abspath(__file__))
PERT = ["r_int_alto", "r_int_muito_alto", "c_mod_menos30", "c_mod_mais30", "eps_baixo", "combinada"]
ENV = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
NOMES = {1: "nominal", 2: "burst", 3: "quant", 4: "pert", 5: "pareto", 6: "ood"}
PARES = [("TEMS_PINN", "R_SINT"), ("TEMS_RC3", "R_SINT"), ("TEMS_PINN", "TEMS_RC3"), ("TEMS_SOLVER", "TEMS_RC3"), ("R_SINT", "S1")]
KEYS = ["bloco", "cenario", "sensor", "quant", "mtotal", "dt_ctrl"]
LIM_PUE = 1e-3     # [HIPOTESE DECLARADA, pendente de aprovacao] |dE_bess| / energia de TI na janela de 3 orbitas


def chunks(xs, n):
    return [xs[i:i + n] for i in range(0, len(xs), n)]


def monta(cfg):
    """Retorna {bloco: [(arquivo_parte, args_do_marco5_tems, tipo)]}; tipo 'res' (resultados) ou 'trace' (series da ultima orbita)."""
    S, S3 = cfg["seeds"], cfg["seeds"][:3]
    base = ["--norb", str(cfg["norb"]), "--ckpt", cfg["ckpt"], "--rtune", cfg["rtune"]]
    B = {k: [] for k in range(1, 8)}

    def add(b, tag, cen, ctrls, sensor, seeds, extra=()):
        for ch in chunks(seeds, cfg["chunk"]):
            B[b].append(("%s_s%d-%d.json" % (tag, ch[0], ch[-1]),
                         ["--scenarios"] + cen + ["--ctrls"] + ctrls + ["--seeds"] + [str(s) for s in ch] + ["--sensor", sensor] + base + list(extra), "res"))
    for sn in cfg["sensores"]:
        add(1, "nom_" + sn, ["nominal"], ["R_SINT", "TEMS_RC3", "TEMS_PINN"], sn, S)
        add(1, "nomdt300_" + sn, ["nominal"], ["TEMS_PINN"], sn, S, ["--dt_ctrl", "300"])     # sem R_SINT: nao ha sintonia para dt = 300 s
        add(1, "nomsolver_" + sn, ["nominal"], ["TEMS_SOLVER"], sn, S3)
        add(2, "burst_" + sn, ["burst"], ["S1", "R_SINT", "TEMS_RC3", "TEMS_PINN"], sn, cfg["seeds_burst"], ["--burst_per_orbit", str(cfg["burst_per_orbit"])])
        for q in cfg["quants"]:
            add(3, "q%s_%s" % (q, sn), ["nominal", "burst"], ["R_SINT", "TEMS_RC3", "TEMS_PINN"], sn, S, ["--quant", q, "--burst_per_orbit", str(cfg["burst_per_orbit"])])
        for m in cfg["mt"]:
            add(5, "pareto_m%s_%s" % (m, sn), ["nominal"], ["TEMS_RC3", "TEMS_PINN"], sn, S, ["--mtotal", str(m)])
    sp = cfg["sensores"][-1]
    add(1, "ref_det", ["nominal"], ["S0", "S1"], sp, [0])        # S0 e S1 nao dependem de sensor nem de semente (confirmado nos JSONs do Colab): 1 execucao basta
    add(4, "pert", cfg["pert"], ["S1", "R_SINT", "TEMS_RC3", "TEMS_PINN"], sp, S)
    add(6, "ood", ["ood"], ["S1", "R_SINT", "TEMS_RC3", "TEMS_PINN"], sp, S)
    for cen in ("nominal", cfg["pert"][0]):
        B[7].append(("traces_%s.json" % cen, ["--scenarios", cen, "--ctrls", "S0", "S1", "R_SINT", "TEMS_RC3", "TEMS_PINN", "--seeds", "0", "--sensor", sp] + base, "trace"))
    return B


def roda_um(parte, args, tipo, cwd):
    if os.path.exists(parte):
        return parte, "pulado", 0.0
    t0 = time.time()
    tmp = parte + ".tmp"
    cmd = list(args) + (["--save_traces", tmp, "--out", parte + ".main.tmp"] if tipo == "trace" else ["--out", tmp])
    r = subprocess.run([sys.executable, "marco5_tems.py"] + cmd, capture_output=True, text=True, env=ENV, cwd=cwd)
    if r.returncode == 0 and os.path.exists(tmp):
        os.replace(tmp, parte)
        if os.path.exists(parte + ".main.tmp"):
            os.remove(parte + ".main.tmp")
        return parte, "ok", time.time() - t0
    return parte, "ERRO: " + ((r.stderr or r.stdout)[-700:]), time.time() - t0


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


def salva_tempos(out, novos):
    p = os.path.join(out, "tempos.json")
    ant = json.load(open(p, encoding="utf-8")) if os.path.exists(p) else []
    with open(p, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(ant + novos, fh, ensure_ascii=False)


def latencia(out, ckpt):
    p = os.path.join(out, "latencia_isolada.json")
    r = subprocess.run([sys.executable, "marco5_tems.py", "--scenarios", "nominal", "--ctrls", "TEMS_RC3", "TEMS_PINN", "--seeds", "0", "--sensor", "sem_calibracao", "--norb", "4",
                        "--ckpt", ckpt, "--out", p], capture_output=True, text=True, env=ENV, cwd=AQUI)
    print("latencia isolada: codigo de saida", r.returncode, (r.stderr or "")[-400:])


def gerar_resumo(out, imprimir=True):
    """Consolida tudo em UM TXT (m5_resumo_para_colar.txt) para voce colar na conversa."""
    import pandas as pd
    txt = []

    def sec(titulo, fn):
        txt.append("\n##### %s #####" % titulo)
        try:
            r = fn()
            txt.append(r if isinstance(r, str) else str(r))
        except Exception as e:
            txt.append("ERRO nesta secao: %s: %s" % (type(e).__name__, str(e)[:300]))

    linhas, hs, cont = [], {"contrato": set(), "ckpt": set(), "rtune": set()}, {}
    for b, nome in NOMES.items():
        arqs = sorted(glob.glob(os.path.join(out, "parts", str(b), "*.json")))
        cont[nome] = len(arqs)
        for f in arqs:
            d = json.load(open(f, encoding="utf-8"))
            for k, kk in (("contrato", "contrato_sha256"), ("ckpt", "ckpt_sha256"), ("rtune", "rtune_sha256")):
                if d.get(kk):
                    hs[k].add(d[kk])
            for r in d["resultados"]:
                r = dict(r)
                r["bloco"] = nome
                linhas.append(r)
    df = pd.DataFrame(linhas)
    if not df.empty:
        for c in ("quant", "mtotal", "dt_ctrl"):
            if c not in df:
                df[c] = None
            df[c] = df[c].astype(str)
        for c in ("pred_err_rms_last_K", "dE_over_ti_last3", "n_infeasible", "n_fallback", "decision_ms_mean", "infeasible_diag"):
            if c not in df:
                df[c] = None
        df["pue_ok"] = pd.to_numeric(df["dE_over_ti_last3"], errors="coerce") < LIM_PUE
        df["pred_last"] = pd.to_numeric(df["pred_err_rms_last_K"], errors="coerce")

    def cabecalho():
        sha = {f: hashlib.sha256(open(os.path.join(AQUI, f), "rb").read()).hexdigest().upper() for f in ("marco5_tems.py", "sintonia_reativo.py", "rodar_campanha_marco5.py", "contrato_parametros.py", "test_marco5.py")
               if os.path.exists(os.path.join(AQUI, f))}
        return json.dumps({"partes_por_bloco": cont, "execucoes": int(len(df)), "contratos(2 esperados: 60s e 300s)": sorted(hs["contrato"]), "ckpt": sorted(hs["ckpt"]), "rtune": sorted(hs["rtune"]),
                           "sha256_arquivos": sha, "LIM_PUE": LIM_PUE}, ensure_ascii=False, indent=1)

    def pytest_txt():
        p = os.path.join(out, "pytest.txt")
        if not os.path.exists(p):
            return "pytest.txt ausente (rode a celula de testes com tee, ou copie a saida do pytest)"
        return "".join(open(p, encoding="utf-8", errors="replace").readlines()[-25:])

    def sintonia():
        p = os.path.join(out, "sintonia_reativo.json")
        if not os.path.exists(p):
            return "sintonia_reativo.json ausente"
        d = json.load(open(p, encoding="utf-8"))
        o = {"seeds_cal": d["seeds_cal"], "norb": d["norb"], "grade": d["grade"],
             "sensores": {k: {x: v for x, v in e.items() if x != "grade_completa"} for k, e in d["sensores"].items()}}
        return json.dumps(o, ensure_ascii=False)

    def tempos():
        p = os.path.join(out, "tempos.json")
        if not os.path.exists(p):
            return "tempos.json ausente"
        g = {}
        for x in json.load(open(p, encoding="utf-8")):
            g.setdefault(x["bloco"], []).append(x)
        L = ["bloco,partes_executadas,tempo_total_min,maior_parte_min,erros"]
        for b, xs in sorted(g.items(), key=lambda kv: kv[0]):
            ex = [x for x in xs if x["status"] != "pulado"]
            L.append("%s,%d,%.1f,%.1f,%d" % (b, len(ex), sum(x["s"] for x in ex) / 60.0, (max([x["s"] for x in ex]) if ex else 0.0) / 60.0, sum(1 for x in xs if x["status"].startswith("ERRO"))))
        return "\n".join(L)

    def agregado():
        if df.empty:
            return "sem resultados"
        g = df.groupby(KEYS + ["ctrl"]).agg(n=("seed", "count"), tj_last_pior=("tj_max_last_C", "max"), runs_viol=("viol_int_Ks", lambda x: int((x > 0).sum())),
                                            t80_total=("t_above_80_s", "sum"), tops=("tops_frac_last", "mean"), soc_min=("soc_min", "min"), infeas=("n_infeasible", "sum"),
                                            fallback=("n_fallback", "sum"), frac_pue_ok=("pue_ok", "mean"), pred_last=("pred_last", "mean"), lat_ms=("decision_ms_mean", "mean")).reset_index()
        return g.round(4).to_csv(index=False)

    def pareado():
        rows = []
        for k, g in df.groupby(KEYS):
            p = g.pivot_table(index="seed", columns="ctrl", values="tops_frac_last")
            for a, b in PARES:
                if a in p and b in p:
                    x = (p[a] - p[b]).dropna()
                    if len(x) >= 2:
                        rows.append(dict(zip(KEYS, k), a=a, b=b, n=len(x), dif=x.mean(), ep=x.std(ddof=1) / math.sqrt(len(x)), frac_a_maior=(x > 1e-6).mean(), frac_iguais=(x.abs() <= 1e-6).mean()))
        return pd.DataFrame(rows).round(5).to_csv(index=False) if rows else "sem pares"

    def inviaveis():
        L = []
        for _, r in df.iterrows():
            dg = r["infeasible_diag"]
            if isinstance(dg, list) and dg:
                for x in dg[:5]:
                    L.append(json.dumps(dict(cenario=r["cenario"], ctrl=r["ctrl"], sensor=r["sensor"], seed=int(r["seed"]), **x), ensure_ascii=False))
        return "\n".join(L[:60]) or "nenhuma decisao inviavel registrada (ou campo ausente)"

    def lat():
        p = os.path.join(out, "latencia_isolada.json")
        if not os.path.exists(p):
            return "latencia_isolada.json ausente (rode --latencia sozinho)"
        return json.dumps([{k: r[k] for k in ("ctrl", "decision_ms_mean", "decision_ms_p99")} for r in json.load(open(p, encoding="utf-8"))["resultados"]])

    def traces():
        p = os.path.join(out, "m5_traces_figuras.json")
        if not os.path.exists(p):
            return "m5_traces_figuras.json ausente"
        t = json.load(open(p, encoding="utf-8"))["traces"]
        return "%d series (grandes: nao impressas): %s" % (len(t), sorted({(x["cenario"], x["ctrl"]) for x in t}))

    for nome, fn in (("CABECALHO", cabecalho), ("PYTEST", pytest_txt), ("SINTONIA DO R (sem a grade completa)", sintonia), ("TEMPOS", tempos), ("AGREGADO (CSV)", agregado),
                     ("DIFERENCA PAREADA DE VAZAO (CSV)", pareado), ("DECISOES INVIAVEIS (diagnostico)", inviaveis), ("LATENCIA ISOLADA", lat), ("SERIES TEMPORAIS", traces)):
        sec(nome, fn)
    txt.append("\nNOTA: veredito P1-P5 NAO calculado (o texto dos criterios nao foi fornecido). Os numeros acima sao a entrada dele.")
    s = "\n".join(txt)
    p = os.path.join(out, "m5_resumo_para_colar.txt")
    with open(p, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(s)
    if imprimir:
        print("##### INICIO #####\n" + s + "\n##### FIM (%d caracteres; tambem gravado em %s) #####" % (len(s), p))
    return s


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
    ap.add_argument("--burst-per-orbit", type=int, default=12)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--resumo", action="store_true")
    ap.add_argument("--latencia", action="store_true")
    a = ap.parse_args(argv)
    os.chdir(AQUI)
    out = os.path.abspath(a.out)
    os.makedirs(out, exist_ok=True)
    if a.resumo:
        gerar_resumo(out, True)
        return 0
    if a.latencia:
        latencia(out, a.ckpt)
        return 0
    rtune = os.path.join(out, "sintonia_reativo.json")
    cfg = {"seeds": [0, 1, 2, 3, 4], "seeds_burst": list(range(20)), "norb": a.norb or 6, "ckpt": a.ckpt, "rtune": rtune, "chunk": a.chunk, "burst_per_orbit": a.burst_per_orbit,
           "sensores": ["calibrado", "sem_calibracao"], "quants": ["trunc", "round"], "mt": [0, 0.5, 1, 2, 3], "pert": PERT}
    sint = ["--sensors", "calibrado", "sem_calibracao", "--quants", "none", "trunc", "round", "--seeds-cal", "100", "101", "102", "--norb", str(a.norb_sint)]
    if a.quick:
        cfg.update({"seeds": [0, 1], "seeds_burst": [0, 1], "norb": a.norb or 4, "sensores": ["sem_calibracao"], "quants": ["trunc"], "mt": [0, 2], "pert": PERT[:1], "chunk": 2})
        sint = ["--quick", "--sensors", "sem_calibracao", "--quants", "none", "trunc", "--seeds-cal", "100", "--norb", "2"]
    B = monta(cfg)
    if a.dry_run:
        print({NOMES.get(b, "traces"): len(v) for b, v in B.items()}, "| norb", cfg["norb"], "| seeds", cfg["seeds"], "| sensores", cfg["sensores"], "| rajadas/orbita", cfg["burst_per_orbit"])
        return 0
    from contrato_parametros import PROJETO_1P0, contract_hash, verify_checkpoint
    print("hash do contrato:", contract_hash(PROJETO_1P0), "| checkpoint:", verify_checkpoint(a.ckpt))
    tempos = []
    if 0 in a.blocos:
        if os.path.exists(rtune):
            print("bloco 0: sintonia_reativo.json ja existe (apague para refazer)")
        else:
            t0 = time.time()
            r = subprocess.run([sys.executable, "sintonia_reativo.py"] + sint + ["--jobs", str(a.jobs), "--out", rtune], capture_output=True, text=True, env=ENV)
            print(r.stdout[-1800:], r.stderr[-600:] if r.returncode else "", "| bloco 0: %.1f min" % ((time.time() - t0) / 60))
            tempos.append({"bloco": 0, "parte": "sintonia", "s": time.time() - t0, "status": "ok" if r.returncode == 0 else "ERRO"})
            if r.returncode:
                salva_tempos(out, tempos)
                raise SystemExit("sintonia falhou")
    if any(b >= 1 for b in a.blocos) and not os.path.exists(rtune):
        raise SystemExit("falta %s: rode antes --blocos 0" % rtune)
    erros = []
    for b in sorted(x for x in a.blocos if x >= 1):
        d = os.path.join(out, "parts", str(b))
        os.makedirs(d, exist_ok=True)
        t0 = time.time()
        with ThreadPoolExecutor(max_workers=a.jobs) as ex:
            fs = [ex.submit(roda_um, os.path.join(d, p), args, tipo, AQUI) for p, args, tipo in B[b]]
            for f in as_completed(fs):
                p, st, dt = f.result()
                print("bloco %d %-40s %6.0fs %s" % (b, os.path.basename(p), dt, st[:300]))
                tempos.append({"bloco": b, "parte": os.path.basename(p), "s": dt, "status": st[:12]})
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
    salva_tempos(out, tempos)
    try:
        gerar_resumo(out, False)
        print("resumo atualizado em", os.path.join(out, "m5_resumo_para_colar.txt"), "(imprima com --resumo)")
    except Exception as e:
        print("resumo nao gerado:", type(e).__name__, str(e)[:200])
    print("partes com erro:", len(erros))
    for e in erros:
        print(e)
    return 0 if not erros else 1


if __name__ == "__main__":
    sys.exit(main())
