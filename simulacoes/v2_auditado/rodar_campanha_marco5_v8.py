"""rodar_campanha_marco5.py (v8 = v7 + correcao do --resumo; nenhuma simulacao muda) -- orquestrador da campanha do Marco 5 (Colab ou local).
v8: no --resumo, valores ausentes em quant/mtotal/dt_ctrl podiam virar NaN (pandas novo) e o groupby descartava essas linhas
(HIPOTESE, nao confirmada: AGREGADO e PAREADA so com o bloco pareto). Agora ausentes viram texto ('none', 'nan', '60.0'),
os groupby usam dropna=False e o AGREGADO imprime '# cobertura: X de Y execucoes' para a falha nao passar despercebida.
ESTADO: ESCRITO SEM TER SIDO EXECUTADO. Rode antes: pytest test_marco5.py test_orquestrador_v7.py, depois --dry-run e --quick.
Ao commitar, salve este arquivo com o nome canonico rodar_campanha_marco5.py (sem o sufixo _v7); o test_orquestrador_v7.py importa por esse nome.

Mudancas do v7 em relacao ao v5 commitado:
  (a) --seeds (sementes de TESTE; default 0 1 2 3 4; minimo 3; proibido usar 100 101 102 = calibracao).
  (b) bloco 4 (pert) e bloco 6 (ood) por controlador x cenario; S1 por SEMENTE (1 processo por semente). Todo processo grava <parte>.err.
  (c) divergencia do integrador ("Newton nao convergiu") NAO e mais ERRO: grava <parte>.divergiu e segue. Uma divergencia de S1 conta como "S1 viola" no P2.
      O P2 traz as DUAS leituras (com e sem divergencia); o status principal e o que conta divergencia.
  (d) clausula de falsificacao ESTRITA (R_SINT sem nenhuma violacao) + sensibilidade secundaria com tolerancia de 0,25 K (declarada POS-HOC).
  (e) P5: fracao de execucoes com PUE periodico por grupo (RC3 e PINN) e flag de ressalva.
  (f) PUE periodico usa |dE_over_ti_last3| < LIM_PUE (valor absoluto, como no criterio registrado). O v5 comparava SEM valor absoluto;
      o resumo imprime tambem a fracao sem abs (comparavel com o piloto).
  (g) P4: o nao-alcance do Pareto e separado em acima (R_SINT mais quente que o ponto m=0) e abaixo.
  (h) campanha_config.json: o cache (partes ja existentes) NAO olha o numero de orbitas; o orquestrador recusa retomar com norb/norb-sint/sementes
      diferentes das gravadas. --nova exige pasta sem partes (campanha final). --ignora-config desliga a recusa.
  (i) --latencia registra a maquina (--maquina "texto") e a carga do sistema; rode SOZINHO.
Uso:
  python rodar_campanha_marco5.py --dry-run
  python rodar_campanha_marco5.py --quick --jobs 5 --out m5_quick_srv
  python rodar_campanha_marco5.py --nova --blocos 0 --norb-sint 6 --norb 6 --seeds 10 11 12 13 14 --jobs 5 --chunk 1 --out m5_final
  python rodar_campanha_marco5.py --blocos 1 2 3 4 5 6 7 --norb 6 --seeds 10 11 12 13 14 --jobs 5 --chunk 1 --out m5_final   (mesmos --norb/--seeds/--norb-sint sempre)
  python rodar_campanha_marco5.py --latencia --maquina "i5-8500T, Debian" --out m5_final
  python rodar_campanha_marco5.py --resumo --out m5_final
"""
import argparse
import glob
import hashlib
import json
import math
import os
import platform
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
LIM_PUE = 1e-3          # |dE_bess| / energia de TI na janela de 3 orbitas; APROVADO pelo usuario (exige >= 4 orbitas simuladas)
LIM_TJ = 80.0
TOL_SENS_K = 0.25       # sensibilidade secundaria da clausula de falsificacao (POS-HOC: definida depois de ver o piloto)
TEMS_CTRLS = ["TEMS_RC3", "TEMS_PINN"]
DIVERG_TXT = "Newton nao convergiu"
SEEDS_CAL = (100, 101, 102)


def chunks(xs, n):
    return [xs[i:i + n] for i in range(0, len(xs), n)]


def valida_seeds(seeds, minimo=3):
    seeds = [int(s) for s in seeds]
    if len(set(seeds)) != len(seeds):
        raise SystemExit("--seeds com repeticao")
    if len(seeds) < minimo:
        raise SystemExit("--seeds: minimo 3 sementes (o bootstrap do P4/P5 exige >= 3 pares)")
    if set(seeds) & set(SEEDS_CAL):
        raise SystemExit("--seeds nao pode conter %s (sementes de calibracao do R_SINT)" % (list(SEEDS_CAL),))
    return seeds


def _lista(args, flag):
    """Valores que seguem 'flag' ate o proximo argumento que comeca com '--'."""
    if flag not in args:
        return []
    i = args.index(flag) + 1
    out = []
    while i < len(args) and not str(args[i]).startswith("--"):
        out.append(args[i])
        i += 1
    return out


def monta(cfg):
    """Retorna {bloco: [(arquivo_parte, args_do_marco5_tems, tipo)]}; tipo 'res' (resultados) ou 'trace' (series da ultima orbita)."""
    S = [int(s) for s in cfg["seeds"]]
    S3 = S[:3]
    base = ["--norb", str(cfg["norb"]), "--ckpt", cfg["ckpt"], "--rtune", cfg["rtune"]]
    B = {k: [] for k in range(1, 8)}

    def add(b, tag, cen, ctrls, sensor, seeds, extra=(), chunk=None):
        for ch in chunks(seeds, chunk or cfg["chunk"]):
            B[b].append(("%s_s%d-%d.json" % (tag, ch[0], ch[-1]),
                         ["--scenarios"] + cen + ["--ctrls"] + ctrls + ["--seeds"] + [str(s) for s in ch] + ["--sensor", sensor] + base + list(extra), "res"))
    for sn in cfg["sensores"]:
        add(1, "nom_" + sn, ["nominal"], ["R_SINT", "TEMS_RC3", "TEMS_PINN"], sn, S)
        add(1, "nomdt300_" + sn, ["nominal"], ["TEMS_PINN"], sn, S, ["--dt_ctrl", "300"])     # sem R_SINT: nao ha sintonia para dt = 300 s
        add(1, "nomsolver_" + sn, ["nominal"], ["TEMS_SOLVER"], sn, S3)
        add(2, "burst_" + sn, ["burst"], ["S1", "R_SINT", "TEMS_RC3", "TEMS_PINN"], sn, S, ["--burst_per_orbit", str(cfg["burst_per_orbit"])])
        for q in cfg["quants"]:
            add(3, "q%s_%s" % (q, sn), ["nominal", "burst"], ["R_SINT", "TEMS_RC3", "TEMS_PINN"], sn, S, ["--quant", q, "--burst_per_orbit", str(cfg["burst_per_orbit"])])
        for m in cfg["mt"]:
            add(5, "pareto_m%s_%s" % (m, sn), ["nominal"], ["TEMS_RC3", "TEMS_PINN"], sn, S, ["--mtotal", str(m)])
    sp = cfg["sensores"][-1]
    add(1, "ref_det", ["nominal"], ["S0", "S1"], sp, [S[0]])        # S0 e S1 nao dependem de sensor nem de semente no cenario nominal (confirmado nos JSONs do Colab): 1 execucao basta
    for cen in cfg["pert"]:
        for ct in ("S1", "R_SINT", "TEMS_RC3", "TEMS_PINN"):
            add(4, "pert_%s_%s" % (ct, cen), [cen], [ct], sp, S, chunk=1 if ct == "S1" else None)      # S1 por semente: uma divergencia nao perde as outras sementes
    for ct in ("S1", "R_SINT", "TEMS_RC3", "TEMS_PINN"):
        add(6, "ood_%s" % ct, ["ood"], [ct], sp, S, chunk=1 if ct == "S1" else None)
    for cen in ("nominal", cfg["pert"][0]):
        cts = ["S0", "S1", "R_SINT", "TEMS_RC3", "TEMS_PINN"]
        grupos = [cts] if cen == "nominal" else [[c] for c in cts]      # perturbacao: 1 processo por controlador (S0 diverge em r_int_alto: esperado)
        for g in grupos:
            B[7].append(("traces_%s%s.json" % (cen, "" if len(g) > 1 else "_" + g[0]), ["--scenarios", cen, "--ctrls"] + g + ["--seeds", str(S[0]), "--sensor", sp] + base, "trace"))
    return B


def roda_um(parte, args, tipo, cwd):
    if os.path.exists(parte):
        return parte, "pulado", 0.0
    if os.path.exists(parte + ".divergiu"):
        return parte, "pulado_divergiu", 0.0
    t0 = time.time()
    tmp = parte + ".tmp"
    cmd = list(args) + (["--save_traces", tmp, "--out", parte + ".main.tmp"] if tipo == "trace" else ["--out", tmp])
    r = subprocess.run([sys.executable, "marco5_tems.py"] + cmd, capture_output=True, text=True, env=ENV, cwd=cwd)
    err = r.stderr or ""
    with open(parte + ".err", "w", encoding="utf-8", newline="\n") as fh:
        fh.write(err[-6000:])
    if r.returncode == 0 and os.path.exists(tmp):
        os.replace(tmp, parte)
        if os.path.exists(parte + ".main.tmp"):
            os.remove(parte + ".main.tmp")
        return parte, "ok", time.time() - t0
    if DIVERG_TXT in err or DIVERG_TXT in (r.stdout or ""):
        meta = {"parte": os.path.basename(parte), "motivo": DIVERG_TXT, "cenarios": _lista(args, "--scenarios"), "ctrls": _lista(args, "--ctrls"),
                "seeds": [int(s) for s in _lista(args, "--seeds")], "sensor": (_lista(args, "--sensor") or [None])[0], "args": [str(a) for a in args]}
        with open(parte + ".divergiu", "w", encoding="utf-8", newline="\n") as fh:
            json.dump(meta, fh, ensure_ascii=False)
        return parte, "DIVERGIU", time.time() - t0
    return parte, "ERRO: " + ((r.stderr or r.stdout or "")[-700:]), time.time() - t0


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


def confere_config(out, cfg_id, nova, ignora):
    """O cache pula qualquer parte existente SEM olhar o numero de orbitas. Aqui: --nova exige pasta sem partes; senao recusa config diferente da gravada."""
    p = os.path.join(out, "campanha_config.json")
    if nova:
        sobras = glob.glob(os.path.join(out, "parts", "*", "*")) + glob.glob(os.path.join(out, "sintonia_reativo.json*"))
        if sobras:
            raise SystemExit("--nova: %s nao esta vazia (%d arquivos de partes/sintonia). Use pasta NOVA e VAZIA." % (out, len(sobras)))
    elif os.path.exists(p) and not ignora:
        ant = json.load(open(p, encoding="utf-8"))
        dif = {k: (ant.get(k), v) for k, v in cfg_id.items() if ant.get(k) != v}
        if dif:
            raise SystemExit("config diferente da gravada em %s (gravado, atual): %s. O cache nao olha o numero de orbitas; use pasta nova ou --ignora-config." % (p, dif))
        return
    with open(p, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(cfg_id, fh, ensure_ascii=False, indent=1)


def latencia(out, ckpt, maquina):
    p = os.path.join(out, "latencia_isolada.json")
    try:
        carga = os.getloadavg()[0]
    except Exception:
        carga = None
    if carga is not None and carga > 0.5:
        print("AVISO: carga do sistema (1 min) = %.2f; a latencia so vale com a maquina ociosa (< 0,5). Considere repetir." % carga)
    if not maquina:
        print("AVISO: --maquina nao informada; a latencia so pode ser citada com a maquina.")
    r = subprocess.run([sys.executable, "marco5_tems.py", "--scenarios", "nominal", "--ctrls", "TEMS_RC3", "TEMS_PINN", "--seeds", "0", "--sensor", "sem_calibracao", "--norb", "4",
                        "--ckpt", ckpt, "--out", p], capture_output=True, text=True, env=ENV, cwd=AQUI)
    meta = {"maquina_declarada": maquina, "platform": platform.platform(), "processor": platform.processor(), "cpu_count": os.cpu_count(), "python": sys.version.split()[0],
            "loadavg_1min_antes": carga, "returncode": r.returncode, "quando": time.strftime("%Y-%m-%d %H:%M:%S")}
    with open(os.path.join(out, "latencia_maquina.json"), "w", encoding="utf-8", newline="\n") as fh:
        json.dump(meta, fh, ensure_ascii=False)
    print("latencia isolada: codigo de saida", r.returncode, (r.stderr or "")[-400:])


def _nat(o):
    return o.item() if hasattr(o, "item") else str(o)


def _boot_ic(d, nboot=4000, seed=12345):
    import numpy as np
    d = np.asarray(d, float)
    rng = np.random.default_rng(seed)
    m = [rng.choice(d, len(d)).mean() for _ in range(nboot)]
    return float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def _estado_clausula(r_seguro, p1_ok, val):
    if not val:
        return "INDETERMINADA"
    if r_seguro and p1_ok and all(c in ("EMPATE", "TEMS_MENOR") for c in val):
        return "ACIONADA"
    if r_seguro and p1_ok and any(c in ("EMPATE", "TEMS_MENOR") for c in val):
        return "PARCIAL"
    return "NAO_ACIONADA"


def calcula_veredito(linhas, lat=None, tol_soc=1e-6, diverg=None):
    """Veredito P1-P5 a partir de uma lista de dicts (uma por execucao, com a chave 'bloco'). Funcao pura: sem leitura de arquivos.
    diverg: lista de dicts (bloco, cenarios, ctrls, seeds) das partes que divergiram (arquivos .divergiu).
    Regras (declaradas antes dos dados): P1 = T_j real maximo da execucao <= 80 C em TEMS (nominal, burst, quant; margem da regra, sem override de m_T, dt_ctrl = 60 s);
    P2 = nos cenarios de perturbacao em que o S1 viola (T_j > 80 C OU divergiu), o TEMS nao viola; P3 = SoC >= 0,70 (tol 1e-6) em R_SINT, R, TEMS;
    P4 = vazao do TEMS (interpolada na curva do Pareto) no T_j real maximo da ultima orbita do R_SINT, comparada por semente com IC95 bootstrap;
    P5 = diferenca PINN-RC3 so e atribuivel se |media| > desvio entre sementes do RC3 e IC95 sem o zero (com ressalva de periodicidade por grupo).
    Clausula: estrita (R_SINT sem nenhuma violacao) + sensibilidade secundaria com tolerancia TOL_SENS_K (pos-hoc)."""
    import numpy as np
    import pandas as pd
    df = pd.DataFrame(linhas)
    if df.empty:
        return {"erro": "sem resultados"}
    for c, v in (("quant", "none"), ("mtotal", None), ("dt_ctrl", 60.0), ("n_infeasible", 0), ("infeasible_t_s", None), ("dE_over_ti_last3", None), ("decision_ms_mean", None)):
        if c not in df:
            df[c] = v
    df["mt"] = pd.to_numeric(df["mtotal"], errors="coerce")
    df["dt"] = pd.to_numeric(df["dt_ctrl"], errors="coerce").fillna(60.0)
    df["quant"] = df["quant"].fillna("none").astype(str)
    df["n_infeasible"] = pd.to_numeric(df["n_infeasible"], errors="coerce").fillna(0)
    base = df[(df["dt"] == 60.0) & df["mt"].isna()]
    out = {"n_execucoes": int(len(df))}
    p1 = base[base["ctrl"].isin(TEMS_CTRLS + ["R_SINT"]) & base["bloco"].isin(["nominal", "burst", "quant"])]
    g1 = p1.groupby(["ctrl", "bloco", "cenario", "sensor", "quant"]).agg(n=("seed", "count"), n_viol=("tj_max_C", lambda x: int((x > LIM_TJ).sum())), tj_pior=("tj_max_C", "max")).reset_index()
    g1t, g1r = g1[g1["ctrl"].isin(TEMS_CTRLS)], g1[g1["ctrl"] == "R_SINT"]
    p1_ok = bool(len(g1t) and (g1t["n_viol"] == 0).all())
    r_seguro = bool(len(g1r) and (g1r["n_viol"] == 0).all())
    r_seguro_tol = bool(len(g1r) and (g1r["tj_pior"] <= LIM_TJ + TOL_SENS_K).all())
    out["P1"] = {"status": "INSUFICIENTE" if not len(g1t) else ("ATENDIDO" if p1_ok else "NAO_ATENDIDO"), "grupos_tems": int(len(g1t)), "execucoes_tems": int(g1t["n"].sum()) if len(g1t) else 0,
                 "pior_tj_tems_C": float(g1t["tj_pior"].max()) if len(g1t) else None, "grupos_tems_com_violacao": g1t[g1t["n_viol"] > 0].round(4).to_dict("records"),
                 "R_SINT_sem_violacao": r_seguro, "R_SINT_pior_tj_C": float(g1r["tj_pior"].max()) if len(g1r) else None,
                 "R_SINT_grupos_com_violacao": g1r[g1r["n_viol"] > 0].round(4).to_dict("records")}
    # ---- P2: duas leituras (com e sem divergencia do S1)
    pert = df[(df["bloco"] == "pert") & (df["dt"] == 60.0)]
    s1 = pert[pert["ctrl"] == "S1"].groupby("cenario")["tj_max_C"].max()
    s1_viola_tj = sorted(c for c, v in s1.items() if v > LIM_TJ)
    div_s1 = {}
    for d in (diverg or []):
        if d.get("bloco") == "pert" and "S1" in d.get("ctrls", []):
            for c in d.get("cenarios", []):
                div_s1.setdefault(c, set()).update(int(s) for s in d.get("seeds", []))
    s1_viola_div = sorted(set(s1_viola_tj) | set(div_s1))

    def _p2(cenarios):
        l2 = []
        for c in cenarios:
            for ct in TEMS_CTRLS + ["R_SINT"]:
                x = pert[(pert["cenario"] == c) & (pert["ctrl"] == ct)]
                l2.append({"cenario": c, "ctrl": ct, "n": int(len(x)), "n_viol": int((x["tj_max_C"] > LIM_TJ).sum()), "tj_pior_C": float(x["tj_max_C"].max()) if len(x) else None})
        t_viol = [r for r in l2 if r["ctrl"] in TEMS_CTRLS and r["n_viol"] > 0]
        t_sem = [r for r in l2 if r["ctrl"] in TEMS_CTRLS and r["n"] == 0]
        st = "NAO_APLICAVEL" if not cenarios else ("INSUFICIENTE" if t_sem else ("ATENDIDO" if not t_viol else "NAO_ATENDIDO"))
        return st, l2
    st_b, l2_b = _p2(s1_viola_div)
    st_a, l2_a = _p2(s1_viola_tj)
    out["P2"] = {"status": st_b, "cenarios_em_que_S1_viola": s1_viola_div, "detalhe": l2_b,
                 "S1_viola_por_T_j": s1_viola_tj, "S1_divergiu_por_cenario": {c: sorted(v) for c, v in sorted(div_s1.items())},
                 "leitura_sem_divergencia": {"status": st_a, "cenarios_em_que_S1_viola": s1_viola_tj, "detalhe": l2_a},
                 "nota": "status principal: S1 que divergiu conta como 'S1 viola' (decisao do usuario). Causa fisica da divergencia NAO confirmada (erro numerico 'Newton nao convergiu'). Se NAO_ATENDIDO: reportar honestamente (clausula do pre-registro)"}
    sup = df[df["ctrl"].isin(TEMS_CTRLS + ["R_SINT", "R", "TEMS_SOLVER"])]
    v3 = sup[sup["soc_min"] < 0.70 - tol_soc]
    out["P3"] = {"status": "INSUFICIENTE" if sup.empty else ("ATENDIDO" if v3.empty else "NAO_ATENDIDO"), "execucoes": int(len(sup)),
                 "soc_min_pior": float(sup["soc_min"].min()) if len(sup) else None,
                 "violacoes_amostra": v3[["bloco", "cenario", "ctrl", "sensor", "seed", "soc_min"]].head(50).to_dict("records")}
    # ---- P4
    ref = base[(base["bloco"] == "nominal") & (base["cenario"] == "nominal") & (base["ctrl"] == "R_SINT")]
    par = df[(df["bloco"] == "pareto") & df["ctrl"].isin(TEMS_CTRLS) & (df["dt"] == 60.0)]
    rows4, classes = [], []
    for sn in sorted(set(par["sensor"]) & set(ref["sensor"])):
        for ct in TEMS_CTRLS:
            diffs, f_ac, f_ab = [], 0, 0
            for sd in sorted(set(ref[ref["sensor"] == sn]["seed"])):
                r = ref[(ref["sensor"] == sn) & (ref["seed"] == sd)]
                c = par[(par["sensor"] == sn) & (par["ctrl"] == ct) & (par["seed"] == sd)].sort_values("tj_max_last_C")
                if len(c) < 2:
                    continue
                tj_r, tp_r = float(r["tj_max_last_C"].iloc[0]), float(r["tops_frac_last"].iloc[0])
                x, y = c["tj_max_last_C"].values.astype(float), c["tops_frac_last"].values.astype(float)
                if tj_r < x.min():
                    f_ab += 1
                    continue
                if tj_r > x.max():
                    f_ac += 1
                    continue
                diffs.append(float(np.interp(tj_r, x, y)) - tp_r)
            if len(diffs) >= 3:
                lo, hi = _boot_ic(diffs)
                cl = "TEMS_MAIOR" if lo > 0 else ("TEMS_MENOR" if hi < 0 else "EMPATE")
            else:
                lo = hi = None
                cl = "INSUFICIENTE"
            rows4.append({"sensor": sn, "ctrl": ct, "n_pares": len(diffs), "fora_do_alcance_do_pareto": f_ac + f_ab, "fora_R_SINT_mais_quente_que_m0": f_ac, "fora_R_SINT_mais_frio_que_o_maior_m": f_ab,
                          "dif_media_vazao": float(np.mean(diffs)) if diffs else None, "ic95": [lo, hi], "classe": cl})
            classes.append(cl)
    val = [c for c in classes if c != "INSUFICIENTE"]
    if not val:
        st4 = "INSUFICIENTE"
    elif any(c == "TEMS_MENOR" for c in val):
        st4 = "NAO_ATENDIDO"
    elif all(c == "TEMS_MAIOR" for c in val):
        st4 = "ATENDIDO"
    else:
        st4 = "ATENDIDO_COM_EMPATE"
    out["P4"] = {"status": st4, "detalhe": rows4, "nota": "EMPATE = IC95 bootstrap sobre sementes inclui 0; com poucas sementes o IC e largo (pouco poder)"}
    # ---- clausula de falsificacao: estrita + sensibilidade secundaria (pos-hoc)
    out["clausula_de_falsificacao"] = {"estado": _estado_clausula(r_seguro, p1_ok, val),
                                       "regra": "ESTRITA: R_SINT sem nenhuma violacao de 80 C e TEMS seguro (P1) e vazao sem TEMS_MAIOR => resultado reportado como empate/derrota do TEMS",
                                       "sensibilidade_tol_0p25K_posthoc": {"estado": _estado_clausula(r_seguro_tol, p1_ok, val), "tolerancia_K": TOL_SENS_K,
                                                                           "R_SINT_pior_tj_C": float(g1r["tj_pior"].max()) if len(g1r) else None,
                                                                           "nota": "analise SECUNDARIA; tolerancia definida DEPOIS de ver o piloto (pos-hoc); nao substitui a regra estrita"}}
    # ---- P5 (com periodicidade por grupo)
    pn = base[base["ctrl"].isin(TEMS_CTRLS) & base["bloco"].isin(["nominal", "burst", "quant", "pert", "ood"])]
    rows5 = []

    def _pue_frac(g, ct):
        v = pd.to_numeric(g[g["ctrl"] == ct]["dE_over_ti_last3"], errors="coerce").abs().dropna()
        return float((v < LIM_PUE).mean()) if len(v) else None
    for k, g in pn.groupby(["bloco", "cenario", "sensor", "quant"]):
        p = g.pivot_table(index="seed", columns="ctrl", values="tops_frac_last")
        if "TEMS_PINN" in p.columns and "TEMS_RC3" in p.columns:
            d = (p["TEMS_PINN"] - p["TEMS_RC3"]).dropna().values
            if len(d) >= 3:
                sd_rc3 = float(p["TEMS_RC3"].std(ddof=1))
                lo, hi = _boot_ic(d)
                m = float(np.mean(d))
                pr, pp = _pue_frac(g, "TEMS_RC3"), _pue_frac(g, "TEMS_PINN")
                rows5.append({"bloco": k[0], "cenario": k[1], "sensor": k[2], "quant": k[3], "n": int(len(d)), "dif_media": m, "desvio_RC3_entre_sementes": sd_rc3, "ic95": [lo, hi],
                              "atribuivel_a_PINN": bool(abs(m) > max(sd_rc3, 1e-6) and (lo > 0 or hi < 0)),
                              "pue_ok_RC3": pr, "pue_ok_PINN": pp, "ressalva_periodicidade": bool(any(x is None or x < 1.0 for x in (pr, pp)))})
    n_at = sum(1 for r in rows5 if r["atribuivel_a_PINN"])
    if lat is None:
        lt = df[df["ctrl"].isin(TEMS_CTRLS)].groupby("ctrl")["decision_ms_mean"].median()
        lat_txt = {k: {"mediana_ms_campanha_paralela": float(v)} for k, v in lt.items()}
    else:
        lat_txt = lat
    out["P5"] = {"status": "INSUFICIENTE" if not rows5 else ("ATRIBUIVEL_EM_ALGUNS_GRUPOS" if n_at else "NAO_ATRIBUIVEL"), "grupos_avaliados": len(rows5), "grupos_atribuiveis": n_at,
                 "detalhe": rows5, "latencia_decisao_ms": lat_txt, "regra": "|dif media| > max(desvio entre sementes do RC3, 1e-6) e IC95 sem o zero; senao o argumento da PINN e a latencia. ressalva_periodicidade = algum controlador do grupo com fracao de PUE periodico < 1 (ou sem dado)"}
    inf = df[df["n_infeasible"] > 0]
    inst = []
    for _, r in inf.iterrows():
        ts = r["infeasible_t_s"] if isinstance(r["infeasible_t_s"], list) else []
        for t in ts:
            inst.append({"bloco": r["bloco"], "cenario": r["cenario"], "ctrl": r["ctrl"], "sensor": r["sensor"], "seed": int(r["seed"]), "t_s": float(t)})
    out["limitacoes_n_infeasible"] = {"por_ctrl": {k: int(v) for k, v in df.groupby("ctrl")["n_infeasible"].sum().items()}, "execucoes_com_infeasible": int(len(inf)), "instantes_amostra": inst[:200],
                                      "nota": "decisao inviavel = pico previsto com f_min acima do limite (a LEITURA de T_j ja passa do limite, por ruido/vies do sensor, mesmo com T_j real abaixo; ou R_int real maior que o do modelo); o controlador manda f_min por 1 passo"}
    z = df[df["ctrl"] != "S0"].copy()
    z["dEti_raw"] = pd.to_numeric(z["dE_over_ti_last3"], errors="coerce")
    z = z[z["dEti_raw"].notna()]
    out["PUE_periodico_fracao_por_ctrl(|dE|<%g)" % LIM_PUE] = {k: float(v) for k, v in (z["dEti_raw"].abs() < LIM_PUE).groupby(z["ctrl"]).mean().items()}
    out["PUE_periodico_fracao_por_ctrl_SEM_abs(comparavel_com_piloto)"] = {k: float(v) for k, v in (z["dEti_raw"] < LIM_PUE).groupby(z["ctrl"]).mean().items()}
    out["divergencias_registradas"] = [{k: d.get(k) for k in ("bloco", "parte", "cenarios", "ctrls", "seeds")} for d in (diverg or [])]
    return out


def _texto_chave(x, padrao):
    """Valor de chave de agrupamento como texto; ausente (None/NaN/NA) vira 'padrao'."""
    try:
        if x is None or x != x:
            return padrao
    except Exception:
        return padrao
    return str(x)


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
    div = []
    for b, nome in list(NOMES.items()) + [(7, "traces")]:
        for f in sorted(glob.glob(os.path.join(out, "parts", str(b), "*.divergiu"))):
            d = json.load(open(f, encoding="utf-8"))
            d["bloco"] = nome
            div.append(d)
    df = pd.DataFrame(linhas)
    if not df.empty:
        for c, pad in (("quant", "none"), ("mtotal", "nan"), ("dt_ctrl", "60.0")):
            if c not in df:
                df[c] = None
            df[c] = [_texto_chave(x, pad) for x in df[c]]
        for c in ("pred_err_rms_last_K", "dE_over_ti_last3", "n_infeasible", "n_fallback", "decision_ms_mean", "infeasible_diag"):
            if c not in df:
                df[c] = None
        de = pd.to_numeric(df["dE_over_ti_last3"], errors="coerce")
        df["pue_ok"] = de.abs() < LIM_PUE
        df["pue_ok_semabs"] = de < LIM_PUE
        df["pred_last"] = pd.to_numeric(df["pred_err_rms_last_K"], errors="coerce")

    def cabecalho():
        sha = {f: hashlib.sha256(open(os.path.join(AQUI, f), "rb").read()).hexdigest().upper() for f in ("marco5_tems.py", "sintonia_reativo.py", "rodar_campanha_marco5.py", "contrato_parametros.py", "test_marco5.py", "test_orquestrador_v7.py")
               if os.path.exists(os.path.join(AQUI, f))}
        pc = os.path.join(out, "campanha_config.json")
        conf = json.load(open(pc, encoding="utf-8")) if os.path.exists(pc) else "campanha_config.json ausente"
        return json.dumps({"partes_por_bloco": cont, "execucoes": int(len(df)), "contratos(2 esperados: 60s e 300s)": sorted(hs["contrato"]), "ckpt": sorted(hs["ckpt"]), "rtune": sorted(hs["rtune"]),
                           "sha256_arquivos": sha, "LIM_PUE(|dE|)": LIM_PUE, "config_da_campanha": conf}, ensure_ascii=False, indent=1)

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
        L = ["bloco,partes_executadas,tempo_total_min,maior_parte_min,erros,divergiu"]
        for b, xs in sorted(g.items(), key=lambda kv: str(kv[0])):
            ex = [x for x in xs if not x["status"].startswith("pulado")]
            L.append("%s,%d,%.1f,%.1f,%d,%d" % (b, len(ex), sum(x["s"] for x in ex) / 60.0, (max([x["s"] for x in ex]) if ex else 0.0) / 60.0,
                                                sum(1 for x in xs if x["status"].startswith("ERRO")), sum(1 for x in xs if x["status"].startswith("DIVERGIU"))))
        return "\n".join(L)

    def diverg_txt():
        if not div:
            return "nenhuma parte divergiu"
        return "bloco,parte,cenarios,ctrls,seeds\n" + "\n".join("%s,%s,%s,%s,%s" % (d["bloco"], d.get("parte"), "|".join(d.get("cenarios", [])), "|".join(d.get("ctrls", [])), "|".join(str(s) for s in d.get("seeds", []))) for d in div)

    def agregado():
        if df.empty:
            return "sem resultados"
        g = df.groupby(KEYS + ["ctrl"], dropna=False).agg(n=("seed", "count"), tj_last_pior=("tj_max_last_C", "max"), runs_viol=("viol_int_Ks", lambda x: int((x > 0).sum())),
                                            t80_total=("t_above_80_s", "sum"), tops=("tops_frac_last", "mean"), soc_min=("soc_min", "min"), infeas=("n_infeasible", "sum"),
                                            fallback=("n_fallback", "sum"), frac_pue_ok=("pue_ok", "mean"), frac_pue_semabs=("pue_ok_semabs", "mean"), pred_last=("pred_last", "mean"), lat_ms=("decision_ms_mean", "mean")).reset_index()
        return "# cobertura: %d de %d execucoes\n" % (int(g["n"].sum()), len(df)) + g.round(4).to_csv(index=False)

    def pareado():
        rows = []
        for k, g in df.groupby(KEYS, dropna=False):
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
        pm = os.path.join(out, "latencia_maquina.json")
        mq = json.load(open(pm, encoding="utf-8")) if os.path.exists(pm) else "latencia_maquina.json ausente (latencia sem registro da maquina)"
        return json.dumps({"maquina": mq, "latencia": [{k: r[k] for k in ("ctrl", "decision_ms_mean", "decision_ms_p99")} for r in json.load(open(p, encoding="utf-8"))["resultados"]]}, ensure_ascii=False)

    def traces():
        p = os.path.join(out, "m5_traces_figuras.json")
        if not os.path.exists(p):
            return "m5_traces_figuras.json ausente"
        t = json.load(open(p, encoding="utf-8"))["traces"]
        return "%d series (grandes: nao impressas): %s" % (len(t), sorted({(x["cenario"], x["ctrl"]) for x in t}))

    def veredito():
        lat_d = None
        pl = os.path.join(out, "latencia_isolada.json")
        if os.path.exists(pl):
            lat_d = {r["ctrl"]: {"media_ms": r["decision_ms_mean"], "p99_ms": r["decision_ms_p99"]} for r in json.load(open(pl, encoding="utf-8"))["resultados"]}
        v = calcula_veredito(linhas, lat_d, diverg=div)
        with open(os.path.join(out, "m5_veredito_P1_P5.json"), "w", encoding="utf-8", newline="\n") as fh:
            json.dump(v, fh, ensure_ascii=False, default=_nat)
        c = {}
        for k, x in v.items():
            if isinstance(x, dict):
                c[k] = {kk: (vv[:12] if isinstance(vv, list) else vv) for kk, vv in x.items()}
            else:
                c[k] = x
        return json.dumps(c, ensure_ascii=False, default=_nat, indent=1)

    for nome, fn in (("CABECALHO", cabecalho), ("PYTEST", pytest_txt), ("SINTONIA DO R (sem a grade completa)", sintonia), ("TEMPOS", tempos), ("DIVERGENCIAS (Newton nao convergiu)", diverg_txt),
                     ("AGREGADO (CSV)", agregado), ("DIFERENCA PAREADA DE VAZAO (CSV)", pareado), ("DECISOES INVIAVEIS (diagnostico)", inviaveis), ("LATENCIA ISOLADA", lat),
                     ("SERIES TEMPORAIS", traces), ("VEREDITO P1-P5 (listas cortadas em 12 itens; completo em m5_veredito_P1_P5.json)", veredito)):
        sec(nome, fn)
    txt.append("\nNOTA: veredito P1-P5 calculado por calcula_veredito (regras no docstring da funcao); arquivo m5_veredito_P1_P5.json.")
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
    ap.add_argument("--seeds", nargs="+", type=int, default=None, help="sementes de TESTE (default 0 1 2 3 4; minimo 3; nao pode ter 100 101 102)")
    ap.add_argument("--chunk", type=int, default=5)
    ap.add_argument("--burst-per-orbit", type=int, default=12)
    ap.add_argument("--nova", action="store_true", help="campanha final: exige --out sem partes/sintonia (sem cache)")
    ap.add_argument("--ignora-config", action="store_true")
    ap.add_argument("--maquina", default="", help="texto citado junto com a latencia, ex.: 'i5-8500T, Debian'")
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
        latencia(out, a.ckpt, a.maquina)
        return 0
    rtune = os.path.join(out, "sintonia_reativo.json")
    seeds = valida_seeds(a.seeds if a.seeds else ([0, 1] if a.quick else [0, 1, 2, 3, 4]), minimo=2 if a.quick else 3)
    cfg = {"seeds": seeds, "seeds_burst": list(seeds), "norb": a.norb or 6, "ckpt": a.ckpt, "rtune": rtune, "chunk": a.chunk, "burst_per_orbit": a.burst_per_orbit,
           "sensores": ["calibrado", "sem_calibracao"], "quants": ["trunc", "round"], "mt": [0, 0.5, 1, 2, 3], "pert": PERT}
    sint = ["--sensors", "calibrado", "sem_calibracao", "--quants", "none", "trunc", "round", "--seeds-cal", "100", "101", "102", "--norb", str(a.norb_sint)]
    if a.quick:
        cfg.update({"norb": a.norb or 4, "sensores": ["sem_calibracao"], "quants": ["trunc"], "mt": [0, 2], "pert": PERT[:1], "chunk": 2})
        sint = ["--quick", "--sensors", "sem_calibracao", "--quants", "none", "trunc", "--seeds-cal", "100", "--norb", "2"]
    B = monta(cfg)
    if a.dry_run:
        print({NOMES.get(b, "traces"): len(v) for b, v in B.items()}, "| norb", cfg["norb"], "| norb-sint", a.norb_sint, "| seeds", cfg["seeds"], "| sensores", cfg["sensores"], "| rajadas/orbita", cfg["burst_per_orbit"])
        return 0
    confere_config(out, {"norb": cfg["norb"], "norb_sint": a.norb_sint, "seeds": cfg["seeds"], "burst_per_orbit": cfg["burst_per_orbit"], "quick": bool(a.quick), "ckpt": a.ckpt}, a.nova, a.ignora_config)
    from contrato_parametros import PROJETO_1P0, contract_hash, verify_checkpoint
    print("hash do contrato:", contract_hash(PROJETO_1P0), "| checkpoint:", verify_checkpoint(a.ckpt))
    print("norb", cfg["norb"], "| norb-sint", a.norb_sint, "| seeds", cfg["seeds"], "| chunk", a.chunk, "| jobs", a.jobs)
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
    erros, divs = [], []
    for b in sorted(x for x in a.blocos if x >= 1):
        d = os.path.join(out, "parts", str(b))
        os.makedirs(d, exist_ok=True)
        t0 = time.time()
        with ThreadPoolExecutor(max_workers=a.jobs) as ex:
            fs = [ex.submit(roda_um, os.path.join(d, p), args, tipo, AQUI) for p, args, tipo in B[b]]
            for f in as_completed(fs):
                p, st, dt = f.result()
                print("bloco %d %-48s %6.0fs %s" % (b, os.path.basename(p), dt, st[:300]))
                tempos.append({"bloco": b, "parte": os.path.basename(p), "s": dt, "status": st[:12]})
                if st.startswith("ERRO"):
                    erros.append((b, p, st))
                if st == "DIVERGIU":
                    divs.append((b, os.path.basename(p)))
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
    print("partes com divergencia (Newton nao convergiu):", len(divs), divs[:20])
    print("partes com erro:", len(erros))
    for e in erros:
        print(e)
    return 0 if not erros else 1


if __name__ == "__main__":
    sys.exit(main())
