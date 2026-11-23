"""junta_campanhas.py -- junta campanhas longas COMPATIVEIS (mesmo horizonte, mesma sintonia, mesmo checkpoint, sementes DISJUNTAS) e recalcula o veredito P1-P5
com a MESMA funcao do orquestrador v8 (rodar_campanha_marco5.calcula_veredito), mais a janela de PUE periodico de k orbitas (janela_pue.py).
SO LEITURA das pastas de entrada: nao simula, nao edita nada congelado. AUDITORIA 22/11/2026: codigo de saida != 0 se o veredito falhar; o hash do contrato agora inclui a planta (c_area), entao campanhas de plantas diferentes nao juntam.

Entradas: pastas de campanha (cada uma com campanha_config.json e parts/[1-6]/*.json OU m5_campanha_*.json). Prefere parts/ (fonte primaria) como o v8.
Verificacoes ANTES de juntar (erro => nao junta, exceto com --forca, que vira aviso; sementes repetidas SEMPRE sao erro):
  - campanha_config.json igual em: norb, norb_sint, burst_per_orbit, quick, ckpt, worker, janela_pue;
  - sementes das execucoes disjuntas entre pastas;
  - hashes de contrato, checkpoint e sintonia (rtune) iguais entre pastas;
  - avisos: numero de execucoes diferente, listas orb_* ausentes, comprimento das listas != norb, divergencias do S1 indisponiveis (sem parts/).
Saidas em --out: junta_info.json, resumo_juntas.txt, comparacao.csv (uma coluna por campanha, uma por --ref e uma 'juntas'), agregado_juntas.csv,
  m5_veredito_P1_P5_juntas.json, m5_pue_janela_juntas.json.
--ref: pastas calculadas SEPARADAMENTE (ex.: m5_final, a campanha de 6 orbitas); entram na comparacao mas NUNCA na juncao.
Os SHA256 dos codigos impressos sao os da maquina onde o script roda (LF normalizado); nao provam qual codigo gerou as simulacoes.
Uso:
  python junta_campanhas.py --autoteste
  python junta_campanhas.py --pastas m5_n20 --out m5_n20_so_B                         (so uma campanha: analise por campanha, nada a juntar)
  python junta_campanhas.py --pastas m5_n20 m5_n20_s10_14 --ref m5_final --lat m5_n20 --janela 9 --out m5_n20_juntas
"""
import argparse
import glob
import hashlib
import json
import os
import sys

import rodar_campanha_marco5 as rc
import janela_pue as jp

CHAVES_CFG = ("norb", "norb_sint", "burst_per_orbit", "quick", "ckpt", "worker", "janela_pue")
NOMES = dict(rc.NOMES)
HASHES = (("contrato", "contrato_sha256"), ("ckpt", "ckpt_sha256"), ("rtune", "rtune_sha256"))


def sha_lf(path):
    try:
        data = open(path, "rb").read().replace(b"\r\n", b"\n")
    except OSError:
        return None
    return hashlib.sha256(data).hexdigest().upper()


def ler_pasta(pasta):
    """Le uma pasta de campanha. Retorna dict com cfg, linhas (com 'bloco'), hashes, divergencias, fonte e contagens."""
    cfgp = os.path.join(pasta, "campanha_config.json")
    cfg = json.load(open(cfgp, encoding="utf-8")) if os.path.exists(cfgp) else None
    linhas, hs = [], {"contrato": set(), "ckpt": set(), "rtune": set()}
    arqs = sorted(glob.glob(os.path.join(pasta, "parts", "[1-6]", "*.json")))
    if arqs:
        fonte = "parts"
        for f in arqs:
            b = int(os.path.basename(os.path.dirname(f)))
            d = json.load(open(f, encoding="utf-8"))
            for k, kk in HASHES:
                if d.get(kk):
                    hs[k].add(d[kk])
            for r in d.get("resultados", []):
                r = dict(r)
                r["bloco"] = NOMES[b]
                linhas.append(r)
    else:
        fonte = "m5_campanha_*.json"
        for f in sorted(glob.glob(os.path.join(pasta, "m5_campanha_*.json"))):
            nome = os.path.basename(f)[len("m5_campanha_"):-len(".json")]
            if nome not in NOMES.values():
                continue
            d = json.load(open(f, encoding="utf-8"))
            for k, kk in (("contrato", "contratos"), ("ckpt", "ckpt_sha256"), ("rtune", "rtune_sha256")):
                for v in d.get(kk) or []:
                    hs[k].add(v)
            for r in d.get("resultados", []):
                r = dict(r)
                r["bloco"] = nome
                linhas.append(r)
    div = []
    for b, nome in list(NOMES.items()) + [(7, "traces")]:
        for f in sorted(glob.glob(os.path.join(pasta, "parts", str(b), "*.divergiu"))):
            d = json.load(open(f, encoding="utf-8"))
            d["bloco"] = nome
            div.append(d)
    por_bloco = {}
    for r in linhas:
        por_bloco[r["bloco"]] = por_bloco.get(r["bloco"], 0) + 1
    seeds = sorted({int(r["seed"]) for r in linhas if r.get("seed") is not None})
    n_orb = [len(r["orb_soc_end"]) for r in linhas if r.get("orb_soc_end")]
    return {"pasta": pasta, "cfg": cfg, "linhas": linhas, "hashes": hs, "diverg": div, "fonte": fonte, "por_bloco": por_bloco,
            "seeds_linhas": seeds, "n_exec": len(linhas), "n_com_orb": len(n_orb), "orb_comprimentos": sorted(set(n_orb))}


def checa_compat(infos, forca=False):
    """infos: lista de (nome, ler_pasta). Retorna (erros, avisos). Sementes repetidas sao erro mesmo com --forca."""
    erros, avisos = [], []

    def reg(msg):
        if forca:
            avisos.append("[FORCADO] " + msg)
        else:
            erros.append(msg)
    for nome, i in infos:
        if i["cfg"] is None:
            reg("%s: sem campanha_config.json" % nome)
        if i["n_exec"] == 0:
            erros.append("%s: nenhuma execucao lida (sem parts/ e sem m5_campanha_*.json)" % nome)
    for k in CHAVES_CFG:
        vals = {nome: (i["cfg"] or {}).get(k) for nome, i in infos}
        if len({json.dumps(v, sort_keys=True) for v in vals.values()}) > 1:
            reg("config diferente em '%s': %s" % (k, vals))
    nomes = [n for n, _ in infos]
    sl = {n: set(i["seeds_linhas"]) for n, i in infos}
    for x in range(len(nomes)):
        for y in range(x + 1, len(nomes)):
            inter = sl[nomes[x]] & sl[nomes[y]]
            if inter:
                erros.append("sementes repetidas entre %s e %s: %s (juntar contaria a mesma semente duas vezes)" % (nomes[x], nomes[y], sorted(inter)))
    for n, i in infos:
        cs = set(int(s) for s in (i["cfg"] or {}).get("seeds", []))
        fora = sorted(sl[n] - cs) if cs else []
        if fora:
            avisos.append("%s: sementes nas execucoes fora do campanha_config.json: %s" % (n, fora))
    for k, _ in HASHES:
        vals = {n: sorted(i["hashes"][k]) for n, i in infos}
        if len({json.dumps(v) for v in vals.values()}) > 1:
            reg("hash '%s' diferente entre pastas: %s" % (k, vals))
        for n, v in vals.items():
            if not v:
                avisos.append("%s: nenhum hash '%s' encontrado" % (n, k))
    cont = {n: i["n_exec"] for n, i in infos}
    if len(set(cont.values())) > 1:
        avisos.append("numero de execucoes diferente entre pastas: %s" % cont)
    for n, i in infos:
        if i["n_com_orb"] < i["n_exec"]:
            avisos.append("%s: %d de %d execucoes sem listas orb_* (nao rodaram pelo marco5_orbitas.py?)" % (n, i["n_exec"] - i["n_com_orb"], i["n_exec"]))
        norb = (i["cfg"] or {}).get("norb")
        if norb is not None and i["orb_comprimentos"] and set(i["orb_comprimentos"]) != {int(norb)}:
            avisos.append("%s: comprimento das listas orb_* %s != norb %s" % (n, i["orb_comprimentos"], norb))
        if i["fonte"] != "parts":
            avisos.append("%s: sem parts/ -> divergencias do S1 (.divergiu) indisponiveis; o P2 nao contara S1 que divergiu" % n)
    return erros, avisos


def le_latencia(pasta):
    p = os.path.join(pasta, "latencia_isolada.json")
    if not os.path.exists(p):
        return None
    return {r["ctrl"]: {"media_ms": r["decision_ms_mean"], "p99_ms": r["decision_ms_p99"]} for r in json.load(open(p, encoding="utf-8"))["resultados"]}


def chaves_resumo(v, pue=None):
    """Achata o veredito (e a janela k) em chave -> texto, para a tabela comparativa."""
    o = {"n_execucoes": v.get("n_execucoes")}
    p1, p2, p3, p4, p5 = (v.get(k, {}) for k in ("P1", "P2", "P3", "P4", "P5"))
    o["P1_status"], o["P1_pior_tj_tems_C"] = p1.get("status"), p1.get("pior_tj_tems_C")
    o["P1_R_SINT_pior_tj_C"], o["P1_R_SINT_grupos_com_violacao"] = p1.get("R_SINT_pior_tj_C"), len(p1.get("R_SINT_grupos_com_violacao", []))
    o["P2_status"] = p2.get("status")
    o["P2_S1_viola(com_divergencia)"] = ";".join(p2.get("cenarios_em_que_S1_viola", []))
    o["P2_S1_divergiu"] = ";".join("%s:%d" % (c, len(s)) for c, s in sorted((p2.get("S1_divergiu_por_cenario") or {}).items()))
    for r in p2.get("detalhe", []):
        o["P2_%s_%s" % (r["cenario"], r["ctrl"])] = "%d/%d tj_pior=%.3f" % (r["n_viol"], r["n"], r["tj_pior_C"]) if r.get("tj_pior_C") is not None else "%d/%d" % (r["n_viol"], r["n"])
    o["P3_status"], o["P3_soc_min_pior"], o["P3_execucoes"] = p3.get("status"), p3.get("soc_min_pior"), p3.get("execucoes")
    o["P4_status"] = p4.get("status")
    for r in p4.get("detalhe", []):
        o["P4_%s_%s" % (r["sensor"], r["ctrl"])] = "%s dif=%s n_pares=%s fora=%s" % (r.get("classe"), None if r.get("dif_media_vazao") is None else "%.5f" % r["dif_media_vazao"], r.get("n_pares"), r.get("fora_do_alcance_do_pareto"))
    cf = v.get("clausula_de_falsificacao", {})
    o["clausula_estrita"] = cf.get("estado")
    o["clausula_sens_0p25K_posthoc"] = (cf.get("sensibilidade_tol_0p25K_posthoc") or {}).get("estado")
    o["P5_status"] = p5.get("status")
    o["P5_atribuiveis/avaliados"] = "%s/%s" % (p5.get("grupos_atribuiveis"), p5.get("grupos_avaliados"))
    o["P5_grupos_atribuiveis"] = ";".join("%s/%s/%s/%s dif=%.5f" % (r["bloco"], r["cenario"], r["sensor"], r["quant"], r["dif_media"]) for r in p5.get("detalhe", []) if r.get("atribuivel_a_PINN"))
    o["P5_grupos_com_ressalva_periodicidade"] = sum(1 for r in p5.get("detalhe", []) if r.get("ressalva_periodicidade"))
    for c, d in (p5.get("latencia_decisao_ms") or {}).items():
        for kk, vv in d.items():
            o["lat_%s_%s" % (c, kk)] = vv
    for c, x in (v.get("limitacoes_n_infeasible", {}).get("por_ctrl") or {}).items():
        o["infeasible_%s" % c] = x
    for c, x in (v.get("PUE_periodico_fracao_por_ctrl(|dE|<%g)" % rc.LIM_PUE) or {}).items():
        o["PUE_janela3_%s" % c] = "%.4f" % x
    if pue:
        for c, e in sorted(pue["por_ctrl"].items()):
            o["PUE_janela%d_%s" % (pue["k"], c)] = ("%d/%d" % (e["final_ok"], e["final_tot"])) if e["final_tot"] else "-"
    return o


def tabela_csv(colunas):
    chaves = []
    for _, d in colunas:
        for k in d:
            if k not in chaves:
                chaves.append(k)
    L = ["chave," + ",".join(n.replace(",", ";") for n, _ in colunas)]
    for k in chaves:
        L.append(k + "," + ",".join(str(d.get(k, "-")).replace(",", ";") for _, d in colunas))
    return "\n".join(L)


def agregado(linhas):
    import pandas as pd
    df = pd.DataFrame(linhas)
    if df.empty:
        return "sem resultados"
    for c, pad in (("quant", "none"), ("mtotal", "nan"), ("dt_ctrl", "60.0")):
        if c not in df:
            df[c] = None
        df[c] = [rc._texto_chave(x, pad) for x in df[c]]
    for c in ("viol_int_Ks", "t_above_80_s", "tops_frac_last", "soc_min", "n_infeasible", "n_fallback", "tj_max_last_C"):
        if c not in df:
            df[c] = 0.0
    g = df.groupby(rc.KEYS + ["ctrl"], dropna=False).agg(n=("seed", "count"), sementes=("seed", "nunique"), tj_last_pior=("tj_max_last_C", "max"),
                                                        runs_viol=("viol_int_Ks", lambda x: int((x > 0).sum())), t80_total=("t_above_80_s", "sum"),
                                                        tops=("tops_frac_last", "mean"), soc_min=("soc_min", "min"), infeas=("n_infeasible", "sum"),
                                                        fallback=("n_fallback", "sum")).reset_index()
    return "# cobertura: %d de %d execucoes\n" % (int(g["n"].sum()), len(df)) + g.round(4).to_csv(index=False)


def _nat(o):
    return o.item() if hasattr(o, "item") else str(o)


def analisa(nome, linhas, div, lat, k):
    try:
        v = rc.calcula_veredito(linhas, lat, diverg=div)
    except Exception as e:
        v = {"erro": "%s: %s" % (type(e).__name__, str(e)[:300]), "n_execucoes": len(linhas)}
    pue = jp.resume(linhas, k)
    return v, pue


def main(argv=None):
    ap = argparse.ArgumentParser(description="Junta campanhas longas compativeis e recalcula P1-P5 e a janela de PUE")
    ap.add_argument("--pastas", nargs="+", help="campanhas a JUNTAR (sementes disjuntas, mesma configuracao)")
    ap.add_argument("--ref", nargs="*", default=[], help="campanhas calculadas SEPARADAMENTE (nunca juntadas)")
    ap.add_argument("--lat", default=None, help="pasta com latencia_isolada.json para a secao de latencia")
    ap.add_argument("--janela", type=int, default=None, help="janela de PUE em orbitas (padrao: janela_pue do campanha_config da 1a pasta, ou 9)")
    ap.add_argument("--forca", action="store_true", help="transforma incompatibilidades de config/hash em avisos (sementes repetidas continuam sendo erro)")
    ap.add_argument("--out", default="m5_juntas")
    ap.add_argument("--autoteste", action="store_true")
    a = ap.parse_args(argv)
    if a.autoteste:
        autoteste()
        return 0
    if not a.pastas:
        raise SystemExit("--pastas exige ao menos 1 pasta")
    for p in a.pastas + a.ref:
        if not os.path.isdir(p):
            raise SystemExit("pasta inexistente: %s" % p)
    infos = [(os.path.basename(os.path.normpath(p)), ler_pasta(p)) for p in a.pastas]
    refs = [(os.path.basename(os.path.normpath(p)) + "(ref)", ler_pasta(p)) for p in a.ref]
    k = a.janela or int(((infos[0][1]["cfg"] or {}).get("janela_pue")) or 9)
    lat = le_latencia(a.lat) if a.lat else None
    erros, avisos = checa_compat(infos, a.forca) if len(infos) > 1 else ([], [])
    for n, i in infos[:1] if len(infos) == 1 else []:
        if i["n_com_orb"] < i["n_exec"]:
            avisos.append("%s: %d de %d execucoes sem listas orb_*" % (n, i["n_exec"] - i["n_com_orb"], i["n_exec"]))
    os.makedirs(a.out, exist_ok=True)
    colunas, vereditos = [], {}
    for n, i in infos + refs:
        v, pue = analisa(n, i["linhas"], i["diverg"], lat, k)
        vereditos[n] = (v, pue)
        colunas.append((n, chaves_resumo(v, pue)))
    juntas = None
    if len(infos) > 1 and not erros:
        linhas = [r for _, i in infos for r in i["linhas"]]
        div = [d for _, i in infos for d in i["diverg"]]
        v, pue = analisa("juntas", linhas, div, lat, k)
        juntas = (linhas, v, pue)
        colunas.append(("JUNTAS", chaves_resumo(v, pue)))
        with open(os.path.join(a.out, "m5_veredito_P1_P5_juntas.json"), "w", encoding="utf-8", newline="\n") as fh:
            json.dump(v, fh, ensure_ascii=False, default=_nat)
        with open(os.path.join(a.out, "m5_pue_janela_juntas.json"), "w", encoding="utf-8", newline="\n") as fh:
            json.dump(pue, fh, ensure_ascii=False, default=_nat)
        try:
            ag = agregado(linhas)
        except Exception as e:
            ag = "ERRO no agregado: %s: %s" % (type(e).__name__, str(e)[:200])
        with open(os.path.join(a.out, "agregado_juntas.csv"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(ag)
    csv = tabela_csv(colunas)
    with open(os.path.join(a.out, "comparacao.csv"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(csv)
    codigos = ("marco5_tems.py", "marco5_orbitas.py", "rodar_campanha_marco5.py", "rodar_campanha_longa.py", "janela_pue.py", "junta_campanhas.py")
    aqui = os.path.dirname(os.path.abspath(__file__))
    info = {"pastas": [i["pasta"] for _, i in infos], "refs": [i["pasta"] for _, i in refs], "janela_pue_orbitas": k, "forca": a.forca, "erros": erros, "avisos": avisos,
            "juntou": juntas is not None, "por_pasta": {n: {"cfg": i["cfg"], "fonte": i["fonte"], "n_exec": i["n_exec"], "por_bloco": i["por_bloco"], "sementes": i["seeds_linhas"],
                                                           "hashes": {kk: sorted(vv) for kk, vv in i["hashes"].items()}, "n_divergiu": len(i["diverg"])} for n, i in infos + refs},
            "sha256_codigos_desta_maquina_LF": {c: sha_lf(os.path.join(aqui, c)) for c in codigos}}
    with open(os.path.join(a.out, "junta_info.json"), "w", encoding="utf-8", newline="\n") as fh:
        json.dump(info, fh, ensure_ascii=False, indent=1, default=_nat)
    txt = ["##### JUNTA DE CAMPANHAS #####", json.dumps({k2: info[k2] for k2 in ("pastas", "refs", "janela_pue_orbitas", "forca", "juntou")}, ensure_ascii=False),
           "##### COMPATIBILIDADE #####", "ERROS: %s" % (erros or "nenhum"), "AVISOS: %s" % (avisos or "nenhum"),
           "##### POR PASTA #####", json.dumps(info["por_pasta"], ensure_ascii=False, default=_nat),
           "##### COMPARACAO (uma coluna por campanha) #####", csv]
    if juntas is not None:
        v = juntas[1]
        cort = {kk: ({k3: (v3[:12] if isinstance(v3, list) else v3) for k3, v3 in vv.items()} if isinstance(vv, dict) else vv) for kk, vv in v.items()}
        txt += ["##### VEREDITO JUNTAS (listas cortadas em 12; completo em m5_veredito_P1_P5_juntas.json) #####", json.dumps(cort, ensure_ascii=False, default=_nat, indent=1)]
    elif len(infos) > 1:
        txt += ["##### NAO JUNTOU: corrija os ERROS acima (ou use --forca, sabendo o que isso significa) #####"]
    s = "\n".join(txt)
    with open(os.path.join(a.out, "resumo_juntas.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(s)
    print(s)
    falhas = [n for n, (v, _) in vereditos.items() if "erro" in v] + (["JUNTAS"] if (juntas is not None and "erro" in juntas[1]) else [])
    if falhas:      # AUDITORIA: antes a excecao de calcula_veredito virava {"erro":...} e o script saia com 0
        print("FALHA no calculo do veredito em: %s" % falhas, file=sys.stderr)
    return 1 if (erros or falhas) else 0


def autoteste():
    base = {"cfg": {"norb": 20, "norb_sint": 6, "burst_per_orbit": 12, "quick": False, "ckpt": "v2c.pt", "worker": "marco5_orbitas.py", "janela_pue": 9, "seeds": [10]},
            "hashes": {"contrato": {"C"}, "ckpt": {"K"}, "rtune": {"R"}}, "seeds_linhas": [10], "n_exec": 2, "n_com_orb": 2, "orb_comprimentos": [20], "fonte": "parts"}
    b = json.loads(json.dumps({k2: v2 for k2, v2 in base.items() if k2 != "hashes"}))
    b["hashes"] = {"contrato": {"C"}, "ckpt": {"K"}, "rtune": {"R"}}
    b["seeds_linhas"] = [15]
    e, av = checa_compat([("A", base), ("B", b)])
    assert e == [], e
    b2 = dict(b, seeds_linhas=[10])
    e, _ = checa_compat([("A", base), ("B", b2)])
    assert any("repetidas" in x for x in e)
    e, _ = checa_compat([("A", base), ("B", b2)], forca=True)
    assert any("repetidas" in x for x in e)
    b3 = dict(b, cfg=dict(b["cfg"], norb=6))
    e, _ = checa_compat([("A", base), ("B", b3)])
    assert any("'norb'" in x for x in e)
    e, av = checa_compat([("A", base), ("B", b3)], forca=True)
    assert not any("'norb'" in x for x in e) and any("FORCADO" in x for x in av)
    print("autoteste OK")


if __name__ == "__main__":
    sys.exit(main())
