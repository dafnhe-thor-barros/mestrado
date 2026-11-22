"""janela_pue.py -- janela de PUE periodico calculada FORA do codigo congelado (opcao b). Python puro; numpy so em por_orbita/enriquece.
ESCRITO SEM TER SIDO EXECUTADO.

Definicao (a mesma de marco5_tems.summarize.janela, L559-571, mas com k livre):
  rel(n, k) = |SoC_fim(n) - SoC_fim(n-k)| * E_cap / soma(E_TI das orbitas n-k+1..n)      (exige n >= k+1)
  'passa' se rel < LIM_PUE (1e-3). Com k = 9 o menor n valido e 10.
Os JSONs de resultado ganham (via marco5_orbitas.py) as listas por orbita orb_soc_end, orb_ti_Wh, orb_tops, orb_tj_max_C; a janela k e calculada depois,
aqui, sem tocar no marco5_tems.py. A janela de 9 orbitas e emenda POS-HOC (escolhida depois do teste de convergencia de 1 semente, RC3, sem_calibracao).
"""
LIM_PUE = 1e-3


def por_orbita(ser, T, fmax, soc0, norb, t0k=273.15):
    """ser: dict de listas (t, h, tj [K], soc, f [MHz], p_ti [W]) da simulacao inteira. Retorna listas por orbita."""
    import numpy as np
    t, h = np.array(ser["t"], float), np.array(ser["h"], float)
    soc, f, p_ti = np.array(ser["soc"], float), np.array(ser["f"], float), np.array(ser["p_ti"], float)
    tj = np.array(ser["tj"], float) - t0k
    out = {"orb_soc_end": [], "orb_ti_Wh": [], "orb_tops": [], "orb_tj_max_C": []}
    for k in range(int(norb)):
        m = (t > k * T + 1e-9) & (t <= (k + 1) * T + 1e-9)
        if not m.any():
            break
        out["orb_soc_end"].append(float(soc[m][-1]))
        out["orb_ti_Wh"].append(float(np.sum(p_ti[m] * h[m]) / 3600.0))
        out["orb_tops"].append(float(np.sum(f[m] * h[m]) / np.sum(h[m]) / fmax))
        out["orb_tj_max_C"].append(float(tj[m].max()))
    return out


def enriquece(r, T, fmax, soc0, norb, t0k=273.15, nd=7):
    """Acrescenta ao dict de resultado r (com r['series'] e r['ecap_Wh']) as listas por orbita. Nao remove 'series'."""
    o = por_orbita(r["series"], T, fmax, soc0, norb, t0k)
    for k, v in o.items():
        r[k] = [round(x, nd) for x in v]
    r["orb_T_s"] = float(T)
    r["orb_soc0"] = float(soc0)
    return r


def rel_ponta(soc, ti, ecap, k, n):
    """Valor relativo da janela de k orbitas terminando na orbita n (1-based). None se n < k+1 ou sem energia de TI."""
    if k < 1 or n < k + 1 or n > len(soc) or n > len(ti):
        return None
    tt = sum(ti[n - k:n])
    if tt <= 0:
        return None
    return abs(soc[n - 1] - soc[n - k - 1]) * ecap / tt


def avalia_run(r, k, lim=LIM_PUE):
    """Resultado da janela k para um registro de execucao. None se o registro nao tem as listas por orbita."""
    soc, ti, ecap = r.get("orb_soc_end"), r.get("orb_ti_Wh"), r.get("ecap_Wh")
    if not soc or not ti or ecap is None or len(soc) != len(ti):
        return None
    N = len(soc)
    vals = []
    for n in range(k + 1, N + 1):
        v = rel_ponta(soc, ti, float(ecap), k, n)
        if v is not None:
            vals.append((n, v))
    if not vals:
        return {"N": N, "final": None, "final_ok": None, "n_pass": 0, "n_tot": 0, "max_rel": None}
    nf, vf = vals[-1]
    return {"N": N, "final": vf if nf == N else None, "final_ok": (vf < lim) if nf == N else None,
            "n_pass": sum(1 for _, v in vals if v < lim), "n_tot": len(vals), "max_rel": max(v for _, v in vals)}


def _txt(x, pad):
    try:
        if x is None or x != x:
            return pad
    except Exception:
        return pad
    return str(x)


def resume(linhas, k, lim=LIM_PUE):
    """Agrega por controlador e por grupo (cenario, sensor, quant, mtotal, dt_ctrl, ctrl)."""
    por_ctrl, grupos = {}, {}
    for r in linhas:
        a = avalia_run(r, k, lim)
        ctrl = _txt(r.get("ctrl"), "?")
        gk = (_txt(r.get("cenario"), "?"), _txt(r.get("sensor"), "?"), _txt(r.get("quant"), "none"), _txt(r.get("mtotal"), "nan"), _txt(r.get("dt_ctrl"), "60.0"), ctrl)
        for d, key in ((por_ctrl, ctrl), (grupos, gk)):
            e = d.setdefault(key, {"n_exec": 0, "sem_orb": 0, "curtas": 0, "final_ok": 0, "final_tot": 0, "frac_n": [], "pior": []})
            e["n_exec"] += 1
            if a is None:
                e["sem_orb"] += 1
            elif a["n_tot"] == 0 or a["final"] is None:
                e["curtas"] += 1
            else:
                e["final_tot"] += 1
                e["final_ok"] += 1 if a["final_ok"] else 0
                e["frac_n"].append(a["n_pass"] / a["n_tot"])
                e["pior"].append(a["max_rel"])
    return {"k": k, "limiar": lim, "por_ctrl": por_ctrl, "grupos": {"|".join(g): v for g, v in grupos.items()}}


def linha_csv(nome, e):
    fo = ("%.4f" % (e["final_ok"] / e["final_tot"])) if e["final_tot"] else "-"
    fn = ("%.4f" % (sum(e["frac_n"]) / len(e["frac_n"]))) if e["frac_n"] else "-"
    pw = ("%.2e" % max(e["pior"])) if e["pior"] else "-"
    return "%s,%d,%d,%d,%s,%s,%s" % (nome, e["n_exec"], e["sem_orb"], e["curtas"], fo, fn, pw)


CABECALHO_CSV = "chave,n_exec,sem_listas_por_orbita,orbitas_insuficientes,frac_passa_no_n_final,frac_media_dos_n_que_passam,pior_valor_relativo"
