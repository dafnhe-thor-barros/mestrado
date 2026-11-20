"""modelo_ia.py -- Marco 2d: carga de IA (DVFS + fuga) -> P_TI(f_clk, V_core, T_j) e termo fonte q_gen.
Modulo de referencia: NVIDIA A100 SXM4 (TDP 400 W). 100 modulos = 40 kW (capacidade de pico, decisao do usuario).
Etiquetas: [FONTE-RESUMO] valor visto so em resumo de busca (PDF NAO aberto por mim) | [CALCULADO] | [HIPOTESE DECLARADA] (com sensibilidade)
| [PENDENTE] sem fonte | [NAO VERIFICADO] memoria | [LEGADO] valor do legado, sem fonte.
Modelo (por modulo):  P = P_fix + P_dyn + P_leak ;  P_dyn = util (1-fix-leak) TDP (V/Vnom)^2 (f/fnom) = alpha C V^2 f normalizado ;
P_leak = V * I_leak(T_j) ;  I_leak = (leak TDP/Vnom) (V/Vnom)^(kv-1) exp(beta (T_j - T_ref)).  Calibrado: P(fnom, Vnom, T_ref, util=1) = TDP.
NAO modela: transitorios de f/V, acoplamento memoria-HBM detalhado, efeitos de radiacao (A100 comercial nao e qualificado para espaco),
curva real V(f) (nao publica; ver ressalva), capacitancia termica (Marco 3). Nao validado experimentalmente.
Dependencias: numpy. Cruzamento opcional com modelo_barramento.py na mesma pasta.
"""
import math
import json
import numpy as np

DEFAULTS = {
    "tdp": 400.0,        # W por modulo [FONTE-RESUMO] datasheet NVIDIA A100 80GB SXM4 (CTS ate 500 W)
    "tops": 624.0,       # TOPS INT8 densos por modulo [FONTE-RESUMO] mesmo datasheet (1248 com esparsidade)
    "f_nom": 1410.0,     # MHz [NAO VERIFICADO] memoria; PDF nao aberto
    "f_min": 300.0,      # MHz [HIPOTESE DECLARADA]
    "v_nom": 0.80,       # V [HIPOTESE DECLARADA] = saida do estagio POL do enunciado (0.8 V)
    "v_min": 0.65,       # V [HIPOTESE DECLARADA] piso de tensao; sem fonte
    "fv": 0.40,          # fracao de f_nom abaixo da qual V = v_min [HIPOTESE DECLARADA]
    "fix": 0.20,         # fracao do TDP que nao escala com DVFS do core (HBM, SerDes...) [HIPOTESE DECLARADA]
    "leak": 0.20,        # fracao do TDP que e fuga em (f_nom, v_nom, t_ref) [HIPOTESE DECLARADA]
    "t_ref": 80.0,       # C, T_j em que o TDP e calibrado [HIPOTESE DECLARADA]
    "beta": 0.02,        # 1/K, inclinacao exponencial da fuga [HIPOTESE DECLARADA]
    "kv": 3.0,           # expoente de tensao da potencia de fuga [HIPOTESE DECLARADA]
    "p_off": 0.0,        # W por modulo desligado [HIPOTESE DECLARADA]
    "eta_dab": 0.9655,   # Marco 2b [HIPOTESE DECLARADA]
    "eta_pol": 0.90,     # Marco 2b [HIPOTESE DECLARADA]
    "q_oisl": 343.2,     # W calor do OISL [LEGADO 4_1]
    "q_aux_other": 560.0,  # W calor de bomba+aviônica+margem (200+350+10) [LEGADO 4_1]
}
N_MOD_REF = 100          # 100 x 400 W = 40 kW [CALCULADO]
TJ_LIMIT_C = 80.0        # limite de controle do projeto


def _h(h=None):
    d = dict(DEFAULTS)
    if h:
        bad = set(h) - set(DEFAULTS)
        if bad:
            raise ValueError("parametros desconhecidos: %s" % sorted(bad))
        d.update(h)
    for k, v in d.items():
        if not (isinstance(v, (int, float)) and math.isfinite(v)):
            raise ValueError("%s deve ser numero finito" % k)
    if not (d["tdp"] > 0 and d["tops"] > 0 and 0.0 < d["f_min"] < d["f_nom"]):
        raise ValueError("tdp, tops > 0 e 0 < f_min < f_nom")
    if not (0.0 < d["v_min"] <= d["v_nom"] and 0.0 < d["fv"] < 1.0 and d["kv"] >= 1.0 and d["beta"] >= 0.0):
        raise ValueError("0 < v_min <= v_nom; 0 < fv < 1; kv >= 1; beta >= 0")
    if not (0.0 <= d["fix"] and 0.0 <= d["leak"] and d["fix"] + d["leak"] < 1.0 and d["p_off"] >= 0.0):
        raise ValueError("fix, leak >= 0 com fix+leak < 1; p_off >= 0")
    if not (0.0 < d["eta_dab"] <= 1.0 and 0.0 < d["eta_pol"] <= 1.0):
        raise ValueError("eficiencias em (0, 1]")
    return d


def v_of_f(f, h=None):
    """V(f): piso v_min ate f_v = fv f_nom; depois linear ate v_nom em f_nom. Contínua. [HIPOTESE DECLARADA]"""
    d = _h(h)
    if not (d["f_min"] - 1e-9 <= f <= d["f_nom"] + 1e-9):
        raise ValueError("f fora de [f_min, f_nom]")
    fv = d["fv"] * d["f_nom"]
    if f <= fv:
        return d["v_min"]
    return d["v_min"] + (d["v_nom"] - d["v_min"]) * (f - fv) / (d["f_nom"] - fv)


def p_components(f, tj, util=1.0, h=None):
    d = _h(h)
    if not (0.0 <= util <= 1.0):
        raise ValueError("util deve estar em [0, 1]")
    if not (-100.0 <= tj <= 300.0):
        raise ValueError("tj fora de [-100, 300] C")
    v = v_of_f(f, d)
    p_fix = d["fix"] * d["tdp"]
    p_dyn = util * (1.0 - d["fix"] - d["leak"]) * d["tdp"] * (v / d["v_nom"]) ** 2 * (f / d["f_nom"])
    i_leak = d["leak"] * d["tdp"] / d["v_nom"] * (v / d["v_nom"]) ** (d["kv"] - 1.0) * math.exp(d["beta"] * (tj - d["t_ref"]))
    p_leak = v * i_leak
    return {"v": v, "p_fix": p_fix, "p_dyn": p_dyn, "i_leak_A": i_leak, "p_leak": p_leak, "p_mod": p_fix + p_dyn + p_leak}


def p_mod(f, tj, util=1.0, h=None):
    return p_components(f, tj, util, h)["p_mod"]


def f_for_power(p_target, tj, util=1.0, h=None):
    """Inverte P(f) em f por bissecao (P e nao-decrescente em f). Retorna (f_MHz, flag) com flag em {'ok','piso','teto'}."""
    d = _h(h)
    lo, hi = d["f_min"], d["f_nom"]
    if p_target <= p_mod(lo, tj, util, d):
        return lo, "piso"
    if p_target >= p_mod(hi, tj, util, d):
        return hi, "teto"
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if p_mod(mid, tj, util, d) < p_target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), "ok"


def thermal_fixed_point(f, t_cold, rth, util=1.0, h=None, tol=1e-9, itmax=1000):
    """T_j = t_cold + rth P(f, T_j). rth em K/W por modulo [PENDENTE: vem do Marco 3]. Detecta fuga termica (runaway)."""
    d = _h(h)
    if rth < 0.0:
        raise ValueError("rth >= 0")
    tj = t_cold
    for k in range(itmax):
        p = p_mod(f, tj, util, d)
        new = t_cold + rth * p
        if new > 300.0:
            return {"tj": new, "p": p, "converged": False, "runaway": True, "iters": k}
        if abs(new - tj) < tol:
            tj = new
            p = p_mod(f, tj, util, d)
            return {"tj": tj, "p": p, "converged": True, "runaway": False, "iters": k, "resid": abs(t_cold + rth * p - tj)}
        tj = new
    return {"tj": tj, "p": p_mod(f, tj, util, d), "converged": False, "runaway": False, "iters": itmax}


def rth_critical(f, tj, util=1.0, h=None):
    """Estabilidade local: rth * dP/dT_j < 1, com dP/dT_j = beta P_leak. Retorna rth critico [K/W]."""
    d = _h(h)
    c = p_components(f, tj, util, d)
    return float("inf") if d["beta"] * c["p_leak"] == 0.0 else 1.0 / (d["beta"] * c["p_leak"])


def heat_chain(p_ti, h=None):
    """P_bus = P_TI/(eta_dab eta_pol); Q_pol + Q_dab = P_bus - P_TI.  q_gen = P_TI + Q_perdas + Q_OISL (definicao do enunciado).
    q_aux_other e calor de bomba/aviônica do legado, FORA da definicao do enunciado; q_node_total = q_gen + q_aux_other."""
    d = _h(h)
    if not (math.isfinite(p_ti) and p_ti >= 0.0):
        raise ValueError("p_ti >= 0")
    p48 = p_ti / d["eta_pol"]
    p400 = p48 / d["eta_dab"]
    q_perdas = p400 - p_ti
    return {"p_bus": p400, "q_pol": p48 - p_ti, "q_dab": p400 - p48, "q_perdas": q_perdas, "q_oisl": d["q_oisl"],
            "q_gen": p_ti + q_perdas + d["q_oisl"], "q_aux_other": d["q_aux_other"], "q_node_total": p_ti + q_perdas + d["q_oisl"] + d["q_aux_other"]}


def throughput_tops(n_act, f, util=1.0, h=None):
    """Escala linear em f (carga limitada por computo) [HIPOTESE DECLARADA]; cargas limitadas por memoria escalam menos."""
    d = _h(h)
    return n_act * d["tops"] * (f / d["f_nom"]) * util


def best_allocation(p_budget, tj, util=1.0, n_total=N_MOD_REF, h=None):
    """Maximiza a vazao com orcamento de potencia nos dies: escolhe n_act modulos ligados, igual f em todos, resto desligado (p_off).
    Busca exaustiva em n_act. Nao e o TEMS: e a referencia estatica de alocacao (Marco 5 faz a versao dinamica)."""
    d = _h(h)
    if not (math.isfinite(p_budget) and p_budget >= 0.0) or n_total < 1:
        raise ValueError("p_budget >= 0 e n_total >= 1")
    pmin, pmax = p_mod(d["f_min"], tj, util, d), p_mod(d["f_nom"], tj, util, d)
    best = None
    for n in range(1, int(n_total) + 1):
        per = (p_budget - (n_total - n) * d["p_off"]) / n
        if per < pmin - 1e-12:
            continue
        f, flag = f_for_power(min(per, pmax), tj, util, d)
        thr = throughput_tops(n, f, util, d)
        if best is None or thr > best["tops"] + 1e-12:
            used = n * p_mod(f, tj, util, d) + (n_total - n) * d["p_off"]
            best = {"feasible": True, "n_act": n, "f_MHz": f, "flag": flag, "tops": thr, "p_used_W": used,
                    "p_per_module_W": p_mod(f, tj, util, d), "tops_per_W": thr / used if used > 0 else float("nan")}
    return best if best else {"feasible": False, "p_min_modulo_W": pmin}


# ---------------------------------------------------------------- testes
def test_calibracao_tdp():
    c = p_components(DEFAULTS["f_nom"], DEFAULTS["t_ref"], 1.0)
    assert abs(c["p_mod"] - DEFAULTS["tdp"]) < 1e-9, c
    assert abs(c["p_leak"] - DEFAULTS["leak"] * DEFAULTS["tdp"]) < 1e-9 and abs(N_MOD_REF * DEFAULTS["tdp"] - 40000.0) < 1e-9


def test_decomposicao_igual_formula_do_enunciado():
    c = p_components(900.0, 70.0, 0.8)
    v = v_of_f(900.0)
    assert abs(c["p_leak"] - v * c["i_leak_A"]) < 1e-12
    assert abs(c["p_mod"] - (c["p_fix"] + c["p_dyn"] + v * c["i_leak_A"])) < 1e-12


def test_escala_dinamica_v2f():
    a, b = p_components(700.0, 80.0), p_components(1200.0, 80.0)
    r = (b["v"] / a["v"]) ** 2 * (1200.0 / 700.0)
    assert abs(b["p_dyn"] / a["p_dyn"] - r) < 1e-12


def test_monotonia_e_continuidade():
    fs = np.linspace(DEFAULTS["f_min"], DEFAULTS["f_nom"], 400)
    v = [v_of_f(float(f)) for f in fs]
    assert all(v[i + 1] >= v[i] - 1e-15 for i in range(399)) and abs(v[0] - DEFAULTS["v_min"]) < 1e-12 and abs(v[-1] - DEFAULTS["v_nom"]) < 1e-12
    assert max(abs(v[i + 1] - v[i]) for i in range(399)) < 2e-3
    p = [p_mod(float(f), 80.0) for f in fs]
    assert all(p[i + 1] >= p[i] - 1e-12 for i in range(399))
    assert p_mod(900.0, 90.0) > p_mod(900.0, 70.0) and p_mod(900.0, 80.0, 1.0) > p_mod(900.0, 80.0, 0.5)


def test_inversao_ida_e_volta():
    for tj in (60.0, 80.0):
        for f in (400.0, 800.0, 1300.0):
            p = p_mod(f, tj)
            f2, flag = f_for_power(p, tj)
            assert abs(f2 - f) < 1e-6 and flag == "ok", (f, f2, flag)
    assert f_for_power(1.0, 80.0)[1] == "piso" and f_for_power(1e6, 80.0)[1] == "teto"


def test_ponto_fixo_termico():
    beta = DEFAULTS["beta"]
    r = thermal_fixed_point(1000.0, 60.0, 0.05)
    assert r["converged"], ("a_convergiu", r)
    assert r["resid"] < 1e-8, ("b_resid", r)
    assert abs(r["tj"] - (60.0 + 0.05 * r["p"])) < 1e-8, ("c_consistencia", r)
    assert 0.05 * beta * p_components(1000.0, r["tj"])["p_leak"] < 1.0, ("d_estabilidade_no_ponto_fixo", r)
    r2 = thermal_fixed_point(1410.0, 40.0, 0.05)
    assert r2["converged"], ("e_fnom_rth_moderado", r2)
    assert 0.05 * beta * p_components(1410.0, r2["tj"])["p_leak"] < 1.0, ("f_estabilidade_fnom", r2)
    r3 = thermal_fixed_point(1410.0, 40.0, 5.0 * rth_critical(1410.0, 80.0))
    assert r3["runaway"], ("g_rth_muito_alto_sem_ponto_fixo", r3)
    r4 = thermal_fixed_point(1410.0, 40.0, 0.5 * rth_critical(1410.0, 80.0))
    assert r4["runaway"] and not r4["converged"], ("h_sem_ponto_fixo_apesar_de_rth<rc_local", r4)


def test_fechamento_do_calor():
    p = 30000.0
    c = heat_chain(p)
    assert abs((c["p_bus"] - p) - (c["q_pol"] + c["q_dab"])) < 1e-9
    assert abs(c["q_gen"] - (p + c["q_perdas"] + c["q_oisl"])) < 1e-9 and abs(c["q_node_total"] - c["q_gen"] - c["q_aux_other"]) < 1e-9


def test_alocacao_respeita_orcamento_e_monotonia():
    prev = -1.0
    for b in (16000.0, 18000.0, 25000.0, 33904.08, 40000.0):
        r = best_allocation(b, 80.0)
        assert r["feasible"] and r["p_used_W"] <= b + 1e-6, (b, r)
        assert r["tops"] >= prev - 1e-9, (b, r["tops"], prev)
        prev = r["tops"]
    assert not best_allocation(1000.0, 80.0, n_total=100)["feasible"] or best_allocation(1000.0, 80.0)["n_act"] < 100


def test_dominio():
    for bad in (lambda: p_mod(100.0, 80.0), lambda: p_mod(900.0, 80.0, 1.5), lambda: p_mod(900.0, 400.0), lambda: p_mod(900.0, 80.0, 1.0, {"x": 1}),
                lambda: p_mod(900.0, 80.0, 1.0, {"fix": 0.6, "leak": 0.5}), lambda: heat_chain(-1.0), lambda: best_allocation(-5.0, 80.0)):
        try:
            bad()
        except ValueError:
            continue
        raise AssertionError("deveria recusar entrada invalida")


def test_reprodutibilidade():
    assert best_allocation(18000.0, 80.0) == best_allocation(18000.0, 80.0)


def _cross_check():
    try:
        import modelo_barramento as mb
    except ImportError:
        return "SKIPPED (modelo_barramento.py ausente)"
    c, o = heat_chain(30000.0), mb.load_chain(30000.0)
    assert abs(c["p_bus"] - o["p_bus"]) < 1e-9 and abs(c["q_pol"] - o["q_pol"]) < 1e-9 and abs(c["q_dab"] - o["q_dab"]) < 1e-9
    return "PASS"


def test_conferencia_cruzada_com_modelo_barramento():
    assert _cross_check().split()[0] in ("PASS", "SKIPPED")


TESTS = [test_calibracao_tdp, test_decomposicao_igual_formula_do_enunciado, test_escala_dinamica_v2f, test_monotonia_e_continuidade, test_inversao_ida_e_volta,
         test_ponto_fixo_termico, test_fechamento_do_calor, test_alocacao_respeita_orcamento_e_monotonia, test_dominio, test_reprodutibilidade]


def run_selftests():
    res = {}
    for f in TESTS:
        try:
            f()
            res[f.__name__] = "PASS"
        except AssertionError as e:
            res[f.__name__] = "FAIL: " + str(e)[:300]
        except Exception as e:
            res[f.__name__] = "ERROR " + type(e).__name__ + ": " + str(e)[:300]
    try:
        res["test_conferencia_cruzada_com_modelo_barramento"] = _cross_check()
    except Exception as e:
        res["test_conferencia_cruzada_com_modelo_barramento"] = "FAIL: " + type(e).__name__ + " " + str(e)[:300]
    return res


def run_report():
    d = _h()
    out = {"calibracao_modulo_W": p_components(d["f_nom"], d["t_ref"], 1.0)}
    rows = []
    for tj in (60.0, 80.0):
        for f in (300.0, 500.0, 700.0, 900.0, 1100.0, 1410.0):
            c = p_components(f, tj, 1.0)
            rows.append([tj, f, c["v"], c["p_dyn"], c["p_leak"], c["p_mod"], throughput_tops(1, f), throughput_tops(1, f) / c["p_mod"]])
    out["modulo_[Tj | f_MHz | V | P_dyn | P_leak | P_mod_W | TOPS | TOPS_por_W]"] = rows
    out["piso_modulo_W_[Tj=60,80]"] = [p_mod(d["f_min"], 60.0), p_mod(d["f_min"], 80.0)]
    budgets = {"sol_sustentavel_2c": 33904.08, "eclipse_18kW": 18000.0, "pico_40kW": 40000.0}
    al = {}
    for nome, b in budgets.items():
        al[nome] = {"tj80": best_allocation(b, 80.0), "tj60": best_allocation(b, 60.0), "heat_80": heat_chain(min(b, 40000.0))}
    out["alocacao"] = al
    out["f_todos_os_modulos_[W_por_modulo -> f_MHz,flag] Tj=80"] = {k: f_for_power(v / N_MOD_REF, 80.0) for k, v in budgets.items()}
    out["rth_critico_K_por_W_modulo_[f=fnom;Tj=80]"] = rth_critical(d["f_nom"], 80.0)
    fp = thermal_fixed_point(d["f_nom"], 40.0, 0.10)
    out["ponto_fixo_exemplo_[t_cold=40,rth=0.10]"] = fp
    sens = []
    for nome, lo, hi in (("fix", 0.10, 0.30), ("leak", 0.10, 0.35), ("beta", 0.01, 0.03), ("v_min", 0.60, 0.72), ("t_ref", 60.0, 100.0),
                         ("kv", 2.0, 4.0), ("fv", 0.20, 0.60), ("f_min", 210.0, 450.0)):
        for val in (lo, hi):
            h = {nome: val}
            try:
                r = best_allocation(18000.0, 80.0, 1.0, N_MOD_REF, h)
                sens.append([nome, val, p_mod(h.get("f_min", d["f_min"]), 80.0, 1.0, h), f_for_power(180.0, 80.0, 1.0, h)[0],
                             r.get("n_act"), r.get("f_MHz"), r.get("tops")])
            except ValueError as e:
                sens.append([nome, val, "ERRO", str(e), None, None, None])
    out["sensibilidade_OAT_[param | valor | piso_modulo_W | f_MHz_p/_180W | n_act_eclipse | f_MHz_eclipse | TOPS_eclipse]"] = sens
    out["legado_80dies_x_250TOPS_TOPS"] = 80 * 250.0
    out["novo_100_modulos_x_624_TOPS_INT8_denso_TOPS"] = N_MOD_REF * d["tops"]
    return out


if __name__ == "__main__":
    tests = run_selftests()
    falhas = {k: v for k, v in tests.items() if not (v == "PASS" or v.startswith("SKIPPED"))}
    if falhas:
        raise RuntimeError("TESTES FALHARAM: " + json.dumps(falhas, ensure_ascii=False))
    payload = {"testes": tests, "relatorio": run_report()}
    print(json.dumps(payload, indent=1, default=float, ensure_ascii=False))
    result = {"json": json.dumps(payload, default=float, ensure_ascii=False)}
