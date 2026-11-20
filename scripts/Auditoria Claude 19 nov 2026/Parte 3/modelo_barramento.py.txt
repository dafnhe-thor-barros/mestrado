"""modelo_barramento.py -- Marco 2b: balanco de potencia do barramento HVDC regulado de 400 V (nivel macro, quase-estatico por fase orbital).
Arquitetura (fronteira sem dupla contagem; cada perda aparece UMA vez):
  arranjo FV (Marco 2a, potencia na SAIDA do arranjo) --Buck/MPPT--> barramento 400 V <--conversor bidirecional--> BESS (terminais)
  barramento 400 V --DAB (400->48 V)--> 48 V --POL/VRM (48->0.8 V)--> dies (P_TI)  ;  cargas auxiliares ligadas ao barramento (P_aux, eficiencia embutida).
Etiquetas: [CALCULADO] | [HIPOTESE DECLARADA] (com sensibilidade) | [PENDENTE] (sem fonte) | [LEGADO] valor do legado, sem fonte lida.
NAO modela: dinamica do BESS (SoC, tensao, resistencia interna -> Marco 2c), curva de eficiencia vs carga, transitorios, SSPC, Middlebrook.
Nao validado experimentalmente.
"""
import math
import json
import numpy as np

V_BUS = 400.0                 # V, barramento regulado [HIPOTESE DECLARADA: decisao do usuario]
ETA_BUCK = 0.91               # [PENDENTE] legado 3_2 (sem fonte); estagio Buck/MPPT do arranjo
ETA_DAB = 0.9655              # [HIPOTESE DECLARADA] 400->48 V; valor CALCULADO a partir do legado 3_4 a 20 kW (entradas [HIP]/[VERIFICAR]); legado 3_1/3_3/4_1 usa 0.95
ETA_POL = 0.90                # [HIPOTESE DECLARADA] 48->0.8 V; nenhum script do legado modela este estagio; sem fonte
ETA_BD = 0.97                 # [HIPOTESE DECLARADA] conversor bidirecional BESS<->400 V, por sentido; nao existe no legado; sem fonte
P_AUX_W = 1000.0              # W [LEGADO 4_1]: bomba 200 + avionica 350 + OISL 440 + margem 10 (eficiencia propria ja embutida)
Q_AUX_HEAT_W = 903.2          # W [LEGADO 4_1]: calor interno de P_AUX_W (OISL: 343.2 W de calor + 96.8 W de potencia optica que sai do veiculo)
SOC_MIN, SOC_MAX = 0.70, 0.95  # janela util do BESS [decisao do usuario]
PV_OUT_DESIGN_W = 50000.0     # W na saida do arranjo, EOL, antes do Buck [decisao do usuario; Marco 2a]
HEAT_IN_NODE = {"buck": True, "bd": True}   # [HIPOTESE DECLARADA] perdas do Buck e do conversor do BESS rejeitadas pelo mesmo radiador do no


def _nn(x, name):
    if not (isinstance(x, (int, float, np.floating)) and math.isfinite(x) and x >= 0.0):
        raise ValueError("%s deve ser finito e >= 0 (recebido %r)" % (name, x))


def _eta(x, name):
    if not (0.0 < x <= 1.0):
        raise ValueError("%s deve estar em (0, 1] (recebido %r)" % (name, x))


def load_chain(p_ti, eta_dab=ETA_DAB, eta_pol=ETA_POL):
    """P_48 = P_TI/eta_pol ; P_400 = P_48/eta_dab ; Q_pol = P_48 - P_TI ; Q_dab = P_400 - P_48.  [W]"""
    _nn(p_ti, "p_ti"); _eta(eta_dab, "eta_dab"); _eta(eta_pol, "eta_pol")
    p48 = p_ti / eta_pol
    p400 = p48 / eta_dab
    return {"p_bus": p400, "p_48": p48, "q_pol": p48 - p_ti, "q_dab": p400 - p48}


def bus_flow(p_pv_out, p_ti, p_aux=P_AUX_W, eta_buck=ETA_BUCK, eta_dab=ETA_DAB, eta_pol=ETA_POL, eta_bd=ETA_BD, p_chg_max_bus=float("inf")):
    """Balanco instantaneo no barramento. Convencao: s = potencia nos TERMINAIS do BESS, >0 carregando, <0 descarregando.
    Excedente vai ao BESS (limitado por p_chg_max_bus); o que sobrar e curtailment do arranjo (MPPT recua; Buck so dissipa o que processa)."""
    _nn(p_pv_out, "p_pv_out"); _nn(p_aux, "p_aux"); _eta(eta_buck, "eta_buck"); _eta(eta_bd, "eta_bd")
    if not (p_chg_max_bus >= 0.0):   # aceita +inf (sem limite); recusa negativo e NaN
        raise ValueError("p_chg_max_bus deve ser >= 0 ou inf (recebido %r)" % (p_chg_max_bus,))
    ch = load_chain(p_ti, eta_dab, eta_pol)
    l_bus = ch["p_bus"] + p_aux
    p_pv_bus = eta_buck * p_pv_out
    net = p_pv_bus - l_bus
    if net >= 0.0:
        to_bat = min(net, p_chg_max_bus)
        used_bus = l_bus + to_bat
        s = eta_bd * to_bat
        q_bd = to_bat - s
    else:
        deficit = -net
        used_bus = p_pv_bus
        p_term = deficit / eta_bd
        s = -p_term
        q_bd = p_term - deficit
    used_out = used_bus / eta_buck
    q_buck = used_out - used_bus
    curtail = p_pv_out - used_out
    heat_aux = p_aux * (Q_AUX_HEAT_W / P_AUX_W)
    opt_aux = p_aux - heat_aux
    q_node = p_ti + ch["q_pol"] + ch["q_dab"] + heat_aux + (q_buck if HEAT_IN_NODE["buck"] else 0.0) + (q_bd if HEAT_IN_NODE["bd"] else 0.0)
    resid = (p_pv_out - curtail) - s - (p_ti + ch["q_pol"] + ch["q_dab"] + p_aux + q_buck + q_bd)
    return {"p_pv_bus_disp": p_pv_bus, "l_bus": l_bus, "s_bat_term": s, "q_pol": ch["q_pol"], "q_dab": ch["q_dab"], "q_buck": q_buck, "q_bd": q_bd,
            "curtail_array": curtail, "heat_aux": heat_aux, "opt_aux": opt_aux, "q_node": q_node, "residuo_W": resid}


def orbit_budget(t_ecl, T, p_pv_out=PV_OUT_DESIGN_W, p_ti_ecl=18000.0, p_ti_sun=None, p_aux=P_AUX_W, eta_buck=ETA_BUCK, eta_dab=ETA_DAB,
                 eta_pol=ETA_POL, eta_bd=ETA_BD, e_nom_j=None):
    """Regime periodico, bateria ideal alem do conversor (eficiencia coulombica 1: [HIPOTESE DECLARADA], refinada no 2c). Duas fases: sol (t_sun) e eclipse (t_ecl).
    Eclipse: L_ecl = P_TI_ecl/(eta_dab eta_pol) + P_aux ; E_term = L_ecl t_ecl/eta_bd.  Sol: carga = eta_bd (eta_buck P_pv - L_sun) t_sun.
    Fechamento periodico: carga = E_term  =>  L_sun = eta_buck P_pv - E_term/(eta_bd t_sun) ;  P_TI_sun = (L_sun - P_aux) eta_dab eta_pol.  [CALCULADO]"""
    _nn(t_ecl, "t_ecl"); _nn(T, "T")
    if t_ecl >= T:
        raise ValueError("t_ecl deve ser menor que T")
    t_sun = T - t_ecl
    l_ecl = load_chain(p_ti_ecl, eta_dab, eta_pol)["p_bus"] + p_aux
    e_term = l_ecl * t_ecl / eta_bd
    out = {"t_sun_s": t_sun, "L_ecl_bus_W": l_ecl, "E_term_ecl_Wh": e_term / 3600.0, "E_nom_req_Wh": e_term / 3600.0 / (SOC_MAX - SOC_MIN)}
    l_sun_max = eta_buck * p_pv_out - e_term / (eta_bd * t_sun)
    out["P_TI_sun_sustentavel_W"] = max(0.0, (l_sun_max - p_aux) * eta_dab * eta_pol)
    out["viavel_sol"] = bool(l_sun_max >= p_aux)
    out["P_TI_sun_sem_carga_W"] = max(0.0, (eta_buck * p_pv_out - p_aux) * eta_dab * eta_pol)
    if p_ti_sun is not None:
        l_sun = load_chain(p_ti_sun, eta_dab, eta_pol)["p_bus"] + p_aux
        out["saldo_energia_orbita_Wh"] = (eta_bd * (eta_buck * p_pv_out - l_sun) * t_sun - e_term) / 3600.0
    if e_nom_j is not None:
        e_win = (SOC_MAX - SOC_MIN) * e_nom_j
        l_ecl_max = e_win * eta_bd / t_ecl
        out["P_TI_ecl_max_para_E_nom_W"] = max(0.0, (l_ecl_max - p_aux) * eta_dab * eta_pol)
    return out


# ---------------------------------------------------------------- testes (compativeis com pytest)
def test_fechamento_aleatorio():
    rng = np.random.default_rng(12345)
    worst = 0.0
    for _ in range(2000):
        a = rng.uniform(0.0, 60e3); b = rng.uniform(0.0, 45e3)
        e = rng.uniform(0.8, 1.0, 4)
        r = bus_flow(a, b, p_aux=rng.uniform(0, 2000), eta_buck=e[0], eta_dab=e[1], eta_pol=e[2], eta_bd=e[3], p_chg_max_bus=rng.choice([float("inf"), 3000.0, 0.0]))
        worst = max(worst, abs(r["residuo_W"]))
        assert r["curtail_array"] >= -1e-9 and r["q_buck"] >= -1e-9 and r["q_bd"] >= -1e-9 and r["q_dab"] >= -1e-9 and r["q_pol"] >= -1e-9
    assert worst < 1e-8, worst


def test_eficiencias_unitarias_sem_perdas():
    r = bus_flow(50000.0, 30000.0, p_aux=1000.0, eta_buck=1.0, eta_dab=1.0, eta_pol=1.0, eta_bd=1.0)
    assert abs(r["q_buck"]) + abs(r["q_bd"]) + abs(r["q_dab"]) + abs(r["q_pol"]) < 1e-9
    assert abs(r["s_bat_term"] - 19000.0) < 1e-9


def test_cadeia_de_carga_composta():
    c = load_chain(10000.0, 0.95, 0.90)
    assert abs(10000.0 / c["p_bus"] - 0.95 * 0.90) < 1e-12
    assert abs(c["q_pol"] + c["q_dab"] - (c["p_bus"] - 10000.0)) < 1e-9


def test_descarga_no_eclipse():
    r = bus_flow(0.0, 18000.0)
    deficit = r["l_bus"]
    assert abs(r["s_bat_term"] + deficit / ETA_BD) < 1e-9 and abs(r["q_bd"] - (deficit / ETA_BD - deficit)) < 1e-9 and r["q_buck"] == 0.0


def test_limite_de_carga_gera_curtailment():
    r = bus_flow(50000.0, 20000.0, p_chg_max_bus=1000.0)
    assert abs(r["s_bat_term"] - ETA_BD * 1000.0) < 1e-9 and r["curtail_array"] > 0.0 and abs(r["residuo_W"]) < 1e-9


def test_fronteira_termica_fecha():
    r = bus_flow(50000.0, 30000.0)
    used_minus_s = 50000.0 - r["curtail_array"] - r["s_bat_term"]
    assert abs(r["q_node"] + r["opt_aux"] - used_minus_s) < 1e-8
    assert abs(r["heat_aux"] + r["opt_aux"] - P_AUX_W) < 1e-9


def test_orbita_fecha_periodicamente():
    o = orbit_budget(1105.39, 5926.379)
    ps = o["P_TI_sun_sustentavel_W"]
    chk = orbit_budget(1105.39, 5926.379, p_ti_sun=ps)
    assert abs(chk["saldo_energia_orbita_Wh"]) < 1e-6


def test_monotonia_em_eficiencia():
    a = orbit_budget(1105.39, 5926.379, eta_pol=0.85)["P_TI_sun_sustentavel_W"]
    b = orbit_budget(1105.39, 5926.379, eta_pol=0.95)["P_TI_sun_sustentavel_W"]
    assert b > a
    assert orbit_budget(2023.33, 5926.379)["P_TI_sun_sustentavel_W"] < orbit_budget(1105.39, 5926.379)["P_TI_sun_sustentavel_W"]


def test_dominio():
    for bad in (lambda: load_chain(-1.0), lambda: bus_flow(1.0, 1.0, eta_buck=0.0), lambda: orbit_budget(6000.0, 5926.0), lambda: bus_flow(float("nan"), 1.0)):
        try:
            bad()
        except ValueError:
            continue
        raise AssertionError("deveria recusar entrada invalida")


def test_reprodutibilidade():
    assert orbit_budget(1105.39, 5926.379) == orbit_budget(1105.39, 5926.379)


TESTS = [test_fechamento_aleatorio, test_eficiencias_unitarias_sem_perdas, test_cadeia_de_carga_composta, test_descarga_no_eclipse, test_limite_de_carga_gera_curtailment,
         test_fronteira_termica_fecha, test_orbita_fecha_periodicamente, test_monotonia_em_eficiencia, test_dominio, test_reprodutibilidade]


def run_selftests():
    res = {}
    for f in TESTS:
        try:
            f(); res[f.__name__] = "PASS"
        except AssertionError as e:
            res[f.__name__] = "FAIL: " + str(e)[:200]
        except Exception as e:
            res[f.__name__] = "ERROR " + type(e).__name__ + ": " + str(e)[:200]
    return res


def run_report():
    T = 5926.379071134441                      # s, saida local do usuario (modelo_orbital, 700 km)
    out = {"parametros": {"eta_buck": ETA_BUCK, "eta_dab": ETA_DAB, "eta_pol": ETA_POL, "eta_bd": ETA_BD, "eta_cadeia_carga": ETA_DAB * ETA_POL, "P_aux_W": P_AUX_W, "P_pv_out_W": PV_OUT_DESIGN_W}}
    out["fluxo_sol_40kW_TI_[W]"] = bus_flow(PV_OUT_DESIGN_W, 40000.0)
    out["fluxo_sol_30kW_TI_[W]"] = bus_flow(PV_OUT_DESIGN_W, 30000.0)
    out["fluxo_eclipse_18kW_TI_[W]"] = bus_flow(0.0, 18000.0)
    out["TI_maximo_sol_sem_carregar_BESS_W"] = (ETA_BUCK * PV_OUT_DESIGN_W - P_AUX_W) * ETA_DAB * ETA_POL
    rows = []
    for nome, te in (("LTAN18_projeto_Marco1_[nao_verificado_local]", 1105.39), ("LTAN10.5_2026-06-01_[verificado_local_oblato]", 2023.3335273230243)):
        o = orbit_budget(te, T, e_nom_j=9.0e3 * 3600.0)
        rows.append([nome, te, o["L_ecl_bus_W"], o["E_term_ecl_Wh"], o["E_nom_req_Wh"], o["P_TI_sun_sustentavel_W"], o["P_TI_ecl_max_para_E_nom_W"]])
        o2 = orbit_budget(te, T, p_ti_sun=40000.0)
        rows[-1].append(o2["saldo_energia_orbita_Wh"])
    out["orbita_[caso | t_ecl_s | L_ecl_bus_W | E_term_ecl_Wh | E_nom_req_Wh | P_TI_sol_sustentavel_W | P_TI_ecl_max_com_9kWh_W | saldo_Wh_se_TI_sol=40kW]"] = rows
    base = orbit_budget(1105.39, T)
    sens = []
    for nome, kw in (("eta_buck 0.88", {"eta_buck": 0.88}), ("eta_buck 0.94", {"eta_buck": 0.94}), ("eta_dab 0.95", {"eta_dab": 0.95}), ("eta_dab 0.98", {"eta_dab": 0.98}),
                     ("eta_pol 0.85", {"eta_pol": 0.85}), ("eta_pol 0.95", {"eta_pol": 0.95}), ("eta_bd 0.95", {"eta_bd": 0.95}), ("eta_bd 0.99", {"eta_bd": 0.99}),
                     ("P_pv -5%", {"p_pv_out": 0.95 * PV_OUT_DESIGN_W}), ("P_pv +5%", {"p_pv_out": 1.05 * PV_OUT_DESIGN_W}), ("P_aux 500 W", {"p_aux": 500.0}),
                     ("P_aux 1500 W", {"p_aux": 1500.0}), ("P_TI_ecl 10 kW", {"p_ti_ecl": 10000.0}), ("P_TI_ecl 25 kW", {"p_ti_ecl": 25000.0})):
        o = orbit_budget(1105.39, T, **kw)
        sens.append([nome, o["P_TI_sun_sustentavel_W"], o["P_TI_sun_sustentavel_W"] - base["P_TI_sun_sustentavel_W"], o["E_nom_req_Wh"], o["E_nom_req_Wh"] - base["E_nom_req_Wh"]])
    out["base_LTAN18"] = {"P_TI_sun_sustentavel_W": base["P_TI_sun_sustentavel_W"], "E_nom_req_Wh": base["E_nom_req_Wh"]}
    out["sensibilidade_OAT_[caso | P_TI_sol_W | dP_W | E_nom_req_Wh | dE_Wh]"] = sens
    # comparacao com o legado 3_1 (eta_conv=0.95 so na TI, P_pv_eol=49733 W ja pos-Buck; eclipse 18 kW TI)
    out["legado_3_1"] = {"L_ecl_bus_W": (18000.0 + 1000.0) / 0.95, "novo_L_ecl_bus_W": orbit_budget(1105.39, T)["L_ecl_bus_W"],
                         "TI_sol_max_sem_carga_W": 49733.0 * 0.95 - 1000.0}
    return out


if __name__ == "__main__":
    payload = {"testes": run_selftests(), "relatorio": run_report()}
    print(json.dumps(payload, indent=1, default=float, ensure_ascii=False))
    result = {"json": json.dumps(payload, default=float, ensure_ascii=False)}
