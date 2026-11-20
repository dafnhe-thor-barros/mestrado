"""modelo_bess.py -- Marco 2c: BESS em nivel macro de energia (d SoC/dt, rendimentos, janela de SoC, limite de taxa C, dimensionamento por eclipse).
Celula de referencia: Saft VES 16 [FONTE: datasheet Doc. n 33026-2-0926, lido em texto e imagem; celula vendida so em modulos/baterias].
Etiquetas: [FONTE] | [CALCULADO] | [HIPOTESE DECLARADA] (com sensibilidade) | [PENDENTE] sem fonte | [LEGADO] valor do legado, sem fonte.
Fronteira (sem dupla contagem): os rendimentos do conversor bidirecional (eta_bd) e dos estagios 400->48->0.8 V pertencem ao Marco 2b e sao
repetidos aqui como parametros (mesmos valores); eta_ch / eta_dis sao SO da bateria (terminais -> energia quimica). Nao ha I2R separado.
NAO modela: OCV(SoC), tensao terminal, resistencia interna, circuito DP, AUKF, temperatura/aquecedor (Marco 3), penumbra (PV=0 no eclipse).
Nao validado experimentalmente. Autocontido (numpy). Se modelo_barramento.py estiver na pasta, run_selftests() faz a conferencia cruzada.
"""
import math
import json
import numpy as np

# ---------------------------------------------------------------- celula VES 16 [FONTE: datasheet Saft, lido]
V_CELL_NOM = 3.6            # V, tensao nominal
C_CELL_AH = 4.5             # Ah, C/1.5, 4.1 V/3 V, 20 C
V_CELL_MAX = 4.1            # V, tensao maxima de carga
V_CELL_MIN = 2.7            # V, limite inferior de descarga continua (0 a +45 C)
M_CELL_KG = 0.115           # kg, massa maxima
SPEC_E_MIN_WH_KG = 140.0    # Wh/kg, energia especifica minima
E_DENS_WH_L = 314.0         # Wh/l
I_DIS_MAX_A = 20.0          # A, "Continuous ~2 s pulse 20 A": texto ambiguo (continua e/ou pulso de ~2 s) -> checagem sob a leitura mais restritiva
LEO_CYCLES_AT_40DOD = 60000  # ciclos LEO a 40% DOD (declaracao do fabricante)
CHG_C_RATE_LEO = 1.0 / 2.5  # corrente continua de carga recomendada a 20 C, LEO (20% DOD): C/2.5
T_CHG_C = (10.0, 35.0)      # temperatura de carga [C]
T_DIS_C = (0.0, 40.0)       # temperatura de descarga [C]
E_CELL_WH = V_CELL_NOM * C_CELL_AH   # [CALCULADO] 16.2 Wh

# ---------------------------------------------------------------- hipoteses de projeto
SOC_MIN, SOC_MAX = 0.70, 0.95   # janela util [decisao do usuario]
FADE_EOL = 0.20     # [HIPOTESE DECLARADA] perda de capacidade no EOL (sem fonte lida); sensibilidade 0.10..0.30
N_SER = 108         # [HIPOTESE DECLARADA] celulas em serie (herdado do legado 108S): V_pack = 388.8 V nominal, 291.6..442.8 V
MISSION_YEARS = 10.0  # [LEGADO] vida de projeto, sem fonte
Q_AUX_FRAC = 903.2 / 1000.0   # [LEGADO 4_1] fracao de P_aux que vira calor interno
PAR = {"eta_buck": 0.91, "eta_dab": 0.9655, "eta_pol": 0.90, "eta_bd": 0.97, "p_aux": 1000.0, "p_pv_out": 50000.0,   # valores do Marco 2b (hipoteses)
       "eta_ch": 0.995, "eta_dis": 0.97}   # [PENDENTE] valores do legado 3_1/3_3 (sem fonte; datasheet nao informa eficiencia)


def _par(par=None):
    p = dict(PAR)
    if par:
        p.update(par)
    for k in ("eta_buck", "eta_dab", "eta_pol", "eta_bd", "eta_ch", "eta_dis"):
        if not (0.0 < p[k] <= 1.0):
            raise ValueError("%s deve estar em (0, 1] (recebido %r)" % (k, p[k]))
    for k in ("p_aux", "p_pv_out"):
        if not (math.isfinite(p[k]) and p[k] >= 0.0):
            raise ValueError("%s deve ser finito e >= 0 (recebido %r)" % (k, p[k]))
    return p


def mission_cycles(T, years=MISSION_YEARS):
    """N = anos x 365.25 x 86400 / T  [CALCULADO]"""
    return years * 365.25 * 86400.0 / T


def analytic_orbit(t_ecl, T, p_ti_ecl, par=None):
    """Regime periodico. Eclipse: L_ecl = P_TI/(eta_dab eta_pol) + P_aux; s_out = L_ecl/eta_bd; E_quim = s_out t_ecl/eta_dis.
    Sol: excedente no barramento S = L_ecl t_ecl / (eta_ch eta_dis eta_bd^2 t_sun)  (carga = descarga);
    L_sun = eta_buck P_pv - S ;  P_TI_sun = (L_sun - P_aux) eta_dab eta_pol.  Com eta_ch = eta_dis = 1 reduz-se ao Marco 2b. [CALCULADO]"""
    p = _par(par)
    if not (0.0 <= t_ecl < T) or not (math.isfinite(p_ti_ecl) and p_ti_ecl >= 0.0):
        raise ValueError("t_ecl deve estar em [0, T) e p_ti_ecl >= 0")
    t_sun = T - t_ecl
    chain = p["eta_dab"] * p["eta_pol"]
    l_ecl = p_ti_ecl / chain + p["p_aux"]
    s_out = l_ecl / p["eta_bd"]
    e_chem_out = s_out * t_ecl / p["eta_dis"] / 3600.0
    surplus = l_ecl * t_ecl / (p["eta_ch"] * p["eta_dis"] * p["eta_bd"] ** 2 * t_sun)
    l_sun = p["eta_buck"] * p["p_pv_out"] - surplus
    return {"L_ecl_bus_W": l_ecl, "s_out_W": s_out, "E_chem_out_Wh": e_chem_out, "S_bus_W": surplus,
            "P_TI_sun_sust_W": max(0.0, (l_sun - p["p_aux"]) * chain), "viavel_sol": bool(l_sun >= p["p_aux"])}


def size_bess(t_ecl, T, p_ti_ecl, par=None, fade=FADE_EOL, n_ser=N_SER, continuous=False):
    """Capacidade no EOL tal que a energia quimica do eclipse ocupe a janela: E_cap_EOL = E_quim/(SOC_MAX-SOC_MIN);
    BOL = EOL/(1-fade); n_par = ceil(E_BOL/(n_ser E_cel)).  Massa = so celulas (massa maxima do datasheet), sem modulo/estrutura."""
    if not (0.0 <= fade < 1.0):
        raise ValueError("fade deve estar em [0, 1)")
    a = analytic_orbit(t_ecl, T, p_ti_ecl, par)
    e_eol = a["E_chem_out_Wh"] / (SOC_MAX - SOC_MIN)
    e_bol = e_eol / (1.0 - fade)
    n_par_real = e_bol / (n_ser * E_CELL_WH)
    n_par = n_par_real if continuous else max(1, math.ceil(n_par_real - 1e-12))
    n_cells = n_ser * n_par
    return {"n_ser": n_ser, "n_par": n_par, "n_cells": n_cells, "fade": fade, "E_req_EOL_Wh": e_eol, "E_req_BOL_Wh": e_bol,
            "E_inst_BOL_Wh": n_cells * E_CELL_WH, "mass_cells_kg": n_cells * M_CELL_KG,
            "dod_eol": a["E_chem_out_Wh"] / (n_cells * E_CELL_WH * (1.0 - fade)), "dod_bol": a["E_chem_out_Wh"] / (n_cells * E_CELL_WH),
            "P_TI_sun_sust_W": a["P_TI_sun_sust_W"], "L_ecl_bus_W": a["L_ecl_bus_W"]}


def max_ti_ecl(t_ecl, n_cells, par=None, fade=FADE_EOL):
    """TI maxima no eclipse para um BESS dado (janela inteira usada). Retorna (P_TI_W, aux_atendido)."""
    p = _par(par)
    e_chem = (SOC_MAX - SOC_MIN) * n_cells * E_CELL_WH * (1.0 - fade)
    l = e_chem * 3600.0 * p["eta_dis"] * p["eta_bd"] / t_ecl
    return max(0.0, (l - p["p_aux"]) * p["eta_dab"] * p["eta_pol"]), bool(l >= p["p_aux"])


def simulate(bess, t_ecl, T, n_orb=10, dt=10.0, p_ti_sun_req=40000.0, p_ti_ecl_req=18000.0, par=None, soc0=SOC_MAX,
             life="eol", c_rate_chg=CHG_C_RATE_LEO, return_records=False):
    """d E/dt = eta_ch s_chg - s_dis/eta_dis  [W], E em Wh, SoC = E/E_cap. Eclipse: PV = 0 (sem penumbra). Politica de despacho (NAO e TEMS, nao otimiza nada):
    atende a TI pedida; se a energia acima de SOC_MIN nao cobre o restante do eclipse, REDUZ a TI de forma uniforme (racionamento: energia disponivel
    / tempo restante de eclipse, que a geometria orbital fornece). Correcao: a versao gulosa anterior (energia / dt) esgotava a bateria no meio do
    eclipse e deixava P_aux sem energia -- o docstring dizia 'P_aux tem prioridade', o que era falso. P_aux so falha se a bateria nao cobre nem P_aux.
    Carga limitada por headroom ate SOC_MAX e por c_rate_chg; o que sobra e curtailment do arranjo. Discretizacao: fase do eclipse por ponto medio do passo."""
    p = _par(par)
    if not (dt > 0.0 and T > 0.0 and 0.0 <= t_ecl < T and n_orb >= 1):
        raise ValueError("dt, T > 0; t_ecl em [0, T); n_orb >= 1")
    if not (SOC_MIN - 1e-12 <= soc0 <= SOC_MAX + 1e-12):
        raise ValueError("soc0 fora da janela")
    for nm, v in (("p_ti_sun_req", p_ti_sun_req), ("p_ti_ecl_req", p_ti_ecl_req)):
        if not (math.isfinite(v) and v >= 0.0):
            raise ValueError("%s deve ser finito e >= 0" % nm)
    n_ser, n_par = bess["n_ser"], bess["n_par"]
    f_life = (1.0 - bess["fade"]) if life == "eol" else 1.0
    e_cap = n_ser * n_par * E_CELL_WH * f_life
    e_min, e_max = SOC_MIN * e_cap, SOC_MAX * e_cap
    v_pack = n_ser * V_CELL_NOM
    p_chg_term_max = n_par * C_CELL_AH * c_rate_chg * v_pack
    chain = p["eta_dab"] * p["eta_pol"]
    n_steps = int(round(n_orb * T / dt))
    spo = int(round(T / dt))
    keys = ("ti", "ti_req", "aux_unmet", "curtail", "q_batt", "q_bd", "i_cell", "resid", "ecl")
    R = {k: np.zeros(n_steps) for k in keys}
    soc = np.zeros(n_steps + 1)
    E = soc0 * e_cap
    soc[0] = soc0
    for k in range(n_steps):
        ecl = ((k + 0.5) * dt) % T < t_ecl
        pv = 0.0 if ecl else p["p_pv_out"]
        ti_req = p_ti_ecl_req if ecl else p_ti_sun_req
        pv_bus = p["eta_buck"] * pv
        l_req = ti_req / chain + p["p_aux"]
        net = pv_bus - l_req
        ti, aux = ti_req, p["p_aux"]
        chem_in = chem_out = 0.0
        if net >= 0.0:
            head = max(0.0, e_max - E) * 3600.0 / dt
            to_bat = min(net, head / (p["eta_ch"] * p["eta_bd"]), p_chg_term_max / p["eta_bd"])
            s_term = p["eta_bd"] * to_bat
            chem_in = p["eta_ch"] * s_term
            used_bus = l_req + to_bat
            q_bd, q_b, i_sign = to_bat - s_term, s_term - chem_in, +1.0
        else:
            t_phase = (k * dt) % T
            t_rem = max(t_ecl - t_phase, dt) if ecl else dt     # tempo restante de eclipse (geometria orbital, deterministica)
            avail = max(0.0, E - e_min) * 3600.0 / t_rem        # racionamento: energia acima do piso espalhada no restante do eclipse
            deficit = -net
            deficit_max = avail * p["eta_dis"] * p["eta_bd"]
            if deficit > deficit_max:
                deficit = deficit_max
                l_allowed = pv_bus + deficit
                if l_allowed >= p["p_aux"]:
                    ti = (l_allowed - p["p_aux"]) * chain
                else:
                    ti, aux = 0.0, l_allowed
            s_term = deficit / p["eta_bd"]
            chem_out = s_term / p["eta_dis"]
            used_bus = pv_bus
            q_bd, q_b, i_sign = s_term - deficit, chem_out - s_term, -1.0
        p48 = ti / p["eta_pol"]
        q_pol, q_dab = p48 - ti, p48 / p["eta_dab"] - p48
        used_out = used_bus / p["eta_buck"]
        q_buck = used_out - used_bus
        d_e = chem_in - chem_out
        R["resid"][k] = used_out - (q_buck + ti + q_pol + q_dab + aux + q_bd + q_b + d_e)
        E += d_e * dt / 3600.0
        soc[k + 1] = E / e_cap
        R["ti"][k], R["ti_req"][k], R["aux_unmet"][k] = ti, ti_req, p["p_aux"] - aux
        R["curtail"][k], R["q_batt"][k], R["q_bd"][k] = pv - used_out, q_b, q_bd
        R["i_cell"][k], R["ecl"][k] = i_sign * s_term / v_pack / n_par, 1.0 if ecl else 0.0
    L = slice(max(0, n_steps - spo), n_steps)
    m = R["ecl"][L] > 0.5
    h = dt / 3600.0
    out = {"e_cap_Wh": e_cap, "soc_min_run": float(soc.min()), "soc_max_run": float(soc.max()),
           "soc_min_last": float(soc[L.start:L.stop + 1].min()), "soc_max_last": float(soc[L.start:L.stop + 1].max()),
           "ti_ecl_mean_W": float(R["ti"][L][m].mean()) if m.any() else float("nan"),
           "ti_sun_mean_W": float(R["ti"][L][~m].mean()) if (~m).any() else float("nan"),
           "ti_orbit_mean_W": float(R["ti"][L].mean()),
           "shed_Wh_orbit": float((R["ti_req"][L] - R["ti"][L]).sum() * h), "aux_unmet_Wh_orbit": float(R["aux_unmet"][L].sum() * h),
           "curtail_Wh_orbit": float(R["curtail"][L].sum() * h), "q_batt_mean_W": float(R["q_batt"][L].mean()), "q_batt_peak_W": float(R["q_batt"][L].max()),
           "q_bd_mean_W": float(R["q_bd"][L].mean()),
           "i_cell_dis_max_A": float(max(0.0, -R["i_cell"].min())), "i_cell_chg_max_A": float(max(0.0, R["i_cell"].max())),
           "dE_last_orbit_Wh": float((soc[L.stop] - soc[L.start]) * e_cap), "ledger_resid_max_W": float(np.abs(R["resid"]).max())}
    if return_records:
        out["_rec"], out["_soc"] = R, soc
    return out


# ---------------------------------------------------------------- testes
def test_celula_consistente_datasheet():
    assert abs(E_CELL_WH - 16.2) < 1e-12
    assert abs(SPEC_E_MIN_WH_KG * M_CELL_KG - E_CELL_WH) / E_CELL_WH < 0.02     # 140 Wh/kg x 0.115 kg = 16.1 Wh
    assert V_CELL_MIN * N_SER < 400.0 < V_CELL_MAX * N_SER                      # barramento de 400 V dentro da faixa do pack 108S


def test_regressao_contra_saida_local_do_2b():
    """Com eta_ch = eta_dis = 1 e fade = 0 deve reproduzir os numeros do Marco 2b CONFIRMADOS LOCALMENTE pelo usuario."""
    T, te = 5926.379071134441, 1105.39
    par = {"eta_ch": 1.0, "eta_dis": 1.0}
    a = analytic_orbit(te, T, 18000.0, par)
    s = size_bess(te, T, 18000.0, par, fade=0.0, continuous=True)
    assert abs(a["P_TI_sun_sust_W"] - 34070.12) < 0.05, a["P_TI_sun_sust_W"]
    assert abs(s["E_req_EOL_Wh"] - 27495.03) < 0.05, s["E_req_EOL_Wh"]


def test_fechamento_periodico_simulacao_x_analitico():
    T, te, dt = 6000.0, 1200.0, 10.0     # t_ecl e T multiplos de dt: sem erro de discretizacao de fase
    s = size_bess(te, T, 18000.0, continuous=True)
    a = analytic_orbit(te, T, 18000.0)
    r = simulate(s, te, T, n_orb=6, dt=dt, p_ti_sun_req=a["P_TI_sun_sust_W"], p_ti_ecl_req=18000.0)
    assert abs(r["soc_min_last"] - SOC_MIN) < 1e-6 and abs(r["soc_max_last"] - SOC_MAX) < 1e-6, (r["soc_min_last"], r["soc_max_last"])
    assert abs(r["dE_last_orbit_Wh"]) < 1e-6 and r["shed_Wh_orbit"] < 1e-6 and r["aux_unmet_Wh_orbit"] == 0.0


def test_subdimensionado_reduz_ti_e_respeita_piso():
    T, te = 6000.0, 1200.0
    s = size_bess(te, T, 18000.0, continuous=True)
    half = dict(s, n_par=0.5 * s["n_par"])
    r = simulate(half, te, T, n_orb=6, dt=10.0, p_ti_sun_req=20000.0, p_ti_ecl_req=18000.0)
    assert r["soc_min_run"] >= SOC_MIN - 1e-9, ("piso", r["soc_min_run"])
    assert r["shed_Wh_orbit"] > 0.0, ("shed", r["shed_Wh_orbit"])
    assert r["ti_ecl_mean_W"] < 18000.0, ("ti_ecl", r["ti_ecl_mean_W"])
    assert r["aux_unmet_Wh_orbit"] == 0.0, ("aux_meia_bateria", r["aux_unmet_Wh_orbit"])
    tiny = dict(s, n_par=0.02 * s["n_par"])
    r2 = simulate(tiny, te, T, n_orb=3, dt=10.0, p_ti_sun_req=20000.0, p_ti_ecl_req=18000.0)
    assert r2["aux_unmet_Wh_orbit"] > 0.0, ("aux_minuscula", r2["aux_unmet_Wh_orbit"])
    assert r2["soc_min_run"] >= SOC_MIN - 1e-9, ("piso_minuscula", r2["soc_min_run"])


def test_invariantes_aleatorios():
    rng = np.random.default_rng(2026)
    for _ in range(25):
        T = rng.uniform(5000.0, 7000.0)
        te = rng.uniform(0.0, 0.45 * T)
        par = {"eta_buck": rng.uniform(0.85, 1.0), "eta_dab": rng.uniform(0.9, 1.0), "eta_pol": rng.uniform(0.85, 1.0), "eta_bd": rng.uniform(0.9, 1.0),
               "eta_ch": rng.uniform(0.95, 1.0), "eta_dis": rng.uniform(0.9, 1.0), "p_aux": rng.uniform(0.0, 2000.0)}
        b = {"n_ser": 108, "n_par": float(rng.integers(2, 40)), "fade": 0.2}
        r = simulate(b, te, T, n_orb=3, dt=25.0, p_ti_sun_req=rng.uniform(0, 45e3), p_ti_ecl_req=rng.uniform(0, 30e3), par=par,
                     soc0=rng.uniform(SOC_MIN, SOC_MAX), return_records=True)
        assert r["soc_min_run"] >= SOC_MIN - 1e-9 and r["soc_max_run"] <= SOC_MAX + 1e-9
        assert r["ledger_resid_max_W"] < 1e-6 and r["shed_Wh_orbit"] >= -1e-9 and r["aux_unmet_Wh_orbit"] >= -1e-9 and r["curtail_Wh_orbit"] >= -1e-9
        assert (r["_rec"]["q_batt"] >= -1e-9).all() and (r["_rec"]["q_bd"] >= -1e-9).all()


def test_convergencia_em_dt():
    T, te = 5926.379071134441, 1105.3908
    s = size_bess(te, T, 18000.0, continuous=True)
    a = analytic_orbit(te, T, 18000.0)
    ref = simulate(s, te, T, n_orb=3, dt=1.0, p_ti_sun_req=a["P_TI_sun_sust_W"])["soc_min_last"]
    err = {dt: abs(simulate(s, te, T, n_orb=3, dt=dt, p_ti_sun_req=a["P_TI_sun_sust_W"])["soc_min_last"] - ref) for dt in (40.0, 20.0, 10.0, 5.0)}
    assert err[5.0] < 2e-3 and err[5.0] <= err[40.0] + 1e-12, err


def test_vida_em_ciclos_cobre_a_missao():
    n = mission_cycles(5926.379071134441)
    assert n < LEO_CYCLES_AT_40DOD and 5.0e4 < n < 5.6e4, n


def test_monotonia_e_dominio():
    T, te = 5926.379071134441, 1105.3908
    assert size_bess(te, T, 18000.0, {"eta_dis": 0.95})["n_cells"] >= size_bess(te, T, 18000.0, {"eta_dis": 0.99})["n_cells"]
    assert size_bess(te, T, 25000.0)["n_cells"] > size_bess(te, T, 10000.0)["n_cells"]
    assert size_bess(te, T, 18000.0, fade=0.3)["n_cells"] >= size_bess(te, T, 18000.0, fade=0.1)["n_cells"]
    for bad in (lambda: _par({"eta_ch": 0.0}), lambda: analytic_orbit(7000.0, 5926.0, 1.0), lambda: simulate({"n_ser": 108, "n_par": 5, "fade": 0.2}, 1000.0, 5926.0, soc0=0.5),
                lambda: size_bess(1000.0, 5926.0, -1.0), lambda: size_bess(1000.0, 5926.0, 1.0, fade=1.0)):
        try:
            bad()
        except ValueError:
            continue
        raise AssertionError("deveria recusar entrada invalida")


def test_reprodutibilidade():
    s = size_bess(1105.3908, 5926.379071134441, 18000.0)
    assert simulate(s, 1105.3908, 5926.379071134441, n_orb=2) == simulate(s, 1105.3908, 5926.379071134441, n_orb=2)


def _cross_check():
    """Retorna 'PASS' ou 'SKIPPED (...)'. Fora do pytest porque funcoes de teste NAO devem retornar valor (pytest recente pode tratar como falha)."""
    try:
        import modelo_barramento as mb
    except ImportError:
        return "SKIPPED (modelo_barramento.py ausente)"
    T, te = 5926.379071134441, 1105.3908
    o = mb.orbit_budget(te, T)
    a = analytic_orbit(te, T, 18000.0, {"eta_ch": 1.0, "eta_dis": 1.0})
    assert abs(o["P_TI_sun_sustentavel_W"] - a["P_TI_sun_sust_W"]) < 1e-6
    f = mb.bus_flow(0.0, 18000.0)
    assert abs(-f["s_bat_term"] - a["s_out_W"]) < 1e-6
    return "PASS"


def test_conferencia_cruzada_com_modelo_barramento():
    assert _cross_check().split()[0] in ("PASS", "SKIPPED")


TESTS = [test_celula_consistente_datasheet, test_regressao_contra_saida_local_do_2b, test_fechamento_periodico_simulacao_x_analitico,
         test_subdimensionado_reduz_ti_e_respeita_piso, test_invariantes_aleatorios, test_convergencia_em_dt, test_vida_em_ciclos_cobre_a_missao,
         test_monotonia_e_dominio, test_reprodutibilidade]


def run_selftests():
    res = {}
    for f in TESTS:
        try:
            f()
            res[f.__name__] = "PASS"
        except AssertionError as e:
            res[f.__name__] = "FAIL: " + str(e)[:200]
        except Exception as e:
            res[f.__name__] = "ERROR " + type(e).__name__ + ": " + str(e)[:200]
    try:
        res["test_conferencia_cruzada_com_modelo_barramento"] = _cross_check()
    except AssertionError as e:
        res["test_conferencia_cruzada_com_modelo_barramento"] = "FAIL: " + str(e)[:200]
    except Exception as e:
        res["test_conferencia_cruzada_com_modelo_barramento"] = "ERROR " + type(e).__name__ + ": " + str(e)[:200]
    return res


def _brief(r):
    keys = ("soc_min_last", "soc_max_last", "ti_ecl_mean_W", "ti_sun_mean_W", "ti_orbit_mean_W", "shed_Wh_orbit", "aux_unmet_Wh_orbit",
            "curtail_Wh_orbit", "q_batt_mean_W", "q_batt_peak_W", "i_cell_dis_max_A", "i_cell_chg_max_A", "dE_last_orbit_Wh")
    return {k: r[k] for k in keys}


def run_report():
    T = 5926.379071134441
    te18, te105 = 1105.3908, 2023.3335273230243     # eclipses confirmados localmente pelo usuario (LTAN 18:00 em 21/06/2026; LTAN 10:30 em 01/06/2026)
    out = {"celula_VES16": {"E_cel_Wh": E_CELL_WH, "E_via_energia_especifica_Wh": SPEC_E_MIN_WH_KG * M_CELL_KG, "V_pack_nom_V": N_SER * V_CELL_NOM,
                            "V_pack_faixa_V": [N_SER * V_CELL_MIN, N_SER * V_CELL_MAX], "ciclos_missao": mission_cycles(T), "ciclos_declarados_40DOD": LEO_CYCLES_AT_40DOD},
           "legado_massa": {"9kWh/140Wh_kg_so_celulas_kg": 9000.0 / SPEC_E_MIN_WH_KG, "legado_C_BAT_comentario_kg": 60.0}}
    rows = []
    for nome, te in (("LTAN18", te18), ("LTAN10.5", te105)):
        for ti in (0.0, 5000.0, 10000.0, 18000.0, 25000.0):
            s = size_bess(te, T, ti)
            rows.append([nome, te, ti, s["n_par"], s["n_cells"], s["mass_cells_kg"], s["E_inst_BOL_Wh"], s["E_req_EOL_Wh"], s["dod_eol"], s["P_TI_sun_sust_W"]])
    out["dimensionamento_[caso | t_ecl_s | P_TI_ecl_W | n_par | n_celulas | massa_celulas_kg | E_BOL_inst_Wh | E_req_EOL_Wh | DoD_EOL | P_TI_sol_sust_W]"] = rows
    inv = []
    for nome, te in (("LTAN18", te18), ("LTAN10.5", te105)):
        for n_par in (5, 10, 20):
            n = N_SER * n_par
            ti, ok = max_ti_ecl(te, n)
            inv.append([nome, n_par, n, n * M_CELL_KG, ti, ok])
    out["inverso_[caso | n_par | n_celulas | massa_kg | P_TI_ecl_max_W | aux_atendido]"] = inv
    base = size_bess(te18, T, 18000.0)
    sens = []
    for nome, kw in (("fade 0.10", {"fade": 0.10}), ("fade 0.30", {"fade": 0.30}), ("eta_dis 0.99", {"par": {"eta_dis": 0.99}}), ("eta_dis 0.95", {"par": {"eta_dis": 0.95}}),
                     ("eta_ch 0.99", {"par": {"eta_ch": 0.99}}), ("eta_bd 0.95", {"par": {"eta_bd": 0.95}}), ("eta_bd 0.99", {"par": {"eta_bd": 0.99}}),
                     ("eta_dab 0.95", {"par": {"eta_dab": 0.95}}), ("eta_pol 0.85", {"par": {"eta_pol": 0.85}}), ("eta_pol 0.95", {"par": {"eta_pol": 0.95}}),
                     ("P_aux 500 W", {"par": {"p_aux": 500.0}}), ("P_aux 1500 W", {"par": {"p_aux": 1500.0}})):
        s = size_bess(te18, T, 18000.0, **kw)
        sens.append([nome, s["n_cells"], s["n_cells"] - base["n_cells"], s["mass_cells_kg"], s["E_req_EOL_Wh"], s["E_req_EOL_Wh"] - base["E_req_EOL_Wh"]])
    out["base_LTAN18_TIecl18kW"] = {k: base[k] for k in ("n_par", "n_cells", "mass_cells_kg", "E_inst_BOL_Wh", "E_req_EOL_Wh", "E_req_BOL_Wh", "dod_eol", "dod_bol", "P_TI_sun_sust_W")}
    out["sensibilidade_OAT_[caso | n_celulas | dN | massa_kg | E_req_EOL_Wh | dE_Wh]"] = sens
    a = analytic_orbit(te18, T, 18000.0)
    ps = a["P_TI_sun_sust_W"]
    out["sim_A_dimensionado_TIsol=sustentavel"] = _brief(simulate(base, te18, T, n_orb=12, dt=10.0, p_ti_sun_req=ps))
    out["sim_B_mesmo_BESS_TIsol=40kW_pedido"] = _brief(simulate(base, te18, T, n_orb=12, dt=10.0, p_ti_sun_req=40000.0))
    leg = {"n_ser": N_SER, "n_par": 5, "fade": FADE_EOL}
    out["sim_C_5P_(~9kWh)_TIecl=18kW_pedido_TIsol=sust"] = _brief(simulate(leg, te18, T, n_orb=12, dt=10.0, p_ti_sun_req=ps))
    out["sim_D_BESS_do_caso_A_em_LTAN10.5"] = _brief(simulate(base, te105, T, n_orb=12, dt=10.0, p_ti_sun_req=ps))
    ref = simulate(base, te18, T, n_orb=3, dt=1.0, p_ti_sun_req=ps)["soc_min_last"]
    errs = {dt: abs(simulate(base, te18, T, n_orb=3, dt=dt, p_ti_sun_req=ps)["soc_min_last"] - ref) for dt in (40.0, 20.0, 10.0, 5.0)}
    out["convergencia_dt_|SoCmin(dt)-SoCmin(1s)|"] = {"ref_SoCmin": ref, "erros": errs}
    return out


if __name__ == "__main__":
    payload = {"testes": run_selftests(), "relatorio": run_report()}
    print(json.dumps(payload, indent=1, default=float, ensure_ascii=False))
    result = {"json": json.dumps(payload, default=float, ensure_ascii=False)}
