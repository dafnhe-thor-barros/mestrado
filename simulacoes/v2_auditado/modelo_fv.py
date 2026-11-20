"""modelo_fv.py -- Marco 2a: gerador fotovoltaico AM0 (Spectrolab XTJ 29,5%) em nivel de potencia de MPP.
Etiquetas: [FONTE] valor lido no texto do datasheet colado pelo usuario (Spectrolab XTJ, rev. 5/20/2010, upd. 9/28/12, p.2)
| [CALCULADO] | [HIPOTESE DECLARADA] (com sensibilidade) | [PENDENTE] | [NAO VERIFICADO] (memoria / sem fonte lida).
Escopo: potencia no MPP da celula (T, G, fluencia), dimensionamento de area para P_req na SAIDA DO ARRANJO em EOL (antes do Buck/MPPT),
equilibrio termico simples da celula (diagnostico), sensibilidade OAT. NAO modela curva I-V (3 diodos do legado descartado), Buck, MPPT, hot spot.
Autocontido: nao importa modelo_orbital. Nao validado experimentalmente.
"""
import math
import json
import numpy as np

# ---------------------------------------------------------------- constantes fisicas
K_B, H_P, C_L = 1.380649e-23, 6.62607015e-34, 299792458.0     # [NAO VERIFICADO] valores exatos do SI 2019, de memoria
SIGMA = 2.0 * math.pi ** 5 * K_B ** 4 / (15.0 * H_P ** 3 * C_L ** 2)   # [CALCULADO] Stefan-Boltzmann
T_SPACE = 3.0                       # K, [FONTE] enunciado do projeto (T_space ~ 3 K)
S_NOM = 1361.0                      # W/m2, [FONTE-RESUMO] IAU 2015 B3 (Kopp&Lean 2011: 1360.8 +- 0.5)
R_AU_APH = 1.00014 + 0.01671 - 0.00014   # [NAO VERIFICADO] USNO R=1.00014-0.01671cos g-0.00014cos 2g em g=180 deg
R_AU_PER = 1.00014 - 0.01671 - 0.00014   # idem em g=0

# ---------------------------------------------------------------- datasheet XTJ (texto colado, p.2) [FONTE]
G_REF = 1353.0                      # W/m2 (AM0 135.3 mW/cm2)
JSC, JMP = 17.76e-3, 17.02e-3       # A/cm2, 28 C, celula nua
VOC, VMP = 2.633, 2.348             # V
EFF_MP, EFF_LOAD, CFF = 0.295, 0.293, 0.85
T_REF = 28.0                        # C
RAD_FLU = np.array([1e14, 5e14, 1e15])      # e-/cm2 (1 MeV)
RAD_I = np.array([1.00, 0.99, 0.95])        # Imp/Imp0
RAD_V = np.array([0.94, 0.91, 0.89])        # Vmp/Vmp0
RAD_P = np.array([0.95, 0.90, 0.85])        # Pmp/Pmp0
DJMP = (6.6e-6, 10.0e-6, 13.2e-6)           # A/cm2/C: BOL, 5e14, 1e15 (faixa 15-75 C)
DVMP = (-6.5e-3, -6.8e-3, -6.9e-3)          # V/C
DJSC = (11.6e-6, 10.9e-6, 11.9e-6)
DVOC = (-5.8e-3, -6.5e-3, -6.6e-3)
ALPHA_S = 0.90                      # absortancia solar (vidro 5 mil CMG/AR) [FONTE]
EPS_FRONT = 0.85                    # emitancia normal (Ceria Doped Microsheet) [FONTE]
CELL_SIZES_CM2 = (26.62, 59.65)     # tamanhos padrao [FONTE]

# ---------------------------------------------------------------- hipoteses de projeto
FLUENCE_EOL = 2.5e14   # e-/cm2 1 MeV, [HIPOTESE DECLARADA] = 2.5e13/ano x 10 anos do legado (sem fonte); faixa de sensibilidade 1e14..5e14
UV_LOSS = 0.052        # [HIPOTESE DECLARADA] legado (5.2% em 10 anos, sem fonte); sensibilidade 0..10%
F_OTHER = 1.0          # [HIPOTESE DECLARADA] perdas de montagem/mismatch/cabeamento NAO incluidas; sensibilidade 0.97
PACKING = 0.90         # [HIPOTESE DECLARADA] fator de empacotamento (area de celulas / area do painel); sensibilidade 0.80..0.95
T_ARRAY_HOT = 65.0     # C, [HIPOTESE DECLARADA] legado; sensibilidade 45..85 (alem de 75 C: extrapolacao)
POINTING_DEG = 0.0     # [HIPOTESE DECLARADA] apontamento perfeito ao Sol (decisao do usuario); sensibilidade 2 deg
EPS_BACK_RANGE = (0.30, 0.85)   # [PENDENTE] emitancia da face traseira (faixa so para o diagnostico termico)
G_DESIGN = S_NOM / R_AU_APH ** 2   # W/m2, irradiancia de projeto no AFELIO [CALCULADO]
P_REQ_W = 50000.0      # W na saida do arranjo, EOL, antes do Buck (decisao do usuario)
CELL_DESIGN_CM2 = 59.65


def _log_interp(f, ys):
    f = float(f)
    if not (RAD_FLU[0] <= f <= RAD_FLU[-1]):
        raise ValueError("fluencia fora de [1e14, 1e15]: sem dado no datasheet")
    return float(np.interp(math.log10(f), np.log10(RAD_FLU), ys))


def _lin_coef(f, vals):
    """Coeficiente termico em funcao da fluencia: interpolacao LINEAR entre BOL (0), 5e14 e 1e15. [HIPOTESE DECLARADA]"""
    return float(np.interp(f, [0.0, 5e14, 1e15], list(vals)))


def cell_mpp(T_c, G, fluence=0.0):
    """Ponto de maxima potencia da celula. j = (Jmp + dJ dT) G/Gref (linear em G); v = Vmp + dV dT (sem termo logaritmico em G: [HIPOTESE DECLARADA]);
    degradacao: P = P(T,G) * (Pmp/Pmp0)(fluencia, interp. log-linear [HIPOTESE DECLARADA]); Vmp com razao Vmp/Vmp0; Imp = P/Vmp."""
    dT = T_c - T_REF
    if fluence and fluence > 0.0:
        rp, rv = _log_interp(fluence, RAD_P), _log_interp(fluence, RAD_V)
        dj, dv = _lin_coef(fluence, DJMP), _lin_coef(fluence, DVMP)
    else:
        rp = rv = 1.0
        dj, dv = DJMP[0], DVMP[0]
    j = (JMP + dj * dT) * G / G_REF
    v = VMP + dv * dT
    p = j * v * rp                       # W/cm2
    vmp = v * rv
    return {"pmp_W_cm2": p, "pmp_W_m2": p * 1e4, "vmp_V": vmp, "imp_A_cm2": p / vmp, "eff": p * 1e4 / G,
            "extrapolado_T": not (15.0 <= T_c <= 75.0)}


def size_array(p_req_w=P_REQ_W, T_c=T_ARRAY_HOT, G=G_DESIGN, fluence=FLUENCE_EOL, uv=UV_LOSS, f_other=F_OTHER,
               pf=PACKING, cell_cm2=CELL_DESIGN_CM2, pointing_deg=POINTING_DEG):
    c = cell_mpp(T_c, G, fluence)
    p_m2 = c["pmp_W_m2"] * math.cos(math.radians(pointing_deg)) * (1.0 - uv) * f_other
    a_cells = p_req_w / p_m2
    return {"p_cell_util_W_m2": p_m2, "A_cells_m2": a_cells, "A_panel_m2": a_cells / pf, "n_cells": math.ceil(a_cells / (cell_cm2 * 1e-4)),
            "eff_cell_hot_eol": c["eff"], "vmp_cell_V": c["vmp_V"], "extrapolado_T": c["extrapolado_T"]}


def pv_power_w(a_cells_m2, T_c, G, sunlit_fraction=1.0, fluence=FLUENCE_EOL, uv=UV_LOSS, f_other=F_OTHER, pointing_deg=POINTING_DEG):
    """Interface para o Marco 2b: potencia instantanea na saida do arranjo [W]. sunlit_fraction vem do eclipse do Marco 1."""
    c = cell_mpp(T_c, G, fluence)
    return a_cells_m2 * c["pmp_W_m2"] * math.cos(math.radians(pointing_deg)) * (1.0 - uv) * f_other * sunlit_fraction


def t_eq_c(G, fluence, eps_back, alpha_s=ALPHA_S, eps_front=EPS_FRONT):
    """Equilibrio termico por area de celula (sem substrato, sem IR/albedo da Terra): alpha G - P_el(T) = (eps_f+eps_b) sigma (T^4 - Ts^4). Bissecao."""
    def f(t):
        return alpha_s * G - cell_mpp(t, G, fluence)["pmp_W_m2"] - (eps_front + eps_back) * SIGMA * ((t + 273.15) ** 4 - T_SPACE ** 4)
    lo, hi = -150.0, 250.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    t = 0.5 * (lo + hi)
    return t, f(t)


def run_selftests():
    res, M = {}, {}

    def chk(name, fn):
        try:
            fn()
            res[name] = "PASS"
        except AssertionError as e:
            res[name] = "FAIL: " + str(e)[:200]
        except Exception as e:
            res[name] = "ERROR " + type(e).__name__ + ": " + str(e)[:200]

    def t_sigma():
        M["sigma"] = SIGMA
        assert abs(SIGMA / 5.670374419e-8 - 1.0) < 1e-8

    def t_eff():
        c = cell_mpp(28.0, G_REF)
        M["eff_ref"] = c["eff"]
        assert abs(c["eff"] - EFF_MP) < 0.002

    def t_ff():
        ff = JMP * VMP / (JSC * VOC)
        M["ff_calc"] = ff
        assert abs(ff - CFF) < 0.01

    def t_load():
        e = 17.14e-3 * 2.310 * 1e4 / G_REF
        M["eff_load_calc"] = e
        assert abs(e - EFF_LOAD) < 0.002

    def t_table():
        iv = RAD_I * RAD_V
        M["IxV_vs_P_datasheet"] = [[float(a), float(b)] for a, b in zip(iv, RAD_P)]
        M["dif_max_IxV_P"] = float(np.max(np.abs(iv - RAD_P)))
        assert np.max(np.abs(iv - RAD_P)) < 0.0125

    def t_mono():
        p = [cell_mpp(t, G_DESIGN, FLUENCE_EOL)["pmp_W_m2"] for t in np.arange(15.0, 90.0, 5.0)]
        assert all(p[i + 1] < p[i] for i in range(len(p) - 1))
        q = [cell_mpp(65.0, G_DESIGN, f)["pmp_W_m2"] for f in (1e14, 2e14, 5e14, 1e15)]
        assert all(q[i + 1] < q[i] for i in range(len(q) - 1))

    def t_edge():
        c = cell_mpp(28.0, G_REF, 1e14)
        assert abs(c["pmp_W_m2"] / cell_mpp(28.0, G_REF, 0.0)["pmp_W_m2"] - 0.95) < 1e-9

    def t_lin():
        a, b = size_array(), size_array(p_req_w=2 * P_REQ_W)
        assert abs(b["A_cells_m2"] / a["A_cells_m2"] - 2.0) < 1e-12

    def t_balance():
        t, r = t_eq_c(G_DESIGN, FLUENCE_EOL, 0.5)
        M["t_eq_eps0.5"] = t
        M["residuo_balanco_W_m2"] = r
        assert abs(r) < 1e-6

    def t_domain():
        try:
            cell_mpp(65.0, G_DESIGN, 5e13)
        except ValueError:
            return
        raise AssertionError("deveria recusar fluencia fora do datasheet")

    def t_repro():
        assert size_array() == size_array()

    for n, f in (("sigma", t_sigma), ("eficiencia_ref", t_eff), ("fill_factor", t_ff), ("eficiencia_carga", t_load), ("tabela_radiacao_IxV_vs_P", t_table),
                 ("monotonia_T_fluencia", t_mono), ("limite_1e14", t_edge), ("linearidade_area", t_lin), ("balanco_termico_residuo", t_balance),
                 ("dominio_fluencia", t_domain), ("reprodutibilidade", t_repro)):
        chk(n, f)
    return res, M


def run_report():
    base = size_array()
    out = {"G_design_Wm2_afelio": G_DESIGN, "R_au_afelio_periélio": [R_AU_APH, R_AU_PER], "base": base}
    out["razoes_degradacao_em_FLUENCE_EOL"] = {"Pmp": _log_interp(FLUENCE_EOL, RAD_P), "Vmp": _log_interp(FLUENCE_EOL, RAD_V), "Imp": _log_interp(FLUENCE_EOL, RAD_I),
                                              "dJmp_uA": _lin_coef(FLUENCE_EOL, DJMP) * 1e6, "dVmp_mV": _lin_coef(FLUENCE_EOL, DVMP) * 1e3}
    # modelo da celula nas condicoes de projeto, por etapa [CALCULADO]
    c28 = cell_mpp(28.0, G_REF)
    out["etapas_W_m2_celula"] = {"datasheet_28C_G1353": c28["pmp_W_m2"], "aphelio_28C_BOL": cell_mpp(28.0, G_DESIGN)["pmp_W_m2"],
                                 "aphelio_65C_BOL": cell_mpp(65.0, G_DESIGN)["pmp_W_m2"], "aphelio_65C_EOL_rad": cell_mpp(65.0, G_DESIGN, FLUENCE_EOL)["pmp_W_m2"],
                                 "aphelio_65C_EOL_rad_UV": base["p_cell_util_W_m2"]}
    sens = []
    def row(nome, **kw):
        d = size_array(**kw)
        sens.append([nome, d["A_cells_m2"], 100.0 * (d["A_cells_m2"] / base["A_cells_m2"] - 1.0), d["A_panel_m2"], d["n_cells"], d["extrapolado_T"]])
    row("base")
    for T in (45.0, 85.0):
        row("T_array=%g C" % T, T_c=T)
    for f in (1e14, 5e14):
        row("fluencia=%.0e" % f, fluence=f)
    row("UV=0", uv=0.0)
    row("UV=10%", uv=0.10)
    row("f_other=0.97", f_other=0.97)
    row("pf=0.80", pf=0.80)
    row("pf=0.95", pf=0.95)
    row("perielio G", G=S_NOM / R_AU_PER ** 2)
    row("G=1361 (sem 1/R2)", G=S_NOM)
    row("apontamento 2 deg", pointing_deg=2.0)
    row("celula 26.62 cm2", cell_cm2=26.62)
    out["sensibilidade_OAT_[caso | A_celulas_m2 | dA_% | A_painel_m2 | n_celulas | T_extrapolada]"] = sens
    # equilibrio termico (diagnostico)
    teq = {}
    for eb in (0.30, 0.50, 0.85):
        t, r = t_eq_c(G_DESIGN, FLUENCE_EOL, eb)
        teq["eps_back=%.2f" % eb] = {"T_eq_C": t, "residuo": r}
    out["T_equilibrio_celula_C_(sem substrato, sem IR/albedo)"] = teq
    # legado
    out["legado_3_2"] = {"A_celulas_m2": 52416 * 30.18e-4, "n_celulas": 52416, "P_pre_buck_65C_EOL_kW_@1361": 49.733 / 0.91,
                         "Pmp_catalogo_legado_mW_cm2": 17.5 * 2.406, "Pmp_catalogo_datasheet_mW_cm2": 17.02 * 2.348}
    # tensoes para 2b (informativo)
    v_cold = {}
    for T in (28.0, 65.0, 75.0):
        c = cell_mpp(T, G_DESIGN, FLUENCE_EOL)
        v_cold["T=%g" % T] = {"Vmp_EOL_V": c["vmp_V"], "Voc_BOL_V": VOC + _lin_coef(0.0, DVOC) * (T - T_REF)}
    out["tensoes_celula_para_2b"] = v_cold
    return out


if __name__ == "__main__":
    tests, meas = run_selftests()
    rep = run_report()
    payload = {"testes": tests, "medidas": meas, "relatorio": rep}
    print(json.dumps(payload, indent=1, default=float, ensure_ascii=False))
    result = {"json": json.dumps(payload, default=float, ensure_ascii=False)}
