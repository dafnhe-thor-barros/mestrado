#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ==============================================================================
# SIMULAÇÃO 3.2: GERADOR FOTOVOLTAICO MULTIJUNÇÃO AM0 (SPECTROLAB XTJ PRIME),
# MODELO DE 3 DIODOS EM SÉRIE (GaInP/GaAs/Ge), DEGRADAÇÃO DDD E BUCK MPPT
# ==============================================================================
# Dissertação de Mestrado Profissional - ProPECAut / IFES
# Autor: Dafnhê Thor Castro de Martins Barros
# Revisão Técnica: Modelo Físico 3 Diodos, Eficiência Buck (nbuck = 0.91) e Hot Spot
# ==============================================================================

import os
import csv
import json
import math
import time
import numpy as np

# Constantes Fundamentais da Física
K_B = 1.380649e-23        # Constante de Boltzmann [J/K]
Q_E = 1.602176634e-19     # Carga elementar do elétron [C]
SIG = 5.670374e-8         # Constante de Stefan-Boltzmann [W/(m^2 K^4)]

OUT_DIR = "saida_3_2"

# Dicionário de Parâmetros de Projeto e Físicos
PAR = {
    # Irradiância solar e geometria da célula (Spectrolab XTJ Prime)
    "g_am0": 1361.0,      # Irradiância solar AM0 na órbita da Terra [W/m^2]
    "g_ds": 1353.0,       # Irradiância padrão de ensaio em datasheet [W/m^2]
    "a_cm2": 30.18,       # Área ativa de cada célula [cm^2]
    "jsc": 18.0e-3,       # Densidade de corrente de curto-circuito nominal [A/cm^2]
    "jmp": 17.5e-3,       # Densidade de corrente no MPP nominal [A/cm^2]
    "vmp_cat": 2.406,     # Tensão de célula no MPP em catálogo a 28°C [V]

    # Parâmetros das 3 subcélulas (1: GaInP Top, 2: GaAs Middle, 3: Ge Bottom)
    "voc_sub": [1.42, 1.01, 0.29],                # Voc individual a 28°C (Soma = 2.720 V)
    "dvoc_sub": [-2.20e-3, -1.90e-3, -1.50e-3],   # Coeficientes térmicos dVoc/dT [V/K]
    "n_sub": [1.15, 1.15, 1.75],                  # Fatores de idealidade de cada junção
    "rsh": [2500.0, 3500.0, 200.0],               # Resistências shunt [Ohm]
    "djsc": 10.0e-6,                              # Coeficiente térmico dJsc/dT [A/(cm^2 K)]
    "w_dph": [0.85, 1.0, 1.8],                    # Ponderação da variação de fotocorrente
    "ratio": [1.0, 18.35 / 18.0, 26.5 / 18.0],   # Razão de fotocorrentes (Top limita corrente)

    # Modo de calibração da resistência série interna
    "rs_modo": "calibrado", # Calibra Rs para atingir exatamente Vmp_cat a 28°C
    "rs_texto": 0.20,       # Valor de referência do texto [Ohm]

    # Topologia da matriz solar do nó de 50 kW
    "ns": 208,              # Células em série por string (garante Vmp > 400V a 65°C)
    "np": 252,              # Strings em paralelo por satélite
    "vbus": 400.0,          # Tensão nominal do barramento primário HVDC [V]
    "v_ov": 442.8,          # Tensão máxima de recarga do pack (4,10 V/célula) [V]
    "d_max": 0.90,          # Razão cíclica máxima contínua do Buck (limite de projeto)

    # AJUSTE: Eficiência de conversão do estágio Buck solar
    "eta_buck": 0.910,      # Eficiência nominal de 91,0% (ajustado de 0,975)
    "p_alvo_bus": 46500.0,  # Potência líquida mínima requerida no barramento [W]

    # Parâmetros dinâmicos do conversor Buck e filtro de entrada
    "L": 45.0e-6,           # Indutância do conversor Buck [H] (45 uH)
    "cin": 250.0e-6,        # Capacitância do filtro de entrada [F] (250 uF)
    "rl": 0.025,            # Resistência série parasita do indutor [Ohm]
    "fsw": 150.0e3,         # Frequência de chaveamento PWM [Hz] (150 kHz)
    "fmppt": 10.0e3,        # Frequência de atualização do algoritmo MPPT [Hz] (10 kHz)
    "m_gain": 1.2e-4,       # Ganho adaptativo do MPPT P&O-AP [V/W]
    "dd_max": 0.002,        # Passo máximo do ciclo de trabalho no modo adaptativo
    "dd_fixo": 0.005,       # Passo fixo de comparação convencional
    "d_lo": 0.20, "d_hi": 0.98,

    # Degradação por Dose de Deslocamento (DDD) e óptica ambiental
    "cd": [0.085, 0.12851], # Coeficientes de degradação (nominal e conservador)
    "phi_x": 1.0e14,        # Fluência de escala de referência [e-/cm^2]
    "dp_env": 0.052,        # Perda acumulada por escurecimento UV em 10 anos (5,2%)
    "fluencia": [1.4e13, 2.5e13, 5.8e13], # Fluência anual (500 km, 700 km, 1000 km)
    "alt_km": [500, 700, 1000],
    "anos": 10.0,

    # Temperaturas de análise da matriz e rampas térmicas
    "temps": [-100.0, 15.0, 28.0, 45.0, 60.0, 65.0, 85.0],
    "rampas_Ks": [1100.0, 330.0, 0.25],

    # Parâmetros térmicos para análise de Hot Spot
    "hs_cth": 3.8,          # Capacitância térmica da célula unitária [J/K]
    "hs_eps": 0.85,         # Emissividade térmica da célula no infravermelho
    "hs_v_rev": 28.0,       # Tensão reversa de breakdown em sombra total [V]
    "hs_v_byp": 0.35,       # Queda direta no diodo de bypass Schottky [V]
    "hs_q_amb": 0.08,       # Carga ambiental residual absorvida [W]
    "hs_t0": 60.0,          # Temperatura inicial [°C]
    "hs_tlim": 180.0,       # Limite de delaminação do coverglass [°C]
}


def sub_params(par):
    """Calcula as correntes e tensões de saturação para cada uma das 3 subcélulas."""
    sc = par["g_am0"] / par["g_ds"]
    isc = par["jsc"] * par["a_cm2"] * sc
    return {
        "iph0": np.array(par["ratio"]) * isc,
        "diph": np.array(par["w_dph"]) * par["djsc"] * par["a_cm2"] * sc,
        "voc0": np.array(par["voc_sub"]),
        "dvoc": np.array(par["dvoc_sub"]),
        "n": np.array(par["n_sub"]),
        "rsh": np.array(par["rsh"]),
        "isc_cat": isc,
        "imp_cat": par["jmp"] * par["a_cm2"] * sc
    }


def v_sub(ic, iph, i0, vth, rsh, voc):
    """Resolve a tensão interna de uma subcélula pelo método de bissecção/Newton."""
    lo = np.zeros_like(ic)
    hi = np.full_like(ic, voc * 1.6)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        f = iph - i0 * (np.exp(np.minimum(mid / vth, 700.0)) - 1.0) - mid / rsh - ic
        pos = f > 0
        lo = np.where(pos, mid, lo)
        hi = np.where(pos, hi, mid)
    return 0.5 * (lo + hi)


def cell_iv(sp, rs, T, n=3000):
    """Calcula a curva I-V completa somando as 3 junções em série: Vcel = V1 + V2 + V3 - I*Rs."""
    dT = T - 28.0
    TK = T + 273.15
    iph = sp["iph0"] + sp["diph"] * dT
    voc = sp["voc0"] + sp["dvoc"] * dT
    vth = sp["n"] * K_B * TK / Q_E
    i0 = (iph - voc / sp["rsh"]) / (np.exp(voc / vth) - 1.0)

    ic = np.linspace(0.0, iph.min() * 0.9999, n)
    vs = np.zeros(n)
    for j in range(3):
        vs = vs + v_sub(ic, iph[j], i0[j], vth[j], sp["rsh"][j], voc[j])
    return ic, np.maximum(0.0, vs - ic * rs), float(voc.sum())


def mpp_cell(sp, rs, T):
    """Encontra o Ponto de Máxima Potência (MPP) unitário da célula a dada temperatura."""
    ic, v, voc = cell_iv(sp, rs, T)
    p = v * ic
    k = int(np.argmax(p))
    return {"vmp": float(v[k]), "imp": float(ic[k]), "pmp": float(p[k]), "voc": voc}


def calib_rs(sp, vmp_cat):
    """Calibra a resistência série interna para casar com a tensão de catálogo a 28°C."""
    lo, hi = 0.0, 0.3
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        if mpp_cell(sp, mid, 28.0)["vmp"] > vmp_cat:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def array_curve(sp, rs, T, ns, npar):
    """Escala a curva da célula unitária para toda a matriz (Ns células em série, Np strings)."""
    ic, v, voc = cell_iv(sp, rs, T)
    return ns * v[::-1], npar * ic[::-1]


def bus_power(par, sp, rs, T, vbus_bat, ret):
    """Calcula a potência EOL entregue ao barramento HVDC pós-Buck com eta_buck = 0.91."""
    ic, v, voc = cell_iv(sp, rs, T)
    p = v * ic
    k = int(np.argmax(p))
    ns = par["ns"]

    # Restrição de saturação de razão cíclica D <= 0.90
    vmin_cell = (vbus_bat / par["d_max"]) / ns
    vref = max(v[k], vmin_cell)
    order = np.argsort(v)
    pat = float(np.interp(vref, v[order], p[order]))

    # Aplicação da eficiência do Buck ajustada para 91,0%
    p_str_bus = ns * pat * ret * par["eta_buck"]
    return {
        "vref_str": ns * vref,
        "p_mpp_str_W": ns * p[k],
        "p_bus_str_W": p_str_bus,
        "fracao_mpp": pat / float(p[k]),
        "np_req": math.ceil(par["p_alvo_bus"] / p_str_bus)
    }


def hotspot(par, pdiss, faces, q_amb=0.0, dt=0.01, tend=240.0):
    """Calcula a evolução térmica no vácuo de uma célula sombreada via Stefan-Boltzmann."""
    a = par["a_cm2"] * 1e-4
    T = par["hs_t0"] + 273.15
    t_lim = None
    tl = par["hs_tlim"] + 273.15

    for k in range(int(tend / dt)):
        q = faces * par["hs_eps"] * SIG * a * (T ** 4 - 3.0 ** 4)
        T = T + (pdiss + q_amb - q) * dt / par["hs_cth"]
        if t_lim is None and T >= tl:
            t_lim = (k + 1) * dt
    teq = ((pdiss + q_amb) / (faces * par["hs_eps"] * SIG * a) + 3.0 ** 4) ** 0.25 - 273.15
    return t_lim, teq, T - 273.15


def ddd_ret(par, cd, rate, yrs=None):
    """Calcula a retenção de potência pelo modelo de Dose de Deslocamento (DDD)."""
    if yrs is None:
        yrs = par["anos"]
    return (1.0 - cd * math.log10(1.0 + yrs * rate / par["phi_x"])) * (1.0 - par["dp_env"] * yrs / 10.0)


def run_32(par):
    out = {}
    sp = sub_params(par)
    rs_cal = calib_rs(sp, par["vmp_cat"])
    rs = rs_cal if par["rs_modo"] == "calibrado" else par["rs_texto"]
    ns, npar = par["ns"], par["np"]
    out["rs"] = rs
    out["rs_cal"] = rs_cal

    # 1. Curvas I-V da matriz
    out["iv"] = []
    for T in par["temps"]:
        m = mpp_cell(sp, rs, T)
        out["iv"].append({
            "T": T, "vmp_cel": m["vmp"], "vmp_str": ns * m["vmp"],
            "voc_str": ns * m["voc"], "imp_arr": npar * m["imp"],
            "pmp_kW": ns * npar * m["pmp"] / 1e3,
            "d_req_400": par["vbus"] / (ns * m["vmp"])
        })

    # 2. Degradação DDD
    out["ddd"] = []
    for cd in par["cd"]:
        for rate, alt in zip(par["fluencia"], par["alt_km"]):
            out["ddd"].append({"cd": cd, "alt_km": alt, "ret_pct": 100.0 * ddd_ret(par, cd, rate)})

    # 3. Potência entregue ao barramento HVDC EOL em SSO 700 km
    ret700 = ddd_ret(par, par["cd"][-1], par["fluencia"][1])
    out["bus"] = []
    for T in (28.0, 65.0, 85.0):
        for vb in (par["vbus"], par["v_ov"]):
            b = bus_power(par, sp, rs, T, vb, ret700)
            out["bus"].append({
                "T": T, "v_bat": vb, "v_ref_str": b["vref_str"],
                "p_mpp_BOL_kW": b["p_mpp_str_W"] * npar / 1e3,
                "p_bus_EOL_kW": b["p_bus_str_W"] * npar / 1e3,
                "fracao_mpp": b["fracao_mpp"], "np_req_alvo": b["np_req"]
            })

    # 4. Hot spot com e sem diodo de bypass
    out["hot"] = []
    imp = sp["imp_cat"]
    for nm, pd, f, qa in (
        ("sem_bypass_2faces", par["hs_v_rev"] * imp, 2, 0.0),
        ("sem_bypass_1face", par["hs_v_rev"] * imp, 1, 0.0),
        ("com_bypass_2faces", par["hs_v_byp"] * imp, 2, 0.0),
        ("com_bypass_2faces_q_amb", par["hs_v_byp"] * imp, 2, par["hs_q_amb"])
    ):
        t_l, teq, t_end = hotspot(par, pd, f, qa)
        out["hot"].append({"caso": nm, "p_diss_W": pd, "t_180_s": t_l, "t_eq_C": teq, "t_240s_C": t_end})

    return out


def write_csv(path, rows):
    if not rows:
        return
    keys = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter=";")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("%.6g" % v if isinstance(v, float) else v) for k, v in r.items()})


def main():
    t0 = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)
    R = run_32(PAR)

    for nome in ("iv", "ddd", "bus", "hot"):
        write_csv(os.path.join(OUT_DIR, nome + ".csv"), R[nome])

    print("=== POTENCIA NO BARRAMENTO (EOL SSO 700 km | nbuck = 0.91) ===")
    for r in R["bus"]:
        print("T=%4.0f C | Vbat=%5.1f V | Vref=%6.1f V | P_BOL=%5.2f kW | P_EOL_bus=%6.2f kW | Np_req=%d" % (
            r["T"], r["v_bat"], r["v_ref_str"], r["p_mpp_BOL_kW"], r["p_bus_EOL_kW"], r["np_req_alvo"]))
    print(f"\nSimulação concluída em {time.time() - t0:.2f} s. Saída em ./{OUT_DIR}/")


if __name__ == "__main__":
    main()