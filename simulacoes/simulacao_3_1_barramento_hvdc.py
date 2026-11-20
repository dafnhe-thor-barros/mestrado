#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ==============================================================================
# SIMULAÇÃO 3.1: DINÂMICA DO BARRAMENTO HVDC (400 V), BALANÇO DE ENERGIA
# MULTI-ÓRBITA DO BESS, CURTO-CIRCUITO E PROTEÇÃO ATIVA POR CHAVE SSPC
# ==============================================================================
# Dissertação de Mestrado Profissional - ProPECAut / IFES
# Autor: Dafnhê Thor Castro de Martins Barros
# Revisão Técnica: Cargas com OISL Contínuo (440 W), ndab = 0.95 e nbuck = 0.91
# ==============================================================================

import os
import csv
import json
import math
import time
import numpy as np
import matplotlib.pyplot as plt

N_ORB = 20                       # Número de órbitas completas simuladas
OUT_DIR = "saida_3_1"

PAR = {
    # Potências da carga e plataforma
    "p_node": 50000.0,           # Potência nominal de pico do nó [W] (50 kW)
    "p_it_nom": 40000.0,         # Potência máxima de TI computacional [W] (40 kW)
    "p_it_sun_nom": 32000.0,     # Carga de cruzeiro de TI em insolação [W]
    "p_it_ecl_falha": 20000.0,   # Carga em contingência BAT-F01 (sem modulação) [W]
    
    # AJUSTE OISL: 560 W (aviônica/bombas) + 440 W (4x terminais laser OISL) = 1000 W
    "p_aux": 1000.0,             # Cargas contínuas essenciais de plataforma e OISL [W]
    "p_oisl_eletrico": 440.0,    # Potência elétrica dedicada aos 4 terminais OISL [W]
    "q_oisl": 343.2,             # Calor residual interno gerado pelos lasers no chassi [W]

    # Barramento primário estabilizado HVDC
    "v_nom": 400.0,              # Tensão nominal de barramento [V]
    "v_uv": 324.0,               # Limiar crítico de subtensão UVLO (3,0 V/célula) [V]
    "v_ov": 442.8,               # Limiar de sobretensão CC-CV (4,10 V/célula) [V]

    # Banco de Baterias BESS (108S1P Li-ion NMC)
    "n_cells": 108,              # Número de células em série
    "ocv0": 3.0, "ocv1": 1.1,    # Relação linear Voc = n*(ocv0 + ocv1*SoC)
    "e_bat_j": 9.0 * 3.6e6,      # Energia nominal armazenada (9,0 kWh) [J]
    "cap_ah": 22.52,             # Capacidade nominal [Ah]
    "c_rate_chg": 1.0,           # Limite de taxa de recarga CC (1,0C = 22,52 A)
    "eta_chg": 0.995,            # Eficiência coulômbica de recarga
    "eta_dis": 0.970,            # Eficiência eletroquímica de descarga
    "r0": 0.085,                 # Resistência interna equivalente EOL a frio [Ohm]
    "soc_min_req": 0.20,         # Piso mínimo regulamentar de Estado de Carga (20%)

    # Conversores estáticos de potência
    "eta_conv": 0.950,           # Eficiência do conversor de distribuição DAB (ndab = 0.95)
    "eta_conv_wca": 0.900,       # Eficiência degradada no pior caso (WCA = 90,0%)
    "p_pv_eol": 49733.0,         # Potência PV EOL em 65°C após o Buck (nbuck = 0.91) [W]

    # Parâmetros de curto-circuito e proteção SSPC (SiC 1200 V)
    "r_line": 14.0e-3,           # Resistência dos cabos de distribuição [Ohm]
    "r_sc": 1.0e-3,              # Resistência de contato do curto-circuito franco [Ohm]
    "r_ds_on": 0.010,            # Resistência de condução do interruptor SSPC [Ohm]
    "v_clamp": 550.0,            # Tensão de corte do diodo de grampeamento TVS [V]
    "i2t_lim": 1500.0,           # Limite adiabático de fusão do SiC MOSFET [A^2 s]
    "v_ds_max": 1200.0,          # Tensão de ruptura máxima da chave de potência [V]
    "v_bat_sc": 399.6, "r0_bat_sc": 0.045, "l_bat": 2.5e-6, "tau_ext": 0.5e-6,

    # Modelo térmico concentrado do chassi (1 nó)
    "c_th_chassi": 18000.0,      # Capacitância térmica da estrutura metálica [J/K]
    "t_chassi0": 25.0,           # Temperatura inicial do chassi [°C]
    "t_lim_chassi": 65.0,        # Limite operacional máximo admissível [°C]
    "t_ref_chassi": 5926.0,      # Período de referência orbital [s]
}

ORBITS = [
    {"name": "LEO_500km",   "T": 5677.0,  "t_ecl": 2148.0, "p_it_ecl": 6500.0},
    {"name": "SSO_700km",   "T": 5926.0,  "t_ecl": 1128.0, "p_it_ecl": 18000.0},
    {"name": "GEO_35786km", "T": 86164.0, "t_ecl": 4320.0, "p_it_ecl": 4500.0},
]

PV_SCEN = [
    ("EOL_85C_nbuck91", 46137.0),
    ("EOL_65C_nbuck91", 49733.0),
    ("EOL_28C_nbuck91", 54183.0),
]

SC_CASES = [
    ("Nominal 8 us",            8e-6,    1.8e-6, 4.7e-3, True,  False),
    ("Nominal 15 us",           15e-6,   1.8e-6, 4.7e-3, True,  False),
    ("WCA 74.54 us",             74.54e-6, 4.0e-6, 1.5e-3, True,  False),
    ("WCA 100 us",              100e-6,  4.0e-6, 1.5e-3, True,  False),
    ("WCA 100 us + ramo bateria", 100e-6,  4.0e-6, 1.5e-3, True,  True),
    ("Falha sem TVS 15 us",     15e-6,   4.0e-6, 1.5e-3, False, False),
]

SINKS = [45.0, 51.94, 55.0]      # Temperatura do radiador sumidouro [°C]


def conv_orbit(series, tol):
    ref = series[-1]
    for j in range(len(series)):
        if all(abs(x - ref) <= tol for x in series[j:]):
            return j + 1
    return len(series)


def bess_run(par, orb, mode="nominal", n_orb=N_ORB, dt=2.0, p_it_sun=None, p_pv=None):
    T = orb["T"]
    te = orb["t_ecl"]

    if mode == "falha":
        p_it_ecl = par["p_it_ecl_falha"]
        p_it_sol = par["p_it_nom"]
    else:
        p_it_ecl = orb["p_it_ecl"]
        p_it_sol = par["p_it_sun_nom"]

    if p_it_sun is not None:
        p_it_sol = p_it_sun
    if p_pv is None:
        p_pv = par["p_pv_eol"]

    eta_c = par["eta_conv"]  # ndab = 0.95
    p_ecl_bus = (p_it_ecl + par["p_aux"]) / eta_c
    surplus = max(0.0, p_pv - (p_it_sol + par["p_aux"]) / eta_c)

    n = int(round(n_orb * T / dt))
    soc = np.empty(n + 1)
    v = np.empty(n + 1)
    ib = np.empty(n + 1)
    soc[0] = 1.0

    E = par["e_bat_j"]
    r0 = par["r0"]
    nc = par["n_cells"]
    o0 = par["ocv0"]
    o1 = par["ocv1"]
    v_ov = par["v_ov"]
    v_uv = par["v_uv"]
    soc_req = par["soc_min_req"]
    i_cc_max = par["c_rate_chg"] * par["cap_ah"]
    eta_chg = par["eta_chg"]
    eta_dis = par["eta_dis"]

    trip_soc = None
    trip_uv = None
    collapse = None
    e_dis = 0.0
    e_chg = 0.0
    clip = 0.0

    for k in range(n):
        t = k * dt
        s = soc[k]
        ocv = nc * (o0 + o1 * s)

        # Descarga a potência constante no eclipse
        if (t % T) < te:
            disc = ocv * ocv - 4.0 * r0 * p_ecl_bus
            if disc <= 0.0:
                if collapse is None:
                    collapse = t / 60.0
                i_k = ocv / (2.0 * r0)
                v_k = 0.5 * ocv
            else:
                i_k = (ocv - math.sqrt(disc)) / (2.0 * r0)
                v_k = ocv - i_k * r0
            ds = -p_ecl_bus / (E * eta_dis) * dt
            e_dis += p_ecl_bus * dt
            ib_k = i_k
        # Recarga CC-CV no período de insolação
        else:
            if s >= 1.0:
                i_c = 0.0
                v_k = ocv
            else:
                i_pow = (-ocv + math.sqrt(ocv * ocv + 4.0 * r0 * surplus)) / (2.0 * r0)
                i_c = min(i_pow, i_cc_max)
                v_k = ocv + i_c * r0
                if v_k > v_ov:
                    i_c = max(0.0, (v_ov - ocv) / r0)
                    v_k = v_ov
            ds = i_c * v_k * eta_chg / E * dt
            e_chg += i_c * v_k * dt
            ib_k = -i_c

        s_new = s + ds
        if s_new > 1.0:
            clip += (s_new - 1.0) * E
            s_new = 1.0
        elif s_new < 0.0:
            clip += s_new * E
            s_new = 0.0

        soc[k + 1] = s_new
        v[k] = v_k
        ib[k] = ib_k

        if trip_soc is None and s_new < soc_req:
            trip_soc = (k + 1) * dt / 60.0
        if trip_uv is None and v_k < v_uv:
            trip_uv = t / 60.0

    v[n] = v[n - 1]
    ib[n] = ib[n - 1]

    closure = (soc[-1] - 1.0) * E - (e_chg * eta_chg - e_dis / eta_dis - clip)
    t_arr = np.arange(n + 1) * dt
    idx = np.minimum((t_arr / T).astype(int), n_orb - 1)
    bounds = np.searchsorted(idx, np.arange(n_orb + 1))
    pw = ib * v
    pw[-1] = 0.0

    orbits = []
    for j in range(n_orb):
        a = int(bounds[j])
        b = int(bounds[j + 1])
        if b - a < 2:
            continue
        sl = slice(a, b)
        p_j = pw[sl]
        orbits.append({
            "j": j + 1,
            "soc_min": float(soc[sl].min() * 100.0),
            "soc_end": float(soc[b - 1] * 100.0),
            "v_min": float(v[sl].min()),
            "v_max": float(v[sl].max()),
            "i_dis_max": float(ib[sl].max()),
            "i_chg_max": float(-ib[sl].min()),
            "e_dis_kwh": float(p_j[p_j > 0].sum() * dt / 3.6e6),
            "e_chg_kwh": float(-p_j[p_j < 0].sum() * dt / 3.6e6)
        })

    st = max(1, int(round(10.0 / dt)))
    return {
        "orbit": orb["name"], "T": T, "mode": mode,
        "p_it_sol": p_it_sol, "p_pv": p_pv, "orbits": orbits,
        "trip_soc_min": trip_soc, "trip_uv_min": trip_uv,
        "collapse_min": collapse, "closure_rel": float(closure / E),
        "t_h": t_arr[::st] / 3600.0, "soc": soc[::st] * 100.0,
        "v": v[::st], "i": ib[::st]
    }


def bess_recovers(res, tol=99.9):
    return (
        res["trip_soc_min"] is None and
        res["trip_uv_min"] is None and
        res["collapse_min"] is None and
        all(o["soc_end"] >= tol for o in res["orbits"])
    )


def p_it_sun_max(par, orb, p_pv, n_orb=N_ORB, dt=2.0):
    hi = min(par["p_it_nom"], p_pv * par["eta_conv"] - par["p_aux"])
    r_hi = bess_run(par, orb, "nominal", 8, 20.0, p_it_sun=hi, p_pv=p_pv)
    if bess_recovers(r_hi):
        lo = hi
    else:
        lo = 0.0
        for _ in range(16):
            mid = 0.5 * (lo + hi)
            r = bess_run(par, orb, "nominal", 8, 20.0, p_it_sun=mid, p_pv=p_pv)
            if bess_recovers(r):
                lo = mid
            else:
                hi = mid
    ps = lo
    r = None
    for _ in range(60):
        r = bess_run(par, orb, "nominal", n_orb, dt, p_it_sun=ps, p_pv=p_pv)
        if bess_recovers(r):
            return ps, r
        ps -= 100.0
        if ps < 0.0:
            break
    return 0.0, r


def short_circuit(par, t_trip, L, C, clamp=True, battery=False, dt=10e-9):
    R = par["r_line"] + par["r_sc"]
    vcl = par["v_clamp"]
    vb = par["v_bat_sc"]
    r0b = par["r0_bat_sc"]
    Lb = par["l_bat"]
    rds = par["r_ds_on"]

    n_max = int((t_trip + 600e-6) / dt)
    iL = 0.0
    ib = 0.0
    vc = par["v_nom"]
    i2t = 0.0
    i2t_pre = None
    i_trip = 0.0
    vc_trip = vc
    e_tvs = 0.0
    ipk = 0.0
    vsw_max = 0.0
    t_end_cond = None
    extinguished = False

    for k in range(n_max):
        t = k * dt
        if t < t_trip:
            vs = iL * rds
            iL = iL + (vc - R * iL) / L * dt
        else:
            if i2t_pre is None:
                i2t_pre = i2t
                i_trip = iL
                vc_trip = vc
            if iL <= 1e-3:
                t_end_cond = t
                extinguished = True
                break
            if clamp:
                vs = vcl
                iL = max(0.0, iL + (vc - vcl - R * iL) / L * dt)
                e_tvs += vcl * iL * dt
            else:
                vs = vc + L * iL / par["tau_ext"] - R * iL
                iL = max(0.0, iL - iL / par["tau_ext"] * dt)

        if battery:
            ib = ib + (vb - r0b * ib - vc) / Lb * dt
            vc = vc + (ib - iL) / C * dt
        else:
            vc = vc - iL / C * dt

        i2t += iL * iL * dt
        if iL > ipk:
            ipk = iL
        if vs > vsw_max:
            vsw_max = vs

    if i2t_pre is None:
        i2t_pre = i2t
    if t_end_cond is None:
        t_end_cond = n_max * dt

    return {
        "t_trip_us": t_trip * 1e6, "ipk_kA": ipk / 1e3, "i_trip_kA": i_trip / 1e3,
        "vc_trip_V": vc_trip, "i2t_total": i2t, "i2t_pre_trip": i2t_pre,
        "e_tvs_J": e_tvs, "v_sw_max_V": vsw_max,
        "t_ext_us": (t_end_cond - t_trip) * 1e6, "extinguished": extinguished
    }


def chassis_run(par, eta, g_th, t_sink, n_orb=N_ORB, dt=2.0):
    t_ref = par["t_ref_chassi"]
    p_conv = par["p_node"] * (1.0 - eta)
    q_o = par["q_oisl"]
    c = par["c_th_chassi"]

    n = int(round(n_orb * t_ref / dt))
    T = np.empty(n + 1)
    T[0] = par["t_chassi0"]

    for k in range(n):
        tc = (k * dt) % 1800.0
        duty = 0.4 + (0.6 if tc < 1200.0 else 0.0)
        q = p_conv * duty + q_o
        T[k + 1] = T[k] + (q - g_th * (T[k] - t_sink)) / c * dt

    idx = np.minimum((np.arange(n + 1) * dt / t_ref).astype(int), n_orb - 1)
    bounds = np.searchsorted(idx, np.arange(n_orb + 1))
    tmax_orb = [float(T[bounds[j]:bounds[j + 1]].max()) for j in range(n_orb)]
    q_avg = p_conv * 0.8 + q_o

    return {
        "p_loss_pk": p_conv + q_o,
        "t_max_last": tmax_orb[-1],
        "t_eq_medio": t_sink + q_avg / g_th,
        "margin": par["t_lim_chassi"] - tmax_orb[-1],
        "status": "OK" if (par["t_lim_chassi"] - tmax_orb[-1]) >= 0 else "VIOLA"
    }


def gerar_graficos(bess_resultados):
    """Gera e salva os gráficos das simulações na pasta saida_3_1 em formato PNG."""
    print(f"\n[INFO] Gerando gráficos e salvando em ./{OUT_DIR}/...")
    
    # 1. Gráfico de Estado de Carga (SoC) por órbita para os casos nominais
    plt.figure(figsize=(10, 5))
    for res in bess_resultados:
        if "nominal" in res["caso"]:
            plt.plot(res["t_h"], res["soc"], label=f"{res['orbita']} - {res['caso']}")
    plt.axhline(20.0, color='r', linestyle='--', label='Piso Mínimo (20%)')
    plt.title("Simulação 3.1: Dinâmica do SoC do BESS (Multi-Órbita)")
    plt.xlabel("Tempo de Simulação [h]")
    plt.ylabel("Estado de Carga - SoC [%]")
    plt.grid(True)
    plt.legend(loc='best', fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "bess_soc_multiorbita.png"), dpi=300)
    plt.close()

    # 2. Gráfico de Tensão nas Células do BESS
    plt.figure(figsize=(10, 5))
    for res in bess_resultados:
        if "nominal" in res["caso"]:
            plt.plot(res["t_h"], res["v"], label=f"{res['orbita']} - {res['caso']}")
    plt.axhline(PAR["v_uv"], color='orange', linestyle='--', label='Limiar Subtensão (UVLO)')
    plt.axhline(PAR["v_ov"], color='purple', linestyle='--', label='Limiar Sobretensão (CC-CV)')
    plt.title("Simulação 3.1: Tensão do Banco de Baterias (BESS)")
    plt.xlabel("Tempo de Simulação [h]")
    plt.ylabel("Tensão do Banco [V]")
    plt.grid(True)
    plt.legend(loc='best', fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "bess_tensao_banco.png"), dpi=300)
    plt.close()

    # 3. Gráfico de Corrente do BESS
    plt.figure(figsize=(10, 5))
    for res in bess_resultados:
        if "nominal" in res["caso"]:
            plt.plot(res["t_h"], res["i"], label=f"{res['orbita']} - {res['caso']}")
    plt.title("Simulação 3.1: Corrente do Banco de Baterias (Descarga/Recarga)")
    plt.xlabel("Tempo de Simulação [h]")
    plt.ylabel("Corrente [A]")
    plt.grid(True)
    plt.legend(loc='best', fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "bess_corrente_banco.png"), dpi=300)
    plt.close()

    print("[INFO] Gráficos salvos com sucesso na pasta especificada.")


def main():
    t0 = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)
    print("=" * 80)
    print("SIMULAÇÃO 3.1: DINÂMICA DO BARRAMENTO E BESS (COM OISL CONTÍNUO = 440 W)")
    print("=" * 80)

    # 1. Simulação BESS
    bess_resumo = []
    bess_resultados_completos = []
    for orb in ORBITS:
        r_sus = p_it_sun_max(par=PAR, orb=orb, p_pv=PAR["p_pv_eol"])[1]
        cases = [
            ("nominal_32kW_script", bess_run(PAR, orb, "nominal", p_it_sun=PAR["p_it_sun_nom"])),
            ("nominal_40kW_cap4", bess_run(PAR, orb, "nominal", p_it_sun=PAR["p_it_nom"])),
            ("nominal_Pmax_sustentavel", r_sus),
            ("falha_BAT_F01", bess_run(PAR, orb, "falha"))
        ]
        for cname, r in cases:
            ob = r["orbits"]
            ok = (
                r["trip_soc_min"] is None and
                r["trip_uv_min"] is None and
                min(o["soc_min"] for o in ob) >= 20.0
            )
            bess_resumo.append({
                "caso": cname, "orbita": orb["name"],
                "p_it_sol_kW": r["p_it_sol"] / 1e3,
                "soc_min_pct": ob[-1]["soc_min"],
                "v_min_V": ob[-1]["v_min"],
                "status": "OK" if ok else "VIOLA"
            })
            bess_resultados_completos.append(r)

    # 2. Resumo da Carga Máxima de TI Sustentável
    pmax_resumo = []
    for orb in ORBITS:
        for pvn, pvv in PV_SCEN:
            ps, rr = p_it_sun_max(PAR, orb, pvv)
            pmax_resumo.append({
                "orbita": orb["name"], "cenario_pv": pvn,
                "p_pv_kW": pvv / 1e3, "p_it_sun_max_kW": ps / 1e3,
                "folga_kW": (ps - PAR["p_it_nom"]) / 1e3
            })

    print("\n--- RESUMO BESS MULTI-ÓRBITA (20 Órbitas | Cargas Contínuas = 1000 W) ---")
    for r in bess_resumo:
        print("%-26s %-12s TI_sol=%4.1f kW | SoCmin=%5.2f%% | Vmin=%5.1f V | %s" % (
            r["caso"], r["orbita"], r["p_it_sol_kW"], r["soc_min_pct"], r["v_min_V"], r["status"]))

    print("\n--- CARGA MÁXIMA DE TI EM SOL (Pmax Sustentável) ---")
    for r in pmax_resumo:
        print("%-12s %-18s PV=%5.2f kW -> TI_sol_max=%5.2f kW (Folga: %5.2f kW)" % (
            r["orbita"], r["cenario_pv"], r["p_pv_kW"], r["p_it_sun_max_kW"], r["folga_kW"]))

    # Chamada da função para salvar os gráficos em PNG na pasta saida_3_1
    gerar_graficos(bess_resumo_completo := bess_resultados_completos)

    print(f"\nSimulação 3.1 executada em {time.time() - t0:.2f} s. Saída em ./{OUT_DIR}/")


if __name__ == "__main__":
    main()