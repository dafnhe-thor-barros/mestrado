#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ==============================================================================
# SIMULAÇÃO 3.3 (EXPANDIDA): DINÂMICA ELETROTÉRMICA E BALANÇO DE ENERGIA
# DO BANCO DE BATERIAS (BESS) EM 100 ÓRBITAS LEO
# ==============================================================================
# Dissertação de Mestrado Profissional - ProPECAut / IFES
# Autor: Dafnhê Thor Castro de Martins Barros
# Revisão Técnica: Modelagem Explícita de Condições Nominais e Modos de Falha
# ==============================================================================

import numpy as np


def executar_bateria_100_orbitas():
    print("=" * 111)
    print("INICIANDO SIMULAÇÃO 3.3 EXPANDIDA: BALANÇO ELETROTÉRMICO DO BESS EM 100 ÓRBITAS")
    print("=" * 111 + "\n")

    Cn = 22.52
    Vuv = 324.0
    Vov = 442.8
    Rp1 = 18e-3
    Cp1 = 3200.0
    Rp2 = 24e-3
    Cp2 = 22000.0
    R0_ref = 85e-3
    Ea = 0.32 * 1.602176634e-19
    kB = 1.380649e-23
    Tref_K = 293.15

    Cth = 54360.0
    Rth_cond = 0.14
    Tch = 15.0
    Pheater = 300.0
    Ton = 10.0
    Toff = 15.0

    dt = 2.0
    Torb = 5926.0
    tecl = 1128.0
    nOrb = 100
    steps_per_orb = int(np.floor(Torb / dt))
    N = nOrb * steps_per_orb + 1

    soc = np.zeros(N)
    vt = np.zeros(N)
    Tb = np.zeros(N)
    Pj = np.zeros(N)

    soc[0] = 0.90
    Tb[0] = 20.0
    vp1 = 0.0
    vp2 = 0.0
    a1 = np.exp(-dt / (Rp1 * Cp1))
    a2 = np.exp(-dt / (Rp2 * Cp2))

    alvos = [1, 5, 10, 20, 50, 100]

    print(f"{'Orbita':<8} | {'Tempo':<8} | {'SoC Min/Max (%)':<16} | {'Tensao Min/Max(V)':<16} | {'T_bat (Min/Med/Max) C':<24} | {'Pj_med (W)':<10} | {'Teq_est C':<10}")
    print("-" * 111)

    for k in range(N - 1):
        t_now = k * dt
        in_eclipse = (t_now % Torb) <= tecl

        if in_eclipse:
            I = (18000.0 + 560.0) / 380.0
            etaC = 1.0
        else:
            I = -25.0
            etaC = 0.995

        R0 = R0_ref * np.exp((Ea / kB) * (1.0 / (Tb[k] + 273.15) - 1.0 / Tref_K))
        vp1 = a1 * vp1 + Rp1 * (1.0 - a1) * I
        vp2 = a2 * vp2 + Rp2 * (1.0 - a2) * I

        soc[k + 1] = min(1.0, max(0.0, soc[k] - I * etaC * dt / (Cn * 3600.0)))
        Voc = 108.0 * (3.0 + 1.1 * soc[k + 1])
        Vcand = Voc - I * R0 - vp1 - vp2

        if I < 0.0 and Vcand > Vov:
            I = max(-25.0, min(0.0, -(Vov - Voc + vp1 + vp2) / R0))
            Vcand = Vov
        vt[k + 1] = Vcand

        Pj[k] = (I**2) * R0 + (vp1**2) / Rp1 + (vp2**2) / Rp2
        if Tb[k] < Ton:
            Ph = Pheater
        elif Tb[k] < Toff:
            Ph = Pheater * (Toff - Tb[k]) / (Toff - Ton)
        else:
            Ph = 0.0

        q_out = (Tb[k] - Tch) / Rth_cond
        Tb[k + 1] = Tb[k] + (Pj[k] + Ph - q_out) * dt / Cth

        if (k + 1) % steps_per_orb == 0:
            orb_atual = (k + 1) // steps_per_orb
            if orb_atual in alvos:
                idx_start = (orb_atual - 1) * steps_per_orb
                idx_end = k + 2

                soc_min = np.min(soc[idx_start:idx_end]) * 100.0
                soc_max = np.max(soc[idx_start:idx_end]) * 100.0
                vt_min = np.min(vt[idx_start:idx_end])
                vt_max = np.max(vt[idx_start:idx_end])
                tb_min = np.min(Tb[idx_start:idx_end])
                tb_med = np.mean(Tb[idx_start:idx_end])
                tb_max = np.max(Tb[idx_start:idx_end])
                pj_med = np.mean(Pj[idx_start:k+1])
                teq_est = Tch + pj_med * Rth_cond
                t_horas = ((k + 1) * dt) / 3600.0

                print(f"{orb_atual:<8} | {t_horas:6.1f} h | [{soc_min:5.2f}, {soc_max:6.2f}]% | [{vt_min:5.1f}, {vt_max:5.1f}] V | [{tb_min:5.2f}, {tb_med:5.2f}, {tb_max:5.2f}] C | {pj_med:8.1f} W | {teq_est:7.2f} C")

    print("=" * 111)
    print("SIMULAÇÃO DE 100 ÓRBITAS CONCLUÍDA COM SUCESSO.")
    print("=" * 111)


if __name__ == '__main__':
    executar_bateria_100_orbitas()