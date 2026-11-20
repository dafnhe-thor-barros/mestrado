#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ==============================================================================
# Simulacao Numerica de Transiente de Curto-Circuito em Barramento HVDC (400 V)
# ==============================================================================
# Dissertação de Mestrado Profissional - ProPECAut / IFES
# Autor: Dafnhê Thor Castro de Martins Barros
# Revisão Técnica: Cargas com OISL Contínuo (440 W), Rth = 0.14 K/W e AUKF
# ==============================================================================

"""
Simulacao Numerica de Transiente de Curto-Circuito em Barramento HVDC (400 V)

Determinacao dos limites criticos de abertura da chave SSPC (SiC MOSFET):
1. Tempo critico de desarme (t_trip) para atingir I^2t = 1500 A^2s.
2. Tempo de desarme que atinge V_spike = 1200 V.
3. Tempo minimo de extincao (soft turn-off) para evitar V_spike > 1200 V.

ATENCAO:
- O tempo de extincao calculado na secao 3 usa uma estimativa analitica:
  V_spike ~= V_bus + L * I / tau_ext.
- Os resultados dependem das hipoteses e limitacoes do modelo numerico.
"""

import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# MODELO TRANSITORIO
# ============================================================

def simular_curto(
    t_trip,
    L_line,
    C_bus,
    R_line=14.0e-3,
    R_sc=1.0e-3,
    V_bus=400.0,
    V_clamp=550.0,
    clamp_ativo=True,
    tau_ext=0.5e-6
):
    dt = 0.01e-6  # 10 ns
    t_end = max(150e-6, t_trip + 40e-6)
    t = np.arange(0.0, t_end, dt)

    i_curr = 0.0
    v_curr = V_bus
    i2t = 0.0
    v_sw_max = 0.0
    i_peak = 0.0

    for tk in t:

        if tk < t_trip:
            di_dt = (
                v_curr - i_curr * (R_line + R_sc)
            ) / L_line

            dv_dt = -i_curr / C_bus

            i_curr = i_curr + di_dt * dt
            v_curr = v_curr + dv_dt * dt

            v_sw = i_curr * 0.010  # Modelo original de Rds_on

        else:
            if clamp_ativo:
                di_dt = -(
                    V_clamp - v_curr + i_curr * R_line
                ) / L_line

                dv_dt = -i_curr / C_bus

                i_curr = max(0.0, i_curr + di_dt * dt)
                v_curr = v_curr + dv_dt * dt

                v_sw = V_clamp if i_curr > 0.0 else v_curr

            else:
                di_dt = -i_curr / tau_ext
                v_indutivo = abs(L_line * di_dt)

                i_curr = max(0.0, i_curr + di_dt * dt)

                dv_dt = -i_curr / C_bus
                v_curr = v_curr + dv_dt * dt

                v_sw = v_curr + v_indutivo

        i2t += (i_curr ** 2) * dt

        if i_curr > i_peak:
            i_peak = i_curr

        if v_sw > v_sw_max:
            v_sw_max = v_sw

        if tk > t_trip and i_curr <= 0.0:
            break

    return i_peak, i2t, v_sw_max


# ============================================================
# PARAMETROS
# ============================================================

I2T_LIMITE = 1500.0       # A^2s
V_SPIKE_LIMITE = 1200.0   # V
V_BUS = 400.0             # V
TAU_EXT_TESTE = 0.5e-6    # 0.5 us

# Cenario WCA
L_wca = 4.0e-6
C_wca = 1.5e-3

# Cenario nominal
L_nom = 1.8e-6
C_nom = 4.7e-3

# Varredura de tempo de desarme
t_trips = np.linspace(1e-6, 120e-6, 600)


# ============================================================
# 1. VARREDURA PARA O LIMITE ADIABATICO I^2t
# ============================================================

print("=" * 75)
print("1. VARREDURA PARA O LIMITE ADIABATICO (I^2t = 1500 A^2s)")
print("=" * 75)

i2t_wca = []
i2t_nom = []
v_spike_unclamp_wca = []
v_spike_unclamp_nom = []

for tt in t_trips:

    _, i2t_w, _ = simular_curto(
        tt, L_wca, C_wca, clamp_ativo=True
    )

    _, i2t_n, _ = simular_curto(
        tt, L_nom, C_nom, clamp_ativo=True
    )

    i2t_wca.append(i2t_w)
    i2t_nom.append(i2t_n)

    # Sem clamp: tau_ext = 0.5 us
    _, _, v_w = simular_curto(
        tt,
        L_wca,
        C_wca,
        clamp_ativo=False,
        tau_ext=TAU_EXT_TESTE
    )

    _, _, v_n = simular_curto(
        tt,
        L_nom,
        C_nom,
        clamp_ativo=False,
        tau_ext=TAU_EXT_TESTE
    )

    v_spike_unclamp_wca.append(v_w)
    v_spike_unclamp_nom.append(v_n)


i2t_wca = np.array(i2t_wca)
i2t_nom = np.array(i2t_nom)
v_spike_unclamp_wca = np.array(v_spike_unclamp_wca)
v_spike_unclamp_nom = np.array(v_spike_unclamp_nom)


# ------------------------------------------------------------
# Encontrar cruzamento do limite I^2t no WCA
# ------------------------------------------------------------

indices = np.where(i2t_wca >= I2T_LIMITE)[0]

if len(indices) == 0:
    t_crit_wca = None
    i_pk_crit_wca = None
    i2t_crit_wca = None

    print(
        "AVISO: o limite I^2t = 1500 A^2s nao foi atingido "
        "na faixa de tempos simulada."
    )

else:
    idx_crit_wca = indices[0]

    if idx_crit_wca == 0:
        t_crit_wca = float(t_trips[0])

    else:
        t_crit_wca = float(
            np.interp(
                I2T_LIMITE,
                [
                    i2t_wca[idx_crit_wca - 1],
                    i2t_wca[idx_crit_wca]
                ],
                [
                    t_trips[idx_crit_wca - 1],
                    t_trips[idx_crit_wca]
                ]
            )
        )

    i_pk_crit_wca, i2t_crit_wca, _ = simular_curto(
        t_crit_wca,
        L_wca,
        C_wca,
        clamp_ativo=True
    )

    print("\nCASO WCA (L = 4.0 uH, C = 1.5 mF):")

    print(
        f" -> Tempo critico interpolado para I^2t = 1500 A^2s: "
        f"t_trip = {t_crit_wca * 1e6:.2f} us"
    )

    print(
        f" -> Corrente de pico atingida: "
        f"{i_pk_crit_wca / 1e3:.2f} kA"
    )

    print(
        f" -> I^2t recalculado na simulacao: "
        f"{i2t_crit_wca:.2f} A^2s"
    )

    print(
        " -> Nota: o I^2t recalculado pode diferir de 1500 A^2s "
        "devido a interpolacao e discretizacao temporal."
    )


# ============================================================
# 2. TEMPO DE RETARDO QUE ATINGE V_SPIKE = 1200 V
# ============================================================

print("\n" + "=" * 75)
print(
    "2. TEMPO DE RETARDO QUE ATINGE V_spike = 1200 V "
    "(SEM CLAMP, tau_ext = 0.5 us)"
)
print("=" * 75)


def primeiro_cruzamento(valores, limite, tempos):
    """Retorna o tempo do primeiro cruzamento, ou None."""
    indices = np.where(valores >= limite)[0]

    if len(indices) == 0:
        return None

    idx = indices[0]

    if idx == 0:
        return float(tempos[0])

    # Interpolacao linear para melhorar a estimativa temporal
    return float(
        np.interp(
            limite,
            [valores[idx - 1], valores[idx]],
            [tempos[idx - 1], tempos[idx]]
        )
    )


t_v1200_wca = primeiro_cruzamento(
    v_spike_unclamp_wca,
    V_SPIKE_LIMITE,
    t_trips
)

t_v1200_nom = primeiro_cruzamento(
    v_spike_unclamp_nom,
    V_SPIKE_LIMITE,
    t_trips
)


if t_v1200_wca is not None:

    i_pk_v1200_wca, _, v_sw_1200_wca = simular_curto(
        t_v1200_wca,
        L_wca,
        C_wca,
        clamp_ativo=False,
        tau_ext=TAU_EXT_TESTE
    )

    print("\nCENARIO WCA (L = 4.0 uH):")

    print(
        f" -> Tempo de retardo t_trip: "
        f"{t_v1200_wca * 1e6:.3f} us "
        f"({t_v1200_wca * 1e9:.1f} ns)"
    )

    print(
        f" -> Corrente de pico na simulacao: "
        f"{i_pk_v1200_wca:.1f} A"
    )

    print(
        f" -> Tensao de pico V_sw: {v_sw_1200_wca:.1f} V"
    )

else:
    print(
        "\nCENARIO WCA: V_spike = 1200 V nao foi atingido "
        "na faixa simulada."
    )


if t_v1200_nom is not None:

    i_pk_v1200_nom, _, v_sw_1200_nom = simular_curto(
        t_v1200_nom,
        L_nom,
        C_nom,
        clamp_ativo=False,
        tau_ext=TAU_EXT_TESTE
    )

    print("\nCENARIO NOMINAL (L = 1.8 uH):")

    print(
        f" -> Tempo de retardo t_trip: "
        f"{t_v1200_nom * 1e6:.3f} us "
        f"({t_v1200_nom * 1e9:.1f} ns)"
    )

    print(
        f" -> Corrente de pico na simulacao: "
        f"{i_pk_v1200_nom:.1f} A"
    )

    print(
        f" -> Tensao de pico V_sw: {v_sw_1200_nom:.1f} V"
    )

else:
    print(
        "\nCENARIO NOMINAL: V_spike = 1200 V nao foi atingido "
        "na faixa simulada."
    )


# ============================================================
# 3. TEMPO MINIMO DE EXTINCAO (SOFT TURN-OFF)
# ============================================================

print("\n" + "=" * 75)
print(
    "3. TEMPO DE EXTINCAO MINIMO (SOFT TURN-OFF) "
    "PARA EVITAR V_spike > 1200 V"
)
print("=" * 75)

# Estimativa analitica:
# V_spike ~= V_bus + L * I / tau_ext
#
# Portanto:
# tau_ext >= L * I / (V_limite - V_bus)

delta_V_max = V_SPIKE_LIMITE - V_BUS


def calcular_tau_min(L_line, i_corrente):
    """Calcula tau minimo estimado, em segundos."""
    if delta_V_max <= 0.0:
        raise ValueError(
            "O limite de tensao deve ser superior a V_BUS."
        )

    return L_line * i_corrente / delta_V_max


# ------------------------------------------------------------
# Caso 1: nominal, desarme em 15 us
# ------------------------------------------------------------

i_pk_15us_nom, _, _ = simular_curto(
    15e-6, L_nom, C_nom, clamp_ativo=True
)

tau_min_15us_nom = calcular_tau_min(
    L_nom, i_pk_15us_nom
)

print(
    "\nCENARIO NOMINAL "
    "(t_trip = 15 us, L = 1.8 uH):"
)

print(
    f" -> Corrente de pico: {i_pk_15us_nom / 1e3:.3f} kA"
)

print(
    f" -> Tempo minimo de extincao estimado: "
    f"tau_ext >= {tau_min_15us_nom * 1e6:.3f} us"
)


# ------------------------------------------------------------
# Caso 2: WCA, desarme em 15 us
# ------------------------------------------------------------

i_pk_15us_wca, _, _ = simular_curto(
    15e-6, L_wca, C_wca, clamp_ativo=True
)

tau_min_15us_wca = calcular_tau_min(
    L_wca, i_pk_15us_wca
)

print(
    "\nCENARIO WCA "
    "(t_trip = 15 us, L = 4.0 uH):"
)

print(
    f" -> Corrente de pico: {i_pk_15us_wca / 1e3:.3f} kA"
)

print(
    f" -> Tempo minimo de extincao estimado: "
    f"tau_ext >= {tau_min_15us_wca * 1e6:.3f} us"
)


# ------------------------------------------------------------
# Caso 3: WCA, desarme em 100 us
# ------------------------------------------------------------

i_pk_100us_wca, _, _ = simular_curto(
    100e-6, L_wca, C_wca, clamp_ativo=True
)

tau_min_100us_wca = calcular_tau_min(
    L_wca, i_pk_100us_wca
)

print(
    "\nCENARIO WCA "
    "(t_trip = 100 us, L = 4.0 uH):"
)

print(
    f" -> Corrente de pico: {i_pk_100us_wca / 1e3:.3f} kA"
)

print(
    f" -> Tempo minimo de extincao estimado: "
    f"tau_ext >= {tau_min_100us_wca * 1e6:.3f} us"
)


# ------------------------------------------------------------
# Caso 4: WCA no tempo critico I^2t = 1500 A^2s
# ------------------------------------------------------------

if t_crit_wca is not None:

    tau_min_crit_wca = calcular_tau_min(
        L_wca, i_pk_crit_wca
    )

    print(
        "\nCASO WCA CRITICO "
        "(I^2t alvo = 1500 A^2s, L = 4.0 uH, C = 1.5 mF):"
    )

    print(
        f" -> Tempo critico de desarme: "
        f"t_trip = {t_crit_wca * 1e6:.3f} us"
    )

    print(
        f" -> Corrente de pico: "
        f"{i_pk_crit_wca / 1e3:.3f} kA"
    )

    print(
        f" -> I^2t recalculado: "
        f"{i2t_crit_wca:.2f} A^2s"
    )

    print(
        f" -> Tempo minimo de extincao estimado: "
        f"tau_ext >= {tau_min_crit_wca * 1e6:.3f} us"
    )

    V_spike_est_crit = (
        V_BUS + L_wca * i_pk_crit_wca / TAU_EXT_TESTE
    )

    print(
        f" -> V_spike estimado para tau_ext = 0.5 us: "
        f"{V_spike_est_crit:.1f} V"
    )

    print(
        f" -> Excesso sobre o limite de 1200 V: "
        f"{max(0.0, V_spike_est_crit - V_SPIKE_LIMITE):.1f} V"
    )

else:
    print(
        "\nCASO WCA CRITICO: nao calculado porque "
        "o limite I^2t nao foi atingido na varredura."
    )


# ------------------------------------------------------------
# Caso 5: WCA, desarme em 100 us e tau_ext = 0.5 us
# ------------------------------------------------------------

V_spike_est_100us = (
    V_BUS + L_wca * i_pk_100us_wca / TAU_EXT_TESTE
)

print(
    "\nWCA COM t_trip = 100 us E tau_ext = 0.5 us:"
)

print(
    f" -> V_spike estimado: {V_spike_est_100us:.1f} V"
)

print(
    f" -> Limite especificado: {V_SPIKE_LIMITE:.1f} V"
)

print(
    f" -> Excesso sobre o limite: "
    f"{max(0.0, V_spike_est_100us - V_SPIKE_LIMITE):.1f} V"
)


# ============================================================
# 4. GRAFICOS
# ============================================================

plt.figure(figsize=(10, 6))

plt.plot(
    t_trips * 1e6,
    i2t_wca,
    label="WCA: L = 4.0 uH, C = 1.5 mF"
)

plt.plot(
    t_trips * 1e6,
    i2t_nom,
    label="Nominal: L = 1.8 uH, C = 4.7 mF"
)

plt.axhline(
    I2T_LIMITE,
    linestyle="--",
    label="Limite I^2t = 1500 A^2s"
)

if t_crit_wca is not None:
    plt.axvline(
        t_crit_wca * 1e6,
        linestyle=":",
        label=f"Tempo critico WCA = {t_crit_wca*1e6:.2f} us"
    )

plt.xlabel("Tempo de desarme t_trip (us)")
plt.ylabel("Integral de Joule I^2t (A^2s)")
plt.title("Curto-circuito HVDC: limite adiabatico")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()


plt.figure(figsize=(10, 6))

plt.plot(
    t_trips * 1e6,
    v_spike_unclamp_wca,
    label="WCA: L = 4.0 uH"
)

plt.plot(
    t_trips * 1e6,
    v_spike_unclamp_nom,
    label="Nominal: L = 1.8 uH"
)

plt.axhline(
    V_SPIKE_LIMITE,
    linestyle="--",
    label="Limite V_spike = 1200 V"
)

if t_v1200_wca is not None:
    plt.axvline(
        t_v1200_wca * 1e6,
        linestyle=":",
        label=f"Cruzamento WCA = {t_v1200_wca*1e6:.2f} us"
    )

if t_v1200_nom is not None:
    plt.axvline(
        t_v1200_nom * 1e6,
        linestyle=":",
        label=f"Cruzamento nominal = {t_v1200_nom*1e6:.2f} us"
    )

plt.xlabel("Tempo de desarme t_trip (us)")
plt.ylabel("Tensao de pico estimada V_sw (V)")
plt.title("Curto-circuito HVDC: sobretensao sem clamp")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
