#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ==============================================================================
# SIMULAÇÃO 4.1: SIMULAÇÃO 4.1 — MODELO ORBITAL, ENERGÉTICO E 
# TÉRMICO CONCENTRADO
# Modelo exploratório in silico. Não é um modelo de voo validado.
# ==============================================================================
# Dissertação de Mestrado Profissional - ProPECAut / IFES
# Autor: Dafnhê Thor Castro de Martins Barros
# Revisão Técnica: Modelo Físico 3 Diodos, Eficiência Buck (nbuck = 0.91) e Hot Spot
# ==============================================================================

"""
SIMULAÇÃO 4.1 — MODELO ORBITAL, ENERGÉTICO E TÉRMICO CONCENTRADO

Modelo exploratório in silico. Não é um modelo de voo validado.

Convenções e limitações:
- A_RAD = 105 m² é a área emissora total adotada e fica fixa.
- O beta é calculado por aproximação a partir de data, LTAN e inclinação.
- A varredura de LTAN é feita numa grade; não identifica a distribuição
  real dos LTANs da constelação.
- A simulação transitória abrange até 60 órbitas, cerca de quatro dias.
- pv_input_kw é potência FV já disponível no barramento após Buck.
  ETA_BUCK não é aplicada novamente.
- P_BASE = 1000 W é carga elétrica agregada; seus componentes não são
  adicionados uma segunda vez.
- O aquecedor da bateria fica fora de P_BASE e é variável.
- 343,2 W de calor interno OISL entram no balanço térmico agregado.
- A redução proporcional da carga em caso de déficit não representa
  controle DVFS, NMPC ou TEMS implementado.
- As capacitâncias, resistências, propriedades radiativas e hipóteses
  de partição térmica precisam de rastreabilidade e validação.
"""

import math
import os

import csv
from dataclasses import dataclass
from datetime import datetime, timedelta

# Importação opcional para geração de gráficos
try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "saida_4_1")

# ==============================================================================
# 1. Parâmetros
# ==============================================================================

# Constantes físicas e geometria
SIGMA = 5.670374419e-8       # W/(m² K⁴)
G_SOL = 1361.0               # W/m²
Q_EARTH_IR = 234.0           # W/m²; hipótese nominal do modelo (ECSS handbook: ~237)
ALBEDO = 0.30                # hipótese de albedo médio

RE_KM = 6378.137             # km
RSUN_KM = 696340.0           # km
D_ES_KM = 1.496e8            # km
MU_EARTH_KM3_S2 = 398600.4418 # km³/s²

# Órbita SSO de referência
ALTITUDE_KM = 700.0
SSO_INCLINATION_DEG = 98.2   # provisória; confirmar pelo dimensionamento J2
T_ORB = 2.0 * math.pi * math.sqrt(
    (RE_KM + ALTITUDE_KM) ** 3 / MU_EARTH_KM3_S2
)

# Radiador: área emissora total de projeto, fixa
A_RAD = 105.0                # m²; não variar nesta configuração
ALPHA_RAD = 0.14             # absortividade solar EOL adotada
EPS_RAD = 0.88               # emissividade IR adotada
ETA_FIN = 0.945              # eficiência adotada; requer confirmação geométrica

# Propriedades térmicas agregadas, provisórias
RHO_AREAL = 4.5              # kg/m²
CP_RAD = 900.0               # J/(kg K)
C_BAT = 54360.0              # J/K
C_CHASSIS = 18000.0          # J/K
R_BAT_CH = 0.20              # K/W
R_CH_RAD = 0.0005            # K/W

# Eficiência dos conversores
ETA_DAB = 0.95
ETA_DAB_WCA = 0.90           # referência WCA; não usada no caso nominal
ETA_BUCK = 0.91              # informativa; não reaplicar à PV pós-Buck

# Cargas úteis
TI_NOM_SUN = 40000.0         # W úteis de TI em insolação
TI_NOM_ECL = 18000.0         # W úteis de TI no cenário de eclipse
N_DIES = 80                  # hipótese de arquitetura

# Potências FV já pós-Buck, disponíveis no barramento
PV_SCENARIOS_KW = {
    "85C": 46.14,
    "65C": 49.73,
    "28C": 54.18,
}
PV_CASE = "65C"

# Carga elétrica agregada P_BASE.
# Composição adotada: 200 W bomba + 350 W aviônica + 440 W OISL
# + 10 W de margem.
P_BASE = 1000.0                # W
P_PUMP_ELECTRICAL = 200.0      # W, incluída em P_BASE
P_AVIONICS_ELECTRICAL = 350.0  # W, incluída em P_BASE
P_OISL_ELECTRICAL = 440.0     # W, incluída em P_BASE
P_BASE_MARGIN = 10.0           # W

# Partição térmica provisória de P_BASE:
# bomba 200 W + aviônica 350 W + calor OISL 343,2 W + margem 10 W.
Q_BASE_THERMAL = 903.2         # W, hipótese agregada de calor interno
Q_OISL_INTERNAL_HEAT = 343.2   # W, parte térmica do OISL
P_OISL_OPTICAL = 96.8          # W, potência óptica que deixa o veículo

# Aquecedor: variável, separado de P_BASE
HEATER_SET_C = 15.0            # °C
HEATER_KP = 120.0              # W/K
HEATER_MAX_W = 600.0           # W

# BESS agregado
V_BAT = 400.0                  # V, aproximação agregada
E_BAT_KWH = 9.008              # kWh
SOC_MAX = 1.0
SOC_MIN = 0.20
R_BAT = 0.085                  # ohm, resistência agregada adotada
CN_AH = 22.52                  # Ah
C_RATE_MAX_CHARGE = 1.0        # C
P_MAX_CHARGE_W = CN_AH * C_RATE_MAX_CHARGE * V_BAT

# Condições iniciais
T_CHASSIS_INIT_C = 21.0
T_RAD_INIT_C = 20.0
T_BAT_INIT_C = 26.0

# Limites de análise adotados; não equivalem automaticamente a limites
# qualificados de hardware.
T_BAT_MIN_C = 10.0
T_BAT_MAX_C = 45.0
T_RAD_MAX_C = 55.0

# Controles numéricos
TOL_T_K = 0.01
TOL_SOC = 1e-5
TOL_ENVELOPE_K = 0.02
MAX_ORBITS = 60
SIM_DT_S = 1.0

# Data inicial usada na aproximação solar
START_DATE = datetime(2026, 1, 1, 0, 0, 0)

# Configuração da varredura LTAN
LTAN_STEP_HOURS = 0.5
ANNUAL_DAY_STEP = 5
PHASE_SAMPLES_PER_ORBIT = 120

# Erro de atitude nominal:
# a variação máxima deve ser avaliada separadamente na FMECA.
ATTITUDE_ERROR_DEG = 0.0


# ==============================================================================
# 2. Estrutura de resultados
# ==============================================================================

@dataclass
class Result:
    converged: bool
    envelope_converged: bool
    periodic_state_converged: bool
    operationally_feasible: bool
    orbits: int
    beta_description: str

    t_rad_min_c: float
    t_rad_max_c: float
    t_ch_min_c: float
    t_ch_max_c: float
    t_bat_min_c: float
    t_bat_max_c: float

    soc_min_last_orbit: float
    soc_min_run: float
    soc_final: float
    unmet_energy_j_last_orbit: float

    e_pv_available_j_last_orbit: float
    e_pv_used_j_last_orbit: float
    e_pv_curtailed_j_last_orbit: float
    e_load_served_j_last_orbit: float

    e_bat_charge_terminal_j_last_orbit: float
    e_bat_discharge_terminal_j_last_orbit: float
    e_bat_chemical_change_j_last_orbit: float
    e_battery_loss_j_last_orbit: float
    e_heater_j_last_orbit: float

    e_dab_loss_j_last_orbit: float
    e_ti_heat_j_last_orbit: float
    e_base_heat_j_last_orbit: float
    e_chassis_heat_j_last_orbit: float

    e_radiated_j_last_orbit: float
    e_ambient_j_last_orbit: float

    electrical_residual_j_last_orbit: float
    battery_energy_residual_j_last_orbit: float
    thermal_residual_j_last_orbit: float

    delta_t_state_k: float
    delta_soc_state: float
    delta_envelope_k: float


# ==============================================================================
# 3. Utilitários
# ==============================================================================

def clamp(value, lower, upper):
    return max(lower, min(upper, value))


def solar_declination_approx(date):
    day = date.timetuple().tm_yday
    hour = date.hour + date.minute / 60.0 + date.second / 3600.0
    fractional_day = day + (hour - 12.0) / 24.0
    gamma = math.radians(
        360.0 * (fractional_day - 81.0) / 365.0
    )
    return 23.44 * math.sin(gamma)


def beta_from_ltan(
    date,
    ltan_hours,
    inclination_deg=SSO_INCLINATION_DEG,
):
    if not 0.0 <= ltan_hours < 24.0:
        raise ValueError("LTAN deve estar em [0,24).")

    delta = math.radians(solar_declination_approx(date))
    inclination = math.radians(inclination_deg)
    delta_omega = math.radians(15.0 * (ltan_hours - 12.0))

    sin_beta = (
        math.sin(inclination)
        * math.cos(delta)
        * math.sin(delta_omega)
        + math.cos(inclination)
        * math.sin(delta)
    )

    return math.degrees(
        math.asin(clamp(sin_beta, -1.0, 1.0))
    )


def beta_for_orbit(beta_profile_deg, orbit_index):
    if isinstance(beta_profile_deg, (int, float)):
        return float(beta_profile_deg)

    if beta_profile_deg is None or len(beta_profile_deg) == 0:
        raise ValueError(
            "beta_profile_deg deve ser escalar ou sequência não vazia."
        )

    return float(
        beta_profile_deg[
            min(orbit_index, len(beta_profile_deg) - 1)
        ]
    )


# ==============================================================================
# 4. Geometria de eclipse
# ==============================================================================

def circle_overlap_area(radius_sun, radius_earth, separation):
    if separation >= radius_sun + radius_earth:
        return 0.0

    if separation <= abs(radius_earth - radius_sun):
        return math.pi * min(radius_sun, radius_earth) ** 2

    arg_sun = clamp(
        (
            separation**2
            + radius_sun**2
            - radius_earth**2
        ) / (2.0 * separation * radius_sun),
        -1.0,
        1.0,
    )

    arg_earth = clamp(
        (
            separation**2
            + radius_earth**2
            - radius_sun**2
        ) / (2.0 * separation * radius_earth),
        -1.0,
        1.0,
    )

    radicand = max(
        0.0,
        (-separation + radius_sun + radius_earth)
        * (separation + radius_sun - radius_earth)
        * (separation - radius_earth + radius_sun)
        * (separation + radius_earth + radius_sun),
    )

    return (
        radius_sun**2 * math.acos(arg_sun)
        + radius_earth**2 * math.acos(arg_earth)
        - 0.5 * math.sqrt(radicand)
    )


def solar_disk_fraction(theta_rad, beta_deg, altitude_km):
    radius_orbit_km = RE_KM + altitude_km

    alpha_earth = math.asin(
        clamp(RE_KM / radius_orbit_km, -1.0, 1.0)
    )
    alpha_sun = math.asin(
        clamp(RSUN_KM / D_ES_KM, -1.0, 1.0)
    )

    beta = math.radians(beta_deg)
    cos_gamma = -math.cos(beta) * math.cos(theta_rad)
    gamma = math.acos(clamp(cos_gamma, -1.0, 1.0))

    overlap = circle_overlap_area(
        alpha_sun,
        alpha_earth,
        gamma,
    )

    fraction = 1.0 - overlap / (math.pi * alpha_sun**2)
    return clamp(fraction, 0.0, 1.0)


# ==============================================================================
# 5. Fator de visão e ambiente radiativo
# ==============================================================================

def earth_view_factor(alt_km, phi_deg):
    H = (RE_KM + alt_km) / RE_KM

    if H <= 1.0:
        raise ValueError("A altitude deve ser positiva.")

    phi = math.radians(phi_deg)
    theta_m = math.asin(1.0 / H)
    root = math.sqrt(H**2 - 1.0)

    if phi <= math.pi / 2.0 - theta_m:
        return clamp(math.cos(phi) / H**2, 0.0, 1.0)

    if phi <= math.pi / 2.0 + theta_m:
        sin_phi = math.sin(phi)
        if abs(sin_phi) < 1e-12:
            return 0.0

        arg1 = clamp(root / (H * sin_phi), -1.0, 1.0)
        arg2 = clamp(
            root * math.cos(phi) / sin_phi,
            -1.0,
            1.0,
        )
        radicand = max(
            0.0,
            1.0 - H**2 * math.cos(phi)**2,
        )

        factor = (
            0.5
            - math.asin(arg1) / math.pi
            + math.cos(phi)
            / (math.pi * H**2)
            * math.acos(arg2)
            - root
            * math.sqrt(radicand)
            / (math.pi * H**2)
        )
        return clamp(factor, 0.0, 1.0)

    return 0.0


def orbital_geometry(theta_rad, beta_deg):
    beta = math.radians(beta_deg)

    radial_hat = (
        math.cos(theta_rad),
        math.sin(theta_rad),
        0.0,
    )

    sun_hat = (
        math.cos(beta),
        0.0,
        math.sin(beta),
    )

    nadir_hat = tuple(-value for value in radial_hat)
    radiator_hat = tuple(-value for value in sun_hat)

    cos_phi = sum(
        radiator_hat[i] * nadir_hat[i]
        for i in range(3)
    )
    phi_deg = math.degrees(
        math.acos(clamp(cos_phi, -1.0, 1.0))
    )

    mu_albedo = max(
        0.0,
        sum(radial_hat[i] * sun_hat[i] for i in range(3)),
    )

    mu_solar_nominal = max(
        0.0,
        sum(radiator_hat[i] * sun_hat[i] for i in range(3)),
    )

    return phi_deg, mu_albedo, mu_solar_nominal


def thermal_environment(
    theta_rad,
    beta_deg,
    altitude_km,
    albedo=ALBEDO,
    q_earth_ir_wm2=Q_EARTH_IR,
    solar_fraction=1.0,
    attitude_error_deg=0.0,
):
    phi_deg, mu_albedo, mu_solar_nominal = orbital_geometry(
        theta_rad,
        beta_deg,
    )

    f_earth = earth_view_factor(
        altitude_km,
        phi_deg,
    )

    mu_solar_error = max(
        0.0,
        math.sin(math.radians(abs(attitude_error_deg))),
    )
    mu_solar = max(mu_solar_nominal, mu_solar_error)

    q_solar_abs = (
        ALPHA_RAD
        * G_SOL
        * mu_solar
        * solar_fraction
    )

    q_albedo_abs = (
        ALPHA_RAD
        * G_SOL
        * albedo
        * f_earth
        * mu_albedo
    )

    q_ir_abs = EPS_RAD * q_earth_ir_wm2 * f_earth

    return {
        "f_earth": f_earth,
        "phi_deg": phi_deg,
        "mu_albedo": mu_albedo,
        "mu_solar": mu_solar,
        "q_solar_abs": q_solar_abs,
        "q_albedo_abs": q_albedo_abs,
        "q_ir_abs": q_ir_abs,
        "q_ambient_abs": (
            q_solar_abs
            + q_albedo_abs
            + q_ir_abs
        ),
    }


# ==============================================================================
# 6. Modelo energético da bateria
# ==============================================================================

def battery_loss(p_terminal_w):
    current_a = abs(p_terminal_w) / V_BAT
    return R_BAT * current_a**2


def max_discharge_for_energy(
    p_requested_w,
    energy_available_j,
    dt_s,
):
    if (
        p_requested_w <= 0.0
        or energy_available_j <= 0.0
        or dt_s <= 0.0
    ):
        return 0.0, 0.0

    def chemical_energy_required(power_w):
        return (power_w + battery_loss(power_w)) * dt_s

    if chemical_energy_required(p_requested_w) <= energy_available_j:
        return p_requested_w, battery_loss(p_requested_w)

    low = 0.0
    high = p_requested_w

    for _ in range(60):
        middle = 0.5 * (low + high)
        if chemical_energy_required(middle) <= energy_available_j:
            low = middle
        else:
            high = middle

    return low, battery_loss(low)


def max_charge_for_headroom(
    p_available_w,
    energy_headroom_j,
    dt_s,
):
    if (
        p_available_w <= 0.0
        or energy_headroom_j <= 0.0
        or dt_s <= 0.0
    ):
        return 0.0, 0.0, 0.0

    p_limit = min(p_available_w, P_MAX_CHARGE_W)

    def chemical_energy_stored(power_w):
        return max(
            0.0,
            power_w - battery_loss(power_w),
        ) * dt_s

    if chemical_energy_stored(p_limit) <= energy_headroom_j:
        loss_w = battery_loss(p_limit)
        return (
            p_limit,
            loss_w,
            max(0.0, p_limit - loss_w),
        )

    low = 0.0
    high = p_limit

    for _ in range(60):
        middle = 0.5 * (low + high)
        if chemical_energy_stored(middle) <= energy_headroom_j:
            low = middle
        else:
            high = middle

    loss_w = battery_loss(low)
    return (
        low,
        loss_w,
        max(0.0, low - loss_w),
    )

# ==============================================================================
# 8. Simulação de três nós (Corrigido para aceitar start_date)
# ==============================================================================

def simulate(
    ti_sun_w=TI_NOM_SUN,
    ti_eclipse_w=TI_NOM_ECL,
    pv_input_kw=PV_SCENARIOS_KW["65C"],
    beta_profile_deg=None,
    dt_s=SIM_DT_S,
    r_bat_ch=R_BAT_CH,
    r_ch_rad=R_CH_RAD,
    albedo=ALBEDO,
    q_earth_ir_wm2=Q_EARTH_IR,
    heater_set_c=HEATER_SET_C,
    heater_kp=HEATER_KP,
    heater_max_w=HEATER_MAX_W,
    attitude_error_deg=0.0,
    eta_dab=ETA_DAB,
    n_orbits_max=MAX_ORBITS,
    start_date=START_DATE,  # <--- Adicionado aqui
    record_history=False,
):
    if beta_profile_deg is None:
        raise ValueError(
            "Informe beta_profile_deg: escalar ou sequência por órbita."
        )

    if min(dt_s, A_RAD, r_bat_ch, r_ch_rad, n_orbits_max) <= 0.0:
        raise ValueError(
            "dt, área, resistências e número de órbitas devem ser positivos."
        )

    if min(ti_sun_w, ti_eclipse_w, pv_input_kw, heater_max_w) < 0.0:
        raise ValueError("Potências não podem ser negativas.")

    if heater_kp < 0.0 or not 0.0 < eta_dab <= 1.0:
        raise ValueError("heater_kp ou eta_dab inválido.")

    if not (0.0 <= SOC_MIN < SOC_MAX <= 1.0):
        raise ValueError("Limites de SoC inválidos.")

    cap_j = E_BAT_KWH * 3.6e6
    c_rad = A_RAD * RHO_AREAL * CP_RAD
    n_steps = max(1, int(round(T_ORB / dt_s)))
    dt = T_ORB / n_steps

    tr = 273.15 + T_RAD_INIT_C
    tc = 273.15 + T_CHASSIS_INIT_C
    tb = 273.15 + T_BAT_INIT_C
    soc = SOC_MAX

    previous_start = None
    previous_envelope = None
    stable_state_count = 0
    stable_envelope_count = 0

    global_soc_min = soc
    final_result = None

    profile_is_constant = isinstance(
        beta_profile_deg,
        (int, float),
    )

    if profile_is_constant:
        beta_description = (
            f"beta constante={float(beta_profile_deg):.4f} deg"
        )
    else:
        beta_description = (
            f"perfil beta com {len(beta_profile_deg)} valores"
        )

    history_records = []
    elapsed_total_s = 0.0

    for orbit in range(1, n_orbits_max + 1):
        beta_now = beta_for_orbit(
            beta_profile_deg,
            orbit - 1,
        )

        start_state = (tr, tc, tb, soc)
        thermal_energy_start = (
            c_rad * tr
            + C_CHASSIS * tc
            + C_BAT * tb
        )
        soc_start = soc

        tr_min = tc_min = tb_min = float("inf")
        tr_max = tc_max = tb_max = float("-inf")
        soc_min_orbit = soc

        e_pv_available = 0.0
        e_pv_used = 0.0
        e_pv_curtailed = 0.0
        e_load_served = 0.0
        e_bat_charge_terminal = 0.0
        e_bat_discharge_terminal = 0.0
        e_bat_loss = 0.0
        e_heater = 0.0
        e_dab_loss = 0.0
        e_ti_heat = 0.0
        e_base_heat = 0.0
        e_chassis_heat = 0.0
        e_radiated = 0.0
        e_ambient = 0.0
        e_unmet = 0.0

        for k in range(n_steps):
            theta = 2.0 * math.pi * (k * dt) / T_ORB

            solar_fraction = solar_disk_fraction(
                theta,
                beta_now,
                ALTITUDE_KM,
            )

            p_ti_requested = (
                ti_eclipse_w
                + (ti_sun_w - ti_eclipse_w) * solar_fraction
            )

            p_pv_available = (
                pv_input_kw * 1000.0 * solar_fraction
            )

            env = thermal_environment(
                theta_rad=theta,
                beta_deg=beta_now,
                altitude_km=ALTITUDE_KM,
                albedo=albedo,
                q_earth_ir_wm2=q_earth_ir_wm2,
                solar_fraction=solar_fraction,
                attitude_error_deg=attitude_error_deg,
            )

            tb_c = tb - 273.15
            p_heater_requested = min(
                heater_max_w,
                max(
                    0.0,
                    heater_kp * (heater_set_c - tb_c),
                ),
            )

            p_ti_bus_requested = p_ti_requested / eta_dab
            p_total_requested = (
                p_ti_bus_requested
                + P_BASE
                + p_heater_requested
            )

            p_pv_to_load = min(
                p_pv_available,
                p_total_requested,
            )

            p_surplus = max(
                0.0,
                p_pv_available - p_pv_to_load,
            )

            p_deficit = max(
                0.0,
                p_total_requested - p_pv_to_load,
            )

            p_bat_charge_terminal = 0.0
            p_bat_discharge_terminal = 0.0
            p_bat_loss_charge = 0.0
            p_bat_loss_discharge = 0.0
            p_bat_chemical_charge = 0.0

            if p_deficit > 0.0:
                energy_available = max(
                    0.0,
                    (soc - SOC_MIN) * cap_j,
                )

                (
                    p_bat_discharge_terminal,
                    p_bat_loss_discharge,
                ) = max_discharge_for_energy(
                    p_deficit,
                    energy_available,
                    dt,
                )

                chemical_draw = (
                    p_bat_discharge_terminal
                    + p_bat_loss_discharge
                ) * dt

                soc = max(
                    SOC_MIN,
                    soc - chemical_draw / cap_j,
                )

            elif p_surplus > 0.0:
                energy_headroom = max(
                    0.0,
                    (SOC_MAX - soc) * cap_j,
                )

                (
                    p_bat_charge_terminal,
                    p_bat_loss_charge,
                    p_bat_chemical_charge,
                ) = max_charge_for_headroom(
                    p_surplus,
                    energy_headroom,
                    dt,
                )

                soc = min(
                    SOC_MAX,
                    soc
                    + p_bat_chemical_charge * dt / cap_j,
                )

            p_bat_loss = (
                p_bat_loss_charge
                + p_bat_loss_discharge
            )

            p_load_served = (
                p_pv_to_load
                + p_bat_discharge_terminal
            )

            p_unmet = max(
                0.0,
                p_total_requested - p_load_served,
            )

            e_unmet += p_unmet * dt

            p_pv_used = (
                p_pv_to_load
                + p_bat_charge_terminal
            )

            p_pv_curtailed = max(
                0.0,
                p_pv_available - p_pv_used,
            )

            served_fraction = clamp(
                p_load_served / max(p_total_requested, 1e-12),
                0.0,
                1.0,
            )

            p_ti_served = p_ti_requested * served_fraction
            p_heater_served = p_heater_requested * served_fraction

            p_ti_bus_served = p_ti_served / eta_dab
            p_dab_loss = max(
                0.0,
                p_ti_bus_served - p_ti_served,
            )

            q_ti = p_ti_served
            q_base = Q_BASE_THERMAL * served_fraction
            q_dab = p_dab_loss

            q_chassis_source = q_ti + q_base + q_dab

            q_bat_to_chassis = (tb - tc) / r_bat_ch
            q_chassis_to_rad = (tc - tr) / r_ch_rad

            q_ambient = A_RAD * env["q_ambient_abs"]

            q_radiated = (
                A_RAD
                * ETA_FIN
                * EPS_RAD
                * SIGMA
                * (max(tr, 1.0) ** 4 - 3.0**4)
            )

            d_tb = (
                p_bat_loss
                + p_heater_served
                - q_bat_to_chassis
            ) / C_BAT

            d_tc = (
                q_chassis_source
                + q_bat_to_chassis
                - q_chassis_to_rad
            ) / C_CHASSIS

            d_tr = (
                q_chassis_to_rad
                + q_ambient
                - q_radiated
            ) / c_rad

            tb += d_tb * dt
            tc += d_tc * dt
            tr += d_tr * dt

            e_pv_available += p_pv_available * dt
            e_pv_used += p_pv_used * dt
            e_pv_curtailed += p_pv_curtailed * dt
            e_load_served += p_load_served * dt

            e_bat_charge_terminal += (
                p_bat_charge_terminal * dt
            )
            e_bat_discharge_terminal += (
                p_bat_discharge_terminal * dt
            )
            e_bat_loss += p_bat_loss * dt
            e_heater += p_heater_served * dt

            e_dab_loss += p_dab_loss * dt
            e_ti_heat += q_ti * dt
            e_base_heat += q_base * dt
            e_chassis_heat += q_chassis_source * dt

            e_radiated += q_radiated * dt
            e_ambient += q_ambient * dt

            tr_c = tr - 273.15
            tc_c = tc - 273.15
            tb_c_new = tb - 273.15

            tr_min = min(tr_min, tr_c)
            tr_max = max(tr_max, tr_c)
            tc_min = min(tc_min, tc_c)
            tc_max = max(tc_max, tc_c)
            tb_min = min(tb_min, tb_c_new)
            tb_max = max(tb_max, tb_c_new)

            soc_min_orbit = min(soc_min_orbit, soc)
            global_soc_min = min(
                global_soc_min if "global_soc_min" in locals() else SOC_MAX,
                soc,
            )

            if record_history:
                history_records.append({
                    "orbit": orbit,
                    "time_s": elapsed_total_s,
                    "t_rad_c": tr_c,
                    "t_chassis_c": tc_c,
                    "t_bat_c": tb_c_new,
                    "soc": soc,
                    "solar_fraction": solar_fraction,
                })

            elapsed_total_s += dt

        thermal_energy_end = (
            c_rad * tr
            + C_CHASSIS * tc
            + C_BAT * tb
        )

        delta_thermal = (
            thermal_energy_end - thermal_energy_start
        )

        electrical_residual = (
            e_pv_used
            + e_bat_discharge_terminal
            - e_load_served
            - e_bat_charge_terminal
        )

        battery_chemical_change = (
            (soc - soc_start) * cap_j
        )

        battery_energy_residual = (
            battery_chemical_change
            - (
                e_bat_charge_terminal
                - e_bat_discharge_terminal
                - e_bat_loss
            )
        )

        thermal_residual = (
            e_chassis_heat
            + e_bat_loss
            + e_heater
            + e_ambient
            - e_radiated
            - delta_thermal
        )

        end_state = (tr, tc, tb, soc)
        envelope = (
            tr_min,
            tr_max,
            tc_min,
            tc_max,
            tb_min,
            tb_max,
        )

        if previous_start is None:
            delta_t_state = float("inf")
            delta_soc_state = float("inf")
            delta_envelope = float("inf")
        else:
            delta_t_state = max(
                abs(end_state[i] - previous_start[i])
                for i in range(3)
            )
            delta_soc_state = abs(
                end_state[3] - previous_start[3]
            )
            delta_envelope = max(
                abs(envelope[i] - previous_envelope[i])
                for i in range(6)
            )

        periodic_candidate = (
            isinstance(beta_profile_deg, (int, float))
            and delta_t_state < TOL_T_K
            and delta_soc_state < TOL_SOC
        )

        operationally_feasible = (
            e_unmet < 1.0
            and soc > SOC_MIN + TOL_SOC
        )

        if periodic_candidate and operationally_feasible:
            stable_state_count += 1
        else:
            stable_state_count = 0

        envelope_candidate = (
            delta_envelope < TOL_ENVELOPE_K
        )

        if envelope_candidate and operationally_feasible:
            stable_envelope_count += 1
        else:
            stable_envelope_count = 0

        periodic_state_converged = (
            stable_state_count >= 2
        )
        envelope_converged = (
            stable_envelope_count >= 2
        )
        converged = (
            operationally_feasible
            and (
                periodic_state_converged
                or envelope_converged
            )
        )

        final_result = Result(
            converged=converged,
            envelope_converged=envelope_converged,
            periodic_state_converged=periodic_state_converged,
            operationally_feasible=operationally_feasible,
            orbits=orbit,
            beta_description=beta_description,

            t_rad_min_c=tr_min,
            t_rad_max_c=tr_max,
            t_ch_min_c=tc_min,
            t_ch_max_c=tc_max,
            t_bat_min_c=tb_min,
            t_bat_max_c=tb_max,

            soc_min_last_orbit=soc_min_orbit,
            soc_min_run=global_soc_min,
            soc_final=soc,
            unmet_energy_j_last_orbit=e_unmet,

            e_pv_available_j_last_orbit=e_pv_available,
            e_pv_used_j_last_orbit=e_pv_used,
            e_pv_curtailed_j_last_orbit=e_pv_curtailed,
            e_load_served_j_last_orbit=e_load_served,

            e_bat_charge_terminal_j_last_orbit=e_bat_charge_terminal,
            e_bat_discharge_terminal_j_last_orbit=e_bat_discharge_terminal,
            e_bat_chemical_change_j_last_orbit=battery_chemical_change,
            e_battery_loss_j_last_orbit=e_bat_loss,
            e_heater_j_last_orbit=e_heater,

            e_dab_loss_j_last_orbit=e_dab_loss,
            e_ti_heat_j_last_orbit=e_ti_heat,
            e_base_heat_j_last_orbit=e_base_heat,
            e_chassis_heat_j_last_orbit=e_chassis_heat,

            e_radiated_j_last_orbit=e_radiated,
            e_ambient_j_last_orbit=e_ambient,

            electrical_residual_j_last_orbit=electrical_residual,
            battery_energy_residual_j_last_orbit=battery_energy_residual,
            thermal_residual_j_last_orbit=thermal_residual,

            delta_t_state_k=delta_t_state,
            delta_soc_state=delta_soc_state,
            delta_envelope_k=delta_envelope,
        )

        if converged:
            break

        previous_start = start_state
        previous_envelope = envelope

    if record_history:
        return final_result, history_records
    return final_result


# ==============================================================================
# 9. Diagnósticos de aleta e temperatura quase-estática de junção
# ==============================================================================

def fin_efficiency_diagnostic():
    t_ref_k = 318.15
    h_rad = 4.0 * EPS_RAD * SIGMA * t_ref_k**3
    k_material = 1275.0  # W/(m K), hipótese
    thickness = 0.002    # m, hipótese
    length = 0.075       # m, hipótese

    m = math.sqrt(
        2.0 * h_rad / (k_material * thickness)
    )
    ml = m * length
    eta_fin_calculated = (
        math.tanh(ml) / ml if ml > 0.0 else 1.0
    )
    return h_rad, m, eta_fin_calculated


# ==============================================================================
# 10. Varredura aproximada de LTAN
# ==============================================================================

def rank_ltan_by_insolation(
    start_date,
    altitude_km=ALTITUDE_KM,
    inclination_deg=SSO_INCLINATION_DEG,
    ltan_step_hours=0.5,
    day_step=5,
    n_phase_samples=120,
):
    if ltan_step_hours <= 0.0:
        raise ValueError("ltan_step_hours deve ser positivo.")
    if day_step <= 0 or n_phase_samples <= 0:
        raise ValueError("day_step e n_phase_samples devem ser positivos.")

    n_ltan = int(round(24.0 / ltan_step_hours))
    candidates = []
    worst_beta_dates = {}

    for j in range(n_ltan):
        ltan = j * ltan_step_hours
        total_fraction = 0.0
        count = 0

        for day in range(0, 365, day_step):
            date = start_date + timedelta(days=day)
            beta = beta_from_ltan(
                date,
                ltan,
                inclination_deg,
            )

            for k in range(n_phase_samples):
                theta = (
                    2.0 * math.pi * k / n_phase_samples
                )
                total_fraction += solar_disk_fraction(
                    theta,
                    beta,
                    altitude_km,
                )
                count += 1

        candidates.append(
            (ltan, total_fraction / count)
        )

        daily_betas = []
        for day in range(365):
            date = start_date + timedelta(days=day)
            beta = beta_from_ltan(
                date,
                ltan,
                inclination_deg,
            )
            daily_betas.append((abs(beta), date))

        worst_beta_dates[ltan] = min(
            daily_betas,
            key=lambda item: item[0],
        )[1]

    worst = min(candidates, key=lambda item: item[1])
    best = max(candidates, key=lambda item: item[1])

    return {
        "worst_ltan_h": worst[0],
        "worst_mean_solar_fraction": worst[1],
        "best_ltan_h": best[0],
        "best_mean_solar_fraction": best[1],
        "worst_beta_dates": worst_beta_dates,
        "ranking": sorted(candidates, key=lambda item: item[1]),
    }


def build_beta_profile(
    start_date,
    ltan_hours,
    n_orbits,
    orbit_period_s_value=T_ORB,
    inclination_deg=SSO_INCLINATION_DEG,
):
    return [
        beta_from_ltan(
            start_date + timedelta(
                seconds=i * orbit_period_s_value
            ),
            ltan_hours,
            inclination_deg,
        )
        for i in range(n_orbits)
    ]


# ==============================================================================
# 11. Relatório
# ==============================================================================

def print_result(label, ltan, beta_profile, result):
    print(f"\n--- {label} ---")
    print(f"LTAN = {ltan:05.2f} h")
    print(
        "Beta do perfil (min/max) = "
        f"{min(beta_profile):.3f} / {max(beta_profile):.3f} deg"
    )
    print(f"Órbitas simuladas = {result.orbits}")
    print(
        "Convergência geral / envelope / estado periódico / "
        "viabilidade = "
        f"{result.converged} / "
        f"{result.envelope_converged} / "
        f"{result.periodic_state_converged} / "
        f"{result.operationally_feasible}"
    )

    print(
        f"Radiador min/max = "
        f"{result.t_rad_min_c:.2f} / {result.t_rad_max_c:.2f} °C"
    )
    print(
        f"Chassi min/max = "
        f"{result.t_ch_min_c:.2f} / {result.t_ch_max_c:.2f} °C"
    )
    print(
        f"Bateria min/max = "
        f"{result.t_bat_min_c:.2f} / {result.t_bat_max_c:.2f} °C"
    )

    print(
        f"SoC mínimo última órbita/global = "
        f"{100.0 * result.soc_min_last_orbit:.2f} / "
        f"{100.0 * result.soc_min_run:.2f} %"
    )
    print(
        f"SoC final = {100.0 * result.soc_final:.2f} %"
    )

    print(
        "Energia não atendida última órbita = "
        f"{result.unmet_energy_j_last_orbit / 3.6e6:.6f} kWh"
    )


# ==============================================================================
# 12. Funções de Salvamento e Plotagem
# ==============================================================================

def save_outputs(case_label, history_records, result):
    output_dir = OUT_DIR
    os.makedirs(output_dir, exist_ok=True)
    
    safe_label = case_label.lower().replace(" ", "_")
    csv_filename = os.path.join(output_dir, f"simulacao_{safe_label}.csv")

    # 1. Salvar CSV com a série temporal
    if history_records:
        keys = history_records[0].keys()
        with open(csv_filename, mode="w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            for row in history_records:
                writer.writerow(row)
        print(f"-> Arquivo CSV salvo com sucesso em: {csv_filename}")

    # 2. Gerar Gráficos PNG (caso matplotlib esteja instalado)
    if HAS_MATPLOTLIB and history_records:
        times_h = [r["time_s"] / 3600.0 for r in history_records]
        t_rad = [r["t_rad_c"] for r in history_records]
        t_ch = [r["t_chassis_c"] for r in history_records]
        t_bat = [r["t_bat_c"] for r in history_records]
        socs = [100.0 * r["soc"] for r in history_records]

        # Gráfico de Temperaturas
        plt.figure(figsize=(10, 5))
        plt.plot(times_h, t_rad, label="Radiador (°C)", color="tab:blue")
        plt.plot(times_h, t_ch, label="Chassi (°C)", color="tab:orange")
        plt.plot(times_h, t_bat, label="Bateria (°C)", color="tab:green")
        plt.title(f"Evolução Térmica — {case_label}")
        plt.xlabel("Tempo de Simulação (horas)")
        plt.ylabel("Temperatura (°C)")
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.legend()
        plt.tight_layout()
        plot_temp_path = os.path.join(output_dir, f"temperaturas_{safe_label}.png")
        plt.savefig(plot_temp_path, dpi=300)
        plt.close()
        print(f"-> Gráfico de temperaturas salvo em: {plot_temp_path}")

        # Gráfico de Estado de Carga (SoC)
        plt.figure(figsize=(10, 4))
        plt.plot(times_h, socs, label="SoC (%)", color="tab:purple")
        plt.axhline(y=SOC_MIN * 100, color="r", linestyle=":", label=f"SoC Min ({SOC_MIN*100}%)")
        plt.title(f"Estado de Carga da Bateria (SoC) — {case_label}")
        plt.xlabel("Tempo de Simulação (horas)")
        plt.ylabel("SoC (%)")
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.legend()
        plt.tight_layout()
        plot_soc_path = os.path.join(output_dir, f"soc_{safe_label}.png")
        plt.savefig(plot_soc_path, dpi=300)
        plt.close()
        print(f"-> Gráfico de SoC salvo em: {plot_soc_path}")
    elif not HAS_MATPLOTLIB:
        print("-> Matplotlib não encontrado. Apenas os dados CSV foram salvos.")


# ==============================================================================
# 13. Execução principal
# ==============================================================================

def main():
    print("SIMULAÇÃO 4.1 — extremos aproximados de LTAN/insolação")
    print(
        f"Altitude={ALTITUDE_KM:.1f} km; "
        f"período={T_ORB:.2f} s; "
        f"inclinação provisória={SSO_INCLINATION_DEG:.2f}°"
    )
    print(
        f"Data inicial={START_DATE.isoformat()}; "
        f"passo={SIM_DT_S:g} s; "
        f"PV pós-Buck={PV_SCENARIOS_KW[PV_CASE]:.2f} kW"
    )
    print(f"Área emissora total fixa: {A_RAD:.1f} m²")

    h_fin, m_fin, eta_fin_calc = fin_efficiency_diagnostic()
    print(
        "Diagnóstico separado de aleta: "
        f"h_rad={h_fin:.3f} W/(m² K); "
        f"m={m_fin:.3f} 1/m; "
        f"eta calculada={eta_fin_calc:.5f}; "
        f"eta adotada={ETA_FIN:.3f}"
    )

    ranking = rank_ltan_by_insolation(START_DATE)

    print("\nVarredura anual aproximada de LTAN:")
    print(
        f"Menor fração solar média: LTAN "
        f"{ranking['worst_ltan_h']:.2f} h; "
        f"score={ranking['worst_mean_solar_fraction']:.6f}"
    )
    print(
        f"Maior fração solar média: LTAN "
        f"{ranking['best_ltan_h']:.2f} h; "
        f"score={ranking['best_mean_solar_fraction']:.6f}"
    )

    cases = (
        ("Pior Insolacao", "worst_ltan_h"),
        ("Melhor Insolacao", "best_ltan_h"),
    )

    for label, key in cases:
        ltan = ranking[key]
        seasonal_start = ranking["worst_beta_dates"][ltan]

        print(
            f"\nÉpoca selecionada para {label.lower()}: "
            f"{seasonal_start.date()}"
        )

        beta_profile = build_beta_profile(
            seasonal_start,
            ltan,
            MAX_ORBITS,
        )

        result, history_records = simulate(
            pv_input_kw=PV_SCENARIOS_KW[PV_CASE],
            beta_profile_deg=beta_profile,
            dt_s=SIM_DT_S,
            n_orbits_max=MAX_ORBITS,
            start_date=seasonal_start,
            record_history=True,
        )

        print_result(
            label,
            ltan,
            beta_profile,
            result,
        )

        # Salva os resultados da simulação atual (CSV e Gráficos PNG) na pasta saida_4_1
        save_outputs(label, history_records, result)

    print("\nProcesso concluído com sucesso. Todos os arquivos foram direcionados para a pasta 'saida_4_1'.")


if __name__ == "__main__":
    main()