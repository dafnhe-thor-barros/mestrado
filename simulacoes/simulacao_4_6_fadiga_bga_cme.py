#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# ==============================================================================
# SIMULAÇÃO 4.6 — COMPARAÇÃO DE FADIGA BGA, SSPC E PUE ORBITAL
# ==============================================================================
# Dissertação de Mestrado Profissional - ProPECAut / IFES
# Autor: Dafnhê Thor Castro de Martins Barros
# Revisão Técnica: Modelo Físico 3 Diodos, Eficiência Buck (nbuck = 0.91) e Hot Spot
# ==============================================================================
"""

Triagem in silico. Não é qualificação de hardware nem previsão
certificada de vida.

O valor TWO_EPS_F_REFERENCE = 2.887 é uma hipótese otimista de
capacidade mecânica teórica superior para o conjunto com underfill,
condicionada a ΔT_j <= 7 °C com TEMS. A retenção por 10 anos requer
aprovação do polímero em desgaseificação (ECSS-Q-ST-70-02C:
TML < 1,0%, CVCM < 0,1%) e validação experimental em TVAC.

O script:
1. Compara Coffin–Manson com Engelmaier usando a mesma estimativa
   de deformação da junta.
2. Compara os cenários térmicos com e sem TEMS.
3. Calcula envelopes de energia para SSPC, distinguindo-os da energia
   real do dispositivo.
4. Calcula um indicador PUE orbital adaptado pela razão entre energias
   integradas, com P_BASE e aquecimento configuráveis por fase orbital.

Os parâmetros marcados como hipótese ou não validados devem ser
substituídos por dados rastreáveis antes de interpretar os resultados
como previsão de missão.
"""

import csv
import json
import math
from pathlib import Path


# ---------------------------------------------------------------------
# Saídas
# ---------------------------------------------------------------------

OUTPUT_DIR = Path("saida_4_6")


# ---------------------------------------------------------------------
# Parâmetros do modelo de fadiga
# ---------------------------------------------------------------------

L_X_M = 28e-3                 # [HIP] Dimensão do encapsulamento
L_Y_M = 28e-3                 # [HIP] Dimensão do encapsulamento
H_S_M = 0.18e-3               # [HIP] Altura efetiva da junta
DELTA_ALPHA_K_INV = 14e-6     # [HIP] Diferença efetiva de CTE
F_CONF = 0.90                 # [HIP] Fator de conformidade
K_ATT = 0.85                  # [HIP] Razão ΔT_solda/ΔT_junção

# Hipótese otimista de limite superior; não validada para a montagem
# SAC305/CoWoS/underfill específica.
TWO_EPS_F_REFERENCE = 2.887

# Varredura exploratória; não é uma faixa qualificada de SAC305.
TWO_EPS_F_SWEEP = [0.65, 1.0, 2.887]

# Coffin–Manson com expoente fixo, apenas como caso de sensibilidade.
C_COFFIN_MANSON = -0.5

# Coeficientes previamente usados na dissertação.
# Tratados como hipótese/calibração não validada.
ENGELMAIER_THESIS = {
    "a": -0.470,
    "b": -7.0e-4,
    "d": 2.0e-2,
    "status": "HIPOTESE_CALIBRACAO_NAO_VALIDADA",
}

# Forma de Engelmaier usada nos scripts anteriores como sensibilidade.
# Confirmar a fonte primária e sua aplicabilidade ao empilhamento real.
ENGELMAIER_REFERENCE = {
    "a": -0.442,
    "b": -6.0e-4,
    "d": 1.74e-2,
    "status": "FONTE_E_APLICABILIDADE_A_CONFIRMAR",
}

LIFE_TARGET_YEARS = 10.0


# ---------------------------------------------------------------------
# Cenários orbitais e térmicos
# ---------------------------------------------------------------------

# Períodos e eclipses reproduzem os valores dos scripts enviados.
# As variações térmicas são entradas herdadas: confirmar sua origem
# quando os resultados orbitais finais forem reconciliados.
ORBIT_CASES = [
    {
        "orbit": "SSO",
        "period_s": 5926.0,
        "eclipse_s": 1128.0,
        "delta_tj_tems_c": 7.0,
        "delta_tj_no_tems_c": 22.5,
        "mean_solder_temp_c": 60.0,
        "ti_sun_w": 40000.0,
        "ti_eclipse_tems_w": 18000.0,
        "ti_eclipse_no_tems_w": 20000.0,
    },
    {
        "orbit": "LEO",
        "period_s": 5677.0,
        "eclipse_s": 2148.0,
        "delta_tj_tems_c": 18.2,
        "delta_tj_no_tems_c": 48.0,
        "mean_solder_temp_c": 50.0,
        "ti_sun_w": 40000.0,
        "ti_eclipse_tems_w": 6500.0,
        "ti_eclipse_no_tems_w": 20000.0,
    },
]


# ---------------------------------------------------------------------
# Parâmetros do PUE orbital adaptado
# ---------------------------------------------------------------------

ETA_DAB = 0.95

# Auxiliares elétricos distintos: não confundir consumo OISL com seu
# calor residual interno. O calor residual é parâmetro térmico, não
# uma segunda parcela de consumo elétrico.
P_PUMP_W = 200.0
P_AVIONICS_W = 350.0
P_OISL_ELECTRICAL_W = 440.0
Q_OISL_INTERNAL_HEAT_W = 343.2  # Usado no balanço térmico, não somado ao PUE.

# P_BASE é provisório e configurável por fase orbital.
# Verificar no pente-fino quais cargas contém antes da versão final.
P_BASE_SUN_W = 1000.0
P_BASE_ECLIPSE_W = 1000.0

# Potência elétrica do aquecedor, também separada de P_BASE.
# Valores atuais são hipóteses configuráveis; usar potência efetiva
# média por fase, não automaticamente a potência máxima instalada.
P_BATTERY_HEATER_SUN_W = 0.0
P_BATTERY_HEATER_ECLIPSE_W = 0.0

# Perdas BESS usadas como hipótese nos scripts anteriores.
# Não são recalculadas por este indicador e devem ser verificadas
# contra o modelo eletrotérmico escolhido.
P_BESS_LOSS_W = 65.0

# Número e desempenho teórico dos dies: hipótese, não throughput medido.
N_DIES = 80
TOPS_PER_DIE = 250.0


# ---------------------------------------------------------------------
# Parâmetros SSPC
# ---------------------------------------------------------------------

V_BUS_V = 400.0
P_NODE_W = 50000.0
I_TRIP_MULTIPLIER = 1.25
T_TRIP_S = 15e-6
L_HARNESS_H = 2.0e-6
V_TVS_CLAMP_V = 644.0

# Limite ilustrativo: substituir pelo componente selecionado e suas
# condições de operação. Não é inferido apenas da tensão nominal.
V_SWITCH_ABS_MAX_V = 1200.0


# ---------------------------------------------------------------------
# Funções comuns
# ---------------------------------------------------------------------

def cycles_per_day(period_s):
    if period_s <= 0:
        raise ValueError("O período orbital precisa ser positivo.")
    return 86400.0 / period_s


def strain_range(delta_tj_c):
    """
    Aproximação geométrica comum aos modelos comparados:

    Δγ = F_conf * L_DNP * Δα * ΔT_s / h_s
    ΔT_s = K_att * ΔT_j

    Esta estimativa não substitui análise mecânica da montagem.
    """
    if delta_tj_c <= 0:
        raise ValueError("ΔT_j precisa ser positivo.")

    l_dnp = 0.5 * math.hypot(L_X_M, L_Y_M)
    delta_ts_c = K_ATT * delta_tj_c
    gamma = (
        F_CONF
        * l_dnp
        * DELTA_ALPHA_K_INV
        * delta_ts_c
        / H_S_M
    )
    return l_dnp, delta_ts_c, gamma


def engelmaier_exponent(temp_c, frequency_cpd, coefficients):
    return (
        coefficients["a"]
        + coefficients["b"] * temp_c
        + coefficients["d"] * math.log1p(frequency_cpd)
    )


def cycles_to_failure(gamma, two_eps_f, exponent):
    """
    Forma comparativa usada nos scripts recebidos:

    Nf = 0.5 * (γ / (2ε'f)) ** (1/c)

    O parâmetro 2ε'f e a forma de estimar γ precisam de validação
    para o material e o empilhamento reais.
    """
    if gamma <= 0 or two_eps_f <= 0:
        raise ValueError("γ e 2ε'f precisam ser positivos.")
    if exponent >= 0:
        raise ValueError("O expoente c precisa ser negativo.")

    return 0.5 * (gamma / two_eps_f) ** (1.0 / exponent)


# ---------------------------------------------------------------------
# Fadiga: Coffin–Manson e Engelmaier
# ---------------------------------------------------------------------

def fatigue_analysis():
    results = []
    sensitivity = []
    gains = []

    for orbit_case in ORBIT_CASES:
        frequency = cycles_per_day(orbit_case["period_s"])
        tmean = orbit_case["mean_solder_temp_c"]

        life_by_control = {}

        control_cases = [
            ("TEMS", orbit_case["delta_tj_tems_c"]),
            ("Sem TEMS", orbit_case["delta_tj_no_tems_c"]),
        ]

        for control_name, delta_tj in control_cases:
            l_dnp, delta_ts, gamma = strain_range(delta_tj)

            models = [
                {
                    "name": "Coffin-Manson, c fixo",
                    "exponent": C_COFFIN_MANSON,
                    "parameter_status": "HIPOTESE_DE_SENSIBILIDADE",
                },
                {
                    "name": "Engelmaier, coeficientes da dissertacao",
                    "exponent": engelmaier_exponent(
                        tmean, frequency, ENGELMAIER_THESIS
                    ),
                    "parameter_status": ENGELMAIER_THESIS["status"],
                },
                {
                    "name": "Engelmaier, forma de referencia",
                    "exponent": engelmaier_exponent(
                        tmean, frequency, ENGELMAIER_REFERENCE
                    ),
                    "parameter_status": ENGELMAIER_REFERENCE["status"],
                },
            ]

            for model in models:
                nf = cycles_to_failure(
                    gamma,
                    TWO_EPS_F_REFERENCE,
                    model["exponent"],
                )
                life_years = nf / (frequency * 365.25)

                row = {
                    "orbit": orbit_case["orbit"],
                    "control": control_name,
                    "model": model["name"],
                    "delta_tj_c": delta_tj,
                    "delta_tsolder_c": delta_ts,
                    "mean_solder_temp_c": tmean,
                    "cycles_per_day": frequency,
                    "l_dnp_m": l_dnp,
                    "gamma_range": gamma,
                    "c_exponent": model["exponent"],
                    "two_eps_f": TWO_EPS_F_REFERENCE,
                    "cycles_to_failure": nf,
                    "life_years": life_years,
                    "margin_vs_10y_percent": (
                        100.0
                        * (life_years - LIFE_TARGET_YEARS)
                        / LIFE_TARGET_YEARS
                    ),
                    "parameter_status": model["parameter_status"],
                    "interpretation": (
                        "Triagem in silico; nao e vida qualificada."
                    ),
                }
                results.append(row)
                life_by_control[(control_name, model["name"])] = life_years

            # Sensibilidade ao parâmetro de ductilidade.
            for two_eps in TWO_EPS_F_SWEEP:
                for model in models:
                    nf_sens = cycles_to_failure(
                        gamma,
                        two_eps,
                        model["exponent"],
                    )
                    sensitivity.append({
                        "orbit": orbit_case["orbit"],
                        "control": control_name,
                        "model": model["name"],
                        "two_eps_f_sensitivity": two_eps,
                        "gamma_range": gamma,
                        "c_exponent": model["exponent"],
                        "cycles_to_failure": nf_sens,
                        "life_years": nf_sens / (frequency * 365.25),
                        "status": "SENSIBILIDADE_NAO_QUALIFICADA",
                    })

        model_names = [
            "Coffin-Manson, c fixo",
            "Engelmaier, coeficientes da dissertacao",
            "Engelmaier, forma de referencia",
        ]

        for model_name in model_names:
            life_tems = life_by_control[("TEMS", model_name)]
            life_no_tems = life_by_control[("Sem TEMS", model_name)]
            gains.append({
                "orbit": orbit_case["orbit"],
                "model": model_name,
                "life_with_tems_years": life_tems,
                "life_without_tems_years": life_no_tems,
                "life_ratio_tems_over_no_tems": (
                    life_tems / life_no_tems
                    if life_no_tems > 0 else None
                ),
                "life_gain_percent": (
                    100.0 * (life_tems / life_no_tems - 1.0)
                    if life_no_tems > 0 else None
                ),
                "warning": (
                    "Ganho condicionado aos ΔT e parâmetros assumidos; "
                    "nao e validacao experimental."
                ),
            })

    return results, sensitivity, gains


# ---------------------------------------------------------------------
# SSPC: envelopes condicionais e energia idealizada de desmagnetização
# ---------------------------------------------------------------------

def sspc_analysis():
    if ETA_DAB <= 0 or V_TVS_CLAMP_V <= V_BUS_V:
        raise ValueError("Eficiência ou tensão de clamp inválida.")

    i_nom = P_NODE_W / V_BUS_V
    i_limit = I_TRIP_MULTIPLIER * i_nom

    # Limite retangular condicional, se tensão e corrente permanecerem
    # limitadas pelos valores assumidos durante todo o tempo de trip.
    e_switch_rectangular_j = (
        V_TVS_CLAMP_V * i_limit * T_TRIP_S
    )

    # Energia inicialmente armazenada na indutância do chicote.
    e_inductor_j = 0.5 * L_HARNESS_H * i_limit**2

    # Modelo idealizado de clamp, assumindo corrente que decresce
    # linearmente por (Vclamp - Vbus) sobre a indutância.
    delta_v = V_TVS_CLAMP_V - V_BUS_V
    t_demag_s = L_HARNESS_H * i_limit / delta_v

    # Energia absorvida pelo TVS no modelo idealizado.
    # Inclui energia que a fonte pode fornecer durante a desmagnetização.
    e_tvs_ideal_j = (
        e_inductor_j * V_TVS_CLAMP_V / delta_v
    )
    e_source_during_demag_j = (
        0.5 * V_BUS_V * i_limit * t_demag_s
    )

    return {
        "i_nom_a": i_nom,
        "i_limit_assumed_a": i_limit,
        "trip_time_us_assumed": T_TRIP_S * 1e6,
        "switch_energy_rectangular_bound_j": e_switch_rectangular_j,
        "inductor_energy_j": e_inductor_j,
        "clamp_voltage_assumed_v": V_TVS_CLAMP_V,
        "ideal_demagnetization_time_us": t_demag_s * 1e6,
        "ideal_tvs_energy_j": e_tvs_ideal_j,
        "source_energy_during_demag_j": e_source_during_demag_j,
        "semiconductor_abs_max_voltage_assumed_v": V_SWITCH_ABS_MAX_V,
        "nominal_voltage_margin_v": V_SWITCH_ABS_MAX_V - V_TVS_CLAMP_V,
        "qualification_status": (
            "NAO QUALIFICADO: exige componente selecionado, curva I-V "
            "do clamp, parasitas, forma de onda e SOA."
        ),
    }


# ---------------------------------------------------------------------
# PUE orbital adaptado: razão entre energias integradas
# ---------------------------------------------------------------------

def pue_for_case(orbit_case, tems_enabled):
    period = orbit_case["period_s"]
    eclipse = orbit_case["eclipse_s"]

    if not (0 <= eclipse < period):
        raise ValueError("Duração de eclipse incompatível com o período.")

    sun_time = period - eclipse
    eclipse_fraction = eclipse / period
    sun_fraction = sun_time / period

    ti_sun = orbit_case["ti_sun_w"]

    if tems_enabled:
        ti_eclipse = orbit_case["ti_eclipse_tems_w"]
        control_name = "TEMS"
    else:
        ti_eclipse = orbit_case["ti_eclipse_no_tems_w"]
        control_name = "Sem TEMS"

    # Energias úteis de TI em uma órbita.
    e_it_sun = ti_sun * sun_time
    e_it_eclipse = ti_eclipse * eclipse
    e_it_total = e_it_sun + e_it_eclipse

    # Perdas DAB: a potência útil de TI é considerada saída do DAB.
    e_dab_loss = (
        ti_sun * (1.0 / ETA_DAB - 1.0) * sun_time
        + ti_eclipse * (1.0 / ETA_DAB - 1.0) * eclipse
    )

    # Cargas variáveis por fase orbital.
    e_base = (
        P_BASE_SUN_W * sun_time
        + P_BASE_ECLIPSE_W * eclipse
    )
    e_heater = (
        P_BATTERY_HEATER_SUN_W * sun_time
        + P_BATTERY_HEATER_ECLIPSE_W * eclipse
    )

    # Cargas elétricas auxiliares e OISL.
    p_aux_electrical = (
        P_PUMP_W + P_AVIONICS_W + P_OISL_ELECTRICAL_W
    )
    e_aux = p_aux_electrical * period

    # Perda BESS como hipótese independente. Revisar a fronteira e
    # evitar somá-la se já estiver embutida em P_BASE ou noutra parcela.
    e_bess_loss = P_BESS_LOSS_W * period

    e_total = (
        e_it_total
        + e_dab_loss
        + e_base
        + e_heater
        + e_aux
        + e_bess_loss
    )

    if e_it_total <= 0:
        raise ValueError("Energia de TI precisa ser positiva.")

    pue = e_total / e_it_total
    tops_peak = N_DIES * TOPS_PER_DIE

    # Indicador de pico nominal; não é throughput medido.
    tops_per_w = tops_peak / (e_total / period)

    return {
        "orbit": orbit_case["orbit"],
        "control": control_name,
        "period_s": period,
        "eclipse_s": eclipse,
        "sun_fraction": sun_fraction,
        "eclipse_fraction": eclipse_fraction,
        "ti_sun_w": ti_sun,
        "ti_eclipse_w": ti_eclipse,
        "p_base_sun_w_provisional": P_BASE_SUN_W,
        "p_base_eclipse_w_provisional": P_BASE_ECLIPSE_W,
        "heater_sun_w_assumed": P_BATTERY_HEATER_SUN_W,
        "heater_eclipse_w_assumed": P_BATTERY_HEATER_ECLIPSE_W,
        "oisl_electrical_w": P_OISL_ELECTRICAL_W,
        "oisl_internal_heat_w_separate": Q_OISL_INTERNAL_HEAT_W,
        "p_aux_electrical_w": p_aux_electrical,
        "bess_loss_w_assumed": P_BESS_LOSS_W,
        "eta_dab_assumed": ETA_DAB,
        "it_energy_j": e_it_total,
        "dab_loss_energy_j": e_dab_loss,
        "base_energy_j": e_base,
        "heater_energy_j": e_heater,
        "auxiliary_energy_j": e_aux,
        "bess_loss_energy_j": e_bess_loss,
        "total_boundary_energy_j": e_total,
        "pue_orbital_adapted": pue,
        "tops_peak_total": tops_peak,
        "tops_per_w_nominal_peak": tops_per_w,
        "warning": (
            "Indicador adaptado, nao alegacao de conformidade ISO. "
            "P_BASE e aquecimento sao configuraveis; verificar "
            "sobreposicao antes de interpretar."
        ),
    }


def pue_analysis():
    rows = []
    gains = []

    for orbit_case in ORBIT_CASES:
        tems = pue_for_case(orbit_case, True)
        no_tems = pue_for_case(orbit_case, False)
        rows.extend([tems, no_tems])

        gains.append({
            "orbit": orbit_case["orbit"],
            "pue_with_tems": tems["pue_orbital_adapted"],
            "pue_without_tems": no_tems["pue_orbital_adapted"],
            "pue_change_percent": (
                100.0
                * (
                    tems["pue_orbital_adapted"]
                    / no_tems["pue_orbital_adapted"]
                    - 1.0
                )
            ),
            "interpretation": (
                "Compara os perfis de carga assumidos; nao isola "
                "causalidade se as cargas auxiliares ou P_BASE mudarem."
            ),
        })

    return rows, gains


# ---------------------------------------------------------------------
# Exportação
# ---------------------------------------------------------------------

def write_csv(path, rows):
    if not rows:
        return

    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    fatigue, fatigue_sensitivity, fatigue_gains = fatigue_analysis()
    pue_rows, pue_gains = pue_analysis()
    sspc = sspc_analysis()

    write_csv(
        OUTPUT_DIR / "simulacao_4_6_fadiga.csv",
        fatigue,
    )
    write_csv(
        OUTPUT_DIR / "simulacao_4_6_fadiga_sensibilidade.csv",
        fatigue_sensitivity,
    )
    write_csv(
        OUTPUT_DIR / "simulacao_4_6_ganho_tems_fadiga.csv",
        fatigue_gains,
    )
    write_csv(
        OUTPUT_DIR / "simulacao_4_6_pue_orbital.csv",
        pue_rows,
    )
    write_csv(
        OUTPUT_DIR / "simulacao_4_6_ganho_tems_pue.csv",
        pue_gains,
    )

    payload = {
        "script": "simulacao_4_6_comparacao.py",
        "status": (
            "triagem_in_silico_nao_qualificada"
        ),
        "fatigue_results": fatigue,
        "fatigue_sensitivity": fatigue_sensitivity,
        "tems_fatigue_gain": fatigue_gains,
        "sspc": sspc,
        "pue_results": pue_rows,
        "tems_pue_comparison": pue_gains,
        "pue_boundary": {
            "included": [
                "potencia util de TI",
                "perdas DAB calculadas a partir de ETA_DAB",
                "P_BASE configuravel por fase orbital",
                "aquecimento da bateria configuravel por fase",
                "bomba",
                "avionica",
                "consumo eletrico integral dos OISL",
                "perda BESS hipotetica",
            ],
            "separate_thermal_quantity": (
                "Q_OISL_INTERNAL_HEAT_W: calor interno OISL; "
                "nao somar novamente como consumo eletrico."
            ),
            "limitations": [
                "P_BASE e provisório e pode conter cargas já contabilizadas.",
                "A perda BESS pode sobrepor-se a parcelas de P_BASE.",
                "O aquecedor deve usar potência média efetivamente acionada.",
                "O indicador é adaptado e não declara conformidade normativa.",
            ],
        },
        "fatigue_limitations": [
            "Coffin-Manson e Engelmaier são comparados, não combinados por média.",
            "O expoente fixo de Coffin-Manson é hipótese de sensibilidade.",
            "A aplicabilidade dos parâmetros de Engelmaier deve ser confirmada.",
            "A geometria e deformação da junta são aproximações não calibradas.",
            "O efeito mecânico do underfill não é representado por simples ajuste de 2ε'f.",
        ],
        "sspc_limitations": [
            "A energia retangular é envelope condicional, não energia medida.",
            "A energia real exige forma de onda e dispositivo selecionado.",
            "O modelo ideal de TVS não substitui SOA e curva de clamp.",
        ],
    }

    json_path = OUTPUT_DIR / "simulacao_4_6_resultados.json"
    with json_path.open("w", encoding="utf-8") as file:
        json.dump(payload, file, ensure_ascii=False, indent=2)

    print("Simulação 4.6 concluída: triagem in silico, não qualificada.")
    print(f"SSPC: I_lim assumida = {sspc['i_limit_assumed_a']:.2f} A")
    print(
        "PUE orbital adaptado por cenário:"
    )
    for row in pue_rows:
        print(
            f"  {row['orbit']} / {row['control']}: "
            f"{row['pue_orbital_adapted']:.4f}"
        )

    print(f"Arquivos gerados em: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
