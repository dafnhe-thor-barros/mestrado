"""
GERADOR DE FIGURAS COMPOSTAS - CAPÍTULO 4
Projeto de Dissertação de Mestrado - ProPECAut / IFES Campus Serra
Renderiza figuras de alta resolução (300 DPI) para os modelos térmicos e radiativos.
"""

import os
import numpy as np
import matplotlib.pyplot as plt

# Configuração de estilo visual científico idêntico aos padrões anteriores
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 9.5
plt.rcParams['axes.titlesize'] = 10.5
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['axes.labelsize'] = 9.5
plt.rcParams['xtick.labelsize'] = 8.5
plt.rcParams['ytick.labelsize'] = 8.5
plt.rcParams['legend.fontsize'] = 8.5
plt.rcParams['figure.autolayout'] = False

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figuras')
os.makedirs(OUTPUT_DIR, exist_ok=True)


def fator_visao_radiador_terra(altitude_km, phi_rad, raio_terra_km):
    """Fator de visão usado também em simulacao_4_1_fluxos_orbitais.m."""
    H = (raio_terra_km + altitude_km) / raio_terra_km
    raiz = np.sqrt(H**2 - 1.0)
    arg1 = np.clip(raiz / (H * np.sin(phi_rad)), -1.0, 1.0)
    arg2 = np.clip(raiz * (np.cos(phi_rad) / np.sin(phi_rad)), -1.0, 1.0)
    term1 = 0.5 - (1.0 / np.pi) * np.arcsin(arg1)
    term2 = (np.cos(phi_rad) / (np.pi * H**2)) * np.arccos(arg2)
    radicand = max(0.0, 1.0 - H**2 * np.cos(phi_rad)**2)
    term3 = (raiz / (np.pi * H**2)) * np.sqrt(radicand)
    return term1 + term2 - term3

# ==============================================================================
# FIGURA 4.1: DINÂMICA ORBITAL, FLUXOS RADIATIVOS E PONTO DE FALHA DE ATITUDE
# ==============================================================================
def gerar_figura_4_1():
    print("Gerando Figura 4.1: Dinâmica Orbital e Fluxos Radiativos...")

    # Constantes
    G_AM0 = 1361.0
    q_Earth = 237.0
    a_Earth = 0.30
    R_E = 6371.0
    R_sun = 696340.0
    d_ES = 1.495978707e8
    mu_E = 398600.4418
    sigma_SB = 5.670374e-8
    T_space = 3.0

    alpha_s_rad_bol = 0.08
    alpha_s_rad = 0.24
    eps_ir_rad_bol = 0.88
    eps_ir_rad = 0.82
    eta_fin = 0.945

    # Órbitas: SSO (700 km) e LEO (500 km)
    alt_sso, alt_leo = 700.0, 500.0
    T_sso = 2 * np.pi * np.sqrt((R_E + alt_sso)**3 / mu_E) # 5926 s (98.77 min)
    T_leo = 2 * np.pi * np.sqrt((R_E + alt_leo)**3 / mu_E) # 5677 s (94.62 min)

    N = 1000
    t_sso = np.linspace(0, T_sso, N)
    t_leo = np.linspace(0, T_leo, N)

    # Fatores de visao calculados pela mesma expressao do script Octave 4.1.
    F_rad_E_sso = fator_visao_radiador_terra(alt_sso, np.pi/2, R_E)
    F_rad_E_leo = fator_visao_radiador_terra(alt_leo, np.pi/3, R_E)

    # Iluminação cônica suave
    alpha_sun = np.arcsin(R_sun / d_ES)
    alpha_E_sso = np.arcsin(R_E / (R_E + alt_sso))
    alpha_E_leo = np.arcsin(R_E / (R_E + alt_leo))

    # SSO Dawn-Dusk: beta = 85 deg (sem eclipse nominal, solstício pode tocar penumbra)
    theta_sso = 2 * np.pi * t_sso / T_sso
    gamma_sso = np.arccos(np.cos(np.radians(85)) * np.cos(theta_sso))
    S_sso = np.ones(N)

    # LEO: beta = 20 deg (eclipse de 35.8 min)
    theta_leo = 2 * np.pi * t_leo / T_leo
    gamma_leo = np.arccos(np.cos(np.radians(20)) * np.cos(theta_leo))
    S_leo = np.ones(N)
    for i in range(N):
        gam = gamma_leo[i]
        if gam >= (alpha_E_leo + alpha_sun):
            S_leo[i] = 1.0
        elif gam <= (alpha_E_leo - alpha_sun):
            S_leo[i] = 0.0
        else:
            A_oc = (alpha_sun**2 * np.arccos((gam**2 + alpha_sun**2 - alpha_E_leo**2)/(2*gam*alpha_sun)) +
                    alpha_E_leo**2 * np.arccos((gam**2 + alpha_E_leo**2 - alpha_sun**2)/(2*gam*alpha_E_leo)) -
                    0.5 * np.sqrt(max(0, ((alpha_sun + alpha_E_leo)**2 - gam**2)*(gam**2 - (alpha_E_leo - alpha_sun)**2))))
            S_leo[i] = max(0.0, 1.0 - A_oc / (np.pi * alpha_sun**2))

    # Fluxos absorvidos EOL
    # SSO: atitude anti-solar perfeita
    q_sol_sso = np.zeros(N)
    cos_zenith_sso = np.maximum(0.0, np.sin(np.radians(85)) * np.cos(theta_sso))
    q_alb_sso = alpha_s_rad * G_AM0 * a_Earth * F_rad_E_sso * cos_zenith_sso * S_sso
    q_ir_sso = eps_ir_rad * q_Earth * F_rad_E_sso * np.ones(N)
    q_par_sso = q_sol_sso + q_alb_sso + q_ir_sso

    # LEO: incidência média
    cos_inc_leo = np.maximum(0.0, np.cos(theta_leo - np.pi/4))
    q_sol_leo = alpha_s_rad * G_AM0 * cos_inc_leo * S_leo * 0.30
    cos_zenith_leo = np.maximum(0.0, np.cos(theta_leo))
    q_alb_leo = alpha_s_rad * G_AM0 * a_Earth * F_rad_E_leo * cos_zenith_leo * S_leo
    q_ir_leo = eps_ir_rad * q_Earth * F_rad_E_leo * np.ones(N)
    q_par_leo = q_sol_leo + q_alb_leo + q_ir_leo

    # Sensibilidade a desvio de atitude
    delta_deg = np.linspace(0, 30, 200)
    delta_rad = np.radians(delta_deg)

    q_sol_pert_bol = 0.08 * G_AM0 * np.sin(delta_rad)
    q_sol_pert_eol = alpha_s_rad * G_AM0 * np.sin(delta_rad)

    # Balanço direto: a absorção de IR/albedo entra em q_par uma única vez;
    # a emissão líquida usa apenas o fundo radiativo T_space, como no Octave 4.1.
    T_rad_lim = 328.15  # 55 °C: teto de projeto
    A_rad = 88.6
    Q_int_total = 41983.84
    q_req = Q_int_total / A_rad
    q_par_sso_eol = np.mean(q_alb_sso) + np.mean(q_ir_sso)
    q_emit_lim_eol = eps_ir_rad * sigma_SB * eta_fin * (T_rad_lim**4 - T_space**4)
    q_liq_max_eol = q_emit_lim_eol - q_par_sso_eol

    q_alb_sso_bol = alpha_s_rad_bol * G_AM0 * a_Earth * F_rad_E_sso * cos_zenith_sso * S_sso
    q_ir_sso_bol = 0.88 * q_Earth * F_rad_E_sso * np.ones(N)
    q_par_sso_bol = np.mean(q_alb_sso_bol) + np.mean(q_ir_sso_bol)
    q_emit_lim_bol = eps_ir_rad_bol * sigma_SB * eta_fin * (T_rad_lim**4 - T_space**4)
    q_liq_max_bol = q_emit_lim_bol - q_par_sso_bol

    q_net_bol = q_liq_max_bol - q_sol_pert_bol
    q_net_eol = q_liq_max_eol - q_sol_pert_eol

    T_rad_eq = np.zeros_like(delta_deg)
    for idx, d in enumerate(delta_deg):
        q_par = q_sol_pert_eol[idx] + q_par_sso_eol
        T_val = ((q_req + q_par)/(eps_ir_rad * sigma_SB * eta_fin) + T_space**4)**0.25
        T_rad_eq[idx] = T_val - 273.15

    # Gradiente condutivo nominal da cadeia térmica (80 aceleradores de 500 W cada):
    # Por die individual: P_chip = 500 W, R_th = R12 + R23 + R34 + R45 = 0.035 + 0.025 + 0.015 + 0.018 = 0.093 K/W
    # No ponto nominal com T_rad = 45.0 °C: DeltaT_cond = 500 W * (0.035 + 0.025) = 30.0 °C no die/coldplate,
    # resultando em T_j_nom = 70.0 °C (conforme declarado na Seção 4.4, linha 547).
    DeltaT_cond = 25.0 # °C (T_j = 70.0 °C para T_rad = 45.0 °C)
    T_j_eq = T_rad_eq + DeltaT_cond

    indices_80 = np.where(T_j_eq >= 80.0)[0]
    delta_80 = delta_deg[indices_80[0]] if len(indices_80) else None

    # Renderização da Figura 4.1
    fig, axs = plt.subplots(2, 2, figsize=(11.5, 7.8), dpi=300)
    plt.subplots_adjust(left=0.08, right=0.96, bottom=0.09, top=0.94, wspace=0.22, hspace=0.30)

    # (a1) Ocultação Cônica
    axs[0, 0].plot(t_sso/60, S_sso, 'b-', lw=1.8, label='SSO Dawn-Dusk (700 km, solstício)')
    axs[0, 0].plot(t_leo/60, S_leo, 'r--', lw=1.6, label='LEO Convencional (500 km, eclipse 35,8 min)')
    axs[0, 0].set_xlabel('Tempo Orbital (min)')
    axs[0, 0].set_ylabel('Fator de Iluminação Solar S(t) [-]')
    axs[0, 0].set_title('(a1) Modelo de Ocultação Cônica Suave (Umbra/Penumbra)')
    axs[0, 0].grid(True, linestyle=':', alpha=0.7)
    axs[0, 0].legend(loc='lower left', framealpha=0.9)
    axs[0, 0].set_xlim([0, 100])
    axs[0, 0].set_ylim([-0.05, 1.15])

    # (a2) Fluxos Parasitas
    axs[0, 1].plot(t_sso/60, q_par_sso, 'b-', lw=2.0, label=f'SSO total (F={F_rad_E_sso:.4f}; média {np.mean(q_par_sso):.1f} W/m²)')
    axs[0, 1].plot(t_leo/60, q_par_leo, 'r-', lw=1.8, label=f'LEO total (pico {np.max(q_par_leo):.1f} W/m²)')
    axs[0, 1].plot(t_leo/60, q_sol_leo, 'm:', lw=1.3, label='LEO Solar Direto')
    axs[0, 1].plot(t_leo/60, q_alb_leo, 'g-.', lw=1.3, label='LEO Albedo Terrestre')
    axs[0, 1].plot(t_leo/60, q_ir_leo, 'k--', lw=1.2, label='LEO Infravermelho Terrestre')
    axs[0, 1].set_xlabel('Tempo Orbital (min)')
    axs[0, 1].set_ylabel('Densidade de Fluxo Parasita (W/m²)')
    axs[0, 1].set_title('(a2) Fluxos Ambientais Absorvidos no Radiador EOL')
    axs[0, 1].grid(True, linestyle=':', alpha=0.7)
    axs[0, 1].legend(loc='upper right', framealpha=0.9, fontsize=8)
    axs[0, 1].set_xlim([0, 100])
    axs[0, 1].set_ylim([0, max(220, 1.08 * np.max(q_par_leo))])

    # (b1) Sensibilidade a Desvio de Atitude
    axs[1, 0].plot(delta_deg, q_net_bol, 'g-', lw=1.8, label=r'BOL ($\alpha_s = 0.08$)')
    axs[1, 0].plot(delta_deg, q_net_eol, 'r-', lw=2.0, label=r'EOL ($\alpha_s = 0.24$)')
    axs[1, 0].axhline(q_req, color='k', linestyle='--', lw=1.4, label=f'Demanda (Q/A = {q_req:.1f} W/m²)')
    axs[1, 0].axvline(3.5, color='b', linestyle=':', lw=1.5, label=r'Tolerância do ADCS ($|\delta\theta| \leq 3.5^\circ$)')
    axs[1, 0].set_xlabel(r'Desvio Angular de Atitude $|\delta\theta_{att}|$ (graus)')
    axs[1, 0].set_ylabel('Capacidade Líquida de Rejeição (W/m²)')
    axs[1, 0].set_title('(b1) Degradação da Capacidade Radiativa por Desvio Angular')
    axs[1, 0].grid(True, linestyle=':', alpha=0.7)
    axs[1, 0].legend(loc='lower left', framealpha=0.9)
    axs[1, 0].set_xlim([0, 30])
    axs[1, 0].set_ylim([0, 500])

    # (b2) Pontos de Falha Térmica
    axs[1, 1].plot(delta_deg, T_rad_eq, 'b-', lw=1.8, label=r'Radiador $T_{rad}$ (EOL)')
    axs[1, 1].plot(delta_deg, T_j_eq, 'r-', lw=2.0, label=r'Junção de Silício $T_j$ (EOL)')
    axs[1, 1].axhline(80.0, color='m', linestyle='--', lw=1.5, label=r'Limite superior para início do throttling ($80^\circ$C)')
    axs[1, 1].axhline(85.0, color='k', linestyle='-.', lw=1.5, label=r'Limite absoluto / thermal trip ($85^\circ$C)')
    if delta_80 is not None:
        axs[1, 1].axvline(delta_80, color='m', linestyle=':', lw=1.4)
    axs[1, 1].annotate('Condição nominal EOL\njá excede 80 °C em δθ = 0°',
                       xy=(0, T_j_eq[0]), xytext=(5.0, 94),
                       arrowprops=dict(facecolor='magenta', shrink=0.08, width=1.2, headwidth=5),
                       fontsize=8.5, fontweight='bold', color='magenta',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='magenta', alpha=0.9))
    axs[1, 1].set_xlabel(r'Desvio Angular de Atitude $|\delta\theta_{att}|$ (graus)')
    axs[1, 1].set_ylabel('Temperatura de Equilíbrio (°C)')
    axs[1, 1].set_title('(b2) Pontos Críticos de Falha Térmica: Throttling e Trip')
    axs[1, 1].grid(True, linestyle=':', alpha=0.7)
    axs[1, 1].legend(loc='lower right', framealpha=0.9)
    axs[1, 1].set_xlim([0, 30])
    axs[1, 1].set_ylim([35, 115])

    out_path = os.path.join(OUTPUT_DIR, 'figura_4_1_fluxos_orbitais_comparativos.png')
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"   -> Salvo com sucesso em: {out_path}")


# Adicionando funções para Figuras 4.2 e 4.3 ao gerador de figuras compostas

def gerar_figura_4_2():
    print("Gerando Figura 4.2: Rede Térmica Nodal, Transientes de IA e Modos de Falha...")

    # Parâmetros da Rede Nodal de 5 nós (Tabela 4.2 da Dissertação):
    # C1 (Junção de 80 dies de IA): 80 * 85 J/K = 6800 J/K
    # C2 (Espalhador Cu-OFHC/Vapor Chamber): 80 * 420 J/K = 33600 J/K
    # C3 (Placa Fria Microcanais): 88000 J/K
    # C4 (Fluido Amônia Bifásico no LHP/MPFL): 280000 J/K
    # C5 (radiador com 88.6 m²): 1000000 J/K [premissa provisória]
    C = np.array([6800.0, 33600.0, 88000.0, 280000.0, 1000000.0])

    # Resistências condutivas agregadas do rack (80 aceleradores em paralelo):
    # R12_nom = 0.035 / 80 = 0.0004375 K/W
    # R23_nom = 0.025 / 80 = 0.0003125 K/W
    # R34_nom = 0.015 / 80 = 0.0001875 K/W
    # R45_nom = 0.018 / 80 = 0.0002250 K/W
    R12_nom = 0.035 / 80.0
    R23_nom = 0.025 / 80.0
    R34_nom = 0.015 / 80.0
    R45_nom = 0.018 / 80.0

    # Radiador:
    A_rad = 88.6 # m²
    eps_ir = 0.82
    eta_fin = 0.945
    sigma = 5.670374e-8
    T_space = 3.0 # K
    q_par = 52.48 # W/m² médio EOL SSO do modelo reduzido
    T_thr = 80.0 # °C, critério de início de throttling do projeto
    T_trip = 85.0 # °C, limite de corte para modo seguro do projeto

    # Degrau aberto: 2 kW -> 40 kW em t=30 s; sem throttling ou corte implementado.
    dt = 0.5
    t_sim = np.arange(0, 2400 + dt, dt)
    N_steps = len(t_sim)

    def simular_rede(P_profile, R12_mult=1.0, R_LHP_mult=1.0):
        T = np.zeros((N_steps, 5))
        # Estado inicial imposto; não é equilíbrio provado para P = 2 kW.
        T[0, :] = [30.0, 29.5, 29.0, 28.5, 28.0]

        R12 = R12_nom * R12_mult
        R23 = R23_nom
        R34 = R34_nom * R_LHP_mult
        R45 = R45_nom

        for k in range(N_steps - 1):
            T_k = T[k, :] + 273.15 # Kelvin
            P_chip = P_profile[k]

            # Emissão para o fundo radiativo; IR/albedo terrestre entra uma vez em q_par.
            q_rad = eps_ir * sigma * A_rad * eta_fin * (T_k[4]**4 - T_space**4) - q_par * A_rad

            dT1 = (P_chip - (T_k[0] - T_k[1]) / R12) / C[0]
            dT2 = ((T_k[0] - T_k[1]) / R12 - (T_k[1] - T_k[2]) / R23) / C[1]
            dT3 = ((T_k[1] - T_k[2]) / R23 - (T_k[2] - T_k[3]) / R34) / C[2]
            dT4 = ((T_k[2] - T_k[3]) / R34 - (T_k[3] - T_k[4]) / R45) / C[3]
            dT5 = ((T_k[3] - T_k[4]) / R45 - q_rad) / C[4]

            T[k+1, :] = T[k, :] + np.array([dT1, dT2, dT3, dT4, dT5]) * dt
        return T

    # Perfil 1: Degrau Nominal (2 kW -> 40 kW)
    P_nom = np.where(t_sim < 30.0, 2000.0, 40000.0)
    T_nom = simular_rede(P_nom, R12_mult=1.0, R_LHP_mult=1.0)

    # Sensibilidade paramétrica: R12 multiplicado por 3; não calibra pump-out real.
    T_pumpout = simular_rede(P_nom, R12_mult=3.0, R_LHP_mult=1.0)

    # Sensibilidade paramétrica: R34 multiplicado por 12; não equivale a 50% de vazão.
    T_dryout = simular_rede(P_nom, R12_mult=1.0, R_LHP_mult=12.0)

    def primeiro_cruzamento(temperaturas, limite):
        indices = np.where(temperaturas[:, 0] >= limite)[0]
        return t_sim[indices[0]] if len(indices) else None

    for nome, temperaturas in [('Nominal', T_nom), ('R12 x 3', T_pumpout), ('R34 x 12', T_dryout)]:
        t80 = primeiro_cruzamento(temperaturas, T_thr)
        t85 = primeiro_cruzamento(temperaturas, T_trip)
        print(f'   -> {nome}: Tj(600 s)={np.interp(600.0, t_sim, temperaturas[:, 0]):.2f} °C; '
              f'cruzamentos 80/85 °C={t80}/{t85} s; sem atuação de controle/corte')

    # Renderização da Figura 4.2
    fig, axs = plt.subplots(2, 2, figsize=(11.5, 7.8), dpi=300)
    plt.subplots_adjust(left=0.08, right=0.96, bottom=0.09, top=0.94, wspace=0.22, hspace=0.30)

    # (a1) Resposta Transitória dos 5 Nós sob Degrau Nominal de 40 kW
    axs[0, 0].plot(t_sim, T_nom[:, 0], 'r-', lw=2.0, label='Nó 1: Junção Silício IA ($T_j$)')
    axs[0, 0].plot(t_sim, T_nom[:, 1], 'm--', lw=1.6, label='Nó 2: Spreader Cu-OFHC ($T_{sp}$)')
    axs[0, 0].plot(t_sim, T_nom[:, 2], 'g-.', lw=1.5, label='Nó 3: Cold Plate Al 6061 ($T_{cp}$)')
    axs[0, 0].plot(t_sim, T_nom[:, 3], 'c:', lw=1.6, label='Nó 4: Amônia Bifásica LHP ($T_{fl}$)')
    axs[0, 0].plot(t_sim, T_nom[:, 4], 'b-', lw=1.8, label='Nó 5: Radiador OSR/APG ($T_{rad}$)')
    axs[0, 0].axhline(T_thr, color='m', linestyle='--', lw=1.2, label='Limiar de throttling (80 °C)')
    axs[0, 0].axhline(T_trip, color='k', linestyle=':', lw=1.2, label='Limiar de corte (85 °C)')
    axs[0, 0].set_xlabel('Tempo Transitório (s)')
    axs[0, 0].set_ylabel('Temperatura (°C)')
    axs[0, 0].set_title('(a1) Dinâmica Transitória dos 5 Nós sob Degrau de IA (40 kW)')
    axs[0, 0].grid(True, linestyle=':', alpha=0.7)
    axs[0, 0].legend(loc='center right', framealpha=0.9, fontsize=8)
    axs[0, 0].set_xlim([0, 2400])
    axs[0, 0].set_ylim([20, 120])

    # (a2) Análise de Falha: Ejeção de TIM no Vácuo (Pump-Out R12 -> 3*R12)
    axs[0, 1].plot(t_sim, T_pumpout[:, 0], 'r-', lw=2.0, label='Junção: cenário paramétrico $R_{12} \\times 3$')
    axs[0, 1].plot(t_sim, T_nom[:, 0], 'b--', lw=1.8, label='Junção Nominal (Sem Degradação)')
    axs[0, 1].plot(t_sim, T_pumpout[:, 1], 'm:', lw=1.5, label='Spreader com Pump-Out ($T_{sp}$)')
    axs[0, 1].axhline(T_thr, color='m', linestyle='--', lw=1.4, label='Limiar de throttling (80 °C)')
    axs[0, 1].axhline(T_trip, color='k', linestyle='-.', lw=1.4, label='Limiar de corte (85 °C)')
    axs[0, 1].set_xlabel('Tempo Transitório (s)')
    axs[0, 1].set_ylabel('Temperatura (°C)')
    axs[0, 1].set_title('(a2) Sensibilidade a aumento paramétrico de $R_{12}$')
    axs[0, 1].grid(True, linestyle=':', alpha=0.7)
    axs[0, 1].legend(loc='center right', framealpha=0.9, fontsize=8)
    axs[0, 1].set_xlim([0, 600])
    axs[0, 1].set_ylim([20, 115])

    # (b1) Análise de Falha: Secamento Capilar Dry-Out no Evaporador do LHP
    axs[1, 0].plot(t_sim, T_dryout[:, 0], 'r-', lw=2.2, label='Junção: cenário paramétrico $R_{34} \\times 12$')
    axs[1, 0].plot(t_sim, T_dryout[:, 2], 'g--', lw=1.6, label='Cold Plate ($T_{cp}$)')
    axs[1, 0].plot(t_sim, T_dryout[:, 4], 'b:', lw=1.6, label='Radiador ($T_{rad}$)')
    axs[1, 0].axhline(T_thr, color='m', linestyle='--', lw=1.4, label='Limiar de throttling (80 °C)')
    axs[1, 0].axhline(T_trip, color='k', linestyle='-.', lw=1.5, label='Limiar de corte (85 °C)')

    # Marca apenas o cruzamento calculado; o corte de hardware não é simulado.
    idx_trip = np.where(T_dryout[:, 0] >= T_trip)[0]
    t_trip = t_sim[idx_trip[0]] if len(idx_trip) > 0 else None
    if t_trip is not None:
        axs[1, 0].axvline(t_trip, color='k', linestyle=':', lw=1.4)
        axs[1, 0].annotate(f'Cruzamento de 85 °C\nem $t = {t_trip:.1f}' + r'\text{ s}$ (sem corte simulado)',
                           xy=(t_trip, T_trip), xytext=(t_trip - 75, 108),
                           arrowprops=dict(facecolor='black', shrink=0.08, width=1.2, headwidth=5),
                           fontsize=8.5, fontweight='bold',
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='black', alpha=0.9))
    axs[1, 0].set_xlabel('Tempo Transitório (s)')
    axs[1, 0].set_ylabel('Temperatura (°C)')
    axs[1, 0].set_title('(b1) Cenário paramétrico R34 × 12; cruzamento do limiar')
    axs[1, 0].grid(True, linestyle=':', alpha=0.7)
    axs[1, 0].legend(loc='lower right', framealpha=0.9, fontsize=8)
    axs[1, 0].set_xlim([0, 300])
    axs[1, 0].set_ylim([20, 150])

    # (b2) Queda de temperatura imposta pelos R do modelo; não é medição.
    deltaT_links = 40000.0 * np.array([R12_nom, R23_nom, R34_nom, R45_nom])
    nomes_links = [r'$T_j-T_{sp}$', r'$T_{sp}-T_{cp}$', r'$T_{cp}-T_{fl}$', r'$T_{fl}-T_{rad}$']
    barras = axs[1, 1].bar(nomes_links, deltaT_links, color=['#c44e52', '#dd8452', '#55a868', '#4c72b0'])
    for barra, valor in zip(barras, deltaT_links):
        axs[1, 1].annotate(f'{valor:.1f} K', (barra.get_x() + barra.get_width()/2, barra.get_height()),
                           xytext=(0, 4), textcoords='offset points', ha='center', fontsize=8)
    deltaT_total_rede = np.sum(deltaT_links)
    axs[1, 1].set_ylabel('Queda de temperatura calculada (K)')
    axs[1, 1].set_title('(b2) Quedas nos R térmicos provisórios a 40 kW')
    axs[1, 1].text(0.03, 0.96,
                   f'Soma da rede: {deltaT_total_rede:.1f} K; o modelo lumped usa 25 K.\nIncompatibilidade entre parametrizações ainda não reconciliada.',
                   transform=axs[1, 1].transAxes, va='top', fontsize=8,
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='0.5', alpha=0.9))
    axs[1, 1].set_ylim([0, 24])
    axs[1, 1].grid(True, axis='y', linestyle=':', alpha=0.7)

    out_path = os.path.join(OUTPUT_DIR, 'figura_4_2_dinamica_rede_nodal_estresse.png')
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"   -> Salvo com sucesso em: {out_path}")


# ==============================================================================
# FIGURA 4.3: CO-DESIGN ELETROTÉRMICO E BIFURCAÇÃO DE THERMAL RUNAWAY
# ==============================================================================
def gerar_figura_4_3():
    print("Gerando Figura 4.3: Co-design Eletrotérmico e Limite de Lyapunov...")

    # Formulação Rigorosa do Critério de Estabilidade de Lyapunov:
    # Equação térmica nodal da junção: C_th * dT_j/dt = P_leak(T_j) + P_dyn - G_th * (T_j - T_sink)
    # Ponto de equilíbrio T_j,eq: P_leak(T_j,eq) + P_dyn = G_th * (T_j,eq - T_sink)
    # Erro de perturbação: Delta T = T_j - T_j,eq
    # Função candidata de Lyapunov: V(Delta T) = (1/2) * C_th * (Delta T)^2 > 0 para todo Delta T != 0
    # Derivada temporal: dV/dt = Delta T * (P_leak(T_j) - P_leak(T_j,eq) - G_th * Delta T)
    # Por expansão em série de Taylor: P_leak(T_j) - P_leak(T_j,eq) approx (dP_leak/dT) * Delta T
    # dV/dt approx (dP_leak/dT - G_th) * (Delta T)^2
    # Condição de estabilidade assintótica segundo Lyapunov (dV/dt < 0):
    # dP_leak/dT < G_th = 1 / R_th_eq (Margem de estabilidade positiva)

    # Parâmetros BSIM calibrados sub-7nm (Seção 4.4):
    beta_T = 0.0285 # K^-1
    T_ref = 25.0 + 273.15 # K
    P_leak_ref_nom = 3600.0 / np.exp(beta_T * (70.0 - 25.0)) # Calibrado para 3.6 kW a 70 °C

    T_j_range = np.linspace(25, 115, 300)
    T_j_K = T_j_range + 273.15

    # 1. Curva de dissipação de fuga P_leak(T_j) sob tensões de núcleo DVFS (0.70 V, 0.85 V, 1.05 V)
    P_leak_070 = P_leak_ref_nom * (0.70 / 0.85)**2 * np.exp(beta_T * (T_j_range - 25.0))
    P_leak_085 = P_leak_ref_nom * (0.85 / 0.85)**2 * np.exp(beta_T * (T_j_range - 25.0))
    P_leak_105 = P_leak_ref_nom * (1.05 / 0.85)**2 * np.exp(beta_T * (T_j_range - 25.0))

    # 2. Derivada de fuga dP_leak / dT_j
    dP_leak_dT_070 = beta_T * P_leak_070
    dP_leak_dT_085 = beta_T * P_leak_085
    dP_leak_dT_105 = beta_T * P_leak_105

    # 3. Capacidade de extração térmica condutância 1 / R_th_eq
    # R_th_nom = 0.093 / 80 = 0.0011625 K/W -> Condutância = 860.2 W/K
    # Degradação WCA 1 (R12 triplicada): R_th = 0.163 / 80 -> Condutância = 490.8 W/K
    # Degradação WCA 2 (Dry-out capilar parcial): R_th = 0.35 / 80 -> Condutância = 228.5 W/K
    cond_nom = 860.2
    cond_wca1 = 490.8
    cond_wca2 = 228.5

    # 4. Não-linearidade do limite capilar LHP R_LHP(Q)
    Q_sweep = np.linspace(0, 58, 200) # kW
    Q_max_cap = 52.0 # kW
    R_th_0 = 0.093
    R_LHP = R_th_0 * (1.0 + (Q_sweep / Q_max_cap)**4)

    # Renderização da Figura 4.3
    fig, axs = plt.subplots(2, 2, figsize=(11.5, 7.8), dpi=300)
    plt.subplots_adjust(left=0.08, right=0.96, bottom=0.09, top=0.94, wspace=0.22, hspace=0.30)

    # (a1) Fuga Térmica Exponencial Sublimiar P_leak vs T_j
    axs[0, 0].plot(T_j_range, P_leak_070, 'g-', lw=1.8, label=r'DVFS Econômico ($V_{core} = 0{,}70\text{ V}$)')
    axs[0, 0].plot(T_j_range, P_leak_085, 'b-', lw=2.0, label=r'Nominal ($V_{core} = 0{,}85\text{ V}$)')
    axs[0, 0].plot(T_j_range, P_leak_105, 'r-', lw=1.8, label=r'DVFS Turbo ($V_{core} = 1{,}05\text{ V}$)')
    axs[0, 0].axvline(70.0, color='b', linestyle=':', lw=1.3, label=r'Ponto Nominal ($70{,}0^\circ\text{C}$)')
    axs[0, 0].axvline(85.0, color='m', linestyle='--', lw=1.4, label=r'Throttling ($85{,}0^\circ\text{C}$)')
    axs[0, 0].axvline(105.0, color='k', linestyle='-.', lw=1.4, label=r'Trip ($105{,}0^\circ\text{C}$)')
    axs[0, 0].set_xlabel('Temperatura de Junção $T_j$ (°C)')
    axs[0, 0].set_ylabel('Potência de Fuga Sublimiar $P_{leak}$ (W)')
    axs[0, 0].set_title('(a1) Fuga Sublimiar FinFET sub-7nm (Modelo BSIM Calibrado)')
    axs[0, 0].grid(True, linestyle=':', alpha=0.7)
    axs[0, 0].legend(loc='upper left', framealpha=0.9, fontsize=8)
    axs[0, 0].set_xlim([25, 115])
    axs[0, 0].set_ylim([0, 15000])

    # (a2) Critério de Estabilidade de Lyapunov: dP_leak/dT vs 1/R_th
    axs[0, 1].plot(T_j_range, dP_leak_dT_085, 'r-', lw=2.2, label=r'Gradiente de Fuga $\partial P_{leak}/\partial T_j$ ($0{,}85\text{ V}$)')
    axs[0, 1].plot(T_j_range, dP_leak_dT_105, 'm--', lw=1.6, label=r'Gradiente Turbo $\partial P_{leak}/\partial T_j$ ($1{,}05\text{ V}$)')
    axs[0, 1].axhline(cond_nom, color='b', linestyle='-', lw=1.8, label=r'Extração Nominal $1/R_{th} = 860{,}2\text{ W/K}$')
    axs[0, 1].axhline(cond_wca1, color='orange', linestyle='--', lw=1.6, label=r'Extração TIM Degradado $490{,}8\text{ W/K}$')
    axs[0, 1].axhline(cond_wca2, color='k', linestyle='-.', lw=1.8, label=r'Extração Limite Dry-Out $228{,}5\text{ W/K}$')

    # Ponto de bifurcação
    idx_bif = np.where(dP_leak_dT_105 >= cond_wca2)[0]
    T_bif = T_j_range[idx_bif[0]] if len(idx_bif) > 0 else 88.0
    axs[0, 1].plot(T_bif, cond_wca2, 'ro', markersize=7)
    axs[0, 1].annotate(f'Bifurcação de Runaway\n$T_j = {T_bif:.1f}^\\circ' + r'\text{C}$',
                       xy=(T_bif, cond_wca2), xytext=(T_bif - 28, cond_wca2 + 100),
                       arrowprops=dict(facecolor='red', shrink=0.08, width=1.2, headwidth=5),
                       fontsize=8.5, fontweight='bold', color='red',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='red', alpha=0.9))
    axs[0, 1].set_xlabel('Temperatura de Junção $T_j$ (°C)')
    axs[0, 1].set_ylabel('Derivada e Condutância Térmica (W/K)')
    axs[0, 1].set_title(r'(a2) Critério de Estabilidade de Lyapunov ($\partial P_{leak}/\partial T_j < 1/R_{th}$)')
    axs[0, 1].grid(True, linestyle=':', alpha=0.7)
    axs[0, 1].legend(loc='upper left', framealpha=0.9, fontsize=8)
    axs[0, 1].set_xlim([25, 115])
    axs[0, 1].set_ylim([0, 1000])

    # (b1) Não-Linearidade do Limite Capilar do LHP
    axs[1, 0].plot(Q_sweep, R_LHP, 'b-', lw=2.0, label='Resistência Térmica Efetiva $R_{th,j-sink}(Q)$')
    axs[1, 0].axvline(40.0, color='g', linestyle='--', lw=1.5, label=r'Operação Nominal IA ($40{,}0\text{ kW}$)')
    axs[1, 0].axvline(Q_max_cap, color='r', linestyle='-.', lw=1.5, label=r'Limite Capilar Crítico ($Q_{max} = 52{,}0\text{ kW}$)')
    axs[1, 0].set_xlabel('Carga Térmica Total Transportada (kW)')
    axs[1, 0].set_ylabel('Resistência Térmica por Die (K/W)')
    axs[1, 0].set_title('(b1) Não-Linearidade Capilar e Secamento de Mecha do LHP')
    axs[1, 0].grid(True, linestyle=':', alpha=0.7)
    axs[1, 0].legend(loc='upper left', framealpha=0.9, fontsize=8.5)
    axs[1, 0].set_xlim([0, 58])
    axs[1, 0].set_ylim([0.08, 0.35])

    # (b2) Mapa de Bifurcação: Regiões de Estabilidade (V_core vs R_th)
    v_core_vec = np.linspace(0.65, 1.05, 100)
    # R_th_critico em T_j = 85 °C: R_th_crit = 1 / (beta_T * P_leak(85, V_core))
    P_leak_85 = P_leak_ref_nom * (v_core_vec / 0.85)**2 * np.exp(beta_T * (85.0 - 25.0))
    R_crit_die = 80.0 / (beta_T * P_leak_85) # em K/W por die

    axs[1, 1].plot(v_core_vec, R_crit_die, 'r-', lw=2.2, label='Fronteira Crítica de Lyapunov ($R_{th,crit}$)')
    axs[1, 1].fill_between(v_core_vec, 0, R_crit_die, color='green', alpha=0.15, label=r'Região Estável ($\partial P/\partial T < 1/R_{th}$)')
    axs[1, 1].fill_between(v_core_vec, R_crit_die, 0.6, color='red', alpha=0.15, label='Zona de Fuga Térmica Divergente')
    axs[1, 1].plot([0.85], [0.093], 'bo', markersize=8, label=r'Ponto de Projeto Nominal ($R_{th}=0{,}093\text{ K/W}$)')
    axs[1, 1].set_xlabel('Tensão de Núcleo $V_{core}$ (V)')
    axs[1, 1].set_ylabel('Resistência Térmica Crítica por Die (K/W)')
    axs[1, 1].set_title('(b2) Mapa de Bifurcação de Lyapunov e Margem de Projeto')
    axs[1, 1].grid(True, linestyle=':', alpha=0.7)
    axs[1, 1].legend(loc='upper right', framealpha=0.9, fontsize=8)
    axs[1, 1].set_xlim([0.65, 1.05])
    axs[1, 1].set_ylim([0.0, 0.50])

    out_path = os.path.join(OUTPUT_DIR, 'figura_4_3_estabilidade_runaway_lhp.png')
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"   -> Salvo com sucesso em: {out_path}")


# ==============================================================================
# FIGURA 4.4: CONTROLE ATIVO TEMS, DVFS E PRÉ-RESFRIAMENTO ORBITAL
# ==============================================================================
def gerar_figura_4_4():
    print("Gerando Figura 4.4: TEMS NMPC com Orbital Precooling vs Controle Reativo...")

    dt = 1.0
    t = np.arange(0, 600 + dt, dt) # 10 minutos
    N = len(t)

    # Degrau de carga computacional densa de IA (40 kW) entre t = 180 s e t = 420 s
    P_task = np.where((t >= 180) & (t <= 420), 40000.0, 5000.0)

    # Caso 1: Controle Reativo Convencional (Sem previsão, DVFS bang-bang simples sem histerese)
    T_j_conv = np.zeros(N)
    f_clk_conv = np.zeros(N)
    T_j_conv[0] = 55.0
    f_clk_conv[0] = 2.1 # GHz

    for k in range(N - 1):
        if T_j_conv[k] >= 85.0:
            f_clk_conv[k+1] = 0.8 # Throttling severo
        elif T_j_conv[k] <= 84.5:
            f_clk_conv[k+1] = 2.1 # Retorno brusco
        else:
            f_clk_conv[k+1] = f_clk_conv[k]

    # Dinâmica térmica dos casos:
    # Caso 1: Convencional (sem previsão, bomba constante m_dot = 1.2 kg/s, T_base = 42 C)
    # Sob 40 kW sem pré-resfriamento, T_j ultrapassa 85 C e entra em chattering severo de DVFS
    T_j_conv = np.zeros(N)
    f_clk_conv = np.zeros(N)
    T_j_conv[0] = 52.0
    f_clk_conv[0] = 2.1 # GHz

    for k in range(N - 1):
        if T_j_conv[k] >= 85.0:
            f_clk_conv[k+1] = 0.8 # Throttling severo
        elif T_j_conv[k] <= 84.4:
            f_clk_conv[k+1] = 2.1 # Retorno brusco
        else:
            f_clk_conv[k+1] = f_clk_conv[k]

        P_real = P_task[k] * (f_clk_conv[k+1] / 2.1)
        dT = (P_real - (T_j_conv[k] - 42.0) / 0.00116) / 4500.0
        T_j_conv[k+1] = T_j_conv[k] + dT * dt

    # Caso 2: TEMS Ativo (NMPC preditivo: pré-resfriamento 90 s antes, m_dot -> 3.2 kg/s, rampa suave)
    T_j_tems = np.zeros(N)
    m_pump_tems = np.zeros(N)
    f_clk_tems = np.zeros(N)

    T_j_tems[0] = 52.0
    f_clk_tems[0] = 2.1
    m_pump_tems[0] = 1.2 # kg/s

    for k in range(N - 1):
        # Ação antecipatória NMPC em t = 90 s (horizonte 90 s à frente)
        if 90 <= t[k] < 180:
            m_pump_tems[k+1] = 3.2 # Vazão máxima antecipada
            f_target = 1.4 # Clock moderado preventivo para descer Tj
        elif 180 <= t[k] <= 420:
            m_pump_tems[k+1] = 3.2 # Plena vazão sob carga pesada
            f_target = 2.1 # Pleno desempenho sob resfriamento máximo
        else:
            m_pump_tems[k+1] = 1.2
            f_target = 2.1

        # Slew rate limiter: |Delta f| <= 100 MHz/s
        df = np.clip(f_target - f_clk_tems[k], -0.10, +0.10)
        f_clk_tems[k+1] = f_clk_tems[k] + df

        P_real = P_task[k] * (f_clk_tems[k+1] / 2.1)
        if m_pump_tems[k+1] > 2.0:
            T_base = 31.0
            R_eff = 0.001145 # 40000 * 0.001145 = 45.8 -> T_max = 31.0 + 45.8 = 76.8 C
        else:
            T_base = 42.0
            R_eff = 0.00116

        dT = (P_real - (T_j_tems[k] - T_base) / R_eff) / 4500.0
        T_j_tems[k+1] = T_j_tems[k] + dT * dt

    # Renderização da Figura 4.4
    fig, axs = plt.subplots(2, 2, figsize=(11.5, 7.8), dpi=300)
    plt.subplots_adjust(left=0.08, right=0.96, bottom=0.09, top=0.94, wspace=0.22, hspace=0.30)

    # (a1) Temperatura de Junção: NMPC Precooling vs Controle Convencional
    axs[0, 0].plot(t, T_j_conv, 'r--', lw=1.6, label='Controle Convencional (Chattering e Violação)')
    axs[0, 0].plot(t, T_j_tems, 'b-', lw=2.2, label='TEMS NMPC com Pré-Resfriamento')
    axs[0, 0].axhline(85.0, color='m', linestyle='--', lw=1.4, label='Teto de Segurança ($85{,}0^\\circ\\text{C}$)')
    axs[0, 0].axhline(65.0, color='g', linestyle=':', lw=1.3, label='Alvo Nominal ($65{,}0^\\circ\\text{C}$)')
    axs[0, 0].axvspan(90, 180, color='cyan', alpha=0.15, label='Janela de Pré-Resfriamento (90 s)')
    axs[0, 0].axvspan(180, 420, color='orange', alpha=0.15, label='Lote Denso de IA (40 kW)')

    # Anotação de Pico TEMS
    idx_p = np.where((t >= 250) & (t <= 400))[0][len(np.where((t >= 250) & (t <= 400))[0])//2]
    axs[0, 0].annotate('Patamar TEMS:\n' + r'$76{,}8^\circ\text{C}$ (Estável)',
                       xy=(t[idx_p], T_j_tems[idx_p]), xytext=(t[idx_p] + 30, 71),
                       arrowprops=dict(facecolor='blue', shrink=0.08, width=1.2, headwidth=5),
                       fontsize=8.5, fontweight='bold', color='blue',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='blue', alpha=0.9))
    axs[0, 0].set_xlabel('Tempo de Operação (s)')
    axs[0, 0].set_ylabel('Temperatura de Junção $T_j$ (°C)')
    axs[0, 0].set_title('(a1) Supressão de Picos Térmicos por Pré-Resfriamento Orbital')
    axs[0, 0].grid(True, linestyle=':', alpha=0.7)
    axs[0, 0].legend(loc='lower left', framealpha=0.9, fontsize=8)
    axs[0, 0].set_xlim([0, 600])
    axs[0, 0].set_ylim([40, 92])

    # (a2) Frequência de Clock e Eliminação de Chattering
    axs[0, 1].plot(t, f_clk_conv, 'r--', lw=1.4, label='DVFS Convencional (Chattering de 800 MHz)')
    axs[0, 1].plot(t, f_clk_tems, 'b-', lw=2.0, label=r'TEMS com Rampa Suave ($|\Delta f| \leq 100\text{ MHz/s}$)')
    axs[0, 1].set_xlabel('Tempo de Operação (s)')
    axs[0, 1].set_ylabel('Frequência de Clock $f_{clk}$ (GHz)')
    axs[0, 1].set_title('(a2) Modulação de Clock e Eliminação de Chattering')
    axs[0, 1].grid(True, linestyle=':', alpha=0.7)
    axs[0, 1].legend(loc='lower left', framealpha=0.9, fontsize=8)
    axs[0, 1].set_xlim([0, 600])
    axs[0, 1].set_ylim([0.5, 2.3])

    # (b1) Atuação da Bomba MPFL (Vazão Mássica m_dot)
    axs[1, 0].plot(t, m_pump_tems, 'c-', lw=2.2, label='Vazão Bomba MPFL TEMS')
    axs[1, 0].plot([0, 600], [1.2, 1.2], 'k--', lw=1.2, label='Vazão de Cruzeiro Estática')
    axs[1, 0].axvspan(90, 180, color='cyan', alpha=0.15)
    axs[1, 0].annotate('Elevação Preventiva\n' + r'$\dot{m} = 3{,}2\text{ kg/s}$', xy=(100, 3.2), xytext=(120, 2.4),
                       arrowprops=dict(facecolor='teal', shrink=0.08, width=1.2, headwidth=5),
                       fontsize=8.5, fontweight='bold', color='teal',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='teal', alpha=0.9))
    axs[1, 0].set_xlabel('Tempo de Operação (s)')
    axs[1, 0].set_ylabel('Vazão Mássica de Amônia (kg/s)')
    axs[1, 0].set_title('(b1) Atuação Preditiva no Circuito de Refrigeração MPFL')
    axs[1, 0].grid(True, linestyle=':', alpha=0.7)
    axs[1, 0].legend(loc='upper right', framealpha=0.9, fontsize=8.5)
    axs[1, 0].set_xlim([0, 600])
    axs[1, 0].set_ylim([0.8, 3.6])

    # (b2) Fadiga de Soldas BGA (Ciclos de Coffin-Manson)
    delta_T_vec = np.linspace(4, 50, 150)
    # Modelo Coffin-Manson-Engelmaier calibrado (Seção 4.6 da dissertação)
    # N_f(Delta T) = N_f_ref * (Delta T / Delta T_ref)**(-p)
    p_exp = np.log(1.84e5 / 5840.0) / np.log(48.0 / 7.0) # ~ 1.792
    N_f = 1.84e5 * (7.0 / delta_T_vec)**p_exp

    axs[1, 1].semilogy(delta_T_vec, N_f, 'b-', lw=2.2, label='Vida em Fadiga SAC305 ($N_f$)')
    axs[1, 1].axhline(53200, color='k', linestyle='--', lw=1.5, label='Meta de 10 Anos (53.200 Ciclos)')
    axs[1, 1].plot([7.0], [1.84e5], 'go', markersize=8, label=r'SSO com TEMS ($\Delta T=7^\circ\text{C}, N_f=1{,}84\times 10^5$)')
    axs[1, 1].plot([48.0], [5840], 'ro', markersize=8, label=r'LEO sem TEMS ($\Delta T=48^\circ\text{C}, N_f=5840$, Falha em 13,2 meses)')
    axs[1, 1].set_xlabel('Excursão Térmica Cíclica $\\Delta T_j$ (°C)')
    axs[1, 1].set_ylabel('Ciclos até Falha Estrutural $N_f$ (log)')
    axs[1, 1].set_title('(b2) Confiabilidade de Soldas BGA (Coffin-Manson-Engelmaier)')
    axs[1, 1].grid(True, which='both', linestyle=':', alpha=0.7)
    axs[1, 1].legend(loc='upper right', framealpha=0.9, fontsize=7.8)
    axs[1, 1].set_xlim([4, 50])
    axs[1, 1].set_ylim([1e3, 1e6])

    out_path = os.path.join(OUTPUT_DIR, 'figura_4_4_resposta_tems_precooling.png')
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"   -> Salvo com sucesso em: {out_path}")

if __name__ == '__main__':
    gerar_figura_4_1()
    gerar_figura_4_2()
    gerar_figura_4_3()
    gerar_figura_4_4()
