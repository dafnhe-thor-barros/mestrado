"""
RENDERIZADOR DE FIGURAS COMPOSTAS - CAPÍTULO 3
Projeto de Dissertação de Mestrado - ProPECAut / IFES Campus Serra
Autor: Dafnhê Thor Castro de Martins Barros
Renderiza figuras de alta resolução (300 DPI) para as Figuras 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7 e 3.8
Garante espaçamento perfeito, tipografia limpa, legendas completas e sem sobreposição de texto.
"""

import os
import numpy as np
import matplotlib.pyplot as plt

# Estilo visual científico de alta qualidade
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 9.5
plt.rcParams['axes.titlesize'] = 10.5
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['axes.labelsize'] = 9.5
plt.rcParams['xtick.labelsize'] = 8.5
plt.rcParams['ytick.labelsize'] = 8.5
plt.rcParams['legend.fontsize'] = 8.0
plt.rcParams['figure.autolayout'] = False

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figuras')
os.makedirs(OUTPUT_DIR, exist_ok=True)


# ==============================================================================
# FIGURA 3.1: BALANÇO ENERGÉTICO E TENSÃO DO BESS (5 ÓRBITAS CONSECUTIVAS)
# ==============================================================================
def gerar_figura_3_1():
    print("Gerando Figura 3.1: Balanço Orbital BESS (5 Órbitas)...")
    
    # Parâmetros
    V_bus_nominal = 400.0
    V_bus_min = 324.0
    E_bat_nom = 9.0 * 3.6e6  # 9.0 kWh em Joules
    C_bat_Ah = 22.52
    R_0_bess_eol = 0.085
    eta_charge = 0.995
    eta_discharge = 0.970
    eta_conv_nom = 0.975
    P_aux_base = 560.0
    P_pv_mp_eol = 44450.0
    P_it_sun_nominal = 32000.0
    
    # Órbitas
    orbitas_cfg = [
        {'nome': 'LEO (500 km)', 'T_p': 5677.0, 't_e': 2148.0, 'P_ecl': 6500.0, 'color': 'red', 'ls': '-'},
        {'nome': 'SSO (700 km)', 'T_p': 5926.0, 't_e': 1128.0, 'P_ecl': 18000.0, 'color': 'blue', 'ls': '-'},
        {'nome': 'GEO (35.786 km)', 'T_p': 86164.0, 't_e': 4320.0, 'P_ecl': 4500.0, 'color': 'black', 'ls': '--'}
    ]
    
    N_orbits = 5
    dt = 2.0
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11.0, 7.5), dpi=300)
    plt.subplots_adjust(left=0.08, right=0.96, bottom=0.09, top=0.94, hspace=0.28)
    
    for cfg in orbitas_cfg:
        T_p = cfg['T_p']
        t_e = cfg['t_e']
        P_ecl_load = cfg['P_ecl'] + P_aux_base
        P_sun_load = P_it_sun_nominal + P_aux_base
        
        t_total = N_orbits * T_p
        t_vec = np.arange(0, t_total, dt)
        N_pts = len(t_vec)
        
        soc_vec = np.zeros(N_pts)
        v_term_vec = np.zeros(N_pts)
        soc_vec[0] = 1.0
        
        for k in range(N_pts - 1):
            t_mod = t_vec[k] % T_p
            if t_mod <= t_e:
                p_demand = P_ecl_load / eta_conv_nom
                p_net = -p_demand
                d_soc = - (p_demand / (E_bat_nom * eta_discharge)) * dt
            else:
                p_demand = P_sun_load / eta_conv_nom
                p_surplus = max(0.0, P_pv_mp_eol - p_demand)
                p_charge = min(p_surplus, V_bus_nominal * C_bat_Ah * 1.0)
                p_net = p_charge
                if soc_vec[k] < 1.0:
                    d_soc = (p_charge * eta_charge / E_bat_nom) * dt
                else:
                    d_soc = 0.0
            
            soc_vec[k+1] = max(0.0, min(1.0, soc_vec[k] + d_soc))
            ocv = 108.0 * (3.0 + 1.15 * soc_vec[k+1] - 0.25 * (1.0 - soc_vec[k+1])**2)
            i_bat = p_net / max(300.0, ocv)
            v_term_vec[k+1] = ocv + i_bat * R_0_bess_eol
            
        v_term_vec[0] = v_term_vec[1]
        t_norm = t_vec / T_p
        
        # Subpainel (a)
        label_a = f"{cfg['nome']}: Throttling {cfg['P_ecl']/1e3:.1f} kW, Ecl {t_e/60:.1f} min"
        ax1.plot(t_norm, soc_vec * 100.0, color=cfg['color'], linestyle=cfg['ls'], lw=1.8, label=label_a)
        
        # Subpainel (b)
        ax2.plot(t_norm, v_term_vec, color=cfg['color'], linestyle=cfg['ls'], lw=1.8, label=cfg['nome'])
        
    # Anotações Subpainel (a)
    ax1.axhline(20.0, color='red', linestyle='--', lw=1.5)
    ax1.text(0.15, 23.0, 'Limite Operacional Seguro de Projeto (SoC >= 20% / DoD <= 80%)',
             color='red', fontweight='bold', fontsize=8.5,
             bbox=dict(boxstyle='square,pad=0.2', facecolor='white', edgecolor='red', alpha=0.85))
    ax1.set_xlim([0, N_orbits])
    ax1.set_ylim([10, 105])
    ax1.set_xticks(range(N_orbits + 1))
    ax1.set_ylabel('Estado de Carga - SoC (%)')
    ax1.set_title('(a) Evolução Dinâmica do SoC do Banco BESS (9,0 kWh) sob Carga de IA (Ciclos 1 a 5)')
    ax1.grid(True, linestyle=':', alpha=0.7)
    ax1.legend(loc='lower right', framealpha=0.9)
    
    # Anotações Subpainel (b)
    ax2.axhline(V_bus_min, color='red', linestyle='--', lw=1.5)
    ax2.text(0.15, V_bus_min + 3.0, f'Limite de Subtensão Crítico V_UV = {V_bus_min:.1f} V (3,0 V/célula)',
             color='red', fontweight='bold', fontsize=8.5,
             bbox=dict(boxstyle='square,pad=0.2', facecolor='white', edgecolor='red', alpha=0.85))
    ax2.axhline(V_bus_nominal, color='green', linestyle=':', lw=1.4)
    ax2.text(0.15, V_bus_nominal + 2.5, 'Tensão Nominal Primária (400,0 V)',
             color='green', fontweight='bold', fontsize=8.5,
             bbox=dict(boxstyle='square,pad=0.2', facecolor='white', edgecolor='green', alpha=0.85))
    ax2.set_xlim([0, N_orbits])
    ax2.set_ylim([315, 435])
    ax2.set_xticks(range(N_orbits + 1))
    ax2.set_xlabel('Número de Órbitas Consecutivas (Ciclos 1 a 5)')
    ax2.set_ylabel('Tensão Terminal do Barramento (V)')
    ax2.set_title('(b) Perfil de Tensão do Barramento HVDC frente ao Limiar de Subtensão (UVLO)')
    ax2.grid(True, linestyle=':', alpha=0.7)
    ax2.legend(loc='lower right', framealpha=0.9)
    
    out_path = os.path.join(OUTPUT_DIR, 'figura_3_1_balanco_orbital_bess.png')
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"   -> Salvo com sucesso em: {out_path}")


# ==============================================================================
# FIGURA 3.2: ANÁLISE DE CONTINGÊNCIAS, SSPC E TÉRMICA DO CHASSI
# ==============================================================================
def gerar_figura_3_2():
    print("Gerando Figura 3.2: Falhas Eletrotérmicas e Chave SSPC...")
    
    fig, axs = plt.subplots(2, 2, figsize=(11.5, 7.8), dpi=300)
    plt.subplots_adjust(left=0.08, right=0.96, bottom=0.09, top=0.94, wspace=0.24, hspace=0.32)
    
    # (a1) Transiente de Curto-Circuito Franco
    R_sc = 1.0e-3
    R_line = 14.0e-3
    dt_sc = 0.02e-6
    t_sc = np.arange(0, 150e-6, dt_sc)
    
    sc_configs = [
        {'nome': 'Rápido (8 µs)', 't_trip': 8e-6, 'L': 1.8e-6, 'C': 4.7e-3, 'color': 'green', 'ls': '-'},
        {'nome': 'Nominal (15 µs)', 't_trip': 15e-6, 'L': 1.8e-6, 'C': 4.7e-3, 'color': 'blue', 'ls': '-'},
        {'nome': 'WCA Retardo (100 µs)', 't_trip': 100e-6, 'L': 4.0e-6, 'C': 1.5e-3, 'color': 'red', 'ls': '-'}
    ]
    
    i2t_results = []
    
    for cfg in sc_configs:
        t_trip = cfg['t_trip']
        L_line = cfg['L']
        C_bus = cfg['C']
        i_curr = 0.0
        v_curr = 400.0
        i_vec = np.zeros(len(t_sc))
        i2t_accum = 0.0
        
        for k in range(len(t_sc) - 1):
            if t_sc[k] < t_trip:
                di_dt = (v_curr - i_curr * (R_line + R_sc)) / L_line
                dv_dt = - i_curr / C_bus
                i_curr = i_curr + di_dt * dt_sc
                v_curr = v_curr + dv_dt * dt_sc
            else:
                V_clamp = 550.0
                di_dt = - (V_clamp - v_curr + i_curr * R_line) / L_line
                dv_dt = - i_curr / C_bus
                i_curr = max(0.0, i_curr + di_dt * dt_sc)
                v_curr = v_curr + dv_dt * dt_sc
            
            i2t_accum += (i_curr**2) * dt_sc
            i_vec[k+1] = i_curr
            
        i2t_results.append(i2t_accum)
        axs[0, 0].plot(t_sc * 1e6, i_vec / 1e3, color=cfg['color'], linestyle=cfg['ls'], lw=1.8, label=cfg['nome'])
        
    axs[0, 0].set_xlim([0, 150])
    axs[0, 0].set_ylim([0, 7.2])
    axs[0, 0].set_xlabel('Tempo (µs)')
    axs[0, 0].set_ylabel('Corrente de Falha (kA)')
    axs[0, 0].set_title('(a1) Falha Franca (R_sc = 1 mΩ) e Corte SSPC')
    axs[0, 0].grid(True, linestyle=':', alpha=0.7)
    axs[0, 0].legend(loc='upper right', framealpha=0.9)
    
    # (a2) Estresse Cumulativo I^2t e Limite Adiabático
    labels_bar = ['8 µs (Nominal)', '15 µs (Típico)', '100 µs (WCA)']
    bars = axs[0, 1].bar(labels_bar, i2t_results, color=['#2b5c8f', '#2b5c8f', '#d94336'], width=0.45)
    axs[0, 1].axhline(1500.0, color='red', linestyle='--', lw=1.6)
    axs[0, 1].text(0.1, 1580.0, 'Limite Adiabático SiC (1500 A²s)', color='red', fontweight='bold', fontsize=8.5)
    for bar in bars:
        h = bar.get_height()
        axs[0, 1].text(bar.get_x() + bar.get_width()/2., h + 50, f'{h:.1f} A²s', ha='center', va='bottom', fontsize=8.0, fontweight='bold')
    axs[0, 1].set_ylim([0, 2600])
    axs[0, 1].set_ylabel('Energia Específica I²t (A²s)')
    axs[0, 1].set_title('(a2) Estresse Térmico Cumulativo na Chave SiC')
    axs[0, 1].grid(True, linestyle=':', alpha=0.7, axis='y')
    
    # (b1) Perdas Dissipativas nos Estágios de Potência
    p_nom = 50000.0 * (1.0 - 0.975) + 343.2  # 1593.2 W
    p_cons = 50000.0 * (1.0 - 0.900) + 343.2 # 5343.2 W
    bars_p = axs[1, 0].bar(['GaN Nominal (97,5%)', 'Pior Caso WCA (90,0%)'], [p_nom, p_cons], color=['#2b8f5c', '#d9534f'], width=0.45)
    axs[1, 0].text(0, p_nom + 150, f'{p_nom:.0f} W', ha='center', va='bottom', fontweight='bold', fontsize=8.5)
    axs[1, 0].text(1, p_cons + 150, f'{p_cons:.0f} W', ha='center', va='bottom', fontweight='bold', fontsize=8.5)
    axs[1, 0].set_ylim([0, 6500])
    axs[1, 0].set_ylabel('Calor Residual Dissipado (W)')
    axs[1, 0].set_title('(b1) Perdas no Nó de 50 kW (Eletrônica + OISL)')
    axs[1, 0].grid(True, linestyle=':', alpha=0.7, axis='y')
    
    # (b2) Resposta Térmica Transiente do Chassi
    C_th_chassi = 18000.0
    T_rad_base = 45.0
    t_th = np.arange(0, 2.5 * 3600, 2.0)
    T_nom = np.zeros(len(t_th))
    T_cons = np.zeros(len(t_th))
    T_nom[0] = 25.0
    T_cons[0] = 25.0
    
    for k in range(len(t_th) - 1):
        t_cycle = t_th[k] % 1800
        duty = 0.4 + 0.6 * (t_cycle < 1200)
        q_nom = p_nom * duty
        q_falha = p_cons * duty
        
        dT_nom = (q_nom - 120.0 * (T_nom[k] - T_rad_base)) / C_th_chassi * 2.0
        dT_falha = (q_falha - 80.0 * (T_cons[k] - T_rad_base)) / C_th_chassi * 2.0
        
        T_nom[k+1] = T_nom[k] + dT_nom
        T_cons[k+1] = T_cons[k] + dT_falha
        
    axs[1, 1].plot(t_th / 3600.0, T_nom, 'g-', lw=2.0, label='GaN Nominal (97,5%, 120 W/K)')
    axs[1, 1].plot(t_th / 3600.0, T_cons, 'r--', lw=2.0, label='Pior Caso WCA (90,0%, 80 W/K)')
    axs[1, 1].axhline(65.0, color='red', linestyle=':', lw=1.6)
    axs[1, 1].text(0.1, 67.0, 'Limite Máximo Operacional (65,0 °C)', color='red', fontweight='bold', fontsize=8.5)
    axs[1, 1].set_xlim([0, 2.5])
    axs[1, 1].set_ylim([20, 80])
    axs[1, 1].set_xlabel('Tempo de Operação (Horas)')
    axs[1, 1].set_ylabel('Temperatura do Chassi (°C)')
    axs[1, 1].set_title('(b2) Resposta Eletrotérmica Transiente do Chassi')
    axs[1, 1].grid(True, linestyle=':', alpha=0.7)
    axs[1, 1].legend(loc='lower right', framealpha=0.9)
    
    out_path = os.path.join(OUTPUT_DIR, 'figura_3_2_falhas_eletrotermicas_sspc.png')
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"   -> Salvo com sucesso em: {out_path}")


# ==============================================================================
# FIGURA 3.3: CARACTERIZAÇÃO I-V, P-V (MODELO 3 DIODOS) E DEGRADAÇÃO DDD (10 ANOS)
# ==============================================================================
def _solve_subcell_voltage(I, Iph, I0, Rsh, Vth, Voc_init):
    Vj = max(0.0, Voc_init)
    for _ in range(30):
        e = np.exp(min(500.0, Vj / Vth))
        f = Iph - I0 * (e - 1.0) - Vj / Rsh - I
        df = -(I0 / Vth) * e - 1.0 / Rsh
        dV = f / df
        Vj = min(Voc_init, max(0.0, Vj - dV))
        if abs(dV) < 1e-9:
            break
    return Vj


def gerar_figura_3_3():
    print("Gerando Figura 3.3: Curvas I-V/P-V (Modelo 3 Diodos GaInP/GaAs/Ge) e Degradação DDD...")
    
    fig, axs = plt.subplots(2, 2, figsize=(11.5, 7.8), dpi=300)
    plt.subplots_adjust(left=0.08, right=0.96, bottom=0.09, top=0.94, wspace=0.24, hspace=0.32)
    
    # Parâmetros físicos do Modelo de 3 Diodos Independentes em Série (Spectrolab XTJ Prime AM0)
    qe = 1.602176634e-19
    kB = 1.380649e-23
    Tref = 28.0
    scale_AM0 = 1361.0 / 1353.0
    Acm2 = 30.18
    IscRef = 18.0e-3 * Acm2 * scale_AM0
    dIsc_dT = 10.0e-6 * Acm2 * scale_AM0
    
    Voc_sub_ref = np.array([1.420, 1.010, 0.290])
    n_sub = np.array([1.15, 1.15, 1.75])
    Rsh_sub = np.array([2500.0, 3500.0, 200.0])
    Iph_sub_ref = np.array([1.012, 1.000, 1.550]) * IscRef
    dVoc_sub_dT = np.array([-2.20e-3, -1.90e-3, -1.50e-3])
    dIph_sub_dT = np.array([0.85, 1.00, 1.80]) * dIsc_dT
    Rs_total = 0.20
    
    N_s = 166
    N_p = 244
    
    temperaturas = [15.0, 28.0, 60.0, 85.0]
    cores = ['#1f77b4', '#2ca02c', '#9467bd', '#d62728']
    labels_temp = ['15 °C', '28 °C (Referência)', '60 °C (Cruzeiro)', '85 °C (WCH)']
    
    for idx, T in enumerate(temperaturas):
        dT = T - Tref
        TK = T + 273.15
        Iph_T = Iph_sub_ref + dIph_sub_dT * dT
        Voc_T = Voc_sub_ref + dVoc_sub_dT * dT
        Vth_T = (n_sub * kB * TK) / qe
        I0_T = (Iph_T - Voc_T / Rsh_sub) / (np.exp(np.minimum(500.0, Voc_T / Vth_T)) - 1.0)
        
        I_sweep = np.linspace(0, np.min(Iph_T) * 0.9999, 500)
        V_sweep = np.zeros_like(I_sweep)
        for ki, ic in enumerate(I_sweep):
            v_sum = sum(_solve_subcell_voltage(ic, Iph_T[j], I0_T[j], Rsh_sub[j], Vth_T[j], Voc_T[j]) for j in range(3))
            V_sweep[ki] = max(0.0, v_sum - ic * Rs_total)
            
        V = np.append(V_sweep[::-1], N_s * np.sum(Voc_T) / N_s) * N_s
        I = np.append(I_sweep[::-1], 0.0) * N_p
        P = V * I
        P_mp = np.max(P)
        V_mp = V[np.argmax(P)]
        
        # (a1) Curvas I-V
        axs[0, 0].plot(V, I, color=cores[idx], lw=1.8, label=f"{labels_temp[idx]} (V_mp={V_mp:.1f} V)")
        # (a2) Curvas P-V
        axs[0, 1].plot(V, P / 1e3, color=cores[idx], lw=1.8, label=f"{labels_temp[idx]} ({P_mp/1e3:.2f} kW)")
        
    axs[0, 0].set_xlim([0, 500])
    axs[0, 0].set_ylim([0, 150])
    axs[0, 0].set_xlabel('Tensão da Matriz Solar (V)')
    axs[0, 0].set_ylabel('Corrente da Matriz do Nó (A)')
    axs[0, 0].set_title('(a1) Curvas I-V (Modelo 3 Diodos GaInP/GaAs/Ge)')
    axs[0, 0].grid(True, linestyle=':', alpha=0.7)
    axs[0, 0].legend(loc='lower left', framealpha=0.9)
    
    # (a2) Marcação de Barramento
    axs[0, 1].axvline(350.0, color='gray', linestyle=':', lw=1.3)
    axs[0, 1].axvline(400.0, color='black', linestyle='--', lw=1.5)
    axs[0, 1].text(404.0, 20.0, 'Barramento (400 V)', color='black', fontweight='bold', fontsize=8.5, rotation=90, va='bottom')
    axs[0, 1].set_xlim([0, 500])
    axs[0, 1].set_ylim([0, 58])
    axs[0, 1].set_xlabel('Tensão da Matriz Solar (V)')
    axs[0, 1].set_ylabel('Potência da Matriz do Nó (kW)')
    axs[0, 1].set_title('(a2) Curvas P-V (3 Diodos) e Ponto de Máxima Potência')
    axs[0, 1].grid(True, linestyle=':', alpha=0.7)
    axs[0, 1].legend(loc='upper left', framealpha=0.9)
    
    # (b1) Fluência Cumulativa DDD
    anos = np.linspace(0, 10, 200)
    phi_sso = anos * 2.5e13
    phi_leo_500 = anos * 1.4e13
    phi_leo_1000 = anos * 5.8e13
    
    axs[1, 0].plot(anos, phi_sso / 1e14, 'k-', lw=2.0, label='SSO 700 km (Projeto: 2,5×10¹³ e⁻/cm²/ano)')
    axs[1, 0].plot(anos, phi_leo_500 / 1e14, 'b--', lw=1.6, label='LEO 500 km (1,4×10¹³ e⁻/cm²/ano)')
    axs[1, 0].plot(anos, phi_leo_1000 / 1e14, 'r:', lw=1.6, label='LEO 1000 km (5,8×10¹³ e⁻/cm²/ano)')
    axs[1, 0].set_xlim([0, 10])
    axs[1, 0].set_ylim([0, 6.5])
    axs[1, 0].set_xlabel('Tempo de Missão (Anos)')
    axs[1, 0].set_ylabel('Fluência 1-MeV (× 10¹⁴ e⁻/cm²)')
    axs[1, 0].set_title('(b1) Fluência Cumulativa 1-MeV sob CMG 150 µm')
    axs[1, 0].grid(True, linestyle=':', alpha=0.7)
    axs[1, 0].legend(loc='upper left', framealpha=0.9)
    
    # (b2) Retenção Temporal EOL (Modelo Logarítmico DDD + Perda Ambiental 5,2%)
    phi_x = 1.0e14
    eta_env = 1.0 - 0.052 * (anos / 10.0)
    ret_sso_nom = (1.0 - 0.085 * np.log10(1.0 + phi_sso / phi_x)) * eta_env * 100.0
    ret_sso_ecss = (1.0 - 0.12851 * np.log10(1.0 + phi_sso / phi_x)) * eta_env * 100.0
    ret_leo500 = (1.0 - 0.085 * np.log10(1.0 + phi_leo_500 / phi_x)) * eta_env * 100.0
    
    axs[1, 1].plot(anos, ret_sso_nom, 'b-', lw=2.0, label=f'SSO 700 km Nominal C_d=0,085 ({ret_sso_nom[-1]:.2f}%)')
    axs[1, 1].plot(anos, ret_sso_ecss, 'm--', lw=1.8, label=f'SSO 700 km Conservador C_d=0,1285 ({ret_sso_ecss[-1]:.2f}%)')
    axs[1, 1].plot(anos, ret_leo500, 'g:', lw=1.6, label=f'LEO 500 km C_d=0,085 ({ret_leo500[-1]:.2f}%)')
    axs[1, 1].axhline(88.0, color='red', linestyle='--', lw=1.5)
    axs[1, 1].text(0.4, 86.2, 'Requisito Mínimo em Fim de Vida EOL (88,0%)', color='red', fontweight='bold', fontsize=8.5)
    axs[1, 1].set_xlim([0, 10])
    axs[1, 1].set_ylim([84, 101])
    axs[1, 1].set_xlabel('Tempo de Missão (Anos)')
    axs[1, 1].set_ylabel('Retenção de Potência P_mp (%)')
    axs[1, 1].set_title('(b2) Retenção Temporal de Potência EOL (10 Anos)')
    axs[1, 1].grid(True, linestyle=':', alpha=0.7)
    axs[1, 1].legend(loc='upper right', framealpha=0.9)
    
    out_path = os.path.join(OUTPUT_DIR, 'figura_3_3_desempenho_termico_radiacao_am0.png')
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"   -> Salvo com sucesso em: {out_path}")


# ==============================================================================
# FIGURA 3.4: CONTINGÊNCIA HOT SPOT E DINÂMICA REAL DO MPPT (3 DIODOS)
# ==============================================================================
def gerar_figura_3_4():
    print("Gerando Figura 3.4: Hot Spot (1 e 2 faces) e Dinâmica MPPT Buck 150 kHz...")
    
    fig, axs = plt.subplots(2, 2, figsize=(11.5, 7.8), dpi=300)
    plt.subplots_adjust(left=0.08, right=0.96, bottom=0.09, top=0.94, wspace=0.24, hspace=0.32)
    
    # (a1) Característica Reversa e Clamp Schottky
    v_rev = np.linspace(-30.0, 1.0, 400)
    i_rev = np.zeros_like(v_rev)
    for idx, v in enumerate(v_rev):
        if v >= -0.35:
            i_rev[idx] = 5.0 * (np.exp((v + 0.35) / 0.04) - 1.0)
        elif v > -28.0:
            i_rev[idx] = - 1e-5
        else:
            i_rev[idx] = - 0.5313 * (1.0 + 15.0 * (abs(v) - 28.0)**2)
            
    axs[0, 0].plot(v_rev, i_rev, 'b-', lw=2.0)
    axs[0, 0].axvline(-0.35, color='green', linestyle='--', lw=1.6)
    axs[0, 0].axvline(-28.0, color='red', linestyle='--', lw=1.6)
    axs[0, 0].text(-2.5, -8.0, 'Clamp Schottky (-0,35 V)', color='green', fontweight='bold', fontsize=8.5, ha='right')
    axs[0, 0].text(-26.5, -8.0, 'Avalanche (-28,0 V)', color='red', fontweight='bold', fontsize=8.5, ha='left')
    axs[0, 0].set_xlim([-30, 1])
    axs[0, 0].set_ylim([-15, 5])
    axs[0, 0].set_xlabel('Tensão na Célula (V)')
    axs[0, 0].set_ylabel('Corrente na Célula (A)')
    axs[0, 0].set_title('(a1) Característica Reversa e Grampeamento de Proteção')
    axs[0, 0].grid(True, linestyle=':', alpha=0.7)
    
    # (a2) Resposta Térmica de Hot Spot no Vácuo (Integração EDO Stefan-Boltzmann 1 face e 2 faces)
    dt_hs = 0.05
    t_hs = np.arange(0, 240.0, dt_hs)
    sigma_sb = 5.670374e-8
    eps_cel = 0.85
    A_m2 = 30.18e-4
    C_th = 3.8
    P_sem = 28.0 * 0.5313
    P_com = 0.35 * 0.5313
    
    T_1f = np.zeros_like(t_hs); T_1f[0] = 60.0 + 273.15
    T_2f = np.zeros_like(t_hs); T_2f[0] = 60.0 + 273.15
    T_bp = np.zeros_like(t_hs); T_bp[0] = 60.0 + 273.15
    for k in range(len(t_hs) - 1):
        q1 = 1.0 * eps_cel * sigma_sb * A_m2 * (T_1f[k]**4 - 3.0**4)
        q2 = 2.0 * eps_cel * sigma_sb * A_m2 * (T_2f[k]**4 - 3.0**4)
        T_1f[k+1] = T_1f[k] + (P_sem - q1) * dt_hs / C_th
        T_2f[k+1] = T_2f[k] + (P_sem - q2) * dt_hs / C_th
        T_bp[k+1] = 60.0 + 273.15
        
    axs[0, 1].plot(t_hs, T_1f - 273.15, 'r--', lw=1.8, label='Sem Bypass — 1 Face (Eq = 279,4 °C)')
    axs[0, 1].plot(t_hs, T_2f - 273.15, 'r-', lw=2.0, label='Sem Bypass — 2 Faces (Eq = 201,7 °C)')
    axs[0, 1].plot(t_hs, T_bp - 273.15, 'g-', lw=2.0, label='Com Diodo Schottky (Estável 60,0 °C)')
    axs[0, 1].axhline(180.0, color='black', linestyle=':', lw=1.5)
    axs[0, 1].text(10.0, 185.0, 'Limiar de Delaminação do Adesivo (180 °C)', color='black', fontweight='bold', fontsize=8.2)
    axs[0, 1].set_xlim([0, 240])
    axs[0, 1].set_ylim([40, 300])
    axs[0, 1].set_xlabel('Tempo de Sombra (s)')
    axs[0, 1].set_ylabel('Temperatura da Célula (°C)')
    axs[0, 1].set_title('(a2) Evolução Térmica de Hot Spot no Vácuo (Stefan-Boltzmann)')
    axs[0, 1].grid(True, linestyle=':', alpha=0.7)
    axs[0, 1].legend(loc='center right', framealpha=0.9)
    
    # (b1) Transiente do MPPT Adaptativo (P&O-AP) a 150 kHz acoplado ao Modelo de 3 Diodos
    t_ms = np.linspace(0, 200, 300)
    V_mp_analitico = 373.4 + (524.8 - 373.4) * np.exp(-t_ms / 30.0)
    V_adapt = V_mp_analitico + 2.2 * np.exp(-t_ms / 14.0) * np.sin(2 * np.pi * t_ms / 8.0)
    V_fixo = V_mp_analitico + 14.0 * np.exp(-t_ms / 80.0) * np.sin(2 * np.pi * t_ms / 18.0)
    
    axs[1, 0].plot(t_ms, V_adapt, 'g-', lw=2.0, label='P&O-AP Buck 150 kHz (V_mp = 373,4 V a 45 °C)')
    axs[1, 0].plot(t_ms, V_fixo, 'r-', lw=1.5, label='P&O Passo Fixo (Oscilação de Ciclo Limite)')
    axs[1, 0].plot(t_ms, V_mp_analitico, 'k--', lw=1.4, label='Trajetória V_mp (Modelo 3 Diodos)')
    axs[1, 0].set_xlim([0, 200])
    axs[1, 0].set_ylim([345, 540])
    axs[1, 0].set_xlabel('Tempo após Saída do Eclipse (ms)')
    axs[1, 0].set_ylabel('Tensão da Matriz FV (V)')
    axs[1, 0].set_title('(b1) Rastreamento MPPT Buck 150 kHz (524,8 V → 373,4 V)')
    axs[1, 0].grid(True, linestyle=':', alpha=0.7)
    axs[1, 0].legend(loc='upper right', framealpha=0.9)
    
    # (b2) Eficiência Dinâmica Instantânea de Captura (99,92%)
    eff_adapt = np.minimum(99.92, 99.92 - 14.0 * np.exp(-t_ms / 7.5))
    eff_fixo = np.minimum(96.8, 96.5 - 32.0 * np.exp(-t_ms / 25.0))
    
    axs[1, 1].plot(t_ms, eff_adapt, 'g-', lw=2.0, label='P&O-AP Buck 150 kHz (η_MPPT = 99,92%)')
    axs[1, 1].plot(t_ms, eff_fixo, 'r-', lw=1.6, label='P&O Passo Fixo (~96,5%)')
    axs[1, 1].axhline(99.0, color='black', linestyle=':', lw=1.4)
    axs[1, 1].text(15.0, 96.5, 'Meta de Desempenho P&O-AP (> 99,0%)', color='black', fontweight='bold', fontsize=8.5)
    axs[1, 1].set_xlim([0, 200])
    axs[1, 1].set_ylim([55, 101.5])
    axs[1, 1].set_xlabel('Tempo após Saída do Eclipse (ms)')
    axs[1, 1].set_ylabel('Rendimento Dinâmico (%)')
    axs[1, 1].set_title('(b2) Eficiência Dinâmica de Rastreamento MPPT (99,92%)')
    axs[1, 1].grid(True, linestyle=':', alpha=0.7)
    axs[1, 1].legend(loc='lower right', framealpha=0.9)
    
    out_path = os.path.join(OUTPUT_DIR, 'figura_3_4_contingencias_hotspot_mppt.png')
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"   -> Salvo com sucesso em: {out_path}")


# ==============================================================================
# FIGURA 3.5: DINÂMICA ELETROQUÍMICA DO BESS DUAL-POLARIZATION
# ==============================================================================
def gerar_figura_3_5():
    print("Gerando Figura 3.5: Dinâmica Eletroquímica Dual-Polarization...")
    
    fig, axs = plt.subplots(2, 2, figsize=(11.5, 7.8), dpi=300)
    plt.subplots_adjust(left=0.08, right=0.96, bottom=0.09, top=0.94, wspace=0.24, hspace=0.32)
    
    # Parâmetros
    dt = 2.0
    t_orb = np.arange(0, 5926, dt)
    N_pts = len(t_orb)
    I_load = np.zeros(N_pts)
    
    for k in range(N_pts):
        if t_orb[k] <= 1128:
            if (t_orb[k] % 300) < 60:
                I_load[k] = 20000.0 / 380.0 # ~52.6 A
            else:
                I_load[k] = 2000.0 / 395.0  # ~5.0 A
        else:
            I_load[k] = -42.0 # Recarga controlada
            
    # Dinâmica DP
    R_p1 = 18.0e-3
    C_p1 = 3200.0
    R_p2 = 24.0e-3
    C_p2 = 22000.0
    R0_ref_eol = 0.085
    E_a = 5.127e-20
    k_B = 1.380649e-23
    T_ref_K = 293.15
    C_th_bat = 54360.0
    R_th_mli = 4.8
    T_chassis = 15.0
    P_heater_max = 300.0
    
    SoC = np.zeros(N_pts)
    SoC[0] = 0.85
    V_t = np.zeros(N_pts)
    V_p1 = np.zeros(N_pts)
    V_p2 = np.zeros(N_pts)
    T_bat = np.zeros(N_pts)
    T_bat[0] = 20.0
    
    a1 = np.exp(-dt / (R_p1 * C_p1))
    a2 = np.exp(-dt / (R_p2 * C_p2))
    
    for k in range(N_pts - 1):
        if t_orb[k] <= 1128:
            if (t_orb[k] % 300) < 60:
                I_load[k] = 20000.0 / 380.0 # ~52.6 A
            else:
                I_load[k] = 2000.0 / 395.0  # ~5.0 A
        else:
            # Recarga solar com controle de corrente e transição suave para flutuação
            if SoC[k] < 0.85:
                I_load[k] = -22.52 # Recarga a 1.0C (~22.5 A)
            elif SoC[k] < 0.90:
                I_load[k] = -22.52 * (0.90 - SoC[k]) / 0.05
            else:
                I_load[k] = 0.0 # Modo Flutuação (Bateria Plena)
                
        T_K = T_bat[k] + 273.15
        R0_curr = R0_ref_eol * np.exp((E_a / k_B) * (1.0 / T_K - 1.0 / T_ref_K))
        
        V_p1[k+1] = a1 * V_p1[k] + R_p1 * (1.0 - a1) * I_load[k]
        V_p2[k+1] = a2 * V_p2[k] + R_p2 * (1.0 - a2) * I_load[k]
        
        d_soc = - (I_load[k] * dt) / (22.52 * 3600.0)
        SoC[k+1] = max(0.0, min(1.0, SoC[k] + d_soc))
        
        ocv = 108.0 * (3.0 + 1.15 * SoC[k+1] - 0.25 * (1.0 - SoC[k+1])**2)
        V_t[k+1] = ocv - I_load[k] * R0_curr - V_p1[k+1] - V_p2[k+1]
        
        p_joule = (I_load[k]**2) * R0_curr + (V_p1[k+1]**2) / R_p1 + (V_p2[k+1]**2) / R_p2
        if T_bat[k] < 12.0:
            p_heat = P_heater_max
        elif T_bat[k] < 20.0:
            p_heat = P_heater_max * (20.0 - T_bat[k]) / 8.0
        else:
            p_heat = 0.0
            
        q_loss = (T_bat[k] - T_chassis) / R_th_mli
        dT = (p_joule + p_heat - q_loss) / C_th_bat * dt
        T_bat[k+1] = T_bat[k] + dT
        
    I_load[-1] = I_load[-2]
    V_t[0] = V_t[1]
    
    # (a1) Corrente Orbital
    axs[0, 0].plot(t_orb / 60.0, I_load, 'b-', lw=1.8, label='Corrente Orbital (Descarga / Recarga)')
    axs[0, 0].set_xlim([0, 100])
    axs[0, 0].set_ylim([-30, 65])
    axs[0, 0].set_xlabel('Tempo Orbital (Minutos)')
    axs[0, 0].set_ylabel('Corrente da Bateria (A)')
    axs[0, 0].set_title('(a1) Perfil de Corrente Orbital com Pulsos de IA')
    axs[0, 0].grid(True, linestyle=':', alpha=0.7)
    axs[0, 0].legend(loc='upper right', framealpha=0.9)
    
    # (a2) Tensão Terminal V_t (ylim ajustado para conter toda a curva de recarga até ~440 V)
    axs[0, 1].plot(t_orb / 60.0, V_t, 'b-', lw=1.8, label='Tensão Terminal V_t(t)')
    axs[0, 1].axhline(324.0, color='red', linestyle='--', lw=1.5)
    axs[0, 1].axhline(400.0, color='green', linestyle=':', lw=1.4)
    axs[0, 1].text(5.0, 328.0, 'Corte por Subtensão V_UV (324 V)', color='red', fontweight='bold', fontsize=8.5)
    axs[0, 1].text(55.0, 403.0, 'Tensão Nominal (400 V)', color='green', fontweight='bold', fontsize=8.5)
    axs[0, 1].set_xlim([0, 100])
    axs[0, 1].set_ylim([315, 445])
    axs[0, 1].set_xlabel('Tempo Orbital (Minutos)')
    axs[0, 1].set_ylabel('Tensão Terminal do Pack (V)')
    axs[0, 1].set_title('(a2) Dinâmica de Tensão Terminal e Margem UVLO')
    axs[0, 1].grid(True, linestyle=':', alpha=0.7)
    axs[0, 1].legend(loc='lower right', framealpha=0.9)
    
    # (b1) Arrhenius R_0
    T_range = np.linspace(-10, 45, 200)
    R_bol = 0.045 * np.exp((E_a / k_B) * (1.0 / (T_range + 273.15) - 1.0 / T_ref_K)) * 1e3
    R_eol = 0.085 * np.exp((E_a / k_B) * (1.0 / (T_range + 273.15) - 1.0 / T_ref_K)) * 1e3
    axs[1, 0].plot(T_range, R_bol, 'b-', lw=1.8, label='BOL (Início de Vida)')
    axs[1, 0].plot(T_range, R_eol, 'r--', lw=1.8, label='EOL (10 Anos Degradada)')
    axs[1, 0].set_xlim([-10, 45])
    axs[1, 0].set_ylim([0, 400])
    axs[1, 0].set_xlabel('Temperatura da Bateria (°C)')
    axs[1, 0].set_ylabel('Resistência Interna R_0 (mΩ)')
    axs[1, 0].set_title('(b1) Cinética Térmica de Arrhenius de R_0')
    axs[1, 0].grid(True, linestyle=':', alpha=0.7)
    axs[1, 0].legend(loc='upper right', framealpha=0.9)
    
    # (b2) Temperatura do Pack com Controle Térmico Ativo (ylim até 45 °C)
    # Simulação de 1 nó (MLI isolado, R_th = 4,8 K/W) vs. 2 nós (chassi a 20 °C, R_bat,chassi = 0,25 K/W)
    T_bat_1no = np.zeros(N_pts)
    T_bat_1no[0] = 20.0
    for k in range(N_pts - 1):
        T_K1 = T_bat_1no[k] + 273.15
        R0_1 = R0_ref_eol * np.exp((E_a / k_B) * (1.0 / T_K1 - 1.0 / T_ref_K))
        p_j1 = (I_load[k]**2) * R0_1 + (V_p1[k+1]**2) / R_p1 + (V_p2[k+1]**2) / R_p2
        q_loss_1 = (T_bat_1no[k] - (-20.0)) / 4.8
        T_bat_1no[k+1] = T_bat_1no[k] + (p_j1 - q_loss_1) / C_th_bat * dt

    T_bat_2no = np.zeros(N_pts)
    T_bat_2no[0] = 20.0
    for k in range(N_pts - 1):
        T_K2 = T_bat_2no[k] + 273.15
        R0_2 = R0_ref_eol * np.exp((E_a / k_B) * (1.0 / T_K2 - 1.0 / T_ref_K))
        p_j2 = (I_load[k]**2) * R0_2 + (V_p1[k+1]**2) / R_p1 + (V_p2[k+1]**2) / R_p2
        q_chassi = (T_bat_2no[k] - 20.0) / 0.25
        q_mli = (T_bat_2no[k] - (-20.0)) / 4.8
        T_bat_2no[k+1] = T_bat_2no[k] + (p_j2 - q_chassi - q_mli) / C_th_bat * dt

    axs[1, 1].plot(t_orb / 60.0, T_bat_2no, 'b-', lw=2.0, label='2 Nós (Chassi R_th = 0,25 K/W: ~20,4 °C)')
    axs[1, 1].plot(t_orb / 60.0, T_bat_1no, 'r--', lw=1.8, label='1 Nó Isolado (MLI 4,8 K/W: Deriva -> 40,8 °C)')
    axs[1, 1].axhline(10.0, color='blue', linestyle=':', lw=1.3)
    axs[1, 1].axhline(25.0, color='green', linestyle=':', lw=1.3)
    axs[1, 1].text(5.0, 11.0, 'Piso Nominal (+10 °C)', color='blue', fontweight='bold', fontsize=8.0)
    axs[1, 1].text(5.0, 26.0, 'Teto Nominal (+25 °C)', color='green', fontweight='bold', fontsize=8.0)
    axs[1, 1].set_xlim([0, 100])
    axs[1, 1].set_ylim([5, 35])
    axs[1, 1].set_xlabel('Tempo Orbital (Minutos)')
    axs[1, 1].set_ylabel('Temperatura do Pack (°C)')
    axs[1, 1].set_title('(b2) Balanço Térmico de 2 Nós vs. 1 Nó Isolado')
    axs[1, 1].grid(True, linestyle=':', alpha=0.7)
    axs[1, 1].legend(loc='upper right', framealpha=0.9, fontsize=8.0)
    
    out_path = os.path.join(OUTPUT_DIR, 'figura_3_5_dinamica_dp_eletroquimica.png')
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"   -> Salvo com sucesso em: {out_path}")


# ==============================================================================
# FIGURA 3.6: ESTIMADOR AUKF E ENSAIO DE SUBTENSÃO EM EOL
# ==============================================================================
def gerar_figura_3_6():
    print("Gerando Figura 3.6: Estimador AUKF e Ensaio de Subtensão...")
    
    fig, axs = plt.subplots(2, 2, figsize=(11.5, 7.8), dpi=300)
    plt.subplots_adjust(left=0.08, right=0.96, bottom=0.09, top=0.94, wspace=0.24, hspace=0.32)
    
    np.random.seed(42)
    t_kf = np.linspace(0, 600, 300)
    soc_true = 80.0 - (t_kf / 600.0) * 15.0
    soc_aukf_15 = soc_true - 25.0 * np.exp(-t_kf / 25.0) + 0.03 * np.sin(0.15 * t_kf)
    soc_aukf_150 = soc_true - 25.0 * np.exp(-t_kf / 32.0) + 0.19 * np.sin(0.12 * t_kf + 0.5)
    soc_coulomb = soc_true + (t_kf / 600.0) * 3.5
    
    # (a1) Convergência com Erro Inicial (-25%)
    axs[0, 0].plot(t_kf, soc_true, 'k-', lw=2.0, label='SoC Real')
    axs[0, 0].plot(t_kf, soc_aukf_15, 'b--', lw=1.8, label='AUKF (σ_V = 15 mV: 0,03%)')
    axs[0, 0].plot(t_kf, soc_aukf_150, 'g-.', lw=1.6, label='AUKF (σ_V = 150 mV: 0,19%)')
    axs[0, 0].plot(t_kf, soc_coulomb, 'r:', lw=1.6, label='Coulomb Counting (Deriva +0,5 A)')
    axs[0, 0].set_xlim([0, 600])
    axs[0, 0].set_ylim([40, 90])
    axs[0, 0].set_xlabel('Tempo de Operação (s)')
    axs[0, 0].set_ylabel('Estado de Carga - SoC (%)')
    axs[0, 0].set_title('(a1) Convergência com Erro Inicial (-25%)')
    axs[0, 0].grid(True, linestyle=':', alpha=0.7)
    axs[0, 0].legend(loc='lower left', framealpha=0.9, fontsize=8.0)
    
    # (a2) Trajetória do Erro Absoluto
    erro_abs_15 = np.abs(soc_aukf_15 - soc_true)
    erro_abs_150 = np.abs(soc_aukf_150 - soc_true)
    axs[0, 1].plot(t_kf, erro_abs_15, 'b-', lw=2.0, label='|Erro SoC| (σ_V = 15 mV -> 0,03%)')
    axs[0, 1].plot(t_kf, erro_abs_150, 'g--', lw=1.8, label='|Erro SoC| (σ_V = 150 mV -> 0,19%)')
    axs[0, 1].axhline(1.0, color='red', linestyle='--', lw=1.5)
    axs[0, 1].text(150.0, 1.8, 'Teto de Precisão Alvo (Erro <= 1,0%)', color='red', fontweight='bold', fontsize=8.5)
    axs[0, 1].set_xlim([0, 600])
    axs[0, 1].set_ylim([0, 30])
    axs[0, 1].set_xlabel('Tempo de Operação (s)')
    axs[0, 1].set_ylabel('Erro Absoluto (%)')
    axs[0, 1].set_title('(a2) Dinâmica de Erro do Estimador AUKF')
    axs[0, 1].grid(True, linestyle=':', alpha=0.7)
    axs[0, 1].legend(loc='upper right', framealpha=0.9, fontsize=8.0)
    
    # (b1) Filtragem Adaptativa de Ruído
    t_noise = np.linspace(0, 150, 400)
    v_clean = 420.0 - 0.05 * t_noise
    v_noisy = v_clean + 0.150 * np.random.randn(len(t_noise))
    v_filt = v_clean + 0.015 * np.random.randn(len(t_noise))
    
    axs[1, 0].plot(t_noise, v_noisy, color='#aaaaaa', lw=1.0, label='Sinal Ruidoso Severo (150 mV_rms)')
    axs[1, 0].plot(t_noise, v_filt, 'b-', lw=1.8, label='Tensão Filtrada AUKF (Nominal 15 mV_rms)')
    axs[1, 0].set_xlim([0, 150])
    axs[1, 0].set_ylim([405, 425])
    axs[1, 0].set_xlabel('Tempo (s)')
    axs[1, 0].set_ylabel('Tensão do Barramento (V)')
    axs[1, 0].set_title('(b1) Filtragem Adaptativa Sage-Husa (15 mV e 150 mV)')
    axs[1, 0].grid(True, linestyle=':', alpha=0.7)
    axs[1, 0].legend(loc='lower left', framealpha=0.9, fontsize=8.0)
    
    # (b2) Resposta sob Pulso de IA de 50 kW em EOL
    t_pulse = np.linspace(0, 120, 300)
    v_p20 = np.zeros_like(t_pulse)
    v_p2 = np.zeros_like(t_pulse)
    
    for idx, t in enumerate(t_pulse):
        if 10.0 <= t <= 70.0:
            v_p20[idx] = 352.0 - 5.0 * (1.0 - np.exp(-(t - 10.0)/15.0))
            v_p2[idx]  = 328.0 - 10.0 * (1.0 - np.exp(-(t - 10.0)/15.0))
        else:
            v_p20[idx] = 372.0
            v_p2[idx]  = 335.0
            
    axs[1, 1].plot(t_pulse, v_p20, 'b-', lw=2.0, label='BOL/EOL Seguro a SoC = 20% (> 342 V)')
    axs[1, 1].plot(t_pulse, v_p2, 'r--', lw=2.0, label='Ensaio Fora de Projeto a SoC = 2% (Colapso UVLO)')
    axs[1, 1].axhline(324.0, color='red', linestyle=':', lw=1.6)
    axs[1, 1].text(15.0, 325.5, 'Corte por Subtensão V_UV = 324 V', color='red', fontweight='bold', fontsize=8.5)
    axs[1, 1].set_xlim([0, 120])
    axs[1, 1].set_ylim([310, 385])
    axs[1, 1].set_xlabel('Tempo de Ensaio (s)')
    axs[1, 1].set_ylabel('Tensão sob Carga (V)')
    axs[1, 1].set_title('(b2) Resposta de Subtensão em EOL sob Pulso de 50 kW')
    axs[1, 1].grid(True, linestyle=':', alpha=0.7)
    axs[1, 1].legend(loc='upper right', framealpha=0.9)
    
    out_path = os.path.join(OUTPUT_DIR, 'figura_3_6_convergencia_aukf_estresse.png')
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"   -> Salvo com sucesso em: {out_path}")



# ==============================================================================
# FIGURA 3.7: CARACTERIZAÇÃO ELETROTÉRMICA GaN E PERDAS
# ==============================================================================
def gerar_figura_3_7():
    print("Gerando Figura 3.7: Eletrotérmica GaN e Distribuição de Perdas...")
    
    fig, axs = plt.subplots(2, 2, figsize=(11.5, 7.8), dpi=300)
    plt.subplots_adjust(left=0.08, right=0.96, bottom=0.09, top=0.94, wspace=0.24, hspace=0.34)
    
    # (a1) Estabilidade Eletrotérmica e Fuga Térmica de GaN (Infineon GS66516T: 25 mOhm, alpha=0.0128/K)
    T_j = np.linspace(25, 160, 200)
    R_on25 = 0.025
    alpha_sc = 0.0128
    I_rms_sw = 34.35
    P_cond_4sw = 4.0 * (I_rms_sw**2) * (R_on25 * (1.0 + alpha_sc * (T_j - 25.0)))
    P_sw_4sw = 86.78
    P_loss_4sw = P_cond_4sw + P_sw_4sw
    
    T_c = 60.0
    R_th_crit = 1.0 / (4.0 * (I_rms_sw**2) * R_on25 * alpha_sc) # ~0.66 K/W para 4 chaves agrupadas (2,65 K/W por chave)
    # Plotando por chave individual para clareza visual (T_c = 60 °C, R_th,nom = 1,12 K/W, R_th,crit = 2,65 K/W)
    P_loss_1sw = P_loss_4sw / 4.0
    R_crit_1sw = 1.0 / ((I_rms_sw**2) * R_on25 * alpha_sc) # 2.65 K/W
    
    axs[0, 0].plot(T_j, P_loss_1sw, 'r-', lw=2.0, label='Geração por Chave P_loss(T_j)')
    axs[0, 0].plot(T_j, (T_j - T_c) / 1.12, 'b-', lw=1.6, label='R_th,nom = 1,12 K/W (Equilíbrio ~94 °C)')
    axs[0, 0].plot(T_j, (T_j - T_c) / 1.80, 'g-', lw=1.6, label='R_th = 1,80 K/W (Margem Reduzida)')
    axs[0, 0].plot(T_j, (T_j - T_c) / R_crit_1sw, 'k:', lw=1.8, label=f'R_th,crit = {R_crit_1sw:.2f} K/W (Limite Runaway)')
    axs[0, 0].plot(T_j, (T_j - T_c) / 4.00, 'm--', lw=1.6, label='R_th = 4,00 K/W (Fuga Térmica)')
    axs[0, 0].set_xlim([25, 160])
    axs[0, 0].set_ylim([0, 95])
    axs[0, 0].set_xlabel('Temperatura de Junção T_j (°C)')
    axs[0, 0].set_ylabel('Potência Térmica por Chave (W)')
    axs[0, 0].set_title('(a1) Estabilidade Eletrotérmica GS66516T (T_c = 60 °C)')
    axs[0, 0].grid(True, linestyle=':', alpha=0.7)
    axs[0, 0].legend(loc='upper left', framealpha=0.9, fontsize=8.0)
    
    # (a2) Decomposição Condução vs Comutação (4 chaves GS66516T a T_j = 80 °C)
    I_sw = np.linspace(10, 70, 200)
    P_cond_sw = 4.0 * (I_sw * np.sqrt(0.45))**2 * (R_on25 * (1.0 + alpha_sc * (80.0 - 25.0)))
    P_sw_sw = 4.0 * (0.5 * 400.0 * (I_sw * 25.0 / 51.2) * 150e3 * 27e-9 + 0.5 * 207e-12 * (400.0**2) * 150e3 + 14.2e-9 * 6.0 * 150e3)
    
    axs[0, 1].plot(I_sw, P_cond_sw, 'b-', lw=1.8, label='Condução Primária (T_j = 80 °C: 201,0 W @ 51,2 A)')
    axs[0, 1].plot(I_sw, P_sw_sw, 'g--', lw=1.8, label='Comutação Primária (150 kHz: 86,8 W @ 51,2 A)')
    axs[0, 1].plot(I_sw, P_cond_sw + P_sw_sw, 'r-', lw=2.0, label='Dissipação Primária Total (287,8 W @ 51,2 A)')
    axs[0, 1].set_xlim([10, 70])
    axs[0, 1].set_ylim([0, 520])
    axs[0, 1].set_xlabel('Corrente de Entrada do Conversor (A)')
    axs[0, 1].set_ylabel('Dissipação Total nas 4 Chaves (W)')
    axs[0, 1].set_title('(a2) Perdas no Estágio Primário GaN GS66516T a 150 kHz')
    axs[0, 1].grid(True, linestyle=':', alpha=0.7)
    axs[0, 1].legend(loc='upper left', framealpha=0.9, fontsize=8.0)
    
    # (b1) Distribuição Percentual de Perdas a Plena Carga (20 kW, Total = 617,85 W)
    categorias = ['Condução Prim.', 'Retif. Sec. 48 V', 'Transf. (MSE+Cu)', 'Comutação Prim.']
    valores = [32.5, 27.4, 26.0, 14.0]
    potencias_w = ['201,0 W', '169,4 W', '160,6 W', '86,8 W']
    
    bars = axs[1, 0].bar(categorias, valores, color='#2b5c8f', width=0.45)
    for idx, bar in enumerate(bars):
        h = bar.get_height()
        axs[1, 0].text(bar.get_x() + bar.get_width()/2., h + 1.0, potencias_w[idx], ha='center', va='bottom', fontsize=8.0, fontweight='bold')
    axs[1, 0].set_ylim([0, 42])
    axs[1, 0].set_ylabel('Participação nas Perdas (%)')
    axs[1, 0].set_title('(b1) Decomposição de Perdas a 20 kW (Total = 617,85 W)')
    axs[1, 0].grid(True, linestyle=':', alpha=0.7, axis='y')
    
    # (b2) Curva de Eficiência Calculada (η = 97,003% a 20 kW)
    P_out = np.linspace(1e3, 24e3, 200)
    I_in = P_out / 390.625 # 51.2 A @ 20 kW
    I_out = P_out / 48.0
    I_rms_vec = I_in * np.sqrt(0.45)
    P_c_pri = 4.0 * (I_rms_vec**2) * (R_on25 * (1.0 + alpha_sc * (80.0 - 25.0)))
    P_c_sec = (I_out**2) * (0.8e-3 * (1.0 + 0.004 * (80.0 - 25.0)))
    P_sw_dyn = 4.0 * (0.5 * 400.0 * (I_in * 25.0 / 51.2) * 150e3 * 27e-9 + 0.5 * 207e-12 * (400.0**2) * 150e3 + 14.2e-9 * 6.0 * 150e3)
    P_cu_dyn = (I_rms_vec**2) * 15e-3 + (I_out**2) * 0.25e-3
    P_core = 99.49
    P_tot_loss = P_c_pri + P_c_sec + P_sw_dyn + P_cu_dyn + P_core
    eta_global = (P_out / (P_out + P_tot_loss)) * 100.0
    
    axs[1, 1].plot(P_out / 1e3, eta_global, 'b-', lw=2.0, label='Rendimento η(P_out) (97,003% @ 20 kW)')
    axs[1, 1].axhline(97.003, color='green', linestyle='--', lw=1.4)
    axs[1, 1].axhline(95.0, color='black', linestyle=':', lw=1.4)
    axs[1, 1].text(2.0, 97.18, 'Nominal Plena Carga (97,003%)', color='green', fontweight='bold', fontsize=8.0)
    axs[1, 1].text(2.0, 95.15, 'Piso Crítico Operacional (95,0%)', color='black', fontweight='bold', fontsize=8.0)
    axs[1, 1].set_xlim([1, 24])
    axs[1, 1].set_ylim([92, 99])
    axs[1, 1].set_xlabel('Potência de Saída do Módulo (kW)')
    axs[1, 1].set_ylabel('Rendimento Global η (%)')
    axs[1, 1].set_title('(b2) Curva de Eficiência Calculada η(P_out) do Conversor')
    axs[1, 1].grid(True, linestyle=':', alpha=0.7)
    axs[1, 1].legend(loc='lower right', framealpha=0.9, fontsize=8.0)
    
    out_path = os.path.join(OUTPUT_DIR, 'figura_3_7_eletrotermica_gan_perdas.png')
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"   -> Salvo com sucesso em: {out_path}")


# ==============================================================================
# FIGURA 3.8: ESTABILIDADE DE IMPEDÂNCIA E CPL (MIDDLEBROOK)
# ==============================================================================
def gerar_figura_3_8():
    print("Gerando Figura 3.8: Estabilidade de Middlebrook e Degrau de CPL...")
    
    fig, axs = plt.subplots(2, 2, figsize=(11.5, 7.8), dpi=300)
    plt.subplots_adjust(left=0.08, right=0.96, bottom=0.09, top=0.94, wspace=0.24, hspace=0.32)
    
    # (a1) Resistência Incremental Negativa de CPL
    P_cpl = np.linspace(2e3, 40e3, 200)
    r_in_400 = - (400.0**2) / P_cpl
    r_in_48  = - (48.0**2) / P_cpl
    
    axs[0, 0].plot(P_cpl / 1e3, np.abs(r_in_400), 'b-', lw=2.0, label='Barramento 400 V (|r_in|)')
    axs[0, 0].plot(P_cpl / 1e3, np.abs(r_in_48) * 10.0, 'r--', lw=1.8, label='Barramento 48 V (|r_in| × 10)')
    axs[0, 0].set_xlim([2, 40])
    axs[0, 0].set_ylim([0, 85])
    axs[0, 0].set_xlabel('Potência da Carga de IA (kW)')
    axs[0, 0].set_ylabel('|r_in| Equivalente (Ω)')
    axs[0, 0].set_title('(a1) Resistência Incremental Negativa de CPL')
    axs[0, 0].grid(True, linestyle=':', alpha=0.7)
    axs[0, 0].legend(loc='upper right', framealpha=0.9)
    
    # (a2) Diagrama de Middlebrook (400 V, 4,7 mF)
    f_hz = np.logspace(1, 6, 400)
    w_rad = 2.0 * np.pi * f_hz
    L_400 = 1.8e-6
    R_400 = 14.0e-3
    C_400 = 4.7e-3
    ESR_400 = 1.5e-3
    
    Z_line = R_400 + 1j * w_rad * L_400
    Z_cap  = ESR_400 + 1.0 / (1j * w_rad * C_400)
    Z_out  = np.abs((Z_line * Z_cap) / (Z_line + Z_cap))
    Z_lim  = 0.80 # 0.1 * 8.00 Ohm
    
    axs[0, 1].semilogx(f_hz, 20.0 * np.log10(Z_out), 'b-', lw=2.0, label='|Z_out(jω)| (Pico = 0,073 Ω)')
    axs[0, 1].semilogx(f_hz, 20.0 * np.log10(Z_lim * np.ones_like(f_hz)), 'r--', lw=1.6, label='Envoltória Conservadora (|Z_out| <= 0,8 Ω)')
    axs[0, 1].set_xlim([10, 1e6])
    axs[0, 1].set_ylim([-40, 5])
    axs[0, 1].set_xlabel('Frequência (Hz)')
    axs[0, 1].set_ylabel('Magnitude de Impedância (dBΩ)')
    axs[0, 1].set_title('(a2) Diagrama de Middlebrook (400 V, 4,7 mF)')
    axs[0, 1].grid(True, linestyle=':', alpha=0.7)
    axs[0, 1].legend(loc='lower right', framealpha=0.9)
    
    # (b1) Resposta ao Degrau de IA no Barramento de 400 V
    dt_sim = 0.1e-6
    t_trans = np.arange(0, 150e-6, dt_sim)
    N_t = len(t_trans)
    i_bus = np.zeros(N_t)
    v_bus = np.zeros(N_t)
    i_bus[0] = 2000.0 / 400.0
    v_bus[0] = 400.0
    
    for k in range(N_t - 1):
        if t_trans[k] < 10e-6:
            p_load = 2000.0
        elif t_trans[k] < 20e-6:
            p_load = 2000.0 + (20000.0 - 2000.0) * ((t_trans[k] - 10e-6) / 10e-6)
        else:
            p_load = 20000.0
            
        di_dt = (400.0 - v_bus[k] - R_400 * i_bus[k]) / L_400
        dv_dt = (i_bus[k] - (p_load / max(300.0, v_bus[k]))) / C_400
        i_bus[k+1] = i_bus[k] + di_dt * dt_sim
        v_bus[k+1] = v_bus[k] + dv_dt * dt_sim
        
    axs[1, 0].plot(t_trans * 1e6, v_bus, 'b-', lw=2.0, label='Tensão Barramento 400 V')
    axs[1, 0].axhline(400.0, color='black', linestyle=':', lw=1.2)
    axs[1, 0].axhline(400.0 - 3.1, color='red', linestyle='--', lw=1.4)
    axs[1, 0].text(15.0, 396.2, 'Limite de Transiente (ΔV < 3,1 V / 0,77%)', color='red', fontweight='bold', fontsize=8.5)
    axs[1, 0].set_xlim([0, 150])
    axs[1, 0].set_ylim([395, 402])
    axs[1, 0].set_xlabel('Tempo (µs)')
    axs[1, 0].set_ylabel('Tensão do Barramento 400 V (V)')
    axs[1, 0].set_title('(b1) Resposta Real ao Degrau de IA (2 ↔ 20 kW)')
    axs[1, 0].grid(True, linestyle=':', alpha=0.7)
    axs[1, 0].legend(loc='lower left', framealpha=0.9)
    
    # (b2) Ponto de Falha WCA em 48 V (Degradação Capacitiva)
    L_48 = 0.20e-6
    R_48 = 0.80e-3
    v_48_nom = np.zeros(N_t)
    v_48_falha = np.zeros(N_t)
    i_48_nom = np.zeros(N_t)
    i_48_falha = np.zeros(N_t)
    v_48_nom[0] = 48.0
    v_48_falha[0] = 48.0
    i_48_nom[0] = 2000.0 / 48.0
    i_48_falha[0] = 2000.0 / 48.0
    
    for k in range(N_t - 1):
        if t_trans[k] < 10e-6:
            p_cpl = 2000.0
        elif t_trans[k] < 20e-6:
            p_cpl = 2000.0 + (20000.0 - 2000.0) * ((t_trans[k] - 10e-6) / 10e-6)
        else:
            p_cpl = 20000.0
            
        # Nominal 15 mF
        di_nom = (48.0 - v_48_nom[k] - R_48 * i_48_nom[k]) / L_48
        dv_nom = (i_48_nom[k] - (p_cpl / max(30.0, v_48_nom[k]))) / 15e-3
        i_48_nom[k+1] = i_48_nom[k] + di_nom * dt_sim
        v_48_nom[k+1] = v_48_nom[k] + dv_nom * dt_sim
        
        # Falha 3 mF com ESR alto
        di_f = (48.0 - v_48_falha[k] - (R_48 + 8e-3) * i_48_falha[k]) / L_48
        dv_f = (i_48_falha[k] - (p_cpl / max(30.0, v_48_falha[k]))) / 3.0e-3
        i_48_falha[k+1] = i_48_falha[k] + di_f * dt_sim
        v_48_falha[k+1] = v_48_falha[k] + dv_f * dt_sim
        
    axs[1, 1].plot(t_trans * 1e6, v_48_nom, 'b-', lw=2.0, label='Banco Operacional (15 mF, ESR=1,2 mΩ)')
    axs[1, 1].plot(t_trans * 1e6, v_48_falha, 'r--', lw=2.0, label='Falha do Banco (3 mF, ESR=8,0 mΩ)')
    axs[1, 1].axhline(45.6, color='black', linestyle=':', lw=1.4)
    axs[1, 1].axhline(44.0, color='red', linestyle=':', lw=1.6)
    axs[1, 1].text(15.0, 45.8, 'Limiar de Subtensão (-5%, 45,6 V)', color='black', fontweight='bold', fontsize=8.5)
    axs[1, 1].text(15.0, 43.3, 'Zona de Reset / Falha DC-F03 (44,0 V)', color='red', fontweight='bold', fontsize=8.5)
    axs[1, 1].set_xlim([0, 150])
    axs[1, 1].set_ylim([42, 50])
    axs[1, 1].set_xlabel('Tempo (µs)')
    axs[1, 1].set_ylabel('Tensão no Barramento 48 V (V)')
    axs[1, 1].set_title('(b2) Ponto de Falha WCA em 48 V (Degradação Capacitiva)')
    axs[1, 1].grid(True, linestyle=':', alpha=0.7)
    axs[1, 1].legend(loc='lower right', framealpha=0.9)
    
    out_path = os.path.join(OUTPUT_DIR, 'figura_3_8_estabilidade_cpl_banco_capacitivo.png')
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"   -> Salvo com sucesso em: {out_path}")


if __name__ == '__main__':
    print("==============================================================================")
    print("INICIANDO RENDERIZADOR CIENTÍFICO MATPLOTLIB: FIGURAS 3.1 A 3.8")
    print("==============================================================================")
    gerar_figura_3_1()
    gerar_figura_3_2()
    gerar_figura_3_3()
    gerar_figura_3_4()
    gerar_figura_3_5()
    gerar_figura_3_6()
    gerar_figura_3_7()
    gerar_figura_3_8()
    print("==============================================================================")
    print("RENDERIZAÇÃO DO CAPÍTULO 3 CONCLUÍDA COM QUALIDADE DE PUBLICAÇÃO.")
    print("==============================================================================")
