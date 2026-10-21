"""
Simulacao Numerica de Transiente de Curto-Circuito em Barramento HVDC (400 V)
Determinacao exata dos limites criticos de abertura da chave SSPC (SiC MOSFET):
1. Tempo critico de desarme (t_trip) para atingir o Limite Adiabatico I^2t = 1500 A^2s
2. Tempo de desarme e tempo de extincao (soft turn-off) para atingir o Limite de Ruptura Dieletrica V_spike = 1200 V
"""

import numpy as np
import matplotlib.pyplot as plt

def simular_curto(t_trip, L_line, C_bus, R_line=14.0e-3, R_sc=1.0e-3, V_bus=400.0, V_clamp=550.0, clamp_ativo=True, tau_ext=0.5e-6):
    dt = 0.01e-6 # 10 ns
    t_end = max(150e-6, t_trip + 40e-6)
    t = np.arange(0, t_end, dt)
    
    i_curr = 0.0
    v_curr = V_bus
    i2t = 0.0
    v_sw_max = 0.0
    i_peak = 0.0
    
    for tk in t:
        if tk < t_trip:
            di_dt = (v_curr - i_curr * (R_line + R_sc)) / L_line
            dv_dt = - i_curr / C_bus
            i_curr = i_curr + di_dt * dt
            v_curr = v_curr + dv_dt * dt
            v_sw = i_curr * 0.010 # Queda direta Rds_on
        else:
            if clamp_ativo:
                di_dt = - (V_clamp - v_curr + i_curr * R_line) / L_line
                dv_dt = - i_curr / C_bus
                i_curr = max(0.0, i_curr + di_dt * dt)
                v_curr = v_curr + dv_dt * dt
                v_sw = V_clamp if i_curr > 0 else v_curr
            else:
                di_dt = - i_curr / tau_ext
                v_indutivo = abs(L_line * di_dt)
                i_curr = max(0.0, i_curr + di_dt * dt)
                dv_dt = - i_curr / C_bus
                v_sw = v_curr + v_indutivo
                
        i2t += (i_curr**2) * dt
        if i_curr > i_peak:
            i_peak = i_curr
        if v_sw > v_sw_max:
            v_sw_max = v_sw
            
        if tk > t_trip and i_curr <= 0.0:
            break
            
    return i_peak, i2t, v_sw_max

print("="*75)
print("1. VARREDURA PARA O LIMITE ADIABATICO (I^2t = 1500 A^2s)")
print("="*75)

# Cenario WCA (L = 4.0 uH, C = 1.5 mF, R_line = 14 mOhm, R_sc = 1 mOhm)
L_wca = 4.0e-6
C_wca = 1.5e-3

# Cenario Nominal (L = 1.8 uH, C = 4.7 mF)
L_nom = 1.8e-6
C_nom = 4.7e-3

t_trips = np.linspace(1e-6, 120e-6, 600)

i2t_wca = []
i2t_nom = []
v_spike_unclamp_wca = []
v_spike_unclamp_nom = []

for tt in t_trips:
    _, i2t_w, _ = simular_curto(tt, L_wca, C_wca, clamp_ativo=True)
    _, i2t_n, _ = simular_curto(tt, L_nom, C_nom, clamp_ativo=True)
    i2t_wca.append(i2t_w)
    i2t_nom.append(i2t_n)
    
    # Sem snubber (tau_ext = 0.5 us)
    _, _, v_w = simular_curto(tt, L_wca, C_wca, clamp_ativo=False, tau_ext=0.5e-6)
    _, _, v_n = simular_curto(tt, L_nom, C_nom, clamp_ativo=False, tau_ext=0.5e-6)
    v_spike_unclamp_wca.append(v_w)
    v_spike_unclamp_nom.append(v_n)

i2t_wca = np.array(i2t_wca)
i2t_nom = np.array(i2t_nom)
v_spike_unclamp_wca = np.array(v_spike_unclamp_wca)
v_spike_unclamp_nom = np.array(v_spike_unclamp_nom)

# Encontrar cruze de 1500 A^2s no WCA
idx_crit_wca = np.where(i2t_wca >= 1500.0)[0][0]
t_crit_wca = t_trips[idx_crit_wca]
i_pk_crit_wca, i2t_crit_wca, _ = simular_curto(t_crit_wca, L_wca, C_wca, clamp_ativo=True)

# Encontrar cruze de 1500 A^2s no Nominal
idx_crit_nom = np.where(i2t_nom >= 1500.0)[0]
if len(idx_crit_nom) > 0:
    t_crit_nom = t_trips[idx_crit_nom[0]]
else:
    t_crit_nom = None

print(f"CENARIO WCA (L = 4.0 uH, C = 1.5 mF):")
print(f"  -> Tempo critico para atingir I^2t = 1500 A^2s: t_trip = {t_crit_wca*1e6:.2f} s")
print(f"  -> Corrente de pico atingida: {i_pk_crit_wca/1e3:.2f} kA")
print(f"  -> I^2t exato calculado: {i2t_crit_wca:.1f} A^2s")

# Encontrar cruze de V_spike = 1200 V sem clamp (tau_ext = 0.5 us)
idx_v1200_wca = np.where(v_spike_unclamp_wca >= 1200.0)[0][0]
t_v1200_wca = t_trips[idx_v1200_wca]
i_pk_v1200_wca, _, v_sw_1200 = simular_curto(t_v1200_wca, L_wca, C_wca, clamp_ativo=False, tau_ext=0.5e-6)

idx_v1200_nom = np.where(v_spike_unclamp_nom >= 1200.0)[0][0]
t_v1200_nom = t_trips[idx_v1200_nom]
i_pk_v1200_nom, _, v_sw_1200_nom = simular_curto(t_v1200_nom, L_nom, C_nom, clamp_ativo=False, tau_ext=0.5e-6)

print("\n" + "="*75)
print("2. TEMPO DE RETARDO QUE ATINGE V_spike = 1200 V (SEM CLAMP, tau_ext = 0.5 us)")
print("="*75)
print(f"CENARIO WCA (L = 4.0 uH):")
print(f"  -> Tempo de retardo t_trip para V_spike = 1200 V: {t_v1200_wca*1e6:.3f} us ({t_v1200_wca*1e9:.1f} ns)")
print(f"  -> Corrente no momento da abertura: {i_pk_v1200_wca:.1f} A")
print(f"  -> Tensao de pico induzida V_sw: {v_sw_1200:.1f} V")

print(f"\nCENARIO NOMINAL (L = 1.8 uH):")
print(f"  -> Tempo de retardo t_trip para V_spike = 1200 V: {t_v1200_nom*1e6:.3f} us ({t_v1200_nom*1e9:.1f} ns)")
print(f"  -> Corrente no momento da abertura: {i_pk_v1200_nom:.1f} A")
print(f"  -> Tensao de pico induzida V_sw: {v_sw_1200_nom:.1f} V")

print("\n" + "="*75)
print("3. TEMPO DE EXTINCAO MINIMO (SOFT TURN-OFF) PARA EVITAR V_spike > 1200 V")
print("="*75)
# Se t_trip = 15 us (nominal) ou 100 us (WCA), qual deve ser o tau_ext minimo sem TVS?
# V_spike = V_bus + L * (i_trip / tau_ext) <= 1200 V
# => L * i_trip / tau_ext <= 1200 - 400 = 800 V
# => tau_ext >= L * i_trip / 800 V

i_pk_15us_nom, _, _ = simular_curto(15e-6, L_nom, C_nom)
tau_min_15us_nom = L_nom * i_pk_15us_nom / (1200.0 - 400.0)

i_pk_15us_wca, _, _ = simular_curto(15e-6, L_wca, C_wca)
tau_min_15us_wca = L_wca * i_pk_15us_wca / (1200.0 - 400.0)

i_pk_100us_wca, _, _ = simular_curto(100e-6, L_wca, C_wca)
tau_min_100us_wca = L_wca * i_pk_100us_wca / (1200.0 - 400.0)

print(f"Para desarme em t_trip = 15 us (Nominal, L = 1.8 uH, i_trip = {i_pk_15us_nom/1e3:.2f} kA):")
print(f"  -> Tempo de extincao suave minimo necessario: tau_ext >= {tau_min_15us_nom*1e6:.2f} s")

print(f"Para desarme em t_trip = 15 us (WCA, L = 4.0 uH, i_trip = {i_pk_15us_wca/1e3:.2f} kA):")
print(f"  -> Tempo de extincao suave minimo necessario: tau_ext >= {tau_min_15us_wca*1e6:.2f} s")

print(f"Para desarme em t_trip = 100 us (WCA, L = 4.0 uH, i_trip = {i_pk_100us_wca/1e3:.2f} kA):")
print(f"  -> Tempo de extincao suave minimo necessario: tau_ext >= {tau_min_100us_wca*1e6:.2f} s")
print(f"  -> Se forcado a desligar em tau_ext = 0.5 us: V_spike = 400 + 4e-6 * ({i_pk_100us_wca:.1f} / 0.5e-6) = {400 + L_wca * i_pk_100us_wca / 0.5e-6:.1f} V (11.99 kV!)")
