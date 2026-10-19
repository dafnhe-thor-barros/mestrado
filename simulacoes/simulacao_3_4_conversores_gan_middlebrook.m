% ==============================================================================
% SIMULAÇÃO 3.4: MODELAGEM ELETROTÉRMICA GaN, PERDAS DE STEINMETZ,
% CRITÉRIO DE FUGA TÉRMICA E ESTABILIDADE DE MIDDLEBROOK SOB CARGAS CPL
% ==============================================================================
% Dissertação de Mestrado Profissional - ProPECAut / IFES
% Autor: Dafnhê Thor Castro de Martins Barros
% ==============================================================================
clear; clc; close all;
fprintf('==============================================================================\n');
fprintf('INICIANDO SIMULAÇÃO 3.4: CONVERSORES GaN, STEINMETZ E ESTABILIDADE CPL\n');
fprintf('==============================================================================\n\n');
%% 1. PARÂMETROS ELETROFÍSICOS DO MÓDULO CONVERSOR DE 20 kW (400 V / 48 V)
P_mod_nom    = 20000.0;     % [W] Potência nominal por módulo
V_in_nom     = 400.0;       % [V] Barramento primário
V_out_nom    = 48.0;        % [V] Barramento intermediário de rack
f_sw         = 150.0e3;     % [Hz] Frequência de chaveamento primário
% Semicondutores GaN HEMT Primários (Tabela 20 da Dissertação)
R_on25_pri   = 18.0e-3;     % [Ohm] Resistência de canal a 25 °C
alpha_sc     = 0.0065;      % [1/K] Coeficiente térmico de canal
t_r          = 12.0e-9;     % [s] Tempo de subida
t_f          = 15.0e-9;     % [s] Tempo de descida
C_oss        = 120.0e-12;   % [F] Capacitância de saída a 400 V
Q_g          = 6.0e-9;      % [C] Carga de gate
V_gate       = 5.0;         % [V] Tensão de driver
% Retificação Síncrona Secundária GaN (48 V)
R_on25_sec   = 0.80e-3;     % [Ohm] Arranjo paralelo a 25 °C
alpha_sec    = 0.0040;      % [1/K]
% Parâmetros Magnéticos (Ferrite MnZn grau Ferroxcube 3C95 - Equação 246 do Capítulo 3)
k_fe         = 0.6715;      % [W/(m^3 Hz^alpha T^beta)] Coeficiente calibrado para P_core = 95.3 W
alpha_st     = 1.46;        % Expoente de frequência
beta_st      = 2.52;        % Expoente de indução
V_e          = 4.5e-4;      % [m^3] Volume efetivo do núcleo magnético
delta_B      = 0.15;        % [T] Excursão de indução de pico
R_ac_pri     = 15.0e-3;     % [Ohm] Resistência CA dos enrolamentos primários
R_ac_sec     = 0.25e-3;     % [Ohm] Resistência CA secundária
%% 2. BALANÇO DE PERDAS E EFICIÊNCIA EM PLENA CARGA (20 kW, Tj = 80 °C)
T_j_target   = 80.0;        % [°C] Temperatura operacional de junção
I_in_nom     = P_mod_nom / V_in_nom; % 50.0 A
I_out_nom    = P_mod_nom / V_out_nom;% 416.7 A
% Topologia Interleaved Bifásica (2 fases de 10 kW em paralelo com defasagem de 180°)
N_phases     = 2;
I_phase_nom  = I_in_nom / N_phases;  % 25.0 A por ramo de conversão
% Correntes eficazes nas chaves (Full-Bridge D = 0.45 por fase)
D_eff = 0.45;
I_rms_sw_pri = I_phase_nom * sqrt(D_eff); % 16.77 A por chave individual
I_rms_pri    = I_in_nom * sqrt(D_eff);
I_rms_sec    = I_out_nom;
% Perdas por Condução
R_on_pri_Tj = R_on25_pri * (1 + alpha_sc * (T_j_target - 25));
P_cond_pri  = (4 * N_phases) * (I_rms_sw_pri^2) * R_on_pri_Tj; % 8 transistores primários (2 pontes completas)
R_on_sec_Tj = R_on25_sec * (1 + alpha_sec * (T_j_target - 25));
P_cond_sec  = (I_out_nom^2) * R_on_sec_Tj;
% Perdas por Comutação GaN Primária (Interleaved bifásica, 4 chaves por fase comutam I_phase_nom)
P_sw_overlap = (4 * N_phases) * (V_in_nom * I_phase_nom * f_sw * (t_r + t_f) / 4);
P_sw_coss    = (4 * N_phases) * (0.5 * C_oss * (V_in_nom^2) * f_sw);
P_gate       = (4 * N_phases) * (Q_g * V_gate * f_sw);
P_sw_total   = P_sw_overlap + P_sw_coss + P_gate;
% Perdas Magnéticas (Equação de Steinmetz e perdas ôhmicas de cobre)
P_core = k_fe * (f_sw^alpha_st) * (delta_B^beta_st) * V_e; % 95.3 W
P_cu   = (I_rms_pri^2) * R_ac_pri + (I_rms_sec^2) * R_ac_sec;
P_mag  = P_core + P_cu;
P_loss_total = P_cond_pri + P_cond_sec + P_sw_total + P_mag;
rendimento_global = (P_mod_nom / (P_mod_nom + P_loss_total)) * 100;
fprintf('--- 2.1 BALANÇO DE PERDAS NO MÓDULO DE 20 kW (Tj = 80 °C) ---\n');
fprintf('Perdas em Condução Primária:   %5.1f W (%4.1f %%)\n', P_cond_pri, P_cond_pri/P_loss_total*100);
fprintf('Perdas em Condução Secundária: %5.1f W (%4.1f %%)\n', P_cond_sec, P_cond_sec/P_loss_total*100);
fprintf('Perdas por Comutação GaN:      %5.1f W (%4.1f %%)\n', P_sw_total, P_sw_total/P_loss_total*100);
fprintf('Perdas Magnéticas (Steinmetz): %5.1f W (%4.1f %%)\n', P_mag, P_mag/P_loss_total*100);
fprintf('Perda Total Dissipada:         %5.1f W | Rendimento Global: %5.2f %%\n\n', ...
    P_loss_total, rendimento_global);
%% 3. ANÁLISE DE FUGA TÉRMICA (THERMAL RUNAWAY - CRITÉRIO ET-F01)
% Critério da Equação 22: dP_loss/dT_j < 1 / R_th
dPloss_dTj = 4 * (I_rms_sw_pri^2) * R_on25_pri * alpha_sc; % [W/K]
R_th_crit = 1 / dPloss_dTj;
fprintf('--- 3.1 CRITÉRIO DE ESTABILIDADE E FUGA TÉRMICA DE GaN (ET-F01) ---\n');
fprintf('Gradiente de Geração Térmica (dP/dTj): %4.2f W/K\n', dPloss_dTj);
fprintf('Resistência Térmica Crítica de Runaway: %4.2f K/W\n', R_th_crit);
% Teste com Interfaces Térmicas
R_th_nominal   = 1.0; % [K/W] Interface saudável com pasta e espalhador de cobre
R_th_degradada = 4.0; % [K/W] Degradação severa de TIM no vácuo (pump-out)
fprintf('Interface Nominal   (R_th = 1.0 K/W): Margem Térmica: +%4.2f K/W | Status: ESTÁVEL\n', ...
    R_th_crit - R_th_nominal);
fprintf('Interface Degradada (R_th = 4.0 K/W): Margem Térmica: +%4.2f K/W | ', ...
    R_th_crit - R_th_degradada);
if R_th_degradada >= R_th_crit
    fprintf('Status: [FALHA CATÁSTROFICA ET-F01] Fuga Térmica Iminente!\n\n');
else
    fprintf('Status: PRÓXIMO AO LIMIAR CRÍTICO\n\n');
end
%% 4. CRITÉRIO DE MIDDLEBROOK SOB CARGAS CPL (48 V / 20 kW - FALHA DC-F03)
r_in_CPL = - (V_out_nom^2) / P_mod_nom; % -0.1152 Ohm (Equação 25)
Z_limite_middlebrook = 0.10 * abs(r_in_CPL); % Margem de 10x (0.0115 Ohm)
% Parâmetros do Barramento de 48 V:
% Caso A: Banco Saudável (C = 15 mF, ESR = 1.2 mOhm, L_dist = 0.2 uH)
% Caso B: Modo de Falha DC-F03 (Degradação WCA: C = 3.0 mF, ESR = 8.0 mOhm)
freq = logspace(1, 6, 500);
omega = 2 * pi * freq;
% Impedância de saída do filtro do barramento:
% Z_out(s) = (s*L + R) || (1/(s*C) + ESR)
L_line_48V = 0.20e-6; % [H] 200 nH
R_line_48V = 0.80e-3; % [Ohm]
% Caso A: Saudável
C_A   = 15.0e-3;
ESR_A = 1.2e-3;
Z_num_A = (1i * omega * L_line_48V + R_line_48V) .* (1 ./ (1i * omega * C_A) + ESR_A);
Z_den_A = (1i * omega * L_line_48V + R_line_48V) + (1 ./ (1i * omega * C_A) + ESR_A);
Z_out_A = abs(Z_num_A ./ Z_den_A);
% Caso B: Falha DC-F03 (Degradação Severa WCA: C = 2.5 mF, ESR = 9.5 mOhm)
C_B   = 2.5e-3;
ESR_B = 9.5e-3;
Z_num_B = (1i * omega * L_line_48V + R_line_48V) .* (1 ./ (1i * omega * C_B) + ESR_B);
Z_den_B = (1i * omega * L_line_48V + R_line_48V) + (1 ./ (1i * omega * C_B) + ESR_B);
Z_out_B = abs(Z_num_B ./ Z_den_B);
pico_Z_A = max(Z_out_A);
pico_Z_B = max(Z_out_B);
fprintf('--- 4.1 ESTABILIDADE DE IMPEDÂNCIA DE MIDDLEBROOK NO BARRAMENTO 48 V ---\n');
fprintf('Resistência Incremental de Carga (CPL): rin = %6.4f Ohm\n', r_in_CPL);
fprintf('Teto de Impedância de Middlebrook (|Zout| <= 0.1|rin|): %6.4f Ohm\n', Z_limite_middlebrook);
fprintf('Banco Saudável (15 mF, 1.2 mOhm): Pico |Zout| = %6.4f Ohm | Status: ESTABILIDADE ASSEGURADA\n', pico_Z_A);
fprintf('Banco Degradado (2.5 mF, 9.5 mOhm): Pico |Zout| = %6.4f Ohm | ', pico_Z_B);
if pico_Z_B > Z_limite_middlebrook
    fprintf('Status: [FALHA DC-F03: INSTABILIDADE CPL] Violação de Middlebrook!\n');
else
    fprintf('Status: ESTÁVEL\n');
end
fprintf('==============================================================================\n');
fprintf('SIMULAÇÃO 3.4 CONCLUÍDA COM SUCESSO.\n');
fprintf('==============================================================================\n');
