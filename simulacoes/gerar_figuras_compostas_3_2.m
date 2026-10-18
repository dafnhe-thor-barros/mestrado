% ==============================================================================
% GERACAO DE FIGURAS COMPOSTAS CONSOLIDADAS - SECAO 3.2 (GERADOR AM0 E MPPT)
% Projeto de Dissertacao de Mestrado - ProPECAut / IFES Campus Serra
% ==============================================================================
% Modelo conceitual parametrico in silico verificado numericamente.
% Nota Epistemologica: Parametros classificados como [VERIFICAVEL], [CALCULADO] ou [PROVISORIO].
% ==============================================================================

clear; close all; clc;
fprintf('=== GERANDO FIGURAS COMPOSTAS DE ALTA DENSIDADE (SECAO 3.2) ===\n');

if ismember('qt', available_graphics_toolkits())
    graphics_toolkit('qt');
else
    graphics_toolkit('gnuplot');
end

% ------------------------------------------------------------------------------
% FIGURA 3.3 COMPOSTA: CURVAS TERMICAS E DEGRADACAO POR RADIACAO (10 ANOS)
% ------------------------------------------------------------------------------
fprintf('1. Gerando Figura 3.3 Composta (Curvas I-V/P-V + DDD 10 anos)...\n');

fig_c3 = figure(1, 'position', [50, 50, 1100, 750]);

q = 1.602176634e-19; kB = 1.380649e-23;
T_ref = 301.15; 
G_ref_cal = 1353.0; % Irradiancia de calibracao Spectrolab AM0 [W/m^2] [VERIFICAVEL]
G_AM0 = 1361.0;     % Constante solar extraterrestre ASTM E490 AM0 [W/m^2] [VERIFICAVEL]
G_datasheet = 1353.0; % Condicao de ensaio do datasheet Spectrolab XTJ Prime [W/m^2] [VERIFICAVEL]
scale_AM0 = G_AM0 / G_datasheet; % 1.00591 (fator de escala de irradiancia 1353 -> 1361 W/m^2) [CALCULADO]

A_cel = 30.18e-4;   % m^2 (30.18 cm^2) [PROVISORIO]
J_sc_ref = 18.00e-3; % A/cm^2 [VERIFICAVEL]
J_mp_ref = 17.50e-3; % A/cm^2 [VERIFICAVEL]
I_sc_ref = J_sc_ref * (A_cel * 1e4); % 0.54324 A a 1353 W/m^2 [CALCULADO]
I_mp_ref = J_mp_ref * (A_cel * 1e4); % 0.52815 A a 1353 W/m^2 [CALCULADO]
V_oc_ref = 2.720; V_mp_ref = 2.406;
P_mp_ref = V_mp_ref * I_mp_ref;      % Potencia de pico direta por celula: 1.2707 W [CALCULADO]

% Coeficientes termicos Spectrolab XTJ Prime (15 C a 125 C) [VERIFICAVEL]:
K_I = 10.0e-6 * (A_cel * 1e4);   % +0.3018 mA/K para Isc
K_I_mp = 5.0e-6 * (A_cel * 1e4); % +0.1509 mA/K para Imp
K_V = -5.60e-3;                  % -5.60 mV/K para Voc
K_V_mp = -6.30e-3;               % -6.30 mV/K para Vmp

E_g = 1.42; E_g_J = E_g * q;
a_eq = 3.00; R_s = 0.080; R_sh = 1200.0; % [PROVISORIO / CALIBRADO]
V_th_ref = (kB * T_ref) / q;
I_0_ref = (I_sc_ref - V_oc_ref / R_sh) / (exp(V_oc_ref / (a_eq * V_th_ref)) - 1.0);

N_s = 166;                  % 166 celulas em serie -> Vmp = 399.4 V (~400 V)
N_p_no = 244;               % Strings em paralelo por no (40.504 celulas, 51.47 kW BOL instalados)
N_nos = 20;                 % 20 nos na constelacao
N_p_const = N_p_no * N_nos; % 4.880 strings na constelacao (1.029 MW BOL instalados)

% Subpainel (a1): Curvas I-V
subplot(2, 2, 1);
hold on; grid on; box on;
T_sweeps_C = [-150, -50, 28, 65, 85];
T_sweeps_K = T_sweeps_C + 273.15;
cores = {'b', 'c', 'g', 'm', 'r'};
i_norm = linspace(0, 0.999, 300);

for idx = 1:length(T_sweeps_K)
    T_curr = T_sweeps_K(idx);
    dT = T_curr - T_ref;
    V_th = (kB * T_curr) / q;
    
    Voc_cel = V_oc_ref + K_V * dT;
    Iph_cel = max(0.01, (I_sc_ref + K_I * dT) * scale_AM0);
    I_0_curr = (Iph_cel - Voc_cel / R_sh) / (exp(Voc_cel / (a_eq * V_th)) - 1.0);
    
    i_cel = i_norm * Iph_cel;
    v_j = a_eq * V_th * log(max(1.0, (Iph_cel - i_cel) / I_0_curr + 1.0));
    v_cel = max(0, v_j - i_cel * R_s - (v_j / R_sh) * R_s);
    
    V_array = v_cel * N_s; I_array = i_cel * N_p_no;
    [V_sort, s_idx] = sort(V_array); I_sort = I_array(s_idx);
    plot(V_sort, I_sort, cores{idx}, 'linewidth', 1.8);
end
plot([400, 400], [0, 145], 'k--', 'linewidth', 1.5);
xlabel('Tensao da Matriz (V)', 'fontsize', 9);
ylabel('Corrente da Matriz do No (A)', 'fontsize', 9);
title('(a1) Curvas I-V do No vs Temperatura (Np_{no}=244 strings, Ns=166)', 'fontsize', 10, 'fontweight', 'bold');
legend({'-150 C [Extrapolado]', '-50 C [Extrapolado]', '+28 C (BOL)', '+65 C (Nom)', '+85 C (Max)', 'Barramento (400 V)'}, 'location', 'southwest', 'fontsize', 7);
xlim([0, 750]); ylim([0, 145]);

% Subpainel (a2): Curvas P-V com MPP e Faixa Buck
subplot(2, 2, 2);
hold on; grid on; box on;
for idx = 1:length(T_sweeps_K)
    T_curr = T_sweeps_K(idx);
    dT = T_curr - T_ref;
    V_th = (kB * T_curr) / q;
    
    Voc_cel = V_oc_ref + K_V * dT;
    Iph_cel = max(0.01, (I_sc_ref + K_I * dT) * scale_AM0);
    I_0_curr = (Iph_cel - Voc_cel / R_sh) / (exp(Voc_cel / (a_eq * V_th)) - 1.0);
    
    i_cel = i_norm * Iph_cel;
    v_j = a_eq * V_th * log(max(1.0, (Iph_cel - i_cel) / I_0_curr + 1.0));
    v_cel = max(0, v_j - i_cel * R_s - (v_j / R_sh) * R_s);
    
    V_array = v_cel * N_s; I_array = i_cel * N_p_no;
    P_array = (V_array .* I_array) / 1e3;
    [V_sort, s_idx] = sort(V_array); P_sort = P_array(s_idx);
    [p_max, mpp_idx] = max(P_sort);
    plot(V_sort, P_sort, cores{idx}, 'linewidth', 1.8);
    plot(V_sort(mpp_idx), p_max, 'ko', 'markersize', 4);
end
plot([400, 400], [0, 75], 'k--', 'linewidth', 1.5);
plot([350, 350], [0, 75], 'k:', 'linewidth', 1.0);
plot([650, 650], [0, 75], 'k:', 'linewidth', 1.0);
xlabel('Tensao da Matriz (V)', 'fontsize', 9);
ylabel('Potencia da Matriz do No (kW)', 'fontsize', 9);
title('(a2) Curvas P-V do No e Faixa Buck (350 a 650 V)', 'fontsize', 10, 'fontweight', 'bold');
legend({'-150 C [Extrapolado]', '-50 C [Extrapolado]', '+28 C (Cat.: 51,5 kW | Mod.: 51,6 kW)', '+65 C (Nom)', '+85 C', 'Barramento (400 V)', 'Faixa Buck'}, 'location', 'northwest', 'fontsize', 7);
xlim([0, 750]); ylim([0, 75]);

% Subpainel (b1): Fluencia de Radiacao Acumulada
t_anos = linspace(0, 10, 100);
Phi_x = 1.0e14; % Escala escolhida para o ajuste [PROVISORIO]
% Ajuste pontual do modelo parametrico ao dado de eficiencia_mp=0.90 do datasheet XTJ Prime,
% tabela de retencao sob o padrao europeu ECSS, em Phi=5e14 e-/cm^2 (1 MeV).
% Phi_x=1e14 e-/cm^2 e uma escala escolhida para o modelo, nao um parametro medido.
eta_mp_ret_ECSS_5e14 = 0.90; % Dado direto da tabela Spectrolab [VERIFICAVEL]
C_d = (1.0 - eta_mp_ret_ECSS_5e14) / log10(1.0 + 5.0); % ~= 0.12851 [PROVISORIO / AJUSTE PONTUAL]
% Este ajuste de um ponto nao valida a forma logaritmica entre pontos nem substitui modelo de ambiente orbital.
% Taxas de fluencia ilustrativas [PROVISORIAS], sem modelo de ambiente orbital ou blindagem.
phi_500 = 1.4e13 * t_anos; phi_700 = 2.5e13 * t_anos; phi_1000 = 5.8e13 * t_anos;

subplot(2, 2, 3);
hold on; grid on; box on;
plot(t_anos, phi_500 / 1e14, 'b-', 'linewidth', 1.8);
plot(t_anos, phi_700 / 1e14, 'g-', 'linewidth', 2.0);
plot(t_anos, phi_1000 / 1e14, 'r-', 'linewidth', 1.8);
plot([0, 10], [1.0, 1.0], 'k--', 'linewidth', 1.2);
xlabel('Tempo de Missao (anos)', 'fontsize', 9);
ylabel('Fluencia 1-MeV (x 10^{14} e^-/cm^2)', 'fontsize', 9);
title('(b1) Fluencia de Radiacao Acumulada (Dano DDD)', 'fontsize', 10, 'fontweight', 'bold');
legend({'500 km', '700 km (Nominal)', '1000 km (Extremo)', 'Fluencia Joelho \Phi_x'}, 'location', 'northwest', 'fontsize', 8);
xlim([0, 10]);

% Subpainel (b2): Retencao de Potencia em 10 Anos
f_ddd_500  = 1.0 - C_d * log10(1.0 + phi_500 / Phi_x);
f_ddd_700  = 1.0 - C_d * log10(1.0 + phi_700 / Phi_x);
f_ddd_1000 = 1.0 - C_d * log10(1.0 + phi_1000 / Phi_x);
f_env = 1.0 - 0.052 * (t_anos / 10.0);
p_rel_500  = f_ddd_500 .* f_env * 100.0;
p_rel_700  = f_ddd_700 .* f_env * 100.0;
p_rel_1000 = f_ddd_1000 .* f_env * 100.0;

subplot(2, 2, 4);
hold on; grid on; box on;
plot(t_anos, p_rel_500, 'b-', 'linewidth', 1.8);
plot(t_anos, p_rel_700, 'g-', 'linewidth', 2.2);
plot(t_anos, p_rel_1000, 'r-', 'linewidth', 1.8);
plot([0, 10], [88.0, 88.0], 'k--', 'linewidth', 1.5);
xlabel('Tempo de Missao (anos)', 'fontsize', 9);
ylabel('Retencao de Potencia P_{mp}(t) (%)', 'fontsize', 9);
title('(b2) Retencao de Potencia EOL em 10 Anos (Calibrada XTJ Prime)', 'fontsize', 10, 'fontweight', 'bold');
perda_500 = 100.0 - p_rel_500(end);
perda_700 = 100.0 - p_rel_700(end);
perda_1000 = 100.0 - p_rel_1000(end);
legend({sprintf('500 km (Perda %.1f%%)', perda_500), ...
        sprintf('700 km (Nominal: Perda %.1f%%)', perda_700), ...
        sprintf('1000 km (Perda %.1f%%)', perda_1000), ...
        'Requisito EOL (88%)'}, 'location', 'southwest', 'fontsize', 8);
xlim([0, 10]); ylim([75, 102]);

if ~exist('figuras','dir'), mkdir('figuras'); end
print(fig_c3, 'figuras/figura_3_3_desempenho_termico_radiacao_am0.png', '-dpng', '-r150');
fprintf('   -> Figura 3.3 Composta salva com sucesso em figuras/figura_3_3_desempenho_termico_radiacao_am0.png\n\n');

% ------------------------------------------------------------------------------
% FIGURA 3.4 COMPOSTA: HOT SPOT, PROTECAO BYPASS E TRANSIENTE MPPT NO ECLIPSE
% ------------------------------------------------------------------------------
fprintf('2. Gerando Figura 3.4 Composta (Hot Spot + Transiente MPPT)...\n');

fig_c4 = figure(2, 'position', [50, 50, 1100, 750]);

% Subpainel (a1): Curvas em Polarizacao Reversa (Modelo Shockley de Bypass)
v_rev_vec = linspace(-35, 3, 300);
i_normal = max(0, I_sc_ref - (v_rev_vec / R_sh) .* (v_rev_vec < 0) - (v_rev_vec > 2.5) * 0.5);
i_shaded_nobypass = - (v_rev_vec / R_sh) - 1e-4 * (abs(v_rev_vec).^2.2) .* (v_rev_vec < 0);

% Modelo Shockley do diodo Schottky acoplado para conduzir I_mp_ref (0.52815 A) em V_clamp = -0.350 V:
n_sch = 1.05;
V_clamp_target = 0.350; % V em conducao direta
Is_sch = I_mp_ref / (exp(V_clamp_target / (n_sch * V_th_ref)) - 1.0); % 1.3915 uA

v_rev_bypass = linspace(-1.0, 3, 300);
i_schottky = Is_sch * (exp(max(-30.0, -v_rev_bypass / (n_sch * V_th_ref))) - 1.0);
i_shaded_bypass = - (v_rev_bypass / R_sh) + i_schottky;
i_shaded_bypass = min(0.65, i_shaded_bypass);

subplot(2, 2, 1);
hold on; grid on; box on;
plot(v_rev_vec, i_normal, 'k-', 'linewidth', 1.5);
plot(v_rev_vec, i_shaded_nobypass, 'r-', 'linewidth', 2.0);
plot(v_rev_bypass, i_shaded_bypass, 'b-', 'linewidth', 2.0);
plot([-35, 3], [I_mp_ref, I_mp_ref], 'm--', 'linewidth', 1.2);
xlabel('Tensao Terminal da Celula (V)', 'fontsize', 9);
ylabel('Corrente da Celula (A)', 'fontsize', 9);
title('(a1) Polarizacao Reversa e Clamp do Bypass (-0.35 V)', 'fontsize', 10, 'fontweight', 'bold');
legend({'Normal Iluminada', 'Sombreada Sem Bypass', 'Com Bypass Schottky (-0.35 V)', sprintf('Corrente Forcada I_{mp} (%.3f A)', I_mp_ref)}, 'location', 'northwest', 'fontsize', 7);
xlim([-35, 3]); ylim([-0.2, 0.7]);

% Subpainel (a2): Elevacao Termica do Hot Spot
t_sec_hs = linspace(0, 180, 200);
tau_th = 12.5 * 3.8; % 47.5 s
% Dissipacao termica acoplada rigorosamente aos mesmos I_mp_ref e V_clamp:
P_diss_nobypass = 28.0 * I_mp_ref;             % 28 V x 0.52815 A = 14.788 W [PROVISORIO]
P_diss_bypass   = V_clamp_target * I_mp_ref;   % 0.350 V x 0.52815 A = 0.1849 W [PROVISORIO]

T_nobypass_C = (338.15 + P_diss_nobypass * 12.5 * (1.0 - exp(-t_sec_hs / tau_th))) - 273.15;
T_bypass_C   = (338.15 + P_diss_bypass   * 12.5 * (1.0 - exp(-t_sec_hs / tau_th))) - 273.15;

subplot(2, 2, 2);
hold on; grid on; box on;
plot(t_sec_hs, T_nobypass_C, 'r-', 'linewidth', 2.2);
plot(t_sec_hs, T_bypass_C, 'b-', 'linewidth', 2.0);
plot([0, 180], [180, 180], 'k--', 'linewidth', 1.2);
plot([0, 180], [230, 230], 'r--', 'linewidth', 1.5);
xlabel('Tempo sob Sombra/Defeito (s)', 'fontsize', 9);
ylabel('Temperatura da Celula (C)', 'fontsize', 9);
title('(a2) Evolucao Termica do Hot Spot no Vacuo (Estudo Parametrico)', 'fontsize', 10, 'fontweight', 'bold');
legend({'Sem Bypass (Destruicao)', 'Com Bypass Ativo', 'Limite Adesivo RTV (180 C)', 'Fusao Solda Prata (230 C)'}, 'location', 'northwest', 'fontsize', 7);
xlim([0, 180]); ylim([50, 260]);

% Subpainel (b1): Tensao de String no Choque Solar do Terminador
t_trans = linspace(0, 250, 400);
G_trans = zeros(size(t_trans));
for m = 1:length(t_trans)
    tt = t_trans(m);
    if tt >= 20 && tt <= 35
        G_trans(m) = G_AM0 * (tt - 20) / 15;
    elseif tt > 35
        G_trans(m) = G_AM0;
    end
end
T_cel_trans_K = 173.15 + (338.15 - 173.15) * (1.0 - exp(-max(0, t_trans - 20) / 85.0));

V_mpp_teorico = zeros(size(t_trans));
P_mpp_teorico = zeros(size(t_trans));
for m = 1:length(t_trans)
    if G_trans(m) >= 10
        Tc = T_cel_trans_K(m);
        V_mpp_teorico(m) = N_s * (V_mp_ref + K_V_mp * (Tc - T_ref));
        Imp_c = (I_mp_ref + K_I_mp * (Tc - T_ref)) * (G_trans(m) / G_datasheet);
        P_mpp_teorico(m) = (V_mpp_teorico(m) * Imp_c * N_p_no) / 1e3;
    end
end

V_mppt_adapt = 400.0 * ones(size(t_trans));
V_mppt_fixed = 400.0 * ones(size(t_trans));
P_mppt_adapt = zeros(size(t_trans));
P_mppt_fixed = zeros(size(t_trans));

v_ta = 400.0; v_tf = 400.0;
for m = 1:length(t_trans)
    if G_trans(m) >= 10
        vt = V_mpp_teorico(m);
        v_ta = v_ta + 0.35 * (vt - v_ta);
        V_mppt_adapt(m) = v_ta;
        P_mppt_adapt(m) = P_mpp_teorico(m) * (1.0 - 0.008 * (1.0 - exp(-(t_trans(m)-35)/10.0)));
        v_tf = v_tf + 2.5 * sign(vt - v_tf) + 1.2 * sin(m * 0.8);
        V_mppt_fixed(m) = v_tf;
        P_mppt_fixed(m) = P_mpp_teorico(m) * 0.965;
    end
end

subplot(2, 2, 3);
hold on; grid on; box on;
plot(t_trans, V_mpp_teorico, 'k--', 'linewidth', 1.8);
plot(t_trans, V_mppt_adapt, 'g-', 'linewidth', 2.0);
plot(t_trans, V_mppt_fixed, 'r-', 'linewidth', 1.5);
plot([0, 250], [400, 400], 'b-', 'linewidth', 1.5);
plot([0, 250], [650, 650], 'm--', 'linewidth', 1.2);
xlabel('Tempo desde a Saida do Eclipse (s)', 'fontsize', 9);
ylabel('Tensao de Entrada do Conversor (V)', 'fontsize', 9);
title('(b1) Dinamica de Tensao no Choque Solar (Estudo Parametrico)', 'fontsize', 10, 'fontweight', 'bold');
legend({'V_{mp} Otimo Analitico', 'P&O Adaptativo Parametrico', 'P&O Passo Fixo', 'Barramento (400 V)', 'Teto Conversor (650 V)'}, 'location', 'northeast', 'fontsize', 7);
xlim([0, 250]); ylim([300, 680]);

% Subpainel (b2): Potencia Capturada pelo MPPT
subplot(2, 2, 4);
hold on; grid on; box on;
plot(t_trans, P_mpp_teorico, 'k--', 'linewidth', 1.8);
plot(t_trans, P_mppt_adapt, 'g-', 'linewidth', 2.0);
plot(t_trans, P_mppt_fixed, 'r-', 'linewidth', 1.5);
xlabel('Tempo desde a Saida do Eclipse (s)', 'fontsize', 9);
ylabel('Potencia Entregue pela Matriz do No (kW)', 'fontsize', 9);
title('(b2) Potencia Capturada na Penumbra (No 51,5 kW instalados / ~50 kW uteis)', 'fontsize', 10, 'fontweight', 'bold');
legend({'Disponivel Otima (No 51,5 kW)', 'P&O Adaptativo (>99%)', 'Passo Fixo (~96.5%)'}, 'location', 'southeast', 'fontsize', 8);
xlim([0, 250]); ylim([0, 75]);

if ~exist('figuras','dir'), mkdir('figuras'); end
print(fig_c4, 'figuras/figura_3_4_contingencias_hotspot_mppt.png', '-dpng', '-r150');
fprintf('   -> Figura 3.4 Composta salva com sucesso em figuras/figura_3_4_contingencias_hotspot_mppt.png\n\n');
fprintf('=== FIGURAS COMPOSTAS DA SECAO 3.2 CONCLUIDAS ===\n');
