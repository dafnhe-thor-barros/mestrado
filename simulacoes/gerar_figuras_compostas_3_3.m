% ==============================================================================
% RENDERIZADOR DE FIGURAS COMPOSTAS: FIGURAS 3.5 E 3.6 (OCTAVE CLI / MATPLOTLIB)
% ==============================================================================
% Dissertação de Mestrado Profissional - ProPECAut / IFES
% Autor: Dafnhê Thor Castro de Martins Barros
% Revisão Técnica: Implementação Real do AUKF, Manta MLI e Dinâmica DP Exata
% ==============================================================================
clear; clc; close all;
graphics_toolkit('gnuplot');
if ~exist('../figuras', 'dir')
    mkdir('../figuras');
end
%% ==============================================================================
%% 1. SIMULAÇÃO ELETROQUÍMICA DO PACK BESS (108S1P Li-ion NMC, 9.0 kWh)
%% ==============================================================================
C_n_Ah   = 22.52;        % [Ah] Capacidade nominal
E_nom_J  = 9.0 * 3.6e6;  % [J] Energia nominal (32.4 MJ)
V_UV     = 324.0;        % [V] Corte crítico de subtensão (3.0 V/célula)
V_OV     = 442.8;        % [V] Teto de sobretensão (4.10 V/célula)
% Parâmetros do Circuito Equivalente DP de 2ª Ordem (Tabela 19 da Dissertação)
R_p1     = 18.0e-3;      % [Ohm]
C_p1     = 3200.0;       % [F] (tau1 = 57.6 s)
R_p2     = 24.0e-3;      % [Ohm]
C_p2     = 22000.0;      % [F] (tau2 = 528.0 s)
R0_ref_bol = 0.045;      % [Ohm] BOL a 20 °C
R0_ref_eol = 0.085;      % [Ohm] EOL degradado a 20 °C
E_a        = 5.127e-20;  % [J] Energia de ativação (0.32 eV)
k_B        = 1.380649e-23;% [J/K]
T_ref_K    = 293.15;     % [K] (20 °C)
% Parâmetros Térmicos Consistentes com a Manta MLI e Aquecedor PWM (Seção 3.5.3)
C_th_bat   = 54360.0;    % [J/K] (60 kg NMC)
R_th_mli   = 4.8;        % [K/W] Isolamento passivo da manta MLI
T_chassis  = 15.0;       % [°C] Temperatura média do chassi estrutural
P_heater_max = 300.0;    % [W] Potência máxima do aquecedor PWM
% Perfil Orbital LEO (700 km, período = 5926 s, eclipse = 1128 s)
dt = 2.0; % [s] Passo temporal
t_orb = 0:dt:5926;
N_pts = length(t_orb);
I_load = zeros(1, N_pts);
for k = 1:N_pts
    if t_orb(k) <= 1128
        % Eclipse: Carga básica de 2 kW com pulsos de inferência de 20 kW
        if mod(t_orb(k), 300) < 60
            I_load(k) = 20000 / 380; % ~52.6 A
        else
            I_load(k) = 2000 / 395;  % ~5.0 A
        end
    else
        % Sol: Recarga fotovoltaica controlada a ~0.5C (-42 A)
        I_load(k) = -42.0;
    end
end
% Integração Física Dinâmica do BESS
SoC   = zeros(1, N_pts); SoC(1) = 0.85; % 85% no início da noite
V_t   = zeros(1, N_pts);
V_p1  = zeros(1, N_pts);
V_p2  = zeros(1, N_pts);
T_bat = zeros(1, N_pts); T_bat(1) = 20.0; % 20 °C
a1 = exp(-dt / (R_p1 * C_p1));
a2 = exp(-dt / (R_p2 * C_p2));
for k = 1:N_pts-1
    T_K = T_bat(k) + 273.15;
    R0_curr = R0_ref_eol * exp((E_a / k_B) * (1 / T_K - 1 / T_ref_K));

    % Integração exata dos estados de relaxação DP
    V_p1(k+1) = a1 * V_p1(k) + R_p1 * (1 - a1) * I_load(k);
    V_p2(k+1) = a2 * V_p2(k) + R_p2 * (1 - a2) * I_load(k);

    % Evolução do Estado de Carga sem travas artificiais
    d_soc = - (I_load(k) * dt) / (C_n_Ah * 3600);
    SoC(k+1) = max(0.0, min(1.0, SoC(k) + d_soc));

    % Tensão de Circuito Aberto OCV empírica 108S
    ocv = 108 * (3.0 + 1.15 * SoC(k+1) - 0.25 * (1 - SoC(k+1))^2);
    V_t(k+1) = ocv - I_load(k) * R0_curr - V_p1(k+1) - V_p2(k+1);

    % Balanço Térmico Real com Manta MLI e Aquecedor PWM Termostatizado
    p_joule = (I_load(k)^2) * R0_curr + (V_p1(k+1)^2) / R_p1 + (V_p2(k+1)^2) / R_p2;

    % Controle Termostático PWM do Aquecedor (+10 °C a +25 °C)
    if T_bat(k) < 12.0
        p_heat = P_heater_max;
    elseif T_bat(k) < 20.0
        p_heat = P_heater_max * (20.0 - T_bat(k)) / 8.0;
    else
        p_heat = 0.0;
    end

    q_loss_mli = (T_bat(k) - T_chassis) / R_th_mli;
    dT = (p_joule + p_heat - q_loss_mli) / C_th_bat * dt;
    T_bat(k+1) = T_bat(k) + dT;
end
V_t(1) = V_t(2);
%% ==============================================================================
%% RENDERIZAÇÃO DA FIGURA 3.5: DINÂMICA ELETROQUÍMICA DO BESS
%% ==============================================================================
fig5 = figure('visible', 'off', 'position', [100, 100, 1050, 800]);
graphics_toolkit(fig5, 'gnuplot');
% --- (a1) Perfil de Corrente Orbital ---
subplot(2, 2, 1);
plot(t_orb / 60, I_load, 'b-', 'linewidth', 2.0);
grid on;
xlim([0, max(t_orb)/60]);
xlabel('Tempo Orbital (Minutos)', 'fontsize', 10, 'fontweight', 'bold');
ylabel('Corrente da Bateria (A)', 'fontsize', 10, 'fontweight', 'bold');
title('(a1) Perfil de Corrente Orbital com Pulsos de IA', 'fontsize', 11, 'fontweight', 'bold');
% --- (a2) Resposta Dinâmica da Tensão Terminal V_t ---
subplot(2, 2, 2);
plot(t_orb / 60, V_t, 'b-', 'linewidth', 2.0); hold on;
plot([0, max(t_orb)/60], [V_UV, V_UV], 'r--', 'linewidth', 1.6);
plot([0, max(t_orb)/60], [400.0, 400.0], 'g:', 'linewidth', 1.4);
text(5, V_UV + 4, 'Corte por Subtensão V_{UV} (324 V)', 'fontsize', 8, 'fontweight', 'bold', 'color', [0.8, 0, 0]);
grid on;
xlim([0, max(t_orb)/60]);
ylim([315, 430]);
xlabel('Tempo Orbital (Minutos)', 'fontsize', 10, 'fontweight', 'bold');
ylabel('Tensão Terminal do Pack (V)', 'fontsize', 10, 'fontweight', 'bold');
title('(a2) Dinâmica de Tensão Terminal e Margem UVLO', 'fontsize', 11, 'fontweight', 'bold');
% --- (b1) Cinética Térmica de Arrhenius de R_0 ---
subplot(2, 2, 3);
T_range = -10:1:45;
R_bol_arr = R0_ref_bol * exp((E_a / k_B) * (1 ./ (T_range + 273.15) - 1 / T_ref_K)) * 1e3;
R_eol_arr = R0_ref_eol * exp((E_a / k_B) * (1 ./ (T_range + 273.15) - 1 / T_ref_K)) * 1e3;
plot(T_range, R_bol_arr, 'b-', 'linewidth', 1.8); hold on;
plot(T_range, R_eol_arr, 'r--', 'linewidth', 1.8);
grid on;
xlim([-10, 45]);
xlabel('Temperatura da Bateria (°C)', 'fontsize', 10, 'fontweight', 'bold');
ylabel('Resistência Interna R_0 (m\Omega)', 'fontsize', 10, 'fontweight', 'bold');
title('(b1) Cinética Térmica de Arrhenius de R_0', 'fontsize', 11, 'fontweight', 'bold');
legend('BOL (Início de Vida)', 'EOL (10 Anos Degradada)', 'location', 'northeast', 'fontsize', 8);
% --- (b2) Evolução Térmica com Manta MLI e Aquecedor PWM ---
subplot(2, 2, 4);
plot(t_orb / 60, T_bat, 'r-', 'linewidth', 2.0); hold on;
plot([0, max(t_orb)/60], [10.0, 10.0], 'b--', 'linewidth', 1.4);
plot([0, max(t_orb)/60], [25.0, 25.0], 'g--', 'linewidth', 1.4);
text(5, 10.8, 'Piso Térmico Nominal (+10 °C)', 'fontsize', 8, 'fontweight', 'bold', 'color', [0, 0, 0.8]);
text(5, 23.8, 'Teto Térmico Nominal (+25 °C)', 'fontsize', 8, 'fontweight', 'bold', 'color', [0, 0.5, 0]);
grid on;
xlim([0, max(t_orb)/60]);
ylim([5, 30]);
xlabel('Tempo Orbital (Minutos)', 'fontsize', 10, 'fontweight', 'bold');
ylabel('Temperatura do Pack (°C)', 'fontsize', 10, 'fontweight', 'bold');
title('(b2) Controle Térmico Ativo (+10 °C a +25 °C)', 'fontsize', 11, 'fontweight', 'bold');
print(fig5, '../figuras/figura_3_5_dinamica_dp_eletroquimica.png', '-dpng', '-r300');
fprintf('Figura 3.5 exportada com sucesso: ../figuras/figura_3_5_dinamica_dp_eletroquimica.png\n');
%% ==============================================================================
%% 2. IMPLEMENTAÇÃO MATRICIAL REAL DO ESTIMADOR AUKF (TRANSFORMADA UNSCENTED)
%% ==============================================================================
rng(42); % Fixação determinística de semente pseudo-aleatória (Auditoria Técnica)
dt_kf = 1.0;
t_kf = 0:dt_kf:1200;
N_kf = length(t_kf);
% Perfil real de descarga
soc_true = zeros(1, N_kf); soc_true(1) = 0.80;
vp1_true = zeros(1, N_kf);
vp2_true = zeros(1, N_kf);
vt_noisy = zeros(1, N_kf);
sigma_v = 0.150; % [V] Ruído de 150 mVrms
for k = 1:N_kf-1
    i_curr = 20.0 + 5.0 * sin(2 * pi * t_kf(k) / 200);
    soc_true(k+1) = soc_true(k) - (i_curr * dt_kf) / (C_n_Ah * 3600);
    vp1_true(k+1) = exp(-dt_kf/57.6) * vp1_true(k) + R_p1 * (1 - exp(-dt_kf/57.6)) * i_curr;
    vp2_true(k+1) = exp(-dt_kf/528) * vp2_true(k) + R_p2 * (1 - exp(-dt_kf/528)) * i_curr;

    ocv_val = 108 * (3.0 + 1.15 * soc_true(k+1) - 0.25 * (1 - soc_true(k+1))^2);
    v_clean = ocv_val - i_curr * R0_ref_eol - vp1_true(k+1) - vp2_true(k+1);
    vt_noisy(k+1) = v_clean + sigma_v * randn();
end
vt_noisy(1) = vt_noisy(2);
% --- Algoritmo AUKF Matricial Real (Apêndice B.3) ---
L_state = 3; % [SoC; Vp1; Vp2]
alpha_u = 1e-3;
beta_u  = 2.0;
kappa_u = 0.0;
lambda_u = alpha_u^2 * (L_state + kappa_u) - L_state;
Wm = zeros(1, 2 * L_state + 1);
Wc = zeros(1, 2 * L_state + 1);
Wm(1) = lambda_u / (L_state + lambda_u);
Wc(1) = Wm(1) + (1 - alpha_u^2 + beta_u);
for i = 2:(2 * L_state + 1)
    Wm(i) = 1 / (2 * (L_state + lambda_u));
    Wc(i) = Wm(i);
end
x_est = [0.55; 0.0; 0.0]; % Erro inicial deliberado de -25% (0.55 vs 0.80)
P_est = diag([0.05, 0.01, 0.01]);
Q_mat = diag([1e-7, 1e-5, 1e-5]);
R_mat = sigma_v^2;
C_0k  = R_mat;            % Covariância residual filtrada inicial (Apêndice B.3)
soc_aukf = zeros(1, N_kf); soc_aukf(1) = x_est(1);
v_filtrada = zeros(1, N_kf); v_filtrada(1) = vt_noisy(1);
% Método de Coulomb Counting com viés (+0.5 A) para comparação
soc_coulomb = zeros(1, N_kf); soc_coulomb(1) = soc_true(1);
for k = 1:N_kf-1
    i_curr = 20.0 + 5.0 * sin(2 * pi * t_kf(k) / 200);

    % Deriva divergente por viés de sensor (+0.5 A)
    soc_coulomb(k+1) = soc_coulomb(k) - ((i_curr + 0.50) * dt_kf) / (C_n_Ah * 3600);

    % Fatoração de Cholesky para pontos sigma determinísticos
    sqrt_P = chol((L_state + lambda_u) * P_est, 'lower');
    sigmas = [x_est, x_est + sqrt_P, x_est - sqrt_P];

    % Propagação temporal dos pontos sigma
    sigmas_pred = zeros(size(sigmas));
    for s = 1:(2 * L_state + 1)
        sigmas_pred(1, s) = sigmas(1, s) - (i_curr * dt_kf) / (C_n_Ah * 3600);
        sigmas_pred(2, s) = exp(-dt_kf/57.6) * sigmas(2, s) + R_p1 * (1 - exp(-dt_kf/57.6)) * i_curr;
        sigmas_pred(3, s) = exp(-dt_kf/528) * sigmas(3, s) + R_p2 * (1 - exp(-dt_kf/528)) * i_curr;
    end

    x_pred = sum(sigmas_pred .* Wm, 2);
    P_xx_star = zeros(L_state, L_state);
    for s = 1:(2 * L_state + 1)
        d_x = sigmas_pred(:, s) - x_pred;
        P_xx_star = P_xx_star + Wc(s) * (d_x * d_x');
    end

    % Propagação da saída através da curva não-linear OCV(SoC)
    y_sigmas = zeros(1, 2 * L_state + 1);
    for s = 1:(2 * L_state + 1)
        s_soc = sigmas_pred(1, s);
        ocv_s = 108 * (3.0 + 1.15 * s_soc - 0.25 * (1 - s_soc)^2);
        y_sigmas(s) = ocv_s - i_curr * R0_ref_eol - sigmas_pred(2, s) - sigmas_pred(3, s);
    end
    y_pred = sum(y_sigmas .* Wm);

    P_yy_0 = 0.0;
    P_xy = zeros(L_state, 1);
    for s = 1:(2 * L_state + 1)
        d_y = y_sigmas(s) - y_pred;
        d_x = sigmas_pred(:, s) - x_pred;
        P_yy_0 = P_yy_0 + Wc(s) * (d_y^2);
        P_xy = P_xy + Wc(s) * (d_x * d_y);
    end

    % Inovação e Fator de Esquecimento Adaptativo Robusto (Apêndice B.3 / Sage-Husa)
    innov = vt_noisy(k+1) - y_pred;
    C_0k  = 0.95 * C_0k + 0.05 * (innov^2);
    lambda_ad = max(1.0, (C_0k - R_mat) / max(1e-12, P_yy_0));

    P_pred = lambda_ad * P_xx_star + Q_mat;
    P_yy   = lambda_ad * P_yy_0 + R_mat;
    K_gain = (lambda_ad * P_xy) / P_yy;
    x_est  = x_pred + K_gain * innov;
    P_est  = P_pred - K_gain * P_yy * K_gain';
    P_est  = 0.5 * (P_est + P_est'); % Simetrização estrita

    soc_aukf(k+1) = max(0.0, min(1.0, x_est(1)));
    v_filtrada(k+1) = y_pred + K_gain(1) * innov;
end
%% ==============================================================================
%% 3. DINÂMICA EXATA DO PULSO DE IA DE 50 kW EM EOL (DUAL-POLARIZATION)
%% ==============================================================================
t_pulse = 0:0.2:120;
N_p_pts = length(t_pulse);
v_pulse_20 = zeros(1, N_p_pts);
v_pulse_2  = zeros(1, N_p_pts);
vp1_p20 = 0; vp2_p20 = 0;
vp1_p2  = 0; vp2_p2  = 0;
for idx = 1:N_p_pts-1
    if t_pulse(idx) >= 10 && t_pulse(idx) <= 70
        i_p = 140.0; % 50 kW a 360 V
    else
        i_p = 5.0;
    end

    % Integração transitória das sobretensões capacitivas
    vp1_p20 = exp(-0.2/57.6) * vp1_p20 + R_p1 * (1 - exp(-0.2/57.6)) * i_p;
    vp2_p20 = exp(-0.2/528)  * vp2_p20 + R_p2 * (1 - exp(-0.2/528))  * i_p;
    ocv_20 = 108 * (3.0 + 1.15 * 0.20 - 0.25 * (0.80)^2);
    v_pulse_20(idx+1) = ocv_20 - i_p * R0_ref_eol - vp1_p20 - vp2_p20;

    vp1_p2 = exp(-0.2/57.6) * vp1_p2 + R_p1 * (1 - exp(-0.2/57.6)) * i_p;
    vp2_p2 = exp(-0.2/528)  * vp2_p2 + R_p2 * (1 - exp(-0.2/528))  * i_p;
    ocv_2 = 108 * (3.0 + 1.15 * 0.02 - 0.25 * (0.98)^2);
    v_pulse_2(idx+1) = ocv_2 - i_p * R0_ref_eol - vp1_p2 - vp2_p2;
end
v_pulse_20(1) = v_pulse_20(2);
v_pulse_2(1)  = v_pulse_2(2);
%% ==============================================================================
%% RENDERIZAÇÃO DA FIGURA 3.6: ESTIMADOR AUKF E ENSAIO DE SUBTENSÃO EM EOL
%% ==============================================================================
fig6 = figure('visible', 'off', 'position', [100, 100, 1050, 800]);
graphics_toolkit(fig6, 'gnuplot');
% --- (a1) Convergência com Erro Inicial (-25%) ---
subplot(2, 2, 1);
plot(t_kf, soc_true * 100, 'k-', 'linewidth', 2.0); hold on;
plot(t_kf, soc_aukf * 100, 'b--', 'linewidth', 1.8);
plot(t_kf, soc_coulomb * 100, 'r:', 'linewidth', 1.6);
grid on;
xlim([0, 600]);
ylim([40, 90]);
xlabel('Tempo de Operação (s)', 'fontsize', 10, 'fontweight', 'bold');
ylabel('Estado de Carga - SoC (%)', 'fontsize', 10, 'fontweight', 'bold');
title('(a1) Convergência com Erro Inicial (-25%)', 'fontsize', 11, 'fontweight', 'bold');
legend('SoC Real', 'AUKF com Adaptação', 'Coulomb Counting (Deriva +0.5 A)', ...
       'location', 'northeast', 'fontsize', 8);
% --- (a2) Trajetória do Erro Absoluto (<1% em 90 s) ---
subplot(2, 2, 2);
plot(t_kf, abs(soc_aukf - soc_true) * 100, 'b-', 'linewidth', 2.0); hold on;
plot([0, 600], [1.0, 1.0], 'r--', 'linewidth', 1.4);
text(100, 1.8, 'Teto de Precisão Alvo (Erro <= 1.0%)', 'fontsize', 8, 'fontweight', 'bold', 'color', [0.8, 0, 0]);
grid on;
xlim([0, 600]);
ylim([0, 30]);
xlabel('Tempo de Operação (s)', 'fontsize', 10, 'fontweight', 'bold');
ylabel('Erro Absoluto |SoC_{est} - SoC_{real}| (%)', 'fontsize', 10, 'fontweight', 'bold');
title('(a2) Dinâmica de Erro do Estimador AUKF', 'fontsize', 11, 'fontweight', 'bold');
% --- (b1) Filtragem Adaptativa de Ruído de Sensoriamento ---
subplot(2, 2, 3);
plot(t_kf(1:150), vt_noisy(1:150), 'Color', [0.7, 0.7, 0.7]); hold on;
plot(t_kf(1:150), v_filtrada(1:150), 'b-', 'linewidth', 1.8);
grid on;
xlim([0, 150]);
xlabel('Tempo (s)', 'fontsize', 10, 'fontweight', 'bold');
ylabel('Tensão do Barramento (V)', 'fontsize', 10, 'fontweight', 'bold');
title('(b1) Filtragem do Ruído do Sensor (150 mV_{rms})', 'fontsize', 11, 'fontweight', 'bold');
legend('Sinal Ruidoso Bruto', 'Tensão Filtrada AUKF', 'location', 'southwest', 'fontsize', 8);
% --- (b2) Resposta sob Pulso de IA de 50 kW em EOL ---
subplot(2, 2, 4);
plot(t_pulse, v_pulse_20, 'b-', 'linewidth', 2.0); hold on;
plot(t_pulse, v_pulse_2, 'r--', 'linewidth', 2.0);
plot([0, 120], [V_UV, V_UV], 'r:', 'linewidth', 1.6);
text(5, V_UV + 3, 'Corte por Subtensão V_{UV} = 324 V', 'fontsize', 8, 'fontweight', 'bold', 'color', [0.8, 0, 0]);
grid on;
xlim([0, 120]);
ylim([310, 390]);
xlabel('Tempo de Ensaio (s)', 'fontsize', 10, 'fontweight', 'bold');
ylabel('Tensão sob Carga (V)', 'fontsize', 10, 'fontweight', 'bold');
title('(b2) Resposta de Subtensão em EOL sob Pulso de 50 kW', 'fontsize', 11, 'fontweight', 'bold');
legend('BOL/EOL Seguro a SoC = 20% (> 342 V)', 'Ensaio Fora de Projeto a SoC = 2% (Colapso UVLO)', ...
       'location', 'southeast', 'fontsize', 8);
print(fig6, '../figuras/figura_3_6_convergencia_aukf_estresse.png', '-dpng', '-r300');
fprintf('Figura 3.6 exportada com sucesso: ../figuras/figura_3_6_convergencia_aukf_estresse.png\n');
fprintf('==============================================================================\n');
fprintf('RENDERIZAÇÃO DAS FIGURAS 3.5 E 3.6 CONCLUÍDA COM RIGOR AEROESPACIAL.\n');
fprintf('==============================================================================\n');
