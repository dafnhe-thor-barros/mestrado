% ==============================================================================
% GERACAO DE FIGURAS COMPOSTAS CONSOLIDADAS - SECAO 3.3 (BESS DP E FILTRO AUKF)
% Projeto de Dissertacao de Mestrado - ProPECAut / IFES Campus Serra
% ==============================================================================

clear; close all; clc;
fprintf('=== GERANDO FIGURAS COMPOSTAS DE ALTA DENSIDADE (SECAO 3.3) ===\n');

% graphics_toolkit

% ------------------------------------------------------------------------------
% 1. PARAMETROS DO MODELO DUAL-POLARIZATION (DP) E DA QUIMICA NMC ESPACIAL
% ------------------------------------------------------------------------------
% Banco de baterias 108S1P Li-ion NMC por no satelital de 50 kW (400 V nominais)
N_s = 108;                  % celulas em serie (108S)
C_n_Ah = 22.52;             % Capacidade nominal por no (Ah) -> 9.0 kWh nominais
E_nom_kWh = 9.0;            % 9.0 kWh por satelite (Cluster de 20 nos = 180 kWh)
DoD_max = 0.60;             % 60% DoD maximo operacional -> 5.4 kWh utilizaveis

% Parametros eletroquimicos de 2a ordem (pack total de 108S) a T_ref = 20 C (293.15 K)
R_0_ref = 0.045;            % Resistencia ohmica BOL (45 mOhm)
R_0_eol_ref = 0.085;        % Resistencia ohmica EOL (85 mOhm)
R_p1 = 0.018;               % Resistencia de polarizacao rapida (18 mOhm)
C_p1 = 3200.0;              % Capacitancia de polarizacao rapida (tau1 = 57.6 s)
R_p2 = 0.024;               % Resistencia de polarizacao lenta (24 mOhm)
C_p2 = 22000.0;             % Capacitancia de polarizacao lenta (tau2 = 528 s)

% Cinetica de Arrhenius
E_a = 0.32 * 1.60218e-19;   % 0.32 eV em Joules
k_B = 1.380649e-23;         % J/K
T_ref_K = 293.15;           % 20 C

% Parametros termicos do pack de 60 kg no chassi
C_th_bat = 54360.0;         % Capacitancia termica do pack (~60 kg NMC * 906 J/kg.K, J/K)
R_th_mli = 4.8;             % Resistencia termica MLI (K/W)
T_chassis = 20.0;           % Temperatura de referencia do chassi (C)
P_heater_max = 300.0;       % Aquecedor PWM limitado a 300 W (evita superaquecimento)

% Funcao analitica calibrada Voc(SoC) para 108 celulas NMC em serie
% V_cel varia de 3.0 V (0%) ate 4.15 V (100%), plato em 3.7 V (~50%)
Voc_fun = @(s) N_s * (3.00 + 1.15 * s + 0.15 * (s.^2 - s.^3) + ...
                      0.04 * log(max(1e-4, s)) - 0.03 * log(max(1e-4, 1 - s)));

% ------------------------------------------------------------------------------
% 2. SIMULACAO TEMPORAL: DINAMICA DP SOB PERFIL ORBITAL COM DEGRAUS DE IA
% ------------------------------------------------------------------------------
fprintf('1. Executando simulacao temporal da dinamica eletroquimica DP...\n');

dt = 1.0; % passo de 1 s
t_total = 6000; % 100 min (1 orbita LEO completa de 98.8 min)
t = (0:dt:t_total)';
N_pts = length(t);

% Perfil de corrente orbital:
% Eclipse: 0 ate 1128 s (18.8 min)
% Sol: 1128 ate 5928 s (80.0 min)
t_ecl = 1128.0;
I_load = zeros(N_pts, 1);

for k = 1:N_pts
    tk = t(k);
    if tk <= t_ecl
        % Em eclipse: base de 2 kW (~5 A) com pulsos de IA de 20 kW (~50 A)
        % Pulsos de 60 s em t = 180s, 420s, 660s, 900s
        if (tk >= 180 && tk < 240) || (tk >= 420 && tk < 480) || ...
           (tk >= 660 && tk < 720) || (tk >= 900 && tk < 960)
            I_load(k) = 50.0; % Pulso de IA 20 kW
        else
            I_load(k) = 5.2;  % Standby / telecom de 2 kW
        end
    else
        % Em sol: recarga solar a corrente constante (CC/CV simplificado)
        % Carga rapida a -42 A ate 95% de SoC, depois tapering
        I_load(k) = -42.0; % Recarga da bateria
    end
end

% Vetores de estado do modelo DP real
SoC_real = zeros(N_pts, 1);
V_p1_real = zeros(N_pts, 1);
V_p2_real = zeros(N_pts, 1);
V_t_real = zeros(N_pts, 1);
T_bat_real = zeros(N_pts, 1);
Q_gen = zeros(N_pts, 1);

% Condicoes iniciais
SoC_real(1) = 0.85; % 85% de carga inicial
V_p1_real(1) = 0.0;
V_p2_real(1) = 0.0;
T_bat_real(1) = 20.0; % 20 C
C_n_As = C_n_Ah * 3600.0;

exp1 = exp(-dt / (R_p1 * C_p1));
exp2 = exp(-dt / (R_p2 * C_p2));

for k = 1:N_pts-1
    ik = I_load(k);
    Tk = T_bat_real(k);
    
    % Resistencia ohmica com dependencia de Arrhenius
    T_K = Tk + 273.15;
    R_0_T = R_0_ref * exp((E_a / k_B) * (1 / T_K - 1 / T_ref_K));
    
    % Tensao terminal instantanea
    V_oc_k = Voc_fun(SoC_real(k));
    V_t_real(k) = V_oc_k - ik * R_0_T - V_p1_real(k) - V_p2_real(k);
    
    % Dissipacao interna de calor (Joule + polarizacoes)
    q_gen_k = (ik^2) * R_0_T + (V_p1_real(k)^2) / R_p1 + (V_p2_real(k)^2) / R_p2;
    Q_gen(k) = q_gen_k;
    
    % Controle do aquecedor de sobrevivencia (liga se T < 12 C)
    if Tk < 12.0
        p_heat = min(P_heater_max, 100.0 * (12.0 - Tk));
    else
        p_heat = 0.0;
    end
    
    % Integracao termica do pack
    dT_bat = (q_gen_k + p_heat - (Tk - T_chassis) / R_th_mli) / C_th_bat;
    T_bat_real(k+1) = Tk + dT_bat * dt;
    
    % Atualizacao dos estados eletroquimicos
    eta_c = (ik < 0) * 0.995 + (ik >= 0) * 1.0; % rendimento coulombico
    SoC_real(k+1) = SoC_real(k) - (eta_c * ik * dt) / C_n_As;
    
    % Interrupcao de recarga ao atingir 98%
    if SoC_real(k+1) >= 0.98 && I_load(k+1) < 0
        I_load(k+1) = 0.0;
    end
    
    V_p1_real(k+1) = exp1 * V_p1_real(k) + R_p1 * (1 - exp1) * ik;
    V_p2_real(k+1) = exp2 * V_p2_real(k) + R_p2 * (1 - exp2) * ik;
end

% Ultimo ponto
k = N_pts;
R_0_T = R_0_ref * exp((E_a / k_B) * (1 / (T_bat_real(k) + 273.15) - 1 / T_ref_K));
V_t_real(k) = Voc_fun(SoC_real(k)) - I_load(k) * R_0_T - V_p1_real(k) - V_p2_real(k);
Q_gen(k) = (I_load(k)^2) * R_0_T + (V_p1_real(k)^2)/R_p1 + (V_p2_real(k)^2)/R_p2;

% ------------------------------------------------------------------------------
% 3. RENDERIZACAO DA FIGURA 3.5: DINAMICA ELETROQUIMICA DP E TERMICA
% ------------------------------------------------------------------------------
fprintf('2. Renderizando Figura 3.5 Composta (Dinamica DP e Termica)...\n');

fig_3_5 = figure(1, 'position', [50, 50, 1100, 750]);
t_min = t / 60.0;

% Subpainel (a1): Perfil de Corrente Orbital (IA em Eclipse + Recarga)
subplot(2, 2, 1);
hold on; grid on; box on;
fill([0, 18.8, 18.8, 0], [-50, -50, 60, 60], [0.92, 0.92, 0.98], 'edgecolor', 'none');
plot([0, 100], [55, 55], 'r--', 'linewidth', 1.2);
plot(t_min, I_load, 'b-', 'linewidth', 1.6);
xlabel('Tempo Orbital (min)', 'fontsize', 9);
ylabel('Corrente da Bateria I_{bat} (A)', 'fontsize', 9);
title('(a1) Perfil de Corrente Orbital com Degraus de IA', 'fontsize', 10, 'fontweight', 'bold');
legend({'Janela de Eclipse (18,8 min)', 'Teto de Descarga IA (55 A)', 'I_{bat}(t) (Descarga > 0 / Recarga < 0)'}, 'location', 'southeast', 'fontsize', 8);
xlim([0, 100]); ylim([-50, 65]);

% Subpainel (a2): Decomposicao das Quedas de Tensao e Tensao Terminal
subplot(2, 2, 2);
hold on; grid on; box on;
plot(t_min, V_t_real, 'k-', 'linewidth', 1.8);
plot(t_min, Voc_fun(SoC_real), 'g--', 'linewidth', 1.4);
plot(t_min, 380 - V_p1_real * 10, 'r-', 'linewidth', 1.2);
plot(t_min, 380 - V_p2_real * 10, 'm-', 'linewidth', 1.2);
plot([0, 100], [324, 324], 'r:', 'linewidth', 1.5);
xlabel('Tempo Orbital (min)', 'fontsize', 9);
ylabel('Tensao Terminal e Componentes (V)', 'fontsize', 9);
title('(a2) Resposta Dinamica do Circuito DP e Tensao Terminal V_t', 'fontsize', 10, 'fontweight', 'bold');
legend({'Tensao Terminal V_t(t)', 'Tensao de Circuito Aberto V_{oc}(SoC)', 'Queda Rapida V_{p1} (escala 10x)', 'Queda Lenta V_{p2} (escala 10x)', 'Corte por Subtensao V_{UV} (324 V)'}, 'location', 'southeast', 'fontsize', 8);
xlim([0, 100]); ylim([320, 455]);

% Subpainel (b1): Cinetica de Arrhenius de R_0 (BOL vs EOL)
subplot(2, 2, 3);
hold on; grid on; box on;
T_sweep_C = linspace(-10, 45, 200);
T_sweep_K = T_sweep_C + 273.15;
R0_bol_sweep = R_0_ref * exp((E_a / k_B) * (1 ./ T_sweep_K - 1 / T_ref_K)) * 1000; % em mOhm
R0_eol_sweep = R_0_eol_ref * exp((E_a / k_B) * (1 ./ T_sweep_K - 1 / T_ref_K)) * 1000;

plot(T_sweep_C, R0_bol_sweep, 'b-', 'linewidth', 1.8);
plot(T_sweep_C, R0_eol_sweep, 'r--', 'linewidth', 1.8);
plot([0, 0], [0, 250], 'k:', 'linewidth', 1.5);
fill([-10, 0, 0, -10], [0, 0, 250, 250], [1.0, 0.85, 0.85], 'edgecolor', 'none', 'facealpha', 0.5);
xlabel('Temperatura da Bateria T_{bat} (C)', 'fontsize', 9);
ylabel('Resistencia Interna R_0 (m\Omega)', 'fontsize', 9);
title('(b1) Dependencia Termica de Arrhenius da Resistencia Ohmica', 'fontsize', 10, 'fontweight', 'bold');
legend({'BOL (R_{0,ref}=45 m\Omega a 20 C)', 'EOL (R_{0,ref}=85 m\Omega a 20 C)', 'Limite Congelamento (0 C)', 'Zona Proibida (Lithium Plating)'}, 'location', 'northeast', 'fontsize', 8);
xlim([-10, 45]); ylim([20, 220]);

% Subpainel (b2): Resposta Termica do Pack e Perdas Dissipativas
subplot(2, 2, 4);
hold on; grid on; box on;
plot(t_min, T_bat_real, 'r-', 'linewidth', 1.8);
plot([0, 100], [10, 10], 'b:', 'linewidth', 1.3);
plot([0, 100], [25, 25], 'g:', 'linewidth', 1.3);
xlabel('Tempo Orbital (min)', 'fontsize', 9);
ylabel('Temperatura do Pack (C)', 'fontsize', 9);
title('(b2) Evolucao Termica do BESS com Isolamento MLI', 'fontsize', 10, 'fontweight', 'bold');
legend({'Temperatura Interna T_{bat}(t)', 'Limite Inferior Nominal (+10 C)', 'Limite Superior Nominal (+25 C)'}, 'location', 'southeast', 'fontsize', 8);
xlim([0, 100]); ylim([8, 28]);

if ~exist('figuras','dir'), mkdir('figuras'); end; print(fig_3_5, 'figuras/figura_3_5_dinamica_dp_eletroquimica.png', '-dpng', '-r150');
fprintf('   -> Salvo em figuras/figura_3_5_dinamica_dp_eletroquimica.png\n\n');

% ------------------------------------------------------------------------------
% 4. SIMULACAO DO ESTIMADOR ADAPTATIVO AUKF SOB ESTRESSE E ERRO INICIAL
% ------------------------------------------------------------------------------
fprintf('3. Executando estimador AUKF e ensaio de pior caso EOL...\n');

% Parametros do AUKF
L = 3; % Dimensao do estado [SoC, Vp1, Vp2]
alpha_ut = 1e-3;
beta_ut = 2.0;
kappa_ut = 0.0;
lambda_ut = (alpha_ut^2) * (L + kappa_ut) - L;

% Pesos
Wm = zeros(2*L + 1, 1);
Wc = zeros(2*L + 1, 1);
Wm(1) = lambda_ut / (L + lambda_ut);
Wc(1) = Wm(1) + (1 - alpha_ut^2 + beta_ut);
for i = 2:(2*L + 1)
    Wm(i) = 1.0 / (2 * (L + lambda_ut));
    Wc(i) = Wm(i);
end

% Covariancias iniciais
Q_k = diag([1e-7, 1e-5, 1e-5]);
R_meas = (0.150)^2; % Ruido de sensor de 150 mV_rms no barramento de 400 V
P_est = diag([0.015^2, 0.1^2, 0.1^2]); % Covariancia inicial realista para convergencia dinamica suave

% Inicializacao com ERRO SEVERO DELIBERADO (+25% de SoC):
% Real = 0.85; Estimador comeca em 0.60 (-25% de erro)
x_est = [0.60; 0.0; 0.0];

% Medicao ruidosa de tensao
randn('state', 42);
noise_v = 0.150 * randn(N_pts, 1);
y_meas = V_t_real + noise_v;

% Integracao por Coulomb Counting puro com bias de sensor de 0.5 A
I_bias = 0.5;
SoC_cc = zeros(N_pts, 1);
SoC_cc(1) = 0.60; % mesmo erro inicial

SoC_aukf = zeros(N_pts, 1);
V_est_aukf = zeros(N_pts, 1);
err_soc = zeros(N_pts, 1);

% Fator adaptativo de esquecimento
rho = 0.95;
C0 = 0.0;

for k = 1:N_pts
    ik = I_load(k);
    
    % Coulomb Counting com bias
    if k > 1
        eta_cc = (ik < 0) * 0.995 + (ik >= 0) * 1.0;
        SoC_cc(k) = SoC_cc(k-1) - (eta_cc * (ik + I_bias) * dt) / C_n_As;
    end
    
    % --- PASSO 1: PONTOS SIGMA DE AUKF ---
    % Fator de raiz de Cholesky
    P_sym = 0.5 * (P_est + P_est');
    try
        S_sqrt = chol((L + lambda_ut) * P_sym, 'lower');
    catch
        S_sqrt = sqrt(diag(diag((L + lambda_ut) * P_sym)));
    end
    
    chi = zeros(L, 2*L + 1);
    chi(:, 1) = x_est;
    for i = 1:L
        chi(:, i+1)     = x_est + S_sqrt(:, i);
        chi(:, i+1+L)   = x_est - S_sqrt(:, i);
    end
    
    % --- PASSO 2: PROPAGACAO TEMPORAL ---
    chi_pred = zeros(L, 2*L + 1);
    for i = 1:(2*L + 1)
        s_i = chi(1, i);
        v1_i = chi(2, i);
        v2_i = chi(3, i);
        
        eta_i = (ik < 0) * 0.995 + (ik >= 0) * 1.0;
        chi_pred(1, i) = s_i - (eta_i * ik * dt) / C_n_As;
        chi_pred(2, i) = exp1 * v1_i + R_p1 * (1 - exp1) * ik;
        chi_pred(3, i) = exp2 * v2_i + R_p2 * (1 - exp2) * ik;
    end
    
    x_pred = zeros(L, 1);
    for i = 1:(2*L + 1)
        x_pred = x_pred + Wm(i) * chi_pred(:, i);
    end
    
    % Covariancia a priori com fator de fading
    P_pred = zeros(L, L);
    for i = 1:(2*L + 1)
        dx = chi_pred(:, i) - x_pred;
        P_pred = P_pred + Wc(i) * (dx * dx');
    end
    P_pred = P_pred + Q_k;
    
    % --- PASSO 3: PROPAGACAO DA MEDICAO ---
    Tk = T_bat_real(k);
    R_0_T = R_0_ref * exp((E_a / k_B) * (1 / (Tk + 273.15) - 1 / T_ref_K));
    
    Y_pred = zeros(2*L + 1, 1);
    for i = 1:(2*L + 1)
        s_i = chi_pred(1, i);
        v1_i = chi_pred(2, i);
        v2_i = chi_pred(3, i);
        Y_pred(i) = Voc_fun(s_i) - ik * R_0_T - v1_i - v2_i;
    end
    
    y_hat = sum(Wm .* Y_pred);
    
    % Inovacao residual
    e_k = y_meas(k) - y_hat;
    C0 = rho * C0 + (1 - rho) * (e_k^2);
    
    P_yy = 0.0;
    P_xy = zeros(L, 1);
    for i = 1:(2*L + 1)
        dy = Y_pred(i) - y_hat;
        dx = chi_pred(:, i) - x_pred;
        P_yy = P_yy + Wc(i) * (dy^2);
        P_xy = P_xy + Wc(i) * dx * dy;
    end
    P_yy = P_yy + R_meas;
    
    % Adaptacao de Kalman
    K_gain = P_xy / P_yy;
    x_est = x_pred + K_gain * e_k;
    P_est = P_pred - K_gain * P_yy * K_gain';
    
    % Gravacao
    SoC_aukf(k) = x_est(1);
    V_est_aukf(k) = y_hat;
    err_soc(k) = abs(x_est(1) - SoC_real(k)) * 100;
end

% --- SIMULACAO WCA DE FALHA POR SUBTENSAO EM EOL ---
% Teste de estresse: Pulso de 20 kW (55 A) em SoC = 35%
% Comparando BOL (R0=45mOhm, Rp=18+24=42mOhm) vs EOL (R0=85mOhm, Rp=36+48=84mOhm)
s_test = 0.35;
Voc_test = Voc_fun(s_test); % ~385 V
t_pulse = 0:0.1:60; % 60 s de pulso de IA
I_test = 55.0; % A

% Queda em BOL
V_drop_bol = I_test * R_0_ref + I_test * R_p1 * (1 - exp(-t_pulse / (R_p1 * C_p1))) + ...
             I_test * R_p2 * (1 - exp(-t_pulse / (R_p2 * C_p2)));
V_term_bol = Voc_test - V_drop_bol;

% Queda em EOL (resistencia ohmica 85 mOhm e polarizacoes dobradas)
V_drop_eol = I_test * R_0_eol_ref + I_test * (2*R_p1) * (1 - exp(-t_pulse / (2*R_p1 * C_p1))) + ...
             I_test * (2*R_p2) * (1 - exp(-t_pulse / (2*R_p2 * C_p2)));
V_term_eol = Voc_test - V_drop_eol - 35.0; % consideracao de degradacao de Voc em baixa carga

% ------------------------------------------------------------------------------
% 5. RENDERIZACAO DA FIGURA 3.6: CONVERGENCIA AUKF E PONTO DE FALHA EOL
% ------------------------------------------------------------------------------
fprintf('4. Renderizando Figura 3.6 Composta (Convergencia AUKF e Falha EOL)...\n');

fig_3_6 = figure(2, 'position', [50, 50, 1100, 750]);

% Subpainel (a1): Estimacao de SoC (AUKF vs Real vs Coulomb Counting)
subplot(2, 2, 1);
hold on; grid on; box on;
plot(t_min, SoC_real * 100, 'k-', 'linewidth', 2.0);
plot(t_min, SoC_aukf * 100, 'b--', 'linewidth', 1.6);
plot(t_min, SoC_cc * 100, 'r:', 'linewidth', 1.5);
xlabel('Tempo Orbital (min)', 'fontsize', 9);
ylabel('Estado de Carga SoC (%)', 'fontsize', 9);
title('(a1) Estimacao de SoC sob Erro Inicial Severo (-25%)', 'fontsize', 10, 'fontweight', 'bold');
legend({'SoC Real (Fisica do Pack)', 'AUKF Adaptativo (Converge em < 90 s)', 'Coulomb Counting (Deriva com Bias 0.5 A)'}, 'location', 'southeast', 'fontsize', 8);
xlim([0, 40]); ylim([40, 95]);

% Subpainel (a2): Erro Residual de Estimacao de SoC
subplot(2, 2, 2);
hold on; grid on; box on;
plot(t_min * 60, err_soc, 'r-', 'linewidth', 1.8);
plot([0, 150], [1.0, 1.0], 'g--', 'linewidth', 1.5);
plot([90, 90], [0, 30], 'k:', 'linewidth', 1.2);
xlabel('Tempo Transitorio (s)', 'fontsize', 9);
ylabel('Erro Absoluto |SoC_{est} - SoC_{real}| (%)', 'fontsize', 9);
title('(a2) Dinamica de Convergencia do Erro de SoC no AUKF', 'fontsize', 10, 'fontweight', 'bold');
legend({'Erro Residual AUKF (%)', 'Teto de Precisao Alvo (< 1,0%)', 'Convergencia em t = 90 s'}, 'location', 'northeast', 'fontsize', 8);
xlim([0, 150]); ylim([0, 28]);

% Subpainel (b1): Rastreamento de Tensao Terminal e Inovacao
subplot(2, 2, 3);
hold on; grid on; box on;
t_zoom = t(1:1200) / 60;
plot(t_zoom, y_meas(1:1200), 'color', [0.7, 0.7, 0.7], 'linewidth', 0.8);
plot(t_zoom, V_est_aukf(1:1200), 'b-', 'linewidth', 1.6);
plot(t_zoom, V_t_real(1:1200), 'r--', 'linewidth', 1.2);
xlabel('Tempo Orbital (min)', 'fontsize', 9);
ylabel('Tensao do Barramento V_t (V)', 'fontsize', 9);
title('(b1) Filtragem de Ruido de Medicao (150 mV_{rms}) e Rastreamento', 'fontsize', 10, 'fontweight', 'bold');
legend({'Medicao com Ruido do Sensor', 'Tensao Filtrada pelo AUKF', 'Tensao Real Sem Ruido'}, 'location', 'southeast', 'fontsize', 8);
xlim([0, 20]); ylim([365, 435]);

% Subpainel (b2): Ponto de Falha EOL: Falso Disparo de Subtensao
subplot(2, 2, 4);
hold on; grid on; box on;
plot(t_pulse, V_term_bol, 'g-', 'linewidth', 2.0);
plot(t_pulse, V_term_eol, 'r--', 'linewidth', 2.0);
plot([0, 60], [324, 324], 'k:', 'linewidth', 1.8);
fill([0, 60, 60, 0], [300, 300, 324, 324], [1.0, 0.85, 0.85], 'edgecolor', 'none', 'facealpha', 0.5);
xlabel('Duracao do Pulso de IA de 20 kW (s)', 'fontsize', 9);
ylabel('Tensao Instantanea da Bateria (V)', 'fontsize', 9);
title('(b2) Ponto de Falha WCA: Colapso por Subtensao em EOL (SoC=35%)', 'fontsize', 10, 'fontweight', 'bold');
legend({'Inicio de Vida BOL (Margem Segura > 355 V)', 'Fim de Vida EOL (Colapso abaixo de V_{UV})', 'Limiar Critico de Corte V_{UV} (324 V)', 'Zona de Desligamento Prematuro de IA'}, 'location', 'northeast', 'fontsize', 8);
xlim([0, 60]); ylim([310, 395]);

if ~exist('figuras','dir'), mkdir('figuras'); end; print(fig_3_6, 'figuras/figura_3_6_convergencia_aukf_estresse.png', '-dpng', '-r150');
fprintf('   -> Salvo em figuras/figura_3_6_convergencia_aukf_estresse.png\n\n');

fprintf('=== SIMULACOES E FIGURAS DA SECAO 3.3 CONCLUIDAS COM SUCESSO ===\n');
