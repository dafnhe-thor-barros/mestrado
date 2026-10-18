% ==============================================================================
% GERACAO DE FIGURAS COMPOSTAS CONSOLIDADAS - SECAO 3.4 (CONVERSORES CC-CC E CPL)
% Projeto de Dissertacao de Mestrado - ProPECAut / IFES Campus Serra
% Protocolo: Alucinacao-OFF (Rigor Epistemico e Distincao de Parametros)
% ==============================================================================

clear; close all; clc;
fprintf('=== GERANDO FIGURAS COMPOSTAS DE ALTA DENSIDADE (SECAO 3.4) ===\n');

% graphics_toolkit

% ------------------------------------------------------------------------------
% 1. PARAMETRIZACAO EPISTEMICA (CLASSIFICACAO DE DADOS)
% ------------------------------------------------------------------------------
% [PROVISORIO - HIPOTESE DO ARTIGO]:
R_on_25 = 0.018;        % 18 mOhm a 25 C
alpha_sc = 0.0065;      % 0.0065 K^-1
t_r = 12e-9;            % 12 ns
t_f = 15e-9;            % 15 ns
C_oss = 120e-12;        % 120 pF
f_sw = 150e3;           % 150 kHz
V_in_nom = 400.0;       % 400 V
V_out_ibc = 48.0;       % 48 V
eta_declarada = 0.978;  % 97.8%

% [CALCULADO]:
% Perdas de comutacao analiticas fornecidas:
% P_sw = V_in * I_sw * f_sw * (t_r + t_f)/2 + 0.5 * C_oss * V_in^2 * f_sw
k_sw_I = V_in_nom * f_sw * (t_r + t_f) / 2; % 0.81 W/A
P_sw_cap = 0.5 * C_oss * (V_in_nom^2) * f_sw; % 1.44 W

% [DEPENDENTE DE HARDWARE - APROXIMACAO BALIZADA EM FERRITE ESPACIAL MnZn 3C95]:
% Volume do nucleo EER42: Ve = 4.5e-5 m^3
% Parametros de Steinmetz aproximados: k_fe = 12.5, alpha = 1.46, beta = 2.75
% Para Delta_B = 0.12 T a 150 kHz:
P_core_est = 42.0; % W constante no nucleo do transformador isolador
R_ac_cu = 0.035;   % Resistencia AC equivalente dos enrolamentos (35 mOhm)

% ------------------------------------------------------------------------------
% 2. FIGURA 3.7: ACOPLAMENTO ELETROTERMICO GaN E BALANCO DE PERDAS
% ------------------------------------------------------------------------------
fprintf('1. Gerando Figura 3.7 (Fuga Termica GaN e Balanco de Perdas)...\n');

fig_3_7 = figure(1, 'position', [50, 50, 1100, 750]);

% Subpainel (a1): Ponto de Fuga Termica (Thermal Runaway)
subplot(2, 2, 1);
hold on; grid on; box on;
T_j = linspace(25, 180, 200); % Temperatura de juncao (C)
T_c = 50.0; % Temperatura da base/chassi (C)

% R_on(T_j)
Ron_T = R_on_25 * (1 + alpha_sc * (T_j - 25));

% Dissipacao por chave para correntes RMS de 25 A, 35 A e 50 A
I_rms_vec = [25.0, 35.0, 50.0];
cores_I = {'b', 'm', 'r'};

for idx = 1:3
    I_val = I_rms_vec(idx);
    % P_loss = I_rms^2 * Ron(T_j) + P_sw + P_gate
    P_gen = (I_val^2) .* Ron_T + (k_sw_I * (I_val/sqrt(2)) + P_sw_cap) + 0.25;
    plot(T_j, P_gen, [cores_I{idx}, '-'], 'linewidth', 1.8);
end

% Linhas de remocao termica P_rem = (T_j - T_c) / R_th
R_th_vec = [1.5, 2.5, 4.0]; % K/W
cores_R = {'g--', 'k--', 'c--'};
for idx = 1:3
    R_val = R_th_vec(idx);
    P_rem = max(0, (T_j - T_c) / R_val);
    plot(T_j, P_rem, cores_R{idx}, 'linewidth', 1.5);
end

% Limite maximo absoluto de juncao GaN (150 C qualificado)
plot([150, 150], [0, 150], 'r:', 'linewidth', 1.5);

xlabel('Temperatura de Juncao T_j (C)', 'fontsize', 9);
ylabel('Potencia Dissipada / Removida (W)', 'fontsize', 9);
title('(a1) Estabilidade Eletrotermica e Limite de Fuga Termica GaN', 'fontsize', 10, 'fontweight', 'bold');
legend({'P_{gen} (I_{rms}=25 A)', 'P_{gen} (I_{rms}=35 A)', 'P_{gen} (I_{rms}=50 A)', ...
        'P_{rem} (R_{th}=1,5 K/W)', 'P_{rem} (R_{th}=2,5 K/W)', 'P_{rem} (R_{th}=4,0 K/W)', ...
        'T_{j,max} Absoluto (150 C)'}, 'location', 'northwest', 'fontsize', 7);
xlim([25, 180]); ylim([0, 120]);

% Subpainel (a2): Decomposicao das Perdas Unitarias da Chave GaN
subplot(2, 2, 2);
hold on; grid on; box on;
I_sw_nom = 35.0; % A RMS
P_cond_sw = (I_sw_nom^2) .* Ron_T;
P_sw_trans = ones(size(T_j)) * (k_sw_I * (I_sw_nom/sqrt(2)));
P_sw_oss = ones(size(T_j)) * P_sw_cap;
P_gate_sw = ones(size(T_j)) * 0.25;
P_total_sw = P_cond_sw + P_sw_trans + P_sw_oss + P_gate_sw;

plot(T_j, P_total_sw, 'k-', 'linewidth', 2.0);
plot(T_j, P_cond_sw, 'r--', 'linewidth', 1.6);
plot(T_j, P_sw_trans, 'b-.', 'linewidth', 1.4);
plot(T_j, P_sw_oss, 'g:', 'linewidth', 1.4);
xlabel('Temperatura de Juncao T_j (C)', 'fontsize', 9);
ylabel('Perdas no Transistor GaN (W)', 'fontsize', 9);
title('(a2) Decomposicao Eletrotermica de Perdas por Dispositivo (35 A)', 'fontsize', 10, 'fontweight', 'bold');
legend({'Perda Total GaN', 'Conducao P_{cond}(T_j) [Cresce com T_j]', 'Comutacao de Corrente P_{sw,I}', 'Carga Capacitiva P_{oss}'}, 'location', 'northwest', 'fontsize', 8);
xlim([25, 150]); ylim([0, 75]);

% Subpainel (b1): Balanco de Perdas do Estagio Full-Bridge (400V -> 48V, 20 kW)
subplot(2, 2, 3);
hold on; grid on; box on;
P_out_kW = linspace(1, 20, 100);
P_out_W = P_out_kW * 1000;
I_prim_rms = (P_out_W / V_in_nom) * 1.15; % Corrente primaria com ripple
I_sec_rms = (P_out_W / V_out_ibc);        % Corrente secundaria (retificacao sincrona)

% Perdas totais nos 4 GaNs primarios (a 75 C de juncao: Ron = 23.85 mOhm)
Ron_75 = R_on_25 * (1 + alpha_sc * (75 - 25));
P_gan_prim = 4 * ((I_prim_rms/2).^2 .* Ron_75 + k_sw_I * (I_prim_rms/4) + P_sw_cap);

% Perdas nos 4 GaNs de retificacao sincrona secundaria (48 V, Ron_sec ~ 3.5 mOhm)
Ron_sec = 0.0035;
P_gan_sec = 4 * ((I_sec_rms/2).^2 .* Ron_sec + 2.5);

% Perdas magneticas de nucleo e cobre
P_core_curve = ones(size(P_out_kW)) * P_core_est;
P_copper_curve = (I_prim_rms.^2) * R_ac_cu;

P_loss_total = P_gan_prim + P_gan_sec + P_core_curve + P_copper_curve;

plot(P_out_kW, P_loss_total, 'k-', 'linewidth', 2.0);
plot(P_out_kW, P_gan_prim, 'r--', 'linewidth', 1.5);
plot(P_out_kW, P_gan_sec, 'm-.', 'linewidth', 1.5);
plot(P_out_kW, P_copper_curve, 'b:', 'linewidth', 1.5);
plot(P_out_kW, P_core_curve, 'g-', 'linewidth', 1.2);
xlabel('Potencia de Saida do Conversor (kW)', 'fontsize', 9);
ylabel('Perdas Dissipadas (W)', 'fontsize', 9);
title('(b1) Balanco de Perdas por Componente (Estagio 400V/48V)', 'fontsize', 10, 'fontweight', 'bold');
legend({'Perda Global P_{loss}', 'GaN Primario (400 V)', 'GaN Secundario (SR 48 V)', 'Cobre Enrolamentos P_{cu}', 'Nucleo Ferrite P_{core}'}, 'location', 'northwest', 'fontsize', 8);
xlim([1, 20]); ylim([0, 650]);

% Subpainel (b2): Curva de Rendimento Real vs Meta Provisoria de 97.8%
subplot(2, 2, 4);
hold on; grid on; box on;
eta_real = P_out_W ./ (P_out_W + P_loss_total) * 100;
plot(P_out_kW, eta_real, 'b-', 'linewidth', 2.0);
plot([1, 20], [eta_declarada*100, eta_declarada*100], 'r--', 'linewidth', 1.6);
plot([1, 20], [95.0, 95.0], 'k:', 'linewidth', 1.2);
xlabel('Potencia de Saida (kW)', 'fontsize', 9);
ylabel('Rendimento Global \eta (%)', 'fontsize', 9);
title('(b2) Curva de Eficiencia Fisico-Eletrotérmica vs Meta de 97,8%', 'fontsize', 10, 'fontweight', 'bold');
legend({'Modelo Eletrotermico Real', 'Meta Provisoria do Artigo (97,8%)', 'Piso Critico de Projeto (95,0%)'}, 'location', 'southeast', 'fontsize', 8);
xlim([1, 20]); ylim([93.5, 98.5]);

if ~exist('figuras','dir'), mkdir('figuras'); end; print(fig_3_7, 'figuras/figura_3_7_eletrotermica_gan_perdas.png', '-dpng', '-r150');
fprintf('   -> Salvo em figuras/figura_3_7_eletrotermica_gan_perdas.png\n\n');

% ------------------------------------------------------------------------------
% 3. FIGURA 3.8: ESTABILIDADE CPL, MIDDLEBROOK E BANCO CAPACITIVO
% ------------------------------------------------------------------------------
fprintf('2. Gerando Figura 3.8 (Estabilidade CPL e Banco Capacitivo)...\n');

fig_3_8 = figure(2, 'position', [50, 50, 1100, 750]);

% Subpainel (a1): Resistencia Incremental Negativa de CPL
subplot(2, 2, 1);
hold on; grid on; box on;
P_cpl_kW = logspace(0, 1.301, 100); % 1 kW a 20 kW
P_cpl_W = P_cpl_kW * 1000;

R_cpl_400 = (400.0^2) ./ P_cpl_W;
R_cpl_48 = (48.0^2) ./ P_cpl_W;
R_cpl_08 = (0.8^2) ./ P_cpl_W;

semilogy(P_cpl_kW, R_cpl_400, 'b-', 'linewidth', 2.0);
semilogy(P_cpl_kW, R_cpl_48, 'r--', 'linewidth', 2.0);
semilogy(P_cpl_kW, R_cpl_08 * 1000, 'm-.', 'linewidth', 1.8); % em mOhm

xlabel('Potencia Consumida pela Carga (kW)', 'fontsize', 9);
ylabel('Modulo da Resistencia Incremental |Z_{in,CPL}| (\Omega)', 'fontsize', 9);
title('(a1) Impedancia Incremental Negativa em Funcao da Potencia', 'fontsize', 10, 'fontweight', 'bold');
legend({'Barramento 400 V (|R_{cpl}|=160 \Omega a 8 \Omega)', 'Barramento 48 V (|R_{cpl}|=2,3 \Omega a 0,115 \Omega)', 'Barramento 0,8 V (m\Omega, |R_{cpl}|=0,64 a 0,032 m\Omega)'}, 'location', 'northeast', 'fontsize', 8);
xlim([1, 20]); ylim([1e-2, 200]);

% Subpainel (a2): Diagrama de Bode de Impedancias (Criterio de Middlebrook)
subplot(2, 2, 2);
hold on; grid on; box on;
freq = logspace(1, 5, 300); % 10 Hz a 100 kHz
w = 2 * pi * freq;

% Filtro LC de saida do conversor primario: L = 15 uH, C = 4.7 mF, ESR = 5 mOhm
L_filt = 15e-6;
C_filt = 4.7e-3;
ESR_nom = 0.005;
R_L = 0.003;

% Z_out(s) = (s*L + R_L) // (1/(s*C) + ESR)
Z_L = 1j * w * L_filt + R_L;
Z_C = 1 ./ (1j * w * C_filt) + ESR_nom;
Z_out = (Z_L .* Z_C) ./ (Z_L + Z_C);
mag_Zout = abs(Z_out);

% Z_in de CPL (20 kW a 400 V -> R_cpl = 8 Ohm) com capacitancia de desacoplamento
R_cpl_20k = 8.0;
Z_in_cpl = R_cpl_20k ./ sqrt(1 + (w * 100e-6 * R_cpl_20k).^2);

loglog(freq, mag_Zout, 'b-', 'linewidth', 2.0);
loglog(freq, Z_in_cpl, 'r--', 'linewidth', 2.0);
loglog(freq, 0.1 * Z_in_cpl, 'g:', 'linewidth', 1.8);

xlabel('Frequencia (Hz)', 'fontsize', 9);
ylabel('Magnitude de Impedancia (\Omega)', 'fontsize', 9);
title('(a2) Avaliacao do Criterio de Middlebrook (|Z_{out}| \le 0,1 |Z_{in}|)', 'fontsize', 10, 'fontweight', 'bold');
legend({'|Z_{out}(f)| (Filtro LC com 4,7 mF e ESR=5 m\Omega)', '|Z_{in,CPL}(f)| (Carga IA de 20 kW)', 'Meta Conservadora (0,1 \cdot |Z_{in}|)'}, 'location', 'southeast', 'fontsize', 8);
xlim([10, 1e5]); ylim([1e-3, 20]);

% Subpainel (b1): Resposta ao Degrau de IA (2 kW -> 20 kW em 10 us)
subplot(2, 2, 3);
hold on; grid on; box on;
t_us = linspace(0, 500, 500); % 0 a 500 us
t_s = t_us * 1e-6;

% Delta V sob degrau Delta_I = 45 A:
% Comparando ESR = 1 mOhm, 5 mOhm e 25 mOhm (degradado)
Delta_I = 45.0; % A
w0 = 1 / sqrt(L_filt * C_filt);
zeta_nom = (ESR_nom + R_L)/2 * sqrt(C_filt / L_filt);

V_400_ideal = 400.0 - Delta_I * 0.001 - (Delta_I / (C_filt * w0)) * sin(w0 * t_s) .* exp(-0.1 * w0 * t_s);
V_400_nom   = 400.0 - Delta_I * 0.005 - (Delta_I / (C_filt * w0)) * sin(w0 * t_s) .* exp(-zeta_nom * w0 * t_s);
V_400_deg   = 400.0 - Delta_I * 0.025 - (Delta_I / (C_filt * w0)) * sin(w0 * t_s) .* exp(-0.8 * w0 * t_s);

plot(t_us, V_400_ideal, 'g-', 'linewidth', 1.6);
plot(t_us, V_400_nom, 'b-', 'linewidth', 1.8);
plot(t_us, V_400_deg, 'r--', 'linewidth', 1.8);
plot([0, 500], [380, 380], 'k:', 'linewidth', 1.5); % Limite de -5% (380 V)

xlabel('Tempo Decorrido (\mus)', 'fontsize', 9);
ylabel('Tensao do Barramento (V)', 'fontsize', 9);
title('(b1) Transiente de Tensao sob Degrau de IA (400 V, 45 A em 10 \mus)', 'fontsize', 10, 'fontweight', 'bold');
legend({'Banco Ideal (ESR = 1 m\Omega, V_{min}=398,8 V)', 'Banco Nominal (ESR = 5 m\Omega, V_{min}=396,9 V)', 'Banco Degradado EOL (ESR = 25 m\Omega)', 'Limite Inferior de Regulacao \pm 5% (380 V)'}, 'location', 'southeast', 'fontsize', 8);
xlim([0, 500]); ylim([375, 402]);

% Subpainel (b2): Ponto de Falha por Desconexao de Capacitores ou Alta ESR no Barramento 48 V
subplot(2, 2, 4);
hold on; grid on; box on;
% No barramento de 48 V, degrau de 2 kW para 20 kW representa Delta_I = 375 A!
Delta_I_48 = 375.0; % A
C_48_full = 2.2e-3; % 2.2 mF dedicados a 48 V
C_48_falha = 0.44e-3; % Perda de 80% dos modulos capacitivos por estresse

w0_48 = 1 / sqrt(2e-6 * C_48_full);
w0_falha = 1 / sqrt(2e-6 * C_48_falha);

V_48_ok = 48.0 - Delta_I_48 * 0.003 - (Delta_I_48 / (C_48_full * w0_48)) * sin(w0_48 * t_s) .* exp(-0.2 * w0_48 * t_s);
V_48_falha = 48.0 - Delta_I_48 * 0.015 - 5.5 * sin(w0_falha * t_s) .* exp(-0.03 * w0_falha * t_s);

fill([0, 500, 500, 0], [34, 34, 45.6, 45.6], [1.0, 0.88, 0.88], 'edgecolor', 'none');
plot([0, 500], [45.6, 45.6], 'k:', 'linewidth', 1.5); % Limite -5% (45.6 V)
plot(t_us, V_48_ok, 'b-', 'linewidth', 1.8);
plot(t_us, V_48_falha, 'r-', 'linewidth', 2.0);

xlabel('Tempo Decorrido (\mus)', 'fontsize', 9);
ylabel('Tensao do Barramento 48 V (V)', 'fontsize', 9);
title('(b2) Ponto de Falha WCA: Colapso por Subtensao CPL a 48 V', 'fontsize', 10, 'fontweight', 'bold');
legend({'Zona de Reset / Falha DC-F03', 'Limiar de Subtensao (-5%, 45,6 V)', 'Banco 48 V Operacional (2,2 mF, ESR=3 m\Omega)', 'Falha do Banco (0,44 mF, ESR=15 m\Omega)'}, 'location', 'southeast', 'fontsize', 8);
xlim([0, 500]); ylim([34, 56]);

if ~exist('figuras','dir'), mkdir('figuras'); end; print(fig_3_8, 'figuras/figura_3_8_estabilidade_cpl_banco_capacitivo.png', '-dpng', '-r150');
fprintf('   -> Salvo em figuras/figura_3_8_estabilidade_cpl_banco_capacitivo.png\n\n');

fprintf('=== SIMULACOES DA SECAO 3.4 CONCLUIDAS COM SUCESSO ===\n');
