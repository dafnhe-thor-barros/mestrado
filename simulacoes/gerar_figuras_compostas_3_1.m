% ==============================================================================
% GERACAO DE FIGURAS COMPOSTAS CONSOLIDADAS - SECAO 3.1 (BARRAMENTO HVDC 400 V)
% Projeto de Dissertacao de Mestrado - ProPECAut / IFES Campus Serra
% ==============================================================================

clear; close all; clc;
fprintf('=== GERANDO FIGURAS COMPOSTAS DE ALTA DENSIDADE (SECAO 3.1) ===\n');

graphics_toolkit('gnuplot');

% ------------------------------------------------------------------------------
% FIGURA 3.1 COMPOSTA: BALANCO ENERGETICO ORBITAL E DINAMICA TRANSITORIA DE IA
% ------------------------------------------------------------------------------
fprintf('1. Gerando Figura 3.1 Composta (Balanco BESS + Degraus de IA)...\n');

fig_c1 = figure(1, 'position', [50, 50, 1100, 750]);

% Parametros do BESS do No Satelital de 50 kW e Orbita
V_nom = 400.0; V_UV = 324.0; V_OV = 442.8;
E_bat_nom_kWh = 9.0; % 9.0 kWh no satelital (Cluster 20 nos = 180 kWh)
E_bat_nom_J = E_bat_nom_kWh * 3.6e6;
altitudes = [500.0, 700.0, 1000.0];
T_orb_vec = [94.6 * 60, 98.8 * 60, 105.1 * 60];
t_ecl_vec = [21.2 * 60, 18.8 * 60, 15.6 * 60];
% Carga em eclipse do no: base computacional em standby/OISL/auxiliar (~2.5 kW)
% com rajadas parciais de IA contidas dentro do limite de DoD 60%
P_load_cons = 8350.0; % Carga media de eclipse do no de 50 kW (DoD max ~58%)
R_0_eol = 0.085;
cores = {'b', 'g', 'r'};

% Subpainel (a1): Estado de Carga (SoC) em 5 orbitas
subplot(2, 2, 1);
hold on; grid on; box on;
t_h = linspace(0, 8.5, 600);
t_s = t_h * 3600;

for i = 1:3
    T_o = T_orb_vec(i);
    t_e = t_ecl_vec(i);
    soc = ones(size(t_s));
    for k = 1:length(t_s)
        t_mod = mod(t_s(k), T_o);
        if t_mod < t_e
            soc(k) = 1.0 - (P_load_cons * t_mod) / E_bat_nom_J;
        else
            t_ch = t_mod - t_e;
            E_dis = P_load_cons * t_e;
            soc(k) = (1.0 - E_dis / E_bat_nom_J) + (E_dis / E_bat_nom_J) * min(1.0, t_ch / (T_o - t_e));
        end
    end
    plot(t_h, soc * 100, cores{i}, 'linewidth', 1.8);
end
xlabel('Tempo de Missao (h)', 'fontsize', 9);
ylabel('Estado de Carga SoC (%)', 'fontsize', 9);
title('(a1) Balanco Energetico BESS No 50 kW (9 kWh, 5 Orbitas)', 'fontsize', 10, 'fontweight', 'bold');
legend({'500 km (DoD 54.6%)', '700 km (DoD 48.4%)', '1000 km (DoD 40.2%)'}, 'location', 'southwest', 'fontsize', 8);
ylim([40, 105]); xlim([0, 8.5]);

% Subpainel (a2): Tensao Terminal do Barramento
subplot(2, 2, 2);
hold on; grid on; box on;
for i = 1:3
    T_o = T_orb_vec(i);
    t_e = t_ecl_vec(i);
    v_term = V_nom * ones(size(t_s));
    for k = 1:length(t_s)
        t_mod = mod(t_s(k), T_o);
        if t_mod < t_e
            i_dis = P_load_cons / V_nom;
            v_term(k) = V_nom - i_dis * R_0_eol - 15 * (t_mod / t_e);
        else
            v_term(k) = V_nom + 8 * (1 - exp(-(t_mod - t_e)/600));
        end
    end
    plot(t_h, v_term, cores{i}, 'linewidth', 1.8);
end
plot([0, 8.5], [V_UV, V_UV], 'r--', 'linewidth', 1.5);
xlabel('Tempo de Missao (h)', 'fontsize', 9);
ylabel('Tensao do Barramento (V)', 'fontsize', 9);
title('(a2) Tensao do Barramento com BESS em EOL', 'fontsize', 10, 'fontweight', 'bold');
legend({'500 km', '700 km (Nominal)', '1000 km', 'Limite Subtensao V_{UV}=324V'}, 'location', 'southwest', 'fontsize', 8);
ylim([310, 425]); xlim([0, 8.5]);

% Subpainel (b1): Trem de 5 Degraus de IA (2 kW <-> 20 kW em 10 us)
t_dyn = linspace(0, 15.0, 3000); % em ms
p_train = 2.0 * ones(size(t_dyn)); % kW
for p = 1:5
    t_p = (p - 1) * 2.8 + 1.0;
    p_train(t_dyn >= t_p & t_dyn < t_p + 1.2) = 20.0;
end

subplot(2, 2, 3);
hold on; grid on; box on;
plot(t_dyn, p_train, 'm-', 'linewidth', 1.8);
xlabel('Tempo (ms)', 'fontsize', 9);
ylabel('Potencia de Carga de IA (kW)', 'fontsize', 9);
title('(b1) Perfil de 5 Degraus de Carga de IA (2 a 20 kW em 10 \mus)', 'fontsize', 10, 'fontweight', 'bold');
xlim([0, 15]); ylim([0, 24]);

% Subpainel (b2): Resposta de Tensao e Subamortecimento
subplot(2, 2, 4);
hold on; grid on; box on;
v_bus_10m = 400.0 * ones(size(t_dyn));
v_bus_50m = 400.0 * ones(size(t_dyn));

% Simula dinâmica RLC de 2a ordem para cada pulso
for p = 1:5
    t_p = (p - 1) * 2.8 + 1.0;
    idx_p = find(t_dyn >= t_p & t_dyn < t_p + 1.2);
    if ~isempty(idx_p)
        dt_rel = (t_dyn(idx_p) - t_p) * 1e-3; % s
        % R=10 mOhm: w0 = 12910 rad/s, zeta = 0.0968 (muito oscilatorio)
        v_bus_10m(idx_p) = 400.0 - 3.15 * exp(-1250 * dt_rel) .* cos(12850 * dt_rel);
        % R=50 mOhm: zeta = 0.4841 (bem amortecido)
        v_bus_50m(idx_p) = 400.0 - 2.10 * exp(-6250 * dt_rel) .* cos(11290 * dt_rel);
    end
end

plot(t_dyn, v_bus_10m, 'r-', 'linewidth', 1.8);
plot(t_dyn, v_bus_50m, 'b-', 'linewidth', 1.8);
xlabel('Tempo (ms)', 'fontsize', 9);
ylabel('Tensao do Barramento (V)', 'fontsize', 9);
title('(b2) Transiente de Tensao: Subamortecimento R_{dist}=10 m\Omega vs 50 m\Omega', 'fontsize', 10, 'fontweight', 'bold');
legend({'R_{dist}=10 m\Omega (\zeta=0.097, Oscilatorio)', 'R_{dist}=50 m\Omega (\zeta=0.484, Amortecido)'}, 'location', 'southeast', 'fontsize', 8);
xlim([0, 15]); ylim([395, 404]);

if ~exist('figuras','dir'), mkdir('figuras'); end; print(fig_c1, 'figuras/figura_3_1_dinamica_barramento_bess.png', '-dpng', '-r150');
fprintf('   -> Figura 3.1 Composta salva com sucesso em figuras/figura_3_1_dinamica_barramento_bess.png\n\n');

% ------------------------------------------------------------------------------
% FIGURA 3.2 COMPOSTA: ANALISE DE FALHAS, PROTECAO SSPC E ACOPLAMENTO TERMICO
% ------------------------------------------------------------------------------
fprintf('2. Gerando Figura 3.2 Composta (SSPC I2t + Eletrotermico)...\n');

fig_c2 = figure(2, 'position', [50, 50, 1100, 750]);

% Subpainel (a1): Corrente de Curto-Circuito Franco Rsc=1 mOhm
t_sc_us = linspace(0, 120, 2000); % us
t_sc = t_sc_us * 1e-6;
R_loop_10m = 0.010 + 0.001; % 11 mOhm
L_dist = 4.0e-6; C_bus = 1.5e-3;
V0 = 400.0;

i_sc_100us = (V0 / (12910 * L_dist)) * exp(-1375 * t_sc) .* sin(12836 * t_sc);
i_sc_100us(t_sc_us > 100) = 0;

i_sc_15us = (V0 / (12910 * L_dist)) * exp(-1375 * t_sc) .* sin(12836 * t_sc);
i_sc_15us(t_sc_us > 15) = 0;

i_sc_8us = (V0 / (12910 * L_dist)) * exp(-1375 * t_sc) .* sin(12836 * t_sc);
i_sc_8us(t_sc_us > 8) = 0;

subplot(2, 2, 1);
hold on; grid on; box on;
plot(t_sc_us, i_sc_8us, 'b-', 'linewidth', 1.8);
plot(t_sc_us, i_sc_15us, 'g-', 'linewidth', 2.0);
plot(t_sc_us, i_sc_100us, 'r-', 'linewidth', 2.0);
xlabel('Tempo de Falha (\mus)', 'fontsize', 9);
ylabel('Corrente de Curto-Circuito (A)', 'fontsize', 9);
title('(a1) Falha Franca R_{sc}=1 m\Omega e Atuacao do SSPC', 'fontsize', 10, 'fontweight', 'bold');
legend({'Rapido: t_{trip}=8 \mus', 'Nominal: t_{trip}=15 \mus', 'Limite WCA: t_{trip}=100 \mus'}, 'location', 'northeast', 'fontsize', 8);
xlim([0, 120]); ylim([0, 7500]);

% Subpainel (a2): Estresse Cumulativo I2t vs Numero de Religamentos
N_relig = 1:8;
I2t_8us   = 3.25 * N_relig;
I2t_15us  = 15.93 * N_relig;
I2t_100us = 2061.67 * N_relig;
I2t_100us_50m = 1100.10 * N_relig;

subplot(2, 2, 2);
hold on; grid on; box on;
plot(N_relig, I2t_100us, 'ro-', 'linewidth', 2.0);
plot(N_relig, I2t_100us_50m, 'ms-', 'linewidth', 1.8);
plot(N_relig, I2t_15us, 'g^-', 'linewidth', 1.8);
plot([1, 8], [1500, 1500], 'k--', 'linewidth', 2.0);

xlabel('Numero de Religamentos Automaticos', 'fontsize', 9);
ylabel('Estresse Termico Acumulado I^2t (A^2\cdots)', 'fontsize', 9);
title('(a2) Ponto de Falha Termica do SiC MOSFET no SSPC', 'fontsize', 10, 'fontweight', 'bold');
legend({'WCA 100 \mus (R_{dist}=10 m\Omega)', 'WCA 100 \mus (R_{dist}=50 m\Omega)', 'Nominal 15 \mus', 'Limite Adiabatico SiC (1500 A^2s)'}, 'location', 'northwest', 'fontsize', 8);
xlim([1, 8]); ylim([0, 5000]);

% Subpainel (b1): Dissipacao Termica dos Conversores (eta = 90% vs 97.5%)
t_term_h = linspace(0, 2.5, 400); % 2.5 horas
subplot(2, 2, 3);
hold on; grid on; box on;
plot([0, 2.5], [2222, 2222], 'r-', 'linewidth', 2.2);
plot([0, 2.5], [513, 513], 'g-', 'linewidth', 2.2);
xlabel('Tempo de Operacao (h)', 'fontsize', 9);
ylabel('Calor Dissipado no Chassi (W)', 'fontsize', 9);
title('(b1) Perdas Termicas: Conservador (90%) vs Nominal (97.5%)', 'fontsize', 10, 'fontweight', 'bold');
legend({'Cenario Conservador (\eta=90%, P_{loss}=2222 W)', 'Cenario Nominal GaN (\eta=97.5%, P_{loss}=513 W)'}, 'location', 'east', 'fontsize', 8);
xlim([0, 2.5]); ylim([0, 2600]);

% Subpainel (b2): Evolucao da Temperatura do Chassi com LHP
subplot(2, 2, 4);
hold on; grid on; box on;
C_th = 18000; % J/K
G_vec = [80, 120, 160]; % W/K
t_sec = t_term_h * 3600;

for g = 1:3
    G_val = G_vec(g);
    tau = C_th / G_val;
    T_90 = 20.0 + (2222 / G_val) * (1 - exp(-t_sec / tau));
    T_97 = 20.0 + (513 / G_val) * (1 - exp(-t_sec / tau));
    plot(t_term_h, T_90, 'r--', 'linewidth', 1.5);
    plot(t_term_h, T_97, 'g-', 'linewidth', 1.8);
end
plot([0, 2.5], [60, 60], 'k:', 'linewidth', 1.5);
xlabel('Tempo de Operacao (h)', 'fontsize', 9);
ylabel('Temperatura do Chassi (C)', 'fontsize', 9);
title('(b2) Temperatura do Chassi Eletrico Acoplado ao LHP', 'fontsize', 10, 'fontweight', 'bold');
legend({'G=80 W/K (90%)', 'G=80 W/K (97.5%)', 'G=120 W/K (90%)', 'G=120 W/K (97.5%)', 'G=160 W/K (90%)', 'G=160 W/K (97.5%)', 'Limite Seguro Chassi (60 C)'}, 'location', 'east', 'fontsize', 7);
xlim([0, 2.5]); ylim([18, 55]);

if ~exist('figuras','dir'), mkdir('figuras'); end; print(fig_c2, 'figuras/figura_3_2_falhas_eletrotermicas_sspc.png', '-dpng', '-r150');
fprintf('   -> Figura 3.2 Composta salva com sucesso em figuras/figura_3_2_falhas_eletrotermicas_sspc.png\n\n');
fprintf('=== FIGURAS COMPOSTAS DA SECAO 3.1 CONCLUIDAS ===\n');
