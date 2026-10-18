% ==============================================================================
% GERACAO DE FIGURAS COMPOSTAS CONSOLIDADAS - SECAO 3.1 (BARRAMENTO HVDC 400 V)
% Projeto de Dissertacao de Mestrado - ProPECAut / IFES Campus Serra
% ==============================================================================

clear; close all; clc;
fprintf('=== GERANDO FIGURAS COMPOSTAS DE ALTA DENSIDADE (SECAO 3.1) ===\n');

if ismember('qt', available_graphics_toolkits())
    graphics_toolkit('qt');
else
    graphics_toolkit('gnuplot');
end

% ------------------------------------------------------------------------------
% FIGURA 3.1 COMPOSTA: BALANCO ENERGETICO ORBITAL E DINAMICA TRANSITORIA DE IA
% ------------------------------------------------------------------------------
fprintf('1. Gerando Figura 3.1 Composta (Balanco BESS + Degraus de IA)...\n');

fig_c1 = figure(1, 'position', [50, 50, 1100, 750]);

% Parametros do BESS do No Satelital de 50 kW e Orbitas Oficiais
V_nom = 400.0; V_UV = 324.0; V_OV = 442.8;
E_bat_nom_kWh = 9.0; % 9.0 kWh no satelital (Cluster 20 nos = 180 kWh) [PROVISORIO]
E_bat_nom_J = E_bat_nom_kWh * 3.6e6;
SoC_min_safe = 0.20; % [VERIFICAVEL / REQUISITO DE PROJETO]
R_0_eol = 0.085;     % Resistencia ohmica em EOL [PROVISORIO]

% Cenarios orbitais oficiais:
altitudes = [500.0, 700.0, 35786.0]; % LEO, SSO, GEO
cores = {'r', 'b', 'k'};
labels_alt = {'500 km (LEO Coplanar: Eclipse 35.8 min, Throttling P_{IT,ecl}=6.5 kW)', ...
              '700 km (SSO Dawn-Dusk Sazonal: Eclipse 18.8 min, P_{IT,ecl}=18.0 kW)', ...
              '35.786 km (GEO Equinocio: Eclipse 72.0 min, Throttling P_{IT,ecl}=4.5 kW)'};

mu_earth = 398600.4418;
R_earth = 6378.14;
N_orbits = 5;

% Execucao da integracao dinamica exata do BESS por cenario orbital:
sim_results = cell(3, 1);

for i = 1:3
    h = altitudes(i);
    r_orb = R_earth + h;
    T_o = 2 * pi * sqrt((r_orb^3) / mu_earth);
    
    if h == 500.0
        t_e = 35.8 * 60.0;
        P_IT = 6500.0;
    elseif h == 700.0
        t_e = 18.8 * 60.0;
        P_IT = 18000.0;
    else
        t_e = 72.0 * 60.0;
        P_IT = 4500.0;
    end
    
    % Livro-razao da carga total no barramento durante eclipse:
    P_load_ecl = (P_IT / 0.975) + 440.0 + 550.0 + 30.0;
    
    t_tot = N_orbits * T_o;
    t_vec = linspace(0, t_tot, 5000);
    dt_s = t_vec(2) - t_vec(1);
    
    soc_vec = zeros(size(t_vec));
    v_vec   = zeros(size(t_vec));
    soc_vec(1) = 0.95;
    Voc_init = 330.0 + 90.0 * soc_vec(1);
    v_vec(1) = Voc_init - (P_load_ecl / Voc_init) * R_0_eol;
    P_solar_avail = 48000.0;
    P_load_total_nom = (40000.0 / 0.975) + 440.0 + 550.0 + 65.0;
    
    for k = 2:length(t_vec)
        t_mod = mod(t_vec(k), T_o);
        is_ecl = (t_mod <= t_e);
        
        if is_ecl
            if soc_vec(k-1) <= SoC_min_safe
                P_b = 440.0 + 550.0; % Throttling de protecao de avionica
            else
                P_b = P_load_ecl;
            end
            dE = P_b * dt_s / 0.98;
            soc_vec(k) = max(SoC_min_safe, soc_vec(k-1) - (dE / E_bat_nom_J));
        else
            P_charge = min(P_solar_avail - P_load_total_nom, 10000.0);
            dE = P_charge * dt_s * 0.98;
            soc_vec(k) = min(0.95, soc_vec(k-1) + (dE / E_bat_nom_J));
            P_b = 0;
        end
        
        Voc_k = 330.0 + 90.0 * soc_vec(k);
        I_b = is_ecl * (P_b / Voc_k);
        v_vec(k) = Voc_k - I_b * R_0_eol;
    end
    
    sim_results{i}.t_norm = t_vec / T_o;
    sim_results{i}.soc = soc_vec;
    sim_results{i}.v = v_vec;
end

% Subpainel (a1): Estado de Carga (SoC) em 5 orbitas (saida da simulacao)
subplot(2, 2, 1);
hold on; grid on; box on;
for i = 1:3
    plot(sim_results{i}.t_norm, sim_results{i}.soc * 100, cores{i}, 'linewidth', 1.8);
end
plot([0, N_orbits], [20, 20], 'r--', 'linewidth', 1.5);
text(0.1, 24, 'Limite Minimo Seguranca SoC = 20% (DoD max = 80%)', 'color', 'r', 'fontsize', 8);
xlabel('Numero de Orbitas Consecutivas (Ciclos 1 a 5)', 'fontsize', 9);
ylabel('Estado de Carga SoC (%)', 'fontsize', 9);
title('(a1) Balanco Energetico BESS No 50 kW (9.0 kWh, 5 Orbitas)', 'fontsize', 10, 'fontweight', 'bold');
legend(labels_alt, 'location', 'southeast', 'fontsize', 7);
xlim([0, N_orbits]); ylim([10, 105]);

% Subpainel (a2): Tensao Terminal do Barramento (saida da simulacao)
subplot(2, 2, 2);
hold on; grid on; box on;
for i = 1:3
    plot(sim_results{i}.t_norm, sim_results{i}.v, cores{i}, 'linewidth', 1.8);
end
plot([0, N_orbits], [V_UV, V_UV], 'r--', 'linewidth', 1.5);
xlabel('Numero de Orbitas Consecutivas (Ciclos 1 a 5)', 'fontsize', 9);
ylabel('Tensao do Barramento (V)', 'fontsize', 9);
title('(a2) Tensao do Barramento com BESS em EOL (Modelo Dinamico)', 'fontsize', 10, 'fontweight', 'bold');
legend(labels_alt, 'location', 'southeast', 'fontsize', 7);
xlim([0, N_orbits]); ylim([310, 430]);

% Subpainel (b1): Trem de 5 Degraus de IA (2 kW <-> 40 kW em 10 us)
t_dyn = linspace(0, 15.0, 3000); % em ms
p_train = 2.0 * ones(size(t_dyn)); % kW
for p = 1:5
    t_p = (p - 1) * 2.8 + 1.0;
    p_train(t_dyn >= t_p & t_dyn < t_p + 1.2) = 40.0;
end

subplot(2, 2, 3);
hold on; grid on; box on;
plot(t_dyn, p_train, 'm-', 'linewidth', 1.8);
xlabel('Tempo (ms)', 'fontsize', 9);
ylabel('Potencia de Carga de IA (kW)', 'fontsize', 9);
title('(b1) Perfil de 5 Degraus de Carga de IA (2 a 40 kW em 10 \mus)', 'fontsize', 10, 'fontweight', 'bold');
xlim([0, 15]); ylim([0, 48]);

% Subpainel (b2): Resposta de Tensao e Subamortecimento (Integracao Numerica Exata)
subplot(2, 2, 4);
hold on; grid on; box on;

t_dyn_s = t_dyn * 1e-3;
dt_dyn = t_dyn_s(2) - t_dyn_s(1);
P_bus_train = (p_train * 1e3 / 0.900) + 440.0 + 550.0;

Cb = 1.5e-3; Lb = 4.0e-6;
R_cases = [0.010, 0.050];
colors_b2 = {'r-', 'b-'};
labels_b2 = {'R_{dist}=10 m\Omega (\zeta=0.097, Oscilatorio)', 'R_{dist}=50 m\Omega (\zeta=0.484, Amortecido)'};

for c = 1:2
    Rb = R_cases(c);
    V_b = zeros(size(t_dyn_s));
    I_b = zeros(size(t_dyn_s));
    V_b(1) = 400.0;
    I_b(1) = P_bus_train(1) / V_b(1);
    V_src = V_b(1) + Rb * I_b(1);
    
    for k = 2:length(t_dyn_s)
        dI = (V_src - V_b(k-1) - Rb * I_b(k-1)) / Lb;
        dV = (I_b(k-1) - (P_bus_train(k-1) / max(V_b(k-1), 50.0))) / Cb;
        I_b(k) = I_b(k-1) + dI * dt_dyn;
        V_b(k) = V_b(k-1) + dV * dt_dyn;
    end
    plot(t_dyn, V_b, colors_b2{c}, 'linewidth', 1.8);
end

xlabel('Tempo (ms)', 'fontsize', 9);
ylabel('Tensao do Barramento (V)', 'fontsize', 9);
title('(b2) Transiente de Tensao: Subamortecimento R_{dist}=10 m\Omega vs 50 m\Omega', 'fontsize', 10, 'fontweight', 'bold');
legend(labels_b2, 'location', 'southeast', 'fontsize', 8);
xlim([0, 15]); ylim([390, 408]);

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

% Subpainel (b1): Dissipacao Termica dos Conversores (eta = 90% vs 97.5% a 40 kW)
P_loss_cons = 40000.0 * (1.0 - 0.900) / 0.900; % 4444.44 W
P_loss_nom  = 40000.0 * (1.0 - 0.975) / 0.975; % 1025.64 W
t_term_h = linspace(0, 2.5, 400); % 2.5 horas
subplot(2, 2, 3);
hold on; grid on; box on;
plot([0, 2.5], [P_loss_cons, P_loss_cons], 'r-', 'linewidth', 2.2);
plot([0, 2.5], [P_loss_nom, P_loss_nom], 'g-', 'linewidth', 2.2);
xlabel('Tempo de Operacao (h)', 'fontsize', 9);
ylabel('Perdas no Chassi (W)', 'fontsize', 9);
title('(b1) Perdas Termicas: Conservador (90%) vs Nominal (97.5%)', 'fontsize', 10, 'fontweight', 'bold');
legend({sprintf('Cenario Conservador (\\eta=90%%, P_{loss}=%.0f W)', P_loss_cons), ...
        sprintf('Cenario Nominal GaN (\\eta=97.5%%, P_{loss}=%.0f W)', P_loss_nom)}, 'location', 'east', 'fontsize', 8);
xlim([0, 2.5]); ylim([0, 5200]);

% Subpainel (b2): Evolucao da Temperatura do Chassi com LHP
subplot(2, 2, 4);
hold on; grid on; box on;
C_th = 18000; % J/K
G_vec = [80, 120, 160]; % W/K
t_sec = t_term_h * 3600;

Q_cons = P_loss_cons + 440.0 * (1 - 0.22) + 50.0; % + perdas laser e chassi
Q_nom  = P_loss_nom  + 440.0 * (1 - 0.22) + 50.0;

for g = 1:3
    G_val = G_vec(g);
    tau = C_th / G_val;
    T_90 = 20.0 + (Q_cons / G_val) * (1 - exp(-t_sec / tau));
    T_97 = 20.0 + (Q_nom  / G_val) * (1 - exp(-t_sec / tau));
    plot(t_term_h, T_90, 'r--', 'linewidth', 1.5);
    plot(t_term_h, T_97, 'g-', 'linewidth', 1.8);
end
plot([0, 2.5], [65, 65], 'r:', 'linewidth', 1.5);
xlabel('Tempo de Operacao (h)', 'fontsize', 9);
ylabel('Temperatura do Chassi (C)', 'fontsize', 9);
title('(b2) Temperatura do Chassi Eletrico Acoplado ao LHP', 'fontsize', 10, 'fontweight', 'bold');
legend({'G=80 W/K (90%)', 'G=80 W/K (97.5%)', 'G=120 W/K (90%)', 'G=120 W/K (97.5%)', 'G=160 W/K (90%)', 'G=160 W/K (97.5%)', 'Limite Maximo (65 C)'}, 'location', 'southeast', 'fontsize', 7);
xlim([0, 2.5]); ylim([18, 85]);

if ~exist('figuras','dir'), mkdir('figuras'); end; print(fig_c2, 'figuras/figura_3_2_falhas_eletrotermicas_sspc.png', '-dpng', '-r150');
fprintf('   -> Figura 3.2 Composta salva com sucesso em figuras/figura_3_2_falhas_eletrotermicas_sspc.png\n\n');
fprintf('=== FIGURAS COMPOSTAS DA SECAO 3.1 CONCLUIDAS ===\n');
