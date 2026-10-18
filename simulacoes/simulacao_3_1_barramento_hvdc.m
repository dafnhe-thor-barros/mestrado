%% ==============================================================================
%% PROPECAUT / IFES CAMPUS SERRA - DISSERTAÇÃO DE MESTRADO
%% MODELO 3.1: BARRAMENTO HVDC DE 400 V - ANÁLISE DE ESTRESSE EXTREMO (WCA)
%% Autor: Dafnhê Thor Castro de Martins Barros
%% Orientador: Prof. Dr. Flávio Barcelos Braz da Silva
%%
%% NOVAS ANÁLISES DE ESTRESSE:
%%   1. Trem de 5 Degraus Extremos de IA com R_dist = 10 mOhms, C_bus = 1.5 mF, L_dist = 4 uH
%%      (Análise de Subamortecimento extremo zeta = 0.0968 e ressonância CPL)
%%   2. Simulação do Número de Religamentos Automáticos até violação do limite I^2t
%%      (Limite térmico adiabático SiC MOSFET I^2t_max = 1500 A^2.s)
%% ==============================================================================

clear; close all; clc;
fprintf('=== INICIANDO SIMULAÇÃO NUMÉRICA: ANÁLISE DE RELIGAMENTOS E SUBAMORTECIMENTO ===\n\n');

if ismember('qt', available_graphics_toolkits())
    graphics_toolkit('qt');
else
    graphics_toolkit('gnuplot');
end

%% ==============================================================================
%% 1. DEFINIÇÃO DA MATRIZ DE PARÂMETROS
%% ==============================================================================

V_nom = 400.0;
V_UV = 324.0; % Limite estrito de subtensao (108s x 3.00 V/cel) [VERIFICAVEL]
V_OV = 442.8; % Limite estrito de sobretensao (108s x 4.10 V/cel) [VERIFICAVEL]

% Livro-Razao Eletrico Nominal do No Satelital (Fronteira do Barramento 400 V):
P_IT_nom = 40000.0;     % [W] Potencia util de processamento neural (80 dies x 500 W) [VERIFICAVEL]
P_IT_idle = 2000.0;     % [W] Modo standby/ocioso [PROVISORIO]
eta_conv_nom  = 0.975;  % Rendimento nominal GaN (perda de 1.025,64 W) [PROVISORIO]
eta_conv_cons = 0.900;  % Rendimento degradado WCH (perda de 4.444,44 W) [PROVISORIO]

P_OISL_total = 440.0;   % [W] 4 transceivers laser no barramento [VERIFICAVEL]
eta_eo = 0.22;          % Rendimento eletro-optico [PROVISORIO]
P_OISL_opt = P_OISL_total * eta_eo;        % 96.8 W irradiados ao espaco [CALCULADO]
Q_OISL_heat = P_OISL_total * (1.0 - eta_eo); % 343.2 W calor interno [CALCULADO]
Q_OISL = Q_OISL_heat;

P_aux = 550.0;          % [W] Auxiliares: 350 W avionica/OBDH + 200 W bomba MPFL [PROVISORIO]
P_loss_BESS_nom = 65.0; % [W] Perdas ohmicas internas em flutuacao [CALCULADO]

% Consumo eletrico total no barramento em operacao nominal plena:
% 40000/0.975 + 440 + 550 + 65 = 42.080,64 W
P_load_total_nom = (P_IT_nom / eta_conv_nom) + P_OISL_total + P_aux + P_loss_BESS_nom;
P_load_total_cons = (P_IT_nom / eta_conv_cons) + P_OISL_total + P_aux + 120.0; % WCH

% BESS do No Satelital de 50 kW
E_bat_nom_kWh = 9.0;    % Capacidade nominal de 9.0 kWh [PROVISORIO]
E_bat_nom_J = E_bat_nom_kWh * 3.6e6; 
SoC_min_safe = 0.20;    % SoC minimo operacional (DoD max 80%) [VERIFICAVEL / ECSS]
E_util_J = E_bat_nom_J * (1.0 - SoC_min_safe); % 7.20 kWh utilizaveis ideais
R_0_vec = [0.032, 0.045, 0.085]; % [Ohm] BOL a EOL [PROVISORIO]

% Parâmetros de distribuição e filtro
C_bus_ext = 1.5e-3;   % 1.5 mF
L_dist_ext = 4.0e-6;  % 4.0 uH
R_dist_low = 0.010;   % 10 mOhms (Condição subamortecida)
R_dist_high = 0.050;  % 50 mOhms (Condição resistiva EOL)
R_dist_nom = 0.014;   % 14 mOhms
C_bus_nom = 4.7e-3;   % 4.7 mF
L_dist_nom = 1.8e-6;  % 1.8 uH

% Curto-circuito franco ultra-severo
R_sc = 0.001; % 1 mOhm
t_trip_vec = [8.0e-6, 15.0e-6, 100.0e-6]; 
I2t_limit_SiC = 1500.0; % Limite térmico adiabático admissível do SiC MOSFET (A^2.s)

% Modelo Térmico
C_th_elec = 18000.0; 
G_th_vec = [80.0, 120.0, 160.0]; 
T_rad = 20.0 + 273.15; 

fprintf('1. Parametros eletricos e orbitais configurados.\n\n');

%% ==============================================================================
%% NÍVEL 1: BALANÇO ENERGÉTICO ORBITAL MULTICÍCLICO (5 ÓRBITAS CONSECUTIVAS)
%% ==============================================================================
fprintf('2. Balanco Energetico Orbital ao longo de 5 ORBITAS CONSECUTIVAS...\n');

mu_earth = 398600.4418;
R_earth = 6378.14;
altitudes = [500.0, 700.0, 35786.0]; % km (LEO 500 km, SSO 700 km, GEO 35.786 km)
cores = {'r', 'b', 'k'};

% Cenarios orbitais oficiais da dissertacao:
% 500 km: LEO convencional, pior caso coplanar (beta=0), eclipse 35.8 min (T_orb = 94.6 min)
% 700 km: SSO Dawn-Dusk, nominal sem eclipse (beta=85 deg); sazonal/equinocio eclipse max 18.8 min (T_orb = 98.8 min)
% 35.786 km: GEO equatorial, eclipse de equinocio maximo 72.0 min (T_orb = 1440.0 min = 24.0 h)
labels_alt = {'500 km (LEO Coplanar: Eclipse 35.8 min, Throttling P_{IT,ecl}=6.5 kW)', ...
              '700 km (SSO Dawn-Dusk Sazonal: Eclipse 18.8 min, P_{IT,ecl}=18.0 kW)', ...
              '35.786 km (GEO Equinocio: Eclipse 72.0 min, Throttling P_{IT,ecl}=4.5 kW)'};

fig1 = figure(1, 'position', [100, 100, 1000, 700]);
N_orbits = 5;

for i = 1:length(altitudes)
    h = altitudes(i);
    r_orb = R_earth + h;
    T_orb_s = 2 * pi * sqrt((r_orb^3) / mu_earth);
    
    if h == 500.0
        t_ecl_min = 35.8;
        P_IT_ecl = 6500.0; % Gestao de carga obrigatoria para respeitar SoC >= 20%
    elseif h == 700.0
        t_ecl_min = 18.8;
        P_IT_ecl = 18000.0;
    else
        % GEO (35.786 km): eclipse de equinocio maximo de 72 min (1.2 h)
        t_ecl_min = 72.0;
        P_IT_ecl = 4500.0; % Carga reduzida para respeitar SoC >= 20% com 9.0 kWh
    end
    t_ecl_s = t_ecl_min * 60.0;
    P_load_ecl = (P_IT_ecl / eta_conv_nom) + P_OISL_total + P_aux + 30.0;
    
    t_total_s = N_orbits * T_orb_s;
    t_sim = linspace(0, t_total_s, 5000);
    dt = t_sim(2) - t_sim(1);
    
    SoC = zeros(size(t_sim));
    V_bat = zeros(size(t_sim));
    SoC(1) = 0.95; 
    V_ocv_init = 330.0 + 90.0 * SoC(1);
    is_ecl_1 = (0 <= t_ecl_s);
    V_bat(1) = V_ocv_init - is_ecl_1 * (P_load_ecl / V_ocv_init) * R_dist_high;
    P_solar_avail = 48000.0; % 48 kW disponiveis no barramento em temperatura de equilibrio
    
    for k = 2:length(t_sim)
        t_mod = mod(t_sim(k), T_orb_s);
        is_eclipse = (t_mod <= t_ecl_s);
        
        if is_eclipse
            % Verificacao e atuacao de protecao do SoC minimo seguro (20% / 80% DoD max)
            if SoC(k-1) <= SoC_min_safe
                P_bat = P_aux + P_OISL_total; % Throttling estrito para manter avionica
            else
                P_bat = P_load_ecl;
            end
            dE = P_bat * dt / 0.98; % Rendimento de descarga de 98%
            SoC(k) = max(SoC_min_safe, SoC(k-1) - (dE / E_bat_nom_J));
        else
            % Recarga com excedente solar
            P_charge = min(P_solar_avail - P_load_total_nom, 10000.0);
            dE = P_charge * dt * 0.98; % Rendimento de carga de 98%
            SoC(k) = min(0.95, SoC(k-1) + (dE / E_bat_nom_J));
            P_bat = 0;
        end
        
        V_ocv = 330.0 + 90.0 * SoC(k);
        I_bat = is_eclipse * (P_bat / V_ocv);
        V_bat(k) = V_ocv - I_bat * R_dist_high;
    end
    
    subplot(2, 1, 1);
    plot(t_sim / T_orb_s, SoC * 100.0, 'linewidth', 1.8, 'color', cores{i}); hold on;
    
    subplot(2, 1, 2);
    plot(t_sim / T_orb_s, V_bat, 'linewidth', 1.8, 'color', cores{i}); hold on;
end

subplot(2, 1, 1);
grid on;
title('Balanco Energetico do BESS (9.0 kWh, 400 V) em 5 Orbitas Consecutivas com Gestao de Carga', 'fontsize', 11);
xlabel('Numero de Orbitas Consecutivas (Ciclos 1 a 5)', 'fontsize', 10);
ylabel('Estado de Carga - SoC (%)', 'fontsize', 10);
plot([0, N_orbits], [20, 20], 'r--', 'linewidth', 1.5);
text(0.1, 23, 'Limite Minimo de Seguranca: SoC = 20% (Requisito de Projeto / DoD = 80%)', 'color', 'r', 'fontsize', 9);
legend(labels_alt, 'location', 'southeast', 'fontsize', 8);
xlim([0, N_orbits]); ylim([10, 105]);

subplot(2, 1, 2);
grid on;
title('Tensao do Barramento Durante 5 Ciclos Consecutivos (R_{0} EOL)', 'fontsize', 11);
xlabel('Numero de Orbitas Consecutivas (Ciclos 1 a 5)', 'fontsize', 10);
ylabel('Tensao Terminal (V)', 'fontsize', 10);
plot([0, N_orbits], [V_UV, V_UV], 'r--', 'linewidth', 1.5);
text(0.1, V_UV + 6, 'Limite de Subtensao V_{UV} = 324 V', 'color', 'r', 'fontsize', 9);
xlim([0, N_orbits]); ylim([310, 430]);
legend(labels_alt, 'location', 'southeast', 'fontsize', 8);

print(fig1, 'figuras/figura_3_1_balanco_orbital_bess.png', '-dpng', '-r300');
fprintf('   -> Figura 3.1 salva com sucesso.\n\n');

%% ==============================================================================
%% NÍVEL 2: DINÂMICA SOB 5 DEGRAUS COM R_dist = 10 mOhm (SUBAMORTECIMENTO)
%% ==============================================================================
fprintf('3. Simulando 5 Degraus com R_dist = 10 mOhm vs 50 mOhm (C=1.5 mF, L=4 uH)...\n');

t_dyn = linspace(0, 15.0e-3, 30000); % 15 ms, 30k pontos para capturar dinâmica rápida
dt_dyn = t_dyn(2) - t_dyn(1);

P_IT_train = P_IT_idle * ones(size(t_dyn));
t_pulses = [1.0e-3, 3.5e-3, 6.0e-3, 8.5e-3, 11.0e-3]; 
t_width = 1.2e-3; 
t_rise = 10.0e-6; 

for p = 1:5
    t_start = t_pulses(p);
    t_end = t_start + t_width;
    for k = 1:length(t_dyn)
        t_now = t_dyn(k);
        if t_now >= t_start && t_now < t_start + t_rise
            P_IT_train(k) = P_IT_idle + (P_IT_nom - P_IT_idle) * ((t_now - t_start) / t_rise);
        elseif t_now >= t_start + t_rise && t_now <= t_end
            P_IT_train(k) = P_IT_nom;
        elseif t_now > t_end && t_now <= t_end + t_rise
            P_IT_train(k) = P_IT_nom - (P_IT_nom - P_IT_idle) * ((t_now - t_end) / t_rise);
        end
    end
end

P_bus_train = (P_IT_train / eta_conv_cons) + P_OISL_total + P_aux;

fig2 = figure(2, 'position', [100, 100, 1000, 700]);

% 3 Casos Comparativos:
% Caso 1: R_dist = 10 mOhms, C_bus = 1.5 mF, L_dist = 4.0 uH (Subamortecido zeta = 0.0968)
% Caso 2: R_dist = 50 mOhms, C_bus = 1.5 mF, L_dist = 4.0 uH (Resistivo zeta = 0.4841)
% Caso 3: R_dist = 14 mOhms, C_bus = 4.7 mF, L_dist = 1.8 uH (Nominal zeta = 0.3578)
C_comp = [1.5e-3, 1.5e-3, 4.7e-3];
L_comp = [4.0e-6, 4.0e-6, 1.8e-6];
R_comp = [0.010,  0.050,  0.014];

labels_comp = {'Pior Caso Subamortecido (C=1.5 mF, L=4 \muH, R=10 m\Omega, \zeta=0.097)', ...
               'Pior Caso Resistivo (C=1.5 mF, L=4 \muH, R=50 m\Omega, \zeta=0.484)', ...
               'Caso Nominal (C=4.7 mF, L=1.8 \muH, R=14 m\Omega, \zeta=0.358)'};
colors_comp = {'m', 'r', 'b'};

for c = 1:3
    Cb = C_comp(c);
    Lb = L_comp(c);
    Rb = R_comp(c);
    
    V_b = zeros(size(t_dyn));
    I_b = zeros(size(t_dyn));
    
    V_b(1) = 400.0;
    I_b(1) = P_bus_train(1) / V_b(1);
    V_src = V_b(1) + Rb * I_b(1);
    
    for k = 2:length(t_dyn)
        dI = (V_src - V_b(k-1) - Rb * I_b(k-1)) / Lb;
        dV = (I_b(k-1) - (P_bus_train(k-1) / max(V_b(k-1), 50.0))) / Cb;
        
        I_b(k) = I_b(k-1) + dI * dt_dyn;
        V_b(k) = V_b(k-1) + dV * dt_dyn;
    end
    
    subplot(2, 1, 1);
    plot(t_dyn * 1e3, V_b, 'linewidth', 1.8, 'color', colors_comp{c}); hold on;
    
    subplot(2, 1, 2);
    plot(t_dyn * 1e3, I_b, 'linewidth', 1.8, 'color', colors_comp{c}); hold on;
    
    undershoot = 400.0 - min(V_b);
    overshoot = max(V_b) - 400.0;
    fprintf('   -> %s:\n      V_min = %.2f V (Undershoot = %.2f V) | V_max = %.2f V (Overshoot = %.2f V)\n', ...
            labels_comp{c}, min(V_b), undershoot, max(V_b), overshoot);
end

subplot(2, 1, 1);
grid on;
title('Resposta Dinâmica a Trem de 5 Degraus de IA: Efeito do Subamortecimento (R_{dist} = 10 m\Omega, \zeta = 0.097)', 'fontsize', 11);
xlabel('Tempo (ms)', 'fontsize', 10);
ylabel('Tensão do Barramento (V)', 'fontsize', 10);
legend(labels_comp, 'location', 'southeast', 'fontsize', 8.5);
ylim([392, 405]);

subplot(2, 1, 2);
grid on;
title('Corrente Transitória sob Degraus de IA: Oscilações Ressonantes de Alta Frequência (f_0 \approx 2055 Hz)', 'fontsize', 11);
xlabel('Tempo (ms)', 'fontsize', 10);
ylabel('Corrente do Barramento (A)', 'fontsize', 10);
legend(labels_comp, 'location', 'southeast', 'fontsize', 8.5);
ylim([0, 80]);

print(fig2, 'figuras/figura_3_2_transiente_degrau_ia.png', '-dpng', '-r300');
fprintf('   -> Figura 3.2 salva com sucesso.\n\n');

%% ==============================================================================
%% NÍVEL 3: QUANTIFICAÇÃO DE RELIGAMENTOS AUTOMÁTICOS ATÉ RUPTURA DO SiC
%% ==============================================================================
fprintf('4. Simulando Religamentos Automáticos até violação do Limite I^2t = 1500 A^2.s...\n');

t_sc = linspace(0, 200.0e-6, 20000); 
dt_sc = t_sc(2) - t_sc(1);

% Cálculo de I^2t por disparo único para as seguintes configurações:
% 1. WCA 100 us com R_dist = 10 mOhm, R_sc = 1 mOhm
% 2. WCA 100 us com R_dist = 50 mOhm, R_sc = 1 mOhm
% 3. Nominal 15 us com R_dist = 10 mOhm, R_sc = 1 mOhm
% 4. Rápido 8 us com R_dist = 10 mOhm, R_sc = 1 mOhm

test_trips = [100.0e-6, 100.0e-6, 15.0e-6, 8.0e-6];
test_Rdist = [0.010,    0.050,    0.010,   0.010];
test_names = {'WCA 100 us (R = 10 mOhm)', ...
              'WCA 100 us (R = 50 mOhm)', ...
              'Nominal 15 us (R = 10 mOhm)', ...
              'Rapido 8 us (R = 10 mOhm)'};
test_colors = {'m', 'r', 'b', 'g'};

I2t_single = zeros(1, 4);
I_peak = zeros(1, 4);
E_arc_single = zeros(1, 4);

for idx = 1:4
    ttrip = test_trips(idx);
    Rdist = test_Rdist(idx);
    
    I_sc = zeros(size(t_sc));
    V_sc = zeros(size(t_sc));
    V_sc(1) = 400.0;
    I_sc(1) = 58.0;
    
    for k = 2:length(t_sc)
        if t_sc(k) < ttrip
            dI = (V_sc(k-1) - (Rdist + R_sc) * I_sc(k-1)) / L_dist_ext;
            dV = -I_sc(k-1) / C_bus_ext;
            I_sc(k) = I_sc(k-1) + dI * dt_sc;
            V_sc(k) = max(0.0, V_sc(k-1) + dV * dt_sc);
        else
            ext_time = 5.0e-6;
            if t_sc(k) <= ttrip + ext_time
                decay = 1.0 - (t_sc(k) - ttrip) / ext_time;
                I_sc(k) = I_sc(round(ttrip / dt_sc)) * decay;
            else
                I_sc(k) = 0.0;
            end
            V_sc(k) = 0.0;
        end
    end
    
    I_peak(idx) = max(I_sc);
    I2t_single(idx) = sum(I_sc.^2) * dt_sc;
    E_arc_single(idx) = sum(I_sc.^2 * R_sc) * dt_sc;
    
    fprintf('   -> %s: I_pico = %.1f A, I^2t único = %.2f A^2.s, E_arc único = %.4f J\n', ...
            test_names{idx}, I_peak(idx), I2t_single(idx), E_arc_single(idx));
end

% Vetor de religamentos automáticos (1 a 10 tentativas)
N_reclosures = 1:10;

fig3 = figure(3, 'position', [100, 100, 1000, 750]);

% Painel 1: Corrente transitória do curto-circuito franco
subplot(2, 1, 1);
for idx = 1:4
    ttrip = test_trips(idx);
    Rdist = test_Rdist(idx);
    I_sc = zeros(size(t_sc));
    V_sc = zeros(size(t_sc));
    V_sc(1) = 400.0;
    I_sc(1) = 58.0;
    for k = 2:length(t_sc)
        if t_sc(k) < ttrip
            dI = (V_sc(k-1) - (Rdist + R_sc) * I_sc(k-1)) / L_dist_ext;
            dV = -I_sc(k-1) / C_bus_ext;
            I_sc(k) = I_sc(k-1) + dI * dt_sc;
            V_sc(k) = max(0.0, V_sc(k-1) + dV * dt_sc);
        else
            ext_time = 5.0e-6;
            if t_sc(k) <= ttrip + ext_time
                decay = 1.0 - (t_sc(k) - ttrip) / ext_time;
                I_sc(k) = I_sc(round(ttrip / dt_sc)) * decay;
            else
                I_sc(k) = 0.0;
            end
            V_sc(k) = 0.0;
        end
    end
    plot(t_sc * 1e6, I_sc, 'linewidth', 1.8, 'color', test_colors{idx}); hold on;
end
grid on;
title('Corrente de Curto-Circuito Franco (R_{sc} = 1 mOhm) sob Diferentes Tempos de Corte e Resistencias', 'fontsize', 11);
xlabel('Tempo (us)', 'fontsize', 10);
ylabel('Corrente I_{sc} (A)', 'fontsize', 10);
legend(test_names, 'location', 'northeast', 'fontsize', 8.5);
xlim([0, 150]);

% Painel 2: Acúmulo de I^2t por religamento automático e superação do limite térmico
subplot(2, 1, 2);
for idx = 1:4
    I2t_accum = N_reclosures * I2t_single(idx);
    plot(N_reclosures, I2t_accum, '-o', 'linewidth', 2, 'color', test_colors{idx}, ...
         'markersize', 6, 'markerfacecolor', test_colors{idx}); hold on;
end

% Linha do limite térmico do semicondutor SiC
plot([1, 10], [I2t_limit_SiC, I2t_limit_SiC], 'k--', 'linewidth', 2.0);
text(1.2, I2t_limit_SiC + 450, 'Limite Termico SiC MOSFET (I^2t = 1500 A^2 s) - Queima Destrutiva', ...
     'color', 'k', 'fontsize', 9.5);

grid on;
title('Estresse Termico Cumulativo (I^2t) vs Numero de Religamentos Automaticos', 'fontsize', 11);
xlabel('Numero de Ciclos de Religamento Consecutivos (N)', 'fontsize', 10);
ylabel('Integral de Joule Acumulada I^2t (A^2 s)', 'fontsize', 10);
legend({'WCA 100 us (R_{dist} = 10 mOhm)', ...
        'WCA 100 us (R_{dist} = 50 mOhm)', ...
        'Nominal 15 us (R_{dist} = 10 mOhm)', ...
        'Rapido 8 us (R_{dist} = 10 mOhm)', ...
        'Limite Termico SiC (1500 A^2 s)'}, 'location', 'northwest', 'fontsize', 8.5);
xlim([1, 10]);
ylim([0, 12000]);

% Determinação exata do número de religamentos tolerados:
fprintf('\n--- QUANTIFICAÇÃO EXATA DE RELIGAMENTOS ATÉ FALHA TÉRMICA (I^2t > 1500 A^2.s) ---\n');
for idx = 1:4
    n_fail = ceil(I2t_limit_SiC / I2t_single(idx));
    n_allowed = n_fail - 1;
    fprintf('   -> %s:\n', test_names{idx});
    fprintf('      I^2t por disparo = %.2f A^2.s\n', I2t_single(idx));
    if n_allowed <= 0
        fprintf('      FALHA NO DISPARO INICIAL (0 religamentos tolerados!)\n');
    else
        fprintf('      Supera limite no ciclo %d (Tolerância máxima: %d religamento(s))\n', n_fail, n_allowed);
    end
end
fprintf('---------------------------------------------------------------------------------\n\n');

print(fig3, 'figuras/figura_3_3_curto_circuito_sspc.png', '-dpng', '-r300');
fprintf('   -> Figura 3.3 salva com sucesso.\n\n');

%% ==============================================================================
%% NÍVEL 4: ACOPLAMENTO TÉRMICO EM 5 CICLOS TÉRMICOS (2.5 HORAS)
%% ==============================================================================
fprintf('5. Acoplamento Eletrotérmico sob 5 CICLOS TÉRMICOS (2.5 h)...\n');

t_th_total = 5 * 1800.0;
t_th = linspace(0, t_th_total, 4500);
dt_th = t_th(2) - t_th(1);

P_loss_cons = P_IT_nom * ((1.0 - eta_conv_cons) / eta_conv_cons); 
P_loss_nom  = P_IT_nom * ((1.0 - eta_conv_nom) / eta_conv_nom);   

Q_chassis_cons = P_loss_cons + Q_OISL + 50.0; 
Q_chassis_nom  = P_loss_nom  + Q_OISL + 50.0; 

fig4 = figure(4, 'position', [100, 100, 1000, 700]);

for g = 1:3
    Gth = G_th_vec(g);
    
    T_cons = zeros(size(t_th));
    T_nom  = zeros(size(t_th));
    T_cons(1) = 25.0;
    T_nom(1)  = 25.0;
    
    for k = 2:length(t_th)
        dT_cons = (Q_chassis_cons - Gth * (T_cons(k-1) - (T_rad - 273.15))) / C_th_elec;
        dT_nom  = (Q_chassis_nom  - Gth * (T_nom(k-1)  - (T_rad - 273.15))) / C_th_elec;
        
        T_cons(k) = T_cons(k-1) + dT_cons * dt_th;
        T_nom(k)  = T_nom(k-1)  + dT_nom  * dt_th;
    end
    
    subplot(2, 1, 1);
    plot(t_th / 3600.0, T_cons, 'linewidth', 1.8, 'color', cores{g}); hold on;
    
    subplot(2, 1, 2);
    plot(t_th / 3600.0, T_nom, 'linewidth', 1.8, 'color', cores{g}); hold on;
end

subplot(2, 1, 1);
grid on;
title(sprintf('Evolucao Termica Continua em 5 Ciclos (2.5 h) - Cenario Conservador (\\eta_{conv} = 90%%, P_{loss} = %.0f W)', P_loss_cons), 'fontsize', 11);
xlabel('Tempo de Operação (horas)', 'fontsize', 10);
ylabel('Temperatura do Chassi (°C)', 'fontsize', 10);
plot([0, t_th_total/3600.0], [65, 65], 'r--', 'linewidth', 1.5);
text(0.1, 67, 'Limite Máximo de Confiabilidade Espacial (65 °C)', 'color', 'r', 'fontsize', 9);
legend({'G_{th} Mín = 80 W/K (Degradado)', 'G_{th} Nom = 120 W/K', 'G_{th} Máx = 160 W/K'}, ...
       'location', 'southeast', 'fontsize', 9);
ylim([20, 75]);

subplot(2, 1, 2);
grid on;
title(sprintf('Evolucao Termica Continua em 5 Ciclos (2.5 h) - Cenario Nominal GaN/SiC (\\eta_{conv} = 97.5%%, P_{loss} = %.0f W)', P_loss_nom), 'fontsize', 11);
xlabel('Tempo de Operação (horas)', 'fontsize', 10);
ylabel('Temperatura do Chassi (°C)', 'fontsize', 10);
plot([0, t_th_total/3600.0], [65, 65], 'r--', 'linewidth', 1.5);
text(0.1, 62, 'Limite Máximo de Confiabilidade Espacial (65 °C)', 'color', 'r', 'fontsize', 9);
legend({'G_{th} Mín = 80 W/K (Degradado)', 'G_{th} Nom = 120 W/K', 'G_{th} Máx = 160 W/K'}, ...
       'location', 'southeast', 'fontsize', 9);
ylim([20, 75]);

print(fig4, 'figuras/figura_3_4_temperatura_chassi_eficiencia.png', '-dpng', '-r300');
fprintf('   -> Figura 3.4 salva com sucesso.\n\n');

fprintf('=== TODAS AS SIMULAÇÕES FORAM CONCLUÍDAS COM SUCESSO! ===\n');
