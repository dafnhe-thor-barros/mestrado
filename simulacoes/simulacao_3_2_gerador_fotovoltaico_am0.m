% ==============================================================================
% SIMULACAO 3.2: MODELAGEM E ANALISE DE ESTRESSE DO GERADOR FOTOVOLTAICO MULTIJUNCAO AM0
% Projeto de Dissertacao de Mestrado - ProPECAut / IFES Campus Serra
% Autor: Dafnhe Thor Castro de Martins Barros
% Orientador: Prof. Dr. Flavio Barcelos Braz da Silva
% ==============================================================================
% Modelo conceitual parametrico in silico verificado numericamente.
% Nota Epistemologica: Ausencia de ensaios de bancada e termovacuo (TVAC);
% parametros classificados como [VERIFICAVEL], [CALCULADO] ou [PROVISORIO].
% ==============================================================================

clear; close all; clc;
fprintf('=== INICIANDO SIMULACAO 3.2: GERADOR FOTOVOLTAICO MULTIJUNCAO AM0 ===\n');

if ismember('qt', available_graphics_toolkits())
    graphics_toolkit('qt');
else
    graphics_toolkit('gnuplot');
end

% Diretorio de saida de figuras
fig_dir = 'figuras';
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

% ------------------------------------------------------------------------------
% 1. PARAMETROS FUNDAMENTAIS E DA CELULA MULTIJUNCAO (SPECTROLAB XTJ PRIME)
% ------------------------------------------------------------------------------
q  = 1.602176634e-19;    % Carga elementar [C]
kB = 1.380649e-23;       % Constante de Boltzmann [J/K]

% Condicoes de Referencia AM0:
% Reconciliacao: 1353 W/m^2 e a calibracao oficial do datasheet Spectrolab AM0;
% 1361 W/m^2 e a constante solar extraterrestre padrao ASTM E490 (1 UA).
T_ref = 301.15;          % Temperatura de referencia [K] (28 degC)
G_ref_cal = 1353.0;      % Irradiancia solar de calibracao Spectrolab AM0 [W/m^2] [VERIFICAVEL]
G_AM0 = 1361.0;              % Irradiancia solar extraterrestre ASTM E490 AM0 [W/m^2] [VERIFICAVEL]
G_datasheet = 1353.0;        % Condicao padrao de ensaio do datasheet Spectrolab XTJ Prime [W/m^2] [VERIFICAVEL]
scale_AM0 = G_AM0 / G_datasheet; % 1.00591 (fator de escala de irradiancia 1353 -> 1361 W/m^2) [CALCULADO]

% Dados eletricos de 1 celula Spectrolab XTJ Prime (AM0, 28 degC, G_datasheet = 1353 W/m^2):
% [PROVISORIO]: Formato geometrico e area ativa adotada dependem do corte e integracao
A_cel    = 30.18e-4;     % Area ativa da celula adotada [m^2] (30.18 cm^2) [PROVISORIO]
J_sc_ref = 18.00e-3;     % Densidade de curto-circuito [A/cm^2] (18.0 mA/cm^2) [VERIFICAVEL]
J_mp_ref = 17.50e-3;     % Densidade no MPP [A/cm^2] (17.5 mA/cm^2) [VERIFICAVEL]
V_oc_ref = 2.720;        % Tensao de circuito aberto [V] [VERIFICAVEL]
V_mp_ref = 2.406;        % Tensao no MPP [V] [VERIFICAVEL]
eta_ref  = 0.307;        % Eficiencia BOL a 28 degC AM0 (30.7%) [VERIFICAVEL]

I_sc_ref = J_sc_ref * (A_cel * 1e4); % Corrente Isc = 0.54324 A a 1353 W/m^2 [CALCULADO]
I_mp_ref = J_mp_ref * (A_cel * 1e4); % Corrente Imp = 0.52815 A a 1353 W/m^2 [CALCULADO]
P_mp_ref = V_mp_ref * I_mp_ref;      % Potencia de pico direta por celula: 1.2707 W [CALCULADO]

% Coeficientes termicos da celula (faixa valida de 15 degC a 125 degC no datasheet XTJ Prime):
% [VERIFICAVEL]: Datasheet Spectrolab XTJ Prime
dJsc_dT = +10.0e-6;      % dJsc/dT = +10 uA/(cm^2.degC) [VERIFICAVEL]
dJmp_dT = +5.0e-6;       % dJmp/dT = +5 uA/(cm^2.degC) [VERIFICAVEL]
dVoc_dT = -5.60e-3;      % dVoc/dT = -5.60 mV/degC [VERIFICAVEL]
dVmp_dT = -6.30e-3;      % dVmp/dT = -6.30 mV/degC [VERIFICAVEL]

K_I = dJsc_dT * (A_cel * 1e4); % +0.3018 mA/K para Isc [CALCULADO]
K_I_mp = dJmp_dT * (A_cel * 1e4); % +0.1509 mA/K para Imp [CALCULADO]
K_V = dVoc_dT;           % -5.60 mV/K para Voc [VERIFICAVEL]
K_V_mp = dVmp_dT;        % -6.30 mV/K para Vmp [VERIFICAVEL]

E_g = 1.42;              % Bandgap da subjuncao limitante InGaAs [eV]
E_g_J = E_g * q;         % Bandgap em Joules [J]

% Parametros do modelo concentrado de 1 diodo equivalente calibrado ao datasheet:
% [PROVISORIO / MODELO PARAMETRICO NUMERICO]:
a_eq = 3.00;             % Ideabilidade equivalente de 3 juncoes em serie calibrada [PROVISORIO]
R_s  = 0.080;            % Resistencia serie por celula [Ohm] [PROVISORIO]
R_sh = 1200.0;           % Resistencia shunt por celula [Ohm] [PROVISORIO]

% Corrente de saturacao de referencia I_0,ref consistente com Voc_ref
V_th_ref = (kB * T_ref) / q;
I_0_ref = (I_sc_ref - V_oc_ref / R_sh) / (exp(V_oc_ref / (a_eq * V_th_ref)) - 1.0);
fprintf('1. I_0,ref calculado analiticamente: %.4e A (a_eq = %.2f)\n', I_0_ref, a_eq);

% ------------------------------------------------------------------------------
% 2. CONFIGURACAO MATRICIAL: NO SATELITAL E CONJUNTO DA CONSTELACAO
% ------------------------------------------------------------------------------
% Arranjo em serie: N_s = 166 celulas gera V_mp_string = 166 * 2.406 = 399.4 V (~400 V)
N_s = 166;

% Strings em paralelo:
% Opcao nominal estrita de 50 kW uteis: N_p_no = 237 -> 39.342 celulas (49.99 kW BOL, 118.7 m^2)
% Opcao adotada com margem de cabeamento/mismatch (2.94% ~ 3%): N_p_no = 244 -> 40.504 celulas (51.47 kW BOL instalados, 122.2 m^2)
N_p_no = 244;                  % Strings em paralelo por no satelital [MARGEM INSTALADA DE 3%]
N_total_no = N_s * N_p_no;     % 40.504 celulas ativas por satelite [CALCULADO]
A_total_no = N_total_no * A_cel; % Area ativa do no = 122.24 m^2 [CALCULADO]
P_nom_BOL_no = N_total_no * P_mp_ref; % 51.47 kW instalados BOL (para ~50 kW uteis no barramento) [CALCULADO]

% Constelacao Agregada (20 nos ativos):
N_nos = 20;
N_p_const = N_p_no * N_nos;    % 4.880 strings em paralelo na constelacao [CALCULADO]
N_total_const = N_total_no * N_nos; % 810.080 celulas ativas [CALCULADO]
A_total_const = A_total_no * N_nos; % 2.444,8 m^2 [CALCULADO]
P_nom_BOL_const = P_nom_BOL_no * N_nos; % 1.029 MW instalados BOL (para 1.0 MW agregados) [CALCULADO]

fprintf('2. Configuracao No Satelital: %d strings de %d celulas (Total: %d celulas)\n', N_p_no, N_s, N_total_no);
fprintf('   Area ativa no: %.2f m^2 | Potencia instalada BOL (margem 3%%): %.2f kW (Nominal util: ~50 kW)\n', A_total_no, P_nom_BOL_no / 1e3);
fprintf('   Constelacao Agregada (20 nos): %d celulas | Area: %.2f m^2 | Potencia instalada: %.3f MW\n\n', N_total_const, A_total_const, P_nom_BOL_const / 1e6);

% ==============================================================================
% FIGURA 3.5: CURVAS CARACTERISTICAS I-V E P-V EM FUNCAO DA TEMPERATURA ORBITAL
% ==============================================================================
fprintf('3. Gerando Figura 3.5: Curvas I-V e P-V (modelo calibrado ao datasheet Spectrolab)...\n');
T_sweeps_C = [-150, -50, 28, 65, 85]; % [degC]
T_sweeps_K = T_sweeps_C + 273.15;
cores = {'b', 'c', 'g', 'm', 'r'};

fig5 = figure(1, 'position', [50, 50, 1000, 480]);

% Subplot 1: Curvas I-V do No Satelital
subplot(1, 2, 1);
hold on; grid on; box on;

i_norm = linspace(0, 0.999, 300);

for idx = 1:length(T_sweeps_K)
    T_curr = T_sweeps_K(idx);
    dT = T_curr - T_ref;
    V_th = (kB * T_curr) / q;
    
    % Coeficientes termicos aplicados rigorosamente:
    % (Faixa valida: 15 C a 125 C; T < 15 C representa extrapolacao analitica linear [PROVISORIO])
    Voc_cel = V_oc_ref + K_V * dT;
    Iph_cel = max(0.01, (I_sc_ref + K_I * dT) * scale_AM0);
    
    % Corrente de saturacao calibrada para garantir Voc(T) rigoroso
    I_0_curr = (Iph_cel - Voc_cel / R_sh) / (exp(Voc_cel / (a_eq * V_th)) - 1.0);
    
    % Corrente da celula
    i_cel = i_norm * Iph_cel;
    
    % Tensao da celula calculada explicitamente via relacao logaritmica com correcao de shunt
    v_j = a_eq * V_th * log(max(1.0, (Iph_cel - i_cel) / I_0_curr + 1.0));
    v_cel = v_j - i_cel * R_s - (v_j / R_sh) * R_s;
    v_cel = max(0, v_cel);
    
    % Escala da Matriz do No Satelital (Ns x Np_no)
    V_array = v_cel * N_s;
    I_array = i_cel * N_p_no;
    
    % Ordenar por tensao crescente para plotagem perfeita
    [V_sort, s_idx] = sort(V_array);
    I_sort = I_array(s_idx);
    
    plot(V_sort, I_sort, cores{idx}, 'linewidth', 2.0);
end

plot([400, 400], [0, 145], 'k--', 'linewidth', 1.8);
xlabel('Tensao da Matriz Solar (V)', 'fontsize', 10);
ylabel('Corrente da Matriz do No (A)', 'fontsize', 10);
title('Curvas I-V do No vs Temperatura (Np_{no}=244 strings, Ns=166)', 'fontsize', 11);
legend({'-150 C [Extrapolado]', '-50 C [Extrapolado]', '+28 C (BOL)', '+65 C (Nominal)', '+85 C (Max)', 'Barramento (400 V)'}, 'location', 'southwest', 'fontsize', 8);
xlim([0, 750]); ylim([0, 145]);

% Subplot 2: Curvas P-V da Matriz do No
subplot(1, 2, 2);
hold on; grid on; box on;

P_mpp_num_vec = zeros(length(T_sweeps_K), 1);
V_mpp_num_vec = zeros(length(T_sweeps_K), 1);

for idx = 1:length(T_sweeps_K)
    T_curr = T_sweeps_K(idx);
    dT = T_curr - T_ref;
    V_th = (kB * T_curr) / q;
    
    Voc_cel = V_oc_ref + K_V * dT;
    Iph_cel = max(0.01, (I_sc_ref + K_I * dT) * scale_AM0);
    I_0_curr = (Iph_cel - Voc_cel / R_sh) / (exp(Voc_cel / (a_eq * V_th)) - 1.0);
    
    i_cel = i_norm * Iph_cel;
    v_j = a_eq * V_th * log(max(1.0, (Iph_cel - i_cel) / I_0_curr + 1.0));
    v_cel = v_j - i_cel * R_s - (v_j / R_sh) * R_s;
    v_cel = max(0, v_cel);
    
    V_array = v_cel * N_s;
    I_array = i_cel * N_p_no;
    P_array = (V_array .* I_array) / 1e3; % em kW
    
    [V_sort, s_idx] = sort(V_array);
    P_sort = P_array(s_idx);
    
    [p_max, mpp_idx] = max(P_sort);
    v_mpp = V_sort(mpp_idx);
    P_mpp_num_vec(idx) = p_max;
    V_mpp_num_vec(idx) = v_mpp;
    
    plot(V_sort, P_sort, cores{idx}, 'linewidth', 2.0);
    plot(v_mpp, p_max, 'ko', 'markersize', 5);
end

P_mpp_num_28 = P_mpp_num_vec(3);
V_mpp_num_28 = V_mpp_num_vec(3);

plot([400, 400], [0, 75], 'k--', 'linewidth', 1.8);
plot([350, 350], [0, 75], 'k:', 'linewidth', 1.2);
plot([650, 650], [0, 75], 'k:', 'linewidth', 1.2);

xlabel('Tensao da Matriz Solar (V)', 'fontsize', 10);
ylabel('Potencia da Matriz do No (kW)', 'fontsize', 10);
title('Curvas P-V do No e Faixa MPPT Buck (350 a 650 V)', 'fontsize', 11);
legend({'-150 C [Extrapolado]', '-50 C [Extrapolado]', ...
        sprintf('+28 C (Catálogo: 51,5 kW | Modelo: %.2f kW)', P_mpp_num_28), ...
        sprintf('+65 C (Modelo: %.2f kW)', P_mpp_num_vec(4)), ...
        sprintf('+85 C (Modelo: %.2f kW)', P_mpp_num_vec(5)), ...
        'Barramento (400 V)', 'Faixa MPPT Buck'}, 'location', 'northwest', 'fontsize', 8);
xlim([0, 750]); ylim([0, 75]);

fprintf('   -> Potencia Nominal Direta de Catalogo (1353 W/m^2): %.2f kW\n', P_nom_BOL_no / 1e3);
fprintf('   -> Potencia de Pico do Modelo 1 Diodo a 28 C (1361 W/m^2): %.2f kW (Vmp = %.1f V)\n', P_mpp_num_28, V_mpp_num_28);

print(fig5, 'figuras/figura_3_5_curvas_iv_pv_temperatura.png', '-dpng', '-r300');
fprintf('   -> Figura 3.5 salva com sucesso.\n\n');

% ==============================================================================
% FIGURA 3.6: DEGRADACAO POR RADIACAO COSMICA (DDD 10 ANOS) E PERDAS AMBIENTAIS
% ==============================================================================
fprintf('4. Gerando Figura 3.6: Degradacao por radiacao e fatores ambientais...\n');
t_anos = linspace(0, 10, 100);

Phi_x = 1.0e14;
% Calibracao semi-empirica ajustada a tabela ECSS (ECSS-E-ST-20-08C) do datasheet XTJ Prime:
% Na fluencia de joelho Phi_x = 5x10^14 e-/cm^2 (1 MeV), P/P0 = 0.89 * 0.99 = 0.8811 (88.1%).
% Com 1 - C_d * log10(1 + 5) = 0.8811 -> C_d = 0.1189 / 0.77815 = 0.1528.
% NOTA EPISTEMICA: O ajuste foi calibrado pontualmente na tabela ECSS no ponto 5x10^14 e-/cm^2,
% constituindo modelo semi-empirico pontual e nao validacao experimental da curva completa.
C_d   = 0.1528; % Calibrado a tabela ECSS da Spectrolab [PROVISORIO / AJUSTE PONTUAL]

rate_500km  = 1.4e13;
rate_700km  = 2.5e13;
rate_1000km = 5.8e13;

loss_env_10y = 0.052; % 5.2% optica e erosao [PROVISORIO]

fig6 = figure(2, 'position', [50, 50, 1000, 450]);

subplot(1, 2, 1);
hold on; grid on; box on;
phi_500  = rate_500km * t_anos;
phi_700  = rate_700km * t_anos;
phi_1000 = rate_1000km * t_anos;

plot(t_anos, phi_500 / 1e14, 'b-', 'linewidth', 2.0);
plot(t_anos, phi_700 / 1e14, 'g-', 'linewidth', 2.2);
plot(t_anos, phi_1000 / 1e14, 'r-', 'linewidth', 2.0);
plot([0, 10], [1.0, 1.0], 'k--', 'linewidth', 1.5);

xlabel('Tempo de Missao (anos)', 'fontsize', 10);
ylabel('Fluencia Acumulada 1-MeV (x 10^{14} e^-/cm^2)', 'fontsize', 10);
title('Fluencia Acumulada de Radiacao (DDD)', 'fontsize', 11);
legend({'h = 500 km', 'h = 700 km (Nominal)', 'h = 1000 km (Extremo)', 'Fluencia Joelho Phi_x'}, 'location', 'northwest', 'fontsize', 8);
xlim([0, 10]);

subplot(1, 2, 2);
hold on; grid on; box on;

f_ddd_500  = 1.0 - C_d * log10(1.0 + phi_500 / Phi_x);
f_ddd_700  = 1.0 - C_d * log10(1.0 + phi_700 / Phi_x);
f_ddd_1000 = 1.0 - C_d * log10(1.0 + phi_1000 / Phi_x);

f_env = 1.0 - loss_env_10y * (t_anos / 10.0);
p_rel_500  = f_ddd_500 .* f_env * 100.0;
p_rel_700  = f_ddd_700 .* f_env * 100.0;
p_rel_1000 = f_ddd_1000 .* f_env * 100.0;

plot(t_anos, p_rel_500, 'b-', 'linewidth', 2.0);
plot(t_anos, p_rel_700, 'g-', 'linewidth', 2.4);
plot(t_anos, p_rel_1000, 'r-', 'linewidth', 2.0);
plot([0, 10], [88.0, 88.0], 'k--', 'linewidth', 1.8);

xlabel('Tempo de Missao (anos)', 'fontsize', 10);
ylabel('Retencao de Potencia P_{mp}(t) (%)', 'fontsize', 10);
title('Degradacao EOL em 10 Anos de Voo Espacial (Tabela ECSS Spectrolab)', 'fontsize', 11);
perda_500 = 100.0 - p_rel_500(end);
perda_700 = 100.0 - p_rel_700(end);
perda_1000 = 100.0 - p_rel_1000(end);
legend({sprintf('500 km (Perda %.1f%%)', perda_500), ...
        sprintf('700 km (Nominal: Perda %.1f%%)', perda_700), ...
        sprintf('1000 km (Perda %.1f%%)', perda_1000), ...
        'Requisito EOL (88% = 12% Perda Max)'}, 'location', 'southwest', 'fontsize', 8);
xlim([0, 10]); ylim([75, 102]);

print(fig6, 'figuras/figura_3_6_degradacao_radiacao_10anos.png', '-dpng', '-r300');
fprintf('   -> Figura 3.6 salva com sucesso.\n\n');

% ==============================================================================
% FIGURA 3.7: HOT SPOT SEVERO, POLARIZACAO REVERSA E RISCO DE ARCO NO VACUO
% ==============================================================================
fprintf('5. Gerando Figura 3.7: Analise de Hot Spot e Polarizacao Reversa...\n');

fig7 = figure(3, 'position', [50, 50, 1000, 450]);

subplot(1, 2, 1);
hold on; grid on; box on;

v_rev_vec = linspace(-35, 3, 300);

% Caso normal iluminado
i_normal = max(0, I_sc_ref - (v_rev_vec / R_sh) .* (v_rev_vec < 0) - (v_rev_vec > 2.5) * 0.5);
plot(v_rev_vec, i_normal, 'k-', 'linewidth', 1.8);

% Caso sombreado sem bypass (celula polarizada reversamente ate conducao de breakdown/shunt)
i_shaded_nobypass = - (v_rev_vec / R_sh) - 1e-4 * (abs(v_rev_vec).^2.2) .* (v_rev_vec < 0);
plot(v_rev_vec, i_shaded_nobypass, 'r-', 'linewidth', 2.2);

% Caso sombreado com diodo de bypass Schottky:
% Modelo Shockley do diodo Schottky acoplado para conduzir I_mp_ref (0.52815 A) em V_clamp = -0.350 V:
n_sch = 1.05;
V_clamp_target = 0.350; % V em conducao direta
Is_sch = I_mp_ref / (exp(V_clamp_target / (n_sch * V_th_ref)) - 1.0); % 1.3915 uA

v_rev_bypass = linspace(-1.0, 3, 300);
i_schottky = Is_sch * (exp(max(-30.0, -v_rev_bypass / (n_sch * V_th_ref))) - 1.0);
i_shaded_bypass = - (v_rev_bypass / R_sh) + i_schottky;
i_shaded_bypass = min(0.65, i_shaded_bypass);
plot(v_rev_bypass, i_shaded_bypass, 'b-', 'linewidth', 2.2);
plot([-35, 3], [I_mp_ref, I_mp_ref], 'm--', 'linewidth', 1.5);

xlabel('Tensao Terminal da Celula (V)', 'fontsize', 10);
ylabel('Corrente da Celula (A)', 'fontsize', 10);
title('Curva no Quadrante Reverso e Clamp do Bypass (-0.35 V)', 'fontsize', 11);
legend({'Celula Normal Iluminada', 'Sombreada SEM Bypass (Perigo)', 'Sombreada COM Bypass Schottky', sprintf('Corrente Forcada I_{mp} (%.3f A)', I_mp_ref)}, 'location', 'northwest', 'fontsize', 8);
xlim([-35, 3]); ylim([-0.2, 0.7]);

subplot(1, 2, 2);
hold on; grid on; box on;

t_sec = linspace(0, 180, 200);
C_th_cel = 3.8;
R_th_rad = 12.5;
T_amb_rad = 338.15; % 65 degC

% Dissipacao termica acoplada rigorosamente a mesma corrente I_mp_ref (0.52815 A) e tensao V_clamp (-0.350 V):
P_diss_nobypass = 28.0 * I_mp_ref;           % 28 V x 0.52815 A = 14.788 W [PROVISORIO]
P_diss_bypass   = V_clamp_target * I_mp_ref; % 0.350 V x 0.52815 A = 0.1849 W [CALCULADO]

tau_th = R_th_rad * C_th_cel;
T_nobypass_C = (T_amb_rad + P_diss_nobypass * R_th_rad * (1.0 - exp(-t_sec / tau_th))) - 273.15;
T_bypass_C   = (T_amb_rad + P_diss_bypass   * R_th_rad * (1.0 - exp(-t_sec / tau_th))) - 273.15;

plot(t_sec, T_nobypass_C, 'r-', 'linewidth', 2.5);
plot(t_sec, T_bypass_C, 'b-', 'linewidth', 2.2);
plot([0, 180], [180, 180], 'k--', 'linewidth', 1.5);
plot([0, 180], [230, 230], 'r--', 'linewidth', 1.8);

xlabel('Tempo sob Sombra/Defeito (s)', 'fontsize', 10);
ylabel('Temperatura Local da Celula (C)', 'fontsize', 10);
title('Evolucao Termica do Hot Spot no Vacuo (Estudo Parametrico)', 'fontsize', 11);
legend({'Sem Bypass (Destruicao Termica)', 'Com Diodo Bypass Ativo', 'Limite Termico Adesivo RTV (180 C)', 'Fusao da Solda de Prata (230 C)'}, 'location', 'northwest', 'fontsize', 8);
xlim([0, 180]); ylim([50, 260]);

print(fig7, 'figuras/figura_3_7_hotspot_polarizacao_reversa.png', '-dpng', '-r300');
fprintf('   -> Figura 3.7 salva com sucesso.\n\n');

% ==============================================================================
% FIGURA 3.8: TRANSIENTE DO TERMINADOR DE ECLIPSE E RASTREAMENTO MPPT
% ==============================================================================
fprintf('6. Gerando Figura 3.8: Resposta dinamica do MPPT no terminador...\n');

t_trans = linspace(0, 250, 400);
G_trans = zeros(size(t_trans));
for m = 1:length(t_trans)
    tt = t_trans(m);
    if tt < 20
        G_trans(m) = 0;
    elseif tt <= 35
        G_trans(m) = G_AM0 * (tt - 20) / 15;
    else
        G_trans(m) = G_AM0;
    end
end

tau_cel_heat = 85.0;
T_cel_trans_K = 173.15 + (338.15 - 173.15) * (1.0 - exp(-max(0, t_trans - 20) / tau_cel_heat));

V_mpp_teorico = zeros(size(t_trans));
P_mpp_teorico = zeros(size(t_trans));

for m = 1:length(t_trans)
    if G_trans(m) < 10
        V_mpp_teorico(m) = 0;
        P_mpp_teorico(m) = 0;
    else
        Tc = T_cel_trans_K(m);
        Vmp_c = V_mp_ref + K_V_mp * (Tc - T_ref);
        V_mpp_teorico(m) = N_s * Vmp_c;
        Imp_c = (I_mp_ref + K_I_mp * (Tc - T_ref)) * (G_trans(m) / G_datasheet);
        % Potencia no nivel do no satelital (N_p_no = 244 strings, margem 3%):
        P_mpp_teorico(m) = (V_mpp_teorico(m) * Imp_c * N_p_no) / 1e3;
    end
end

V_mppt_adapt = zeros(size(t_trans));
V_mppt_fixed = zeros(size(t_trans));
P_mppt_adapt = zeros(size(t_trans));
P_mppt_fixed = zeros(size(t_trans));

v_track_a = 400.0;
v_track_f = 400.0;

for m = 1:length(t_trans)
    if G_trans(m) < 10
        V_mppt_adapt(m) = 400.0;
        V_mppt_fixed(m) = 400.0;
        P_mppt_adapt(m) = 0.0;
        P_mppt_fixed(m) = 0.0;
    else
        v_target = V_mpp_teorico(m);
        v_track_a = v_track_a + 0.35 * (v_target - v_track_a);
        V_mppt_adapt(m) = v_track_a;
        P_mppt_adapt(m) = P_mpp_teorico(m) * (1.0 - 0.008 * (1.0 - exp(-(t_trans(m)-35)/10.0)));
        
        step_f = 2.5 * sign(v_target - v_track_f);
        v_track_f = v_track_f + step_f + 1.2 * sin(m * 0.8);
        V_mppt_fixed(m) = v_track_f;
        P_mppt_fixed(m) = P_mpp_teorico(m) * 0.965;
    end
end

fig8 = figure(4, 'position', [50, 50, 1000, 450]);

subplot(1, 2, 1);
hold on; grid on; box on;

plot(t_trans, V_mpp_teorico, 'k--', 'linewidth', 2.0);
plot(t_trans, V_mppt_adapt, 'g-', 'linewidth', 2.2);
plot(t_trans, V_mppt_fixed, 'r-', 'linewidth', 1.8);
plot([0, 250], [400, 400], 'b-', 'linewidth', 1.8);
plot([0, 250], [650, 650], 'm--', 'linewidth', 1.4);

xlabel('Tempo desde a Saida do Eclipse (s)', 'fontsize', 10);
ylabel('Tensao de Entrada do Conversor (V)', 'fontsize', 10);
title('Tensao no Choque Solar e Saida do Eclipse (Modelo Parametrico)', 'fontsize', 11);
legend({'V_{mp} Otimo Analitico', 'Rastreamento Parametrico Adaptativo', 'Rastreamento Degrau Fixo', 'Barramento (400 V)', 'Teto Conversor (650 V)'}, 'location', 'northeast', 'fontsize', 8);
xlim([0, 250]); ylim([300, 680]);

subplot(1, 2, 2);
hold on; grid on; box on;

plot(t_trans, P_mpp_teorico, 'k--', 'linewidth', 2.0);
plot(t_trans, P_mppt_adapt, 'g-', 'linewidth', 2.2);
plot(t_trans, P_mppt_fixed, 'r-', 'linewidth', 1.8);

xlabel('Tempo desde a Saida do Eclipse (s)', 'fontsize', 10);
ylabel('Potencia Entregue pela Matriz do No (kW)', 'fontsize', 10);
title('Potencia Capturada na Penumbra (Modelo Parametrico Sintetico)', 'fontsize', 11);
legend({'Disponivel Otima (No 51,5 kW)', 'Rastreamento Adaptativo (>99%)', 'Passo Fixo (~96.5%)'}, 'location', 'southeast', 'fontsize', 8);
xlim([0, 250]); ylim([0, 75]);

print(fig8, 'figuras/figura_3_8_transiente_mppt_penumbra.png', '-dpng', '-r300');
fprintf('   -> Figura 3.8 salva com sucesso.\n\n');

% ==============================================================================
% 5. SINTESE DOS PONTOS DE FALHA E RESULTADOS
% ==============================================================================
fprintf('=== RESUMO EXECUTIVO DOS ENSAIOS DE ESTRESSE E PONTOS DE FALHA (ITEM 3.2) ===\n');
fprintf('1. Compatibilidade com 400V: Strings de Ns = 166 celulas geram Vmp = 399.4 V a 28 C (BOL).\n');
fprintf('   A frio (-150 C), Voc atinge 617.0 V (extrapolacao analitica linear [PROVISORIO]). Exige conversor Buck com V_in,max >= 650 V.\n');
fprintf('2. Potencia e Area: 40.504 celulas (No: Np_no=244 strings, margem de 2.94%%) ocupam 122.2 m^2, gerando 51.47 kW BOL instalados (~50 kW uteis).\n');
fprintf('   Constelacao (20 nos: Np_const=4880): 810.080 celulas, 2.444,8 m^2, gerando 1.029 MW BOL instalados (1.0 MW agregados).\n');
fprintf('3. Degradacao em 10 Anos (700 km): Perda DDD de %.1f%% somada a perdas ambientais resulta em retencao de %.1f%% (Perda total %.1f%%).\n', ...
        (1.0 - f_ddd_700(end))*100, p_rel_700(end), perda_700);
fprintf('4. Ponto de Falha de Hot Spot: Sem diodo de bypass, a celula dissipa calor severo (%.2f W, T > 230 C);\n', P_diss_nobypass);
fprintf('   Com diodo Schottky, o clamp em -0.35 V limita a dissipacao a %.3f W e temperatura a 67 C.\n', P_diss_bypass);
fprintf('5. Ponto de Falha de Arco: Risco de arco sustentado no vacuo se V_string > 450 V em presenca de plasma.\n');
fprintf('   Mitigacao: Conversores MPPT distribuidos descentralizados limitando secoes a 40 kW.\n');
fprintf('=== SIMULACAO 3.2 CONCLUIDA COM SUCESSO ===\n');
