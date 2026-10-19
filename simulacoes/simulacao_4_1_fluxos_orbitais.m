% ==============================================================================
% SIMULAÇÃO E RENDERIZADOR 4.1 (CORRIGIDA E RECONCILIADA COM A DISSERTAÇÃO)
% FLUXOS RADIATIVOS AMBIENTAIS ORBITAIS (SSO, LEO E GEO) E BALANÇO TÉRMICO EOL
% ==============================================================================
% Programa de Pós-Graduação em Engenharia de Controle e Automação - ProPECAut / IFES
% Autor: Dafnhê Thor Castro de Martins Barros
% Auditoria Técnica: Dedução ECSS-E-HB-31-01, Ocultação Cônica e Análise de Falhas
% ==============================================================================
clear; clc; close all;
graphics_toolkit('gnuplot');
fprintf('==============================================================================\n');
fprintf('INICIANDO SIMULAÇÃO 4.1: FLUXOS AMBIENTAIS ORBITAIS E BALANÇO DE RADIADORES\n');
fprintf('==============================================================================\n\n');
% Diretório de destino
if ~exist('../figuras', 'dir')
    mkdir('../figuras');
end
%% 1. CONSTANTES FUNDAMENTAIS E PARÂMETROS AMBIENTAIS DE PROJETO
sigma_SB = 5.670374e-8;      % [W/(m^2 K^4)] Constante fundamental de Stefan-Boltzmann
G_sc     = 1361.0;           % [W/m^2] Constante solar extraterrestre AM0 a 1 UA
q_Earth  = 237.0;            % [W/m^2] Emissão térmica infravermelha média da Terra
albedo_E = 0.30;             % Coeficiente de albedo planetário médio terrestre
R_E      = 6371.0e3;         % [m] Raio volumétrico médio da Terra
R_Sun    = 696340.0e3;       % [m] Raio do Sol
d_ES     = 1.496e11;         % [m] Distância média Terra-Sol (1 UA)
% Propriedades Termo-Ópticas em Fim de Vida Útil (EOL 10 Anos - Tabela 10)
alpha_s_eol  = 0.24;         % Absorvidade solar OSR degradada (UV e contaminação)
eps_ir_eol   = 0.82;         % Emissividade infravermelha OSR degradada
A_rad_adot   = 88.6;         % [m^2] Área efetiva adotada no nó satelital
eta_fin      = 0.945;        % Eficiência da aleta composta APG/K13D2U (Seção 4.3.1)
Q_int_diss   = 41983.84;     % [W] Calor interno total dissipado no nó (Seção 4.1.7)
Delta_T_cond = 25.0;         % [K] Salto térmico condutivo junção-radiador (LHP)
% Limiares Operacionais Rígidos da Dissertação (Tabela 13)
T_rad_max_lim = 55.0;        % [°C] Limite operacional máximo de projeto do radiador
T_j_target    = 75.0;        % [°C] Temperatura alvo em regime normal
T_j_throt_lim = 80.0;        % [°C] Limite de início de thermal throttling
T_j_trip_lim  = 85.0;        % [°C] Corte térmico de emergência (Thermal Trip)
%% 2. MODELAGEM ANALÍTICA DOS FATORES DE VISÃO PLANETÁRIOS (NORMA ECSS-E-HB-31-01)
% Função analítica fechada para placa plana vertical (phi = 90° em relação ao Nadir)
% Equação 115 do Apêndice B.6 da dissertação
% Fatores de visão planetários tabelados pela norma ECSS-E-HB-31-01
h_sso = 700.0e3;             % [m] SSO 700 km
h_leo = 500.0e3;             % [m] LEO 500 km
h_geo = 35786.0e3;           % [m] GEO 35.786 km
F_rad_E_sso = 0.2328;        % Fator de visão ECSS para radiador anti-solar em SSO 700 km
F_rad_E_leo = 0.4473;        % Fator de visão ECSS em LEO 500 km
F_rad_E_geo = 0.0035;        % Fator de visão ECSS em GEO 35.786 km
%% 3. BALANÇO DE FLUXOS ABSORVIDOS E TEMPERATURAS DE EQUILÍBRIO (EOL)
% --- Órbita SSO Dawn-Dusk (700 km, sem eclipse, atitude anti-solar) ---
q_alb_sso = alpha_s_eol * G_sc * albedo_E * F_rad_E_sso * 0.32; % 7.25 W/m^2 médio
q_ir_sso  = eps_ir_eol * q_Earth * F_rad_E_sso;                  % 45.23 W/m^2
q_amb_sso = q_alb_sso + q_ir_sso;                                % 52.48 W/m^2
% --- Órbita LEO Convencional (500 km) ---
q_alb_leo = alpha_s_eol * G_sc * albedo_E * F_rad_E_leo * 0.45; % 19.45 W/m^2 médio
q_ir_leo  = eps_ir_eol * q_Earth * F_rad_E_leo;                  % 78.50 W/m^2
q_amb_leo = q_alb_leo + q_ir_leo;                                % 97.95 W/m^2
% --- Órbita GEO Equatorial (35.786 km) ---
q_alb_geo = alpha_s_eol * G_sc * albedo_E * F_rad_E_geo * 0.10; % 0.03 W/m^2
q_ir_geo  = eps_ir_eol * q_Earth * F_rad_E_geo;                  % 0.68 W/m^2
q_amb_geo = q_alb_geo + q_ir_geo;                                % 0.71 W/m^2
% Cálculo das Temperaturas de Equilíbrio Térmico (A_rad = 88.6 m^2)
calc_T_rad_C = @(q_amb, A) (((Q_int_diss / A + q_amb) / (eta_fin * eps_ir_eol * sigma_SB))^0.25) - 273.15;
T_rad_sso = calc_T_rad_C(q_amb_sso, A_rad_adot);
T_j_sso   = T_rad_sso + Delta_T_cond;
T_rad_leo = calc_T_rad_C(q_amb_leo, A_rad_adot);
T_j_leo   = T_rad_leo + Delta_T_cond;
T_rad_geo = calc_T_rad_C(q_amb_geo, A_rad_adot);
T_j_geo   = T_rad_geo + Delta_T_cond;
% Capacidade Líquida Máxima a 55 °C e Área de Equilíbrio Requerida
q_emit_55C = eta_fin * eps_ir_eol * sigma_SB * ((55.0 + 273.15)^4); % 509.50 W/m^2
q_liq_sso  = q_emit_55C - q_amb_sso; % 457.02 W/m^2
q_liq_leo  = q_emit_55C - q_amb_leo; % 411.55 W/m^2
q_liq_geo  = q_emit_55C - q_amb_geo; % 508.79 W/m^2
A_rad_req_sso = Q_int_diss / q_liq_sso; % 91.86 m^2
A_rad_req_leo = Q_int_diss / q_liq_leo; % 102.02 m^2
A_rad_req_geo = Q_int_diss / q_liq_geo; % 82.52 m^2
fprintf('--- 3.1 BALANÇO RADIATIVO EOL (A_rad = 88.6 m^2) ---\n');
fprintf('SSO 700 km: Fluxo Parasita = %5.2f W/m^2 | T_rad = %5.2f °C | T_j = %5.2f °C\n', ...
    q_amb_sso, T_rad_sso, T_j_sso);
if T_j_sso >= T_j_throt_lim
    fprintf('  -> Status: [FALHA DE DIMENSIONAMENTO] Excede teto de throttling (%5.2f °C >= 80.0 °C)!\n', T_j_sso);
end
fprintf('  -> Área Requerida para manter T_rad = 55 °C em SSO: %5.2f m^2 (Déficit: %4.2f m^2)\n\n', ...
    A_rad_req_sso, A_rad_req_sso - A_rad_adot);
fprintf('LEO 500 km: Fluxo Parasita = %5.2f W/m^2 | T_rad = %5.2f °C | T_j = %5.2f °C\n', ...
    q_amb_leo, T_rad_leo, T_j_leo);
fprintf('  -> Área Requerida para manter T_rad = 55 °C em LEO: %5.2f m^2\n\n', A_rad_req_leo);
fprintf('GEO 35786 km: Fluxo Parasita = %5.2f W/m^2 | T_rad = %5.2f °C | T_j = %5.2f °C\n', ...
    q_amb_geo, T_rad_geo, T_j_geo);
fprintf('  -> Área Requerida para manter T_rad = 55 °C em GEO: %5.2f m^2 (Superdimensionado)\n\n', A_rad_req_geo);
%% 4. SENSIBILIDADE A DESVIOS ANGULARES DE ATITUDE (VAZAMENTO SOLAR DIRETO)
delta_theta_vec = 0:0.1:10; % [graus]
T_j_sens_sso = zeros(size(delta_theta_vec));
T_rad_sens_sso = zeros(size(delta_theta_vec));
for idx = 1:length(delta_theta_vec)
    d_th = delta_theta_vec(idx);
    % Radiação solar direta incidente absorvida por desapontamento anti-solar
    q_sol_pert = alpha_s_eol * G_sc * sind(d_th);
    q_tot_pert = q_amb_sso + q_sol_pert;
    
    t_rad_k = calc_T_rad_C(q_tot_pert, A_rad_adot);
    T_rad_sens_sso(idx) = t_rad_k;
    T_j_sens_sso(idx)   = t_rad_k + Delta_T_cond;
end
% Localização dos pontos de bifurcação e cruzamento
idx_throt = find(T_j_sens_sso >= T_j_throt_lim, 1);
idx_trip  = find(T_j_sens_sso >= T_j_trip_lim, 1);
if ~isempty(idx_trip)
    theta_trip = delta_theta_vec(idx_trip);
    trip_str = sprintf('%4.2f graus', theta_trip);
else
    trip_str = '> 10.0 graus (Margem Segura)';
end
fprintf('--- 4.1 SENSIBILIDADE TÉRMICA AO CONTROLE DE ATITUDE (ADCS) ---\n');
fprintf('Ponto Nominal (0.0°): T_j = %5.2f °C\n', T_j_sens_sso(1));
fprintf('Disparo de Thermal Trip (Corte a 85 °C): Ocorre em erro de atitude de %s!\n\n', trip_str);
%% 5. MODELO DE OCULTAÇÃO CÔNICA SUAVE (EQUAÇÕES 45 E 114 DA DISSERTAÇÃO)
% Propagação temporal de uma órbita LEO para fator de iluminação S(t) contínuo
t_orb_leo = linspace(0, 5677, 500); % Período LEO: 94.6 min
S_conico_leo = zeros(size(t_orb_leo));
r_sat_leo = R_E + h_leo;
alpha_Earth = asin(R_E / r_sat_leo); % Semi-ângulo aparente da Terra
alpha_Solar = asin(R_Sun / d_ES);    % Semi-ângulo aparente do Sol
for k = 1:length(t_orb_leo)
    % Posição angular em órbita coplanar com eclipse central
    theta_orb = (2 * pi / 5677) * t_orb_leo(k);
    gamma_sep = acos(cos(theta_orb)); % Separação angular aparente Sol-Terra
    
    if gamma_sep >= (alpha_Earth + alpha_Solar)
        S_conico_leo(k) = 1.0; % Iluminação Plena
    elseif gamma_sep <= (alpha_Earth - alpha_Solar)
        S_conico_leo(k) = 0.0; % Umbra Total
    else
        % Penumbra gradual contínua (Equação 114 do Apêndice B.5)
        d_val = gamma_sep;
        r1 = alpha_Solar; r2 = alpha_Earth;
        term1 = r1^2 * acos((d_val^2 + r1^2 - r2^2) / (2 * d_val * r1));
        term2 = r2^2 * acos((d_val^2 + r2^2 - r1^2) / (2 * d_val * r2));
        term3 = 0.5 * sqrt(max(0, (-d_val + r1 + r2) * (d_val + r1 - r2) * (d_val - r1 + r2) * (d_val + r1 + r2)));
        A_overlap = term1 + term2 - term3;
        S_conico_leo(k) = max(0.0, min(1.0, 1.0 - (A_overlap / (pi * r1^2))));
    end
end
%% ==============================================================================
%% RENDERIZAÇÃO DA FIGURA 4.1: MODELO RADIATIVO E FLUXOS COMPARATIVOS
%% ==============================================================================
fig_4_1 = figure('visible', 'off', 'position', [100, 100, 1050, 800]);
graphics_toolkit(fig_4_1, 'gnuplot');
% % --- (a1) Fator de Iluminação Solar S(t) (Ocultação Cônica Exata) ---
subplot(2, 2, 1);
ha1 = plot(t_orb_leo / 60, S_conico_leo, 'r-', 'linewidth', 2.0); hold on;
ha2 = plot([0, max(t_orb_leo)/60], [1.0, 1.0], 'b-', 'linewidth', 1.8);
grid on;
xlim([0, max(t_orb_leo)/60]);
ylim([-0.05, 1.10]);
xlabel('Tempo Orbital (Minutos)', 'fontsize', 10, 'fontweight', 'bold');
ylabel('Fator de Iluminação Solar S(t) [-]', 'fontsize', 10, 'fontweight', 'bold');
title('(a1) Ocultação Cônica Suave (Umbra e Penumbra)', 'fontsize', 11, 'fontweight', 'bold');
legend([ha1, ha2], {'LEO Convencional (Eclipse Cônico 35.8 min)', 'SSO Dawn-Dusk (\beta = 85°, Solstício Contínuo)'}, ...
       'location', 'southwest', 'fontsize', 8);
% --- (a2) Decomposição dos Fluxos Parasitas Absorvidos em EOL ---
subplot(2, 2, 2);
flux_matrix = [q_alb_sso, q_ir_sso, q_amb_sso; ...
               q_alb_leo, q_ir_leo, q_amb_leo; ...
               q_alb_geo, q_ir_geo, q_amb_geo];
b = bar(1:3, flux_matrix, 0.65);
set(b(1), 'facecolor', [0.2, 0.6, 0.3]);
set(b(2), 'facecolor', [0.85, 0.35, 0.2]);
set(b(3), 'facecolor', [0.2, 0.4, 0.8]);
set(gca, 'xtick', 1:3, 'xticklabel', {'SSO 700 km', 'LEO 500 km', 'GEO 35786 km'}, 'fontsize', 9);
grid on;
ylim([0, 115]);
ylabel('Densidade de Fluxo Absorvido (W/m^2)', 'fontsize', 10, 'fontweight', 'bold');
title('(a2) Fluxos Parasitas em Fim de Vida Útil (EOL)', 'fontsize', 11, 'fontweight', 'bold');
legend(b, {'Albedo Terrestre', 'Infravermelho Planetário', 'Fluxo Total Médio'}, ...
       'location', 'northwest', 'fontsize', 8);
% --- (b1) Sensibilidade Térmica ao Apontamento de Atitude ---
subplot(2, 2, 3);
p1 = plot(delta_theta_vec, T_j_sens_sso, 'b-', 'linewidth', 2.0); hold on;
p2 = plot(delta_theta_vec, T_rad_sens_sso, 'm--', 'linewidth', 1.6);
% Linhas guias dos limiares de controle TEMS
plot([0, 10], [T_j_target, T_j_target], 'g:', 'linewidth', 1.4);
plot([0, 10], [T_j_throt_lim, T_j_throt_lim], 'k--', 'linewidth', 1.4);
plot([0, 10], [T_j_trip_lim, T_j_trip_lim], 'r--', 'linewidth', 1.8);
plot([0, 10], [T_rad_max_lim, T_rad_max_lim], 'm:', 'linewidth', 1.4);
text(4.0, T_j_target + 1.5, 'Alvo Normal (75 °C)', 'fontsize', 8, 'color', [0, 0.5, 0], 'fontweight', 'bold');
text(4.0, T_j_throt_lim - 1.8, 'Início Throttling (80 °C)', 'fontsize', 8, 'fontweight', 'bold');
text(4.0, T_j_trip_lim + 1.5, 'Corte Trip (85 °C)', 'fontsize', 8, 'color', [0.8, 0, 0], 'fontweight', 'bold');
text(0.2, T_rad_max_lim + 1.5, 'T_{rad,max} (55 °C)', 'fontsize', 8, 'color', [0.6, 0, 0.6], 'fontweight', 'bold');
grid on;
xlim([0, 10]);
ylim([50, 100]);
xlabel('Desvio de Atitude do Radiador |\delta\theta_{att}| (Graus)', 'fontsize', 10, 'fontweight', 'bold');
ylabel('Temperatura de Equilíbrio (°C)', 'fontsize', 10, 'fontweight', 'bold');
title('(b1) Impacto Térmico de Erro de Apontamento ADCS', 'fontsize', 11, 'fontweight', 'bold');
legend([p1, p2], {'Temperatura de Junção T_j (Salto 25 K)', 'Temperatura do Radiador T_{rad}'}, ...
       'location', 'southeast', 'fontsize', 8);
% --- (b2) Dimensionamento de Área vs. Capacidade Líquida em EOL ---
subplot(2, 2, 4);
areas_eval = 75:1:115;
q_dem_eval = Q_int_diss ./ areas_eval;
h1 = plot(areas_eval, q_dem_eval, 'b-', 'linewidth', 2.0); hold on;
% Capacidades líquidas a 55 °C em EOL
h2 = plot([75, 115], [q_liq_sso, q_liq_sso], 'r--', 'linewidth', 1.6);
h3 = plot([75, 115], [q_liq_leo, q_liq_leo], 'm:',  'linewidth', 1.6);
% Marcação das áreas
plot([A_rad_adot, A_rad_adot], [350, 560], 'k:', 'linewidth', 1.6);
plot([A_rad_req_sso, A_rad_req_sso], [350, 560], 'g--', 'linewidth', 1.8);
text(A_rad_adot + 0.5, 520, sprintf('Adotada: %4.1f m^2 (T_j = 82.7 °C)', A_rad_adot), ...
    'fontsize', 8, 'fontweight', 'bold', 'color', [0.8, 0, 0]);
text(A_rad_req_sso + 0.5, 470, sprintf('Requerida SSO: %4.1f m^2 (T_{rad} = 55 °C)', A_rad_req_sso), ...
    'fontsize', 8, 'fontweight', 'bold', 'color', [0, 0.5, 0]);
grid on;
xlim([75, 115]);
ylim([360, 560]);
xlabel('Área Efetiva de Radiadores A_{rad} (m^2)', 'fontsize', 10, 'fontweight', 'bold');
ylabel('Densidade de Fluxo Térmico (W/m^2)', 'fontsize', 10, 'fontweight', 'bold');
title('(b2) Balanço de Dimensionamento Térmico no Vácuo', 'fontsize', 11, 'fontweight', 'bold');
legend([h1, h2, h3], {'Demanda Térmica Q_{int}/A_{rad}', 'Capacidade Líquida SSO (457 W/m^2)', ...
       'Capacidade Líquida LEO (411 W/m^2)'}, 'location', 'northeast', 'fontsize', 8);
print(fig_4_1, '../figuras/figura_4_1_fluxos_orbitais_comparativos.png', '-dpng', '-r300');
fprintf('Figura 4.1 exportada com sucesso: ../figuras/figura_4_1_fluxos_orbitais_comparativos.png\n');
fprintf('==============================================================================\n');
fprintf('SIMULAÇÃO 4.1 CONCLUÍDA COM RIGOR AEROESPACIAL.\n');
fprintf('==============================================================================\n');
