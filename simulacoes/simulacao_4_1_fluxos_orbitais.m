% ==============================================================================
% SIMULACAO 4.1: DINAMICA ORBITAL E FLUXOS RADIATIVOS COMPARATIVOS
% Projeto de Dissertacao de Mestrado - ProPECAut / IFES Campus Serra
% ==============================================================================
% Modela a propagacao kepleriana e a geometria de iluminacao solar com:
% - Ocultacao conica suave (Umbra, Penumbra e Iluminacao Plena)
% - Fator de visao satelite-Terra exato (ECSS-E-HB-31-01 Parte 1)
% - Comparacao entre regimes orbitais: SSO Dawn-Dusk (700 km), LEO (500 km) e GEO (35.786 km)
% - Analise de sensibilidade e ponto de falha sob perturbacao angular de atitude (delta_theta)
% ==============================================================================

clear; close all; clc;
fprintf('=== EXECUTANDO SIMULACAO 4.1: DINAMICA ORBITAL E FLUXOS RADIATIVOS ===\n');

% ------------------------------------------------------------------------------
% 1. CONSTANTES FUNDAMENTAIS E PARAMETROS PLANETARIOS
% ------------------------------------------------------------------------------
G_AM0 = 1361.0;            % Irradiancia solar extraterrestre [W/m^2]
q_Earth_nom = 237.0;       % Emissao infravermelha terrestre media [W/m^2]
a_Earth_nom = 0.30;        % Albedo medio planetario [-]
R_E = 6371.0;              % Raio volumetrico medio da Terra [km]
R_sun = 696340.0;          % Raio medio do Sol [km]
d_ES = 1.495978707e8;      % Distancia media Terra-Sol (1 UA) [km]
mu_E = 398600.4418;        % Parametro gravitacional da Terra [km^3/s^2]
sigma_SB = 5.670374e-8;    % Constante de Stefan-Boltzmann [W/(m^2 K^4)]
T_space = 3.0;             % Fundo cosmico [K]

% Propriedades termo-opticas em Fim de Vida (EOL - 10 anos)
alpha_s_rad = 0.14;        % OSR Quartzo/Prata EOL
eps_ir_rad  = 0.88;        % Emissividade OSR
alpha_s_sc  = 0.26;        % Tinta branca AZ-93 EOL
eps_ir_sc   = 0.89;        % Emissividade AZ-93
eta_fin     = 0.945;       % Eficiencia da aleta de APG composito

% ------------------------------------------------------------------------------
% 2. CONFIGURACAO DOS REGIMES ORBITAIS
% ------------------------------------------------------------------------------
altitudes = [700.0, 500.0, 35786.0]; % km
nomes_orbitas = {'SSO Dawn-Dusk (700 km)', 'LEO Convencional (500 km)', 'GEO (35.786 km)'};

T_orb = zeros(3, 1);
for i = 1:3
    a_orb = R_E + altitudes(i);
    T_orb(i) = 2 * pi * sqrt(a_orb^3 / mu_E); % segundos
end
fprintf('Periodos Orbitais: SSO = %.2f min | LEO = %.2f min | GEO = %.2f h\n', ...
        T_orb(1)/60, T_orb(2)/60, T_orb(3)/3600);

% ------------------------------------------------------------------------------
% 3. SIMULACAO DE 1 ORBITA COMPLETA PARA SSO E LEO
% ------------------------------------------------------------------------------
N_pts = 1000;

% SSO Dawn-Dusk:
t_sso = linspace(0, T_orb(1), N_pts)';
theta_orb_sso = 2 * pi * t_sso / T_orb(1);

H_sso = (R_E + altitudes(1)) / R_E;
theta_m_sso = asin(1 / H_sso);
phi_rad_sso = pi / 2; % 90 deg (anti-solar, perpendicular ao Nadir)
term1 = 0.5 - (1/pi)*asin(sqrt(H_sso^2 - 1)/(H_sso * sin(phi_rad_sso)));
term2 = (cos(phi_rad_sso)/(pi*H_sso^2)) * acos(sqrt(H_sso^2 - 1)*cot(phi_rad_sso));
term3 = (sqrt(H_sso^2 - 1)/(pi*H_sso^2)) * sqrt(1 - H_sso^2 * cos(phi_rad_sso)^2);
F_rad_E_sso = term1 + term2 - term3;
fprintf('Fator de visao Radiador-Terra SSO Dawn-Dusk (phi = 90 deg): F_rad_E = %.4f\n', F_rad_E_sso);

S_sso = ones(N_pts, 1);
gamma_sso = acos(cos(85 * pi/180) * cos(theta_orb_sso));
alpha_sun = asin(R_sun / d_ES);
alpha_E_sso = asin(R_E / (R_E + altitudes(1)));
for k = 1:N_pts
    gam = gamma_sso(k);
    if gam >= (alpha_E_sso + alpha_sun)
        S_sso(k) = 1.0;
    elseif gam <= (alpha_E_sso - alpha_sun)
        S_sso(k) = 0.0;
    else
        A_oc = alpha_sun^2 * acos((gam^2 + alpha_sun^2 - alpha_E_sso^2)/(2*gam*alpha_sun)) + ...
               alpha_E_sso^2 * acos((gam^2 + alpha_E_sso^2 - alpha_sun^2)/(2*gam*alpha_E_sso)) - ...
               0.5 * sqrt(max(0, ((alpha_sun + alpha_E_sso)^2 - gam^2)*(gam^2 - (alpha_E_sso - alpha_sun)^2)));
        S_sso(k) = max(0.0, 1.0 - A_oc / (pi * alpha_sun^2));
    end
end

q_sol_abs_sso = zeros(N_pts, 1);
cos_zenith_sso = max(0.0, sin(85*pi/180)*cos(theta_orb_sso));
q_alb_abs_sso = alpha_s_rad * G_AM0 * a_Earth_nom * F_rad_E_sso .* cos_zenith_sso .* S_sso;
q_ir_abs_sso = eps_ir_rad * q_Earth_nom * F_rad_E_sso * ones(N_pts, 1);
q_parasita_sso = q_sol_abs_sso + q_alb_abs_sso + q_ir_abs_sso;

% LEO Convencional (500 km, beta = 20 deg):
t_leo = linspace(0, T_orb(2), N_pts)';
theta_orb_leo = 2 * pi * t_leo / T_orb(2);
H_leo = (R_E + altitudes(2)) / R_E;
theta_m_leo = asin(1 / H_leo);
phi_rad_leo = pi / 3; % 60 deg
term1_l = 0.5 - (1/pi)*asin(sqrt(H_leo^2 - 1)/(H_leo * sin(phi_rad_leo)));
term2_l = (cos(phi_rad_leo)/(pi*H_leo^2)) * acos(sqrt(H_leo^2 - 1)*cot(phi_rad_leo));
term3_l = (sqrt(H_leo^2 - 1)/(pi*H_leo^2)) * sqrt(max(0, 1 - H_leo^2 * cos(phi_rad_leo)^2));
F_rad_E_leo = term1_l + term2_l - term3_l;
fprintf('Fator de visao Radiador-Terra LEO Convencional (phi = 60 deg): F_rad_E = %.4f\n', F_rad_E_leo);

S_leo = ones(N_pts, 1);
gamma_leo = acos(cos(20 * pi/180) * cos(theta_orb_leo));
alpha_E_leo = asin(R_E / (R_E + altitudes(2)));
for k = 1:N_pts
    gam = gamma_leo(k);
    if gam >= (alpha_E_leo + alpha_sun)
        S_leo(k) = 1.0;
    elseif gam <= (alpha_E_leo - alpha_sun)
        S_leo(k) = 0.0;
    else
        A_oc = alpha_sun^2 * acos((gam^2 + alpha_sun^2 - alpha_E_leo^2)/(2*gam*alpha_sun)) + ...
               alpha_E_leo^2 * acos((gam^2 + alpha_E_leo^2 - alpha_sun^2)/(2*gam*alpha_E_leo)) - ...
               0.5 * sqrt(max(0, ((alpha_sun + alpha_E_leo)^2 - gam^2)*(gam^2 - (alpha_E_leo - alpha_sun)^2)));
        S_leo(k) = max(0.0, 1.0 - A_oc / (pi * alpha_sun^2));
    end
end

cos_inc_leo = max(0.0, cos(theta_orb_leo - pi/4));
q_sol_abs_leo = alpha_s_rad * G_AM0 * cos_inc_leo .* S_leo * 0.3;
cos_zenith_leo = max(0.0, cos(theta_orb_leo));
q_alb_abs_leo = alpha_s_rad * G_AM0 * a_Earth_nom * F_rad_E_leo .* cos_zenith_leo .* S_leo;
q_ir_abs_leo = eps_ir_rad * q_Earth_nom * F_rad_E_leo * ones(N_pts, 1);
q_parasita_leo = q_sol_abs_leo + q_alb_abs_leo + q_ir_abs_leo;

% ------------------------------------------------------------------------------
% 4. ENSAIO DE ESTRESSE E PONTO DE FALHA: DESVIO ANGULAR DE ATITUDE
% ------------------------------------------------------------------------------
fprintf('2. Analisando ponto de falha sob desvio angular de atitude (0 a 30 deg)...\n');
delta_theta_deg = linspace(0, 30, 200)';
delta_theta_rad = delta_theta_deg * pi / 180;

q_sol_pert_bol = 0.08 * G_AM0 * sin(delta_theta_rad);
q_sol_pert_eol = alpha_s_rad * G_AM0 * sin(delta_theta_rad);

T_rad_nom = 318.15;
T_sink_sso = ((1 - F_rad_E_sso)*T_space^4 + F_rad_E_sso*(q_Earth_nom/sigma_SB))^0.25;
q_emit_bruta = eps_ir_rad * sigma_SB * (T_rad_nom^4 - T_sink_sso^4);
q_liq_max = q_emit_bruta * eta_fin - (mean(q_alb_abs_sso) + mean(q_ir_abs_sso));

q_net_pert_bol = q_liq_max - q_sol_pert_bol;
q_net_pert_eol = q_liq_max - q_sol_pert_eol;

q_req = 472.5;
T_rad_eq_eol = zeros(size(delta_theta_deg));
T_j_eq_eol = zeros(size(delta_theta_deg));

for idx = 1:length(delta_theta_deg)
    q_par = q_sol_pert_eol(idx) + mean(q_alb_abs_sso) + mean(q_ir_abs_sso);
    T_val = ((q_req + q_par)/(eps_ir_rad * sigma_SB * eta_fin) + T_sink_sso^4)^0.25;
    T_rad_eq_eol(idx) = T_val - 273.15;
    T_j_eq_eol(idx) = T_rad_eq_eol(idx) + 25.0;
end

idx_fail_85 = find(T_j_eq_eol >= 85.0, 1);
delta_fail_85 = delta_theta_deg(idx_fail_85);
fprintf('Ponto de Infracao Termica (T_j >= 85 C): delta_theta = %.2f graus\n', delta_fail_85);

idx_trip_105 = find(T_j_eq_eol >= 105.0, 1);
if ~isempty(idx_trip_105)
    delta_trip_105 = delta_theta_deg(idx_trip_105);
else
    delta_trip_105 = 30.0;
end

% ------------------------------------------------------------------------------
% 5. CONSTRUCAO DA FIGURA 4.1 COMPOSTA MULTIPAINEL (GNUPLOT DE ALTA QUALIDADE)
% ------------------------------------------------------------------------------
fprintf('3. Gerando Figura 4.1 Composta (Dinamica Orbital e Ponto de Falha de Atitude)...\n');

graphics_toolkit('gnuplot');

fig_4_1 = figure(1, 'position', [50, 50, 1100, 750]);

% Subpainel (a1): Perfil de Iluminacao Conica Suave S(t)
subplot(2, 2, 1);
hold on; grid on; box on;
t_norm_sso = t_sso / 60;
t_norm_leo = t_leo / 60;
plot(t_norm_sso, S_sso, 'b-', 'linewidth', 1.8);
plot(t_norm_leo, S_leo, 'r--', 'linewidth', 1.6);
xlabel('Tempo Orbital (min)', 'fontsize', 9);
ylabel('Fator Solar S(t) [-]', 'fontsize', 9);
title('(a1) Ocultacao Conica Suave (Umbra/Penumbra)', 'fontsize', 10, 'fontweight', 'bold');
legend({'SSO Dawn-Dusk 700 km', 'LEO Convencional 500 km'}, 'location', 'southwest', 'fontsize', 8);
xlim([0, 100]); ylim([-0.05, 1.15]);

% Subpainel (a2): Decomposicao dos Fluxos Radiativos Parasitas Incidentes
subplot(2, 2, 2);
hold on; grid on; box on;
plot(t_norm_sso, q_parasita_sso, 'b-', 'linewidth', 2.0);
plot(t_norm_leo, q_parasita_leo, 'r-', 'linewidth', 1.8);
plot(t_norm_leo, q_sol_abs_leo, 'm:', 'linewidth', 1.3);
plot(t_norm_leo, q_alb_abs_leo, 'g-.', 'linewidth', 1.3);
plot(t_norm_leo, q_ir_abs_leo, 'k--', 'linewidth', 1.2);
xlabel('Tempo Orbital (min)', 'fontsize', 9);
ylabel('Fluxo Parasita (W/m^2)', 'fontsize', 9);
title('(a2) Fluxos Parasitas Absorvidos no Radiador EOL', 'fontsize', 10, 'fontweight', 'bold');
legend({'SSO Total (F_{rad-E}=0.0314)', 'LEO Total', 'LEO Solar', 'LEO Albedo', 'LEO IR'}, ...
       'location', 'northeast', 'fontsize', 8);
xlim([0, 100]); ylim([0, 220]);

% Subpainel (b1): Sensibilidade a Perturbacao Angular de Atitude (BOL vs. EOL)
subplot(2, 2, 3);
hold on; grid on; box on;
plot(delta_theta_deg, q_net_pert_bol, 'g-', 'linewidth', 1.8);
plot(delta_theta_deg, q_net_pert_eol, 'r-', 'linewidth', 2.0);
plot([0, 30], [q_req, q_req], 'k--', 'linewidth', 1.4);
plot([3.5, 3.5], [0, 520], 'b:', 'linewidth', 1.4);
xlabel('Desvio de Atitude |\delta\theta_{att}| (graus)', 'fontsize', 9);
ylabel('Capacidade Liquida (W/m^2)', 'fontsize', 9);
title('(b1) Degradacao da Rejeicao por Desvio de Atitude', 'fontsize', 10, 'fontweight', 'bold');
legend({'BOL (\alpha_s = 0.08)', 'EOL (\alpha_s = 0.14)', 'Requisito (472.5 W/m^2)', 'Margem (\le 3.5^\circ)'}, ...
       'location', 'southwest', 'fontsize', 8);
xlim([0, 25]); ylim([250, 500]);

% Subpainel (b2): Curva de Ponto de Falha Eletrotermica (T_rad e T_juncao)
subplot(2, 2, 4);
hold on; grid on; box on;
plot(delta_theta_deg, T_rad_eq_eol, 'b-', 'linewidth', 1.8);
plot(delta_theta_deg, T_j_eq_eol, 'r-', 'linewidth', 2.0);
plot([0, 30], [85, 85], 'm--', 'linewidth', 1.4);
plot([0, 30], [105, 105], 'k-.', 'linewidth', 1.4);
plot([delta_fail_85, delta_fail_85], [20, 115], 'm:', 'linewidth', 1.2);
plot([delta_trip_105, delta_trip_105], [20, 115], 'k:', 'linewidth', 1.2);
xlabel('Desvio de Atitude |\delta\theta_{att}| (graus)', 'fontsize', 9);
ylabel('Temperatura de Equilibrio (^oC)', 'fontsize', 9);
title('(b2) Ponto de Falha: Throttling (85 ^oC) e Trip (105 ^oC)', 'fontsize', 10, 'fontweight', 'bold');
legend({'T_{rad} EOL', 'T_j EOL', 'Throttling (85 ^oC)', 'Trip (105 ^oC)'}, ...
       'location', 'northwest', 'fontsize', 8);
xlim([0, 25]); ylim([35, 115]);

% Exportacao em arquivo final
root_dir = fileparts(fileparts(mfilename('fullpath')));
fig_dir = fullfile(root_dir, 'figuras');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
output_path = fullfile(fig_dir, 'figura_4_1_fluxos_orbitais_comparativos.png');
print(fig_4_1, output_path, '-dpng', '-r150');
fprintf('=== FIGURA 4.1 SALVA COM SUCESSO EM: %s ===\n\n', output_path);
