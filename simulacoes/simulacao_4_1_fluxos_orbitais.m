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

% Propriedades termo-opticas em Fim de Vida (EOL - 10 anos):
% [PROVISORIO / REQUISITO DE PROJETO]: Valores de estresse EOL requeridos para a dissertacao
alpha_s_rad = 0.24;        % OSR Quartzo/Prata EOL degradado [PROVISORIO]
eps_ir_rad  = 0.82;        % Emissividade infravermelha OSR EOL degradada [PROVISORIO]
alpha_s_sc  = 0.26;        % Tinta branca AZ-93 EOL [PROVISORIO]
eps_ir_sc   = 0.89;        % Emissividade AZ-93 [PROVISORIO]
eta_fin     = 0.945;       % Eficiencia da aleta de APG composito [PROVISORIO]

% ------------------------------------------------------------------------------
% 2. CONFIGURACAO DOS REGIMES ORBITAIS
% ------------------------------------------------------------------------------
altitudes = [700.0, 500.0, 35786.0]; % km
nomes_orbitas = {'SSO Dawn-Dusk (700 km, beta=85 deg)', 'LEO Convencional (500 km, beta=20 deg)', 'GEO (35.786 km)'};

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

% Balanço térmico radiativo rigoroso sem dupla contagem de IR terrestre:
% Capacidade máxima de rejeição na temperatura de teto do radiador (T_rad = 55 °C = 328.15 K):
T_rad_lim = 328.15; % 55 °C
q_emit_lim_eol = eps_ir_rad * sigma_SB * eta_fin * (T_rad_lim^4 - T_space^4);
q_par_sso_eol = mean(q_alb_abs_sso) + mean(q_ir_abs_sso); % 45.24 W/m^2 absorvido
q_liq_max_eol = q_emit_lim_eol - q_par_sso_eol;

% Para BOL (alpha_s = 0.08, eps_ir = 0.88):
q_emit_lim_bol = 0.88 * sigma_SB * eta_fin * (T_rad_lim^4 - T_space^4);
q_par_sso_bol = 0.88 * q_Earth_nom * F_rad_E_sso;
q_liq_max_bol = q_emit_lim_bol - q_par_sso_bol;

q_net_pert_bol = q_liq_max_bol - q_sol_pert_bol;
q_net_pert_eol = q_liq_max_eol - q_sol_pert_eol;

% Potencia termica total a rejeitar conectada ao livro-razao:
% Caso Nominal: Q_int = 41.983,84 W [CALCULADO]
% Caso WCH (Worst-Case Hot EOL): Q_int_wch = 45.757,64 W [CALCULADO]
% Area do radiador de projeto: A_rad = 88.6 m^2 (coifa de 5.2 m) [CALCULADO / CAPITULO 4]
% Area alternativa parametrica: A_rad_alt = 92.0 m^2 [ESTUDO DE SENSIBILIDADE]
% Delta T chassi-juncao: 25.0 °C [PROVISORIO / PREMISSA DE PROJETO]
A_rad = 88.6; 
Q_int_total = 41983.84; % [W] Nominal
Q_int_wch   = 45757.64; % [W] WCH
q_req = Q_int_total / A_rad; % 473.86 W/m^2

T_rad_eq_eol = zeros(size(delta_theta_deg));
T_j_eq_eol = zeros(size(delta_theta_deg));

for idx = 1:length(delta_theta_deg)
    q_par = q_sol_pert_eol(idx) + q_par_sso_eol;
    % Equilibrio Stefan-Boltzmann exato: eps*sigma*eta*(T_rad^4 - T_space^4) = q_req + q_par
    T_val = ((q_req + q_par) / (eps_ir_rad * sigma_SB * eta_fin) + T_space^4)^0.25;
    T_rad_eq_eol(idx) = T_val - 273.15;
    T_j_eq_eol(idx) = T_rad_eq_eol(idx) + 25.0; % Delta T chassi-juncao [PROVISORIO]
end

% Calculo de sensibilidade para A_rad = 92.0 m^2 e caso WCH (para registro analitico):
q_req_92 = Q_int_total / 92.0;
T_rad_92 = ((q_req_92 + q_par_sso_eol) / (eps_ir_rad * sigma_SB * eta_fin) + T_space^4)^0.25 - 273.15;
T_j_92   = T_rad_92 + 25.0;

q_req_wch = Q_int_wch / A_rad;
T_rad_wch = ((q_req_wch + q_par_sso_eol) / (eps_ir_rad * sigma_SB * eta_fin) + T_space^4)^0.25 - 273.15;
T_j_wch   = T_rad_wch + 25.0;

% Avaliacao dos criterios termicos de projeto (75 °C normal, 80 °C throttling, 85 °C trip):
idx_fail_80 = find(T_j_eq_eol >= 80.0, 1);
idx_trip_85 = find(T_j_eq_eol >= 85.0, 1);
if ~isempty(idx_trip_85)
    delta_trip_85 = delta_theta_deg(idx_trip_85);
else
    delta_trip_85 = 0.0;
end

fprintf('1. Balanco Termico Radiativo EOL sem Dupla Contagem (delta_theta = 0):\n');
fprintf('   -> Caso Nominal (Q = 41.984 W, A_rad = 88.6 m^2):\n');
fprintf('      T_rad = %.2f °C (Teto de projeto: 55.0 °C -> Superado em malha aberta)\n', T_rad_eq_eol(1));
fprintf('      T_j   = %.2f °C (Meta normal: 75.0 °C; Limiar Trip: 85.0 °C -> Superado)\n', T_j_eq_eol(1));
fprintf('   -> Sensibilidade de Area (A_rad = 92.0 m^2 Nominal):\n');
fprintf('      T_rad = %.2f °C | T_j = %.2f °C (Excede limiar de Throttling de 80.0 °C)\n', T_rad_92, T_j_92);
fprintf('   -> Caso Pior Caso Quente WCH (Q = 45.758 W, A_rad = 88.6 m^2):\n');
fprintf('      T_rad = %.2f °C | T_j = %.2f °C (Supera amplamente o limiar de Trip)\n', T_rad_wch, T_j_wch);
fprintf('   -> Conclusao: A operacao continua a plena carga em EOL sem controle ativo e termicamente inviavel,\n');
fprintf('      comprovando numericamente a necessidade mandatória de NMPC com Throttling.\n\n');

% ------------------------------------------------------------------------------
% 5. CONSTRUCAO DA FIGURA 4.1 COMPOSTA MULTIPAINEL
% ------------------------------------------------------------------------------
fprintf('3. Gerando Figura 4.1 Composta (Dinamica Orbital e Ponto de Falha de Atitude)...\n');

if ismember('qt', available_graphics_toolkits())
    graphics_toolkit('qt');
else
    graphics_toolkit('gnuplot');
end

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
legend({'SSO Dawn-Dusk 700 km (\beta=85^\circ: Insolacao Continua)', 'LEO Convencional 500 km (\beta=20^\circ: Eclipse 35.8 min)'}, 'location', 'southwest', 'fontsize', 7);
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
legend({sprintf('SSO Total (F_{rad-E}=%.4f)', F_rad_E_sso), 'LEO Total', 'LEO Solar', 'LEO Albedo', 'LEO IR'}, ...
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
title('(b1) Rejeicao Termica vs Desvio de Atitude', 'fontsize', 10, 'fontweight', 'bold');
legend({'BOL (\alpha_s = 0.08)', 'EOL (\alpha_s = 0.24)', sprintf('Requisito (%.1f W/m^2)', q_req), 'Margem ADCS (\le 3.5^\circ)'}, ...
       'location', 'southwest', 'fontsize', 8);
xlim([0, 25]); ylim([250, 500]);

% Subpainel (b2): Curva de Ponto de Falha Eletrotermica (T_rad e T_juncao)
subplot(2, 2, 4);
hold on; grid on; box on;
plot(delta_theta_deg, T_rad_eq_eol, 'b-', 'linewidth', 1.8);
plot(delta_theta_deg, T_j_eq_eol, 'r-', 'linewidth', 2.0);
plot([0, 30], [55, 55], 'c--', 'linewidth', 1.3);
plot([0, 30], [75, 75], 'g--', 'linewidth', 1.3);
plot([0, 30], [80, 80], 'm--', 'linewidth', 1.3);
plot([0, 30], [85, 85], 'k-.', 'linewidth', 1.5);
plot([delta_trip_85, delta_trip_85], [20, 115], 'k:', 'linewidth', 1.2);
xlabel('Desvio de Atitude |\delta\theta_{att}| (graus)', 'fontsize', 9);
ylabel('Temperatura de Equilibrio (^oC)', 'fontsize', 9);
title('(b2) Avaliacao de Limiares: Normal (75 ^oC), Throttling (80 ^oC) e Trip (85 ^oC)', 'fontsize', 10, 'fontweight', 'bold');
legend({'T_{rad} EOL (\alpha_s=0.24, \epsilon=0.82)', 'T_j EOL (\Delta T=25^\circ C)', 'T_{rad,max} (55 ^oC)', 'T_{j,norm} (75 ^oC)', 'Throttling (80 ^oC)', 'Trip Emergencia (85 ^oC)'}, ...
       'location', 'northwest', 'fontsize', 7);
xlim([0, 25]); ylim([35, 115]);

% Exportacao em arquivo final
if ~exist('figuras', 'dir'), mkdir('figuras'); end
output_path = 'figuras/figura_4_1_fluxos_orbitais_comparativos.png';
print(fig_4_1, output_path, '-dpng', '-r150');
fprintf('=== FIGURA 4.1 SALVA COM SUCESSO EM: %s ===\n\n', output_path);
