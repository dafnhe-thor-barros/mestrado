% ==============================================================================
% RENDERIZADOR DE FIGURAS COMPOSTAS: FIGURAS 3.1 E 3.2
% DINÂMICA DO BARRAMENTO HVDC E PROTEÇÃO POR CHAVE SSPC
% ==============================================================================
% Dissertação de Mestrado Profissional - ProPECAut / IFES
% Autor: Dafnhê Thor Castro de Martins Barros
% Revisão Técnica: Modelagem Explícita de Condições Nominais e Modos de Falha
% ==============================================================================
clear; clc; close all;
graphics_toolkit('gnuplot');
% 1. Executa a simulação base para carregar os vetores numéricos reconciliados
simulacao_3_1_barramento_hvdc;
% Cria o diretório de destino se não existir
if ~exist('../figuras', 'dir')
    mkdir('../figuras');
end
%% ==============================================================================
%% RENDERIZAÇÃO DA FIGURA 3.1: BALANÇO ENERGÉTICO E TENSÃO DO BESS (5 ÓRBITAS)
%% ==============================================================================
fig1 = figure('visible', 'off', 'position', [100, 100, 1050, 800]);
graphics_toolkit(fig1, 'gnuplot');
% Normalização dos eixos para "Número de Órbitas Consecutivas (Ciclos 1 a 5)"
t_norm_leo = sim_bess.LEO_500km.t * 3600 / T_orb(1);
t_norm_sso = sim_bess.SSO_700km.t * 3600 / T_orb(2);
t_norm_geo = sim_bess.GEO_35786km.t * 3600 / T_orb(3);
% --- Subpainel (a): Evolução do Estado de Carga (SoC) ---
subplot(2, 1, 1);
plot(t_norm_leo, sim_bess.LEO_500km.soc, 'r-', 'linewidth', 2.0); hold on;
plot(t_norm_sso, sim_bess.SSO_700km.soc, 'b-', 'linewidth', 2.0);
plot(t_norm_geo, sim_bess.GEO_35786km.soc, 'k--', 'linewidth', 1.8);
% Linha de limite operacional seguro de projeto (SoC >= 20%)
plot([0, N_orbits], [20, 20], 'r--', 'linewidth', 1.5);
text(0.15, 23, 'Limite Operacional Seguro de Projeto (SoC >= 20% / DoD <= 80%)', ...
    'fontsize', 9, 'fontweight', 'bold', 'color', [0.8, 0, 0]);
grid on;
xlim([0, N_orbits]);
ylim([10, 105]);
set(gca, 'xtick', 0:1:N_orbits, 'fontsize', 10);
ylabel('Estado de Carga - SoC (%)', 'fontsize', 11, 'fontweight', 'bold');
title('(a) Evolução Dinâmica do SoC do Banco BESS (9,0 kWh) sob Carga de IA (Ciclos 1 a 5)', ...
    'fontsize', 12, 'fontweight', 'bold');
legend('LEO Coplanar (Eclipse 35.8 min, Throttling P_{IT,ecl}=6.5 kW)', ...
       'SSO Dawn-Dusk Sazonal (Eclipse 18.8 min, P_{IT,ecl}=18.0 kW)', ...
       'GEO Equinócio (Eclipse 72.0 min, Throttling P_{IT,ecl}=4.5 kW)', ...
       'location', 'southwest', 'fontsize', 9);
% --- Subpainel (b): Perfil de Tensão Terminal do Barramento HVDC ---
subplot(2, 1, 2);
plot(t_norm_leo, sim_bess.LEO_500km.v_term, 'r-', 'linewidth', 2.0); hold on;
plot(t_norm_sso, sim_bess.SSO_700km.v_term, 'b-', 'linewidth', 2.0);
plot(t_norm_geo, sim_bess.GEO_35786km.v_term, 'k--', 'linewidth', 1.8);
% Linhas guias de UVLO crítico e tensão nominal
plot([0, N_orbits], [V_bus_min, V_bus_min], 'r--', 'linewidth', 1.6);
text(0.15, V_bus_min + 3, sprintf('Limite de Subtensão Crítico V_{UV} = %3.1f V (3.0 V/célula)', V_bus_min), ...
    'fontsize', 9, 'fontweight', 'bold', 'color', [0.8, 0, 0]);
plot([0, N_orbits], [V_bus_nominal, V_bus_nominal], 'g:', 'linewidth', 1.4);
text(0.15, V_bus_nominal + 2, 'Tensão Nominal Primária (400.0 V)', ...
    'fontsize', 9, 'fontweight', 'bold', 'color', [0, 0.5, 0]);
grid on;
xlim([0, N_orbits]);
ylim([315, 435]);
set(gca, 'xtick', 0:1:N_orbits, 'fontsize', 10);
xlabel('Número de Órbitas Consecutivas (Ciclos 1 a 5)', 'fontsize', 11, 'fontweight', 'bold');
ylabel('Tensão Terminal do Barramento (V)', 'fontsize', 11, 'fontweight', 'bold');
title('(b) Perfil de Tensão do Barramento HVDC frente ao Limiar de Subtensão (UVLO)', ...
    'fontsize', 12, 'fontweight', 'bold');
print(fig1, '../figuras/figura_3_1_balanco_orbital_bess.png', '-dpng', '-r300');
fprintf('Figura 3.1 exportada com sucesso: ../figuras/figura_3_1_balanco_orbital_bess.png\n');
%% ==============================================================================
%% RENDERIZAÇÃO DA FIGURA 3.2: ANÁLISE DE CONTINGÊNCIAS, SSPC E TÉRMICA DO CHASSI
%% ==============================================================================
fig2 = figure('visible', 'off', 'position', [100, 100, 1050, 800]);
graphics_toolkit(fig2, 'gnuplot');
% --- Subpainel (a1): Transiente de Curto-Circuito Franco ---
subplot(2, 2, 1);
plot(sim_sc(1).t_us, sim_sc(1).i_kA, 'g-', 'linewidth', 1.8); hold on;
plot(sim_sc(2).t_us, sim_sc(2).i_kA, 'b-', 'linewidth', 1.8);
plot(sim_sc(3).t_us, sim_sc(3).i_kA, 'r-', 'linewidth', 1.8);
grid on;
xlim([0, 150]);
xlabel('Tempo (\mus)', 'fontsize', 10, 'fontweight', 'bold');
ylabel('Corrente de Falha (kA)', 'fontsize', 10, 'fontweight', 'bold');
title('(a1) Falha Franca R_{sc}=1 m\Omega e Abertura SSPC', 'fontsize', 11, 'fontweight', 'bold');
legend('Rápido: t_{trip}=8 \mus', 'Nominal: t_{trip}=15 \mus', 'WCA Retardo: t_{trip}=100 \mus', ...
       'location', 'northeast', 'fontsize', 8);
% --- Subpainel (a2): Estresse Cumulativo I^2t e Limite Adiabático SiC ---
subplot(2, 2, 2);
trips_bar = [8, 15, 100];
i2t_vals  = [sim_sc(1).i2t, sim_sc(2).i2t, sim_sc(3).i2t];
b = bar(1:3, i2t_vals, 0.45);
set(b, 'facecolor', [0.25, 0.45, 0.85]);
set(gca, 'xtick', 1:3, 'xticklabel', {'8 \mus (N=1)', '15 \mus (N=1)', '100 \mus (WCA)'}, 'fontsize', 9);
hold on;
% Linha do limite térmico de queima adiabática (1500 A^2s)
plot([0.5, 3.5], [1500, 1500], 'r--', 'linewidth', 1.8);
text(0.6, 1600, 'Limite Adiabático SiC MOSFET (1500 A^2s)', ...
    'fontsize', 9, 'fontweight', 'bold', 'color', [0.8, 0, 0]);
grid on;
ylim([0, 2400]);
ylabel('Energia Específica I^2t (A^2s)', 'fontsize', 10, 'fontweight', 'bold');
title('(a2) Estresse Térmico Cumulativo na Chave SiC', 'fontsize', 11, 'fontweight', 'bold');
% --- Subpainel (b1): Perdas Dissipativas nos Estágios de Potência ---
subplot(2, 2, 3);
p_loss_display = [P_loss_nom, P_loss_cons];
b2 = bar([1, 2], p_loss_display, 0.45);
set(b2, 'facecolor', [0.85, 0.35, 0.25]);
set(gca, 'xtick', [1, 2], 'xticklabel', {'GaN Nominal (97,5%)', 'Pior Caso WCA (90,0%)'}, 'fontsize', 9);
hold on;
grid on;
ylabel('Calor Residual Dissipado (W)', 'fontsize', 10, 'fontweight', 'bold');
title('(b1) Perdas no Nó de 50 kW (Eletrônica + OISL)', 'fontsize', 11, 'fontweight', 'bold');
text(0.75, P_loss_nom + 200, sprintf('%4.0f W', P_loss_nom), 'fontweight', 'bold', 'fontsize', 9);
text(1.75, P_loss_cons + 200, sprintf('%4.0f W', P_loss_cons), 'fontweight', 'bold', 'fontsize', 9);
ylim([0, 6000]);
% --- Subpainel (b2): Resposta Térmica Transiente da Carcaça Estrutural ---
subplot(2, 2, 4);
plot(t_th / 3600, T_chassi_nom, 'g-', 'linewidth', 2.0); hold on;
plot(t_th / 3600, T_chassi_cons, 'r--', 'linewidth', 2.0);
% Linha de limite operacional máximo permitido (65 °C - Tabela 3 da Dissertação)
plot([0, max(t_th)/3600], [65.0, 65.0], 'r:', 'linewidth', 1.8);
text(0.1, 67.5, 'Limite Máximo Operacional (65,0 °C)', ...
    'fontsize', 9, 'fontweight', 'bold', 'color', [0.8, 0, 0]);
grid on;
xlim([0, max(t_th)/3600]);
ylim([20, 125]);
set(gca, 'ytick', 20:20:120);
xlabel('Tempo de Operação (Horas)', 'fontsize', 10, 'fontweight', 'bold');
ylabel('Temperatura do Chassi (°C)', 'fontsize', 10, 'fontweight', 'bold');
title('(b2) Resposta Eletrotérmica Transiente do Chassi', 'fontsize', 11, 'fontweight', 'bold');
legend('GaN Nominal (97,5%, 120 W/K)', 'Pior Caso WCA (90,0%, 80 W/K)', ...
       'location', 'southeast', 'fontsize', 8);
print(fig2, '../figuras/figura_3_2_falhas_eletrotermicas_sspc.png', '-dpng', '-r300');
fprintf('Figura 3.2 exportada com sucesso: ../figuras/figura_3_2_falhas_eletrotermicas_sspc.png\n');
fprintf('==============================================================================\n');
fprintf('RENDERIZAÇÃO DAS FIGURAS 3.1 E 3.2 CONCLUÍDA COM TOTAL ADERÊNCIA AO TEXTO.\n');
fprintf('==============================================================================\n');
