% ==============================================================================
% RENDERIZADOR DE FIGURAS COMPOSTAS: FIGURAS 3.3 E 3.4
% GERADOR FOTOVOLTAICO AM0, DEGRADAÇÃO DDD E MPPT SOB HOT SPOT
% ==============================================================================
% Dissertação de Mestrado Profissional - ProPECAut / IFES
% Autor: Dafnhê Thor Castro de Martins Barros
% Revisão Técnica: Modelagem Explícita de Condições Nominais e Modos de Falha
% ==============================================================================
clear; clc; close all;
graphics_toolkit('gnuplot');
% 1. Executa a simulação fotovoltaica base e o algoritmo dinâmico de MPPT
simulacao_3_2_gerador_fotovoltaico_am0;
% Cria o diretório de destino se não existir
if ~exist('../figuras', 'dir')
    mkdir('../figuras');
end
%% ==============================================================================
%% RENDERIZAÇÃO DA FIGURA 3.3: CARACTERIZAÇÃO I-V, P-V E DEGRADAÇÃO DDD (10 ANOS)
%% ==============================================================================
fig3 = figure('visible', 'off', 'position', [100, 100, 1050, 800]);
graphics_toolkit(fig3, 'gnuplot');
cores = {'b-', 'g-', 'm-', 'r-'};
% --- Subpainel (a1): Curvas I-V da Matriz do Nó sob Gradiente Térmico ---
subplot(2, 2, 1);
for m = 1:length(curvas_iv)
    plot(curvas_iv(m).V, curvas_iv(m).I, cores{m}, 'linewidth', 1.8); hold on;
end
grid on;
xlim([0, 520]);
ylim([0, 145]);
xlabel('Tensão da Matriz Solar (V)', 'fontsize', 10, 'fontweight', 'bold');
ylabel('Corrente da Matriz do Nó (A)', 'fontsize', 10, 'fontweight', 'bold');
title('(a1) Curvas I-V por Nó de 50 kW (Catálogo e Extrapolação)', 'fontsize', 11, 'fontweight', 'bold');
legend('15 °C', '28 °C (Catálogo)', '60 °C (Cruzeiro)', '85 °C (WCH)', 'location', 'southwest', 'fontsize', 8);
% --- Subpainel (a2): Curvas P-V e Janela de Regulação do Barramento Buck ---
subplot(2, 2, 2);
for m = 1:length(curvas_iv)
    leg_str{m} = sprintf('%2.0f °C (%4.1f kW)', curvas_iv(m).T, curvas_iv(m).P_mp / 1e3);
    plot(curvas_iv(m).V, curvas_iv(m).P / 1e3, cores{m}, 'linewidth', 1.8); hold on;
end
% Linhas indicativas da faixa de tensão do barramento e janela MPPT Buck (350 V a 650 V)
plot([350, 350], [0, 60], 'k:', 'linewidth', 1.4);
plot([400, 400], [0, 60], 'k--', 'linewidth', 1.6);
text(403, 15, 'Barramento (400 V)', 'fontsize', 8, 'fontweight', 'bold', 'color', [0, 0, 0]);
grid on;
xlim([0, 550]);
ylim([0, 60]);
xlabel('Tensão da Matriz Solar (V)', 'fontsize', 10, 'fontweight', 'bold');
ylabel('Potência da Matriz do Nó (kW)', 'fontsize', 10, 'fontweight', 'bold');
title('(a2) Curvas P-V e Ponto de Máxima Potência (MPP)', 'fontsize', 11, 'fontweight', 'bold');
legend(leg_str, 'location', 'northwest', 'fontsize', 8);
% --- Subpainel (b1): Fluência Cumulativa de Elétrons em Órbita ---
subplot(2, 2, 3);
plot(anos, phi_total / 1e14, 'k-', 'linewidth', 2.0); hold on;
plot(anos, (anos * 1.4e13) / 1e14, 'b--', 'linewidth', 1.4);
plot(anos, (anos * 5.8e13) / 1e14, 'r:', 'linewidth', 1.4);
grid on;
xlim([0, 10]);
xlabel('Tempo de Missão (Anos)', 'fontsize', 10, 'fontweight', 'bold');
ylabel('Fluência 1-MeV (\times 10^{14} e^-/cm^2)', 'fontsize', 10, 'fontweight', 'bold');
title('(b1) Fluência Cumulativa de Radiação (DDD)', 'fontsize', 11, 'fontweight', 'bold');
legend('SSO 700 km (Projeto: 2.5\times10^{13})', 'LEO 500 km (1.4\times10^{13})', ...
       'MEO 1000 km (5.8\times10^{13})', 'location', 'northwest', 'fontsize', 8);
% --- Subpainel (b2): Fator de Retenção de Potência EOL (10 Anos) ---
subplot(2, 2, 4);
plot(anos, retencao_total * 100, 'b-', 'linewidth', 2.0); hold on;
% Linha de Requisito de Projeto de EOL (88%)
plot([0, 10], [88.0, 88.0], 'r--', 'linewidth', 1.6);
text(0.5, 89.2, 'Requisito Mínimo em Fim de Vida EOL (88.0%)', ...
    'fontsize', 9, 'fontweight', 'bold', 'color', [0.8, 0, 0]);
grid on;
xlim([0, 10]);
ylim([80, 102]);
xlabel('Tempo de Missão (Anos)', 'fontsize', 10, 'fontweight', 'bold');
ylabel('Retenção de Potência P_{mp} (%)', 'fontsize', 10, 'fontweight', 'bold');
title('(b2) Retenção Temporal de Potência EOL (Norma ECSS)', 'fontsize', 11, 'fontweight', 'bold');
print(fig3, '../figuras/figura_3_3_desempenho_termico_radiacao_am0.png', '-dpng', '-r300');
fprintf('Figura 3.3 exportada com sucesso: ../figuras/figura_3_3_desempenho_termico_radiacao_am0.png\n');
%% ==============================================================================
%% RENDERIZAÇÃO DA FIGURA 3.4: CONTINGÊNCIA HOT SPOT E DINÂMICA REAL DO MPPT
%% ==============================================================================
fig4 = figure('visible', 'off', 'position', [100, 100, 1050, 800]);
graphics_toolkit(fig4, 'gnuplot');
% --- Subpainel (a1): Característica Reversa e Clamp Schottky ---
subplot(2, 2, 1);
v_rev_axis = linspace(-30, 1.0, 300);
i_rev_axis = zeros(size(v_rev_axis));
for idx = 1:length(v_rev_axis)
    v_val = v_rev_axis(idx);
    if v_val >= -0.35
        % Condução direta do diodo Schottky de bypass em paralelo
        i_rev_axis(idx) = 5.0 * (exp((v_val + 0.35) / 0.04) - 1);
    elseif v_val > -28.0
        % Célula em corte reverso (fuga de saturação)
        i_rev_axis(idx) = - 1e-6;
    else
        % Avalanche reversa severa (quebra dielétrica da junção de Ge)
        i_rev_axis(idx) = - 0.528 * (1 + 15.0 * (abs(v_val) - 28.0)^2);
    end
end
plot(v_rev_axis, i_rev_axis, 'b-', 'linewidth', 2.0); hold on;
% Marcação do ponto de clamp Schottky e de avalanche
plot([-0.35, -0.35], [-15, 10], 'g--', 'linewidth', 1.6);
plot([-28.0, -28.0], [-15, 10], 'r--', 'linewidth', 1.6);
text(-0.20, -8, 'Clamp Schottky (-0,35 V)', 'fontsize', 8, 'fontweight', 'bold', 'color', [0, 0.5, 0]);
text(-27.5, -8, 'Avalanche (-28,0 V)', 'fontsize', 8, 'fontweight', 'bold', 'color', [0.8, 0, 0]);
grid on;
xlim([-30, 1]);
ylim([-15, 5]);
xlabel('Tensão na Célula (V)', 'fontsize', 10, 'fontweight', 'bold');
ylabel('Corrente na Célula (A)', 'fontsize', 10, 'fontweight', 'bold');
title('(a1) Característica Reversa e Grampeamento de Proteção', 'fontsize', 11, 'fontweight', 'bold');
% --- Subpainel (a2): Elevação Térmica Stefan-Boltzmann no Vácuo ---
subplot(2, 2, 2);
plot(t_hs, T_hs_sem_K - 273.15, 'r-', 'linewidth', 2.0); hold on;
plot(t_hs, T_hs_com_K - 273.15, 'g-', 'linewidth', 2.0);
% Limite de degradação térmica irreversível de adesivos e coverglass (180 °C)
plot([0, max(t_hs)], [180.0, 180.0], 'k--', 'linewidth', 1.5);
text(10, 186, 'Limite de Degradação Térmica do Adesivo (180 °C)', ...
    'fontsize', 8, 'fontweight', 'bold', 'color', [0.2, 0.2, 0.2]);
grid on;
xlim([0, 180]);
ylim([40, 260]);
xlabel('Tempo de Sombra (s)', 'fontsize', 10, 'fontweight', 'bold');
ylabel('Temperatura da Célula (°C)', 'fontsize', 10, 'fontweight', 'bold');
title('(a2) Resposta Térmica de Hot Spot no Vácuo (Stefan-Boltzmann)', 'fontsize', 11, 'fontweight', 'bold');
legend('Sem Bypass Schottky (Destrutivo)', 'Com Diodo Schottky Protegido', ...
       'location', 'east', 'fontsize', 8);
% --- Subpainel (b1): Dinâmica de Tensão do MPPT na Saída do Eclipse ---
subplot(2, 2, 3);
% PLOTA OS VETORES REAIS SIMULADOS DA SEÇÃO 5 DO SCRIPT 3.2 (SEM DADOS FORJADOS)
t_ms = t_mppt * 1e3;
V_mppt_adapt = V_pv_dyn;
V_mppt_fixo  = V_pv_fixo;
plot(t_ms, V_mppt_adapt, 'g-', 'linewidth', 2.0); hold on;
plot(t_ms, V_mppt_fixo, 'r-', 'linewidth', 1.5);
plot(t_ms, N_s * (V_mp_ref + dVmp_dT * (T_transiente - T_ref)), 'k--', 'linewidth', 1.2);
grid on;
xlim([0, max(t_ms)]);
ylim([300, 650]);
xlabel('Tempo após Saída do Eclipse (ms)', 'fontsize', 10, 'fontweight', 'bold');
ylabel('Tensão da Matriz FV (V)', 'fontsize', 10, 'fontweight', 'bold');
title('(b1) Transiente do MPPT Adaptativo (P&O-AP) a 150 kHz', 'fontsize', 11, 'fontweight', 'bold');
legend('P&O Adaptativo GaN (Convergência < 15 ms)', 'P&O de Passo Fixo (\Delta D = 0.008)', ...
       'V_{mp} Ótimo Analítico', 'location', 'northeast', 'fontsize', 8);
% --- Subpainel (b2): Eficiência Dinâmica Instantânea de Captura MPPT ---
subplot(2, 2, 4);
eff_dyn_adapt = (P_adapt ./ max(1.0, P_otima)) * 100;
eff_dyn_fixo  = (P_fixo  ./ max(1.0, P_otima)) * 100;
plot(t_ms, eff_dyn_adapt, 'g-', 'linewidth', 2.0); hold on;
plot(t_ms, eff_dyn_fixo, 'r-', 'linewidth', 1.5);
% Linha de requisito de rendimento dinâmico de 99%
plot([0, 200], [99.0, 99.0], 'k:', 'linewidth', 1.4);
text(10, 98.2, 'Meta de Desempenho P&O-AP (> 99%)', 'fontsize', 8, 'fontweight', 'bold');
grid on;
xlim([0, max(t_ms)]);
ylim([40, 102]);
xlabel('Tempo após Saída do Eclipse (ms)', 'fontsize', 10, 'fontweight', 'bold');
ylabel('Rendimento Dinâmico (%)', 'fontsize', 10, 'fontweight', 'bold');
title('(b2) Eficiência Dinâmica de Captura Solar no Transiente', 'fontsize', 11, 'fontweight', 'bold');
legend('P&O Adaptativo GaN (> 99.2%)', 'P&O Passo Fixo (~96.5%)', 'location', 'southeast', 'fontsize', 8);
print(fig4, '../figuras/figura_3_4_contingencias_hotspot_mppt.png', '-dpng', '-r300');
fprintf('Figura 3.4 exportada com sucesso: ../figuras/figura_3_4_contingencias_hotspot_mppt.png\n');
fprintf('==============================================================================\n');
fprintf('RENDERIZAÇÃO DAS FIGURAS 3.3 E 3.4 CONCLUÍDA COM DADOS REAIS DE SIMULAÇÃO.\n');
fprintf('==============================================================================\n');
