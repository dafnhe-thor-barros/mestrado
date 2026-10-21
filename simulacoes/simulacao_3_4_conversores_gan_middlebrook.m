% ==============================================================================
% SIMULACAO 3.4: CONVERSORES ESTATICOS CC-CC GaN DE ALTA FREQUENCIA,
% BALANCO DE PERDAS TERMICAS E ESTABILIDADE CPL POR CRITERIO DE MIDDLEBROOK
% ==============================================================================
% Dissertacao de Mestrado Profissional - ProPECAut / IFES
% Autor: Dafnhe Thor Castro de Martins Barros
% Revisao Tecnica: Modelagem Explicita de Condicoes Nominais e Modos de Falha
% ==============================================================================
function simulacao_3_4_conversores_gan_middlebrook()
  close all; clc;
  P_mod = 20000.0; V_in = 400.0; V_out = 48.0; f_sw = 150e3; N_sw = 4;
  I_in = 51.2; D = 0.45;  % I_in arredondado de P/(V*eta) = 51,28 A
  I_out = P_mod / V_out;
  I_rms_sw = I_in * sqrt(D);
  Ron25 = 25e-3; alpha = (65e-3/Ron25 - 1) / (150 - 25);
  Coss = 207e-12; Qg = 14.2e-9; Vg = 6.0; tr = 12e-9; tf = 15e-9; I_turn = 25.0;
  P_ov = V_in * I_turn * f_sw * (tr + tf) / 2;
  P_oss = 0.5 * Coss * V_in^2 * f_sw; P_gate = Qg * Vg * f_sw;
  P_sw_chave = P_ov + P_oss + P_gate;
  k_fe = 0.70; a_st = 1.46; b_st = 2.52; Ve = 4.5e-4; dB = 0.15;
  P_core = k_fe * (f_sw^a_st) * (dB^b_st) * Ve;
  P_cu = 17.70 + 43.40;
  Ron_sec = 0.80e-3 * (1 + 0.004 * (80 - 25));
  P_cond_sec = I_out^2 * Ron_sec;
  Rth_tot = 0.27 + 0.85;
  % Cenario de referencia: Tj = 80 C (Ron = 42,6 mOhm) e cenario com Tc = 60 C
  Ron80 = Ron25 * (1 + alpha * (80 - 25));
  P_cond_80 = I_rms_sw^2 * Ron80;
  Tc_para_80 = 80 - Rth_tot * (P_cond_80 + P_sw_chave);
  Ploss_ref = N_sw * (P_cond_80 + P_sw_chave) + P_cond_sec + P_core + P_cu;
  eta_ref = 100 * P_mod / (P_mod + Ploss_ref);
  Tc = 60.0; Tj = Tc;
  for it = 1:500
    Ron_T = Ron25 * (1 + alpha * (Tj - 25));
    Tj_new = Tc + Rth_tot * (I_rms_sw^2 * Ron_T + P_sw_chave);
    if abs(Tj_new - Tj) < 1e-6, break; end
    Tj = Tj_new;
  end
  P_cond_Tc = I_rms_sw^2 * Ron25 * (1 + alpha * (Tj - 25));
  Ploss_Tc = N_sw * (P_cond_Tc + P_sw_chave) + P_cond_sec + P_core + P_cu;
  eta_Tc = 100 * P_mod / (P_mod + Ploss_Tc);
  dP_dT = I_rms_sw^2 * Ron25 * alpha; Rth_crit = 1 / dP_dT;
  fprintf('=== 1. Balanco Eletrotermico do Conversor GaN ===\n');
  fprintf('I_in = %.2f A | D = %.2f | I_rms por chave = %.2f A\n', I_in, D, I_rms_sw);
  fprintf('Ron(80C) = %.2f mOhm | P_cond/chave(80C) = %.2f W | P_cond_pri = %.2f W\n', 1e3*Ron80, P_cond_80, N_sw*P_cond_80);
  fprintf('P_ov = %.2f W | P_oss = %.3f W | P_gate = %.4f W | P_sw/chave = %.2f W | P_sw_tot = %.2f W\n', P_ov, P_oss, P_gate, P_sw_chave, N_sw*P_sw_chave);
  fprintf('P_cond_sec = %.2f W | P_core = %.2f W | P_cu = %.2f W\n', P_cond_sec, P_core, P_cu);
  fprintf('[Ref. Tj = 80 C] Perda Total: %.2f W | Rendimento: %.3f %% | Tc requerido: %.1f C\n', Ploss_ref, eta_ref, Tc_para_80);
  fprintf('[Tc = %.0f C] Tj convergido = %.1f C | Perda Total: %.2f W | Rendimento: %.3f %%\n', Tc, Tj, Ploss_Tc, eta_Tc);
  if Tj > 150, fprintf('ALERTA: Tj acima do maximo absoluto de 150 C do GS66516T\n'); else fprintf('Tj dentro do limite de 150 C do GS66516T (margem %.1f K)\n', 150 - Tj); end
  % --- A10: CLARIFICACAO DO CIRCUITO TERMICO EQUIVALENTE INDIVIDUAL ---
  % Rth_tot representa o acoplamento termico INDIVIDUAL da chave ao dissipador
  fprintf('=== 1. Estabilidade Eletrotermica Individual por Chave GaN ===\n');
  fprintf('dP/dTj por chave = %.4f W/K | Limite Rth_crit individual = %.2f K/W\n', dP_dT, Rth_crit);
  fprintf('Interface aplicada por chave (AlN + PGS): %.2f K/W\n', Rth_tot);
  if Rth_tot < Rth_crit
    fprintf('Status: ESTAVEL (Margem por chave contra runaway: +%.2f K/W)\n\n', Rth_crit - Rth_tot);
  else
    fprintf('Status: ALERTA DE FUGA TERMICA REGENERATIVA\n\n');
  end

% --- A11: EXPORTACAO DE DADOS PARA A FIGURA 3.8 (Abordagem 2) ---
save('-v7', 'simulacoes/dados_simulacao_3_4.mat', ...
     'f', 's', 'casos', 'r_in', 'Z_lim', 'P_mod', 'Ploss_ref', 'eta_ref');
fprintf('Dados salvos em simulacoes/dados_simulacao_3_4.mat (Abordagem 2)\n');% 4. Criterio de Impedancia de Middlebrook (Barramento 48 V)
  fprintf('=== 2. Analise de Impedancia de Middlebrook ===\n');
  r_in = -V_out^2 / P_mod;
  Z_lim = 0.10 * abs(r_in);
  
  f = logspace(1, 6, 5000); s = 1i * 2 * pi * f;
  L_line = 200e-9; R_line = 0.8e-3;

  casos = {'Banco Saudavel (15 mF / 1.2 mOhm)', 15e-3, 1.2e-3;
           'WCA Texto Tese (3 mF / 8 mOhm)',     3e-3,  8e-3;
           'Script Antigo (2.5 mF / 9.5 mOhm)',  2.5e-3, 9.5e-3;
           'Fig. 6 Operacional (2.2 mF / 3 mOhm)', 2.2e-3, 3e-3;
           'Fig. 6 Falha (0.44 mF / 15 mOhm)',   0.44e-3, 15e-3};

  fprintf('r_in da CPL = %.4f Ohm | Limite Middlebrook: %.5f Ohm\n', r_in, Z_lim);
  for c = 1:size(casos, 1)
    C_b = casos{c,2}; ESR_b = casos{c,3};
    Z1 = s * L_line + R_line;
    Z2 = 1 ./ (s * C_b) + ESR_b;
    Z = abs(Z1 .* Z2 ./ (Z1 + Z2));
    [Zpk, idx_f] = max(Z);
    if Zpk <= Z_lim, st = 'ESTAVEL'; else, st = 'VIOLA'; end
    fprintf('%-36s Pico: %.5f Ohm @ %.0f Hz -> %s\n', ...
            casos{c,1}, Zpk, f(idx_f), st);
  end
end
