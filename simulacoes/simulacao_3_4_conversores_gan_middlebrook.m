function simulacao_3_4_conversores_gan_middlebrook()
  % ============================================================================
  % SIMULACAO 3.4: CONVERSORES GaN (GS66516T), STEINMETZ MODIFICADO (MSE,
  % NUCLEO FERROXCUBE U93/152/30 3C95), INTERFACE AlN+PGS E MIDDLEBROOK
  % ============================================================================
  % Dissertacao de Mestrado Profissional - ProPECAut / IFES
  % Autor: Dafnhe Thor Castro de Martins Barros
  % ============================================================================
 function simulacao_3_4_conversores_gan_middlebrook()
  close all; clc;

  % =========================================================================
  % 1. Parametros do Conversor e do Transformador (400 V -> 48 V, 20 kW)
  % =========================================================================
  P_mod = 20000.0;            % Potencia nominal do estagio [W]
  V_in = 400.0;               % Tensao nominal de entrada [V]
  V_out = 48.0;               % Tensao de saida [V]
  f_sw = 150e3;               % Frequencia de chaveamento [Hz]
  N_sw = 4;                   % Ponte completa primária (4 chaves)
  n_tr = 400.0 / 55.0;        % Relacao de transformacao

  I_out = P_mod / V_out;      % 416.7 A
  D_eff = (V_out * n_tr) / V_in; % Razao ciclica efetiva (~0.44)
  I_pri_pk = I_out / n_tr;    % Corrente de pico refletida no primario (~57.3 A)

  % Corrente eficaz real por chave primária (forma de onda trapezoidal simetrica)
  I_rms_sw = I_pri_pk * sqrt(D_eff / 2);

  % =========================================================================
  % 2. Parametros Reais do GS66516T e Fisica de Comutacao GaN
  % =========================================================================
  Ron25 = 25e-3;              % 25 mOhm a 25 C
  alpha = (65e-3/Ron25 - 1) / (150 - 25); % Coeficiente real = 0.0128 K^-1
  Coss = 207e-12; Qg = 14.2e-9; Vg = 6.0;

  % Energias de comutacao sob VDS = 400 V (Zero Reverse Recovery, Qrr = 0)
  E_on_ref = 134.1e-6; E_off_ref = 17.0e-6; I_ref = 20.0;
  E_sw_chave = (E_on_ref + E_off_ref) * (I_pri_pk / I_ref);
  P_sw_chave = E_sw_chave * f_sw + (Qg * Vg * f_sw);

  % Perdas no nucleo de ferrite 3C95 e condutores
  k_fe = 0.70; a_st = 1.46; b_st = 2.52; Ve = 4.5e-4; dB = 0.15;
  P_core = k_fe * (f_sw^a_st) * (dB^b_st) * Ve;
  Rac_p = 15e-3; Rac_s = 0.25e-3;
  P_cu = (I_pri_pk * sqrt(D_eff))^2 * Rac_p + I_out^2 * Rac_s;

  % =========================================================================
  % 3. Resolucao Eletrotermica Iterativa e Limiar de Fuga Termica
  % =========================================================================
  Tc = 60.0;                  % Temperatura da placa fria [degC]
  Rth_jc = 0.27; Rth_tim = 0.85; % Resistencia de interface
  Rth_tot = Rth_jc + Rth_tim;

  Tj = Tc;
  for it = 1:500
    Ron_T = Ron25 * (1 + alpha * (Tj - 25));
    P_cond_sw = I_rms_sw^2 * Ron_T;
    P_tot_sw = P_cond_sw + P_sw_chave;
    Tj_new = Tc + Rth_tot * P_tot_sw;
    if abs(Tj_new - Tj) < 1e-4, break; end
    Tj = Tj_new;
  end

  % Condicao analitica de fuga termica: dP/dT < 1/Rth
  dP_dT = I_rms_sw^2 * Ron25 * alpha;
  Rth_crit = 1 / dP_dT;

  P_cond_tot = N_sw * P_cond_sw;
  P_sw_tot = N_sw * P_sw_chave;
  Ploss_tot = P_cond_tot + P_sw_tot + P_core + P_cu;
  eta = 100 * P_mod / (P_mod + Ploss_tot);

  fprintf('=== 1. Balanco Eletrotermico do Conversor GaN ===\n');
  fprintf('Corrente Eficaz por Chave: %.2f A | Tj Estavel: %.1f C (Tc = %.1f C)\n', ...
          I_rms_sw, Tj, Tc);
  fprintf('Perda Conducao: %.1f W | Chaveamento: %.1f W | Nucleo: %.1f W | Cobre: %.1f W\n', ...
          P_cond_tot, P_sw_tot, P_core, P_cu);
  fprintf('Perda Total: %.1f W | Rendimento do Estagio: %.2f %%\n', Ploss_tot, eta);
  fprintf('Resistencia Termica Critica (Rth_crit): %.2f K/W (Aplicada: %.2f K/W)\n', ...
          Rth_crit, Rth_tot);
  if Rth_tot < Rth_crit
    fprintf('Status: ESTAVEL CONTRA FUGA TERMICA (Margem: +%.2f K/W)\n\n', Rth_crit - Rth_tot);
  else
    fprintf('Status: ALERTA DE FUGA TERMICA\n\n');
  end

  % =========================================================================
  % 4. Criterio de Estabilidade de Impedancia de Middlebrook
  % =========================================================================
  fprintf('=== 2. Analise de Impedancia de Middlebrook (Barramento 48 V) ===\n');
  r_in = -V_out^2 / P_mod;    % Resistencia incremental negativa de CPL
  Z_lim = 0.10 * abs(r_in);  % Criterio com margem de 20 dB (10%)

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
