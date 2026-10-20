function simulacao_3_2_gerador_fotovoltaico_am0()
  % ============================================================================
  % SIMULACAO 3.2: GERADOR FOTOVOLTAICO MULTIJUNCAO AM0 (SPECTROLAB XTJ PRIME)
  % MODELO FISICO DE 3 DIODOS INDEPENDENTES EM SERIE (GaInP / GaAs / Ge),
  % DEGRADACAO DDD, HOT SPOT E CONVERSOR BUCK CHAVEADO (150 kHz)
  % ============================================================================
  % Dissertacao de Mestrado Profissional - ProPECAut / IFES
  % Autor: Dafnhe Thor Castro de Martins Barros
  % Referencias Primarias:
  %   [1] Spectrolab XTJ Prime Space Solar Cell Datasheet (AM0, 1353 W/m^2, 28 C)
  %   [2] Cotal et al. (2000), Triple-Junction Solar Cell Efficiencies Above 32%
  % ============================================================================
  close all; clc;

  fprintf('==============================================================================\n');
  fprintf('INICIANDO SIMULACAO 3.2: MODELO DE 3 DIODOS EM SERIE (GaInP/GaAs/Ge) E BUCK\n');
  fprintf('==============================================================================\n\n');

  % ----------------------------------------------------------------------------
  % 1. Constantes fisicas e parametros de datasheet (Spectrolab XTJ Prime)
  % ----------------------------------------------------------------------------
  qe = 1.602176634e-19;      % [C] Carga elementar
  kB = 1.380649e-23;         % [J/K] Constante de Boltzmann
  sigmaSB = 5.670374e-8;     % [W/(m^2 K^4)] Constante de Stefan-Boltzmann

  Tref = 28.0;               % [degC] Temperatura padrao de ensaio
  G_datasheet = 1353.0;      % [W/m^2] Irradiancia de ensaio Spectrolab AM0
  G_AM0 = 1361.0;            % [W/m^2] Constante solar ASTM E490 (1 UA)
  scale_AM0 = G_AM0 / G_datasheet;

  Acm2 = 30.18;              % [cm^2] Area ativa por celula
  Am2  = Acm2 * 1e-4;        % [m^2]
  VocRef_tot = 2.720;        % [V] Voc total do datasheet (1.42 + 1.01 + 0.29 V)
  JscRef     = 18.0e-3;      % [A/cm^2] Densidade de corrente de curto-circuito
  VmpRef_cat = 2.406;        % [V] Tensao MPP de catalogo
  JmpRef     = 17.5e-3;      % [A/cm^2] Densidade de corrente MPP de catalogo

  IscRef = JscRef * Acm2 * scale_AM0; % [A] 0.5465 A sob 1361 W/m^2
  ImpRef = JmpRef * Acm2 * scale_AM0; % [A] 0.5313 A sob 1361 W/m^2

  % Coeficientes termicos oficiais do datasheet Spectrolab XTJ Prime
  dVoc_dT_tot = -5.60e-3;    % [V/K] Coeficiente termico global de Voc
  dJsc_dT     = +10.0e-6;    % [A/(cm^2 K)] Coeficiente termico de Jsc
  dIsc_dT     = dJsc_dT * Acm2 * scale_AM0; % [A/K] +0.3036 mA/K

  % Parametros fisicos das 3 subcelulas em serie (GaInP / GaAs / Ge) a 28 C:
  % [1: Top GaInP, 2: Middle GaAs, 3: Bottom Ge]
  Voc_sub_ref = [1.420, 1.010, 0.290];        % [V] Soma = 2.720 V
  n_sub       = [1.15,  1.15,  1.75];         % Fatores de idealidade por subcelula
  Rsh_sub     = [2500.0, 3500.0, 200.0];      % [Ohm] Resistencias shunt individuais
  % Fotocorrentes individuais a 28 C (GaAs limitante, GaInP levemente acima, Ge com excesso IR)
  Iph_sub_ref = [1.012 * IscRef, 1.000 * IscRef, 1.550 * IscRef]; % [A]
  % Reparticao dos coeficientes termicos de Voc e Isc pelas 3 subcelulas
  dVoc_sub_dT = [-2.20e-3, -1.90e-3, -1.50e-3]; % [V/K] Soma = -5.60 mV/K
  dIph_sub_dT = [0.85 * dIsc_dT, 1.00 * dIsc_dT, 1.80 * dIsc_dT]; % [A/K]

  % Resistencia serie total (contatos metalicos + juncoes tunel)
  Rs_total = 0.20;           % [Ohm/celula]

  % Arquitetura da matriz por no satelital (166S x 244P = 40.504 celulas)
  Ns = 166;
  Np = 244;

  Pmp_cell_cat = VmpRef_cat * ImpRef;
  eta_cell_cat = Pmp_cell_cat / (G_AM0 * Am2);
  fprintf('--- 1.1 MODELO DE 3 DIODOS EM SERIE (GaInP / GaAs / Ge) ---\n');
  fprintf('Subcelulas Voc(28 C): GaInP=%.2f V | GaAs=%.2f V | Ge=%.2f V (Total=%.3f V)\n', ...
          Voc_sub_ref(1), Voc_sub_ref(2), Voc_sub_ref(3), sum(Voc_sub_ref));
  fprintf('Fatores de idealidade: n = [%.2f, %.2f, %.2f] | Rsh = [%.0f, %.0f, %.0f] Ohm | Rs = %.2f Ohm\n', ...
          n_sub(1), n_sub(2), n_sub(3), Rsh_sub(1), Rsh_sub(2), Rsh_sub(3), Rs_total);
  fprintf('Matriz: %dS x %dP (%d celulas) | Pmp catalogo (28 C) = %.2f kW (eta = %.2f %%)\n\n', ...
          Ns, Np, Ns * Np, (Ns * Np * Pmp_cell_cat) / 1e3, 100 * eta_cell_cat);

  % ----------------------------------------------------------------------------
  % 2. Resolucao das curvas I-V e P-V do modelo de 3 diodos por temperatura
  % ----------------------------------------------------------------------------
  fprintf('--- 2.1 CURVAS I-V E P-V DO MODELO DE 3 DIODOS EM SERIE ---\n');
  Tgrid = [15, 28, 45, 60, 85];
  VcellGrid = linspace(0, 2.85, 600);
  IcellGrid = zeros(numel(Tgrid), numel(VcellGrid));
  PmpNode = zeros(size(Tgrid));
  VmpNode = zeros(size(Tgrid));
  ImpNode = zeros(size(Tgrid));

  for m = 1:numel(Tgrid)
    T = Tgrid(m);
    dT = T - Tref;
    TK = T + 273.15;

    Iph_T = Iph_sub_ref + dIph_sub_dT * dT;
    Voc_T = Voc_sub_ref + dVoc_sub_dT * dT;
    Vth_T = (n_sub * kB * TK) / qe;
    I0_T  = (Iph_T - Voc_T ./ Rsh_sub) ./ (exp(Voc_T ./ Vth_T) - 1);

    % Parametrizacao por corrente I_cel de 0 ate min(Iph_T) para obter V_cel(I)
    I_sweep = linspace(0, min(Iph_T) * 0.9999, 1000);
    V_sweep = zeros(size(I_sweep));
    for k_i = 1:numel(I_sweep)
      ic = I_sweep(k_i);
      v_sum = 0.0;
      for j_sub = 1:3
        v_j = solve_subcell_voltage(ic, Iph_T(j_sub), I0_T(j_sub), Rsh_sub(j_sub), Vth_T(j_sub), Voc_T(j_sub));
        v_sum = v_sum + v_j;
      end
      V_sweep(k_i) = max(0.0, v_sum - ic * Rs_total);
    end

    % Interpolacao sobre a malha uniforme de tensao VcellGrid
    [V_uniq, idx_u] = unique(V_sweep, 'stable');
    I_uniq = I_sweep(idx_u);
    IcellGrid(m, :) = interp1(V_uniq, I_uniq, VcellGrid, 'linear', 0.0);

    Pcell = VcellGrid .* IcellGrid(m, :);
    [pmax, idx] = max(Pcell);
    PmpNode(m) = Ns * Np * pmax;
    VmpNode(m) = Ns * VcellGrid(idx);
    ImpNode(m) = Np * IcellGrid(m, idx);
    fprintf('T = %3.0f C | Vmp_cel = %.3f V | Vmp_str = %5.1f V | Imp_no = %6.2f A | Pmp_no = %5.2f kW\n', ...
            T, VcellGrid(idx), VmpNode(m), ImpNode(m), PmpNode(m) / 1e3);
  end
  fprintf('\n');

  % ----------------------------------------------------------------------------
  % 3. Cenario radiativo nominal (1 MeV eq., CMG 150 um, SSO 700 km)
  % ----------------------------------------------------------------------------
  fprintf('--- 3.1 DEGRADACAO DDD (EL. 1 MeV SOB COVERGLASS CMG 150 um EM SSO 700 km) ---\n');
  years = 0:0.1:10;
  Cd = 0.085;                % Eq. 34 da tese
  PhiX = 1.0e14;             % [e-/cm^2]
  PhiAnnual = 2.5e13;        % [e-/cm^2/ano] Fluencia eq. 1 MeV atras de CMG 150 um
  lossEnv10 = 0.052;         % Escurecimento UV + micrometeoroides (5.2% em 10 anos)
  Phi = years * PhiAnnual;
  etaRad = 1 - Cd * log10(1 + Phi / PhiX);
  etaEnv = 1 - lossEnv10 * years / 10;
  etaTotal = etaRad .* etaEnv;
  fprintf('Retencao DDD (10 anos) = %.2f %% | Retencao Ambiental = %.2f %% | Retencao Total EOL = %.2f %%\n', ...
          100 * etaRad(end), 100 * etaEnv(end), 100 * etaTotal(end));
  fprintf('Potencia EOL do no a 28 C (Modelo 3 Diodos) = %.2f kW\n\n', (PmpNode(2) / 1e3) * etaTotal(end));

  % ----------------------------------------------------------------------------
  % 4. Contingencia de Hot Spot no vacuo (Stefan-Boltzmann)
  % ----------------------------------------------------------------------------
  fprintf('--- 4.1 CONTINGENCIA DE HOT SPOT NO VACUO (STEFAN-BOLTZMANN) ---\n');
  Vbreak = 28.0; Vsch = 0.35; Istring = ImpRef;
  Psem = Vbreak * Istring; Pcom = Vsch * Istring;
  epsCell = 0.85; Cth = 3.8; Tspace = 3.0;
  for nFaces = [1, 2]
    dt_hs = 0.01; t_hs = 0:dt_hs:240;
    Tsem = zeros(size(t_hs)); Tsem(1) = 60 + 273.15;
    for k = 1:numel(t_hs)-1
      qsem = nFaces * epsCell * sigmaSB * Am2 * (Tsem(k)^4 - Tspace^4);
      Tsem(k+1) = Tsem(k) + (Psem - qsem) * dt_hs / Cth;
    end
    idx180 = find(Tsem >= 180 + 273.15, 1);
    TeqSem = (Psem / (nFaces * epsCell * sigmaSB * Am2) + Tspace^4)^0.25 - 273.15;
    TeqCom = (Pcom / (nFaces * epsCell * sigmaSB * Am2) + Tspace^4)^0.25 - 273.15;
    fprintf('nFaces = %d: P_sem = %.2f W | t(180 C) = %.2f s | Teq_sem = %.2f C | Teq_com = %.2f C\n', ...
            nFaces, Psem, t_hs(idx180), TeqSem, TeqCom);
  end
  fprintf('\n');

  % ----------------------------------------------------------------------------
  % 5. MALHA DINAMICA DO BUCK (150 kHz) E MPPT ADAPTATIVO P&O-AP (10 kHz)
  % ----------------------------------------------------------------------------
  fprintf('--- 5.1 DINAMICA DO CONVERSOR BUCK E MPPT ADAPTATIVO P&O-AP (150 kHz, T = 45 C) ---\n');
  V_bus_nom = 345.0;         % [V] Barramento pos-eclipse (janela 324-442.8 V do BESS 108S)
  L_buck    = 45.0e-6;       % [H] Indutancia do Buck (45 uH)
  C_in_buck = 250.0e-6;      % [F] Capacitancia de entrada (250 uF)
  R_L_buck  = 0.008;         % [Ohm] Resistencia parasita serie do indutor

  f_sw   = 150e3;            % [Hz] Frequencia de chaveamento PWM (150 kHz)
  dt_sim = 1 / f_sw;         % [s] Passo de integracao no periodo medio PWM (~6.667 us)
  t_end  = 0.050;            % [s] Janela transitoria de 50 ms
  t_mppt = 0:dt_sim:t_end;
  N_steps = numel(t_mppt);

  f_mppt = 10e3;             % [Hz] Frequencia do algoritmo MPPT (10 kHz)
  N_sync = round(f_sw / f_mppt); % 15 ciclos PWM

  M_adapt    = 1.2e-4;       % [V/W] Ganho adaptativo (Eq. eq:po_adaptativo)
  dD_max     = 0.002;        % Limite de saturacao de passo por amostragem
  delta_fixo = 0.005;        % Passo fixo comparativo

  % Curva de 3 Diodos em T = 45 C (indice 3: Vmp = 373.4 V, Pmp = 48.94 kW)
  V_curve_45 = Ns * VcellGrid;
  I_curve_45 = Np * IcellGrid(3, :);

  V_pv_a = zeros(1, N_steps); I_L_a = zeros(1, N_steps); D_a = zeros(1, N_steps); P_pv_a = zeros(1, N_steps);
  V_pv_f = zeros(1, N_steps); I_L_f = zeros(1, N_steps); D_f = zeros(1, N_steps); P_pv_f = zeros(1, N_steps);

  % Condicao inicial fora do MPP (partida em tensao elevada pos-eclipse: 410 V)
  V_pv_a(1) = 410.0; I_L_a(1) = 125.0; D_a(1) = V_bus_nom / V_pv_a(1);
  V_pv_f(1) = 410.0; I_L_f(1) = 125.0; D_f(1) = D_a(1);

  v_prev_mppt_a = V_pv_a(1);
  p_prev_mppt_a = V_pv_a(1) * interp1(V_curve_45, I_curve_45, V_pv_a(1), 'linear', 0.0);
  v_prev_mppt_f = V_pv_f(1);
  p_prev_mppt_f = p_prev_mppt_a;

  for k = 1:N_steps-1
    % Corrente exata interpolada da curva de 3 Diodos (GaInP/GaAs/Ge)
    I_pv_now_a = interp1(V_curve_45, I_curve_45, V_pv_a(k), 'linear', 0.0);
    I_pv_now_f = interp1(V_curve_45, I_curve_45, V_pv_f(k), 'linear', 0.0);
    P_pv_a(k)  = V_pv_a(k) * I_pv_now_a;
    P_pv_f(k)  = V_pv_f(k) * I_pv_now_f;

    if mod(k, N_sync) == 0
      % 1. MPPT Adaptativo P&O-AP: D[k] = D[k-1] - M * (dP/dV)
      dV_a = V_pv_a(k) - v_prev_mppt_a;
      dP_a = P_pv_a(k) - p_prev_mppt_a;
      if abs(dV_a) > 1e-4
        slope_a = dP_a / dV_a;
        step_adapt = - min(dD_max, max(-dD_max, M_adapt * slope_a));
        D_a(k) = max(0.50, min(0.96, D_a(k) + step_adapt));
      end
      v_prev_mppt_a = V_pv_a(k);
      p_prev_mppt_a = P_pv_a(k);

      % 2. MPPT de Passo Fixo Convencional
      dV_f = V_pv_f(k) - v_prev_mppt_f;
      dP_f = P_pv_f(k) - p_prev_mppt_f;
      if abs(dV_f) > 1e-4
        D_f(k) = max(0.50, min(0.96, D_f(k) - delta_fixo * sign(dP_f / dV_f)));
      end
      v_prev_mppt_f = V_pv_f(k);
      p_prev_mppt_f = P_pv_f(k);
    end

    D_a(k+1) = D_a(k);
    D_f(k+1) = D_f(k);

    % Integracao Euler-Cromer das EDOs do filtro LC do Buck
    dVpv_dt_a = (I_pv_now_a - D_a(k) * I_L_a(k)) / C_in_buck;
    V_pv_a(k+1) = max(0, V_pv_a(k) + dVpv_dt_a * dt_sim);
    dIL_dt_a = (D_a(k) * V_pv_a(k+1) - V_bus_nom - R_L_buck * I_L_a(k)) / L_buck;
    I_L_a(k+1) = max(0, I_L_a(k) + dIL_dt_a * dt_sim);

    dVpv_dt_f = (I_pv_now_f - D_f(k) * I_L_f(k)) / C_in_buck;
    V_pv_f(k+1) = max(0, V_pv_f(k) + dVpv_dt_f * dt_sim);
    dIL_dt_f = (D_f(k) * V_pv_f(k+1) - V_bus_nom - R_L_buck * I_L_f(k)) / L_buck;
    I_L_f(k+1) = max(0, I_L_f(k) + dIL_dt_f * dt_sim);
  end
  P_pv_a(end) = V_pv_a(end) * interp1(V_curve_45, I_curve_45, V_pv_a(end), 'linear', 0.0);
  P_pv_f(end) = V_pv_f(end) * interp1(V_curve_45, I_curve_45, V_pv_f(end), 'linear', 0.0);

  idx_reg = (N_steps - 2000):N_steps;
  fprintf('Tensao media no painel P&O-AP (Vin): %.2f V (Vmp 3-Diodos a 45 C = %.1f V)\n', ...
          mean(V_pv_a(idx_reg)), VmpNode(3));
  fprintf('Corrente media no indutor (IL): %.2f A\n', mean(I_L_a(idx_reg)));
  fprintf('Potencia extraida no MPP (45 C): %.2f kW (Eficiencia MPPT P&O-AP: %.2f %% vs Fixo: %.2f %%)\n', ...
          mean(P_pv_a(idx_reg)) / 1e3, 100 * mean(P_pv_a(idx_reg)) / PmpNode(3), 100 * mean(P_pv_f(idx_reg)) / PmpNode(3));
  fprintf('==============================================================================\n');
end

function Vj = solve_subcell_voltage(I, Iph, I0, Rsh, Vth, Voc_init)
  Vj = max(0.0, Voc_init);
  for it = 1:30
    e = exp(min(500.0, Vj / Vth));
    f = Iph - I0 * (e - 1.0) - Vj / Rsh - I;
    df = -(I0 / Vth) * e - 1.0 / Rsh;
    dV = f / df;
    Vj = min(Voc_init, max(0.0, Vj - dV));
    if abs(dV) < 1e-9, break; end
  end
end
