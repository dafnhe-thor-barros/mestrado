function simulacao_3_2_gerador_fotovoltaico_am0()
  % SIMULACAO 3.2 - GERADOR FOTOVOLTAICO AM0, DEGRADACAO E BUCK MPPT
  % Resolucao de modelo de 1 diodo implicito e integracao dinamica em 150 kHz
  close all; clc;

  % 1. Constantes Fisicas e Parametros de Referencia da Celula (Anexo A)
  qe = 1.602176634e-19;      % Carga elementar [C]
  kB = 1.380649e-23;          % Constante de Boltzmann [J/K]
  sigmaSB = 5.670374e-8;      % Stefan-Boltzmann [W/(m^2 K^4)]
  G0 = 1361.0;                % Irradiancia AM0 [W/m^2]

  Tref = 28.0;                % Temperatura de referencia [degC]
  Acm2 = 30.18;               % Area da celula [cm^2]
  Am2 = Acm2 * 1e-4;          % Area da celula [m^2]
  VocRef = 2.720;             % Tensao de circuito aberto [V]
  IscRef = 0.538;             % Corrente de curto-circuito [A]
  VmpRef = 2.410;             % Tensao no MPP [V]
  ImpRef = 0.512;             % Corrente no MPP [A]
  dVoc_dT = -6.2e-3;          % Coeficiente termico Voc [V/K]
  dIsc_dT = 0.32e-3;          % Coeficiente termico Isc [A/K]

  aeq = 3.45; Rs = 0.085; Rsh = 1200.0;
  Ns = 180; Np = 225;

  % Eficiencia real calculada em AMO
  Pmp_cell = VmpRef * ImpRef;
  eta_cell = Pmp_cell / (G0 * Am2);
  fprintf('=== 1. Caracterizacao do Painel Solar ===\n');
  fprintf('Pmp celula nominal: %.4f W | Eficiencia real AM0: %.2f %%\n', Pmp_cell, 100*eta_cell);
  fprintf('Matriz: %dS x %dP (%d celulas no total)\n\n', Ns, Np, Ns*Np);

  % 2. Resolucao das Curvas I-V por Newton-Raphson
  Tgrid = [15 28 60 85];
  VcellGrid = linspace(0, 2.85, 600);
  for m = 1:numel(Tgrid)
    T = Tgrid(m);
    dT = T - Tref;
    Iph = IscRef + dIsc_dT * dT;
    Voc = VocRef + dVoc_dT * dT;
    Vth = aeq * kB * (T + 273.15) / qe;
    I0 = (Iph - Voc/Rsh) / (exp(Voc/Vth) - 1);

    Icell = zeros(size(VcellGrid));
    for j = 1:numel(VcellGrid)
      Icell(j) = solve_cell_i(VcellGrid(j), Iph, I0, Rs, Rsh, Vth);
    end
    [pmax, idx] = max(VcellGrid .* Icell);
    fprintf('T = %2d C | Vmp_matriz = %6.1f V | Imp_matriz = %6.2f A | Pmp = %5.2f kW\n', ...
            T, Ns*VcellGrid(idx), Np*Icell(idx), (Ns*Np*pmax)/1e3);
  end

  % 3. Degradacao por Dano de Deslocamento (DDD) e Perdas Ambientais
  years = 0:0.1:10;
  Cd = 0.085; PhiX = 1.0e14; PhiAnnual = 2.5e13; lossEnv10 = 0.052;
  Phi = years * PhiAnnual;
  etaRad = 1 - Cd * log10(1 + Phi / PhiX);
  etaEnv = 1 - lossEnv10 * (years / 10);
  etaTotal = etaRad .* etaEnv;

  fprintf('\n=== 2. Analise de Degradacao EOL (10 Anos) ===\n');
  fprintf('Retencao por radiacao (DDD): %.2f %%\n', 100*etaRad(end));
  fprintf('Retencao global EOL: %.2f %% (Perda total: %.2f %%)\n\n', ...
          100*etaTotal(end), 100*(1-etaTotal(end)));

  % 4. Hot Spot: Modelagem Analitica
  Vbreak = 28.0; Vsch = 0.35; Istring = ImpRef;
  P_hot_sem = Vbreak * Istring; P_hot_com = Vsch * Istring;
  epsCell = 0.85; Tspace = 3.0;

  for nFaces = [1 2]
    Teq_sem = (P_hot_sem / (nFaces * epsCell * sigmaSB * Am2) + Tspace^4)^0.25 - 273.15;
    Teq_com = (P_hot_com / (nFaces * epsCell * sigmaSB * Am2) + Tspace^4)^0.25 - 273.15;
    fprintf('Hot Spot (n_faces = %d): Teq sem bypass = %.1f C | com bypass = %.1f C\n', ...
            nFaces, Teq_sem, Teq_com);
  end

  % 5. Dinamica Chaveada do Conversor Buck e MPPT (150 kHz)
  fprintf('\n=== 3. Dinamica do Conversor Buck e MPPT (150 kHz) ===\n');
  fsw = 150e3; Tsw = 1/fsw; Vbus = 345.0;
  L_buck = 60e-6; RL_buck = 12e-3; Cin = 220e-6;
  dt = Tsw / 20; t_sim = 0:dt:0.025; N_steps = numel(t_sim);
  
  v_in = zeros(1, N_steps); i_L = zeros(1, N_steps);
  v_in(1) = 430.0; i_L(1) = 15.0;
  D = 0.80; P_prev = 0; V_prev = v_in(1);
  mppt_div = 15; mppt_cnt = 0; step_D = 0.004;

  T_op = 45.0;
  Iph_matriz = (IscRef + dIsc_dT*(T_op-Tref)) * Np;
  Voc_matriz = (VocRef + dVoc_dT*(T_op-Tref)) * Ns;
  Vth_matriz = (aeq * kB * (T_op+273.15) / qe) * Ns;
  Rs_m = Rs * (Ns/Np); Rsh_m = Rsh * (Ns/Np);
  I0_m = (Iph_matriz - Voc_matriz/Rsh_m) / (exp(Voc_matriz/Vth_matriz) - 1);

  for k = 1:N_steps-1
    ipv = solve_cell_i(v_in(k), Iph_matriz, I0_m, Rs_m, Rsh_m, Vth_matriz);
    sw = (mod(t_sim(k), Tsw) / Tsw) < D;

    diL = (sw * v_in(k) - Vbus - i_L(k)*RL_buck) / L_buck;
    dVin = (ipv - sw * i_L(k)) / Cin;

    i_L(k+1) = max(0, i_L(k) + diL * dt);
    v_in(k+1) = max(0, v_in(k) + dVin * dt);

    mppt_cnt = mppt_cnt + 1;
    if mppt_cnt >= mppt_div
      mppt_cnt = 0;
      P_now = v_in(k+1) * ipv;
      dP = P_now - P_prev; dV = v_in(k+1) - V_prev;
      if dP > 0
        if dV < 0, D = D + step_D; else, D = D - step_D; end
      else
        if dV < 0, D = D - step_D; else, D = D + step_D; end
      end
      D = min(0.95, max(0.45, D));
      P_prev = P_now; V_prev = v_in(k+1);
    end
  end

  P_final = mean(v_in(end-200:end)) * mean(i_L(end-200:end)) * (Vbus/mean(v_in(end-200:end)));
  fprintf('Tensao final no painel (Vin): %.2f V\n', mean(v_in(end-200:end)));
  fprintf('Corrente media no indutor (IL): %.2f A (Ripple: %.2f A)\n', ...
          mean(i_L(end-200:end)), max(i_L(end-200:end)) - min(i_L(end-200:end)));
  fprintf('Potencia gerada em regime: %.2f kW\n', P_final/1e3);
end

function I = solve_cell_i(V, Iph, I0, Rs, Rsh, Vth)
  I = Iph;
  for it = 1:40
    e = exp(min(60, (V + I*Rs)/Vth));
    f = Iph - I0*(e - 1) - (V + I*Rs)/Rsh - I;
    df = -I0*(Rs/Vth)*e - Rs/Rsh - 1;
    step = f / df;
    I = I - step;
    if abs(step) < 1e-7, break; end
  end
  I = max(0, I);
end
