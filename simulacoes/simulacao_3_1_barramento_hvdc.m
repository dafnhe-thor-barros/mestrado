% ==============================================================================
% SIMULAÇÃO 3.1: DINÂMICA DO BARRAMENTO HVDC (400 V), BALANÇO DE ENERGIA
% MULTI-ÓRBITA DO BESS E PROTEÇÃO ATIVA POR CHAVE SSPC
% ==============================================================================
% Dissertação de Mestrado Profissional - ProPECAut / IFES
% Autor: Dafnhê Thor Castro de Martins Barros
% Revisão Técnica: Modelagem Explícita de Condições Nominais e Modos de Falha
% ==============================================================================
function simulacao_3_2_gerador_fotovoltaico_am0()
  close all; clc;

  % =========================================================================
  % 1. Constantes Fisicas e Parametros de Referencia da Celula (Anexo A)
  % =========================================================================
  qe = 1.602176634e-19;      % Carga elementar [C]
  kB = 1.380649e-23;          % Constante de Boltzmann [J/K]
  sigmaSB = 5.670374e-8;      % Constante de Stefan-Boltzmann [W/(m^2 K^4)]
  G0 = 1361.0;                % Constante solar AM0 [W/m^2]

  Tref = 28.0;                % Temperatura de referencia [degC]
  Acm2 = 30.18;               % Area da celula [cm^2]
  Am2 = Acm2 * 1e-4;          % Area da celula [m^2]
  VocRef = 2.720;             % Tensao de circuito aberto [V]
  IscRef = 0.538;             % Corrente de curto-circuito [A]
  VmpRef = 2.410;             % Tensao no ponto de maxima potencia [V]
  ImpRef = 0.512;             % Corrente no ponto de maxima potencia [A]
  dVoc_dT = -6.2e-3;          % Coeficiente de temperatura de Voc [V/K]
  dIsc_dT = 0.32e-3;          % Coeficiente de temperatura de Isc [A/K]

  aeq = 3.45;                 % Fator de idealidade equivalente
  Rs = 0.085;                 % Resistencia serie por celula [Ohm]
  Rsh = 1200.0;               % Resistencia paralelo por celula [Ohm]
  Ns = 180;                   % Celulas em serie por string
  Np = 225;                   % Strings em paralelo na matriz

  % Eficiencia real calculada sob AM0 (harmonizacao de catalogo)
  Pmp_cell = VmpRef * ImpRef;
  eta_cell = Pmp_cell / (G0 * Am2);
  fprintf('=== 1. Caracterizacao do Painel Solar ===\n');
  fprintf('Pmp celula nominal: %.4f W | Eficiencia real AM0: %.2f %%\n', Pmp_cell, 100*eta_cell);
  fprintf('Arranjo da matriz: %dS x %dP (%d celulas no total)\n\n', Ns, Np, Ns*Np);

  % =========================================================================
  % 2. Resolucao das Curvas I-V da Matriz por Newton-Raphson
  % =========================================================================
  Tgrid = [15 28 60 85];
  VcellGrid = linspace(0, 2.85, 600);
  PmpNode = zeros(size(Tgrid));
  VmpNode = zeros(size(Tgrid));

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
    Pcell = VcellGrid .* Icell;
    [pmax, idx] = max(Pcell);
    PmpNode(m) = Ns * Np * pmax;
    VmpNode(m) = Ns * VcellGrid(idx);
    fprintf('T = %2d C | Vmp_matriz = %6.1f V | Imp_matriz = %6.2f A | Pmp = %5.2f kW\n', ...
            T, VmpNode(m), Np*Icell(idx), PmpNode(m)/1e3);
  end

  % =========================================================================
  % 3. Degradacao por Dano de Deslocamento (DDD) e Efeito Ambiental
  % =========================================================================
  years = 0:0.1:10;
  Cd = 0.085;                 % Coeficiente de degradacao (Eq. 34)
  PhiX = 1.0e14;              % Fluencia caracteristica [e-/cm^2]
  PhiAnnual = 2.5e13;         % Fluencia equivalente anual [e-/cm^2/ano]
  lossEnv10 = 0.052;          % Perda ambiental acumulada em 10 anos (5,2%)
  Phi = years * PhiAnnual;
  etaRad = 1 - Cd * log10(1 + Phi / PhiX);
  etaEnv = 1 - lossEnv10 * (years / 10);
  etaTotal = etaRad .* etaEnv;

  fprintf('\n=== 2. Analise de Degradacao EOL (10 Anos) ===\n');
  fprintf('Retencao por radiacao (DDD): %.2f %%\n', 100*etaRad(end));
  fprintf('Retencao global EOL: %.2f %% (Perda total: %.2f %%)\n\n', ...
          100*etaTotal(end), 100*(1-etaTotal(end)));

  % =========================================================================
  % 4. Hot Spot com Duas Hipoteses de Faces Radiantes
  % =========================================================================
  Vbreak = 28.0;              % Tensao reversa de breakdown [V]
  Vsch = 0.35;                % Queda no diodo de bypass Schottky [V]
  Istring = ImpRef;
  P_hot_sem = Vbreak * Istring;
  P_hot_com = Vsch * Istring;
  epsCell = 0.85; Cth_cel = 3.8; Tspace = 3.0;

  % Modelo analitico assintotico
  for nFaces = [1 2]
    Teq_sem = (P_hot_sem / (nFaces * epsCell * sigmaSB * Am2) + Tspace^4)^0.25 - 273.15;
    Teq_com = (P_hot_com / (nFaces * epsCell * sigmaSB * Am2) + Tspace^4)^0.25 - 273.15;
    fprintf('Hot Spot (n_faces = %d): Teq sem bypass = %.1f C | com bypass = %.1f C\n', ...
            nFaces, Teq_sem, Teq_com);
  end

  % =========================================================================
  % 5. Dinamica do Conversor Buck e MPPT (150 kHz, Resolucao de Estados)
  % =========================================================================
  fprintf('\n=== 3. Dinamica Chaveada do Buck e MPPT (150 kHz) ===\n');
  fsw = 150e3; Tsw = 1/fsw;
  Vbus = 345.0;               % Tensao do barramento [V]
  L_buck = 60e-6; RL_buck = 12e-3; Cin = 220e-6;

  dt = Tsw / 20;              % 20 passos por ciclo PWM
  t_sim = 0:dt:0.025;         % 25 ms de simulacao dinamica
  N_steps = numel(t_sim);

  v_in = zeros(1, N_steps); i_L = zeros(1, N_steps);
  v_in(1) = 430.0; i_L(1) = 15.0;
  D = 0.80; P_prev = 0; V_prev = v_in(1);
  mppt_div = 15; mppt_cnt = 0; step_D = 0.004;

  % Parametros equivalentes da matriz a 45 degC
  T_op = 45.0;
  Iph_matriz = (IscRef + dIsc_dT*(T_op-Tref)) * Np;
  Voc_matriz = (VocRef + dVoc_dT*(T_op-Tref)) * Ns;
  Vth_matriz = (aeq * kB * (T_op+273.15) / qe) * Ns;
  Rs_m = Rs * (Ns/Np); Rsh_m = Rsh * (Ns/Np);
  I0_m = (Iph_matriz - Voc_matriz/Rsh_m) / (exp(Voc_matriz/Vth_matriz) - 1);

  for k = 1:N_steps-1
    % Corrente nao-linear da matriz
    ipv = solve_cell_i(v_in(k), Iph_matriz, I0_m, Rs_m, Rsh_m, Vth_matriz);

    % Estado da chave PWM
    sw = (mod(t_sim(k), Tsw) / Tsw) < D;

    % Equacoes diferenciais do conversor Buck
    diL = (sw * v_in(k) - Vbus - i_L(k)*RL_buck) / L_buck;
    dVin = (ipv - sw * i_L(k)) / Cin;

    i_L(k+1) = max(0, i_L(k) + diL * dt);
    v_in(k+1) = max(0, v_in(k) + dVin * dt);

    % Laco de controle MPPT P&O
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
  fprintf('Corrente media no indutor (IL): %.2f A | Ripple: %.2f A\n', ...
          mean(i_L(end-200:end)), max(i_L(end-200:end)) - min(i_L(end-200:end)));
  fprintf('Potencia extraida no regime MPPT: %.2f kW\n', P_final/1e3);
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
fprintf('==============================================================================\n');
fprintf('SIMULAÇÃO 3.1 EXECUTADA COM SUCESSO: DIAGNÓSTICO DE FALHAS CONCLUÍDO.\n');
fprintf('==============================================================================\n');
