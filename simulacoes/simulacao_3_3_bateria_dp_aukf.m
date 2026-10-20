function simulacao_3_3_bateria_dp_aukf()
  % ============================================================================
  % SIMULACAO 3.3: BESS, MODELO DUAL-POLARIZATION (2a ORDEM), REDE TERMICA DE
  % DOIS NOS (BATERIA-CHASSI-RADIADOR, NASA NTRS) E AUKF SAGE-HUSA
  % ============================================================================
  % Dissertacao de Mestrado Profissional - ProPECAut / IFES
  % Autor: Dafnhe Thor Castro de Martins Barros
  % Referencias Verificaveis:
  %   [1] NASA NTRS 20205006303: Practical Battery Thermal Modeling Techniques
  %   [2] NASA NTRS 20090015884: ORS Standard Bus Battery Thermal Balance
  % ============================================================================
  close all; clc;

  fprintf('==============================================================================\n');
  fprintf('INICIANDO SIMULACAO 3.3: MODELO DP, REDE TERMICA DE 2 NOS E AUKF SAGE-HUSA\n');
  fprintf('==============================================================================\n\n');

  % ----------------------------------------------------------------------------
  % 1. Parametros Eletroquimicos Sinteticos Parametrizados (108S Li-ion NMC)
  % ----------------------------------------------------------------------------
  Cn = 22.52;                  % [Ah] Capacidade equivalente [PROVISORIO]
  Vuv = 324.0; Vov = 442.8;    % [V] Limites UVLO (3.0 V/cel) e OV (4.1 V/cel)
  Rp1 = 18e-3; Cp1 = 3200.0;   % [ohm, F] tau1 = 57.6 s [PROVISORIO]
  Rp2 = 24e-3; Cp2 = 22000.0;  % [ohm, F] tau2 = 528.0 s [PROVISORIO]
  R0bol20 = 45e-3;             % [ohm] BOL a 20 C [PROVISORIO]
  R0eol20 = 85e-3;             % [ohm] EOL a 20 C [PROVISORIO]
  Ea = 0.32 * 1.602176634e-19; % [J] Energia de ativacao (0.32 eV) [PROVISORIO]
  kB = 1.380649e-23;
  sigmaSB = 5.670374e-8;
  Tref = 293.15;

  % ----------------------------------------------------------------------------
  % 2. Parametros da Rede Termica de 2 Nos (Bateria + Chassi Radiativo)
  % ----------------------------------------------------------------------------
  C_bat    = 54360.0;          % [J/K] Capacidade termica do pack (60 kg NMC)
  C_chassi = 38000.0;          % [J/K] Capacidade termica do chassi de montagem [HIPOTESE]
  R_bat_mli    = 4.80;         % [K/W] Isolamento radiativo da manta MLI [HIPOTESE]
  R_bat_chassi = 0.25;         % [K/W] Interface condutiva pack-chassi [HIPOTESE]
  T_amb_space  = -40.0;        % [degC] Temperatura equivalente do entorno MLI
  T_space_K    = 3.0;          % [K] Fundo cosmico
  A_rad_bat    = 0.22;         % [m^2] Area radiativa dedicada do chassi da bateria
  eps_rad_bat  = 0.88;         % Emissividade infravermelha da face radiativa
  q_amb_bat    = 52.8;         % [W/m^2] Fluxo ambiental absorvido em SSO

  Pheater = 300.0;             % [W] Aquecedor de sobrevivencia (banda 10 C - 15 C)
  Ton = 10.0; Toff = 15.0;

  dt = 2.0; Torb = 5926.0; tecl = 1128.0; nOrb = 5;
  N = floor(nOrb * Torb / dt) + 1;

  % Simulacao comparativa:
  % Caso 1: 1 No isolado apenas por R_th,MLI = 4.8 K/W (mostra aquecimento continuo)
  % Caso 2: 2 Nos acoplados (Bateria + Chassi Radiativo com R_bat,chassi = 0.25 K/W)
  for modo_termico = [1, 2]
    soc = zeros(1, N); vt = zeros(1, N);
    Tb = zeros(1, N); Tch = zeros(1, N);
    soc(1) = 0.90; Tb(1) = 20.0; Tch(1) = 15.0;
    vp1 = 0.0; vp2 = 0.0; heaterOn = false;
    a1 = exp(-dt / (Rp1 * Cp1)); a2 = exp(-dt / (Rp2 * Cp2));
    Itrace = zeros(1, N); Pjtrace = zeros(1, N);

    for k = 1:N-1
      tk = (k - 1) * dt;
      phase = mod(tk, Torb);
      if phase <= tecl
        I = (18000.0 + 560.0) / 380.0; % [A] Descarga no eclipse (48.84 A)
      else
        I = -25.0;                     % [A] Recarga solar (-25.0 A)
      end
      R0 = R0eol20 * exp(Ea / kB * (1.0 / (Tb(k) + 273.15) - 1.0 / Tref));
      vp1 = a1 * vp1 + Rp1 * (1.0 - a1) * I;
      vp2 = a2 * vp2 + Rp2 * (1.0 - a2) * I;
      etaC = 1.0; if I < 0, etaC = 0.995; end
      soc(k+1) = min(1.0, max(0.0, soc(k) - I * etaC * dt / (Cn * 3600.0)));

      Voc = 108.0 * (3.0 + 1.1 * soc(k+1));
      Vcand = Voc - I * R0 - vp1 - vp2;
      if I < 0 && Vcand > Vov
        Icv = -(Vov - Voc + vp1 + vp2) / R0;
        I = max(-25.0, min(0.0, Icv));
        Vcand = Vov;
      end
      vt(k+1) = Vcand;

      % Geracao termica irreversivel (Joule + polarizacao eletroquimica)
      Pj = (I^2) * R0 + (vp1^2) / Rp1 + (vp2^2) / Rp2;
      if Tb(k) <= Ton, heaterOn = true; elseif Tb(k) >= Toff, heaterOn = false; end
      Ph = Pheater * double(heaterOn);

      if modo_termico == 1
        % Modelo simplificado de 1 no (apenas R_th,MLI = 4.8 K/W para chassi fixo a 15 C)
        qout = (Tb(k) - 15.0) / R_bat_mli;
        Tb(k+1) = Tb(k) + (Pj + Ph - qout) * dt / C_bat;
        Tch(k+1) = 15.0;
      else
        % Modelo de 2 nos (Bateria <-> Chassi <-> Radiador Espacial)
        q_bat_ch  = (Tb(k) - Tch(k)) / R_bat_chassi;
        q_bat_mli = (Tb(k) - T_amb_space) / R_bat_mli;
        Q_ch_rad  = A_rad_bat * (eps_rad_bat * sigmaSB * ((Tch(k) + 273.15)^4 - T_space_K^4) - q_amb_bat);

        Tb(k+1)  = Tb(k)  + (Pj + Ph - q_bat_ch - q_bat_mli) * dt / C_bat;
        Tch(k+1) = Tch(k) + (q_bat_ch - Q_ch_rad) * dt / C_chassi;
      end
      Itrace(k) = I; Pjtrace(k) = Pj;
    end
    vt(1) = vt(2); Itrace(end) = Itrace(end-1); Pjtrace(end) = Pjtrace(end-1);

    idx_ecl = find(mod((0:N-1)*dt, Torb) <= tecl);
    idx_sol = find(mod((0:N-1)*dt, Torb) > tecl);
    fprintf('--- 1.%d BALANCO MULTI-ORBITA (5 ORBITAS) - MODO TERMICO %d ---\n', modo_termico, modo_termico);
    fprintf('SoC min = %.2f %% | Vt min/max = %.2f / %.2f V\n', 100 * min(soc), min(vt), max(vt));
    fprintf('Dissipacao Joule: Eclipse = %.1f W | Sol (CC-CV) = %.1f W | Media Orbital = %.2f W\n', ...
            mean(Pjtrace(idx_ecl)), mean(Pjtrace(idx_sol)), mean(Pjtrace));
    fprintf('T_bat min/max (5 orbitas) = %.2f / %.2f C | T_chassi final = %.2f C\n\n', ...
            min(Tb), max(Tb), Tch(end));
  end

  % ----------------------------------------------------------------------------
  % 3. Estimador AUKF com Adaptacao de Sage-Husa (15 mVrms vs 150 mVrms)
  % ----------------------------------------------------------------------------
  fprintf('--- 2.1 VALIDACAO IN SILICO DO AUKF SAGE-HUSA (15 mV vs 150 mV) ---\n');
  for sigmaV = [0.015, 0.150]
    m = aukf_demo(Cn, Rp1, Rp2, R0eol20, sigmaV, true);
    fprintf('Ruido sigmaV = %3.0f mVrms | tconv = %2g s | Erro SoC medio/max = %.4f / %.4f %% | RMSE V = %.2f mV\n', ...
            sigmaV * 1e3, m.tconv, m.errMean, m.errMax, m.rmseClean * 1e3);
  end
  fprintf('==============================================================================\n');
end

function M = aukf_demo(Cn, Rp1, Rp2, R0, sigmaV, adaptR)
  rng(42, 'twister');
  dt = 1; N = 1201; L = 3; b = 0.95;
  socTrue = zeros(1, N); v1True = zeros(1, N); v2True = zeros(1, N);
  vClean = zeros(1, N); vMeas = zeros(1, N); socTrue(1) = 0.85;
  current = @(t) 25 + 15 * sin(2 * pi * t / 180) + 10 * (mod(t, 120) < 30);
  for k = 1:N-1
    I = current(k-1);
    socTrue(k+1) = socTrue(k) - I * dt / (Cn * 3600);
    v1True(k+1) = exp(-dt / 57.6) * v1True(k) + Rp1 * (1 - exp(-dt / 57.6)) * I;
    v2True(k+1) = exp(-dt / 528.0) * v2True(k) + Rp2 * (1 - exp(-dt / 528.0)) * I;
    vClean(k+1) = 108 * (3 + 1.1 * socTrue(k+1)) - I * R0 - v1True(k+1) - v2True(k+1);
    vMeas(k+1) = vClean(k+1) + sigmaV * randn();
  end
  vClean(1) = vClean(2); vMeas(1) = vMeas(2);

  alpha = 1e-3; beta = 2; kappa = 0;
  lambda = alpha^2 * (L + kappa) - L;
  Wm = ones(1, 2 * L + 1) / (2 * (L + lambda));
  Wc = Wm; Wm(1) = lambda / (L + lambda);
  Wc(1) = Wm(1) + (1 - alpha^2 + beta);

  x = [0.60; 0; 0]; P = diag([0.05, 0.01, 0.01]);
  Q = diag([1e-7, 1e-5, 1e-5]); R = sigmaV^2;
  socHat = zeros(1, N); vHat = zeros(1, N); socHat(1) = x(1);
  for k = 1:N-1
    I = current(k-1);
    Ps = (P + P') / 2;
    [V, D] = eig(Ps);
    D = diag(max(diag(D), 1e-12));
    Ps = V * D * V';
    S = chol((L + lambda) * Ps, 'lower');
    X = [x, x + S, x - S];
    Xp = zeros(size(X));
    for j = 1:2*L+1
      Xp(1, j) = X(1, j) - I * dt / (Cn * 3600);
      Xp(2, j) = exp(-dt / 57.6) * X(2, j) + Rp1 * (1 - exp(-dt / 57.6)) * I;
      Xp(3, j) = exp(-dt / 528.0) * X(3, j) + Rp2 * (1 - exp(-dt / 528.0)) * I;
    end
    xp = Xp * Wm';
    Pp = Q;
    for j = 1:2*L+1
      dx = Xp(:, j) - xp; Pp = Pp + Wc(j) * (dx * dx');
    end
    Y = 108 * (3 + 1.1 * Xp(1, :)) - I * R0 - Xp(2, :) - Xp(3, :);
    yp = Y * Wm'; Pyy0 = 0; Pxy = zeros(L, 1);
    for j = 1:2*L+1
      dy = Y(j) - yp; dx = Xp(:, j) - xp;
      Pyy0 = Pyy0 + Wc(j) * dy^2; Pxy = Pxy + Wc(j) * dx * dy;
    end
    innov = vMeas(k+1) - yp;
    if adaptR
      dk = (1 - b) / (1 - b^(k+1));
      R = max(1e-8, (1 - dk) * R + dk * (innov^2 - Pyy0));
    end
    Pyy = Pyy0 + R; K = Pxy / Pyy;
    x = xp + K * innov; P = Pp - K * Pyy * K';
    x(1) = min(1, max(0, x(1)));
    socHat(k+1) = x(1);
    vHat(k+1) = 108 * (3 + 1.1 * x(1)) - I * R0 - x(2) - x(3);
  end
  err = abs(socHat - socTrue) * 100;
  idx = find(err > 1, 1, 'last');
  if isempty(idx), M.tconv = 0; else M.tconv = idx - 1; end
  M.errMean = mean(err(61:end)); M.errMax = max(err(61:end));
  M.rmseClean = sqrt(mean((vHat(61:end) - vClean(61:end)).^2));
end
