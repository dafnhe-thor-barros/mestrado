% ==============================================================================
% SIMULACAO 3.3: MODELAGEM ELETROTERMICA DO BANCO DE BATERIAS (BESS),
% CIRCUITO DUAL POLARIZATION (DP 2a ORDEM) E ESTIMACAO DE SOC VIA AUKF
% ==============================================================================
% Dissertacao de Mestrado Profissional - ProPECAut / IFES
% Autor: Dafnhe Thor Castro de Martins Barros
% Revisao Tecnica: Modelagem Explicita de Condicoes Nominais e Modos de Falha
% ==============================================================================
function simulacao_3_3_bateria_dp_aukf()
  close all; clc;

  % 1. Parametros Eletroquimicos e Termicos da Bateria
  Cn = 22.52;                 % Capacidade nominal equivalente [Ah]
  Vuv = 324.0; Vov = 442.8;   % Limites de tensao [V]
  Rp1 = 18e-3; Cp1 = 3200;    % R-C difusao rapida
  Rp2 = 24e-3; Cp2 = 22000;   % R-C difusao lenta
  R0_ref = 85e-3;             % Resistencia EOL a 20 C [Ohm]
  Ea = 0.32 * 1.602176634e-19;% Energia de ativacao [J]
  kB = 1.380649e-23;
  Tref_K = 293.15;            % 20 C em K

  % Acoplamento condutivo corrigido via barras de aluminio e heat pipes
  Cth = 54360;                % Capacitancia termica [J/K]
  Rth_cond = 0.14;            % Resistencia termica condutiva [K/W] (corrigida de 4.8)
  Tch = 15.0;                 % Temperatura de chassi [degC]
  Pheater = 300.0; Ton = 10.0; Toff = 15.0;

  dt = 2.0; Torb = 5926; tecl = 1128; nOrb = 5;
  N = floor(nOrb * Torb / dt) + 1;

  soc = zeros(1, N); vt = zeros(1, N); Tb = zeros(1, N); Pj = zeros(1, N);
  soc(1) = 0.90; Tb(1) = 20.0;
  vp1 = 0; vp2 = 0;
  a1 = exp(-dt / (Rp1 * Cp1)); a2 = exp(-dt / (Rp2 * Cp2));

  % 2. Simulacao Multiorbital Fechada
  for k = 1:N-1
    t_now = (k-1) * dt;
    in_eclipse = mod(t_now, Torb) <= tecl;
    
    if in_eclipse
      I = 18500 / 400;         % Descarga no eclipse (~46.3 A)
      etaC = 1.0;
    else
      I = -11.8;               % Recarga solar CC-CV [A]
      etaC = 0.995;
    end

    R0 = R0_ref * exp((Ea/kB) * (1/(Tb(k) + 273.15) - 1/Tref_K));
    vp1 = a1 * vp1 + Rp1 * (1 - a1) * I;
    vp2 = a2 * vp2 + Rp2 * (1 - a2) * I;

    soc(k+1) = min(1.0, max(0.0, soc(k) - I * etaC * dt / (Cn * 3600)));
    Voc = 108 * (3.0 + 1.1 * soc(k+1));
    Vcand = Voc - I * R0 - vp1 - vp2;

    if I < 0 && Vcand > Vov
      I = max(-11.8, min(0.0, -(Vov - Voc + vp1 + vp2) / R0));
      Vcand = Vov;
    end
    vt(k+1) = Vcand;

    Pj(k) = I^2 * R0 + vp1^2/Rp1 + vp2^2/Rp2;
    if Tb(k) < Ton
      Ph = Pheater;
    elseif Tb(k) < Toff
      Ph = Pheater * (Toff - Tb(k)) / (Toff - Ton);
    else
      Ph = 0;
    end
    q_out = (Tb(k) - Tch) / Rth_cond;
    Tb(k+1) = Tb(k) + (Pj(k) + Ph - q_out) * dt / Cth;
  end
  vt(1) = vt(2); Pj(end) = Pj(end-1);

  fprintf('=== 1. Balanco Eletrotermico Multiorbital (5 Orbitas) ===\n');
  fprintf('SoC Minimo: %.2f %% | Tensao Min/Max: %.1f V / %.1f V\n', ...
          100*min(soc), min(vt), max(vt));
  fprintf('Temperatura Inicial/Final: %.1f C -> %.1f C (Pico: %.1f C)\n', ...
          Tb(1), Tb(end), max(Tb));
  fprintf('Perda Joule Media: %.1f W | Equilibrio Assintotico Teq: %.1f C\n', ...
          mean(Pj), Tch + mean(Pj)*Rth_cond);
  fprintf('Energia Eclipse: %.2f kWh (Util a 60%% DoD: %.2f kWh)\n\n', ...
          (18560*tecl/3.6e6), 9.0*0.60);

  % 3. Estimador AUKF com Sage-Husa (sigma_v = 150 mV da tese)
  sigma_v = 0.150;
  res_aukf = executa_aukf(Cn, Rp1, Rp2, R0_ref, sigma_v);
  fprintf('=== 2. Desempenho do Estimador AUKF ===\n');
  fprintf('Tempo de convergencia (erro < 1%%): %d s\n', res_aukf.tconv);
  fprintf('Erro medio de SoC: %.4f %% | Erro maximo: %.4f %%\n', ...
          res_aukf.err_mean, res_aukf.err_max);
  fprintf('RMSE contra tensao limpa: %.2f mV\n', res_aukf.rmse_mV);
end

function res = executa_aukf(Cn, Rp1, Rp2, R0, sigmaV)
  randn('seed', 42);
  dt = 1.0; N = 1201; L = 3; b = 0.95;
  socT = zeros(1, N); v1T = zeros(1, N); v2T = zeros(1, N);
  vClean = zeros(1, N); vMeas = zeros(1, N); socT(1) = 0.85;

  prof_i = @(t) 25 + 15*sin(2*pi*t/180) + 10*(mod(t, 120) < 30);
  for k = 1:N-1
    i_curr = prof_i(k-1);
    socT(k+1) = socT(k) - i_curr * dt / (Cn * 3600);
    v1T(k+1) = exp(-dt/57.6)*v1T(k) + Rp1*(1 - exp(-dt/57.6))*i_curr;
    v2T(k+1) = exp(-dt/528)*v2T(k) + Rp2*(1 - exp(-dt/528))*i_curr;
    vClean(k+1) = 108*(3 + 1.1*socT(k+1)) - i_curr*R0 - v1T(k+1) - v2T(k+1);
    vMeas(k+1) = vClean(k+1) + sigmaV * randn();
  end
  vClean(1) = vClean(2); vMeas(1) = vMeas(2);

  alpha = 1e-3; beta = 2; kappa = 0;
  lambda = alpha^2 * (L + kappa) - L;
  Wm = ones(1, 2*L+1) / (2*(L + lambda)); Wc = Wm;
  Wm(1) = lambda / (L + lambda); Wc(1) = Wm(1) + (1 - alpha^2 + beta);

  x = [0.60; 0; 0];
  P = diag([0.05, 0.01, 0.01]); Q = diag([1e-7, 1e-5, 1e-5]); R = sigmaV^2;
  socHat = zeros(1, N); vHat = zeros(1, N); socHat(1) = x(1);

  for k = 1:N-1
    i_curr = prof_i(k-1);
    Ps = 0.5 * (P + P');
    [V, D_mat] = eig(Ps);
    Ps = V * diag(max(diag(D_mat), 1e-12)) * V';
    S = chol((L + lambda) * Ps, 'lower');
    X = [x, x + S, x - S];

    Xp = zeros(size(X));
    for j = 1:2*L+1
      Xp(1,j) = X(1,j) - i_curr * dt / (Cn * 3600);
      Xp(2,j) = exp(-dt/57.6)*X(2,j) + Rp1*(1 - exp(-dt/57.6))*i_curr;
      Xp(3,j) = exp(-dt/528)*X(3,j) + Rp2*(1 - exp(-dt/528))*i_curr;
    end
    xp = Xp * Wm';
    Pp = Q;
    for j = 1:2*L+1, dx = Xp(:,j) - xp; Pp = Pp + Wc(j)*(dx * dx'); end

    Y = 108*(3 + 1.1*Xp(1,:)) - i_curr*R0 - Xp(2,:) - Xp(3,:);
    yp = Y * Wm';
    Pyy0 = 0; Pxy = zeros(L,1);
    for j = 1:2*L+1
      dy = Y(j) - yp; dx = Xp(:,j) - xp;
      Pyy0 = Pyy0 + Wc(j)*dy^2; Pxy = Pxy + Wc(j)*dx*dy;
    end

    innov = vMeas(k+1) - yp;
    dk = (1 - b) / (1 - b^(k+1));
    R = max(1e-8, (1 - dk)*R + dk*(innov^2 - Pyy0));
    Pyy = Pyy0 + R;
    K = Pxy / Pyy;
    x = xp + K * innov;
    P = Pp - K * Pyy * K';

    x(1) = min(1.0, max(0.0, x(1)));
    socHat(k+1) = x(1);
    vHat(k+1) = 108*(3 + 1.1*x(1)) - i_curr*R0 - x(2) - x(3);
  end

  err = abs(socHat - socT) * 100;
  idx = find(err > 1.0, 1, 'last');
  if isempty(idx), res.tconv = 0; else res.tconv = idx; end
  res.err_mean = mean(err(60:end));
  res.err_max = max(err(60:end));
  res.rmse_mV = 1e3 * sqrt(mean((vHat(60:end) - vClean(60:end)).^2));
end
