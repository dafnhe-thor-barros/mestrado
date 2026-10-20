function simulacao_3_4_conversores_gan_middlebrook()
  % ============================================================================
  % SIMULACAO 3.4: CONVERSORES GaN (GS66516T), STEINMETZ MODIFICADO (MSE,
  % NUCLEO FERROXCUBE U93/152/30 3C95), INTERFACE AlN+PGS E MIDDLEBROOK
  % ============================================================================
  % Dissertacao de Mestrado Profissional - ProPECAut / IFES
  % Autor: Dafnhe Thor Castro de Martins Barros
  % ============================================================================
  close all; clc;

  % Parametros de Operacao do Conversor
  Pmod = 20000;      % Potencia modular (20 kW)
  Vin = 400;         % Tensao de entrada (400 V)
  Vout = 48;         % Tensao de saida (48 V)
  fsw = 150e3;       % Frequencia de chaveamento (150 kHz)
  Nsw = 4;           % Quatro chaves no caso de ponte unica
  D = 0.45;          % Razao de ciclo (Duty cycle)
  Iin = 51.2;        % Corrente de entrada (A)
  Iout = Pmod / Vout;% Corrente de saida (A)
  Tc = 60;           % Temperatura do case/chassi (C)

  % GS66516T: Dados nominais extraidos do datasheet (Infineon/GaN Systems)
  Ron25 = 25e-3;
  alpha = (65e-3 / Ron25 - 1) / (150 - 25); % Coeficiente termico (0.0128 1/K)
  Coss = 207e-12;    % Capacitancia equivalente de saida
  Qg = 14.2e-9;      % Carga de gate total
  Vgate = 6;         % Tensao de acionamento do gate (V)
  tr = 12e-9;        % Tempo de subida (Rise time)
  tf = 15e-9;        % Tempo de descida (Fall time)
  Isw = 25;          % Corrente de comutacao considerada
  Rth_jc  = 0.27;    % Resistencia juncao-case do GS66516T (K/W)
  Rth_tim = 0.85;    % Isolador ceramico AlN + folha grafite pirolitico PGS a vacuo (K/W)
  Rth_total_hw = Rth_jc + Rth_tim; % 1.12 K/W

  % Parametros do Secundario e Elementos Magneticos (Par U93/152/30, Ferrite 3C95)
  RonSec25 = 0.80e-3;
  alphaSec = 0.004;
  kfe = 0.70;        % Coeficiente efetivo MSE (excitacao trapezoidal)
  alphaSt = 1.46;
  betaSt = 2.52;
  Ve = 4.5e-4;       % Volume efetivo U93/152/30 = 450.000 mm^3 (Ae=840 mm^2, le=536 mm)
  dB = 0.15;
  RacPri = 15e-3;
  RacSec = 0.25e-3;

  % Calculo das Correntes e Perdas por Conducao
  Irms = Iin * sqrt(D);
  Tj_est = 80; % Temperatura de juncao estimada inicial para perdas nominais
  RonT = Ron25 * (1 + alpha * (Tj_est - 25));
  PcondPri = Nsw * (Irms^2) * RonT;

  RonSecT = RonSec25 * (1 + alphaSec * (Tj_est - 25));
  PcondSec = (Iout^2) * RonSecT;

  % Perdas por Chaveamento (Por chave e total)
  PoverlapOne = Vin * Isw * fsw * (tr + tf) / 2;
  PossOne = 0.5 * Coss * (Vin^2) * fsw;
  PgateOne = Qg * Vgate * fsw;
  Psw = Nsw * (PoverlapOne + PossOne + PgateOne);

  % Perdas Magneticas (Nucleo + Cobre corrigido por Irms no primario)
  Pcore = kfe * (fsw^alphaSt) * (dB^betaSt) * Ve;
  Pcu = (Irms^2) * RacPri + (Iout^2) * RacSec;
  Pmag = Pcore + Pcu;

  % Perda Total e Eficiencia
  Ploss = PcondPri + PcondSec + Psw + Pmag;
  eta = 100 * Pmod / (Pmod + Ploss);

  % Derivada termica e limites de estabilidade intrinseca
  dPdT_one = (Irms^2) * Ron25 * alpha;
  dPdT = Nsw * dPdT_one;
  Rcrit_one = 1 / dPdT_one;
  Rcrit = 1 / dPdT;
  Icrit4 = sqrt(1 / (4 * Ron25 * alpha));

  fprintf('--- ANALISE DE PERDAS E EFICIENCIA NOMINAL ---\n');
  fprintf('Ron(%.0f C) = %.2f mOhm; perda total = %.2f W; eta = %.3f %%\n', Tj_est, 1e3*RonT, Ploss, eta);
  fprintf('Conducao prim. = %.2f W; comutacao = %.2f W; secundaria = %.2f W\n', PcondPri, Psw, PcondSec);
  fprintf('Nucleo (U93/152/30 3C95, MSE) = %.2f W; cobre = %.2f W\n', Pcore, Pcu);
  fprintf('dPloss/dT (1 chave) = %.4f W/K (Rcrit_1 = %.3f K/W)\n', dPdT_one, Rcrit_one);
  fprintf('dPloss/dT (4 chaves) = %.4f W/K (Rcrit_4 = %.3f K/W); Icrit(Rth=4) = %.2f A\n\n', dPdT, Rcrit, Icrit4);

  % Equilibrio Termico Iterativo da Chave Primaria (incluindo Rth_total_hw = 1.12 K/W)
  fprintf('--- ANALISE DE EQUILIBRIO TERMICO DA CHAVE (AlN + PGS: Rth_total = %.2f K/W) ---\n', Rth_total_hw);
  PswOne = PoverlapOne + PossOne + PgateOne;
  for Rth = [0.27, 1.00, Rth_total_hw, 2.00, 4.00]
    Tje = Tc;
    stable = true;
    for n = 1:500
      Pj = (Irms^2) * Ron25 * (1 + alpha * (Tje - 25)) + PswOne;
      Tnew = Tc + Rth * Pj;
      if Tnew > 400
        stable = false;
        break;
      end
      if abs(Tnew - Tje) < 1e-5
        Tje = Tnew;
        break;
      end
      Tje = Tnew;
    end
    if stable
      fprintf('Rth = %.2f K/W: Tj estimada = %.2f C\n', Rth, Tje);
    else
      fprintf('Rth = %.2f K/W: SEM EQUILIBRIO (fuga termica / limite do modelo)\n', Rth);
    end
  end
  fprintf('\n');

  % Criterio de Estabilidade de Middlebrook no Barramento de 48 V
  fprintf('--- ANALISE DE ESTABILIDADE DE MIDDLEBROOK ---\n');
  rin = -(Vout^2) / Pmod;
  Zlim = 0.10 * abs(rin);
  freq = logspace(1, 6, 5000);
  s = 1i * 2 * pi * freq;
  Lline = 200e-9;
  Rline = 0.8e-3;

  cases = {15e-3, 1.2e-3; 3e-3, 8e-3; 2.5e-3, 9.5e-3; 2.2e-3, 3e-3; 0.44e-3, 15e-3};
  names = {'Saudavel', 'WCA texto tese', 'Script antigo', 'Figura operacional', 'Figura falha'};

  fprintf('rin = %.5f ohm; limite (Middlebrook 10%%) = %.5f ohm\n', rin, Zlim);
  for j = 1:size(cases, 1)
    C = cases{j, 1};
    ESR = cases{j, 2};
    Z1 = s * Lline + Rline;
    Z2 = 1 ./ (s * C) + ESR;
    Z = abs((Z1 .* Z2) ./ (Z1 + Z2));
    [Zpeak, k] = max(Z);
    if Zpeak <= Zlim
      status = 'ESTAVEL';
    else
      status = 'VIOLA';
    end
    fprintf('%-18s: C = %6.3g F  ESR = %6.3g ohm | pico = %.5f ohm em %7.1f Hz | %s\n', ...
            names{j}, C, ESR, Zpeak, freq(k), status);
  end
end
