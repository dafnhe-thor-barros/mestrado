function simulacao_4_6_fadiga_bga_pue()
  % ============================================================================
  % SIMULACAO 4.6: FADIGA TERMOMECANICA DE SOLDAS BGA SAC305 (ENGELMAIER COM
  % LDNP = 19.80 mm) E INDICADORES DE EFICACIA ENERGETICA PUE ORBITAL
  % ============================================================================
  % Dissertacao de Mestrado Profissional - ProPECAut / IFES
  % Autor: Dafnhe Thor Castro de Martins Barros
  % ============================================================================
  close all; clc;

  fprintf('==============================================================================\n');
  fprintf('INICIANDO SIMULACAO 4.6: FADIGA BGA SAC305 (Ldnp = 19.80 mm) E PUE ORBITAL\n');
  fprintf('==============================================================================\n\n');

  Lx = 28.0e-3;
  Ly = 28.0e-3;
  Ldnp = 0.5 * sqrt(Lx^2 + Ly^2); % Geometria 28x28 mm: 19.799 mm (~19.80 mm)
  hs = 0.18e-3;
  deltaAlpha = 14.0e-6;
  Fconf = 0.90;
  katt = 0.85;
  twoEpsF = 2.887;                % [A CONFERIR: fonte primaria]
  c0 = -0.470;
  alphaT = -7.0e-4;
  betaF = 0.020;

  names = {'SSO + TEMS', 'SSO sem TEMS', 'LEO + TEMS', 'LEO sem TEMS'};
  dTj = [7.0, 22.5, 18.2, 48.0];
  dTcp = [0.0, 0.0, 0.0, 0.0];
  Tmean = [60.0, 60.0, 50.0, 50.0];
  fcpd = [14.58, 14.58, 15.22, 15.22];

  fprintf('Ldnp calculado (canto do die 28x28 mm) = %.3f mm\n', Ldnp * 1e3);
  fprintf('2epsF = %.3f [A CONFERIR]\n\n', twoEpsF);

  for j = 1:4
    dTs = katt * dTj(j) + (1 - katt) * dTcp(j);
    gamma = Fconf * Ldnp * deltaAlpha * dTs / hs;
    c = c0 + alphaT * Tmean(j) + betaF * log(1 + fcpd(j));
    Nf = 0.5 * (gamma / twoEpsF)^(1 / c);
    years = Nf / (fcpd(j) * 365.25);
    margin = 100 * (years - 10.0) / 10.0;
    fprintf('%-15s | dTs=%5.2f C | gamma=%.4e | c=%.4f | Nf=%7.0f | vida=%5.2f anos | margem=%+6.1f%%\n', ...
            names{j}, dTs, gamma, c, Nf, years, margin);
  end
  fprintf('\n');

  % PUE nominal e pior caso quente, usando os valores da tese
  PIT = 40000;
  PcoolNom = 200; PcoolWCH = 400;
  PlossGaNNom = 1025; PlossBESSNom = 65;
  PlossGaNWCH = 4444; PlossBESSWCH = 120;
  PauxNom = 350 + 343; PauxWCH = 450 + 343;
  PtotalNom = PIT + PcoolNom + PlossGaNNom + PlossBESSNom + PauxNom;
  PtotalWCH = PIT + PcoolWCH + PlossGaNWCH + PlossBESSWCH + PauxWCH;
  tops = 80 * 250;
  fprintf('Nominal: Ptotal=%.0f W | PUE=%.4f | etaComp=%.4f TOPS/W\n', ...
          PtotalNom, PtotalNom / PIT, tops / PtotalNom);
  fprintf('WCH:     Ptotal=%.0f W | PUE=%.4f | etaComp=%.4f TOPS/W\n', ...
          PtotalWCH, PtotalWCH / PIT, tops / PtotalWCH);
  fprintf('==============================================================================\n');
end
