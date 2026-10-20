function simulacao_4_1_fluxos_orbitais()
  % ============================================================================
  % SIMULACAO 4.1: FLUXOS RADIATIVOS ORBITAIS (SSO, LEO E GEO), FATOR DE VISAO
  % ANALITICO (EQ. 75), OCULTACAO CONICA E BALANCO TERMICO EOL (TABELA 10)
  % ============================================================================
  % Dissertacao de Mestrado Profissional - ProPECAut / IFES
  % Autor: Dafnhe Thor Castro de Martins Barros
  % ============================================================================
  close all; clc;

  fprintf('==============================================================================\n');
  fprintf('INICIANDO SIMULACAO 4.1: FLUXOS ORBITAIS E DIMENSIONAMENTO TERMICO (TABELA 10)\n');
  fprintf('==============================================================================\n\n');

  sigma = 5.670374e-8;       % [W/(m^2 K^4)] Constante de Stefan-Boltzmann
  G = 1361.0;                % [W/m^2] Constante solar AM0 (ASTM E490)
  qEarth = 237.0;            % [W/m^2] Emissao infravermelha media terrestre
  albedo = 0.30;             % Albedo planetario medio
  RE = 6371.0e3;             % [m] Raio volumetrico medio da Terra
  RSun = 696340.0e3;         % [m] Raio do Sol
  dES = 1.496e11;            % [m] Distancia Terra-Sol (1 UA)

  % Propriedades termo-opticas EOL (10 anos) do espelho OSR segundo a Tabela 10
  alphaEOL = 0.14;           % Absortividade solar EOL (Tabela 10)
  epsEOL   = 0.88;           % Emissividade infravermelha EOL (Tabela 10)
  etaFin   = 0.945;          % Eficiencia da aleta APG/K13D2U
  Qint     = 41983.0;        % [W] Calor interno total dissipado no no
  Arad     = 88.6;           % [m^2] Area efetiva instalada de radiadores
  dTcond   = 25.0;           % [K] Salto termico condutivo juncao-radiador (LHP)
  phiRad   = 90.0;           % [deg] Placa ortogonal ao nadir [HIPOTESE geometrica]

  TjThrot = 80.0;            % [degC] Limiar de inicio de thermal throttling
  TjTrip  = 85.0;            % [degC] Limiar de corte termico (Thermal Trip)
  TradLim = 55.0;            % [degC] Teto de projeto do radiador

  hkm = [700, 500, 35786];
  labels = {'SSO 700 km', 'LEO 500 km', 'GEO 35786 km'};
  Kalb = [0.32, 0.45, 0.10]; % [HIPOTESE: fator medio orbital de albedo]

  qemit = etaFin * epsEOL * sigma * (TradLim + 273.15)^4;
  fprintf('Capacidade de emissao bruta a 55 C (alpha=0.14, eps=0.88): %.2f W/m^2\n\n', qemit);

  for j = 1:numel(hkm)
    F = view_factor(hkm(j), phiRad);
    qAlb = alphaEOL * G * albedo * F * Kalb(j);
    qIR  = epsEOL * qEarth * F;
    qAmb = qAlb + qIR;
    Trad = ((Qint / Arad + qAmb) / (etaFin * epsEOL * sigma))^0.25 - 273.15;
    Tj   = Trad + dTcond;
    Areq = Qint / (qemit - qAmb);
    fprintf('%s (Eq. 75 phi=90 deg): F=%.5f | qalb=%.2f | qir=%.2f | qamb=%.2f W/m^2\n', ...
            labels{j}, F, qAlb, qIR, qAmb);
    fprintf('  Trad=%.2f C | Tj=%.2f C | Areq(Trad<=55 C)=%.2f m^2 (Folga: %+.2f m^2)\n', ...
            Trad, Tj, Areq, Arad - Areq);
    print_angle(qAmb, alphaEOL, G, Qint, Arad, etaFin, epsEOL, sigma, dTcond, TjThrot, TjTrip);
  end

  % Comparacao adicional para LEO 500 km sob fator de visao tabelado F = 0.4473
  F_leo_tab = 0.4473;
  qAlb_lt = alphaEOL * G * albedo * F_leo_tab * Kalb(2);
  qIR_lt  = epsEOL * qEarth * F_leo_tab;
  qAmb_lt = qAlb_lt + qIR_lt;
  Trad_lt = ((Qint / Arad + qAmb_lt) / (etaFin * epsEOL * sigma))^0.25 - 273.15;
  Tj_lt   = Trad_lt + dTcond;
  Areq_lt = Qint / (qemit - qAmb_lt);
  fprintf('\nLEO 500 km (Fator tabelado F=0.4473): qalb=%.2f | qir=%.2f | qamb=%.2f W/m^2\n', ...
          qAlb_lt, qIR_lt, qAmb_lt);
  fprintf('  Trad=%.2f C | Tj=%.2f C | Areq(Trad<=55 C)=%.2f m^2 (Deficit: %.2f m^2)\n\n', ...
          Trad_lt, Tj_lt, Areq_lt, Areq_lt - Arad);

  % Eclipse conico em LEO circular 500 km
  h = 500e3;
  mu = 398600.4418e9;
  Torb = 2 * pi * sqrt((RE + h)^3 / mu);
  t = linspace(0, Torb, 2000);
  aEarth = asin(RE / (RE + h));
  aSun = asin(RSun / dES);
  S = ones(size(t));
  for k = 1:numel(t)
    theta = 2 * pi * t(k) / Torb;
    gamma = acos(max(-1, min(1, -cos(theta))));
    if gamma <= aEarth - aSun
      S(k) = 0;
    elseif gamma < aEarth + aSun
      S(k) = sun_fraction(gamma, aEarth, aSun);
    else
      S(k) = 1;
    end
  end
  eclipseMinutes = sum(S < 1) * Torb / (numel(t) - 1) / 60;
  fprintf('Orbita LEO 500 km = %.2f min; eclipse geometrico aproximado = %.2f min\n', ...
          Torb / 60, eclipseMinutes);
  fprintf('==============================================================================\n');
end

function F = view_factor(hkm, phiDeg)
  RE = 6371.0;
  H = (RE + hkm) / RE;
  phi = phiDeg * pi / 180;
  thetaM = asin(1 / H);
  r = sqrt(H^2 - 1);
  if phi <= pi / 2 - thetaM
    F = cos(phi) / H^2;
  elseif phi <= pi / 2 + thetaM
    arg1 = max(-1, min(1, r / (H * sin(phi))));
    arg2 = max(-1, min(1, r * cot(phi)));
    rad = max(0, 1 - H^2 * cos(phi)^2);
    F = 0.5 - asin(arg1) / pi + cos(phi) / (pi * H^2) * acos(arg2) ...
        - r * sqrt(rad) / (pi * H^2);
  else
    F = 0;
  end
  F = max(0, min(1, F));
end

function S = sun_fraction(gamma, aEarth, aSun)
  r1 = aSun;
  r2 = aEarth;
  d = gamma;
  if d == 0
    S = 0;
    return;
  end
  x1 = max(-1, min(1, (d^2 + r1^2 - r2^2) / (2 * d * r1)));
  x2 = max(-1, min(1, (d^2 + r2^2 - r1^2) / (2 * d * r2)));
  rad = max(0, (-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2));
  overlap = r1^2 * acos(x1) + r2^2 * acos(x2) - 0.5 * sqrt(rad);
  S = max(0, min(1, 1 - overlap / (pi * r1^2)));
end

function print_angle(qAmb, alpha, G, Q, A, eta, epsi, sigma, dTc, Tthrot, Ttrip)
  theta = 0:0.01:90;
  Tj = zeros(size(theta));
  for k = 1:numel(theta)
    qsol = alpha * G * sind(theta(k));
    Tj(k) = ((Q / A + qAmb + qsol) / (eta * epsi * sigma))^0.25 - 273.15 + dTc;
  end
  i1 = find(Tj >= Tthrot, 1);
  i2 = find(Tj >= Ttrip, 1);
  if isempty(i1), x1 = NaN; else x1 = theta(i1); end
  if isempty(i2), x2 = NaN; else x2 = theta(i2); end
  fprintf('  Angulos estimados: throttling=%.2f deg; trip=%.2f deg\n', x1, x2);
end
