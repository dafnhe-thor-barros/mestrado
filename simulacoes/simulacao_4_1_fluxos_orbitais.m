function simulacao_4_1_fluxos_orbitais()
  % ============================================================================
  % SIMULACAO 4.1: FLUXOS RADIATIVOS ORBITAIS (SSO, LEO E GEO), FATOR DE VISAO
  % ANALITICO (EQ. 75), OCULTACAO CONICA E BALANCO TERMICO EOL (TABELA 10)
  % ============================================================================
  % Dissertacao de Mestrado Profissional - ProPECAut / IFES
  % Autor: Dafnhe Thor Castro de Martins Barros
  % ============================================================================
function simulacao_4_1_fluxos_orbitais()
  close all; clc;

  % =========================================================================
  % 1. Constantes e Propriedades Termo-Opticas EOL (Tabela 10)
  % =========================================================================
  sigma = 5.670374e-8; G0 = 1361.0; qEarth_IR = 237.0; albedo = 0.30;
  RE = 6371e3; RSun = 696340e3; dES = 1.496e11;

  alphaEOL = 0.14; epsEOL = 0.88; % OSR de fim de vida (EOL)
  etaFin = 0.945; Qint = 41983.0; Arad = 88.6; dT_cond = 25.0;
  phiRad = 90.0;              % Radiadores normais ao plano orbital (graus)
  Trad_lim = 55.0;            % Limite para Tj <= 80 C

  hkm = [700 500 35786];
  labels = {'SSO 700 km', 'LEO 500 km', 'GEO'};
  Kalb = [0.32 0.45 0.10];     % Coeficientes de incidencia de albedo

  fprintf('=== 1. Balanco Termico dos Radiadores em Regime Permanente ===\n');
  q_emit_max = etaFin * epsEOL * sigma * (Trad_lim + 273.15)^4;

  for j = 1:numel(hkm)
    F = calcula_fator_visao(hkm(j), phiRad);
    qAlb = alphaEOL * G0 * albedo * F * Kalb(j);
    qIR = epsEOL * qEarth_IR * F;
    qAmb = qAlb + qIR;

    Trad = ((Qint/Arad + qAmb) / (etaFin * epsEOL * sigma))^0.25 - 273.15;
    Tj = Trad + dT_cond;
    Areq = Qint / (q_emit_max - qAmb);

    fprintf('%s:\n', labels{j});
    fprintf('  F_rad_E: %.5f | Carga Ambiental: %.2f W/m2 (Albedo: %.1f, IR: %.1f)\n', ...
            F, qAmb, qAlb, qIR);
    fprintf('  Trad: %.2f C | Tj: %.2f C | Area Requerida (Trad <= 55 C): %.2f m2\n', ...
            Trad, Tj, Areq);
    avalia_margem_atitude(qAmb, alphaEOL, G0, Qint, Arad, etaFin, epsEOL, sigma, dT_cond);
  end

  % =========================================================================
  % 2. Geometria de Eclipse Conico em LEO 500 km
  % =========================================================================
  h_leo = 500e3; mu = 398600.4418e9;
  Torb_leo = 2 * pi * sqrt((RE + h_leo)^3 / mu);
  t_vec = linspace(0, Torb_leo, 3000);
  aEarth = asin(RE / (RE + h_leo));
  aSun = asin(RSun / dES);

  S_frac = ones(size(t_vec));
  for k = 1:numel(t_vec)
    theta = 2 * pi * t_vec(k) / Torb_leo;
    gamma = acos(max(-1, min(1, -cos(theta))));
    if gamma <= aEarth - aSun
      S_frac(k) = 0.0; % Umbra total
    elseif gamma < aEarth + aSun
      % Penumbra
      r1 = aSun; r2 = aEarth; d = gamma;
      x1 = (d^2 + r1^2 - r2^2) / (2*d*r1);
      x2 = (d^2 + r2^2 - r1^2) / (2*d*r2);
      rad = max(0, (-d + r1 + r2)*(d + r1 - r2)*(d - r1 + r2)*(d + r1 + r2));
      overlap = r1^2 * acos(max(-1,min(1,x1))) + r2^2 * acos(max(-1,min(1,x2))) - 0.5*sqrt(rad);
      S_frac(k) = max(0, min(1, 1 - overlap / (pi * r1^2)));
    else
      S_frac(k) = 1.0; % Sol pleno
    end
  end

  dt_step = t_vec(2) - t_vec(1);
  t_umbra = sum(S_frac == 0) * dt_step / 60;
  t_ecl_tot = sum(S_frac < 1) * dt_step / 60;
  fprintf('\n=== 2. Geometria de Eclipse em LEO 500 km ===\n');
  fprintf('Periodo Orbital: %.2f min | Umbra: %.2f min | Eclipse Total: %.2f min\n', ...
          Torb_leo/60, t_umbra, t_ecl_tot);
end

function F = calcula_fator_visao(hkm, phiDeg)
  RE = 6371.0; H = (RE + hkm) / RE; phi = phiDeg * pi / 180;
  thetaM = asin(1 / H); r = sqrt(H^2 - 1);
  if phi <= pi/2 - thetaM
    F = cos(phi) / H^2;
  elseif phi <= pi/2 + thetaM
    arg1 = max(-1, min(1, r / (H * sin(phi))));
    arg2 = max(-1, min(1, r * cot(phi)));
    rad = max(0, 1 - H^2 * cos(phi)^2);
    F = 0.5 - asin(arg1)/pi + cos(phi)/(pi*H^2)*acos(arg2) - r*sqrt(rad)/(pi*H^2);
  else
    F = 0;
  end
  F = max(0, min(1, F));
end

function avalia_margem_atitude(qAmb, alpha, G0, Q, A, eta, epsi, sigma, dTc)
  theta = 0:0.01:90;
  Tj = zeros(size(theta));
  for k = 1:numel(theta)
    qsol = alpha * G0 * sind(theta(k));
    Tj(k) = ((Q/A + qAmb + qsol) / (eta * epsi * sigma))^0.25 - 273.15 + dTc;
  end
  i_throt = find(Tj >= 80.0, 1);
  i_trip = find(Tj >= 85.0, 1);
  if isempty(i_throt), a1 = NaN; else, a1 = theta(i_throt); end
  if isempty(i_trip), a2 = NaN; else, a2 = theta(i_trip); end
  fprintf('  Margem de Atitude Solar: Throttling (80 C) = %.2f deg | Trip (85 C) = %.2f deg\n', a1, a2);
end
