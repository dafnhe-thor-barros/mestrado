% ==============================================================================
% SIMULAÇÃO 4.1 (EXPANDIDA): PROPAGAÇÃO DE FLUXOS RADIATIVOS AMBIENTAIS
% E DINÂMICA TÉRMICA MULTI-ÓRBITA EM 100 ÓRBITAS LEO
% ==============================================================================
% Dissertação de Mestrado Profissional - ProPECAut / IFES
% Autor: Dafnhê Thor Castro de Martins Barros
% Revisão Técnica: Modelagem Explícita de Condições Nominais e Modos de Falha
% ==============================================================================
function executar_fluxos_orbitais_100_orbitas()
  close all; clc;

  sigma = 5.670374e-8; G0 = 1361.0; qEarth_IR = 237.0; albedo = 0.30;
  RE = 6371e3; RSun = 696340e3; dES = 1.496e11; mu = 398600.4418e9;
  
  alphaEOL = 0.14; epsEOL = 0.88; % OSR EOL
  etaFin = 0.945; Qint = 41983.0; Arad = 88.6; dT_cond = 25.0;
  phiRad = 90.0;
  Trad_lim = 55.0; T_space = 3.0;

  % Capacitancia termica total dos paineis radiadores + condutores (APG / Al / OSR)
  % 88.6 m2 * 4.5 kg/m2 * 900 J/(kg*K) = 358.830 J/K
  Crad = 358830.0;

  regimes = {
    'SSO 700 km', 700, 5926.0, 0.32, 58.3, 0.0;
    'LEO 500 km', 500, 5677.0, 0.45, 20.0, 2148.0
  };

  alvos = [1, 5, 10, 20, 50, 100];
  nOrb = 100;
  dt = 2.0;

  for r = 1:size(regimes, 1)
    nome_regime = regimes{r, 1};
    hkm = regimes{r, 2};
    Torb = regimes{r, 3};
    Kalb = regimes{r, 4};
    beta_deg = regimes{r, 5};
    tecl = regimes{r, 6};

    F = calcula_fator_visao(hkm, phiRad);
    qAlb_base = alphaEOL * G0 * albedo * F * Kalb;
    qIR_base = epsEOL * qEarth_IR * F;
    qAmb_nominal = qAlb_base + qIR_base;

    % Margens de apontamento analiticas
    [a_throt, a_trip] = calcula_margens_atitude(qAmb_nominal, alphaEOL, G0, Qint, Arad, etaFin, epsEOL, sigma, dT_cond);

    steps_per_orb = floor(Torb / dt);
    N = nOrb * steps_per_orb + 1;
    Trad = zeros(1, N);
    Trad(1) = 20.0 + 273.15; % Condicao inicial de implantacao orbital: 20 C

    fprintf('\n===============================================================================================================\n');
    fprintf('REGIME ORBITAL: %s (Altitude = %d km | Periodo = %.1f min | F_rad_E = %.5f)\n', ...
            nome_regime, hkm, Torb/60, F);
    fprintf('Carga Ambiental Media: %.2f W/m2 (Albedo: %.1f, IR: %.1f) | Margem Atitude: Throt=%.2f deg, Trip=%.2f deg\n', ...
            qAmb_nominal, qAlb_base, qIR_base, a_throt, a_trip);
    fprintf('===============================================================================================================\n');
    fprintf('%-8s | %-8s | %-24s | %-24s | %-12s | %-12s\n', ...
            'Orbita', 'Tempo', 'T_rad (Min/Med/Max) C', 'T_j (Min/Med/Max) C', 'q_amb (W/m2)', 'Status Tj');
    fprintf('---------------------------------------------------------------------------------------------------------------\n');

    for k = 1:N-1
      t_curr = (k-1) * dt;
      t_in_orb = mod(t_curr, Torb);

      % Perfil de fluxo ambiental incidente com fase orbital
      if tecl > 0
        % LEO com eclipse periodico
        if t_in_orb < tecl
          q_amb_now = qIR_base; % Em eclipse apenas IR terrestre
        else
          t_sun = Torb - tecl;
          fase_sun = (t_in_orb - tecl) / t_sun;
          q_amb_now = qIR_base + (qAlb_base / (t_sun/Torb)) * max(0, sin(pi * fase_sun));
        end
      else
        % SSO Dawn-Dusk com variacao suave de irradiacao terrestre
        q_amb_now = qAmb_nominal + 4.3 * cos(2 * pi * t_in_orb / Torb);
      end

      % Equacao diferencial do balanco termico do radiador no vacuo
      q_emit = etaFin * epsEOL * sigma * Arad * (Trad(k)^4 - T_space^4);
      dTrad = (Qint + q_amb_now * Arad - q_emit) / Crad;
      if k < N, Trad(k+1) = Trad(k) + dTrad * dt; end

      if mod(k, steps_per_orb) == 0
        orb_atual = k / steps_per_orb;
        if ismember(orb_atual, alvos)
          idx_s = (orb_atual - 1) * steps_per_orb + 1;
          idx_e = k + 1;
          t_slice = Trad(idx_s:idx_e) - 273.15;
          tj_slice = t_slice + dT_cond;
          t_horas = (k * dt) / 3600;

          if max(tj_slice) <= 80.0
            st = 'CONFORME';
          elseif max(tj_slice) <= 85.0
            st = 'THROTTLING';
          else
            st = 'TRIP CRIT';
          end

          fprintf('%-8d | %6.1f h | [%5.2f, %5.2f, %5.2f] C | [%5.2f, %5.2f, %5.2f] C | %10.2f   | %-12s\n', ...
                  orb_atual, t_horas, min(t_slice), mean(t_slice), max(t_slice), ...
                  min(tj_slice), mean(tj_slice), max(tj_slice), qAmb_nominal, st);
        end
      end
    end
  end
  fprintf('\n===============================================================================================================\n');
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

function [a1, a2] = calcula_margens_atitude(qAmb, alpha, G0, Q, A, eta, epsi, sigma, dTc)
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
end
