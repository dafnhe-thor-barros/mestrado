function executar_bateria_100_orbitas()
  close all; clc;

  Cn = 22.52;
  Vuv = 324.0; Vov = 442.8;
  Rp1 = 18e-3; Cp1 = 3200;
  Rp2 = 24e-3; Cp2 = 22000;
  R0_ref = 85e-3;
  Ea = 0.32 * 1.602176634e-19;
  kB = 1.380649e-23;
  Tref_K = 293.15;

  Cth = 54360;
  Rth_cond = 0.14;
  Tch = 15.0;
  Pheater = 300.0; Ton = 10.0; Toff = 15.0;

  dt = 2.0; Torb = 5926; tecl = 1128; nOrb = 100;
  steps_per_orb = floor(Torb / dt);
  N = nOrb * steps_per_orb + 1;

  soc = zeros(1, N); vt = zeros(1, N); Tb = zeros(1, N); Pj = zeros(1, N);
  soc(1) = 0.90; Tb(1) = 20.0;
  vp1 = 0; vp2 = 0;
  a1 = exp(-dt / (Rp1 * Cp1)); a2 = exp(-dt / (Rp2 * Cp2));

  alvos = [1, 5, 10, 20, 50, 100];

  fprintf('%-8s | %-8s | %-16s | %-16s | %-24s | %-10s | %-10s\n', ...
          'Orbita', 'Tempo', 'SoC Min/Max (%)', 'Tensao Min/Max(V)', 'T_bat (Min/Med/Max) C', 'Pj_med (W)', 'Teq_est C');
  fprintf('---------------------------------------------------------------------------------------------------------------\n');

  for k = 1:N-1
    t_now = (k-1) * dt;
    in_eclipse = mod(t_now, Torb) <= tecl;
    
    if in_eclipse
      I = (18000 + 560) / 380;
      etaC = 1.0;
    else
      I = -25.0;
      etaC = 0.995;
    end

    R0 = R0_ref * exp((Ea/kB) * (1/(Tb(k) + 273.15) - 1/Tref_K));
    vp1 = a1 * vp1 + Rp1 * (1 - a1) * I;
    vp2 = a2 * vp2 + Rp2 * (1 - a2) * I;

    soc(k+1) = min(1.0, max(0.0, soc(k) - I * etaC * dt / (Cn * 3600)));
    Voc = 108 * (3.0 + 1.1 * soc(k+1));
    Vcand = Voc - I * R0 - vp1 - vp2;

    if I < 0 && Vcand > Vov
      I = max(-25.0, min(0.0, -(Vov - Voc + vp1 + vp2) / R0));
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

    if mod(k, steps_per_orb) == 0
      orb_atual = k / steps_per_orb;
      if ismember(orb_atual, alvos)
        idx_start = (orb_atual - 1) * steps_per_orb + 1;
        idx_end = k + 1;
        
        soc_min = min(soc(idx_start:idx_end)) * 100;
        soc_max = max(soc(idx_start:idx_end)) * 100;
        vt_min = min(vt(idx_start:idx_end));
        vt_max = max(vt(idx_start:idx_end));
        tb_min = min(Tb(idx_start:idx_end));
        tb_med = mean(Tb(idx_start:idx_end));
        tb_max = max(Tb(idx_start:idx_end));
        pj_med = mean(Pj(idx_start:k));
        teq_est = Tch + pj_med * Rth_cond;
        t_horas = (k * dt) / 3600;

        fprintf('%-8d | %6.1f h | [%5.2f, %6.2f]%% | [%5.1f, %5.1f] V | [%5.2f, %5.2f, %5.2f] C | %8.1f W | %7.2f C\n', ...
                orb_atual, t_horas, soc_min, soc_max, vt_min, vt_max, tb_min, tb_med, tb_max, pj_med, teq_est);
      end
    end
  end
end
