function simulacao_4_6_fadiga_bga_pue()
  % ============================================================================
  % SIMULACAO 4.6: FADIGA TERMOMECANICA DE SOLDAS BGA SAC305 (ENGELMAIER COM
  % LDNP = 19.80 mm) E INDICADORES DE EFICACIA ENERGETICA PUE ORBITAL
  % ============================================================================
  % Dissertacao de Mestrado Profissional - ProPECAut / IFES
  % Autor: Dafnhe Thor Castro de Martins Barros
  % ============================================================================
function simulacao_4_6_fadiga_bga_protecao()
  close all; clc;

  % =========================================================================
  % 1. Fadiga Termomecanica de Solda BGA SAC305 (Modelo de Engelmaier)
  % =========================================================================
  % Geometria real: encapsulamento quadrado de 28 mm x 28 mm
  Lx = 28e-3; Ly = 28e-3;
  Ldnp = 0.5 * sqrt(Lx^2 + Ly^2); % Meia-diagonal real = 19.80 mm (corrigida de 14 mm)
  hs = 0.18e-3;                   % Altura da esfera de solda [m]
  deltaAlpha = 14e-6;             % Descasamento CTE (pastilha vs PCB) [1/K]
  Fconf = 0.90; katt = 0.85;
  twoEpsF = 2.887;                % Ductilidade em fadiga SAC305
  c0 = -0.470; alphaT = -7e-4; betaF = 0.020;

  casos = {'SSO + TEMS', 7.0, 60.0, 14.58;
           'SSO sem TEMS', 22.5, 60.0, 14.58;
           'LEO + TEMS', 18.2, 50.0, 15.22;
           'LEO sem TEMS', 48.0, 50.0, 15.22};

  fprintf('=== 1. Analise de Fadiga BGA SAC305 (Engelmaier) ===\n');
  fprintf('Distancia ao Ponto Neutro (L_DNP): %.2f mm\n', Ldnp*1e3);

  for j = 1:size(casos, 1)
    nome = casos{j,1}; dTj = casos{j,2}; Tmean = casos{j,3}; fcpd = casos{j,4};

    dTs = katt * dTj;
    gamma = Fconf * Ldnp * deltaAlpha * dTs / hs;
    c = c0 + alphaT * Tmean + betaF * log(1 + fcpd);
    Nf = 0.5 * (gamma / twoEpsF)^(1 / c);
    vida_anos = Nf / (fcpd * 365.25);
    margem = 100 * (vida_anos - 10.0) / 10.0;

    fprintf('%-15s: dTs = %4.1f C | gamma = %.3e | Nf = %8.0f | Vida = %5.1f anos (Margem: %+5.1f %%)\n', ...
            nome, dTs, gamma, Nf, vida_anos, margem);
  end

  % =========================================================================
  % 2. Protecao LCL em Barramento de 400 V / 50 kW (Tailoring Funcional)
  % =========================================================================
  fprintf('\n=== 2. Dimensionamento do LCL de 400 V / 50 kW ===\n');
  V_bus = 400.0; P_total = 50e3;
  I_nom = P_total / V_bus;        % 125 A
  I_lim = 1.25 * I_nom;           % Ponto medio de limitacao ativa (156.25 A)
  L_harness = 2.0e-6;             % Indutancia de linha [H]

  % Tempo maximo de trip para suportar curto sem exceder limite termico SiC
  t_trip_us = 15.0;               % 15 us (Janela rapida SSPC, contra ms da ECSS 50 V)
  E_dissipada_chave = V_bus * I_lim * (t_trip_us * 1e-6);

  % Desmagnetizacao indutiva absorvida por TVS (ex: 2x 15KPA200A ou SMCJ400A)
  V_clamp_tvs = 644.0;
  E_indutiva = 0.5 * L_harness * I_lim^2;
  t_desmag = (L_harness * I_lim) / (V_clamp_tvs - V_bus);

  fprintf('Corrente Nominal: %.1f A | Limite de Corte Ativo: %.1f A\n', I_nom, I_lim);
  fprintf('Janela de Desarme (Trip-off): %.1f us\n', t_trip_us);
  fprintf('Energia Dissipada no Semicondutor no Trip: %.2f J\n', E_dissipada_chave);
  fprintf('Energia Indutiva no TVS: %.3f J | Tempo de Desmagnetizacao: %.2f us\n', ...
          E_indutiva, t_desmag * 1e6);

  % =========================================================================
  % 3. Metricas PUE da Carga Util de Processamento
  % =========================================================================
  PIT = 40000.0; tops_total = 80 * 250;
  Ptotal_nom = PIT + 200 + 1025 + 65 + (350 + 343);
  Ptotal_wch = PIT + 400 + 4444 + 120 + (450 + 343);

  fprintf('\n=== 3. Metricas PUE do Sistema ===\n');
  fprintf('Nominal: Ptotal = %.0f W | PUE = %.4f | Eficiencia: %.4f TOPS/W\n', ...
          Ptotal_nom, Ptotal_nom/PIT, tops_total/Ptotal_nom);
  fprintf('Pior Caso: Ptotal = %.0f W | PUE = %.4f | Eficiencia: %.4f TOPS/W\n', ...
          Ptotal_wch, Ptotal_wch/PIT, tops_total/Ptotal_wch);
end
