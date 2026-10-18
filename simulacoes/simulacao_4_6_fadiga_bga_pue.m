% ==============================================================================
% SIMULAÇÃO NUMÉRICA DO MODELO 4.6: FADIGA DE SOLDAS BGA E EFICÁCIA PUE
% Dissertação de Mestrado ProPECAut / IFES Campus Serra
% ==============================================================================
%
% Este script implementa:
% 1. O modelo analítico de Coffin-Manson generalizado por Engelmaier (IPC-SM-785 / IPC-9701A)
%    calibrado com coeficientes da literatura para a liga de solda SAC305.
%    NOTA EPISTEMICA: Analise parametrica de sensibilidade in silico, sem ensaios fisicos
%    ou de bancada TVAC no componente real.
% 2. O modelo térmico de desacoplamento dinâmico entre a junção de silício (die)
%    e a esfera crítica BGA (DNP).
% 3. A avaliação da Matriz Fatorial 2x2 isolando os efeitos da órbita (SSO vs. LEO)
%    e do sistema de controle ativo TEMS (com TEMS vs. sem TEMS).
% 4. A avaliacao parametrica in silico da Eficacia Energetica Orbital (PUE_orbital, adaptacao
%    conceitual da ISO/IEC 30134-2 para plataformas espaciais) e rendimento eta_comp.
% ==============================================================================

clear all;
close all;
clc;

fprintf('==============================================================================\n');
fprintf('  SIMULAÇÃO NUMÉRICA 4.6: FADIGA TERMOMECÂNICA BGA SAC305 E PUE ORBITAL\n');
fprintf('==============================================================================\n\n');

% ------------------------------------------------------------------------------
% 1. PARÂMETROS MECÂNICOS E ESTRUTURAIS DO ENCAPSULAMENTO BGA (COWOS)
% ------------------------------------------------------------------------------
L_die = 28.0;              % [mm] Dimensão linear do substrato do chip
L_DNP = 0.5 * sqrt(2) * L_die; % [mm] Distância ao centro neutro da esfera no canto (DNP)
Delta_alpha = 14.0e-6;     % [1/K] Diferencial de expansão térmica (alpha_sub - alpha_PCB)
h_solder = 0.18;           % [mm] Altura nominal pós-refluxo da esfera BGA
F_conformidade = 0.90;     % [-] Fator de conformidade estrutural e flexão da placa
k_att = 0.85;              % [-] Fator de atenuação térmica entre junção e solda

% Coeficientes de fadiga oligocíclica para SAC305 (Syed 2004, Gao et al. 2022)
eps_f_2 = 2.8866;          % [-] Ductilidade cíclica efetiva de cisalhamento BGA (2*eps_f_eff)
c_0 = -0.470;              % [-] Termo base do expoente de Engelmaier
alpha_T = -7.0e-4;         % [1/degC] Sensibilidade à temperatura média
beta_f = 2.0e-2;           % [-] Sensibilidade à frequência cíclica diária

% Requisito da missão espacial
anos_missao = 10.0;
T_orb_sso_min = 98.77;     % [min] Período em órbita SSO Dawn-Dusk (700 km)
T_orb_leo_min = 94.62;     % [min] Período em órbita LEO convencional (500 km)

f_cpd_sso = (24.0 * 60.0) / T_orb_sso_min; % 14.58 ciclos/dia
f_cpd_leo = (24.0 * 60.0) / T_orb_leo_min; % 15.22 ciclos/dia

N_ciclos_10anos_sso = f_cpd_sso * 365.25 * anos_missao; % ~ 53.245 ciclos
N_ciclos_10anos_leo = f_cpd_leo * 365.25 * anos_missao; % ~ 55.580 ciclos

fprintf('1. Requisitos Orbitais e Ciclagem em 10 Anos:\n');
fprintf('   - Órbita SSO Dawn-Dusk: %6.2f ciclos/dia -> %7.0f ciclos em 10 anos\n', f_cpd_sso, N_ciclos_10anos_sso);
fprintf('   - Órbita LEO Convencional: %6.2f ciclos/dia -> %7.0f ciclos em 10 anos\n\n', f_cpd_leo, N_ciclos_10anos_leo);

% ------------------------------------------------------------------------------
% 2. MATRIZ FATORIAL 2x2: AVALIAÇÃO DOS CENÁRIOS
% ------------------------------------------------------------------------------
% Cenário 1: SSO com TEMS (Ponto de Projeto)
% Cenário 2: SSO sem TEMS (Variação de carga não amortecida)
% Cenário 3: LEO com TEMS (Mitigação do eclipse pelo TEMS)
% Cenário 4: LEO sem TEMS (Comutação reativa sob eclipse severo)

Delta_Tj_vec     = [ 7.0,  22.5,  18.2,  48.0]; % [degC]
Delta_Tcp_vec    = [ 0.0,   0.0,   0.0,   0.0]; % [degC] Variação da placa fria
T_mean_solder    = [60.0,  60.0,  50.0,  50.0]; % [degC]
f_cpd_vec        = [f_cpd_sso, f_cpd_sso, f_cpd_leo, f_cpd_leo];
cenarios_nomes   = {'SSO + TEMS (Projeto)', 'SSO sem TEMS', 'LEO + TEMS', 'LEO sem TEMS (Falha)'};

fprintf('2. Avaliação dos Cenários Termomecânicos:\n');
fprintf('------------------------------------------------------------------------------------------------------\n');
fprintf('%-25s | Delta Tj | Delta Tsold | Gamma_eff   | Exp c    | Nf (ciclos)  | Vida (anos) | Margem\n', 'Cenário');
fprintf('------------------------------------------------------------------------------------------------------\n');

Nf_results = zeros(1, 4);
vida_anos_results = zeros(1, 4);

for i = 1:4
    % Atenuação térmica até a solda
    Delta_T_solder = k_att * Delta_Tj_vec(i) + (1.0 - k_att) * Delta_Tcp_vec(i);
    
    % Deformação cisalhante inelástica efetiva (Eq. 4.79)
    gamma_eff = F_conformidade * (L_DNP * Delta_alpha * Delta_T_solder) / h_solder;
    
    % Expoente c de Engelmaier calibrado para SAC305 (Eq. 4.81)
    c_exp = c_0 + alpha_T * T_mean_solder(i) + beta_f * log(1.0 + f_cpd_vec(i));
    
    % Vida útil em fadiga oligocíclica (Eq. 4.78)
    Nf = 0.5 * (gamma_eff / eps_f_2)^(1.0 / c_exp);
    Nf_results(i) = Nf;
    
    % Tempo até falha física
    vida_anos = Nf / (f_cpd_vec(i) * 365.25);
    vida_anos_results(i) = vida_anos;
    
    margem_10anos = ((vida_anos - 10.0) / 10.0) * 100.0;
    
    fprintf('%-25s | %5.1f °C |   %5.2f °C  | %9.3e | %8.4f | %12.2e | %7.2f anos | %+6.1f%%\n', ...
            cenarios_nomes{i}, Delta_Tj_vec(i), Delta_T_solder, gamma_eff, c_exp, Nf, vida_anos, margem_10anos);
end
fprintf('------------------------------------------------------------------------------------------------------\n\n');

% ------------------------------------------------------------------------------
% 3. BALANCO ELETRICO, BALANCO TERMICO E ADAPTACAO CONCEITUAL DO PUE ORBITAL
% ------------------------------------------------------------------------------
fprintf('3. Balanco Eletrico, Balanco Termico e Adaptacao Conceitual do PUE Orbital:\n');

% Livro-razao eletrico no barramento de 400 V (Regime Nominal Pleno):
P_IT_nom = 40000.00;          % [W] Potencia util de processamento neural (80 dies x 500 W) [VERIFICAVEL]
P_pump_nom = 200.00;          % [W] Bomba centrifuga MPFL amonia [PROVISORIO]
P_loss_GaN_nom = 1025.64;     % [W] Perdas conversor GaN (eta = 97.5%: 40000/0.975 - 40000) [PROVISORIO]
P_loss_BESS_nom = 65.00;      % [W] Resistencia interna em flutuacao [CALCULADO]
P_avionics_nom = 350.00;      % [W] Computador de bordo, ADCS e telemetria [PROVISORIO]
P_OISL_el_nom = 440.00;       % [W] 4 transceivers laser no barramento [VERIFICAVEL]

% Total eletrico drenado do barramento:
P_total_el_nom = P_IT_nom + P_pump_nom + P_loss_GaN_nom + P_loss_BESS_nom + P_avionics_nom + P_OISL_el_nom;
% P_total_el_nom = 42.080,64 W

% Balanco termico interno (dissipacao de calor a ser rejeitada pelo radiador):
eta_eo = 0.22;                % Rendimento eletro-optico do transmissor laser [PROVISORIO]
P_opt_irradiada = P_OISL_el_nom * eta_eo; % 96.80 W deixam a fronteira do satelite via feixe optico [CALCULADO]
Q_OISL_calor = P_OISL_el_nom * (1.0 - eta_eo); % 343.20 W retidos como calor [CALCULADO]
Q_int_total = P_total_el_nom - P_opt_irradiada;
% Q_int_total = 42.080,64 - 96,80 = 41.983,84 W (~41.984 W)

% Adaptacao Conceitual do PUE para Data Center Orbital (sem conformidade normativa):
PUE_orbital_nom = P_total_el_nom / P_IT_nom; % 42080.64 / 40000 = 1.0520

% Pior Caso Quente Degradado EOL (WCH):
P_pump_wch = 400.00;
P_loss_GaN_wch = 4444.44;     % eta = 90.0% (40000/0.90 - 40000)
P_loss_BESS_wch = 120.00;
P_avionics_wch = 450.00;
P_OISL_el_wch = 440.00;

P_total_el_wch = P_IT_nom + P_pump_wch + P_loss_GaN_wch + P_loss_BESS_wch + P_avionics_wch + P_OISL_el_wch;
PUE_orbital_wch = P_total_el_wch / P_IT_nom; % 45854.44 / 40000 = 1.1464

% Rendimento computacional especifico eta_comp (TOPS/W eletrico):
TOPS_total = 80 * 250.0; % 20.000 TOPS
eta_comp_nom = TOPS_total / P_total_el_nom;
eta_comp_wch = TOPS_total / P_total_el_wch;

fprintf('   - Livro-Razao Eletrico Nominal (400 V):\n');
fprintf('     P_IT: %.2f W | Auxiliares+Perdas: %.2f W | Total Eletrico: %.2f W\n', ...
        P_IT_nom, P_total_el_nom - P_IT_nom, P_total_el_nom);
fprintf('     -> PUE_orbital nominal: %.4f (Adaptacao conceitual proposta)\n', PUE_orbital_nom);
fprintf('     -> Rendimento eta_comp: %.3f TOPS/W\n\n', eta_comp_nom);

fprintf('   - Balanco Termico Interno Liquido:\n');
fprintf('     Total Eletrico: %.2f W - Feixe Laser Emitido: %.2f W = Calor Interno: %.2f W (~41.984 W)\n\n', ...
        P_total_el_nom, P_opt_irradiada, Q_int_total);

fprintf('   - Pior Caso Quente Degradado EOL (WCH):\n');
fprintf('     P_IT: %.2f W | Auxiliares+Perdas: %.2f W | Total Eletrico: %.2f W\n', ...
        P_IT_nom, P_total_el_wch - P_IT_nom, P_total_el_wch);
fprintf('     -> PUE_orbital WCH:     %.4f\n', PUE_orbital_wch);
fprintf('     -> Rendimento eta_comp: %.3f TOPS/W\n\n', eta_comp_wch);

fprintf('==============================================================================\n');
fprintf('  SIMULACAO 4.6 CONCLUIDA COM SUCESSO (ANALISE PARAMETRICA)\n');
fprintf('==============================================================================\n');
