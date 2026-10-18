% ==============================================================================
% SIMULAÇÃO 4.6 (CORRIGIDA E RECONCILIADA COM A DISSERTAÇÃO DE MESTRADO)
% FADIGA TERMOMECÂNICA DE SOLDAS BGA SAC305 (ENGELMAIER) E MÉTRICA PUE ORBITAL
% ==============================================================================
% Programa de Pós-Graduação em Engenharia de Controle e Automação - ProPECAut / IFES
% Autor da Dissertação: Dafnhê Thor Castro de Martins Barros
% Auditoria Técnica: Reconciliação Numérica com as Tabelas 14 e 15 e Equações 88-93
% ==============================================================================
clear; clc; close all;
fprintf('==============================================================================\n');
fprintf('INICIANDO SIMULAÇÃO 4.6 CORRIGIDA: FADIGA BGA SAC305 E INDICADOR PUE ORBITAL\n');
fprintf('==============================================================================\n\n');
%% 1. PARÂMETROS TERMOMECÂNICOS DAS JUNTAS BGA SAC305 (TABELA 14 DA DISSERTAÇÃO)
L_DNP        = 14.0e-3;      % [m] Distância ao ponto neutro crítico (die 28x28 mm)
h_solder     = 0.18e-3;      % [m] Altura da esfera de solda SAC305 padrão espacial
delta_alpha  = 14.0e-6;      % [K^-1] Diferencial de expansão térmica (substrato BT vs silício)
F_conform    = 0.90;         % Fator de conformidade mecânica da placa (IPC-SM-785)
k_att        = 0.85;         % Fator de atenuação térmica da junção até a esfera (Eq. 88)
two_eps_f    = 2.887;        % Ductilidade cíclica efetiva SAC305 (Syed, 2004; Gao et al., 2022)
% Coeficientes de Calibração de Engelmaier (Norma IPC-9701A / Syed, 2004)
c_0          = -0.470;       % Coeficiente base de Engelmaier
alpha_T      = -7.0e-4;      % [°C^-1] Sensibilidade térmica a creep/relaxação
beta_f       = +2.0e-2;      % Sensibilidade à frequência cíclica diária
%% 2. MATRIZ FATORIAL 2x2: ISOLAMENTO DO EFEITO ORBITAL E CONTROLE TEMS (TABELA 15)
% Cenários Analisados:
% 1: SSO + TEMS (Ponto de Projeto Nominal)
% 2: SSO sem TEMS (Controle Bang-Bang / Transientes Abruptos)
% 3: LEO + TEMS (Órbita de Estresse com Mitigação Ativa)
% 4: LEO sem TEMS (Colapso Prematuro por Ciclagem de Eclipse)
cenarios_nomes = { ...
    'SSO + TEMS (Projeto)', ...
    'SSO sem TEMS', ...
    'LEO + TEMS', ...
    'LEO sem TEMS (Falha)' ...
};
% Entradas Estritamente Reconciliadas com a Tabela 15 da Dissertação
delta_T_j       = [7.0, 22.5, 18.2, 48.0];  % [°C] Excursão térmica na junção do chip
delta_T_cp      = [0.0, 0.0, 0.0, 0.0];     % [°C] Excursão na placa fria (isolada pelo LHP)
T_mean_solder_T = [60.0, 60.0, 50.0, 50.0]; % [°C] Temperatura média na junta de solda
f_cpd_vec       = [14.58, 14.58, 15.22, 15.22]; % [ciclos/dia] Frequência orbital diária
N_f_calc    = zeros(1, 4);
vida_anos   = zeros(1, 4);
gamma_eff   = zeros(1, 4);
c_engel     = zeros(1, 4);
margem_10a  = zeros(1, 4);
fprintf('--- 2.1 AVALIAÇÃO ANALÍTICA DE FADIGA DE SOLDAS BGA SAC305 (ENGELMAIER) ---\n');
for idx = 1:4
    % Atenuação Térmica da Junção até a Esfera BGA (Equação 88)
    dT_solder = k_att * delta_T_j(idx) + (1 - k_att) * delta_T_cp(idx);
    T_mean_s = T_mean_solder_T(idx);
    f_cpd    = f_cpd_vec(idx);
    
    % Deformação Efetiva de Cisalhamento por Termofluência (Equação 89)
    gamma_e = F_conform * (L_DNP * delta_alpha * dT_solder) / h_solder;
    gamma_eff(idx) = gamma_e;
    
    % Expoente Visco-Plástico de Engelmaier Dependente de Temperatura e Frequência
    c_val = c_0 + alpha_T * T_mean_s + beta_f * log(1 + f_cpd);
    c_engel(idx) = c_val;
    
    % Número de Ciclos Térmicos até a Falha (N_f)
    N_f = 0.5 * (gamma_e / two_eps_f)^(1 / c_val);
    N_f_calc(idx) = N_f;
    
    % Vida Útil Equivalente em Anos de Voo Espacial
    ciclos_por_ano = f_cpd * 365.25;
    t_anos = N_f / ciclos_por_ano;
    vida_anos(idx) = t_anos;
    
    % Margem de Segurança Frente ao Requisito de Missão de 10 Anos
    margem = ((t_anos - 10.0) / 10.0) * 100;
    margem_10a(idx) = margem;
    
    fprintf('Cenário %d: %-22s\n', idx, cenarios_nomes{idx});
    fprintf('  Delta_T Junção: %4.1f °C | Delta_T Solda: %4.2f °C | Deformação gamma_eff: %7.2e\n', ...
        delta_T_j(idx), dT_solder, gamma_e);
    fprintf('  Expoente c: %6.3f | Ciclos até Falha: %9.0f | Vida Útil: %5.2f anos\n', ...
        c_val, N_f, t_anos);
    
    if t_anos >= 10.0
        fprintf('  -> Status: [CONFORME] Margem de Segurança: +%5.1f %% (Atende aos 10 Anos)\n\n', margem);
    else
        fprintf('  -> Status: [FALHA PREMATURA] Margem: %5.1f %% (Colapso por Fadiga Termomecânica)\n\n', margem);
    end
end
%% 3. INDICADOR PUE ORBITAL E RENDIMENTO COMPUTACIONAL (EQUAÇÕES 90 A 93)
% Potência Computacional Útil dos Aceleradores de IA
P_IT = 40000.0; % [W] 40 kW nominal
% --- Cenário 1: Regime Nominal de Cruzeiro em Órbita SSO (Equação 91) ---
P_pump_nom    = 200.0;  % [W] Bomba mecânica MPFL / LHP
P_heaters_nom = 0.0;    % [W] Aquecedores em repouso durante cruzeiro iluminado
P_cooling_nom = P_pump_nom + P_heaters_nom;
P_loss_GaN_nom  = 1025.0; % [W] Perdas nos conversores GaN a 97.5% de rendimento
P_loss_BESS_nom = 65.0;   % [W] Perdas ôhmicas de flutuação no pack de baterias
P_elec_loss_nom = P_loss_GaN_nom + P_loss_BESS_nom;
P_avionics_nom = 350.0; % [W] Computador de bordo, TT&C e subsistema ADCS
P_OISL_nom     = 343.0; % [W] Consumo dos 4 transceivers ópticos laser
P_aux_nom      = P_avionics_nom + P_OISL_nom;
P_total_nom = P_IT + P_cooling_nom + P_elec_loss_nom + P_aux_nom;
PUE_orbital_nom = P_total_nom / P_IT;
% --- Cenário 2: Pior Caso Quente Degradado EOL (WCH - Equação 92) ---
P_pump_wch    = 400.0;  % [W] Bomba em vazão máxima de contingência
P_heaters_wch = 0.0;
P_cooling_wch = P_pump_wch + P_heaters_wch;
P_loss_GaN_wch  = 4444.0; % [W] Conversores degradados a 90.0% de rendimento
P_loss_BESS_wch = 120.0;  % [W] Ciclagem em alta corrente sob R0 degradada EOL
P_elec_loss_wch = P_loss_GaN_wch + P_loss_BESS_wch;
P_avionics_wch = 450.0; % [W] ADCS em atuação plena de correção de atitude
P_OISL_wch     = 343.0;
P_aux_wch      = P_avionics_wch + P_OISL_wch;
P_total_wch = P_IT + P_cooling_wch + P_elec_loss_wch + P_aux_wch;
PUE_orbital_wch = P_total_wch / P_IT;
% --- Rendimento Computacional Específico (Equação 93) ---
% Arranjo do Nó: 80 dies operando a 2.1 GHz com 250 TOPS nominais por die
N_dies = 80;
TOPS_per_die = 250.0;
TOPS_total_no = N_dies * TOPS_per_die; % 20.000 TOPS por satélite
eta_comp_nom = TOPS_total_no / P_total_nom; % [TOPS/W]
eta_comp_wch = TOPS_total_no / P_total_wch; % [TOPS/W]
fprintf('--- 3.1 INDICADORES DE EFICÁCIA ENERGÉTICA E RENDIMENTO COMPUTACIONAL ---\n');
fprintf('Potência Total Nominal do Nó: %5.0f W | Carga Útil de IA: %5.0f W\n', P_total_nom, P_IT);
fprintf('PUE Orbital Nominal (Cruzeiro SSO): %6.4f (Aderência com Eq. 91: 1.050)\n', PUE_orbital_nom);
fprintf('PUE Orbital Pior Caso Quente (WCH): %6.4f (Aderência com Eq. 92: 1.144)\n\n', PUE_orbital_wch);
fprintf('Rendimento Computacional Específico Nominal:   %5.3f TOPS/W (20.000 TOPS / 41.983 W)\n', eta_comp_nom);
fprintf('Rendimento Computacional Específico Pior Caso: %5.3f TOPS/W (20.000 TOPS / 45.757 W)\n', eta_comp_wch);
fprintf('==============================================================================\n');
fprintf('SIMULAÇÃO 4.6 EXECUTADA COM SUCESSO: RASTREABILIDADE TOTAL COM A DISSERTAÇÃO.\n');
fprintf('==============================================================================\n');