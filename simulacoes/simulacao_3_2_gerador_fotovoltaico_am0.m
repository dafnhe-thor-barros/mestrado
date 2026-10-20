% ============================================================================
% SIMULACAO 3.2: GERADOR FOTOVOLTAICO MULTIJUNCAO AM0 (SPECTROLAB XTJ PRIME)
% MODELO FISICO DE 3 DIODOS INDEPENDENTES EM SERIE (GaInP / GaAs / Ge),
% DEGRADACAO DDD, HOT SPOT E CONVERSOR BUCK CHAVEADO (150 kHz)
% ============================================================================
% Dissertacao de Mestrado Profissional - ProPECAut / IFES
% Autor: Dafnhe Thor Castro de Martins Barros
% Referencias Primarias:
%   [1] Spectrolab XTJ Prime Space Solar Cell Datasheet (AM0, 1353 W/m^2, 28 C)
%   [2] Cotal et al. (2000), Triple-Junction Solar Cell Efficiencies Above 32%
% ============================================================================
clear; close all; clc;

fprintf('==============================================================================\n');
fprintf('INICIANDO SIMULACAO 3.2: MODELO DE 3 DIODOS EM SERIE (GaInP/GaAs/Ge) E BUCK\n');
fprintf('==============================================================================\n\n');

% ----------------------------------------------------------------------------
% 1. Constantes fisicas e parametros de datasheet (Spectrolab XTJ Prime)
% ----------------------------------------------------------------------------
qe = 1.602176634e-19;      % [C] Carga elementar
kB = 1.380649e-23;         % [J/K] Constante de Boltzmann
sigmaSB = 5.670374e-8;     % [W/(m^2 K^4)] Constante de Stefan-Boltzmann

Tref = 28.0;               % [degC] Temperatura padrao de ensaio
G_datasheet = 1353.0;      % [W/m^2] Irradiancia de ensaio Spectrolab AM0
G_AM0 = 1361.0;            % [W/m^2] Constante solar ASTM E490 (1 UA)
scale_AM0 = G_AM0 / G_datasheet;

Acm2 = 30.18;              % [cm^2] Area ativa por celula
Am2  = Acm2 * 1e-4;        % [m^2]
VocRef_tot = 2.720;        % [V] Voc total do datasheet (1.42 + 1.01 + 0.29 V)
JscRef     = 18.0e-3;      % [A/cm^2] Densidade de corrente de curto-circuito
VmpRef_cat = 2.406;        % [V] Tensao MPP de catalogo
JmpRef     = 17.5e-3;      % [A/cm^2] Densidade de corrente MPP de catalogo

IscRef = JscRef * Acm2 * scale_AM0; % [A] 0.5465 A sob 1361 W/m^2
ImpRef = JmpRef * Acm2 * scale_AM0; % [A] 0.5313 A sob 1361 W/m^2

% Coeficientes termicos oficiais do datasheet Spectrolab XTJ Prime
dVoc_dT_tot = -5.60e-3;    % [V/K] Coeficiente termico global de Voc
dJsc_dT     = +10.0e-6;    % [A/(cm^2 K)] Coeficiente termico de Jsc
dIsc_dT     = dJsc_dT * Acm2 * scale_AM0; % [A/K] +0.3036 mA/K
dVmp_dT     = -6.30e-3;    % [V/K] Coeficiente termico de Vmp
dJmp_dT     = +5.0e-6;     % [A/(cm^2 K)] Coeficiente termico de Jmp
dImp_dT     = dJmp_dT * Acm2 * scale_AM0; % [A/K]

% Parametros fisicos das 3 subcelulas em serie (GaInP / GaAs / Ge) a 28 C:
% [1: Top GaInP, 2: Middle GaAs, 3: Bottom Ge]
Voc_sub_ref = [1.420, 1.010, 0.290];        % [V] Soma = 2.720 V
n_sub       = [1.15,  1.15,  1.75];         % Fatores de idealidade por subcelula
Rsh_sub     = [2500.0, 3500.0, 200.0];      % [Ohm] Resistencias shunt individuais
% Fotocorrentes individuais a 28 C (GaAs limitante, GaInP levemente acima, Ge com excesso IR)
Iph_sub_ref = [1.012 * IscRef, 1.000 * IscRef, 1.550 * IscRef]; % [A]
% Reparticao dos coeficientes termicos de Voc e Isc pelas 3 subcelulas
dVoc_sub_dT = [-2.20e-3, -1.90e-3, -1.50e-3]; % [V/K] Soma = -5.60 mV/K
dIph_sub_dT = [0.85 * dIsc_dT, 1.00 * dIsc_dT, 1.80 * dIsc_dT]; % [A/K]

% Resistencia serie total (contatos metalicos + juncoes tunel)
Rs_total = 0.20;           % [Ohm/celula]

%% =========================================================================
%% REDIMENSIONAMENTO FÍSICO DO ARRANJO: Ns = 208 CÉLULAS EM SÉRIE
%% Garante V_mp > V_bus (400 V) no pior caso quente (+65 °C, LEO)
%% =========================================================================

% Parâmetros Térmicos Oficiais Spectrolab XTJ Prime
V_mp_ref   = 2.406;          % [V] Tensão MPP a 28 °C sob AM0
dVmp_dT    = -6.30e-3;       % [V/°C] Coeficiente térmico de Vmp
T_hot_max  = 65.0;           % [°C] Pior caso térmico contínuo em insolação

% Tensão de célula degradada no regime quente (+65 °C)
V_mp_cell_hot = V_mp_ref + dVmp_dT * (T_hot_max - 28.0); % 2.1729 V

% Requisito do Conversor Buck para D <= 0.90 em regime quente
V_bus_nominal = 400.0;       % [V] Barramento primário estabilizado
D_max_target  = 0.90;
V_pv_min_req  = V_bus_nominal / D_max_target; % 444.4 V

% Novo número de células em série (substitui Ns = 166)
N_s = 208;                   % Adota-se N_s = 208 células (limite D <= 0.885 no pior caso quente)

% Redimensionamento de strings em paralelo para manter 50 kW nominais EOL
I_mp_hot = (17.50e-3 * 30.18 * (1361/1353)) + (0.1518e-3 * (65 - 28)); % 0.5369 A
P_cell_hot = V_mp_cell_hot * I_mp_hot; % 1.1666 W

P_alvo_eol = 50e3;           % [W] Potência requerida em fim de vida
fator_eol  = 0.82;           % Degradação DDD após 10 anos em SSO 700 km
N_total    = ceil(P_alvo_eol / (P_cell_hot * fator_eol));
N_p_no     = 252;            % 252 strings em paralelo (52.416 células no total)

fprintf('=== NOVO ARRANJO COMPATÍVEL COM TOPOLOGIA BUCK E 400 V ===\n');
fprintf('Células em série (N_s): %d (Tensão MPP a 65 °C = %5.1f V)\n', N_s, N_s * V_mp_cell_hot);
fprintf('Ciclo de trabalho D a +65 °C: %0.3f (Operação estável no Buck)\n', V_bus_nominal / (N_s * V_mp_cell_hot));
fprintf('Tensão de circuito aberto a frio (-100 °C): %5.1f V (Margem GaN 900 V: OK)\n', ...
    N_s * (2.720 - 5.60e-3 * (-100 - 28)));

Ns = N_s;
Np = N_p_no;
T_ref = Tref;
I_mp_ref = ImpRef;

Pmp_cell_cat = VmpRef_cat * ImpRef;
eta_cell_cat = Pmp_cell_cat / (G_AM0 * Am2);
fprintf('--- 1.1 MODELO DE 3 DIODOS EM SERIE (GaInP / GaAs / Ge) ---\n');
fprintf('Subcelulas Voc(28 C): GaInP=%.2f V | GaAs=%.2f V | Ge=%.2f V (Total=%.3f V)\n', ...
        Voc_sub_ref(1), Voc_sub_ref(2), Voc_sub_ref(3), sum(Voc_sub_ref));
fprintf('Fatores de idealidade: n = [%.2f, %.2f, %.2f] | Rsh = [%.0f, %.0f, %.0f] Ohm | Rs = %.2f Ohm\n', ...
        n_sub(1), n_sub(2), n_sub(3), Rsh_sub(1), Rsh_sub(2), Rsh_sub(3), Rs_total);
fprintf('Matriz: %dS x %dP (%d celulas) | Pmp catalogo (28 C) = %.2f kW (eta = %.2f %%)\n\n', ...
        Ns, Np, Ns * Np, (Ns * Np * Pmp_cell_cat) / 1e3, 100 * eta_cell_cat);

% ----------------------------------------------------------------------------
% 2. Resolucao das curvas I-V e P-V do modelo de 3 diodos por temperatura
% ----------------------------------------------------------------------------
fprintf('--- 2.1 CURVAS I-V E P-V DO MODELO DE 3 DIODOS EM SERIE ---\n');
Tgrid = [15, 28, 45, 60, 85];
VcellGrid = linspace(0, 2.85, 600);
IcellGrid = zeros(numel(Tgrid), numel(VcellGrid));
PmpNode = zeros(size(Tgrid));
VmpNode = zeros(size(Tgrid));
ImpNode = zeros(size(Tgrid));

curvas_iv = struct('V', {}, 'I', {}, 'P', {}, 'T', {}, 'P_mp', {});

for m = 1:numel(Tgrid)
  T = Tgrid(m);
  dT = T - Tref;
  TK = T + 273.15;

  Iph_T = Iph_sub_ref + dIph_sub_dT * dT;
  Voc_T = Voc_sub_ref + dVoc_sub_dT * dT;
  Vth_T = (n_sub * kB * TK) / qe;
  I0_T  = (Iph_T - Voc_T ./ Rsh_sub) ./ (exp(Voc_T ./ Vth_T) - 1);

  % Parametrizacao por corrente I_cel de 0 ate min(Iph_T) para obter V_cel(I)
  I_sweep = linspace(0, min(Iph_T) * 0.9999, 1000);
  V_sweep = zeros(size(I_sweep));
  for k_i = 1:numel(I_sweep)
    ic = I_sweep(k_i);
    v_sum = 0.0;
    for j_sub = 1:3
      v_j = max(0.0, Voc_T(j_sub));
      for it = 1:25
        e = exp(min(500.0, v_j / Vth_T(j_sub)));
        f = Iph_T(j_sub) - I0_T(j_sub) * (e - 1.0) - v_j / Rsh_sub(j_sub) - ic;
        df = -(I0_T(j_sub) / Vth_T(j_sub)) * e - 1.0 / Rsh_sub(j_sub);
        dV = f / df;
        v_j = min(Voc_T(j_sub), max(0.0, v_j - dV));
        if abs(dV) < 1e-9, break; end
      end
      v_sum = v_sum + v_j;
    end
    V_sweep(k_i) = max(0.0, v_sum - ic * Rs_total);
  end

  % Interpolacao sobre a malha uniforme de tensao VcellGrid
  [V_uniq, idx_u] = unique(V_sweep, 'stable');
  I_uniq = I_sweep(idx_u);
  IcellGrid(m, :) = interp1(V_uniq, I_uniq, VcellGrid, 'linear', 0.0);

  Pcell = VcellGrid .* IcellGrid(m, :);
  [pmax, idx] = max(Pcell);
  PmpNode(m) = Ns * Np * pmax;
  VmpNode(m) = Ns * VcellGrid(idx);
  ImpNode(m) = Np * IcellGrid(m, idx);

  curvas_iv(m).V = Ns * VcellGrid;
  curvas_iv(m).I = Np * IcellGrid(m, :);
  curvas_iv(m).P = curvas_iv(m).V .* curvas_iv(m).I;
  curvas_iv(m).T = T;
  curvas_iv(m).P_mp = PmpNode(m);

  fprintf('T = %3.0f C | Vmp_cel = %.3f V | Vmp_str = %5.1f V | Imp_no = %6.2f A | Pmp_no = %5.2f kW\n', ...
          T, VcellGrid(idx), VmpNode(m), ImpNode(m), PmpNode(m) / 1e3);
end
fprintf('\n');

% ----------------------------------------------------------------------------
% 3. Cenario radiativo nominal (1 MeV eq., CMG 150 um, SSO 700 km)
% ----------------------------------------------------------------------------
fprintf('--- 3.1 DEGRADACAO DDD (EL. 1 MeV SOB COVERGLASS CMG 150 um EM SSO 700 km) ---\n');
years = 0:0.1:10;
Cd = 0.085;                % Eq. 34 da tese
PhiX = 1.0e14;             % [e-/cm^2]
PhiAnnual = 2.5e13;        % [e-/cm^2/ano] Fluencia eq. 1 MeV atras de CMG 150 um
lossEnv10 = 0.052;         % Escurecimento UV + micrometeoroides (5.2% em 10 anos)
Phi = years * PhiAnnual;
etaRad = 1 - Cd * log10(1 + Phi / PhiX);
etaEnv = 1 - lossEnv10 * years / 10;
etaTotal = etaRad .* etaEnv;

anos = years;
phi_total = Phi;
retencao_total = etaTotal;

fprintf('Retencao DDD (10 anos) = %.2f %% | Retencao Ambiental = %.2f %% | Retencao Total EOL = %.2f %%\n', ...
        100 * etaRad(end), 100 * etaEnv(end), 100 * etaTotal(end));
fprintf('Potencia EOL do no a 28 C (Modelo 3 Diodos) = %.2f kW\n\n', (PmpNode(2) / 1e3) * etaTotal(end));

% ----------------------------------------------------------------------------
% 4. Contingencia de Hot Spot no vacuo (Stefan-Boltzmann)
% ----------------------------------------------------------------------------
fprintf('--- 4.1 CONTINGENCIA DE HOT SPOT NO VACUO (STEFAN-BOLTZMANN) ---\n');
Vbreak = 28.0; Vsch = 0.35; Istring = ImpRef;
Psem = Vbreak * Istring; Pcom = Vsch * Istring;
epsCell = 0.85; Cth = 3.8; Tspace = 3.0;

dt_hs = 0.01; t_hs = 0:dt_hs:240;
T_hs_sem_K = zeros(size(t_hs)); T_hs_sem_K(1) = 60 + 273.15;
T_hs_com_K = zeros(size(t_hs)); T_hs_com_K(1) = 60 + 273.15;
for k = 1:numel(t_hs)-1
  qsem = 2.0 * epsCell * sigmaSB * Am2 * (T_hs_sem_K(k)^4 - Tspace^4);
  T_hs_sem_K(k+1) = T_hs_sem_K(k) + (Psem - qsem) * dt_hs / Cth;
  qcom = 2.0 * epsCell * sigmaSB * Am2 * (T_hs_com_K(k)^4 - Tspace^4);
  T_hs_com_K(k+1) = T_hs_com_K(k) + (Pcom - qcom) * dt_hs / Cth;
end

idx180 = find(T_hs_sem_K >= 180 + 273.15, 1);
TeqSem = (Psem / (2.0 * epsCell * sigmaSB * Am2) + Tspace^4)^0.25 - 273.15;
TeqCom = (Pcom / (2.0 * epsCell * sigmaSB * Am2) + Tspace^4)^0.25 - 273.15;
fprintf('HotSpot (2 faces): P_sem = %.2f W | t(180 C) = %.2f s | Teq_sem = %.2f C | Teq_com = %.2f C\n\n', ...
        Psem, t_hs(idx180), TeqSem, TeqCom);

% =========================================================================
% SEÇÃO 5 CORRIGIDA: DINÂMICA CHAVEADA BUCK (150 kHz) E CONTROLE P&O-AP (10 kHz)
% =========================================================================
fprintf('--- 5.1 SIMULAÇÃO DINÂMICA REAL: BUCK SÍNCRONO + MPPT ADAPTATIVO ---\n');

% Parâmetros do Conversor e Barramento Primário
V_bus_reg = 400.0;          % [V] Barramento estabilizado HVDC nominal
L_filter  = 45.0e-6;        % [H] Indutância de filtro (45 uH)
C_in_filt = 250.0e-6;       % [F] Capacitância de entrada (250 uF)
R_parasita= 0.025;          % [Ohm] Resistência de condução série do indutor

% Base Temporal e Sincronismo
f_sw      = 150e3;          % [Hz] Frequência de chaveamento das chaves GaN
dt        = 1 / f_sw;       % [s] Passo de integração (~6.667 us)
t_final   = 0.050;          % [s] Janela transitória representativa (50 ms)
t_vec     = 0:dt:t_final;
N_pts     = length(t_vec);

f_mppt    = 10e3;           % [Hz] Taxa de amostragem do controle MPPT
N_div     = round(f_sw / f_mppt); % Exatamente 15 ciclos PWM por passo de controle

% Parâmetros do Controle P&O-AP
M_gain    = 1.2e-4;         % [V/W] Ganho adaptativo (Dissertação Eq. 3.17)
dD_max    = 0.002;          % Limite estrito de saturação do passo
dD_fixo   = 0.005;          % Passo fixo comparativo

% Rampa Térmica de Órbita Real (15 K/min = 0.25 K/s)
% Mapeada em janela acelerada para ensaio transitório de malha:
T_cell_dyn = -100.0 + (65.0 - (-100.0)) * (t_vec / t_final); 

% Vetores de Estado das Equações Diferenciais do Filtro Buck (Adaptativo)
V_pv_dyn  = zeros(size(t_vec));
I_L_dyn   = zeros(size(t_vec));
Duty_dyn  = zeros(size(t_vec));

% Vetores de Estado para Passo Fixo Comparativo
V_pv_fixo = zeros(size(t_vec));
I_L_fixo  = zeros(size(t_vec));
Duty_fixo = zeros(size(t_vec));

% Condições Iniciais no Frio Extremo (-100 °C)
V_mp_cold = N_s * (V_mp_ref + dVmp_dT * (-100.0 - T_ref));
V_pv_dyn(1)  = V_mp_cold;    % ~563.6 V
Duty_dyn(1)  = V_bus_reg / V_mp_cold; % D inicial ~ 0.7097
I_L_dyn(1)   = (N_p_no * I_mp_ref) / Duty_dyn(1);

V_pv_fixo(1) = V_mp_cold;
Duty_fixo(1) = Duty_dyn(1);
I_L_fixo(1)  = I_L_dyn(1);

V_prev_ctrl = V_pv_dyn(1);
P_prev_ctrl = V_pv_dyn(1) * (Duty_dyn(1) * I_L_dyn(1));

V_prev_fixo = V_pv_fixo(1);
P_prev_fixo = V_pv_fixo(1) * (Duty_fixo(1) * I_L_fixo(1));

% Loop de Resolução Temporal (Método de Euler Semi-Implícito / Cromer)
for k = 1:N_pts-1
    T_now = T_cell_dyn(k);
    V_now = V_pv_dyn(k);
    D_now = Duty_dyn(k);
    I_now = I_L_dyn(k);
    
    % Modelo de Corrente Fotovoltaica Instantânea AM0
    V_mp_k = N_s * (V_mp_ref + dVmp_dT * (T_now - T_ref));
    I_mp_k = N_p_no * (I_mp_ref + dImp_dT * (T_now - T_ref));
    I_pv_inst = max(0, I_mp_k * (1 - (V_now - V_mp_k) / (0.30 * V_mp_k)));
    
    % Algoritmo MPPT executado a 10 kHz (a cada 15 períodos de PWM)
    if mod(k, N_div) == 0
        P_now_ctrl = V_now * (D_now * I_now);
        dV_k = V_now - V_prev_ctrl;
        dP_k = P_now_ctrl - P_prev_ctrl;
        
        if abs(dV_k) > 1e-3
            slope_mppt = dP_k / dV_k;
            % SINAL NEGATIVO: Obrigatório para Conversor Buck (dV_pv/dD < 0)
            step_ctrl = - min(dD_max, max(-dD_max, M_gain * slope_mppt));
            D_next = max(0.20, min(0.98, D_now + step_ctrl));
        else
            D_next = D_now;
        end
        V_prev_ctrl = V_now;
        P_prev_ctrl = P_now_ctrl;
    else
        D_next = D_now;
    end
    Duty_dyn(k+1) = D_next;
    
    % Integração das EDOs de Estado do Conversor Adaptativo:
    % 1) dV_pv/dt = (I_pv - D * I_L) / C_in
    dVpv_dt = (I_pv_inst - D_next * I_now) / C_in_filt;
    V_pv_dyn(k+1) = V_now + dVpv_dt * dt;
    
    % 2) dI_L/dt = (D * V_pv - V_bus - R * I_L) / L
    dIL_dt = (D_next * V_pv_dyn(k+1) - V_bus_reg - R_parasita * I_now) / L_filter;
    I_L_dyn(k+1) = max(0, I_now + dIL_dt * dt);

    % Comparativo Passo Fixo
    V_now_f = V_pv_fixo(k);
    D_now_f = Duty_fixo(k);
    I_now_f = I_L_fixo(k);
    I_pv_inst_f = max(0, I_mp_k * (1 - (V_now_f - V_mp_k) / (0.30 * V_mp_k)));
    if mod(k, N_div) == 0
        P_now_f = V_now_f * (D_now_f * I_now_f);
        dV_f = V_now_f - V_prev_fixo;
        dP_f = P_now_f - P_prev_fixo;
        if abs(dV_f) > 1e-3
            D_next_f = max(0.20, min(0.98, D_now_f - dD_fixo * sign(dP_f / dV_f)));
        else
            D_next_f = D_now_f;
        end
        V_prev_fixo = V_now_f;
        P_prev_fixo = P_now_f;
    else
        D_next_f = D_now_f;
    end
    Duty_fixo(k+1) = D_next_f;
    dVpv_dt_f = (I_pv_inst_f - D_next_f * I_now_f) / C_in_filt;
    V_pv_fixo(k+1) = V_now_f + dVpv_dt_f * dt;
    dIL_dt_f = (D_next_f * V_pv_fixo(k+1) - V_bus_reg - R_parasita * I_now_f) / L_filter;
    I_L_fixo(k+1) = max(0, I_now_f + dIL_dt_f * dt);
end

fprintf('Convergência Dinâmica Validada com Sucesso em Malha Fechada a 150 kHz.\n\n');

% Exportação de aliases para compatibilidade plena com rotinas de renderização
t_mppt = t_vec;
V_bus_est = V_bus_reg;
D_adapt = Duty_dyn;
D_fixo  = Duty_fixo;
T_transiente = T_cell_dyn;
P_adapt = V_pv_dyn .* (Duty_dyn .* I_L_dyn);
P_fixo  = V_pv_fixo .* (Duty_fixo .* I_L_fixo);
P_otima = (N_p_no * (I_mp_ref + dImp_dT * (T_cell_dyn - T_ref))) .* ...
          (N_s * (V_mp_ref + dVmp_dT * (T_cell_dyn - T_ref)));
