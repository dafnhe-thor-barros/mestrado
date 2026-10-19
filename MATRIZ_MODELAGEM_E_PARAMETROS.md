# MATRIZ CANÔNICA DE MODELAGEM MATEMÁTICA E CATÁLOGO DE PARÂMETROS
## Dissertação de Mestrado — ProPECAut / IFES Campus Serra
**Título da Pesquisa:** Modelagem Computacional de Gerenciamento Térmico-Energético para Centros de Dados Orbitais em Escala Megawatt  
**Autor:** Dafnhê Thor Castro de Martins Barros  
**Orientador:** Prof. Dr. Flávio Barcelos Braz da Silva  
**Data de Consolidação:** Outubro de 2026  

---

### SUMÁRIO DA MATRIZ
1. [MÓDULO 1: Geração Fotovoltaica Espacial AM0 (EPS-PV)](#módulo-1-geração-fotovoltaica-espacial-am0-eps-pv)
2. [MÓDULO 2: Algoritmos de Rastreamento de Máxima Potência (MPPT)](#módulo-2-algoritmos-de-rastreamento-de-máxima-potência-mppt)
3. [MÓDULO 3: Conversão e Distribuição de Potência DC-DC de Alta Tensão (HVDC)](#módulo-3-conversão-e-distribuição-de-potência-dc-dc-de-alta-tensão-hvdc)
4. [MÓDULO 4: Subsistema de Armazenamento Eletroquímico de Baterias (BESS)](#módulo-4-subsistema-de-armazenamento-eletroquímico-de-baterias-bess)
5. [MÓDULO 5: Algoritmo de Estimação de Estado de Carga (AUKF)](#módulo-5-algoritmo-de-estimação-de-estado-de-carga-aukf)
6. [MÓDULO 6: Telecomunicações e Enlaces Ópticos Inter-Satélite (OISL)](#módulo-6-telecomunicações-e-enlaces-ópticos-inter-satélite-oisl)
7. [MÓDULO 7: Mecânica Orbital e Fluxos Radiativos Ambientais (TCS-Orbital)](#módulo-7-mecânica-orbital-e-fluxos-radiativos-ambientais-tcs-orbital)
8. [MÓDULO 8: Dissipação Eletrotérmica dos Processadores Neurais de IA](#módulo-8-dissipação-eletrotérmica-dos-processadores-neurais-de-ia)
9. [MÓDULO 9: Rede Térmica Nodal Multinodal de Resfriamento no Vácuo](#módulo-9-rede-térmica-nodal-multinodal-de-resfriamento-no-vácuo)
10. [MÓDULO 10: Sistema de Gerenciamento Térmico-Energético Ativo (TEMS / NMPC)](#módulo-10-sistema-de-gerenciamento-térmico-energético-ativo-tems--nmpc)
11. [MÓDULO 11: Confiabilidade, Degradação Cíclica e Fadiga Térmica (10 Anos)](#módulo-11-confiabilidade-degradação-cíclica-e-fadiga-térmica-10-anos)

---

### MÓDULO 1: Geração Fotovoltaica Espacial AM0 (EPS-PV)

#### 1.1 Modelos Matemáticos Implementados
* **Modelo Elétrico de 1 Diodo / 5 Parâmetros:**
  $$I_{pv} = I_{ph} - I_0 \left[ \exp\left( \frac{q (V_{pv} + I_{pv} R_s)}{a k_B T_{cel}} \right) - 1 \right] - \frac{V_{pv} + I_{pv} R_s}{R_{sh}}$$
* **Modelo Comparativo de 2 Diodos (Recombinação de Junção):**
  $$I_{pv} = I_{ph} - I_{01} \left[ \exp\left( \frac{q (V_{pv} + I_{pv} R_s)}{a_1 k_B T_{cel}} \right) - 1 \right] - I_{02} \left[ \exp\left( \frac{q (V_{pv} + I_{pv} R_s)}{a_2 k_B T_{cel}} \right) - 1 \right] - \frac{V_{pv} + I_{pv} R_s}{R_{sh}}$$
* **Dependência Térmica e de Irradiância AM0:**
  $$I_{ph}(T_{cel}, G) = \max\left(0{,}01, [I_{sc,ref} + K_I(T_{cel}-T_{ref})]\frac{G}{G_{datasheet}}\right)$$
  $$V_{oc}(T_{cel}) = V_{oc,ref} + K_V (T_{cel} - T_{ref})$$
  $$I_0(T_{cel},G) = \frac{I_{ph}(T_{cel},G)-V_{oc}(T_{cel})/R_{sh}}{\exp[V_{oc}(T_{cel})/(a_{eq}V_T)]-1}, \quad V_T=k_BT_{cel}/q$$
  A expressão de $I_0$ é a recalibração algébrica usada por `simulacao_3_2_gerador_fotovoltaico_am0.m` para reproduzir o $V_{oc}(T)$ escolhido; não é uma lei física de temperatura validada independentemente.
* **Degradação por Radiação (ajuste paramétrico DDD):**
  $$P_{mp,total}(t) = P_{mp,0} \left[ 1 - C_d \log_{10}\left( 1 + \frac{\Phi_{eq}(t)}{\Phi_x} \right) \right]\left[1-\Delta P_{env}(t/10)\right]$$
  O coeficiente $C_d$ é ajustado a um único ponto de retenção do datasheet Spectrolab XTJ Prime (padrão ECSS); $Phi_x$, as taxas de fluência orbitais e a perda ambiental são premissas provisórias. Esse ajuste pontual não valida a forma da curva entre pontos nem substitui a análise do ambiente radiativo e da blindagem.

#### 1.2 Catálogo de Parâmetros do Painel Solar
| Símbolo | Descrição do Parâmetro | Valor Nominal | Unidade | Referência / Justificativa |
| :--- | :--- | :---: | :---: | :--- |
| $G_{AM0}$ | Irradiância solar de referência a 1 UA | $1361{,}0$ | $\text{W/m}^2$ | ASTM E490; usada como cenário do modelo |
| $G_{datasheet}$ | Irradiância de referência do datasheet XTJ Prime | $1353{,}0$ | $\text{W/m}^2$ | Spectrolab, AM0, 28°C; https://www.spectrolab.com/photovoltaics/XTJ-Prime_Data_Sheet.pdf |
| $\eta_{mp}$ | Eficiência MPP BOL declarada | $30{,}7$ | $\%$ | Catálogo XTJ Prime; não usada no sizing devido à discrepância com $V_{mp}J_{mp}$ arredondados |
| $A_{cel}$ | Área ativa adotada por célula | $30{,}18$ | $\text{cm}^2$ | Premissa geométrica do modelo [PROVISÓRIA] |
| $V_{mp,cel}$ | Tensão de máxima potência por célula | $2{,}406$ | $\text{V}$ | Datasheet Spectrolab |
| $J_{mp,cel}$ | Densidade de corrente de MPP | $17{,}5$ | $\text{mA/cm}^2$ | Datasheet Spectrolab |
| $I_{mp,cel}$ | Corrente de MPP por célula adotada | $0{,}52815$ | $\text{A}$ | Calculada de $J_{mp}A_{cel}$ |
| $V_{oc,cel}$ | Tensão de circuito aberto por célula | $2{,}720$ | $\text{V}$ | Datasheet Spectrolab |
| $J_{sc,cel}$ | Densidade de corrente de curto-circuito | $18{,}0$ | $\text{mA/cm}^2$ | Datasheet Spectrolab |
| $K_V$ | Coeficiente térmico de $V_{oc}$ | $-5{,}60$ | $\text{mV/K}$ | Datasheet Spectrolab |
| $K_I$ | Coeficiente térmico de $I_{sc}$ na área adotada | $+0{,}3018$ | $\text{mA/K}$ | Calculado de $dJ_{sc}/dT$ e $A_{cel}$ |
| $a_{eq}$ | Fator de idealidade equivalente | $3{,}00$ | -- | Parâmetro concentrado [PROVISÓRIO] |
| $R_s$ | Resistência série equivalente por célula | $0{,}080$ | $\Omega$ | Parâmetro do modelo [PROVISÓRIO] |
| $R_{sh}$ | Resistência shunt por célula | $1200$ | $\Omega$ | Parâmetro do modelo [PROVISÓRIO] |
| $N_{s,no},N_{p,no}$ | Série/paralelo por nó | $166,244$ | -- | Script de geração AM0; ponto de catálogo em 400 V |
| $N_{cel,no}$ | Células ativas por nó | $40504$ | -- | $166\times244$ |
| $P_{BOL,no}$ | Potência MPP ideal BOL por nó | $51{,}47$ | $\text{kW}$ | Antes de perdas de conversão |
| $A_{array,cluster}$ | Área ativa total | $2444{,}82$ | $\text{m}^2$ | 20 nós; área por célula adotada |
| $\Phi_x$ | Escala do modelo DDD | $1{,}0\times10^{14}$ | $\text{e}^-/\text{cm}^2$ | Escolha paramétrica [PROVISÓRIA] |
| $C_d$ | Coeficiente de ajuste DDD | $0{,}12851$ | -- | Ajuste pontual à retenção ECSS $0{,}90$ a $5\times10^{14}$; forma contínua provisória |
| $\Phi_x$ | Escala do ajuste DDD | $1{,}0\times10^{14}$ | $\text{e}^-/\text{cm}^2$ | Escolha paramétrica [PROVISÓRIA] |
| $\Delta P_{env}$ | Perda ambiental em dez anos | $5{,}2$ | $\%$ | Premissa sem ensaio/calibração [PROVISÓRIA] |

---

### MÓDULO 2: Algoritmos de Rastreamento de Máxima Potência (MPPT)

#### 2.1 Modelos Matemáticos Implementados
* **Perturbação e Observação Adaptativo (P&O-AP):**
  $$D[k] = D[k-1] + M \cdot \text{sgn}\left( \frac{\Delta P_{pv}}{\Delta V_{pv}} \right) \cdot \left| \frac{\Delta P_{pv}}{\Delta V_{pv}} \right|$$
* **Condutância Incremental (IncCond):**
  $$\frac{dI_{pv}}{dV_{pv}} = -\frac{I_{pv}}{V_{pv}} \quad \text{(Condição no MPP)}$$
  $$\text{Erro de rastreamento: } e_{mpp}[k] = \frac{I_{pv}[k] - I_{pv}[k-1]}{V_{pv}[k] - V_{pv}[k-1]} + \frac{I_{pv}[k]}{V_{pv}[k]}$$

#### 2.2 Catálogo de Parâmetros do MPPT
| Símbolo | Descrição do Parâmetro | Valor Nominal | Unidade | Referência / Justificativa |
| :--- | :--- | :---: | :---: | :--- |
| $f_{mppt}$ | Frequência de amostragem e controle do MPPT | $1{,}0$ | $\text{kHz}$ | Rastreamento dinâmico de gradientes de 10 K/min |
| $\Delta D_{base}$ | Passo de razão cíclica base | $0{,}005$ | -- | Compromisso entre ripple no MPP e velocidade |
| $M$ | Ganho de escala adaptativo de passo | $1{,}2 \times 10^{-4}$ | $\text{V/W}$ | Normalização pelo gradiente de potência |
| $\eta_{mppt}$ | Eficiência estática de rastreamento no MPP | $> 99{,}4$ | $\%$ | Simulação em regime permanente |
| $\tau_{conv}$ | Tempo de convergência na saída de eclipse | $< 18$ | $\text{ms}$ | Resposta a degrau térmico de penumbra |

---

### MÓDULO 3: Conversão e Distribuição de Potência DC-DC de Alta Tensão (HVDC)

#### 3.1 Modelos Matemáticos Implementados
* **Topologia Primária:** Conversor Full-Bridge Isolado com Retificação Síncrona baseada em semicondutores GaN (Nitreto de Gálio) e SiC (Carbeto de Silício).
* **Balanço Analítico de Perdas Eletrotérmicas:**
  $$P_{loss,conv} = P_{cond} + P_{sw} + P_{gate} + P_{core} + P_{cap}$$
  * *Perdas de Condução:* $P_{cond} = 2 I_{rms}^2 R_{on}(T_{j,sw})$ com $R_{on}(T) = R_{on,25} \left( 1 + \alpha_{sc} (T - 25) \right)$
  * *Perdas de Comutação:* $P_{sw} = V_{in} I_{out} f_{sw} \frac{t_r + t_f}{2} + \frac{1}{2} C_{oss} V_{in}^2 f_{sw}$
  * *Perdas Magnéticas (Steinmetz Modificada):* $P_{core} = k_{fe} f_{sw}^\alpha \Delta B^\beta V_e$
* **Critério de Estabilidade de Impedância de Middlebrook:**
  $$\| Z_{out,source}(j\omega) \| \ll \| Z_{in,load}(j\omega) \|, \quad \forall \omega$$
  $$\text{Impedância de entrada da Carga de Potência Constante (CPL): } Z_{in,CPL} = -\frac{V_{bus}^2}{P_{comp}}$$

#### 3.2 Catálogo de Parâmetros dos Conversores DC-DC
| Símbolo | Descrição do Parâmetro | Valor Nominal | Unidade | Referência / Justificativa |
| :--- | :--- | :---: | :---: | :--- |
| $V_{bus,prim}$ | Tensão nominal do barramento primário HVDC | $400$ | $\text{V}$ | Faixa regulada de $324{,}0\text{ V}$ (UVLO) a $442{,}8\text{ V}$ (OV) |
| $V_{bus,sec}$ | Tensão do barramento intermediário do rack | $48$ | $\text{V}$ | Barramento padrão de data centers modernos (Open Compute) |
| $V_{core}$ | Tensão de alimentação lógica das GPUs/TPUs | $0{,}75 - 1{,}10$ | $\text{V}$ | Ajustável via DVFS nos conversores PoL |
| $f_{sw}$ | Frequência de chaveamento dos conversores | $150$ | $\text{kHz}$ | Otimização de volume de magnéticos vs. perdas $P_{sw}$ |
| $R_{on,25}$ | Resistência de condução do chaveador GaN a 25°C | $18$ | $\text{m}\Omega$ | Transistores GaN de alta tensão (650V) |
| $\alpha_{sc}$ | Coeficiente térmico de resistência do GaN | $0{,}0065$ | $\text{K}^{-1}$ | Crescimento de $R_{on}$ sob aquecimento |
| $t_r, t_f$ | Tempos de subida e descida de chaveamento | $12, \, 15$ | $\text{ns}$ | Comutação ultra-rápida |
| $C_{oss}$ | Capacitância de saída parasita do semicondutor | $120$ | $\text{pF}$ | Carga capacitiva comutada |
| $\eta_{conv,nom}$ | Rendimento nominal do conversor primário | $97{,}7$ | $\%$ | Medido em plena carga ($20\text{ kW}$) |
| $C_{bus}$ | Capacitância de amortecimento do barramento 400 V | $4{,}7$ | $\text{mF}$ | Bancos cerâmicos/filme espaciais ($1{,}5\text{ mF}$ em WCA) |
| $I^2t_{lim,sic}$ | Limite adiabático de fusão do SiC MOSFET (SSPC) | $1500$ | $\text{A}^2\text{s}$ | Norma ECSS-E-ST-20-20C / MIL-HDBK-338B |
| $t_{trip,crit}$ | Retardo crítico WCA para atingir limite adiabático | $74{,}8$ | $\mu\text{s}$ | **Conclusão Técnica 1** ($L=4\,\mu\text{H}, C=1{,}5\,\text{mF}$) |
| $V_{spike,max}$ | Tensão máxima de ruptura dielétrica $V_{ds}$ | $1200$ | $\text{V}$ | Tensão de breakdown do semicondutor SiC |
| $t_{trip,v1200}$ | Retardo de abertura que rompe 1200 V sem clamp | $1{,}00$ | $\mu\text{s}$ | **Conclusão Técnica 2** ($1000\text{ ns}$ sob $\tau_{ext}=0{,}5\,\mu\text{s}$) |
| $\tau_{ext,min}$ | Tempo mínimo de extinção suave (*soft turn-off*) | $31{,}1$ | $\mu\text{s}$ | Evita $V_{ds} > 1200\text{ V}$ para $I_{peak} = 6{,}21\text{ kA}$ |

---

### MÓDULO 4: Subsistema de Armazenamento Eletroquímico de Baterias (BESS)

#### 4.1 Modelos Matemáticos Implementados
* **Circuito Equivalente de Thévenin de 2ª Ordem (Dual Polarization - DP):**
  $$V_t(t) = V_{oc}(SoC) - I_L(t) R_0(T_{bat}) - V_{p1}(t) - V_{p2}(t)$$
  $$\frac{dV_{p1}(t)}{dt} = -\frac{1}{R_{p1} C_{p1}} V_{p1}(t) + \frac{1}{C_{p1}} I_L(t)$$
  $$\frac{dV_{p2}(t)}{dt} = -\frac{1}{R_{p2} C_{p2}} V_{p2}(t) + \frac{1}{C_{p2}} I_L(t)$$
* **Dinâmica do Estado de Carga (Coulomb Counting com Eficiência Coulombiana):**
  $$\frac{dSoC(t)}{dt} = -\frac{\eta_c(I_L)}{C_n} I_L(t)$$
* **Cinética Térmica de Arrhenius da Resistência Interna:**
  $$R_0(T_{bat}) = R_{0,ref} \exp\left[ \frac{E_a}{k_B} \left( \frac{1}{T_{bat}} - \frac{1}{T_{ref}} \right) \right]$$
* **Balanço Térmico Interno no Chassi com Manta MLI e Heaters:**
  $$C_{th,bat} \frac{dT_{bat}(t)}{dt} = I_L^2 R_0 + \frac{V_{p1}^2}{R_{p1}} + \frac{V_{p2}^2}{R_{p2}} + P_{heater} - \frac{T_{bat} - T_{chassis}}{R_{th,MLI}}$$

#### 4.2 Catálogo de Parâmetros da Bateria (Química Li-ion NMC Espacial)
| Símbolo | Descrição do Parâmetro | Valor Nominal | Unidade | Referência / Justificativa |
| :--- | :--- | :---: | :---: | :--- |
| $E_{bat,nom}$ | Energia nominal instalada por nó satelital | $9{,}0$ | $\text{kWh}$ | $32{,}4\text{ MJ}$ por satélite ($180\text{ kWh}$ no cluster de 20 nós) |
| $C_{n,Ah}$ | Capacidade nominal equivalente do pack | $22{,}52$ | $\text{Ah}$ | Células NMC espaciais de alta energia específica |
| $V_{bat,nom}$ | Tensão nominal do banco de baterias | $400$ | $\text{V}$ | Associação 108S1P de células NMC ($324\text{ V}$ a $442{,}8\text{ V}$) |
| $R_{0,ref}$ | Resistência ôhmica interna total a 20°C (BOL) | $45$ | $\text{m}\Omega$ | Equivale a $32\text{ m}\Omega$ a 30°C sob $E_a = 0{,}32\text{ eV}$ ($85\text{ m}\Omega$ em EOL) |
| $R_{p1}$ | Resistência de polarização rápida (transferência) | $18$ | $\text{m}\Omega$ | Dinâmica de dupla camada eletroquímica |
| $C_{p1}$ | Capacitância de polarização rápida | $3200$ | $\text{F}$ | Constante de tempo $\tau_1 \approx 57{,}6\text{ s}$ |
| $R_{p2}$ | Resistência de polarização lenta (difusão de fase) | $24$ | $\text{m}\Omega$ | Difusão de íons de lítio no estado sólido |
| $C_{p2}$ | Capacitância de polarização lenta | $22000$ | $\text{F}$ | Constante de tempo $\tau_2 \approx 528\text{ s}$ |
| $E_a$ | Energia de ativação eletroquímica | $0{,}32$ | $\text{eV}$ | Calibração experimental para LiNiMnCoO2 ($5{,}127\times 10^{-20}\text{ J}$) |
| $T_{bat,min}$ | Limite inferior de controle térmico nominal | $+10$ | $^\circ\text{C}$ | Prevenção de deposição de lítio (*lithium plating*) |
| $T_{bat,max}$ | Limite superior de controle térmico nominal | $+25$ | $^\circ\text{C}$ | Minimização de degradação da camada SEI (máximo tolerado: $+35^\circ\text{C}$) |
| $P_{heater,max}$ | Potência máxima dos aquecedores de sobrevivência | $300$ | $\text{W}$ | 3 resistores de 100 W PWM termostatizados na placa fria do pack |
| $R_{th,MLI}$ | Resistência térmica equivalente do isolamento MLI | $4{,}8$ | $\text{K/W}$ | 20 camadas de Mylar aluminizado a vácuo |
| $\eta_c$ | Rendimento coulométrico de recarga | $0{,}995$ | -- | Fator de reversibilidade eletroquímica |

---

### MÓDULO 5: Algoritmo de Estimação de Estado de Carga (AUKF)

#### 5.1 Modelos Matemáticos Implementados
* **Espaço de Estados Discretizado:**
  $$\mathbf{x}_k = \begin{bmatrix} SoC_k & V_{p1,k} & V_{p2,k} \end{bmatrix}^T, \quad u_k = I_{L,k}$$
  $$\mathbf{x}_{k} = f(\mathbf{x}_{k-1}, u_{k-1}) + \mathbf{w}_{k-1}$$
  $$y_k = V_{oc}(SoC_k) - u_k R_0(T_{bat,k}) - V_{p1,k} - V_{p2,k} + v_k$$
* **Geração de $2L+1$ Pontos Sigma Determinísticos:**
  $$\boldsymbol{\chi}_{0,k-1} = \hat{\mathbf{x}}_{k-1}$$
  $$\boldsymbol{\chi}_{i,k-1} = \hat{\mathbf{x}}_{k-1} + \left( \sqrt{(L + \gamma) \mathbf{P}_{k-1}} \right)_i, \quad i = 1, \dots, L$$
  $$\boldsymbol{\chi}_{i+L,k-1} = \hat{\mathbf{x}}_{k-1} - \left( \sqrt{(L + \gamma) \mathbf{P}_{k-1}} \right)_i, \quad i = 1, \dots, L$$
* **Fator de Esquecimento Adaptativo e Atualização de Covariância:**
  $$\mathbf{P}_{k|k-1} = \lambda_k \sum_{i=0}^{2L} W_i^{(c)} \left[ \boldsymbol{\chi}_{i, k|k-1} - \hat{\mathbf{x}}_{k|k-1} \right] \left[ \boldsymbol{\chi}_{i, k|k-1} - \hat{\mathbf{x}}_{k|k-1} \right]^T + \mathbf{Q}_k$$
  $$\lambda_k = \max\left( 1, \, \frac{\text{tr}\left( \mathbf{N}_k \right)}{\text{tr}\left( \mathbf{M}_k \right)} \right)$$

#### 5.2 Catálogo de Parâmetros do AUKF
| Símbolo | Descrição do Parâmetro | Valor Nominal | Unidade | Referência / Justificativa |
| :--- | :--- | :---: | :---: | :--- |
| $L$ | Dimensão do vetor de estado | $3$ | -- | $[SoC, V_{p1}, V_{p2}]^T$ |
| $\alpha_{ut}$ | Parâmetro de espalhamento dos pontos sigma | $10^{-3}$ | -- | Controla a dispersão de amostragem |
| $\beta_{ut}$ | Parâmetro de incorporação de distribuição prévia | $2{,}0$ | -- | Valor ótimo para variáveis gaussianas |
| $\kappa_{ut}$ | Parâmetro de ajuste secundário | $0$ | -- | Garante semidefinição positiva |
| $\mathbf{Q}_0$ | Covariância inicial do ruído de processo | $\text{diag}(10^{-7}, 10^{-5}, 10^{-5})$ | -- | Incerteza no modelo do circuito |
| $\mathbf{R}_0$ | Covariância inicial do ruído de medição | $2{,}5 \times 10^{-4}$ | $\text{V}^2$ | Ruído de sensor de tensão ($15\text{ mV}_{rms}$) |
| $\Delta t_{ukf}$ | Período de execução do filtro de Kalman | $100$ | $\text{ms}$ | Frequência de 10 Hz no computador de bordo (OBC) |
| $RMSE_{meta}$ | Erro quadrático médio máximo tolerado de SoC | $< 1{,}2$ | $\%$ | Garantia de precisão sob degraus de corrente |

---

### MÓDULO 6: Telecomunicações e Enlaces Ópticos Inter-Satélite (OISL)

#### 6.1 Modelos Matemáticos Implementados
* **Consumo Elétrico e Dissipação Térmica dos Terminais Laser:**
  $$P_{OISL,total}(t) = N_{term} \cdot \left[ P_{laser,elec}(t) + P_{fsm} + P_{edfa} + P_{proc,dsp} \right]$$
  $$P_{th,OISL}(t) = (1 - \eta_{eo}) \cdot P_{OISL,total}(t)$$
* **Balanço de Potência Global do Barramento DC:**
  $$P_{total}(t) = P_{compute}(t) + P_{OISL,total}(t) + P_{aux}(t) + P_{loss,conv}(t)$$

#### 6.2 Catálogo de Parâmetros do OISL
| Símbolo | Descrição do Parâmetro | Valor Nominal | Unidade | Referência / Justificativa |
| :--- | :--- | :---: | :---: | :--- |
| $N_{term}$ | Número de terminais ópticos por nó do cluster | $4$ | -- | Topologia de rede em malha e anel redundante |
| $\lambda_{opt}$ | Comprimento de onda de emissão óptica | $1550$ | $\text{nm}$ | Padrão COTS espacial (baixa atenuação e alta potência) |
| $R_{data}$ | Taxa líquida de transmissão por enlace | $100$ | $\text{Gbps}$ | Sincronização de tensores de redes neurais |
| $P_{elec,term}$ | Potência elétrica consumida por terminal ativo | $110$ | $\text{W}$ | Diodo laser + EDFA + FSM + DSP de banda base |
| $\eta_{eo}$ | Eficiência de conversão eletro-óptica | $22$ | $\%$ | Rendimento de emissores semicondutores espaciais |
| $P_{th,term}$ | Calor residual injetado no chassi por terminal | $85{,}8$ | $\text{W}$ | $(1 - \eta_{eo}) \times P_{elec} \rightarrow$ onera o radiador |

---

### MÓDULO 7: Mecânica Orbital e Fluxos Radiativos Ambientais (TCS-Orbital)

#### 7.1 Modelos Matemáticos Implementados
* **Geometria Orbital da Trajetória de Referência (SSO Dawn-Dusk):**
  $$h = 700\text{ km}, \quad i = 98{,}2^\circ, \quad \text{LTAN} = 18\text{h:00min}, \quad \beta \approx 90^\circ$$
* **Incidência Radiativa Superficial Total sobre o Painel $k$:**
  $$q_{in,k}(t) = q_{solar,k}(t) + q_{albedo,k}(t) + q_{IR,Earth,k}(t)$$
  * *Radiação Solar Direta:* $q_{solar,k} = S_{AM0} \cdot \cos(\theta_{sun,k}) \cdot \mathbb{I}_{visivel}$
  * *Albedo Terrestre:* $q_{albedo,k} = S_{AM0} \cdot a_E \cdot F_{k-Earth} \cdot \cos(\theta_{sun,Earth})$
  * *Infravermelho Planetário:* $q_{IR,k} = q_{Earth} \cdot F_{k-Earth}$
* **Fator de Forma Satélite-Terra para Placa Plana em Atitude Nadir/Anti-Nadir:**
  $$F_{k-Earth} = \left( \frac{R_E}{R_E + h} \right)^2 \cos(\phi_{nadir})$$

#### 7.2 Catálogo de Parâmetros Ambientais e Órbita
| Símbolo | Descrição do Parâmetro | Valor Nominal | Unidade | Referência / Justificativa |
| :--- | :--- | :---: | :---: | :--- |
| $h_{sso}$ | Altitude de voo na órbita heliossíncrona | $700$ | $\text{km}$ | Compromisso arrasto residual vs. radiação de Van Allen |
| $i_{sso}$ | Inclinação da órbita heliossíncrona | $98{,}2$ | $^\circ$ | Precessão nodal síncrona com o movimento aparente do Sol |
| $T_{orb}$ | Período de revolução orbital | $98{,}7$ | $\text{min}$ | $\approx 14{,}6$ órbitas por dia terrestre |
| $f_{ecl,sso}$ | Fração anual de eclipse na SSO Dawn-Dusk | $< 0{,}8$ | $\%$ | Eclipses ocorrem apenas em breves semanas de solstício |
| $t_{ecl,max}$ | Duração máxima de eclipse sazonal | $18{,}8$ | $\text{min}$ | Barros e Fernandes (2026b) |
| $R_E$ | Raio médio volumétrico da Terra | $6371$ | $\text{km}$ | Constante geodésica WGS-84 |
| $a_E$ | Coeficiente médio de albedo da Terra | $0{,}30$ | -- | Média global anual refletida |
| $q_{Earth}$ | Densidade de fluxo médio de radiação infravermelha | $237$ | $\text{W/m}^2$ | Emissão de corpo cinzento da Terra para o espaço |
| $T_{space}$ | Temperatura da radiação cósmica de fundo | $3{,}0$ | $\text{K}$ | Fundo térmico do espaço profundo |
| $\sigma$ | Constante universal de Stefan-Boltzmann | $5{,}67037 \times 10^{-8}$ | $\text{W/m}^2\text{K}^4$ | Constante fundamental da termodinâmica |

---

### MÓDULO 8: Dissipação Eletrotérmica dos Processadores Neurais de IA

#### 8.1 Modelos Matemáticos Implementados
* **Potência Dissipada por Nó de Processamento Neural (GPU/TPU/ASIC):**
  $$P_{chip}(t) = P_{dyn}(V_{core}, f_{clk}) + P_{leak}(T_j, V_{core})$$
* **Potência Dinâmica de Comutação Lógica:**
  $$P_{dyn} = \alpha_{act} \cdot C_{eff} \cdot V_{core}^2 \cdot f_{clk}$$
* **Potência Estática de Fuga Sublimiar (*Subthreshold Leakage*):**
  $$P_{leak}(T_j) = V_{core} \cdot I_{leak,0} \cdot \left( \frac{T_j}{T_{ref}} \right)^2 \exp\left[ \frac{-q E_g}{2 k_B T_j} \right]$$
* **Instabilidade de Realimentação Positiva (*Thermal Runaway*):**
  $$\frac{dP_{chip}}{dT_j} = \frac{dP_{leak}}{dT_j} > 0$$
  Se $\frac{dP_{leak}}{dT_j} > \frac{1}{R_{th,j-rad}}$, ocorre fuga térmica descontrolada e queima do silício.

#### 8.2 Catálogo de Parâmetros dos Chips de IA (Nó Acelerador Espacial)
| Símbolo | Descrição do Parâmetro | Valor Nominal | Unidade | Referência / Justificativa |
| :--- | :--- | :---: | :---: | :--- |
| $N_{nodes}$ | Número total de aceleradores neurais no data center | $1200$ | -- | Escala de processamento agregada $\approx 600\text{ kW}$ IT |
| $V_{core,nom}$ | Tensão nominal de alimentação do silício | $0{,}85$ | $\text{V}$ | Nós litográficos FinFET/GAA sub-7nm |
| $V_{core,min}, V_{core,max}$ | Faixa operacional dinâmica de tensão (DVFS) | $0{,}65 - 1{,}05$ | $\text{V}$ | Ajuste dinâmico coordenado pelo TEMS |
| $f_{clk,nom}$ | Frequência nominal de clock | $1{,}50$ | $\text{GHz}$ | Operação balanceada desempenho/eficiência |
| $f_{clk,min}, f_{clk,max}$ | Faixa de modulação de frequência de clock | $0{,}60 - 2{,}10$ | $\text{GHz}$ | Margem de redução sob estresse térmico |
| $C_{eff}$ | Capacitância lógica efetiva comutada | $380$ | $\text{nF}$ | Representativo de nós de tensores densos de IA |
| $\alpha_{act,inf}$ | Fator de atividade em inferência neural | $0{,}45$ | -- | Cargas de ativação esparsa de modelos |
| $\alpha_{act,train}$ | Fator de atividade em treinamento denso | $0{,}88$ | -- | Operação contínua de matrizes FP16/BF16 |
| $P_{chip,nom}$ | Potência nominal dissipada por acelerador | $450$ | $\text{W}$ | Consumo por GPU de IA de alta performance |
| $T_{j,max}$ | Limite estrito de temperatura de junção de segurança | $85{,}0$ | $^\circ\text{C}$ | Prevenção de quebra de dielétrico e eletromigração |
| $T_{j,crit}$ | Limite de corte por hardware (*thermal trip*) | $105{,}0$ | $^\circ\text{C}$ | Fusível térmico de emergência |

---

### MÓDULO 9: Rede Térmica Nodal Multinodal de Resfriamento no Vácuo

#### 9.1 Modelos Matemáticos Implementados
* **Rede de Parâmetros Concentrados de 5 Nós Acoplados:**
  1. *Nó 1 (Junção de Silício):*
     $$C_1 \frac{dT_1(t)}{dt} = P_{chip}(t) - \frac{T_1(t) - T_2(t)}{R_{12}}$$
  2. *Nó 2 (TIM e Espalhador de Calor de Cobre/Vapor Chamber):*
     $$C_2 \frac{dT_2(t)}{dt} = \frac{T_1(t) - T_2(t)}{R_{12}} - \frac{T_2(t) - T_3(t)}{R_{23}}$$
  3. *Nó 3 (Placa Fria com Microcanais de Fluido):*
     $$C_3 \frac{dT_3(t)}{dt} = \frac{T_2(t) - T_3(t)}{R_{23}} - \dot{m} c_p \left( T_{out} - T_{in} \right) - \frac{T_3(t) - T_4(t)}{R_{34}}$$
  4. *Nó 4 (Loop Bifásico Bombeado Mecanicamente - MPFL / Tubos de Calor):*
     $$C_4 \frac{dT_4(t)}{dt} = \dot{m} c_p \left( T_{out} - T_{in} \right) - \frac{T_4(t) - T_5(t)}{R_{45}}$$
  5. *Nó 5 (Painel Radiador Externo Desdobrável - Stefan-Boltzmann):*
     $$C_5 \frac{dT_5(t)}{dt} = \frac{T_4(t) - T_5(t)}{R_{45}} + q_{in,ext}(t) A_{rad} - \varepsilon \sigma A_{rad} \left( T_5(t)^4 - T_{space}^4 \right)$$

#### 9.2 Catálogo de Parâmetros Térmicos e Radiadores
| Símbolo | Descrição do Parâmetro | Valor Nominal | Unidade | Referência / Justificativa |
| :--- | :--- | :---: | :---: | :--- |
| $C_1$ | Capacitância térmica do die de silício | $85$ | $\text{J/K}$ | Massa do chip semicondutor $\approx 50\text{ g}$ |
| $C_2$ | Capacitância térmica do espalhador de calor (spreader) | $420$ | $\text{J/K}$ | Cobre de alta condutividade ($c_p = 385\text{ J/kg}\cdot\text{K}$) |
| $C_3$ | Capacitância da placa fria (*cold plate*) | $1100$ | $\text{J/K}$ | Alumínio estrutural com fluido interno |
| $C_4$ | Capacitância do fluido circulante no loop (MPFL) | $3500$ | $\text{J/K}$ | Amônia ou refrigerante dielétrico fluorado |
| $C_5$ | Capacitância térmica do painel radiador composto | $12500$ | $\text{J/K}$ | Painéis sanduíche de fibra de carbono/alumínio |
| $R_{12}$ | Resistência térmica condutiva junção-spreader (TIM) | $0{,}035$ | $\text{K/W}$ | Solda de índio ou metal líquido estável a vácuo |
| $R_{23}$ | Resistência térmica spreader-cold plate | $0{,}025$ | $\text{K/W}$ | Fixação mecânica sob torque e pasta de alto rendimento |
| $R_{45}$ | Resistência de transferência loop-painel radiador | $0{,}018$ | $\text{K/W}$ | Soldagem por difusão de tubos de condensador |
| $A_{rad}$ | Área total efetiva de radiadores externos desdobráveis | $1450$ | $\text{m}^2$ | Painéis bifaciais articulados orientados no vetor anti-solar |
| $\varepsilon$ | Emissividade no infravermelho do radiador | $0{,}91$ | -- | Revestimento de espelho solar óptico (OSR) de quartzo |
| $\alpha_s$ | Absortância solar da superfície externa do radiador | $0{,}08$ | -- | Baixa absorção de radiação incidente |
| $\dot{m}_{fluid}$ | Vazão mássica nominal de bombeamento do fluido | $2{,}8$ | $\text{kg/s}$ | Bomba mecânica centrífuga redundante (200 W de potência) |

---

### MÓDULO 10: Sistema de Gerenciamento Térmico-Energético Ativo (TEMS / NMPC)

#### 10.1 Modelos Matemáticos Implementados
* **Problema de Controle Ótimo em Horizonte Deslizante ($N$ passos):**
  $$\min_{\mathbf{u}_0, \dots, \mathbf{u}_{N-1}} J = \sum_{k=0}^{N-1} \left( w_T \left[ T_j[k] - T_{ref} \right]^2 + w_P \left[ P_{total}[k] - P_{avail}[k] \right]^2 - w_Q \Psi_{ops}[k] + \Phi_{barrier}(T_j[k]) \right)$$
* **Vetor de Variáveis de Controle ($\mathbf{u}[k]$):**
  $$\mathbf{u}[k] = \begin{bmatrix} f_{clk,1}[k] & \dots & f_{clk,M}[k] & V_{core,1}[k] & \dots & V_{core,M}[k] & \dot{m}_{pump}[k] & \theta_{gimbal}[k] \end{bmatrix}^T$$
* **Restrições Físicas Hard e Soft:**
  $$T_j[k] \le T_{j,max} = 85^\circ\text{C}, \quad \forall k$$
  $$P_{total}[k] \le P_{pv,mpp}[k] + P_{bat,disch,max}[k]$$
  $$SoC[k] \ge SoC_{min} = 30\% \quad \text{(Preservação de vida útil)}$$
  $$f_{clk,min} \le f_{clk}[k] \le f_{clk,max}, \quad |\Delta f_{clk}| \le \Delta f_{slew,max}$$

#### 10.2 Catálogo de Parâmetros do Controlador TEMS
| Símbolo | Descrição do Parâmetro | Valor Nominal | Unidade | Referência / Justificativa |
| :--- | :--- | :---: | :---: | :--- |
| $N$ | Horizonte de predição do controlador MPC | $20$ | passos | Cobre a constante de tempo térmica principal ($\approx 60\text{ s}$) |
| $\Delta t_{ctrl}$ | Intervalo de amostragem do controle térmico | $3{,}0$ | $\text{s}$ | Tempo compatível com a dinâmica térmica lenta dos radiadores |
| $w_T$ | Fator de ponderação do erro de temperatura | $10{,}0$ | $\text{K}^{-2}$ | Prioridade máxima para estabilidade térmica |
| $w_P$ | Fator de ponderação do balanço de potência | $5{,}0$ | $\text{kW}^{-2}$ | Prevenção de sobrecorrente no barramento DC |
| $w_Q$ | Fator de ponderação da vazão de computação | $2{,}5$ | $\text{TFLOPS}^{-1}$ | Incentivo econômico à entrega de processamento útil |
| $\rho_{pen}$ | Coeficiente da barreira para $T_j > 80^\circ\text{C}$ | $10^{4}$ | -- | Função penalidade exponencial para evitar violações |

---

### MÓDULO 11: Confiabilidade, Degradação Cíclica e Fadiga Térmica (10 Anos)

#### 11.1 Modelos Matemáticos Implementados
* **Fadiga Térmica Mecânica das Soldas de Chips (Equação de Coffin-Manson Modificada por Engelmaier):**
  $$N_f = \frac{1}{2} \left( \frac{\Delta \gamma}{2 \varepsilon_f'} \right)^{1/c}$$
  $$\Delta \gamma = \frac{L_{die} \cdot \Delta \alpha_{CTE} \cdot \Delta T_j}{2 h_{solder}}$$
  $$c = -0{,}442 - 6 \times 10^{-4} T_{mean} + 1{,}74 \times 10^{-2} \ln\left( 1 + f_{cycle} \right)$$
* **Eficácia Energética Equivalente Orbital ($PUE_{orbital}$):**
  $$PUE_{orbital} = \frac{P_{total,gerada}}{P_{compute,liquida}} = \frac{P_{IT} + P_{cooling,pumps} + P_{OISL} + P_{losses,DC} + P_{battery,eff}}{P_{IT}}$$

#### 11.2 Catálogo de Parâmetros de Confiabilidade
| Símbolo | Descrição do Parâmetro | Valor Nominal | Unidade | Referência / Justificativa |
| :--- | :--- | :---: | :---: | :--- |
| $L_{die}$ | Dimensão linear diagonal do chip de silício | $28$ | $\text{mm}$ | Aceleradores avançados com encapsulamento CoWoS |
| $\Delta \alpha_{CTE}$ | Diferença de expansão térmica silício/substrato | $14 \times 10^{-6}$ | $\text{K}^{-1}$ | $\alpha_{solder} - \alpha_{silicon}$ |
| $h_{solder}$ | Altura das esferas de solda BGA (Ball Grid Array) | $0{,}18$ | $\text{mm}$ | Encapsulamento microeletrônico avançado |
| $2 \varepsilon_f'$ | Coeficiente de ductilidade de fadiga de solda SAC305 | $0{,}65$ | -- | Liga livre de chumbo qualificada para espaço |
| $t_{mission}$ | Vida útil nominal de projeto em órbita | $10$ | $\text{anos}$ | Operação sem manutenção tripulada |
| $PUE_{meta}$ | Meta de eficácia energética do sistema em SSO | $< 1{,}12$ | -- | Rendimento global frente a data centers terrestres |

---

### TABELA SÍNTESE DE ESTADOS E ENTRADAS DO FRAMEWORK INTEGRADO

```
VETOR DE ESTADOS CONTÍNUOS MULTIFÍSICOS (16 ESTADOS):
  x = [ SoC,  V_p1,  V_p2,         -> Dinâmica Eletroquímica do BESS (3 estados)
        T_j,  T_sp,  T_cp,  T_fl,  T_rad, -> Rede Térmica Nodal de Resfriamento (5 estados)
        V_bus, I_bus,              -> Dinâmica do Barramento DC-DC (2 estados)
        q_0, q_1, q_2, q_3,        -> Quatérnions de Atitude do Satélite (4 estados)
        omega_x, omega_y ]         -> Velocidades Angulares do Corpo Rígido (2 estados)

VETOR DE ENTRADAS DE CONTROLE ATIVO (7 ENTRADAS):
  u = [ f_clk (Frequência DVFS), 
        V_core (Tensão de Núcleo), 
        m_pump (Vazão da Bomba MPFL), 
        D_mppt (Razão Cíclica do MPPT), 
        tau_x, tau_y, tau_z (Torques de Rodas de Reação ADCS) ]
```

---
*Este documento constitui a especificação paramétrica oficial e fechada que rege o equacionamento analítico dos Capítulos 3 e 4, alimentando diretamente os modelos de simulação numérica em MATLAB/Python do Capítulo 5 da Dissertação de Mestrado.*
