# MATRIZ CANÔNICA DE MODELAGEM MATEMÁTICA E CATÁLOGO DE PARÂMETROS
## Dissertação de Mestrado — ProPECAut / IFES Campus Serra
**Título da Pesquisa:** Modelagem Computacional de Gerenciamento Térmico-Energético para Centros de Dados Orbitais em Escala Megawatt  
**Autor:** Dafnhê Thor Castro de Martins Barros  
**Orientador:** Prof. Dr. Flávio Barcelos Braz da Silva  
**Data de Consolidação e Sincronização:** Outubro de 2026  

---

### SUMÁRIO DA MATRIZ
1. [MÓDULO 1: Geração Fotovoltaica Espacial AM0 — Modelo de 3 Diodos em Série (EPS-PV)](#módulo-1-geração-fotovoltaica-espacial-am0--modelo-de-3-diodos-em-série-eps-pv)
2. [MÓDULO 2: Conversor Buck Síncrono e Rastreamento de Máxima Potência (MPPT)](#módulo-2-conversor-buck-síncrono-e-rastreamento-de-máxima-potência-mppt)
3. [MÓDULO 3: Conversão e Distribuição de Potência DC-DC de Alta Tensão (HVDC)](#módulo-3-conversão-e-distribuição-de-potência-dc-dc-de-alta-tensão-hvdc)
4. [MÓDULO 4: Subsistema de Armazenamento Eletroquímico de Baterias (BESS)](#módulo-4-subsistema-de-armazenamento-eletroquímico-de-baterias-bess)
5. [MÓDULO 5: Algoritmo de Estimação de Estado de Carga (AUKF Sage–Husa Livre de Jacobianas)](#módulo-5-algoritmo-de-estimação-de-estado-de-carga-aukf-sagehusa-livre-de-jacobianas)
6. [MÓDULO 6: Telecomunicações e Enlaces Ópticos Inter-Satélite (OISL)](#módulo-6-telecomunicações-e-enlaces-ópticos-inter-satélite-oisl)
7. [MÓDULO 7: Mecânica Orbital e Fluxos Radiativos Ambientais (TCS-Orbital)](#módulo-7-mecânica-orbital-e-fluxos-radiativos-ambientais-tcs-orbital)
8. [MÓDULO 8: Dissipação Eletrotérmica dos Processadores Neurais de IA e Estabilidade FinFET](#módulo-8-dissipação-eletrotérmica-dos-processadores-neurais-de-ia-e-estabilidade-finfet)
9. [MÓDULO 9: Rede Térmica Nodal Multinodal de Resfriamento no Vácuo](#módulo-9-rede-térmica-nodal-multinodal-de-resfriamento-no-vácuo)
10. [MÓDULO 10: Sistema de Gerenciamento Térmico-Energético Ativo (TEMS / NMPC)](#módulo-10-sistema-de-gerenciamento-térmico-energético-ativo-tems--nmpc)
11. [MÓDULO 11: Confiabilidade, Degradação Cíclica e Fadiga Térmica (10 Anos)](#módulo-11-confiabilidade-degradação-cíclica-e-fadiga-térmica-10-anos)

---

### MÓDULO 1: Geração Fotovoltaica Espacial AM0 — Modelo de 3 Diodos em Série (EPS-PV)

#### 1.1 Modelos Matemáticos Implementados
* **Modelo Físico de Três Diodos Independentes em Série ($\text{GaInP}/\text{GaAs}/\text{Ge}$):**
  $$V_{cel}(I) = \sum_{i=1}^{3} V_i(I) - I R_s, \quad \text{com} \quad I = I_{ph,i}(T,G) - I_{0,i}(T) \left[ \exp\!\left( \frac{q V_i}{n_i k_B T_{cel}} \right) - 1 \right] - \frac{V_i}{R_{sh,i}}$$
  onde $i \in \{1: \text{Top GaInP}, \, 2: \text{Middle GaAs}, \, 3: \text{Bottom Ge}\}$ representa as três subcélulas monolíticas interligadas por junções túnel (COTAL et al., 2000; SPECTROLAB, s.d.), resolvidas numericamente por Newton-Raphson para cada nível de corrente $I$.
* **Dependência Térmica e de Irradiância AM0 por Subcélula:**
  $$I_{ph,i}(T_{cel}, G) = \left[ I_{ph,i,ref} + K_{I,i}(T_{cel} - T_{ref}) \right] \frac{G}{G_{datasheet}}$$
  $$V_{oc,i}(T_{cel}) = V_{oc,i,ref} + K_{V,i} (T_{cel} - T_{ref}), \qquad \sum_{i=1}^{3} V_{oc,i,ref} = 1{,}42 + 1{,}01 + 0{,}29 = 2{,}720\text{ V}$$
  $$I_{0,i}(T_{cel},G) = \frac{I_{ph,i}(T_{cel},G) - V_{oc,i}(T_{cel})/R_{sh,i}}{\exp\!\left[\frac{q V_{oc,i}(T_{cel})}{n_i k_B T_{cel}}\right] - 1}$$
* **Dimensionamento Físico do Arranjo Fotovoltaico ($N_s = 208$ e $N_{p,no} = 252$):**
  * *Tensão MPP por célula no regime quente ($+65^\circ\text{C}$):*
    $$V_{mp,cell}(65^\circ\text{C}) = V_{mp,ref} + \frac{dV_{mp}}{dT}(65 - 28) = 2{,}406 - 6{,}30\times 10^{-3}\times 37 = 2{,}1729\text{ V}$$
  * *Requisito do Conversor Buck para $D \le 0{,}90$ em regime quente:*
    $$V_{pv,min} \ge \frac{V_{bus}}{D_{max}} = \frac{400\text{ V}}{0{,}90} = 444{,}4\text{ V} \implies N_s \ge \frac{444{,}4\text{ V}}{2{,}1729\text{ V}} = 204{,}5 \implies \mathbf{N_s = 208\text{ células em série}}$$
    garantindo $V_{mp,str}(65^\circ\text{C}) = 208 \times 2{,}1729\text{ V} = 452{,}0\text{ V} > 400\text{ V}$ e $D = 0{,}885 \le 0{,}90$ (eliminando a inviabilidade física de $D = 1{,}11 > 1$ que ocorreria com $N_s = 166$).
  * *Dimensionamento de strings para entrega de $50\text{ kW}$ em EOL (10 anos, SSO 700 km, $f_{EOL} = 0{,}82$):*
    $$I_{mp,cell,hot} = 0{,}5369\text{ A}, \quad P_{cell,hot} = 1{,}1666\text{ W}$$
    $$N_{total} = \left\lceil \frac{50.000\text{ W}}{1{,}1666\text{ W} \times 0{,}82} \right\rceil = 52.268 \implies \mathbf{N_{p,no} = \left\lceil \frac{52.268}{208} \right\rceil = 252\text{ strings em paralelo}}$$
  * *Total por satélite:* $N_{cel,no} = 208 \times 252 = 52.416\text{ células}$ ($A_{ativa,no} = 158{,}19\text{ m}^2$).
  * *Total na constelação (20 nós):* $1.048.320\text{ células}$, $3163{,}8\text{ m}^2$ de área ativa, $P_{BOL,tot} = 1{,}332\text{ MW}$ (catálogo) / $1{,}303\text{ MW}$ (3 diodos), garantindo $> 1{,}0\text{ MW}$ útil em EOL.
* **Degradação por Radiação em 10 Anos (Método da Dose de Deslocamento — DDD):**
  $$P_{mp,total}(t) = P_{mp,0} \left[ 1 - C_d \log_{10}\left( 1 + \frac{\Phi_{eq}(t)}{\Phi_x} \right) \right] \left[ 1 - \Delta P_{env}\left(\frac{t}{10}\right) \right]$$
  As taxas anuais de fluência são enquadradas como o **cenário radiativo nominal de projeto equivalente a elétrons de $1\text{ MeV}$ sob blindagem de coverglass CMG de $150\ \mu\text{m}$ em SSO $700\text{ km}$**. Avaliam-se $C_d = 0{,}085$ (nominal analítico, retenção EOL de $90{,}42\%$ em 700 km) e $C_d = 0{,}12851$ (ajuste pontual conservador ECSS, retenção EOL de $88{,}17\%$ em 700 km).
* **Contingência de Hot Spot no Vácuo (Stefan-Boltzmann):**
  $$C_{th,cel}\frac{dT_{cel}(t)}{dt} = P_{rev} - n_{faces} \varepsilon_{cel} \sigma A_{cel} \left( T_{cel}^4(t) - T_{space}^4 \right)$$
  onde 1 célula sombreada é forçada à polarização reversa pelas demais 207 células da string.

#### 1.2 Catálogo de Parâmetros do Painel Solar (Spectrolab XTJ Prime — 3 Diodos)
| Símbolo | Descrição do Parâmetro | Valor Nominal | Unidade | Classificação e Referência / Justificativa |
| :--- | :--- | :---: | :---: | :--- |
| $G_{AM0}$ | Irradiância solar extraterrestre a 1 UA | $1361{,}0$ | $\text{W/m}^2$ | *(Especificação de catálogo COTS / Referência validada)* ASTM E490-22; cenário operacional AM0 |
| $G_{datasheet}$ | Irradiância de ensaio de catálogo XTJ Prime | $1353{,}0$ | $\text{W/m}^2$ | *(Especificação de catálogo COTS / Referência validada)* Spectrolab XTJ Prime Datasheet ($28^\circ\text{C}$) |
| $A_{cel}$ | Área ativa adotada por célula | $30{,}18$ | $\text{cm}^2$ | **Premissa de geometria:** Célula espacial de grande área |
| $V_{oc,1..3}$ | Tensões de circuito aberto por subcélula ($28^\circ\text{C}$) | $1{,}42 / 1{,}01 / 0{,}29$ | $\text{V}$ | *(Especificação de catálogo COTS / Referência validada)* GaInP / GaAs / Ge; soma $V_{oc} = 2{,}720\text{ V}$ (COTAL et al., 2000) |
| $n_{1..3}$ | Fatores de idealidade por subcélula | $1{,}15 / 1{,}15 / 1{,}75$ | -- | *(Especificação de catálogo COTS / Referência validada)* Faixas físicas $1{,}1\text{--}1{,}2$ (GaInP/GaAs) e $1{,}5\text{--}2{,}0$ (Ge) |
| $R_{sh,1..3}$ | Resistências shunt por subcélula | $2500 / 3500 / 200$ | $\Omega$ | *(Especificação de catálogo COTS / Referência validada)* GaInP ($>2000\,\Omega$), GaAs ($>3000\,\Omega$), Ge ($\approx 200\,\Omega$) |
| $R_s$ | Resistência série total (contatos + junções túnel) | $0{,}20$ | $\Omega$ | *(Especificação de catálogo COTS / Referência validada)* Arquitetura monolítica em série (COTAL et al., 2000) |
| $J_{sc,cel}, J_{mp,cel}$ | Densidades de corrente de catálogo ($1353\text{ W/m}^2$) | $18{,}0 / 17{,}5$ | $\text{mA/cm}^2$ | *(Especificação de catálogo COTS / Referência validada)* Datasheet Spectrolab XTJ Prime ($V_{mp,cat}=2{,}406\text{ V}$) |
| $I_{sc,ref}, I_{mp,ref}$ | Correntes escaladas para $1361\text{ W/m}^2$ ($30{,}18\text{ cm}^2$) | $0{,}5465 / 0{,}5313$ | $\text{A}$ | **(dimensionado analiticamente)** $J \times A_{cel} \times (1361/1353)$ |
| $K_V, K_I$ | Coeficientes térmicos globais de $V_{oc}$ e $I_{sc}$ | $-5{,}60\text{ mV/K}, +0{,}3036\text{ mA/K}$ | -- | *(Especificação de catálogo COTS / Referência validada)* Datasheet XTJ Prime ($dJ_{sc}/dT = +10\,\mu\text{A/cm}^2\text{K}$) |
| $N_{s,no}, N_{p,no}$ | Células em série / strings em paralelo por satélite | $208 \times 252$ | -- | **(dimensionado analiticamente)** $52.416\text{ células/nó}$ ($A_{ativa,no} = 158{,}19\text{ m}^2$) |
| $A_{array,cluster}$ | Área ativa total na constelação (20 nós) | $3163{,}8$ | $\text{m}^2$ | **(dimensionado analiticamente)** $1.048.320\text{ células}$ no cluster de $1{,}0\text{ MW}$ EOL |
| $V_{mp,str}(28^\circ\text{C})$ | Tensão MPP por string a $28^\circ\text{C}$ (3 diodos / catálogo) | $489{,}6 / 500{,}45$ | $\text{V}$ | **(dimensionado analiticamente)** $489{,}6\text{ V}$ ($2{,}354\text{ V/cel}$ com $R_s=0{,}20\,\Omega$) / $500{,}45\text{ V}$ catálogo |
| $P_{BOL,no}, P_{EOL,no}$ | Potência MPP BOL e EOL (10 anos) por nó | $65{,}55\text{--}66{,}61 / \ge 50{,}0$ | $\text{kW}$ | **(dimensionado analiticamente)** $65{,}55\text{ kW}$ (3 diodos) / $66{,}61\text{ kW}$ (catálogo); $\ge 50{,}0\text{ kW}$ garantido em EOL |
| $V_{mp,str}(T)$ | Faixa térmica de tensão MPP ($65^\circ\text{C}$ a $-100^\circ\text{C}$) | $452{,}0\text{ a }657{,}2$ | $\text{V}$ | **(dimensionado analiticamente)** Modelo de 3 diodos ($425{,}8\text{ V}$ em WCH a $85^\circ\text{C}$; $V_{oc,-100}=714{,}9\text{ V} \le 900\text{ V}$) |
| $\Phi_{anual,700}$ | Taxa de fluência anual equivalente (1 MeV, CMG $150\,\mu\text{m}$) | $2{,}5 \times 10^{13}$ | $\text{e}^-/\text{cm}^2/\text{ano}$ | **[CENÁRIO NOMINAL]** SSO 700 km ($\Phi_x = 1{,}0\times 10^{14}\text{ e}^-/\text{cm}^2$) |
| $C_d$ | Coeficiente de degradação logarítmica DDD | $0{,}085 / 0{,}12851$ | -- | **(dimensionado analiticamente)** Nominal ($90{,}42\%$ EOL) / Conservador ECSS ($88{,}17\%$ EOL) |
| $\Delta P_{env}$ | Degradação óptica ambiental acumulada em 10 anos | $5{,}2$ | $\%$ | **(premissa de modelagem)** Escurecimento UV e micrometeoroides |
| $T_{eq,hotspot}$ | Temperatura de equilíbrio em hot spot ($V_{rev}=-28\text{ V}$) | $201{,}7\text{ (2f)} / 279{,}4\text{ (1f)}$ | $^\circ\text{C}$ | **(dimensionado analiticamente)** Sem bypass ($14{,}88\text{ W}$) vs. $60{,}0^\circ\text{C}$ com Schottky ($0{,}186\text{ W}$) |

---

### MÓDULO 2: Conversor Buck Síncrono e Rastreamento de Máxima Potência (MPPT)

#### 2.1 Modelos Matemáticos Implementados
* **Dinâmica Chaveada do Conversor Buck Síncrono Acoplado ao Arranjo de 3 Diodos:**
  $$C_{in}\frac{dV_{pv}(t)}{dt} = I_{pv,3D}(V_{pv}, T_{cel}) - D(t) I_L(t), \qquad L\frac{dI_L(t)}{dt} = D(t) V_{pv}(t) - V_{bus} - R_L I_L(t)$$
* **Algoritmo Perturba e Observa com Razão Cíclica Adaptativa (P&O-AP):**
  $$D[k] = D[k-1] - M \cdot \left( \frac{\Delta P_{pv}[k]}{\Delta V_{pv}[k]} \right), \qquad |\Delta D| \le \Delta D_{max} = 0{,}002$$
  onde o sinal negativo assegura realimentação estável de malha fechada frente à característica estritamente redutora do Buck ($\partial V_{pv}/\partial D < 0$).

#### 2.2 Catálogo de Parâmetros do Estágio Buck / MPPT
| Símbolo | Descrição do Parâmetro | Valor Nominal | Unidade | Classificação e Referência / Justificativa |
| :--- | :--- | :---: | :---: | :--- |
| $f_{sw,buck}$ | Frequência de chaveamento PWM do conversor Buck | $150{,}0$ | $\text{kHz}$ | **[PROJETO]** Semicondutores GaN com bloqueio $\ge 650\text{--}900\text{ V}$ |
| $L_{buck}, C_{in}$ | Filtro indutivo-capacitivo do estágio Buck | $45{,}0\,\mu\text{H}, \, 250{,}0\,\mu\text{F}$ | -- | **[PROJETO]** $R_L = 25{,}0\text{ m}\Omega$; ripple indutivo $\Delta i_{L,pp} \approx 4{,}72\text{ A}$ a $+65^\circ\text{C}$ |
| $f_{mppt}$ | Frequência de atualização da malha MPPT | $10{,}0$ | $\text{kHz}$ | **[PROJETO]** Sincronizada a cada 15 ciclos PWM de $150\text{ kHz}$ |
| $M, \Delta D_{max}$ | Ganho adaptativo e limitador de passo P&O-AP | $1{,}2\times 10^{-4}\text{ V/W}, \, 0{,}002$ | -- | **[PROJETO]** Janela de razão cíclica $D \in [0{,}55, \, 0{,}90]$ sob $V_{bus}=400\text{ V}$ |
| $\eta_{MPPT}$ | Eficiência média de rastreamento dinâmico no MPP | $99{,}92$ | $\%$ | **(dimensionado analiticamente)** Rastreamento dinâmico contínuo sob rampa térmica pós-eclipse ($15\text{ K/min}$) |

---

### MÓDULO 3: Conversão e Distribuição de Potência DC-DC de Alta Tensão (HVDC)

#### 3.1 Modelos Matemáticos Implementados
* **Arquitetura Modular Redundante $N+1$:** 3 módulos conversores Full-Bridge GaN de $20\text{ kW}$ em paralelo compartilhando os $40\text{ kW}$ de TI ($13{,}33\text{ kW}$ ou $66{,}7\%$ por módulo em regime nominal).
* **Balanço Analítico de Perdas Eletrotérmicas (GaN GS66516T + Ferrite 3C95):**
  $$P_{loss,total} = P_{cond,pri}(T_j) + P_{sw,total} + P_{cond,sec}(T_j) + P_{core,MSE} + P_{cu}$$
  * *Resistência de Canal GaN:* $R_{on}(T_j) = R_{on,25} \left[ 1 + \alpha_{sc} (T_j - 25) \right]$
  * *Critério de Estabilidade contra Fuga Térmica de GaN:* $\frac{\partial P_{loss,sw}}{\partial T_j} = I_{rms,sw}^2 R_{on,25} \alpha_{sc} < \frac{1}{R_{th,total}}$
  * *Perdas Magnéticas (Steinmetz Modificada — MSE):* $P_{core} = k_{fe} f_{sw}^\alpha \Delta B^\beta V_e$
* **Diagnóstico de Impedância de Middlebrook sob Cargas CPL (Formulações Retificadas):**
  $$r_{in,CPL} = -\frac{V_{bus}^2}{P_{CPL}}, \qquad Z_{out}(s) = \frac{(s L_{dist} + R_{dist})(s C_{bus} ESR + 1)}{s^2 L_{dist} C_{bus} + s C_{bus}(R_{dist} + ESR) + 1}$$
  $$C_{bus,min} \ge \max\left( \frac{4 L_{dist}}{R_{dist}^2}, \, \frac{L_{dist} \Delta I^2}{V_{nom}^2 - V_{min}^2}, \, \frac{L_{dist} P_{CPL}^2 \cdot GM^2}{V_{bus}^4} \right)$$
  com coeficiente de amortecimento $\zeta = \frac{R_{dist} + ESR}{2}\sqrt{\frac{C_{bus}}{L_{dist}}}$ e margem de ganho $GM = 10$.
* **Proteção SSPC e Grampeamento Ativo TVS/Snubber:**
  Extinção ultrarrápida ($I_{cut} = 120\text{ A}$ em $t_{off} = 200\text{ ns}$ sob $L_{dist} = 1{,}8\ \mu\text{H}$) gerando sobretensão $V_{spike} = L_{dist} \frac{di}{dt} = 1080\text{ V}$, totalizando $V_{max} = 1480\text{ V}$. Impõe rede de supressão ativa TVS/snubber para absorver $\frac{1}{2} L_{dist} I_{cut}^2 = 12{,}96\text{ mJ}$ e proteger chaves GaN ($650\text{--}900\text{ V}$).

#### 3.2 Catálogo de Parâmetros dos Conversores DC-DC e Proteção SSPC
| Símbolo | Descrição do Parâmetro | Valor Nominal | Unidade | Classificação e Referência / Justificativa |
| :--- | :--- | :---: | :---: | :--- |
| $V_{bus,prim}$ | Tensão nominal do barramento primário HVDC | $400{,}0$ | $\text{V}$ | **(dimensionado analiticamente)** Faixa $324{,}0\text{ V}$ (UVLO) a $442{,}8\text{ V}$ (OV) (BEVER et al., 2006; FERGUSON, 2006) |
| $V_{bus,sec}$ | Tensão do barramento intermediário de rack | $48{,}0$ | $\text{V}$ | *(Especificação de catálogo COTS / Referência validada)* Padrão Open Compute / conversores PoL sub-1 V |
| $f_{sw}$ | Frequência de chaveamento do conversor isolado | $150{,}0$ | $\text{kHz}$ | **[PROJETO]** Ponte completa GaN + transformador planar 3C95 |
| $R_{on,25}, \alpha_{sc}$ | Resistência a $25^\circ\text{C}$ e coef. térmico (GS66516T) | $25{,}0\text{ m}\Omega, \, 0{,}0128\text{ K}^{-1}$ | -- | *(Especificação de catálogo COTS / Referência validada)* Datasheet Infineon GS66516T ($R_{on,150}=65\text{ m}\Omega$; $42{,}60\text{ m}\Omega$ a $80^\circ\text{C}$) |
| $C_{oss}, Q_g$ | Capacitância de saída e carga de porta ($Q_{rr}=0$) | $207\text{ pF}, \, 14{,}2\text{ nC}$ | -- | *(Especificação de catálogo COTS / Referência validada)* Datasheet Infineon GS66516T ($V_{gate}=6{,}0\text{ V}$) |
| $t_r, t_f$ | Tempos de subida e descida de comutação | $12, \, 15$ | $\text{ns}$ | **[PREMISSA DE GATE]** Encapsulamento GaNPX top-side cooled |
| $R_{th,total}$ | Resistência térmica junção-dissipador (GaNPX + AlN/PGS) | $1{,}12$ | $\text{K/W}$ | *(Especificação de catálogo COTS / Referência validada)* $R_{th,jc}=0{,}27\text{ K/W} + R_{th,TIM}=0{,}85\text{ K/W}$ |
| $R_{th,crit,1}, R_{th,crit,4}$ | Resistência térmica crítica de fuga térmica (1 ch / 4 ch) | $2{,}65 / 0{,}662$ | $\text{K/W}$ | **(dimensionado analiticamente)** $\partial P/\partial T_j = 0{,}3775\text{ W/K}$ por chave a $I_{rms}=34{,}35\text{ A}$ |
| $P_{core}, P_{cu}$ | Perdas no núcleo (Ferroxcube U93/152/30 3C95) e cobre | $95{,}31 / 61{,}10$ | $\text{W}$ | **(dimensionado analiticamente)** $V_e=450\text{ cm}^3, \Delta B=0{,}15\text{ T}, k_{fe}=0{,}70, \alpha=1{,}46, \beta=2{,}52$ |
| $P_{loss,total}, \eta_{20kW}$ | Dissipação total e rendimento em plena carga ($20\text{ kW}$) | $617{,}85\text{ W} / 97{,}003\%$ | -- | **(dimensionado analiticamente)** $201{,}01\text{ W}$ pri + $90{,}99\text{ W}$ sw + $169{,}44\text{ W}$ sec + $156{,}41\text{ W}$ mag |
| $|Z_{out,400}|_{max}$ | Pico de impedância passiva no barramento de $400\text{ V}$ | $0{,}0732$ | $\Omega$ | **(dimensionado analiticamente)** $C_{bus}=4{,}7\text{ mF}, L_{dist}=1{,}8\,\mu\text{H}$ (limite Middlebrook $0{,}800\,\Omega$) |
| $|Z_{out,48}|_{max}$ | Pico de impedância no barramento laminado de $48\text{ V}$ | $7{,}23\text{ (BOL)} / 12{,}14\text{ (EOL)}$ | $\text{m}\Omega$ | **[DIAGNÓSTICO]** $L_{48}=0{,}15\,\mu\text{H}$; EOL ($2{,}5\text{ mF}, 9{,}5\text{ m}\Omega$) viola $11{,}52\text{ m}\Omega$ (DC-F03) |
| $t_{trip,crit}, \tau_{ext,min}$ | Retardo adiabático ($1500\text{ A}^2\text{s}$) e extinção suave SSPC | $74{,}8\,\mu\text{s}, \, 31{,}1\,\mu\text{s}$ | -- | **(dimensionado analiticamente)** Proteção contra fusão SiC e ruptura dielétrica ($V_{ds} \le 1200\text{ V}$) |
| $V_{spike,GaN}, E_{L}$ | Sobretensão indutiva ultrarrápida GaN ($200\text{ ns}$) e energia | $1080\text{ V}, \, 12{,}96\text{ mJ}$ | -- | **(dimensionado analiticamente)** Exige grampeamento ativo TVS para conter $V_{max}=1480\text{ V}$ |

---

### MÓDULO 4: Subsistema de Armazenamento Eletroquímico de Baterias (BESS)

#### 4.1 Modelos Matemáticos Implementados
* **Circuito Equivalente de Thévenin de 2ª Ordem (Dual Polarization — DP Sintético Parametrizado):**
  $$V_{bat}(t) = V_{oc}(SoC, T_{bat}) - I_{bat}(t) R_0(T_{bat}) - V_{p1}(t) - V_{p2}(t)$$
  $$V_{pi,k+1} = a_i V_{pi,k} + R_{pi}(1 - a_i)I_{bat,k}, \qquad a_i = e^{-\Delta t / (R_{pi}C_{pi})}, \quad i \in \{1, 2\}$$
* **Cinética Térmica de Arrhenius da Resistência Interna:**
  $$R_0(T_{bat}) = R_{0,ref} \exp\left[ \frac{E_a}{k_B} \left( \frac{1}{T_{bat}} - \frac{1}{T_{ref}} \right) \right]$$
* **Geração Térmica Irreversível e Rede Térmica de Dois Nós (Bateria + Chassi Radiativo):**
  $$P_{Joule}(t) = I_{bat}^2(t) R_0(T_{bat}) + \frac{V_{p1}^2(t)}{R_{p1}} + \frac{V_{p2}^2(t)}{R_{p2}}$$
  $$C_{bat}\frac{dT_{bat}(t)}{dt} = P_{Joule}(t) + P_{heat}(t) - \frac{T_{bat}(t) - T_{chassi}(t)}{R_{bat,chassi}} - \frac{T_{bat}(t) - T_{amb}}{R_{bat,MLI}}$$
  $$C_{chassi}\frac{dT_{chassi}(t)}{dt} = \frac{T_{bat}(t) - T_{chassi}(t)}{R_{bat,chassi}} - \varepsilon_{rad}\sigma A_{rad,bat}\left(T_{chassi}^4(t) - T_{space}^4\right) + \dot{Q}_{amb,chassi}$$
* **Restrição Normativa de Vida Útil e Ciclagem LEO (ECSS-E-ST-20-32C):**
  Para órbitas LEO sujeitas a $> 5.500$ ciclos anuais de carga e descarga, a profundidade de descarga contínua ($DoD$) de regime nominal deve ser limitada entre $20\%\text{ e }30\%$ ($SoC \ge 70\%\text{--}80\%$) para garantir retenção de capacidade $> 80\%$ ao final de 10 anos de missão.

#### 4.2 Catálogo de Parâmetros da Bateria (Arquitetura A — Conexão Direta 108S1P)
| Símbolo | Descrição do Parâmetro | Valor Nominal | Unidade | Classificação e Referência / Justificativa |
| :--- | :--- | :---: | :---: | :--- |
| $E_{bat,nom}, C_{n}$ | Energia e capacidade nominal por nó satelital | $9{,}0\text{ kWh}, \, 22{,}52\text{ Ah}$ | -- | **(dimensionado analiticamente)** Pack 108S1P (6 submódulos 18S; $60\text{ kg}$ a $150\text{ Wh/kg}$; $180\text{ kWh}$ no cluster) |
| $V_{bat,nom}$ | Tensão nominal e limites operacionais do pack | $399{,}6$ ($324{,}0\text{--}442{,}8$) | $\text{V}$ | **(dimensionado analiticamente)** Conexão direta ao barramento (NASA, 2021; ECSS-E-HB-20-02A) |
| $R_{0,ref}$ | Resistência ôhmica interna a $20^\circ\text{C}$ (BOL / EOL) | $45{,}0 / 85{,}0$ | $\text{m}\Omega$ | **[PARAMETRIZAÇÃO FENOMENOLÓGICA]** Parametrização representativa 108S ($32\text{ m}\Omega$ a $30^\circ\text{C}$ BOL) |
| $R_{p1}, C_{p1}$ | Polarização rápida ($\tau_1 = 57{,}6\text{ s}$) | $18{,}0\text{ m}\Omega, \, 3200\text{ F}$ | -- | **[PARAMETRIZAÇÃO FENOMENOLÓGICA]** Dupla camada eletroquímica ($a_1 = 0{,}998266$ para $\Delta t=0{,}1\text{ s}$) |
| $R_{p2}, C_{p2}$ | Polarização lenta ($\tau_2 = 528{,}0\text{ s}$) | $24{,}0\text{ m}\Omega, \, 22000\text{ F}$ | -- | **[PARAMETRIZAÇÃO FENOMENOLÓGICA]** Difusão em fase sólida ($a_2 = 0{,}999811$ para $\Delta t=0{,}1\text{ s}$) |
| $E_a$ | Energia de ativação eletroquímica de Arrhenius | $0{,}32$ ($5{,}127\times 10^{-20}$) | $\text{eV}$ ($\text{J}$) | **[PARAMETRIZAÇÃO FENOMENOLÓGICA]** Calibração representativa Li-ion NMC |
| $\bar{P}_{Joule,orb}$ | Dissipação Joule média orbital em SSO ($T_{orb}=5926\text{ s}$) | $52{,}82$ | $\text{W}$ | **(dimensionado analiticamente)** $212{,}4\text{ W}$ em eclipse ($1128\text{ s}$) e $15{,}3\text{ W}$ em sol ($4798\text{ s}$) |
| $C_{bat}, C_{chassi}$ | Capacidades térmicas do pack e do chassi | $54.360, \, 18.000\text{--}38.000$ | $\text{J/K}$ | **Valor derivado:** $60\text{ kg}$ NMC ($906\text{ J/kgK}$) + estrutura de alumínio |
| $R_{bat,chassi}, R_{bat,MLI}$ | Resistências térmicas ao chassi e através do MLI | $0{,}25, \, 4{,}80$ | $\text{K/W}$ | **(premissa de modelagem)** 2 nós estabiliza em $28{,}4^\circ\text{C}$; 1 nó isolado em MLI atinge $40{,}82^\circ\text{C}$ |
| $DoD_{LEO,cont}$ | Faixa de profundidade de descarga contínua LEO | $20\%\text{--}30\%$ | -- | **[NORMA ESPACIAL]** ECSS-E-ST-20-32C para retenção $> 80\%$ em 10 anos |

---

### MÓDULO 5: Algoritmo de Estimação de Estado de Carga (AUKF Sage–Husa Livre de Jacobianas)

#### 5.1 Modelos Matemáticos Implementados
* **Espaço de Estados Discretizado ($L = 3$, $2L+1 = 7$ pontos sigma):**
  $$\mathbf{x}_k = \begin{bmatrix} SoC_k & V_{p1,k} & V_{p2,k} \end{bmatrix}^T, \qquad y_k = V_{oc}(SoC_k) - I_{bat,k} R_0(T_{bat,k}) - V_{p1,k} - V_{p2,k} + v_k$$
* **Estimação Adaptativa de Ruído de Sage–Husa Derivative-Free:**
  $$\mathbf{R}_k = (1 - d_k) \mathbf{R}_{k-1} + d_k \left( \mathbf{e}_k \mathbf{e}_k^T - \mathbf{P}_{yy,k|k-1}^0 \right)$$
  $$\mathbf{Q}_k = (1 - d_k)\mathbf{Q}_{k-1} + d_k \left[ \mathbf{K}_k \mathbf{e}_k \mathbf{e}_k^T \mathbf{K}_k^T + \mathbf{P}_{k|k} - \mathbf{P}_{xx,k|k-1}^* \right]$$
  onde a dispersão estatística dos pontos sigma substitui rigorosamente qualquer aproximação de Taylor de 1ª ordem ou matriz Jacobiana $\mathbf{F}_k = \partial f / \partial \mathbf{x}$.
* **Fator de Esquecimento Adaptativo para Garantia de Semidefinição Positiva ($\mathbf{P}_{k|k-1} \succ 0$):**
  $$\mathbf{P}_{k|k-1} = \lambda_{ad,k}\mathbf{P}_{xx,k|k-1}^* + \mathbf{Q}_0, \qquad \lambda_{ad,k} = \max\left(1, \, \frac{\operatorname{tr}(\mathbf{C}_{0,k} - \mathbf{R}_k)}{\operatorname{tr}(\mathbf{P}_{yy,k|k-1}^0)}\right), \qquad d_k = \frac{1 - b}{1 - b^{k+1}}$$
  com $\mathbf{P}_{xx,k|k-1}^* = \sum_{i=0}^{2L} W_i^{(c)} (\boldsymbol{\chi}_{i,k|k-1}^* - \hat{\mathbf{x}}_{k|k-1})(\boldsymbol{\chi}_{i,k|k-1}^* - \hat{\mathbf{x}}_{k|k-1})^T$, fator de memória $b = 0{,}95$, simetrização $\mathbf{P} \leftarrow \frac{1}{2}(\mathbf{P} + \mathbf{P}^T)$ e piso de autovalores para assegurar estabilidade estrita na fatoração de Cholesky.

#### 5.2 Catálogo de Parâmetros do Estimador AUKF
| Símbolo | Descrição do Parâmetro | Valor Nominal | Unidade | Classificação e Referência / Justificativa |
| :--- | :--- | :---: | :---: | :--- |
| $\alpha_{ut}, \beta_{ut}, \kappa_{ut}$ | Parâmetros da Transformada Unscented | $10^{-3}, \, 2{,}0, \, 0$ | -- | *(Especificação de catálogo COTS / Referência validada)* Fator de memória de Sage–Husa $b = 0{,}95$ |
| $\sigma_{V,meas}$ | Cenários de ruído do sensor de tensão avaliados | $15{,}0\text{ e }150{,}0$ | $\text{mV}_{rms}$ | **Ensaio de sensibilidade:** Precisão ($R_0=2{,}25\times 10^{-4}\text{ V}^2$) vs. EMI ($0{,}0225\text{ V}^2$) |
| $e_{SoC,15mV}$ | Erro médio de regime de SoC sob $\sigma_V = 15\text{ mV}_{rms}$ | $0{,}028$ ($\approx 0{,}03$) | $\%$ | **(dimensionado analiticamente)** $\text{RMSE}_V = 1{,}18\text{ mV}$, convergência $< 90\text{ s}$ (\texttt{rng(42)}) |
| $e_{SoC,150mV}$ | Erro médio de regime de SoC sob $\sigma_V = 150\text{ mV}_{rms}$ | $0{,}194$ ($\approx 0{,}19$) | $\%$ | **(dimensionado analiticamente)** $\text{RMSE}_V = 11{,}60\text{ mV}$, atende folgadamente à meta $< 1{,}20\%$ |

---

### MÓDULO 6: Telecomunicações e Enlaces Ópticos Inter-Satélite (OISL)

#### 6.1 Modelos Matemáticos Implementados
* **Consumo Elétrico e Dissipação Térmica dos Terminais Laser:**
  $$P_{OISL,total} = N_{term} P_{elec,term} = 4 \times 110\text{ W} = 440\text{ W}, \qquad Q_{OISL,heat} = (1 - \eta_{eo}) P_{OISL,total} = 343{,}2\text{ W}$$

#### 6.2 Catálogo de Parâmetros do OISL
| Símbolo | Descrição do Parâmetro | Valor Nominal | Unidade | Classificação e Referência / Justificativa |
| :--- | :--- | :---: | :---: | :--- |
| $N_{term}, \lambda_{opt}$ | Número de terminais por satélite e comprimento de onda | $4, \, 1550\text{ nm}$ | -- | **(premissa de projeto)** Interconexão em anel/malha ($100\text{ Gbps}$ a $1\text{ Tbps}$) |
| $P_{elec,term}, \eta_{eo}$ | Potência elétrica por terminal e eficiência eletro-óptica | $110\text{ W}, \, 22\%$ | -- | **(especificação COTS)** Gera $85{,}8\text{ W/terminal}$ ($343{,}2\text{ W}$ totais) de calor interno |

---

### MÓDULO 7: Mecânica Orbital e Fluxos Radiativos Ambientais (TCS-Orbital)

#### 7.1 Modelos Matemáticos Implementados
* **Ocultação Solar Cônica Exata (Umbra e Penumbra Contínua $C^1$) e Fator de Forma Analítico ECSS:**
  $$q_{abs,k}(t) = q_{sol,abs,k}(t) + q_{alb,abs,k}(t) + q_{ir,abs,k}(t)$$
  $$q_{alb,abs} = \alpha_{s,EOL} G_{AM0} a_E F_{rad-E}(h, \phi) K_{alb}, \qquad q_{ir,abs} = \varepsilon_{ir,EOL} q_{Earth} F_{rad-E}(h, \phi)$$
* **Prevenção de Dupla Contagem Radiativa (Norma ECSS-E-ST-31C):**
  Em modelos térmicos nodais com injeção direta de fluxos ambientais absorvidos ($Q_{sol}, Q_{alb}, Q_{ir}$), a dissipação radiativa descarrega unicamente contra a temperatura de fundo cósmico ($T_{space} \approx 2{,}7\text{--}3{,}0\text{ K}$):
  $$\dot{Q}_{rad,i} = \varepsilon_i \sigma A_i (T_i^4 - T_{space}^4) \approx \varepsilon_i \sigma A_i T_i^4$$
  eliminando a sobreposição de calor ambiental que ocorreria se $T_{sink}$ fosse acoplada no termo emissivo.
* **Sensibilidade a Desvios de Atitude do ADCS ($\delta\theta_{att}$):**
  $$q_{sol,abs,pert}(\delta\theta_{att}) = \alpha_{s,EOL} G_{AM0} \sin|\delta\theta_{att}| = 190{,}54 \sin|\delta\theta_{att}|\text{ W/m}^2$$

#### 7.2 Catálogo de Parâmetros Orbitais e Balanço Radiativo EOL (Tabela 10 / Script 4.1)
| Símbolo | Descrição do Parâmetro | Valor Nominal | Unidade | Classificação e Referência / Justificativa |
| :--- | :--- | :---: | :---: | :--- |
| $h_{orb}$ | Altitudes avaliadas (SSO / LEO / GEO) | $700 / 500 / 35.786$ | $\text{km}$ | **[PARÂMETRO ORBITAL]** SSO Dawn-Dusk ($i=98{,}19^\circ$) como projeto nominal |
| $T_{orb}, t_{ecl,max}$ | Período orbital e eclipse máximo em SSO / LEO / GEO | $98{,}8 / 94{,}6 / 1436\text{ min}$ | -- | **(dimensionado analiticamente)** Eclipses máximos: $18{,}8\text{ min}$ (SSO), $35{,}8\text{ min}$ (LEO), $72{,}0\text{ min}$ (GEO) |
| $\alpha_{s,EOL}, \varepsilon_{ir,EOL}$ | Propriedades termo-ópticas EOL do espelho OSR (10 anos) | $0{,}14, \, 0{,}88$ | -- | *(Especificação de catálogo COTS / Referência validada)* Tabela 10 ($\alpha_{s,BOL}=0{,}08$; ECSS-E-HB-31-01 Parte 6; GILMORE, 2002) |
| $F_{rad-E}(90^\circ)$ | Fator de visão analítico placa--Terra ($\phi_{rad}=90^\circ$) | $0{,}23282 / 0{,}27125 / 0{,}00729$ | -- | **(dimensionado analiticamente)** SSO ($700\text{ km}$), LEO ($500\text{ km}$; ou $0{,}4473$ inclinado), GEO ($35.786\text{ km}$) |
| $K_{alb}$ | Fator médio orbital de iluminação de albedo | $0{,}32 / 0{,}45 / 0{,}10$ | -- | **Parâmetro de referência:** Curvas de Lambert NASA TM-4527 para SSO / LEO / GEO |
| $q_{amb,EOL}$ | Fluxo ambiental absorvido em EOL (SSO / LEO / GEO) | $52{,}82 / 63{,}55\text{ (104,79)} / 1{,}56$ | $\text{W/m}^2$ | **(dimensionado analiticamente)** Em SSO: $4{,}26\text{ W/m}^2$ (albedo) $+ 48{,}56\text{ W/m}^2$ (IR) |
| $A_{rad}, A_{req,55}$ | Área instalada por nó vs. área requerida a $T_{rad}=55^\circ\text{C}$ | $88{,}60\text{ vs. }84{,}99\text{ (SSO)}$ | $\text{m}^2$ | **(dimensionado analiticamente)** Folga de $+3{,}61\text{ m}^2$ em SSO ($T_{rad}=51{,}40^\circ\text{C}, T_j=76{,}40^\circ\text{C}$) |
| $|\delta\theta_{thr}|, |\delta\theta_{trip}|$ | Desvios de atitude para início de throttling e trip em SSO | $5{,}87^\circ, \, 15{,}24^\circ$ | graus | **(dimensionado analiticamente)** Throttling a $T_j=80^\circ\text{C}$ ($T_{rad}=55^\circ\text{C}$) e trip a $T_j=85^\circ\text{C}$ ($60^\circ\text{C}$) |

---

### MÓDULO 8: Dissipação Eletrotérmica dos Processadores Neurais de IA e Estabilidade FinFET

#### 8.1 Modelos Matemáticos Implementados
* **Decomposição da Potência Dissipada e Corrente de Fuga Sublimiar:**
  $$P_{chip}(V_{core}, f_{clk}, T_j) = (\alpha_{act} C_{eff}) V_{core}^2 f_{clk} + P_{leak,ref}(V_{core}) \exp\!\left[ \beta_T (T_j - T_{ref}) \right]$$
  (SKADRON et al., 2004; MURALI et al., 2007; LE SUEUR; HEISER, 2010; YUAN et al., 2025).
* **Critério Analítico de Estabilidade Eletrotérmica de Lyapunov:**
  $$\left.\frac{\partial P_{chip}}{\partial T_j}\right|_{T_j^*} = \beta_T P_{leak}(T_j^*, V_{core}) < \left.\frac{d\dot{Q}_{out}}{dT_j}\right|_{T_j^*}$$
* **Critério de Estabilidade Térmica em Semicondutores FinFET sub-7nm ($\theta_{ja} < 0{,}182\text{ K/W}$):**
  $$\frac{\partial P_{leak}}{\partial T_j} < \frac{1}{\theta_{ja}}$$
  Em nível de die individual ($P_0 = 45\text{ W}$ a $25^\circ\text{C}$, $\gamma = 0{,}025\text{ K}^{-1}$): a $T_j = 88{,}4^\circ\text{C}$, $\partial P_{leak}/\partial T_j = 0{,}025 \times 45 \times e^{0{,}025(88{,}4 - 25)} = 5{,}49\text{ W/K}$.
  Para interface convencional degradada ($\theta_{ja} = 0{,}85\text{ K/W}$, condutância $1/\theta_{ja} = 1{,}176\text{ W/K}$), $\partial P_{leak}/\partial T_j > 1/\theta_{ja}$ deflagrando colapso térmico irreversível (*thermal runaway*).
  A estabilidade assintótica estática exige estritamente:
  $$\theta_{ja} < \frac{1}{5{,}49\text{ W/K}} \approx 0{,}182\text{ K/W}$$
  o que comprova a necessidade mandatória de evaporadores de placas frias com microcanais bifásicos de circulação forçada ($R_{12}+R_{23}+R_{34} \le 0{,}075\text{ K/W}$).

#### 8.2 Catálogo de Parâmetros dos Aceleradores Neurais (Nó de 40 kW TI)
| Símbolo | Descrição do Parâmetro | Valor Nominal | Unidade | Classificação e Referência / Justificativa |
| :--- | :--- | :---: | :---: | :--- |
| $P_{IT,no}, Q_{int,no}$ | Potência útil de TI e carga térmica interna total por nó | $40{,}0\text{ kW}, \, 41{,}984\text{ kW}$ | -- | **(dimensionado analiticamente)** 80 dies por satélite ($800\text{ kW}$ TI em 20 satélites operacionais) |
| $V_{core}, f_{clk}$ | Faixas operacionais de tensão e frequência (DVFS) | $0{,}65\text{--}1{,}05\text{ V}, \, 1{,}0\text{--}1{,}8\text{ GHz}$ | -- | **[PROJETO]** Nominal $V_{core}=0{,}85\text{ V}$; inferência ($171\text{ nF}$) e treino ($334{,}4\text{ nF}$) |
| $\beta_T$ | Coeficiente térmico exponencial de fuga sublimiar | $0{,}0285$ | $\text{K}^{-1}$ | **(premissa de projeto)** A $70^\circ\text{C}$, $\partial P_{leak}/\partial T_j \approx 102{,}6\text{ W/K} \ll 860{,}2\text{ W/K}$ |
| $T_{target}, T_{thr}, T_{crit}$ | Alvo nominal, limiar de throttling e corte de emergência | $75{,}0 / <80{,}0 / 85{,}0$ | $^\circ\text{C}$ | **[REQUISITO DE PROJETO]** Prevenção de runaway e fadiga termomecânica |
| $\theta_{ja,crit}$ | Resistência térmica máxima para estabilidade FinFET a $88{,}4^\circ\text{C}$ | $0{,}182$ | $\text{K/W}$ | **(dimensionado analiticamente)** $\theta_{ja} < 1 / (5{,}49\text{ W/K})$; exige placa fria microcanalada |

---

### MÓDULO 9: Rede Térmica Nodal Multinodal de Resfriamento no Vácuo

#### 9.1 Modelos Matemáticos Implementados
* **Rede de Parâmetros Concentrados de 5 Nós Acoplados ($\mathbf{T} = [T_1, T_2, T_3, T_4, T_5]^T$):**
  $$\mathbf{C} \frac{d\mathbf{T}(t)}{dt} = \mathbf{Q}_{in}(t) - \mathbf{K}_{cond} \mathbf{T}(t) - \mathbf{Q}_{rad}(\mathbf{T}(t))$$
* **Aleta Compósita Anisotrópica de Grafite Pirolítico Encapsulado (APG + K13D2U):**
  $$k_{eff} = \frac{k_{APG} t_{APG} + 2 k_c t_c}{t_{APG} + 2 t_c} = 1275{,}0\text{ W/(m}\cdot\text{K)}, \qquad \eta_{fin,ideal} = \frac{\tanh(m L_{fin})}{m L_{fin}} \approx 0{,}9978 \quad (\eta_{fin,adotado} = 0{,}945)$$

#### 9.2 Catálogo de Parâmetros Térmicos da Rede de 5 Nós
| Símbolo | Descrição do Parâmetro | Valor Nominal | Unidade | Classificação e Referência / Justificativa |
| :--- | :--- | :---: | :---: | :--- |
| $C_1 \dots C_5$ | Capacitâncias térmicas dos nós 1 a 5 | $85 / 420 / 1100 / 3500 / 12500$ | $\text{J/K}$ | **(premissa de projeto)** Die, spreader, cold plate, loop LHP/MPFL e painel radiador |
| $R_{12}, R_{23}, R_{34}, R_{45}$ | Resistências térmicas inter-nodais | $0{,}035 / 0{,}025 / 0{,}015 / 0{,}018$ | $\text{K/W}$ | **(premissa de projeto)** TIM1 (GaInSn), TIM2, evaporador-loop e condensador-radiador |
| $k_{eff}, \eta_{fin}$ | Condutividade planar APG/K13D2U e eficiência efetiva | $1275\text{ W/(mK)}, \, 0{,}945$ | -- | **[VERIFICÁVEL / PROVISÓRIO]** Inclui resistência de colagem tubo-aleta |

---

### MÓDULO 10: Sistema de Gerenciamento Térmico-Energético Ativo (TEMS / NMPC)

#### 10.1 Modelos Matemáticos Implementados
* **Camada Rápida Descentralizada ($\Delta t_{fast} \le 1{,}0\text{ ms}$):** DVFS com banda de histerese anti-chattering ($\Delta T_{hyst} = 3{,}0^\circ\text{C}$) e limitador de taxa ($\Delta f_{slew,max} = 100\text{ MHz/s}$).
* **Camada Supervisória Lenta NMPC ($\Delta t_{ctrl} = 3{,}0\text{ s}$, $N_p = 20$, $N_u = 10$):**
  $$\min_{\mathbf{u}} J = \sum_{i=0}^{N_p-1} \left( w_T [T_j - T_{target}]^2 + w_P [P_{total} - P_{avail}]^2 - w_\Psi \Psi_{ops} + \Phi_{barrier}(T_j) \right) + \sum_{i=0}^{N_u-1} \|\Delta \mathbf{u}\|_{\mathbf{W}_u}^2$$
  (FERREAU et al., 2014; HEGRENÆS; GRAVDAHL; TØNDEL, 2014; BARROS; FERNANDES, 2026a).

---

### MÓDULO 11: Confiabilidade, Degradação Cíclica e Fadiga Térmica (10 Anos)

#### 11.1 Modelos Matemáticos Implementados
* **Fadiga Termomecânica das Soldas BGA SAC305 (Coffin-Manson-Engelmaier com Underfill):**
  $$\Delta T_{solder} = k_{att} \Delta T_j + (1 - k_{att}) \Delta T_{cp}, \qquad \Delta \gamma_{eff} = F \frac{L_{DNP} \Delta \alpha_{CTE} \Delta T_{solder}}{h_{solder}}$$
  $$N_f = \frac{1}{2} \left( \frac{\Delta \gamma_{eff}}{2 \varepsilon_f'} \right)^{1/c}, \qquad c = -0{,}470 - 7{,}0 \times 10^{-4} T_{mean,solder} + 2{,}0 \times 10^{-2} \ln(1 + f_{cpd})$$
* **Eficácia Energética Orbital ($PUE_{orbital}$) e Rendimento Computacional ($\eta_{comp}$):**
  $$PUE_{orbital} = \frac{P_{IT} + P_{cooling} + P_{elec,loss} + P_{aux}}{P_{IT}}, \qquad \eta_{comp} = \frac{\text{TOPS}_{agregado}}{P_{total,no}}$$

#### 11.2 Catálogo de Parâmetros de Confiabilidade e Desempenho (Script 4.6)
| Símbolo | Descrição do Parâmetro | Valor Nominal | Unidade | Classificação e Referência / Justificativa |
| :--- | :--- | :---: | :---: | :--- |
| $L_{DNP}$ | Distância do centro neutro à esfera de canto ($28\times 28\text{ mm}$) | $19{,}80$ | $\text{mm}$ | **(dimensionado analiticamente)** Semidiagonal exata $\frac{1}{2}\sqrt{28^2 + 28^2} = 19{,}799\text{ mm}$ |
| $h_{solder}, \Delta\alpha_{CTE}$ | Altura da esfera BGA e descasamento de expansão térmica | $0{,}18\text{ mm}, \, 14{,}0\times 10^{-6}\text{ K}^{-1}$ | -- | *(Especificação de catálogo COTS / Referência validada)* IPC-9701A / IPC-SM-785 ($F = 0{,}90$, $k_{att} = 0{,}85$ (premissa de modelagem)) |
| $2\varepsilon_f'$ | Coeficiente de ductilidade equivalente com *underfill* | $2{,}887$ ($0{,}65$ sem underfill) | -- | **Calibração empírica equivalente:** Redistribuição por *underfill* epóxi (OSTERMAN, 2018; SYED, 2004) |
| $N_f, \text{Vida}_{SSO+TEMS}$ | Vida em fadiga no cenário nominal (`SSO + TEMS`, $\Delta T_j \le 7^\circ\text{C}$) | $73.063\text{ ciclos} / 13{,}72\text{ anos}$ | -- | **(dimensionado analiticamente)** $L_{DNP}=19{,}80\text{ mm}$; margem $+37{,}2\%$ sobre 10 anos (LEO sem TEMS: $3.165\text{ ciclos} / 6{,}8\text{ meses}$) |
| $PUE_{orbital}$ | Eficácia energética orbital (Nominal / Pior Caso WCH) | $1{,}050 / 1{,}144$ | -- | **(dimensionado analiticamente)** $41.983\text{ W}/40.000\text{ W}$ (nom.) e $45.757\text{ W}/40.000\text{ W}$ (WCH); $\eta_{comp}=0{,}476\text{ TOPS/W}$ |

---
*Este documento constitui a especificação paramétrica oficial, fechada e sincronizada que rege o equacionamento analítico dos Capítulos 3 e 4, dos Apêndices e das rotinas de simulação numérica da Dissertação de Mestrado.*
