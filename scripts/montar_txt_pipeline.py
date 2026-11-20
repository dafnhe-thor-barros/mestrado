# -*- coding: utf-8 -*-
import os
import datetime

out_path_1 = r"c:\ZZ - Antigravity\Mestrado\projeto-mestrado\simulacoes\PIPELINE_SCRIPTS_SIMULACAO_ENCADEADOS.txt"
out_path_2 = r"c:\ZZ - Antigravity\Mestrado\01 - Artigos\Artigo 02\scripts\PIPELINE_SCRIPTS_SIMULACAO_ENCADEADOS.txt"

scripts = [
    {
        "id": "1",
        "name": "simulacao_4_1_fluxos_orbitais.py",
        "path": r"c:\ZZ - Antigravity\Mestrado\projeto-mestrado\simulacoes\simulacao_4_1_fluxos_orbitais.py",
        "category": "CAMADA 1: AMBIENTE ESPACIAL E TERMO-ORBITAL",
        "subsystem": "Fluxos Orbitais, Irradiância AM0, Albedo, Radiação Infravermelha e Balanço Térmico",
        "inputs": "Altitude orbital (h = 600 km, SSO 97,8°), constante solar AM0 (S = 1361 W/m²), coeficiente de albedo terrestre (Ka = 0,30), irradiância infravermelha terrestre (q_IR = 237 W/m²), emissividade e absortividade dos materiais, massas concentradas.",
        "outputs": "Históricos temporais de irradiância incidente total e por componente, temperatura de equilíbrio orbital, duração precisa das fases de insolação (61,5 min) e eclipse (35,2 min) ao longo do período orbital de 96,7 min.",
        "downstream": "Define o perfil térmico e de iluminância que alimenta o modelo do gerador fotovoltaico (simulacao_3_2), o controle térmico de carga (simulacao_4_5) e a co-simulação multifísica (simula_eps_completo)."
    },
    {
        "id": "2",
        "name": "simulacao_3_2_gerador_fotovoltaico_am0.py",
        "path": r"c:\ZZ - Antigravity\Mestrado\projeto-mestrado\simulacoes\simulacao_3_2_gerador_fotovoltaico_am0.py",
        "category": "CAMADA 2: GERAÇÃO FOTOVOLTAICA ESPACIAL",
        "subsystem": "Gerador Fotovoltaico Multijunção AM0 (Spectrolab XTJ Prime), 3 Diodos e Degradação DDD",
        "inputs": "Irradiância solar AM0 incidente e perfil de temperatura Tpv(t) (oriundos de 4.1), parâmetros de 3 diodos em série (subcélulas GaInP, GaAs, Ge com modelo de Varshni para bandgap), taxa de dano acumulado por radiação espacial (DDD - Displacement Damage Dose equivalente a 1 MeV elétrons).",
        "outputs": "Curvas características I-V e P-V em início de vida (BOL) e fim de vida (EOL), tensão de circuito aberto (Voc), corrente de curto-circuito (Isc), ponto de máxima potência (MPP), fator de forma (FF), eficiência global de conversão e área líquida calibrada (4008 m²).",
        "downstream": "Fornece a curva estática e dinâmica de geração e o ponto de operação para os conversores MPPT e para a injeção de potência no barramento (simulacao_3_1 e simula_eps_completo)."
    },
    {
        "id": "3",
        "name": "simulacao_3_4_conversores_gan_middlebrook.py",
        "path": r"c:\ZZ - Antigravity\Mestrado\projeto-mestrado\simulacoes\simulacao_3_4_conversores_gan_middlebrook.py",
        "category": "CAMADA 3: ELETRÔNICA DE POTÊNCIA E ESTABILIDADE DE IMPEDÂNCIA",
        "subsystem": "Conversores CC-CC GaN de Alta Frequência, Balanço de Perdas e Critério de Middlebrook",
        "inputs": "Tensão nominal do barramento (400 V), potência total de 1 MW fracionada em topologia modular paralela (módulos DAB e Buck multifásico), frequência de chaveamento (100 kHz), parâmetros de comutação e condução de transistores GaN (Rds_on, Qg, Coss), impedância negativa incremental de cargas CPL (Z_cpl = -Vbus² / P).",
        "outputs": "Mapa de perdas (condução, comutação e magnéticos), curvas de rendimento elétrico global, matriz de impedância de pequenos sinais de saída Zout(s), e síntese analítica da rede amortecedora passiva paralela Rd - Cd que assegura margem de ganho de 6 dB pelo critério de Middlebrook.",
        "downstream": "Fornece a parametrização do amortecedor passivo e as restrições de impedância para a estabilidade do barramento em simulacao_3_1 e simula_eps_completo."
    },
    {
        "id": "4",
        "name": "simulacao_3_3_bateria_dp_aukf.py",
        "path": r"c:\ZZ - Antigravity\Mestrado\projeto-mestrado\simulacoes\simulacao_3_3_bateria_dp_aukf.py",
        "category": "CAMADA 4: ARMAZENAMENTO ELETROQUÍMICO E ESTIMATIVA DE ESTADO",
        "subsystem": "Banco de Baterias (BESS) Eletrotérmico Dual Polarization (2-RC) e Filtragem Adaptativa AUKF",
        "inputs": "Perfil de corrente de carga e descarga orbital (descarga contínua de 1 MW durante os 35,2 min de eclipse; recarga controlada durante os 61,5 min de sol), parâmetros do circuito equivalente Dual Polarization (OCV(SoC), resistência ôhmica R0, redes RC de transferência de carga R1-C1 e difusão R2-C2), capacidade nominal do pack (38,88 kAh) e modelo térmico acoplado (Rth = 0,14 K/W).",
        "outputs": "Perfis contínuos de Estado de Carga (SoC), tensão terminal Vbatt(t), temperatura interna do pack, estimativa em tempo real de resistência ôhmica interna R0 e validação de convergência do filtro AUKF.",
        "downstream": "Define as características elétricas e limites de operação das baterias que sustentam o barramento nos eclipses (simulacao_3_1 e simula_eps_completo)."
    },
    {
        "id": "5",
        "name": "simulacao_3_1_barramento_hvdc.py",
        "path": r"c:\ZZ - Antigravity\Mestrado\projeto-mestrado\simulacoes\simulacao_3_1_barramento_hvdc.py",
        "category": "CAMADA 5: DINÂMICA DO BARRAMENTO E PROTEÇÃO ATIVA",
        "subsystem": "Barramento HVDC 400 V, Balanço Multi-Órbita, Curto-Circuito e Chave de Estado Sólido SSPC",
        "inputs": "Injeção de corrente dos conversores solares (3.2 / 3.4), suporte bidirecional do BESS (3.3), carga pulsada e contínua de 1 MW a 2,2 MW pico com dinâmica CPL, parâmetros parasitas da linha de 50 m (indutância e capacitância de barramentos planares laminares).",
        "outputs": "Dinâmica da tensão Vbus(t), resposta a degraus violentos de carga, análise de sobretensão/subtensão transiente, caracterização de curto-circuito e perfil de disparo de proteção ultra-rápida por chave semicondutora SSPC.",
        "downstream": "Valida o fechamento de estabilidade global em malha fechada e parametriza os testes de estresse em simula_eps_completo."
    },
    {
        "id": "6",
        "name": "calcular_tempo_critico_sspc.py",
        "path": r"c:\ZZ - Antigravity\Mestrado\projeto-mestrado\simulacoes\calcular_tempo_critico_sspc.py",
        "category": "CAMADA 5: DINÂMICA DO BARRAMENTO E PROTEÇÃO ATIVA (ANEXO)",
        "subsystem": "Tempo Crítico de Desarme da Proteção Ativa SSPC sob Cargas CPL",
        "inputs": "Capacitância efetiva total do barramento (incluindo banco de hold-up), indutância do barramento laminar, resistência negativa da carga no ponto de operação e taxa dv/dt durante curto-circuito.",
        "outputs": "Tempo máximo admissível de detecção e desarme da chave SSPC (tempo crítico em microssegundos) antes que o barramento caia irreversivelmente abaixo de 250 V.",
        "downstream": "Especifica o tempo de resposta máximo para os controladores digitais do barramento em simulacao_3_1 e no artigo."
    },
    {
        "id": "7",
        "name": "simulacao_4_5_tems_nmpc.py",
        "path": r"c:\ZZ - Antigravity\Mestrado\projeto-mestrado\simulacoes\simulacao_4_5_tems_nmpc.py",
        "category": "CAMADA 6: CARGA COMPUTACIONAL E GERENCIAMENTO TÉRMICO PREDITIVO",
        "subsystem": "Controle Preditivo Não Linear de Temperatura (TEMS-NMPC) de Cargas de 50 kW",
        "inputs": "Perfil orbital periódico de temperatura ambiente e dissipação radiativa (4.1), capacidade térmica concentrada dos servidores de processamento de IA, limite térmico de segurança de junção de silício (Tj,max = 85 °C), horizonte de predição do NMPC.",
        "outputs": "Trajetória de controle preditivo combinando modulação de vazão de bomba de fluido refrigerante com escalonamento dinâmico de frequência e tensão (DVFS), pré-resfriamento térmico orbital antecipado e achatamento de picos térmicos em comparação a controladores histeréticos/reativos.",
        "downstream": "Fornece os ciclos térmicos reais de operação dos chips para alimentar a modelagem de degradação termomecânica (simulacao_4_6)."
    },
    {
        "id": "8",
        "name": "simulacao_4_6_fadiga_bga_cme.py",
        "path": r"c:\ZZ - Antigravity\Mestrado\projeto-mestrado\simulacoes\simulacao_4_6_fadiga_bga_cme.py",
        "category": "CAMADA 7: CONFIABILIDADE ESTRUTURAL E EFICIÊNCIA ENERGÉTICA",
        "subsystem": "Fadiga Termomecânica de Juntas de Solda BGA e Métricas de Eficiência Energética (PUE Orbital)",
        "inputs": "Perfis térmicos de ciclagem gerados pelo TEMS-NMPC e controle reativo (4.5), propriedades mecânicas da liga de solda livre de chumbo SAC305 (coeficiente de expansão térmica CTE, módulo de elasticidade, expoentes de fadiga), geometria de esferas BGA e consumo dos auxiliares térmicos.",
        "outputs": "Variação de deformação plástica cisalhante por ciclo (Δγ), previsão de ciclos de vida até a falha (Nf) por modelos de Coffin-Manson e Engelmaier, e cálculo do PUE (Power Usage Effectiveness) anualizado do data center orbital.",
        "downstream": "Consolida os indicadores de durabilidade e viabilidade computacional para a conclusão da dissertação de mestrado."
    },
    {
        "id": "9",
        "name": "simula_eps_completo.py",
        "path": r"C:\ZZ - Antigravity\Mestrado\01 - Artigos\Artigo 02\scripts\simula_eps_completo.py",
        "category": "CAMADA 8: MOTOR DE CO-SIMULAÇÃO MULTIFÍSICA INTEGRADA (ARTIGO 02)",
        "subsystem": "Co-Simulador Multifísico Integrado do Artigo 02 (PV + Térmica + MPPT + Bateria + Barramento)",
        "inputs": "Acoplamento direto e simultâneo em passo de tempo contínuo/discreto dos modelos desenvolvidos nas camadas 1 a 5: solar 3J com deriva térmica, rastreamento MPPT (P&O, IncCond, VP&O e TC-PMPPT), bateria Thevenin 2-RC com estimador SR-UKF adaptativo seletivo em Q(R0) avaliado em 5 sementes, e barramento 400 V com degrau de carga de 1 MW a 2,2 MW, amortecedor de Middlebrook e teto de saturação do conversor de 6 kA.",
        "outputs": "Geração do banco de dados completo de resultados numéricos validados (resultados_eps.json), interface visual interativa (painel_eps.html), e matriz de dados para preenchimento das Tabelas 1 a 8 e Figuras 2 a 5.",
        "downstream": "Alimenta gerar_figuras_pdf.py e o construtor automatizado de manuscrito gerar_artigo_remunom.py."
    },
    {
        "id": "10",
        "name": "gerar_figuras_pdf.py",
        "path": r"C:\ZZ - Antigravity\Mestrado\01 - Artigos\Artigo 02\scripts\gerar_figuras_pdf.py",
        "category": "CAMADA 9: PÓS-PROCESSAMENTO GRÁFICO VETORIAL",
        "subsystem": "Gerador de Figuras Vetoriais em Alta Definição (PDF) para o Artigo 02",
        "inputs": "Séries temporais e dados tabulados calculados pelo motor de co-simulação integrada simula_eps_completo.py.",
        "outputs": "Arquivos vetoriais em formato PDF prontos para submissão: fig_termica.pdf (temperatura e rampa térmica da manta), fig_mppt.pdf (trajetória transitória de rastreamento de MPPT), fig_bateria.pdf (tensão terminal, SoC e estimativa de R0) e fig_barramento.pdf (resposta dinâmica a degrau com e sem amortecimento).",
        "downstream": "Inseridos diretamente no documento LaTeX (main.tex) e convertidos para imagens de alta resolução no Word."
    },
    {
        "id": "11",
        "name": "gerar_artigo_remunom.py",
        "path": r"C:\ZZ - Antigravity\Mestrado\01 - Artigos\Artigo 02\scripts\gerar_artigo_remunom.py",
        "category": "CAMADA 10: CONSTRUÇÃO DOCUMENTAL E SÍNTESE EDITORIAL",
        "subsystem": "Gerador Paramétrico Oficial do Manuscrito Word (.docx Completo, Cego e Revisado)",
        "inputs": "Valores numéricos exatos de todas as simulações, textos em inglês técnico acadêmico, 17 equações matemáticas OMML completas, 8 tabelas de dados formatadas no padrão visual da REMUNOM, e lista canônica das 26 referências bibliográficas no padrão ABNT.",
        "outputs": "Três documentos Word oficiais sincronizados: Artigo_02_REMUNOM_Versao_Completa.docx, Artigo_02_REMUNOM_Versao_Avaliacao_Cega.docx e Artigo_02_REMUNOM_Versao_Revisada_Destacada.docx.",
        "downstream": "Submissão formal na plataforma da Revista Multidisciplinar do Nordeste Mineiro (REMUNOM)."
    }
]

def generate_pipeline_txt():
    global out_path_1, out_path_2, scripts
    lines = []
    
    def bar(c="="):
        return c * 80
        
    lines.append(bar("="))
    lines.append("PIPELINE ESTRUTURADO DE SIMULAÇÃO E MODELAGEM MULTIFÍSICA - EPS 1 MW LEO")
    lines.append("Subsistema Elétrico e Térmico de Potência para Data Centers Orbitais")
    lines.append("Dissertação de Mestrado Profissional (ProPECAut / IFES) & Artigo 02 (REMUNOM)")
    lines.append(bar("="))
    lines.append("Autor: Dafnhe Thor Castro de Martins Barros (IFES / Petrobras)")
    lines.append("Orientador: Prof. Dr. Flavio Barros Braz da Silva (IFES Serra)")
    lines.append(f"Data de Consolidacao: {datetime.datetime.now().strftime('%d de %B de %Y')}")
    lines.append("Ambiente: Python 3.12+ (NumPy, SciPy, Matplotlib, python-docx)")
    lines.append(bar("="))
    lines.append("")
    lines.append(bar("-"))
    lines.append("SUMÁRIO DO DOCUMENTO")
    lines.append(bar("-"))
    lines.append("PARTE 1: ARQUITETURA GERAL DO SISTEMA E FLUXO DE DADOS (DAG)")
    lines.append("PARTE 2: TABELA DE ENCADEAMENTO DOS SUBSISTEMAS E MATRIZ DE INTERDEPENDÊNCIAS")
    lines.append("PARTE 3: SCRIPT MESTRE DE ORQUESTRAÇÃO / EXECUÇÃO DO PIPELINE COMPLETO")
    lines.append("PARTE 4: CÓDIGO-FONTE INTEGRAL E ESTRUTURADO DOS 11 SCRIPTS ENCADEADOS:")
    for s in scripts:
        lines.append(f"   [{s['id'].rjust(2)}] {s['name']:<46} | {s['subsystem'][:45]}...")
    lines.append(bar("="))
    lines.append("")
    lines.append("")

    # PARTE 1
    lines.append(bar("="))
    lines.append("PARTE 1: ARQUITETURA GERAL DO SISTEMA E FLUXO DE DADOS (DAG)")
    lines.append(bar("="))
    lines.append("""
O subsistema de suprimento de energia eletroquimica e fotovoltaica de 1 MW para data
centers orbitais em orbita baixa terrestre (SSO a 600 km) requer o acoplamento de
fenomenos multifisicos operando em multiplas escalas de tempo:
  - Escala lenta (minutos/horas): Ciclos orbitais sol/eclipse (96,7 min), balanco termico
    radiativo da manta solar e do radiador, integracao de estado de carga da bateria (SoC).
  - Escala media (segundos): Resposta transiente termica na saida da penumbra (144,3 K/min),
    deriva termica de tensao (-2,317 V/K), controle de ponto de maxima potencia (MPPT) e
    adaptacao estocastica de ruido do filtro de bateria (SR-UKF).
  - Escala rapida (micro a milissegundos): Dinamica de chaveamento GaN (100 kHz), transitorios
    do barramento de 400 V frente a degraus de carga CPL, atencao a saturacao de corrente
    do conversor (6 kA) e tempo critico de desarme de protecoes SSPC.

DIAGRAMA DO FLUXO DE DADOS E DEPENDÊNCIAS (DAG):

 [ CAMADA 1: simulacao_4_1_fluxos_orbitais.py ]
   │
   ├──────> Irradiancia AM0 (1361 W/m2), Albedo, IR, Perfil Sol (61,5 min) / Eclipse (35,2 min)
   │        │
   │        ▼
   ├─> [ CAMADA 2: simulacao_3_2_gerador_fotovoltaico_am0.py ]
   │     │  Modelo 3 Diodos (Spectrolab XTJ Prime: GaInP/GaAs/Ge) + Degradacao DDD
   │     │  Curvas I-V/P-V, MPP, Voc, Isc, Deriva Termica (-2,317 V/K) e Area (4008 m2)
   │     │  │
   │     │  ▼
   │     ├─> Curva de Potencia Gerada & MPP
   │     │
   ├─> [ CAMADA 3: simulacao_3_4_conversores_gan_middlebrook.py ]
   │     │  Eletronica de Potencia GaN (100 kHz), Topologia Modular DAB + Buck
   │     │  Mapa de Perdas, Impedancia Zout(s) e Amortecedor Passivo Middlebrook (Rd - Cd)
   │     │  │
   │     ▼  ▼
   ├──> [ CAMADA 4: simulacao_3_3_bateria_dp_aukf.py ]
   │     │  Banco de Baterias BESS NMC Thevenin 2-RC (38,88 kAh) + Eletrotermico (Rth=0,14 K/W)
   │     │  Estimacao de SoC/R0 por AUKF / SR-UKF com Adaptacao Seletiva em Q(R0)
   │     │  │
   │     ▼  ▼
   └──> [ CAMADA 5: simulacao_3_1_barramento_hvdc.py & calcular_tempo_critico_sspc.py ]
         │  Barramento 400 V, 1 MW continuo / 2,2 MW pico, Cargas CPL (-R_dyn)
         │  Barramento Laminar Planar (0,314 uH), Hold-up e Protecao SSPC (tempo critico)
         │
         ▼
 [ CAMADA 6: simulacao_4_5_tems_nmpc.py ]  ──> Gestao Termica Preditiva NMPC (50 kW / no)
         │
         ▼
 [ CAMADA 7: simulacao_4_6_fadiga_bga_cme.py ] ──> Fadiga Termomecanica BGA & PUE Orbital
         │
         ▼
 ═══════════════════════════════════════════════════════════════════════════════
 [ CAMADA 8: simula_eps_completo.py ]
   Motor Integrado de Co-Simulacao Multifisica (Artigo 02 - REMUNOM / Elsevier)
   Loop RK45 unificado, 5 sementes estocasticas, fechamento de energia (-3,9% / 4008 m2)
   Saida: resultados_eps.json, painel_eps.html, Tabelas 1 a 8
 ═══════════════════════════════════════════════════════════════════════════════
         │
         ├─────────────────────────────────────────┐
         ▼                                         ▼
 [ CAMADA 9: gerar_figuras_pdf.py ]       [ CAMADA 10: gerar_artigo_remunom.py ]
   Exporta figuras vetoriais PDF            Construtor do Manuscrito Word (.docx)
   (fig_termica, fig_mppt,                  Injecao parametrica exata, 17 equacoes OMML,
    fig_bateria, fig_barramento)            8 tabelas ABNT e 26 referencias consolidadas.
""")
    lines.append(bar("="))
    lines.append("")
    lines.append("")

    # PARTE 2
    lines.append(bar("="))
    lines.append("PARTE 2: TABELA DE ENCADEAMENTO DOS SUBSISTEMAS E MATRIZ DE INTERDEPENDÊNCIAS")
    lines.append(bar("="))
    lines.append("")
    for s in scripts:
        lines.append(f"ETAPA {s['id'].zfill(2)}: {s['name']}")
        lines.append(f"  Categoria   : {s['category']}")
        lines.append(f"  Subsistema  : {s['subsystem']}")
        lines.append(f"  Entradas    : {s['inputs']}")
        lines.append(f"  Saidas      : {s['outputs']}")
        lines.append(f"  Consumidores: {s['downstream']}")
        lines.append(bar("-"))
    lines.append("")
    lines.append("")

    # PARTE 3
    lines.append(bar("="))
    lines.append("PARTE 3: SCRIPT MESTRE DE ORQUESTRAÇÃO / EXECUÇÃO DO PIPELINE COMPLETO")
    lines.append(bar("="))
    lines.append("""
O script Python a seguir pode ser executado diretamente para orquestrar e rodar
todo o pipeline em sequencia estrita, verificando o status de saida de cada etapa.

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ==============================================================================
# ORQUESTRADOR MESTRE DO PIPELINE DE SIMULAÇÃO - EPS 1 MW LEO
# ==============================================================================
import os
import sys
import subprocess
import time

PIPELINE_STEPS = [
    {
        "step": 1,
        "name": "Fluxos Orbitais e Balanco Termico",
        "cwd": r"c:\\ZZ - Antigravity\\Mestrado\\projeto-mestrado\\simulacoes",
        "cmd": ["python", "simulacao_4_1_fluxos_orbitais.py"],
        "desc": "Calcula irradiancia AM0, albedo, IR e temperaturas da orbita."
    },
    {
        "step": 2,
        "name": "Gerador Solar 3J e Degradacao DDD",
        "cwd": r"c:\\ZZ - Antigravity\\Mestrado\\projeto-mestrado\\simulacoes",
        "cmd": ["python", "simulacao_3_2_gerador_fotovoltaico_am0.py"],
        "desc": "Calcula curva I-V de 3 diodos, MPP, Voc, Isc e area do arranjo."
    },
    {
        "step": 3,
        "name": "Conversores GaN e Estabilidade Middlebrook",
        "cwd": r"c:\\ZZ - Antigravity\\Mestrado\\projeto-mestrado\\simulacoes",
        "cmd": ["python", "simulacao_3_4_conversores_gan_middlebrook.py"],
        "desc": "Calcula mapa de perdas GaN e sintese do amortecedor Rd-Cd."
    },
    {
        "step": 4,
        "name": "Banco de Baterias Thevenin 2-RC e AUKF",
        "cwd": r"c:\\ZZ - Antigravity\\Mestrado\\projeto-mestrado\\simulacoes",
        "cmd": ["python", "simulacao_3_3_bateria_dp_aukf.py"],
        "desc": "Simula dinamica eletroquimica e estimacao de SoC e R0."
    },
    {
        "step": 5,
        "name": "Dinamica do Barramento 400 V e Protecao SSPC",
        "cwd": r"c:\\ZZ - Antigravity\\Mestrado\\projeto-mestrado\\simulacoes",
        "cmd": ["python", "simulacao_3_1_barramento_hvdc.py"],
        "desc": "Valida estabilidade com cargas CPL e desarme SSPC."
    },
    {
        "step": 6,
        "name": "Tempo Critico de Desarme SSPC",
        "cwd": r"c:\\ZZ - Antigravity\\Mestrado\\projeto-mestrado\\simulacoes",
        "cmd": ["python", "calcular_tempo_critico_sspc.py"],
        "desc": "Calcula tempo limite para corte de falhas no barramento."
    },
    {
        "step": 7,
        "name": "Gestao Termica Preditiva TEMS-NMPC",
        "cwd": r"c:\\ZZ - Antigravity\\Mestrado\\projeto-mestrado\\simulacoes",
        "cmd": ["python", "simulacao_4_5_tems_nmpc.py"],
        "desc": "Controla preditivamente temperatura de nos de processamento."
    },
    {
        "step": 8,
        "name": "Fadiga Termomecanica BGA e PUE",
        "cwd": r"c:\\ZZ - Antigravity\\Mestrado\\projeto-mestrado\\simulacoes",
        "cmd": ["python", "simulacao_4_6_fadiga_bga_cme.py"],
        "desc": "Estima confiabilidade BGA por Coffin-Manson/Engelmaier e PUE."
    },
    {
        "step": 9,
        "name": "Co-Simulador Multifisico Integrado do Artigo 02",
        "cwd": r"c:\\ZZ - Antigravity\\Mestrado\\01 - Artigos\\Artigo 02\\scripts",
        "cmd": ["python", "simula_eps_completo.py"],
        "desc": "Gera resultados_eps.json e valida os 11 pontos de checagem."
    },
    {
        "step": 10,
        "name": "Geracao de Figuras Vetoriais PDF",
        "cwd": r"c:\\ZZ - Antigravity\\Mestrado\\01 - Artigos\\Artigo 02\\scripts",
        "cmd": ["python", "gerar_figuras_pdf.py"],
        "desc": "Exporta graficos vetoriais de alta resolucao do artigo."
    },
    {
        "step": 11,
        "name": "Construcao dos Manuscritos Word (.docx)",
        "cwd": r"c:\\ZZ - Antigravity\\Mestrado\\01 - Artigos\\Artigo 02\\scripts",
        "cmd": ["python", "gerar_artigo_remunom.py"],
        "desc": "Compila os 3 arquivos Word oficiais (.docx) para submissao."
    }
]

def main():
    print("=" * 80)
    print("INICIANDO EXECUÇÃO SEQUENCIAL DO PIPELINE DE SIMULAÇÃO")
    print("=" * 80)
    t_start_total = time.time()
    
    for item in PIPELINE_STEPS:
        step_num = item["step"]
        name = item["name"]
        cwd = item["cwd"]
        cmd = item["cmd"]
        desc = item["desc"]
        
        print(f"\\n>>> [PASSO {step_num}/11] {name}")
        print(f"    Descricao : {desc}")
        print(f"    Diretorio : {cwd}")
        print(f"    Comando   : {' '.join(cmd)}")
        
        t0 = time.time()
        try:
            res = subprocess.run(cmd, cwd=cwd, check=True, capture_output=True, text=True)
            elapsed = time.time() - t0
            print(f"    STATUS: SUCESSO ({elapsed:.2f} s)")
        except subprocess.CalledProcessError as e:
            elapsed = time.time() - t0
            print(f"    STATUS: FALHA apos {elapsed:.2f} s")
            print("    ERRO STDERR:")
            print(e.stderr)
            sys.exit(1)
            
    total_elapsed = time.time() - t_start_total
    print("\\n" + "=" * 80)
    print(f"PIPELINE COMPLETO CONCLUÍDO COM SUCESSO EM {total_elapsed:.2f} s!")
    print("=" * 80)

if __name__ == '__main__':
    main()
```
""")
    lines.append(bar("="))
    lines.append("")
    lines.append("")

    # PARTE 4
    lines.append(bar("="))
    lines.append("PARTE 4: CÓDIGO-FONTE INTEGRAL E ESTRUTURADO DOS 11 SCRIPTS ENCADEADOS")
    lines.append(bar("="))
    lines.append("")

    for s in scripts:
        lines.append(bar("#"))
        lines.append(f"# SCRIPT {s['id'].zfill(2)} DE 11: {s['name']}")
        lines.append(f"# CAMINHO: {s['path']}")
        lines.append(f"# CATEGORIA: {s['category']}")
        lines.append(f"# SUBSISTEMA: {s['subsystem']}")
        lines.append(bar("#"))
        lines.append("")
        
        if os.path.exists(s['path']):
            with open(s['path'], 'r', encoding='utf-8', errors='ignore') as fp:
                code_content = fp.read()
            lines.append(code_content)
        else:
            lines.append(f"# ERRO: Arquivo nao encontrado em {s['path']}")
            
        lines.append("")
        lines.append(bar("#"))
        lines.append(f"# FIM DO SCRIPT {s['id'].zfill(2)}: {s['name']}")
        lines.append(bar("#"))
        lines.append("")
        lines.append("")

    full_text = "\n".join(lines)
    
    # Save to both paths
    for p in [out_path_1, out_path_2]:
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "w", encoding="utf-8") as fp:
            fp.write(full_text)
        print(f"Arquivo gerado com sucesso em: {p} ({len(full_text.encode('utf-8')):,} bytes)")

if __name__ == '__main__':
    generate_pipeline_txt()
