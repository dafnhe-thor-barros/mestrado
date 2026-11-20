#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ==============================================================================
# SIMULAÇÃO 3.3 (v2): MODELAGEM ELETROTÉRMICA DO BANCO DE BATERIAS (BESS),
# CIRCUITO DUAL POLARIZATION (DP 2ª ORDEM) E ESTIMAÇÃO DE SOC VIA AUKF
# ==============================================================================
# Dissertação de Mestrado Profissional - ProPECAut / IFES
# Autor: Dafnhê Thor Castro de Martins Barros
# Revisão Técnica: Cargas com OISL Contínuo (440 W), Rth = 0.14 K/W e AUKF
# ==============================================================================

import os
import math
import csv
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt

# Configuração de estilo para gráficos científicos
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.titlesize'] = 14

# Constantes Físicas Fundamentais
Q_E, K_B = 1.602176634e-19, 1.380649e-23
EA_K = 0.32 * Q_E / K_B                 # Energia de ativação de Arrhenius [K] (0,32 eV)
NS = 108                                # Células em série no pack Li-ion NMC
V_UV, V_OV = 324.0, 442.8               # Limiares de corte: UVLO (3,0 V) e CV (4,10 V)
CN_AH = 22.52                           # Capacidade nominal [Ah]
E_NOM_J = 9.0e3 * 3600.0                # Energia nominal armazenada (9,0 kWh) [J]

# Parâmetros Eletroquímicos do Modelo DP de 2ª Ordem
RP1, CP1 = 18e-3, 3200.0                # Ramo de difusão rápida: tau1 = 57,6 s
RP2, CP2 = 24e-3, 22000.0               # Ramo de difusão lenta: tau2 = 528 s
R0_REF = 85e-3                          # Resistência ôhmica degradada em EOL a 20 °C [Ohm]
R_DIST = 14e-3                          # Resistência da cablagem de distribuição [Ohm]
T_REF_K = 293.15                        # 20 °C em Kelvin

# Acoplamento Térmico Condutivo Corrigido
C_BAT = 54360.0                         # Capacitância térmica do pack (60 kg NMC) [J/K]
RTH_COND = 0.14                         # Resistência condutiva estrutural [K/W]
T_CH = 15.0                             # Temperatura média do chassi de ancoragem [°C]
P_HEATER, T_ON, T_OFF = 300.0, 10.0, 15.0

# Dinâmica Orbital SSO (700 km)
T_ORB, T_ECL = 5926.0, 1128.0            # Período = 98,77 min; Eclipse sazonal máx = 18,8 min
ETA_DAB = 0.95                          # Eficiência do conversor de distribuição (ndab = 0.95)
ETA_DIS = 0.97                          # Eficiência de descarga eletroquímica
ETA_CH = 0.995                          # Eficiência coulômbica de recarga

# Potências de Consumo Incluindo OISL Contínuo (440 W)
P_AUX_KW = 0.56                         # Aviônica e bombeamento essencial [kW]
P_OISL_KW = 0.44                        # 4 terminais ópticos inter-satélite 24/7 [kW]
P_IT_ECL_KW = 18.0                      # Carga computacional prioritária no eclipse [kW]

# Cenários de Geração Solar Pós-Buck (nbuck = 0.91) e Carga de TI em Sol Sustentável
CENARIOS_FV = {                         # nome: (P_pv líquido [kW], P_IT sustentável em sol [kW])
    "EOL_85C": (46.14, 38.10),          # Ajustado para 38.10 kW devido ao OISL contínuo
    "EOL_65C": (49.73, 40.00),
    "EOL_28C": (54.18, 40.00),
}


def ocv(soc):
    """Relação linear unificada de circuito aberto: Voc = Ns * (3.0 + 1.1 * SoC)."""
    return NS * (3.0 + 1.1 * soc)


def r0_arrhenius(t_c):
    """Variação térmica da resistência interna de condução via lei de Arrhenius."""
    return R0_REF * math.exp(EA_K * (1.0 / (t_c + 273.15) - 1.0 / T_REF_K))


def simula_orbital(p_pv_kw, p_it_sol_kw, rth=RTH_COND, t_ch=T_CH, soc0=1.0,
                   t0=20.0, n_orb=10, dt=2.0, e_nom=E_NOM_J, cn=CN_AH):
    """
    Simula a dinâmica eletrotérmica multiorbital do BESS com carga de potência
    constante (CPL) no eclipse e recarga CC-CV no período iluminado.
    """
    p_aux_tot = P_AUX_KW + P_OISL_KW  # 1.00 kW total (plataforma + OISL)
    p_chg_w = max(0.0, p_pv_kw - (p_it_sol_kw + p_aux_tot) / ETA_DAB) * 1e3
    p_ecl_w = (P_IT_ECL_KW + p_aux_tot) * 1e3 / ETA_DAB  # 20.00 kW no barramento

    n = int(math.floor(n_orb * T_ORB / dt)) + 1
    soc, vt, tb, pj, ib = (np.zeros(n) for _ in range(5))
    ecl = np.zeros(n, dtype=bool)
    soc[0], tb[0] = soc0, t0
    a1 = math.exp(-dt / (RP1 * CP1))
    a2 = math.exp(-dt / (RP2 * CP2))
    vp1 = vp2 = 0.0

    for k in range(n - 1):
        ecl[k] = (k * dt) % T_ORB < T_ECL
        r0 = r0_arrhenius(tb[k])
        rs = r0 + R_DIST
        ve = ocv(soc[k]) - vp1 - vp2

        # Eclipse: Descarga como Carga de Potência Constante (Eq. 3.10)
        if ecl[k]:
            i = (ve - math.sqrt(max(ve * ve - 4.0 * rs * p_ecl_w, 0.0))) / (2.0 * rs)
            v = ve - i * rs
            dsoc = -(v * i * dt) / (ETA_DIS * e_nom)  # Formulação energética Eq. 3.7
        # Sol: Recarga controlada em dois estágios (CC -> CV)
        else:
            ich = (-ve + math.sqrt(ve * ve + 4.0 * rs * p_chg_w)) / (2.0 * rs)
            ich = min(ich, cn)  # Teto de 1,0C
            v = ve + ich * rs
            if v > V_OV:
                ich = max(0.0, (V_OV - ve) / rs)
                v = ve + ich * rs
            i = -ich
            dsoc = ETA_CH * (v * ich * dt) / e_nom

        soc[k + 1] = min(1.0, max(0.0, soc[k] + dsoc))
        vt[k], ib[k] = v, i
        vp1 = a1 * vp1 + RP1 * (1.0 - a1) * i
        vp2 = a2 * vp2 + RP2 * (1.0 - a2) * i
        pj[k] = i * i * r0 + vp1 ** 2 / RP1 + vp2 ** 2 / RP2

        # Controle ativo dos aquecedores de sobrevivência
        if tb[k] < T_ON:
            ph = P_HEATER
        elif tb[k] < T_OFF:
            ph = P_HEATER * (T_OFF - tb[k]) / (T_OFF - T_ON)
        else:
            ph = 0.0

        # Balanço térmico condutivo corrigido com o chassi estrutural
        tb[k + 1] = tb[k] + (pj[k] + ph - (tb[k] - t_ch) / rth) * dt / C_BAT

    ecl[-1], vt[-1], pj[-1], ib[-1] = ecl[-2], vt[-2], pj[-2], ib[-2]
    i0 = int((n_orb - 1) * T_ORB / dt)  # Amostragem da última órbita estabilizada
    sl = slice(i0, n - 1)
    m = ecl[sl]
    pj_med = float(np.mean(pj[sl]))

    resumo = dict(
        p_chg_kw=p_chg_w / 1e3,
        soc_min=float(soc[i0:].min()),
        v_min=float(vt[sl].min()),
        soc_fim=float(soc[-1]),
        i_ecl_med=float(ib[sl][m].mean()),
        i_ecl_pico=float(ib[sl][m].max()),
        pj_ecl=float(pj[sl][m].mean()),
        pj_sol=float(pj[sl][~m].mean()),
        pj_med=pj_med,
        tb_min=float(tb[i0:].min()),
        tb_pico=float(tb.max()),
        t_eq=t_ch + pj_med * rth,
        tau_h=C_BAT * rth / 3600.0
    )
    
    serie = dict(t=np.arange(n) * dt, soc=soc, vt=vt, tb=tb, pj=pj, ib=ib, ecl=ecl)
    return resumo, serie


def _psd(m, piso):
    """Simetriza a matriz e impõe piso positivo nos autovalores para garantir Cholesky."""
    s = 0.5 * (m + m.T)
    w, v = np.linalg.eigh(s)
    return v @ np.diag(np.maximum(w, piso)) @ v.T


def _perfil_i(t):
    """Perfil dinâmico com rajadas periódicas simulando cargas de tensores de IA."""
    return 25.0 + 15.0 * math.sin(2.0 * math.pi * t / 180.0) + (10.0 if (t % 120.0) < 30.0 else 0.0)


def executa_aukf(sigma_v, soc_real0=0.80, soc_est0=0.60, r0_real=R0_REF, cn_real=CN_AH,
                 i_bias=0.0, adapta_q=True, lim_r=(0.1, 10.0), piso=1e-8,
                 alpha=1e-3, beta=2.0, kappa=0.0, dt=1.0, n=1201, seed=42):
    """
    Executa o estimador AUKF de Sage-Husa completo (Q e R adaptativos) sob
    incertezas paramétricas da planta (R0 degradado, Cn reduzido e bias de sensor).
    """
    ruido = np.random.RandomState(seed).randn(n - 1)
    a1, a2 = math.exp(-dt / 57.6), math.exp(-dt / 528.0)
    s_t, v1, v2, v_cl, v_ms = (np.zeros(n) for _ in range(5))
    s_t[0] = soc_real0

    for k in range(n - 1):
        i = _perfil_i(k * dt)
        s_t[k + 1] = s_t[k] - i * dt / (cn_real * 3600.0)
        v1[k + 1] = a1 * v1[k] + RP1 * (1 - a1) * i
        v2[k + 1] = a2 * v2[k] + RP2 * (1 - a2) * i
        v_cl[k + 1] = ocv(s_t[k + 1]) - i * r0_real - v1[k + 1] - v2[k + 1]
        v_ms[k + 1] = v_cl[k + 1] + sigma_v * ruido[k]
    v_cl[0], v_ms[0] = v_cl[1], v_ms[1]

    L = 3
    lam = alpha ** 2 * (L + kappa) - L
    wm = np.full(2 * L + 1, 1.0 / (2.0 * (L + lam)))
    wc = wm.copy()
    wm[0] = lam / (L + lam)
    wc[0] = wm[0] + (1.0 - alpha ** 2 + beta)

    x = np.array([soc_est0, 0.0, 0.0])
    P = np.diag([0.05, 0.01, 0.01])
    Q = np.diag([1e-7, 1e-5, 1e-5])
    R, b = sigma_v ** 2, 0.95
    s_h, v_h = np.zeros(n), np.zeros(n)
    s_h[0] = x[0]

    for k in range(n - 1):
        i = _perfil_i(k * dt) + i_bias
        S = np.linalg.cholesky((L + lam) * _psd(P, piso))
        X = np.zeros((L, 2 * L + 1))
        X[:, 0] = x
        for j in range(L):
            X[:, 1 + j], X[:, 1 + L + j] = x + S[:, j], x - S[:, j]

        Xp = np.zeros_like(X)
        Xp[0] = X[0] - i * dt / (CN_AH * 3600.0)
        Xp[1] = a1 * X[1] + RP1 * (1 - a1) * i
        Xp[2] = a2 * X[2] + RP2 * (1 - a2) * i

        xp = Xp @ wm
        p_star = sum(wc[j] * np.outer(Xp[:, j] - xp, Xp[:, j] - xp) for j in range(2 * L + 1))
        Pp = p_star + Q

        Y = ocv(Xp[0]) - i * R0_REF - Xp[1] - Xp[2]
        yp = float(Y @ wm)
        pyy0 = sum(wc[j] * (Y[j] - yp) ** 2 for j in range(2 * L + 1))
        pxy = sum(wc[j] * (Xp[:, j] - xp) * (Y[j] - yp) for j in range(2 * L + 1))

        inov = v_ms[k + 1] - yp
        dk = (1.0 - b) / (1.0 - b ** (k + 1))
        R = min(lim_r[1] * sigma_v ** 2, max(lim_r[0] * sigma_v ** 2,
                (1.0 - dk) * R + dk * (inov ** 2 - pyy0)))

        pyy = pyy0 + R
        K = pxy / pyy
        x = xp + K * inov
        P = Pp - np.outer(K, K) * pyy

        if adapta_q:
            Q = _psd((1 - dk) * Q + dk * (inov ** 2 * np.outer(K, K) + P - p_star), 1e-12)

        x[0] = min(1.0, max(0.0, x[0]))
        s_h[k + 1] = x[0]
        v_h[k + 1] = ocv(x[0]) - i * R0_REF - x[1] - x[2]

    err = np.abs(s_h - s_t) * 100.0
    idx = np.where(err > 1.0)[0]
    
    resumo = dict(
        t_conv=(int(idx[-1]) + 1) * dt if len(idx) else 0.0,
        err_med=float(err[59:].mean()),
        err_max=float(err[59:].max()),
        rmse_mV=float(1e3 * np.sqrt(np.mean((v_h[59:] - v_cl[59:]) ** 2)))
    )
    serie = dict(t=np.arange(n) * dt, s_t=s_t, s_h=s_h, err=err, v_cl=v_cl, v_h=v_h)
    return resumo, serie


def main():
    # Criação do diretório de saída
    out_dir = "saida_3_3"
    os.makedirs(out_dir, exist_ok=True)

    print("=" * 80)
    print("SIMULAÇÃO 3.3 (v2): BESS 108S1P - SSO 700 km, 10 ÓRBITAS (COM OISL CONTÍNUO)")
    print(f"Salvando resultados e gráficos em: ./{out_dir}/")
    print("=" * 80)

    # 1. Simulação Multiorbital por Cenário Solar (10 Órbitas Consecutivas)
    resultados_fv_csv = []
    
    for nome, (p_pv, p_it) in CENARIOS_FV.items():
        r, serie = simula_orbital(p_pv, p_it, n_orb=10)
        print(f"[{nome}] P_chg={r['p_chg_kw']:.3f} kW | SoC_min={100*r['soc_min']:.2f} % | "
              f"V_min={r['v_min']:.1f} V | SoC_fim={100*r['soc_fim']:.2f} %")
        print(f"   I_ecl={r['i_ecl_med']:.1f} A (pico {r['i_ecl_pico']:.1f} A) | "
              f"Pj média={r['pj_med']:.1f} W | Tb pico={r['tb_pico']:.1f} °C | Teq={r['t_eq']:.1f} °C")

        # Salvar série temporal detalhada em CSV
        csv_path_serie = os.path.join(out_dir, f"serie_temporal_orbital_{nome}.csv")
        with open(csv_path_serie, mode="w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["Tempo_s", "SoC", "Tensao_V", "Temp_Bat_C", "Perdas_J_W", "Corrente_A", "Eclipse"])
            for k in range(len(serie['t'])):
                writer.writerow([serie['t'][k], serie['soc'][k], serie['vt'][k], serie['tb'][k], serie['pj'][k], serie['ib'][k], int(serie['ecl'][k])])

        # Gerar e salvar gráfico eletrotérmico da órbita
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
        t_horas = serie['t'] / 3600.0
        
        ax1.plot(t_horas, serie['soc'] * 100.0, color='#1f77b4', lw=1.5)
        ax1.set_ylabel("SoC [%]")
        ax1.set_title(f"Dinâmica Eletrotérmica - Cenário {nome} (10 Órbitas SSO)")
        ax1.grid(True)

        ax2.plot(t_horas, serie['vt'], color='#ff7f0e', lw=1.5, label="Tensão Terminal (V)")
        ax2.set_ylabel("Tensão [V]")
        ax2.legend(loc="upper right")
        ax2.grid(True)

        ax3.plot(t_horas, serie['tb'], color='#2ca02c', lw=1.5, label="Temp. Pack (°C)")
        ax3.axhline(T_CH, color='gray', linestyle='--', label=f"Chassi ({T_CH}°C)")
        ax3.set_xlabel("Tempo [horas]")
        ax3.set_ylabel("Temperatura [°C]")
        ax3.legend(loc="upper right")
        ax3.grid(True)

        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"grafico_orbital_{nome}.png"), dpi=300)
        plt.close()

        # Acumular para CSV de resumo
        row = {"Cenario": nome}
        row.update(r)
        resultados_fv_csv.append(row)

    # Salvar CSV de resumo dos cenários solares
    if resultados_fv_csv:
        csv_resumo_fv = os.path.join(out_dir, "resumo_cenarios_fv.csv")
        with open(csv_resumo_fv, mode="w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=resultados_fv_csv[0].keys())
            writer.writeheader()
            for row in resultados_fv_csv:
                writer.writerow(row)

    # 2. Testes de Sensibilidade do Estimador AUKF Sage-Husa
    print("\n=== AVALIAÇÃO DE ROBUSTEZ DO ESTIMADOR AUKF SAGE-HUSA ===")
    casos_aukf = {
        "Modelo exato (Q e R adapt.)": {},
        "R0 real +30 % (Envelhecimento)": dict(r0_real=1.3 * R0_REF),
        "Cn real -5 % (Perda de Cap.)": dict(cn_real=0.95 * CN_AH),
        "Offset de corrente +0,5 A": dict(i_bias=0.5),
        "Caso Combinado Severo": dict(r0_real=1.3 * R0_REF, cn_real=0.95 * CN_AH, i_bias=0.5),
    }

    resultados_aukf_csv = []

    for nome, kw in casos_aukf.items():
        for sig in (0.015, 0.150):
            a, serie_aukf = executa_aukf(sig, **kw)
            sig_mV_str = f"{1e3*sig:.0f}mV"
            print(f"{nome:32s} sigma={sig_mV_str:>5s} | e_med={a['err_med']:.4f} % | "
                  f"e_max={a['err_max']:.3f} % | RMSE_V={a['rmse_mV']:.2f} mV | t_conv={a['t_conv']:.0f} s")

            # Salvar CSV detalhado do AUKF
            nome_arq_seg = nome.replace(" ", "_").replace("%", "").replace("+", "").replace("(", "").replace(")", "").replace(".", "")
            csv_path_aukf = os.path.join(out_dir, f"serie_aukf_{nome_arq_seg}_{sig_mV_str}.csv")
            with open(csv_path_aukf, mode="w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                writer.writerow(["Tempo_s", "SoC_Real", "SoC_Estimado", "Erro_Percentual", "Tensao_Planta_V", "Tensao_Estimada_V"])
                for k in range(len(serie_aukf['t'])):
                    writer.writerow([serie_aukf['t'][k], serie_aukf['s_t'][k], serie_aukf['s_h'][k], serie_aukf['err'][k], serie_aukf['v_cl'][k], serie_aukf['v_h'][k]])

            # Gerar gráfico de estimação do AUKF para o ruído de 150 mV (ou ambos)
            if sig == 0.150:
                fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
                
                ax1.plot(serie_aukf['t'], serie_aukf['s_t'] * 100.0, 'k-', lw=1.5, label="SoC Real (Planta)")
                ax1.plot(serie_aukf['t'], serie_aukf['s_h'] * 100.0, 'r--', lw=1.5, label="SoC Estimado (AUKF)")
                ax1.set_ylabel("SoC [%]")
                ax1.set_title(f"AUKF Estimação de SoC - {nome} ($\sigma_v$ = 150 mV)")
                ax1.legend(loc="upper right")
                ax1.grid(True)

                ax2.plot(serie_aukf['t'], serie_aukf['err'], color='purple', lw=1.2)
                ax2.set_xlabel("Tempo [s]")
                ax2.set_ylabel("Erro Absoluto [%]")
                ax2.grid(True)

                plt.tight_layout()
                plt.savefig(os.path.join(out_dir, f"grafico_aukf_{nome_arq_seg}.png"), dpi=300)
                plt.close()

            # Acumular para resumo
            row_a = {"Caso": nome, "Sigma_mV": 1e3*sig}
            row_a.update(a)
            resultados_aukf_csv.append(row_a)

    # Salvar CSV de resumo do AUKF
    if resultados_aukf_csv:
        csv_resumo_aukf = os.path.join(out_dir, "resumo_aukf.csv")
        with open(csv_resumo_aukf, mode="w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=resultados_aukf_csv[0].keys())
            writer.writeheader()
            for row in resultados_aukf_csv:
                writer.writerow(row)

    print("=" * 80)
    print(f"Arquivos CSV e gráficos salvos com sucesso na pasta ./{out_dir}/")
    print("=" * 80)


if __name__ == "__main__":
    main()