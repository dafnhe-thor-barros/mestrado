#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ==============================================================================
# SIMULAÇÃO 3.4: CONVERSORES ESTÁTICOS CC-CC GaN DE ALTA FREQUÊNCIA,
# BALANÇO DE PERDAS TÉRMICAS E ESTABILIDADE CPL POR CRITÉRIO DE MIDDLEBROOK
# ==============================================================================
# Dissertação de Mestrado Profissional - ProPECAut / IFES
# Autor: Dafnhê Thor Castro de Martins Barros
# Revisão Técnica: Modelagem Explícita de Condições Nominais e Modos de Falha
# ==============================================================================

import math
import os
import csv
from collections import namedtuple
import numpy as np
import matplotlib.pyplot as plt

# ---- primário GS66516T (8 chaves, 2 por braço)
RON25, ALPHA_PRI = 25e-3, 0.0128                 # [DS] / [VERIFICAR curva Ron x Tj]
RTH_JC, RTH_TIM = 0.27, 0.85                     # [DS] / [HIP]
RTH = RTH_JC + RTH_TIM
I_IN, D_EFF, N_PAR_PRI = 51.2, 0.45, 2
V_IN, FSW = 400.0, 150e3
TR, TF, ISW_ON = 12e-9, 15e-9, 25.0              # [HIP]
COSS, QG, VGATE = 207e-12, 14.2e-9, 6.0          # [VERIFICAR] Qg (12,1 nC máx em listagem comercial)
IRMS = I_IN * math.sqrt(D_EFF)                   # 34,35 A por chave (arranjo de 4 chaves)

# ---- secundário IG1NT052N10R
R_SEC25, ALPHA_SEC, I_D_MAX = 6e-3, 0.004, 52.0  # [DS-web] / [EST] / [DS-web]
N_PAR_SEC, R_INTERLIG, I_OUT = 12, 0.732e-3, 416.7   # derating adotado / [HIP]

# ---- magnéticos
K_FE, A_F, B_F, V_E, DB = 0.70, 1.46, 2.52, 4.5e-4, 0.15
RAC_PRI, RAC_SEC = 15e-3, 0.25e-3

# ---- capacitores 48 V
C_NOM_UNIT, FRAC_BIAS, DF_MAX = 22e-6, 0.50, 0.025   # [DS] / [EST] / [DS]
C_BANCO_BOL, FATOR_EOL = 15e-3, 0.70
ESR_BOL, ESR_EOL = 1.2e-3, 2.4e-3                    # [HIP]
LIM_48 = 11.52                                       # mΩ

Perdas = namedtuple("Perdas", "p_cp p_sw p_cs_disp p_cs_int p_core p_cu total q_primario eta")


def ron_pri(tj):
    return RON25 * (1.0 + ALPHA_PRI * (tj - 25.0))


def r_sec(tj=80.0):
    return R_SEC25 * (1.0 + ALPHA_SEC * (tj - 25.0))


def p_sw_dev(n):
    p_ov = V_IN * ISW_ON / n * FSW * (TR + TF) / 2.0
    return p_ov + 0.5 * COSS * V_IN ** 2 * FSW + QG * VGATE * FSW


def tj_eq(tb, n=N_PAR_PRI):
    a = (IRMS / n) ** 2 * RON25
    c0 = a * (1.0 - ALPHA_PRI * 25.0) + p_sw_dev(n)
    return (tb + RTH * c0) / (1.0 - RTH * a * ALPHA_PRI)


def tb_max(tj, n=N_PAR_PRI):
    a = (IRMS / n) ** 2 * RON25
    c0 = a * (1.0 - ALPHA_PRI * 25.0) + p_sw_dev(n)
    return tj * (1.0 - RTH * a * ALPHA_PRI) - RTH * c0


def perdas(n_pri=N_PAR_PRI, n_sec=N_PAR_SEC, r_int=R_INTERLIG):
    nsw = 4 * n_pri
    r_disp = 2.0 * r_sec() / n_sec
    p_cp = nsw * (IRMS / n_pri) ** 2 * ron_pri(80.0)
    p_sw = nsw * p_sw_dev(n_pri)
    p_cs_disp = I_OUT ** 2 * r_disp
    p_cs_int = I_OUT ** 2 * r_int
    p_core = K_FE * FSW ** A_F * DB ** B_F * V_E
    p_cu = 2.0 * IRMS ** 2 * RAC_PRI + I_OUT ** 2 * RAC_SEC   # enrolamento = sqrt(2) * I_chave
    tot = p_cp + p_sw + p_cs_disp + p_cs_int + p_core + p_cu
    eta = 100.0 * 20000.0 / (20000.0 + tot)
    return Perdas(p_cp, p_sw, p_cs_disp, p_cs_int, p_core, p_cu, tot, p_cp + p_sw, eta)


def zout_pico(L, R, C, ESR, n=8001):
    f = np.logspace(1, 5.3, n)
    s = 2j * np.pi * f
    z = np.abs((s * L + R) * (s * C * ESR + 1.0) / (s * s * L * C + s * C * (R + ESR) + 1.0))
    k = int(np.argmax(z))
    return z[k] * 1e3, f[k]


def zout_curva(L, R, C, ESR, n=8001):
    f = np.logspace(1, 5.3, n)
    s = 2j * np.pi * f
    z = np.abs((s * L + R) * (s * C * ESR + 1.0) / (s * s * L * C + s * C * (R + ESR) + 1.0))
    return f, z * 1e3


def regioes_ok(L, R, esr, limite, c_min=1.0e-3, c_max=15.0e-3, n=141):
    cs = np.linspace(c_min, c_max, n)
    ok = [zout_pico(L, R, c, esr)[0] <= limite for c in cs]
    regs, ini = [], None
    for i in range(len(ok)):
        if ok[i] and ini is None:
            ini = i
        if (not ok[i]) and ini is not None:
            regs.append((cs[ini], cs[i - 1]))
            ini = None
    if ini is not None:
        regs.append((cs[ini], cs[-1]))
    return regs


def degrau_cpl(c, esr, vs=400.0, L=1.8e-6, R=14e-3, p1=2000.0, p2=20000.0, dt=5e-7, t_fim=0.01):
    vb = vs
    for _ in range(60):
        vb = vs - R * p1 / vb
    i, vc, v0, vmin = p1 / vb, vb, vb, 1e9
    
    t_arr = []
    v_arr = []
    for k in range(int(t_fim / dt)):
        t = k * dt
        p = p1 if t < 0.5e-3 else p2
        vbb = vc
        for _ in range(8):
            vbb = vc + esr * (i - p / vbb)
        il = p / vbb
        i += (vs - vbb - R * i) / L * dt
        vc += (i - il) / c * dt
        vmin = min(vmin, vbb)
        t_arr.append(t)
        v_arr.append(vbb)
    return v0 - vmin, t_arr, v_arr


def main():
    out_dir = "saida_3_4"
    os.makedirs(out_dir, exist_ok=True)

    linha = "=" * 78
    print(linha)
    print("SIMULAÇÃO 3.4 (v4) - Com exportação de CSV e Gráficos")
    print(linha)
    assert abs(tj_eq(60.0, 1) - 187.0) < 0.3
    print("8 chaves: Tj(55/65)=%.1f/%.1f °C | Tc_max(Tj<=105)=%.1f °C"
          % (tj_eq(55.0), tj_eq(65.0), tb_max(105.0)))

    p = perdas()
    print("\nPerdas (W): cond_pri=%.1f sw=%.1f sec_disp=%.1f sec_interlig=%.1f nucleo=%.1f cobre=%.1f"
          % (p.p_cp, p.p_sw, p.p_cs_disp, p.p_cs_int, p.p_core, p.p_cu))
    print("Total=%.1f W | eta=%.2f %% | Q_primario=%.1f W" % (p.total, p.eta, p.q_primario))
    
    # Salvar CSV de perdas por R_interlig
    csv_perdas_path = os.path.join(out_dir, "perdas_interligacao.csv")
    with open(csv_perdas_path, mode="w", newline="", encoding="utf-8") as f_csv:
        writer = csv.writer(f_csv)
        writer.writerow(["R_interlig_mOhm", "Total_W", "Eficiencia_pct", "Secundario_W"])
        for r_int in (0.0, 0.3e-3, 0.732e-3, 1.0e-3):
            q = perdas(r_int=r_int)
            writer.writerow([f"{r_int*1e3:.3f}", f"{q.total:.2f}", f"{q.eta:.2f}", f"{q.p_cs_disp + q.p_cs_int:.2f}"])
            print("  R_interlig=%.3f mΩ -> total=%.1f W, eta=%.2f %%, secundário=%.1f W"
                  % (r_int * 1e3, q.total, q.eta, q.p_cs_disp + q.p_cs_int))

    print("Cold plate (estimativa; exige rede térmica por componente): primário R_th<=%.3f/%.3f K/W (45/55 °C); "
          "módulo inteiro R_th<=%.3f K/W (45 °C)"
          % ((tb_max(105.0) - 45.0) / p.q_primario, (tb_max(105.0) - 55.0) / p.q_primario,
             (tb_max(105.0) - 45.0) / p.total))

    c_ef = C_NOM_UNIT * FRAC_BIAS
    n_comp = math.ceil(C_BANCO_BOL / c_ef)
    esr_cap = DF_MAX / (2.0 * math.pi * 1e3 * c_ef) / n_comp
    c_eol = C_BANCO_BOL * FATOR_EOL
    print("\nBanco 48 V: C_ef=%.1f uF | N=%d | ESR_cap=%.3f mΩ | C_EOL=%.1f mF"
          % (c_ef * 1e6, n_comp, esr_cap * 1e3, c_eol * 1e3))
    
    casos = [("Saudável", C_BANCO_BOL, ESR_BOL), ("EOL de projeto", c_eol, ESR_EOL),
             ("Controle só C", c_eol, ESR_BOL), ("ESR só capacitor", c_eol, esr_cap),
             ("Severo 2,5 mF", 2.5e-3, ESR_BOL), ("Adverso 2,5 mF/9,5 mΩ", 2.5e-3, 9.5e-3)]
    
    csv_banco48_path = os.path.join(out_dir, "banco_48v_casos.csv")
    with open(csv_banco48_path, mode="w", newline="", encoding="utf-8") as f_csv:
        writer = csv.writer(f_csv)
        writer.writerow(["Caso", "Capacitancia_F", "ESR_Ohm", "Z_pico_mOhm", "Freq_Hz", "Status"])
        for nome, cc, ee in casos:
            z, f = zout_pico(0.15e-6, 1.5e-3, cc, ee)
            status = "OK" if z <= LIM_48 else "VIOLA"
            writer.writerow([nome, f"{cc}", f"{ee}", f"{z:.2f}", f"{f:.1f}", status])
            print("  %-24s %6.2f mΩ @ %6.0f Hz -> %s" % (nome, z, f, status))

    print("Regiões de C que atendem (varredura 1-15 mF):")
    for esr in (esr_cap, 1.2e-3, 2.4e-3, 4.0e-3, 9.5e-3):
        regs = regioes_ok(0.15e-6, 1.5e-3, esr, LIM_48)
        txt = "; ".join(["%.1f-%.1f mF" % (a * 1e3, b * 1e3) for a, b in regs]) if len(regs) > 0 else "nenhuma"
        print("  ESR=%.3f mΩ: %s" % (esr * 1e3, txt))

    print("\nBarramento 400 V (limite 800 mΩ; razão = limite/pico, não é margem de estabilidade):")
    csv_barramento_path = os.path.join(out_dir, "barramento_400v.csv")
    with open(csv_barramento_path, mode="w", newline="", encoding="utf-8") as f_csv:
        writer = csv.writer(f_csv)
        writer.writerow(["C_mF", "Z_pico_mOhm", "Freq_Hz", "Razao_Limite", "dV_degrau_V"])
        for c in (4.7e-3, 2.0e-3, 1.0e-3, 0.5e-3):
            z_fix, f_fix = zout_pico(1.8e-6, 14e-3, c, 1.5e-3)
            dv = degrau_cpl(c, 1.5e-3)[0]
            writer.writerow([f"{c*1e3:.1f}", f"{z_fix:.2f}", f"{f_fix:.1f}", f"{800.0/z_fix:.2f}", f"{dv:.2f}"])
            print("  C=%.1f mF | ESR fixa 1,5 mΩ: Zpk=%.1f mΩ @%.0f Hz (%.1fx) | dV=%.2f V"
                  % (c * 1e3, z_fix, f_fix, 800.0 / z_fix, dv))

    # ---- GERAÇÃO DE GRÁFICOS ----
    print(f"\nGerando gráficos em '{out_dir}/'...")

    # Gráfico 1: Impedância de Saída do Barramento 400V
    plt.figure(figsize=(9, 5))
    for c in (4.7e-3, 2.0e-3, 1.0e-3, 0.5e-3):
        f_vec, z_vec = zout_curva(1.8e-6, 14e-3, c, 1.5e-3)
        plt.loglog(f_vec, z_vec, label=f"C = {c*1e3:.1f} mF")
    plt.axhline(800.0, color='r', linestyle='--', label='Limite (800 mΩ)')
    plt.title("Impedância de Saída - Barramento 400V (CPL Stability)")
    plt.xlabel("Frequência (Hz)")
    plt.ylabel("Impedância Zout (mΩ)")
    plt.grid(True, which="both", ls="--")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "zout_barramento_400v.png"), dpi=300)
    plt.close()

    # Gráfico 2: Resposta transitória CPL (C = 2.0 mF)
    _, t_arr, v_arr = degrau_cpl(2.0e-3, 1.5e-3)
    plt.figure(figsize=(9, 5))
    plt.plot([t * 1e3 for t in t_arr], v_arr, color='b', lw=1.5)
    plt.axhline(400.0, color='gray', linestyle=':', label='Tensão nominal (400V)')
    plt.title("Resposta Transitória a Degrau de Carga CPL (C = 2.0 mF)")
    plt.xlabel("Tempo (ms)")
    plt.ylabel("Tensão do Barramento (V)")
    plt.grid(True, ls="--")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "transit_cpl_2mf.png"), dpi=300)
    plt.close()

    print("Processo concluído com sucesso! Arquivos salvos em:", out_dir)
    print(linha)


if __name__ == "__main__":
    main()