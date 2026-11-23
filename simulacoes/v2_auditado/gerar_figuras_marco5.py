"""gerar_figuras_marco5.py -- Gera figuras 2x2 (300 DPI) do Marco 5 a partir de m5_n20_4860_s10_29/."""
import csv
import json
import os
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(DIR, "..", ".."))
FIG_DIR = os.path.join(ROOT, "figuras")
CAMP_DIR = os.path.join(DIR, "m5_n20_4860_s10_29")
TRACES_PATH = os.path.join(CAMP_DIR, "m5_traces_figuras.json")

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 9.5,
    "axes.titlesize": 10.5,
    "axes.labelsize": 9.5,
    "legend.fontsize": 8.0,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "axes.grid": True,
    "grid.alpha": 0.3,
})

STYLE = {
    "S0":        {"color": "#7f7f7f", "ls": ":",  "lw": 1.4, "label": "S0 (sem controle)"},
    "S1":        {"color": "#2ca02c", "ls": "--", "lw": 1.5, "label": "S1 (estático)"},
    "R_SINT":    {"color": "#d62728", "ls": "-.", "lw": 1.6, "label": "R_SINT (PI)"},
    "TEMS_RC3":  {"color": "#1f77b4", "ls": "-",  "lw": 1.8, "label": "TEMS_RC3"},
    "TEMS_PINN": {"color": "#ff7f0e", "ls": "--", "lw": 1.8, "label": "TEMS_PINN"},
}

T_ORB = 5926.379071134441
T_ECL = 1105.3908


def shade_eclipse(ax):
    t_ecl_start = (0.5 * T_ORB - 0.5 * T_ECL) / 60.0
    t_ecl_end = (0.5 * T_ORB + 0.5 * T_ECL) / 60.0
    ax.axvspan(t_ecl_start, t_ecl_end, color="#d9d9d9", alpha=0.45, label="_nolegend_")


def main():
    with open(TRACES_PATH, encoding="utf-8") as f:
        traces = json.load(f)["traces"]

    by_key = {(r["cenario"], r["ctrl"]): r for r in traces}

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.6), sharex=True)
    ax_tj, ax_f = axes[0]
    ax_soc, ax_pert = axes[1]

    # (a) Nominal: Tj
    shade_eclipse(ax_tj)
    for ctrl in ("S0", "S1", "R_SINT", "TEMS_RC3", "TEMS_PINN"):
        r = by_key[("nominal", ctrl)]
        t_min = (np.array(r["t_s"]) - r["t_s"][0]) / 60.0
        ax_tj.plot(t_min, r["tj_C"], **STYLE[ctrl])
    ax_tj.axhline(80.0, color="black", ls="--", lw=1.0, label=r"Limite $80^\circ\mathrm{C}$")
    ax_tj.set_ylabel(r"Temperatura de junção $T_j$ ($^\circ\mathrm{C}$)")
    ax_tj.set_title("(a) Cenário nominal: $T_j(t)$ na 20ª órbita (semente 10)")
    ax_tj.legend(loc="lower right", frameon=True, ncol=2)

    # (b) Nominal: f(t)
    shade_eclipse(ax_f)
    for ctrl in ("S1", "R_SINT", "TEMS_RC3", "TEMS_PINN"):
        r = by_key[("nominal", ctrl)]
        t_min = (np.array(r["t_s"]) - r["t_s"][0]) / 60.0
        ax_f.plot(t_min, r["f_MHz"], **STYLE[ctrl])
    ax_f.set_ylabel("Frequência de relógio $f$ (MHz)")
    ax_f.set_title("(b) Cenário nominal: atuação DVFS $f(t)$ na 20ª órbita")
    ax_f.legend(loc="lower right", frameon=True, ncol=2)

    # (c) Nominal: SoC(t)
    shade_eclipse(ax_soc)
    for ctrl in ("S0", "S1", "R_SINT", "TEMS_RC3", "TEMS_PINN"):
        r = by_key[("nominal", ctrl)]
        t_min = (np.array(r["t_s"]) - r["t_s"][0]) / 60.0
        ax_soc.plot(t_min, r["soc"], **STYLE[ctrl])
    ax_soc.axhline(0.70, color="black", ls="--", lw=1.0, label=r"$SoC_{\min} = 0{,}70$")
    ax_soc.axhline(0.95, color="black", ls=":", lw=0.9, label=r"$SoC_{\max} = 0{,}95$")
    ax_soc.set_xlabel("Tempo na 20ª órbita (min)")
    ax_soc.set_ylabel("Estado de carga $SoC$ (—)")
    ax_soc.set_title("(c) Cenário nominal: $SoC(t)$ na 20ª órbita")
    ax_soc.legend(loc="lower right", frameon=True, ncol=2)

    # (d) Perturbacao r_int_alto: Tj
    shade_eclipse(ax_pert)
    for ctrl in ("S1", "R_SINT", "TEMS_RC3", "TEMS_PINN"):
        r = by_key[("r_int_alto", ctrl)]
        t_min = (np.array(r["t_s"]) - r["t_s"][0]) / 60.0
        ax_pert.plot(t_min, r["tj_C"], **STYLE[ctrl])
    ax_pert.axhline(80.0, color="black", ls="--", lw=1.0, label=r"Limite $80^\circ\mathrm{C}$")
    ax_pert.set_xlabel("Tempo na 20ª órbita (min)")
    ax_pert.set_ylabel(r"Temperatura de junção $T_j$ ($^\circ\mathrm{C}$)")
    ax_pert.set_title(r"(d) Perturbação $R_{int}$ $+50\%$ (\texttt{r\_int\_alto}): $T_j(t)$")
    ax_pert.legend(loc="lower right", frameon=True, ncol=2)

    fig.tight_layout()
    out = os.path.join(FIG_DIR, "fig_4_2_marco5_malha_fechada.png")
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print("Saved:", out, os.path.getsize(out), "bytes")

    # ---- Figura 4.3: Sintese quantitativa 2x2 (P2, Latencia, Convergencia SoC, PUE 3 vs 9) ----
    with open(os.path.join(CAMP_DIR, "m5_veredito_P1_P5.json"), encoding="utf-8") as f:
        ver = json.load(f)
    with open(os.path.join(CAMP_DIR, "latencia_isolada.json"), encoding="utf-8") as f:
        lat_raw = json.load(f)
    with open(os.path.join(CAMP_DIR, "m5_pue_janela.json"), encoding="utf-8") as f:
        pue9 = json.load(f)
    with open(os.path.join(DIR, "conv_curto", "conv_TEMS_RC3_sem_calibracao_s10_soc020_n20.json"), encoding="utf-8") as f:
        c20 = json.load(f)
    with open(os.path.join(DIR, "conv_curto", "conv_TEMS_RC3_sem_calibracao_s10_soc095_n20.json"), encoding="utf-8") as f:
        c95 = json.load(f)

    # Calcula PUE janela 3 a partir da tabela agregada do resumo
    with open(os.path.join(CAMP_DIR, "m5_resumo_para_colar.txt"), encoding="utf-8") as f:
        txt_res = f.read()
    blk = txt_res.split("bloco,cenario,")[1].split("\n\n")[0]
    lines = ["bloco,cenario," + blk.splitlines()[0]] + blk.splitlines()[1:]
    tot_j3, ok_j3 = defaultdict(int), defaultdict(float)
    for row in csv.DictReader(lines):
        tot_j3[row["ctrl"]] += int(row["n"])
        ok_j3[row["ctrl"]] += int(row["n"]) * float(row["frac_pue_ok"])

    fig2, axes2 = plt.subplots(2, 2, figsize=(11.0, 7.6))
    ax_p2, ax_lat = axes2[0]
    ax_conv, ax_pue = axes2[1]

    # (a) P2: pior Tj por cenario
    cens = ["eps_baixo", "r_int_alto", "r_int_muito_alto", "combinada"]
    p2_map = {(d["cenario"], d["ctrl"]): d["tj_pior_C"] for d in ver["P2"]["detalhe"]}
    x = np.arange(len(cens))
    w = 0.25
    ax_p2.bar(x - w, [p2_map[(c, "TEMS_RC3")] for c in cens], w, color="#1f77b4", label="TEMS_RC3")
    ax_p2.bar(x,     [p2_map[(c, "TEMS_PINN")] for c in cens], w, color="#ff7f0e", label="TEMS_PINN")
    ax_p2.bar(x + w, [p2_map[(c, "R_SINT")] for c in cens], w, color="#d62728", label="R_SINT (PI)")
    ax_p2.axhline(80.0, color="black", ls="--", lw=1.1, label=r"Limite $80^\circ\mathrm{C}$")
    ax_p2.set_xticks(x)
    ax_p2.set_xticklabels(["eps_baixo", "r_int_alto", "r_int_muito_alto", "combinada"], fontsize=8.5)
    ax_p2.set_ylim(75, 90)
    ax_p2.set_ylabel(r"Pior $T_j$ em 20 sementes ($^\circ\mathrm{C}$)")
    ax_p2.set_title("(a) Critério P2: pior $T_j$ sob perturbações (20 órbitas, $n=20$)")
    ax_p2.legend(loc="upper left", frameon=True, ncol=2)

    # (b) Latencia isolada (m5_n20_4860_s10_29)
    lat_map = {r["ctrl"]: r for r in lat_raw["resultados"]}
    cats = ["TEMS_RC3", "TEMS_PINN"]
    means = [lat_map[c]["decision_ms_mean"] for c in cats]
    p99s = [lat_map[c]["decision_ms_p99"] for c in cats]
    speedup = means[0] / means[1]
    xl = np.arange(len(cats))
    wl = 0.32
    b1 = ax_lat.bar(xl - wl/2, means, wl, color="#1f77b4", label="Média (ms)")
    b2 = ax_lat.bar(xl + wl/2, p99s,  wl, color="#aec7e8", edgecolor="#1f77b4", label=r"$p_{99}$ (ms)")
    for bars in (b1, b2):
        for bar in bars:
            h = bar.get_height()
            ax_lat.text(bar.get_x() + bar.get_width()/2, h + 1.2, f"{h:.2f}", ha="center", va="bottom", fontsize=8.5)
    ax_lat.set_xticks(xl)
    ax_lat.set_xticklabels(cats)
    ax_lat.set_ylim(0, 78)
    ax_lat.set_ylabel("Tempo de decisão completa (ms)")
    ax_lat.set_title(rf"(b) Latência isolada em \texttt{{i5-8500T}} (PINN ${speedup:.2f}\times$ mais rápida)".replace(".", "{,}"))
    ax_lat.legend(loc="upper right", frameon=True)

    # (c) Convergencia orbital (soc_end por orbita, 20 orbitas)
    orbs20 = np.arange(1, len(c20["por_orbita"]) + 1)
    soc20 = [o["soc_end"] for o in c20["por_orbita"]]
    soc95 = [o["soc_end"] for o in c95["por_orbita"]]
    ax_conv.plot(orbs20, soc20, "o-", color="#1f77b4", lw=1.6, ms=4.5, label=r"$SoC_0 = 0{,}20$")
    ax_conv.plot(orbs20, soc95, "s--", color="#ff7f0e", lw=1.6, ms=4.5, label=r"$SoC_0 = 0{,}95$")
    ax_conv.axvline(4, color="gray", ls=":", lw=1.1, label="Convergência (órbita 4)")
    ax_conv.set_xlabel("Órbita")
    ax_conv.set_ylabel(r"$SoC$ ao final da órbita (—)")
    ax_conv.set_title("(c) Convergência de $SoC_0$ e batimento de fase (~4,4 órbitas)")
    ax_conv.set_xticks(range(2, 21, 2))
    ax_conv.legend(loc="lower right", frameon=True)

    # (d) PUE periodico: janela 3 vs janela 9
    ctrls_pue = ["R_SINT", "TEMS_RC3", "TEMS_PINN", "TEMS_SOLVER", "S1"]
    j3 = [100.0 * ok_j3[c] / tot_j3[c] for c in ctrls_pue]
    j9 = [100.0 * pue9["por_ctrl"][c]["final_ok"] / pue9["por_ctrl"][c]["final_tot"] for c in ctrls_pue]
    print("PUE J3:", {c: round(v, 2) for c, v in zip(ctrls_pue, j3)})
    print("PUE J9:", {c: round(v, 2) for c, v in zip(ctrls_pue, j9)})
    print("P5 atribuiveis:", [(d["bloco"], d["cenario"], d["sensor"], d["quant"], round(100*d["dif_media"], 3)) for d in ver["P5"]["detalhe"] if d["atribuivel_a_PINN"]])
    xp = np.arange(len(ctrls_pue))
    wp = 0.34
    ax_pue.bar(xp - wp/2, j3, wp, color="#9edae5", edgecolor="#17becf", label="Janela 3 órbitas")
    ax_pue.bar(xp + wp/2, j9, wp, color="#2ca02c", label="Janela 9 órbitas (emenda)")
    ax_pue.set_xticks(xp)
    ax_pue.set_xticklabels(ctrls_pue, fontsize=8.0)
    ax_pue.set_ylim(0, 115)
    ax_pue.set_ylabel("Execuções com PUE periódico (%)")
    ax_pue.set_title(r"(d) Periodicidade do PUE ($|\Delta E|/E_{TI} < 10^{-3}$, 1.728 exec.)")
    ax_pue.legend(loc="lower right", frameon=True)

    fig2.tight_layout()
    out2 = os.path.join(FIG_DIR, "fig_4_3_marco5_sintese_quantitativa.png")
    fig2.savefig(out2, dpi=300)
    plt.close(fig2)
    print("Saved:", out2, os.path.getsize(out2), "bytes")


if __name__ == "__main__":
    main()
