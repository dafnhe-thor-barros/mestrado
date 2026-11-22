"""gerar_figuras_marco5.py -- Gera figura 2x2 (300 DPI) do Marco 5 a partir de m5_n20/m5_traces_figuras.json."""
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(DIR, "..", ".."))
FIG_DIR = os.path.join(ROOT, "figuras")
TRACES_PATH = os.path.join(DIR, "m5_n20", "m5_traces_figuras.json")

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
    ax_tj.set_title("(a) Cenário nominal: $T_j(t)$ na 20ª órbita (semente 15)")
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


if __name__ == "__main__":
    main()
