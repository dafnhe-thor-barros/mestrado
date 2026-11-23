"""gerar_esquematico_planta.py -- Gera esquematico fisico 300 DPI do acoplamento eletrotermico do no orbital (PROJETO_1P0)."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(DIR, "..", ".."))
OUT = os.path.join(ROOT, "figuras", "fig_3_5_esquematico_planta_fisica.png")


def box(ax, x, y, w, h, text, fc="#ffffff", ec="#2b4c7e", lw=1.4, fontsize=8.3):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.018",
                                facecolor=fc, edgecolor=ec, linewidth=lw))
    lines = text.split("\n")
    ax.text(x + w / 2, y + h * 0.72, lines[0], ha="center", va="center",
            fontsize=fontsize, fontweight="bold", color="#1a1a1a", family="serif")
    if len(lines) > 1:
        ax.text(x + w / 2, y + h * 0.33, "\n".join(lines[1:]), ha="center", va="center",
                fontsize=fontsize - 0.8, color="#333333", family="serif")


def arrow(ax, x1, y1, x2, y2, label="", lx=None, ly=None, color="#333333"):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=1.35))
    if label:
        tx = (x1 + x2) / 2 if lx is None else lx
        ty = (y1 + y2) / 2 if ly is None else ly
        ax.text(tx, ty, label, ha="center", va="center", fontsize=7.5,
                color="#222222", family="serif",
                bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.88))


def main():
    fig, ax = plt.subplots(figsize=(11.5, 5.8))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Subsistema Eletrico (topo)
    ax.add_patch(FancyBboxPatch((0.02, 0.54), 0.96, 0.43, boxstyle="round,pad=0.01,rounding_size=0.02",
                                facecolor="#eef5fb", edgecolor="#3b6ea5", lw=1.2, ls="--"))
    ax.text(0.035, 0.945, "SUBSISTEMA ELÉTRICO-ENERGÉTICO (PROJETO_1P0 — 50 kW EOL)",
            fontsize=9.0, fontweight="bold", color="#234973", family="serif")

    box(ax, 0.04, 0.72, 0.20, 0.18,
        "Arranjo FV AM0\nSpectrolab XTJ (59,65 cm²)\n27.073 cél. | 161,49 m²\n50,0 kW EOL a 65 °C",
        fc="#fffde7", ec="#c68400")
    box(ax, 0.29, 0.72, 0.16, 0.18,
        "Conversor Buck\n$\\eta_{buck} = 0{,}91$\n$Q_{buck} = 4{,}50$ kW",
        fc="#ffffff", ec="#3b6ea5")
    box(ax, 0.50, 0.72, 0.21, 0.18,
        "Barramento HVDC 400 V\nSol: 45,5 kW | Ecl.: BESS\nAuxiliar: 1,0 kW\n(903,2 W calor + 96,8 W ópt.)",
        fc="#e3f2fd", ec="#1565c0")
    box(ax, 0.77, 0.72, 0.19, 0.18,
        "BESS Saft VES 16\n19P108S (2.052 cél., 236 kg)\n$E_{EOL} = 26{,}59$ kWh\n$SoC \\in [0{,}70;\\; 0{,}95]$",
        fc="#e8f5e9", ec="#2e7d32")

    box(ax, 0.50, 0.56, 0.21, 0.11,
        "Conversor DAB (400→48 V)\n$\\eta_{dab} = 0{,}9655$ | $Q_{dab} = 1{,}36$ kW",
        fc="#ffffff", ec="#3b6ea5", fontsize=8.0)

    # Subsistema Termico (base)
    ax.add_patch(FancyBboxPatch((0.02, 0.04), 0.96, 0.44, boxstyle="round,pad=0.01,rounding_size=0.02",
                                facecolor="#fdf3ef", edgecolor="#c0392b", lw=1.2, ls="--"))
    ax.text(0.035, 0.455, "PLANTA TÉRMICA ACOPLADA (EDP 1D + EDO — Conservação Estrita de Energia)",
            fontsize=9.0, fontweight="bold", color="#922b21", family="serif")

    box(ax, 0.04, 0.14, 0.24, 0.24,
        "Nó de Junção ($T_j$)\n100x NVIDIA A100 SXM4\nDVFS $f \\in [300;\\; 1410]$ MHz\n$C_{mod}=500$ J/K | $R_{int}=0{,}08$ K/W",
        fc="#ffebee", ec="#c62828")
    box(ax, 0.34, 0.14, 0.24, 0.24,
        "Coletor / Fluido ($T_{raiz}$)\n$C_{fluid} = 50.000$ J/K\nRecebe $P_{TI}/R_{int}$ + $Q_{rest}$\nSensor RTD Pt100 (1 Hz)",
        fc="#fff3e0", ec="#e65100")
    box(ax, 0.64, 0.14, 0.32, 0.24,
        "Aleta Radiativa 1D (Al 6061-T6)\n$L=75$ mm | $\\delta=2$ mm | $c_a=4860$ J/(m²K)\n$A_{pl}=49{,}90$ m² ($A_{emit}=99{,}80$ m² bifacial)\nZ-93 EOL ($\\alpha_s=0{,}36, \\varepsilon=0{,}90$) $\\rightarrow T_s=3$ K",
        fc="#fce4ec", ec="#ad1457")

    # Arrows Eletrico
    arrow(ax, 0.24, 0.81, 0.29, 0.81, "50 kW", ly=0.84)
    arrow(ax, 0.45, 0.81, 0.50, 0.81, "45,5 kW", ly=0.84)
    arrow(ax, 0.71, 0.83, 0.77, 0.83, "Carga", ly=0.86, color="#2e7d32")
    arrow(ax, 0.77, 0.78, 0.71, 0.78, "Desc.", ly=0.75, color="#c62828")
    arrow(ax, 0.605, 0.72, 0.605, 0.67, "")

    # Coupling Eletrico -> Termico
    arrow(ax, 0.50, 0.61, 0.16, 0.38, "$P_{TI}$ (38,18 kW sol / 18 kW ecl.)", lx=0.25, ly=0.505, color="#1565c0")
    arrow(ax, 0.41, 0.72, 0.46, 0.38, "$Q_{rest}$ (perdas Buck, DAB, BD + $Q_{aux}$)", lx=0.52, ly=0.505, color="#c0392b")

    # Arrows Termico
    arrow(ax, 0.28, 0.26, 0.34, 0.26, "$\\Delta T/R_{int}$", ly=0.30, color="#c0392b")
    arrow(ax, 0.58, 0.26, 0.64, 0.26, "Condução\n$x=0$", ly=0.31, color="#c0392b")
    arrow(ax, 0.80, 0.14, 0.80, 0.055, "Radiação bifacial $2\\varepsilon\\sigma(T^4-T_s^4)$", lx=0.80, ly=0.095, color="#ad1457")

    fig.tight_layout()
    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved:", OUT, os.path.getsize(OUT), "bytes")


if __name__ == "__main__":
    main()
