"""gerar_diagrama_controle_tems.py -- Gera diagrama de blocos 300 DPI da arquitetura de controle TEMS em malha fechada."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(DIR, "..", ".."))
OUT = os.path.join(ROOT, "figuras", "fig_3_4_arquitetura_controle_tems.png")


def box(ax, x, y, w, h, text, fc="#f7f9fc", ec="#2b4c7e", lw=1.4, fontsize=8.5, bold_first=True):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.015,rounding_size=0.018",
                       facecolor=fc, edgecolor=ec, linewidth=lw)
    ax.add_patch(p)
    lines = text.split("\n")
    if bold_first and len(lines) > 1:
        ax.text(x + w / 2, y + h * 0.72, lines[0], ha="center", va="center",
                fontsize=fontsize, fontweight="bold", color="#1a1a1a", family="serif")
        ax.text(x + w / 2, y + h * 0.34, "\n".join(lines[1:]), ha="center", va="center",
                fontsize=fontsize - 0.8, color="#333333", family="serif")
    else:
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
                fontsize=fontsize, fontweight="bold" if bold_first else "normal",
                color="#1a1a1a", family="serif")


def arrow(ax, x1, y1, x2, y2, label="", lx=None, ly=None, color="#333333", rad=0.0):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=1.3,
                                connectionstyle=f"arc3,rad={rad}"))
    if label:
        tx = (x1 + x2) / 2 if lx is None else lx
        ty = (y1 + y2) / 2 if ly is None else ly
        ax.text(tx, ty, label, ha="center", va="center", fontsize=7.6,
                color="#222222", family="serif",
                bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85))


def main():
    fig, ax = plt.subplots(figsize=(11.5, 6.2))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Background group: Planta Fisica
    ax.add_patch(FancyBboxPatch((0.02, 0.56), 0.96, 0.41, boxstyle="round,pad=0.01,rounding_size=0.02",
                                facecolor="#eef4fa", edgecolor="#4a76a8", lw=1.2, ls="--"))
    ax.text(0.035, 0.945, "PLANTA FÍSICA E AMBIENTE ORBITAL (passo de integração $\\Delta t_{sim} = 5$ s)",
            fontsize=9.0, fontweight="bold", color="#2b4c7e", family="serif")

    # Background group: Controlador TEMS
    ax.add_patch(FancyBboxPatch((0.02, 0.03), 0.96, 0.45, boxstyle="round,pad=0.01,rounding_size=0.02",
                                facecolor="#fdf7ee", edgecolor="#c77d20", lw=1.2, ls="--"))
    ax.text(0.035, 0.455, "CONTROLADOR PREDITIVO TEMS / DVFS ($\\Delta t_{ctrl} = 60$ s, horizonte $W = 300$ s)",
            fontsize=9.0, fontweight="bold", color="#9a5b0c", family="serif")

    # Top row (Planta)
    box(ax, 0.04, 0.72, 0.21, 0.18,
        "Ambiente Orbital\nSSO 700 km, $T_{orb}=5926$ s\nSol/Eclipse, Albedo, IR",
        fc="#ffffff", ec="#4a76a8")
    box(ax, 0.31, 0.72, 0.25, 0.18,
        "Barramento HVDC + BESS\nPV 50 kW, Buck, DAB, BD\n$P_{aux}=1$ kW, $SoC\\in[0{,}70;0{,}95]$",
        fc="#ffffff", ec="#4a76a8")
    box(ax, 0.62, 0.72, 0.24, 0.18,
        "Planta Térmica 1D + EDO\n100x A100 DVFS + Coletor\nAleta Al 6061-T6 (Z-93 EOL)",
        fc="#ffffff", ec="#4a76a8")
    box(ax, 0.89, 0.72, 0.08, 0.18,
        "Sensores\n5 s\nRTD / Diodo\n(NVML)",
        fc="#e8f5e9", ec="#2e7d32", fontsize=8.0)

    # Bottom row (TEMS)
    box(ax, 0.04, 0.24, 0.23, 0.17,
        "Supervisor Energético\nSol: $P_{TI,sol,sust}=38{,}18$ kW\nEclipse: $(E-0{,}70E_{EOL})/t_{rem}$",
        fc="#ffffff", ec="#c77d20")
    box(ax, 0.04, 0.05, 0.23, 0.14,
        "Teto de Frequência $f_{hi}$\n$\\min(f_{max},\\; f_{dem},\\; f_{cap})$",
        fc="#fff3e0", ec="#c77d20")
    box(ax, 0.33, 0.12, 0.23, 0.24,
        "Otimizador por Bisseção\n7 iterações ($\\approx 9$ MHz)\n$\\max f \\in [300, f_{hi}]$ t.q.\n$\\max \\hat{T}_j \\leq 80^\\circ\\mathrm{C} - m_T$",
        fc="#fff3e0", ec="#b25900")
    box(ax, 0.61, 0.25, 0.17, 0.16,
        "Preditor PINN (v2c)\nMLP 4x64 (7,75 ms)\nBarreira OOD",
        fc="#e3f2fd", ec="#1565c0")
    box(ax, 0.61, 0.05, 0.17, 0.15,
        "Preditor / Fallback RC3\n3 nós (48,75 ms)",
        fc="#f3e5f5", ec="#6a1b9a")
    box(ax, 0.82, 0.15, 0.15, 0.21,
        "Estimador $\\hat{\\Delta}_{ponta}$\nMín. Quadrados\nJanela $w_{obs}=60$ s",
        fc="#e8f5e9", ec="#2e7d32")

    # Arrows Planta
    arrow(ax, 0.25, 0.81, 0.31, 0.81, "Sol/Ecl.")
    arrow(ax, 0.145, 0.72, 0.74, 0.72, "$q_{abs}(t)$", lx=0.45, ly=0.635, rad=0.14)
    arrow(ax, 0.56, 0.81, 0.62, 0.81, "$Q_{rest}, f_{eff}$", ly=0.845)
    arrow(ax, 0.86, 0.81, 0.89, 0.81, "")

    # Feedback Planta -> TEMS
    arrow(ax, 0.93, 0.72, 0.93, 0.36, "$\\tilde{T}_j, \\tilde{T}_{raiz}$", lx=0.955, ly=0.52)
    arrow(ax, 0.43, 0.72, 0.155, 0.41, "$E(t)$", lx=0.25, ly=0.52)

    # Inside TEMS
    arrow(ax, 0.155, 0.24, 0.155, 0.19, "$f_{cap}$", lx=0.185)
    arrow(ax, 0.27, 0.12, 0.33, 0.20, "$f_{hi}$", ly=0.14)
    arrow(ax, 0.82, 0.29, 0.78, 0.33, "")
    arrow(ax, 0.82, 0.21, 0.78, 0.13, "")
    arrow(ax, 0.695, 0.25, 0.695, 0.20, "OOD", lx=0.725, color="#b71c1c")
    arrow(ax, 0.61, 0.33, 0.56, 0.28, "$\\hat{T}_j(\\tau;f)$", lx=0.585, ly=0.35)
    arrow(ax, 0.61, 0.13, 0.56, 0.18, "")

    # Command TEMS -> Planta
    arrow(ax, 0.445, 0.36, 0.445, 0.72, "Comando $f^\\star(t)$ (DVFS único)", lx=0.525, ly=0.52, color="#b25900")
    arrow(ax, 0.445, 0.58, 0.68, 0.72, "", color="#b25900")

    fig.tight_layout()
    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved:", OUT, os.path.getsize(OUT), "bytes")


if __name__ == "__main__":
    main()
