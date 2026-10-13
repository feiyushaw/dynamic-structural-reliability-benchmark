"""
Minimalist 2D Plane-Stress Cantilever Beam Schematic Generator.
Clean, uncluttered, publication-grade layout:
  - 2D Cantilever beam geometry (L = 300 mm, H = 40 mm, B = 10 mm)
  - Clamped boundary condition at x = 0
  - Tip load F(t) at x = L
  - Core analytical closed-form equations (stress & deflection)
  - Clean load history profiles
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle, FancyArrowPatch

# Minimalist scientific styling
mpl.rcParams["font.family"] = ["Noto Sans CJK SC", "DejaVu Sans", "sans-serif"]
mpl.rcParams["axes.unicode_minus"] = False
mpl.rcParams["figure.dpi"] = 300
mpl.rcParams["savefig.dpi"] = 300
mpl.rcParams["font.size"] = 10
mpl.rcParams["lines.linewidth"] = 1.8


def plot_experimental_schematic(out_dir: str = "results"):
    os.makedirs(out_dir, exist_ok=True)
    fig, (ax_beam, ax_wave) = plt.subplots(2, 1, figsize=(11, 6.2), height_ratios=[1.2, 1.0])

    # -------------------------------------------------------------
    # 1. Clean 2D Plane-Stress Cantilever Beam Schematic
    # -------------------------------------------------------------
    ax_beam.set_xlim(-45, 345)
    ax_beam.set_ylim(-20, 60)
    ax_beam.set_aspect("equal")
    ax_beam.axis("off")

    L = 300.0  # mm
    H = 40.0   # mm

    # Clamped Wall (X = 0)
    ax_beam.plot([0, 0], [-10, H + 10], color="#2c3e50", lw=3.0)
    for y_hatch in np.arange(-10, H + 10, 6):
        ax_beam.plot([-12, 0], [y_hatch - 6, y_hatch + 4], color="#7f8c8d", lw=1.5)

    # Beam Rectangle (2D Plane Stress)
    beam_rect = Rectangle((0, 0), L, H, facecolor="#ebf5fb", edgecolor="#2980b9", lw=2.0)
    ax_beam.add_patch(beam_rect)

    # Neutral Axis (Y = 20 mm)
    ax_beam.plot([0, L], [20, 20], color="#8e44ad", lw=1.5, linestyle="-.")
    ax_beam.text(150, 21.5, "中性轴 (y = 20 mm)", color="#8e44ad", fontsize=9, ha="center")

    # Downward Force F(t) at Free Tip
    f_arrow = FancyArrowPatch((L, H + 20), (L, H + 1), arrowstyle="-|>", mutation_scale=18, color="#c0392b", lw=2.8)
    ax_beam.add_patch(f_arrow)
    ax_beam.text(L, H + 24, "F(t)", color="#c0392b", fontsize=11, fontweight="bold", ha="center")

    # Dimension lines: Length L
    ax_beam.annotate("", xy=(0, -8), xytext=(L, -8), arrowprops=dict(arrowstyle="<->", color="#34495e", lw=1.2))
    ax_beam.text(150, -16, "L = 300 mm", color="#34495e", fontsize=9.5, ha="center")

    # Dimension lines: Height H
    ax_beam.annotate("", xy=(-6, 0), xytext=(-6, H), arrowprops=dict(arrowstyle="<->", color="#34495e", lw=1.2))
    ax_beam.text(-10, 20, "H = 40 mm", color="#34495e", fontsize=9.5, ha="right", va="center")

    # Analytical Formula Box
    formula_text = (
        "【二维平面应力经典解析解】\n"
        "• 根部最大正应力: $\\sigma_{max} = \\frac{6FL}{BH^2} = 168.75$ MPa\n"
        "• 自由端静挠度: $w_{max} = \\frac{FL^3}{3EI} = 1.205$ mm\n"
        "• 一阶自振频率: $f_1 = \\frac{1.875^2}{2\\pi}\\sqrt{\\frac{EI}{\\rho A L^4}} = 371.3$ Hz ($T_1=2.69$ ms)"
    )
    ax_beam.text(15, 47, formula_text, fontsize=8.8, color="#1a5276",
                 bbox=dict(boxstyle="square,pad=0.3", fc="#f4f6f7", ec="#bdc3c7", lw=0.8))

    ax_beam.text(0, -16, "固定端 (x=0)", color="#2c3e50", fontsize=9.5, ha="center")
    ax_beam.text(L, -16, "自由端 (x=L)", color="#2c3e50", fontsize=9.5, ha="center")
    ax_beam.set_title("二维平面应力悬臂梁力学模型与解析解标定 (厚度 B = 10 mm)", fontsize=11, fontweight="bold", pad=10)

    # -------------------------------------------------------------
    # 2. Clean Load Waveforms
    # -------------------------------------------------------------
    t_ms = np.linspace(0, 150, 400)
    f_c1 = np.full_like(t_ms, 1500.0)
    f_c2 = np.zeros_like(t_ms)
    f_c3 = np.zeros_like(t_ms)

    for i, t in enumerate(t_ms):
        # Case 2: Ramp 50ms
        if t <= 50:
            f_c2[i] = 1500.0 * (t / 50.0)
        elif t <= 100:
            f_c2[i] = 1500.0
        else:
            f_c2[i] = 1500.0 * max(0.0, 1.0 - (t - 100.0) / 50.0)
        # Case 3: Impact 10ms spike
        if t < 10:
            f_c3[i] = 2250.0 * (t / 10.0)
        elif t < 25:
            ratio = (t - 10.0) / 15.0
            f_c3[i] = 2250.0 * (1.0 - ratio) + 1500.0 * ratio
        elif t <= 100:
            f_c3[i] = 1500.0
        else:
            f_c3[i] = 1500.0 * max(0.0, 1.0 - (t - 100.0) / 50.0)

    ax_wave.plot(t_ms, f_c1, "b--", lw=1.8, label="Case 1: 准静态 (1500 N 恒定)")
    ax_wave.plot(t_ms, f_c2, color="#27ae60", lw=2.0, label="Case 2: 动态急踩 (50 ms 斜坡)")
    ax_wave.plot(t_ms, f_c3, color="#c0392b", lw=2.2, label="Case 3: 冲击跺脚 (10 ms 尖峰 2250 N)")

    ax_wave.set_title("三大加载工况时间曲线 F(t)", fontsize=10.5, fontweight="bold")
    ax_wave.set_xlabel("时间 t (ms)", fontsize=9.5)
    ax_wave.set_ylabel("载荷 (N)", fontsize=9.5)
    ax_wave.set_xlim(0, 150)
    ax_wave.set_ylim(0, 2550)
    ax_wave.grid(True, linestyle="--", alpha=0.5)
    ax_wave.legend(loc="upper right", fontsize=8.5)

    plt.tight_layout()
    schematic_path = os.path.join(out_dir, "experimental_setup_schematic.png")
    plt.savefig(schematic_path)
    plt.close()
    print(f"Saved Clean Schematic: {schematic_path}")


if __name__ == "__main__":
    plot_experimental_schematic()
