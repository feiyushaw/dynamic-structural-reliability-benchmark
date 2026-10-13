"""
Figure 1: Minimalist 2D Model Schematic & Load Waveforms.
Zero clutter, clean scientific style.
Saves to results/fig1_model_and_loading.png.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle, FancyArrowPatch

mpl.rcParams["font.family"] = ["Noto Sans CJK SC", "DejaVu Sans", "sans-serif"]
mpl.rcParams["axes.unicode_minus"] = False
mpl.rcParams["figure.dpi"] = 300
mpl.rcParams["savefig.dpi"] = 300
mpl.rcParams["font.size"] = 9.5
mpl.rcParams["lines.linewidth"] = 1.8


def plot_model_and_loading(out_dir: str = "results"):
    os.makedirs(out_dir, exist_ok=True)
    fig, (ax_beam, ax_load) = plt.subplots(2, 1, figsize=(10, 5.2), height_ratios=[1.1, 1.0])

    # 1. 2D Plane-Stress Beam Layout
    L, H = 300.0, 40.0
    ax_beam.set_xlim(-40, 340)
    ax_beam.set_ylim(-18, 55)
    ax_beam.set_aspect("equal")
    ax_beam.axis("off")

    # Clamped Wall (x=0)
    ax_beam.plot([0, 0], [-8, H + 8], color="#2c3e50", lw=2.5)
    for y_hatch in np.arange(-8, H + 8, 6):
        ax_beam.plot([-10, 0], [y_hatch - 5, y_hatch + 4], color="#7f8c8d", lw=1.2)

    # Beam Body
    ax_beam.add_patch(Rectangle((0, 0), L, H, facecolor="#ebf5fb", edgecolor="#2980b9", lw=1.8))
    # Neutral Axis
    ax_beam.plot([0, L], [20, 20], color="#8e44ad", lw=1.2, linestyle="-.")
    ax_beam.text(150, 21.5, "中性轴 (y = 20 mm)", color="#8e44ad", fontsize=8.5, ha="center")

    # Force Arrow
    ax_beam.add_patch(FancyArrowPatch((L, H + 18), (L, H + 1), arrowstyle="-|>", mutation_scale=16, color="#c0392b", lw=2.5))
    ax_beam.text(L, H + 21, "F(t)", color="#c0392b", fontsize=10.5, fontweight="bold", ha="center")

    # Dimensions
    ax_beam.annotate("", xy=(0, -6), xytext=(L, -6), arrowprops=dict(arrowstyle="<->", color="#34495e", lw=1.0))
    ax_beam.text(150, -14, "L = 300 mm", color="#34495e", fontsize=9, ha="center")
    ax_beam.annotate("", xy=(-5, 0), xytext=(-5, H), arrowprops=dict(arrowstyle="<->", color="#34495e", lw=1.0))
    ax_beam.text(-8, 20, "H = 40 mm", color="#34495e", fontsize=9, ha="right", va="center")

    ax_beam.text(10, 44, "解析解: $\\sigma_{max}=\\frac{6FL}{BH^2}=168.8$ MPa, $w=\\frac{FL^3}{3EI}=1.21$ mm (B = 10 mm)",
                 fontsize=8.5, color="#1b4f72", bbox=dict(boxstyle="square,pad=0.25", fc="#f8f9f9", ec="#ccd1d1", lw=0.6))
    ax_beam.text(0, -14, "x = 0 (固定端)", color="#2c3e50", fontsize=8.5, ha="center")
    ax_beam.text(L, -14, "x = L (载荷端)", color="#2c3e50", fontsize=8.5, ha="center")
    ax_beam.set_title("(a) 二维悬臂梁力学模型与几何尺寸", fontsize=10.5, fontweight="bold", pad=8)

    # 2. Three Load Cases Waveforms
    t = np.linspace(0, 150, 400)
    f1 = np.full_like(t, 1500.0)
    f2 = np.array([1500.0 * (ti / 50.0) if ti <= 50 else (1500.0 if ti <= 100 else 1500.0 * max(0.0, 1.0 - (ti - 100.0) / 50.0)) for ti in t])
    f3 = np.array([2250.0 * (ti / 10.0) if ti < 10 else (2250.0 * (1.0 - (ti - 10.0) / 15.0) + 1500.0 * ((ti - 10.0) / 15.0) if ti < 25 else (1500.0 if ti <= 100 else 1500.0 * max(0.0, 1.0 - (ti - 100.0) / 50.0))) for ti in t])

    ax_load.plot(t, f1, "b--", lw=1.6, label="Case 1: 准静态 (1500 N)")
    ax_load.plot(t, f2, color="#27ae60", lw=1.8, label="Case 2: 动态急踩 (50 ms 斜坡)")
    ax_load.plot(t, f3, color="#c0392b", lw=2.0, label="Case 3: 冲击跺脚 (10 ms 尖峰 2250 N)")
    ax_load.set_title("(b) 三大工况载荷时间历程 F(t)", fontsize=10.5, fontweight="bold", pad=6)
    ax_load.set_xlabel("时间 t (ms)", fontsize=9.5)
    ax_load.set_ylabel("载荷 (N)", fontsize=9.5)
    ax_load.set_xlim(0, 150)
    ax_load.set_ylim(0, 2500)
    ax_load.grid(True, linestyle="--", alpha=0.5)
    ax_load.legend(loc="upper right", fontsize=8.5)

    plt.tight_layout()
    fig_path = os.path.join(out_dir, "fig1_model_and_loading.png")
    plt.savefig(fig_path)
    plt.close()
    print(f"Saved: {fig_path}")


if __name__ == "__main__":
    plot_model_and_loading()
