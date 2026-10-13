"""
Experimental Design & Boundary Conditions Schematic Diagram Generator (Intuitive Edition).
Renders publication-grade, self-explanatory Figure 1:
  - 3D Mechanical model of the Cantilever Beam (Pedal Arm / Bracket)
  - Plain Chinese annotations for boundaries, loads, and failure danger zones
  - Side-by-side comparison of load input profiles F(t) and vehicle deceleration a_x(t)
Saves to results/experimental_setup_schematic.png.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle, FancyArrowPatch, Polygon

# Publication styling
mpl.rcParams["font.family"] = ["Noto Sans CJK SC", "DejaVu Sans", "sans-serif"]
mpl.rcParams["axes.unicode_minus"] = False
mpl.rcParams["figure.dpi"] = 300
mpl.rcParams["savefig.dpi"] = 300
mpl.rcParams["font.size"] = 10.5
mpl.rcParams["lines.linewidth"] = 2.0


def plot_experimental_schematic(out_dir: str = "results"):
    os.makedirs(out_dir, exist_ok=True)
    fig = plt.figure(figsize=(15.5, 8.8))

    gs = fig.add_gridspec(2, 2, width_ratios=[1.38, 1.0], height_ratios=[1.0, 1.0], hspace=0.35, wspace=0.22)
    ax_mech = fig.add_subplot(gs[:, 0])
    ax_force = fig.add_subplot(gs[0, 1])
    ax_accel = fig.add_subplot(gs[1, 1])

    # -------------------------------------------------------------
    # 1. Mechanical Schematic (Isometric/Perspective 3D Beam)
    # -------------------------------------------------------------
    ax_mech.set_xlim(-65, 365)
    ax_mech.set_ylim(-38, 98)
    ax_mech.set_aspect("equal")
    ax_mech.axis("off")

    dx_z, dy_z = 18.0, 14.0
    L = 300.0  # mm
    H = 40.0   # mm

    p1 = np.array([0.0, 0.0])
    p2 = np.array([L, 0.0])
    p3 = np.array([L, H])
    p4 = np.array([0.0, H])

    # Draw Clamped Wall at X = 0 (Fixed Support)
    wall_y_bot = -25.0
    wall_y_top = H + dy_z + 25.0
    wall_poly = Polygon(
        [[-25, wall_y_bot], [0, wall_y_bot], [dx_z, wall_y_bot + dy_z],
         [dx_z, wall_y_top], [0, wall_y_top - dy_z], [-25, wall_y_top - dy_z]],
        closed=True, facecolor="#bdc3c7", edgecolor="#7f8c8d", lw=1.5, zorder=1
    )
    ax_mech.add_patch(wall_poly)

    for y_hatch in np.arange(wall_y_bot + 5, wall_y_top - 5, 8):
        ax_mech.plot([-22, -2], [y_hatch - 8, y_hatch + 6], color="#5d6d7e", lw=1.2, zorder=2)

    ax_mech.text(-12, (wall_y_top + wall_y_bot) / 2.0, "【固定端】\n模拟车身防火墙\n与踏板支架安装板\n(全固支 $u_x=u_y=u_z=0$)",
                 ha="right", va="center", fontsize=9.5, fontweight="bold", color="#2c3e50")

    # Top Face (Y = H)
    top_poly = Polygon([p4, p3, p3 + [dx_z, dy_z], p4 + [dx_z, dy_z]], closed=True, facecolor="#d6eaf8", edgecolor="#2980b9", lw=1.8, zorder=3)
    ax_mech.add_patch(top_poly)

    # Free End Tip Face (X = L)
    tip_poly = Polygon([p2, p2 + [dx_z, dy_z], p3 + [dx_z, dy_z], p3], closed=True, facecolor="#aed6f1", edgecolor="#2980b9", lw=1.8, zorder=3)
    ax_mech.add_patch(tip_poly)

    # Front Face (Z = 0)
    front_poly = Polygon([p1, p2, p3, p4], closed=True, facecolor="#ebf5fb", edgecolor="#2980b9", lw=2.0, zorder=3)
    ax_mech.add_patch(front_poly)

    # Hex8 Mesh Indicative Lines
    for x_line in np.linspace(0, L, 13):
        ax_mech.plot([x_line, x_line], [0, H], color="#aed6f1", lw=0.9, linestyle="--", zorder=4)
        ax_mech.plot([x_line, x_line + dx_z], [H, H + dy_z], color="#85c1e9", lw=0.8, linestyle="--", zorder=4)
    for y_line in np.linspace(0, H, 5):
        ax_mech.plot([0, L], [y_line, y_line], color="#aed6f1", lw=0.9, linestyle="--", zorder=4)
        ax_mech.plot([L, L + dx_z], [y_line, y_line + dy_z], color="#85c1e9", lw=0.8, linestyle="--", zorder=4)

    # Neutral Axis (Y = 20mm)
    ax_mech.plot([0, L], [20.0, 20.0], color="#7d3c98", lw=2.0, linestyle="-.", zorder=5)
    ax_mech.text(L / 2.0, 22.0, "弯曲中性轴 (正应力 $\\sigma_{xx}=0$, 无拉压变形)",
                 color="#7d3c98", fontsize=9, fontweight="bold", ha="center")

    # Critical Failure Zone (Root Top Tension Fiber)
    crack_zone = Rectangle((0, 32), 25, 8, facecolor="crimson", alpha=0.35, edgecolor="crimson", lw=1.8, linestyle="--", zorder=6)
    ax_mech.add_patch(crack_zone)
    ax_mech.annotate("⚠️ 最容易被踩断的危险区！\n(承受最大弯曲拉应力)\n动态猛踩下位错被锁死，极易在此诱发瞬间脆断",
                     xy=(10, 40), xytext=(35, 75),
                     arrowprops=dict(facecolor="crimson", shrink=0.08, width=1.5, headwidth=6),
                     fontsize=9.5, fontweight="bold", color="crimson",
                     bbox=dict(boxstyle="round,pad=0.3", fc="#fce4e4", ec="crimson", alpha=0.95))

    # Applied Downward Pedal Force F(t) at Free Tip
    f_arrow = FancyArrowPatch((L + dx_z / 2.0, H + dy_z / 2.0 + 35), (L + dx_z / 2.0, H + dy_z / 2.0 + 2),
                              arrowstyle="-|>", mutation_scale=22, color="#c0392b", lw=3.2, zorder=7)
    ax_mech.add_patch(f_arrow)
    ax_mech.text(L + dx_z / 2.0 + 4, H + dy_z / 2.0 + 26, "【脚踩受力端】\n驾驶员向下踩踏力 $F(t)$\n(垂直踏板板面加载)",
                 color="#c0392b", fontsize=9.5, fontweight="bold", va="center")

    # Body Deceleration Force (D'Alembert Inertia)
    for y_b in [10, 30]:
        b_arrow = FancyArrowPatch((120, y_b), (170, y_b),
                                  arrowstyle="-|>", mutation_scale=15, color="#27ae60", lw=2.0, zorder=7)
        ax_mech.add_patch(b_arrow)
    ax_mech.text(145, -13, "急刹车身减速时踏板自重惯性力\n$F_{inertia} = 9.24$ N (仅占脚踩力的 0.61%, 极微小)",
                 color="#1e8449", fontsize=9.5, fontweight="bold", ha="center")

    # Dimension Annotations: Length L = 300 mm
    ax_mech.annotate("", xy=(0, -8), xytext=(L, -8),
                     arrowprops=dict(arrowstyle="<->", color="#34495e", lw=1.5))
    ax_mech.text(L / 2.0, -18, "梁有效力臂长度 L = 300 mm", color="#34495e", fontsize=9.5, fontweight="bold", ha="center")

    # Dimension: Height H = 40 mm
    ax_mech.annotate("", xy=(-6, 0), xytext=(-6, H),
                     arrowprops=dict(arrowstyle="<->", color="#34495e", lw=1.5))
    ax_mech.text(-12, H / 2.0, "抗弯高度\nH = 40 mm", color="#34495e", fontsize=9, fontweight="bold", ha="right", va="center")

    # Dimension: Thickness B = 10 mm
    ax_mech.annotate("", xy=(L + 4, 0), xytext=(L + dx_z + 4, dy_z),
                     arrowprops=dict(arrowstyle="<->", color="#34495e", lw=1.2))
    ax_mech.text(L + dx_z / 2.0 + 12, dy_z / 2.0 - 5, "板厚 B = 10 mm", color="#34495e", fontsize=8.5, fontweight="bold")

    # Mesh Specification Tag
    ax_mech.text(L / 2.0, 93, "【有限元模型】结构钢悬臂梁 (长300×高40×厚10mm, 共960个六面体单元)\n材料: 结构钢 ($E=210$ GPa, $\\nu=0.3$, $\\rho=7850$ kg/m$^3$, 标称屈服 250 MPa)",
                 color="#2c3e50", fontsize=10, fontweight="bold", ha="center",
                 bbox=dict(boxstyle="round,pad=0.4", fc="#eaeded", ec="#bdc3c7", alpha=0.9))

    ax_mech.set_title("(a) 悬臂钢构件（代表踏板臂/支架）物理力学模型与边界条件", fontsize=12, fontweight="bold", pad=15)

    # -------------------------------------------------------------
    # 2. Applied Force Curves (Panel b)
    # -------------------------------------------------------------
    t_ms = np.linspace(0, 200, 500)
    f_c1 = np.full_like(t_ms, 1500.0)
    f_c2 = np.zeros_like(t_ms)
    for i, t in enumerate(t_ms):
        if t <= 50:
            f_c2[i] = 1500.0 * (t / 50.0)
        elif t <= 100:
            f_c2[i] = 1500.0
        else:
            f_c2[i] = 1500.0 * max(0.0, 1.0 - (t - 100.0) / 100.0)

    f_c3 = np.zeros_like(t_ms)
    for i, t in enumerate(t_ms):
        if t < 10:
            f_c3[i] = 2250.0 * (t / 10.0)
        elif t < 25:
            ratio = (t - 10.0) / 15.0
            f_c3[i] = 2250.0 * (1.0 - ratio) + 1500.0 * ratio
        elif t <= 100:
            f_c3[i] = 1500.0
        else:
            f_c3[i] = 1500.0 * max(0.0, 1.0 - (t - 100.0) / 100.0)

    ax_force.plot(t_ms, f_c1, "b--", lw=2.0, label="Case 1: 准静态台架 (QC/T 788 标称 1500N 恒定慢压)")
    ax_force.plot(t_ms, f_c2, color="#27ae60", lw=2.2, label="Case 2: 动态平缓急踩 (50ms 快速踩下)")
    ax_force.plot(t_ms, f_c3, color="#c0392b", lw=2.5, label="Case 3: 猛踩跺脚冲击 (10ms 接触尖峰 2250N)")
    ax_force.axvspan(0, 10, color="crimson", alpha=0.12, label="初瞬接触碰撞区 (0-10ms)")

    ax_force.set_title("(b) 三大对比实验：脚踏力输入随时间变化曲线 F(t)", fontsize=11, fontweight="bold")
    ax_force.set_ylabel("脚踏板端部作用力 (N)", fontsize=10)
    ax_force.set_xlim(0, 200)
    ax_force.set_ylim(0, 2600)
    ax_force.grid(True, linestyle="--", alpha=0.5)
    ax_force.legend(loc="upper right", fontsize=8.5, framealpha=0.95)

    # -------------------------------------------------------------
    # 3. Vehicle Deceleration Acceleration Curves (Panel c)
    # -------------------------------------------------------------
    a_x = np.zeros_like(t_ms)
    for i, t in enumerate(t_ms):
        if t < 30:
            a_x[i] = 0.0
        elif t <= 60:
            a_x[i] = -1.0 * ((t - 30.0) / 30.0)
        elif t <= 120:
            a_x[i] = -1.0
        else:
            a_x[i] = -1.0 * max(0.0, 1.0 - (t - 120.0) / 80.0)

    ax_accel.plot(t_ms, a_x, color="#1e8449", lw=2.5, label="Case 3: 急刹车身减速度 $a_x(t)$ (峰值 -1.0g)")
    ax_accel.axhline(0, color="gray", linestyle=":", lw=1.2)
    ax_accel.text(90, -0.65, "整车极限减速度 $\\approx 1.0g$\n踏板受自身惯性力 $F_i = 9.24$ N (仅占 0.61%)",
                  fontsize=8.5, color="#145a32", fontweight="bold", ha="center")

    ax_accel.set_title("(c) 车辆运动环境：整车急刹减速度时程 a(t)", fontsize=11, fontweight="bold")
    ax_accel.set_xlabel("时间 t (ms)", fontsize=10)
    ax_accel.set_ylabel("车身减速度 (g)", fontsize=10)
    ax_accel.set_xlim(0, 200)
    ax_accel.set_ylim(-1.3, 0.3)
    ax_accel.grid(True, linestyle="--", alpha=0.5)
    ax_accel.legend(loc="lower left", fontsize=8.5, framealpha=0.95)

    plt.tight_layout()
    schematic_path = os.path.join(out_dir, "experimental_setup_schematic.png")
    plt.savefig(schematic_path)
    plt.close()
    print(f"Saved Publication Schematic Diagram: {schematic_path}")


if __name__ == "__main__":
    plot_experimental_schematic()
