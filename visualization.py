"""
Publication-Quality Visualization Suite for Brake Pedal Structural Dynamic Reliability.
Generates:
  - Fig 1: Load Profiles & Working Condition Time Histories
  - Fig 2: Transient Dynamic Response Curves (Displacement & Stress)
  - Fig 3: Max Stress Comparison Bar Chart & Safety Margin
  - Fig 4: Dynamic Amplification Factor (DAF) & Rate Sensitivity Curve
  - Fig 5: Force Flow & Vehicle Weight Decoupling Breakdown (2t vs 3t)
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Set publication-grade styling and Chinese font
mpl.rcParams["font.family"] = ["Noto Sans CJK SC", "DejaVu Sans", "sans-serif"]
mpl.rcParams["axes.unicode_minus"] = False
mpl.rcParams["figure.dpi"] = 300
mpl.rcParams["savefig.dpi"] = 300
mpl.rcParams["font.size"] = 11
mpl.rcParams["axes.labelsize"] = 12
mpl.rcParams["axes.titlesize"] = 13
mpl.rcParams["xtick.labelsize"] = 10
mpl.rcParams["ytick.labelsize"] = 10
mpl.rcParams["legend.fontsize"] = 10
mpl.rcParams["lines.linewidth"] = 2.0
mpl.rcParams["grid.alpha"] = 0.4
mpl.rcParams["grid.linestyle"] = "--"


def load_data(json_path: str = "results/simulation_results.json"):
    with open(json_path, "r", encoding="utf-8") as f:
        return json.load(f)


def plot_fig1_load_profiles(data, out_dir: str = "figures"):
    """Fig 1: Load profiles for static, transient ramp, panic stamping, and accelerations."""
    os.makedirs(out_dir, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.2))

    t_d1 = np.array(data["Case_D1"]["time"]) * 1000.0  # ms
    f_d1 = np.array(data["Case_D1"]["force"])
    f_imp = np.array(data["Case_D1_Impact"]["force"])

    # Subplot 1: Pedal Force Input
    ax1.plot([0, 200], [1500, 1500], "b--", label="Case S1: 1500N 静态标准 (QC/T 788)", alpha=0.8)
    ax1.plot([0, 200], [3000, 3000], "m-.", label="Case S2: 3000N 静态加严 (2倍冗余)", alpha=0.8)
    ax1.plot(t_d1, f_d1, color="#1f77b4", lw=2.5, label="Case D1: 快速踩踏 (50ms 上升沿)")
    ax1.plot(t_d1, f_imp, color="#d62728", lw=2.5, label="Case D1-Imp: 紧急跺脚冲击 (2250N 峰值)")
    ax1.axvspan(0, 10, color="#d62728", alpha=0.1, label="初瞬接触碰撞区 (0-10ms)")

    ax1.set_title("(a) 踏板作用力输入曲线 F(t)", fontweight="bold")
    ax1.set_xlabel("时间 t (ms)")
    ax1.set_ylabel("踏板端部载荷 (N)")
    ax1.set_xlim(0, 200)
    ax1.set_ylim(0, 3400)
    ax1.grid(True)
    ax1.legend(loc="upper right", framealpha=0.9)

    # Subplot 2: Acceleration Environments
    ax_d2 = np.array(data["Case_D2"]["accel_x"]) / 9.80665  # in g
    ay_d3 = np.array(data["Case_D3"]["accel_y"]) / 9.80665  # in g

    ax2.plot(t_d1, ax_d2, color="#2ca02c", lw=2.5, label="Case D2: 车体制动减速度 $a_x$ (峰值 -1.0g)")
    ax2.plot(t_d1, ay_d3, color="#ff7f0e", lw=2.5, label="Case D3: 道路冲击加速度 $a_y$ (峰值 3.0g, 30ms)")
    ax2.axhline(0, color="gray", lw=1.0, linestyle=":")
    ax2.axvspan(50, 80, color="#ff7f0e", alpha=0.12, label="路面冲击激励窗口 (50-80ms)")

    ax2.set_title("(b) 车体基座与路面动态加速度激励 a(t)", fontweight="bold")
    ax2.set_xlabel("时间 t (ms)")
    ax2.set_ylabel("加速度 (g)")
    ax2.set_xlim(0, 200)
    ax2.set_ylim(-1.5, 3.6)
    ax2.grid(True)
    ax2.legend(loc="upper right", framealpha=0.9)

    plt.tight_layout()
    save_path = os.path.join(out_dir, "fig1_load_profiles.png")
    plt.savefig(save_path)
    plt.close()
    print(f"Saved: {save_path}")


def plot_fig2_dynamic_responses(data, out_dir: str = "figures"):
    """Fig 2: Free-end deflection and max von Mises stress histories."""
    os.makedirs(out_dir, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 9.5), sharex=True)

    t_ms = np.array(data["Case_D1"]["time"]) * 1000.0

    # Deflections
    d_s1 = data["Case_S1"]["max_deflection_mm"]
    d_s2 = data["Case_S2"]["max_deflection_mm"]
    d_d1 = data["Case_D1"]["deflection_mm"]
    d_imp = data["Case_D1_Impact"]["deflection_mm"]
    d_d2 = data["Case_D2"]["deflection_mm"]
    d_d3 = data["Case_D3"]["deflection_mm"]

    ax1.axhline(d_s1, color="blue", linestyle="--", alpha=0.7, label=f"Case S1 (1500N 静态): {d_s1:.2f} mm")
    ax1.axhline(d_s2, color="purple", linestyle="-.", alpha=0.7, label=f"Case S2 (3000N 静态): {d_s2:.2f} mm")
    ax1.plot(t_ms, d_d1, color="#1f77b4", label=f"Case D1 (快速平滑急踩): 峰值 {max(d_d1):.2f} mm")
    ax1.plot(t_ms, d_imp, color="#d62728", lw=2.2, label=f"Case D1-Imp (紧急冲击跺脚): 峰值 {max(d_imp):.2f} mm")
    ax1.plot(t_ms, d_d2, color="#2ca02c", linestyle=":", label=f"Case D2 (减速惯性耦合): 峰值 {max(d_d2):.2f} mm")
    ax1.plot(t_ms, d_d3, color="#ff7f0e", linestyle="--", label=f"Case D3 (3g 道路冲击): 峰值 {max(d_d3):.2f} mm")

    ax1.set_title("(a) 踏板自由端垂直挠度响应时程 $u_y(t)$", fontweight="bold")
    ax1.set_ylabel("自由端挠度 (mm)")
    ax1.set_ylim(0, 2.7)
    ax1.grid(True)
    ax1.legend(loc="upper right", framealpha=0.9, ncol=2)

    # Stresses
    s_s1 = data["Case_S1"]["max_stress_MPa"]
    s_s2 = data["Case_S2"]["max_stress_MPa"]
    s_d1 = data["Case_D1"]["stress_MPa"]
    s_imp = data["Case_D1_Impact"]["stress_MPa"]
    s_d2 = data["Case_D2"]["stress_MPa"]
    s_d3 = data["Case_D3"]["stress_MPa"]

    ax2.axhline(s_s1, color="blue", linestyle="--", alpha=0.7, label=f"Case S1 (1500N 静态): {s_s1:.1f} MPa")
    ax2.axhline(s_s2, color="purple", linestyle="-.", alpha=0.7, label=f"Case S2 (3000N 静态): {s_s2:.1f} MPa")
    ax2.axhline(250.0, color="crimson", linewidth=2.0, linestyle="-", label="材料屈服极限 $\\sigma_y = 250$ MPa")
    ax2.axhspan(250.0, 320.0, color="crimson", alpha=0.08, label="屈服/塑性变形区")

    ax2.plot(t_ms, s_d1, color="#1f77b4", label=f"Case D1 (平滑急踩): 峰值 {max(s_d1):.1f} MPa (DAF=1.01)")
    ax2.plot(t_ms, s_imp, color="#d62728", lw=2.2, label=f"Case D1-Imp (冲击跺脚): 峰值 {max(s_imp):.1f} MPa (DAF=1.54)")
    ax2.plot(t_ms, s_d2, color="#2ca02c", linestyle=":", label=f"Case D2 (惯性耦合): 峰值 {max(s_d2):.1f} MPa (DAF=1.01)")
    ax2.plot(t_ms, s_d3, color="#ff7f0e", linestyle="--", label=f"Case D3 (道路冲击): 峰值 {max(s_d3):.1f} MPa (DAF=1.01)")

    ax2.set_title("(b) 结构最大 von Mises 等效应力响应时程 $\\sigma_{vM}(t)$", fontweight="bold")
    ax2.set_xlabel("时间 t (ms)")
    ax2.set_ylabel("最大 von Mises 应力 (MPa)")
    ax2.set_xlim(0, 200)
    ax2.set_ylim(0, 310)
    ax2.grid(True)
    ax2.legend(loc="upper right", framealpha=0.9, ncol=2)

    plt.tight_layout()
    save_path = os.path.join(out_dir, "fig2_dynamic_response_history.png")
    plt.savefig(save_path)
    plt.close()
    print(f"Saved: {save_path}")


def plot_fig3_stress_bar_comparison(data, out_dir: str = "figures"):
    """Fig 3: Bar chart comparing peak stresses against yield limit."""
    os.makedirs(out_dir, exist_ok=True)
    plt.figure(figsize=(11, 6.2))

    cases = [
        "Case S1\n(1500N 静态)",
        "Case D1\n(50ms 平滑急踩)",
        "Case D2\n(制动惯性耦合)",
        "Case D3\n(3g 道路冲击)",
        "Case D1-Imp\n(冲击跺脚 2250N)",
        "Case S2\n(3000N 极限静态)",
        "Case Step\n(突加阶跃上限)",
    ]

    stresses = [
        data["Case_S1"]["max_stress_MPa"],
        data["Case_D1"]["peak_stress_MPa"],
        data["Case_D2"]["peak_stress_MPa"],
        data["Case_D3"]["peak_stress_MPa"],
        data["Case_D1_Impact"]["peak_stress_MPa"],
        data["Case_S2"]["max_stress_MPa"],
        data["Case_Step"]["peak_stress_MPa"],
    ]

    colors = [
        "#4a90e2",  # S1: blue
        "#50e3c2",  # D1: teal
        "#27ae60",  # D2: green
        "#f39c12",  # D3: amber
        "#e67e22",  # D1-Imp: orange
        "#8e44ad",  # S2: purple
        "#c0392b",  # Step: dark red
    ]

    bars = plt.bar(cases, stresses, color=colors, width=0.55, edgecolor="black", linewidth=1.1, alpha=0.9)
    plt.axhline(250.0, color="crimson", linewidth=2.2, linestyle="--", label="材料屈服极限 $\\sigma_y = 250$ MPa")
    plt.axhspan(250.0, 320.0, color="crimson", alpha=0.1, label="进入塑性屈服失效区")

    # Add value annotations
    for bar in bars:
        h = bar.get_height()
        status = " (屈服)" if h >= 250.0 else " (安全)"
        plt.text(
            bar.get_x() + bar.get_width() / 2.0,
            h + 4.0,
            f"{h:.1f} MPa{status}",
            ha="center",
            va="bottom",
            fontsize=9.5,
            fontweight="bold",
        )

    plt.title("各工况结构最大 von Mises 等效应力对比与安全裕度评估", fontsize=14, fontweight="bold", pad=12)
    plt.ylabel("最大 von Mises 应力 (MPa)", fontsize=12)
    plt.ylim(0, 330)
    plt.grid(axis="y", linestyle="--", alpha=0.5)
    plt.legend(loc="upper left", framealpha=0.95)

    plt.tight_layout()
    save_path = os.path.join(out_dir, "fig3_stress_comparison_bar.png")
    plt.savefig(save_path)
    plt.close()
    print(f"Saved: {save_path}")


def plot_fig4_daf_and_sweep(data, out_dir: str = "figures"):
    """Fig 4: Dynamic Amplification Factor bar chart and rise-time sensitivity curve."""
    os.makedirs(out_dir, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    # Subplot 1: DAF across operating cases
    daf_cases = [
        "S1 (基准)",
        "D1 (平滑急踩)",
        "D2 (惯性耦合)",
        "D3 (路面冲击)",
        "D1-Imp (冲击跺脚)",
        "Step (阶跃上限)",
    ]
    daf_values = [
        1.000,
        data["Case_D1"]["daf"],
        data["Case_D2"]["daf"],
        data["Case_D3"]["daf"],
        data["Case_D1_Impact"]["daf"],
        data["Case_Step"]["daf"],
    ]

    colors = ["#3498db", "#1abc9c", "#2ecc71", "#f1c40f", "#e67e22", "#e74c3c"]
    bars = ax1.bar(daf_cases, daf_values, color=colors, width=0.52, edgecolor="black", linewidth=1.0)
    ax1.axhline(1.0, color="gray", linestyle=":", lw=1.2)

    for bar in bars:
        h = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width() / 2.0, h + 0.03, f"{h:.2f}x", ha="center", va="bottom", fontweight="bold")

    ax1.set_title("(a) 各动态工况相对静态标称应力放大系数 (DAF)", fontweight="bold")
    ax1.set_ylabel("动态放大系数 DAF = $\\sigma_{dyn} / \\sigma_{S1}$")
    ax1.set_ylim(0, 2.35)
    ax1.grid(axis="y", linestyle="--", alpha=0.5)

    # Subplot 2: Rise-time sensitivity sweep
    sweep = data["Parametric_Sweep"]
    r_times = [item["rise_time_ms"] for item in sweep]
    dafs = [item["daf"] for item in sweep]

    ax2.plot(r_times, dafs, "o-", color="#8e44ad", lw=2.5, markersize=7, label="FEM 瞬态动力学仿真曲线")
    ax2.axvline(2.69, color="#d62728", linestyle="--", lw=1.8, label="悬臂梁一阶自振周期 $T_1 = 2.69$ ms")
    ax2.axhline(1.0, color="gray", linestyle=":", lw=1.2, label="准静态极限 (DAF = 1.00)")
    ax2.axhline(2.034, color="crimson", linestyle="-.", lw=1.5, label="突加阶跃理论上限 (DAF = 2.03)")

    ax2.annotate(
        "快速冲击放大区\n$T_{rise} < 5 T_1$",
        xy=(2.69, 1.3),
        xytext=(0.8, 1.6),
        arrowprops=dict(facecolor="black", shrink=0.08, width=1, headwidth=6),
        fontweight="bold",
        color="#d62728",
    )
    ax2.annotate(
        "准静态等效区\n$T_{rise} \\gg T_1$\n(平缓踩踏/台架试验)",
        xy=(100, 1.01),
        xytext=(40, 1.35),
        arrowprops=dict(facecolor="black", shrink=0.08, width=1, headwidth=6),
        fontweight="bold",
        color="#2980b9",
    )

    ax2.set_xscale("log")
    ax2.set_title("(b) 踏板加载上升时间与动态放大效应敏感度扫描", fontweight="bold")
    ax2.set_xlabel("载荷上升时间 $T_{rise}$ (ms, 对数坐标)")
    ax2.set_ylabel("动态放大系数 DAF")
    ax2.set_xlim(0.8, 1200)
    ax2.set_ylim(0.95, 2.15)
    ax2.grid(True, which="both", linestyle="--", alpha=0.5)
    ax2.legend(loc="upper right", framealpha=0.9)

    plt.tight_layout()
    save_path = os.path.join(out_dir, "fig4_daf_and_rise_time_sensitivity.png")
    plt.savefig(save_path)
    plt.close()
    print(f"Saved: {save_path}")


def plot_fig5_mechanics_decoupling(data, out_dir: str = "figures"):
    """Fig 5: Mechanics force breakdown & vehicle weight decoupling comparison."""
    os.makedirs(out_dir, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    # Subplot 1: Force comparison on pedal component
    categories = ["2.0 吨整车\n(1.0g 制动)", "3.0 吨整车\n(1.0g 制动)"]
    foot_force = [1500.0, 1500.0]
    inertia_force = [9.24, 9.24]  # m * a = 0.942 * 9.81 N

    x = np.arange(len(categories))
    width = 0.35

    b1 = ax1.bar(x - width / 2, foot_force, width, label="驾驶员脚踏力输入", color="#3498db", edgecolor="black")
    b2 = ax1.bar(x + width / 2, inertia_force, width, label="踏板自身构件惯性力 ($m_{pedal} \\cdot a$)", color="#e74c3c", edgecolor="black")

    for bar in b1:
        ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 30, f"{bar.get_height():.0f} N\n(99.39%)", ha="center", fontsize=9.5, fontweight="bold")
    for bar in b2:
        ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 30, f"{bar.get_height():.1f} N\n(0.61%)", ha="center", fontsize=9.5, fontweight="bold", color="#c0392b")

    ax1.set_title("(a) 踏板构件直接受载定量分解 (2吨 vs 3吨)", fontweight="bold")
    ax1.set_ylabel("受力大小 (N)")
    ax1.set_ylim(0, 1800)
    ax1.set_xticks(x)
    ax1.set_xticklabels(categories)
    ax1.grid(axis="y", linestyle="--", alpha=0.5)
    ax1.legend(loc="upper right", framealpha=0.9)

    # Subplot 2: Force transmission path explanation
    labels = ["地面制动附着力\n(闭环作用于轮胎-悬架-车身)", "踏板承受力\n(脚踏力 + 踏板自身微小惯性)"]
    forces_2t = [19.62, 1.509]  # kN
    forces_3t = [29.43, 1.509]  # kN

    x2 = np.arange(len(labels))
    b3 = ax2.bar(x2 - width / 2, forces_2t, width, label="2.0 吨车型 (地面附着力 19.6 kN)", color="#2ecc71", edgecolor="black")
    b4 = ax2.bar(x2 + width / 2, forces_3t, width, label="3.0 吨车型 (地面附着力 29.4 kN)", color="#9b59b6", edgecolor="black")

    ax2.text(x2[0] - width / 2, forces_2t[0] + 0.5, f"{forces_2t[0]:.1f} kN", ha="center", fontweight="bold")
    ax2.text(x2[0] + width / 2, forces_3t[0] + 0.5, f"{forces_3t[0]:.1f} kN\n(+50% 底盘附着力)", ha="center", fontweight="bold", color="#6c3483")
    ax2.text(x2[1] - width / 2, forces_2t[1] + 0.5, f"{forces_2t[1]:.3f} kN\n(几乎完全相同)", ha="center", fontweight="bold")
    ax2.text(x2[1] + width / 2, forces_3t[1] + 0.5, f"{forces_3t[1]:.3f} kN\n(几乎完全相同)", ha="center", fontweight="bold")

    ax2.set_title("(b) 底盘力流与踏板力流的隔离解耦对比", fontweight="bold")
    ax2.set_ylabel("总受载力 (kN)")
    ax2.set_ylim(0, 34)
    ax2.set_xticks(x2)
    ax2.set_xticklabels(labels)
    ax2.grid(axis="y", linestyle="--", alpha=0.5)
    ax2.legend(loc="upper right", framealpha=0.9)

    plt.tight_layout()
    save_path = os.path.join(out_dir, "fig5_force_flow_and_inertia_breakdown.png")
    plt.savefig(save_path)
    plt.close()
    print(f"Saved: {save_path}")


def main():
    data = load_data("results/simulation_results.json")
    out_dir = "figures"
    print("Generating publication-quality figures...")
    plot_fig1_load_profiles(data, out_dir)
    plot_fig2_dynamic_responses(data, out_dir)
    plot_fig3_stress_bar_comparison(data, out_dir)
    plot_fig4_daf_and_sweep(data, out_dir)
    plot_fig5_mechanics_decoupling(data, out_dir)
    print("All figures successfully created in figures/ directory!")


if __name__ == "__main__":
    main()
