"""
Figure 2: Three Cases Systematic Comparison Suite.
Generates results/fig2_three_cases_comparison.png:
  - (a) Tip deflection histories u(t)
  - (b) Max von Mises stress histories sigma(t)
  - (c) Work stress vs Failure limit comparison
  - (d) Structural safety factor SF comparison
Clean scientific style with zero text clutter.
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams["font.family"] = ["Noto Sans CJK SC", "DejaVu Sans", "sans-serif"]
mpl.rcParams["axes.unicode_minus"] = False
mpl.rcParams["figure.dpi"] = 300
mpl.rcParams["savefig.dpi"] = 300
mpl.rcParams["font.size"] = 9.5
mpl.rcParams["lines.linewidth"] = 1.8


def postprocess_all(results_dir: str = "results"):
    c1_path = os.path.join(results_dir, "case1_static.json")
    c2_path = os.path.join(results_dir, "case2_dynamic.json")
    c3_path = os.path.join(results_dir, "case3_acceleration.json")

    with open(c1_path, "r", encoding="utf-8") as f: c1 = json.load(f)
    with open(c2_path, "r", encoding="utf-8") as f: c2 = json.load(f)
    with open(c3_path, "r", encoding="utf-8") as f: c3 = json.load(f)

    t_ms = np.array(c2["time_s"]) * 1000.0
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5))

    # -------------------------------------------------------------
    # (a) Deflection History
    # -------------------------------------------------------------
    ax1 = axes[0, 0]
    ax1.axhline(c1["max_deflection_mm"], color="blue", linestyle="--", lw=1.5, label=f"Case 1 (静态): {c1['max_deflection_mm']:.2f} mm")
    ax1.plot(t_ms, c2["deflection_mm"], color="#27ae60", lw=1.8, label=f"Case 2 (急踩): 峰值 {c2['peak_deflection_mm']:.2f} mm")
    ax1.plot(t_ms, c3["deflection_mm"], color="#c0392b", lw=2.0, label=f"Case 3 (冲击): 峰值 {c3['peak_deflection_mm']:.2f} mm")
    ax1.set_title("(a) 自由端挠度时程 u(t)", fontsize=10.5, fontweight="bold")
    ax1.set_xlabel("时间 t (ms)")
    ax1.set_ylabel("挠度 (mm)")
    ax1.set_xlim(0, 150)
    ax1.set_ylim(0, 2.2)
    ax1.grid(True, linestyle="--", alpha=0.5)
    ax1.legend(loc="upper right", fontsize=8.5)

    # -------------------------------------------------------------
    # (b) Stress History
    # -------------------------------------------------------------
    ax2 = axes[0, 1]
    ax2.axhline(c1["max_stress_MPa"], color="blue", linestyle="--", lw=1.5, label=f"Case 1 (静态): {c1['max_stress_MPa']:.1f} MPa")
    ax2.plot(t_ms, c2["stress_MPa"], color="#27ae60", lw=1.8, label=f"Case 2 (急踩): 峰值 {c2['peak_stress_MPa']:.1f} MPa")
    ax2.plot(t_ms, c3["stress_MPa"], color="#c0392b", lw=2.0, label=f"Case 3 (冲击): 峰值 {c3['peak_stress_MPa']:.1f} MPa")
    ax2.set_title("(b) 最大等效应力时程 $\\sigma_{vM}(t)$", fontsize=10.5, fontweight="bold")
    ax2.set_xlabel("时间 t (ms)")
    ax2.set_ylabel("应力 (MPa)")
    ax2.set_xlim(0, 150)
    ax2.set_ylim(0, 250)
    ax2.grid(True, linestyle="--", alpha=0.5)
    ax2.legend(loc="upper right", fontsize=8.5)

    # -------------------------------------------------------------
    # (c) Work Stress vs Failure Limit
    # -------------------------------------------------------------
    ax3 = axes[1, 0]
    cases = ["Case 1\n准静态", "Case 2\n动态急踩", "Case 3\n高率冲击"]
    work_s = [c1["max_stress_MPa"], c2["peak_stress_MPa"], c3["peak_stress_MPa"]]
    fail_s = [400.0, 320.0, 245.0]
    x = np.arange(len(cases))
    w = 0.35

    b1 = ax3.bar(x - w/2, work_s, w, label="工作应力 $\\sigma_{work}$", color="#3498db", edgecolor="black", lw=0.8)
    b2 = ax3.bar(x + w/2, fail_s, w, label="破坏极限 $\\sigma_{fail}$", color="#e74c3c", edgecolor="black", lw=0.8)

    for bar in b1:
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5, f"{bar.get_height():.1f}", ha="center", fontsize=8.5, fontweight="bold")
    for bar in b2:
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5, f"{bar.get_height():.1f}", ha="center", fontsize=8.5, fontweight="bold", color="#c0392b")

    ax3.set_title("(c) 峰值工作应力 vs 真实破坏极限", fontsize=10.5, fontweight="bold")
    ax3.set_ylabel("应力 (MPa)")
    ax3.set_ylim(0, 470)
    ax3.set_xticks(x)
    ax3.set_xticklabels(cases)
    ax3.grid(axis="y", linestyle="--", alpha=0.5)
    ax3.legend(loc="upper left", fontsize=8.5)

    # -------------------------------------------------------------
    # (d) Safety Factor
    # -------------------------------------------------------------
    ax4 = axes[1, 1]
    sf = [fail_s[i] / work_s[i] for i in range(3)]
    colors_sf = ["#27ae60", "#f39c12", "#c0392b"]
    bars_sf = ax4.bar(cases, sf, color=colors_sf, width=0.45, edgecolor="black", lw=0.8)
    ax4.axhline(1.0, color="crimson", linestyle="--", lw=1.2, label="破坏临界线 (SF=1.0)")

    for bar in bars_sf:
        h = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2, h + 0.06, f"{h:.2f}x", ha="center", fontsize=9, fontweight="bold")

    ax4.set_title("(d) 破坏安全系数 $SF = \\sigma_{fail} / \\sigma_{work}$", fontsize=10.5, fontweight="bold")
    ax4.set_ylabel("安全系数 (倍)")
    ax4.set_ylim(0, 3.3)
    ax4.grid(axis="y", linestyle="--", alpha=0.5)
    ax4.legend(loc="upper right", fontsize=8.5)

    plt.tight_layout()
    fig2_path = os.path.join(results_dir, "fig2_three_cases_comparison.png")
    plt.savefig(fig2_path)
    plt.close()
    print(f"Saved Core Comparison Figure: {fig2_path}")

    # Print Clean Summary Table
    print("\n" + "=" * 70)
    print(f"{'工况编号':<16} | {'等效应变率':<13} | {'工作应力':<10} | {'破坏极限':<10} | {'安全系数'}")
    print("-" * 70)
    print(f"{'Case 1: 准静态':<16} | {c1['strain_rate_s1']:>9.1e} s^-1 | {c1['max_stress_MPa']:>7.1f} MPa | {fail_s[0]:>7.1f} MPa | {sf[0]:>6.2f} x")
    print(f"{'Case 2: 动态急踩':<16} | {c2['peak_strain_rate_s1']:>9.1e} s^-1 | {c2['peak_stress_MPa']:>7.1f} MPa | {fail_s[1]:>7.1f} MPa | {sf[1]:>6.2f} x")
    print(f"{'Case 3: 高率冲击':<16} | {c3['peak_strain_rate_s1']:>9.1e} s^-1 | {c3['peak_stress_MPa']:>7.1f} MPa | {fail_s[2]:>7.1f} MPa | {sf[2]:>6.2f} x")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    postprocess_all()
