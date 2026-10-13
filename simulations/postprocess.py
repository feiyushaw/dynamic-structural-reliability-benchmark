"""
Post-processing script for Cantilever Beam FEM Validation Study.
Generates the 3 essential publication-grade figures:
  1. results/load_cases.png
  2. results/displacement.png
  3. results/stress_comparison.png
and outputs the core summary comparison table.
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Publication-grade styling
mpl.rcParams["font.family"] = ["Noto Sans CJK SC", "DejaVu Sans", "sans-serif"]
mpl.rcParams["axes.unicode_minus"] = False
mpl.rcParams["figure.dpi"] = 300
mpl.rcParams["savefig.dpi"] = 300
mpl.rcParams["font.size"] = 10.5
mpl.rcParams["lines.linewidth"] = 2.0


def postprocess_all(results_dir: str = "results"):
    c1_path = os.path.join(results_dir, "case1_static.json")
    c2_path = os.path.join(results_dir, "case2_dynamic.json")
    c3_path = os.path.join(results_dir, "case3_acceleration.json")

    with open(c1_path, "r", encoding="utf-8") as f:
        c1 = json.load(f)
    with open(c2_path, "r", encoding="utf-8") as f:
        c2 = json.load(f)
    with open(c3_path, "r", encoding="utf-8") as f:
        c3 = json.load(f)

    t_ms = np.array(c2["time_s"]) * 1000.0

    # -------------------------------------------------------------
    # Figure 1: Load Cases (Force and Base Acceleration)
    # -------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

    # (a) Force curves
    ax1.plot([0, 200], [1500, 1500], "b--", label="Case 1: Static (1500 N)")
    ax1.plot(t_ms, c2["force_N"], color="#d62728", label="Case 2 & 3: Dynamic Force (50 ms Ramp)")
    ax1.set_title("(a) Applied Tip Force Profile F(t)", fontweight="bold")
    ax1.set_xlabel("Time t (ms)")
    ax1.set_ylabel("Force (N)")
    ax1.set_xlim(0, 200)
    ax1.set_ylim(0, 1750)
    ax1.grid(True, linestyle="--", alpha=0.5)
    ax1.legend(loc="lower right")

    # (b) Base acceleration curve
    a_x_g = np.array(c3["accel_x_mps2"]) / 9.80665
    ax2.plot(t_ms, a_x_g, color="#2ca02c", label="Case 3: Base Deceleration $a_x(t)$ (Peak -1.0g)")
    ax2.axhline(0, color="gray", linestyle=":", lw=1.2)
    ax2.set_title("(b) Vehicle Base Deceleration Profile", fontweight="bold")
    ax2.set_xlabel("Time t (ms)")
    ax2.set_ylabel("Base Acceleration (g)")
    ax2.set_xlim(0, 200)
    ax2.set_ylim(-1.3, 0.4)
    ax2.grid(True, linestyle="--", alpha=0.5)
    ax2.legend(loc="lower left")

    plt.tight_layout()
    fig1_path = os.path.join(results_dir, "load_cases.png")
    plt.savefig(fig1_path)
    plt.close()
    print(f"Saved: {fig1_path}")

    # -------------------------------------------------------------
    # Figure 2: Free-end Displacement Comparison
    # -------------------------------------------------------------
    plt.figure(figsize=(8, 4.8))
    plt.axhline(c1["max_deflection_mm"], color="blue", linestyle="--", label=f"Case 1 (Static): {c1['max_deflection_mm']:.3f} mm")
    plt.plot(t_ms, c2["deflection_mm"], color="#d62728", label=f"Case 2 (Dynamic Force): Peak {c2['peak_deflection_mm']:.3f} mm")
    plt.plot(t_ms, c3["deflection_mm"], color="#2ca02c", linestyle=":", label=f"Case 3 (Dynamic + 1g Accel): Peak {c3['peak_deflection_mm']:.3f} mm")

    plt.title("Free-End Deflection Response Comparison $u_y(t)$", fontweight="bold")
    plt.xlabel("Time t (ms)")
    plt.ylabel("Tip Deflection (mm)")
    plt.xlim(0, 200)
    plt.ylim(0, 1.45)
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.legend(loc="upper right")

    plt.tight_layout()
    fig2_path = os.path.join(results_dir, "displacement.png")
    plt.savefig(fig2_path)
    plt.close()
    print(f"Saved: {fig2_path}")

    # -------------------------------------------------------------
    # Figure 3: Peak Stress Comparison & Amplification
    # -------------------------------------------------------------
    plt.figure(figsize=(7.5, 5.0))
    cases = ["Case 1\nStatic (1500N)", "Case 2\nDynamic Force", "Case 3\nDynamic + 1g Accel"]
    stresses = [c1["max_stress_MPa"], c2["peak_stress_MPa"], c3["peak_stress_MPa"]]
    colors = ["#3498db", "#e74c3c", "#2ecc71"]

    bars = plt.bar(cases, stresses, color=colors, width=0.48, edgecolor="black", linewidth=1.1)
    plt.axhline(250.0, color="crimson", linestyle="--", lw=1.8, label="Yield Limit $\\sigma_y = 250$ MPa")

    for bar in bars:
        h = bar.get_height()
        amp = h / c1["max_stress_MPa"]
        plt.text(bar.get_x() + bar.get_width() / 2.0, h + 3.0, f"{h:.2f} MPa\n({amp:.3f}x)", ha="center", fontsize=9.5, fontweight="bold")

    plt.title("Peak von Mises Stress Comparison", fontweight="bold")
    plt.ylabel("Peak Stress (MPa)")
    plt.ylim(0, 280)
    plt.grid(axis="y", linestyle="--", alpha=0.5)
    plt.legend(loc="upper left")

    plt.tight_layout()
    fig3_path = os.path.join(results_dir, "stress_comparison.png")
    plt.savefig(fig3_path)
    plt.close()
    print(f"Saved: {fig3_path}")

    # -------------------------------------------------------------
    # Print Core Summary Table
    # -------------------------------------------------------------
    print("\n" + "=" * 72)
    print("                FEM EXPERIMENTAL VALIDATION SUMMARY TABLE                ")
    print("=" * 72)
    print(f"{'Case':<20} | {'Load Description':<24} | {'Peak Disp':<10} | {'Peak Stress':<12}")
    print("-" * 72)
    print(f"{'Case 1: Static':<20} | {'1500 N constant':<24} | {c1['max_deflection_mm']:>6.3f} mm  | {c1['max_stress_MPa']:>8.2f} MPa")
    print(f"{'Case 2: Dynamic':<20} | {'1500 N (50ms ramp)':<24} | {c2['peak_deflection_mm']:>6.3f} mm  | {c2['peak_stress_MPa']:>8.2f} MPa")
    print(f"{'Case 3: Dynamic+1g':<20} | {'Force + 1g Base Decel':<24} | {c3['peak_deflection_mm']:>6.3f} mm  | {c3['peak_stress_MPa']:>8.2f} MPa")
    print("=" * 72)
    print("> The static test underestimates structural response under transient vehicle operating conditions.\n")


if __name__ == "__main__":
    postprocess_all()
