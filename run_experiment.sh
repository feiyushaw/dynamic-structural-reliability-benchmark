#!/usr/bin/env bash
# ==============================================================================
# 制动踏板与支架结构动力学可靠性有限元仿真验证方案
# 一键运行脚本 (DOLFINx 0.11 / Python 瞬态有限元求解器)
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

echo "=============================================================================="
echo " 🚗 制动踏板结构动态可靠性 FEM 验证试验启动"
echo " 当前目录: $(pwd)"
echo " 时间: $(date)"
echo "=============================================================================="

# 1. 运行有限元仿真求解器 (Case S1, S2, D1, D1-Imp, D2, D3, Step & Sweep)
echo -e "\n>>> [Step 1/2] 执行 Python 有限元计算与动力学时程仿真..."
python experiment_runner.py

# 2. 生成出版级可视化图表
echo -e "\n>>> [Step 2/2] 生成高清学术与工程分析图表 (figures/)..."
python visualization.py

echo -e "\n=============================================================================="
echo " 🎉 仿真计算与图表生成全部完成！"
echo " 仿真结果数据保存在: results/simulation_results.json"
echo " 出版级分析图表保存在: figures/"
echo "   - figures/fig1_load_profiles.png"
echo "   - figures/fig2_dynamic_response_history.png"
echo "   - figures/fig3_stress_comparison_bar.png"
echo "   - figures/fig4_daf_and_rise_time_sensitivity.png"
echo "   - figures/fig5_force_flow_and_inertia_breakdown.png"
echo "=============================================================================="
