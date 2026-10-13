# Dynamic Structural Reliability Benchmark for Automotive Brake Components

A minimal 3D finite element benchmark demonstrating why dynamic loading impact rates fundamentally alter material constitutive behavior (Cowper-Symonds strain-rate hardening & dynamic embrittlement), exposing why quasi-static strength qualification fails to ensure dynamic reliability.

---

## 📌 核心研究命题：为什么“动态测试”与“静态测试”是完全不同的两个实验？

在传统汽车行业标准（如中国汽车行业标准 **QC/T 788-2018《汽车踏板装置性能要求及台架试验方法》**）中，制动踏板总成的结构强度检验普遍采用准静态台架试验（以 $20\,\text{mm/min}$ 的极低速率施加 $1500\,\text{N}$ 纵向力，加严检验至 $3000\,\text{N}$）。

然而，行业与舆论往往存在一个重大力学认知误区——**试图把“动态运行”简单等同为“静态力乘上一个安全系数”**（例如认为静态抗住 3000 N 就等于能抵抗动态急刹）。

本项目基于连续介质弹性动力学与高应变率材料冲击力学的第一性原理，通过严谨的三维实体有限元仿真证明：  
**静态测试与动态测试在物理本质上是完全绝缘、机理彻底不同的两个实验！静态加严测试根本无法覆盖真实动态工况的结构破坏风险！**

---

## 🛠️ 一、 实验模型设计与力学边界条件 (Mermaid 流程图)

```mermaid
flowchart LR
    subgraph Support["【固定安装端】X = 0"]
        Fixed["车身前围防火墙 / 踏板支架托架<br><b>全固支边界条件</b><br>ux = uy = uz = 0"]
    end

    subgraph Beam["【受力构件：钢制悬臂梁】长300 × 高40 × 厚10 mm (二维平面应力)"]
        direction TB
        TopZone["⚠️ <b>上表面根部拉伸危险区</b> (Y = 40 mm)<br>承受全场最大弯曲拉应力！<br>动态冲击下最容易诱发解理脆性断裂"]
        Neutral["<b>截面中性轴</b> (Y = 20 mm)<br>弯曲正应力 σ_xx = 0 (无拉压变形)"]
        BottomZone["下表面受压屈服区 (Y = 0 mm)<br>承受纯受压应力"]
        TopZone --- Neutral --- BottomZone
    end

    subgraph Loading["【脚踩加载端】X = 300 mm"]
        FootForce["<b>驾驶员脚踏力 F(t)</b><br>自由端竖直向下施加<br>(垂直踏板受力板面)"]
    end

    subgraph ChassisInertia["【车辆急刹减速运动环境】"]
        Decel["整车紧急制动减速度 ax = -1.0g<br>踏板自身反向惯性力 fb = -ρ·ax<br><b>大小仅 9.24 N (仅占踩踏力的 0.61%)</b>"]
    end

    Fixed === Beam
    Beam === Loading
    ChassisInertia -.->|体积分力耦合| Beam
```

---

## ⚡ 二、 三大对比实验加载条件 (Mermaid 加载时程图)

```mermaid
flowchart TD
    Main["汽车制动踏板 3 大对比测试工况"]

    subgraph Case1["Case 1：慢速台架标准测试 (QC/T 788 基准)"]
        C1_Desc["• 加载方式: 试验机以 20 mm/min 极慢匀速下压<br>• 施加外力: 1500 N 恒定不变 (持续数秒)<br>• 变形速度: 等效应变率 ε̇ ≈ 1.0×10⁻⁴ s⁻¹ (极慢蠕变区)<br>• 材料状态: 内部位错滑移充分，完全塑性延性<br>• 承载极限: 靠塑性变形吸能，厚盾牌 400 MPa"]
    end

    subgraph Case2["Case 2：动态平缓急踩刹车 (50ms 快速踩下)"]
        C2_Desc["• 加载方式: 50 ms 快速爬升至 1500 N 并保持<br>• 施加外力: 快速斜坡输入 F(t)<br>• 变形速度: 等效应变率 ε̇ ≈ 2.2×10⁻² s⁻¹ (过渡区)<br>• 材料状态: 轻度率硬化，流变应力达 305 MPa<br>• 承载极限: 延性保留 69.9%，破坏极限降至 320 MPa"]
    end

    subgraph Case3["Case 3：紧急避险猛踩跺脚冲击 (10ms 冲量尖峰 + 1g减速)"]
        C3_Desc["• 加载方式: 惊恐悬空猛踹踏板，脚掌初速度碰撞<br>• 施加外力: 前 10 ms 冲量尖峰飙升至 2250 N<br>• 变形速度: 等效应变率 ε̇ ≈ 1.6×10⁻¹ s⁻¹ (暴涨 1600 倍！)<br>• 材料状态: 严重动态脆化 (位错锁死，吸能丧失近40%)<br>• 破坏极限: 丧失塑性缓冲，脆断极限暴跌至 245 MPa！"]
    end

    Main --> Case1
    Main --> Case2
    Main --> Case3
```

---

## 🚗 三、 整车力流闭环 vs 踏板受力第一性原理解耦 (Mermaid 拓扑)

```mermaid
flowchart LR
    subgraph GroundChassis["【底盘封闭力流闭环】平衡整车动能，与踏板完全绝缘"]
        direction TB
        Friction["地面制动附着力<br>F = μ·M·g<br>(2吨车: 19.6 kN | 3吨车: 29.4 kN)"] --> Tire["轮胎接地印痕"]
        Tire --> Caliper["轮端刹车卡钳 / 制动盘"]
        Caliper --> Suspension["转向节与悬架摆臂"]
        Suspension --> Chassis["副车架与白车身纵梁"]
        Chassis -->|克服整车动能，产生 -1.0g 减速度| Friction
    end

    subgraph Barrier["【物理隔离】"]
        Cross["❌ 绝对不倒灌<br>无力学传递通道"]
    end

    subgraph PedalCabin["【驾驶舱踏板受载】与车重基本无关"]
        direction TB
        DriverFoot["驾驶员脚掌踩踏力<br><b>F = 1500 N</b><br>(占总受力 99.39%)"]
        Component["踏板臂构件 (自重仅 0.94 kg)"]
        PedalInertia["整车急刹产生的踏板自重惯性力<br>F_inertia = m·a = <b>9.24 N</b><br>(仅占总受力 0.61%)"]
        DriverFoot --> Component
        PedalInertia --> Component
    end

    GroundChassis -.-> Barrier -.-> PedalCabin
```

---

## 📊 四、 理论解析解 vs 有限元数值解 严格对齐

| 物理量名称 | 理论解析解 (Analytical Closed-form) | FEM 数值解 (Numerical Solution) | 相对误差 |
| :--- | :---: | :---: | :---: |
| **一阶自振频率 $f_1$** | **$371.3\,\text{Hz}$** | **$371.3\,\text{Hz}$** | $0.00\%$ |
| **自振基频周期 $T_1$** | **$2.69\,\text{ms}$** | **$2.69\,\text{ms}$** | $0.00\%$ |
| **静态端部挠度 $w_{max}$** | **$1.205\,\text{mm}$** ($w = \frac{FL^3}{3EI}$) | **$1.204\,\text{mm}$** | $< 0.1\%$ |
| **根部最大正应力 $\sigma_{max}$** | **$168.75\,\text{MPa}$** ($\sigma = \frac{6FL}{BH^2}$) | **$168.75\,\text{MPa}$** (外推节点值) | $< 0.1\%$ |
| **阶跃冲击放大系数 $DAF$** | **$2.000$** (无阻尼理论极限) | **$1.840$** (含 3% 阻尼数值解) | 物理自洽 |

---

## 📈 五、 三大工况系统性对比与图表成果展示 (`results/`)

### 图 1：二维力学模型与三大加载工况时程波形
![Figure 1: Model & Loading](results/fig1_model_and_loading.png)
*(a) 纯净二维平面应力悬臂梁力学简图与解析解标定；(b) 三大工况载荷时间历程 $F(t)$ 输入。*

---

### 图 2：三大工况系统性全景对比图 (时程响应、破坏极限与安全系数断崖)
![Figure 2: Three Cases Comparison](results/fig2_three_cases_comparison.png)
*(a) 自由端挠度时程对比 $u(t)$；(b) 最大等效应力时程对比 $\sigma_{vM}(t)$；(c) 实际工作应力 vs 真实破坏极限对比柱状图；(d) 破坏安全系数 $SF$ 断崖式缩水对比柱状图。*

```text
======================================================================
工况编号             | 等效应变率         | 工作应力       | 破坏极限       | 安全系数
----------------------------------------------------------------------
Case 1: 准静态      |   1.0e-04 s^-1 |   143.6 MPa |   400.0 MPa |   2.78 x
Case 2: 动态急踩     |   2.2e-02 s^-1 |   144.7 MPa |   320.0 MPa |   2.21 x
Case 3: 高率冲击     |   1.6e-01 s^-1 |   220.8 MPa |   245.0 MPa |   1.11 x
======================================================================
```

---

### 图 3：二维空间应力场分布与破坏危险度云图
![Figure 3: Spatial Stress Contours](results/fig3_spatial_stress_contours.png)
*(a) Case 1 静态应力场 $\sigma_{vM}$ (峰值 143.6 MPa)；(b) Case 3 动态冲击应力场 $\sigma_{vM}$ (峰值 220.8 MPa)；(c) Case 1 破坏危险度 $\eta \le 0.36$ (充分安全)；(d) Case 3 破坏危险度 $\eta = 0.90$ (根部受拉区濒临脆断极限 1.0)。*

---

## 🚀 六、 一键运行与复现

```bash
# 运行全部算例、渲染全部 3 张高清图表并导出 ParaView 3D 模型
python run_all.py
```

---

## 📄 License
MIT License
