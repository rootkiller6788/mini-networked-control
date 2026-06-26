# Mini Networked Control（迷你网络控制系统）

一套**从零构建、零外部依赖的 C 语言实现**，覆盖网络控制系统（NCS）理论——当反馈回路经由通信网络闭合时，控制性能、带宽、量化、延迟与丢包之间会产生根本性的折衷。每个子模块对应 MIT、Stanford 等顶尖大学的课程，涵盖数据速率定理、事件触发控制、一致性协议、信息物理安全等方向。

## 子模块

| 子模块 | 主题 | 核心课程 |
|--------|------|----------|
| [mini-bandwidth-limited-control](mini-bandwidth-limited-control/) | 数据速率定理、Shannon-Hartley 信道容量、编码与量化策略、最小注意力控制、带宽调度 | MIT 6.241J, Stanford AA 203 |
| [mini-cloud-control-system](mini-cloud-control-system/) | 云-边协同、分层控制、网络延迟模型、资源分配（RM/EDF 调度）、Smith 预估器 | MIT 6.824, Stanford CS 244B |
| [mini-consensus-over-network](mini-consensus-over-network/) | 图 Laplacian 矩阵、连续/离散一致性协议、Gossip 算法、有限时间一致性、事件触发一致性、收敛速率分析 | MIT 6.241J, Stanford AA 203 |
| [mini-cyber-physical-security](mini-cyber-physical-security/) | 虚假数据注入（FDI）检测、卡方/CUSUM 检测器、博弈论安全分析、ℓ₀/ℓ₁ 安全状态估计、物理水印 | MIT 6.241J, Stanford AA 273 |
| [mini-event-based-control](mini-event-based-control/) | 事件触发控制（ETC）、自触发控制（STC）、周期事件触发控制（PETC）、ISS-Lyapunov 稳定性、触发条件设计、最小事件间隔 | MIT 6.241J, Stanford AA 203 |
| [mini-packet-loss-control](mini-packet-loss-control/) | 间歇观测 Kalman 滤波、丢包网络 LQG 最优控制（TCP/UDP 类）、Markov 跳变线性系统、分组预测控制、临界丢包概率 | MIT 6.241J, Stanford AA 273 |
| [mini-quantized-control](mini-quantized-control/) | 数据速率定理、均匀/对数/矢量/抖动量化器、Huffman/算术编码、增量调制、差分脉冲编码、游程编码 | MIT 6.241J, Stanford AA 203 |
| [mini-time-delay-system](mini-time-delay-system/) | 延迟微分方程数值求解、Lyapunov-Krasovskii 泛函、延迟系统 Nyquist 判据、Smith 预估器、网络延迟模型（UDP/TCP/CAN/以太网） | MIT 6.241J, Stanford AA 203 |

## 设计哲学

- **零外部依赖** — 纯 C99/C11，仅使用标准库头文件
- **子模块自包含** — 每个子模块拥有独立的 `include/`、`src/`、`CMakeLists.txt` 与冒烟测试
- **理论-代码一一对应** — 每个模块内联引用奠基性论文（Wong & Brockett 1997, Nair & Evans 2004, Tabuada 2007, Schenato et al. 2007 等）
- **NCS 优先视角** — 所有控制设计均显式考虑通信约束：比特率限制、丢包、延迟、量化与安全

## 构建

每个子模块均可独立构建。使用 CMake：

```bash
cd mini-bandwidth-limited-control
mkdir build && cd build
cmake ..
make
./smoke_test
```

需要 **C99 兼容编译器**与 **CMake ≥ 3.14**。

## 项目结构

```
20. mini-networked-control/
├── mini-bandwidth-limited-control/   # 数据速率定理、Shannon 容量、编码、带宽调度
├── mini-cloud-control-system/        # 云-边协同、网络延迟、资源分配
├── mini-consensus-over-network/      # 图 Laplacian、一致性协议、Gossip、有限时间一致性
├── mini-cyber-physical-security/     # FDI 检测、博弈论安全、安全估计、水印
├── mini-event-based-control/         # ETC、STC、PETC、ISS-Lyapunov 稳定性、触发条件
├── mini-packet-loss-control/         # 间歇观测 Kalman、丢包 LQG、预测控制
├── mini-quantized-control/           # 数据速率定理、均匀/对数/矢量/抖动量化器、编码
├── mini-time-delay-system/           # DDE 求解器、Lyapunov-Krasovskii、Smith 预估器、网络延迟
├── .gitignore
├── README.md
└── README-CN.md
```

## 许可证

MIT
