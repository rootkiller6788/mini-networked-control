# Knowledge Graph — Bandwidth-Limited Control

## L1: Definitions
- **Bandwidth (B)**: Frequency range available for communication (Hz)
- **Bit Rate (R)**: Information transmitted per unit time (bps)
- **Channel Capacity (C)**: C = B·log₂(1+SNR) (Shannon, 1948)
- **Quantization**: Mapping continuous values to discrete levels
- **Quantization Error**: ε = x - Q(x), bounded by ±Δ/2
- **Data Rate Theorem**: Minimum bit rate for stabilizability
- **Packet**: Unit of communication with header, payload, checksum
- **MATI**: Maximum Allowable Transfer Interval
- **Send-on-Delta**: Event-triggered transmission on significant change
- **Zoom Quantizer**: Dynamic range adjustment for quantization

## L2: Core Concepts
- **Bandwidth-limited stabilization**: Control under finite bit rate
- **Data Rate Theorem (DRT)**: R_min = Σ log₂|λᵢ| for |λᵢ|>1
- **Quantization-induced instability**: Error acts as disturbance
- **Separation failure**: Estimation/control separation fails under rate constraints
- **Encoding horizon**: Prediction to reduce required bit rate
- **Minimum attention**: Minimize sensor access rate
- **Event-triggered communication**: Transmit only on significant events
- **Bandwidth scheduling**: Allocate shared bandwidth among control loops
- **Try-Once-Discard**: Optimal contention resolution for NCS
- **Water-filling**: Adaptive bandwidth allocation

## L3: Mathematical Structures
- **Channel**: Shannon-Hartley capacity C = B·log₂(1+SNR)
- **Plant model**: ẋ = Ax + Bu, with quantized measurements
- **Quantizer**: Q: R → {q₁, ..., q_N} with N discrete levels
- **Uniform quantization**: Δ = 2U/N, max error ≤ Δ/2
- **Logarithmic quantization**: Geometric spacing, optimal for stabilization
- **Lyapunov equation**: A'P + PA = -Q (for stability analysis)
- **ARE**: A'P + PA - PBR^{-1}B'P + Q = 0 (LQR optimal control)
- **Rate-distortion**: R(D) = 0.5·log₂(σ²/D) + Σλ/ln2

## L4: Fundamental Theorems
- **Shannon-Hartley Theorem**: C = B·log₂(1+SNR) > 0
- **Data Rate Theorem (Nair-Evans, 2004)**: R > Σ log₂|λᵢ| for |λᵢ|>1
- **Quantizer Error Bound**: |ε| ≤ Δ/2 = U/N
- **Zoom Quantizer Convergence**: Geometric convergence with ρ_in > 1
- **Huffman Optimality**: H ≤ L < H+1
- **Event-Triggered Min Interval**: τ_min > 0 (no Zeno)
- **TOD Optimality**: Maximizes stability region
- **MATI Theorem (Walsh et al., 2002)**: τ_MATI ≤ 1/(16·||A_c||·√n·(1+α)/(1-α))

## L5: Algorithms
- **QR eigenvalue algorithm** (Francis, 1961)
- **Newton-Kleinman ARE iteration** (Kleinman, 1968)
- **Lloyd-Max quantizer design** (Lloyd, 1957; Max, 1960)
- **LBG vector quantization** (Linde-Buzo-Gray, 1980)
- **Huffman entropy coding** (Huffman, 1952)
- **Arithmetic coding** (Witten-Neal-Cleary, 1987)
- **Zoom quantizer with dynamic range**
- **Delta modulation** (1-bit encoding)
- **Predictive encoding for control**
- **Run-length encoding for sparse control**
- **LQR with quantization penalty analysis**
- **Circle Criterion for quantization robustness**
- **Projected gradient MPC with L1 proxy**

## L6: Canonical Problems
- **Double integrator stabilization under bit-rate constraint**
- **DC motor servo with event-triggered CAN bus**
- **Smart grid multi-loop SCADA scheduling**
- **Unstable pole stabilization with quantized feedback**
- **Bandwidth allocation among heterogeneous control loops**
- **TOD contention resolution for NCS**
- **Adaptive bandwidth scheduling under varying conditions**

## L7: Applications
- **Smart grid SCADA** (IEC 61850) — multi-loop bandwidth scheduling
- **DC motor control** (CAN bus, ISO 15745) — event-triggered comms
- **Aerospace telemetry** (Boeing 787 AFDX) — rate-constrained control
- **Tesla vehicle control** — event-triggered CAN messaging
- **NASA Deep Space Network** — multi-spacecraft bandwidth allocation
- **Fukushima monitoring** — post-accident degraded comms

## L8: Advanced Topics
- **Adaptive bandwidth scheduling** with water-filling
- **Self-triggered control** — pre-computed next transmission time
- **Lebesgue sampling** — level-crossing detection
- **TOD protocol with weighted error**
- **Rate-distortion optimization for control**
- **Robust control via Circle Criterion** for quantization

## L9: Research Frontiers
- **Semantic compression** for networked control
- **Learning-based bandwidth allocation**
- **Information-theoretic control limits**
- **Quantum control with limited classical communication**
- **Age-of-Information (AoI) in control loops**
