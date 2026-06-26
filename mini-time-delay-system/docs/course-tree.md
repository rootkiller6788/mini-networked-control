# Course Tree — mini-time-delay-system

Prerequisite dependency tree for knowledge acquisition.

## Prerequisites (Dependencies)

```
mini-time-delay-system
├── mini-general-system-theory (L1-L2)
│   └── System concepts, state-space models
├── mini-system-dynamics (L1-L3)
│   └── ODE integration, RK4 method
├── mini-optimal-control (L4-L5)
│   └── LQR, LQG, Kalman filtering
├── mini-robust-control (L4)
│   └── Stability margins, Nyquist criterion
└── mini-networked-control (L6-L7)
    └── NCS architecture, packet loss
```

## Internal Dependency Tree

```
L1 Definitions ──────────────────────────────────────────┐
  ├── Delay types, DDE types, stability classes           │
  └── Characteristic quasipolynomial                      │
      ↓                                                    │
L2 Core Concepts ──────────────────────────────────────┐  │
  ├── State norm, delay rate, spectral abscissa         │  │
  └── Stability classification                          │  │
      ↓                                                 │  │
L3 Mathematical Structures ←──────────── L1, L2 ───────┘  │
  ├── Matrix measure μ₂(A)                                 │
  ├── Padé approximation (e^{-τs} → rational)              │
  ├── Chebyshev spectral discretization                    │
  └── Routh-Hurwitz via Rekasius substitution              │
      ↓                                                    │
L4 Fundamental Laws ←──────────── L1-L3 ──────────────────┘
  ├── Lyapunov-Krasovskii Theorem
  │   └── LK functional evaluation + derivative
  ├── Razumikhin Theorem
  │   └── History condition check
  ├── Halanay Inequality
  │   └── Exponential decay bound
  ├── Walton-Marshall RHP Counting
  │   └── Substitution → polynomial → Routh array
  └── Delay Margin Formula
      └── ω² = a_d² - a², τ* = φ/ω
      ↓
L5 Algorithms ←──────────── L1-L4 ─────────────────────────
  ├── Method of Steps (DDE → ODE on intervals)
  ├── RK4/RKF45 DDE Integration
  ├── Smith Predictor Implementation
  ├── LMI Stability Checks
  ├── Timestamp-Based LQG
  ├── TCP/AQM Fluid Model
  └── Predictor Feedback (Krstic backstepping)
      ↓
L6 Canonical Problems ←──────────── L1-L5 ─────────────────
  ├── FOPDT stability regions
  ├── NCS benchmark (delay + loss)
  ├── Smith predictor vs PI comparison
  └── TCP Reno congestion control
      ↓
L7 Applications ←──────────── L1-L6 ───────────────────────
  ├── Boeing 747 roll control (aerospace)
  ├── Mars rover teleoperation (space)
  ├── Smart grid WAMS (energy)
  ├── Automotive brake-by-wire (transportation)
  ├── Chemical process control (process industry)
  └── Telesurgery haptics (medical)
      ↓
L8 Advanced Topics ←──────────── L1-L7 ────────────────────
  ├── Stochastic delay systems (mean-square stability)
  ├── Time-varying delay (MATI, piecewise LK)
  ├── Markovian jump (LQG with random loss)
  └── Event-triggered + delay
      ↓
L9 Research Frontiers ←──────────── All ───────────────────
  ├── Cyber-physical delay attacks
  ├── Multi-agent consensus with heterogeneous delays
  ├── Learning-based delay compensation
  └── Quantum networked control
```

## Study Path (Suggested Order)

1. **Week 1**: L1 — Delay types, DDE classification
2. **Week 2**: L2 — Stability concepts, characteristic equation
3. **Week 3**: L3 — Padé, matrix measure, spectral discretization
4. **Week 4**: L4 — LK theorem (2 lectures), Razumikhin (1 lecture)
5. **Week 5**: L4 — Delay margin, Halanay, RHP counting
6. **Week 6**: L5 — DDE solvers, method of steps
7. **Week 7**: L5 — Smith predictor, LMI methods
8. **Week 8**: L5 — NCS, timestamp LQG, TCP/AQM
9. **Week 9**: L6 — Canonical problems, end-to-end examples
10. **Week 10**: L7-L8 — Applications + advanced topics
11. **Week 11**: L9 — Research frontiers seminar

Total: ~11 weeks for complete module mastery.
