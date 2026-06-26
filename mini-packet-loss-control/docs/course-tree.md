# Course Tree — Packet Loss Control

## Prerequisites

```
Linear Algebra (matrices, eigenvalues, SVD)
    └── Probability Theory (random variables, Markov chains)
            └── Linear Systems Theory (state-space, stability, controllability)
                    ├── Optimal Control (LQR, DARE, dynamic programming)
                    │       └── Packet Loss Control (THIS MODULE)
                    ├── Estimation Theory (Kalman filter, observers)
                    │       └── Packet Loss Control (THIS MODULE)
                    └── Stochastic Processes (Markov chains, MJLS)
                            └── Packet Loss Control (THIS MODULE)
```

## Internal Dependency Graph

```
packet_loss_core.h/.c          (L1: definitions, channel models, PRNG)
    ├── packet_loss_controller.h/.c  (L3-L5: LTI, DARE, TCP/UDP control)
    │       ├── packet_loss_estimator.h/.c   (L5: KF, intermittent KF, IMM)
    │       └── packet_loss_predictor.h/.c   (L5: PPC, hold selector, RG)
    └── packet_loss_analysis.h/.c    (L4: stability, MJLS, critical prob)
            └── packet_loss_predictor.h/.c   (uses stability results)
```

## Knowledge Flow

1. **Channel Models** (L1) → Understand the loss process
2. **LTI System** (L3) → Model the plant dynamics
3. **Stability Analysis** (L4) → Determine tolerable loss rates
4. **Estimation** (L5) → Recover state from intermittent measurements
5. **Control Design** (L5) → Compute control under loss
6. **Compensation** (L5) → Predictive/proactive strategies
7. **Application** (L7) → Real-world deployment considerations
