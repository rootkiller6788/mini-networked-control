# Course Tree — Bandwidth-Limited Control

## Prerequisites
```
Linear Algebra → Linear Systems Theory → Bandwidth-Limited Control
                     ↓
              Optimal Control (LQR) ────┘
                     ↓
Information Theory (Shannon) ──────────┘
                     ↓
Networked Control Systems ─────────────┘
```

## Dependency Graph (Internal)
```
blc_core.h ──────► blc_core.c (system lifecycle, quantizer, packet, simulation)
     │
     ├──► blc_datarate.h ──► blc_datarate.c (eigenvalues, zoom Q, delta mod, predictive)
     │         │
     │         ├──► blc_control.h ──► blc_control.c (LQR, quantized FB, MPC, robust)
     │         │
     │         ├──► blc_event.h ────► blc_event.c (SOD, Lebesgue, ETF, self-triggered)
     │         │
     │         └──► blc_scheduling.h ► blc_scheduling.c (TDMA, alloc, adaptive, TOD)
     │
     └──► blc_encoding.h ─► blc_encoding.c (log Q, Lloyd-Max, VQ, Huffman, arithmetic)
```

## Knowledge Prerequisites
1. **Linear Algebra**: Eigenvalues, SVD, matrix exponential, Lyapunov equation
2. **Control Theory**: State-space models, LQR, ARE, stability
3. **Information Theory**: Entropy, channel capacity, rate-distortion
4. **Communication**: Quantization, encoding, packets, scheduling
5. **Real-Time Systems**: TDMA, MATI, priority scheduling

## Postrequisites (What depends on this module)
- mini-cloud-control-system (bandwidth-aware cloud control)
- mini-cyber-physical-security (secure bandwidth allocation)
- mini-consensus-over-network (bandwidth-limited consensus)
- mini-packet-loss-control (related: lossy channels)
- mini-event-based-control (event-triggering)
- mini-quantized-control (quantization)
