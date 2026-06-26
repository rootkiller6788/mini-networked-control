# Course Tree — mini-quantized-control

## Prerequisites (Knowledge Dependencies)

```
Linear Systems Theory
    ↓
Nonlinear Systems (sector boundedness)
    ↓
Information Theory (entropy, channel capacity)
    ↓
┌──────────────────────────────────────┐
│     Quantized Control Theory         │
│  (this module)                       │
└──────────────────────────────────────┘
    ↓                    ↓
Networked Control    Robust Control
(communication      (sector-bound
 constraints)        approach)
```

## Internal Dependencies

```
L1: Definitions
  ↓
L2: Core Concepts (quantization error, data rate, encoding)
  ↓
L3: Math Structures (quantizer types, channel models)
  ↓
L4: Fundamental Laws (Data Rate Theorem, Sector Bound)
  ↓
L5: Algorithms (Lloyd-Max, LBG, Huffman, water-filling)
  ↓
L6: Canonical Problems (quantized LQR, min rate computation)
  ↓
L7: Applications (NCS, embedded control, sensor networks)
  ↓
L8: Advanced Topics (anytime capacity, adaptive quantization)
  ↓
L9: Research Frontiers (info-theoretic limits)
```

## Relation to Sibling Modules

- `mini-bandwidth-limited-control`: Data rate theorem directly applicable
- `mini-packet-loss-control`: Erasure channel model shared
- `mini-event-based-control`: Quantization as trigger condition
- `mini-time-delay-system`: Combined quantization + delay analysis
- `mini-consensus-over-network`: Quantized consensus algorithms
- `mini-cyber-physical-security`: Quantization as security primitive
