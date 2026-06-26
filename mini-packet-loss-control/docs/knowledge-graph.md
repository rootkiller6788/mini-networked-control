# Knowledge Graph — Packet Loss in Networked Control

## L1: Definitions

- Packet loss process (Bernoulli, Gilbert-Elliott, K-state Markov, Burst)
- Packet arrival indicator γ_k ∈ {0,1}
- Network transport protocols: TCP-like vs UDP-like
- Hold strategies: zero-input, zero-order, predictive, LQG-optimal
- Channel state (Good/Bad for GE model)
- Stationary distribution π of Markov channel
- Burstiness metric: P(loss|prev_loss) / P(loss)

## L2: Core Concepts

- Memoryless property of Bernoulli loss
- Temporal correlation in Gilbert-Elliott (Markov memory)
- Separation principle under TCP-like protocols
- Non-separation under UDP-like protocols
- Buffer exhaustion in PPC
- Belief state in UDP-like control
- Mode-dependent estimation (IMM)

## L3: Math Structures

- Markov chain transition matrix P (K×K, row-stochastic)
- Stationary distribution: π = π·P
- LTI system: (A, B, C, Q, R) state-space model
- Controllability matrix: C = [B, AB, ..., A^{n-1}B]
- Observability matrix: O = [C; CA; ...; CA^{n-1}]
- Kronecker product: A ⊗ B for MJLS stability
- Modified Riccati equation for intermittent observations

## L4: Fundamental Laws

- Sinopoli et al. (2004) Theorem 2: γ_c ≥ 1 - 1/ρ(A)²
- MJLS mean-square stability: ρ(A) < 1 (Costa et al. 2005)
- Separation holds for TCP-like (Schenato et al. 2007)
- Separation fails for UDP-like (Gupta et al. 2007)
- PPC stability: H ≥ max consecutive losses
- Bounded-error estimation: ellipsoid intersection

## L5: Algorithms/Methods

- Bernoulli channel simulation (inverse CDF)
- Gilbert-Elliott channel simulation (Markov step)
- Power iteration for steady-state distribution
- Value iteration for DARE (LQR gain)
- Standard Kalman filter (predict + update)
- Intermittent Kalman filter (conditional update)
- Mode-dependent Kalman filter (IMM mixing)
- Packetized predictive control buffer generation
- Hold strategy adaptive selection
- Reference governor for constraint enforcement
- Loss-constrained optimal controller (backward DP)

## L6: Canonical Problems

- Inverted pendulum over lossy network (stabilization under loss)
- Vehicle platoon with inter-vehicle packet loss (string stability)
- Bernoulli channel statistical validation
- Gilbert-Elliott burst detection and analysis

## L7: Applications

- Remote control over wireless networks (WiFi, 5G)
- Teleoperation with communication latency/loss
- Autonomous vehicle platooning (V2V communication)
- UAV formation control with lossy links
- Smart grid demand response over lossy networks
- Martian rover remote control (delay-tolerant networking)

## L8: Advanced Topics

- Fading channel models (Rayleigh, Rician) for wireless NCS
- Set-membership estimation with bounded noise and loss
- MJLS coupled Lyapunov equations for stability
- Predictive compensation with buffer management
- Optimal UDP-like control structure (nonlinear policy)

## L9: Research Frontiers

- Information-theoretic limits: Nair-Evans data rate theorem
- Event-triggered control with packet loss (minimum inter-event time)
- Co-design of communication and control
- Learning-based packet loss compensation
- Age of Information (AoI) for NCS performance
