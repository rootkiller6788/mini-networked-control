# Knowledge Graph — mini-cyber-physical-security

## L1: Definitions (Complete)
- CPSAttackType: DoS, FDI, Replay, Bias, Covert, Surge, Zero-Dynamics
- CPSAttackTarget: Sensor, Actuator, Controller, Network, Estimator, Reference
- CPSSecurityState: Normal, Suspicious, Attacked, Degraded, Recovering, Compromised
- CPSDetectionMethod: Chi2, CUSUM, Watermark, MMD, KL, Interval, Set-Membership
- CPSResilientStrategy: Hold, Fallback, Reconfig, Game, MPC, Redundant
- CPSDynamicalSystem: (A, B, C, D) state-space model
- CPSAttackModel: attack type, target, magnitude, stealthiness, Gamma matrices
- CPSDetector: threshold, statistic, alarm state, history buffer
- CPSWatermark: type, signal, expected response, detection statistics
- CPSSecuritySystem: unified security model encapsulating all subcomponents

## L2: Core Concepts (Complete)
- Attack detection via residual monitoring
- Attack identification via sensor subset analysis
- Secure state estimation under sensor attacks
- Resilient control: maintaining stability under attack
- Physical watermarking: proactive defense mechanism
- Game-theoretic security: optimal defense/attack strategies
- Moving target defense: time-varying security parameters
- Sensor diversity and redundancy for resilience
- Kalman filtering as the foundation for residual-based detection
- Security state machine: transitions between security levels

## L3: Mathematical Structures (Complete)
- CPSMatrix: dense matrix with dimensions, row-major storage
- CPSVector: vector abstraction for state/measurement/control
- CPSDynamicalSystem: (A,B,C,D) with noise covariances Q,R
- Householder QR decomposition for numerical rank
- Observability Gramian Wo = Σ(A')^k C' C A^k
- Controllability Gramian Wc = Σ A^k B B' (A')^k
- Gaussian noise model: w~N(0,Q), v~N(0,R)
- Attack signal generation models for each attack type
- Residual computation: r[k] = y[k] - C*ẋ[k]

## L4: Fundamental Laws (Complete)
- Kalman Observability Theorem: rank(O)=n iff (A,C) is observable
- Kalman Controllability Theorem: rank(Co)=n iff (A,B) is controllable
- Mo-Sinopoli Detectability Theorem: watermarking makes replay detectable
- Fawzi-Tabuada-Diggavi Theorem: 2s<p is necessary for secure estimation
- Pasqualetti Attack Detectability: undetectable iff in unobservable subspace
- Gramian-based security analysis: reachable sets under attack
- Chi-squared distribution property: r'Σ⁻¹r ~ χ²(p) under H0
- CUSUM optimality: minimizes worst-case detection delay (Lorden 1971)
- Wald's SPRT optimality: minimizes expected sample size
- Zero-dynamics vulnerability: invariant zeros outside unit circle

## L5: Algorithms/Methods (Complete)
- Chi-squared attack detector with Wilson-Hilferty p-value approximation
- CUSUM (Cumulative Sum) change detection algorithm
- Wald's Sequential Probability Ratio Test (SPRT)
- Kalman filter predict/update cycle
- DARE (Discrete Algebraic Riccati Equation) iterative solver
- L0 secure state estimation via combinatorial branch-and-bound
- L1 secure state estimation via ISTA (Iterative Soft-Thresholding)
- Resilient Kalman filter with attack-aware covariance inflation
- Game-theoretic saddle-point solving via fictitious play
- Watermark signal generation: Gaussian, binary, sinusoidal, chaotic, adaptive
- Encrypted watermark generation via stream cipher
- Sensor selection via greedy norm-maximization

## L6: Canonical Problems (Complete)
- FDI attack on smart grid state estimation (Mo & Sinopoli 2010)
- GPS spoofing attack on autonomous vehicles
- Replay attack on networked control system (with watermark defense)
- Stuxnet-class covert attack on industrial control systems
- Bias injection attack on water distribution SCADA
- Chi-squared evasion by stealthy FDI (attack in column space of C)
- Zero-dynamics attack exploiting invariant zeros
- DoS attack modeling as intermittent packet loss

## L7: Applications (Complete — 4 real-world scenarios)
1. Smart Grid: Power system state estimation with FDI detection
   - IEEE 5-bus system simulation
   - Chi-squared tests on power flow residuals
2. Autonomous Vehicle: GPS spoofing detection
   - Kalman filter fusion of IMU + GPS
   - Innovation-based spoofing detection
3. Industrial Control System (ICS): Replay attack defense
   - Water tank level control with watermarking
   - Replay recording and detection
4. Water Distribution SCADA: Chlorine sensor attack
   - CUSUM detection of sustained bias
   - Contamination scenario simulation

## L8: Advanced Topics (Complete — 5+ topics)
- Stackelberg security game: defender-first commitment with SSE
- Dynamic multistage game via backward induction
- Bayesian attacker type belief update
- Zero-sum game LP solving via fictitious play (Brown 1951)
- Moving target defense via watermark parameter rotation
- Encrypted control with homomorphic watermark verification
- Adaptive watermarking: energy proportional to security risk
- GLRT (Generalized Likelihood Ratio Test) for unknown attack parameters
- Attack-reachable set computation via controllability Gramian
- Sensor-optimal detection threshold under cost trade-off

## L9: Research Frontiers (Partial — documented)
- Quantum-safe cryptographic watermarking for CPS
- AI-driven adaptive attack strategies using adversarial ML
- Formal verification of CPS security protocols (Lean 4)
- Resilient control under coordinated multi-vector attacks
- Privacy-preserving attack detection via homomorphic encryption
