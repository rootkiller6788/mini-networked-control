# mini-cyber-physical-security

Cyber-Physical System Security — Attack Detection, Resilient Control & Game-Theoretic Defense

## Module Status: COMPLETE ✅

- **L1-L6**: Complete
- **L7**: Complete (4 applications: Smart Grid, Autonomous Vehicle, ICS, SCADA)
- **L8**: Complete (10 advanced topics: Stackelberg games, Bayesian fusion, MTD, encrypted watermarking, etc.)
- **L9**: Partial (5 research frontiers documented)

## Nine-Layer Knowledge Coverage

| Level | Name | Rating | Key Content |
|-------|------|--------|-------------|
| L1 | Definitions | Complete | 10 struct/enum types, attack taxonomy, security states |
| L2 | Core Concepts | Complete | Detection, resilience, watermarking, game theory |
| L3 | Math Structures | Complete | Matrix/Vector ops, QR decomposition, Gramians, state-space |
| L4 | Fundamental Laws | Complete | Kalman observability/controllability, Mo-Sinopoli, Fawzi-Tabuada-Diggavi, Pasqualetti, CUSUM optimality, SPRT optimality |
| L5 | Algorithms | Complete | Chi2, CUSUM, SPRT, Kalman, DARE, L0/L1 estimation, fictitious play, watermark generation |
| L6 | Canonical Problems | Complete | FDI, replay, GPS spoofing, bias injection, zero-dynamics, DoS |
| L7 | Applications | Complete | Smart grid, autonomous vehicle (Tesla), ICS (Stuxnet-class), SCADA (water) |
| L8 | Advanced Topics | Complete | Stackelberg SSE, dynamic games, Bayesian belief, MTD, encrypted watermarking, GLRT, adaptive watermarking |
| L9 | Research Frontiers | Partial | Quantum-safe CPS, adversarial ML, formal verification (documented) |

## Core Definitions (L1)

- **CPSAttackType**: DoS, FDI, Replay, Bias, Covert, Surge, Zero-Dynamics
- **CPSAttackTarget**: Sensor, Actuator, Controller, Network, Estimator, Reference
- **CPSSecurityState**: Normal, Suspicious, Attacked, Degraded, Recovering, Compromised
- **CPSDynamicalSystem**: (A, B, C, D) with noise covariances Q, R
- **CPSAttackModel**: type, target, magnitude, stealthiness, Gamma matrices

## Core Theorems (L4)

1. **Kalman Observability**: rank([C; CA; ...; CA^{n-1}]) = n ⇔ observable
2. **Kalman Controllability**: rank([B, AB, ..., A^{n-1}B]) = n ⇔ controllable
3. **Mo-Sinopoli Watermarking**: K-L divergence ∝ watermark energy → replay detectable
4. **Fawzi-Tabuada-Diggavi**: Secure estimation possible iff 2s < p (sensors > 2×attacks)
5. **Pasqualetti Detectability**: Attack undetectable ⇔ lies in unobservable subspace
6. **CUSUM Optimality** (Lorden 1971): Minimizes worst-case detection delay
7. **SPRT Optimality** (Wald 1945): Minimizes expected sample size
8. **Gramian Bounds**: Attack reachable set bounded by controllability Gramian

## Core Algorithms (L5)

- **Chi-squared detector**: g = r'Σ⁻¹r ~ χ²(p) under H0
- **CUSUM**: S[k] = max(0, S[k-1] + LLR[k] - drift)
- **SPRT**: Sequential testing with Wald thresholds A = log(β/(1-α)), B = log((1-β)/α)
- **Kalman filter**: Predict x̂_{k|k-1} = Ax̂_{k-1} + Bu_{k-1}, Update with innovation
- **DARE**: Iterative Riccati solver for steady-state Kalman gain
- **L0 secure estimation**: Combinatorial search over sensor subsets
- **L1 secure estimation**: ISTA with soft-thresholding prox
- **Fictitious play**: Zero-sum game LP solver
- **Watermark generation**: Gaussian, binary, sinusoidal, chaotic, PN sequence

## Canonical Problems (L6)

1. FDI attack on smart grid state estimation
2. Replay attack on NCS (defeated by watermarking)
3. GPS spoofing on autonomous vehicles
4. Stuxnet-class covert attack on ICS
5. Bias injection on water SCADA
6. Zero-dynamics attack exploiting invariant zeros
7. DoS as intermittent packet loss

## Course Alignment

| MIT | Stanford | Berkeley | CMU | Princeton | Caltech | Cambridge | Oxford | ETH |
|-----|----------|----------|-----|-----------|---------|-----------|--------|-----|
| 6.241J | AA203 | EE221A | 18-771 | ELE 530 | CDS110 | 4F3 | B4 | 227-0216 |
| 16.323 | EE363 | EE222 | 24-677 | MAE 546 | CDS140 | 4F2 | C20 | 252-0400 |

## Build & Test

```bash
make          # Build static library
make test     # Run test suite (~40 tests)
make examples # Build all examples
make run-examples  # Run all examples
make demo     # Run integrated demo
make bench    # Run benchmarks
```

## References

- Mo & Sinopoli (2010) — False Data Injection Attacks in Control Systems
- Pasqualetti, Dorfler & Bullo (2013) — Attack Detection and Identification in CPS
- Fawzi, Tabuada & Diggavi (2014) — Secure Estimation and Control for CPS
- Zhu & Basar (2015) — Game-Theoretic Methods for CPS Security
- Basseville & Nikiforov (1993) — Detection of Abrupt Changes
- Kay (1998) — Fundamentals of Statistical Signal Processing, Vol. II
- Golub & Van Loan (2013) — Matrix Computations, 4th ed.
- Smith (2015) — Covert Misappropriation of Networked Control Systems
