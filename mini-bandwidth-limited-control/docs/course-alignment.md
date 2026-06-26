# Course Alignment — Bandwidth-Limited Control

## Nine-School Curriculum Mapping

| School | Course | Topic | Our Implementation |
|--------|--------|-------|--------------------|
| **MIT** | 6.241J Dynamic Systems | Lyapunov stability, state feedback | `blc_lqr_solve`, closed-loop pole analysis |
| **MIT** | 16.323 Optimal Control | LQR, ARE, Riccati equation | `BLCLQRController`, Newton-Kleinman iteration |
| **MIT** | 6.832 Underactuated Robotics | Event-triggered control | `BLCEventTriggeredFeedback`, Tabuada trigger |
| **Stanford** | AA203 Optimal Control | Optimal control with constraints | `BLCMPCController`, L1-penalized MPC |
| **Stanford** | EE363 Convex Optimization | Rate-distortion optimization | `blc_rate_distortion`, `blc_distortion_rate` |
| **Berkeley** | EE221A Linear Systems | State-space, eigenvalue analysis | `blc_eigenvalues` (QR algorithm) |
| **Berkeley** | EE222 Nonlinear Systems | Quantization nonlinearity, Circle Criterion | `BLCRobustQuantCtrl`, sector analysis |
| **CMU** | 18-771 Linear Systems | Controllability, observability | Plant model, observer gain matrix |
| **Princeton** | MAE 546 Optimal Control | Dynamic programming | MPC with rate constraints |
| **Princeton** | ELE 530 Estimation | Kalman filtering, prediction | `BLCEncoder`, predictive encoding |
| **Caltech** | CDS110 Intro Control | PID, sampling, quantization | Quantizer model, zoom quantizer |
| **Caltech** | CDS140 Nonlinear Dynamics | Bifurcation, chaos | Stability under quantization |
| **Cambridge** | 4F3 Nonlinear Ctrl | Lyapunov methods, passivity | Absolute stability, Circle Criterion |
| **Cambridge** | 4F2 Robust Ctrl | H∞ control, uncertainty | `blc_robust_hinf_norm` |
| **Oxford** | B4 Predictive Ctrl | MPC, constraints | `BLCMPCController` |
| **Oxford** | C20 Adaptive Ctrl | Self-tuning, adaptation | `BLCAdaptiveScheduler` |
| **ETH** | 227-0216 Sys Identification | System modeling, SVD | `blc_singular_values` |
| **ETH** | 227-0220 Model Reduction | Balanced truncation | `blc_lyapunov_solve` |

## Key Textbook Alignment

| Textbook | Chapter | Content | Implementation |
|----------|---------|---------|----------------|
| Khalil (2002) | Ch. 10 | Absolute stability, Circle Criterion | `blc_robust_circle_criterion` |
| Bertsekas (2012) | Ch. 4 | LQR, ARE | `blc_lqr_solve` |
| Åström-Wittenmark (2013) | Ch. 9 | Adaptive control | `BLCAdaptiveScheduler` |
| Rawlings-Mayne-Diehl (2017) | Ch. 2 | MPC formulation | `BLCMPCController` |
| Simon (2006) | Ch. 5 | Kalman filter | `BLCEncoder` (predictive) |
| Cover-Thomas (2006) | Ch. 5 | Huffman, entropy coding | `BLCHuffmanCoder` |
| Golub-Van Loan (2013) | Ch. 7 | QR algorithm | `blc_eigenvalues` |
