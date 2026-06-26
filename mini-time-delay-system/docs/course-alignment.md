# Course Alignment — mini-time-delay-system

Mapping time-delay system topics to nine world-class university curricula.

## MIT

| Course | Topic | Module Coverage |
|--------|-------|----------------|
| 6.241J Dynamic Systems & Control | LTI system stability with delays | `delay_stability.c` — Nyquist, delay margin |
| 16.323 Optimal Control | LQR with delayed actuation | `TimestampLQG` |
| 6.832 Underactuated Robotics | Control of systems with input delay | Smith predictor, predictor feedback |
| 2.151 Nonlinear Control | Lyapunov-Krasovskii for nonlinear DDEs | `lyapunov_krasovskii.c` — Razumikhin |

## Stanford

| Course | Topic | Module Coverage |
|--------|-------|----------------|
| AA203 Optimal Control | NCS with delays (Hespanha et al.) | `networked_delay.c` — NCS step, packet loss |
| EE363 Convex Optimization | LMI for delay-dependent stability | `lmi_delay_dependent_check` |
| CS358 Circuit Complexity | N/A (different domain) | — |

## Berkeley

| Course | Topic | Module Coverage |
|--------|-------|----------------|
| EE221A Linear Systems | Matrix measure, stability margins | `matrix_measure_l2`, `delay_gain_margin` |
| EE222 Nonlinear Systems | Krasovskii method for DDEs | LK functional evaluation |
| ME233 Advanced Control | Dead-time compensation | Smith predictor implementation |

## CMU

| Course | Topic | Module Coverage |
|--------|-------|----------------|
| 18-771 Linear Systems | Pole placement with delay | Characteristic root computation |
| 24-677 Nonlinear Control | Razumikhin theorem | `razumikhin_condition_check` |
| 24-654 Systems Thinking | System dynamics with delays (Sterman) | Delay type taxonomy |

## Princeton

| Course | Topic | Module Coverage |
|--------|-------|----------------|
| MAE 546 Optimal Control | NCS with lossy networks (Schenato) | TimestampLQG |
| ELE 530 Estimation | Kalman filter with delayed measurements | `tslqg_update` |
| COS 522 Computational Complexity | N/A | — |

## Caltech

| Course | Topic | Module Coverage |
|--------|-------|----------------|
| CDS110 Intro to Control | Delay stability fundamentals | Delay margin analysis |
| CDS140 Nonlinear Dynamics | Bifurcations in DDEs | Characteristic root Newton solver |

## Cambridge

| Course | Topic | Module Coverage |
|--------|-------|----------------|
| 4F3 Nonlinear & Predictive Control | Smith predictor, dead-time | `smith_predictor.c` |
| 4F2 Robust Control | Delay robustness margins | `delay_margin_frequency_sweep` |

## Oxford

| Course | Topic | Module Coverage |
|--------|-------|----------------|
| B4 Predictive Control | Delay in MPC | Documented in gap-report |
| C20 Adaptive Control | Adaptive delay compensation | Predictor feedback (Krstic) |

## ETH

| Course | Topic | Module Coverage |
|--------|-------|----------------|
| 227-0216 System Identification | Identification of delay systems | Delay type taxonomy |
| 227-0220 Model Reduction | Padé approximation of delays | `pade_create`, `pade_augment_system` |
| 227-0690 Advanced Control | LMI for time-delay systems | LMI-based stability checks |

## Key References

1. Hale & Verduyn Lunel (1993) — *Introduction to Functional Differential Equations* [MIT, Cambridge]
2. Gu, Kharitonov, Chen (2003) — *Stability of Time-Delay Systems* [Berkeley, ETH]
3. Niculescu (2001) — *Delay Effects on Stability* [Stanford, Caltech]
4. Krstic (2009) — *Delay Compensation* [MIT, Princeton]
5. Hespanha et al. (2007) — *NCS Survey* [Stanford, CMU]
6. Normey-Rico & Camacho (2007) — *Control of Dead-time Processes* [Oxford, Cambridge]
