# Knowledge Graph: Event-Based Control

## L1: Definitions
| Concept | Definition | File |
|---------|-----------|------|
| Event-Triggered Control (ETC) | Control updates triggered by events, not a fixed clock | ebc_core.h: EBC_Paradigm |
| Self-Triggered Control (STC) | Next update time computed from current state + model | ebc_self.h |
| Periodic ETC (PETC) | State sampled periodically, control updated only on event | ebc_periodic.h |
| Inter-Event Time (IET) | Time between consecutive control updates | ebc_core.h: EBC_IET_Class |
| Measurement Error | e(t) = x(t_k) - x(t) for t in [t_k, t_{k+1}) | ebc_core.h: EBC_System.e |
| Event Condition | Gamma(x, e) >= 0 triggers update | ebc_trigger.h |
| Zeno Behavior | Infinite events in finite time | ebc_core.h: ebc_detect_zeno |

## L2: Core Concepts
| Concept | Description | File |
|---------|------------|------|
| Sample-and-Hold | u(t) = k(x(t_k)) constant between events | ebc_core.c |
| Minimum Inter-Event Time | Lower bound tau_min > 0 guarantees no Zeno | ebc_stability.c |
| Communication Reduction | Events saved vs periodic baseline | ebc_performance.c |
| Trigger Threshold | sigma*abs(x) + epsilon determines firing | ebc_trigger.c |

## L3: Mathematical Structures
| Structure | Type | File |
|-----------|------|------|
| EBC_System | State-space system with event tracking | ebc_core.h |
| EBC_TriggerParams | Threshold parameters (sigma, epsilon) | ebc_core.h |
| EBC_StabilityCert | ISS-Lyapunov characterization | ebc_stability.h |
| EBC_DynamicTrigger | Internal dynamic variable for advanced triggering | ebc_trigger.h |
| EBC_PETC_System | Periodic ETC state machine | ebc_periodic.h |
| LyapunovFunction (Lean) | Formal Lyapunov function type | event_control.lean |

## L4: Fundamental Laws
| Theorem | Statement | Source |
|---------|-----------|--------|
| Tabuada ISS-ETC | If ISS-Lyapunov function exists and sigma < sigma_crit, system is asymptotically stable | ebc_stability.c: ebc_stability_certify_linear |
| Linear ETC Stability | For A+BK Hurwitz, there exists sigma in (0,1) guaranteeing GES | ebc_stability.c: ebc_critical_sigma |
| Positive MIET | With epsilon > 0, inter-event time has positive lower bound | ebc_stability.c: ebc_minimum_iet_linear |
| Zeno-Free | Mixed threshold with epsilon > 0 rules out Zeno | ebc_stability.c: ebc_zeno_free_proof |

## L5: Algorithms/Methods
| Algorithm | Description | File |
|-----------|------------|------|
| ETC Event Check | Evaluate mixed threshold condition | ebc_trigger.c |
| Lyapunov Equation Solver | Bartels-Stewart via Schur decomposition | ebc_stability.c |
| Self-Triggered Next Time | Bisection search with matrix exponential | ebc_self.c |
| Matrix Exponential | Scaling-and-squaring (Higham 2005) | ebc_self.c |
| Dynamic Trigger Update | Euler integration of eta dynamics | ebc_trigger.c |
| PETC Step | Sample-evaluate-transmit state machine | ebc_periodic.c |
| RK4 Integration | 4th-order Runge-Kutta for ETC simulation | ebc_core.c |

## L6: Canonical Problems
| Problem | Example | File |
|---------|---------|------|
| Event-Triggered Stabilization | Double integrator with ETC | examples/example1_etc_stabilization.c |
| Trigger Type Comparison | Absolute vs relative vs mixed | examples/example2_trigger_comparison.c |
| Performance Benchmarking | ISE/IAE/ITAE vs comm reduction | examples/example3_performance_analysis.c |

## L7: Applications
| Application | Description | File |
|------------|-------------|------|
| Networked Control | Bandwidth-limited communication reduction | ebc_performance.h: EBC_RobustnessResult |
| Wireless Sensor Networks | Send-on-delta for energy efficiency | ebc_trigger.c: ebc_trigger_send_on_delta |
| Cyber-Physical Systems | Robustness under disturbances | ebc_performance.c: ebc_robustness_analysis |

## L8: Advanced Topics
| Topic | Description | File |
|-------|-------------|------|
| Dynamic Event-Triggering | Girard (2015) with internal variable | ebc_trigger.c: EBC_DynamicTrigger |
| Nonlinear STC | Lipschitz-bound-based next-time computation | ebc_self.c: ebc_self_next_time_nonlinear |
| Pareto Frontier | Communication-performance trade-off analysis | ebc_performance.c: ebc_pareto_frontier |

## L9: Research Frontiers
| Topic | Status | File |
|-------|--------|------|
| Meta-complexity of scheduling | Documented | event_control.lean |
| Distributed ETC | Documented | docs/course-alignment.md |
| Event-triggered learning | Documented | docs/gap-report.md |
