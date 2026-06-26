# Course Dependency Tree: mini-event-based-control

## Prerequisites
```
mini-event-based-control
├── mini-lyapunov-direct-method  [ISS theory, Lyapunov functions]
│   ├── Nonlinear systems theory
│   ├── Quadratic forms
│   └── Positive definiteness
├── mini-linear-system-theory   [State-space models, A,B,K]
│   ├── Controllability/observability
│   ├── Pole placement
│   └── Matrix algebra
├── mini-stability-theory       [Stability definitions, ISS]
│   ├── Asymptotic stability
│   ├── Exponential stability
│   └── Input-to-state stability
├── mini-networked-control:[others]
│   ├── mini-time-delay-system
│   ├── mini-packet-loss-control
│   └── mini-quantized-control
└── Numerical methods
    ├── Matrix exponential
    ├── Schur decomposition
    └── Runge-Kutta integration
```

## Knowledge Dependencies
| Concept | Source Module | Used In |
|---------|--------------|---------|
| Lyapunov function V(x) | mini-lyapunov-direct-method | ebc_stability.c |
| ISS stability | mini-stability-theory | ebc_stability.c |
| State-space model | mini-linear-system-theory | ebc_core.c |
| Matrix exponential | Numerical methods | ebc_self.c |
| Schur decomposition | Numerical linear algebra | ebc_stability.c |

## Downstream Dependencies
Modules that depend on this module:
- mini-consensus-over-network: Event-triggered consensus
- mini-cyber-physical-security: Secure event-based control
- mini-cloud-control-system: Cloud-based event-triggered updates
