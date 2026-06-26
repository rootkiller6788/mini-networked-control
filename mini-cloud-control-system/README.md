# Mini Cloud Control System

Cloud-based control system simulation with delay compensation,
resource management, and network topology modeling.

## Core Concept

Cloud control moves computation from local controllers to cloud infrastructure, enabling
scalable, elastic, and globally optimized control. The key challenge is network-induced
delay, which threatens stability and performance.

## Key Equations

### Cloud Control Loop with Delay
```
x_dot = A x(t) + B u(t - tau)
y = C x(t - tau)
```
where tau is the round-trip network delay.

### Smith Predictor Compensation
```
y_compensated = y_model(t) + [y_measured(t) - y_model(t - tau)]
```

### Rate Monotonic Schedulability (Liu & Layland, 1973)
```
sum(U_i) <= n * (2^(1/n) - 1),  U_i = C_i / T_i
```

### Lyapunov-Krasovskii Functional
```
V(x_t) = x^T(t) P x(t) + integral_{t-tau}^{t} x^T(s) Q x(s) ds
```

## Architecture

```
[Plant] <---> [Edge Node] <---> [Cloud Controller]
                  |                    |
            Fast loop (us)       Slow loop (ms)
            Local K, observer    Optimization, learning
```

## API Reference

### Core Control (cloud_control_core.h)
- `ccs_create`, `ccs_free`, `ccs_set_plant_model`, `ccs_set_controller`
- `ccs_compute_control`, `ccs_update_observer`, `ccs_apply_control`, `ccs_step`
- `ccs_is_stable`, `ccs_max_allowable_delay`, `ccs_compute_performance`
- `ccs_switch_mode`, `ccs_random_init`, `ccs_compare`, `ccs_reset_metrics`

### Delay Models (cloud_control_delay.h)
- `delay_stats_compute`, `delay_is_stable`, `delay_mati_compute`
- `delay_lyapunov_krasovskii`, `delay_generate`, `delay_trace_generate`
- `delay_estimate_online`, `delay_jitter_compute`, `delay_packet_loss_correlate`
- `smith_create`, `smith_predict`, `smith_adapt_delay`, `smith_handle_oos_measurement`
- `delay_compensated_step`, `delay_fit_exponential`, `delay_fit_normal`

### Resource Management (cloud_control_resource.h)
- `rpool_create`, `rpool_allocate`, `rpool_release`, `rpool_rebalance`
- `sched_rate_monotonic_bound`, `sched_edf_bound`, `sched_deadline_monotonic`
- `elastic_create`, `elastic_evaluate`, `elastic_scale_up`, `elastic_scale_down`
- `wp_create`, `wp_forecast_holtwinters`, `wp_forecast_ema`, `wp_forecast_ar`
- `qos_resource_map`, `qos_isolation_check`

### Network Layer (cloud_control_network.h)
- `topology_create`, `topology_add_link`, `topology_find_shortest_path`
- `pktbuf_enqueue`, `pktbuf_dequeue`, `pkt_serialize`, `pkt_deserialize`
- `bw_alloc_reserve`, `bw_alloc_release`, `bw_alloc_get_available`
- `edge_to_cloud_send`, `cloud_to_edge_send`, `edge_holdover_check`
- `cloud_latency_estimate`, `net_detect_condition`

## Lean 4 Formal Verification

Seven theorems formalized:
- **zero_delay_stable_implies_mati_exists**: Stability at zero delay => MATI > 0
- **lyapunov_krasovskii_stability**: LMI-based stability condition
- **smith_predictor_perfect_compensation**: Exact model => delay-free feedback
- **rate_monotonic_bound**: RM schedulability (Liu & Layland)
- **edf_optimal_bound**: EDF necessary and sufficient condition
- **edge_holdover_stability**: Finite holdover preserves stability
- **qos_isolation_guarantee**: Reserved capacity prevents deadline misses

## Build & Test

```bash
make          # Build static library libcloud_control.a
make test     # Build and run test suite (32+ asserts)
make examples # Build all 3 examples
make demo     # Build and run demos
make bench    # Build and run benchmarks
make clean    # Clean build artifacts
```

## Directory Structure

```
mini-cloud-control-system/
  Makefile                        # Build system
  README.md                       # This file
  include/
    cloud_control_core.h          # Core types, control loop (218 lines)
    cloud_control_delay.h         # Delay models, Smith predictor (417 lines)
    cloud_control_resource.h      # Resource management, scheduling (171 lines)
    cloud_control_network.h       # Topology, packets, edge-cloud (150 lines)
  src/
    cloud_control_core.c          # Lifecycle, QR algorithm, control loop (501 lines)
    cloud_control_delay.c         # Delay stats, Smith predictor, generation (503 lines)
    cloud_control_resource.c      # Pool, scheduling, elastic, workload (541 lines)
    cloud_control_network.c       # Topology, buffer, serialization (554 lines)
    cloud_control.lean            # Lean 4 formal verification (7 theorems)
  tests/
    test_cloud_control.c          # 38 test functions, 32+ mathematical asserts
  examples/
    demo_cloud_control_basic.c    # Double integrator via cloud with delay
    demo_delay_compensation.c     # Smith predictor vs uncompensated
    demo_resource_scaling.c       # Elastic scaling with scheduling
  docs/
    knowledge-graph.md, coverage-report.md, gap-report.md
    course-alignment.md, course-tree.md
```

## Quality Metrics

| Metric | Status |
|--------|--------|
| include/ .h files | 4 (req >= 4) |
| src/ .c files | 4 (req >= 4) |
| include/ + src/ lines | 3055 (req >= 3000) |
| Exported functions | 100+ |
| Core structs | 20+ |
| Test asserts | 32+ |
| Examples | 3 (req >= 3) |
| Lean theorems | 7 |
| Docs | 5/5 |

## Key References

- Zhang et al., "Stability of Networked Control Systems", IEEE CSM, 2001
- Smith, "A controller to overcome dead time", ISA Journal, 1959
- Liu & Layland, "Scheduling Algorithms for Multiprogramming", JACM, 1973
- Schenato, "Optimal estimation in networked control...", Automatica, 2008
- Nesic & Teel, "Input-output stability of networked control systems", Automatica, 2004
- Richard, "Time-delay systems: an overview of some recent advances", Automatica, 2003
- Walsh et al., "Stability analysis of networked control systems", IEEE CST, 2002
- Xia et al., "Cloud Control Systems", IEEE/CAA JAS, 2015
- Barroso & Holzle, "The Datacenter as a Computer", 2019
- Abeni & Buttazzo, "QoS guarantee using probabilistic deadlines", ECRTS, 1999

## Course Alignment

| Course | Topic | Implementation |
|--------|-------|----------------|
| MIT 6.241J | State-space, stability | `ccs_is_stable`, eigenvalue analysis |
| Stanford AA203 | Optimal control | State-feedback LQR, observer |
| Berkeley EE221A | Linear systems | Controllability, state estimation |
| Berkeley EE222 | Nonlinear systems | Lyapunov-Krasovskii analysis |
| CMU 18-771 | Linear systems | Observer theory, estimation |
| Princeton MAE 546 | Optimal control | Smith predictor, delay handling |
| Caltech CDS110 | Dynamical systems | Cloud-based control demo |
| Cambridge 4F3 | Nonlinear control | Delay stability analysis |
| Oxford C20 | Adaptive control | Adaptive delay estimation |
| ETH 227-0216 | System identification | Online delay identification |

## Module Status: COMPLETE

- **L1** Definitions: Complete -- 13+ typedefs, 4 headers, all structs with mathematical backing
- **L2** Core Concepts: Complete -- Cloud-edge, delay compensation, elastic scaling, QoS
- **L3** Math Structures: Complete -- State-space, LK functional, RM/EDF bounds, Dijkstra
- **L4** Fundamental Laws: Complete -- 7 Lean theorems, 10+ stability asserts
- **L5** Algorithms: Complete -- QR eigenvalue, Smith predictor, Holt-Winters, Dijkstra
- **L6** Canonical Problems: Complete -- 3 examples (basic control, delay comp, scaling)
- **L7** Applications: Partial -- 3 (industrial cloud, multi-region, scheduling)
- **L8** Advanced Topics: Partial -- 4 (stochastic workload, LK stability, multi-path, elastic)
- **L9** Research Frontiers: Partial -- 3 topics documented (federated learning, quantum-safe, edge AI)
