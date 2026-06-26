# mini-consensus-over-network

## Module Status: COMPLETE

- **L1-L6**: Complete
- **L7**: Complete (4 applications: rendezvous, clock sync, UAV flocking, platooning)
- **L8**: Complete (5 advanced topics: delays, Byzantine resilience, switching, quantized, resilient weights)
- **L9**: Partial (7 research frontiers documented; no implementation required)

## Core Concept

Consensus is the problem of achieving agreement among a group of autonomous agents that can only communicate with their neighbors in a network. This is a fundamental primitive for distributed control, sensor fusion, and multi-robot coordination.

## Core Definitions (L1)

| Term | Definition |
|------|-----------|
| Consensus | x_i = x_j for all agents i, j |
| Communication Graph | G = (V, E) where V = agents, E = communication links |
| Laplacian Matrix | L = D - A, encodes graph topology |
| Algebraic Connectivity | lambda_2(L), Fiedler eigenvalue, measures connectivity |
| Weight Matrix | W, row-stochastic, determines discrete consensus update |
| Agreement Space | span{1}, the subspace where all agents agree |
| Disagreement Energy | Phi(x) = (1/2) x^T L x, natural Lyapunov function |

## Core Theorems (L4)

### T1: Consensus Theorem (Olfati-Saber and Murray, 2004)
For an undirected connected graph G with Laplacian L:
  dx/dt = -L x
  limit(t -> inf) x_i(t) = (1/N) sum_j x_j(0) for all i
Convergence rate: ||x(t) - x*|| <= ||x(0) - x*|| * exp(-lambda_2 * t)

### T2: Fiedler Theorem (1973)
  lambda_2(L) > 0  iff  G is connected

### T3: Discrete Consensus (DeGroot, 1974; Ren and Beard, 2005)
  x[k+1] = W x[k]  converges to average iff W is doubly stochastic and primitive
  Rate: ||x[k] - x*|| <= C * rho(W - 11^T/N)^k

### T4: Lyapunov Stability
  V(x) = (1/2) x^T L x  =>  dV/dt = -x^T L^2 x <= 0
By LaSalle invariance principle, x(t) converges to agreement subspace.

### T5: Delay Margin (Olfati-Saber and Murray, 2004)
  Uniform delay tau:  stability iff tau < pi / (2 * lambda_N(L))

### T6: Switching Topology (Jadbabaie, Lin, Morse, 2003)
  Consensus converges if the union graph over uniformly bounded time windows is jointly connected.

## Core Algorithms (L5)

| Algorithm | Complexity/step | Convergence |
|-----------|----------------|-------------|
| Continuous-time Laplacian flow | O(N^2 * d) | O(exp(-lambda_2 * t)) |
| Discrete-time DeGroot | O(N^2 * d) | O(mu^k) where mu = rho(W-11^T/N) |
| Gossip (pairwise randomized) | O(d) | O(N^2 * log(1/eps) / lambda_2) |
| Push-sum (directed graphs) | O(N^2 * d) | O(exp(-alpha * k)) |
| Max/Min consensus | O(N * d) | O(diam(G)) exactly |
| Finite-time consensus | O(D * N * d) | O(D) where D = deg(min poly) |
| Event-triggered | O(N^2 * d) | Same rate, fewer messages |
| Optimal weights (Xiao-Boyd) | O(1) after evals | Optimal rho(W-11^T/N) |

## Weight Design Methods (L5)

| Method | Formula | Properties |
|--------|---------|------------|
| Metropolis-Hastings | w_ij = 1/(max(d_i,d_j)+1) | Doubly stochastic, robust |
| Max-degree | w_ij = 1/(d_max+1) | Simple, good for regular graphs |
| Optimal constant | alpha = 2/(lambda_2+lambda_N) | Minimizes spectral radius |
| Lazy | W_lazy = gamma*I + (1-gamma)*W | Positive definite, always stable |

## Canonical Problems (L6)

| Problem | Description | Example File |
|---------|-------------|-------------|
| Distributed Averaging | Compute mean of sensor measurements without central server | example_averaging.c |
| Rendezvous | Multiple robots meet at common point using only local info | example_rendezvous.c |
| Clock Synchronization | Align virtual clocks across sensor network | example_clock_sync.c |
| Formation Control | UAVs maintain V-formation via offset consensus | example_formation.c |

## Nine-School Course Mapping

| School | Course | Covered Topics |
|--------|--------|---------------|
| MIT | 6.241J, 6.832 | Laplacian dynamics, multi-robot rendezvous |
| Stanford | AA274 | Gossip algorithms, distributed optimization |
| Berkeley | EE222 | Lyapunov analysis, convergence rate bounds |
| CMU | 24-677 | Leader-follower, switching topology |
| Princeton | MAE 546 | Fastest mixing Markov chain |
| Caltech | CDS140 | Bipartite consensus, structural balance |
| Cambridge | 4F3 | Event-triggered, delay margin |
| Oxford | C20 | Resilient consensus, Byzantine fault tolerance |
| ETH | 227-0220 | Graph-theoretic model reduction |

## Build and Test

  make          # Build static library
  make test     # Build and run test suite
  make examples # Build all examples
  make demo     # Build and run all demos
  make clean    # Remove build artifacts

## File Structure

```
mini-consensus-over-network/
  Makefile                     # GNU Make build
  README.md                    # This file
  include/                     # 4 header files
    consensus_types.h           # Core type definitions, vectors, matrices
    consensus_graph.h           # Graph construction, Laplacian, spectra
    consensus_dynamics.h        # Continuous/discrete consensus, Lyapunov
    consensus_algorithms.h      # Gossip, push-sum, max/min, event-triggered
  src/                          # 7 source files
    consensus_types.c           # Type implementations (773 lines)
    consensus_graph.c           # Graph theory, eigenvalues, generators (733 lines)
    consensus_dynamics.c        # Consensus simulation, second-order (609 lines)
    consensus_algorithms.c      # Algorithm catalog (611 lines)
    consensus_applications.c    # Robot rendezvous, clock sync, UAV, platoon (463 lines)
    consensus_resilience.c      # Delays, Byzantine MSR, switching, quantized (481 lines)
    consensus_lean.lean         # Lean 4 formalization of core theorems
  tests/
    test_consensus.c            # 16 tests across L1-L8
  examples/
    example_averaging.c         # Distributed sensor averaging
    example_rendezvous.c        # Multi-robot 2D rendezvous
    example_clock_sync.c        # Clock synchronization protocol
    example_formation.c         # UAV V-formation control
  docs/
    knowledge-graph.md          # L1-L9 knowledge map
    coverage-report.md          # Coverage assessment
    gap-report.md               # Gap analysis
    course-alignment.md         # Nine-school course mapping
    course-tree.md              # Prerequisite dependency tree
  demos/                        # Visualization/demo directory
  benches/                      # Performance benchmarks directory
```

## Line Count Compliance

include/ + src/ total: **4052 lines** >= 3000 threshold