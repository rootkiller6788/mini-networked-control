# Knowledge Graph — Consensus over Network

## L1: Definitions
- Consensus: agreement of all agent states (x_i = x_j for all i,j)
- Communication graph G = (V,E), Laplacian L = D - A
- Agreement space: span{1}, Disagreement space: orthogonal complement
- Algebraic connectivity: lambda_2(L), Fiedler eigenvalue
- Consensus protocol types: continuous-time, discrete-time, gossip, event-triggered
- Agent dynamics: single/double integrator, general linear, unicycle
- Weight matrix: stochastic, doubly stochastic, Perron-Frobenius eigenvalue
- Convergence criteria: absolute, relative, disagreement-based

## L2: Core Concepts
- Distributed averaging: compute mean without central coordination
- Connectivity requirement: graph must be connected (undirected) or have spanning tree (directed)
- Convergence rate vs. spectral gap: larger lambda_2 => faster convergence
- Leader-follower consensus: one agent anchors, others follow
- Average preservation: doubly stochastic W preserves sum
- Disagreement energy as Lyapunov function
- Second-order consensus: position + velocity agreement
- Proportional interaction: control = weighted sum of relative differences

## L3: Mathematical Structures
- Graph Laplacian L = D - A (N x N symmetric PSD)
- Normalized Laplacian: I - D^{-1/2} A D^{-1/2}
- Perron matrix W: row-stochastic, compatible with graph
- Incidence matrix B: L = B B^T
- Vector space R^{N x d}: agent states as matrix
- Agreement subspace: null(L) = span{1}
- Spectral decomposition: L = U Lambda U^T
- Graph-theoretic metrics: diameter, average path length, clustering coefficient

## L4: Fundamental Theorems
- T1: Consensus Theorem (Olfati-Saber 2004): dx/dt = -Lx converges to average iff G connected
- T2: Fiedler Theorem (1973): lambda_2(L) > 0 iff G connected
- T3: Perron-Frobenius: stochastic W has eigenvalue 1, convergence determined by second largest
- T4: Convergence rate: ||x(t)-x*|| <= exp(-lambda_2 t) ||x(0)-x*||
- T5: Average invariance: 1^T x(t) constant for doubly stochastic W
- T6: Lyapunov stability: V = (1/2)x^T L x, dV/dt = -x^T L^2 x <= 0
- T7: LaSalle invariance: limit set = agreement subspace
- T8: Switching topology (Jadbabaie 2003): consensus if union graph jointly connected
- T9: Delay margin: tau < pi/(2*lambda_N) for uniform delay

## L5: Algorithms/Methods
- A1: Continuous-time Laplacian flow (Euler/RK4 integration)
- A2: Discrete-time weighted averaging (DeGroot iteration)
- A3: Gossip algorithm (pairwise randomized updates, Boyd 2006)
- A4: Push-sum for directed graphs (Kempe 2003)
- A5: Max/Min consensus (converge in diameter steps)
- A6: Finite-time consensus (minimal polynomial, Sundaram 2007)
- A7: Event-triggered consensus (Dimarogonas 2012)
- A8: Weight design: Metropolis-Hastings, max-degree, optimal (Xiao-Boyd 2004)
- A9: Fastest mixing Markov chain design

## L6: Canonical Problems
- P1: Rendezvous problem (multi-agent meeting)
- P2: Distributed averaging (sensor fusion)
- P3: Clock synchronization (sensor networks)
- P4: Formation control (UAV swarms)
- P5: Flocking/alignment (Reynolds rules)
- P6: Autonomous vehicle platooning
- P7: Load balancing (smart grid)

## L7: Applications
- Multi-robot rendezvous (2D/3D position consensus)
- Distributed clock synchronization (TimeSynch protocol)
- UAV formation control (V-formation, offset consensus)
- Autonomous vehicle platooning (spacing consensus)
- Distributed sensor fusion / estimation
- Smart grid load balancing (power consensus)
- Multi-camera calibration (parameter consensus)

## L8: Advanced Topics
- Consensus with communication delays (delay margin analysis)
- Resilient consensus under Byzantine/malicious agents (MSR algorithm)
- Consensus under switching topologies (joint connectivity)
- Quantized consensus (probabilistic/uniform quantization)
- Resilient weight design against link failures
- Bipartite consensus / signed graphs (structural balance)
- Finite-time exact consensus
- Zeno-free event-triggered consensus
- Second-order consensus with damping
- Push-sum for directed graphs

## L9: Research Frontiers
- Byzantine fault tolerance in dynamic networks
- Secure consensus under adversarial attacks
- Learning-based consensus (RL for weight adaptation)
- Quantum consensus protocols
- Consensus in time-varying directed graphs
- Privacy-preserving consensus (differential privacy)
- Resilient consensus with trust metrics
