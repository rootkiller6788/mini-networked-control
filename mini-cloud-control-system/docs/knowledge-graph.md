# Knowledge Graph -- mini-cloud-control-system

## L1: Definitions
- CloudControlSystem, CloudNode, EdgeNode, ControlTask, QoSProfile
- SystemState, CloudControlMode, CloudControlState
- NetworkDelay, DelayStatistics, DelayModelType
- SmithPredictor, CompensationStrategy
- ResourcePool, ComputeNode, ElasticController, WorkloadPredictor
- EdgeCloudTopology, NetworkLink, PacketBuffer, BandwidthAllocator
- ControlPacket, PacketHeader, MultiPathRoute, RouteType

## L2: Core Concepts
- Cloud-edge collaborative control, Hierarchical control
- Operating modes: EDGE_ONLY,CLOUD_ONLY,COLLABORATIVE,HIERARCHICAL,HYBRID,REDUNDANT
- Delay compensation (Smith predictor principle)
- Elastic resource scaling (MAPE-K loop), QoS-aware allocation
- Real-time scheduling (RM, EDF, DM)
- Edge holdover capability, Out-of-order measurement handling
- Multi-tenant isolation, Bumpless mode transfer

## L3: Mathematical Structures
- State-space: x_dot = A x + B u, y = C x + D u
- Luenberger observer: x_hat_dot = A x_hat + B u + L(y - C x_hat)
- Delayed closed-loop: x_dot = A x + B u(t-tau)
- Lyapunov-Krasovskii functional for delay stability
- RM bound: sum(U_i) <= n(2^(1/n)-1), EDF bound: sum(U_i) <= 1.0
- Holt-Winters triple exponential smoothing
- Dijkstra shortest path, Bandwidth-delay product

## L4: Fundamental Theorems
- Zero-delay stability => MATI exists (continuity of eigenvalues)
- Lyapunov-Krasovskii stability (LMI conditions)
- Smith predictor perfect compensation (exact model => delay-free)
- RM schedulability bound (Liu & Layland 1973)
- EDF optimality (necessary and sufficient)
- Edge holdover stability preservation
- QoS isolation guarantee under reserved capacity

## L5: Algorithms
- QR algorithm with Wilkinson shift (eigenvalues)
- Power iteration for dominant eigenvalue
- Luenberger observer (predictor-corrector)
- Smith predictor delay compensation
- Online delay estimation (EMA), Delay trace generation
- Holt-Winters workload forecasting
- Elastic scaling (MAPE-K loop)
- Dijkstra shortest path routing
- Audsley response time analysis

## L6: Canonical Problems
- Cloud-controlled double integrator with delay
- Smith predictor vs. uncompensated control comparison
- Elastic resource scaling under varying workload

## L7: Applications
- Industrial cloud control (edge-to-cloud topology)
- Multi-region cloud deployment (latency modeling)
- Real-time scheduling for control tasks

## L8: Advanced Topics
- Stochastic workload prediction (Holt-Winters)
- Lyapunov-Krasovskii delay stability
- Multi-path routing with reliability
- Elastic/adaptive resource scaling

## L9: Research Frontiers
- Federated learning for cloud control (documented)
- Quantum-safe cloud control (documented)
- Edge AI for ultra-low latency control (documented)
