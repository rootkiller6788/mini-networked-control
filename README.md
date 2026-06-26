# Mini Networked Control

A collection of **from-scratch, zero-dependency C implementations** of networked control systems (NCS) theory — when feedback loops close over communication networks, fundamental tradeoffs emerge between control performance, bandwidth, quantization, delay, and packet loss. Each sub-module maps to MIT, Stanford, and other top-tier university courses, covering data rate theorems, event-triggered control, consensus, cyber-physical security, and more.

## Sub-Modules

| Sub-Module | Topics | Key Courses |
|------------|--------|-------------|
| [mini-bandwidth-limited-control](mini-bandwidth-limited-control/) | Data rate theorem, Shannon-Hartley channel capacity, encoding/quantization strategies, minimum attention control, bandwidth scheduling | MIT 6.241J, Stanford AA 203 |
| [mini-cloud-control-system](mini-cloud-control-system/) | Cloud-edge collaboration, hierarchical control, network delay models, resource allocation (RM/EDF scheduling), Smith predictor | MIT 6.824, Stanford CS 244B |
| [mini-consensus-over-network](mini-consensus-over-network/) | Graph Laplacians, continuous/discrete consensus, gossip algorithms, finite-time consensus, event-triggered consensus, convergence rate analysis | MIT 6.241J, Stanford AA 203 |
| [mini-cyber-physical-security](mini-cyber-physical-security/) | FDI (false data injection) detection, chi-squared/CUSUM detectors, game-theoretic security, ℓ₀/ℓ₁ secure estimation, physical watermarking | MIT 6.241J, Stanford AA 273 |
| [mini-event-based-control](mini-event-based-control/) | Event-triggered (ETC), self-triggered (STC), periodic event-triggered (PETC), ISS-Lyapunov stability, trigger condition design, minimum inter-event time | MIT 6.241J, Stanford AA 203 |
| [mini-packet-loss-control](mini-packet-loss-control/) | Kalman filtering with intermittent observations, LQG over lossy networks (TCP/UDP-like), Markovian jump linear systems, packetized predictive control, critical loss probability | MIT 6.241J, Stanford AA 273 |
| [mini-quantized-control](mini-quantized-control/) | Data rate theorem, uniform/logarithmic/vector/dither quantizers, Huffman/arithmetic coding, delta modulation, DPCM, run-length encoding | MIT 6.241J, Stanford AA 203 |
| [mini-time-delay-system](mini-time-delay-system/) | DDE numerical solvers, Lyapunov-Krasovskii functionals, Nyquist criterion for delay systems, Smith predictor, networked delay models (UDP/TCP/CAN/Ethernet) | MIT 6.241J, Stanford AA 203 |

## Design Philosophy

- **Zero external dependencies** — pure C99/C11, only standard library headers
- **Self-contained sub-modules** — each has its own `include/`, `src/`, `CMakeLists.txt`, and smoke tests
- **Theory-to-code mapping** — every module includes inline references to seminal papers (Wong & Brockett 1997, Nair & Evans 2004, Tabuada 2007, Schenato et al. 2007, etc.)
- **NCS-first perspective** — all control design explicitly accounts for communication constraints: bit-rate limits, packet loss, latency, quantization, and security

## Building

Each sub-module is standalone. Build with CMake:

```bash
cd mini-bandwidth-limited-control
mkdir build && cd build
cmake ..
make
./smoke_test
```

Requires a **C99-compliant compiler** and **CMake ≥ 3.14**.

## Project Structure

```
20. mini-networked-control/
├── mini-bandwidth-limited-control/   # Data rate theorem, Shannon capacity, encoding, bandwidth scheduling
├── mini-cloud-control-system/        # Cloud-edge collaboration, network delay, resource allocation
├── mini-consensus-over-network/      # Graph Laplacians, consensus protocols, gossip, finite-time consensus
├── mini-cyber-physical-security/     # FDI detection, game-theoretic security, secure estimation, watermarking
├── mini-event-based-control/         # ETC, STC, PETC, ISS-Lyapunov stability, trigger conditions
├── mini-packet-loss-control/         # Kalman with intermittent obs., LQG over lossy links, predictive control
├── mini-quantized-control/           # Data rate theorem, uniform/log/vector/dither quantizers, encoding
├── mini-time-delay-system/           # DDE solvers, Lyapunov-Krasovskii, Smith predictor, networked delay
├── .gitignore
├── README.md
└── README-CN.md
```

## License

MIT
