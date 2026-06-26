# Knowledge Graph — mini-quantized-control

## L1: Definitions (Complete)

| # | Concept | C Implementation | Lean Definition |
|---|---------|-----------------|-----------------|
| 1 | Quantizer | `QCQuantizer` struct | `QuantLevel` inductive |
| 2 | Quantized Control System | `QCSystem` struct | `QCSystemModel` structure |
| 3 | Uniform Quantizer | `qc_uniform_quantize_midtread` | - |
| 4 | Logarithmic Quantizer | `QCLogQuantizer` struct | `LogQuantizer` structure |
| 5 | Dynamic Quantizer | `QCDynamicQuantizer` struct | - |
| 6 | Encoder | `QCEncoder` struct | - |
| 7 | Decoder | `QCDecoder` struct | - |
| 8 | Data Rate | `QCDataRate` struct | - |
| 9 | Quantization Error | `qc_quantization_error` | Theorem 1 |
| 10 | Sector Bound | `QCSectorBoundResult` struct | Theorem 2 |
| 11 | Stability Status | `QCStabilityStatus` enum | `StabilityStatus` inductive |
| 12 | Quantization Cell | `QCQuantizationCell` struct | - |

## L2: Core Concepts (Complete)

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Quantization error analysis | `qc_quantization_error`, `qc_quantization_noise_variance` |
| 2 | Data rate → stability relationship | `qc_data_rate_is_stabilizable` |
| 3 | Logarithmic relative error | `qc_log_quantize`, sector delta |
| 4 | Zoom-in/zoom-out strategy | `qc_dyn_zoom_in`, `qc_dyn_zoom_out` |
| 5 | Saturation handling | `qc_quantize_scalar` with range clamping |
| 6 | Encoding schemes | Fixed-length, Huffman, arithmetic, delta, DPCM |
| 7 | Bitstream I/O | `QCBitWriter`, `QCBitReader` |
| 8 | Quantizer overload strategies | Saturate, zoom-out, modulo, extend |

## L3: Mathematical Structures (Complete)

| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | Uniform quantizer partition | `qc_uniform_quantize_midtread/midrise` |
| 2 | Logarithmic partition | `qc_log_quantizer_build_levels` |
| 3 | Lattice vector quantizer | `QCVectorCodebook`, LBG training |
| 4 | Huffman prefix code tree | `QCHuffmanTree`, priority queue build |
| 5 | Arithmetic coding interval | `QCArithEncoder`, `QCArithDecoder` |
| 6 | Delta modulation | `QCDeltaModulator` with adaptive step |
| 7 | DPCM predictor | `QCDPCMEncoder` with linear predictor |
| 8 | Communication channel models | `QCChannel` (AWGN, BSC, erasure) |
| 9 | Quantized control loop | Euler integration with quantization |
| 10 | Observer-based output feedback | Luenberger observer with quantized I/O |

## L4: Fundamental Laws / Theorems (Complete)

| # | Theorem | C Verification | Lean Statement |
|---|---------|---------------|----------------|
| 1 | Data Rate Theorem | `test_data_rate_theorem`, `qc_data_rate_theoretical_min` | Theorem 3 |
| 2 | Sector Bound Theorem | `test_log_quantizer_sector`, `qc_sector_bound_analyze` | Theorem 2 |
| 3 | Zoom Convergence Theorem | `test_zoom_verification`, `qc_zoom_verify_conclusion` | Theorem 4 |
| 4 | Ultimate Bound Theorem | `qc_ultimate_bound_estimate` | Theorem 5 |
| 5 | Quantization Noise: σ²_q = Δ²/12 | `test_quantization_noise` | - |
| 6 | SQNR: 6.02N + 1.76 dB | `test_snr_enob` | - |
| 7 | Shannon-Hartley: C = B·log₂(1+SNR) | `test_channel_capacity` | - |

## L5: Algorithms / Methods (Complete)

| # | Algorithm | Implementation |
|---|-----------|---------------|
| 1 | Uniform quantization | `qc_uniform_quantize_midtread/midrise` |
| 2 | Logarithmic quantization | `qc_log_quantize` |
| 3 | Dynamic zoom quantization | `qc_dyn_quantize` |
| 4 | Lloyd-Max optimal quantizer | `qc_lloyd_max_gaussian` |
| 5 | LBG vector quantizer training | `qc_vec_codebook_train_lbg` |
| 6 | Huffman tree construction | `qc_huffman_build` |
| 7 | Huffman encode/decode | `qc_huffman_encode/decode` |
| 8 | Arithmetic encode/decode | `qc_arith_encode_symbol/decode_symbol` |
| 9 | Delta modulation | `qc_delta_mod_encode/decode` |
| 10 | DPCM encode/decode | `qc_dpcm_encode/decode` |
| 11 | Reverse water-filling | `qc_reverse_waterfill` |
| 12 | QR eigenvalue algorithm | `qc_matrix_eigenvalues` |
| 13 | Lyapunov equation solver | `qc_solve_lyapunov` |
| 14 | Run-length encoding | `qc_runlength_encode/decode` |

## L6: Canonical Problems (Complete)

| # | Problem | Example/Test |
|---|---------|-------------|
| 1 | Uniform quantizer performance analysis | `example1_uniform_quantizer.c` |
| 2 | Logarithmic quantizer sector bound analysis | `example2_log_quantizer.c` |
| 3 | Data rate theorem verification | `example3_data_rate_theorem.c` |
| 4 | Quantized LQR stabilization | `test_quantized_lqr` |
| 5 | Encoder/decoder roundtrip | `test_encoder_decoder` |
| 6 | Huffman coding optimality | `test_huffman_coding` |
| 7 | Dynamic zoom quantizer | `test_dynamic_zoom` |

## L7: Applications (Partial+)

| # | Application | Evidence |
|---|------------|----------|
| 1 | Networked control over bandlimited channels | Data rate theorem integration |
| 2 | Embedded control with finite wordlength | Uniform/log quantizer implementations |
| 3 | Sensor networks with quantized measurements | Output feedback simulation |
| 4 | Teleoperation with limited feedback | Encoder/decoder infrastructure |

## L8: Advanced Topics (Partial+)

| # | Topic | Implementation |
|---|-------|---------------|
| 1 | Anytime capacity | `qc_anytime_capacity` |
| 2 | Entropy-constrained quantization | `qc_entropy_constrained_step` |
| 3 | Adaptive delta modulation | `qc_delta_mod_adapt` |
| 4 | Reverse water-filling | `qc_reverse_waterfill` |
| 5 | Circle criterion for sector-bounded nonlinearity | `qc_circle_criterion_check` |

## L9: Research Frontiers (Partial)

| # | Topic | Documentation |
|---|-------|--------------|
| 1 | Information-theoretic control limits | Data rate theorem extensions |
| 2 | Meta-complexity of quantized stabilization | Zoom strategy complexity analysis |
| 3 | Quantum-limited measurement in control | Future extension |
