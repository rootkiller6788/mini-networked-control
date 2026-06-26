/**
 * blc_encoding.h — Encoding and Quantization Strategies for Control
 *
 * This header implements practical encoding schemes for bandwidth-limited
 * control systems, building on the theoretical foundations established
 * by the Data Rate Theorem.
 *
 * Encoder-decoder pairs must be co-designed with the controller because
 * the classical separation principle fails under communication constraints.
 * The encoder must allocate bits to convey state information that is most
 * critical for stabilization — typically the unstable subspace.
 *
 * Key encoding strategies covered:
 *
 * 1. Uniform Scalar Quantization — simplest, each component independently
 *    quantized with equal bit allocation.
 *
 * 2. Logarithmic Quantization — quantizer levels geometrically spaced,
 *    optimal for stabilizing unstable systems because it provides finer
 *    resolution near the origin (Elia & Mitter, 2001).
 *
 * 3. Lloyd-Max Quantization — minimum mean-square error quantizer for a
 *    given source distribution (Lloyd, 1957; Max, 1960).
 *
 * 4. Vector Quantization — jointly quantize multiple state components,
 *    achieving the rate-distortion bound for vector sources.
 *
 * 5. Differential / Predictive Coding — transmit quantized innovation
 *    ε = x - ẍ̂_pred, exploiting temporal correlation to reduce entropy.
 *
 * 6. Entropy Coding — Huffman or arithmetic coding to approach entropy
 *    limit after quantization.
 *
 * 7. Unequal Error Protection — allocate more bits (or redundancy) to
 *    unstable subspace components.
 *
 * Knowledge coverage: L5 (Algorithms), L7 (Applications)
 */

#ifndef BLC_ENCODING_H
#define BLC_ENCODING_H

#include "blc_core.h"
#include "blc_datarate.h"

/* ================================================================
 * Logarithmic Quantizer
 *
 * For stabilizing an unstable scalar system with pole λ > 0:
 *   The coarsest stabilizing quantizer is logarithmic with
 *   density δ = e^{λ·T_s} (Elia & Mitter, 2001).
 *
 * Levels:  q_k = ρ^k · q_0   for k = 0, ±1, ±2, ...
 * where ρ = δ = e^{λ·T_s} is the optimal density.
 *
 * A logarithmic quantizer with N levels can stabilize a system
 * with maximum unstable eigenvalue:
 *   λ_max = (1/T_s) · ln(ρ)  with ρ = (1+Δ)/(1-Δ), Δ = 1/N
 * ================================================================ */
typedef struct {
    int     num_levels_pos;    /** Number of positive levels */
    int     num_levels_neg;    /** Number of negative levels */
    double  density;           /** Quantizer density δ (e.g., e^{λT_s}) */
    double  deadzone;          /** Deadzone around zero (smallest step) */
    double* level_table;       /** Pre-computed quantization levels */
    double* recon_table;       /** Pre-computed reconstruction values */
    double  overload_value;    /** Saturation value */
} BLCLogQuantizer;

/** Lloyd-Max optimal quantizer (iterative design) */
typedef struct {
    int      num_levels;       /** Number of reconstruction levels */
    double*  boundaries;       /** Decision boundaries (N+1 values) */
    double*  reconstructions;  /** Reconstruction levels (N values) */
    double   mse;              /** Achieved mean squared error */
    int      iterations;       /** Iterations to converge */
} BLCLloydQuantizer;

/** Vector quantizer using LBG algorithm (Linde-Buzo-Gray, 1980) */
typedef struct {
    int      dimension;        /** Vector dimension k */
    int      codebook_size;    /** Number of codewords N */
    double** codebook;         /** N × k matrix of codewords */
    double   mse;              /** Quantization distortion */
    double   bits_per_sample;  /** = log₂(N) / k */
    int*     cell_counts;      /** Number of vectors in each Voronoi cell */
} BLCVectorQuantizer;

/** Entropy coder (Huffman) for quantized indices */
typedef struct {
    int      alphabet_size;    /** Number of symbols */
    double*  probabilities;    /** Symbol probability distribution */
    double   entropy;          /** Source entropy H = -Σ pᵢ log₂ pᵢ */
    int*     code_lengths;     /** Code word lengths in bits */
    int      max_code_length;  /** Maximum code word length */
    double   avg_code_length;  /** Average code word length */
    bool     is_optimal;       /** Huffman optimality flag */
} BLCHuffmanCoder;

/** Arithmetic encoder for near-entropy coding */
typedef struct {
    uint64_t low;              /** Current interval lower bound */
    uint64_t high;             /** Current interval upper bound */
    uint64_t range;            /** Current range = high - low */
    int      pending_bits;     /** Bits pending for carry resolution */
    int      total_bits;       /** Total bits emitted */
    int      symbol_count[256];/** Histogram of encoded symbols */
} BLCArithmeticEncoder;

/* ================================================================
 * Logarithmic Quantizer API
 *
 * The logarithmic quantizer is optimal for stabilization:
 * it provides the coarsest quantization that still achieves
 * asymptotic stability (Elia & Mitter, 2001).
 * ================================================================ */

/** Initialize logarithmic quantizer.
 *  @param lq Quantizer to initialize
 *  @param density Quantizer density δ (must be > 1)
 *  @param num_pos Number of positive levels
 *  @param deadzone Smallest quantization step near zero
 *  @return 0 on success */
int     blc_log_quant_init(BLCLogQuantizer* lq, double density,
                            int num_pos, double deadzone);

/** Free logarithmic quantizer resources */
void    blc_log_quant_free(BLCLogQuantizer* lq);

/** Quantize a value with logarithmic quantizer */
int     blc_log_quant_encode(const BLCLogQuantizer* lq, double value);

/** Reconstruct a value from logarithmic quantizer index */
double  blc_log_quant_decode(const BLCLogQuantizer* lq, int index);

/** Compute optimal density for a given unstable eigenvalue.
 *  δ_opt = e^{λ_max · T_s}
 *  @param lambda_max Maximum unstable eigenvalue (real part)
 *  @param Ts Sampling period
 *  @return Optimal quantizer density */
double  blc_log_quant_optimal_density(double lambda_max, double Ts);

/** Compute the minimum number of levels for stabilization.
 *  N_min = ceil( (δ+1) / (δ-1) )
 *  @param density Quantizer density
 *  @return Minimum number of positive levels */
int     blc_log_quant_min_levels(double density);

/* ================================================================
 * Lloyd-Max Quantizer API
 * ================================================================ */

/** Design Lloyd-Max quantizer for a given probability density.
 *  Uses the iterative Lloyd-Max algorithm:
 *    b_i = (r_i + r_{i+1}) / 2  (boundaries = midpoints)
 *    r_i = E[x | x ∈ (b_{i-1}, b_i]] (centroids)
 *  @param q Quantizer structure (pre-allocated with num_levels set)
 *  @param pdf_samples Array of sample values from the source
 *  @param n_samples Number of samples
 *  @param range_lo Lower bound of signal range
 *  @param range_hi Upper bound of signal range
 *  @param max_iter Maximum iterations
 *  @param tol Convergence tolerance
 *  @return Number of iterations performed */
int     blc_lloyd_max_design(BLCLloydQuantizer* q, const double* pdf_samples,
                              int n_samples, double range_lo, double range_hi,
                              int max_iter, double tol);

/** Free Lloyd-Max quantizer resources */
void    blc_lloyd_max_free(BLCLloydQuantizer* q);

/** Encode value with Lloyd-Max quantizer */
int     blc_lloyd_max_encode(const BLCLloydQuantizer* q, double value);

/** Decode from Lloyd-Max index */
double  blc_lloyd_max_decode(const BLCLloydQuantizer* q, int index);

/** Get Lloyd-Max quantizer MSE */
double  blc_lloyd_max_get_mse(const BLCLloydQuantizer* q);

/* ================================================================
 * Vector Quantizer (LBG) API
 * ================================================================ */

/** Initialize vector quantizer */
int     blc_vq_init(BLCVectorQuantizer* vq, int dimension, int codebook_size);

/** Free vector quantizer */
void    blc_vq_free(BLCVectorQuantizer* vq);

/** Train vector quantizer using LBG algorithm.
 *  @param vq Vector quantizer (codebook_size and dimension pre-set)
 *  @param training_data Array of training vectors (n_vectors × dimension)
 *  @param n_vectors Number of training vectors
 *  @param epsilon Convergence threshold (relative distortion change)
 *  @param max_iter Maximum iterations
 *  @return Number of iterations */
int     blc_vq_train(BLCVectorQuantizer* vq, const double* training_data,
                      int n_vectors, double epsilon, int max_iter);

/** Encode vector: find nearest codeword */
int     blc_vq_encode(const BLCVectorQuantizer* vq, const double* x);

/** Decode: return codeword for index */
void    blc_vq_decode(const BLCVectorQuantizer* vq, int index, double* x_hat);

/** Compute quantization distortion for data set */
double  blc_vq_distortion(const BLCVectorQuantizer* vq,
                           const double* data, int n_vectors);

/** Split codebook (double the size) — used in LBG splitting phase */
int     blc_vq_split_codebook(BLCVectorQuantizer* vq, double epsilon);

/** Compute bits per dimension = log₂(codebook_size) / dimension */
double  blc_vq_rate(const BLCVectorQuantizer* vq);

/* ================================================================
 * Entropy Coding API
 * ================================================================ */

/** Build Huffman code from symbol probabilities.
 *  Uses the Huffman algorithm (1952): greedy bottom-up tree construction.
 *  @param huff Coder structure (alphabet_size must be set)
 *  @param probabilities Symbol probabilities (must sum to 1)
 *  @return 0 on success */
int     blc_huffman_build(BLCHuffmanCoder* huff, const double* probabilities);

/** Free Huffman coder resources */
void    blc_huffman_free(BLCHuffmanCoder* huff);

/** Get code word length for a symbol */
int     blc_huffman_get_length(const BLCHuffmanCoder* huff, int symbol);

/** Encode a symbol, returns bit string length (actual encoding
 *  requires bit packing — this returns the theoretical length) */
int     blc_huffman_encode_symbol(const BLCHuffmanCoder* huff,
                                   int symbol, uint8_t* bits, int* offset);

/** Decode symbol from bit stream */
int     blc_huffman_decode_symbol(const BLCHuffmanCoder* huff,
                                   const uint8_t* bits, int* offset,
                                   int total_bits);

/** Compute source entropy: H = -Σ p_i log₂(p_i) */
double  blc_entropy(const double* probabilities, int n);

/** Compute coding efficiency: H / L_avg */
double  blc_huffman_efficiency(const BLCHuffmanCoder* huff);

/* ================================================================
 * Arithmetic Coding API
 * ================================================================ */

/** Initialize arithmetic encoder */
void    blc_arithmetic_init(BLCArithmeticEncoder* ae);

/** Encode a symbol with given cumulative probability.
 *  @param ae Encoder state
 *  @param cum_low Cumulative probability up to (not including) this symbol
 *  @param cum_high Cumulative probability up to and including this symbol
 *  @param total Total probability mass (typically 1.0 or a power of 2) */
int     blc_arithmetic_encode(BLCArithmeticEncoder* ae,
                               double cum_low, double cum_high, double total);

/** Finalize arithmetic coding stream */
int     blc_arithmetic_finalize(BLCArithmeticEncoder* ae);

/** Get number of bits emitted */
int     blc_arithmetic_get_bits(const BLCArithmeticEncoder* ae);

/* ================================================================
 * Bandwidth-Optimal Encoding
 * ================================================================ */

/** Allocate bit budget across state components based on
 *  unstable eigenvalues. States corresponding to larger |λ|
 *  get more bits.
 *  @param eigenvalues Real parts of eigenvalues
 *  @param n Number of states
 *  @param total_bits Total bit budget
 *  @param allocation Output bit allocation per state */
void    blc_allocate_bits_eigenvalue(const double* eigenvalues, int n,
                                      int total_bits, int* allocation);

/** Compute encoding distortion lower bound for given bit rate.
 *  D ≥ σ² · 2^{-2R}  (for Gaussian source with variance σ²)
 *  @param variance Source variance
 *  @param bits_per_sample Bit allocation
 *  @return Lower bound on MSE distortion */
double  blc_encoding_distortion_bound(double variance, double bits_per_sample);

/** Compute required bit rate for target distortion.
 *  R ≥ 0.5 · log₂(σ² / D)  (Shannon lower bound)
 *  @param variance Source variance
 *  @param target_distortion Target MSE
 *  @return Minimum bits per sample */
double  blc_encoding_bitrate_bound(double variance, double target_distortion);

/** Run-length encode for sparse control signals.
 *  When control is sparse (zero most of the time), run-length
 *  encoding can dramatically reduce bandwidth.
 *  @param signal Input signal array
 *  @param n Samples
 *  @param runs Output (value, count) pairs
 *  @param max_runs Maximum number of run entries
 *  @return Number of runs found */
int     blc_runlength_encode(const double* signal, int n,
                              double* run_values, int* run_lengths,
                              int max_runs);

#endif /* BLC_ENCODING_H */