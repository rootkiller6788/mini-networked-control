/**
 * @file    qc_encoder.h
 * @brief   Encoding/decoding structures for quantized control data transmission
 *
 * Implements encoding schemes for transmitting quantized values over
 * digital channels: fixed-length, Huffman, arithmetic coding,
 * delta modulation, DPCM, and run-length encoding.
 *
 * Key references:
 *   - Huffman, D.A. (1952). Minimum-redundancy codes. Proc. IRE.
 *   - Witten, Neal, Cleary (1987). Arithmetic coding. CACM.
 *   - Sayood, K. (2017). Introduction to Data Compression. 5th ed.
 */

#ifndef QC_ENCODER_H
#define QC_ENCODER_H

#include "quantized_control.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct QC_HuffNode {
    int             symbol;
    double          probability;
    struct QC_HuffNode *left;
    struct QC_HuffNode *right;
    int             is_leaf;
} QCHuffNode;

typedef struct {
    QCHuffNode     *root;
    int             num_symbols;
    int            *codes;
    int            *code_lengths;
    double          avg_code_length;
    double          entropy;
    double          efficiency;
} QCHuffmanTree;

typedef struct {
    uint8_t        *data;
    int             capacity;
    int             byte_pos;
    int             bit_pos;
    int             total_bits;
} QCBitWriter;

typedef struct {
    const uint8_t  *data;
    int             size;
    int             byte_pos;
    int             bit_pos;
    int             total_bits_read;
} QCBitReader;

typedef struct {
    uint32_t        low;
    uint32_t        high;
    uint32_t        range;
    uint32_t        pending_bits;
    QCBitWriter    *writer;
    int            *cumulative_counts;
    int             total_count;
    int             num_symbols;
} QCArithEncoder;

typedef struct {
    uint32_t        low;
    uint32_t        high;
    uint32_t        code;
    QCBitReader    *reader;
    int            *cumulative_counts;
    int             total_count;
    int             num_symbols;
} QCArithDecoder;

typedef struct {
    double          step_size;
    double          step_min;
    double          step_max;
    double          adaptation_factor;
    double          predicted_value;
    int             initialized;
} QCDeltaModulator;

typedef struct {
    int             order;
    double         *coeffs;
    double         *history;
    int             history_pos;
    QCQuantizer     residual_quantizer;
    int             initialized;
    int             history_filled;
} QCDPCMEncoder;

/* Bitstream operations */
void qc_bit_writer_init(QCBitWriter *bw, int initial_capacity);
int qc_bit_writer_write_bit(QCBitWriter *bw, int bit);
int qc_bit_writer_write_bits(QCBitWriter *bw, uint64_t value, int n);
int qc_bit_writer_flush(QCBitWriter *bw);
void qc_bit_writer_free(QCBitWriter *bw);

void qc_bit_reader_init(QCBitReader *br, const uint8_t *data, int size);
int qc_bit_reader_read_bit(QCBitReader *br);
int qc_bit_reader_read_bits(QCBitReader *br, uint64_t *value, int n);
int qc_bit_reader_has_more(const QCBitReader *br);

/* Huffman coding */
int qc_huffman_build(QCHuffmanTree *tree, const double *probs, const int *symbols, int num_symbols);
int qc_huffman_encode(const QCHuffmanTree *tree, int symbol, QCBitWriter *bw);
int qc_huffman_decode(const QCHuffmanTree *tree, QCBitReader *br, int *symbol);
double qc_huffman_avg_length(const QCHuffmanTree *tree);
double qc_huffman_efficiency(QCHuffmanTree *tree);
void qc_huffman_free(QCHuffmanTree *tree);

/* Arithmetic coding */
int qc_arith_encoder_init(QCArithEncoder *ae, QCBitWriter *bw, const int *freqs, int num_symbols);
int qc_arith_encode_symbol(QCArithEncoder *ae, int symbol);
int qc_arith_encoder_finish(QCArithEncoder *ae);
int qc_arith_decoder_init(QCArithDecoder *ad, QCBitReader *br, const int *freqs, int num_symbols);
int qc_arith_decode_symbol(QCArithDecoder *ad, int *symbol);

/* Delta modulation */
void qc_delta_mod_init(QCDeltaModulator *dm, double step_init, double step_min, double step_max);
int qc_delta_mod_encode(QCDeltaModulator *dm, double sample);
double qc_delta_mod_decode(QCDeltaModulator *dm, int bit);
void qc_delta_mod_adapt(QCDeltaModulator *dm, int bit);
void qc_delta_mod_reset(QCDeltaModulator *dm);

/* DPCM */
int qc_dpcm_encoder_init(QCDPCMEncoder *dpcm, int order, const double *coeffs, int quant_bits);
int qc_dpcm_encode(QCDPCMEncoder *dpcm, double sample, QCBitWriter *bw);
double qc_dpcm_decode(QCDPCMEncoder *dpcm, QCBitReader *br);
void qc_dpcm_encoder_free(QCDPCMEncoder *dpcm);

/* Utility */
double qc_entropy_discrete(const double *probs, int num_symbols);
double qc_differential_entropy_gaussian(double sigma);
double qc_entropy_constrained_step(double sigma, double rate);
int qc_runlength_encode(const int *quantized_indices, int len, QCBitWriter *bw);
int qc_runlength_decode(QCBitReader *br, int *indices, int max_len, int *decoded_len);

#ifdef __cplusplus
}
#endif

#endif /* QC_ENCODER_H */
