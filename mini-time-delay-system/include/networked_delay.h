#ifndef NETWORKED_DELAY_H
#define NETWORKED_DELAY_H

#include "time_delay_system.h"

/* ============================================================================
 * Networked Control with Time Delays
 *
 * Reference:
 *   J. P. Hespanha, P. Naghshtabrizi, Y. Xu, "A Survey of Recent Results
 *     in Networked Control Systems" (Proc. IEEE, 2007)
 *   W. Zhang, M. S. Branicky, S. M. Phillips, "Stability of Networked
 *     Control Systems" (IEEE CSM, 2001)
 *   G. C. Walsh, H. Ye, L. G. Bushnell, "Stability analysis of networked
 *     control systems" (IEEE CST, 2002)
 *   L. Schenato, B. Sinopoli, M. Franceschetti, K. Poolla, S. Sastry,
 *     "Foundations of Control and Estimation Over Lossy Networks" (Proc. IEEE, 2007)
 *   Y. Tipsuwan & M.-Y. Chow, "Control methodologies in networked control
 *     systems" (CEP, 2003)
 *
 * Level 6 — Canonical Problem: NCS with Delay
 * Level 7 — Applications: networked, cyber-physical
 * ============================================================================ */

/* ============================================================================
 * Network Delay Model
 * ============================================================================ */

/* Network protocol model */
typedef enum {
    NET_PROTO_UDP = 0,        /* No retransmission, variable delay */
    NET_PROTO_TCP = 1,        /* Retransmission, larger delay variance */
    NET_PROTO_CAN = 2,        /* Automotive CAN bus (deterministic) */
    NET_PROTO_ETHERNET = 3,   /* Standard Ethernet (CSMA/CD) */
    NET_PROTO_WIFI = 4,       /* Wireless (random access) */
    NET_PROTO_FLEXRAY = 5     /* Time-triggered automotive */
} NetworkProtocol;

/* Network Quality of Service */
typedef struct {
    double bandwidth_bps;         /* Available bandwidth (bits/sec) */
    double packet_loss_rate;      /* p ∈ [0, 1] */
    double mean_delay;            /* E[τ] — average round-trip delay */
    double delay_jitter;          /* σ_τ — delay standard deviation */
    double max_delay;             /* τ_max — worst-case delay */
    double min_delay;             /* τ_min — best-case delay */
    double dropout_rate;          /* Consecutive packet dropout probability */
    int max_consecutive_dropouts; /* Maximum consecutive packet losses */
    bool is_time_triggered;       /* Time-triggered vs event-triggered */
} NetworkQoS;

/* ============================================================================
 * Networked Control System (NCS)
 * ============================================================================ */

typedef struct {
    /* Plant */
    TimeDelaySystem* plant;

    /* Controller (discrete-time, with sampling period h) */
    double h;                     /* Sampling period */
    double* K;                    /* Control gain matrix (m×n) */
    double* x_hat;                /* Estimated state (if using observer) */
    double* u;                    /* Current control input */
    double* y;                    /* Latest measurement */

    /* Network parameters — sensor-to-controller */
    NetworkQoS sc_channel;

    /* Network parameters — controller-to-actuator */
    NetworkQoS ca_channel;

    /* Total round-trip delay τ_total = τ_sc + τ_ca + τ_c */
    double computational_delay;   /* τ_c */

    /* Delay compensation method */
    enum {
        NCS_NO_COMPENSATION = 0,        /* Ignore delay (may go unstable) */
        NCS_SMITH_PREDICTOR = 1,        /* Smith predictor compensation */
        NCS_QUEUE_BUFFER = 2,           /* Buffer at actuator, use latest */
        NCS_TIME_STAMP = 3,             /* Timestamp-based LQR */
        NCS_PREDICTOR_FEEDBACK = 4      /* Krstic predictor feedback */
    } compensation_method;

    /* Timestamp-based control variables */
    double* timestamp_buffer;     /* Buffered (u, timestamp) pairs */
    int ts_buffer_size;
    int ts_write_idx;

    /* Performance metrics */
    double settling_time;
    double overshoot;
    double ise;
    double max_control_effort;

    /* NCS-specific */
    int n_packets_sent;
    int n_packets_lost;
    int n_packets_delayed;        /* Exceeded max allowable delay */
    double avg_actual_delay;

} NetworkedControlSystem;

/* ============================================================================
 * NCS API
 * ============================================================================ */

/* Create NCS with given plant and sampling period */
NetworkedControlSystem* ncs_create(TimeDelaySystem* plant,
                                    double h, double comp_delay);

/* Set the controller gain */
void ncs_set_gain(NetworkedControlSystem* ncs, const double* K);

/* Configure sensor-to-controller channel QoS */
void ncs_set_sc_qos(NetworkedControlSystem* ncs,
                     double bandwidth, double loss_rate,
                     double mean_delay, double jitter,
                     double min_delay, double max_delay);

/* Configure controller-to-actuator channel QoS */
void ncs_set_ca_qos(NetworkedControlSystem* ncs,
                     double bandwidth, double loss_rate,
                     double mean_delay, double jitter,
                     double min_delay, double max_delay);

/* Set delay compensation method */
void ncs_set_compensation(NetworkedControlSystem* ncs, int method);

/* Execute one NCS control cycle:
 *  1. Sample plant → y[k]
 *  2. Transmit y[k] over SC network (adds delay)
 *  3. Compute control u[k] at controller
 *  4. Transmit u[k] over CA network (adds delay)
 *  5. Apply u[k] to plant
 *  6. Simulate plant one sampling period */
void ncs_step(NetworkedControlSystem* ncs);

/* Run NCS for N steps */
void ncs_run(NetworkedControlSystem* ncs, int N_steps,
             double* out_t, double* out_y, double* out_u);

/* Get current plant state */
const double* ncs_get_state(const NetworkedControlSystem* ncs);

/* Get performance statistics */
void ncs_get_stats(const NetworkedControlSystem* ncs,
                   double* ise, double* overshoot,
                   double* settling_time);

/* Free NCS */
void ncs_free(NetworkedControlSystem* ncs);

/* Print NCS status */
void ncs_print(const NetworkedControlSystem* ncs);

/* ============================================================================
 * Queue-Induced Delay Model
 * ============================================================================ */

/* M/M/1 queue model for network-induced delay:
 *
 * τ_queue(t) follows from queue dynamics:
 *   q̇(t) = λ - μ·C(q)  (if queue nonempty)
 * where λ = packet arrival rate, μ = service rate.
 *
 * In steady state: E[τ_queue] = 1/(μ - λ) for M/M/1. */

typedef struct {
    double lambda;          /* Packet arrival rate (packets/sec) */
    double mu;              /* Service rate (packets/sec) */
    double queue_length;    /* Current queue length */
    double max_queue;       /* Buffer size */
    double util;            /* ρ = λ/μ (current utilization) */
    double mean_delay;      /* E[τ] = 1/(μ - λ) */
    double delay_variance;  /* Var[τ] = 1/(μ - λ)² */
} MM1Queue;

/* Create M/M/1 queue model */
MM1Queue* mm1_create(double lambda, double mu, double max_queue);

/* Step the queue forward by dt seconds.
 * Returns the delay experienced by a packet arriving now. */
double mm1_step(MM1Queue* q, double dt, int packets_arrived);

/* Get current queue statistics */
double mm1_mean_delay(const MM1Queue* q);
double mm1_loss_probability(const MM1Queue* q);   /* P(queue full) */
double mm1_utilization(const MM1Queue* q);

/* Free queue model */
void mm1_free(MM1Queue* q);

/* ============================================================================
 * Timestamp-Based LQG Control over Networks
 * ============================================================================ */

/* Optimal LQG controller for NCS with random delays modeled by
 * a Markov chain (Schenato et al., 2007).
 *
 * Uses Kalman filter with intermittent observations and
 * LQR with delayed actuation. */

typedef struct {
    /* System matrices */
    double* A; double* B; double* C;
    double* Q_kf;          /* Process noise covariance */
    double* R_kf;          /* Measurement noise covariance */
    double* Q_lqr;         /* State penalty */
    double* R_lqr;         /* Control penalty */

    /* Kalman filter state */
    double* x_hat;         /* Estimated state */
    double* P;             /* Error covariance */
    int n, m, p;

    /* Delay handling */
    double* u_buffer;      /* Control buffer for delayed actuation */
    int buffer_size;
    int buffer_head;

    /* Arrival indicators */
    bool measurement_arrived;
    bool control_applied;
    double current_delay;

    double dt;
} TimestampLQG;

TimestampLQG* tslqg_create(int n, int m, int p,
                            const double* A, const double* B, const double* C,
                            const double* Q_kf, const double* R_kf,
                            const double* Q_lqr, const double* R_lqr,
                            double dt, int delay_buffer_size);

void tslqg_predict(TimestampLQG* lqg);
void tslqg_update(TimestampLQG* lqg, const double* y, bool packet_arrived);
void tslqg_compute_control(TimestampLQG* lqg, double* u);
const double* tslqg_get_estimate(const TimestampLQG* lqg);
void tslqg_free(TimestampLQG* lqg);

/* ============================================================================
 * TCP/AQM Congestion Control Delay Model (Misra, Gong, Towsley, 2000)
 * ============================================================================ */

/* TCP Reno + RED AQM fluid model:
 *   Ẇ(t) = 1/τ(t) - W(t)·W(t-τ(t)) / (2τ(t-τ(t))) · p(t-τ(t))
 *   q̇(t) = N(t)·W(t)/τ(t) - C
 * where W = TCP window size, q = queue length, τ = q/C + T_p,
 * N = number of flows, C = link capacity, p = drop probability. */

typedef struct {
    double W;                /* TCP congestion window (packets) */
    double q;                /* Queue length (packets) */
    double N;                /* Number of TCP flows */
    double C;                /* Link capacity (packets/sec) */
    double T_p;              /* Propagation delay (sec) */
    double p;                /* Drop/mark probability */
    double tau;              /* RTT = q/C + T_p */
    /* History for delayed variables */
    double* W_history;
    double* p_history;
    int hist_size;
    int hist_idx;
    double dt;
} TCPAQMModel;

TCPAQMModel* tcp_aqm_create(double N, double C, double T_p, double W0);
void tcp_aqm_step(TCPAQMModel* model, double* out_W, double* out_q);
void tcp_aqm_set_drop_prob(TCPAQMModel* model, double p);
double tcp_aqm_rtt(const TCPAQMModel* model);
void tcp_aqm_free(TCPAQMModel* model);

#endif /* NETWORKED_DELAY_H */
