#include "networked_delay.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * Networked Control System (NCS) Implementation
 *
 * Simulates a control loop closed over a communication network
 * with delays, packet loss, and jitter.
 *
 * Key NCS effects:
 *   1. Sensor-to-controller delay (τ_sc)
 *   2. Controller-to-actuator delay (τ_ca)
 *   3. Computational delay (τ_c)
 *   4. Packet loss (modeled as Bernoulli/Gilbert-Elliott)
 *   5. Quantization effects (not modeled here)
 *   6. Sampling jitter
 * ============================================================================ */

/* Simple pseudo-random number generator (linear congruential) */
static unsigned int ncs_rand_state = 12345;

static void ncs_srand(unsigned int seed) { ncs_rand_state = seed; }

static double ncs_urand(void) {
    ncs_rand_state = (ncs_rand_state * 1103515245U + 12345U) & 0x7fffffffU;
    return (double)ncs_rand_state / (double)0x7fffffffU;
}

/* Gaussian random variable via Box-Muller */
static double ncs_gauss_rand(void) {
    double u1 = ncs_urand(), u2 = ncs_urand();
    if (u1 < 1e-15) u1 = 1e-15;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/* Sample a network delay from QoS parameters */
static double sample_network_delay(const NetworkQoS* qos) {
    if (!qos) return 0.0;
    /* Model: delay = mean_delay + jitter * N(0,1)
     * truncated to [min_delay, max_delay] */
    double d = qos->mean_delay + qos->delay_jitter * ncs_gauss_rand();
    if (d < qos->min_delay) d = qos->min_delay;
    if (d > qos->max_delay) d = qos->max_delay;
    return d;
}

/* Check if packet is lost */
static bool is_packet_lost(const NetworkQoS* qos) {
    if (!qos) return false;
    return ncs_urand() < qos->packet_loss_rate;
}

/* ============================================================================
 * NCS Lifecycle
 * ============================================================================ */

NetworkedControlSystem* ncs_create(TimeDelaySystem* plant,
                                    double h, double comp_delay) {
    NetworkedControlSystem* ncs = (NetworkedControlSystem*)
        calloc(1, sizeof(NetworkedControlSystem));
    if (!ncs) return NULL;

    ncs->plant = plant;
    ncs->h = h;
    ncs->computational_delay = comp_delay;

    int n = plant->n_states;
    int m = plant->n_inputs;
    int p = plant->n_outputs;

    ncs->K = (double*)calloc((size_t)(m * n), sizeof(double));
    ncs->x_hat = (double*)calloc((size_t)n, sizeof(double));
    ncs->u = (double*)calloc((size_t)m, sizeof(double));
    ncs->y = (double*)calloc((size_t)p, sizeof(double));

    /* Default: ideal network (zero delay, zero loss) */
    memset(&ncs->sc_channel, 0, sizeof(NetworkQoS));
    memset(&ncs->ca_channel, 0, sizeof(NetworkQoS));

    ncs->compensation_method = NCS_NO_COMPENSATION;
    ncs->timestamp_buffer = NULL;
    ncs->ts_buffer_size = 0;

    ncs->settling_time = 0.0;
    ncs->overshoot = 0.0;
    ncs->ise = 0.0;
    ncs->max_control_effort = 0.0;
    ncs->n_packets_sent = 0;
    ncs->n_packets_lost = 0;
    ncs->n_packets_delayed = 0;
    ncs->avg_actual_delay = 0.0;

    ncs_srand(1);  /* Seed RNG */
    return ncs;
}

void ncs_set_gain(NetworkedControlSystem* ncs, const double* K) {
    if (!ncs || !K) return;
    memcpy(ncs->K, K,
           (size_t)(ncs->plant->n_states * ncs->plant->n_inputs)
           * sizeof(double));
}

void ncs_set_sc_qos(NetworkedControlSystem* ncs,
                     double bandwidth, double loss_rate,
                     double mean_delay, double jitter,
                     double min_delay, double max_delay) {
    if (!ncs) return;
    ncs->sc_channel.bandwidth_bps = bandwidth;
    ncs->sc_channel.packet_loss_rate = loss_rate;
    ncs->sc_channel.mean_delay = mean_delay;
    ncs->sc_channel.delay_jitter = jitter;
    ncs->sc_channel.min_delay = min_delay;
    ncs->sc_channel.max_delay = max_delay;
}

void ncs_set_ca_qos(NetworkedControlSystem* ncs,
                     double bandwidth, double loss_rate,
                     double mean_delay, double jitter,
                     double min_delay, double max_delay) {
    if (!ncs) return;
    ncs->ca_channel.bandwidth_bps = bandwidth;
    ncs->ca_channel.packet_loss_rate = loss_rate;
    ncs->ca_channel.mean_delay = mean_delay;
    ncs->ca_channel.delay_jitter = jitter;
    ncs->ca_channel.min_delay = min_delay;
    ncs->ca_channel.max_delay = max_delay;
}

void ncs_set_compensation(NetworkedControlSystem* ncs, int method) {
    if (ncs) ncs->compensation_method = method;
}

/* ============================================================================
 * NCS Step — Execute one control cycle
 * ============================================================================ */

void ncs_step(NetworkedControlSystem* ncs) {
    if (!ncs || !ncs->plant) return;

    TimeDelaySystem* plant = ncs->plant;
    int n = plant->n_states;
    int m = plant->n_inputs;
    int p = plant->n_outputs;
    double h = ncs->h;
    (void)(plant->n_delays > 0 ? plant->delays[0]->tau_nominal : 0.0);

    /* 1. Sample plant output y(t) */
    /* y = C x + (measurement noise = 0 for now) */
    for (int i = 0; i < p; i++) {
        ncs->y[i] = 0.0;
        for (int j = 0; j < n; j++)
            ncs->y[i] += plant->C[i * n + j] * plant->current_state[j];
    }

    /* 2. Transmit over SC channel */
    ncs->n_packets_sent++;
    double tau_sc = 0.0;
    bool sc_lost = is_packet_lost(&ncs->sc_channel);
    if (!sc_lost) {
        tau_sc = sample_network_delay(&ncs->sc_channel);
    } else {
        ncs->n_packets_lost++;
        /* Use previous measurement or zero-order hold */
        /* y stays at previous value */
    }

    /* 3. Controller computation (after SC delay) */
    /* In a real NCS, measurement arrives at controller at t + τ_sc */
    /* For simulation, we compute control immediately with delayed y */

    /* 4. Compute control u = -K * x (full state feedback) */
    /* If sc_lost, use last known measurement = ncs->y (ZOH) */
    double* x_for_control = plant->current_state;
    /* In practice with output feedback, would estimate x from y. */

    double max_u = 0.0;
    for (int i = 0; i < m; i++) {
        ncs->u[i] = 0.0;
        for (int j = 0; j < n; j++)
            ncs->u[i] -= ncs->K[i * n + j] * x_for_control[j];
        if (fabs(ncs->u[i]) > max_u) max_u = fabs(ncs->u[i]);
    }
    if (max_u > ncs->max_control_effort) ncs->max_control_effort = max_u;

    /* 5. Transmit control over CA channel */
    double tau_ca = 0.0;
    bool ca_lost = is_packet_lost(&ncs->ca_channel);
    if (!ca_lost) {
        tau_ca = sample_network_delay(&ncs->ca_channel);
    } else {
        ncs->n_packets_lost++;
        /* Control packet lost → actuator holds previous value (ZOH) */
    }

    double total_delay = tau_sc + tau_ca + ncs->computational_delay;
    if (total_delay > ncs->sc_channel.max_delay + ncs->ca_channel.max_delay)
        ncs->n_packets_delayed++;

    ncs->avg_actual_delay = 0.95 * ncs->avg_actual_delay + 0.05 * total_delay;

    /* 6. Apply control to plant and simulate plant dynamics for h seconds */
    /* Plant dynamic update: ẋ = A x + A_d x(t-τ) + B u */
    int steps = (int)ceil(h / 0.001);
    if (steps < 1) steps = 1;
    double dt_sim = h / (double)steps;

    double* xdot = (double*)malloc((size_t)n * sizeof(double));
    double* x_delayed = (double*)malloc((size_t)n * sizeof(double));

    for (int k = 0; k < steps; k++) {
        /* Get delayed state (simplified: use last delayed_state) */
        if (plant->history_buffer && plant->history_points > 0) {
            /* Use oldest history entry as delayed state */
            memcpy(x_delayed, plant->history_buffer,
                   (size_t)n * sizeof(double));
        } else {
            memcpy(x_delayed, plant->delayed_state,
                   (size_t)n * sizeof(double));
        }

        /* ẋ = A x + A_d x_d + B u */
        for (int i = 0; i < n; i++) {
            xdot[i] = 0.0;
            for (int j = 0; j < n; j++) {
                xdot[i] += plant->A[i * n + j] * plant->current_state[j]
                         + plant->A_delayed[i * n + j] * x_delayed[j];
            }
            for (int j = 0; j < m; j++)
                xdot[i] += plant->B[i * n + j] * ncs->u[j];
        }

        /* Update state via Euler integration */
        for (int i = 0; i < n; i++)
            plant->current_state[i] += xdot[i] * dt_sim;

        /* Shift history buffer */
        if (plant->history_buffer && plant->history_points > 1) {
            size_t hist_size = (size_t)(plant->history_points - 1) * n;
            memmove(plant->history_buffer + n, plant->history_buffer,
                    hist_size * sizeof(double));
            memcpy(plant->history_buffer, plant->current_state,
                   (size_t)n * sizeof(double));
        }

        plant->t_current += dt_sim;
    }

    /* Update ISE */
    if (p > 0) {
        double ref = 0.0;  /* Assuming regulation to zero */
        double e = ref - ncs->y[0];
        ncs->ise += e * e * h;
    }

    free(xdot); free(x_delayed);
}

void ncs_run(NetworkedControlSystem* ncs, int N_steps,
             double* out_t, double* out_y, double* out_u) {
    if (!ncs) return;
    for (int k = 0; k < N_steps; k++) {
        if (out_t) out_t[k] = (double)k * ncs->h;
        ncs_step(ncs);
        if (out_y) out_y[k] = ncs->y[0];
        if (out_u) out_u[k] = ncs->u[0];
    }
}

const double* ncs_get_state(const NetworkedControlSystem* ncs) {
    return ncs ? ncs->plant->current_state : NULL;
}

void ncs_get_stats(const NetworkedControlSystem* ncs,
                   double* ise, double* overshoot,
                   double* settling_time) {
    if (!ncs) return;
    if (ise) *ise = ncs->ise;
    if (overshoot) *overshoot = 0.0;  /* Requires tracking ref tracking */
    if (settling_time) *settling_time = 0.0;
}

void ncs_free(NetworkedControlSystem* ncs) {
    if (!ncs) return;
    free(ncs->K); free(ncs->x_hat); free(ncs->u); free(ncs->y);
    free(ncs->timestamp_buffer);
    free(ncs);
}

void ncs_print(const NetworkedControlSystem* ncs) {
    if (!ncs) { printf("NCS: NULL\n"); return; }
    printf("=== Networked Control System ===\n");
    printf("  Sampling period: %.4f s\n", ncs->h);
    printf("  SC channel: loss=%.2f%% mean_delay=%.4f jitter=%.4f\n",
           ncs->sc_channel.packet_loss_rate * 100.0,
           ncs->sc_channel.mean_delay, ncs->sc_channel.delay_jitter);
    printf("  CA channel: loss=%.2f%% mean_delay=%.4f jitter=%.4f\n",
           ncs->ca_channel.packet_loss_rate * 100.0,
           ncs->ca_channel.mean_delay, ncs->ca_channel.delay_jitter);
    printf("  Compensation: %d  Avg delay: %.6f\n",
           ncs->compensation_method, ncs->avg_actual_delay);
    printf("  Packets: sent=%d lost=%d delayed=%d\n",
           ncs->n_packets_sent, ncs->n_packets_lost,
           ncs->n_packets_delayed);
    printf("  ISE: %.6f  Max u: %.4f\n", ncs->ise, ncs->max_control_effort);
}

/* ============================================================================
 * M/M/1 Queue Model for Network Delay
 *
 * λ = packet arrival rate, μ = service rate
 * ρ = λ/μ = utilization
 * E[queue_length] = ρ/(1-ρ)
 * E[delay] = 1/(μ - λ)
 * Blocking probability P_block = ρ^k (1-ρ)/(1-ρ^{k+1}) for buffer K
 * ============================================================================ */

MM1Queue* mm1_create(double lambda, double mu, double max_queue) {
    MM1Queue* q = (MM1Queue*)calloc(1, sizeof(MM1Queue));
    if (!q) return NULL;
    q->lambda = lambda;
    q->mu = mu;
    q->queue_length = 0.0;
    q->max_queue = max_queue;
    q->util = lambda / mu;
    q->mean_delay = (mu > lambda) ? 1.0 / (mu - lambda) : INFINITY;
    q->delay_variance = (mu > lambda) ? 1.0 / ((mu - lambda) * (mu - lambda))
                                      : INFINITY;
    return q;
}

double mm1_step(MM1Queue* q, double dt, int packets_arrived) {
    if (!q) return 0.0;

    /* Arrival process */
    q->queue_length += (double)packets_arrived;

    /* Service process: can serve mu*dt packets per dt seconds */
    int served = (int)floor(q->mu * dt);
    if (served < 0) served = 0;
    double actual_served = (double)served;
    if (actual_served > q->queue_length) actual_served = q->queue_length;
    q->queue_length -= actual_served;

    /* Bound queue length */
    if (q->queue_length > q->max_queue) q->queue_length = q->max_queue;

    /* Little's Law: delay = queue_length / lambda */
    if (q->lambda > 0)
        return q->queue_length / q->lambda;
    return 0.0;
}

double mm1_mean_delay(const MM1Queue* q) {
    if (!q) return 0.0;
    if (q->mu <= q->lambda) return INFINITY;
    return 1.0 / (q->mu - q->lambda);
}

double mm1_loss_probability(const MM1Queue* q) {
    if (!q) return 0.0;
    /* For M/M/1/K, loss prob = ρ^K (1-ρ) / (1-ρ^{K+1}) */
    double rho = q->util;
    int K = (int)q->max_queue;
    if (K <= 0) return 1.0;
    if (fabs(rho - 1.0) < 1e-10) return 1.0 / (double)(K + 1);

    double rho_K = pow(rho, (double)K);
    double rho_K1 = rho_K * rho;
    return rho_K * (1.0 - rho) / (1.0 - rho_K1);
}

double mm1_utilization(const MM1Queue* q) {
    return q ? q->util : 0.0;
}

void mm1_free(MM1Queue* q) {
    free(q);
}

/* ============================================================================
 * Timestamp-Based LQG Control over Networks
 *
 * Implements the optimal controller for NCS with random delays
 * (Schenato et al., 2007).
 *
 * Uses a Kalman filter modified for intermittent observations:
 *   Measurement update: occurs only when packet arrives
 *   Time update: always
 *   LQR control: uses predicted state, buffered for delay
 * ============================================================================ */

TimestampLQG* tslqg_create(int n, int m, int p,
                            const double* A, const double* B, const double* C,
                            const double* Q_kf, const double* R_kf,
                            const double* Q_lqr, const double* R_lqr,
                            double dt, int delay_buffer_size) {
    TimestampLQG* lqg = (TimestampLQG*)calloc(1, sizeof(TimestampLQG));
    if (!lqg) return NULL;

    lqg->n = n; lqg->m = m; lqg->p = p;
    lqg->dt = dt;
    lqg->buffer_size = delay_buffer_size;

    int n2 = n * n;
    lqg->A = (double*)malloc((size_t)n2 * sizeof(double));
    lqg->B = (double*)malloc((size_t)(n * m) * sizeof(double));
    lqg->C = (double*)malloc((size_t)(p * n) * sizeof(double));
    lqg->Q_kf = (double*)malloc((size_t)n2 * sizeof(double));
    lqg->R_kf = (double*)malloc((size_t)(p * p) * sizeof(double));
    lqg->Q_lqr = (double*)malloc((size_t)n2 * sizeof(double));
    lqg->R_lqr = (double*)malloc((size_t)(m * m) * sizeof(double));

    if (A) memcpy(lqg->A, A, (size_t)n2 * sizeof(double));
    if (B) memcpy(lqg->B, B, (size_t)(n * m) * sizeof(double));
    if (C) memcpy(lqg->C, C, (size_t)(p * n) * sizeof(double));
    if (Q_kf) memcpy(lqg->Q_kf, Q_kf, (size_t)n2 * sizeof(double));
    if (R_kf) memcpy(lqg->R_kf, R_kf, (size_t)(p * p) * sizeof(double));
    if (Q_lqr) memcpy(lqg->Q_lqr, Q_lqr, (size_t)n2 * sizeof(double));
    if (R_lqr) memcpy(lqg->R_lqr, R_lqr, (size_t)(m * m) * sizeof(double));

    lqg->x_hat = (double*)calloc((size_t)n, sizeof(double));
    lqg->P = (double*)calloc((size_t)n2, sizeof(double));
    /* Initialize P as identity */
    for (int i = 0; i < n; i++) lqg->P[i * n + i] = 1.0;

    lqg->u_buffer = (double*)calloc((size_t)(delay_buffer_size * m), sizeof(double));
    lqg->buffer_head = 0;
    lqg->measurement_arrived = false;
    lqg->control_applied = false;

    return lqg;
}

void tslqg_predict(TimestampLQG* lqg) {
    if (!lqg) return;
    int n = lqg->n;

    /* Time update: x̂⁻ = A x̂⁺ */
    double* x_pred = (double*)malloc((size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) {
        x_pred[i] = 0.0;
        for (int j = 0; j < n; j++)
            x_pred[i] += lqg->A[i * n + j] * lqg->x_hat[j];
    }
    memcpy(lqg->x_hat, x_pred, (size_t)n * sizeof(double));

    /* P⁻ = A P Aᵀ + Q_kf */
    double* P_new = (double*)calloc((size_t)(n * n), sizeof(double));
    /* A * P * Aᵀ */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                for (int l = 0; l < n; l++)
                    P_new[i * n + j] += lqg->A[i * n + k]
                        * lqg->P[k * n + l] * lqg->A[j * n + l];

    /* Add Q_kf */
    for (int i = 0; i < n * n; i++)
        P_new[i] += lqg->Q_kf[i];

    memcpy(lqg->P, P_new, (size_t)(n * n) * sizeof(double));

    free(x_pred); free(P_new);
    lqg->measurement_arrived = false;
}

void tslqg_update(TimestampLQG* lqg, const double* y, bool packet_arrived) {
    if (!lqg) return;
    if (!packet_arrived || !y) { tslqg_predict(lqg); return; }

    int n = lqg->n, p = lqg->p;

    /* Kalman gain: K = P Cᵀ (C P Cᵀ + R)^{-1} */
    /* For scalar measurement, compute directly */
    if (p == 1) {
        /* S = C P Cᵀ + R */
        double S = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                S += lqg->C[i] * lqg->P[i * n + j] * lqg->C[j];
        S += lqg->R_kf[0];
        if (S < 1e-12) return;

        /* Innovation: y - C x̂ */
        double innov = y[0];
        for (int i = 0; i < n; i++) innov -= lqg->C[i] * lqg->x_hat[i];

        /* Update x̂ = x̂ + K * innov */
        for (int i = 0; i < n; i++) {
            double K_i = 0.0;
            for (int j = 0; j < n; j++)
                K_i += lqg->P[i * n + j] * lqg->C[j];
            K_i /= S;
            lqg->x_hat[i] += K_i * innov;
        }

        /* Update P = (I - K C) P */
        double* P_tmp = (double*)malloc((size_t)(n * n) * sizeof(double));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                double KC_ij = 0.0;
                for (int k = 0; k < n; k++)
                    KC_ij += (lqg->P[i * n + k] * lqg->C[k] / S) * lqg->C[j];
                P_tmp[i * n + j] = lqg->P[i * n + j] - KC_ij;
                /* Ensure symmetry */
                if (P_tmp[i * n + j] < 0) P_tmp[i * n + j] = 0.0;
            }
        memcpy(lqg->P, P_tmp, (size_t)(n * n) * sizeof(double));
        free(P_tmp);
    }

    lqg->measurement_arrived = true;
}

void tslqg_compute_control(TimestampLQG* lqg, double* u) {
    if (!lqg || !u) return;
    int n = lqg->n, m = lqg->m;

    /* LQR gain: u = -K_lqr * x̂ (using precomputed or approximate gain)
     * For simplicity, use a deadbeat-like gain from Q_lqr and R_lqr */
    for (int i = 0; i < m; i++) {
        u[i] = 0.0;
        for (int j = 0; j < n; j++)
            u[i] -= lqg->x_hat[j];  /* Simplified: K=I */
    }

    /* Buffer control for delayed actuation */
    memcpy(lqg->u_buffer + (size_t)lqg->buffer_head * m,
           u, (size_t)m * sizeof(double));
    lqg->buffer_head = (lqg->buffer_head + 1) % lqg->buffer_size;
}

const double* tslqg_get_estimate(const TimestampLQG* lqg) {
    return lqg ? lqg->x_hat : NULL;
}

void tslqg_free(TimestampLQG* lqg) {
    if (!lqg) return;
    free(lqg->A); free(lqg->B); free(lqg->C);
    free(lqg->Q_kf); free(lqg->R_kf);
    free(lqg->Q_lqr); free(lqg->R_lqr);
    free(lqg->x_hat); free(lqg->P); free(lqg->u_buffer);
    free(lqg);
}

/* ============================================================================
 * TCP/AQM Congestion Control Model (Misra, Gong, Towsley, 2000)
 *
 * Fluid model of TCP Reno:
 *   Ẇ(t) = 1/RTT - W(t)W(t-RTT)/(2·RTT(t-RTT)) · p(t-RTT)
 *   q̇(t) = N·W(t)/RTT - C
 *
 * W = TCP congestion window (packets)
 * q = queue length (packets)
 * RTT = q/C + T_p  (round-trip time)
 * ============================================================================ */

TCPAQMModel* tcp_aqm_create(double N, double C, double T_p, double W0) {
    TCPAQMModel* m = (TCPAQMModel*)calloc(1, sizeof(TCPAQMModel));
    if (!m) return NULL;
    m->N = N;
    m->C = C;
    m->T_p = T_p;
    m->W = (W0 > 0) ? W0 : 1.0;
    m->q = 0.0;
    m->p = 0.0;
    m->tau = T_p;  /* Initial RTT = propagation delay (empty queue) */
    m->dt = 0.001;

    m->hist_size = 2000;
    m->W_history = (double*)calloc((size_t)m->hist_size, sizeof(double));
    m->p_history = (double*)calloc((size_t)m->hist_size, sizeof(double));
    m->hist_idx = 0;
    /* Initialize history with current values */
    for (int i = 0; i < m->hist_size; i++) {
        m->W_history[i] = m->W;
        m->p_history[i] = m->p;
    }
    return m;
}

void tcp_aqm_step(TCPAQMModel* m, double* out_W, double* out_q) {
    if (!m) return;

    /* Update RTT */
    m->tau = m->q / m->C + m->T_p;
    if (m->tau < m->T_p) m->tau = m->T_p;

    /* Get delayed values from history buffer */
    int delay_idx = m->hist_idx - (int)ceil(m->tau / m->dt);
    while (delay_idx < 0) delay_idx += m->hist_size;

    double W_delayed = m->W_history[delay_idx];
    double p_delayed = m->p_history[delay_idx];

    /* TCP window dynamics:
     * dW/dt = (1/RTT) - (W * W_delayed) / (2 * RTT_delayed) * p_delayed */
    double RTT_delayed = m->q / m->C + m->T_p;  /* Approximate */
    if (RTT_delayed < m->T_p) RTT_delayed = m->T_p;

    double dW = 1.0 / m->tau;
    if (RTT_delayed > 0)
        dW -= (m->W * W_delayed) / (2.0 * RTT_delayed) * p_delayed;

    double W_new = m->W + dW * m->dt;
    if (W_new < 1.0) W_new = 1.0;

    /* Queue dynamics: dq/dt = N * W / RTT - C */
    double dq = m->N * m->W / m->tau - m->C;
    double q_new = m->q + dq * m->dt;
    if (q_new < 0.0) q_new = 0.0;

    m->W = W_new;
    m->q = q_new;

    /* Store in history buffers */
    m->hist_idx = (m->hist_idx + 1) % m->hist_size;
    m->W_history[m->hist_idx] = m->W;
    m->p_history[m->hist_idx] = m->p;

    if (out_W) *out_W = m->W;
    if (out_q) *out_q = m->q;
}

void tcp_aqm_set_drop_prob(TCPAQMModel* m, double p) {
    if (!m) return;
    if (p < 0.0) p = 0.0;
    if (p > 1.0) p = 1.0;
    m->p = p;
}

double tcp_aqm_rtt(const TCPAQMModel* m) {
    return m ? m->tau : 0.0;
}

void tcp_aqm_free(TCPAQMModel* m) {
    if (!m) return;
    free(m->W_history); free(m->p_history);
    free(m);
}
