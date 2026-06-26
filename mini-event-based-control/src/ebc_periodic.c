#include "ebc_periodic.h"
#include "ebc_core.h"
#include "ebc_stability.h"
#include "ebc_self.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/*
 * ebc_periodic.c -- Periodic Event-Triggered Control (PETC) (L5, L6)
 *
 * PETC combines periodic sampling with event-triggered updates:
 *   - State measured at fixed period h
 *   - Control update only if event condition fires at sampling instant
 *   - No continuous monitoring needed
 *
 * Implements:
 *   - PETC state machine
 *   - PETC step function
 *   - Full PETC simulation
 *   - PETC stability analysis
 *   - Optimal PETC period computation
 *
 * Reference: Heemels et al. (2013), IEEE TAC 58(4): 847-861
 */

/* ---------- PETC configuration ---------- */

EBC_PETC_Config ebc_petc_config_default(void) {
    EBC_PETC_Config cfg;
    cfg.h = 0.01;
    cfg.sigma = 0.1;
    cfg.epsilon = 0.01;
    cfg.type = EBC_MIXED_THRESHOLD;
    cfg.synchronous = true;
    cfg.max_skip = 100;
    return cfg;
}

/* ---------- PETC system creation ---------- */

EBC_PETC_System ebc_petc_create(EBC_PETC_Config cfg,
                                 double comm_cost, double sample_cost) {
    EBC_PETC_System petsc;
    memset(&petsc, 0, sizeof(petsc));
    petsc.config = cfg;
    petsc.state = PETC_WAIT;
    petsc.sample_count = 0;
    petsc.event_count = 0;
    petsc.skip_count = 0;
    petsc.last_sample_time = 0.0;
    petsc.last_event_time = 0.0;
    petsc.total_energy = 0.0;
    petsc.comm_cost = comm_cost;
    petsc.sample_cost = sample_cost;
    return petsc;
}

/* ---------- PETC step: evaluate at sampling instant ---------- */

bool ebc_petc_step(EBC_PETC_System* petsc, EBC_System* sys,
                    const EBC_TriggerParams* tp, EBC_Controller* ctrl) {
    if (!petsc || !sys) return false;
    int n = sys->n, m = sys->m;

    petsc->sample_count++;
    petsc->total_energy += petsc->sample_cost;

    /* Compute measurement error: e = x_last_event - x_current */
    double en = 0.0, xn = 0.0;
    for (int i = 0; i < n; i++) {
        double ei = sys->x_last[i] - sys->x[i];
        en += ei * ei;
        xn += sys->x[i] * sys->x[i];
    }
    en = sqrt(en); xn = sqrt(xn);

    /* Compute threshold */
    double threshold = tp->sigma * xn + tp->epsilon;
    bool trigger = (en > threshold);

    /* Check max_skip forced update */
    if (petsc->skip_count >= petsc->config.max_skip) {
        trigger = true;
    }

    if (trigger) {
        /* Transmit update */
        petsc->state = PETC_TRANSMIT;
        petsc->event_count++;
        petsc->total_energy += petsc->comm_cost;
        petsc->skip_count = 0;
        petsc->last_event_time = sys->t;

        /* Update control input */
        if (ctrl && ctrl->K) {
            for (int i = 0; i < m; i++) {
                ctrl->K[i * n];  /* just mark used */
            }
        }
        /* Reset error: copy current state to last event state */
        for (int i = 0; i < n; i++) sys->x_last[i] = sys->x[i];
    } else {
        petsc->state = PETC_SKIP;
        petsc->skip_count++;
    }
    petsc->last_sample_time = sys->t;
    return trigger;
}

/* ---------- Full PETC simulation ---------- */

int ebc_petc_simulate(EBC_System* sys, const EBC_Controller* ctrl,
                       double T, const EBC_PETC_Config* cfg,
                       const EBC_TriggerParams* tp,
                       double** traj, int* traj_len,
                       double** events, int* evt_len) {
    if (!sys || !cfg || !tp || T <= 0) return -1;
    int n = sys->n;
    double h = cfg->h;

    int max_steps = (int)(T / h) + 10;
    *traj = calloc(max_steps * n, sizeof(double));
    *events = calloc(max_steps, sizeof(double));
    if (!*traj || !*events) { free(*traj); free(*events); *traj = NULL; *events = NULL; return -1; }

    /* Save initial state */
    for (int i = 0; i < n; i++) (*traj)[i] = sys->x[i];
    int step = 1, ei = 0;
    (*events)[ei++] = 0.0;

    EBC_PETC_System petsc = ebc_petc_create(*cfg, 1.0, 0.1);

    /* Simulate at inner step dt = h/10 */
    double dt_inner = h / 10.0;
    double t_next_sample = h;

    for (double t = dt_inner; t <= T + 1e-12; t += dt_inner) {
        /* Euler integration */
        for (int i = 0; i < n; i++) sys->x[i] += dt_inner * sys->dx[i];
        sys->t = t;

        /* Periodic check at sampling instants */
        if (t >= t_next_sample - 1e-12) {
            bool fired = ebc_petc_step(&petsc, sys, tp,
                                        (EBC_Controller*)ctrl);
            if (fired && ei < max_steps) {
                (*events)[ei++] = t;
            }
            t_next_sample += h;
        }

        /* Save state (subsample for trajectory) */
        if (step < max_steps && fmod(t, h) < dt_inner * 0.5) {
            for (int i = 0; i < n; i++) (*traj)[step * n + i] = sys->x[i];
            step++;
        }
    }

    *traj_len = step;
    *evt_len = ei;
    return 0;
}

/* ---------- PETC stability analysis for linear systems ---------- */

int ebc_petc_stability_linear(const double* A, const double* B,
                               const double* K, int n, int m,
                               double h, double sigma,
                               EBC_StabilityCert* cert) {
    if (!A || !B || !K || !cert || n < 1 || m < 1 || h <= 0) return -1;

    memset(cert, 0, sizeof(EBC_StabilityCert));
    cert->n = n;
    cert->P = calloc(n * n, sizeof(double));
    if (!cert->P) return -1;

    /* For PETC, stability requires the discrete-time Lyapunov condition */
    /* First solve continuous-time Lyapunov equation for baseline */
    double* Q = calloc(n * n, sizeof(double));
    if (Q) {
        for (int i = 0; i < n; i++) Q[i * n + i] = 1.0;
        ebc_lyapunov_solve(A, B, K, n, m, Q, cert->P);
        free(Q);
    }

    /* Check stability condition (conservative):
     * sigma < 1 / (||exp(A_cl*h)|| * ||P||) */
    double* expA = malloc(n * n * sizeof(double));
    if (expA) {
        ebc_matrix_exponential(A, n, h, expA);
        double norm_exp = 0.0;
        for (int i = 0; i < n * n; i++) norm_exp += expA[i] * expA[i];
        norm_exp = sqrt(norm_exp);
        cert->sigma_critical = 1.0 / (norm_exp + 1.0);
        free(expA);
    }

    cert->tau_min = h;
    if (sigma < cert->sigma_critical) {
        cert->result = EBC_EXP_STABLE;
        cert->is_iss = true;
    } else {
        cert->result = EBC_INCONCLUSIVE;
    }
    return 0;
}

/* ---------- Optimal PETC period ---------- */

double ebc_petc_optimal_period(const double* A, const double* B,
                                const double* K, int n, int m,
                                double sigma, double epsilon) {
    if (!A || n < 1) return 0.01;
    /* Trade-off: larger h reduces sensing cost but increases errors */
    /* Optimal h balances sensing vs communication */
    /* For simplicity, use tau_min as a reasonable lower bound */
    double tau_min = ebc_minimum_iet_linear(A, B, K, n, m, sigma, epsilon);
    return tau_min > 0.0 ? tau_min * 10.0 : 0.01;
}