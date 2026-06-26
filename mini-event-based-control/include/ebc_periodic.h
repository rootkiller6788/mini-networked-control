#ifndef EBC_PERIODIC_H
#define EBC_PERIODIC_H
#include "ebc_core.h"

/*
 * ebc_periodic.h -- Periodic Event-Triggered Control (PETC) (L5, L6)
 *
 * PETC (Heemels et al. 2013) combines periodic sampling with
 * event-triggered updates. The state is measured at a fixed
 * period h, but a control update only occurs if the event
 * condition is satisfied at the sampling instant.
 *
 * PETC advantages:
 *   - Implementable on digital hardware (periodic measurement)
 *   - Communication reduction compared to periodic
 *   - No continuous monitoring needed (unlike continuous ETC)
 *   - Easier scheduling analysis
 *
 * PETC event condition (at t = k*h, k = 0,1,2,...):
 *   Gamma(x(kh), e(kh)) >= 0  =>  update control
 *   where e(kh) = x(t_last) - x(kh)
 *
 * References:
 *   Heemels et al. (2013) -- "Periodic event-triggered control
 *     for linear systems", IEEE TAC 58(4): 847-861
 *   Postoyan et al. (2015) -- Periodic ETC for nonlinear systems
 *   Borgers & Heemels (2014) -- Periodic ETC stability analysis
 */

/* ---- PETC system and configuration ---- */

/**
 * PETC configuration parameters.
 */
typedef struct {
    double h;             /* sampling period (measurement interval) */
    double sigma;         /* relative threshold */
    double epsilon;       /* absolute tolerance */
    EBC_TriggerType type; /* trigger condition type */
    bool synchronous;     /* true: all nodes synchronized to same clock */
    double max_skip;      /* max consecutive non-event sampling instants */
} EBC_PETC_Config;

/** Create default PETC configuration */
EBC_PETC_Config ebc_petc_config_default(void);

/**
 * PETC state machine states.
 */
typedef enum {
    PETC_WAIT,           /* waiting for next sampling instant */
    PETC_SAMPLE,         /* just sampled: evaluate condition */
    PETC_TRANSMIT,       /* event triggered: transmit update */
    PETC_SKIP,           /* no event: skip transmission */
    PETC_MAX_SKIP_REACHED /* forced update due to max_skip */
} EBC_PETC_State;

/**
 * PETC system state.
 */
typedef struct {
    EBC_PETC_Config config;
    EBC_PETC_State  state;
    int    sample_count;    /* total samples taken */
    int    event_count;     /* events triggered */
    int    skip_count;      /* consecutive non-event samples */
    double last_sample_time;
    double last_event_time;
    double total_energy;    /* cumulative energy consumption */
    double comm_cost;       /* cost per transmission */
    double sample_cost;     /* cost per sample (sensing) */
} EBC_PETC_System;

/* ---- PETC API ---- */

/** Initialise PETC system */
EBC_PETC_System ebc_petc_create(EBC_PETC_Config cfg,
                                 double comm_cost, double sample_cost);

/**
 * PETC step function: called at each sampling instant k*h.
 *
 * 1. Sample state: x_k = x(k*h)
 * 2. Compute error: e_k = x_last_event - x_k
 * 3. Evaluate condition: Gamma(x_k, e_k) >= 0 ?
 * 4. If yes: transmit update, reset x_last_event = x_k
 * 5. If no: skip transmission, increment skip counter
 * 6. If skip_count > max_skip: force update
 *
 * @param petsc     PETC system state
 * @param sys       EBC system
 * @param tp        Trigger parameters
 * @param ctrl      Controller (used on event)
 * @return          true if update transmitted at this step
 */
bool ebc_petc_step(EBC_PETC_System* petsc, EBC_System* sys,
                    const EBC_TriggerParams* tp, EBC_Controller* ctrl);

/**
 * Run full PETC simulation.
 *
 * Simulates from t=0 to T, sampling at period h and
 * transmitting only when the event condition fires.
 */
int ebc_petc_simulate(EBC_System* sys, const EBC_Controller* ctrl,
                       double T, const EBC_PETC_Config* cfg,
                       const EBC_TriggerParams* tp,
                       double** traj, int* traj_len,
                       double** events, int* evt_len);

/**
 * Stability analysis for PETC.
 *
 * For linear systems with quadratic event condition:
 * The PETC scheme is globally exponentially stable if
 * there exists P > 0 and sigma in (0,1) such that
 * the discrete-time Lyapunov condition holds:
 *   V(x((k+1)h)) - V(x(kh)) <= -beta * V(x(kh))
 * at all sampling instants.
 *
 * @param A, B, K, n, m  System matrices
 * @param h              Sampling period
 * @param sigma          Event threshold
 * @param cert           Output stability certificate
 * @return               0 if stable, -1 if unstable or inconclusive
 */
int ebc_petc_stability_linear(const double* A, const double* B,
                               const double* K, int n, int m,
                               double h, double sigma,
                               EBC_StabilityCert* cert);

/**
 * Compute optimal PETC sampling period h_star.
 *
 * Trade-off: smaller h -> more samples but more potential events.
 *            larger h -> fewer samples but larger errors.
 *
 * Finds h that maximizes communication reduction while
 * maintaining stability (sigma < sigma_crit).
 */
double ebc_petc_optimal_period(const double* A, const double* B,
                                const double* K, int n, int m,
                                double sigma, double epsilon);

#endif /* EBC_PERIODIC_H */
