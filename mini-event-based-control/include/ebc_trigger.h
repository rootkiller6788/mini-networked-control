#ifndef EBC_TRIGGER_H
#define EBC_TRIGGER_H
#include "ebc_core.h"

/*
 * ebc_trigger.h -- Event-Triggering Conditions & Mechanisms (L2, L3, L5)
 *
 * Provides various event-triggering condition functions that
 * determine when a control update should occur. These are the
 * heart of event-based control.
 *
 * The general form (Tabuada 2007):
 *   Gamma(t, x(t), e(t)) = |e(t)| - (sigma * |x(t)| + epsilon)
 *   Event at t_k when Gamma(t_k) >= 0
 *
 * Trigger types covered:
 *   1. Absolute threshold:  |e| > epsilon
 *   2. Relative threshold:  |e| > sigma * |x|
 *   3. Mixed threshold:     |e| > sigma * |x| + epsilon
 *   4. Lyapunov-based:      |e| > sigma * sqrt(V(x)/lambda_min(P))
 *   5. Send-on-delta:       |x_k - x_prev| > delta
 *   6. Dynamic threshold:   threshold evolves with dynamics
 *   7. Integral-based:      integrate error, trigger on accumulated
 *
 * References:
 *   Miskowicz (2006) -- Send-on-delta concept
 *   Tabuada (2007) -- ISS-based event condition
 *   Girard (2015) -- Dynamic triggering mechanisms
 *   Dolk et al. (2017) -- Output-based event-triggered control
 */

/* ---- Trigger initialisation and parameter setting ---- */

/** Create default trigger parameters (mixed threshold) */
EBC_TriggerParams ebc_trigger_default(void);

/** Create relative-only trigger parameters */
EBC_TriggerParams ebc_trigger_make_relative(double sigma);

/** Create mixed threshold trigger parameters */
EBC_TriggerParams ebc_trigger_mixed(double sigma, double epsilon);

/** Create Lyapunov-based trigger parameters */
EBC_TriggerParams ebc_trigger_lyapunov(double sigma, const double* P, int n);

/** Create send-on-delta trigger parameters */
EBC_TriggerParams ebc_trigger_make_send_on_delta(double delta, int n);

/** Free internal allocations in trigger params */
void ebc_trigger_free(EBC_TriggerParams* tp);

/* ---- Trigger condition evaluation functions (L5: Algorithms) ---- */

/**
 * Evaluate the trigger condition Gamma(t, x, e) for given type.
 * Returns true if event should fire (Gamma >= 0).
 */
bool ebc_trigger_evaluate(const EBC_System* sys, const EBC_TriggerParams* tp);

/**
 * Evaluate absolute threshold: |e| > epsilon
 * Simplest form -- essentially a deadband on measurement.
 * Complexity: O(n) where n is state dimension.
 */
bool ebc_trigger_absolute(const EBC_System* sys, const EBC_TriggerParams* tp);

/**
 * Evaluate relative threshold: |e| > sigma * |x|
 * Common in networked control -- scales threshold with state magnitude.
 * Complexity: O(n).
 */
bool ebc_trigger_relative(const EBC_System* sys, const EBC_TriggerParams* tp);

/**
 * Evaluate mixed threshold: |e| > sigma * |x| + epsilon
 * Tabuada (2007) standard form.
 * The epsilon term ensures positive minimum inter-event time.
 * Complexity: O(n).
 */
bool ebc_trigger_mixed_threshold(const EBC_System* sys,
                                  const EBC_TriggerParams* tp);

/**
 * Evaluate Lyapunov-based threshold.
 * Uses the quadratic Lyapunov function V(x) = x'Px.
 * Event fires when |e| > sigma * sqrt(V(x) / lambda_min(P)).
 * Complexity: O(n^2) due to quadratic form computation.
 */
bool ebc_trigger_lyapunov_based(const EBC_System* sys,
                                 const EBC_TriggerParams* tp);

/**
 * Evaluate send-on-delta condition.
 * Event fires when norm of change since last send exceeds delta.
 * Complexity: O(n).
 */
bool ebc_trigger_send_on_delta(const EBC_System* sys,
                                const EBC_TriggerParams* tp);

/**
 * Evaluate dual-threshold hysteresis trigger.
 * Upper threshold triggers event; lower threshold resets.
 * Prevents Zeno behavior via hysteresis band.
 * Complexity: O(n).
 */
bool ebc_trigger_hysteresis(const EBC_System* sys,
                             double sigma_upper, double sigma_lower,
                             double epsilon);

/* ---- Dynamic trigger (L8: Advanced) ---- */

/**
 * Dynamic event-triggering (Girard 2015).
 * Uses an internal dynamic variable eta(t) that evolves:
 *   d(eta)/dt = -beta * eta + (sigma * |x| - |e|)
 *   eta(0) = eta_0 > 0
 * Event fires when eta(t) + theta * (sigma*|x| - |e|) <= 0.
 *
 * This enlarges inter-event times compared to static triggering.
 */
typedef struct {
    double eta;         /* current value of dynamic variable */
    double eta0;        /* initial value */
    double beta;        /* decay rate of eta */
    double theta;       /* weighting parameter, theta > 0 */
    double sigma;       /* relative threshold (same as static) */
} EBC_DynamicTrigger;

/** Create a dynamic trigger state */
EBC_DynamicTrigger ebc_dynamic_trigger_create(double eta0, double beta,
                                               double theta, double sigma);

/** Update dynamic variable eta via Euler integration */
void ebc_dynamic_trigger_update(EBC_DynamicTrigger* dt,
                                 const double* x, const double* e,
                                 int n, double h);

/** Evaluate dynamic trigger condition */
bool ebc_dynamic_trigger_evaluate(const EBC_DynamicTrigger* dt,
                                   const double* x, const double* e, int n);

/* ---- Trigger statistics and analysis ---- */

/**
 * Compute the evolution of the trigger margin over a trajectory.
 * Margin(t) = sigma*|x(t)| + epsilon - |e(t)|
 * Positive margin means no event yet; negative means event fires.
 */
void ebc_trigger_margin_trace(const double* x_traj,
                               const double* e_traj,
                               int n, int len, double sigma, double epsilon,
                               double* margins);

/**
 * Find the optimal sigma for a given inter-event time target.
 * Binary search over sigma in (0, 1) to achieve target average IET.
 */
double ebc_trigger_optimize_sigma(
    void (*system_dynamics)(double, const double*, const double*, int,
                             double*, void*),
    void* ctx, int n, int m,
    const double* K, const double* x0,
    double T, double dt,
    double target_iet, double epsilon);

#endif /* EBC_TRIGGER_H */
