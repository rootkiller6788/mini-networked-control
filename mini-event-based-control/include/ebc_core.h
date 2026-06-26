#ifndef EBC_CORE_H
#define EBC_CORE_H
#include <stdbool.h>
#include <stddef.h>

/*
 * ebc_core.h -- Event-Based Control: Core Definitions & Types
 *
 * Event-based (event-triggered) control is a paradigm where
 * control updates are triggered by events (typically state
 * deviation crossing a threshold) rather than a fixed clock.
 * This reduces communication/computation while maintaining
 * desired closed-loop properties.
 *
 * Pioneering work:
 *   Astrom & Bernhardsson (1999) -- periodic vs event-based sampling
 *   Arzen (1999) -- Simple event-based PID controller
 *   Tabuada (2007) -- ISS-based event-triggering for nonlinear sys
 *     "To sample or not to sample", IEEE TAC 52(9): 1682-1691
 *   Heemels, Johansson & Tabuada (2012) -- Tutorial survey
 *     "Intro to event-triggered and self-triggered control",
 *      IEEE TAC 57(3): 609-626
 *
 * Key concepts:
 *   ETC  = Event-Triggered Control (continuous measurement,
 *          discrete updates)
 *   STC  = Self-Triggered Control (next update time computed
 *          based on current state + model)
 *   PETC = Periodic Event-Triggered Control (measure at fixed
 *          period, update only if event detected)
 *
 * References (L4: ISS-Lyapunov for ETC):
 *   Tabuada (2007) Theorem III.1: For a nonlinear system with
 *     ISS-Lyapunov function V satisfying:
 *       alpha1(|x|) <= V(x) <= alpha2(|x|)
 *       gradV * f(x,k(x+e)) <= -alpha(|x|) + gamma(|e|)
 *     the event condition |e| <= sigma * alpha^{-1}(gamma^{-1}(...))
 *     guarantees asymptotic stability with positive MIET.
 *
 *   Heemels et al. (2012) Theorem V.1: For linear system
 *     dx/dt = Ax + BK(x+e) with A+BK Hurwitz, exists
 *     P > 0, sigma in (0,1) such that
 *       |e| <= sigma * |x|  ==>  dV/dt < 0
 *     guaranteeing global exponential stability.
 */

/* ---- Core enums (L1: Definitions) ---- */

/** Event-triggered control paradigm */
typedef enum {
    EBC_CONTINUOUS_ETC = 0,  /* Continuous ETC: check at every instant */
    EBC_PERIODIC_ETC   = 1,  /* Periodic ETC: check at multiples of h */
    EBC_SELF_TRIGGERED = 2,  /* Self-triggered: compute next time from model */
    EBC_TIME_TRIGGERED = 3,  /* Traditional periodic (baseline) */
    EBC_MIXED_TRIGGER  = 4   /* Mixed time/event-triggered */
} EBC_Paradigm;

/** Event-triggering condition type */
typedef enum {
    EBC_ABSOLUTE_ERROR   = 0, /* |e(t)| > epsilon */
    EBC_RELATIVE_ERROR   = 1, /* |e(t)| > sigma * |x(t)| */
    EBC_MIXED_THRESHOLD  = 2, /* |e(t)| > sigma*|x(t)| + epsilon */
    EBC_LYAPUNOV_BASED   = 3, /* V(x) related triggering */
    EBC_PERFORMANCE_BASED = 4, /* Trigger on performance degradation */
    EBC_DEADBAND         = 5, /* Simple deadband on measurement */
    EBC_SENDBOX_DELTA    = 6  /* Send-on-delta: trigger when change > delta */
} EBC_TriggerType;

/** Stability status for event-based system */
typedef enum {
    EBC_STABLE            = 0,
    EBC_ASYM_STABLE       = 1,
    EBC_EXP_STABLE        = 2,
    EBC_ISS_STABLE        = 3,  /* Input-to-State Stable */
    EBC_PRACTICAL_STABLE  = 4,
    EBC_UNSTABLE          = 5,
    EBC_ZENO              = 6,  /* Zeno behavior detected */
    EBC_INCONCLUSIVE      = 7
} EBC_StabilityResult;

/** Inter-event time classification */
typedef enum {
    EBC_IET_POSITIVE    = 0, /* Positive lower bound guaranteed */
    EBC_IET_ZERO_LIMIT  = 1, /* IET approaches zero (Zeno risk) */
    EBC_IET_MIN_GUAR    = 2, /* Minimum IET guaranteed > tau_min */
    EBC_IET_AUTO        = 3  /* Automatically determined lower bound */
} EBC_IET_Class;

/* ---- Core type definitions (L3: Mathematical Structures) ---- */

/**
 * EBC_System: Event-triggered control system
 *
 * State evolves as: dx/dt = f(t, x, u)
 * With sample-and-hold: u(t) = k(x(t_k)) for t in [t_k, t_{k+1})
 * Measurement error: e(t) = x(t_k) - x(t) for t in [t_k, t_{k+1})
 * Event condition: Gamma(x(t), e(t)) >= 0 triggers next update
 */
typedef struct {
    double* x;           /* current state vector (length n) */
    double* x_last;      /* last sampled state */
    double* e;           /* measurement error e = x_last - x */
    double* u;           /* control input (length m) */
    double* dx;          /* workspace for state derivative */
    int     n;           /* state dimension */
    int     m;           /* input dimension */
    double  t;           /* current time */
    double  t_last;      /* last event time */
    double  t_next;      /* next scheduled event (STC) */
    int     event_count; /* number of events triggered */
    double  last_iet;    /* most recent inter-event time */
    double  min_iet;     /* minimum inter-event time observed */
    double  max_iet;     /* maximum inter-event time observed */
    EBC_Paradigm   paradigm;
    EBC_TriggerType trigger;
    void*   context;     /* user-defined system context */
} EBC_System;

/**
 * EBC_TriggerParams: Event-triggering condition parameters
 *
 * Standard mixed threshold (Tabuada 2007):
 *   Gamma(t,x,e) = |e| - (sigma * |x| + epsilon)
 *   Event occurs when Gamma >= 0, i.e., |e| > sigma*|x| + epsilon
 *
 * For Lyapunov-based (Wang & Lemmon 2009):
 *   Gamma(t,x,e) = |e| - sigma * sqrt(V(x)/lambda_min(P))
 *   where V(x) = x'Px
 */
typedef struct {
    EBC_TriggerType type;
    double sigma;         /* relative threshold coefficient, 0 < sigma < 1 */
    double epsilon;       /* absolute tolerance (guarantees positive MIET) */
    double eta;           /* dynamic threshold decay/growth rate */
    double* P;            /* Lyapunov matrix P (for Lyapunov-based trigger) */
    int     n;            /* dimension of P (n x n) */
    double  V_ref;        /* reference Lyapunov value */
    double  (*ext_gamma)(double t, const double* x, const double* e,
                         int n, void* ctx); /* user-defined trigger function */
    void*   ctx;          /* context for ext_gamma */
} EBC_TriggerParams;

/**
 * EBC_Controller: Controller for event-based system
 *
 * Linear:  u = K * x(t_k)
 * Nonlinear: u = kappa(x(t_k))
 */
typedef struct {
    double* K;            /* feedback gain matrix (m x n) */
    int     n;            /* state dimension */
    int     m;            /* input dimension */
    bool    is_nonlinear; /* true if nonlinear controller */
    void    (*ctrl_law)(const double* x, int n, double* u, int m,
                        void* ctx);
    void*   ctrl_ctx;
} EBC_Controller;

/**
 * EBC_Performance: Performance metrics for event-based control
 *
 * Compares event-triggered with periodic implementation.
 * Communication reduction = 1 - events / periodic_equiv
 */
typedef struct {
    double total_events;        /* total number of control updates */
    double periodic_equiv;      /* equivalent periodic updates (T/h) */
    double comm_reduction;      /* comm reduction ratio */
    double ise;                 /* Integral of Squared Error */
    double iae;                 /* Integral of Absolute Error */
    double itae;                /* Integral of Time-weighted Absolute Error */
    double isci;                /* Integral of Squared Control Input */
    double energy_cost;         /* Total energy (comm + computation) */
    double avg_iet;             /* Average inter-event time */
    double max_state_dev;       /* Maximum state deviation between events */
    double settling_time;       /* Settling time to 2% band */
    double overshoot;           /* Maximum overshoot */
} EBC_Performance;

/**
 * EBC_StabilityCert: Stability certificate for event-based system
 *
 * ISS-Lyapunov characterization (Tabuada 2007, Theorem III.1):
 *   alpha1(|x|) <= V(x) <= alpha2(|x|)
 *   dV/dt <= -alpha3 * V(x) + gamma * |e|^2
 *
 * Exponential stability guaranteed when:
 *   sigma < sqrt(alpha3 / (gamma * norm(PBK)))
 */
typedef struct {
    EBC_StabilityResult result;
    double* P;                 /* Lyapunov matrix */
    int     n;
    double  alpha1;            /* c1*|x|^2 <= V(x) */
    double  alpha2;            /* V(x) <= c2*|x|^2 */
    double  alpha3;            /* dV/dt <= -c3*|x|^2 + gamma*|e|^2 */
    double  gamma;             /* ISS gain from e to x */
    double  sigma_critical;    /* max sigma for guaranteed stability */
    double  tau_min;           /* minimum inter-event time */
    double  convergence_rate;  /* exponential convergence rate */
    bool    is_iss;
} EBC_StabilityCert;

/* ---- Core API: System lifecycle ---- */
EBC_System*      ebc_system_create(int n, int m, EBC_Paradigm p);
void             ebc_system_set_state(EBC_System* sys, const double* x0);
void             ebc_system_set_dynamics(EBC_System* sys,
                     void (*f)(double, const double*, const double*, int,
                               double*, void*), void* ctx);
void             ebc_system_reset(EBC_System* sys);
void             ebc_system_free(EBC_System* sys);

/* ---- Core API: Event detection (L2) ---- */
bool             ebc_check_event(const EBC_System* sys,
                                 const EBC_TriggerParams* tp);
double           ebc_compute_error_norm(const EBC_System* sys);
double           ebc_compute_threshold(const EBC_System* sys,
                                       const EBC_TriggerParams* tp);
int              ebc_mark_event(EBC_System* sys, double t);

/* ---- Core API: Simulation ---- */
int              ebc_step_euler(EBC_System* sys, const EBC_Controller* ctrl,
                                double dt, const EBC_TriggerParams* tp);
int              ebc_step_rk4(EBC_System* sys, const EBC_Controller* ctrl,
                              double dt, const EBC_TriggerParams* tp);
int              ebc_simulate(EBC_System* sys, const EBC_Controller* ctrl,
                              double T, double dt,
                              const EBC_TriggerParams* tp,
                              double** traj, int* traj_len,
                              double** events, int* evt_len);

/* ---- Core API: Diagnostic ---- */
EBC_IET_Class    ebc_classify_iet(const double* event_times, int n_events);
bool             ebc_detect_zeno(const double* event_times, int n_events,
                                 double tol);

/* ---- Utility API ---- */
double           ebc_vector_norm(const double* v, int n);
double           ebc_matrix_norm(const double* M, int n);
void             ebc_matrix_multiply(const double* A, const double* B,
                                     double* C, int n);
void             ebc_matrix_vec_mul(const double* A, const double* x,
                                    double* y, int n);
double           ebc_integral_square(const double* x, int n);
void             ebc_copy_vector(const double* src, double* dst, int n);

#endif /* EBC_CORE_H */
