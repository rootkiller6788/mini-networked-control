#ifndef TIME_DELAY_SYSTEM_H
#define TIME_DELAY_SYSTEM_H

#include <stdbool.h>
#include <stddef.h>

/* ============================================================================
 * Time-Delay System — Core Definitions and Types
 *
 * Reference Works:
 *   J. K. Hale & S. M. Verduyn Lunel, "Introduction to Functional
 *     Differential Equations" (1993) — foundational DDE theory
 *   K. Gu, V. L. Kharitonov, J. Chen, "Stability of Time-Delay Systems" (2003)
 *     — LMI-based stability, delay margin computation
 *   S.-I. Niculescu, "Delay Effects on Stability" (2001) — comprehensive survey
 *   J.-P. Richard, "Time-delay systems: an overview of some recent advances
 *     and open problems" (Automatica, 2003)
 *   M. Krstic, "Delay Compensation for Nonlinear, Adaptive, and PDE Systems"
 *     (2009) — predictor feedback methodology
 *
 * Level 1 — Core Definitions
 * ============================================================================ */

/* --- Delay Type Taxonomy --- */

typedef enum {
    DELAY_CONSTANT = 0,       /* τ(t) = τ₀, fixed known delay */
    DELAY_TIME_VARYING = 1,   /* τ(t) varies within [τ_min, τ_max] */
    DELAY_STOCHASTIC = 2,     /* τ ~ distribution (e.g., exponential, Gaussian) */
    DELAY_DISTRIBUTED = 3,    /* ∫₀^∞ K(θ) x(t-θ) dθ, distributed kernel */
    DELAY_STATE_DEPENDENT = 4 /* τ = τ(x(t)), delay depends on system state */
} DelayType;

/* --- DDE Classification --- */

typedef enum {
    DDE_RETARDED = 0,         /* ẋ(t) = f(t, x(t), x(t-τ))  — most common */
    DDE_NEUTRAL = 1,          /* ẋ(t) = f(t, x(t), x(t-τ), ẋ(t-τ)) */
    DDE_ADVANCED = 2,         /* ẋ(t) = f(t, x(t), x(t+τ)) — future-dependent */
    DDE_INTEGRO = 3           /* ẋ(t) = ∫ f(t,s,x(s)) ds — Volterra-type */
} DDEType;

/* --- Delay System State Classification --- */

typedef enum {
    DELAY_STABLE = 0,              /* All characteristic roots in LHP */
    DELAY_MARGINALLY_STABLE = 1,   /* Roots on imaginary axis */
    DELAY_UNSTABLE = 2,            /* At least one root in RHP */
    DELAY_INDEPENDENT_STABLE = 3,  /* Stable for ALL τ ≥ 0 */
    DELAY_DEPENDENT_STABLE = 4     /* Stable for τ ∈ [0, τ*) */
} DelayStabilityClass;

/* --- Network-Induced Delay Sources --- */

typedef enum {
    NET_DELAY_SENSOR_CONTROLLER = 0,  /* τ_sc: sensor-to-controller delay */
    NET_DELAY_CONTROLLER_ACTUATOR = 1,/* τ_ca: controller-to-actuator delay */
    NET_DELAY_COMPUTATION = 2,        /* τ_c: computational delay */
    NET_DELAY_QUEUING = 3,            /* τ_q: network queuing delay */
    NET_DELAY_PROPAGATION = 4,        /* τ_p: physical propagation delay */
    NET_DELAY_SERIALIZATION = 5       /* τ_s: packet serialization delay */
} NetworkDelaySource;

/* ============================================================================
 * Core Data Structures
 * ============================================================================ */

/* --- Delay Descriptor --- */
typedef struct {
    DelayType type;           /* Classification of delay */
    double tau_nominal;       /* Nominal/average delay value (seconds) */
    double tau_min;           /* Minimum delay bound */
    double tau_max;           /* Maximum delay bound */
    double tau_variance;      /* Variance for stochastic delays */
    double* kernel;           /* Kernel K(θ) for distributed delays */
    int kernel_size;          /* Kernel discretization points */
    double derivative_bound;  /* Upper bound on |dτ/dt| < 1 typically */
    bool is_bounded;          /* True if τ(t) has known bounds */
} DelayDescriptor;

/* --- DDE Right-Hand Side Function Signature ---
 * f(t, x, y, dx_delayed, n, params) -> ẋ(t)
 *   where y = x(t-τ) and dx_delayed = ẋ(t-τ) for neutral DDEs */
typedef void (*DDERHSFunc)(double t,
                           const double* x,
                           const double* x_delayed,
                           const double* dx_delayed,
                           int n,
                           const double* params,
                           double* dxdt);

/* --- History Function: φ(t) for t ∈ [-τ_max, 0] --- */
typedef void (*HistoryFunc)(double t, int n, double* x);

/* --- Time-Delay System (primary structure) --- */
typedef struct {
    /* System identity */
    char* name;                    /* System name */
    int n_states;                  /* Dimension of state vector x ∈ Rⁿ */
    int n_outputs;                 /* Number of outputs */
    int n_inputs;                  /* Number of control inputs */

    /* Nominal linear model: ẋ = A x + A_d x(t-τ) + B u */
    double* A;                     /* n×n system matrix */
    double* A_delayed;             /* n×n delayed-state matrix */
    double* B;                     /* n×m input matrix */
    double* C;                     /* p×n output matrix */

    /* Delay specification */
    int n_delays;                  /* Number of distinct delays */
    DelayDescriptor** delays;      /* Array of delay descriptors */

    /* DDE function pointers */
    DDERHSFunc f_rhs;              /* Nonlinear RHS (NULL for linear systems) */
    HistoryFunc history;           /* History function φ(·) */
    double* rhs_params;            /* Additional parameters for f_rhs */
    int n_params;                  /* Number of RHS parameters */

    /* DDE type */
    DDEType dde_type;

    /* Characteristic roots (for analysis) */
    int n_roots;                   /* Number of computed characteristic roots */
    double* roots_real;            /* Real parts */
    double* roots_imag;            /* Imaginary parts */

    /* Stability metadata */
    DelayStabilityClass stability_class;
    double delay_margin;           /* τ* — critical delay for stability */
    double spectral_abscissa;      /* max Re(λ) of characteristic eqn */

    /* Solution storage */
    double* current_state;         /* x(t) — current state */
    double* delayed_state;         /* x(t-τ) — delayed state */
    double* history_buffer;        /* φ(θ) for θ ∈ [-τ_max, 0] */
    int history_points;            /* Number of history samples */
    double t_current;              /* Current time */

    /* Performance metrics */
    double settling_time;          /* Time to reach ±2% of steady state */
    double overshoot_percent;      /* Maximum overshoot percentage */
    double ise;                    /* Integral of Squared Error */
    double iae;                    /* Integral of Absolute Error */
} TimeDelaySystem;

/* ============================================================================
 * L2 — Core Concept Implementations
 * ============================================================================ */

/* --- Characteristic Quasipolynomial ---
 * Δ(s) = det(sI - A - A_d e^{-τs}) = 0
 * For scalar case: Δ(s) = s + a + b e^{-τs} = 0
 * The presence of e^{-τs} makes this a quasipolynomial with ∞ many roots. */

/* Compute Δ(s) for a given complex s = σ + jω */
double time_delay_characteristic_eqn(const TimeDelaySystem* sys,
                                     double sigma, double omega);

/* Compute the real and imaginary parts of Δ(s) separately */
void time_delay_char_eqn_parts(const TimeDelaySystem* sys,
                                double sigma, double omega,
                                double* out_real, double* out_imag);

/* --- State Norm for DDE ---
 * ||x_t|| = sup_{θ∈[-τ,0]} ||x(t+θ)||  (supremum norm on C([-τ,0], Rⁿ)) */
double time_delay_state_norm(const TimeDelaySystem* sys);

/* --- Delay Rate Condition ---
 * Check if |dτ/dt| ≤ d < 1 (necessary for some stability criteria) */
bool time_delay_rate_check(const DelayDescriptor* delay);

/* ============================================================================
 * System Lifecycle
 * ============================================================================ */

/* Allocate and initialize a time-delay system */
TimeDelaySystem* tds_create(const char* name, int n_states,
                             int n_inputs, int n_outputs);

/* Set the linear model matrices (deep copy) */
void tds_set_linear_model(TimeDelaySystem* sys,
                          const double* A, const double* A_d,
                          const double* B, const double* C);

/* Add a delay to the system */
void tds_add_delay(TimeDelaySystem* sys, DelayType type,
                   double tau_nominal, double tau_min, double tau_max);

/* Set the history function */
void tds_set_history(TimeDelaySystem* sys, HistoryFunc phi);

/* Set the nonlinear RHS function */
void tds_set_nonlinear_rhs(TimeDelaySystem* sys, DDERHSFunc f,
                            const double* params, int n_params);

/* Set DDE type */
void tds_set_dde_type(TimeDelaySystem* sys, DDEType dtype);

/* Free all resources */
void tds_free(TimeDelaySystem* sys);

/* Print system summary to stdout */
void tds_print(const TimeDelaySystem* sys);

/* --- Constant delay helper --- */
DelayDescriptor* delay_create_constant(double tau);
DelayDescriptor* delay_create_time_varying(double tau_min, double tau_max,
                                            double deriv_bound);
DelayDescriptor* delay_create_stochastic(double tau_mean, double tau_var);
void delay_free(DelayDescriptor* delay);
void delay_print(const DelayDescriptor* delay);

/* ============================================================================
 * Characteristic Root Computation
 * ============================================================================ */

/* Compute the rightmost characteristic roots using
 * spectral discretization of the solution operator.
 * Returns number of roots found. */
int tds_compute_characteristic_roots(TimeDelaySystem* sys, int max_roots);

/* Check if a given s = σ + jω satisfies the characteristic equation
 * within tolerance ε. */
bool tds_is_characteristic_root(const TimeDelaySystem* sys,
                                 double sigma, double omega, double eps);

/* Compute the spectral abscissa α = max{Re(λ) : Δ(λ) = 0}. */
double tds_spectral_abscissa(TimeDelaySystem* sys);

#endif /* TIME_DELAY_SYSTEM_H */
