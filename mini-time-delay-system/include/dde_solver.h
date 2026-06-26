#ifndef DDE_SOLVER_H
#define DDE_SOLVER_H

#include "time_delay_system.h"

/* ============================================================================
 * Delay Differential Equation (DDE) Numerical Integration
 *
 * Reference:
 *   L. F. Shampine & S. Thompson, "Solving DDEs in MATLAB" (2001)
 *     — dde23 algorithm description
 *   A. Bellen & M. Zennaro, "Numerical Methods for Delay Differential
 *     Equations" (2003) — comprehensive DDE numerics
 *   C. T. H. Baker, C. A. H. Paul, D. R. Willé, "Issues in the numerical
 *     solution of evolutionary delay differential equations" (1995)
 *
 * Key challenge: DDE solvers must handle discontinuities that propagate
 * along characteristic lines t → t+τ → t+2τ → ...
 *
 * Level 5 — Numerical Algorithms
 * ============================================================================ */

/* ============================================================================
 * DDE Solver Methods
 * ============================================================================ */

typedef enum {
    DDE_METHOD_STEPS = 0,       /* Method of steps (exact for constant delay) */
    DDE_METHOD_EULER = 1,       /* Forward Euler (1st order, fast) */
    DDE_METHOD_RK2 = 2,         /* Runge-Kutta 2 / Heun's method */
    DDE_METHOD_RK4 = 3,         /* Classic Runge-Kutta 4 (4th order) */
    DDE_METHOD_RKF45 = 4,       /* Runge-Kutta-Fehlberg 4(5) adaptive */
    DDE_METHOD_DDE23 = 5        /* MATLAB-style continuous Runge-Kutta */
} DDEMethod;

/* History interpolation method */
typedef enum {
    DDE_INTERP_LINEAR = 0,      /* Linear interpolation between stored points */
    DDE_INTERP_HERMITE = 1,     /* Hermite cubic (uses derivative info) */
    DDE_INTERP_SPLINE = 2       /* Cubic spline interpolation */
} DDEInterpMethod;

/* ============================================================================
 * DDE Solver Configuration
 * ============================================================================ */

typedef struct {
    DDEMethod method;
    DDEInterpMethod interp;
    double dt;                    /* Nominal time step */
    double dt_min;                /* Minimum step (for adaptive methods) */
    double dt_max;                /* Maximum step (for adaptive methods) */
    double t_start;               /* Start time */
    double t_end;                 /* End time */
    double rel_tol;               /* Relative tolerance (for adaptive) */
    double abs_tol;               /* Absolute tolerance (for adaptive) */
    int max_steps;                /* Maximum number of steps */
    bool adaptive;                /* Use adaptive step sizing */
    bool track_discontinuities;   /* Track and handle derivative breaks */
} DDESolverConfig;

/* ============================================================================
 * Solution Output
 * ============================================================================ */

typedef struct {
    int n_steps;                  /* Number of steps taken */
    int n_states;                 /* State dimension */
    double* t;                    /* Time vector (n_steps×1) */
    double* x;                    /* State history (n_steps×n_states) */
    double* xdot;                 /* Derivative history */
    bool success;                 /* Integration succeeded */
    double final_time;            /* Actual final time reached */
    const char* error_msg;        /* Error message if failed */
} DDESolution;

/* ============================================================================
 * Solver API
 * ============================================================================ */

/* Create default solver configuration */
DDESolverConfig dde_config_default(void);

/* Create solver config for method of steps */
DDESolverConfig dde_config_steps(double dt, double t_end);

/* Create solver config for fixed-step RK4 */
DDESolverConfig dde_config_rk4(double dt, double t_end);

/* Create solver config for adaptive RKF45 */
DDESolverConfig dde_config_rkf45(double dt_min, double dt_max,
                                  double t_end, double rel_tol, double abs_tol);

/* --- Main solve --- */

/* Solve a time-delay system from t0 to t_end.
 * Requires:
 *   - sys: fully configured TimeDelaySystem with DDE RHS or linear model
 *   - config: solver configuration
 *   - history_func: φ(t) for t ∈ [-τ_max, 0]
 *   - x0: initial state at t=0 (used if history_func is NULL) */
DDESolution* dde_solve(const TimeDelaySystem* sys,
                       const DDESolverConfig* config,
                       HistoryFunc history_func,
                       const double* x0);

/* Free solution memory */
void dde_solution_free(DDESolution* sol);

/* Print solution summary */
void dde_solution_print(const DDESolution* sol);

/* Extract state at a given time index */
const double* dde_solution_state_at(const DDESolution* sol, int idx);

/* Extract time at a given index */
double dde_solution_time_at(const DDESolution* sol, int idx);

/* ============================================================================
 * Method of Steps
 * ============================================================================ */

/* The method of steps solves DDEs by integrating over intervals of
 * length τ. On each interval [kτ, (k+1)τ], the delayed term x(t-τ)
 * is KNOWN from the previous interval, reducing the DDE to an ODE.
 *
 * For ẋ(t) = f(t, x(t), x(t-τ)):
 *   On [0, τ]:   ẋ(t) = f(t, x(t), φ(t-τ))  — ODE with known forcing
 *   On [τ, 2τ]:  ẋ(t) = f(t, x(t), x_prev(t-τ)) — etc. */

/* Solve one interval [t0, t0+dt_total] given that x(t-τ) is known.
 * Used internally by the method of steps. */
void dde_solve_one_interval(const TimeDelaySystem* sys,
                            const double* delayed_values,
                            int n_delayed, double dt, double dt_total,
                            double* x_start, double* x_end);

/* ============================================================================
 * History Management
 * ============================================================================ */

/* Interpolate the history buffer to get x(t) for t ≤ t_current.
 * For t > t_current, returns 0 and sets x to NaN. */
bool dde_interpolate_history(const DDESolution* sol, double t,
                              double* x_out, int n);

/* Standard constant history: φ(t) = x0 for all t ∈ [-τ_max, 0] */
void history_constant(double t, int n, double* x);

/* Zero history: φ(t) = 0 */
void history_zero(double t, int n, double* x);

/* Step history: φ(t) = x_a for t < -τ/2, φ(t) = x_b for t ≥ -τ/2 */
typedef struct { double* x_a; double* x_b; int n; double split_t; } StepHistoryData;
void history_step(double t, int n, double* x);

/* ============================================================================
 * Discontinuity Tracking
 * ============================================================================ */

/* Derivative discontinuities occur at t = kτ (from initial function φ)
 * and propagate: t = τ, 2τ, 3τ, ...
 * The solver must restart at each discontinuity for accurate results. */

/* Find the next discontinuity time after t_current given delay τ */
double dde_next_discontinuity(double t_current, double delay);

/* List all discontinuity times up to t_end */
double* dde_discontinuity_times(double delay, double t_end, int* out_count);

/* Check if t is at a discontinuity point (within tolerance) */
bool dde_is_discontinuity(double t, double delay, double tol);

#endif /* DDE_SOLVER_H */
