#ifndef CPS_RESILIENCE_H
#define CPS_RESILIENCE_H

#include "cps_security_core.h"

/* ============================================================================
 * CPS Resilient Control and Secure Estimation (L4: Theorems, L5: Algorithms)
 *
 * Key theorems:
 *   Fawzi, Tabuada, Diggavi (2014) — Secure state estimation is possible
 *     iff the number of attacked sensors < p/2 (observability preservation)
 *   Shoukry & Tabuada (2016) — Event-triggered secure estimation
 *   Pajic et al. (2017) — Resilient control via receding horizon under attack
 * ============================================================================ */

/* --- Secure State Estimation (L5: Algorithms) --- */

/* ℓ₀-norm secure estimation (combinatorial search)
 * Theorem (Fawzi et al. 2014): If at most s sensors are attacked,
 * the true state can be recovered exactly if 2s < p and the system
 * is 2s-sparse observable.
 * Complexity: O(choose(p,s) * n^3) — NP-hard in general, but tractable
 * for small s with branch-and-bound */
int cps_secure_l0_estimation(double* x_est, const double* y,
                              const double* C, const double* A,
                              const double* B, const double* u,
                              const double* x_pred,
                              int n, int p, int m, int s_max);

/* ℓ₁-relaxation: min ||z||₁ s.t. z = y - C*x
 * Convex relaxation of the secure estimation problem.
 * Uses iterative soft-thresholding (ISTA)
 * Complexity: O(n^2 * p) per iteration */
int cps_secure_l1_estimation(double* x_est, const double* y,
                              const double* C,
                              const double* x0, int n, int p,
                              double lambda, int max_iter, double tol);

/* Resilient Kalman filter with attack-aware covariance inflation
 * When attack is suspected, inflate measurement covariance to
 * reduce attacker influence on state estimate */
void cps_resilient_kalman_step(double* x_est, double* P_est,
                                const double* y, const double* C,
                                const double* A, const double* B,
                                const double* u, double Q_scale,
                                double R_scale, double R_inflation,
                                int n, int p, int m);

/* Attack identification — which sensors are attacked?
 * Returns bitmask: bit i set = sensor i is attacked */
unsigned int cps_identify_attacked_sensors(const double* residuals,
                                            int p, double threshold,
                                            int min_consistent);

/* --- Resilient Control Strategies (L5: Algorithms) --- */

/* Hold-last-safe: u = u_safe (simplest strategy, no computation) */
void cps_resilient_hold(CPSResilientController* rc, double* u);

/* Fallback control: pre-computed stabilizing controller
 * u_fallback = -K_fallback * x_est */
void cps_resilient_fallback(CPSResilientController* rc,
                             const double* x_est, const double* K_fb,
                             double* u);

/* Game-theoretic resilient control: min_u max_a J(u,a)
 * Solves saddle-point: u* = argmin_u max_a (x'Qx + u'Ru - gamma*a'a)
 * Reference: Zhu & Basar (2015) — Game-Theoretic Methods for Robustness */
void cps_resilient_gametheoretic(CPSResilientController* rc,
                                  const double* x_est,
                                  const double* A, const double* B,
                                  double gamma, double* u);

/* Constraint-tightening MPC for resilience
 * u* = argmin sum_{k=0}^{N-1} (x_k'Q x_k + u_k'R u_k)
 * s.t. x_{k+1} = A x_k + B u_k
 *      x_k in X_tightened (original constraints shrunk by attack margin) */
void cps_resilient_tightened_mpc(CPSResilientController* rc,
                                  const double* x_est,
                                  const double* A, const double* B,
                                  const double* Q, const double* R,
                                  int horizon, double* u_mpc);

/* Redundant sensor switching
 * Select subset of trusted sensors for estimation */
int cps_resilient_sensor_selection(CPSResilientController* rc,
                                    const double* C, int n, int p,
                                    int min_sensors);

/* --- Reachability-Based Security Analysis (L4: Theorems) --- */

/* Compute attack-reachable set: all states an attacker can drive
 * the system to while remaining undetected.
 * Theorem: The attack-reachable set is bounded iff the system
 * under attack is detectable from the attacked measurements. */
void cps_attack_reachable_set(const double* A, const double* B,
                               const double* C, int n, int m, int p,
                               double attack_budget, int horizon,
                               double* R_min, double* R_max);

/* Invariant set under attack: states from which the system
 * cannot be driven to unsafe region despite attack */
double cps_attack_invariant_set_volume(const double* A,
                                        const double* B,
                                        const double* Gamma_a,
                                        int n, int m, int n_a);

/* --- Redundancy and Diversity Analysis --- */

/* Minimum sensor redundancy needed for attack resilience
 * Theorem: Need at least 2s+1 sensors to tolerate s attacks
 * while maintaining observability (Fawzi et al. 2014) */
int cps_min_sensors_for_resilience(int s_attack);

/* Sensor diversity score: higher diversity → harder to
 * simultaneously compromise multiple sensors */
double cps_sensor_diversity_score(const double* C, int n, int p);

/* Actuator redundancy: backup actuators that can substitute
 * for compromised ones */
void cps_actuator_redundancy_map(const double* B, int n, int m,
                                  int* redundancy_matrix);

/* --- System-Theoretic Attack Detectability (L4: Theorems) --- */

/* Check if an attack of given type is detectable
 * Theorem (Pasqualetti et al. 2013): An attack is undetectable
 * iff the attack signal lies in the unobservable subspace of
 * the augmented system (A, [C; detection_filter]). */
bool cps_is_attack_detectable(const double* A, const double* C,
                               const double* Gamma_a, int n, int p,
                               int n_a);

/* Zero-dynamics attack: exploit invariant zeros of the system
 * y = C*(sI-A)^{-1}*B_a * a. If the system has an invariant zero
 * in the unstable region, zero-dynamics attack can be undetectable. */
bool cps_has_zero_dynamics_vulnerability(const double* A,
                                          const double* B,
                                          const double* C,
                                          int n, int m, int p);

#endif /* CPS_RESILIENCE_H */
