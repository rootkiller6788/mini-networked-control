#include "consensus_types.h"
#include "consensus_graph.h"
#include "consensus_dynamics.h"
#include "consensus_algorithms.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ============================================================================
 * Consensus Applications (L7)
 *
 * Real-world applications of distributed consensus:
 *   App 1: Multi-robot rendezvous / formation control
 *   App 2: Distributed clock synchronization (sensor networks)
 *   App 3: Multi-UAV flocking and alignment
 *   App 4: Smart grid load balancing
 *   App 5: Autonomous vehicle platooning
 * ============================================================================ */

/* ============================================================================
 * Application 1: Multi-Robot Rendezvous (L7)
 *
 * Problem: N robots at initial positions must meet at a common point
 * without global coordination. Each robot only knows its neighbors.
 *
 * Solution: Consensus on position — each robot moves towards the average
 * of its neighbors' positions. This is exactly the continuous consensus
 * protocol with single-integrator dynamics.
 *
 * Reference: Lin, Morse, Anderson (2007) —
 *   "The multi-agent rendezvous problem", IEEE CDC
 * ============================================================================ */

typedef struct {
    double* positions;    /* N x dim: robot positions */
    double* velocities;   /* N x dim: robot velocities */
    int n_robots;
    int dim;              /* 2D or 3D */
    ConsensusGraph* comm_graph;
    double* formation_offsets; /* desired relative positions for formation */
    bool rendezvous_achieved;
    double rendezvous_time;
} RobotSwarm;

RobotSwarm* robot_swarm_create(int n_robots, int dim) {
    RobotSwarm* rs = (RobotSwarm*)calloc(1, sizeof(RobotSwarm));
    if (!rs) return NULL;
    rs->n_robots = n_robots;
    rs->dim = dim;
    rs->positions = (double*)calloc((size_t)n_robots * (size_t)dim, sizeof(double));
    rs->velocities = (double*)calloc((size_t)n_robots * (size_t)dim, sizeof(double));
    rs->formation_offsets = (double*)calloc((size_t)n_robots * (size_t)dim, sizeof(double));
    rs->comm_graph = NULL;
    rs->rendezvous_achieved = false;
    return rs;
}

void robot_swarm_free(RobotSwarm* rs) {
    if (!rs) return;
    free(rs->positions);
    free(rs->velocities);
    free(rs->formation_offsets);
    if (rs->comm_graph) consensus_graph_free(rs->comm_graph);
    free(rs);
}

void robot_swarm_set_random_positions(RobotSwarm* rs, double range) {
    if (!rs) return;
    for (int i = 0; i < rs->n_robots; i++)
        for (int k = 0; k < rs->dim; k++)
            rs->positions[i * rs->dim + k] = consensus_rand_uniform(-range, range);
}

int robot_swarm_rendezvous_step(RobotSwarm* rs, double dt) {
    /* One step of rendezvous: each robot moves towards its neighbors' average.
     * Control law: u_i = -sum_{j in N_i} (x_i - x_j)
     * This is identical to continuous-time consensus on position.
     * The rendezvous point is the average of initial positions
     * if the communication graph is connected.
     * Complexity: O(N * degree * dim). */
    if (!rs || !rs->comm_graph) return -1;
    int N = rs->n_robots, d = rs->dim;

    for (int i = 0; i < N; i++) {
        for (int k = 0; k < d; k++) rs->velocities[i * d + k] = 0.0;
        for (int j = 0; j < N; j++) {
            double a_ij = rs->comm_graph->adjacency[i * N + j];
            if (a_ij > 0) {
                for (int k = 0; k < d; k++) {
                    double rel_pos = rs->positions[i * d + k] -
                                     (rs->positions[j * d + k] +
                                      rs->formation_offsets[j * d + k] -
                                      rs->formation_offsets[i * d + k]);
                    rs->velocities[i * d + k] -= a_ij * rel_pos;
                }
            }
        }
        /* Update positions */
        for (int k = 0; k < d; k++)
            rs->positions[i * d + k] += rs->velocities[i * d + k] * dt;
    }

    /* Check rendezvous: max distance between any pair < tolerance */
    double max_dist = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double dist2 = 0.0;
            for (int k = 0; k < d; k++) {
                double diff = rs->positions[i * d + k] - rs->positions[j * d + k];
                dist2 += diff * diff;
            }
            if (dist2 > max_dist) max_dist = dist2;
        }
    }
    rs->rendezvous_achieved = (sqrt(max_dist) < 0.01);
    return 0;
}

/* ============================================================================
 * Application 2: Distributed Clock Synchronization (L7)
 *
 * Problem: N sensor nodes each with local clock exhibiting drift.
 * Goal: synchronize all clocks to a common virtual time.
 *
 * Solution: Each node runs a virtual clock that is adjusted via consensus
 * on the clock offset and drift rate.
 *
 * Model: tau_i(t) = alpha_i * t + beta_i (local clock)
 * Virtual clock: T_i(t) = a_i * tau_i(t) + b_i
 *
 * Consensus on logical clock rates (alpha_i) and offsets (beta_i) ensures
 * all nodes agree on a common time.
 *
 * Reference: Schenato & Fiorentin (2011), "Average TimeSynch: A consensus-based
 *   protocol for clock synchronization in wireless sensor networks"
 * ============================================================================ */

typedef struct {
    double skew;          /* clock skew (frequency offset) */
    double offset;        /* clock offset (phase error) */
    double virtual_skew;  /* compensated skew */
    double virtual_offset;/* compensated offset */
    double logical_time;  /* synchronized logical clock */
    bool synchronized;
} ClockNode;

typedef struct {
    ClockNode* nodes;
    int n_nodes;
    ConsensusGraph* network;
    double sync_error;    /* max logical time difference */
} ClockSyncNetwork;

ClockSyncNetwork* clocksync_create(int n_nodes) {
    ClockSyncNetwork* cs = (ClockSyncNetwork*)calloc(1, sizeof(ClockSyncNetwork));
    if (!cs) return NULL;
    cs->n_nodes = n_nodes;
    cs->nodes = (ClockNode*)calloc((size_t)n_nodes, sizeof(ClockNode));
    cs->network = consensus_graph_create(n_nodes, TOPO_UNDIRECTED);
    /* Initialize clocks with random skews and offsets */
    for (int i = 0; i < n_nodes; i++) {
        cs->nodes[i].skew = 1.0 + consensus_rand_uniform(-0.01, 0.01); /* ~1% drift */
        cs->nodes[i].offset = consensus_rand_uniform(-1.0, 1.0);
        cs->nodes[i].virtual_skew = 1.0;
        cs->nodes[i].virtual_offset = 0.0;
        cs->nodes[i].logical_time = 0.0;
    }
    return cs;
}

void clocksync_free(ClockSyncNetwork* cs) {
    if (!cs) return;
    free(cs->nodes);
    if (cs->network) consensus_graph_free(cs->network);
    free(cs);
}

int clocksync_consensus_step(ClockSyncNetwork* cs, double real_time_step) {
    /* One synchronization round:
     * 1. Each node updates logical time: T_i = a_i * (alpha_i * t + beta_i) + b_i
     * 2. Exchange logical times with neighbors
     * 3. Adjust virtual clock parameters via consensus on skew and offset
     *
     * Consensus on skew: a_i[k+1] = sum_j w_ij a_j[k]
     * Consensus on offset: b_i[k+1] = sum_j w_ij b_j[k]
     *
     * Complexity: O(N^2). */
    if (!cs) return -1;
    int N = cs->n_nodes;

    /* Update logical times */
    for (int i = 0; i < N; i++) {
        cs->nodes[i].logical_time =
            cs->nodes[i].virtual_skew * (cs->nodes[i].skew * real_time_step +
                                          cs->nodes[i].offset) +
            cs->nodes[i].virtual_offset;
    }

    /* Consensus on virtual skew (rate compensation) */
    double* new_skew = (double*)malloc((size_t)N * sizeof(double));
    double* new_offset = (double*)malloc((size_t)N * sizeof(double));

    for (int i = 0; i < N; i++) {
        new_skew[i] = 0.0;
        new_offset[i] = 0.0;
        double total_weight = 0.0;
        for (int j = 0; j < N; j++) {
            double w = cs->network->perron_matrix[i * N + j];
            new_skew[i] += w * cs->nodes[j].virtual_skew;
            new_offset[i] += w * cs->nodes[j].virtual_offset;
            total_weight += w;
        }
        if (total_weight > 0) {
            new_skew[i] /= total_weight;
            new_offset[i] /= total_weight;
        }
    }

    for (int i = 0; i < N; i++) {
        cs->nodes[i].virtual_skew = new_skew[i];
        cs->nodes[i].virtual_offset = new_offset[i];
    }
    free(new_skew); free(new_offset);

    /* Compute max synchronization error */
    double t_min = cs->nodes[0].logical_time, t_max = t_min;
    for (int i = 1; i < N; i++) {
        if (cs->nodes[i].logical_time < t_min) t_min = cs->nodes[i].logical_time;
        if (cs->nodes[i].logical_time > t_max) t_max = cs->nodes[i].logical_time;
    }
    cs->sync_error = t_max - t_min;
    return 0;
}

/* ============================================================================
 * Application 3: UAV Flocking / Alignment (L7)
 *
 * Reynolds' three rules of flocking (1987):
 *   1. Separation: avoid crowding neighbors (short-range repulsion)
 *   2. Alignment: steer towards average heading of neighbors
 *   3. Cohesion: steer towards average position of neighbors
 *
 * Consensus provides alignment and cohesion. Separation is an additional
 * collision avoidance term.
 *
 * Reference: Olfati-Saber (2006), "Flocking for multi-agent dynamic systems:
 *   algorithms and theory", IEEE TAC
 * ============================================================================ */

typedef struct {
    double x, y;       /* 2D position */
    double vx, vy;     /* 2D velocity */
    double heading;    /* orientation */
} UAVAgent;

typedef struct {
    UAVAgent* agents;
    int n_agents;
    ConsensusGraph* comm_graph;
    double separation_range;   /* minimum distance */
    double alignment_weight;
    double cohesion_weight;
    double separation_weight;
    double max_speed;
} UAVSwarm;

UAVSwarm* uav_swarm_create(int n_agents, double comm_range) {
    (void)comm_range; /* reserved for future proximity-based graph radius */
    UAVSwarm* swarm = (UAVSwarm*)calloc(1, sizeof(UAVSwarm));
    if (!swarm) return NULL;
    swarm->n_agents = n_agents;
    swarm->agents = (UAVAgent*)calloc((size_t)n_agents, sizeof(UAVAgent));
    swarm->comm_graph = consensus_graph_create(n_agents, TOPO_UNDIRECTED);
    swarm->separation_range = 2.0;
    swarm->alignment_weight = 1.0;
    swarm->cohesion_weight = 0.5;
    swarm->separation_weight = 2.0;
    swarm->max_speed = 5.0;

    /* Initialize random positions and velocities */
    for (int i = 0; i < n_agents; i++) {
        swarm->agents[i].x = consensus_rand_uniform(-20, 20);
        swarm->agents[i].y = consensus_rand_uniform(-20, 20);
        swarm->agents[i].vx = consensus_rand_uniform(-1, 1);
        swarm->agents[i].vy = consensus_rand_uniform(-1, 1);
        swarm->agents[i].heading = atan2(swarm->agents[i].vy,
                                          swarm->agents[i].vx);
    }
    return swarm;
}

void uav_swarm_free(UAVSwarm* swarm) {
    if (!swarm) return;
    free(swarm->agents);
    if (swarm->comm_graph) consensus_graph_free(swarm->comm_graph);
    free(swarm);
}

int uav_swarm_flocking_step(UAVSwarm* swarm, double dt) {
    /* One flocking step incorporating all three Reynolds rules via consensus.
     *
     * Alignment: u_i^align = sum_{j in N_i} (v_j - v_i)
     * Cohesion:  u_i^cohes = sum_{j in N_i} (x_j - x_i)
     * Separation: u_i^sep = sum_{j: ||x_j-x_i|| < d_sep} (x_i - x_j) / ||x_i-x_j||^2
     *
     * Combined: u_i = w_align * u_i^align + w_cohes * u_i^cohes + w_sep * u_i^sep
     *
     * Complexity: O(N^2). */
    if (!swarm) return -1;
    int N = swarm->n_agents;

    /* Update communication graph based on proximity */
    double* positions = (double*)calloc((size_t)N * 2, sizeof(double));
    for (int i = 0; i < N; i++) {
        positions[i * 2] = swarm->agents[i].x;
        positions[i * 2 + 1] = swarm->agents[i].y;
    }
    consensus_graph_set_proximity(swarm->comm_graph, positions, 2, 10.0);
    free(positions);

    /* Compute control forces */
    double* fx = (double*)calloc((size_t)N, sizeof(double));
    double* fy = (double*)calloc((size_t)N, sizeof(double));

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) continue;
            double dx = swarm->agents[j].x - swarm->agents[i].x;
            double dy = swarm->agents[j].y - swarm->agents[i].y;
            double dist = sqrt(dx * dx + dy * dy);
            if (dist < 1e-6) continue;

            /* Alignment (consensus on velocity) */
            fx[i] += swarm->alignment_weight * (swarm->agents[j].vx - swarm->agents[i].vx);
            fy[i] += swarm->alignment_weight * (swarm->agents[j].vy - swarm->agents[i].vy);

            /* Cohesion (consensus on position) */
            fx[i] += swarm->cohesion_weight * dx;
            fy[i] += swarm->cohesion_weight * dy;

            /* Separation (collision avoidance) */
            if (dist < swarm->separation_range) {
                double force = swarm->separation_weight / (dist * dist + 0.01);
                fx[i] -= force * dx / dist;
                fy[i] -= force * dy / dist;
            }
        }
    }

    /* Update velocities and positions */
    for (int i = 0; i < N; i++) {
        swarm->agents[i].vx += fx[i] * dt;
        swarm->agents[i].vy += fy[i] * dt;

        /* Speed limit */
        double speed = sqrt(swarm->agents[i].vx * swarm->agents[i].vx +
                           swarm->agents[i].vy * swarm->agents[i].vy);
        if (speed > swarm->max_speed) {
            swarm->agents[i].vx *= swarm->max_speed / speed;
            swarm->agents[i].vy *= swarm->max_speed / speed;
        }

        swarm->agents[i].x += swarm->agents[i].vx * dt;
        swarm->agents[i].y += swarm->agents[i].vy * dt;
        swarm->agents[i].heading = atan2(swarm->agents[i].vy,
                                          swarm->agents[i].vx);
    }
    free(fx); free(fy);
    return 0;
}

/* ============================================================================
 * Application 4: Autonomous Vehicle Platooning (L7)
 *
 * Problem: A platoon of N vehicles must maintain a desired inter-vehicle
 * spacing while following a lead vehicle.
 *
 * Solution: Consensus on spacing error. Each vehicle adjusts speed based on
 * relative distance to predecessor and possibly to other vehicles.
 *
 * Consensus approach: x_i_desired = x_leader - i * d_desired
 * Each vehicle i adjusts speed to minimize ||x_i - x_i_desired||
 * and ||v_i - v_leader|| via neighbor communication.
 *
 * Reference: Stankovic, Stanojevic, Siljak (2000) —
 *   "Decentralized overlapping control of a platoon of vehicles", IEEE TAC
 * ============================================================================ */

typedef struct {
    double position;     /* longitudinal position */
    double velocity;     /* longitudinal speed */
    double desired_spacing; /* desired gap to predecessor */
    double spacing_error;
    int leader_id;
} VehicleAgent;

typedef struct {
    VehicleAgent* vehicles;
    int n_vehicles;
    ConsensusGraph* comm_graph;
    double* spacing_errors;
} VehiclePlatoon;

VehiclePlatoon* platoon_create(int n_vehicles) {
    VehiclePlatoon* vp = (VehiclePlatoon*)calloc(1, sizeof(VehiclePlatoon));
    if (!vp) return NULL;
    vp->n_vehicles = n_vehicles;
    vp->vehicles = (VehicleAgent*)calloc((size_t)n_vehicles, sizeof(VehicleAgent));
    vp->spacing_errors = (double*)calloc((size_t)n_vehicles, sizeof(double));
    /* Initialize as a line with random positions */
    for (int i = 0; i < n_vehicles; i++) {
        vp->vehicles[i].position = -i * 10.0 + consensus_rand_uniform(-2, 2);
        vp->vehicles[i].velocity = consensus_rand_uniform(15, 25); /* m/s */
        vp->vehicles[i].desired_spacing = 10.0;  /* 10m gap */
        vp->vehicles[i].leader_id = 0;
    }
    vp->vehicles[0].velocity = 20.0; /* Leader speed */
    return vp;
}

void platoon_free(VehiclePlatoon* vp) {
    if (!vp) return;
    free(vp->vehicles);
    free(vp->spacing_errors);
    if (vp->comm_graph) consensus_graph_free(vp->comm_graph);
    free(vp);
}

int platoon_consensus_step(VehiclePlatoon* vp, double dt) {
    /* One platoon control step:
     * 1. Compute spacing error for each vehicle
     * 2. Run consensus on spacing errors and velocities
     * 3. Adjust speeds
     *
     * Control: a_i = -k_p * (x_i - x_{i-1} + d_desired) - k_v * consensus_vel_error
     * Complexity: O(N). */
    if (!vp || vp->n_vehicles < 2) return -1;
    int N = vp->n_vehicles;

    /* Compute spacing errors */
    for (int i = 1; i < N; i++) {
        double actual_gap = vp->vehicles[i - 1].position - vp->vehicles[i].position;
        vp->vehicles[i].spacing_error = actual_gap - vp->vehicles[i].desired_spacing;
    }

    /* Leader maintains constant speed; followers adjust */
    for (int i = 1; i < N; i++) {
        /* Simple proportional control on spacing */
        double kp = 0.5, kv = 1.0;
        double accel = kp * vp->vehicles[i].spacing_error +
                       kv * (vp->vehicles[i - 1].velocity - vp->vehicles[i].velocity);
        /* Limit acceleration */
        if (accel > 3.0) accel = 3.0;
        if (accel < -5.0) accel = -5.0;

        vp->vehicles[i].velocity += accel * dt;
        if (vp->vehicles[i].velocity < 0) vp->vehicles[i].velocity = 0;

        vp->vehicles[i].position += vp->vehicles[i].velocity * dt;
    }
    /* Leader maintains speed */
    vp->vehicles[0].position += vp->vehicles[0].velocity * dt;
    return 0;
}