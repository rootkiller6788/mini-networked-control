/**
 * cloud_control_core.h - Cloud Control System Core Definitions
 *
 * Defines the fundamental data types, control loop architecture, and core
 * operations for cloud-based control systems.
 *
 * Domain: Networked Control Systems / Cloud Control Systems
 * References:
 *   - Xia et al., "Cloud Control Systems", IEEE/CAA JAS, 2015
 *   - He et al., "Cloud-Based Control: A Survey", Automatica, 2021
 *
 * Knowledge Coverage:
 *   L1: CloudControlSystem, CloudNode, EdgeNode, ControlTask, QoSProfile
 *   L2: Cloud-edge collaboration, hierarchical control, mode switching
 *   L3: State-space model with network delay
 */

#ifndef CLOUD_CONTROL_CORE_H
#define CLOUD_CONTROL_CORE_H

#include <stddef.h>
#include <stdint.h>
#include <time.h>

#define CCS_MAX_STATES          64
#define CCS_MAX_INPUTS          32
#define CCS_MAX_OUTPUTS         32
#define CCS_MAX_EDGE_NODES      256
#define CCS_MAX_CLOUD_INSTANCES 128
#define CCS_MAX_TASK_QUEUE      1024
#define CCS_MAX_HISTORY         4096
#define CCS_DEFAULT_PERIOD_US   10000
#define CCS_DEFAULT_MATI_US     50000
#define CCS_DEFAULT_RTT_US      5000

typedef enum {
    CC_MODE_EDGE_ONLY      = 0,
    CC_MODE_CLOUD_ONLY     = 1,
    CC_MODE_COLLABORATIVE  = 2,
    CC_MODE_HIERARCHICAL   = 3,
    CC_MODE_HYBRID         = 4,
    CC_MODE_REDUNDANT      = 5
} CloudControlMode;

typedef enum {
    CTASK_REGULATION    = 0,
    CTASK_OPTIMIZATION  = 1,
    CTASK_ESTIMATION    = 2,
    CTASK_DIAGNOSTICS   = 3,
    CTASK_PLANNING      = 4,
    CTASK_LEARNING      = 5
} ControlTaskType;

typedef enum {
    QOS_CRITICAL   = 0,
    QOS_HIGH       = 1,
    QOS_MEDIUM     = 2,
    QOS_LOW        = 3,
    QOS_BACKGROUND = 4
} QoSLevel;

typedef enum {
    CC_STATE_UNINITIALIZED = 0,
    CC_STATE_INITIALIZING  = 1,
    CC_STATE_RUNNING       = 2,
    CC_STATE_DEGRADED      = 3,
    CC_STATE_FALLBACK      = 4,
    CC_STATE_RECOVERING    = 5,
    CC_STATE_MAINTENANCE   = 6,
    CC_STATE_ERROR         = 7
} CloudControlState;

typedef struct {
    double   x[CCS_MAX_STATES];
    double   dx[CCS_MAX_STATES];
    double   u[CCS_MAX_INPUTS];
    double   y[CCS_MAX_OUTPUTS];
    int      dim;
    int      input_dim;
    int      output_dim;
    double   t;
    uint64_t seq;
} SystemState;

typedef struct {
    char     id[64];
    char     region[32];
    double   cpu_allocated;
    double   memory_allocated_mb;
    double   bandwidth_mbps;
    double   cpu_utilization;
    double   memory_utilization;
    double   avg_response_time_us;
    double   max_response_time_us;
    double   virtualization_overhead;
    int      task_queue_length;
    int      is_virtual;
    uint64_t uptime_cycles;
    uint64_t missed_deadlines;
} CloudNode;

typedef struct {
    char     id[64];
    double   cpu_allocated;
    double   memory_allocated_mb;
    double   local_loop_period_us;
    double   cloud_sync_period_us;
    double   max_holdover_time_us;
    double   avg_rtt_to_cloud_us;
    double   holdover_remaining_us;
    int      is_active;
    int      cloud_connected;
    uint64_t local_cycles;
    uint64_t cloud_cycles;
} EdgeNode;

typedef struct {
    char            id[64];
    ControlTaskType type;
    QoSLevel        qos;
    double          period_us;
    double          deadline_us;
    double          wcet_us;
    double          avg_et_us;
    double          release_time;
    double          completion_time;
    int             priority;
    int             state_dim;
    int             completed;
    uint64_t        instance_id;
} ControlTask;

typedef struct {
    QoSLevel level;
    double   max_e2e_latency_us;
    double   max_jitter_us;
    double   min_update_rate_hz;
    double   max_packet_loss_rate;
    double   control_accuracy;
    double   availability_target;
    int      require_encryption;
    int      require_redundancy;
} QoSProfile;

typedef struct {
    char              id[128];
    char              description[256];
    CloudControlMode  mode;
    CloudControlState state;
    QoSProfile        qos;
    double A[CCS_MAX_STATES][CCS_MAX_STATES];
    double B[CCS_MAX_STATES][CCS_MAX_INPUTS];
    double C[CCS_MAX_OUTPUTS][CCS_MAX_STATES];
    double D[CCS_MAX_OUTPUTS][CCS_MAX_INPUTS];
    int    n, m, p;
    double K[CCS_MAX_INPUTS][CCS_MAX_STATES];
    double L[CCS_MAX_STATES][CCS_MAX_OUTPUTS];
    double x_hat[CCS_MAX_STATES];
    SystemState  plant_state;
    SystemState  reference;
    CloudNode    cloud_node;
    EdgeNode     edge_node;
    double settling_time;
    double overshoot_percent;
    double steady_state_error;
    double ise;
    double iae;
    double control_effort;
    uint64_t total_cycles;
    double avg_rtt_us;
    double max_rtt_us;
    double avg_jitter_us;
    double packet_loss_rate;
    uint64_t packets_sent;
    uint64_t packets_lost;
    double history_t[CCS_MAX_HISTORY];
    double history_y[CCS_MAX_HISTORY];
    int    history_count;
    struct timespec created_at;
    struct timespec last_update;
} CloudControlSystem;

CloudControlSystem* ccs_create(const char *id, int n, int m, int p,
                                CloudControlMode mode);
void ccs_free(CloudControlSystem *ccs);
int  ccs_set_plant_model(CloudControlSystem *ccs,
                          const double *A, const double *B,
                          const double *C, const double *D);
int  ccs_set_controller(CloudControlSystem *ccs,
                         const double *K, const double *L);
int  ccs_set_qos(CloudControlSystem *ccs, const QoSProfile *qos);
int  ccs_set_reference(CloudControlSystem *ccs, const double *ref);
int  ccs_compute_control(CloudControlSystem *ccs, double *u_out);
int  ccs_update_observer(CloudControlSystem *ccs,
                          const double *y, double timestamp);
int  ccs_apply_control(CloudControlSystem *ccs, const double *u, double dt);
int  ccs_step(CloudControlSystem *ccs, const double *measurement,
               double ts, double dt, double delay_us);
int  ccs_compute_performance(CloudControlSystem *ccs);
int  ccs_is_stable(const CloudControlSystem *ccs, double delay_us);
double ccs_max_allowable_delay(const CloudControlSystem *ccs);
void ccs_print_state(const CloudControlSystem *ccs);
void ccs_print_metrics(const CloudControlSystem *ccs);
int  ccs_switch_mode(CloudControlSystem *ccs, CloudControlMode new_mode);
const char* ccs_mode_string(CloudControlMode mode);
const char* ccs_state_string(CloudControlState state);
int  ccs_register_cloud_node(CloudControlSystem *ccs, const CloudNode *node);
int  ccs_register_edge_node(CloudControlSystem *ccs, const EdgeNode *node);
int  ccs_update_cloud_metrics(CloudControlSystem *ccs,
                               double cpu_util, double mem_util,
                               double resp_time_us, int queue_len);
double ccs_get_cloud_load(const CloudControlSystem *ccs);
int  ccs_random_init(CloudControlSystem *ccs, int n, int m, int p,
                      unsigned int seed);
int  ccs_compare(const CloudControlSystem *a, const CloudControlSystem *b);
void ccs_reset_metrics(CloudControlSystem *ccs);

#endif
