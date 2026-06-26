#ifndef CLOUD_CONTROL_RESOURCE_H
#define CLOUD_CONTROL_RESOURCE_H

#include "cloud_control_core.h"
#include <stdint.h>

#define CCS_MAX_RESOURCE_TYPES   8
#define CCS_MAX_COMPUTE_NODES    256
#define CCS_MAX_ALLOCATIONS      512
#define CCS_DEFAULT_CPU_QUOTA    1.0
#define CCS_DEFAULT_MEM_MB       1024.0
#define CCS_DEFAULT_BW_MBPS      100.0

typedef enum {
    RESOURCE_CPU      = 0,
    RESOURCE_MEMORY   = 1,
    RESOURCE_NETWORK  = 2,
    RESOURCE_STORAGE  = 3,
    RESOURCE_GPU      = 4,
    RESOURCE_IO       = 5,
    RESOURCE_CACHE    = 6,
    RESOURCE_POWER    = 7
} ResourceType;

typedef enum {
    SCHED_RM    = 0,
    SCHED_EDF   = 1,
    SCHED_DM    = 2,
    SCHED_FIFO  = 3,
    SCHED_RR    = 4,
    SCHED_CBS   = 5
} SchedulingPolicy;

typedef enum {
    SCALE_NONE         = 0,
    SCALE_REACTIVE     = 1,
    SCALE_PREDICTIVE   = 2,
    SCALE_SCHEDULE     = 3,
    SCALE_ADAPTIVE     = 4
} ScalingPolicy;

typedef struct {
    ResourceType type;
    double       total_capacity;
    double       reserved_capacity;
    double       used_capacity;
    double       utilization;
    double       overload_threshold;
} ResourceCapacity;

typedef struct {
    char         task_id[64];
    ResourceType resources[CCS_MAX_RESOURCE_TYPES];
    double       amounts[CCS_MAX_RESOURCE_TYPES];
    int          resource_count;
    double       priority_weight;
    double       min_amount;
    double       optimal_amount;
    QoSLevel     qos;
    int          is_guaranteed;
    double       allocated_at;
} ResourceAllocation;

typedef struct {
    char              id[64];
    ResourceCapacity  capacities[CCS_MAX_RESOURCE_TYPES];
    int               capacity_count;
    ResourceAllocation allocations[CCS_MAX_ALLOCATIONS];
    int               allocation_count;
    SchedulingPolicy  sched_policy;
    double            total_utilization;
    double            slack;
    int               overload_count;
    uint64_t          allocation_cycles;
} ResourcePool;

typedef struct {
    char     id[64];
    char     ip_addr[64];
    double   cpu_cores;
    double   cpu_freq_ghz;
    double   memory_mb;
    double   storage_gb;
    double   network_mbps;
    int      gpu_count;
    double   power_watts;
    double   cpu_util;
    double   mem_util;
    double   temperature_c;
    int      is_active;
    int      is_virtual;
    char     hypervisor[32];
    uint64_t tasks_completed;
    uint64_t tasks_missed;
    double   avg_latency_us;
} ComputeNode;

typedef struct {
    double   history[CCS_MAX_HISTORY];
    double   history_t[CCS_MAX_HISTORY];
    int      history_count;
    double   level;
    double   trend;
    double   seasonal[64];
    int      seasonal_period;
    double   alpha;
    double   beta;
    double   gamma;
    double   forecast_horizon;
    double   last_forecast;
    double   forecast_error;
    double   mape;
} WorkloadPredictor;

typedef struct {
    ScalingPolicy     policy;
    ResourcePool     *pool;
    WorkloadPredictor predictor;
    double            scale_up_threshold;
    double            scale_down_threshold;
    double            scale_factor;
    double            cooldown_period_s;
    double            last_scale_time;
    int               current_instances;
    int               min_instances;
    int               max_instances;
    double            scale_up_latency_s;
    double            scale_down_latency_s;
    uint64_t          scale_up_count;
    uint64_t          scale_down_count;
} ElasticController;

ResourcePool* rpool_create(const char *id);
void rpool_free(ResourcePool *pool);
int  rpool_add_capacity(ResourcePool *pool, ResourceType type,
                         double total, double reserved);
int  rpool_allocate(ResourcePool *pool, const char *task_id,
                     const ResourceType *types, const double *amounts,
                     int count, QoSLevel qos);
int  rpool_release(ResourcePool *pool, const char *task_id);
int  rpool_rebalance(ResourcePool *pool);
double rpool_get_utilization(const ResourcePool *pool, ResourceType type);
double rpool_get_slack(const ResourcePool *pool);
int  rpool_is_schedulable(const ResourcePool *pool, SchedulingPolicy policy);
int  sched_rate_monotonic_bound(const double *utilizations, int n);
int  sched_edf_bound(const double *utilizations, int n);
void sched_deadline_monotonic(const double *deadlines, int *priorities, int n);
double sched_response_time(const double *C, const double *T,
                            const int *priorities, int n, int task_idx);
ElasticController* elastic_create(ResourcePool *pool, ScalingPolicy policy);
void elastic_free(ElasticController *ec);
int  elastic_set_thresholds(ElasticController *ec,
                             double up_threshold, double down_threshold);
int  elastic_evaluate(ElasticController *ec, double current_time);
int  elastic_scale_up(ElasticController *ec);
int  elastic_scale_down(ElasticController *ec);
int  elastic_get_instance_count(const ElasticController *ec);
void elastic_print_state(const ElasticController *ec);
WorkloadPredictor* wp_create(void);
void wp_free(WorkloadPredictor *wp);
int  wp_add_sample(WorkloadPredictor *wp, double workload, double timestamp);
double wp_forecast_holtwinters(WorkloadPredictor *wp, double horizon);
double wp_forecast_ema(const WorkloadPredictor *wp, double alpha);
double wp_forecast_ar(const WorkloadPredictor *wp, int order);
double wp_get_mape(const WorkloadPredictor *wp);
void wp_reset(WorkloadPredictor *wp);
void qos_resource_map(QoSLevel qos, double *cpu_out,
                       double *mem_out, double *bw_out);
int  qos_isolation_check(const ResourcePool *pool, const char *task_id);


/* Extended resource functions */
double rpool_estimate_required_capacity(const ResourcePool *pool, double target_utilization);
int rpool_get_overload_count(const ResourcePool *pool);
int rpool_get_allocation_count(const ResourcePool *pool);
double rpool_fragmentation_index(const ResourcePool *pool);
int rpool_admission_test(const ResourcePool *pool, const double *requested_amounts, const ResourceType *requested_types, int count);
int cnode_init(ComputeNode *node, const char *id, double cpu, double mem, double storage, double bw);
void cnode_update_utilization(ComputeNode *node, double cpu, double mem);
int cnode_can_accept_task(const ComputeNode *node, double cpu_req, double mem_req);
void cnode_print(const ComputeNode *node);
double elastic_cost_estimate(const ElasticController *ec, double cost_per_instance_hour);
double elastic_efficiency_ratio(const ElasticController *ec);
double wp_mean_workload(const WorkloadPredictor *wp);
double wp_variance_workload(const WorkloadPredictor *wp);
double wp_peak_to_mean_ratio(const WorkloadPredictor *wp);
int wp_set_parameters(WorkloadPredictor *wp, double alpha, double beta, double gamma, int seasonal_period);

#endif /* CLOUD_CONTROL_RESOURCE_H */
