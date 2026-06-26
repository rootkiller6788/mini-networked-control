/* cloud_control_resource.c - Cloud Resource Management for Control Systems */
#include "cloud_control_resource.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

/* -------- Resource Pool -------- */

ResourcePool* rpool_create(const char *id) {
    ResourcePool *pool = (ResourcePool*)calloc(1, sizeof(ResourcePool));
    if (!pool) return NULL;
    if (id) strncpy(pool->id, id, sizeof(pool->id)-1);
    pool->sched_policy = SCHED_RM;
    pool->total_utilization = 0.0;
    pool->slack = 1.0;
    return pool;
}

void rpool_free(ResourcePool *pool) { free(pool); }

int rpool_add_capacity(ResourcePool *pool, ResourceType type,
                        double total, double reserved) {
    if (!pool || pool->capacity_count >= CCS_MAX_RESOURCE_TYPES) return -1;
    int idx = pool->capacity_count++;
    pool->capacities[idx].type = type;
    pool->capacities[idx].total_capacity = total;
    pool->capacities[idx].reserved_capacity = reserved;
    pool->capacities[idx].used_capacity = 0.0;
    pool->capacities[idx].utilization = 0.0;
    pool->capacities[idx].overload_threshold = 0.85;
    return 0;
}

int rpool_allocate(ResourcePool *pool, const char *task_id,
                    const ResourceType *types, const double *amounts,
                    int count, QoSLevel qos) {
    if (!pool || !task_id || !types || !amounts || count <= 0) return -1;
    if (pool->allocation_count >= CCS_MAX_ALLOCATIONS) return -1;
    /* Check feasibility */
    for (int i = 0; i < count; i++) {
        int found = 0;
        for (int j = 0; j < pool->capacity_count; j++) {
            if (pool->capacities[j].type == types[i]) {
                double avail = pool->capacities[j].total_capacity
                             - pool->capacities[j].reserved_capacity
                             - pool->capacities[j].used_capacity;
                if (amounts[i] > avail) return -1;
                found = 1;
                break;
            }
        }
        if (!found) return -1;
    }
    /* Perform allocation */
    int idx = pool->allocation_count++;
    strncpy(pool->allocations[idx].task_id, task_id, 63);
    for (int i = 0; i < count && i < CCS_MAX_RESOURCE_TYPES; i++) {
        pool->allocations[idx].resources[i] = types[i];
        pool->allocations[idx].amounts[i] = amounts[i];
    }
    pool->allocations[idx].resource_count = count;
    pool->allocations[idx].qos = qos;
    pool->allocations[idx].is_guaranteed = (qos <= QOS_HIGH) ? 1 : 0;
    /* Update used capacity */
    for (int i = 0; i < count; i++) {
        for (int j = 0; j < pool->capacity_count; j++) {
            if (pool->capacities[j].type == types[i]) {
                pool->capacities[j].used_capacity += amounts[i];
                pool->capacities[j].utilization =
                    pool->capacities[j].used_capacity /
                    (pool->capacities[j].total_capacity -
                     pool->capacities[j].reserved_capacity + 1e-12);
                break;
            }
        }
    }
    pool->allocation_cycles++;
    return 0;
}

int rpool_release(ResourcePool *pool, const char *task_id) {
    if (!pool || !task_id) return -1;
    for (int i = 0; i < pool->allocation_count; i++) {
        if (strcmp(pool->allocations[i].task_id, task_id) == 0) {
            /* Return resources */
            for (int k = 0; k < pool->allocations[i].resource_count; k++) {
                for (int j = 0; j < pool->capacity_count; j++) {
                    if (pool->capacities[j].type == pool->allocations[i].resources[k]) {
                        pool->capacities[j].used_capacity -= pool->allocations[i].amounts[k];
                        if (pool->capacities[j].used_capacity < 0)
                            pool->capacities[j].used_capacity = 0;
                        break;
                    }
                }
            }
            /* Remove allocation by shifting */
            for (int j = i; j < pool->allocation_count - 1; j++)
                pool->allocations[j] = pool->allocations[j+1];
            pool->allocation_count--;
            return 0;
        }
    }
    return -1;
}

int rpool_rebalance(ResourcePool *pool) {
    if (!pool) return -1;
    double total_avail = 0.0, total_used = 0.0;
    for (int i = 0; i < pool->capacity_count; i++) {
        total_avail += pool->capacities[i].total_capacity -
                       pool->capacities[i].reserved_capacity;
        total_used += pool->capacities[i].used_capacity;
    }
    pool->total_utilization = (total_avail > 0) ? total_used / total_avail : 0;
    pool->slack = 1.0 - pool->total_utilization;
    if (pool->slack < 0) { pool->slack = 0; pool->overload_count++; }
    return 0;
}

double rpool_get_utilization(const ResourcePool *pool, ResourceType type) {
    if (!pool) return 0.0;
    for (int i = 0; i < pool->capacity_count; i++)
        if (pool->capacities[i].type == type)
            return pool->capacities[i].utilization;
    return 0.0;
}

double rpool_get_slack(const ResourcePool *pool) {
    return pool ? pool->slack : 0.0;
}

int rpool_is_schedulable(const ResourcePool *pool, SchedulingPolicy policy) {
    if (!pool) return 0;
    double U[CCS_MAX_ALLOCATIONS];
    int n = pool->allocation_count;
    if (n == 0) return 1;
    for (int i = 0; i < n; i++) {
        U[i] = 0.0;
        for (int j = 0; j < pool->capacity_count; j++)
            if (pool->capacities[j].type == RESOURCE_CPU)
                U[i] = pool->allocations[i].amounts[j] /
                       pool->capacities[j].total_capacity;
    }
    if (policy == SCHED_EDF) return sched_edf_bound(U, n);
    return sched_rate_monotonic_bound(U, n);
}
/* -------- Scheduling Analysis (Liu & Layland, 1973) -------- */

int sched_rate_monotonic_bound(const double *utilizations, int n) {
    if (!utilizations || n <= 0) return 0;
    double sum_u = 0.0;
    for (int i = 0; i < n; i++) sum_u += utilizations[i];
    double bound = n * (pow(2.0, 1.0/n) - 1.0);
    return (sum_u <= bound) ? 1 : 0;
}

int sched_edf_bound(const double *utilizations, int n) {
    if (!utilizations || n <= 0) return 0;
    double sum_u = 0.0;
    for (int i = 0; i < n; i++) sum_u += utilizations[i];
    return (sum_u <= 1.0 + 1e-12) ? 1 : 0;
}

void sched_deadline_monotonic(const double *deadlines, int *priorities, int n) {
    if (!deadlines || !priorities || n <= 0) return;
    /* Assign priority by sorting deadlines (shorter = higher priority = smaller number) */
    typedef struct { double d; int idx; } pair_t;
    pair_t *pairs = (pair_t*)malloc((size_t)n * sizeof(pair_t));
    if (!pairs) return;
    for (int i = 0; i < n; i++) { pairs[i].d = deadlines[i]; pairs[i].idx = i; }
    for (int i = 0; i < n-1; i++)
        for (int j = i+1; j < n; j++)
            if (pairs[i].d > pairs[j].d) { pair_t t = pairs[i]; pairs[i] = pairs[j]; pairs[j] = t; }
    for (int i = 0; i < n; i++) priorities[pairs[i].idx] = i;
    free(pairs);
}

double sched_response_time(const double *C, const double *T,
                            const int *priorities, int n, int task_idx) {
    if (!C || !T || !priorities || n <= 0 || task_idx < 0 || task_idx >= n)
        return -1.0;
    double R = C[task_idx];
    double R_prev;
    int max_iter = 1000;
    do {
        R_prev = R;
        double interference = 0.0;
        for (int j = 0; j < n; j++) {
            if (j != task_idx && priorities[j] < priorities[task_idx])
                interference += ceil(R_prev / T[j]) * C[j];
        }
        R = C[task_idx] + interference;
        if (R > T[task_idx] * 10.0) return -1.0;
        max_iter--;
    } while (fabs(R - R_prev) > 1e-9 && max_iter > 0);
    return R;
}
/* -------- Elastic Controller (MAPE-K loop) -------- */

ElasticController* elastic_create(ResourcePool *pool, ScalingPolicy policy) {
    if (!pool) return NULL;
    ElasticController *ec = (ElasticController*)calloc(1, sizeof(ElasticController));
    if (!ec) return NULL;
    ec->pool = pool;
    ec->policy = policy;
    ec->scale_up_threshold = 0.80;
    ec->scale_down_threshold = 0.30;
    ec->scale_factor = 2.0;
    ec->cooldown_period_s = 60.0;
    ec->current_instances = 1;
    ec->min_instances = 1;
    ec->max_instances = 16;
    ec->scale_up_latency_s = 30.0;
    ec->scale_down_latency_s = 10.0;
    return ec;
}

void elastic_free(ElasticController *ec) { free(ec); }

int elastic_set_thresholds(ElasticController *ec,
                            double up_threshold, double down_threshold) {
    if (!ec) return -1;
    if (up_threshold <= down_threshold || up_threshold > 1.0 || down_threshold < 0) return -1;
    ec->scale_up_threshold = up_threshold;
    ec->scale_down_threshold = down_threshold;
    return 0;
}

int elastic_evaluate(ElasticController *ec, double current_time) {
    if (!ec) return -1;
    if (current_time - ec->last_scale_time < ec->cooldown_period_s) return 0;
    rpool_rebalance(ec->pool);
    double util = ec->pool->total_utilization;
    if (util > ec->scale_up_threshold && ec->current_instances < ec->max_instances) {
        elastic_scale_up(ec);
        ec->last_scale_time = current_time;
    } else if (util < ec->scale_down_threshold && ec->current_instances > ec->min_instances) {
        elastic_scale_down(ec);
        ec->last_scale_time = current_time;
    }
    return 0;
}

int elastic_scale_up(ElasticController *ec) {
    if (!ec) return -1;
    int new_count = (int)(ec->current_instances * ec->scale_factor);
    if (new_count <= ec->current_instances) new_count = ec->current_instances + 1;
    if (new_count > ec->max_instances) new_count = ec->max_instances;
    ec->current_instances = new_count;
    ec->scale_up_count++;
    return 0;
}

int elastic_scale_down(ElasticController *ec) {
    if (!ec) return -1;
    int new_count = ec->current_instances / 2;
    if (new_count < ec->min_instances) new_count = ec->min_instances;
    ec->current_instances = new_count;
    ec->scale_down_count++;
    return 0;
}

int elastic_get_instance_count(const ElasticController *ec) {
    return ec ? ec->current_instances : -1;
}

void elastic_print_state(const ElasticController *ec) {
    if (!ec) return;
    printf("ElasticController: instances=%d policy=%d\n",
           ec->current_instances, (int)ec->policy);
    printf("  Thresholds: up=%.2f down=%.2f cooldown=%.1fs\n",
           ec->scale_up_threshold, ec->scale_down_threshold, ec->cooldown_period_s);
    printf("  Scale events: up=%llu down=%llu\n",
           (unsigned long long)ec->scale_up_count, (unsigned long long)ec->scale_down_count);
}

/* -------- Workload Prediction (Holt-Winters Exponential Smoothing) -------- */

WorkloadPredictor* wp_create(void) {
    WorkloadPredictor *wp = (WorkloadPredictor*)calloc(1, sizeof(WorkloadPredictor));
    if (!wp) return NULL;
    wp->alpha = 0.3; wp->beta = 0.1; wp->gamma = 0.05;
    wp->seasonal_period = 24;
    wp->forecast_horizon = 1.0;
    for (int i = 0; i < 64; i++) wp->seasonal[i] = 1.0;
    return wp;
}

void wp_free(WorkloadPredictor *wp) { free(wp); }

int wp_add_sample(WorkloadPredictor *wp, double workload, double timestamp) {
    if (!wp || wp->history_count >= CCS_MAX_HISTORY) return -1;
    int idx = wp->history_count;
    wp->history[idx] = workload;
    wp->history_t[idx] = timestamp;
    wp->history_count++;
    /* Update Holt-Winters components */
    if (idx == 0) {
        wp->level = workload;
        wp->trend = 0.0;
    } else {
        double old_level = wp->level;
        int s_idx = idx % wp->seasonal_period;
        wp->level = wp->alpha * workload / (wp->seasonal[s_idx] + 1e-12)
                    + (1.0 - wp->alpha) * (old_level + wp->trend);
        wp->trend = wp->beta * (wp->level - old_level)
                    + (1.0 - wp->beta) * wp->trend;
        wp->seasonal[s_idx] = wp->gamma * workload / (wp->level + 1e-12)
                              + (1.0 - wp->gamma) * wp->seasonal[s_idx];
    }
    /* Compute forecast error */
    double forecast = wp->level + wp->trend * wp->forecast_horizon;
    wp->forecast_error = workload - forecast;
    wp->mape = 0.95 * wp->mape + 0.05 * fabs(wp->forecast_error / (workload + 1e-12));
    wp->last_forecast = forecast;
    return 0;
}

double wp_forecast_holtwinters(WorkloadPredictor *wp, double horizon) {
    if (!wp || wp->history_count < 2) return 0.0;
    int s_idx = (wp->history_count) % wp->seasonal_period;
    return (wp->level + wp->trend * horizon) * wp->seasonal[s_idx];
}

double wp_forecast_ema(const WorkloadPredictor *wp, double alpha) {
    if (!wp || wp->history_count == 0) return 0.0;
    double ema = wp->history[0];
    for (int i = 1; i < wp->history_count; i++)
        ema = alpha * wp->history[i] + (1.0 - alpha) * ema;
    return ema;
}

double wp_forecast_ar(const WorkloadPredictor *wp, int order) {
    if (!wp || wp->history_count < order + 1 || order <= 0) return 0.0;
    /* Simple AR(order): compute phi via Yule-Walker (first-order only for simplicity) */
    int n = wp->history_count;
    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += wp->history[i];
    mean /= n;
    double num = 0.0, den = 0.0;
    for (int i = order; i < n; i++) {
        num += (wp->history[i] - mean) * (wp->history[i-order] - mean);
        den += (wp->history[i-order] - mean) * (wp->history[i-order] - mean);
    }
    double phi = (den > 1e-12) ? num/den : 0.0;
    return mean + phi * (wp->history[n-1] - mean);
}

double wp_get_mape(const WorkloadPredictor *wp) {
    return wp ? wp->mape : -1.0;
}

void wp_reset(WorkloadPredictor *wp) {
    if (!wp) return;
    wp->history_count = 0; wp->level = 0; wp->trend = 0; wp->mape = 0;
    for (int i = 0; i < 64; i++) wp->seasonal[i] = 1.0;
}

/* -------- QoS-Aware Resource Mapping -------- */

void qos_resource_map(QoSLevel qos, double *cpu_out,
                       double *mem_out, double *bw_out) {
    if (!cpu_out || !mem_out || !bw_out) return;
    switch (qos) {
    case QOS_CRITICAL: *cpu_out = 2.0; *mem_out = 2048.0; *bw_out = 1000.0; break;
    case QOS_HIGH:     *cpu_out = 1.0; *mem_out = 1024.0; *bw_out = 500.0;  break;
    case QOS_MEDIUM:   *cpu_out = 0.5; *mem_out = 512.0;  *bw_out = 100.0;  break;
    case QOS_LOW:      *cpu_out = 0.25; *mem_out = 256.0; *bw_out = 50.0;  break;
    case QOS_BACKGROUND: *cpu_out = 0.1; *mem_out = 128.0; *bw_out = 10.0; break;
    default:           *cpu_out = 0.5; *mem_out = 512.0;  *bw_out = 100.0;  break;
    }
}

int qos_isolation_check(const ResourcePool *pool, const char *task_id) {
    if (!pool || !task_id) return 0;
    /* Find the task allocation */
    int alloc_idx = -1;
    for (int i = 0; i < pool->allocation_count; i++) {
        if (strcmp(pool->allocations[i].task_id, task_id) == 0) { alloc_idx = i; break; }
    }
    if (alloc_idx < 0) return 0;
    /* Check if any other high-priority task can preempt */
    QoSLevel my_qos = pool->allocations[alloc_idx].qos;
    for (int i = 0; i < pool->allocation_count; i++) {
        if (i != alloc_idx && pool->allocations[i].qos < my_qos)
            return 0; /* Higher priority task exists, isolation not guaranteed */
    }
    return 1;
}
/* ============================================================================
 * Additional Utility Functions for Resource Management
 * ============================================================================ */

/* -------- Capacity Planning -------- */

double rpool_estimate_required_capacity(const ResourcePool *pool,
                                         double target_utilization) {
    if (!pool || target_utilization <= 0 || target_utilization >= 1.0) return -1.0;
    double total_used = 0.0;
    for (int i = 0; i < pool->capacity_count; i++)
        total_used += pool->capacities[i].used_capacity;
    return total_used / target_utilization;
}

int rpool_get_overload_count(const ResourcePool *pool) {
    return pool ? pool->overload_count : -1;
}

int rpool_get_allocation_count(const ResourcePool *pool) {
    return pool ? pool->allocation_count : -1;
}

/* -------- Resource Fragmentation Analysis -------- */

double rpool_fragmentation_index(const ResourcePool *pool) {
    if (!pool || pool->allocation_count <= 1) return 0.0;
    double max_alloc = 0.0, sum_alloc = 0.0;
    for (int i = 0; i < pool->allocation_count; i++) {
        double a = 0.0;
        for (int j = 0; j < pool->allocations[i].resource_count; j++)
            a += pool->allocations[i].amounts[j];
        sum_alloc += a;
        if (a > max_alloc) max_alloc = a;
    }
    return 1.0 - (max_alloc / (sum_alloc + 1e-12));
}

/* -------- Admission Control -------- */

int rpool_admission_test(const ResourcePool *pool,
                          const double *requested_amounts,
                          const ResourceType *requested_types,
                          int count) {
    if (!pool || !requested_amounts || !requested_types || count <= 0) return 0;
    for (int i = 0; i < count; i++) {
        int found = 0;
        for (int j = 0; j < pool->capacity_count; j++) {
            if (pool->capacities[j].type == requested_types[i]) {
                double slack = pool->capacities[j].total_capacity
                             - pool->capacities[j].reserved_capacity
                             - pool->capacities[j].used_capacity;
                if (requested_amounts[i] > slack) return 0;
                found = 1; break;
            }
        }
        if (!found) return 0;
    }
    return 1;
}

/* -------- ComputeNode Management -------- */

int cnode_init(ComputeNode *node, const char *id, double cpu, double mem,
               double storage, double bw) {
    if (!node || !id) return -1;
    memset(node, 0, sizeof(*node));
    strncpy(node->id, id, sizeof(node->id)-1);
    node->cpu_cores = cpu; node->memory_mb = mem;
    node->storage_gb = storage; node->network_mbps = bw;
    node->is_active = 1;
    return 0;
}

void cnode_update_utilization(ComputeNode *node, double cpu, double mem) {
    if (!node) return;
    node->cpu_util = cpu; node->mem_util = mem;
}

int cnode_can_accept_task(const ComputeNode *node, double cpu_req,
                           double mem_req) {
    if (!node || !node->is_active) return 0;
    double cpu_avail = node->cpu_cores * (1.0 - node->cpu_util);
    double mem_avail = node->memory_mb * (1.0 - node->mem_util);
    return (cpu_avail >= cpu_req && mem_avail >= mem_req) ? 1 : 0;
}

void cnode_print(const ComputeNode *node) {
    if (!node) return;
    printf("ComputeNode[%s]: cpu=%.1f/%.1f mem=%.0f/%.0fMB temp=%.1fC active=%d\n",
           node->id, node->cpu_util*node->cpu_cores, node->cpu_cores,
           node->mem_util*node->memory_mb, node->memory_mb,
           node->temperature_c, node->is_active);
}


/* -------- Elastic Cost Analysis -------- */

double elastic_cost_estimate(const ElasticController *ec,
                              double cost_per_instance_hour) {
    if (!ec) return -1.0;
    return ec->current_instances * cost_per_instance_hour;
}

double elastic_efficiency_ratio(const ElasticController *ec) {
    if (!ec || ec->scale_up_count + ec->scale_down_count == 0) return 1.0;
    return (double)ec->scale_down_count /
           (double)(ec->scale_up_count + ec->scale_down_count + 1);
}

/* -------- Statistical Workload Analysis -------- */

double wp_mean_workload(const WorkloadPredictor *wp) {
    if (!wp || wp->history_count == 0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < wp->history_count; i++) sum += wp->history[i];
    return sum / wp->history_count;
}

double wp_variance_workload(const WorkloadPredictor *wp) {
    if (!wp || wp->history_count < 2) return 0.0;
    double mean = wp_mean_workload(wp), sum = 0.0;
    for (int i = 0; i < wp->history_count; i++) {
        double d = wp->history[i] - mean; sum += d*d;
    }
    return sum / (wp->history_count - 1);
}

double wp_peak_to_mean_ratio(const WorkloadPredictor *wp) {
    if (!wp || wp->history_count == 0) return 0.0;
    double mean = wp_mean_workload(wp);
    if (mean < 1e-12) return 0.0;
    double peak = wp->history[0];
    for (int i = 1; i < wp->history_count; i++)
        if (wp->history[i] > peak) peak = wp->history[i];
    return peak / mean;
}

int wp_set_parameters(WorkloadPredictor *wp, double alpha, double beta,
                       double gamma, int seasonal_period) {
    if (!wp) return -1;
    if (alpha < 0 || alpha > 1 || beta < 0 || beta > 1 || gamma < 0 || gamma > 1)
        return -1;
    wp->alpha = alpha; wp->beta = beta; wp->gamma = gamma;
    wp->seasonal_period = seasonal_period > 0 ? seasonal_period : 1;
    return 0;
}
