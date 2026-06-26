#ifndef CLOUD_CONTROL_NETWORK_H
#define CLOUD_CONTROL_NETWORK_H
#include "cloud_control_core.h"
#include "cloud_control_resource.h"
#include <stdint.h>

#define CCS_MAX_NETWORK_PATHS     16
#define CCS_MAX_PACKET_SIZE       1500
#define CCS_MAX_TOPOLOGY_NODES    512

typedef enum {
    NET_CONDITION_NORMAL      = 0,
    NET_CONDITION_CONGESTED   = 1,
    NET_CONDITION_HIGH_LOSS   = 2,
    NET_CONDITION_PARTITIONED = 3,
    NET_CONDITION_DDOS        = 4,
    NET_CONDITION_DEGRADED    = 5
} NetworkCondition;

typedef enum {
    ROUTE_PRIMARY     = 0,
    ROUTE_BACKUP      = 1,
    ROUTE_LOAD_BALANCE = 2,
    ROUTE_LOW_LATENCY = 3,
    ROUTE_HIGH_BW     = 4
} RouteType;

typedef struct {
    char     src_id[64], dst_id[64];
    double   bandwidth_mbps, latency_us, jitter_us, packet_loss_rate, utilization;
    int      is_wireless, mtu, buffer_size_pkts, buffer_occupancy;
    uint64_t bytes_sent, bytes_received, packets_dropped;
} NetworkLink;

typedef struct {
    uint64_t seq_num, timestamp_us;
    char     src_id[64], dst_id[64];
    int      payload_type, payload_length, priority, ttl;
    double   sent_at, received_at;
} PacketHeader;

typedef struct {
    PacketHeader header;
    uint8_t      payload[CCS_MAX_PACKET_SIZE];
    int          is_control;
    double       value[CCS_MAX_STATES];
    int          value_dim;
} ControlPacket;

typedef struct {
    NetworkLink links[CCS_MAX_NETWORK_PATHS];
    int         link_count;
    RouteType   type;
    double      total_latency_us, total_bandwidth_mbps, reliability;
    int         hop_count;
} MultiPathRoute;

typedef struct {
    char              id[64];
    int               node_count;
    char              node_ids[CCS_MAX_TOPOLOGY_NODES][64];
    NetworkLink       adjacency[CCS_MAX_TOPOLOGY_NODES][CCS_MAX_NETWORK_PATHS];
    int               adj_count[CCS_MAX_TOPOLOGY_NODES];
    double            avg_path_length, diameter, clustering_coefficient, connectivity;
    NetworkCondition  condition;
} EdgeCloudTopology;

typedef struct {
    char          id[64];
    ControlPacket buffer[CCS_MAX_TASK_QUEUE];
    int           head, tail, count, max_size;
    uint64_t      packets_enqueued, packets_dequeued, packets_dropped;
} PacketBuffer;

typedef struct {
    char       link_id[64];
    double     total_bandwidth_mbps, allocated_bandwidth_mbps, available_bandwidth_mbps;
    double     allocation[CCS_MAX_ALLOCATIONS];
    char       alloc_ids[CCS_MAX_ALLOCATIONS][64];
    int        alloc_count;
    QoSLevel   min_qos;
} BandwidthAllocator;

EdgeCloudTopology* topology_create(const char *id);
void topology_free(EdgeCloudTopology *topo);
int  topology_add_node(EdgeCloudTopology *topo, const char *node_id);
int  topology_add_link(EdgeCloudTopology *topo, const char *src,
                        const char *dst, double bw_mbps, double lat_us);
int  topology_find_shortest_path(const EdgeCloudTopology *topo,
                                  const char *src, const char *dst,
                                  MultiPathRoute *route);
double topology_reliability(const EdgeCloudTopology *topo,
                             const char *src, const char *dst);
int  topology_is_connected(const EdgeCloudTopology *topo);
int  topology_detect_partition(EdgeCloudTopology *topo);
void topology_print(const EdgeCloudTopology *topo);
int  topology_import_from_matrix(EdgeCloudTopology *topo,
                                  const double *latency_matrix, int n);

PacketBuffer* pktbuf_create(const char *id, int max_size);
void pktbuf_free(PacketBuffer *buf);
int  pktbuf_enqueue(PacketBuffer *buf, const ControlPacket *pkt);
int  pktbuf_dequeue(PacketBuffer *buf, ControlPacket *pkt);
int  pktbuf_peek(const PacketBuffer *buf, ControlPacket *pkt);
int  pktbuf_is_full(const PacketBuffer *buf);
int  pktbuf_is_empty(const PacketBuffer *buf);
int  pktbuf_count(const PacketBuffer *buf);
void pktbuf_flush(PacketBuffer *buf);

int  netlink_update_metrics(NetworkLink *link, double latency, double jitter,
                              double loss_rate, double utilization);
int  netlink_is_congested(const NetworkLink *link);
double netlink_effective_bandwidth(const NetworkLink *link);
double netlink_bandwidth_delay_product(const NetworkLink *link);
void netlink_print(const NetworkLink *link);

BandwidthAllocator* bw_alloc_create(const char *link_id, double total_bw);
void bw_alloc_free(BandwidthAllocator *ba);
int  bw_alloc_reserve(BandwidthAllocator *ba, const char *consumer_id,
                       double amount, QoSLevel qos);
int  bw_alloc_release(BandwidthAllocator *ba, const char *consumer_id);
double bw_alloc_get_available(const BandwidthAllocator *ba);
int  bw_alloc_is_feasible(const BandwidthAllocator *ba, double amount);

int  edge_to_cloud_send(CloudControlSystem *ccs, const double *state,
                          int dim, double timestamp);
int  cloud_to_edge_send(CloudControlSystem *ccs, const double *control,
                          int dim, double timestamp);
int  edge_holdover_check(const EdgeNode *edge);
int  edge_holdover_enter(EdgeNode *edge);
int  edge_holdover_exit(EdgeNode *edge);
double edge_holdover_remaining(const EdgeNode *edge);
double net_condition_delay_factor(NetworkCondition cond);
double net_condition_loss_factor(NetworkCondition cond);
const char* net_condition_string(NetworkCondition cond);
NetworkCondition net_detect_condition(double loss_rate, double jitter, double rtt);
int  mproute_add_link(MultiPathRoute *route, const NetworkLink *link);
double mproute_total_latency(const MultiPathRoute *route);
double mproute_effective_bw(const MultiPathRoute *route);
int  mproute_compare(const MultiPathRoute *a, const MultiPathRoute *b);
double cloud_latency_estimate(const char *region_from, const char *region_to);
double edge_to_cloud_rtt_model(double edge_distance_km, const char *cloud_region);
double cloud_inter_region_rtt(const char *region_a, const char *region_b);
int  pkt_serialize(const ControlPacket *pkt, uint8_t *buffer, int buf_size);
int  pkt_deserialize(const uint8_t *buffer, int buf_size, ControlPacket *pkt);
int  pkt_create_measurement(ControlPacket *pkt, const char *src,
                              const double *y, int dim, uint64_t seq);
int  pkt_create_control(ControlPacket *pkt, const char *src,
                          const double *u, int dim, uint64_t seq);

#endif

/* ============================================================================
 * L7-L9: Extended API - Time Sync, AoI, Reorder, Token Bucket, TSN
 * ============================================================================ */

/* --- Age of Information (AoI) --- */
void* ccs_aoi_create(double threshold_s);
void  ccs_aoi_free(void *aoi);
void  ccs_aoi_generate(void *aoi, double t);
double ccs_aoi_deliver(void *aoi, double gen_t, double del_t);
double ccs_aoi_current(void *aoi, double now);
double ccs_aoi_peak(void *aoi);
double ccs_aoi_average(void *aoi);
double ccs_aoi_delivery_ratio(void *aoi);
int    ccs_aoi_is_violated(void *aoi, double now);
void   ccs_aoi_reset(void *aoi);
double ccs_aoi_optimal_rate(double service_rate);

/* --- NTP Clock Synchronization --- */
void*  ccs_ntp_create(void);
void   ccs_ntp_free(void *ncs);
int    ccs_ntp_process(void *ncs, double t1, double t2, double t3, double t4);
double ccs_ntp_get_offset_us(void *ncs);
double ccs_ntp_get_delay_us(void *ncs);
double ccs_ntp_get_uncertainty_us(void *ncs);
int    ccs_ntp_is_locked(void *ncs);
double ccs_ntp_skew_ppb(void *ncs);
double ccs_ntp_compensate(void *ncs, double remote_ts);
double ccs_ntp_holdover_time(void *ncs, double max_offset_us);

/* --- Packet Reorder Buffer --- */
void*  ccs_reorder_create(int win_size, double max_wait_s);
void   ccs_reorder_free(void *rb);
int    ccs_reorder_insert(void *rb, uint64_t seq, const double *data, int dim,
                           double arr_t, double gen_t);
int    ccs_reorder_check_timeout(void *rb, double now);
double ccs_reorder_density(void *rb);
double ccs_reorder_avg_wait(void *rb);
void   ccs_reorder_stats(void *rb, uint64_t *rx, uint64_t *inorder,
                          uint64_t *late, uint64_t *dups, uint64_t *drops);

/* --- Token Bucket Traffic Shaper --- */
void*  ccs_token_bucket_create(double rate, double burst);
void   ccs_token_bucket_free(void *tb);
int    ccs_token_bucket_test(void *tb, double size, double now);
double ccs_token_bucket_conformance(void *tb);

/* --- Delay-Gradient Congestion Detection --- */
void*  ccs_congestion_create(int threshold);
void   ccs_congestion_free(void *dc);
int    ccs_congestion_update(void *dc, double rtt_us);
double ccs_congestion_severity(void *dc);

/* --- Weighted Round-Robin Multi-Path Scheduler --- */
void*  ccs_wrr_create(void);
void   ccs_wrr_free(void *ws);
int    ccs_wrr_add_path(void *ws, int id, double bw_mbps);
int    ccs_wrr_schedule(void *ws, int pkt_size);

/* --- One-Way Delay / TSN Guard Band --- */
double ccs_owd_minrtt(const double *rtt, int n);
double ccs_owd_asymmetry(const double *rtt, int n);
double net_aoi_theoretical_avg(double arrival_rate, double service_rate);
double net_tsn_guard_band(int max_frame_bytes, double link_rate_bps,
                           double sync_uncertainty_us);
