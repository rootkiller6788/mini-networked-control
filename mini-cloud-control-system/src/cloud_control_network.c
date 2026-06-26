/* cloud_control_network.c - Network Layer for Cloud Control Systems */
#include "cloud_control_network.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>

/* -------- EdgeCloudTopology -------- */

EdgeCloudTopology* topology_create(const char *id) {
    EdgeCloudTopology *topo = (EdgeCloudTopology*)calloc(1, sizeof(EdgeCloudTopology));
    if (!topo) return NULL;
    if (id) strncpy(topo->id, id, sizeof(topo->id)-1);
    topo->condition = NET_CONDITION_NORMAL;
    return topo;
}

void topology_free(EdgeCloudTopology *topo) { free(topo); }

int topology_add_node(EdgeCloudTopology *topo, const char *node_id) {
    if (!topo || !node_id || topo->node_count >= CCS_MAX_TOPOLOGY_NODES) return -1;
    for (int i = 0; i < topo->node_count; i++)
        if (strcmp(topo->node_ids[i], node_id) == 0) return -1;
    strncpy(topo->node_ids[topo->node_count], node_id, 63);
    topo->adj_count[topo->node_count] = 0;
    topo->node_count++;
    return 0;
}

static int find_node_idx(const EdgeCloudTopology *topo, const char *id) {
    for (int i = 0; i < topo->node_count; i++)
        if (strcmp(topo->node_ids[i], id) == 0) return i;
    return -1;
}

int topology_add_link(EdgeCloudTopology *topo, const char *src,
                       const char *dst, double bw_mbps, double lat_us) {
    if (!topo || !src || !dst) return -1;
    int si = find_node_idx(topo, src), di = find_node_idx(topo, dst);
    if (si < 0 || di < 0) return -1;
    if (topo->adj_count[si] >= CCS_MAX_NETWORK_PATHS) return -1;
    int idx = topo->adj_count[si]++;
    strncpy(topo->adjacency[si][idx].src_id, src, 63);
    strncpy(topo->adjacency[si][idx].dst_id, dst, 63);
    topo->adjacency[si][idx].bandwidth_mbps = bw_mbps;
    topo->adjacency[si][idx].latency_us = lat_us;
    topo->adjacency[si][idx].jitter_us = lat_us * 0.1;
    topo->adjacency[si][idx].packet_loss_rate = 0.001;
    return 0;
}

int topology_remove_node(EdgeCloudTopology *topo, const char *node_id) {
    if (!topo || !node_id) return -1;
    int idx = find_node_idx(topo, node_id);
    if (idx < 0) return -1;
    for (int i = idx; i < topo->node_count - 1; i++) {
        strncpy(topo->node_ids[i], topo->node_ids[i+1], 63);
        topo->adj_count[i] = topo->adj_count[i+1];
        for (int j = 0; j < topo->adj_count[i]; j++)
            topo->adjacency[i][j] = topo->adjacency[i+1][j];
    }
    topo->node_count--;
    return 0;
}

/* Dijkstra shortest path for latency */
int topology_find_shortest_path(const EdgeCloudTopology *topo,
                                 const char *src, const char *dst,
                                 MultiPathRoute *route) {
    if (!topo || !src || !dst || !route) return -1;
    int n = topo->node_count;
    int si = find_node_idx(topo, src), di = find_node_idx(topo, dst);
    if (si < 0 || di < 0) return -1;
    double *dist = (double*)malloc((size_t)n * sizeof(double));
    int *prev = (int*)malloc((size_t)n * sizeof(int));
    int *visited = (int*)calloc((size_t)n, sizeof(int));
    if (!dist || !prev || !visited) { free(dist); free(prev); free(visited); return -1; }
    for (int i = 0; i < n; i++) { dist[i] = 1e300; prev[i] = -1; }
    dist[si] = 0;
    for (int iter = 0; iter < n; iter++) {
        int u = -1; double min_d = 1e300;
        for (int i = 0; i < n; i++)
            if (!visited[i] && dist[i] < min_d) { min_d = dist[i]; u = i; }
        if (u < 0 || u == di) break;
        visited[u] = 1;
        for (int j = 0; j < topo->adj_count[u]; j++) {
            int v = find_node_idx(topo, topo->adjacency[u][j].dst_id);
            if (v < 0) continue;
            double nd = dist[u] + topo->adjacency[u][j].latency_us;
            if (nd < dist[v]) { dist[v] = nd; prev[v] = u; }
        }
    }
    if (dist[di] >= 1e300) { free(dist); free(prev); free(visited); return -1; }
    /* Trace path */
    memset(route, 0, sizeof(*route));
    route->total_latency_us = dist[di];
    int path[CCS_MAX_TOPOLOGY_NODES], path_len = 0;
    for (int cur = di; cur != -1; cur = prev[cur]) path[path_len++] = cur;
    route->hop_count = path_len - 1;
    /* Compute effective bandwidth (min along path) */
    double min_bw = 1e300;
    for (int i = path_len-1; i > 0; i--) {
        int u = path[i], v_idx = path[i-1];
        for (int j = 0; j < topo->adj_count[u]; j++) {
            if (find_node_idx(topo, topo->adjacency[u][j].dst_id) == v_idx) {
                if (topo->adjacency[u][j].bandwidth_mbps < min_bw)
                    min_bw = topo->adjacency[u][j].bandwidth_mbps;
                if (route->link_count < CCS_MAX_NETWORK_PATHS)
                    route->links[route->link_count++] = topo->adjacency[u][j];
                break;
            }
        }
    }
    route->total_bandwidth_mbps = (min_bw < 1e300) ? min_bw : 0;
    route->reliability = 1.0;
    for (int i = 0; i < route->link_count; i++)
        route->reliability *= (1.0 - route->links[i].packet_loss_rate);
    free(dist); free(prev); free(visited);
    return 0;
}

double topology_reliability(const EdgeCloudTopology *topo,
                             const char *src, const char *dst) {
    MultiPathRoute route;
    if (topology_find_shortest_path(topo, src, dst, &route) < 0) return 0.0;
    return route.reliability;
}

int topology_is_connected(const EdgeCloudTopology *topo) {
    if (!topo || topo->node_count <= 1) return 1;
    /* BFS from node 0 */
    int *visited = (int*)calloc((size_t)topo->node_count, sizeof(int));
    if (!visited) return 0;
    int *queue = (int*)malloc((size_t)topo->node_count * sizeof(int));
    if (!queue) { free(visited); return 0; }
    int qh = 0, qt = 0;
    visited[0] = 1; queue[qt++] = 0;
    while (qh < qt) {
        int u = queue[qh++];
        for (int j = 0; j < topo->adj_count[u]; j++) {
            int v = find_node_idx(topo, topo->adjacency[u][j].dst_id);
            if (v >= 0 && !visited[v]) { visited[v] = 1; queue[qt++] = v; }
        }
    }
    int connected = 1;
    for (int i = 0; i < topo->node_count; i++)
        if (!visited[i]) { connected = 0; break; }
    free(visited); free(queue);
    return connected;
}

int topology_detect_partition(EdgeCloudTopology *topo) {
    if (!topo) return 0;
    return topology_is_connected(topo) ? 0 : 1;
}

void topology_print(const EdgeCloudTopology *topo) {
    if (!topo) return;
    printf("Topology[%s]: %d nodes, condition=%s\n",
           topo->id, topo->node_count, net_condition_string(topo->condition));
    for (int i = 0; i < topo->node_count; i++) {
        printf("  Node[%s]: %d links\n", topo->node_ids[i], topo->adj_count[i]);
        for (int j = 0; j < topo->adj_count[i]; j++)
            printf("    -> %s bw=%.1fMbps lat=%.0fus\n",
                   topo->adjacency[i][j].dst_id,
                   topo->adjacency[i][j].bandwidth_mbps,
                   topo->adjacency[i][j].latency_us);
    }
}


int topology_import_from_matrix(EdgeCloudTopology *topo,
                                 const double *latency_matrix, int n) {
    if (!topo || !latency_matrix || n <= 0 || n > CCS_MAX_TOPOLOGY_NODES) return -1;
    for (int i = 0; i < n; i++) {
        char nid[64]; snprintf(nid, sizeof(nid), "N%d", i);
        topology_add_node(topo, nid);
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (i != j && latency_matrix[i*n+j] > 0)
                topology_add_link(topo, topo->node_ids[i], topo->node_ids[j],
                                   100.0, latency_matrix[i*n+j]);
    return 0;
}
/* -------- PacketBuffer (Ring Buffer) -------- */

PacketBuffer* pktbuf_create(const char *id, int max_size) {
    if (max_size <= 0 || max_size > CCS_MAX_TASK_QUEUE) return NULL;
    PacketBuffer *buf = (PacketBuffer*)calloc(1, sizeof(PacketBuffer));
    if (!buf) return NULL;
    if (id) strncpy(buf->id, id, sizeof(buf->id)-1);
    buf->max_size = max_size;
    buf->head = 0; buf->tail = 0; buf->count = 0;
    return buf;
}

void pktbuf_free(PacketBuffer *buf) { free(buf); }

int pktbuf_enqueue(PacketBuffer *buf, const ControlPacket *pkt) {
    if (!buf || !pkt || pktbuf_is_full(buf)) return -1;
    buf->buffer[buf->tail] = *pkt;
    buf->tail = (buf->tail + 1) % buf->max_size;
    buf->count++;
    buf->packets_enqueued++;
    return 0;
}

int pktbuf_dequeue(PacketBuffer *buf, ControlPacket *pkt) {
    if (!buf || !pkt || pktbuf_is_empty(buf)) return -1;
    *pkt = buf->buffer[buf->head];
    buf->head = (buf->head + 1) % buf->max_size;
    buf->count--;
    buf->packets_dequeued++;
    return 0;
}

int pktbuf_peek(const PacketBuffer *buf, ControlPacket *pkt) {
    if (!buf || !pkt || pktbuf_is_empty(buf)) return -1;
    *pkt = buf->buffer[buf->head]; return 0;
}

int pktbuf_is_full(const PacketBuffer *buf) { return buf ? (buf->count >= buf->max_size) : 1; }
int pktbuf_is_empty(const PacketBuffer *buf) { return buf ? (buf->count == 0) : 1; }
int pktbuf_count(const PacketBuffer *buf) { return buf ? buf->count : -1; }

void pktbuf_flush(PacketBuffer *buf) {
    if (!buf) return;
    buf->head = 0; buf->tail = 0; buf->count = 0;
}

/* -------- NetworkLink -------- */

int netlink_update_metrics(NetworkLink *link, double latency, double jitter,
                            double loss_rate, double utilization) {
    if (!link) return -1;
    link->latency_us = 0.7 * link->latency_us + 0.3 * latency;
    link->jitter_us = 0.7 * link->jitter_us + 0.3 * jitter;
    link->packet_loss_rate = 0.7 * link->packet_loss_rate + 0.3 * loss_rate;
    link->utilization = utilization;
    return 0;
}

int netlink_is_congested(const NetworkLink *link) {
    return (link && link->utilization > 0.85) ? 1 : 0;
}

double netlink_effective_bandwidth(const NetworkLink *link) {
    if (!link) return 0.0;
    return link->bandwidth_mbps * (1.0 - link->packet_loss_rate);
}

double netlink_bandwidth_delay_product(const NetworkLink *link) {
    if (!link) return 0.0;
    return link->bandwidth_mbps * 1e6 / 8.0 * (link->latency_us / 1e6);
}

void netlink_print(const NetworkLink *link) {
    if (!link) return;
    printf("Link[%s->%s]: bw=%.1fMbps lat=%.0fus jit=%.0fus loss=%.4f util=%.2f\n",
           link->src_id, link->dst_id, link->bandwidth_mbps, link->latency_us,
           link->jitter_us, link->packet_loss_rate, link->utilization);
}


/* -------- BandwidthAllocator -------- */

BandwidthAllocator* bw_alloc_create(const char *link_id, double total_bw) {
    BandwidthAllocator *ba = (BandwidthAllocator*)calloc(1, sizeof(BandwidthAllocator));
    if (!ba) return NULL;
    if (link_id) strncpy(ba->link_id, link_id, sizeof(ba->link_id)-1);
    ba->total_bandwidth_mbps = total_bw;
    ba->available_bandwidth_mbps = total_bw;
    return ba;
}

void bw_alloc_free(BandwidthAllocator *ba) { free(ba); }

int bw_alloc_reserve(BandwidthAllocator *ba, const char *consumer_id,
                      double amount, QoSLevel qos) {
    if (!ba || !consumer_id || amount <= 0) return -1;
    if (ba->alloc_count >= CCS_MAX_ALLOCATIONS) return -1;
    if (amount > ba->available_bandwidth_mbps) return -1;
    int idx = ba->alloc_count;
    ba->alloc_count++;
    strncpy(ba->alloc_ids[idx], consumer_id, 63);
    ba->allocation[idx] = amount;
    ba->allocated_bandwidth_mbps += amount;
    ba->available_bandwidth_mbps -= amount;
    if ((int)qos < (int)ba->min_qos) ba->min_qos = qos;
    return 0;
}

int bw_alloc_release(BandwidthAllocator *ba, const char *consumer_id) {
    if (!ba || !consumer_id) return -1;
    for (int i = 0; i < ba->alloc_count; i++) {
        if (strcmp(ba->alloc_ids[i], consumer_id) == 0) {
            ba->allocated_bandwidth_mbps -= ba->allocation[i];
            ba->available_bandwidth_mbps += ba->allocation[i];
            for (int j = i; j < ba->alloc_count - 1; j++) {
                ba->allocation[j] = ba->allocation[j+1];
                strncpy(ba->alloc_ids[j], ba->alloc_ids[j+1], 63);
            }
            ba->alloc_count--;
            return 0;
        }
    }
    return -1;
}

double bw_alloc_get_available(const BandwidthAllocator *ba) {
    return ba ? ba->available_bandwidth_mbps : 0.0;
}

int bw_alloc_is_feasible(const BandwidthAllocator *ba, double amount) {
    return (ba && amount <= ba->available_bandwidth_mbps) ? 1 : 0;
}
/* -------- Edge-Cloud Communication -------- */

int edge_to_cloud_send(CloudControlSystem *ccs, const double *state,
                        int dim, double timestamp) {
    (void)timestamp;
    if (!ccs || !state || dim <= 0) return -1;
    int n = ccs->n;
    if (dim > n) dim = n;
    for (int i = 0; i < dim; i++) ccs->plant_state.x[i] = state[i];
    ccs->edge_node.cloud_connected = 1;
    ccs->edge_node.cloud_cycles++;
    ccs->packets_sent++;
    return 0;
}

int cloud_to_edge_send(CloudControlSystem *ccs, const double *control,
                        int dim, double timestamp) {
    (void)timestamp;
    if (!ccs || !control || dim <= 0) return -1;
    int m = ccs->m;
    if (dim > m) dim = m;
    for (int i = 0; i < dim; i++) ccs->plant_state.u[i] = control[i];
    ccs->edge_node.local_cycles++;
    return 0;
}

int edge_holdover_check(const EdgeNode *edge) {
    if (!edge) return 0;
    return (edge->holdover_remaining_us > 0) ? 1 : 0;
}

int edge_holdover_enter(EdgeNode *edge) {
    if (!edge) return -1;
    edge->holdover_remaining_us = edge->max_holdover_time_us;
    edge->cloud_connected = 0;
    return 0;
}

int edge_holdover_exit(EdgeNode *edge) {
    if (!edge) return -1;
    edge->holdover_remaining_us = 0;
    edge->cloud_connected = 1;
    return 0;
}

double edge_holdover_remaining(const EdgeNode *edge) {
    return edge ? edge->holdover_remaining_us : 0.0;
}

/* -------- Network Condition Simulation -------- */

double net_condition_delay_factor(NetworkCondition cond) {
    switch (cond) {
    case NET_CONDITION_NORMAL:      return 1.0;
    case NET_CONDITION_CONGESTED:   return 3.0;
    case NET_CONDITION_HIGH_LOSS:   return 1.5;
    case NET_CONDITION_PARTITIONED: return 1e9;
    case NET_CONDITION_DDOS:        return 10.0;
    case NET_CONDITION_DEGRADED:    return 5.0;
    default: return 1.0;
    }
}

double net_condition_loss_factor(NetworkCondition cond) {
    switch (cond) {
    case NET_CONDITION_NORMAL:      return 1.0;
    case NET_CONDITION_CONGESTED:   return 5.0;
    case NET_CONDITION_HIGH_LOSS:   return 20.0;
    case NET_CONDITION_DDOS:        return 50.0;
    case NET_CONDITION_PARTITIONED: return 1e6;
    case NET_CONDITION_DEGRADED:    return 10.0;
    default: return 1.0;
    }
}

const char* net_condition_string(NetworkCondition cond) {
    switch (cond) {
    case NET_CONDITION_NORMAL:      return "NORMAL";
    case NET_CONDITION_CONGESTED:   return "CONGESTED";
    case NET_CONDITION_HIGH_LOSS:   return "HIGH_LOSS";
    case NET_CONDITION_PARTITIONED: return "PARTITIONED";
    case NET_CONDITION_DDOS:        return "DDOS";
    case NET_CONDITION_DEGRADED:    return "DEGRADED";
    default: return "UNKNOWN";
    }
}

NetworkCondition net_detect_condition(double loss_rate, double jitter,
                                       double rtt) {
    if (loss_rate > 0.5) return NET_CONDITION_DDOS;
    if (loss_rate > 0.1) return NET_CONDITION_HIGH_LOSS;
    if (jitter > rtt) return NET_CONDITION_DEGRADED;
    if (loss_rate > 0.05) return NET_CONDITION_CONGESTED;
    return NET_CONDITION_NORMAL;
}

/* -------- Multi-Path Routing -------- */

int mproute_add_link(MultiPathRoute *route, const NetworkLink *link) {
    if (!route || !link || route->link_count >= CCS_MAX_NETWORK_PATHS) return -1;
    route->links[route->link_count++] = *link;
    route->total_latency_us += link->latency_us;
    if (link->bandwidth_mbps < route->total_bandwidth_mbps || route->total_bandwidth_mbps == 0)
        route->total_bandwidth_mbps = link->bandwidth_mbps;
    route->reliability *= (1.0 - link->packet_loss_rate);
    route->hop_count++;
    return 0;
}

double mproute_total_latency(const MultiPathRoute *route) {
    return route ? route->total_latency_us : -1.0;
}

double mproute_effective_bw(const MultiPathRoute *route) {
    return route ? route->total_bandwidth_mbps : 0.0;
}

int mproute_compare(const MultiPathRoute *a, const MultiPathRoute *b) {
    if (!a || !b) return 0;
    /* Prefer lower latency, then higher bandwidth */
    if (a->total_latency_us < b->total_latency_us) return 1;
    if (a->total_latency_us > b->total_latency_us) return -1;
    if (a->total_bandwidth_mbps > b->total_bandwidth_mbps) return 1;
    if (a->total_bandwidth_mbps < b->total_bandwidth_mbps) return -1;
    return 0;
}

/* -------- Cloud Region Latency Model -------- */

double cloud_latency_estimate(const char *region_from, const char *region_to) {
    if (!region_from || !region_to) return 100000.0;
    /* Simplified cloud region latency model based on geographic distance approximation */
    /* In practice, this would use a latency map (e.g., cloudping.co data) */
    if (strcmp(region_from, region_to) == 0) return 1000.0; /* Same region */
    /* Cross-region latencies (approximate RTTs in us) */
    struct { const char *a; const char *b; double rtt; } pairs[] = {
        {"us-east", "us-west", 60000}, {"us-east", "eu-west", 80000},
        {"us-west", "eu-west", 130000}, {"eu-west", "ap-southeast", 180000},
        {"us-east", "ap-northeast", 150000}, {"us-west", "ap-northeast", 110000}
    };
    int n_pairs = (int)(sizeof(pairs)/sizeof(pairs[0]));
    for (int i = 0; i < n_pairs; i++) {
        if ((strstr(region_from, pairs[i].a) && strstr(region_to, pairs[i].b)) ||
            (strstr(region_from, pairs[i].b) && strstr(region_to, pairs[i].a)))
            return pairs[i].rtt;
    }
    return 100000.0; /* Default inter-continental estimate */
}

double edge_to_cloud_rtt_model(double edge_distance_km,
                                 const char *cloud_region) {
    /* Speed of light in fiber: ~200 km/ms, plus switching/processing overhead */
    (void)cloud_region;
    double propagation_us = edge_distance_km / 200.0 * 1000.0;
    double overhead_us = 5000.0; /* 5ms overhead for processing and queueing */
    return propagation_us * 2.0 + overhead_us;
}

double cloud_inter_region_rtt(const char *region_a, const char *region_b) {
    return cloud_latency_estimate(region_a, region_b);
}
/* -------- Packet Serialization (Simple Binary Protocol) -------- */

int pkt_serialize(const ControlPacket *pkt, uint8_t *buffer, int buf_size) {
    if (!pkt || !buffer || buf_size < 128) return -1;
    int offset = 0;
    memcpy(buffer + offset, &pkt->header.seq_num, 8); offset += 8;
    memcpy(buffer + offset, &pkt->header.timestamp_us, 8); offset += 8;
    memcpy(buffer + offset, pkt->header.src_id, 64); offset += 64;
    memcpy(buffer + offset, pkt->header.dst_id, 64); offset += 64;
    memcpy(buffer + offset, &pkt->header.payload_type, 4); offset += 4;
    memcpy(buffer + offset, &pkt->header.payload_length, 4); offset += 4;
    memcpy(buffer + offset, &pkt->header.priority, 4); offset += 4;
    memcpy(buffer + offset, &pkt->header.ttl, 4); offset += 4;
    memcpy(buffer + offset, &pkt->header.sent_at, 8); offset += 8;
    memcpy(buffer + offset, &pkt->header.received_at, 8); offset += 8;
    memcpy(buffer + offset, &pkt->is_control, 4); offset += 4;
    memcpy(buffer + offset, &pkt->value_dim, 4); offset += 4;
    int val_bytes = pkt->value_dim * (int)sizeof(double);
    if (offset + val_bytes > buf_size) return -1;
    memcpy(buffer + offset, pkt->value, (size_t)val_bytes); offset += val_bytes;
    return offset;
}

int pkt_deserialize(const uint8_t *buffer, int buf_size, ControlPacket *pkt) {
    if (!buffer || buf_size < 128 || !pkt) return -1;
    int offset = 0;
    memcpy(&pkt->header.seq_num, buffer + offset, 8); offset += 8;
    memcpy(&pkt->header.timestamp_us, buffer + offset, 8); offset += 8;
    memcpy(pkt->header.src_id, buffer + offset, 64); offset += 64;
    memcpy(pkt->header.dst_id, buffer + offset, 64); offset += 64;
    memcpy(&pkt->header.payload_type, buffer + offset, 4); offset += 4;
    memcpy(&pkt->header.payload_length, buffer + offset, 4); offset += 4;
    memcpy(&pkt->header.priority, buffer + offset, 4); offset += 4;
    memcpy(&pkt->header.ttl, buffer + offset, 4); offset += 4;
    memcpy(&pkt->header.sent_at, buffer + offset, 8); offset += 8;
    memcpy(&pkt->header.received_at, buffer + offset, 8); offset += 8;
    memcpy(&pkt->is_control, buffer + offset, 4); offset += 4;
    memcpy(&pkt->value_dim, buffer + offset, 4); offset += 4;
    int val_bytes = pkt->value_dim * (int)sizeof(double);
    if (offset + val_bytes > buf_size) return -1;
    memcpy(pkt->value, buffer + offset, (size_t)val_bytes); offset += val_bytes;
    return offset;
}

int pkt_create_measurement(ControlPacket *pkt, const char *src,
                            const double *y, int dim, uint64_t seq) {
    if (!pkt || !src || !y || dim <= 0) return -1;
    memset(pkt, 0, sizeof(*pkt));
    pkt->header.seq_num = seq;
    pkt->header.timestamp_us = (uint64_t)(time(NULL) * 1000000);
    strncpy(pkt->header.src_id, src, 63);
    pkt->header.payload_type = 1;
    pkt->header.priority = 0;
    pkt->header.ttl = 16;
    pkt->is_control = 0;
    pkt->value_dim = dim;
    for (int i = 0; i < dim; i++) pkt->value[i] = y[i];
    return 0;
}

int pkt_create_control(ControlPacket *pkt, const char *src,
                        const double *u, int dim, uint64_t seq) {
    if (!pkt || !src || !u || dim <= 0) return -1;
    memset(pkt, 0, sizeof(*pkt));
    pkt->header.seq_num = seq;
    pkt->header.timestamp_us = (uint64_t)(time(NULL) * 1000000);
    strncpy(pkt->header.src_id, src, 63);
    pkt->header.payload_type = 2;
    pkt->header.priority = 1;
    pkt->header.ttl = 16;
    pkt->is_control = 1;
    pkt->value_dim = dim;
    for (int i = 0; i < dim; i++) pkt->value[i] = u[i];
    return 0;
}
