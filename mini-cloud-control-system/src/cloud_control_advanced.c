/* cloud_control_advanced.c - Advanced Cloud Control
 *
 * Time Synchronization, Age of Information, Packet Reordering,
 * Token Bucket, Congestion Detection, WRR Scheduling, TSN Guard Band.
 *
 * Knowledge Coverage:
 *   L7: NTP clock sync, AoI tracking, WRR scheduling, OWD estimation
 *   L8: Reorder buffer, Token bucket, Congestion detection
 *   L9: TSN guard band (IEEE 802.1Qbv)
 *
 * References:
 *   - Mills, "NTP", RFC 5905
 *   - Kaul et al., "Age of Information", IEEE INFOCOM 2012
 *   - Bennett et al., "Packet Reordering", IEEE/ACM ToN 1999
 *   - Parekh & Gallager, "GPS Flow Control", IEEE/ACM ToN 1993
 *   - Jacobson, "Congestion Avoidance and Control", SIGCOMM 1988
 *   - Shreedhar & Varghese, "Deficit Round-Robin", IEEE/ACM ToN 1996
 *   - IEEE 802.1Qbv-2015, "Scheduled Traffic"
 */

#include "cloud_control_core.h"
#include "cloud_control_network.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ===========================================================================
 * Section 1: NTP Clock Synchronization (L7)
 *
 * NTP four-timestamp exchange (RFC 5905):
 *   Client                Server
 *     |---- t1 (local) ---->|
 *     |                     | (t2 remote)
 *     |                     |
 *     |<--- t4 (local) -----|
 *     |                     (t3 remote)
 *
 *   offset = ((t2 - t1) - (t4 - t3)) / 2
 *   delay  = (t4 - t1) - (t3 - t2)
 *
 * Minimum-delay principle (Mills): the sample with smallest RTT
 * gives the best offset estimate because queueing delay adds
 * positive bias to both directions.
 * =========================================================================== */

#define NTP_MAX_SAMPLES  64
#define ACP_EWMA_ALPHA   0.1
#define REORDER_WMAX     256
#define REORDER_DEF_TO   0.5
#define WRR_MAX_PATHS    8
#define DCD_WINDOW       8

/* --- NTP Data Structures --- */
typedef struct {
    double t1, t2, t3, t4;
    double offset_us;
    double delay_us;
} NtpSample;

typedef struct {
    NtpSample samples[NTP_MAX_SAMPLES];
    int       count;
    int       widx;
    double    best_offset_us;
    double    best_delay_us;
    double    skew_ppb;
    double    uncertainty_us;
    double    last_sync;
    int       locked;
    int       rounds;
} NtpClock;

static void ntp_init(NtpClock *nc) {
    if (!nc) return;
    memset(nc, 0, sizeof(*nc));
    nc->uncertainty_us = 1e6;
}

static int ntp_process(NtpClock *nc, double t1, double t2,
                        double t3, double t4) {
    if (!nc || t4 <= t1 || t3 < t2) return -1;
    double o = ((t2 - t1) - (t4 - t3)) / 2.0;
    double d = (t4 - t1) - (t3 - t2);
    if (d < 0.0) d = 0.0;

    NtpSample *s = &nc->samples[nc->widx];
    s->t1 = t1; s->t2 = t2; s->t3 = t3; s->t4 = t4;
    s->offset_us = o * 1e6;
    s->delay_us  = d * 1e6;
    nc->widx = (nc->widx + 1) % NTP_MAX_SAMPLES;
    if (nc->count < NTP_MAX_SAMPLES) nc->count++;

    /* Minimum-delay best offset selection */
    double md = 1e300;
    for (int i = 0; i < nc->count; i++) {
        if (nc->samples[i].delay_us < md) {
            md = nc->samples[i].delay_us;
            nc->best_offset_us = nc->samples[i].offset_us;
            nc->best_delay_us  = nc->samples[i].delay_us;
        }
    }

    /* Uncertainty: sliding std dev of last 8 */
    if (nc->count >= 3) {
        int w = (nc->count < 8) ? nc->count : 8;
        double su = 0.0, sq = 0.0;
        for (int i = 0; i < w; i++) {
            int idx = (nc->widx - 1 - i + NTP_MAX_SAMPLES) % NTP_MAX_SAMPLES;
            double off = nc->samples[idx].offset_us;
            su += off; sq += off * off;
        }
        double m = su / w;
        double v = sq / w - m * m;
        if (v > 0) nc->uncertainty_us = sqrt(v);
    }

    nc->last_sync = t4;
    nc->rounds++;
    nc->locked = (nc->count >= 4 && nc->uncertainty_us < 500.0);
    return 0;
}

/* L5: Linear regression skew estimation */
static double ntp_skew(const NtpClock *nc) {
    if (!nc || nc->count < 4) return 0.0;
    int n = nc->count;
    double st = 0, so = 0, stt = 0, sto = 0;
    double t0 = nc->samples[0].t1;
    for (int i = 0; i < n; i++) {
        double dt = nc->samples[i].t1 - t0;
        double off = nc->samples[i].offset_us;
        st  += dt;
        so  += off;
        stt += dt * dt;
        sto += dt * off;
    }
    double den = n * stt - st * st;
    if (fabs(den) < 1e-15) return 0.0;
    return (n * sto - st * so) / den * 1e3;  /* us/s -> ppb */
}

/* L3: Timestamp compensation with skew */
static double ntp_compensate(const NtpClock *nc, double remote) {
    if (!nc || !nc->locked) return remote;
    double dt = remote - nc->last_sync;
    return remote + nc->best_offset_us * 1e-6 + nc->skew_ppb * 1e-9 * dt;
}

/* L3: Holdover time estimation */
static double ntp_holdover(const NtpClock *nc, double max_off_us) {
    if (!nc || !nc->locked) return 0.0;
    double skew_abs = fabs(nc->skew_ppb) * 1e-9;
    if (skew_abs < 1e-15) return 1e100;
    double remaining = max_off_us - fabs(nc->best_offset_us);
    if (remaining <= 0) return 0.0;
    return remaining * 1e-6 / skew_abs;
}

/* ===========================================================================
 * Section 2: Age of Information (L7)
 *
 * AoI(t) = t - u(t) where u(t) = generation time of latest delivered update.
 *
 * Theorem (Kaul et al., M/D/1 AoI):
 *   E[AoI] = 1/mu + 1/lambda + lambda/(2*mu*(mu-lambda))
 *
 * Peak AoI = max instantaneous AoI across all update cycles.
 * =========================================================================== */

typedef struct {
    double last_gen;
    double last_del;
    double current_aoi;
    double peak_aoi;
    double avg_aoi;
    double min_aoi;
    uint64_t gens;
    uint64_t dels;
    uint64_t violations;
    double vthreshold;
    double alpha;
} AoiState;

static void aoi_init(AoiState *as, double thr) {
    if (!as) return;
    memset(as, 0, sizeof(*as));
    as->vthreshold = (thr > 0) ? thr : 1.0;
    as->alpha = ACP_EWMA_ALPHA;
    as->min_aoi = 1e300;
}

static void aoi_generate(AoiState *as, double t) {
    if (!as) return;
    as->last_gen = t;
    as->gens++;
}

static double aoi_deliver(AoiState *as, double gen_t, double del_t) {
    if (!as) return -1.0;
    double ia = del_t - gen_t;
    if (ia < 0) ia = 0;
    as->current_aoi = ia;
    as->last_del = del_t;
    as->dels++;
    if (ia > as->peak_aoi) as->peak_aoi = ia;
    if (ia < as->min_aoi) as->min_aoi = ia;
    if (as->dels == 1) as->avg_aoi = ia;
    else as->avg_aoi = as->alpha * ia + (1.0 - as->alpha) * as->avg_aoi;
    if (ia > as->vthreshold) as->violations++;
    return ia;
}

static double aoi_now(const AoiState *as, double now) {
    if (!as || as->dels == 0) return now;
    return now - as->last_gen;
}

static double aoi_delratio(const AoiState *as) {
    if (!as || as->gens == 0) return 0.0;
    return (double)as->dels / (double)as->gens;
}

/* Exported: M/D/1 theoretical average AoI */
double net_aoi_theoretical_avg(double lambda, double mu) {
    if (lambda <= 0 || mu <= 0 || lambda >= mu) return -1.0;
    return 1.0/mu + 1.0/lambda + lambda/(2.0 * mu * (mu - lambda));
}

/* L5: Golden-section search for optimal update rate */
static double aoi_optimal_rate(double mu) {
    if (mu <= 0) return 0.0;
    double lo = 0.001 * mu, hi = 0.999 * mu;
    double phi = (sqrt(5.0) - 1.0) / 2.0;
    for (int i = 0; i < 50; i++) {
        double c = hi - phi * (hi - lo);
        double d = lo + phi * (hi - lo);
        if (net_aoi_theoretical_avg(c, mu) < net_aoi_theoretical_avg(d, mu))
            hi = d;
        else
            lo = c;
    }
    return (lo + hi) / 2.0;
}

/* ===========================================================================
 * Section 3: Packet Reordering Buffer (L8)
 *
 * Sliding-window reorder buffer with bounded delay.
 *   - seq < expected: LATE (discard)
 *   - seq == expected: IN-ORDER (deliver, flush consecutive)
 *   - seq > expected: EARLY (buffer)
 *   - seq > expected + W - 1: DROP (window overflow)
 *
 * Reference: Bennett et al., IEEE/ACM ToN, 1999
 * =========================================================================== */

typedef struct {
    int      occ;
    uint64_t seq;
    double   data[CCS_MAX_STATES];
    int      dim;
    double   arr_t;
    double   gen_t;
} RSlot;

typedef struct {
    RSlot    w[REORDER_WMAX];
    int      wsz;
    uint64_t next;
    uint64_t max_rx;
    double   mwait;
    double   gap_start;
    uint64_t total_rx;
    uint64_t total_inorder;
    uint64_t total_late;
    uint64_t total_dups;
    uint64_t total_drops;
    uint64_t total_timeouts;
    double   sum_wait;
    int      num_delayed;
} RBuffer;

static void rb_init(RBuffer *rb, int wsz, double mw) {
    if (!rb) return;
    memset(rb, 0, sizeof(*rb));
    rb->wsz = (wsz > 0 && wsz <= REORDER_WMAX) ? wsz : REORDER_WMAX;
    rb->mwait = (mw > 0) ? mw : REORDER_DEF_TO;
    rb->next = 1;
}

static int rb_insert(RBuffer *rb, uint64_t seq, const double *data,
                      int dim, double arr_t, double gen_t) {
    if (!rb || !data || dim <= 0) return -1;
    rb->total_rx++;

    if (seq < rb->next) { rb->total_late++; return -1; }  /* Late */

    if (seq == rb->next) {  /* In-order */
        rb->next++;
        rb->total_inorder++;
        rb->gap_start = 0;
        /* Flush buffered consecutive packets */
        for (;;) {
            uint64_t n = rb->next;
            int s = (int)(n % (uint64_t)rb->wsz);
            if (!rb->w[s].occ || rb->w[s].seq != n) break;
            rb->w[s].occ = 0;
            rb->sum_wait += arr_t - rb->w[s].arr_t;
            rb->num_delayed++;
            rb->next++;
            rb->total_inorder++;
        }
        return 1;
    }

    /* Early arrival - buffer */
    if (seq > rb->next + (uint64_t)(rb->wsz - 1)) { rb->total_drops++; return -1; }
    int s = (int)(seq % (uint64_t)rb->wsz);
    if (rb->w[s].occ && rb->w[s].seq == seq) { rb->total_dups++; return -1; }
    rb->w[s].occ = 1;
    rb->w[s].seq = seq;
    rb->w[s].dim = dim;
    rb->w[s].arr_t = arr_t;
    rb->w[s].gen_t = gen_t;
    for (int i = 0; i < dim && i < CCS_MAX_STATES; i++)
        rb->w[s].data[i] = data[i];
    if (seq > rb->max_rx) rb->max_rx = seq;
    if (rb->gap_start == 0) rb->gap_start = arr_t;
    return 0;
}

/* L4: Bounded-delay timeout - skip missing after max_wait */
static int rb_timeout(RBuffer *rb, double now) {
    if (!rb || rb->gap_start == 0) return 0;
    if (now - rb->gap_start < rb->mwait) return 0;
    int skipped = 0;
    while (rb->next <= rb->max_rx) {
        int s = (int)(rb->next % (uint64_t)rb->wsz);
        if (rb->w[s].occ && rb->w[s].seq == rb->next) break;
        rb->next++;
        skipped++;
        rb->total_timeouts++;
    }
    rb->gap_start = 0;
    return skipped;
}

static double rb_density(const RBuffer *rb) {
    if (!rb || rb->total_rx == 0) return 0.0;
    return 1.0 - (double)rb->total_inorder / (double)rb->total_rx;
}

static double rb_avgwait(const RBuffer *rb) {
    if (!rb || rb->num_delayed == 0) return 0.0;
    return rb->sum_wait / (double)rb->num_delayed;
}

/* ===========================================================================
 * Section 4: Token Bucket Traffic Shaper (L8)
 *
 * Enforces long-term rate r with burst tolerance b.
 * Tokens accumulate at rate r (tokens/sec), capped at b.
 * Each packet consumes size tokens. Non-conformant if insufficient.
 *
 * Reference: Parekh & Gallager, IEEE/ACM ToN, 1993
 * =========================================================================== */

typedef struct {
    double   rate;
    double   burst;
    double   tokens;
    double   last_upd;
    uint64_t passed;
    uint64_t dropped;
    uint64_t bpassed;
    uint64_t bdropped;
} TBucket;

static void tb_init(TBucket *tb, double rate, double burst) {
    if (!tb) return;
    memset(tb, 0, sizeof(*tb));
    tb->rate = rate;
    tb->burst = burst;
    tb->tokens = burst;
}

static int tb_consume(TBucket *tb, double sz, double now) {
    if (!tb) return 0;
    if (tb->last_upd > 0 && now > tb->last_upd) {
        double elap = now - tb->last_upd;
        tb->tokens += tb->rate * elap;
        if (tb->tokens > tb->burst) tb->tokens = tb->burst;
    }
    tb->last_upd = now;
    if (tb->tokens >= sz) {
        tb->tokens -= sz;
        tb->passed++;
        tb->bpassed += (uint64_t)sz;
        return 1;
    }
    tb->dropped++;
    tb->bdropped += (uint64_t)sz;
    return 0;
}

static double tb_rate(const TBucket *tb) {
    if (!tb) return 0.0;
    uint64_t tot = tb->passed + tb->dropped;
    return tot > 0 ? (double)tb->passed / (double)tot : 1.0;
}

/* ===========================================================================
 * Section 5: Delay-Gradient Congestion Detection (L8)
 *
 * Congestion signal when RTT increases for N consecutive samples.
 * Jacobson's principle: monotonic RTT growth indicates queue buildup.
 *
 * Reference: Jacobson, SIGCOMM 1988
 * =========================================================================== */

typedef struct {
    double rtt[DCD_WINDOW];
    int    idx;
    int    cnt;
    int    consec;
    int    thresh;
    int    congested;
    double severity;
} DCongestion;

static void dc_init(DCongestion *dc, int thr) {
    if (!dc) return;
    memset(dc, 0, sizeof(*dc));
    dc->thresh = (thr > 0) ? thr : 4;
}

static int dc_update(DCongestion *dc, double rtt_us) {
    if (!dc) return 0;
    dc->rtt[dc->idx] = rtt_us;
    dc->idx = (dc->idx + 1) % DCD_WINDOW;
    if (dc->cnt < DCD_WINDOW) dc->cnt++;

    if (dc->cnt >= 2) {
        int p = (dc->idx - 2 + DCD_WINDOW) % DCD_WINDOW;
        int c = (dc->idx - 1 + DCD_WINDOW) % DCD_WINDOW;
        if (dc->rtt[c] > dc->rtt[p]) dc->consec++;
        else dc->consec = 0;
    }
    int was = dc->congested;
    dc->congested = (dc->consec >= dc->thresh);

    if (dc->congested && dc->cnt >= 2) {
        int l = (dc->idx - 1 + DCD_WINDOW) % DCD_WINDOW;
        int f = (dc->idx - dc->thresh + DCD_WINDOW) % DCD_WINDOW;
        double ratio = dc->rtt[l] / (dc->rtt[f] + 1e-6);
        dc->severity = (ratio - 1.0) / 2.0;
        if (dc->severity > 1.0) dc->severity = 1.0;
        if (dc->severity < 0.0) dc->severity = 0.0;
    } else {
        dc->severity *= 0.9;
    }
    return (was != dc->congested) ? 1 : 0;
}

/* ===========================================================================
 * Section 6: Deficit Round-Robin Multi-Path Scheduler (L7)
 *
 * WRR distributes packets across paths proportional to bandwidth.
 * DRR handles variable packet sizes via deficit counters.
 *
 * Reference: Shreedhar & Varghese, IEEE/ACM ToN, 1996
 * =========================================================================== */

typedef struct {
    int    id;
    double bw;
    double weight;
    int    deficit;
    int    quantum;
    uint64_t sent;
} WPath;

typedef struct {
    WPath paths[WRR_MAX_PATHS];
    int   n;
    int   cur;
} WRR;

static void wrr_init(WRR *w) { if (w) memset(w, 0, sizeof(*w)); }

static int wrr_add(WRR *w, int id, double bw) {
    if (!w || w->n >= WRR_MAX_PATHS || bw <= 0) return -1;
    int idx = w->n++;
    w->paths[idx].id = id;
    w->paths[idx].bw = bw;
    double tbw = 0;
    for (int i = 0; i < w->n; i++) tbw += w->paths[i].bw;
    int bq = 1500;
    for (int i = 0; i < w->n; i++) {
        w->paths[i].weight = w->paths[i].bw / tbw;
        w->paths[i].quantum = (int)(bq * w->paths[i].weight);
        if (w->paths[i].quantum < 1) w->paths[i].quantum = 1;
    }
    return 0;
}

static int wrr_sched(WRR *w, int pkt_sz) {
    if (!w || w->n == 0 || pkt_sz <= 0) return -1;
    int start = w->cur;
    do {
        w->paths[w->cur].deficit += w->paths[w->cur].quantum;
        if (w->paths[w->cur].deficit >= pkt_sz) {
            w->paths[w->cur].deficit -= pkt_sz;
            w->paths[w->cur].sent++;
            int sel = w->paths[w->cur].id;
            w->cur = (w->cur + 1) % w->n;
            return sel;
        }
        w->cur = (w->cur + 1) % w->n;
    } while (w->cur != start);
    return -1;
}

/* ===========================================================================
 * Section 7: One-Way Delay & TSN Guard Band (L7/L9)
 *
 * OWD via min-RTT/2: tightest symmetric lower bound without sync clocks.
 * TSN guard band (IEEE 802.1Qbv): protection before scheduled window.
 *
 * Reference: Choi & Yoo, Computer Comm. 2005; IEEE 802.1Qbv-2015
 * =========================================================================== */

static double owd_minrtt(const double *rtt, int n) {
    if (!rtt || n <= 0) return 0.0;
    double mr = rtt[0];
    for (int i = 1; i < n; i++) if (rtt[i] < mr) mr = rtt[i];
    return mr / 2.0;
}

static double owd_asym(const double *rtt, int n) {
    if (!rtt || n < 2) return 0.0;
    double s = 0, sq = 0;
    for (int i = 0; i < n; i++) { s += rtt[i]; sq += rtt[i] * rtt[i]; }
    double m = s / n;
    double v = sq / n - m * m;
    if (v < 0 || m < 1e-9) return 0.0;
    return sqrt(v) / m;
}

/* L4: TSN guard band = tx_time + 2*sync_uncertainty */
double net_tsn_guard_band(int max_frame_bytes, double link_rate_bps,
                           double sync_unc_us) {
    if (max_frame_bytes <= 0 || link_rate_bps <= 0) return -1.0;
    double bits = (double)(max_frame_bytes + 20) * 8.0;  /* +preamble+IFG */
    double tx_us = bits / link_rate_bps * 1e6;
    return tx_us + 2.0 * sync_unc_us;
}

/* ===========================================================================
 * PUBLIC API
 * =========================================================================== */

/* --- AoI --- */
void* ccs_aoi_create(double thr) {
    AoiState *a = (AoiState*)calloc(1, sizeof(AoiState));
    if (a) aoi_init(a, thr);
    return a;
}
void ccs_aoi_free(void *a) { free(a); }
void ccs_aoi_generate(void *a, double t) { aoi_generate((AoiState*)a, t); }
double ccs_aoi_deliver(void *a, double gt, double dt) {
    return aoi_deliver((AoiState*)a, gt, dt);
}
double ccs_aoi_current(void *a, double now) {
    return aoi_now((AoiState*)a, now);
}
double ccs_aoi_peak(void *a) {
    return a ? ((AoiState*)a)->peak_aoi : 0.0;
}
double ccs_aoi_average(void *a) {
    return a ? ((AoiState*)a)->avg_aoi : 0.0;
}
double ccs_aoi_delivery_ratio(void *a) {
    return aoi_delratio((AoiState*)a);
}
int ccs_aoi_is_violated(void *a, double now) {
    AoiState *as = (AoiState*)a;
    return (as && aoi_now(as, now) > as->vthreshold) ? 1 : 0;
}
void ccs_aoi_reset(void *a) { aoi_init((AoiState*)a, 1.0); }
double ccs_aoi_optimal_rate(double mu) { return aoi_optimal_rate(mu); }

/* --- NTP --- */
void* ccs_ntp_create(void) {
    NtpClock *n = (NtpClock*)calloc(1, sizeof(NtpClock));
    if (n) ntp_init(n);
    return n;
}
void ccs_ntp_free(void *n) { free(n); }
int ccs_ntp_process(void *n, double t1, double t2, double t3, double t4) {
    return ntp_process((NtpClock*)n, t1, t2, t3, t4);
}
double ccs_ntp_get_offset_us(void *n) {
    return n ? ((NtpClock*)n)->best_offset_us : 0.0;
}
double ccs_ntp_get_delay_us(void *n) {
    return n ? ((NtpClock*)n)->best_delay_us : 0.0;
}
double ccs_ntp_get_uncertainty_us(void *n) {
    return n ? ((NtpClock*)n)->uncertainty_us : 1e6;
}
int ccs_ntp_is_locked(void *n) {
    return n ? ((NtpClock*)n)->locked : 0;
}
double ccs_ntp_skew_ppb(void *n) { return ntp_skew((NtpClock*)n); }
double ccs_ntp_compensate(void *n, double rt) {
    return ntp_compensate((NtpClock*)n, rt);
}
double ccs_ntp_holdover_time(void *n, double mo) {
    return ntp_holdover((NtpClock*)n, mo);
}

/* --- Reorder --- */
void* ccs_reorder_create(int wsz, double mw) {
    RBuffer *r = (RBuffer*)calloc(1, sizeof(RBuffer));
    if (r) rb_init(r, wsz, mw);
    return r;
}
void ccs_reorder_free(void *r) { free(r); }
int ccs_reorder_insert(void *r, uint64_t seq, const double *d, int dim,
                        double at, double gt) {
    return rb_insert((RBuffer*)r, seq, d, dim, at, gt);
}
int ccs_reorder_check_timeout(void *r, double now) {
    return rb_timeout((RBuffer*)r, now);
}
double ccs_reorder_density(void *r) { return rb_density((RBuffer*)r); }
double ccs_reorder_avg_wait(void *r) { return rb_avgwait((RBuffer*)r); }
void ccs_reorder_stats(void *r, uint64_t *rx, uint64_t *in,
                        uint64_t *la, uint64_t *du, uint64_t *dr) {
    RBuffer *rb = (RBuffer*)r;
    if (!rb) return;
    if (rx) *rx = rb->total_rx;
    if (in) *in = rb->total_inorder;
    if (la) *la = rb->total_late;
    if (du) *du = rb->total_dups;
    if (dr) *dr = rb->total_drops;
}

/* --- Token Bucket --- */
void* ccs_token_bucket_create(double rate, double burst) {
    TBucket *t = (TBucket*)calloc(1, sizeof(TBucket));
    if (t) tb_init(t, rate, burst);
    return t;
}
void ccs_token_bucket_free(void *t) { free(t); }
int ccs_token_bucket_test(void *t, double sz, double now) {
    return tb_consume((TBucket*)t, sz, now);
}
double ccs_token_bucket_conformance(void *t) { return tb_rate((TBucket*)t); }

/* --- Congestion --- */
void* ccs_congestion_create(int thr) {
    DCongestion *d = (DCongestion*)calloc(1, sizeof(DCongestion));
    if (d) dc_init(d, thr);
    return d;
}
void ccs_congestion_free(void *d) { free(d); }
int ccs_congestion_update(void *d, double rtt_us) {
    return d ? dc_update((DCongestion*)d, rtt_us) : 0;
}
double ccs_congestion_severity(void *d) {
    return d ? ((DCongestion*)d)->severity : 0.0;
}

/* --- WRR --- */
void* ccs_wrr_create(void) {
    WRR *w = (WRR*)calloc(1, sizeof(WRR));
    if (w) wrr_init(w);
    return w;
}
void ccs_wrr_free(void *w) { free(w); }
int ccs_wrr_add_path(void *w, int id, double bw) {
    return wrr_add((WRR*)w, id, bw);
}
int ccs_wrr_schedule(void *w, int sz) { return wrr_sched((WRR*)w, sz); }

/* --- OWD / TSN --- */
double ccs_owd_minrtt(const double *rtt, int n) { return owd_minrtt(rtt, n); }
double ccs_owd_asymmetry(const double *rtt, int n) { return owd_asym(rtt, n); }
