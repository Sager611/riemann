/* isieve2.c (ImqClass2): certify h(-d) > N for every candidate d in [LO, HI] (d = |disc|; candidates: every
 * d = 3, 7, 11, 15, 4, 8 mod 16, fundamental or not), except the fundamental d written to OUT.t4.
 *
 * Bound (Watkins' minima, as isieve.c of ImqClass): for fundamental -d, h(-d) >= M(d) = sum_{a <= X} r(a),
 * X = floor(sqrt(d/4)), r multiplicative: r(p^k) = 2 (p split), 0 (p inert), r(p) = 1, r(p^k>1) = 0 (p ramified).
 * Residue classes c mod QS = 16*3*5*7*11*13 fix the symbols of p <= 13; G_c = integers m built from the split and
 * ramified p <= 13 (ramified ones to the first power), F_c(t) = sum_{m in G_c, m <= t} r(m). Every a = m*n with
 * m in G_c and n built from primes > 13 has r(a) = r(m) r(n), and distinct (m, n) give distinct a. With
 * Xs = X(segment start) <= X(d), summing over a family of distinct n with r(n) > 0:
 *   M(d) >= F_c(Xs)                                                      n = 1
 *         + sum_{17 <= p <= Xs/17, p split} W1(p),  W1(p) = 2 sum_{k>=1} F_c(Xs/p^k)   n = p^k   (phase A)
 *         + sum_{17 <= p < q, pq <= Xs, p, q split} 4 F_c(Xs/(pq))                     n = p*q   (pairs)
 *         + sum_{Xs/17 < p <= min(Xs, PT), p split} 2 F_c(Xs/p)                        n = p     (phase B)
 * (all divisions floor; p^k, pq <= Xs; q <= Xs/17 in every pair since p >= 17). Every weight w is then rounded
 * down to w' = the largest number 2^k or 3*2^(k-1) that is <= w (w' >= 2w/3; isieve.c rounded 2F_c(Xs/p) down to a
 * power of 2, w' >= w/2). New against isieve.c: the 2-bit rounding, the prime powers p^k in W1 and the pairs
 * n = pq. Option -1 drops the pairs (singles only: the bound of isieve.c with the finer rounding).
 * Along each phase the weights are non-increasing (primes ascending; pairs sorted by pq ascending), so the items of
 * one rounded weight form runs. Tier 1 (bit-parallel, 64 lanes d = c + QS*j, consecutive j, per word): the split
 * bits of p > 13 come from per-prime tables T_p (bit x set iff -QS*x is a nonzero square mod p); a run is counted
 * with carry-save (Harley-Seal) adders in chunks, each chunk count times w' is added to 16 bit-sliced planes (exact
 * integer addition; a carry out of plane 15 means a sum >= 2^16 > N+1: certified), and a lane is certified when
 * its sum reaches N+1 (bit-sliced exact comparison). Words are processed in groups of VW consecutive words of one
 * class (VW*64 lanes, one vector per plane) through phase A, pairs and phase B until all lanes are certified.
 * Tier 3 (each uncertified lane): exact M(d), capped at N+1, by depth-first enumeration of the a <= X with r(a) > 0
 * (Kronecker symbols from the tables, checked against Jacobi symbols at start and in test mode).
 * Tier 4: fundamental d with M(d) <= N go to OUT.t4 ("d M"), to be decided exactly (t4.gp).
 * Usage: isieve2 N LO HI OUT [-1] [-P PT] [-T K]
 *   -P PT: tables only for primes <= PT (memory; default 250000); phase B stops at min(Xs, PT).
 *   -T K: test mode: every K-th certified lane (hash of d) gets its exact M(d) by an independent Jacobi-symbol
 *         enumeration (cap 2^20): it must be >= the lane's tier-1 sum (when that sum did not overflow) and >= N+1;
 *         every tier-3 value is recomputed with Jacobi symbols and must agree.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <mach/mach_time.h>
typedef unsigned long long u64;
typedef unsigned int u32;
typedef unsigned short u16;
#define QS 240240ULL
#define NPAT 729
#define WPC 128       /* words per class per segment */
#define NPL 16        /* accumulator planes */
#define WCAP 32767    /* weights are capped (a smaller weight keeps the bound valid) */
#ifndef VW
#define VW 8          /* words per vector group */
#endif
#ifndef PFD
#define PFD 2         /* prefetch distance in groups */
#endif
typedef unsigned long long vec __attribute__((ext_vector_type(VW)));
static inline int vnone(vec v) { for (int e = 0; e < VW; e++) if (v[e]) return 0; return 1; }
static const int SPR[6] = {2, 3, 5, 7, 11, 13};

static int N; static u64 LO, HI, PT = 250000; static long testK = 0; static int nopairs = 0;
static u32 *pr; static int npr;           /* all primes <= PALL = isqrt(HI)+2 */
static int i17;                           /* index of 17 */
static int nmed;                          /* table primes pr[i17 .. i17+nmed-1]: 17 <= p <= min(XMAX, PT) */
static u32 *pmed;                         /* pmed[i] = pr[i17+i] */
static u64 **T; static u32 *invQ, *rG;     /* rG[i] = 64*VW mod p */

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static u64 isqrt(u64 x) { u64 r = (u64)sqrtl((long double)x); while (r * r > x) r--; while ((r + 1) * (r + 1) <= x) r++; return r; }
static u64 Xof(u64 d) { return isqrt(d / 4); } /* floor(sqrt(d/4)) = floor(sqrt(floor(d/4))) */

static int jacobi(long long a, long long n) { int t = 1; a %= n; if (a < 0) a += n;
  while (a) { while (!(a & 1)) { a >>= 1; long long r = n & 7; if (r == 3 || r == 5) t = -t; }
    long long x = a; a = n; n = x; if ((a & 3) == 3 && (n & 3) == 3) t = -t; a %= n; }
  return n == 1 ? t : 0; }
/* Kronecker (-d | p) for prime p, independent of the tables */
static int kron(u64 d, u32 p) {
  if (p == 2) { if (!(d & 1)) return 0; return ((d & 7) == 7) ? 1 : -1; }
  u64 r = d % p; if (!r) return 0; return jacobi((long long)(p - r), p); }
/* the same from the tables (p = pmed[i]) */
static inline int kron_tab(u64 d, int i) {
  u32 p = pmed[i]; u64 r = d % p; if (!r) return 0;
  u64 x = (r * invQ[i]) % p; return ((T[i][x >> 6] >> (x & 63)) & 1) ? 1 : -1; }

/* ---------- patterns of the primes <= 13 ---------- */
typedef struct { int n; u64 *m; u64 *pre; } pat_t;
static pat_t PAT[NPAT];
static u64 XMAX;
static int st_of(int pat, int i) { for (int t = 0; t < i; t++) pat /= 3; return pat % 3; } /* 0 inert 1 split 2 ram */
static int cmpu(const void *a, const void *b) { u64 x = *(const u64 *)a, y = *(const u64 *)b; return x < y ? -1 : x > y; }
static void build_patterns(void) {
  size_t cap = 4000000;
  u64 *tmpm = malloc(sizeof(u64) * cap); u64 *tmpw = malloc(sizeof(u64) * cap);
  for (int pat = 0; pat < NPAT; pat++) {
    size_t n = 1; tmpm[0] = 1; tmpw[0] = 1;
    for (int i = 0; i < 6; i++) {
      int s = st_of(pat, i); if (s == 0) continue;
      size_t cur = n; u64 p = SPR[i];
      for (size_t j = 0; j < cur; j++) {
        u64 m = tmpm[j] * p; u64 w = tmpw[j] * (s == 1 ? 2 : 1);
        while (m <= XMAX) { if (n >= cap) { fprintf(stderr, "pattern overflow\n"); exit(5); }
          tmpm[n] = m; tmpw[n] = w; n++; if (s == 2) break; m *= p; }
      }
    }
    u64 *pack = malloc(sizeof(u64) * n);
    for (size_t j = 0; j < n; j++) pack[j] = (tmpm[j] << 8) | tmpw[j];  /* w <= 2^6 */
    qsort(pack, n, 8, cmpu);
    PAT[pat].n = (int)n; PAT[pat].m = malloc(sizeof(u64) * n); PAT[pat].pre = malloc(sizeof(u64) * n);
    u64 acc = 0;
    for (size_t j = 0; j < n; j++) { PAT[pat].m[j] = pack[j] >> 8; acc += pack[j] & 255; PAT[pat].pre[j] = acc; }
    free(pack);
  }
  free(tmpm); free(tmpw);
}
static u64 Fp(int pat, u64 t) { /* sum of r(m), m in G, m <= t */
  pat_t *P = &PAT[pat]; int lo = 0, hi = P->n;
  while (lo < hi) { int mid = (lo + hi) / 2; if (P->m[mid] <= t) lo = mid + 1; else hi = mid; }
  return lo ? P->pre[lo - 1] : 0; }

/* ---------- tier 3 ---------- */
static u64 *st_n; static long long *st_w; static int *st_i; static signed char *chi3;
static long long dfs_M(u64 X, int np, long long cap) {
  int top = 0; long long tot = 0;
  st_n[0] = 1; st_w[0] = 1; st_i[0] = 0; top = 1;
  while (top) {
    top--; u64 n = st_n[top]; long long w = st_w[top]; int idx = st_i[top];
    tot += w; if (tot >= cap) return cap;
    for (int i = idx; i < np; i++) {
      u64 p = pr[i]; if (n * p > X) break;
      int s = chi3[i]; if (s < 0) continue;
      if (s > 0) { u64 m = n * p; while (m <= X) { st_n[top] = m; st_w[top] = 2 * w; st_i[top] = i + 1; top++; if (m > X / p) break; m *= p; } }
      else { st_n[top] = n * p; st_w[top] = w; st_i[top] = i + 1; top++; }
      if (top > 3900000) { fprintf(stderr, "dfs stack\n"); exit(6); }
    }
  }
  return tot;
}
static long long tier3_tab(u64 d, long long cap) { /* exact M(d) capped at cap; symbols from the tables */
  u64 X = Xof(d); int np = 0;
  while (np < npr && pr[np] <= X) {
    int i = np - i17;
    chi3[np] = (signed char)((np < i17 || i >= nmed) ? kron(d, pr[np]) : kron_tab(d, i)); np++; }
  return dfs_M(X, np, cap);
}
static long long tier3_jac(u64 d, long long cap) { /* the same with Jacobi symbols only (independent of the tables) */
  u64 X = Xof(d); int np = 0;
  while (np < npr && pr[np] <= X) { chi3[np] = (signed char)kron(d, pr[np]); np++; }
  return dfs_M(X, np, cap);
}
static int fundamental(u64 d) { int r16 = d & 15; if (!((d & 3) == 3 || r16 == 4 || r16 == 8)) return 0;
  u64 u = (d & 1) ? d : d / 4;
  for (int i = 1; i < npr; i++) { u64 p = pr[i]; if (p * p > u) break; if (u % (p * p) == 0) return 0; }
  return 1; }

/* ---------- bit-sliced accumulator ---------- */
static inline void addsh(vec *acc, vec *ovf, const vec *C, int nc, int b) { /* acc += C << b (per lane, exact) */
  vec carry = 0;
  for (int z = b; z < NPL; z++) {
    int q = z - b; vec c = q < nc ? C[q] : (vec)0, a = acc[z];
    acc[z] = a ^ c ^ carry; carry = (a & c) | (carry & (a ^ c));
  }
  for (int q = NPL - b; q < nc; q++) if (q >= 0) carry |= C[q];   /* planes beyond 2^15: sum >= 2^16 */
  *ovf |= carry; }
static inline vec geq(const vec *acc, u64 needv) { /* lanes with acc >= needv (needv >= 1), exact */
  u64 thr = needv - 1; vec gt = 0, eq = ~(vec)0;
  for (int z = NPL - 1; z >= 0; z--) { vec a = acc[z]; if ((thr >> z) & 1) eq &= a; else { gt |= eq & a; eq &= ~a; } }
  return gt; }
static inline vec fetch(int i, u32 o) { /* bits o .. o + 64*VW - 1 of T_i, lane t of word e = bit o + 64e + t */
  const u64 *Tp = T[i] + (o >> 6); u32 sh = o & 63; vec lo, hi;
#ifndef NOPF
  __builtin_prefetch(Tp + PFD * VW + 8);   /* the window of group g + PFD (a hint; never faults) */
#endif
  memcpy(&lo, Tp, sizeof lo); memcpy(&hi, Tp + 1, sizeof hi);
  return (lo >> (u64)sh) | ((hi << 1ULL) << (u64)(63 - sh)); }
static inline u32 round2(u64 w) { if (w > WCAP) w = WCAP; if (w == 0) return 0; int k = 63 - __builtin_clzll(w);
  if (k >= 1 && ((w >> (k - 1)) & 1)) return (3u << (k - 1)); return 1u << k; }

/* per-class offsets: cur[i] = offset of group wl[i] (its first lane) for table prime i; primes [0, imax) initialized */
static u32 *cur; static u16 *wl; static int imax; static u64 cls_c, cls_jlo;
__attribute__((unused)) static inline u32 offset_of(int i, int g) {
  u32 p = pmed[i];
  if (i >= imax) { cur[i] = (u32)(((cls_jlo % p) + (u64)(cls_c % p) * invQ[i]) % p); wl[i] = 0; imax = i + 1; }
  int dl = g - (int)wl[i]; u32 o = cur[i];
  if (dl == 1) { o += rG[i]; if (o >= p) o -= p; }
  else if (dl > 1) { u64 s = 64ULL * VW * (u64)dl; if (s >= p) s %= p; o += (u32)s; if (o >= p) o -= p; }
  cur[i] = o; wl[i] = (u16)g; return o; }

static vec *xs; static u16 *ppi, *ppj;
typedef struct { int end; u32 w; } run_t;
#define MAXRUN 40
#define CSA(h, l, a, b_, cc) { vec u_ = (a) ^ (b_); h = ((a) & (b_)) | (u_ & (cc)); l = u_ ^ (cc); }
/* one run [s, e) of items of rounded weight wq; GETX(k) gives the split bits of item k. On exit (all lanes of the
   word certified) jumps to the label done. */
#define RUN(s, e, wq, GETX, WORK) { \
  int chunk_ = (wq) < 128 ? 256 : ((wq) < 1024 ? 64 : 16); \
  for (int k0_ = (s); k0_ < (e); k0_ += chunk_) { \
    int k1_ = k0_ + chunk_ < (e) ? k0_ + chunk_ : (e); \
    vec o_ = 0, t_ = 0, f_ = 0, e_ = 0, c16[5] = {0, 0, 0, 0, 0}; \
    for (int k_ = k0_; k_ < k1_; k_ += 16) { \
      vec x[16]; int nb_ = k1_ - k_ < 16 ? k1_ - k_ : 16; \
      for (int t_i = 0; t_i < nb_; t_i++) x[t_i] = GETX(k_ + t_i) & live; \
      for (int t_i = nb_; t_i < 16; t_i++) x[t_i] = 0; \
      vec tA, tB, fA, fB, eA, eB, s16; \
      CSA(tA, o_, o_, x[0], x[1]); CSA(tB, o_, o_, x[2], x[3]); CSA(fA, t_, t_, tA, tB); \
      CSA(tA, o_, o_, x[4], x[5]); CSA(tB, o_, o_, x[6], x[7]); CSA(fB, t_, t_, tA, tB); \
      CSA(eA, f_, f_, fA, fB); \
      CSA(tA, o_, o_, x[8], x[9]); CSA(tB, o_, o_, x[10], x[11]); CSA(fA, t_, t_, tA, tB); \
      CSA(tA, o_, o_, x[12], x[13]); CSA(tB, o_, o_, x[14], x[15]); CSA(fB, t_, t_, tA, tB); \
      CSA(eB, f_, f_, fA, fB); \
      CSA(s16, e_, e_, eA, eB); \
      { vec cy_ = s16; for (int q_ = 0; q_ < 5; q_++) { vec tt_ = c16[q_] & cy_; c16[q_] ^= cy_; cy_ = tt_; } } \
    } \
    WORK += (u64)(k1_ - k0_); \
    { vec C_[9] = {o_, t_, f_, e_, c16[0], c16[1], c16[2], c16[3], c16[4]}; u32 ww_ = (wq); \
      while (ww_) { int b_ = __builtin_ctz(ww_); ww_ &= ww_ - 1; addsh(acc, &ovf, C_, 9, b_); } } \
    pass = ovf | geq(acc, needv); live = vm & ~pass; if (vnone(live)) goto done; \
  } }

/* group state of the current class (at most MAXG groups of VW words) */
#define MAXG (WPC / VW)
static vec GACC[MAXG][NPL], GOVF[MAXG], GPASS[MAXG], GLIVE[MAXG], GVM[MAXG];
static vec GO[MAXG], GT[MAXG], GF[MAXG], GE[MAXG], GC16[MAXG][5];
static u32 *m512t; static int nsmall;   /* m512t[MAXG*i + g] = 64*VW*g mod p_i for the small table primes i < nsmall */
static inline u32 goff(int i, int g) { u32 p = pmed[i]; u32 o = cur[i] + (i < nsmall ? m512t[MAXG * i + g] : (u32)(64 * VW) * (u32)g); if (o >= p) o -= p; return o; }
static inline void base_init(int i) { while (imax <= i) { u32 p = pmed[imax]; cur[imax] = (u32)(((cls_jlo % p) + (u64)(cls_c % p) * invQ[imax]) % p); imax++; } }
/* items [s, e) of rounded weight wq, prime-major over the live groups gl[0..ngl); STORE: keep the bits in xs (phase A).
   Groups whose lanes are all certified leave the list with stage ST; if none is left, jump to class_done. */
/* items kb_+t, kb_+t+1 of group g_ (zero beyond nb_) into a0_, a1_ */
#define FX1(t, a) { if ((t) < nb_) { int i_ = kb_ + (t); vec y_ = fetch(i_, goff(i_, g_)); if (STORE_) xsg_[i_] = y_; a = y_ & lv_; } else a = 0; }
#define FX2(t) { FX1((t), a0_); FX1((t) + 1, a1_); }
#define PHASEPM(s, e, wq, STORE, WORK, ST) { const int STORE_ = (STORE); \
  int chunk_ = (wq) < 128 ? 256 : ((wq) < 1024 ? 64 : 16); \
  for (int k0_ = (s); k0_ < (e); k0_ += chunk_) { \
    int k1_ = k0_ + chunk_ < (e) ? k0_ + chunk_ : (e); \
    for (int gi_ = 0; gi_ < ngl; gi_++) { int g_ = gl[gi_]; GO[g_] = GT[g_] = GF[g_] = GE[g_] = 0; for (int q_ = 0; q_ < 5; q_++) GC16[g_][q_] = 0; } \
    for (int kb_ = k0_; kb_ < k1_; kb_ += 16) { \
      int nb_ = k1_ - kb_ < 16 ? k1_ - kb_ : 16; \
      base_init(kb_ + nb_ - 1); \
      { int pe_ = kb_ + 32 < (e) ? kb_ + 32 : (e); if (pe_ > kb_ + 16) { base_init(pe_ - 1); \
          for (int i_ = kb_ + 16; i_ < pe_; i_++) { const char *q_ = (const char *)(T[i_] + (goff(i_, gl[0]) >> 6)); \
            __builtin_prefetch(q_); __builtin_prefetch(q_ + 64); __builtin_prefetch(q_ + 64 * VW); } } } \
      for (int gi_ = 0; gi_ < ngl; gi_++) { \
        int g_ = gl[gi_]; vec lv_ = GLIVE[g_]; vec *xsg_ = xs + (size_t)g_ * (nAmax + 1); \
        vec o_ = GO[g_], t_ = GT[g_], f_ = GF[g_], e_ = GE[g_]; \
        vec tA, tB, fA, fB, eA, eB, s16, a0_, a1_; \
        FX2(0); CSA(tA, o_, o_, a0_, a1_); FX2(2); CSA(tB, o_, o_, a0_, a1_); CSA(fA, t_, t_, tA, tB); \
        FX2(4); CSA(tA, o_, o_, a0_, a1_); FX2(6); CSA(tB, o_, o_, a0_, a1_); CSA(fB, t_, t_, tA, tB); \
        CSA(eA, f_, f_, fA, fB); \
        FX2(8); CSA(tA, o_, o_, a0_, a1_); FX2(10); CSA(tB, o_, o_, a0_, a1_); CSA(fA, t_, t_, tA, tB); \
        FX2(12); CSA(tA, o_, o_, a0_, a1_); FX2(14); CSA(tB, o_, o_, a0_, a1_); CSA(fB, t_, t_, tA, tB); \
        CSA(eB, f_, f_, fA, fB); \
        CSA(s16, e_, e_, eA, eB); \
        GO[g_] = o_; GT[g_] = t_; GF[g_] = f_; GE[g_] = e_; \
        { vec cy_ = s16; for (int q_ = 0; q_ < 5; q_++) { vec tt_ = GC16[g_][q_] & cy_; GC16[g_][q_] ^= cy_; cy_ = tt_; } } \
      } \
    } \
    WORK += (u64)(k1_ - k0_) * (u64)ngl; \
    for (int gi_ = 0; gi_ < ngl; ) { \
      int g_ = gl[gi_]; \
      vec C_[9] = {GO[g_], GT[g_], GF[g_], GE[g_], GC16[g_][0], GC16[g_][1], GC16[g_][2], GC16[g_][3], GC16[g_][4]}; u32 ww_ = (wq); \
      while (ww_) { int b_ = __builtin_ctz(ww_); ww_ &= ww_ - 1; addsh(GACC[g_], &GOVF[g_], C_, 9, b_); } \
      GPASS[g_] = GOVF[g_] | geq(GACC[g_], needv); GLIVE[g_] = GVM[g_] & ~GPASS[g_]; \
      if (vnone(GLIVE[g_])) { stage[g_] = (ST); gl[gi_] = gl[--ngl]; } else gi_++; \
    } \
    if (!ngl) goto class_done; \
  } }

/* ---------- main ---------- */
int main(int argc, char **argv) {
  if (argc < 5) { fprintf(stderr, "usage: isieve2 N LO HI OUT [-1] [-P PT] [-T K]\n"); return 1; }
  long Nl = atol(argv[1]); LO = strtoull(argv[2], 0, 10); HI = strtoull(argv[3], 0, 10);
  const char *out = argv[4];
  for (int a = 5; a < argc; a++) {
    if (!strcmp(argv[a], "-1")) nopairs = 1;
    else if (!strcmp(argv[a], "-P") && a + 1 < argc) PT = strtoull(argv[++a], 0, 10);
    else if (!strcmp(argv[a], "-T") && a + 1 < argc) testK = atol(argv[++a]);
    else { fprintf(stderr, "bad option %s\n", argv[a]); return 1; }
  }
  if (Nl < 1 || Nl > 30000 || LO < 1000000 || HI < LO || HI > 4000000000000ULL || PT < 1000) { fprintf(stderr, "bad args\n"); return 1; }
  N = (int)Nl;
  double t0 = now();
  u64 PALL = isqrt(HI) + 2; XMAX = Xof(HI) + 1;
  char *comp = calloc(PALL + 1, 1); pr = malloc(sizeof(u32) * (PALL / 2 + 100)); npr = 0;
  for (u64 i = 2; i <= PALL; i++) { if (comp[i]) continue; pr[npr++] = (u32)i; for (u64 j = i * i; j <= PALL; j += i) comp[j] = 1; }
  free(comp);
  i17 = 0; while (pr[i17] < 17) i17++;
  u64 plim = XMAX < PT ? XMAX : PT;
  nmed = 0; while (i17 + nmed < npr && pr[i17 + nmed] <= plim) nmed++;
  pmed = pr + i17;
  T = malloc(sizeof(u64 *) * nmed); invQ = malloc(sizeof(u32) * nmed); rG = malloc(sizeof(u32) * nmed);
  size_t tbytes = 0;
  unsigned char *qr = malloc(plim + 2);
  for (int i = 0; i < nmed; i++) {
    u32 p = pmed[i];
    memset(qr, 0, p + 1);
    { u64 sq = 0; for (u64 x = 1; x <= (p - 1) / 2; x++) { sq += 2 * x - 1; if (sq >= p) sq -= p; qr[sq] = 1; } } /* sq = x^2 mod p; qr[0] = qr[p] = 0 */
    u64 q = QS % p, inv = 0; /* QS^{-1} mod p */
    { long long a = (long long)q, m = p, x0 = 0, x1 = 1; long long mm = m; while (a > 1) { long long qq = a / mm, t = mm; mm = a % mm; a = t; t = x0; x0 = x1 - qq * x0; x1 = t; } if (x1 < 0) x1 += m; inv = (u64)x1; }
    if ((q * inv) % p != 1) { fprintf(stderr, "inv fail %u\n", p); return 3; }
    invQ[i] = (u32)inv; rG[i] = (u32)((64ULL * VW) % p);
    u64 nbits = (u64)p + 64ULL * VW + 64;                    /* any window [o, o + 64*VW + 64), o < p */
    size_t nw = (nbits + 63) / 64 + 1; T[i] = calloc(nw, 8); tbytes += nw * 8;
    { u64 dm = 0, step = q, word = 0; u64 j;               /* dm = QS*j mod p; bit j = qr[p - dm] (dm = 0: qr[p] = 0) */
      for (j = 0; j < nbits; j++) {
        word |= (u64)qr[p - dm] << (j & 63);
        if ((j & 63) == 63) { T[i][j >> 6] = word; word = 0; }
        dm += step; if (dm >= p) dm -= p;
      }
      if (j & 63) T[i][j >> 6] = word; }
  }
  free(qr);
  double tT = now() - t0;
  build_patterns();
  double tP = now() - t0;
  st_n = malloc(sizeof(u64) * 4000000); st_w = malloc(sizeof(long long) * 4000000); st_i = malloc(sizeof(int) * 4000000);
  chi3 = malloc(npr);
  /* table self-check: T against kron on random d */
  { u64 s = 12345; long bad = 0; for (int t = 0; t < 400000; t++) { s ^= s << 13; s ^= s >> 7; s ^= s << 17; u64 d = LO + s % (HI - LO + 1);
      int i = (int)((s >> 20) % (u64)nmed); u32 p = pmed[i]; u64 j = ((d % p) * invQ[i]) % p;
      int split = (T[i][j >> 6] >> (j & 63)) & 1; if (split != (kron(d, p) == 1)) bad++;
      if (kron_tab(d, i) != kron(d, p)) bad++; }
    if (bad) { fprintf(stderr, "TABLE SELF-CHECK FAIL %ld\n", bad); return 4; } }
  /* classes */
  static int cls_pat[QS]; static unsigned char cls_ok[QS];
  for (u64 c = 0; c < QS; c++) {
    int r16 = (int)(c & 15); int s2;
    if (r16 == 3 || r16 == 11) s2 = 0; else if (r16 == 7 || r16 == 15) s2 = 1; else if (r16 == 4 || r16 == 8) s2 = 2; else { cls_ok[c] = 0; continue; }
    int pat = s2, mul = 3;
    for (int i = 1; i < 6; i++) { int k = kron(c, SPR[i]); int s = (k == 0) ? 2 : (k == 1 ? 1 : 0); pat += s * mul; mul *= 3; }
    cls_ok[c] = 1; cls_pat[c] = pat;
  }
  /* pairs (i < j table primes, p_i p_j <= XMAX), sorted by p_i p_j ascending */
  int nAmax = 0; while (nAmax < nmed && pmed[nAmax] <= XMAX / 17) nAmax++;
  size_t npairmax = 0;
  for (int i = 0; i < nAmax; i++) for (int j = i + 1; j < nAmax && (u64)pmed[i] * pmed[j] <= XMAX; j++) npairmax++;
  u64 *pk = malloc(sizeof(u64) * (npairmax + 1)); size_t npk = 0;
  for (int i = 0; i < nAmax; i++) for (int j = i + 1; j < nAmax && (u64)pmed[i] * pmed[j] <= XMAX; j++) pk[npk++] = (((u64)pmed[i] * pmed[j]) << 28) | ((u64)i << 14) | (u64)j;
  if (nAmax >= (1 << 14)) { fprintf(stderr, "too many phase-A primes\n"); return 7; }
  qsort(pk, npk, 8, cmpu);
  ppi = malloc(sizeof(u16) * (npk + 1)); ppj = malloc(sizeof(u16) * (npk + 1)); u64 *pprod = malloc(sizeof(u64) * (npk + 1));
  for (size_t k = 0; k < npk; k++) { pprod[k] = pk[k] >> 28; ppi[k] = (u16)((pk[k] >> 14) & 16383); ppj[k] = (u16)(pk[k] & 16383); }
  free(pk);
  int tmaxmax = (int)(XMAX / (17 * 19)) + 2;
  int *posT = malloc(sizeof(int) * (tmaxmax + 2));
  run_t *runA = malloc(sizeof(run_t) * NPAT * MAXRUN), *runP = malloc(sizeof(run_t) * NPAT * MAXRUN), *runB = malloc(sizeof(run_t) * NPAT * MAXRUN);
  int *nrA = malloc(sizeof(int) * NPAT), *nrP = malloc(sizeof(int) * NPAT), *nrB = malloc(sizeof(int) * NPAT);
  xs = malloc(sizeof(vec) * (size_t)MAXG * (nAmax + 1));
  nsmall = 0; while (nsmall < nmed && pmed[nsmall] <= 64 * VW * (MAXG - 1)) nsmall++;
  m512t = malloc(sizeof(u32) * MAXG * (nsmall + 1));
  for (int i = 0; i < nsmall; i++) for (int g = 0; g < MAXG; g++) m512t[MAXG * i + g] = (u32)((64ULL * VW * g) % pmed[i]);
  static u64 tpat[NPAT], lpat[NPAT];
  cur = malloc(sizeof(u32) * (nmed + 1)); wl = malloc(sizeof(u16) * (nmed + 1));
  char path[4096]; snprintf(path, sizeof path, "%s.t4", out); FILE *f4 = fopen(path, "w");
  snprintf(path, sizeof path, "%s.log", out); FILE *flog = fopen(path, "w");
  if (!f4 || !flog) { fprintf(stderr, "cannot open outputs\n"); return 1; }
  fprintf(flog, "isieve2 N=%d LO=%llu HI=%llu pairs %s PT=%llu: table primes %d (<= %u), T tables %.1f MB, phase-A primes <= %d, pairs <= %zu, setup %.1f s (tables %.1f, patterns %.1f)\n",
          N, LO, HI, nopairs ? "off" : "on", PT, nmed, pmed[nmed - 1], tbytes / 1e6, nAmax, npk, now() - t0, tT, tP);
  double tsetup = now();
  fflush(flog);
  u64 lanes = 0, cert1 = 0, certF = 0, t3calls = 0, t3pass = 0, nonfund = 0, t4 = 0, tested = 0, testbad = 0, t3chk = 0;
  u64 cA = 0, cP = 0, cB = 0, wA = 0, wP = 0, wB = 0, words = 0, tmA = 0, tmP = 0, tmB = 0, tmR = 0;
  mach_timebase_info_data_t tbi; mach_timebase_info(&tbi);
  u64 LSEG = QS * 64 * WPC;
  int nseg = 0;
  for (u64 seg_lo = LO; seg_lo <= HI; seg_lo += LSEG) {
    u64 seg_hi = seg_lo + LSEG - 1; if (seg_hi > HI) seg_hi = HI;   /* inclusive */
    u64 Xs = Xof(seg_lo);
    /* segment tables */
    static long long need[NPAT];
    int nA = 0; while (nA < nmed && pmed[nA] <= Xs / 17) nA++;
    int nBend = 0; while (nBend < nmed && pmed[nBend] <= Xs) nBend++;
    int nP = 0; while ((size_t)nP < npk && pprod[nP] <= Xs) nP++;
    if (nopairs) nP = 0;
    int tmax = (int)(Xs / (17 * 19)); if (tmax + 1 > tmaxmax) { fprintf(stderr, "tmax\n"); return 7; }
    for (int k = 0; k < nP; k++) if (ppj[k] >= nA || Xs / pprod[k] > (u64)tmax || Xs / pprod[k] < 1) { fprintf(stderr, "pair index\n"); return 7; }
    /* posT[t] = number of pairs k < nP with Xs/pprod[k] >= t (pairs sorted by product: a prefix), t = 1..tmax+1 */
    { int k = nP; for (int t = 1; t <= tmax + 1; t++) { while (k > 0 && Xs / pprod[k - 1] < (u64)t) k--; posT[t] = k; } }
    for (int pat = 0; pat < NPAT; pat++) {
      need[pat] = (long long)N + 1 - (long long)Fp(pat, Xs);
      /* phase A runs: W1(p_i) non-increasing in i */
      run_t *ra = runA + pat * MAXRUN; int nr = 0; u32 prevw = 0;
      for (int i = 0; i < nA; i++) { u64 p = pmed[i], s = 0;
        for (u64 q = p; q <= Xs; ) { s += Fp(pat, Xs / q); if (q > Xs / p) break; q *= p; }
        u32 w = round2(2 * s);
        if (nr && w == prevw) ra[nr - 1].end = i + 1;
        else { if (nr && w > prevw) { fprintf(stderr, "A weights not monotone\n"); return 7; } if (nr >= MAXRUN) { fprintf(stderr, "runs A\n"); return 7; } ra[nr].end = i + 1; ra[nr].w = w; nr++; prevw = w; } }
      nrA[pat] = nr;
      /* pair runs: weight round2(4 F(t)), t = Xs/(pq) non-increasing along the sorted list */
      run_t *rp = runP + pat * MAXRUN; nr = 0; prevw = 0;
      { pat_t *P = &PAT[pat]; u32 wt_prev = 0; int tb = tmax; u64 Ft = Fp(pat, (u64)tmax);
        for (int t = tmax; t >= 1; t--) {
          u64 F = (t == tmax) ? Ft : Fp(pat, (u64)t); (void)P;
          u32 w = round2(4 * F);
          if (t < tmax && w != wt_prev) { /* run of weight wt_prev: t in [t+1, tb] */
            int s0 = posT[tb + 1], e0 = posT[t + 1];
            if (e0 > s0) { if (nr >= MAXRUN) { fprintf(stderr, "runs P\n"); return 7; } rp[nr].end = e0; rp[nr].w = wt_prev; nr++; }
            tb = t; }
          wt_prev = w; }
        { int s0 = posT[tb + 1], e0 = posT[1]; if (e0 > s0) { if (nr >= MAXRUN) { fprintf(stderr, "runs P\n"); return 7; } rp[nr].end = e0; rp[nr].w = wt_prev; nr++; } } }
      nrP[pat] = nr;
      /* phase B runs: weight round2(2 F(Xs/p)), Xs/p in [1, 16] non-increasing */
      run_t *rb = runB + pat * MAXRUN; nr = 0; prevw = 0;
      for (int i = nA; i < nBend; i++) { u32 w = round2(2 * Fp(pat, Xs / pmed[i]));
        if (nr && w == prevw) rb[nr - 1].end = i + 1;
        else { if (nr >= MAXRUN) { fprintf(stderr, "runs B\n"); return 7; } rb[nr].end = i + 1; rb[nr].w = w; nr++; prevw = w; } }
      nrB[pat] = nr;
    }
    for (u64 c = 0; c < QS; c++) {
      if (!cls_ok[c]) continue;
      u64 jlo = (seg_lo > c) ? (seg_lo - c + QS - 1) / QS : 0;
      u64 jhi = (seg_hi >= c) ? (seg_hi - c) / QS + 1 : 0;  /* exclusive */
      if (jhi <= jlo) continue;
      u64 nl = jhi - jlo; int nw = (int)((nl + 63) / 64); int ng = (nw + VW - 1) / VW;
      lanes += nl;
      int pat = cls_pat[c];
      if (need[pat] <= 0) { certF += nl; continue; }
      u64 needv = (u64)need[pat];
      const run_t *ra = runA + pat * MAXRUN, *rp = runP + pat * MAXRUN, *rb = runB + pat * MAXRUN;
      int nra = nrA[pat], nrp = nrP[pat], nrb = nrB[pat];
      cls_c = c; cls_jlo = jlo; imax = 0;
      u64 tmc0 = mach_absolute_time(), tmc1 = tmc0, tmc2 = tmc0; int phaseB = 0;
      int gl[MAXG], ngl = ng; int stage[MAXG];
      for (int g = 0; g < ng; g++) {
        vec v;
        for (int e = 0; e < VW; e++) { long long wi = (long long)VW * g + e; long long cnt = (long long)nl - 64LL * wi;
          v[e] = cnt >= 64 ? ~0ULL : (cnt <= 0 ? 0ULL : ((1ULL << cnt) - 1)); }
        GVM[g] = v; for (int z = 0; z < NPL; z++) GACC[g][z] = 0; GOVF[g] = 0; GPASS[g] = 0; GLIVE[g] = v; stage[g] = 3; gl[g] = g; }
      words += (u64)ng;
      /* phase A, prime-major over the live groups (each prime's window is streamed across the groups) */
      { int s = 0; for (int r = 0; r < nra; r++) { int e = ra[r].end; u32 wq = ra[r].w; PHASEPM(s, e, wq, 1, wA, 0); s = e; } }
      tmc1 = mach_absolute_time(); tmA += tmc1 - tmc0;
      /* pairs, per live group */
      for (int gi = 0; gi < ngl; ) {
        int g = gl[gi]; vec *acc = GACC[g]; vec ovf = GOVF[g], pass = GPASS[g], live = GLIVE[g], vm = GVM[g]; const vec *xg = xs + (size_t)g * (nAmax + 1);
        int fin = 0;
#define GETP(k) (xg[ppi[k]] & xg[ppj[k]])
        { int s = 0; for (int r = 0; r < nrp; r++) { RUN(s, rp[r].end, rp[r].w, GETP, wP); s = rp[r].end; } }
        goto pairs_end;
      done:
        fin = 1;
      pairs_end:
        GOVF[g] = ovf; GPASS[g] = pass; GLIVE[g] = live;
        if (fin) { stage[g] = 1; gl[gi] = gl[--ngl]; } else gi++;
      }
      tmc2 = mach_absolute_time(); tmP += tmc2 - tmc1; phaseB = 1;
      /* phase B, prime-major over the live groups */
      if (ngl) { int s = nA; for (int r = 0; r < nrb; r++) { int e = rb[r].end; u32 wq = rb[r].w; PHASEPM(s, e, wq, 0, wB, 2); s = e; } }
    class_done:
      if (phaseB) tmB += mach_absolute_time() - tmc2; else if (tmc1 == tmc0) tmA += mach_absolute_time() - tmc0; else tmP += mach_absolute_time() - tmc1;
      tpat[pat] += mach_absolute_time() - tmc0; lpat[pat] += nl;
      { u64 tmr0 = mach_absolute_time();
      for (int g = 0; g < ng; g++) {
        vec vm = GVM[g], pass = GPASS[g], ovf = GOVF[g]; vec *acc = GACC[g];
        u64 pc = 0; for (int e = 0; e < VW; e++) pc += (u64)__builtin_popcountll(pass[e] & vm[e]);
        if (stage[g] == 0) cA += pc; else if (stage[g] == 1) cP += pc; else cB += pc;
        cert1 += pc;
        for (int e = 0; e < VW; e++) {
        int w = VW * g + e;
        if (testK) {
          for (int t = 0; t < 64; t++) {
            if (!(((pass[e] & vm[e]) >> t) & 1)) continue;
            u64 d = c + QS * (jlo + 64ULL * w + t);
            if ((d * 2654435761ULL) % (u64)testK) continue;
            long long M = tier3_jac(d, 1LL << 20);
            tested++;
            if (M < (long long)N + 1) { testbad++; fprintf(flog, "TEST FAIL d=%llu M=%lld < N+1\n", d, M); }
            if (!((ovf[e] >> t) & 1)) {
              long long v = 0; for (int z = 0; z < NPL; z++) v |= (long long)((acc[z][e] >> t) & 1) << z;
              if (v > M) { testbad++; fprintf(flog, "TEST FAIL d=%llu bound %lld > M %lld\n", d, v, M); }
            }
          }
        }
        u64 sv = vm[e] & ~pass[e];
        while (sv) {
          int t = __builtin_ctzll(sv); sv &= sv - 1;
          u64 d = c + QS * (jlo + 64ULL * w + t);
          t3calls++;
          long long M = tier3_tab(d, (long long)N + 1);
          if (testK) { long long M2 = tier3_jac(d, (long long)N + 1); t3chk++; if (M2 != M) { testbad++; fprintf(flog, "TEST FAIL tier3 d=%llu tab %lld jac %lld\n", d, M, M2); } }
          if (M >= (long long)N + 1) { t3pass++; continue; }
          if (!fundamental(d)) { nonfund++; continue; }
          t4++; fprintf(f4, "%llu %lld\n", d, M);
        }
        }
      }
      tmR += mach_absolute_time() - tmr0; }
    }
    nseg++;
    if (nseg % 4 == 0 || seg_hi == HI) {
      fprintf(flog, "seg %d [%llu,%llu] X=%llu nA %d nP %d lanes %llu certF %llu cert1 %llu (A %llu, pairs %llu, B %llu) t3 %llu t3pass %llu nonfund %llu t4 %llu; per group (VW=%d): A %.1f pairs %.1f B %.1f; time A %.1f pairs %.1f B %.1f rest %.1f s; elapsed %.1f s (sieve %.1f s, %.2f ns/lane)\n",
              nseg, seg_lo, seg_hi, Xs, nA, nP, lanes, certF, cert1, cA, cP, cB, t3calls, t3pass, nonfund, t4, VW,
              words ? (double)wA / words : 0, words ? (double)wP / words : 0, words ? (double)wB / words : 0, tmA * 1e-9 * tbi.numer / tbi.denom, tmP * 1e-9 * tbi.numer / tbi.denom, tmB * 1e-9 * tbi.numer / tbi.denom, tmR * 1e-9 * tbi.numer / tbi.denom, now() - t0, now() - tsetup, lanes ? (now() - tsetup) * 1e9 / lanes : 0);
      fflush(flog); fflush(f4);
    }
    if (seg_hi == HI) break;
  }
  fprintf(flog, "DONE N=%d [%llu,%llu]: lanes %llu = certF %llu + cert1 %llu + t3 %llu (t3pass %llu, nonfund %llu, t4 %llu); tested %llu, tier3 rechecked %llu, bad %llu; %.1f s\n",
          N, LO, HI, lanes, certF, cert1, t3calls, t3pass, nonfund, t4, tested, t3chk, testbad, now() - t0);
  if (lanes != certF + cert1 + t3calls) fprintf(flog, "LANE COUNT MISMATCH\n");
  { u64 tt = 0; for (int p = 0; p < NPAT; p++) tt += tpat[p];
    fprintf(flog, "patterns by tier-1 time (pattern: states of 2,3,5,7,11,13 as i/s/r, share of time, share of lanes, ns/lane):");
    for (int r = 0; r < 12; r++) { int b = -1; for (int p = 0; p < NPAT; p++) if (tpat[p] && (b < 0 || tpat[p] > tpat[b])) b = p; if (b < 0) break;
      char nm[8]; for (int i = 0; i < 6; i++) nm[i] = "isr"[st_of(b, i)]; nm[6] = 0;
      fprintf(flog, " %s %.1f%% %.2f%% %.0f;", nm, 100.0 * tpat[b] / (tt ? tt : 1), 100.0 * lpat[b] / (lanes ? lanes : 1), lpat[b] ? tpat[b] * 1.0 * tbi.numer / tbi.denom / lpat[b] : 0); tpat[b] = 0; }
    fprintf(flog, "\n"); }
  fclose(flog); fclose(f4);
  return 0;
}
