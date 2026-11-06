/* isieve.c: certify h(-d) > N for every imaginary quadratic field Q(sqrt(-d)), d in [LO, HI] (d = |disc|),
 * except the d written to the tier-4 list.
 *
 * Bound (Watkins' minima): for fundamental -d, every primitive ideal of norm a <= X = floor(sqrt(d/4)) is the
 * first coefficient of a reduced form (a,b,c), b in (-a,a], and distinct ideals give distinct reduced forms, so
 * h(-d) >= M(d) = sum_{a<=X} r(a), r multiplicative: r(p^k) = 2 (p split), 0 (p inert), r(p) = 1, r(p^k>1) = 0
 * (p ramified).
 * Tier 1 (bit-parallel, per residue class c mod QS = 16*3*5*7*11*13): the symbols of p <= 13 are fixed by c;
 *   G_c = integers m composed of split/ramified p <= 13 (ramified to the first power), F_c(t) = sum_{m in G_c, m<=t} r(m).
 *   n = m*p (m in G_c, p > 13 split, m <= X/p) are distinct and r(n) = 2 r(m), so
 *   M(d) >= F_c(X) + sum_{13 < p <= X, p split} 2 F_c(X/p) >= F_c(Xs) + sum_p 2^{k_p} [p split],
 *   with Xs = X(segment start) <= X(d) and 2^{k_p} <= 2 F_c(Xs/p) (rounded down to a power of 2).
 *   Lanes (64 d of one class per word) are certified when this sum reaches N+1.
 * Tier 3 (per uncertified d, independent code: Jacobi symbols): exact M(d) by DFS, capped at N+1.
 * Tier 4: fundamental d with M(d) <= N go to OUT.t4 (exact class number later, PARI).
 * Usage: isieve N LO HI OUT [-T K]   (-T K: test mode, every K-th certified lane is re-checked by tier 3:
 *        exact M(d) must be >= the tier-1 bound of that lane)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
typedef unsigned long long u64;
typedef unsigned int u32;
#define QS 240240ULL
#define NPAT 729
#define WPC 128       /* words per class per segment */
#define NPL 16        /* accumulator planes */
static const int SPR[6] = {2, 3, 5, 7, 11, 13};

static int N; static u64 LO, HI; static long testK = 0;
static u32 *pr; static int npr;           /* all primes <= PALL = isqrt(HI)+1 */
static int i17;                           /* index of 17 */
static int nmed;                          /* medium primes: pr[i17 .. i17+nmed-1] <= PMED = isqrt(HI/4)+1 */
static u64 **T; static u32 *invQ, *r64;   /* per medium prime */

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

/* ---------- patterns of the primes <= 13 ---------- */
typedef struct { int n; u64 *m; u64 *pre; } pat_t;
static pat_t PAT[NPAT];
static u64 XMAX;
static int st_of(int pat, int i) { for (int t = 0; t < i; t++) pat /= 3; return pat % 3; } /* 0 inert 1 split 2 ram */
static int cmpu(const void *a, const void *b) { u64 x = *(const u64 *)a, y = *(const u64 *)b; return x < y ? -1 : x > y; }
static void build_patterns(void) {
  u64 *tmpm = malloc(sizeof(u64) * 20000000); u64 *tmpw = malloc(sizeof(u64) * 20000000);
  for (int pat = 0; pat < NPAT; pat++) {
    int n = 1; tmpm[0] = 1; tmpw[0] = 1;
    for (int i = 0; i < 6; i++) {
      int s = st_of(pat, i); if (s == 0) continue;
      int cur = n; u64 p = SPR[i];
      for (int j = 0; j < cur; j++) {
        u64 m = tmpm[j] * p; u64 w = tmpw[j] * (s == 1 ? 2 : 1);
        while (m <= XMAX) { tmpm[n] = m; tmpw[n] = w; n++; if (s == 2) break; m *= p; }
      }
    }
    u64 *pack = malloc(sizeof(u64) * n);
    for (int j = 0; j < n; j++) pack[j] = (tmpm[j] << 8) | tmpw[j];  /* w <= 2^6 */
    qsort(pack, n, 8, cmpu);
    PAT[pat].n = n; PAT[pat].m = malloc(sizeof(u64) * n); PAT[pat].pre = malloc(sizeof(u64) * n);
    u64 acc = 0;
    for (int j = 0; j < n; j++) { PAT[pat].m[j] = pack[j] >> 8; acc += pack[j] & 255; PAT[pat].pre[j] = acc; }
    free(pack);
  }
  free(tmpm); free(tmpw);
}
static u64 Fp(int pat, u64 t) { /* sum of r(m), m <= t */
  pat_t *P = &PAT[pat]; int lo = 0, hi = P->n;
  while (lo < hi) { int mid = (lo + hi) / 2; if (P->m[mid] <= t) lo = mid + 1; else hi = mid; }
  return lo ? P->pre[lo - 1] : 0; }
static u64 mfirst(int pat, u64 v) { /* smallest m with F(m) >= v, or 0 if none */
  pat_t *P = &PAT[pat]; int lo = 0, hi = P->n;
  while (lo < hi) { int mid = (lo + hi) / 2; if (P->pre[mid] >= v) hi = mid; else lo = mid + 1; }
  return lo < P->n ? P->m[lo] : 0; }

/* ---------- tier 3 ---------- */
static u64 *st_n; static long long *st_w; static int *st_i; static signed char *chi3;
static long long tier3_M(u64 d, long long cap) { /* exact M(d) capped at cap; Jacobi symbols */
  u64 X = Xof(d); int np = 0;
  while (np < npr && pr[np] <= X) { chi3[np] = (signed char)kron(d, pr[np]); np++; }
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
      if (top > 3900000) return cap; /* cannot happen for cap <= 2^20 */
    }
  }
  return tot;
}
static int fundamental(u64 d) { int r16 = d & 15; if (!((d & 3) == 3 || r16 == 4 || r16 == 8)) return 0;
  u64 u = (d & 1) ? d : d / 4;
  for (int i = 1; i < npr; i++) { u64 p = pr[i]; if (p * p > u) break; if (u % (p * p) == 0) return 0; }
  return 1; }

/* ---------- main ---------- */
int main(int argc, char **argv) {
  if (argc < 5) { fprintf(stderr, "usage: isieve N LO HI OUT [-T K]\n"); return 1; }
  N = atoi(argv[1]); LO = strtoull(argv[2], 0, 10); HI = strtoull(argv[3], 0, 10);
  const char *out = argv[4];
  if (argc >= 7 && !strcmp(argv[5], "-T")) testK = atol(argv[6]);
  if (LO < 1000000 || HI < LO || N + 1 >= (1 << 15)) { fprintf(stderr, "bad args\n"); return 1; }
  double t0 = now();
  u64 PALL = isqrt(HI) + 2; XMAX = Xof(HI) + 1;
  char *comp = calloc(PALL + 1, 1); pr = malloc(sizeof(u32) * (PALL / 2 + 100)); npr = 0;
  for (u64 i = 2; i <= PALL; i++) { if (comp[i]) continue; pr[npr++] = (u32)i; for (u64 j = i * i; j <= PALL; j += i) comp[j] = 1; }
  free(comp);
  i17 = 0; while (pr[i17] < 17) i17++;
  nmed = 0; while (i17 + nmed < npr && pr[i17 + nmed] <= XMAX) nmed++;
  T = malloc(sizeof(u64 *) * nmed); invQ = malloc(sizeof(u32) * nmed); r64 = malloc(sizeof(u32) * nmed);
  size_t tbytes = 0;
  char *qr = malloc(PALL + 1);
  for (int i = 0; i < nmed; i++) {
    u32 p = pr[i17 + i];
    memset(qr, 0, p); { u64 sq = 0; for (u64 x = 1; x < p; x++) { sq += 2 * x - 1; if (sq >= p) sq -= p; if (sq >= p) sq -= p; qr[sq] = 1; } } /* sq = x^2 mod p */
    u64 q = QS % p, inv = 0; /* QS^{-1} mod p */
    { long long a = (long long)q, m = p, x0 = 0, x1 = 1; long long mm = m; while (a > 1) { long long qq = a / mm, t = mm; mm = a % mm; a = t; t = x0; x0 = x1 - qq * x0; x1 = t; } if (x1 < 0) x1 += m; inv = (u64)x1; }
    if ((q * inv) % p != 1) { fprintf(stderr, "inv fail %u\n", p); return 3; }
    invQ[i] = (u32)inv; r64[i] = 64 % p;
    size_t nw = (p + 64 + 63) / 64 + 1; T[i] = calloc(nw, 8); tbytes += nw * 8;
    { u64 dm = 0, step = QS % p;                   /* dm = QS*j mod p */
      for (u64 j = 0; j < (u64)p + 64; j++) {
        int split = (dm != 0) && qr[p - dm];
        if (split) T[i][j >> 6] |= 1ULL << (j & 63);
        dm += step; if (dm >= p) dm -= p;
      } }
  }
  free(qr);
  double tT = now() - t0;
  build_patterns();
  double tP = now() - t0;
  st_n = malloc(sizeof(u64) * 4000000); st_w = malloc(sizeof(long long) * 4000000); st_i = malloc(sizeof(int) * 4000000);
  chi3 = malloc(npr);
  /* table self-check: T against kron on random d */
  { u64 s = 12345; long bad = 0; for (int t = 0; t < 200000; t++) { s ^= s << 13; s ^= s >> 7; s ^= s << 17; u64 d = LO + s % (HI - LO + 1);
      int i = (int)(s % (u64)nmed); u32 p = pr[i17 + i]; u64 j = ((d % p) * invQ[i]) % p;
      int split = (T[i][j >> 6] >> (j & 63)) & 1; if (split != (kron(d, p) == 1)) bad++; }
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
  char path[4096]; snprintf(path, sizeof path, "%s.t4", out); FILE *f4 = fopen(path, "w");
  snprintf(path, sizeof path, "%s.log", out); FILE *flog = fopen(path, "w");
  fprintf(flog, "isieve N=%d LO=%llu HI=%llu medium primes %d (<= %u), T tables %.1f MB, setup %.1f s (tables %.1f, patterns %.1f)\n", N, LO, HI, nmed, pr[i17 + nmed - 1], tbytes / 1e6, now() - t0, tT, tP);
  double tsetup = now();
  fflush(flog);
  u64 lanes = 0, cert1 = 0, certF = 0, t3calls = 0, t3pass = 0, nonfund = 0, t4 = 0, wordsdone = 0, primework = 0, tested = 0, testbad = 0;
  u64 LSEG = QS * 64 * WPC;
  static u64 acc[WPC][NPL], passm[WPC], validm[WPC], ones_[WPC], twos_[WPC], fours_[WPC], eights_[WPC], c16[WPC][NPL];
  static int active[WPC];
  int nseg = 0;
  for (u64 seg_lo = LO; seg_lo <= HI; seg_lo += LSEG) {
    u64 seg_hi = seg_lo + LSEG - 1; if (seg_hi > HI) seg_hi = HI;   /* inclusive */
    u64 Xs = Xof(seg_lo);
    /* per pattern: need and level boundaries */
    static long long need[NPAT]; static u64 Pk[NPAT][NPL + 2]; static int kmax[NPAT];
    for (int pat = 0; pat < NPAT; pat++) {
      long long F = (long long)Fp(pat, Xs); need[pat] = (long long)N + 1 - F;
      int km = 0;
      for (int k = 1; k <= NPL; k++) { u64 m = mfirst(pat, 1ULL << (k - 1)); Pk[pat][k] = m ? Xs / m : 0; if (Pk[pat][k] >= 17) km = k; }
      kmax[pat] = km; Pk[pat][NPL + 1] = 0;
    }
    for (u64 c = 0; c < QS; c++) {
      if (!cls_ok[c]) continue;
      u64 jlo = (seg_lo > c) ? (seg_lo - c + QS - 1) / QS : 0;
      u64 jhi = (seg_hi >= c) ? (seg_hi - c) / QS + 1 : 0;  /* exclusive */
      if (jhi <= jlo) continue;
      u64 nl = jhi - jlo; int nw = (int)((nl + 63) / 64);
      lanes += nl;
      int pat = cls_pat[c];
      if (need[pat] <= 0) { certF += nl; continue; }
      int nact = nw;
      for (int w = 0; w < nw; w++) {
        u64 cnt = (w == nw - 1) ? nl - 64ULL * w : 64; validm[w] = (cnt == 64) ? ~0ULL : ((1ULL << cnt) - 1);
        passm[w] = 0; active[w] = 1; memset(acc[w], 0, sizeof acc[w]);
      }
      u64 needv = (u64)need[pat];
      /* levels descending; primes ascending */
      int ip = 0; /* medium prime index */
      for (int k = kmax[pat]; k >= 1 && nact > 0; k--) {
        u64 plim = Pk[pat][k]; /* primes p <= plim have weight >= 2^k; this level: (Pk[k+1], Pk[k]] */
        int chunk = 256; if (k >= 8) { chunk = 1 << (15 - k); if (chunk < 16) chunk = 16; }
        while (ip < nmed && pr[i17 + ip] <= plim && nact > 0) {
          int iend = ip; while (iend < nmed && pr[i17 + iend] <= plim && iend - ip < chunk) iend++;
          /* reset chunk counters */
          for (int w = 0; w < nw; w++) if (active[w]) { ones_[w] = twos_[w] = fours_[w] = eights_[w] = 0; memset(c16[w], 0, sizeof c16[w]); }
          for (int b = ip; b < iend; b += 16) {
            int nb = iend - b; if (nb > 16) nb = 16;
            u32 off[16], pp[16], rr[16]; const u64 *TT[16];
            for (int t = 0; t < 16; t++) {
              if (t < nb) { int i = b + t; u32 p = pr[i17 + i]; pp[t] = p; rr[t] = r64[i]; TT[t] = T[i];
                off[t] = (u32)(((jlo % p) + (u64)(c % p) * invQ[i]) % p); }
              else { pp[t] = 1; rr[t] = 0; TT[t] = 0; off[t] = 0; }
            }
            for (int w = 0; w < nw; w++) {
              u64 x[16];
              if (active[w]) {
                for (int t = 0; t < 16; t++) {
                  if (TT[t]) { u32 o = off[t]; const u64 *Tp = TT[t]; u32 ix = o >> 6, sh = o & 63;
                    x[t] = (Tp[ix] >> sh) | ((Tp[ix + 1] << 1) << (63 - sh)); }
                  else x[t] = 0;
                }
                u64 o_ = ones_[w], t_ = twos_[w], f_ = fours_[w], e_ = eights_[w];
                #define CSA(h, l, a, b_, cc) { u64 u_ = (a) ^ (b_); h = ((a) & (b_)) | (u_ & (cc)); l = u_ ^ (cc); }
                u64 tA, tB, fA, fB, eA, eB, s16;
                CSA(tA, o_, o_, x[0], x[1]); CSA(tB, o_, o_, x[2], x[3]); CSA(fA, t_, t_, tA, tB);
                CSA(tA, o_, o_, x[4], x[5]); CSA(tB, o_, o_, x[6], x[7]); CSA(fB, t_, t_, tA, tB);
                CSA(eA, f_, f_, fA, fB);
                CSA(tA, o_, o_, x[8], x[9]); CSA(tB, o_, o_, x[10], x[11]); CSA(fA, t_, t_, tA, tB);
                CSA(tA, o_, o_, x[12], x[13]); CSA(tB, o_, o_, x[14], x[15]); CSA(fB, t_, t_, tA, tB);
                CSA(eB, f_, f_, fA, fB);
                CSA(s16, e_, e_, eA, eB);
                ones_[w] = o_; twos_[w] = t_; fours_[w] = f_; eights_[w] = e_;
                { u64 carry = s16; for (int q = 0; q < NPL && carry; q++) { u64 tt = c16[w][q] & carry; c16[w][q] ^= carry; carry = tt; } }
                primework += (u64)nb;
              }
              for (int t = 0; t < nb; t++) { u32 o = off[t] + rr[t]; if (o >= pp[t]) o -= pp[t]; off[t] = o; }
            }
          }
          /* add chunk counts << k into acc, check */
          for (int w = 0; w < nw; w++) {
            if (!active[w]) continue;
            u64 vec[4] = {ones_[w], twos_[w], fours_[w], eights_[w]};
            for (int q = 0; q < 4 + NPL; q++) {
              u64 v = (q < 4) ? vec[q] : c16[w][q - 4];
              int pos = k + q; if (!v) continue; if (pos >= NPL) { passm[w] |= v; continue; } /* weight >= 2^16 > need: certified */
              u64 carry = v; for (int z = pos; z < NPL && carry; z++) { u64 tt = acc[w][z] & carry; acc[w][z] ^= carry; carry = tt; }
              passm[w] |= carry; /* overflow past 2^16 means >= 2^16 > need */
            }
            /* acc >= needv ? bit-sliced compare: gt = acc > needv-1 */
            u64 thr = needv - 1, gt = 0, eq = ~0ULL;
            for (int z = NPL - 1; z >= 0; z--) {
              u64 a = acc[w][z];
              if ((thr >> z) & 1) { eq &= a; } else { gt |= eq & a; eq &= ~a; }
            }
            passm[w] |= gt;
            if (((passm[w] | ~validm[w]) == ~0ULL)) { active[w] = 0; nact--; wordsdone++; }
          }
          ip = iend;
        }
        if (ip < nmed && pr[i17 + ip] <= plim) {} /* all words done */
      }
      /* survivors, and test-mode checks */
      for (int w = 0; w < nw; w++) {
        u64 sv = validm[w] & ~passm[w];
        cert1 += (u64)__builtin_popcountll(validm[w] & passm[w]);
        if (testK) {
          for (int t = 0; t < 64; t++) {
            if (!((validm[w] & passm[w]) >> t & 1)) continue;
            u64 d = c + QS * (jlo + 64ULL * w + t);
            if ((d * 2654435761ULL) % (u64)testK) continue;
            /* lane value (may be wrapped/overflowed only if passed earlier: compare against need only) */
            long long M = tier3_M(d, (long long)N + 1);
            tested++; if (M < (long long)N + 1) { testbad++; fprintf(flog, "TEST FAIL d=%llu M=%lld\n", d, M); }
          }
        }
        while (sv) {
          int t = __builtin_ctzll(sv); sv &= sv - 1;
          u64 d = c + QS * (jlo + 64ULL * w + t);
          t3calls++;
          long long M = tier3_M(d, (long long)N + 1);
          if (M >= (long long)N + 1) { t3pass++; continue; }
          if (!fundamental(d)) { nonfund++; continue; }
          t4++; fprintf(f4, "%llu %lld\n", d, M);
        }
      }
    }
    nseg++;
    if (nseg % 8 == 0 || seg_hi == HI) {
      fprintf(flog, "seg %d [%llu,%llu] X=%llu lanes %llu certF %llu cert1 %llu t3 %llu t3pass %llu nonfund %llu t4 %llu primework/lane %.1f elapsed %.1f s (sieve %.1f s)\n",
              nseg, seg_lo, seg_hi, Xs, lanes, certF, cert1, t3calls, t3pass, nonfund, t4, lanes ? (double)primework * 64 / lanes : 0, now() - t0, now() - tsetup);
      fflush(flog); fflush(f4);
    }
    if (seg_hi == HI) break;
  }
  fprintf(flog, "DONE N=%d [%llu,%llu]: lanes %llu = certF %llu + cert1 %llu + t3 %llu (t3pass %llu, nonfund %llu, t4 %llu); tested %llu bad %llu; %.1f s\n",
          N, LO, HI, lanes, certF, cert1, t3calls, t3pass, nonfund, t4, tested, testbad, now() - t0);
  if (lanes != certF + cert1 + t3calls) fprintf(flog, "LANE COUNT MISMATCH\n");
  fclose(flog); fclose(f4);
  return 0;
}
