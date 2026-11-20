// dbnrec3 bb_rec3.c (= dbnrec2 bb_rec2.c + exact two-group anchor prefilter ANY2[], performance only, see build_any2):
// dbnrec2 bb_rec2.c (= dbnrec bb_rec.c + exact anchor prefilter ANY[], performance only, see build_any): 
// dbnrec bb_rec.c (= dbnhead bb.c + float margin FM = 1e-9 on every pruning test: Gneed = U0C - thr + FM):
// dbnhead bb.c: rigorous branch-and-bound upper bound of sup over x' in [xs, xs + K*D] of
//   U(x') = sum_g g_g(theta_g(x')),  g_g(th) = sum_h max(lo_h cos(k_h th), hi_h cos(k_h th)),  theta_g = (x'/2) log b_g,
// i.e. of the finite head sum_n c_n cos(tau log n) over all coefficient values c_n in [lo_n, hi_n] (prep.py).
// Phases: exact uint64 turns (increments rounded to 1 unit; total drift < 2^-34 turns, covered by EPSU = 2^-24 turns).
// Cells: anchors x_k = xs + k D (D = 4pi/log 2), level-1 cells D/M1, level-2 D/(M1 M2), level-3 D/(M1 M2 M3), then Lipschitz bisection.
// Tables: T_L[g][b] >= max of g_g over phases [b/B, (b+1)/B + arc_L + 4 EPSU] (sampled, midpoint + Lipschitz cell bound).
// usage: bb data m1 m2 [Kmax] [thr] [hmin] [dumpfile]
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>
#define B 1024
#define SPB 8
#define NS (B*SPB)
#define EPSU ((uint64_t)1 << 40)
typedef struct { int b, nh; int k[16]; double lo[16], hi[16], ac[16]; uint64_t p0, dK, d1, d2, d3; double om, lip, mx; double *T[3]; } G;
static G *g; static int NG, M1, M2, M3, m1, m2; static FILE *dumpf;
static double D, thr, Gneed, U0C, LIPX, HMIN, h1, h2, h3, PHERR;
static long long n1, s1, s2, s3, nfin, nunres, curk; /* dbnrec: curk = current anchor (exact), dumped with the in-anchor offset */ static double maxUB3 = -1e300, maxUBun = -1e300, maxQ = -1e300, maxQx = 0, maxUBunx = 0;
static double gev(const G *q, double th) { double s = 0; for (int h = 0; h < q->nh; h++) { double c = cos(q->k[h]*th); s += c >= 0 ? q->hi[h]*c : q->lo[h]*c; } return s; }
static double gact(const G *q, double th) { double s = 0; for (int h = 0; h < q->nh; h++) s += q->ac[h]*cos(q->k[h]*th); return s; }
static double U[NS], V[NS + 1], BL[B];
// dbnrec2 prefilter: P1[q][b] = (the original level-1 loop prunes candidate q at gi = 1 when group 1's bin is b), evaluated with the same
// float operations (loss = ll[q]; loss += mx - T; loss >= Gneed). For the 16-bit top f of base = p0 + k dK - EPSU (group 1) the true bin of
// base + j d1 is lo >> 54 or hi >> 54 (lo = (f << 48) + j d1, hi = lo + 2^48 - 1, mod 2^64). ANY[f] = 0 iff every candidate is pruned at gi = 1
// for both bins; then the anchor is skipped with n1 += nl, exactly what the original loop does for it (each candidate: n1++, break at gi = 1, continue).
static unsigned char P1[64][B], ANY[1 << 16];
// dbnrec3 two-group prefilter: NP[q][b1] = bitset of the group-2 bins b2 for which the original level-1 loop does NOT prune candidate q at gi = 1 or
// gi = 2 when group 1's bin is b1 and group 2's bin is b2, with the same float operations (loss = ll[q]; loss += mx1 - T1[b1]; loss >= Gneed;
// loss += mx2 - T2[b2]; loss >= Gneed). For the 11-bit tops a, c of base_g = p0_g + k dK_g - EPSU (g = 1, 2) the true bin of base_g + j d1_g is
// lo >> 54 or hi >> 54 (lo = (a << 53) + j d1_g, hi = lo + 2^53 - 1, mod 2^64; span 2^53 - 1 < 2^54). ANY2 bit (a, c) = 0 iff every candidate is pruned
// at gi = 1 or 2 for all four bin pairs; then the anchor is skipped with n1 += nl, exactly what the original loop does (n1++, break at gi <= 2, continue).
static uint64_t NP[64][B][B/64]; static unsigned char ANY2[1 << 19];
static void build_any2(int nl, const int *jl, const double *ll) {
  for (int q = 0; q < nl; q++) for (int b1 = 0; b1 < B; b1++) { double loss = ll[q]; loss += g[1].mx - g[1].T[0][b1]; memset(NP[q][b1], 0, sizeof(NP[q][b1]));
    if (loss >= Gneed) continue;
    for (int b2 = 0; b2 < B; b2++) { double l2 = loss; l2 += g[2].mx - g[2].T[0][b2]; if (!(l2 >= Gneed)) NP[q][b1][b2 >> 6] |= (uint64_t)1 << (b2 & 63); } }
  memset(ANY2, 0, sizeof(ANY2));
  for (uint64_t a = 0; a < 2048; a++) for (int q = 0; q < nl; q++) {
    uint64_t lo1 = (a << 53) + (uint64_t)jl[q]*g[1].d1, hi1 = lo1 + (((uint64_t)1 << 53) - 1); uint64_t u[B/64];
    for (int w = 0; w < B/64; w++) u[w] = NP[q][lo1 >> 54][w] | NP[q][hi1 >> 54][w];
    for (uint64_t c = 0; c < 2048; c++) { uint64_t lo2 = (c << 53) + (uint64_t)jl[q]*g[2].d1, hi2 = lo2 + (((uint64_t)1 << 53) - 1);
      uint64_t b2a = lo2 >> 54, b2b = hi2 >> 54;
      if (((u[b2a >> 6] >> (b2a & 63)) | (u[b2b >> 6] >> (b2b & 63))) & 1) ANY2[(a << 8) | (c >> 3)] |= (unsigned char)(1 << (c & 7)); } }
}
static void build_any(int nl, const int *jl, const double *ll) {
  for (int q = 0; q < nl; q++) for (int b = 0; b < B; b++) { double loss = ll[q]; loss += g[1].mx - g[1].T[0][b]; P1[q][b] = loss >= Gneed; }
  for (uint64_t f = 0; f < (1 << 16); f++) { int any = 0;
    for (int q = 0; q < nl && !any; q++) { uint64_t lo = (f << 48) + (uint64_t)jl[q]*g[1].d1, hi = lo + (((uint64_t)1 << 48) - 1);
      if (!(P1[q][lo >> 54] && P1[q][hi >> 54])) any = 1; }
    ANY[f] = (unsigned char)any; }
}
static void prepg(G *q) {
  double dth = 2*M_PI/NS; q->lip = 0; for (int h = 0; h < q->nh; h++) q->lip += q->k[h]*fmax(fabs(q->lo[h]), fabs(q->hi[h]));
  for (int i = 0; i <= NS; i++) V[i] = gev(q, i*dth);
  q->mx = -1e300;
  for (int i = 0; i < NS; i++) { U[i] = 0.5*(V[i] + V[i + 1]) + 0.5*q->lip*dth*1.0000001 + 1e-15; if (U[i] > q->mx) q->mx = U[i]; }
  for (int b = 0; b < B; b++) { double m = -1e300; for (int i = b*SPB; i < (b + 1)*SPB; i++) if (U[i] > m) m = U[i]; BL[b] = m; }
}
static void build(G *q, double arc, double *T) {   // arc in turns; uses BL of q (call prepg first)
  double span = 1.0/B + arc + 4.0*ldexp(1.0, -24);
  int nb = (int)ceil(span*B) + 1;
  if (nb >= B) { for (int b = 0; b < B; b++) T[b] = q->mx; return; }
  static int dq[3*B]; int hd = 0, tl = 0;
  for (int i = 0; i < B + nb - 1; i++) {
    double x = BL[i % B];
    while (tl > hd && BL[dq[tl - 1] % B] <= x) tl--;
    dq[tl++] = i;
    if (dq[hd] <= i - nb) hd++;
    if (i >= nb - 1) T[i - nb + 1] = BL[dq[hd] % B];
  }
}
static void fin(const double *th, double len, double xoff, int depth) {
  double ub = 0, half = 0.5*len;
  for (int i = 0; i < NG; i++) ub += gev(&g[i], th[i] + g[i].om*half);
  ub += LIPX*half + 1e-9 + PHERR; nfin++;   // PHERR: start-phase error < 2^-30 turn (review fix, 08:25)
  if (depth == 0) { double qa = 0; for (int i = 0; i < NG; i++) qa += gact(&g[i], th[i] + g[i].om*half); if (qa > maxQ) { maxQ = qa; maxQx = (double)curk*D + xoff + half; } }
  if (ub <= thr) return;
  if (len < HMIN) { nunres++; if (ub > maxUBun) { maxUBun = ub; maxUBunx = (double)curk*D + xoff + half; }
    if (dumpf) fprintf(dumpf, "%lld %.15f %.6e %.9f\n", curk, xoff + half, len, ub);   // dbnrec: anchor k, in-anchor centre offset, length, ub
    if (nunres <= 20) { double qa = 0; for (int i = 0; i < NG; i++) qa += gact(&g[i], th[i] + g[i].om*half); printf("  unresolved x'-xs = %.6f ub %.9f qact %.9f\n", (double)curk*D + xoff + half, ub, qa); }
    if (nunres > 200000) { printf("ABORT: too many unresolved\n"); fflush(stdout); exit(3); }
    return; }
  double *t2 = malloc(sizeof(double)*NG);
  fin(th, half, xoff, depth + 1);
  for (int i = 0; i < NG; i++) t2[i] = th[i] + g[i].om*half;
  fin(t2, half, xoff + half, depth + 1);
  free(t2);
}
static void level3(uint64_t k, int j, int j2) {
  for (int j3 = 0; j3 < M3; j3++) {
    double loss = 0; int gi;
    for (gi = 0; gi < NG; gi++) { uint64_t ph = g[gi].p0 + k*g[gi].dK + (uint64_t)j*g[gi].d1 + (uint64_t)j2*g[gi].d2 + (uint64_t)j3*g[gi].d3 - EPSU; loss += g[gi].mx - g[gi].T[2][ph >> 54]; if (loss >= Gneed) break; }
    if (gi < NG) continue;
    s3++; if (U0C - loss > maxUB3) maxUB3 = U0C - loss;
    double *th = malloc(sizeof(double)*NG);
    for (int i = 0; i < NG; i++) { uint64_t ph = g[i].p0 + k*g[i].dK + (uint64_t)j*g[i].d1 + (uint64_t)j2*g[i].d2 + (uint64_t)j3*g[i].d3; th[i] = ldexp((double)ph, -64)*2*M_PI; }
    fin(th, h3, j*h1 + j2*h2 + j3*h3, 0);   // dbnrec: in-anchor offset; anchor in curk
    free(th);
  }
}
static void level2(uint64_t k, int j) {
  for (int j2 = 0; j2 < M2; j2++) {
    double loss = 0; int gi;
    for (gi = 0; gi < m2; gi++) { uint64_t ph = g[gi].p0 + k*g[gi].dK + (uint64_t)j*g[gi].d1 + (uint64_t)j2*g[gi].d2 - EPSU; loss += g[gi].mx - g[gi].T[1][ph >> 54]; if (loss >= Gneed) break; }
    if (gi < m2) continue;
    s2++; level3(k, j, j2);
  }
}
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r"); m1 = atoi(argv[2]); m2 = atoi(argv[3]);
  long long K; double mxpp, tri, comp;
  if (fscanf(f, "%d %lld %lf %d %d %d %lf %lf %lf %lf", &NG, &K, &D, &M1, &M2, &M3, &thr, &mxpp, &tri, &comp) != 10) return 1;
  if (M1 > 64) { printf("ABORT: M1 > 64 (jl, ll, P1 capacity)\n"); return 4; }   // dbnrec2: nl <= M1 <= 64
  if (argc > 4 && atoll(argv[4]) > 0) K = atoll(argv[4]);
  if (argc > 5 && atof(argv[5]) > 0) thr = atof(argv[5]);
  HMIN = argc > 6 ? atof(argv[6]) : 1e-6;
  dumpf = argc > 7 ? fopen(argv[7], "w") : NULL;   // optional: every unresolved cell (x' - xs, ub)
  h1 = D/M1; h2 = h1/M2; h3 = h2/M3;
  g = calloc(NG, sizeof(G));
  for (int i = 0; i < NG; i++) {
    unsigned long long a, b2, c, d, e;
    if (fscanf(f, "%d %d %llu %llu %llu %llu %llu %lf", &g[i].b, &g[i].nh, &a, &b2, &c, &d, &e, &g[i].om) != 8) return 2;
    g[i].p0 = a; g[i].dK = b2; g[i].d1 = c; g[i].d2 = d; g[i].d3 = e;
    for (int h = 0; h < g[i].nh; h++) if (fscanf(f, "%d %lf %lf %lf", &g[i].k[h], &g[i].lo[h], &g[i].hi[h], &g[i].ac[h]) != 4) return 2;
  }
  if (m1 > NG) m1 = NG; if (m2 > NG) m2 = NG;
  clock_t c0 = clock(); U0C = 0; LIPX = 0;
  for (int i = 0; i < NG; i++) {
    prepg(&g[i]); U0C += g[i].mx; LIPX += g[i].lip*g[i].om;
    double arcs[3] = {h1*g[i].om/(2*M_PI), h2*g[i].om/(2*M_PI), h3*g[i].om/(2*M_PI)};
    for (int L = 0; L < 3; L++) { if ((L == 0 && i >= m1) || (L == 1 && i >= m2)) { g[i].T[L] = NULL; continue; } g[i].T[L] = malloc(sizeof(double)*B); build(&g[i], arcs[L], g[i].T[L]); }
  }
  LIPX *= 1.0000001;
  { double sl = 0; for (int i = 0; i < NG; i++) sl += g[i].lip; PHERR = 2*M_PI*ldexp(1.0, -30)*sl*1.01; }
  Gneed = U0C - thr + 1e-9;   // dbnrec: prune only if bound <= thr - 1e-9 (covers float summation error < 1e-11)
  printf("NG %d K %lld D %.12f h1 %.6f h2 %.6f h3 %.6f | U0C %.9f (mx_pp %.9f comp %.9f tri %.9f) thr %.9f Gneed %.9f LIPX %.4f | tables %.1f s\n", NG, K, D, h1, h2, h3, U0C, mxpp, comp, tri, thr, Gneed, LIPX, (double)(clock() - c0)/CLOCKS_PER_SEC);
  fflush(stdout);
  int jl[64]; double ll[64]; int nl = 0; uint64_t lastkey = ~(uint64_t)0;
  c0 = clock();
  for (long long kk = 0; kk < K; kk++) {
    uint64_t k = (uint64_t)kk; curk = kk;
    uint64_t key = (g[0].p0 + k*g[0].dK) >> 40;
    if (key != lastkey) { nl = 0; lastkey = key;
      for (int j = 0; j < M1; j++) { uint64_t ph = (key << 40) - EPSU + (uint64_t)j*g[0].d1; double l2 = g[0].mx - g[0].T[0][ph >> 54]; if (l2 < Gneed) { jl[nl] = j; ll[nl] = l2; nl++; } }
      if (m1 > 1) build_any(nl, jl, ll);
      if (m1 > 2) build_any2(nl, jl, ll); }
    int skip = 0;   // dbnrec2 prefilter
    if (m1 > 1) { uint64_t base = g[1].p0 + k*g[1].dK - EPSU; if (!ANY[base >> 48]) { skip = 1; n1 += nl; }
      else if (m1 > 2) { uint64_t a = base >> 53, c = (g[2].p0 + k*g[2].dK - EPSU) >> 53; if (!((ANY2[(a << 8) | (c >> 3)] >> (c & 7)) & 1)) { skip = 1; n1 += nl; } } }   // dbnrec3 two-group prefilter
    if (!skip) for (int q = 0; q < nl; q++) {
      int j = jl[q]; double loss = ll[q]; int gi; n1++;
      for (gi = 1; gi < m1; gi++) { uint64_t ph = g[gi].p0 + k*g[gi].dK + (uint64_t)j*g[gi].d1 - EPSU; loss += g[gi].mx - g[gi].T[0][ph >> 54]; if (loss >= Gneed) break; }
      if (gi < m1) continue;
      s1++; level2(k, j);
    }
    if ((kk & ((1LL << 27) - 1)) == 0 && kk) { printf("  k %lld (%.3f) n1 %lld s1 %lld s2 %lld s3 %lld fin %lld unres %lld maxUB3 %.6f maxQ %.6f | %.0f s\n", kk, (double)kk/K, n1, s1, s2, s3, nfin, nunres, maxUB3, maxQ, (double)(clock() - c0)/CLOCKS_PER_SEC); fflush(stdout); }
  }
  printf("DONE K %lld window %.6g | n1 %lld s1 %lld s2 %lld s3 %lld fin %lld unres %lld | maxUB3 %.9f maxUBunres %.9f at %.4f | maxQact %.9f at x'-xs = %.4f | thr %.9f | %.0f s\n",
         K, (double)K*D, n1, s1, s2, s3, nfin, nunres, maxUB3, maxUBun, maxUBunx, maxQ, maxQx, thr, (double)(clock() - c0)/CLOCKS_PER_SEC);
  printf("CERT sup <= %.9f\n", nunres ? fmax(thr, maxUBun) : thr);
  if (dumpf) fclose(dumpf);
  return 0;
}
