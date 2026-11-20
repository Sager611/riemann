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
      for (int j = 0; j < M1; j++) { uint64_t ph = (key << 40) - EPSU + (uint64_t)j*g[0].d1; double l2 = g[0].mx - g[0].T[0][ph >> 54]; if (l2 < Gneed) { jl[nl] = j; ll[nl] = l2; nl++; } } }
    for (int q = 0; q < nl; q++) {
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
