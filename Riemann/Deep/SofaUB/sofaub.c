/* sofaub.c -- rigorous branch-and-bound / verifier for V_ii(u) (3D snaky-corridor sofa, upper bound).
 *
 * Problem (arXiv 2610.14050, Prop. "event-forced tube family"): V*(K_M) <= max{u, V_ii(u)} for every u in (sqrt2, M),
 * V_ii(u) = sup vol(T ^ (S_X u Stub(c))) over unit-square tubes T meeting W_1 and W_c with |w.yhat| >= mu0(u).
 * Slice form (NOTE.md): a tube is Sigma(y) = X + (y-ya)k + G[-1/2,1/2]^2, G = A_k R(t), A_k = I + k k^T/(1+s),
 * s = sqrt(1+|k|^2), R(t) = rotation with cos = (1-t^2)/(1+t^2), sin = 2t/(1+t^2), t in [0,1];
 * V = int_0^c area(Sigma(y) ^ R_y) dy, R_y = (-inf,1]x[0,1] (y<=1), [0,1]^2 (1<y<=c).
 * Box bound: Sigma_p(y) is contained in D(y) = {p : |cross(p-X(y),g2)| <= A(y), |cross(g1,p-X(y))| <= B(y)},
 * X(y) = Xc + (y-ya)kc, A(y) = det/2 + PA0 + PA1|y-ya| (B likewise), for every p in the box; so
 * V(p) <= int_0^c area(D(y) ^ R_y) dy, computed exactly (piecewise quadratic) with interval arithmetic.
 * Arithmetic: IEEE double, round to nearest; every interval operation is widened outward by |x|*2^-50 + 2^-1074
 * (>= 4 ulp), so intervals are rigorous enclosures. Only + - * / sqrt are used (no libm transcendental).
 * Build: cc -O2 -ffp-contract=off -fno-fast-math -o sofaub sofaub.c -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

typedef struct { double lo, hi; } I;
static int g_err = 0;   /* set when an enclosure could not be computed (degenerate division etc.) */

static inline double DN(double x) { return x - (fabs(x)*0x1p-50 + 0x1p-1074); }
static inline double UP(double x) { return x + (fabs(x)*0x1p-50 + 0x1p-1074); }
static inline I mk(double a) { I r = {a, a}; return r; }
static inline I mk2(double a, double b) { I r = {a, b}; return r; }
static inline I add(I a, I b) { I r = {DN(a.lo + b.lo), UP(a.hi + b.hi)}; return r; }
static inline I sub(I a, I b) { I r = {DN(a.lo - b.hi), UP(a.hi - b.lo)}; return r; }
static inline I neg(I a) { I r = {-a.hi, -a.lo}; return r; }
static inline I mul(I a, I b) {
  double p0 = a.lo*b.lo, p1 = a.lo*b.hi, p2 = a.hi*b.lo, p3 = a.hi*b.hi;
  double lo = fmin(fmin(DN(p0), DN(p1)), fmin(DN(p2), DN(p3)));
  double hi = fmax(fmax(UP(p0), UP(p1)), fmax(UP(p2), UP(p3)));
  I r = {lo, hi}; return r;
}
static inline I divi(I a, I b) {
  if (!(b.lo > 0 || b.hi < 0)) { g_err = 1; I r = {-INFINITY, INFINITY}; return r; }
  double p0 = a.lo/b.lo, p1 = a.lo/b.hi, p2 = a.hi/b.lo, p3 = a.hi/b.hi;
  double lo = fmin(fmin(DN(p0), DN(p1)), fmin(DN(p2), DN(p3)));
  double hi = fmax(fmax(UP(p0), UP(p1)), fmax(UP(p2), UP(p3)));
  I r = {lo, hi}; return r;
}
static inline I sqri(I a) {
  I r;
  if (a.lo >= 0) { r.lo = DN(a.lo*a.lo); r.hi = UP(a.hi*a.hi); }
  else if (a.hi <= 0) { r.lo = DN(a.hi*a.hi); r.hi = UP(a.lo*a.lo); }
  else { r.lo = 0; r.hi = UP(fmax(a.lo*a.lo, a.hi*a.hi)); }
  if (r.lo < 0) r.lo = 0;
  return r;
}
static inline I sqrti(I a) { I r; r.lo = DN(sqrt(fmax(a.lo, 0.0))); if (r.lo < 0) r.lo = 0; r.hi = UP(sqrt(fmax(a.hi, 0.0))); return r; }
static inline I maxi(I a, I b) { I r = {fmax(a.lo, b.lo), fmax(a.hi, b.hi)}; return r; }
static inline I mini(I a, I b) { I r = {fmin(a.lo, b.lo), fmin(a.hi, b.hi)}; return r; }
static inline I absi(I a) { I r; if (a.lo >= 0) r = a; else if (a.hi <= 0) r = neg(a); else { r.lo = 0; r.hi = fmax(-a.lo, a.hi); } return r; }
static inline I scal(double s, I a) { return mul(mk(s), a); }
static inline I half(I a) { I r = {a.lo*0.5, a.hi*0.5}; return r; }   /* exact (no underflow at our scales) */
static inline double upv(I a) { return a.hi; }

/* ---------------- global problem data ---------------- */
static double U, C, YA, TARGET;      /* u, c = 1+u (exact), anchor height, target */
static I K2;                          /* K^2 = 8u^2/(u^2-2)^2 enclosure */
static double KHI;                    /* upper bound of K */

/* ---------------- box / slice data ---------------- */
enum { DKX = 0, DKZ, DT, DX1, DX2, ND };
typedef struct { double lo[ND], hi[ND]; } Box;

typedef struct {
  double xc[2], kc[2];   /* center X(ya), center slope (exact floats) */
  I g1x, g1z, g2x, g2z;   /* edge vectors of the center slice */
  I det;                  /* cross(g1,g2) = s(kc) */
  double PA0, PA1, PB0, PB1; /* push-outs (upper bounds) */
  double kmax2, kmin2;    /* bounds on |k|^2 over the box */
  double smax;            /* upper bound of sqrt(1+|k|^2) over the box */
  double cost[ND];        /* heuristic split costs */
} Slice;

static inline I cross2(I ax, I az, I bx, I bz) { return sub(mul(ax, bz), mul(az, bx)); }

static void center_and_rad(double lo, double hi, double *c, double *r) {
  double m = 0.5*(lo + hi);
  if (m < lo) m = lo; if (m > hi) m = hi;
  *c = m;
  *r = UP(fmax(m - lo, hi - m));  /* m - lo, hi - m computed in RN then widened: >= true half-widths */
}

static void setup(const Box *b, Slice *S) {
  double kc[2], rk[2], tc, rt, xc[2], rx[2];
  center_and_rad(b->lo[DKX], b->hi[DKX], &kc[0], &rk[0]);
  center_and_rad(b->lo[DKZ], b->hi[DKZ], &kc[1], &rk[1]);
  center_and_rad(b->lo[DT], b->hi[DT], &tc, &rt);
  center_and_rad(b->lo[DX1], b->hi[DX1], &xc[0], &rx[0]);
  center_and_rad(b->lo[DX2], b->hi[DX2], &xc[1], &rx[1]);
  S->xc[0] = xc[0]; S->xc[1] = xc[1]; S->kc[0] = kc[0]; S->kc[1] = kc[1];
  I kx = mk(kc[0]), kz = mk(kc[1]);
  I s = sqrti(add(add(mk(1.0), sqri(kx)), sqri(kz)));
  I f = divi(mk(1.0), add(mk(1.0), s));
  I A11 = add(mk(1.0), mul(f, sqri(kx))), A12 = mul(f, mul(kx, kz)), A22 = add(mk(1.0), mul(f, sqri(kz)));
  I t = mk(tc), t2 = sqri(t), den = add(mk(1.0), t2);
  I ct = divi(sub(mk(1.0), t2), den), st = divi(scal(2.0, t), den);
  S->g1x = add(mul(A11, ct), mul(A12, st)); S->g1z = add(mul(A12, ct), mul(A22, st));
  S->g2x = sub(mul(A12, ct), mul(A11, st)); S->g2z = sub(mul(A22, ct), mul(A12, st));
  /* det(A R) = det(A) = s exactly; intersect with the computed cross product */
  I det = cross2(S->g1x, S->g1z, S->g2x, S->g2z);
  det.lo = fmax(det.lo, s.lo); det.hi = fmin(det.hi, s.hi);
  if (det.lo > det.hi) { det = s; }
  S->det = det;
  /* |k|^2 bounds over the box */
  double ax_hi = fmax(fabs(b->lo[DKX]), fabs(b->hi[DKX])), az_hi = fmax(fabs(b->lo[DKZ]), fabs(b->hi[DKZ]));
  double ax_lo = (b->lo[DKX] > 0) ? b->lo[DKX] : ((b->hi[DKX] < 0) ? -b->hi[DKX] : 0.0);
  double az_lo = (b->lo[DKZ] > 0) ? b->lo[DKZ] : ((b->hi[DKZ] < 0) ? -b->hi[DKZ] : 0.0);
  I kmax2 = add(sqri(mk(ax_hi)), sqri(mk(az_hi))), kmin2 = add(sqri(mk(ax_lo)), sqri(mk(az_lo)));
  S->kmax2 = kmax2.hi; S->kmin2 = kmin2.lo;
  I kmax = sqrti(mk(kmax2.hi));
  I smaxI = sqrti(add(mk(1.0), mk(kmax2.hi)));
  S->smax = smaxI.hi;
  /* derivative bounds of A over the box */
  I kmI = mk(kmax.hi);
  I sK = sqrti(add(mk(1.0), sqri(kmI)));
  I LA = add(divi(scal(2.0, kmI), add(mk(1.0), sK)),
             divi(mul(sqri(kmI), kmI), mul(sK, sqri(add(mk(1.0), sK)))));
  I smin = sqrti(add(mk(1.0), mk(kmin2.lo)));
  I fmaxI = divi(mk(1.0), add(mk(1.0), smin));
  I km2 = mk(kmax2.hi);
  I D2 = add(add(divi(scal(4.0, mul(sqri(fmaxI), km2)), smin), scal(2.0, fmaxI)),
             add(divi(scal(2.0, mul(mul(sqri(fmaxI), fmaxI), sqri(km2))), sqri(smin)), divi(mul(sqri(fmaxI), km2), smin)));
  I dk = sqrti(add(sqri(mk(rk[0])), sqri(mk(rk[1]))));
  /* rotation half-angle bound: |psi(t)-psi(tc)| <= 2 rt/(1+tmin^2) */
  double tmin = fmax(0.0, DN(tc - rt));
  I delta = divi(scal(2.0, mk(rt)), add(mk(1.0), sqri(mk(tmin))));
  I dl = mk(delta.hi);
  I ROT = half(mul(det, add(dl, half(sqri(dl)))));
  /* dA_a at kc: F (e_a k^T + k e_a^T) + dF_a k k^T, dF_a = -F^2 k_a / s */
  I dA[2][2][2];
  for (int a = 0; a < 2; a++) {
    I ka = a ? kz : kx;
    I dF = neg(divi(mul(sqri(f), ka), s));
    I kk[2] = {kx, kz};
    for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) {
      I e = mk(0.0);
      if (i == a) e = add(e, kk[j]);
      if (j == a) e = add(e, kk[i]);
      dA[a][i][j] = add(mul(f, e), mul(dF, mul(kk[i], kk[j])));
    }
  }
  I v[2][2] = {{ct, st}, {neg(st), ct}};
  I g1n = sqrti(add(sqri(S->g1x), sqri(S->g1z))), g2n = sqrti(add(sqri(S->g2x), sqri(S->g2z)));
  I rem2 = half(mul(D2, sqri(dk)));        /* (1/2) D2 |dk|^2 */
  I remr = mul(mul(LA, dk), dl);           /* L_A |dk| delta */
  I shA = mk(0.0), shB = mk(0.0);
  double cA[2] = {0, 0}, cB[2] = {0, 0};  /* per-k-dim shape cost (heuristic) */
  for (int j = 0; j < 2; j++) {
    I tA = mk(0.0), tB = mk(0.0);
    for (int a = 0; a < 2; a++) {
      I wx = add(mul(dA[a][0][0], v[j][0]), mul(dA[a][0][1], v[j][1]));
      I wz = add(mul(dA[a][1][0], v[j][0]), mul(dA[a][1][1], v[j][1]));
      I cA_ = mul(mk(rk[a]), absi(cross2(wx, wz, S->g2x, S->g2z)));
      I cB_ = mul(mk(rk[a]), absi(cross2(S->g1x, S->g1z, wx, wz)));
      tA = add(tA, cA_); tB = add(tB, cB_);
      cA[a] += 0.5*cA_.hi; cB[a] += 0.5*cB_.hi;
    }
    tA = add(tA, add(mul(rem2, g2n), mul(remr, g2n)));
    tB = add(tB, add(mul(rem2, g1n), mul(remr, g1n)));
    shA = add(shA, half(tA)); shB = add(shB, half(tB));
  }
  I rX1 = mk(rx[0]), rX2 = mk(rx[1]), rK1 = mk(rk[0]), rK2 = mk(rk[1]);
  I PXA = add(mul(rX1, absi(S->g2z)), mul(rX2, absi(S->g2x)));
  I PXB = add(mul(rX1, absi(S->g1z)), mul(rX2, absi(S->g1x)));
  I PKA = add(mul(rK1, absi(S->g2z)), mul(rK2, absi(S->g2x)));
  I PKB = add(mul(rK1, absi(S->g1z)), mul(rK2, absi(S->g1x)));
  S->PA0 = add(add(PXA, ROT), shA).hi;
  S->PB0 = add(add(PXB, ROT), shB).hi;
  S->PA1 = PKA.hi; S->PB1 = PKB.hi;
  /* heuristic split costs (deterministic; soundness does not depend on them) */
  double ybar = (YA*YA + (C - YA)*(C - YA))/(2.0*C);
  S->cost[DX1] = rx[0]*(fabs(S->g2z.hi) + fabs(S->g1z.hi));
  S->cost[DX2] = rx[1]*(fabs(S->g2x.hi) + fabs(S->g1x.hi));
  S->cost[DKX] = rk[0]*(fabs(S->g2z.hi) + fabs(S->g1z.hi))*ybar + cA[0] + cB[0];
  S->cost[DKZ] = rk[1]*(fabs(S->g2x.hi) + fabs(S->g1x.hi))*ybar + cA[1] + cB[1];
  S->cost[DT] = 2.0*ROT.hi;
}

/* ---------------- slice geometry at height y ---------------- */
typedef struct { I Xx, Xz, al, be, A, B; } SliceY;

static void slice_y(const Slice *S, double y, SliceY *q) {
  I dy = sub(mk(y), mk(YA));
  I ady = absi(dy);
  q->Xx = add(mk(S->xc[0]), mul(dy, mk(S->kc[0])));
  q->Xz = add(mk(S->xc[1]), mul(dy, mk(S->kc[1])));
  I pa = add(mk(S->PA0), mul(mk(S->PA1), ady)), pb = add(mk(S->PB0), mul(mk(S->PB1), ady));
  q->A = add(half(S->det), pa); q->B = add(half(S->det), pb);
  q->al = add(mk(0.5), divi(pa, S->det)); q->be = add(mk(0.5), divi(pb, S->det));
}

/* area of D(y) ^ R, R = S=[0,1]^2 (reg 0) or H=(-inf,1]x[0,1] (reg 1), via Green: area = oint x dz */
static I area_D(const Slice *S, const SliceY *q, int reg) {
  I ax = mul(q->al, S->g1x), az = mul(q->al, S->g1z), bx = mul(q->be, S->g2x), bz = mul(q->be, S->g2z);
  I Vx[4], Vz[4], dx[4], dz[4];
  Vx[0] = sub(sub(q->Xx, ax), bx); Vz[0] = sub(sub(q->Xz, az), bz);
  Vx[1] = sub(add(q->Xx, ax), bx); Vz[1] = sub(add(q->Xz, az), bz);
  Vx[2] = add(add(q->Xx, ax), bx); Vz[2] = add(add(q->Xz, az), bz);
  Vx[3] = add(sub(q->Xx, ax), bx); Vz[3] = add(sub(q->Xz, az), bz);
  dx[0] = scal(2.0, ax); dz[0] = scal(2.0, az);
  dx[1] = scal(2.0, bx); dz[1] = scal(2.0, bz);
  dx[2] = scal(-2.0, ax); dz[2] = scal(-2.0, az);
  dx[3] = scal(-2.0, bx); dz[3] = scal(-2.0, bz);
  I tot = mk(0.0);
  for (int e = 0; e < 4; e++) {
    I t0 = mk(0.0), t1 = mk(1.0);
    /* x <= 1, and x >= 0 for reg 0 */
    if (dx[e].lo > 0) {
      t1 = mini(t1, divi(sub(mk(1.0), Vx[e]), dx[e]));
      if (reg == 0) t0 = maxi(t0, divi(neg(Vx[e]), dx[e]));
    } else if (dx[e].hi < 0) {
      t0 = maxi(t0, divi(sub(mk(1.0), Vx[e]), dx[e]));
      if (reg == 0) t1 = mini(t1, divi(neg(Vx[e]), dx[e]));
    } else { g_err = 1; return mk2(-INFINITY, INFINITY); }
    if (dz[e].lo > 0) {
      t0 = maxi(t0, divi(neg(Vz[e]), dz[e]));
      t1 = mini(t1, divi(sub(mk(1.0), Vz[e]), dz[e]));
    } else if (dz[e].hi < 0) {
      t0 = maxi(t0, divi(sub(mk(1.0), Vz[e]), dz[e]));
      t1 = mini(t1, divi(neg(Vz[e]), dz[e]));
    } else { g_err = 1; return mk2(-INFINITY, INFINITY); }
    t1 = maxi(t1, t0);
    I len = sub(t1, t0);
    if (len.lo < 0) len.lo = 0;
    I mid = half(add(t0, t1));
    I term = mul(mul(dz[e], len), add(Vx[e], mul(dx[e], mid)));
    tot = add(tot, term);
  }
  /* right edge x = 1, z in [0,1], points inside D */
  {
    I ox = sub(mk(1.0), q->Xx);
    I w0 = neg(q->Xz), w1 = sub(mk(1.0), q->Xz);   /* w = z - Xz in [w0, w1] */
    /* constraint A: w*g2x in [ox*g2z - A, ox*g2z + A] */
    I cA = mul(ox, S->g2z);
    I la = sub(cA, q->A), ha = add(cA, q->A);
    if (S->g2x.lo > 0) { w0 = maxi(w0, divi(la, S->g2x)); w1 = mini(w1, divi(ha, S->g2x)); }
    else if (S->g2x.hi < 0) { w0 = maxi(w0, divi(ha, S->g2x)); w1 = mini(w1, divi(la, S->g2x)); }
    else { g_err = 1; return mk2(-INFINITY, INFINITY); }
    /* constraint B: w*g1x in [g1z*ox - B, g1z*ox + B] */
    I cB = mul(S->g1z, ox);
    I lb = sub(cB, q->B), hb = add(cB, q->B);
    if (S->g1x.lo > 0) { w0 = maxi(w0, divi(lb, S->g1x)); w1 = mini(w1, divi(hb, S->g1x)); }
    else if (S->g1x.hi < 0) { w0 = maxi(w0, divi(hb, S->g1x)); w1 = mini(w1, divi(lb, S->g1x)); }
    else { g_err = 1; return mk2(-INFINITY, INFINITY); }
    w1 = maxi(w1, w0);
    I len = sub(w1, w0); if (len.lo < 0) len.lo = 0;
    tot = add(tot, len);
  }
  if (tot.lo < 0) tot.lo = 0;
  return tot;
}

/* ---------------- exact piecewise-quadratic integration over a sub-range ---------------- */
typedef struct { double lo, hi; } Ev;
static int cmpev(const void *a, const void *b) { double x = ((const Ev*)a)->lo, y = ((const Ev*)b)->lo; return x < y ? -1 : x > y; }

/* add event y* = ya + num/den (den must exclude 0) if it can lie in (y0,y1); returns 0 if degenerate */
static int add_event(Ev *ev, int *ne, I num, I den, double y0, double y1) {
  if (!(den.lo > 0 || den.hi < 0)) {
    /* den may vanish: check whether phi(y) = num - (y-ya)*den can vanish on [y0,y1] */
    I dyr = mk2(DN(y0 - YA), UP(y1 - YA));
    I ph = sub(num, mul(dyr, den));
    if (ph.lo > 0 || ph.hi < 0) return 1;
    return 0;
  }
  I yy = add(mk(YA), divi(num, den));
  if (yy.hi <= y0 || yy.lo >= y1) return 1;
  ev[*ne].lo = fmax(yy.lo, y0); ev[*ne].hi = fmin(yy.hi, y1); (*ne)++;
  return 1;
}

static long n_area_evals = 0;
static I f_at(const Slice *S, double y, int reg) { SliceY q; slice_y(S, y, &q); n_area_evals++; return area_D(S, &q, reg); }

/* upper bound of int_{y0}^{y1} area(D(y) ^ R_reg) dy; sigma = sign of (y - ya) on the sub-range */
static double integ_sub(const Slice *S, double y0, double y1, int reg, int sigma) {
  if (!(y1 > y0)) return 0.0;
  /* vertex/corner motion on the sub-range: everything affine in (y - ya) */
  I dt = S->det;
  I pa1 = mk(sigma*S->PA1), pb1 = mk(sigma*S->PB1);
  I al0 = add(mk(0.5), divi(mk(S->PA0), dt)), be0 = add(mk(0.5), divi(mk(S->PB0), dt));
  I al1 = divi(pa1, dt), be1 = divi(pb1, dt);
  I A0 = add(half(dt), mk(S->PA0)), B0 = add(half(dt), mk(S->PB0));
  I kx = mk(S->kc[0]), kz = mk(S->kc[1]);
  Ev ev[96]; int ne = 0; int ok = 1;
  double sgn[2] = {-1.0, 1.0};
  /* (i) vertices of D crossing region lines */
  for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) {
    double s1 = sgn[i], s2 = sgn[j];
    I v0x = add(add(mk(S->xc[0]), mul(mk(s1), mul(al0, S->g1x))), mul(mk(s2), mul(be0, S->g2x)));
    I v0z = add(add(mk(S->xc[1]), mul(mk(s1), mul(al0, S->g1z))), mul(mk(s2), mul(be0, S->g2z)));
    I vpx = add(add(kx, mul(mk(s1), mul(al1, S->g1x))), mul(mk(s2), mul(be1, S->g2x)));
    I vpz = add(add(kz, mul(mk(s1), mul(al1, S->g1z))), mul(mk(s2), mul(be1, S->g2z)));
    double xl[2] = {1.0, 0.0}; int nx = reg == 0 ? 2 : 1;
    for (int a = 0; a < nx; a++) {
      int ne0 = ne;
      ok &= add_event(ev, &ne, sub(mk(xl[a]), v0x), vpx, y0, y1);
      if (ne > ne0) { /* filter: crossing point z must be possibly in [0,1] */
        I yy = mk2(ev[ne-1].lo, ev[ne-1].hi); I zz = add(v0z, mul(sub(yy, mk(YA)), vpz));
        if (zz.lo > 1.0 || zz.hi < 0.0) ne--;
      }
    }
    double zl[2] = {0.0, 1.0};
    for (int a = 0; a < 2; a++) {
      int ne0 = ne;
      ok &= add_event(ev, &ne, sub(mk(zl[a]), v0z), vpz, y0, y1);
      if (ne > ne0) {
        I yy = mk2(ev[ne-1].lo, ev[ne-1].hi); I xx = add(v0x, mul(sub(yy, mk(YA)), vpx));
        if (xx.lo > 1.0 || (reg == 0 && xx.hi < 0.0)) ne--;
      }
    }
  }
  /* (ii) region corners crossing D's edge lines */
  double cx[4] = {0, 1, 1, 0}, cz[4] = {0, 0, 1, 1};
  I cKA = cross2(kx, kz, S->g2x, S->g2z);       /* cross(kc, g2) */
  I cKB = cross2(S->g1x, S->g1z, kx, kz);       /* cross(g1, kc) */
  for (int j = 0; j < 4; j++) {
    if (reg == 1 && (j == 0 || j == 3)) continue;
    I px = sub(mk(cx[j]), mk(S->xc[0])), pz = sub(mk(cz[j]), mk(S->xc[1]));
    I cpA = cross2(px, pz, S->g2x, S->g2z);    /* cross(p - Xc, g2) */
    I cpB = cross2(S->g1x, S->g1z, px, pz);    /* cross(g1, p - Xc) */
    for (int sg = 0; sg < 2; sg++) {
      double sv = sgn[sg];
      /* A-lines: (y-ya)[cross(kc,g2) + sv*sigma*PA1] = cross(p-Xc,g2) - sv*A0 */
      int ne0 = ne;
      ok &= add_event(ev, &ne, sub(cpA, mul(mk(sv), A0)), add(cKA, mul(mk(sv), pa1)), y0, y1);
      if (ne > ne0) {  /* filter: p within the B-range at that time */
        I yy = mk2(ev[ne-1].lo, ev[ne-1].hi); I d = sub(yy, mk(YA));
        I bb = absi(sub(cpB, mul(d, cKB))); I Bv = add(B0, mul(d, pb1));
        if (bb.lo > Bv.hi) ne--;
      }
      ne0 = ne;
      /* B-lines: (y-ya)[cross(g1,kc) + sv*sigma*PB1] = cross(g1,p-Xc) - sv*B0 */
      ok &= add_event(ev, &ne, sub(cpB, mul(mk(sv), B0)), add(cKB, mul(mk(sv), pb1)), y0, y1);
      if (ne > ne0) {
        I yy = mk2(ev[ne-1].lo, ev[ne-1].hi); I d = sub(yy, mk(YA));
        I aa = absi(sub(cpA, mul(d, cKA))); I Av = add(A0, mul(d, pa1));
        if (aa.lo > Av.hi) ne--;
      }
    }
  }
  /* f_max on the sub-range: area(D) = 4 al be det, max at max |y-ya| */
  double dmax = fmax(fabs(y0 - YA), fabs(y1 - YA));
  I alm = add(al0, mul(mk(S->PA1/1.0), divi(mk(UP(dmax)), dt)));
  I bem = add(be0, mul(mk(S->PB1/1.0), divi(mk(UP(dmax)), dt)));
  double fmx = mul(mul(scal(4.0, alm), bem), dt).hi;
  if (!ok) { return UP(UP(y1 - y0)*fmx); }
  qsort(ev, ne, sizeof(Ev), cmpev);
  /* merge */
  int m = 0;
  for (int i = 0; i < ne; i++) {
    if (m > 0 && ev[i].lo <= ev[m-1].hi) { if (ev[i].hi > ev[m-1].hi) ev[m-1].hi = ev[i].hi; }
    else ev[m++] = ev[i];
  }
  I total = mk(0.0);
  double a = y0;
  for (int i = 0; i <= m; i++) {
    double b = (i < m) ? ev[i].lo : y1;
    if (b > a) {
      double mid = 0.5*(a + b);
      if (mid > a && mid < b) {
        I fa = f_at(S, a, reg), fm = f_at(S, mid, reg), fb = f_at(S, b, reg);
        I h = sub(mk(b), mk(a));
        I sI = divi(sub(mk(mid), mk(a)), h);
        I w0 = sub(mk(0.5), divi(mk(1.0), scal(6.0, sI)));
        I w1 = divi(mk(1.0), scal(6.0, mul(sI, sub(mk(1.0), sI))));
        I w2 = divi(sub(mk(1.0/3.0), half(sI)), sub(mk(1.0), sI));
        /* 1/3 is not exact: use an enclosure */
        I third = mk2(DN(1.0/3.0), UP(1.0/3.0));
        w2 = divi(sub(third, half(sI)), sub(mk(1.0), sI));
        I sixth = mk2(DN(1.0/6.0), UP(1.0/6.0));
        w0 = sub(mk(0.5), divi(sixth, sI));
        w1 = divi(sixth, mul(sI, sub(mk(1.0), sI)));
        I q = mul(h, add(add(mul(w0, fa), mul(w1, fm)), mul(w2, fb)));
        total = add(total, q);
      } else {
        total = add(total, mk2(0.0, UP(UP(b - a)*fmx)));
      }
    }
    if (i < m) {
      total = add(total, mk2(0.0, UP(UP(ev[i].hi - ev[i].lo)*fmx)));
      a = ev[i].hi;
    }
  }
  return total.hi;
}

static double integral_bound(const Slice *S) {
  double pts[4]; int np = 0;
  pts[np++] = 0.0; pts[np++] = 1.0; pts[np++] = C;
  if (YA > 0.0 && YA < C && YA != 1.0) pts[np++] = YA;
  for (int i = 1; i < np; i++) { double t = pts[i]; int j = i - 1; while (j >= 0 && pts[j] > t) { pts[j+1] = pts[j]; j--; } pts[j+1] = t; }
  I tot = mk(0.0);
  for (int i = 0; i + 1 < np; i++) {
    double y0 = pts[i], y1 = pts[i+1];
    int reg = (y1 <= 1.0) ? 1 : 0;
    int sigma = (y0 >= YA) ? 1 : -1;
    tot = add(tot, mk(integ_sub(S, y0, y1, reg, sigma)));
  }
  return tot.hi;
}

/* certain separation of D(y) from the unit square [0,1]^2 */
static int separated(const Slice *S, double y) {
  SliceY q; slice_y(S, y, &q);
  I ex = add(mul(q.al, absi(S->g1x)), mul(q.be, absi(S->g2x)));
  I ez = add(mul(q.al, absi(S->g1z)), mul(q.be, absi(S->g2z)));
  if (sub(q.Xx, ex).lo > 1.0 || add(q.Xx, ex).hi < 0.0) return 1;
  if (sub(q.Xz, ez).lo > 1.0 || add(q.Xz, ez).hi < 0.0) return 1;
  double cx[4] = {0, 1, 1, 0}, cz[4] = {0, 0, 1, 1};
  double amin = INFINITY, amax = -INFINITY, bmin = INFINITY, bmax = -INFINITY;
  for (int j = 0; j < 4; j++) {
    I px = sub(mk(cx[j]), q.Xx), pz = sub(mk(cz[j]), q.Xz);
    I a = cross2(px, pz, S->g2x, S->g2z), b = cross2(S->g1x, S->g1z, px, pz);
    if (a.lo < amin) amin = a.lo; if (a.hi > amax) amax = a.hi;
    if (b.lo < bmin) bmin = b.lo; if (b.hi > bmax) bmax = b.hi;
  }
  if (amin > q.A.hi || amax < -q.A.hi) return 1;
  if (bmin > q.B.hi || bmax < -q.B.hi) return 1;
  return 0;
}

/* status: 0 = bound by integral, 1 = by fallback, 2 = infeasible (no admissible tube in the box), 3 = degenerate */
static double box_bound(const Box *b, int *status, Slice *Sout) {
  Slice S; g_err = 0;
  setup(b, &S);
  if (Sout) *Sout = S;
  if (S.kmin2 > K2.hi) { *status = 2; return -INFINITY; }
  if (separated(&S, C) || separated(&S, 1.0)) { if (!g_err) { *status = 2; return -INFINITY; } }
  double fb = UP(U + S.smax);
  if (fb <= TARGET) { *status = 1; return fb; }
  g_err = 0;
  double v = integral_bound(&S);
  if (g_err || !(v == v)) { *status = 3; return fb; }
  if (fb <= v) { *status = 1; return fb; }
  *status = 0; return v;
}

/* ---------------- roots ---------------- */
typedef struct { int nkx, nkz, nt, nx; } Grid;
static int make_root(const Grid *G, long idx, Box *r) {
  long id = idx;
  int ix = id % G->nx; id /= G->nx;
  int iz = id % G->nx; id /= G->nx;
  int it = id % G->nt; id /= G->nt;
  int jz = id % G->nkz; id /= G->nkz;
  int jx = (int)id;
  if (jx >= G->nkx) return 0;
  /* k cells: kx in [-KHI, KHI], kz in [0, KHI] (z-mirror symmetry) */
  r->lo[DKX] = -KHI + 2.0*KHI*jx/G->nkx; r->hi[DKX] = (jx + 1 == G->nkx) ? KHI : -KHI + 2.0*KHI*(jx+1)/G->nkx;
  r->lo[DKZ] = KHI*jz/G->nkz; r->hi[DKZ] = (jz + 1 == G->nkz) ? KHI : KHI*(jz+1)/G->nkz;
  if (jx == 0) r->lo[DKX] = -KHI;
  if (jz == 0) r->lo[DKZ] = 0.0;
  r->lo[DT] = (double)it/G->nt; r->hi[DT] = (it + 1 == G->nt) ? 1.0 : (double)(it+1)/G->nt;
  /* X(ya) range: Sigma(c) and Sigma(1) must meet [0,1]^2; half-extent of a slice <= s/sqrt2 */
  double ax = fmax(fabs(r->lo[DKX]), fabs(r->hi[DKX])), az = fmax(fabs(r->lo[DKZ]), fabs(r->hi[DKZ]));
  I km2 = add(sqri(mk(ax)), sqri(mk(az)));
  I ext = mul(sqrti(add(mk(1.0), km2)), sqrti(mk(0.5)));
  double e = ext.hi;
  double lo[2], hi[2];
  for (int a = 0; a < 2; a++) {
    I kr = mk2(r->lo[a == 0 ? DKX : DKZ], r->hi[a == 0 ? DKX : DKZ]);
    I s1 = mul(sub(mk(YA), mk(C)), kr);   /* (ya - c) k */
    I s2 = mul(sub(mk(YA), mk(1.0)), kr); /* (ya - 1) k */
    double l1 = DN(DN(-e) + s1.lo), h1 = UP(UP(1.0 + e) + s1.hi);
    double l2 = DN(DN(-e) + s2.lo), h2 = UP(UP(1.0 + e) + s2.hi);
    lo[a] = fmax(l1, l2); hi[a] = fmin(h1, h2);
  }
  if (!(lo[0] < hi[0] && lo[1] < hi[1])) return -1;   /* empty: no feasible tube */
  r->lo[DX1] = lo[0] + (hi[0] - lo[0])*ix/G->nx; r->hi[DX1] = (ix + 1 == G->nx) ? hi[0] : lo[0] + (hi[0] - lo[0])*(ix+1)/G->nx;
  r->lo[DX2] = lo[1] + (hi[1] - lo[1])*iz/G->nx; r->hi[DX2] = (iz + 1 == G->nx) ? hi[1] : lo[1] + (hi[1] - lo[1])*(iz+1)/G->nx;
  if (ix == 0) r->lo[DX1] = lo[0];
  if (iz == 0) r->lo[DX2] = lo[1];
  return 1;
}

/* ---------------- search / verify ---------------- */
typedef struct { long nodes, leaves, infeas, fallb, integ, degen, stuck; double maxb; int maxdepth; } Stats;
static Stats ST;
static long hist[40];

static int choose_dim(const Box *b, const Slice *S, int status_degen, double bound) {
  (void)bound;
  double kv2 = (TARGET - U)*(TARGET - U) - 1.0;   /* fallback works if |k|^2 <= kv2 */
  if (status_degen) { /* split the widest dimension (relative) */
    int d = 0; double w = -1;
    for (int i = 0; i < ND; i++) { double wi = b->hi[i] - b->lo[i]; if (wi > w) { w = wi; d = i; } }
    return d;
  }
  if (kv2 > 0 && S->kmin2 < kv2) {
    /* near-vertical box: shrink k so the fallback applies */
    double wx = b->hi[DKX] - b->lo[DKX], wz = b->hi[DKZ] - b->lo[DKZ];
    double kxm = fmax(fabs(b->lo[DKX]), fabs(b->hi[DKX])), kzm = fmax(fabs(b->lo[DKZ]), fabs(b->hi[DKZ]));
    if (kxm*kxm + kzm*kzm > kv2 && S->kmin2 < 0.25*kv2) return (wx >= wz) ? DKX : DKZ;
  }
  int d = 0;
  for (int i = 1; i < ND; i++) if (S->cost[i] > S->cost[d]) d = i;
  return d;
}

static Box stk[4096]; static int stkd[4096];
static FILE *fcodes = NULL;   /* search: output codes; verify: input codes */
static int MODE = 0;          /* 0 search, 1 verify, 2 dump (verify + print sampled leaves) */
static long DUMP_EVERY = 0; static FILE *fdump = NULL;
static int DEPTHMAX = 70;
static int g_fail = 0;

static void run_root(const Box *root) {
  int sp = 0; stk[sp] = *root; stkd[sp] = 0; sp++;
  while (sp > 0) {
    sp--; Box b = stk[sp]; int dep = stkd[sp];
    if (dep > ST.maxdepth) ST.maxdepth = dep;
    ST.nodes++;
    int code;
    int st; Slice S; double v = 0;
    if (MODE >= 1) {
      int ch = fgetc(fcodes);
      if (ch == EOF) { fprintf(stderr, "code stream ended early\n"); g_fail = 1; return; }
      code = ch;
      if (code == 0) {
        v = box_bound(&b, &st, &S);
      }
    } else {
      v = box_bound(&b, &st, &S);
      if (st == 2 || v <= TARGET) code = 0;
      else if (dep >= DEPTHMAX) { code = 0; ST.stuck++; }
      else code = 1 + choose_dim(&b, &S, st == 3, v);
      fputc(code, fcodes);
    }
    if (code == 0) {
      ST.leaves++;
      if (MODE == 2 && fdump && DUMP_EVERY > 0 && (ST.leaves % DUMP_EVERY) == 0) {
        fprintf(fdump, "%ld %d %.17g", ST.leaves, st, v);
        for (int d = 0; d < ND; d++) fprintf(fdump, " %.17g %.17g", b.lo[d], b.hi[d]);
        fprintf(fdump, " %.17g %.17g %.17g %.17g", S.kc[0], S.kc[1], S.xc[0], S.xc[1]);
        I gg[4] = {S.g1x, S.g1z, S.g2x, S.g2z};
        for (int i = 0; i < 4; i++) fprintf(fdump, " %.17g %.17g", gg[i].lo, gg[i].hi);
        fprintf(fdump, " %.17g %.17g %.17g %.17g %.17g %.17g\n", S.det.lo, S.det.hi, S.PA0, S.PA1, S.PB0, S.PB1);
      }
      if (st == 2) ST.infeas++;
      else {
        if (st == 1) ST.fallb++; else if (st == 3) ST.degen++; else ST.integ++;
        if (v > ST.maxb) ST.maxb = v;
        if (v > TARGET) { g_fail = 1; if (ST.stuck < 20 || MODE == 1) fprintf(stderr, "leaf above target: v=%.9f st=%d k=[%.6g,%.6g]x[%.6g,%.6g] t=[%.6g,%.6g] X=[%.6g,%.6g]x[%.6g,%.6g]\n", v, st, b.lo[0], b.hi[0], b.lo[1], b.hi[1], b.lo[2], b.hi[2], b.lo[3], b.hi[3], b.lo[4], b.hi[4]); }
        int hb = (int)(sqrt(0.5*(S.kmax2 + S.kmin2))/0.1); if (hb > 39) hb = 39; if (hb < 0) hb = 0; hist[hb]++;
      }
      continue;
    }
    int d = code - 1;
    if (d < 0 || d >= ND) { fprintf(stderr, "bad code %d\n", code); g_fail = 1; return; }
    double m = 0.5*(b.lo[d] + b.hi[d]);
    if (!(m > b.lo[d] && m < b.hi[d])) { fprintf(stderr, "cannot split\n"); g_fail = 1; return; }
    if (sp + 2 >= 4096) { fprintf(stderr, "stack overflow\n"); g_fail = 1; return; }
    Box c1 = b, c2 = b; c1.hi[d] = m; c2.lo[d] = m;
    stk[sp] = c2; stkd[sp] = dep + 1; sp++;
    stk[sp] = c1; stkd[sp] = dep + 1; sp++;
  }
}

int main(int argc, char **argv) {
  if (argc >= 2 && strcmp(argv[1], "roots") == 0) {
    /* print the root k-cells and their X(ya) ranges: roots u ya nkx nkz nt nx */
    U = atof(argv[2]); YA = atof(argv[3]); C = 1.0 + U;
    Grid G = {atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7])};
    I u = mk(U), u2 = sqri(u); K2 = divi(scal(8.0, u2), sqri(sub(u2, mk(2.0)))); KHI = sqrti(mk(K2.hi)).hi;
    printf("KHI %.17g\n", KHI);
    for (int jx = 0; jx < G.nkx; jx++) for (int jz = 0; jz < G.nkz; jz++) {
      long idx = ((long)jx*G.nkz + jz)*G.nt*G.nx*G.nx;
      Box r0, r1; int a = make_root(&G, idx, &r0); int b = make_root(&G, idx + (long)G.nx*G.nx - 1, &r1);
      if (a < 0 || b < 0) { printf("CELL %d %d EMPTY\n", jx, jz); continue; }
      printf("CELL %d %d %.17g %.17g %.17g %.17g X %.17g %.17g %.17g %.17g\n", jx, jz, r0.lo[DKX], r0.hi[DKX], r0.lo[DKZ], r0.hi[DKZ], r0.lo[DX1], r1.hi[DX1], r0.lo[DX2], r1.hi[DX2]);
    }
    for (int it = 0; it < G.nt; it++) { Box r; make_root(&G, (long)it*G.nx*G.nx, &r); printf("TCELL %d %.17g %.17g\n", it, r.lo[DT], r.hi[DT]); }
    return 0;
  }
  if (argc < 10) { fprintf(stderr, "usage: sofaub search|verify|dump u target ya nkx nkz nt nx codesfile [root_from root_to [byte_offset]] | roots u ya nkx nkz nt nx\n"); return 2; }
  MODE = (strcmp(argv[1], "verify") == 0) ? 1 : (strcmp(argv[1], "dump") == 0 ? 2 : 0);
  if (MODE == 2) { DUMP_EVERY = atol(getenv("DUMP_EVERY") ? getenv("DUMP_EVERY") : "1000"); fdump = fopen(getenv("DUMP_FILE") ? getenv("DUMP_FILE") : "dump.txt", "w"); }
  U = atof(argv[2]); TARGET = atof(argv[3]); YA = atof(argv[4]);
  Grid G = {atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8])};
  C = 1.0 + U;
  if (C - 1.0 != U) { fprintf(stderr, "1+u not exact\n"); return 2; }
  I u = mk(U), u2 = sqri(u);
  K2 = divi(scal(8.0, u2), sqri(sub(u2, mk(2.0))));
  KHI = sqrti(mk(K2.hi)).hi;
  long r0 = 0, r1 = (long)G.nkx*G.nkz*G.nt*G.nx*G.nx;
  if (argc >= 12) { r0 = atol(argv[10]); r1 = atol(argv[11]); }
  fcodes = fopen(argv[9], MODE ? "rb" : "wb");
  if (MODE >= 1 && argc >= 13) { long skip = atol(argv[12]); if (fseek(fcodes, skip, SEEK_SET) != 0) { fprintf(stderr, "seek failed\n"); return 2; } }
  if (!fcodes) { fprintf(stderr, "cannot open %s\n", argv[9]); return 2; }
  setvbuf(stdout, NULL, _IOLBF, 0);
  printf("%s u=%.17g c=%.17g target=%.17g ya=%.17g K^2 in [%.17g, %.17g] Khi=%.17g grid=%d,%d,%d,%d roots [%ld,%ld)\n",
         MODE == 2 ? "DUMP" : (MODE ? "VERIFY" : "SEARCH"), U, C, TARGET, YA, K2.lo, K2.hi, KHI, G.nkx, G.nkz, G.nt, G.nx, r0, r1);
  clock_t t0 = clock();
  long nroot = 0, nempty = 0;
  for (long i = r0; i < r1 && !g_fail; i++) {
    Box r; int ok = make_root(&G, i, &r);
    if (ok == 0) break;
    if (ok < 0) { nempty++; continue; }
    run_root(&r); nroot++;
    if ((i & 1023) == 1023) printf("  root %ld nodes=%ld leaves=%ld maxb=%.6f t=%.0fs\n", i+1, ST.nodes, ST.leaves, ST.maxb, (double)(clock()-t0)/CLOCKS_PER_SEC);
  }
  if (MODE >= 1 && r0 == 0 && r1 == (long)G.nkx*G.nkz*G.nt*G.nx*G.nx) { int ch = fgetc(fcodes); if (ch != EOF) { fprintf(stderr, "code stream has extra data\n"); g_fail = 1; } }
  if (fdump) fclose(fdump);
  fclose(fcodes);
  double secs = (double)(clock() - t0)/CLOCKS_PER_SEC;
  printf("roots=%ld empty_roots=%ld nodes=%ld leaves=%ld infeasible=%ld fallback=%ld integral=%ld degenerate=%ld stuck=%ld maxdepth=%d max_leaf_bound=%.12f area_evals=%ld time=%.1fs\n",
         nroot, nempty, ST.nodes, ST.leaves, ST.infeas, ST.fallb, ST.integ, ST.degen, ST.stuck, ST.maxdepth, ST.maxb, n_area_evals, secs);
  printf("RESULT %s\n", g_fail ? "FAIL" : "PASS");
  printf("leaves by |k| (0.1 bins):"); for (int i = 0; i < 40; i++) if (hist[i]) printf(" %d:%ld", i, hist[i]); printf("\n");
  return g_fail ? 1 : 0;
}
