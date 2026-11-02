/* lowcert.c -- rigorous (FLINT/arb ball arithmetic) lower bound for
 *   inf_{x in [X0,X1]} A(x),  A(x) = U_mu(x) - sum_j c_j log|Q_j(x)|,
 * mu = piecewise-constant density: mass m_k uniformly on [e_k, e_{k+1}].
 * U_mu(x) = sum_k beta_k F(x - e_k),  F(u) = u log|u| - u,  beta_k = rho_k - rho_{k-1}.
 * Bisection; on I = [t-H, t+H]:
 *   near terms (edge or real root within 2H of t): exact monotone bounds
 *     (F decreasing on [-1,1]; log|x-r| <= log(|t-r|+H));
 *   far terms: second-order Taylor bound  f(t) - |f'(t)| H + min(0, inf_I f'') H^2/2,
 *     inf_I f'' by ball evaluation on I.
 * Input (stdin, text): K / e_0..e_K / m_0..m_{K-1} / NP / {c deg a_0..a_deg} x NP / X0 X1 V MAXDEPTH
 * Doubles are read with strtod (decimal %.17g -> exact same double).
 * Output: OK/FAIL, #intervals, min leaf lower bound.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/arb_fmpz_poly.h>

#define PREC 80

static int K; static double *ed; static arb_ptr ea, beta; static int *bnz; static int nbnz;
static int NR; static double *rre, *rim; static arb_ptr ra, rb, rc; static int *risreal;
static arb_t constpart; static arb_t one;

static void Ffun(arb_t res, const arb_t u, slong prec) {
    /* F(u) = u log|u| - u ; F(0)=0 ; if the ball u contains 0: |F(u)| <= r|log r| + r, r=|u|max<=1/e */
    if (arb_is_zero(u)) { arb_zero(res); return; }
    if (arb_contains_zero(u)) {
        arb_t r, t; arb_init(r); arb_init(t);
        arb_get_abs_ubound_arf(arb_midref(r), u, prec); mag_zero(arb_radref(r));
        if (arf_cmp_d(arb_midref(r), 0.36) < 0) {
            arb_log(t, r, prec); arb_mul(t, t, r, prec); arb_abs(t, t);
            arb_zero(res); arb_add_error(res, t); arb_add_error(res, r);
        } else { arb_zero(res); arb_add_error_2exp_si(res, 2); }
        arb_clear(r); arb_clear(t); return;
    }
    arb_t a; arb_init(a);
    arb_abs(a, u); arb_log(a, a, prec); arb_mul(a, a, u, prec); arb_sub(res, a, u, prec);
    arb_clear(a);
}

static long nint = 0; static double minleaf = 1e300; static int failed = 0; static int maxd = 60;
static double Vt;

/* lower bound of A on [lo,hi] (exact doubles). returns double lower bound (rounded down) */
static double lbound(double lo, double hi) {
    double t = 0.5 * (lo + hi);   /* exact for doubles in [0,1] unless underflow */
    double H = 0.5 * (hi - lo);   /* exact: hi-lo exact? not always -> use arb below */
    arb_t T, Hb, I, val, der, curv, tmp, u, w, lb, sq, bb;
    arb_init(T); arb_init(Hb); arb_init(I); arb_init(val); arb_init(der); arb_init(curv);
    arb_init(tmp); arb_init(u); arb_init(w); arb_init(lb); arb_init(sq); arb_init(bb);
    arb_set_d(T, t);
    /* H = max(hi - t, t - lo) computed exactly in arb */
    arb_set_d(tmp, hi); arb_sub(Hb, tmp, T, PREC); arb_set_d(tmp, lo); arb_sub(tmp, T, tmp, PREC);
    arb_max(Hb, Hb, tmp, PREC);
    /* I as ball: [lo,hi] */
    arb_set_d(I, lo); arb_set_d(tmp, hi); arb_union(I, I, tmp, PREC);
    arb_set(val, constpart); arb_zero(der); arb_zero(curv);
    double Hd = H * (1 + 1e-15) + 1e-300;
    for (int q = 0; q < nbnz; q++) {
        int k = bnz[q];
        double dist = fabs(t - ed[k]);
        if (dist <= 2.0 * Hd) {
            /* min over I of beta F(x - e): F decreasing on [-1,1] */
            if (arb_is_positive(beta + k)) { arb_set_d(tmp, hi); arb_sub(u, tmp, ea + k, PREC); Ffun(w, u, PREC); arb_mul(w, w, beta + k, PREC); }
            else if (arb_is_negative(beta + k)) { arb_set_d(tmp, lo); arb_sub(u, tmp, ea + k, PREC); Ffun(w, u, PREC); arb_mul(w, w, beta + k, PREC); }
            else { arb_t w2; arb_init(w2);
                arb_set_d(tmp, hi); arb_sub(u, tmp, ea + k, PREC); Ffun(w, u, PREC); arb_mul(w, w, beta + k, PREC);
                arb_set_d(tmp, lo); arb_sub(u, tmp, ea + k, PREC); Ffun(w2, u, PREC); arb_mul(w2, w2, beta + k, PREC);
                arb_union(w, w, w2, PREC); arb_clear(w2); }
            arb_add(val, val, w, PREC);
        } else {
            arb_sub(u, T, ea + k, PREC);
            arb_abs(w, u); arb_log(w, w, PREC); arb_addmul(der, beta + k, w, PREC);   /* F' = log|u| */
            arb_sub(w, w, one, PREC); arb_mul(w, w, u, PREC); arb_addmul(val, beta + k, w, PREC); /* F = u(log|u|-1) */
            arb_sub(w, I, ea + k, PREC); arb_div(w, beta + k, w, PREC); arb_add(curv, curv, w, PREC);
        }
    }
    for (int r = 0; r < NR; r++) {
        if (risreal[r]) {
            double dist = fabs(t - rre[r]);
            if (dist <= 2.0 * Hd + 1e-12) {
                arb_sub(u, T, ra + r, PREC); arb_abs(u, u); arb_add(u, u, Hb, PREC);
                arb_log(w, u, PREC); arb_submul(val, rc + r, w, PREC);
            } else {
                arb_sub(u, T, ra + r, PREC);
                arb_abs(w, u); arb_log(w, w, PREC); arb_submul(val, rc + r, w, PREC);
                arb_inv(w, u, PREC); arb_submul(der, rc + r, w, PREC);
                arb_sub(w, I, ra + r, PREC); arb_sqr(w, w, PREC); arb_div(w, rc + r, w, PREC);
                arb_add(curv, curv, w, PREC);
            }
        } else {
            /* pair a +- ib : -c log((x-a)^2+b^2) */
            arb_sub(u, T, ra + r, PREC); arb_sqr(sq, u, PREC); arb_sqr(bb, rb + r, PREC);
            arb_add(sq, sq, bb, PREC);           /* (t-a)^2+b^2 */
            arb_log(w, sq, PREC); arb_submul(val, rc + r, w, PREC);
            arb_mul_2exp_si(w, u, 1); arb_div(w, w, sq, PREC); arb_submul(der, rc + r, w, PREC);
            /* curvature of -c log(u^2+b^2): -c*2(b^2-u^2)/(u^2+b^2)^2 on I */
            arb_sub(u, I, ra + r, PREC); arb_sqr(u, u, PREC);   /* u^2 on I */
            arb_sub(w, bb, u, PREC); arb_add(sq, u, bb, PREC); arb_sqr(sq, sq, PREC);
            arb_div(w, w, sq, PREC); arb_mul_2exp_si(w, w, 1); arb_submul(curv, rc + r, w, PREC);
        }
    }
    /* lb = val - |der| H + min(0, curv_low) H^2/2 */
    arb_abs(w, der); arb_mul(w, w, Hb, PREC); arb_sub(lb, val, w, PREC);
    arf_t cl; arf_init(cl); arb_get_lbound_arf(cl, curv, PREC);
    if (arf_sgn(cl) < 0) {
        arb_set_arf(w, cl); arb_mul(w, w, Hb, PREC); arb_mul(w, w, Hb, PREC); arb_mul_2exp_si(w, w, -1);
        arb_add(lb, lb, w, PREC);
    }
    arb_get_lbound_arf(cl, lb, PREC);
    double res = arf_get_d(cl, ARF_RND_FLOOR);
    if (!arb_is_finite(lb)) res = -1e300;
    arf_clear(cl);
    arb_clear(T); arb_clear(Hb); arb_clear(I); arb_clear(val); arb_clear(der); arb_clear(curv);
    arb_clear(tmp); arb_clear(u); arb_clear(w); arb_clear(lb); arb_clear(sq); arb_clear(bb);
    return res;
}

static void rec(double lo, double hi, int depth) {
    nint++;
    double b = lbound(lo, hi);
    if (b >= Vt) { if (b < minleaf) minleaf = b; return; }
    if (depth >= maxd) { failed++; if (failed < 20) fprintf(stderr, "FAIL at [%.17g, %.17g] lb=%.12g\n", lo, hi, b); if (b < minleaf) minleaf = b; return; }
    double t = 0.5 * (lo + hi);
    rec(lo, t, depth + 1); rec(t, hi, depth + 1);
}

int main(int argc, char **argv) {
    arb_init(one); arb_one(one);
    if (scanf("%d", &K) != 1) return 2;
    ed = malloc(sizeof(double) * (K + 1));
    double *md = malloc(sizeof(double) * K);
    char buf[4096];
    for (int k = 0; k <= K; k++) { scanf("%s", buf); ed[k] = strtod(buf, NULL); }
    for (int k = 0; k < K; k++) { scanf("%s", buf); md[k] = strtod(buf, NULL); }
    ea = _arb_vec_init(K + 1); beta = _arb_vec_init(K + 1);
    for (int k = 0; k <= K; k++) arb_set_d(ea + k, ed[k]);
    arb_t rho, rhoprev, h; arb_init(rho); arb_init(rhoprev); arb_init(h);
    arb_zero(rhoprev);
    bnz = malloc(sizeof(int) * (K + 1)); nbnz = 0;
    for (int k = 0; k <= K; k++) {
        if (k < K) { arb_sub(h, ea + k + 1, ea + k, PREC); arb_set_d(rho, md[k]); arb_div(rho, rho, h, PREC); }
        else arb_zero(rho);
        arb_sub(beta + k, rho, rhoprev, PREC);
        arb_set(rhoprev, rho);
        if (!arb_is_zero(beta + k)) bnz[nbnz++] = k;
    }
    int NP; scanf("%d", &NP);
    arb_init(constpart); arb_zero(constpart);
    int cap = 4096; NR = 0;
    rre = malloc(sizeof(double) * cap); rim = malloc(sizeof(double) * cap); risreal = malloc(sizeof(int) * cap);
    ra = _arb_vec_init(cap); rb = _arb_vec_init(cap); rc = _arb_vec_init(cap);
    for (int p = 0; p < NP; p++) {
        double c; int deg; scanf("%s", buf); c = strtod(buf, NULL); scanf("%d", &deg);
        fmpz_poly_t P; fmpz_poly_init(P);
        for (int i = 0; i <= deg; i++) { fmpz_t z; fmpz_init(z); scanf("%s", buf); fmpz_set_str(z, buf, 10); fmpz_poly_set_coeff_fmpz(P, i, z); fmpz_clear(z); }
        acb_ptr roots = _acb_vec_init(deg);
        arb_fmpz_poly_complex_roots(roots, P, 0, 2 * PREC);
        arb_t cb, la; arb_init(cb); arb_init(la);
        arb_set_d(cb, c);
        fmpz_t lead; fmpz_init(lead); fmpz_poly_get_coeff_fmpz(lead, P, deg); fmpz_abs(lead, lead);
        arb_set_fmpz(la, lead); arb_log(la, la, PREC); arb_submul(constpart, cb, la, PREC);
        for (int i = 0; i < deg; i++) {
            acb_ptr z = roots + i;
            if (arb_is_zero(acb_imagref(z))) {
                arb_set(ra + NR, acb_realref(z)); arb_zero(rb + NR); arb_set(rc + NR, cb); risreal[NR] = 1;
                rre[NR] = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR); rim[NR] = 0; NR++;
            } else if (arf_sgn(arb_midref(acb_imagref(z))) > 0) {
                arb_set(ra + NR, acb_realref(z)); arb_abs(rb + NR, acb_imagref(z)); arb_set(rc + NR, cb); risreal[NR] = 0;
                rre[NR] = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
                rim[NR] = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR); NR++;
            }
        }
        _acb_vec_clear(roots, deg); arb_clear(cb); arb_clear(la); fmpz_clear(lead); fmpz_poly_clear(P);
    }
    double X0, X1; scanf("%s", buf); X0 = strtod(buf, NULL); scanf("%s", buf); X1 = strtod(buf, NULL);
    scanf("%s", buf); Vt = strtod(buf, NULL); scanf("%d", &maxd);
    fprintf(stderr, "K=%d nonzero jumps=%d roots(real+pairs)=%d  [%g,%g] V=%.15g maxdepth=%d\n", K, nbnz, NR, X0, X1, Vt, maxd);
    /* start from a uniform split into 4096 pieces */
    int N0 = 4096;
    for (int i = 0; i < N0; i++) {
        double lo = X0 + (X1 - X0) * i / N0, hi = X0 + (X1 - X0) * (i + 1) / N0;
        if (i == N0 - 1) hi = X1;
        rec(lo, hi, 0);
    }
    printf("%s intervals=%ld minleaf=%.15g failed=%d\n", failed ? "FAIL" : "OK", nint, minleaf, failed);
    return failed ? 1 : 0;
}
