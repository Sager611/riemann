/* upcert.c -- rigorous (arb) upper bound for sup_{w in [W0,W1]} f(w),
 *   f(w) = sum_j g_j log|R_j(w)|, g_j = e_j / N (exact rationals).
 * On I=[t-H,t+H]: near real roots (within 2H): log|w-r| <= log(|t-r|+H);
 * far terms: f(t) + |f'(t)| H + max(0, sup_I f'') H^2/2 (sup by ball evaluation).
 * Input: NF / {e deg a_0..a_deg} x NF / N / W0 W1 T MAXDEPTH
 * Output: OK/FAIL intervals maxleaf
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/arb_fmpz_poly.h>
#define PREC 96

static int NR; static double *rre; static arb_ptr ra, rb, rg; static int *risreal;
static arb_t constpart; static long nint = 0; static double maxleaf = -1e300; static int failed = 0;
static int maxd; static double Tt;

static double ubound(double lo, double hi) {
    double t = 0.5 * (lo + hi), H = 0.5 * (hi - lo);
    arb_t T, Hb, I, val, der, curv, tmp, u, w, ub, sq, bb;
    arb_init(T); arb_init(Hb); arb_init(I); arb_init(val); arb_init(der); arb_init(curv);
    arb_init(tmp); arb_init(u); arb_init(w); arb_init(ub); arb_init(sq); arb_init(bb);
    arb_set_d(T, t);
    arb_set_d(tmp, hi); arb_sub(Hb, tmp, T, PREC); arb_set_d(tmp, lo); arb_sub(tmp, T, tmp, PREC); arb_max(Hb, Hb, tmp, PREC);
    arb_set_d(I, lo); arb_set_d(tmp, hi); arb_union(I, I, tmp, PREC);
    arb_set(val, constpart); arb_zero(der); arb_zero(curv);
    double Hd = H * (1 + 1e-15) + 1e-300;
    for (int r = 0; r < NR; r++) {
        if (risreal[r]) {
            if (fabs(t - rre[r]) <= 2.0 * Hd + 1e-12) {
                arb_sub(u, T, ra + r, PREC); arb_abs(u, u); arb_add(u, u, Hb, PREC);
                arb_log(w, u, PREC); arb_addmul(val, rg + r, w, PREC);
            } else {
                arb_sub(u, T, ra + r, PREC);
                arb_abs(w, u); arb_log(w, w, PREC); arb_addmul(val, rg + r, w, PREC);
                arb_inv(w, u, PREC); arb_addmul(der, rg + r, w, PREC);
                arb_sub(w, I, ra + r, PREC); arb_sqr(w, w, PREC); arb_div(w, rg + r, w, PREC);
                arb_sub(curv, curv, w, PREC);           /* g log|w-r| has f'' = -g/(w-r)^2 */
            }
        } else {
            arb_sub(u, T, ra + r, PREC); arb_sqr(sq, u, PREC); arb_sqr(bb, rb + r, PREC); arb_add(sq, sq, bb, PREC);
            arb_log(w, sq, PREC); arb_addmul(val, rg + r, w, PREC);
            arb_mul_2exp_si(w, u, 1); arb_div(w, w, sq, PREC); arb_addmul(der, rg + r, w, PREC);
            arb_sub(u, I, ra + r, PREC); arb_sqr(u, u, PREC);
            arb_sub(w, bb, u, PREC); arb_add(sq, u, bb, PREC); arb_sqr(sq, sq, PREC);
            arb_div(w, w, sq, PREC); arb_mul_2exp_si(w, w, 1); arb_addmul(curv, rg + r, w, PREC);
        }
    }
    arb_abs(w, der); arb_mul(w, w, Hb, PREC); arb_add(ub, val, w, PREC);
    arf_t cu; arf_init(cu); arb_get_ubound_arf(cu, curv, PREC);
    if (arf_sgn(cu) > 0) { arb_set_arf(w, cu); arb_mul(w, w, Hb, PREC); arb_mul(w, w, Hb, PREC); arb_mul_2exp_si(w, w, -1); arb_add(ub, ub, w, PREC); }
    arb_get_ubound_arf(cu, ub, PREC);
    double res = arf_get_d(cu, ARF_RND_CEIL);
    if (!arb_is_finite(ub)) res = 1e300;
    arf_clear(cu);
    arb_clear(T); arb_clear(Hb); arb_clear(I); arb_clear(val); arb_clear(der); arb_clear(curv);
    arb_clear(tmp); arb_clear(u); arb_clear(w); arb_clear(ub); arb_clear(sq); arb_clear(bb);
    return res;
}

static void rec(double lo, double hi, int depth) {
    nint++;
    double b = ubound(lo, hi);
    if (b <= Tt) { if (b > maxleaf) maxleaf = b; return; }
    if (depth >= maxd) { failed++; if (failed < 20) fprintf(stderr, "FAIL [%.17g,%.17g] ub=%.15g\n", lo, hi, b); if (b > maxleaf) maxleaf = b; return; }
    double t = 0.5 * (lo + hi);
    rec(lo, t, depth + 1); rec(t, hi, depth + 1);
}

int main(void) {
    char buf[8192]; int NF;
    if (scanf("%d", &NF) != 1) return 2;
    int cap = 8192; NR = 0;
    rre = malloc(sizeof(double) * cap); risreal = malloc(sizeof(int) * cap);
    ra = _arb_vec_init(cap); rb = _arb_vec_init(cap); rg = _arb_vec_init(cap);
    fmpz_t *ev = malloc(sizeof(fmpz_t) * NF); int *dg = malloc(sizeof(int) * NF);
    fmpz_poly_struct *P = malloc(sizeof(fmpz_poly_struct) * NF);
    for (int j = 0; j < NF; j++) {
        fmpz_init(ev[j]); scanf("%s", buf); fmpz_set_str(ev[j], buf, 10); scanf("%d", dg + j);
        fmpz_poly_init(P + j);
        for (int i = 0; i <= dg[j]; i++) { fmpz_t z; fmpz_init(z); scanf("%s", buf); fmpz_set_str(z, buf, 10); fmpz_poly_set_coeff_fmpz(P + j, i, z); fmpz_clear(z); }
    }
    fmpz_t N, Nsum; fmpz_init(N); fmpz_init(Nsum); scanf("%s", buf); fmpz_set_str(N, buf, 10);
    for (int j = 0; j < NF; j++) fmpz_addmul_ui(Nsum, ev[j], dg[j]);
    if (!fmpz_equal(N, Nsum)) { printf("DEGREE MISMATCH\n"); return 3; }
    arb_init(constpart); arb_zero(constpart);
    arb_t g, la; arb_init(g); arb_init(la);
    for (int j = 0; j < NF; j++) {
        arb_set_fmpz(g, ev[j]); arb_div_fmpz(g, g, N, PREC);
        fmpz_t lead; fmpz_init(lead); fmpz_poly_get_coeff_fmpz(lead, P + j, dg[j]); fmpz_abs(lead, lead);
        arb_set_fmpz(la, lead); arb_log(la, la, PREC); arb_addmul(constpart, g, la, PREC); fmpz_clear(lead);
        acb_ptr roots = _acb_vec_init(dg[j]);
        arb_fmpz_poly_complex_roots(roots, P + j, 0, 2 * PREC);
        for (int i = 0; i < dg[j]; i++) {
            acb_ptr z = roots + i;
            if (arb_is_zero(acb_imagref(z))) {
                arb_set(ra + NR, acb_realref(z)); arb_zero(rb + NR); arb_set(rg + NR, g); risreal[NR] = 1;
                rre[NR] = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR); NR++;
            } else if (arf_sgn(arb_midref(acb_imagref(z))) > 0) {
                arb_set(ra + NR, acb_realref(z)); arb_abs(rb + NR, acb_imagref(z)); arb_set(rg + NR, g); risreal[NR] = 0;
                rre[NR] = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR); NR++;
            }
        }
        _acb_vec_clear(roots, dg[j]);
    }
    double W0, W1; scanf("%s", buf); W0 = strtod(buf, NULL); scanf("%s", buf); W1 = strtod(buf, NULL);
    scanf("%s", buf); Tt = strtod(buf, NULL); scanf("%d", &maxd);
    fprintf(stderr, "NF=%d roots(real+pairs)=%d [%g,%g] T=%.17g\n", NF, NR, W0, W1, Tt);
    int N0 = 4096;
    for (int i = 0; i < N0; i++) {
        double lo = W0 + (W1 - W0) * i / N0, hi = W0 + (W1 - W0) * (i + 1) / N0;
        if (i == N0 - 1) hi = W1;
        rec(lo, hi, 0);
    }
    printf("%s intervals=%ld maxleaf=%.17g failed=%d\n", failed ? "FAIL" : "OK", nint, maxleaf, failed);
    return failed ? 1 : 0;
}
