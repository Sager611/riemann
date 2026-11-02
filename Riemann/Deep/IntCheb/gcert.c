/* gcert.c -- rigorous (arb) lower bounds for the divisibility terms
 *   G_i = (1/d_i) [ int log|Q_i| dmu - sum_{j != i} c_j log|Res(Q_i,Q_j)| ]
 * (equal to the paper's G_i when mu(total)=1), with
 *   int log|Q_i| dmu = M log a_i + sum_{alpha: Q_i(alpha)=0} U_mu(alpha),
 *   U_mu(z) = sum_k beta_k Fc(z - e_k),  Fc(u) = Re(u log u - u).
 * Input (stdin): K / edges / masses / J / {c deg a_0..a_deg} x J
 * Output: one line per i: index lowerbound(G_i); then "MIN <value> <index>".
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

static void Freal(arb_t res, const arb_t u, slong prec) {
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

int main(void) {
    int K; char buf[8192];
    if (scanf("%d", &K) != 1) return 2;
    arb_ptr e = _arb_vec_init(K + 1), beta = _arb_vec_init(K + 1);
    double *md = malloc(sizeof(double) * K);
    for (int k = 0; k <= K; k++) { scanf("%s", buf); arb_set_d(e + k, strtod(buf, NULL)); }
    for (int k = 0; k < K; k++) { scanf("%s", buf); md[k] = strtod(buf, NULL); }
    arb_t rho, rp, h, M; arb_init(rho); arb_init(rp); arb_init(h); arb_init(M);
    arb_zero(rp); arb_zero(M);
    for (int k = 0; k <= K; k++) {
        if (k < K) { arb_sub(h, e + k + 1, e + k, PREC); arb_set_d(rho, md[k]); arb_div(rho, rho, h, PREC);
                     arb_set_d(h, md[k]); arb_add(M, M, h, PREC); }
        else arb_zero(rho);
        arb_sub(beta + k, rho, rp, PREC); arb_set(rp, rho);
    }
    int J; scanf("%d", &J);
    fmpz_poly_struct *P = malloc(sizeof(fmpz_poly_struct) * J);
    arb_ptr c = _arb_vec_init(J); int *deg = malloc(sizeof(int) * J);
    for (int j = 0; j < J; j++) {
        scanf("%s", buf); arb_set_d(c + j, strtod(buf, NULL)); scanf("%d", deg + j);
        fmpz_poly_init(P + j);
        for (int i = 0; i <= deg[j]; i++) { fmpz_t z; fmpz_init(z); scanf("%s", buf); fmpz_set_str(z, buf, 10); fmpz_poly_set_coeff_fmpz(P + j, i, z); fmpz_clear(z); }
    }
    arb_t S, val, t1, t2, re, im, lg, ar; acb_t u;
    arb_init(S); arb_init(val); arb_init(t1); arb_init(t2); arb_init(re); arb_init(im); arb_init(lg); arb_init(ar); acb_init(u);
    fmpz_t R; fmpz_init(R);
    double gmin = 1e300; int imin = -1;
    for (int i = 0; i < J; i++) {
        acb_ptr roots = _acb_vec_init(deg[i]);
        arb_fmpz_poly_complex_roots(roots, P + i, 0, 2 * PREC);
        arb_zero(S);
        for (int r = 0; r < deg[i]; r++) {
            for (int k = 0; k <= K; k++) {
                if (arb_is_zero(beta + k)) continue;
                acb_sub_arb(u, roots + r, e + k, PREC);
                if (arb_is_zero(acb_imagref(u))) {
                    Freal(val, acb_realref(u), PREC);
                } else {
                    /* Re(u log u - u) = x log|u| - y arg(u) - x */
                    acb_abs(lg, u, PREC); arb_log(lg, lg, PREC);
                    arb_atan2(ar, acb_imagref(u), acb_realref(u), PREC);
                    arb_mul(val, acb_realref(u), lg, PREC);
                    arb_submul(val, acb_imagref(u), ar, PREC);
                    arb_sub(val, val, acb_realref(u), PREC);
                }
                arb_addmul(S, beta + k, val, PREC);
            }
        }
        /* + M log a_i */
        fmpz_t lead; fmpz_init(lead); fmpz_poly_get_coeff_fmpz(lead, P + i, deg[i]); fmpz_abs(lead, lead);
        arb_set_fmpz(t1, lead); arb_log(t1, t1, PREC); arb_addmul(S, M, t1, PREC); fmpz_clear(lead);
        /* - sum_{j != i} c_j log|Res| */
        for (int j = 0; j < J; j++) {
            if (j == i || arb_is_zero(c + j)) continue;
            fmpz_poly_resultant(R, P + i, P + j);
            if (fmpz_is_zero(R)) { printf("ZERO RESULTANT %d %d\n", i, j); return 3; }
            fmpz_abs(R, R); arb_set_fmpz(t1, R); arb_log(t1, t1, PREC);
            arb_submul(S, c + j, t1, PREC);
        }
        arb_div_si(S, S, deg[i], PREC);
        arf_t lb; arf_init(lb); arb_get_lbound_arf(lb, S, PREC);
        double g = arf_get_d(lb, ARF_RND_FLOOR); arf_clear(lb);
        printf("%d %.17g %.3g\n", i, g, mag_get_d(arb_radref(S)));
        if (g < gmin) { gmin = g; imin = i; }
        _acb_vec_clear(roots, deg[i]);
    }
    printf("MIN %.17g %d\n", gmin, imin);
    return 0;
}
