/* deig.c (bothaaudit): lowest k eigenvalues of a truncated Weil block given by gen.py's generators (orientation only, no rigour).
   T_ij = ps p_i p_j + a_i a_j (u_j - u_i)/(x_i - x_j) (i != j), T_ii = d_i. Dense LDL^T in place (packed lower rows, FLINT arb_approx_dot,
   mid-points only), then block inverse iteration (k vectors, its steps, Gram-Schmidt) and the Rayleigh-Ritz matrix H = Y^T T Y = (L^T Y)^T D (L^T Y).
   Output (stdout): n, #negative pivots, min pivot, per-iteration log10 diag(H), and H (k x k) as decimal strings for the Python eigensolve.
   Usage: deig <genfile> <prec> <k> <its> [vecfile]   (vecfile: the k final orthonormalised vectors, hex mantissa + exponent per entry)      cc -O2 -o deig deig.c -lflint -lgmp -lmpfr -lm -I/opt/homebrew/include -L/opt/homebrew/lib */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <flint/arb.h>
#include <flint/fmpz.h>

static void rd(FILE *f, arb_t v) {
    char ms[8192]; long e; fmpz_t m; fmpz_init(m);
    if (fscanf(f, "%8191s %ld", ms, &e) != 2) { fprintf(stderr, "read error\n"); exit(2); }
    fmpz_set_str(m, ms, 16); arb_set_fmpz(v, m); arb_mul_2exp_si(v, v, e); fmpz_clear(m);
}
static double l10(const arb_t v) { arb_t t; arb_init(t); arb_abs(t, v); double r;
    if (arb_is_zero(t)) r = -1e9; else { arb_log_base_ui(t, t, 10, 64); r = arf_get_d(arb_midref(t), ARF_RND_NEAR); } arb_clear(t); return r; }

int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "r"); slong prec = atol(argv[2]); int k = atoi(argv[3]), its = atoi(argv[4]);
    long n, ps; if (fscanf(f, "%ld %ld", &n, &ps) != 2) return 2;
    arb_ptr d = _arb_vec_init(n), p = _arb_vec_init(n), a = _arb_vec_init(n), u = _arb_vec_init(n), x = _arb_vec_init(n);
    for (long i = 0; i < n; i++) { rd(f, d + i); rd(f, p + i); rd(f, a + i); rd(f, u + i); rd(f, x + i); }
    fclose(f);
    clock_t t0 = clock();
    slong tot = n*(n + 1)/2; arb_ptr A = _arb_vec_init(tot); arb_ptr D = _arb_vec_init(n);
    arb_t s, t, q; arb_init(s); arb_init(t); arb_init(q);
    #define ROW(i) (A + (slong)(i)*((i) + 1)/2)
    for (long i = 0; i < n; i++) {
        arb_ptr r = ROW(i);
        for (long j = 0; j < i; j++) {
            arb_mul(s, p + i, p + j, prec); arb_mul_si(s, s, ps, prec);
            arb_sub(t, u + j, u + i, prec); arb_sub(q, x + i, x + j, prec); arb_div(t, t, q, prec);
            arb_mul(t, t, a + i, prec); arb_mul(t, t, a + j, prec); arb_add(r + j, s, t, prec); mag_zero(arb_radref(r + j));
        }
        arb_set(r + i, d + i); mag_zero(arb_radref(r + i));
    }
    clock_t t1 = clock(); long neg = 0;
    for (long i = 0; i < n; i++) {
        arb_ptr r = ROW(i);
        for (long j = 1; j < i; j++) { arb_approx_dot(s, r + j, 1, r, 1, ROW(j), 1, j, prec); arb_swap(r + j, s); }
        arb_set(s, r + i);
        for (long kk = 0; kk < i; kk++) {      /* r[kk] holds W_ik = L_ik D_k; convert to L_ik, accumulate D_i */
            arb_set(t, r + kk); arb_div(r + kk, t, D + kk, prec); mag_zero(arb_radref(r + kk));
            arb_mul(q, t, r + kk, prec); arb_sub(s, s, q, prec);
        }
        arb_set(D + i, s); mag_zero(arb_radref(D + i));
        if (arf_sgn(arb_midref(D + i)) <= 0) neg++;
    }
    clock_t t2 = clock();
    arb_t mn; arb_init(mn); arb_set(mn, D + 0); for (long i = 1; i < n; i++) if (arf_cmp(arb_midref(D + i), arb_midref(mn)) < 0) arb_set(mn, D + i);
    printf("n %ld neg %ld minpiv_log10 %.3f fill_s %.1f ldl_s %.1f\n", n, neg, l10(mn), (double)(t1 - t0)/CLOCKS_PER_SEC, (double)(t2 - t1)/CLOCKS_PER_SEC);
    arb_ptr *V = malloc(k*sizeof(arb_ptr)), *Z = malloc(k*sizeof(arb_ptr));
    flint_rand_t st; flint_rand_init(st);
    for (int c = 0; c < k; c++) { V[c] = _arb_vec_init(n); Z[c] = _arb_vec_init(n);
        for (long i = 0; i < n; i++) arb_set_si(V[c] + i, (long)(n_randint(st, 2000001)) - 1000000); }
    for (int it = 0; it < its; it++) {
        for (int c = 0; c < k; c++) {      /* V[c] <- T^{-1} V[c] */
            arb_ptr v = V[c];
            for (long i = 1; i < n; i++) { arb_approx_dot(s, v + i, 1, ROW(i), 1, v, 1, i, prec); arb_swap(v + i, s); }
            for (long i = 0; i < n; i++) { arb_div(v + i, v + i, D + i, prec); mag_zero(arb_radref(v + i)); }
            for (long kk = n - 1; kk >= 1; kk--) { arb_ptr r = ROW(kk);
                for (long i = 0; i < kk; i++) { arb_submul(v + i, r + i, v + kk, prec); mag_zero(arb_radref(v + i)); } }
        }
        for (int c = 0; c < k; c++) {      /* Gram-Schmidt */
            for (int b = 0; b < c; b++) { arb_approx_dot(s, NULL, 0, V[c], 1, V[b], 1, n, prec);
                for (long i = 0; i < n; i++) { arb_submul(V[c] + i, s, V[b] + i, prec); mag_zero(arb_radref(V[c] + i)); } }
            arb_approx_dot(s, NULL, 0, V[c], 1, V[c], 1, n, prec); arb_sqrt(s, s, prec);
            for (long i = 0; i < n; i++) { arb_div(V[c] + i, V[c] + i, s, prec); mag_zero(arb_radref(V[c] + i)); }
        }
        for (int c = 0; c < k; c++) {      /* Z = L^T V */
            _arb_vec_zero(Z[c], n);
            for (long i = 0; i < n; i++) { arb_ptr r = ROW(i);
                for (long kk = 0; kk < i; kk++) arb_addmul(Z[c] + kk, r + kk, V[c] + i, prec);
                arb_add(Z[c] + i, Z[c] + i, V[c] + i, prec); }
            for (long i = 0; i < n; i++) mag_zero(arb_radref(Z[c] + i));
        }
        printf("it %d", it);
        for (int c = 0; c < k; c++) { arb_zero(s); for (long i = 0; i < n; i++) { arb_mul(t, Z[c] + i, Z[c] + i, prec); arb_addmul(s, t, D + i, prec); }
            printf(" %.9f", l10(s)); }
        printf("\n");
    }
    printf("H");
    for (int b = 0; b < k; b++) for (int c = 0; c < k; c++) {
        arb_zero(s); for (long i = 0; i < n; i++) { arb_mul(t, Z[b] + i, Z[c] + i, prec); arb_addmul(s, t, D + i, prec); }
        char *str = arb_get_str(s, 40, ARB_STR_NO_RADIUS); printf(" %s", str); flint_free(str); }
    printf("\nsec %.1f\n", (double)(clock() - t0)/CLOCKS_PER_SEC);
    if (argc > 5) { FILE *g = fopen(argv[5], "w"); fmpz_t m, e; fmpz_init(m); fmpz_init(e);
        fprintf(g, "%ld %d\n", n, k);
        for (long i = 0; i < n; i++) { for (int cc = 0; cc < k; cc++) {
            arf_get_fmpz_2exp(m, e, arb_midref(V[cc] + i)); char *ms = fmpz_get_str(NULL, 16, m); fprintf(g, " %s %ld", ms, fmpz_get_si(e)); flint_free(ms); }
            fprintf(g, "\n"); }
        fclose(g); fmpz_clear(m); fmpz_clear(e); }
    return 0;
}
