/* fxkern.c -- dense kernels on fixed-point blocks (weilooc, 2026-11-05), called from Python via ctypes.
   Blocks are row-major arrays of ND uint16 little-endian digits per entry, offset binary (X = Y - 2^(T-1)), value X/2^F.
   Midpoint ("approx") arithmetic only: these kernels produce DATA (a Cholesky factor, a triangular solve); the proof
   never trusts them -- the certificate recomputes the residual H - sigma I - L L^T exactly (oocert.py).
   fx_chol : L lower triangular with L L^T ~ A (Cholesky-Banachiewicz, arb_approx_dot at prec bits); upper part 0.
             returns 0, or -(k+1) if pivot k is not positive.
   fx_trsm : X with X L^T ~ S (rows of S, L lower triangular).
   fx_trsm_left : X with L X ~ B (forward substitution, B n x m).
   fx_trsm_left_t : X with L^T X ~ B (back substitution).
   Build: cc -O2 -shared -fPIC -o libfxkern.dylib fxkern.c -I/opt/homebrew/include -L/opt/homebrew/lib -lflint -lgmp -lmpfr */
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/arb.h>

static void dig_to_arb(arb_t x, const uint16_t *d, int ND, long F, long T) {
    int nl = (ND + 3)/4; ulong limbs[64];
    memset(limbs, 0, sizeof limbs);
    for (int s = 0; s < ND; s++) limbs[s/4] |= ((ulong) d[s]) << (16*(s % 4));
    fmpz_t X, off; fmpz_init(X); fmpz_init(off);
    fmpz_set_ui_array(X, limbs, nl);
    fmpz_one(off); fmpz_mul_2exp(off, off, T - 1); fmpz_sub(X, X, off);
    arb_set_fmpz(x, X); arb_mul_2exp_si(x, x, -F);
    fmpz_clear(X); fmpz_clear(off);
}
static int arb_to_dig(uint16_t *d, const arb_t x, int ND, long F, long T) {
    arf_t t; arf_init(t); fmpz_t X, off; fmpz_init(X); fmpz_init(off);
    arf_mul_2exp_si(t, arb_midref(x), F); arf_get_fmpz(X, t, ARF_RND_NEAR);
    fmpz_one(off); fmpz_mul_2exp(off, off, T - 1); fmpz_add(X, X, off);
    int bad = (fmpz_sgn(X) < 0 || fmpz_bits(X) > (ulong) T);
    if (!bad) {
        ulong limbs[64]; memset(limbs, 0, sizeof limbs); int nl = (ND + 3)/4;
        fmpz_get_ui_array(limbs, nl, X);
        for (int s = 0; s < ND; s++) d[s] = (uint16_t)((limbs[s/4] >> (16*(s % 4))) & 0xFFFF);
    }
    arf_clear(t); fmpz_clear(X); fmpz_clear(off);
    return bad ? -1 : 0;
}
static arb_ptr load(const uint16_t *src, long r, long c, int ND, long F, long T) {
    arb_ptr a = _arb_vec_init(r*c);
    for (long i = 0; i < r*c; i++) dig_to_arb(a + i, src + (size_t) i*ND, ND, F, T);
    return a;
}
static int store(uint16_t *dst, arb_srcptr a, long r, long c, int ND, long F, long T) {
    for (long i = 0; i < r*c; i++) if (arb_to_dig(dst + (size_t) i*ND, a + i, ND, F, T)) return -1;
    return 0;
}

int fx_chol(long n, int ND, long F, long T, const uint16_t *in, uint16_t *out, long prec) {
    arb_ptr A = load(in, n, n, ND, F, T);
    arb_ptr L = _arb_vec_init(n*n);
    arb_t s; arb_init(s);
    int ret = 0;
    for (long i = 0; i < n && !ret; i++) {
        for (long j = 0; j <= i; j++) {
            arb_approx_dot(s, A + i*n + j, 1, L + i*n, 1, L + j*n, 1, j, prec);
            if (i == j) {
                if (!arf_sgn(arb_midref(s)) || arf_sgn(arb_midref(s)) < 0) { ret = -(int)(i + 1); break; }
                arf_sqrt(arb_midref(L + i*n + i), arb_midref(s), prec, ARF_RND_NEAR);
                mag_zero(arb_radref(L + i*n + i));
            } else {
                arf_div(arb_midref(L + i*n + j), arb_midref(s), arb_midref(L + j*n + j), prec, ARF_RND_NEAR);
                mag_zero(arb_radref(L + i*n + j));
            }
        }
    }
    if (!ret && store(out, L, n, n, ND, F, T)) ret = -1000000;
    arb_clear(s); _arb_vec_clear(A, n*n); _arb_vec_clear(L, n*n);
    return ret;
}

int fx_trsm(long m, long n, int ND, long F, long T, const uint16_t *S, const uint16_t *Ld, uint16_t *out, long prec) {
    arb_ptr Sa = load(S, m, n, ND, F, T), L = load(Ld, n, n, ND, F, T), X = _arb_vec_init(m*n);
    arb_t s; arb_init(s);
    for (long r = 0; r < m; r++)
        for (long j = 0; j < n; j++) {
            arb_approx_dot(s, Sa + r*n + j, 1, X + r*n, 1, L + j*n, 1, j, prec);
            arf_div(arb_midref(X + r*n + j), arb_midref(s), arb_midref(L + j*n + j), prec, ARF_RND_NEAR);
            mag_zero(arb_radref(X + r*n + j));
        }
    int ret = store(out, X, m, n, ND, F, T);
    arb_clear(s); _arb_vec_clear(Sa, m*n); _arb_vec_clear(L, n*n); _arb_vec_clear(X, m*n);
    return ret;
}

/* L X = B, L n x n lower, B n x m (row-major); column-oriented forward substitution on X^T rows */
int fx_trsm_left(long n, long m, int ND, long F, long T, const uint16_t *Ld, const uint16_t *B, uint16_t *out, long prec) {
    arb_ptr L = load(Ld, n, n, ND, F, T), Ba = load(B, n, m, ND, F, T);
    arb_ptr Xt = _arb_vec_init(m*n), Bt = _arb_vec_init(m*n);
    for (long i = 0; i < n; i++) for (long c = 0; c < m; c++) arb_swap(Bt + c*n + i, Ba + i*m + c);
    arb_t s; arb_init(s);
    for (long c = 0; c < m; c++)
        for (long i = 0; i < n; i++) {
            arb_approx_dot(s, Bt + c*n + i, 1, L + i*n, 1, Xt + c*n, 1, i, prec);
            arf_div(arb_midref(Xt + c*n + i), arb_midref(s), arb_midref(L + i*n + i), prec, ARF_RND_NEAR);
            mag_zero(arb_radref(Xt + c*n + i));
        }
    for (long i = 0; i < n; i++) for (long c = 0; c < m; c++) arb_swap(Ba + i*m + c, Xt + c*n + i);
    int ret = store(out, Ba, n, m, ND, F, T);
    arb_clear(s); _arb_vec_clear(L, n*n); _arb_vec_clear(Ba, n*m); _arb_vec_clear(Xt, m*n); _arb_vec_clear(Bt, m*n);
    return ret;
}

/* L^T X = B (back substitution), L n x n lower, B n x m */
int fx_trsm_left_t(long n, long m, int ND, long F, long T, const uint16_t *Ld, const uint16_t *B, uint16_t *out, long prec) {
    arb_ptr L = load(Ld, n, n, ND, F, T), Ba = load(B, n, m, ND, F, T);
    arb_ptr U = _arb_vec_init(n*n);                         /* U = L^T stored row-major, reversed index order */
    for (long i = 0; i < n; i++) for (long j = 0; j <= i; j++) arb_set(U + (n - 1 - j)*n + (n - 1 - i), L + i*n + j);
    /* with reversed indices, L^T becomes lower triangular: reuse forward substitution */
    arb_ptr Xt = _arb_vec_init(m*n), Bt = _arb_vec_init(m*n);
    for (long i = 0; i < n; i++) for (long c = 0; c < m; c++) arb_set(Bt + c*n + (n - 1 - i), Ba + i*m + c);
    arb_t s; arb_init(s);
    for (long c = 0; c < m; c++)
        for (long i = 0; i < n; i++) {
            arb_approx_dot(s, Bt + c*n + i, 1, U + i*n, 1, Xt + c*n, 1, i, prec);
            arf_div(arb_midref(Xt + c*n + i), arb_midref(s), arb_midref(U + i*n + i), prec, ARF_RND_NEAR);
            mag_zero(arb_radref(Xt + c*n + i));
        }
    for (long i = 0; i < n; i++) for (long c = 0; c < m; c++) arb_set(Ba + i*m + c, Xt + c*n + (n - 1 - i));
    int ret = store(out, Ba, n, m, ND, F, T);
    arb_clear(s); _arb_vec_clear(L, n*n); _arb_vec_clear(U, n*n); _arb_vec_clear(Ba, n*m); _arb_vec_clear(Xt, m*n); _arb_vec_clear(Bt, m*n);
    return ret;
}
