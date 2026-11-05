// cscan.c: float64 pre-scan of mu(k) for cert_chi.py (flags candidates; every flagged node is re-checked in Arb).
// Input (stdin, text): L a cst ntail {bn cbn}*ntail natoms {B(uint64) W}*natoms kstart kend thr
// Scans k = kstart, kstart+2, ..., <= kend; prints k with mu(k) < thr + 1e-6.  Phases exact mod 2^64; rotator recurrence with resync.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>
int main(void) {
    double L, cst, thr; int a, nt, na; unsigned long long ks, ke;
    if (scanf("%lf %d %lf %d", &L, &a, &cst, &nt) != 4) return 1;
    double bn[16], cbn[16];
    for (int i = 0; i < nt; i++) if (scanf("%lf %lf", &bn[i], &cbn[i]) != 2) return 1;
    if (scanf("%d", &na) != 1) return 1;
    uint64_t B[64]; double W[64];
    for (int i = 0; i < na; i++) if (scanf("%llu %lf", (unsigned long long*)&B[i], &W[i]) != 2) return 1;
    if (scanf("%llu %llu %lf", &ks, &ke, &thr) != 3) return 1;
    const double TWO_PI = 6.283185307179586476925286766559;
    double complex z[64], rho[64];
    for (int i = 0; i < na; i++) { uint64_t st = 2*B[i]; rho[i] = cexp(I*TWO_PI*((double)st*0x1p-64)); }
    long long nflag = 0; int cnt = 0;
    for (unsigned long long k = ks; k <= ke; k += 2) {
        if (cnt == 0) { for (int i = 0; i < na; i++) { uint64_t ph = (uint64_t)k*B[i]; z[i] = cexp(I*TWO_PI*((double)ph*0x1p-64)); } cnt = 1024; }
        double r = M_PI*(double)k/L;
        double complex w = (0.25 + 0.5*a) + I*0.5*r, u = 1.0/w, u2 = u*u;
        double complex ps = clog(w) - 0.5*u - u2*(1.0/12 - u2*(1.0/120 - u2*(1.0/252 - u2*(1.0/240 - u2*(1.0/132)))));
        double tl = 0; for (int i = 0; i < nt; i++) tl += cbn[i]*bn[i]/(bn[i]*bn[i] + r*r);
        double sg = (k % 2 == 0) ? 1.0 : -1.0;
        double s = 0; for (int i = 0; i < na; i++) { s += W[i]*creal(z[i]); z[i] *= rho[i]; }
        cnt--;
        double mu = creal(ps) + cst + 2*sg*tl - 2*s;
        if (mu < thr + 1e-6) { printf("%llu\n", k); nflag++; }
    }
    fprintf(stderr, "flagged %lld\n", nflag);
    return 0;
}
