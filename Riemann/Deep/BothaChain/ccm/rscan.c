// rscan.c: rigorous pre-scan of mu(k) for cert_chi.py (no libm in the loop).
// For every node k in each block (k0, k0+2, ..., <= k1) it computes mu_c(k) = H - T - 2 sum_n W_n cos_c(theta_n(k)), where
//   theta_n(k) = 2 pi frac(k B_n / 2^64) (exact 64-bit integer phase), cos_c from a 2^16-entry table (C_a, S_a ~ cos, sin of 2 pi a/2^16,
//   correctly rounded from Arb) and the Taylor pair 1 - beta^2/2, beta (1 - beta^2/6) for the low 48 bits (beta < 2 pi 2^-16).
// H <= h(r_k) and T >= |tail(r_k)| for every k of the block are Arb bounds (h increasing, tail decreasing in r).
// Only IEEE-754 +, -, * (correctly rounded) are used; compile with -ffp-contract=off. The a-priori error bound E (RIGOR.md) gives
// mu(k) >= mu_c(k) - E. Prints every k with mu_c(k) < thr + margin (margin >= E is checked by the caller); those are re-checked in Arb.
// Input (stdin, hex floats): 65536 x "C S"; K; natoms; natoms x "B W"; thr; margin; nblocks; nblocks x "k0 k1 H T".
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
static double C[65536], S[65536];
int main(void) {
    for (int i = 0; i < 65536; i++) if (scanf("%la %la", &C[i], &S[i]) != 2) return 1;
    double K; if (scanf("%la", &K) != 1) return 1;
    int na; if (scanf("%d", &na) != 1 || na > 64) return 1;
    uint64_t B[64]; double W[64];
    for (int i = 0; i < na; i++) { unsigned long long b; if (scanf("%llu %la", &b, &W[i]) != 2) return 1; B[i] = (uint64_t)b; }
    double thr, mg; if (scanf("%la %la", &thr, &mg) != 2) return 1;
    long nb; if (scanf("%ld", &nb) != 1) return 1;
    const double c6 = 1.0/6.0;
    long long nflag = 0, nnode = 0;
    for (long ib = 0; ib < nb; ib++) {
        unsigned long long k0, k1; double H, T;
        if (scanf("%llu %llu %la %la", &k0, &k1, &H, &T) != 4) return 1;
        double HT = H - T;
        for (uint64_t k = k0; k <= k1; k += 2) {
            double s = 0.0;
            for (int n = 0; n < na; n++) {
                uint64_t ph = k*B[n];
                unsigned a = (unsigned)(ph >> 48);
                double beta = (double)(ph & 0xFFFFFFFFFFFFULL) * K;
                double t = beta*beta;
                double cb = 1.0 - 0.5*t;
                double sb = beta*(1.0 - t*c6);
                double cs = C[a]*cb - S[a]*sb;
                s += W[n]*cs;
            }
            double mu = HT - 2.0*s;
            nnode++;
            if (mu < thr + mg) { printf("%llu\n", (unsigned long long)k); nflag++; }
        }
    }
    fprintf(stderr, "nodes %lld flagged %lld\n", nnode, nflag);
    return 0;
}
