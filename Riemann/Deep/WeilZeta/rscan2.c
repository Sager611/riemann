// rscan2.c: rscan.c (RIGOR.md) with a proved early exit, for large scans (zeta at c = 47: 7e10 nodes).
// Atoms are processed in order of decreasing |W|. After atom n the code knows the partial sum s_n and R_{n+1} >= sum_{m>n} |W_m|
// (suffix sums, rounded upward by a relative 1e-12 safety). Since cos <= 1, mu(k) >= H - T - 2 (sum_{m<=n} w_m cos theta_m + sum_{m>n} |w_m|),
// so if mu_p = fl(HT - 2 fl(s_n + R_{n+1})) >= thr + margin the node is rejected (provably mu(k) >= thr: the error of mu_p is bounded by the
// same E as the full sum plus one extra rounding 2 u sigma_W, added by the caller). Nodes never rejected get exactly the full rscan.c
// evaluation and test. Only IEEE-754 +, -, * (and one exact comparison per atom); compile with -O2 -ffp-contract=off.
// Input format identical to rscan.c. Output: flagged k (stdout); "nodes N flagged F atoms_evaluated A" (stderr).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
static double C[65536], S[65536];
int main(void) {
    for (int i = 0; i < 65536; i++) if (scanf("%la %la", &C[i], &S[i]) != 2) return 1;
    double K; if (scanf("%la", &K) != 1) return 1;
    int na; if (scanf("%d", &na) != 1 || na > 64) return 1;
    uint64_t B[64]; double W[64];
    for (int i = 0; i < na; i++) { unsigned long long b; if (scanf("%llu %la", &b, &W[i]) != 2) return 1; B[i] = (uint64_t)b; }
    for (int i = 0; i < na; i++) for (int j = i + 1; j < na; j++) if (fabs(W[j]) > fabs(W[i])) {   // sort by |W| descending
        double tw = W[i]; W[i] = W[j]; W[j] = tw; uint64_t tb = B[i]; B[i] = B[j]; B[j] = tb; }
    double R[65]; R[na] = 0.0;
    for (int i = na - 1; i >= 0; i--) R[i] = R[i+1] + fabs(W[i]);
    for (int i = 0; i < na; i++) R[i] = R[i] * (1.0 + 1e-12) + 1e-300;   // upward safety (summation error <= na u << 1e-12)
    double thr, mg; if (scanf("%la %la", &thr, &mg) != 2) return 1;
    long nb; if (scanf("%ld", &nb) != 1) return 1;
    const double c6 = 1.0/6.0;
    const double lim = thr + mg;
    long long nflag = 0, nnode = 0, natom = 0;
    for (long ib = 0; ib < nb; ib++) {
        unsigned long long k0, k1; double H, T;
        if (scanf("%llu %llu %la %la", &k0, &k1, &H, &T) != 4) return 1;
        double HT = H - T;
        for (uint64_t k = k0; k <= k1; k += 2) {
            double s = 0.0; int n = 0, rej = 0;
            for (n = 0; n < na; n++) {
                uint64_t ph = k*B[n];
                unsigned a = (unsigned)(ph >> 48);
                double beta = (double)(ph & 0xFFFFFFFFFFFFULL) * K;
                double t = beta*beta;
                double cb = 1.0 - 0.5*t;
                double sb = beta*(1.0 - t*c6);
                double cs = C[a]*cb - S[a]*sb;
                s += W[n]*cs;
                if (n + 1 < na && HT - 2.0*(s + R[n+1]) >= lim) { rej = 1; n++; break; }
            }
            natom += n;
            nnode++;
            if (rej) continue;
            double mu = HT - 2.0*s;
            if (mu < lim) { printf("%llu\n", (unsigned long long)k); nflag++; }
        }
    }
    fprintf(stderr, "nodes %lld flagged %lld atoms_evaluated %lld\n", nnode, nflag, natom);
    return 0;
}
