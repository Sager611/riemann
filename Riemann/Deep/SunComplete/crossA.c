
/* crossA.c : Hamilton cycles in H_q (x~y iff x!=y, x+y primitive) for all NON-PRIME prime powers 8 <= q < QMAX,
   by a seeded Posa rotation-extension search; every cycle is re-checked in full. Field arithmetic as in crossB.c. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int P, K, Q, N; static int *EXPT, *LOGT; static unsigned char *ISPRIM; static int PW[32];
static int gcdi(int a, int b){ while (b){ int t = a % b; a = b; b = t; } return a < 0 ? -a : a; }
static int addf(int a, int b){
  if (K == 1){ int s = a + b; return s >= P ? s - P : s; }
  int r = 0; for (int i = 0; i < K; i++){ int s = a % P + b % P; if (s >= P) s -= P; r += s*PW[i]; a /= P; b /= P; } return r; }
static int mulx(int a, const int *f){ int d[32], nd[32]; for (int i = 0; i < K; i++){ d[i] = a % P; a /= P; }
  int top = d[K-1]; nd[0] = 0; for (int i = 1; i < K; i++) nd[i] = d[i-1];
  int r = 0; for (int i = 0; i < K; i++){ int v = (nd[i] - top*f[i]) % P; if (v < 0) v += P; r += v*PW[i]; } return r; }
static int setup_field(int q){
  int p = 0; for (int d = 2; d <= q; d++) if (q % d == 0){ p = d; break; }
  int k = 0, x = q; while (x % p == 0){ x /= p; k++; } if (x != 1) return 0;
  P = p; K = k; Q = q; N = q - 1; PW[0] = 1; for (int i = 1; i <= K; i++) PW[i] = PW[i-1]*P;
  EXPT = realloc(EXPT, sizeof(int)*N); LOGT = realloc(LOGT, sizeof(int)*q); ISPRIM = realloc(ISPRIM, q);
  int f[32]; int found = 0; long total = 1; for (int i = 0; i < K; i++) total *= P;
  for (long c = 0; c < total && !found; c++){ long cc = c; for (int i = 0; i < K; i++){ f[i] = cc % P; cc /= P; }
    if (f[0] == 0) continue; for (int z = 0; z < q; z++) LOGT[z] = -1;
    int a = 1, ok = 1; for (int i = 0; i < N; i++){ if (a == 0 || LOGT[a] != -1){ ok = 0; break; } LOGT[a] = i; EXPT[i] = a; a = mulx(a, f); }
    if (ok && a == 1) found = 1; }
  if (!found) return 0; LOGT[0] = -1;
  for (int z = 0; z < q; z++) ISPRIM[z] = (z != 0 && gcdi(LOGT[z], N) == 1); return 1; }
static uint64_t rs = 88172645463325252ULL;
static uint64_t rnd(void){ rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static int adj(int a, int b){ return a != b && ISPRIM[addf(a, b)]; }
int main(int argc, char **argv){
  int QMAX = argc > 1 ? atoi(argv[1]) : 30140;
  int *path = malloc(sizeof(int)*QMAX), *unv = malloc(sizeof(int)*QMAX);
  for (int q = 8; q < QMAX; q++){
    int p = 0; for (int d = 2; d*d <= q; d++) if (q % d == 0){ p = d; break; } if (!p) continue; /* skip primes */
    int x = q; while (x % p == 0) x /= p; if (x != 1) continue;
    if (!setup_field(q)){ printf("SETUP FAIL %d\n", q); return 1; }
    rs = 88172645463325252ULL ^ (uint64_t)q * 0x9E3779B97F4A7C15ULL; int ok = 0, attempts = 0; long steps = 0;
    while (!ok && attempts < 1000){
      attempts++;
      int len = 1; path[0] = (int)(rnd() % q); int nu = 0; for (int v = 0; v < q; v++) if (v != path[0]) unv[nu++] = v;
      long lim = 200L*q + 10000; long st = 0;
      while (st++ < lim){
        int end = path[len-1];
        if (nu > 0){
          int off = (int)(rnd() % nu), hit = -1;
          for (int t = 0; t < nu; t++){ int idx = off + t; if (idx >= nu) idx -= nu; if (adj(end, unv[idx])){ hit = idx; break; } }
          if (hit >= 0){ path[len++] = unv[hit]; unv[hit] = unv[--nu]; continue; }
        } else if (adj(end, path[0])){ ok = 1; break; }
        if (len < 3){ break; }
        int off = (int)(rnd() % (len-2)), piv = -1;
        for (int t = 0; t < len-2; t++){ int i = off + t; if (i >= len-2) i -= len-2; if (adj(end, path[i])){ piv = i; break; } }
        if (piv < 0 || (rnd() % 10) == 0){ for (int i = 0, j = len-1; i < j; i++, j--){ int tmp = path[i]; path[i] = path[j]; path[j] = tmp; } continue; }
        for (int i = piv+1, j = len-1; i < j; i++, j--){ int tmp = path[i]; path[i] = path[j]; path[j] = tmp; }
      }
      steps += st;
    }
    if (!ok){ printf("q=%d NOT FOUND\n", q); continue; }
    /* full re-check */
    unsigned char *seen = calloc(q, 1); int good = 1;
    for (int i = 0; i < q; i++){ if (seen[path[i]]) good = 0; seen[path[i]] = 1; if (!ISPRIM[addf(path[i], path[(i+1)%q])]) good = 0; }
    free(seen);
    uint64_t h = 1469598103934665603ULL; for (int i = 0; i < q; i++){ h ^= (uint64_t)path[i]; h *= 1099511628211ULL; }
    printf("q=%d p=%d k=%d cycle_ok=%d attempts=%d steps=%ld hash=%016llx\n", q, P, K, good, attempts, steps, (unsigned long long)h);
    fflush(stdout);
  }
  return 0;
}

