
/* crossB.c : independent cross-check for ALL prime powers 8 <= q < QMAX.
   For every a0 in F_q (incl. 0) search a certificate for the constructions of NOTE.md sec. 5:
     q odd : (u,r)  u non-square, r generator of the squares, a0+u, a0+ur primitive
     q even: (u,s)  s primitive, a0+u, a0+us, a0+us^2 primitive
   and for (A), q prime: smallest x with x, x+1, x+2 primitive (She's construction).
   Field arithmetic: F_p[x]/(f), f monic of degree k with x of multiplicative order q-1
   (checked by building the full power table), elements coded base p.
   Output per q: counts, failures, FNV-1a checksum of the certificate list. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int P, K, Q, N;           /* q = P^K, N = q-1 */
static int *EXPT, *LOGT;          /* EXPT[i] = x^i, LOGT[z] (LOGT[0] = -1) */
static unsigned char *ISPRIM;
static int PW[32];

static int gcdi(int a, int b){ while (b){ int t = a % b; a = b; b = t; } return a < 0 ? -a : a; }

static int addf(int a, int b){
  if (K == 1){ int s = a + b; return s >= P ? s - P : s; }
  int r = 0;
  for (int i = 0; i < K; i++){ int da = a % P, db = b % P; int s = da + db; if (s >= P) s -= P; r += s*PW[i]; a /= P; b /= P; }
  return r;
}
static int negf(int a){
  if (K == 1) return a ? P - a : 0;
  int r = 0;
  for (int i = 0; i < K; i++){ int d = a % P; r += (d ? P - d : 0)*PW[i]; a /= P; }
  return r;
}
static int mulf(int a, int b){ if (!a || !b) return 0; int e = LOGT[a] + LOGT[b]; if (e >= N) e -= N; return EXPT[e]; }

static int mulx(int a, const int *f){ /* a*x mod f */
  int d[32]; for (int i = 0; i < K; i++){ d[i] = a % P; a /= P; }
  int top = d[K-1]; int nd[32]; nd[0] = 0; for (int i = 1; i < K; i++) nd[i] = d[i-1];
  int r = 0; for (int i = 0; i < K; i++){ int v = (nd[i] - top*f[i]) % P; if (v < 0) v += P; r += v*PW[i]; }
  return r;
}

static int setup_field(int q){
  /* returns 1 on success; finds p,k */
  int p = 0; for (int d = 2; d <= q; d++) if (q % d == 0){ p = d; break; }
  int k = 0, x = q; while (x % p == 0){ x /= p; k++; } if (x != 1) return 0;
  P = p; K = k; Q = q; N = q - 1;
  PW[0] = 1; for (int i = 1; i <= K; i++) PW[i] = PW[i-1]*P;
  EXPT = realloc(EXPT, sizeof(int)*N); LOGT = realloc(LOGT, sizeof(int)*q); ISPRIM = realloc(ISPRIM, q);
  if (K == 1){
    for (int g = 1; g < P; g++){
      for (int z = 0; z < q; z++) LOGT[z] = -1;
      int a = 1, ok = 1;
      for (int i = 0; i < N; i++){ if (LOGT[a] != -1){ ok = 0; break; } LOGT[a] = i; EXPT[i] = a; a = (int)((long)a*g % P); }
      if (ok && a == 1) break;
    }
  } else {
    int f[32]; int found = 0;
    long total = 1; for (int i = 0; i < K; i++) total *= P;
    for (long c = 0; c < total && !found; c++){
      long cc = c; for (int i = 0; i < K; i++){ f[i] = cc % P; cc /= P; }
      if (f[0] == 0) continue;
      for (int z = 0; z < q; z++) LOGT[z] = -1;
      int a = 1, ok = 1;
      for (int i = 0; i < N; i++){ if (LOGT[a] != -1 || a == 0){ ok = 0; break; } LOGT[a] = i; EXPT[i] = a; a = mulx(a, f); }
      if (ok && a == 1) found = 1;
    }
    if (!found) return 0;
  }
  LOGT[0] = -1;
  for (int z = 0; z < q; z++) ISPRIM[z] = (z != 0 && gcdi(LOGT[z], N) == 1);
  return 1;
}

static uint64_t fnv(uint64_t h, uint64_t v){ for (int i = 0; i < 8; i++){ h ^= (v >> (8*i)) & 0xff; h *= 1099511628211ULL; } return h; }

int main(int argc, char **argv){
  int QMAX = argc > 1 ? atoi(argv[1]) : 30140;
  int *prims = malloc(sizeof(int)*QMAX);
  long tot_inst = 0, tot_fail = 0;
  for (int q = 8; q < QMAX; q++){
    /* prime power? */
    int p = 0; for (int d = 2; d*d <= q; d++) if (q % d == 0){ p = d; break; } if (!p) p = q;
    int x = q; while (x % p == 0) x /= p; if (x != 1) continue;
    if (!setup_field(q)){ printf("SETUP FAIL %d\n", q); return 1; }
    int np = 0; for (int z = 0; z < q; z++) if (ISPRIM[z]) prims[np++] = z;
    /* (A) She's x for primes */
    int sheX = -1;
    if (K == 1 && (q & 1)) for (int z = 0; z < q; z++) if (ISPRIM[z] && ISPRIM[(z+1)%q] && ISPRIM[(z+2)%q]){ sheX = z; break; }
    uint64_t h = 1469598103934665603ULL; int fails = 0; int failA0[64]; int nfa = 0;
    for (int a0 = 0; a0 < q; a0++){
      int cu = -1, cr = -1;
      if (q & 1){
        int m = N/2;
        for (int u = 1; u < q && cu < 0; u++){
          if (!(LOGT[u] & 1) || !ISPRIM[addf(a0, u)]) continue;
          int lu = LOGT[u];
          int na0 = negf(a0);
          for (int i = 0; i < np; i++){
            int v = addf(prims[i], na0);
            if (!v || !(LOGT[v] & 1)) continue;
            int d = LOGT[v] - lu; if (d < 0) d += N;
            if (gcdi(d/2, m) == 1){ cu = u; cr = EXPT[d]; break; }
          }
        }
      } else {
        for (int u = 1; u < q && cu < 0; u++){
          if (!ISPRIM[addf(a0, u)]) continue;
          for (int i = 0; i < np; i++){ int s = prims[i];
            int us = mulf(u, s), us2 = mulf(us, s);
            if (ISPRIM[addf(a0, us)] && ISPRIM[addf(a0, us2)]){ cu = u; cr = s; break; }
          }
        }
      }
      tot_inst++;
      if (cu < 0){ fails++; tot_fail++; if (nfa < 64) failA0[nfa++] = a0; h = fnv(h, (uint64_t)a0 << 40); }
      else h = fnv(h, ((uint64_t)a0 << 40) | ((uint64_t)cu << 20) | (uint64_t)cr);
    }
    printf("q=%d p=%d k=%d sheX=%d failsB=%d hash=%016llx", q, P, K, sheX, fails, (unsigned long long)h);
    if (fails){ printf(" fail_a0="); for (int i = 0; i < nfa; i++) printf("%d,", failA0[i]); }
    printf("\n");
  }
  fprintf(stderr, "instances=%ld failures=%ld\n", tot_inst, tot_fail);
  return 0;
}

