/* Least prime quadratic, cubic, quintic non-residues for all primes p <= N (default 1e7).
   For r in {2,3,5}: over primes p with r | p-1, the least n with n^((p-1)/r) != 1 mod p (n is then prime).
   Prints the maxima, the worst p, and the largest ratio log n_r(p) / log((log p)^{7.8978}).
   Build: cc -O2 -o /tmp/nr nr.c -lm ; run: /tmp/nr 10000000 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
static uint64_t pw(uint64_t b, uint64_t e, uint64_t m){ uint64_t r=1%m; b%=m; while(e){ if(e&1) r=(r*b)%m; b=(b*b)%m; e>>=1;} return r; }
int main(int argc, char**argv){
  long N = argc>1 ? atol(argv[1]) : 10000000L;
  unsigned char *c = calloc(N+1,1); if(!c) return 1;
  for(long i=2;(long long)i*i<=N;i++) if(!c[i]) for(long j=(long)i*i;j<=N;j+=i) c[j]=1;
  int rs[3]={2,3,5}; const double E=7.8978;
  for(int k=0;k<3;k++){
    int r=rs[k]; long cnt=0, worstp=0, maxn=0; double maxratio=-1; long ratiop=0, ration=0;
    for(long p=3;p<=N;p++){ if(c[p]) continue; if((p-1)%r) continue; if(p==r) continue;
      cnt++; uint64_t e=(p-1)/r; long n=2;
      for(;;n++){ if(pw(n,e,p)!=1) break; }
      if(n>maxn){maxn=n; worstp=p;}
      double bound = E*log(log((double)p));
      double ratio = (bound>0)? log((double)n)/bound : 99;
      if(ratio>maxratio && p>=17){maxratio=ratio; ratiop=p; ration=n;}
    }
    printf("r=%d primes p<=%ld with r|p-1: %ld; max least r-th non-residue %ld at p=%ld; bound (log p)^7.8978 there = %.3g; worst log n/log bound (p>=17) = %.4f at p=%ld (n=%ld)\n",
      r, N, cnt, maxn, worstp, pow(log((double)worstp),E), maxratio, ratiop, ration);
  }
  free(c); return 0;
}
