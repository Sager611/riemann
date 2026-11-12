
// fwht.c — exact Walsh spectrum of knapsack down-sets, n <= 24, int32, single-threaded.
// usage: fwht <system> <x> <seed> <eps>
// system: 0 primes (weights log p, threshold log x, exact integer test)
//         1 equal weights (ball radius r + generic initial segment of level r+1, |D| = Q(x))
//         2 iid perturbed: log p + U[-eps,eps]
//         3 Poisson points of intensity 1/log t on [2,inf): first n=pi(x) points, threshold log x
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
static int n; static double w[32]; static double th; static int32_t *F; static long cnt;
static uint64_t rs=88172645463325252ULL; static double urand(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return (rs>>11)*(1.0/9007199254740992.0);}
static int isprime(int k){ if(k<2) return 0; for(int d=2;d*d<=k;d++) if(k%d==0) return 0; return 1; }
static int ip[32]; static int sysid; static long X;
static void dfs(int start, double sum, long prod, uint32_t mask){
  F[mask]=1; cnt++;
  for(int i=start;i<n;i++){
    if(sysid==0){ if(prod*ip[i] > X) break; dfs(i+1,0,prod*ip[i],mask|(1u<<i)); }
    else { if(sum+w[i] > th+1e-12) break; dfs(i+1,sum+w[i],0,mask|(1u<<i)); }
  }
}
static int cmpd(const void*a,const void*b){ double x=*(double*)a,y=*(double*)b; return x<y?-1:x>y; }
typedef struct {double k; uint32_t m;} KP;
static int cmpP(const void*a,const void*b){ double x=((const KP*)a)->k,y=((const KP*)b)->k; return x<y?-1:x>y; }
static long Qcount(long x){ long c=0; for(long m=1;m<=x;m++){ int sf=1; for(long d=2;d*d<=m;d++) if(m%(d*d)==0){sf=0;break;} c+=sf;} return c; }
static long binom(int a,int b){ if(b<0||b>a) return 0; long r=1; for(int i=1;i<=b;i++) r=r*(a-b+i)/i; return r; }
int main(int argc,char**argv){
  sysid=atoi(argv[1]); X=atol(argv[2]); rs ^= (uint64_t)atol(argv[3])*0x9E3779B97F4A7C15ULL; for(int i=0;i<5;i++) urand();
  double eps=atof(argv[4]);
  n=0; for(int k=2;k<=X;k++) if(isprime(k)) ip[n++]=k;
  if(n>24){ fprintf(stderr,"n>24\n"); return 1; }
  long Q=Qcount(X);
  th=log((double)X);
  if(sysid==2){ for(int i=0;i<n;i++) w[i]=log((double)ip[i])+eps*(2*urand()-1); qsort(w,n,sizeof(double),cmpd); }
  if(sysid==3){ double t=2.0; int i=0; while(i<n){ double lam=1.0/log(t+0.5); if(urand()<lam*0.01){ w[i++]=log(t);} t+=0.01; } }
  size_t N=(size_t)1<<n; F=calloc(N,sizeof(int32_t)); if(!F){fprintf(stderr,"oom\n");return 1;}
  if(sysid==1){ // equal weights: levels 0..r full, then first rem sets of level r+1 by generic order
    long tot=0; int r=0; while(tot+binom(n,r)<=Q){ tot+=binom(n,r); r++; } long rem=Q-tot; r--; // levels 0..r full
    double e[32]; for(int i=0;i<n;i++) e[i]=urand();
    // collect level r+1 sets sorted by sum e: brute force over masks with popcount r+1
    long L=binom(n,r+1); double *key=malloc(sizeof(double)*L); uint32_t *ms=malloc(sizeof(uint32_t)*L); long j=0;
    for(uint32_t m=0;m<N;m++){ int pc=__builtin_popcount(m); if(pc<=r) {F[m]=1; cnt++;} else if(pc==r+1){ double s=0; for(int i=0;i<n;i++) if(m>>i&1) s+=e[i]; key[j]=s; ms[j]=m; j++; } }
    // select rem smallest keys (partial selection by sorting indices)
    { KP*pp=malloc(sizeof(KP)*L); for(long a=0;a<L;a++){pp[a].k=key[a];pp[a].m=ms[a];}
      qsort(pp,L,sizeof(KP),cmpP); for(long a=0;a<rem;a++){ F[pp[a].m]=1; cnt++; } free(pp);
    }
    free(key); free(ms);
    printf("# equal-weights: full levels 0..%d, %ld sets of level %d\n",r,rem,r+1);
  } else dfs(0,0.0,1,0);
  // Bjorner/Morse critical count: S in D, coord 0 (smallest weight) not in S, S+{0} not in D
  long betti=0; for(uint32_t m=0;m<N;m++) if(F[m] && !(m&1) && !F[m|1]) betti++;
  // A_k profile and symmetric top
  long A[33]={0}; for(uint32_t m=0;m<N;m++) if(F[m]) A[__builtin_popcount(m)]++;
  long symtop=0; for(int k=0;k<=n;k++) symtop += (k&1?-1:1)*A[k];
  // FWHT
  for(size_t h=1;h<N;h<<=1) for(size_t i=0;i<N;i+=h<<1) for(size_t j=i;j<i+h;j++){ int32_t a=F[j],b=F[j+h]; F[j]=a+b; F[j+h]=a-b; }
  long Mtop=F[N-1]; long N0=F[0];
  double Wk[33]={0}; long maxabs=0; int maxlev=-1; long ge=0, gt=0; long am=labs(Mtop);
  for(size_t m=0;m<N;m++){ long v=F[m]; int pc=__builtin_popcount((unsigned)m); Wk[pc]+=(double)v*v; if(m){ if(labs(v)>maxabs){maxabs=labs(v);maxlev=pc;} } if(labs(v)>=am) ge++; if(labs(v)>am) gt++; }
  double alpha=(double)N0/N, la=log(1.0/alpha);
  // level-k inequality ratios
  double cum=0, worst=0; int kworst=0; double ratio_n=NAN;
  for(int k=0;k<=n;k++){ cum+=Wk[k]; if(k>=1 && k<=2*la){ double b=pow(2*M_E/k*la,k)*(double)N0*N0; double r=cum/b; if(r>worst){worst=r;kworst=k;} if(k==n) ratio_n=r; } }
  double lk_top = (n<=2*la)? (double)N0*pow(2*M_E/n*la, n/2.0) : NAN;  // |M| <= N((2e/n)ln(1/a))^{n/2}
  double hyp_top = sqrt((double)N*N0);                     // min over rho at rho=1
  double peres_top = 0.5*sqrt(3.0)*pow(2.0,n)*pow((double)n,-0.25); // W^{>=n}[1-2f] <= 3 NS_{1/n} <= 3/sqrt(n) (Peres, O'Donnell's proof)
  double W1=Wk[1], lev1=2.0*N0*(double)N0*la; // in F^2 units: sum_{|T|=1}F^2 vs 2 N^2 ln(1/a)
  printf("sys=%d x=%ld n=%d Q=%ld N=%ld M=%ld symtop=%ld rms=%.3f |M|/rms=%.3f frac_ge=%.4f frac_gt=%.4f maxabs(T!=0)=%ld@lev%d betti=%ld\n",
    sysid,X,n,Q,N0,Mtop,symtop,sqrt((double)N0),am/sqrt((double)N0),(double)ge/N,(double)gt/N,maxabs,maxlev,betti);
  printf("  bounds on |M|: trivial=%ld bjorner=%ld levelk(k=n)=%.3g hypercontr=%.3g peres=%.3g | level-k worst ratio=%.3f@k=%d ratio(k=n)=%.3g | level1 W1/(2N^2 ln(1/a))=%.3f\n",
    N0,betti,lk_top,hyp_top,peres_top,worst,kworst,ratio_n,W1/lev1);
  printf("  A_k:"); for(int k=0;k<=n;k++) if(A[k]) printf(" %ld",A[k]); printf("\n");
  printf("  W_k/(2^n N) vs binom/2^n:"); for(int k=0;k<=n;k++) printf(" %d:%.4f/%.4f",k,Wk[k]/((double)N*N0),(double)binom(n,k)/N); printf("\n");
  free(F); return 0;
}
