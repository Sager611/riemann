/* nblp2 numerics: ||1 + sum_{k<=N} mu(k)chi(k)(1-k/N) lambda_k||_{L^p(0,1)}, p = 1.1 and 2.
   lambda_k = kappa(1/(kx)) - delta k^{-1} kappa(1/x), kappa(y) = delta (phi(q)/q) y - S_chi(y), delta = [chi principal].
   On I_n = (1/(n+1),1/n]: lambda_k = delta S_chi(n)/k - S_chi(floor(n/k)), so the approximant equals
   v(n) = 1 + delta S_chi(n) P_N(1) - sum_{j<=n} chi(j) A_N(j),  A_N(j) = sum_{k|j,k<=N} mu(k)(1-k/N),
   P_N(1) = sum_{k<=N} mu(k)chi(k)(1-k/N)/k.  Exact integer sieve for N*A_N; norm^p = sum_n |v(n)|^p/(n(n+1)) + tail,
   0 <= tail <= Vmax^p/(M+1), Vmax = a priori bound on |v|.
   chi types: 1 trivial (zeta), 3 chi_{-3}, 4 chi_{-4}, 30 principal mod 3, 40 principal mod 4.
   usage: nbchi M type N1 [N2 ...]   (single-threaded) */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
static int mu[100001];
static void mobius(int n){ int *pr=calloc(n+1,sizeof(int)); char *c=calloc(n+1,1); int np=0; mu[1]=1;
  for(int i=2;i<=n;i++){ if(!c[i]){pr[np++]=i; mu[i]=-1;} for(int j=0;j<np && (long)i*pr[j]<=n;j++){ c[i*pr[j]]=1; if(i%pr[j]==0){mu[i*pr[j]]=0;break;} mu[i*pr[j]]=-mu[i]; } } free(pr); free(c);}
static int chi(int type,long long n){ int r;
  switch(type){ case 1: return 1;
    case 3: r=(int)(n%3); return r==0?0:(r==1?1:-1);
    case 4: r=(int)(n%4); return r==1?1:(r==3?-1:0);
    case 30: return (n%3)?1:0;
    case 40: return (n%2)?1:0; }
  return 0; }
int main(int argc,char**argv){
  long long M=atoll(argv[1]); int type=atoi(argv[2]); int nN=argc-3; mobius(100000);
  int principal=(type==1||type==30||type==40);
  /* sup|lambda_k| on the I_n: nonprincipal mod 3,4: |S_chi|<=1; trivial: {n/k}<1; principal q=3: c+2*max|kappa(n)| = 2/3+4/3; q=4: 1/2+1 */
  long double lb = (type==30)?2.0L:((type==40)?1.5L:1.0L);
  const long long B=1LL<<22; long long *d=malloc(B*sizeof(long long));
  for(int a=0;a<nN;a++){ int N=atoi(argv[3+a]);
    long double NA=0; for(int k=1;k<=N;k++) NA+=(long double)mu[k]*chi(type,k)*(N-k)/k;   /* N*P_N(1) */
    long double Vmax=1; for(int k=1;k<=N;k++){ if(principal&&k==1) continue; Vmax+=(long double)abs(mu[k]*chi(type,k))*(N-k)/N*lb; }
    long double S11=0,S2=0,vobs=0; long long ND=0,S0=0,next_ck=1000000;
    for(long long st=1;st<=M;st+=B){ long long en=st+B-1; if(en>M) en=M; long long len=en-st+1;
      memset(d,0,len*sizeof(long long));
      for(int k=1;k<=N;k++){ if(!mu[k]||k==N) continue; long long w=(long long)mu[k]*(N-k);
        long long j0=((st+k-1)/k)*k; for(long long j=j0;j<=en;j+=k) d[j-st]+=w; }
      for(long long i=0;i<len;i++){ long long m=st+i; int cm=chi(type,m); ND+=cm*d[i]; if(principal) S0+=(cm!=0);
        long double v = principal ? ((long double)N+(long double)S0*NA-(long double)ND)/N : ((long double)N-(long double)ND)/N;
        long double wm=1.0L/((long double)m*(m+1)); double av=fabs((double)v); if(av>vobs) vobs=av;
        S11+=(long double)pow(av,1.1)*wm; S2+=(long double)av*av*wm;
        if(m==next_ck){ printf("type=%d N=%d M=%lld norm1.1=%.10Lf norm2=%.10Lf max|v|=%.4Lf\n",type,N,m,powl(S11,1/1.1L),sqrtl(S2),vobs); fflush(stdout); next_ck*=10; }
      } }
    long double tb11=powl(Vmax,1.1L)/(M+1), tb2=Vmax*Vmax/(M+1);
    printf("FINAL type=%d N=%d M=%lld Vmax=%.4Lf norm1.1 in [%.10Lf, %.10Lf] norm2 in [%.10Lf, %.10Lf]\n",type,N,M,Vmax,powl(S11,1/1.1L),powl(S11+tb11,1/1.1L),sqrtl(S2),sqrtl(S2+tb2)); fflush(stdout);
  }
  free(d); return 0; }
