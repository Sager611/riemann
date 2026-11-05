/* nblp numerics: ||1 + sum_{k<=N} mu(k)(1-k/N) f_k||_{L^p(0,1)}, f_k = {1/(kx)} - {1/x}/k.
   On I_m=(1/(m+1),1/m], f_k = {m/k}; value v(m) = 1 + sum_k a_k {m/k}, a_k = mu(k)(1-k/N).
   Exact integer sieve: N*D(m) = sum_{j<=m} sum_{k|j,k<=N} mu(k)(N-k); v(m) = (N + m*NA - N*D(m))/N,
   NA = sum_k mu(k)(N-k)/k.  norm^p = sum_m |v(m)|^p /(m(m+1)) + tail, tail <= Vmax^p/(M+1). */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
static int mu[100001];
static void mobius(int n){ int *pr=calloc(n+1,sizeof(int)); char *c=calloc(n+1,1); int np=0; mu[1]=1;
  for(int i=2;i<=n;i++){ if(!c[i]){pr[np++]=i; mu[i]=-1;} for(int j=0;j<np && (long)i*pr[j]<=n;j++){ c[i*pr[j]]=1; if(i%pr[j]==0){mu[i*pr[j]]=0;break;} mu[i*pr[j]]=-mu[i]; } } free(pr); free(c);}
int main(int argc,char**argv){
  long long M = atoll(argv[1]); int nN=argc-2; mobius(100000);
  const long long B = 1LL<<23; long long *d = malloc(B*sizeof(long long));
  for(int a=0;a<nN;a++){ int N = atoi(argv[2+a]);
    long double NA=0; long long sumabs=0; for(int k=1;k<=N;k++){ NA += (long double)mu[k]*(N-k)/k; }
    long double Vmax=1; for(int k=2;k<=N;k++) Vmax += (long double)abs(mu[k])*(N-k)/N*(k-1)/k;
    long double S11=0,S2=0; long long ND=0; long long next_ck=1000000;
    for(long long st=1; st<=M; st+=B){ long long en = st+B-1; if(en>M) en=M; long long len=en-st+1;
      memset(d,0,len*sizeof(long long));
      for(int k=1;k<=N;k++){ if(!mu[k]||k==N) continue; long long w=(long long)mu[k]*(N-k);
        long long j0 = ((st+k-1)/k)*k; for(long long j=j0;j<=en;j+=k) d[j-st]+=w; }
      for(long long i=0;i<len;i++){ long long m=st+i; ND += d[i];
        long double v = ((long double)N + (long double)m*NA - (long double)ND)/N;
        long double wm = 1.0L/((long double)m*(m+1)); long double av=fabsl(v);
        S11 += powl(av,1.1L)*wm; S2 += av*av*wm;
        if(m==next_ck){ printf("N=%d M=%lld S1.1=%.12Lf norm1.1=%.10Lf S2=%.12Lf norm2=%.10Lf\n",N,m,S11,powl(S11,1/1.1L),S2,sqrtl(S2)); fflush(stdout); next_ck*=10; }
      } }
    long double tb11 = powl(Vmax,1.1L)/(M+1), tb2 = Vmax*Vmax/(M+1);
    printf("FINAL N=%d M=%lld Vmax=%.4Lf norm1.1 in [%.10Lf, %.10Lf]  norm2 in [%.10Lf, %.10Lf]\n",N,M,Vmax,powl(S11,1/1.1L),powl(S11+tb11,1/1.1L),sqrtl(S2),sqrtl(S2+tb2)); fflush(stdout);
  }
  free(d); return 0; }
