/* nblpaudit: ||F_N||_p, F_N = 1 + sum_{k=2}^N mu(k)(1-k/N) f_k, on I_m value
   v(m) = 1 + sum_k mu(k)(1-k/N) (m mod k)/k.  Integer sieve: N v(m) = N + m*T - A(m),
   T = sum_k mu(k)(N-k)/k, A(m) = sum_{j<=m} sum_{k|j,2<=k<=N} mu(k)(N-k).
   Norm^p = sum_m |v(m)|^p/(m(m+1)); tail rigorous <= Vmax^p/(M+1), heuristic = mean_{(M/2,M]}|v|^p/(M+1).
   usage: lpnorm N M ; prints p=1.1,1.3,2 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
int main(int argc,char**argv){
  long N=atol(argv[1]), M=atol(argv[2]); const long B=1000000;
  int *mu=calloc(N+1,sizeof(int)); char *comp=calloc(N+1,1);
  for(long i=1;i<=N;i++) mu[i]=1;
  for(long q=2;q<=N;q++) if(!comp[q]){ for(long j=q;j<=N;j+=q){ if(j>q) comp[j]=1; mu[j]=-mu[j]; }
      for(long j=q*q;j<=N;j+=q*q) mu[j]=0; }
  double T=0, Vmax=1; for(long k=2;k<=N;k++){ T+=(double)mu[k]*(N-k)/k; Vmax+=fabs((double)mu[k]*(N-k)/N)*(k-1)/k; }
  long long *a=malloc(B*sizeof(long long)); long long A=0;
  const double ps[3]={1.1,1.3,2.0}; double S[3]={0,0,0}, tailm[3]={0,0,0}; long cnt=0; double vmaxobs=0;
  long nchk=0; double maxdev=0;
  for(long L=1;L<=M;L+=B){ long R=L+B-1; if(R>M) R=M; long len=R-L+1;
    memset(a,0,len*sizeof(long long));
    for(long k=2;k<=N;k++){ if(!mu[k]) continue; long long w=(long long)mu[k]*(N-k);
      long j0=((L+k-1)/k)*k; for(long j=j0;j<=R;j+=k) a[j-L]+=w; }
    for(long i=0;i<len;i++){ long m=L+i; A+=a[i];
      double v=(N + (double)m*T - (double)A)/(double)N; double av=fabs(v); if(av>vmaxobs) vmaxobs=av;
      double w=1.0/((double)m*(double)(m+1));
      for(int t=0;t<3;t++){ double vp=pow(av,ps[t]); S[t]+=vp*w; if(m>M/2) tailm[t]+=vp; }
      if(m>M/2) cnt++;
      if(m%99991==1 && nchk<200){ /* direct check */ double d=1; for(long k=2;k<=N;k++) if(mu[k]) d+=(double)mu[k]*(N-k)/N*(double)(m%k)/k;
        double dev=fabs(d-v); if(dev>maxdev) maxdev=dev; nchk++; }
    }
  }
  for(int t=0;t<3;t++){ double p=ps[t]; double lo=S[t], hi=S[t]+pow(Vmax,p)/(M+1.0), he=S[t]+tailm[t]/cnt/(M+1.0);
    printf("N=%ld M=%ld p=%.1f norm_lo=%.10f norm_hi=%.10f norm_heur=%.10f\n",N,M,p,pow(lo,1/p),pow(hi,1/p),pow(he,1/p)); }
  printf("N=%ld Vmax_bound=%.4f Vmax_obs=%.4f directcheck n=%ld maxdev=%.2e\n",N,Vmax,vmaxobs,nchk,maxdev);
  return 0; }
