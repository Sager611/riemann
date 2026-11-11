
// sieve2.c: psi and sum Lambda(n)/n extrema (lem:pernic checks) for n<=N.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
int main(int argc,char**argv){
  uint64_t N=strtoull(argv[1],0,10);
  uint64_t r=(uint64_t)sqrt((double)N)+2; char*isc=calloc(r+1,1); uint32_t*P=malloc(sizeof(uint32_t)*(r+1)); int np=0;
  for(uint64_t i=2;i<=r;i++) if(!isc[i]){P[np++]=i; for(uint64_t j=i*i;j<=r;j+=i) isc[j]=1;}
  const uint64_t SEG=1<<20; uint64_t*pr=malloc(SEG*8); double*lam=malloc(SEG*8);
  long double psi=0, S1=0; const long double EG=0.5772156649015328606065120900824024L;
  double pmax=-1e9,pmaxX=0,pmin4=1e9,pmin4X=0,fmax=-1e9,fmaxX=0,fmin=1e9,fminX=0;
  for(uint64_t L=1;L<=N;L+=SEG){
    uint64_t H=L+SEG-1; if(H>N)H=N; uint64_t len=H-L+1;
    for(uint64_t i=0;i<len;i++){pr[i]=1;lam[i]=0;}
    for(int k=0;k<np;k++){ uint64_t p=P[k]; if(p*p>H) break;
      uint64_t s=((L+p-1)/p)*p; for(uint64_t j=s;j<=H;j+=p) pr[j-L]*=p;
      for(uint64_t q=p;q<=H;q*=p){ if(q>=L) lam[q-L]=log((double)p); if(q>H/p) break; } }
    for(uint64_t i=0;i<len;i++){ uint64_t n=L+i;
      if(n>1 && pr[i]==1) lam[i]=logl((long double)n);
      psi+=lam[i]; S1+=(long double)lam[i]/(long double)n;
      if(n>=2){ double a=(double)((psi-(long double)n)/sqrtl((long double)n)); if(a>pmax){pmax=a;pmaxX=n;} }
      if(n>=10000){ double b=(double)((psi-(long double)(n+1))/sqrtl((long double)(n+1))); if(b<pmin4){pmin4=b;pmin4X=n+1;} }
      { double fa=(double)(S1-logl((long double)n)+EG), fb=(double)(S1-logl((long double)(n+1))+EG);
        if(fa>fmax){fmax=fa;fmaxX=n;} if(fb<fmin){fmin=fb;fminX=n+1;} }
    }
  }
  printf("N=%llu\n",(unsigned long long)N);
  printf("max (psi(x)-x)/sqrt x on [1,N]: %.10f at %.0f   (paper: 0.79059275 at 110102617)\n",pmax,pmaxX);
  printf("min (psi(x)-x)/sqrt x on [1e4,N]: %.10f at %.0f^-   (paper on [1e4,1e13]: -0.7509024438 at 36917099^-)\n",pmin4,pmin4X);
  printf("max sum_{n<=x}Lambda(n)/n-(log x-gamma) on [1,N]: %.10f at %.0f   (paper: 0.787.. at 110102617)\n",fmax,fmaxX);
  printf("min sum_{n<=x}Lambda(n)/n-(log x-gamma) on [1,N]: %.10f at %.0f^-   (paper: -0.7585825520 at 1423^-)\n",fmin,fminX);
  return 0;}
