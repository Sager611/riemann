
// sieve.c: exact M(n), m(n)=sum mu(k)/k, psi(n) for n<=N; segmented. Sample outputs + global bound checks.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
int main(int argc,char**argv){
  uint64_t N=strtoull(argv[1],0,10); const char* sampf=argv[2]; // sample file: sorted integers
  FILE*fs=fopen(sampf,"r"); size_t ns=0,cap=1<<20; uint64_t*S=malloc(cap*8); uint64_t v;
  while(fscanf(fs,"%llu",(unsigned long long*)&v)==1){ if(ns==cap){cap*=2;S=realloc(S,cap*8);} S[ns++]=v;} fclose(fs);
  uint64_t r=(uint64_t)sqrt((double)N)+2; char*isc=calloc(r+1,1); uint32_t*P=malloc(sizeof(uint32_t)*(r+1)); int np=0;
  for(uint64_t i=2;i<=r;i++) if(!isc[i]){P[np++]=i; for(uint64_t j=i*i;j<=r;j+=i) isc[j]=1;}
  const uint64_t SEG=1<<20; int8_t*mu=malloc(SEG); uint64_t*pr=malloc(SEG*8); double*lam=malloc(SEG*8);
  long long M=0; long double m=0, psi=0; size_t si=0;
  double Mmax=0,MmaxX=0, Mrat=0,MratX=0, mmax=0,mmaxX=0, mrat=0,mratX=0, pmax=-1e9,pmaxX=0,pmin=1e9,pminX=0, prat=0,pratX=0;
  const double eps1=3.0/(M_PI*1e10), eps2=M_PI/(3e12);
  for(uint64_t L=1;L<=N;L+=SEG){
    uint64_t H=L+SEG-1; if(H>N)H=N; uint64_t len=H-L+1;
    for(uint64_t i=0;i<len;i++){mu[i]=1;pr[i]=1;lam[i]=0;}
    for(int k=0;k<np;k++){ uint64_t p=P[k]; if(p*p>H) break;
      uint64_t s=((L+p-1)/p)*p; for(uint64_t j=s;j<=H;j+=p){mu[j-L]=-mu[j-L]; pr[j-L]*=p;}
      uint64_t p2=p*p; s=((L+p2-1)/p2)*p2; for(uint64_t j=s;j<=H;j+=p2) mu[j-L]=0;
      // prime powers p^k (k>=1) in segment for Lambda
      for(uint64_t q=p;q<=H;q*=p){ if(q>=L) lam[q-L]=log((double)p); if(q>H/p) break; }
    }
    for(uint64_t i=0;i<len;i++){ uint64_t n=L+i;
      if(n==1){mu[i]=1;}
      else if(pr[i]<n){ if(mu[i]) mu[i]=-mu[i]; if(pr[i]==1) lam[i]=log((double)n); }
      M+=mu[i]; m+=(long double)mu[i]/(long double)n; psi+=lam[i];
      double sx=sqrt((double)n), sx1=sqrt((double)n+1.0);
      // M on [n,n+1): value M; ratio to bound at x=n (bound increasing)
      if(n>=33){ double a=fabs((double)M)/sx; if(a>Mmax){Mmax=a;MmaxX=n;} }
      { double b=eps1*n+11.39*sx; double a=fabs((double)M)/b; if(a>Mrat){Mrat=a;MratX=n;} }
      // m on [n,n+1): bound decreasing, worst x->n+1
      if(n>=3){ double a=fabsl(m)*sx1; if(a>mmax){mmax=a;mmaxX=n;} }
      { double b=eps1+11.39/sx1; double a=fabsl(m)/b; if(a>mrat){mrat=a;mratX=n;} }
      // psi on [n,n+1): psi(x)-x in (psi-n-1, psi-n]
      { double d1=(double)(psi-(long double)n), d0=(double)(psi-(long double)n-1.0L);
        if(n>=2){ double a=d1/sx; if(a>pmax){pmax=a;pmaxX=n;} double b=d0/sx1; if(b<pmin){pmin=b;pminX=n;} }
        double r1=fabs(d1)/(eps2*n+113.67*sx), r2=fabs(d0)/(eps2*(n+1.0)+113.67*sx1); if(r1>prat){prat=r1;pratX=n;} if(r2>prat){prat=r2;pratX=n;} }
      while(si<ns && S[si]==n){ printf("S %llu %lld %.21Le %.12Lf %d\n",(unsigned long long)n,M,m,psi,(int)mu[i]); si++; }
    }
  }
  printf("RANGE N=%llu\n",(unsigned long long)N);
  printf("M: max|M(x)|/sqrt(x) for 33<=x<=N: %.6f at %.0f\n",Mmax,MmaxX);
  printf("M: max |M(x)|/(3x/(pi 1e10)+11.39 sqrt x): %.6f at %.0f\n",Mrat,MratX);
  printf("m: max|m(x)|sqrt(x) 3<=x<=N: %.6f at %.0f\n",mmax,mmaxX);
  printf("m: max |m(x)|/(3/(pi 1e10)+11.39/sqrt x): %.6f at %.0f\n",mrat,mratX);
  printf("psi: max (psi(x)-x)/sqrt x: %.6f at %.0f ; min: %.6f at %.0f^-\n",pmax,pmaxX,pmin,pminX+1);
  printf("psi: max |psi(x)-x|/(pi x/3e12+113.67 sqrt x): %.6f at %.0f\n",prat,pratX);
  printf("final M=%lld m=%.15Le psi=%.6Lf\n",M,m,psi);
  return 0;}
