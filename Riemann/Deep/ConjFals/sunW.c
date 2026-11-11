// W(x) = sum_{n<=x} (-2)^{Omega(n)}; segmented sieve; reports octave extremes of W(n)/n and any n>=3078 with |W(n)|>=n
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
int main(int argc,char**argv){
  uint64_t N=strtoull(argv[1],0,10); uint64_t S=1<<23;
  uint64_t R=(uint64_t)sqrtl((long double)N)+2;
  char*isc=calloc(R+1,1); uint32_t*pr=malloc(sizeof(uint32_t)*(R/2+10)); uint64_t np=0;
  for(uint64_t i=2;i<=R;i++){ if(!isc[i]){pr[np++]=i; for(uint64_t j=i*i;j<=R;j+=i) isc[j]=1;} }
  uint64_t*prod=malloc(8*S); uint8_t*cnt=malloc(S);
  int64_t W=0; int oct=0; double mx=-1e9,mn=1e9; uint64_t amx=0,amn=0; long nviol=0;
  for(uint64_t L=1;L<=N;L+=S){
    uint64_t H=L+S-1; if(H>N) H=N; uint64_t len=H-L+1;
    for(uint64_t i=0;i<len;i++){prod[i]=1;cnt[i]=0;}
    for(uint64_t k=0;k<np;k++){ uint64_t p=pr[k]; if(p*p>H) break;
      for(uint64_t q=p; q<=H; q*=p){ uint64_t st=((L+q-1)/q)*q;
        for(uint64_t m=st;m<=H;m+=q){ prod[m-L]*=p; cnt[m-L]++; }
        if(q>H/p) break; }
    }
    for(uint64_t i=0;i<len;i++){ uint64_t n=L+i; int c=cnt[i]; if(prod[i]!=n) c++;
      int64_t t=((int64_t)1)<<c; W+= (c&1)? -t : t;
      int o=63-__builtin_clzll(n);
      if(o!=oct){ if(oct>=11) printf("oct %d [2^%d,2^%d): max W/n=%.9f at %llu, min W/n=%.9f at %llu\n",oct,oct,oct+1,mx,(unsigned long long)amx,mn,(unsigned long long)amn); oct=o; mx=-1e9; mn=1e9; fflush(stdout);}
      double r=(double)W/(double)n; if(r>mx){mx=r;amx=n;} if(r<mn){mn=r;amn=n;}
      if(n>=3078 && (W>=(int64_t)n || -W>=(int64_t)n)){ nviol++; if(nviol<=20) printf("VIOL n=%llu W=%lld\n",(unsigned long long)n,(long long)W);}
      if(n<=3100 && n>=3070) printf("n=%llu W=%lld\n",(unsigned long long)n,(long long)W);
    }
  }
  printf("oct %d partial: max W/n=%.9f at %llu, min W/n=%.9f at %llu\n",oct,mx,(unsigned long long)amx,mn,(unsigned long long)amn);
  printf("N=%llu W(N)=%lld violations(n>=3078)=%ld\n",(unsigned long long)N,(long long)W,nviol);
  return 0;
}
