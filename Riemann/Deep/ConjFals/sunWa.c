// W_a(x) = sum_{n<=x} (-a)^{Omega(n)} for several real a >= 2; per-octave sup of |W_a(n)|/n^{log2 a}
// (2408.04143 Conjecture newcon2: |W_a(x)| < x^{log2 a} for x >= x0(a); s_a = limsup < 1?)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#define NA 64
int main(int argc,char**argv){
  uint64_t N=strtoull(argv[1],0,10); int na=argc-2; double a[NA], al[NA], pw[NA][64]; long double W[NA];
  for(int j=0;j<na;j++){ a[j]=atof(argv[j+2]); al[j]=log(a[j])/log(2.0); W[j]=0; pw[j][0]=1; for(int c=1;c<64;c++) pw[j][c]=pw[j][c-1]*(-a[j]); }
  uint64_t S=1<<22; uint64_t R=(uint64_t)sqrtl((long double)N)+2;
  char*isc=calloc(R+1,1); uint32_t*pr=malloc(sizeof(uint32_t)*(R/2+10)); uint64_t np=0;
  for(uint64_t i=2;i<=R;i++){ if(!isc[i]){pr[np++]=i; for(uint64_t j=i*i;j<=R;j+=i) isc[j]=1;} }
  uint64_t*prod=malloc(8*S); uint8_t*cnt=malloc(S);
  int oct=0; double mx[NA]; uint64_t am[NA]; for(int j=0;j<na;j++){mx[j]=0;am[j]=0;}
  for(uint64_t L=1;L<=N;L+=S){
    uint64_t H=L+S-1; if(H>N) H=N; uint64_t len=H-L+1;
    for(uint64_t i=0;i<len;i++){prod[i]=1;cnt[i]=0;}
    for(uint64_t k=0;k<np;k++){ uint64_t p=pr[k]; if(p*p>H) break;
      for(uint64_t q=p; q<=H; q*=p){ uint64_t st=((L+q-1)/q)*q;
        for(uint64_t m=st;m<=H;m+=q){ prod[m-L]*=p; cnt[m-L]++; }
        if(q>H/p) break; } }
    for(uint64_t i=0;i<len;i++){ uint64_t n=L+i; int c=cnt[i]; if(prod[i]!=n) c++;
      int o=63-__builtin_clzll(n);
      if(o!=oct){ if(oct>=8){ printf("oct %2d",oct); for(int j=0;j<na;j++) printf(" %.5f",mx[j]); printf("\n"); fflush(stdout);} oct=o; for(int j=0;j<na;j++){mx[j]=0;} }
      double ln=log((double)n);
      for(int j=0;j<na;j++){ W[j]+=pw[j][c]; double r=fabs((double)W[j])*exp(-al[j]*ln); if(r>mx[j]){mx[j]=r;am[j]=n;} }
    }
  }
  printf("a   "); for(int j=0;j<na;j++) printf(" %.5g",a[j]); printf("\n");
  return 0;
}
