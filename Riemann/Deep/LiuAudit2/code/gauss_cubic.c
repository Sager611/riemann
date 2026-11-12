/* liuaudit2: normalized cubic Gauss sums gamma_2(P) = l^{-1/2} sum_{x in F_l^*} zeta3^{m(x)} e(-x w / l)
   for degree-1 primes P of Q(zeta12) above l = 1 mod 12, m(x) defined by x^((l-1)/3) = (r^4)^{m(x)} mod l (r = image of zeta12).
   Input lines: l r w (sorted by l). Output lines: re im. Single-threaded. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
static uint64_t pw(uint64_t b, uint64_t e, uint64_t m){ uint64_t r=1; b%=m; while(e){ if(e&1) r=r*b%m; b=b*b%m; e>>=1;} return r; }
static uint64_t primroot(uint64_t l){
  uint64_t q=l-1, f[64]; int nf=0;
  for(uint64_t p=2;p*p<=q;p++){ if(q%p==0){ f[nf++]=p; while(q%p==0) q/=p; } }
  if(q>1) f[nf++]=q;
  for(uint64_t g=2; g<l; g++){ int ok=1; for(int i=0;i<nf;i++) if(pw(g,(l-1)/f[i],l)==1){ok=0;break;} if(ok) return g; }
  return 0;
}
int main(int argc, char**argv){
  FILE*in=fopen(argv[1],"r"); FILE*out=fopen(argv[2],"w");
  if(!in||!out){ fprintf(stderr,"io\n"); return 1; }
  uint64_t l,r,w, curl=0, g=0; double *Er=NULL,*Ei=NULL; size_t cap=0;
  const double c3=-0.5, s3=0.86602540378443864676;
  while(fscanf(in,"%llu %llu %llu",(unsigned long long*)&l,(unsigned long long*)&r,(unsigned long long*)&w)==3){
    if(l!=curl){
      curl=l; g=primroot(l);
      if(l>cap){ free(Er); free(Ei); cap=l; Er=malloc(sizeof(double)*cap); Ei=malloc(sizeof(double)*cap); }
      for(uint64_t j=0;j<l;j++){ double a=-2.0*M_PI*(double)j/(double)l; Er[j]=cos(a); Ei[j]=sin(a); }
    }
    uint64_t r4=pw(r,4,l), r8=r4*r4%l, h=pw(g,(l-1)/3,l); int s = (h==r4)?1:((h==r8)?2:0);
    if(!s){ fprintf(stderr,"bad cube root l=%llu\n",(unsigned long long)l); return 2; }
    double A[3][2]={{0,0},{0,0},{0,0}}; uint64_t y=1; 
    for(uint64_t k=0;k<l-1;k++){ uint64_t idx=y*w%l; int m=(int)((s*k)%3); A[m][0]+=Er[idx]; A[m][1]+=Ei[idx]; y=y*g%l; }
    double re=A[0][0]+c3*A[1][0]-s3*A[1][1]+c3*A[2][0]+s3*A[2][1];
    double im=A[0][1]+c3*A[1][1]+s3*A[1][0]+c3*A[2][1]-s3*A[2][0];
    double sc=1.0/sqrt((double)l);
    fprintf(out,"%.17g %.17g\n",re*sc,im*sc);
  }
  fclose(in); fclose(out); return 0;
}
