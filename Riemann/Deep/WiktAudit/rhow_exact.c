// rhow_exact.c: exhaustive Schur-test bound for ||B||, B = positive prime operator on L^2(0,a), a = 23029/10000,
// with ALL order decisions exact: breakpoint labels (e,P,Q) mean e*a + log(P/Q), gcd(P,Q)=1, e in {-1,0,1};
// equal e: compare P1*Q2 vs P2*Q1 (exact 128-bit); e2-e1 = j != 0: compare R = P1Q2/(Q1P2) with e^{j a} through
// the integer enclosure EA/E2A of 2^56 e^{a}, 2^56 e^{2a} (exact 192-bit products). No libm call decides anything;
// doubles are used only for the values h_k (positive sums) and for diagnostics.
// (B v)(x) = sum_m w_m [1{x+d<a} v(x+d) + 1{x>d} v(x-d) + 1{d-a<x<d} v(d-x)], d = log m.
// h_k = B^k 1 is constant on the cells of BP_k = {0,a} U (U_phi phi^{-1}(BP_{k-1}) n (0,a)).  Usage: rhow_exact K
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "consts.h"
typedef unsigned __int128 u128;
typedef struct { uint64_t P, Q; int e; } Lab;
static long ambiguous = 0, flags = 0;
static uint64_t gcd64(uint64_t a, uint64_t b){ while(b){uint64_t t=a%b;a=b;b=t;} return a; }
static Lab mk(int e, uint64_t P, uint64_t Q){ uint64_t g=gcd64(P,Q); Lab L={P/g,Q/g,e}; return L; }
static int same(const Lab*a,const Lab*b){ return a->e==b->e && a->P==b->P && a->Q==b->Q; }
// 192-bit compare of X*2^S (X < 2^128) with N*Y (N < 2^64, Y < 2^128): returns -1,0,1
static int cmp_shift_vs_prod(u128 X, int S, uint64_t N, u128 Y){
  // left = X << S as (hi64, lo128)
  u128 lo_l = X << S; uint64_t hi_l = (uint64_t)(X >> (128 - S));
  // right = N*Y = N*Yhi*2^64 + N*Ylo
  uint64_t Ylo = (uint64_t)Y, Yhi = (uint64_t)(Y >> 64);
  u128 p0 = (u128)N * Ylo, p1 = (u128)N * Yhi;      // right = p1*2^64 + p0
  u128 lo_r = p0 + (p1 << 64); uint64_t carry = (lo_r < p0) ? 1 : 0;
  uint64_t hi_r = (uint64_t)(p1 >> 64) + carry;
  if(hi_l != hi_r) return hi_l < hi_r ? -1 : 1;
  if(lo_l != lo_r) return lo_l < lo_r ? -1 : 1;
  return 0;
}
// sign of (X/Y - c) where 2^S c in (LO, HI): -1 if X/Y < c, +1 if X/Y > c, 0 if undecided
static int cmp_ratio_const(u128 X, u128 Y, uint64_t LO, uint64_t HI){
  if(cmp_shift_vs_prod(X, ESHIFT, LO, Y) <= 0) return -1;   // X/Y <= LO/2^S < c
  if(cmp_shift_vs_prod(X, ESHIFT, HI, Y) >= 0) return 1;    // X/Y >= HI/2^S > c
  return 0;
}
static int cmp(const Lab*A,const Lab*B){
  if(same(A,B)) return 0;
  u128 X=(u128)A->P*B->Q, Y=(u128)A->Q*B->P;     // log(P1/Q1) - log(P2/Q2) = log(X/Y)
  int j=B->e-A->e;                                // A < B  <=>  log(X/Y) < j a
  if(j==0) return X<Y?-1:1;
  int s;
  if(j==1) s=cmp_ratio_const(X,Y,EA_LO,EA_HI);
  else if(j==2) s=cmp_ratio_const(X,Y,E2A_LO,E2A_HI);
  else if(j==-1) s=-cmp_ratio_const(Y,X,EA_LO,EA_HI);     // X/Y vs e^{-a} <=> Y/X vs e^{a}, reversed
  else if(j==-2) s=-cmp_ratio_const(Y,X,E2A_LO,E2A_HI);
  else { flags++; fprintf(stderr,"bad j\n"); return 0; }
  if(s==0){ ambiguous++; flags++; fprintf(stderr,"ambiguous comparison\n"); return 0; }
  return s;   // s<0: X/Y < e^{ja}  => A < B
}
static Lab ZERO, AEND;
static Lab Bp_img(int m,int t,const Lab*b){ uint64_t mm=MP[m];
  if(t==0) return mk(b->e,b->P,b->Q*mm);      // b - d
  if(t==1) return mk(b->e,b->P*mm,b->Q);      // b + d
  return mk(-b->e,b->Q*mm,b->P);              // d - b
}
static int inside(const Lab*L){ return cmp(&ZERO,L)<0 && cmp(L,&AEND)<0; }
typedef struct { int m,t; long i; Lab cur; int alive; } Stream;
static Lab *Bp; static long nB; static double *Hprev;
static void advance(Stream*s){
  if(s->t==3){ if(s->i==0){s->cur=ZERO;s->i=1;return;} if(s->i==1){s->cur=AEND;s->i=2;return;} s->alive=0; return; }
  while(1){
    if(s->t<2){ if(s->i>=nB){s->alive=0;return;} Lab c=Bp_img(s->m,s->t,&Bp[s->i]); s->i++; if(inside(&c)){s->cur=c;return;} }
    else { if(s->i<0){s->alive=0;return;} Lab c=Bp_img(s->m,s->t,&Bp[s->i]); s->i--; if(inside(&c)){s->cur=c;return;} }
  }
}
static Stream *S; static int nS; static int *heap; static int hn;
static int hless(int a,int b){ return cmp(&S[a].cur,&S[b].cur)<0; }
static void hpush(int s){ int i=hn++; heap[i]=s; while(i>0){int p=(i-1)/2; if(hless(heap[i],heap[p])){int t=heap[i];heap[i]=heap[p];heap[p]=t;i=p;} else break;} }
static int hpop(void){ int top=heap[0]; heap[0]=heap[--hn]; int i=0; while(1){int l=2*i+1,r=l+1,s=i; if(l<hn&&hless(heap[l],heap[s]))s=l; if(r<hn&&hless(heap[r],heap[s]))s=r; if(s==i)break; int t=heap[i];heap[i]=heap[s];heap[s]=t;i=s;} return top; }
// index i with Bp[i] <= L < Bp[i+1] (cells 0..nB-2), monotone pointer
static long *ptr;
static long locate(int key,const Lab*L){
  long i=ptr[key];
  if(i<0){ long lo=0,hi=nB-1; while(hi-lo>1){long md=(lo+hi)/2; if(cmp(&Bp[md],L)<=0) lo=md; else hi=md;} i=lo; }
  while(i+1<nB-1 && cmp(&Bp[i+1],L)<=0) i++;
  while(i>0 && cmp(&Bp[i],L)>0) i--;
  if(!(cmp(&Bp[i],L)<=0 && cmp(L,&Bp[i+1])<0)){ flags++; fprintf(stderr,"locate failed\n"); }
  ptr[key]=i; return i;
}
int main(int argc,char**argv){
  int KF=argc>1?atoi(argv[1]):6;
  ZERO=mk(0,1,1); AEND=mk(1,1,1);
  nB=2; Bp=malloc(2*sizeof(Lab)); Bp[0]=ZERO; Bp[1]=AEND; Hprev=malloc(sizeof(double)); Hprev[0]=1.0;
  nS=3*NMP+1; S=malloc(nS*sizeof(Stream)); heap=malloc(nS*sizeof(int)); ptr=malloc((3*NMP+1)*sizeof(long));
  Lab LAM[64], LD[64], LDMA[64];
  for(int m=0;m<NMP;m++){ LAM[m]=mk(1,1,MP[m]); LD[m]=mk(0,MP[m],1); LDMA[m]=mk(-1,MP[m],1); }
  for(int k=1;k<=KF;k++){
    int final=(k==KF);
    for(int m=0;m<NMP;m++) for(int t=0;t<3;t++){ Stream*s=&S[3*m+t]; s->m=m; s->t=t; s->alive=1; s->i=(t<2)?0:nB-1; }
    S[nS-1].m=0; S[nS-1].t=3; S[nS-1].i=0; S[nS-1].alive=1;
    hn=0; for(int s=0;s<nS;s++){ advance(&S[s]); if(S[s].alive) hpush(s); }
    for(int q=0;q<3*NMP+1;q++) ptr[q]=-1;
    long cap=final?0:1024, n=0; Lab *Nb=final?NULL:malloc(cap*sizeof(Lab)); double *Hn=NULL; long hcap=0;
    Lab last; int have=0; long cells=0; double maxr=0, minh=1e300, maxh=0; int maxterms=0; long maxcell=-1;
    while(hn>0){
      int s=hpop(); Lab c=S[s].cur; advance(&S[s]); if(S[s].alive) hpush(s);
      if(have && same(&c,&last)) continue;
      if(have){
        if(cmp(&last,&c)>=0){ flags++; fprintf(stderr,"order violation\n"); }
        // cell (last, c): exact domain tests and exact image lookups
        double val=0; int terms=0;
        for(int m=0;m<NMP;m++){
          if(cmp(&c,&LAM[m])<=0){            // whole cell has x + d < a
            Lab Li=Bp_img(m,1,&last); long i=locate(3*m,&Li); val+=WUP[m]*Hprev[i]; terms++; }
          if(cmp(&last,&LD[m])>=0){           // whole cell has x > d
            Lab Li=Bp_img(m,0,&last); long i=locate(3*m+1,&Li); val+=WUP[m]*Hprev[i]; terms++; }
          if(cmp(&last,&LDMA[m])>=0 && cmp(&c,&LD[m])<=0){   // whole cell has d - a < x < d
            Lab Li=Bp_img(m,2,&c); long i=locate(3*m+2,&Li); val+=WUP[m]*Hprev[i]; terms++; }
        }
        if(terms>maxterms) maxterms=terms;
        if(final){
          long i=locate(3*NMP,&last); double h=Hprev[i];
          if(h<minh) minh=h; if(h>maxh) maxh=h;
          double r=val/h; if(r>maxr){maxr=r;maxcell=cells;}
        } else {
          if(cells>=hcap){ hcap=hcap?2*hcap:1024; Hn=realloc(Hn,hcap*sizeof(double)); }
          Hn[cells]=val; if(val<minh) minh=val; if(val>maxh) maxh=val;
        }
        cells++;
      }
      if(!final){ if(n>=cap){cap*=2; Nb=realloc(Nb,cap*sizeof(Lab));} Nb[n++]=c; }
      last=c; have=1;
    }
    if(final){
      printf("level %d: cells %ld, max terms %d; test function h = h_%d min %.17g max %.17g\n",k,cells,maxterms,k-1,minh,maxh);
      printf("  computed max ratio (B h)/h = %.17g at cell %ld; ambiguous %ld flags %ld\n",maxr,maxcell,ambiguous,flags);
      printf("RESULT maxr_hex %a minh_hex %a maxterms %d flags %ld\n",maxr,minh,maxterms,flags);
    } else {
      printf("level %d: cells %ld  h_k min %.17g max %.17g  maxterms %d  ambiguous %ld flags %ld\n",k,cells,minh,maxh,maxterms,ambiguous,flags);
      fflush(stdout); free(Bp); free(Hprev); Bp=Nb; nB=n; Hprev=Hn;
    }
  }
  return flags?3:0;
}
