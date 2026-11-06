// Independent Schur-test bound for || B || , B = positive prime operator on L^2(0,a), a = 23029/10000.
// (B v)(x) = sum_m w_m [ 1{x+d<a} v(x+d) + 1{x>d} v(x-d) + 1{0<d-x<a} v(d-x) ],  d = log m, m prime power < e^{2a}.
// h_k = B^k 1 is piecewise constant on the partition by BP_k, BP_k = {0,a} U (U_phi phi^{-1}(BP_{k-1}) n (0,a)).
// Breakpoints carry exact labels (e,P,Q): x = e*a + log(P/Q), e in {-1,0,1}, gcd(P,Q)=1.
// Final level: Schur test ||B|| <= ess sup (B h)/h with h = h_{K-1} (any positive test function).
// Usage: rhow K   (K = final level, e.g. 6). Single-threaded.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "weights.h"

typedef struct { double x; uint64_t P, Q; int e; } Lab;
static double A;            // a as double
static double DLOG[64];     // log m
static const double TOL = 1e-13;
static long near_ties = 0, flags = 0;

static uint64_t gcd64(uint64_t a, uint64_t b){ while(b){uint64_t t=a%b;a=b;b=t;} return a; }
static Lab mk(int e, uint64_t P, uint64_t Q){ uint64_t g=gcd64(P,Q); Lab L; L.e=e; L.P=P/g; L.Q=Q/g; L.x=(double)e*A+log((double)L.P/(double)L.Q); return L; }
static int same(const Lab*a,const Lab*b){ return a->e==b->e && a->P==b->P && a->Q==b->Q; }
static int is0(const Lab*L){ return L->e==0 && L->P==1 && L->Q==1; }
static int isA(const Lab*L){ return L->e==1 && L->P==1 && L->Q==1; }
static int cmp(const Lab*a,const Lab*b){
  if(same(a,b)) return 0;
  double dx=a->x-b->x;
  if(fabs(dx)>TOL) return dx<0?-1:1;
  near_ties++;
  if(a->e==b->e){ unsigned __int128 l=(unsigned __int128)a->P*b->Q, r=(unsigned __int128)b->P*a->Q; return l<r?-1:1; }
  flags++; fprintf(stderr,"UNRESOLVED near tie (%d,%llu,%llu) (%d,%llu,%llu) dx=%g\n",a->e,(unsigned long long)a->P,(unsigned long long)a->Q,b->e,(unsigned long long)b->P,(unsigned long long)b->Q,dx);
  return dx<0?-1:1;
}
// stream: map index m, type t: 0: beta-d (ascending), 1: beta+d (ascending), 2: d-beta (descending beta); 3: endpoints
typedef struct { int m, t; long i; Lab cur; int alive; } Stream;
static Lab *Bp; static long nB;
static int image(int m,int t,const Lab*b,Lab*out){
  uint64_t mm=MP[m];
  if(t==0) *out=mk(b->e,b->P,b->Q*mm);
  else if(t==1) *out=mk(b->e,b->P*mm,b->Q);
  else *out=mk(-b->e,b->Q*mm,b->P);
  if(is0(out)||isA(out)) return 0;
  if(out->x<=TOL||out->x>=A-TOL){
    if(out->x>-TOL && out->x<A+TOL){ flags++; fprintf(stderr,"candidate within TOL of endpoint\n"); }
    return 0;
  }
  return 1;
}
static void advance(Stream*s){
  if(s->t==3){ if(s->i==0){ s->cur=mk(0,1,1); s->i=1; return;} if(s->i==1){ s->cur=mk(1,1,1); s->i=2; return;} s->alive=0; return; }
  while(1){
    if(s->t<2){ if(s->i>=nB){s->alive=0;return;} Lab c; int ok=image(s->m,s->t,&Bp[s->i],&c); s->i++; if(ok){s->cur=c;return;} }
    else { if(s->i<0){s->alive=0;return;} Lab c; int ok=image(s->m,s->t,&Bp[s->i],&c); s->i--; if(ok){s->cur=c;return;} }
  }
}
static Stream *S; static int nS; static int *heap; static int hn;
static int hless(int a,int b){ return cmp(&S[a].cur,&S[b].cur)<0; }
static void hpush(int s){ int i=hn++; heap[i]=s; while(i>0){int p=(i-1)/2; if(hless(heap[i],heap[p])){int t=heap[i];heap[i]=heap[p];heap[p]=t;i=p;} else break;} }
static int hpop(void){ int top=heap[0]; heap[0]=heap[--hn]; int i=0; while(1){int l=2*i+1,r=l+1,s=i; if(l<hn&&hless(heap[l],heap[s]))s=l; if(r<hn&&hless(heap[r],heap[s]))s=r; if(s==i)break; int t=heap[i];heap[i]=heap[s];heap[s]=t;i=s;} return top; }

// lookup of previous-level cell containing y, with monotone pointer
static double *Hprev; static long *ptr;
static double hprev_at(int key,double y,int dir){
  long i=ptr[key];
  if(i<0){ long lo=0,hi=nB-1; while(hi-lo>1){long md=(lo+hi)/2; if(Bp[md].x<y) lo=md; else hi=md;} i=lo; }
  if(dir>0){ while(i+1<nB-1 && Bp[i+1].x<y) i++; } else { while(i>0 && Bp[i].x>y) i--; }
  if(!(Bp[i].x<y && y<Bp[i+1].x)){ while(i+1<nB-1 && Bp[i+1].x<y) i++; while(i>0 && Bp[i].x>y) i--; }
  if(!(y-Bp[i].x>TOL && Bp[i+1].x-y>TOL)){ flags++; if(flags<20) fprintf(stderr,"lookup within TOL: y=%.17g cell [%.17g,%.17g]\n",y,Bp[i].x,Bp[i+1].x); }
  ptr[key]=i; return Hprev[i];
}

int main(int argc,char**argv){
  int KF=argc>1?atoi(argv[1]):6;
  A=23029.0/10000.0;
  for(int m=0;m<NMP;m++) DLOG[m]=log((double)MP[m]);
  // level 0
  nB=2; Bp=malloc(2*sizeof(Lab)); Bp[0]=mk(0,1,1); Bp[1]=mk(1,1,1);
  Hprev=malloc(sizeof(double)); Hprev[0]=1.0;
  nS=3*NMP+1; S=malloc(nS*sizeof(Stream)); heap=malloc(nS*sizeof(int)); ptr=malloc((3*NMP+1)*sizeof(long));
  long total_pos=0;
  for(int k=1;k<=KF;k++){
    int final=(k==KF);
    for(int m=0;m<NMP;m++) for(int t=0;t<3;t++){ Stream*s=&S[3*m+t]; s->m=m; s->t=t; s->alive=1; s->i=(t<2)?0:nB-1; }
    S[nS-1].m=0; S[nS-1].t=3; S[nS-1].i=0; S[nS-1].alive=1;
    hn=0; for(int s=0;s<nS;s++){ advance(&S[s]); if(S[s].alive) hpush(s); }
    for(int q=0;q<3*NMP+1;q++) ptr[q]=-1;
    long cap=final?0:1024, n=0; Lab *Nb=final?NULL:malloc(cap*sizeof(Lab));
    double *Hn=NULL; long hcap=0;
    Lab last; int have=0; long cells=0; double minw=1e9, maxr=0, minr=1e300, minh=1e300, maxh=0; long maxcell=-1; double maxmid=0;
    int maxterms=0;
    while(hn>0){
      int s=hpop(); Lab c=S[s].cur; advance(&S[s]); if(S[s].alive) hpush(s);
      if(have && same(&c,&last)) continue;
      if(have){
        if(!(c.x>last.x)){ flags++; fprintf(stderr,"order violation\n"); }
        // cell (last, c)
        double w=c.x-last.x; if(w<minw) minw=w;
        double mid=0.5*(last.x+c.x);
        double val=0; int terms=0;
        for(int m=0;m<NMP;m++){
          double d=DLOG[m];
          double y1=mid+d; if(fabs(y1-A)<=TOL){flags++;}
          if(y1<A){ val+=WUP[m]*hprev_at(3*m,y1,+1); terms++; }
          double y2=mid-d; if(fabs(y2)<=TOL){flags++;}
          if(y2>0){ val+=WUP[m]*hprev_at(3*m+1,y2,+1); terms++; }
          double y3=d-mid; if(fabs(y3)<=TOL||fabs(y3-A)<=TOL){flags++;}
          if(y3>0 && y3<A){ val+=WUP[m]*hprev_at(3*m+2,y3,-1); terms++; }
        }
        if(terms>maxterms) maxterms=terms;
        if(final){
          double h=hprev_at(3*NMP,mid,+1);
          if(h<minh) minh=h; if(h>maxh) maxh=h;
          double r=val/h; if(r>maxr){maxr=r;maxcell=cells;maxmid=mid;} if(r<minr) minr=r;
        } else {
          if(cells>=hcap){ hcap=hcap?2*hcap:1024; Hn=realloc(Hn,hcap*sizeof(double)); }
          Hn[cells]=val; if(val<minh) minh=val; if(val>maxh) maxh=val;
        }
        cells++;
      }
      if(!final){ if(n>=cap){cap*=2; Nb=realloc(Nb,cap*sizeof(Lab));} Nb[n++]=c; }
      last=c; have=1;
    }
    total_pos+=cells+1;
    if(final){
      printf("level %d: cells %ld (breakpoints incl. ends %ld) min width %.6e max terms %d\n",k,cells,cells+1,minw,maxterms);
      printf("  test function h = h_%d: min %.17g max %.17g\n",k-1,minh,maxh);
      printf("  computed max ratio (B h)/h = %.17g at cell %ld mid %.17g; min ratio %.17g\n",maxr,maxcell,maxmid,minr);
      printf("  near ties resolved exactly %ld, flags %ld\n",near_ties,flags);
      printf("RESULT maxr_hex %a minh_hex %a maxterms %d flags %ld\n",maxr,minh,maxterms,flags);
    } else {
      printf("level %d: cells %ld min width %.6e  h_k min %.17g max %.17g  maxterms %d  near ties %ld flags %ld\n",k,cells,minw,minh,maxh,maxterms,near_ties,flags);
      fflush(stdout);
      free(Bp); free(Hprev); Bp=Nb; nB=n; Hprev=Hn;
    }
  }
  return flags?3:0;
}

