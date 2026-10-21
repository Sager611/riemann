/* floor_arb.c — Arb (ball arithmetic) certificate of the Polymath improved-
   triangle (SHARP) or plain-triangle (T) Euler-mollified floor at ONE window N,
   at fixed exact t = tn/td, y^2 = yn/yd. Mirrors floor.c and Gomila's 'n' mode:
     sigma = (1+y)/2 + (t/4)log(x_N/4pi) - t/(2x_N^2)(1-3y+4y(1+y)/x_N^2)_+
     g = exp(0.02y - (y/2)log(x_N/4pi)),  x_N = 4 pi N^2 - pi t/4
     b_n = sum_{d|n,d|D,n/d<=N} lam_d bt(n/d),  a_n = same * (n/d)^y,  lam_d = prod_{p|d}(-bt(p))
     S = sum_{n>=2} W n^{-sigma},  W = |b|+g|a|  or  max(|b-ga|, c|b+ga|), c=(1-g)/(1+g)
     floor = (1-g-S)/modmoll - corr,  modmoll = prod(1+bt(p)p^{-sigma}),
     corr = g * sum_{n=2}^N bt(n)(e^{K ln n}-1) n^{-(sigma-y)},  K = t y/(2(x_N-6)).
   Prints the certified ball for floor.  usage: floor_arb N tn td yn yd T|S prec p1 p2 ... */
#include <stdio.h>
#include <stdlib.h>
#include <flint/arb.h>
typedef unsigned long long u64;
static int K; static u64 P[16]; static u64 NN; static slong prec; static char wt;
static arb_ptr bt, ky, rpow; static arb_t psig[16], g, c, S, tmp1, tmp2, bb, aa, W, wgt;
static arb_t lamd[1<<12]; static u64 dval[1<<12];
static u64 cnt=0;
static void visit(u64 r, u64 s, const arb_t spow, int radmask){
  u64 n=s*r; if(n==1) return;
  arb_zero(bb); arb_zero(aa);
  int m=radmask;
  for(;;){ u64 d=dval[m]; if(n%d==0){ u64 k=n/d; if(k<=NN){ arb_mul(tmp1, lamd[m], bt+k, prec); arb_add(bb,bb,tmp1,prec); arb_addmul(aa,tmp1,ky+k,prec);} }
           if(m==0) break; m=(m-1)&radmask; }
  if(wt=='T'){ arb_abs(tmp1,bb); arb_abs(tmp2,aa); arb_mul(tmp2,tmp2,g,prec); arb_add(W,tmp1,tmp2,prec); }
  else { arb_mul(tmp2,g,aa,prec); arb_sub(tmp1,bb,tmp2,prec); arb_abs(tmp1,tmp1);
         arb_add(tmp2,bb,tmp2,prec); arb_abs(tmp2,tmp2); arb_mul(tmp2,tmp2,c,prec); arb_max(W,tmp1,tmp2,prec); }
  arb_mul(wgt, rpow+r, spow, prec); arb_addmul(S, W, wgt, prec); cnt++;
}
static void dfs(int i, u64 r, u64 s, const arb_t spow, int radmask, u64 sp, u64 lim){
  if(i==K){ visit(r,s,spow,radmask); return; }
  dfs(i+1,r,s,spow,radmask,sp,lim);
  arb_t sp2w; arb_init(sp2w);
  u64 s2=s*P[i]; u64 sp2=sp; arb_mul(sp2w, spow, psig[i], prec);
  for(;;){
    dfs(i+1,r,s2,sp2w,radmask|(1<<i),sp2,lim);
    if(sp2 > lim/P[i]) break;
    sp2*=P[i]; s2*=P[i]; arb_mul(sp2w, sp2w, psig[i], prec);
    if(sp2>lim) break;
  }
  arb_clear(sp2w);
}
int main(int argc,char**argv){
  NN=strtoull(argv[1],0,10); u64 tn=strtoull(argv[2],0,10), td=strtoull(argv[3],0,10), yn=strtoull(argv[4],0,10), yd=strtoull(argv[5],0,10);
  wt=argv[6][0]; prec=atol(argv[7]); K=0; for(int i=8;i<argc;i++) P[K++]=strtoull(argv[i],0,10);
  arb_t t,y,pi,xN,L,sig,inner,lk,one,modmoll,corr,fl,modK,smy;
  arb_init(t);arb_init(y);arb_init(pi);arb_init(xN);arb_init(L);arb_init(sig);arb_init(inner);arb_init(lk);arb_init(one);arb_init(modmoll);arb_init(corr);arb_init(fl);arb_init(modK);arb_init(smy);
  arb_init(g);arb_init(c);arb_init(S);arb_init(tmp1);arb_init(tmp2);arb_init(bb);arb_init(aa);arb_init(W);arb_init(wgt);
  arb_set_ui(t,tn); arb_div_ui(t,t,td,prec); arb_set_ui(y,yn); arb_div_ui(y,y,yd,prec); arb_sqrt(y,y,prec);
  arb_const_pi(pi,prec);
  arb_set_ui(xN,NN); arb_mul(xN,xN,xN,prec); arb_mul(xN,xN,pi,prec); arb_mul_2exp_si(xN,xN,2);
  arb_mul(tmp1,pi,t,prec); arb_mul_2exp_si(tmp1,tmp1,-2); arb_sub(xN,xN,tmp1,prec);
  arb_div(L,xN,pi,prec); arb_mul_2exp_si(L,L,-2); arb_log(L,L,prec);
  /* inner = 1-3y+4y(1+y)/xN^2 */
  arb_mul(tmp2,xN,xN,prec); arb_add_ui(tmp1,y,1,prec); arb_mul(tmp1,tmp1,y,prec); arb_mul_2exp_si(tmp1,tmp1,2); arb_div(tmp1,tmp1,tmp2,prec);
  arb_mul_ui(inner,y,3,prec); arb_neg(inner,inner); arb_add_ui(inner,inner,1,prec); arb_add(inner,inner,tmp1,prec);
  if(arb_is_negative(inner)) arb_zero(inner); else if(!arb_is_positive(inner)){ printf("inner sign uncertain\n"); return 3; }
  arb_mul(tmp1,t,inner,prec); arb_mul_2exp_si(tmp2,tmp2,1); arb_div(tmp1,tmp1,tmp2,prec);
  arb_mul(sig,t,L,prec); arb_mul_2exp_si(sig,sig,-2); arb_sub(sig,sig,tmp1,prec);
  arb_add_ui(tmp1,y,1,prec); arb_mul_2exp_si(tmp1,tmp1,-1); arb_add(sig,sig,tmp1,prec);
  /* g = exp(0.02y - y/2 L) */
  arb_mul(tmp1,y,L,prec); arb_mul_2exp_si(tmp1,tmp1,-1); arb_set_ui(tmp2,2); arb_div_ui(tmp2,tmp2,100,prec); arb_mul(tmp2,tmp2,y,prec);
  arb_sub(tmp1,tmp2,tmp1,prec); arb_exp(g,tmp1,prec);
  arb_sub_ui(tmp1,g,1,prec); arb_neg(tmp1,tmp1); arb_add_ui(tmp2,g,1,prec); arb_div(c,tmp1,tmp2,prec);
  /* tables */
  bt=_arb_vec_init(NN+1); ky=_arb_vec_init(NN+1); rpow=_arb_vec_init(NN+1);
  for(u64 k=1;k<=NN;k++){ arb_set_ui(lk,k); arb_log(lk,lk,prec);
    arb_mul(tmp1,lk,lk,prec); arb_mul(tmp1,tmp1,t,prec); arb_mul_2exp_si(tmp1,tmp1,-2); arb_exp(bt+k,tmp1,prec);
    arb_mul(tmp1,lk,y,prec); arb_exp(ky+k,tmp1,prec);
    arb_mul(tmp1,lk,sig,prec); arb_neg(tmp1,tmp1); arb_exp(rpow+k,tmp1,prec); }
  arb_one(modmoll);
  arb_t lamp[16];
  for(int i=0;i<K;i++){ arb_init(lamp[i]); arb_init(psig[i]); arb_set_ui(lk,P[i]); arb_log(lk,lk,prec);
    arb_mul(tmp1,lk,lk,prec); arb_mul(tmp1,tmp1,t,prec); arb_mul_2exp_si(tmp1,tmp1,-2); arb_exp(lamp[i],tmp1,prec);
    arb_mul(tmp1,lk,sig,prec); arb_neg(tmp1,tmp1); arb_exp(psig[i],tmp1,prec);
    arb_mul(tmp1,lamp[i],psig[i],prec); arb_add_ui(tmp1,tmp1,1,prec); arb_mul(modmoll,modmoll,tmp1,prec); }
  for(int m=0;m<(1<<K);m++){ arb_init(lamd[m]); arb_one(lamd[m]); u64 d=1;
    for(int i=0;i<K;i++) if(m&(1<<i)){ d*=P[i]; arb_mul(lamd[m],lamd[m],lamp[i],prec); arb_neg(lamd[m],lamd[m]); }
    dval[m]=d; }
  arb_zero(S);
  arb_t one_; arb_init(one_); arb_one(one_);
  for(u64 r=1;r<=NN;r++){ int ok=1; for(int i=0;i<K;i++) if(r%P[i]==0){ok=0;break;} if(!ok) continue;
    dfs(0,r,1,one_,0,1,NN/r); }
  /* corr */
  arb_mul(modK,t,y,prec); arb_sub_ui(tmp1,xN,6,prec); arb_mul_2exp_si(tmp1,tmp1,1); arb_div(modK,modK,tmp1,prec);
  arb_sub(smy,sig,y,prec); arb_zero(corr);
  for(u64 n=2;n<=NN;n++){ arb_set_ui(lk,n); arb_log(lk,lk,prec);
    arb_mul(tmp1,lk,modK,prec); arb_expm1(tmp1,tmp1,prec); arb_mul(tmp1,tmp1,bt+n,prec);
    arb_mul(tmp2,lk,smy,prec); arb_neg(tmp2,tmp2); arb_exp(tmp2,tmp2,prec); arb_addmul(corr,tmp1,tmp2,prec); }
  arb_mul(corr,corr,g,prec);
  arb_one(fl); arb_sub(fl,fl,g,prec); arb_sub(fl,fl,S,prec); arb_div(fl,fl,modmoll,prec); arb_sub(fl,fl,corr,prec);
  printf("N=%llu t=%llu/%llu y2=%llu/%llu W=%c prec=%ld primes=%d terms=%llu\n",NN,tn,td,yn,yd,wt,(long)prec,K,cnt);
  printf("sigma="); arb_printn(sig,20,0); printf("\ng="); arb_printn(g,20,0); printf("\nmodmoll="); arb_printn(modmoll,20,0);
  printf("\nS="); arb_printn(S,20,0); printf("\ncorr="); arb_printn(corr,10,0); printf("\nFLOOR="); arb_printn(fl,20,0);
  arf_t lo; arf_init(lo); arb_get_lbound_arf(lo,fl,prec); printf("\nFLOOR_LOWER="); arf_printd(lo,15); printf("\npositive=%d\n", arb_is_positive(fl));
  return 0;
}
