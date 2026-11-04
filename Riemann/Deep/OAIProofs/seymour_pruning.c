/* Independent test of the Pruning Lemma (Lemma 2.1) of OpenAI's
   "A proof of Seymour's second-neighborhood conjecture" (Sep 23 2026).
   For a finite set X, ANY binary relation ->, and families R, C of ordered
   pairs, the lemma claims subfamilies R', C' with no conflicts, keeping every
   initially conflict-free member, and
      |R\R'| + |C\C'| <= |H| + #{z in Z : z covered by no member of R' u C'}.
   Checks per instance:
     (1) matching bound nu(Gamma) <= |Z| + |H|   (the paper's key claim)
     (2) the Konig cover of a maximum matching gives a valid witness (R',C')
     (3) [mode exact] brute-force optimum over all deletion sets.
   usage: seymour_pruning exh N            exhaustive over all relations,R,C (N<=3)
          seymour_pruning rnd N COUNT SEED random instances
          seymour_pruning exact N COUNT SEED random + exact optimum (N<=3)  */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

static int n, NN;               /* |X|, |X|^2 */
static int rel[8][8];
typedef uint64_t mask;

static int popc(mask m){ return __builtin_popcountll(m); }

/* augmented bipartite graph: left = R members (pair idx) + Z_L copies,
   right = C members + Z_R copies. Vertex ids: left 0..NN-1 = R pair idx,
   NN..2NN-1 = Z_L copy of pair idx; right likewise C / Z_R. */
static mask adjL[128];          /* adjacency of left vertex -> right vertices (128 bits split) */
static mask adjL2[128];         /* high part (right ids >= 64) */
static int matchR[128], matchL[128], seen[128];
static int nleft, nright;

static int hasedge(int l, int r){ return r<64 ? (adjL[l]>>r)&1 : (adjL2[l]>>(r-64))&1; }

static int try_kuhn(int l){
  for(int r=0;r<2*NN;r++) if(hasedge(l,r) && !seen[r]){
    seen[r]=1;
    if(matchR[r]<0 || try_kuhn(matchR[r])){ matchR[r]=l; matchL[l]=r; return 1; }
  }
  return 0;
}

static unsigned long long cnt_inst=0, cnt_conf=0, fail_bound=0, fail_konig=0, fail_exact=0, tight_bound=0, tight_exact=0;

static mask Rm, Cm, Zm, Hm, conflictedR, conflictedC;
static mask confR[64];          /* for R member a: set of C members conflicting with it */

static void build(void){
  Zm=Hm=conflictedR=conflictedC=0;
  for(int a=0;a<NN;a++){ confR[a]=0; }
  for(int a=0;a<NN;a++) if((Rm>>a)&1){
    int p=a/n, j=a%n;
    for(int b=0;b<NN;b++) if((Cm>>b)&1){
      int i=b/n, s=b%n;
      if(rel[p][i] && rel[s][j]){
        confR[a] |= (mask)1<<b;
        conflictedR |= (mask)1<<a; conflictedC |= (mask)1<<b;
        Zm |= (mask)1<<(i*n+j); Hm |= (mask)1<<(p*n+s);
      }
    }
  }
}

/* points covered by members R', C' */
static mask covered(mask Rp, mask Cp){
  mask cov=0;
  for(int a=0;a<NN;a++) if((Rp>>a)&1){ int p=a/n,j=a%n; for(int i=0;i<n;i++) if(rel[p][i]) cov|=(mask)1<<(i*n+j); }
  for(int b=0;b<NN;b++) if((Cp>>b)&1){ int i=b/n,s=b%n; for(int j=0;j<n;j++) if(rel[s][j]) cov|=(mask)1<<(i*n+j); }
  return cov;
}

static int conflict_free(mask Rp, mask Cp){
  for(int a=0;a<NN;a++) if((Rp>>a)&1) if(confR[a] & Cp) return 0;
  return 1;
}

static void print_instance(const char*tag){
  printf("%s n=%d rel:", tag, n);
  for(int a=0;a<n;a++) for(int b=0;b<n;b++) if(rel[a][b]) printf(" %d>%d",a,b);
  printf(" R:"); for(int a=0;a<NN;a++) if((Rm>>a)&1) printf(" (%d,%d)",a/n,a%n);
  printf(" C:"); for(int a=0;a<NN;a++) if((Cm>>a)&1) printf(" (%d,%d)",a/n,a%n);
  printf(" |Z|=%d |H|=%d\n", popc(Zm), popc(Hm));
}

/* exact optimum: min over K subset of conflicted members, K vertex cover of conflict graph,
   of |K| - #(Z uncovered by survivors). returns min value */
static int exact_opt(void){
  int ids[64], side[64], m=0;
  for(int a=0;a<NN;a++) if((conflictedR>>a)&1){ ids[m]=a; side[m]=0; m++; }
  for(int b=0;b<NN;b++) if((conflictedC>>b)&1){ ids[m]=b; side[m]=1; m++; }
  if(m>24){ return -1000; }
  int best=1<<30;
  for(uint32_t K=0; K < (1u<<m); K++){
    mask Rp=Rm, Cp=Cm; int k=0;
    for(int t=0;t<m;t++) if((K>>t)&1){ k++; if(side[t]==0) Rp&=~((mask)1<<ids[t]); else Cp&=~((mask)1<<ids[t]); }
    if(k - 0 >= best + popc(Zm)) continue; /* cannot improve: uncovered <= |Z| */
    if(!conflict_free(Rp,Cp)) continue;
    int unc = popc(Zm & ~covered(Rp,Cp));
    int val = k - unc;
    if(val<best) best=val;
  }
  return best;
}

static void check_instance(int do_exact){
  cnt_inst++;
  build();
  if(!conflictedR){ return; }
  cnt_conf++;
  /* build augmented graph */
  for(int l=0;l<2*NN;l++){ adjL[l]=adjL2[l]=0; }
  for(int a=0;a<NN;a++) if((Rm>>a)&1){
    /* conflicts to C members (right ids 0..NN-1) */
    for(int b=0;b<NN;b++) if((confR[a]>>b)&1){ int r=b; if(r<64) adjL[a]|=(mask)1<<r; else adjL2[a]|=(mask)1<<(r-64); }
    /* covers Z points -> Z_R copies (right ids NN + idx) */
    int p=a/n, j=a%n;
    for(int i=0;i<n;i++) if(rel[p][i]){ int z=i*n+j; if((Zm>>z)&1){ int r=NN+z; if(r<64) adjL[a]|=(mask)1<<r; else adjL2[a]|=(mask)1<<(r-64);} }
  }
  for(int b=0;b<NN;b++) if((Cm>>b)&1){
    int i=b/n, s=b%n;
    for(int j=0;j<n;j++) if(rel[s][j]){ int z=i*n+j; if((Zm>>z)&1){ int l=NN+z; int r=b; if(r<64) adjL[l]|=(mask)1<<r; else adjL2[l]|=(mask)1<<(r-64);} }
  }
  /* sanity: initially conflict-free members isolated */
  for(int a=0;a<NN;a++) if(((Rm>>a)&1) && !((conflictedR>>a)&1)) if(adjL[a]||adjL2[a]){ print_instance("ISOLATION-FAIL"); }
  for(int t=0;t<2*NN;t++){ matchR[t]=-1; matchL[t]=-1; }
  int nu=0;
  for(int l=0;l<2*NN;l++){
    int present = (l<NN) ? (int)((Rm>>l)&1) : (int)((Zm>>(l-NN))&1);
    if(!present) continue;
    memset(seen,0,sizeof(seen));
    if(try_kuhn(l)) nu++;
  }
  int zh = popc(Zm)+popc(Hm);
  if(nu>zh){ fail_bound++; if(fail_bound<=5) { print_instance("BOUND-FAIL"); printf("  nu=%d\n",nu);} }
  if(nu==zh) tight_bound++;
  /* Konig cover: alternating reachability from unmatched left vertices */
  int reachL[128]={0}, reachR[128]={0}, stack[256], sp=0;
  for(int l=0;l<2*NN;l++){
    int present = (l<NN) ? (int)((Rm>>l)&1) : (int)((Zm>>(l-NN))&1);
    if(present && matchL[l]<0){ reachL[l]=1; stack[sp++]=l; }
  }
  while(sp){ int l=stack[--sp];
    for(int r=0;r<2*NN;r++) if(hasedge(l,r) && !reachR[r] && matchL[l]!=r){ reachR[r]=1; int l2=matchR[r]; if(l2>=0 && !reachL[l2]){ reachL[l2]=1; stack[sp++]=l2; } }
  }
  mask delR=0, delC=0; int Ksize=0;
  for(int l=0;l<2*NN;l++){
    int present = (l<NN) ? (int)((Rm>>l)&1) : (int)((Zm>>(l-NN))&1);
    if(present && !reachL[l]){ Ksize++; if(l<NN) delR|=(mask)1<<l; }
  }
  for(int r=0;r<2*NN;r++){
    int present = (r<NN) ? (int)((Cm>>r)&1) : (int)((Zm>>(r-NN))&1);
    if(present && reachR[r]){ Ksize++; if(r<NN) delC|=(mask)1<<r; }
  }
  if(Ksize!=nu){ print_instance("KONIG-SIZE-FAIL"); }
  mask Rp=Rm&~delR, Cp=Cm&~delC;
  int ok = conflict_free(Rp,Cp);
  if((delR & ~conflictedR) || (delC & ~conflictedC)) ok=0;
  int unc = popc(Zm & ~covered(Rp,Cp));
  int deleted = popc(delR)+popc(delC);
  if(deleted > popc(Hm)+unc) ok=0;
  if(!ok){ fail_konig++; if(fail_konig<=5){ print_instance("KONIG-WITNESS-FAIL"); printf("  deleted=%d unc=%d\n",deleted,unc);} }
  if(do_exact){
    int best = exact_opt();
    if(best==-1000) return;
    if(best > popc(Hm)){ fail_exact++; if(fail_exact<=5){ print_instance("LEMMA-FAIL"); printf("  best=%d\n",best);} }
    if(best == popc(Hm)) tight_exact++;
  }
}

static uint64_t rs=88172645463325252ULL;
static uint64_t xr(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return rs; }

int main(int argc, char**argv){
  if(argc<3){ fprintf(stderr,"usage\n"); return 2; }
  n=atoi(argv[2]); NN=n*n;
  if(!strcmp(argv[1],"exh")){
    int nrel = 1<<NN;
    for(int rm=0; rm<nrel; rm++){
      for(int a=0;a<n;a++) for(int b=0;b<n;b++) rel[a][b]=(rm>>(a*n+b))&1;
      for(mask R=0; R<((mask)1<<NN); R++) for(mask C=0; C<((mask)1<<NN); C++){ Rm=R; Cm=C; check_instance(n<=2); }
    }
  } else {
    long count=atol(argv[3]); rs ^= (uint64_t)atoll(argv[4])*0x9E3779B97F4A7C15ULL; if(!rs) rs=1;
    int ex = !strcmp(argv[1],"exact");
    for(long t=0;t<count;t++){
      /* random density for relation and families */
      int dr = 1 + xr()%7, df = 1 + xr()%7; /* density in eighths */
      for(int a=0;a<n;a++) for(int b=0;b<n;b++) rel[a][b] = (int)(xr()%8) < dr;
      Rm=Cm=0;
      for(int a=0;a<NN;a++){ if((int)(xr()%8)<df) Rm|=(mask)1<<a; if((int)(xr()%8)<df) Cm|=(mask)1<<a; }
      check_instance(ex);
    }
  }
  printf("mode=%s n=%d instances=%llu with_conflict=%llu bound_fail=%llu bound_tight=%llu konig_fail=%llu exact_fail=%llu exact_tight=%llu\n",
     argv[1], n, cnt_inst, cnt_conf, fail_bound, tight_bound, fail_konig, fail_exact, tight_exact);
  return 0;
}
