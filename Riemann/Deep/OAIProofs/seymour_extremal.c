/* Test of the extremal argument (Section 3, Proposition 3.1) of OpenAI's
   Seymour paper on small oriented graphs with all indegrees positive.
   For compatible pairs (P,Q) containing the diagonal we check, with our own code:
    C1 every column of P and every row of Q is a proper subset of X
    C2 the diagonal lies in R = X^2 \ F2^2 Q and in C = X^2 \ F1^2 P
    C3 diagonal members of R and C are conflict-free
    C4 H is disjoint from F1 P u F2 Q   (so |H| <= d)
    C5 pruning (R,C) by the Konig cover of a maximum matching of the augmented
       graph gives compatible (P',Q') containing the diagonal with
       |P'|+|Q'|+d' >= |R|+|C|-d, and strictly more than M+d whenever the
       summed concavity defect  sum_j (|P_j|+|F^2 P_j|-2|F P_j|) + rows(Q) < 0
    C6 at a global maximiser of M+d (exhaustive for n<=4) the summed defect is >= 0
   usage: geng -q N | directg -q -o -T | seymour_extremal MODE SAMPLES SEED
          MODE = exh (all pairs, n<=3), max (maximiser + samples, n=4), rnd (samples) */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef uint64_t mask;
static int n, NN; static unsigned out[16];
static int rel(int a,int b){ return (out[a]>>b)&1; }
static mask img1[64], img2[64], DIAG, ALL;
static int pc(mask m){ return __builtin_popcountll(m); }
static mask F1(mask S){ mask r=0; while(S){ int c=__builtin_ctzll(S); S&=S-1; r|=img1[c]; } return r; }
static mask F2(mask S){ mask r=0; while(S){ int c=__builtin_ctzll(S); S&=S-1; r|=img2[c]; } return r; }
static unsigned long long cnt_pairs=0, f1=0,f2=0,f3=0,f4=0,f5=0,f5s=0,f6=0, negS=0, maxchecked=0;
static uint64_t rs=0x9E3779B97F4A7C15ULL; static uint64_t xr(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return rs; }

/* ---- pruning via augmented graph and Konig ---- */
static mask confR[64]; static mask Zm,Hm,conflictedR,conflictedC;
static mask adjlo[128], adjhi[128]; static int matchL[128], matchR[128], seen[128];
static int hasedge(int l,int r){ return r<64 ? (adjlo[l]>>r)&1 : (adjhi[l]>>(r-64))&1; }
static void addedge(int l,int r){ if(r<64) adjlo[l]|=(mask)1<<r; else adjhi[l]|=(mask)1<<(r-64); }
static int kuhn(int l){ for(int r=0;r<2*NN;r++) if(hasedge(l,r)&&!seen[r]){ seen[r]=1; if(matchR[r]<0||kuhn(matchR[r])){ matchR[r]=l; matchL[l]=r; return 1; } } return 0; }
static mask cover_pts(mask Rp, mask Cp){ return F1(Rp)|F2(Cp); }
static int conflict_free(mask Rp, mask Cp){ for(int a=0;a<NN;a++) if((Rp>>a)&1) if(confR[a]&Cp) return 0; return 1; }
/* returns deleted masks */
static void prune(mask Rm, mask Cm, mask*Rp, mask*Cp){
  Zm=Hm=conflictedR=conflictedC=0;
  for(int a=0;a<NN;a++){ confR[a]=0; if(!((Rm>>a)&1)) continue; int p=a/n,j=a%n;
    for(int b=0;b<NN;b++) if((Cm>>b)&1){ int i=b/n,s=b%n; if(rel(p,i)&&rel(s,j)){ confR[a]|=(mask)1<<b; conflictedR|=(mask)1<<a; conflictedC|=(mask)1<<b; Zm|=(mask)1<<(i*n+j); Hm|=(mask)1<<(p*n+s);} } }
  if(!conflictedR){ *Rp=Rm; *Cp=Cm; return; }
  for(int l=0;l<2*NN;l++){ adjlo[l]=adjhi[l]=0; matchL[l]=matchR[l]=-1; }
  for(int a=0;a<NN;a++) if((Rm>>a)&1){ for(int b=0;b<NN;b++) if((confR[a]>>b)&1) addedge(a,b);
    int p=a/n,j=a%n; for(int i=0;i<n;i++) if(rel(p,i)){ int z=i*n+j; if((Zm>>z)&1) addedge(a,NN+z);} }
  for(int b=0;b<NN;b++) if((Cm>>b)&1){ int i=b/n,s=b%n; for(int j=0;j<n;j++) if(rel(s,j)){ int z=i*n+j; if((Zm>>z)&1) addedge(NN+z,b);} }
  for(int l=0;l<2*NN;l++){ int pres = l<NN ? (int)((Rm>>l)&1) : (int)((Zm>>(l-NN))&1); if(!pres) continue; memset(seen,0,sizeof(seen)); kuhn(l); }
  int reachL[128]={0}, reachR[128]={0}, st[256], sp=0;
  for(int l=0;l<2*NN;l++){ int pres = l<NN ? (int)((Rm>>l)&1) : (int)((Zm>>(l-NN))&1); if(pres&&matchL[l]<0){ reachL[l]=1; st[sp++]=l; } }
  while(sp){ int l=st[--sp]; for(int r=0;r<2*NN;r++) if(hasedge(l,r)&&!reachR[r]&&matchL[l]!=r){ reachR[r]=1; int l2=matchR[r]; if(l2>=0&&!reachL[l2]){ reachL[l2]=1; st[sp++]=l2; } } }
  mask dR=0,dC=0;
  for(int l=0;l<NN;l++) if(((Rm>>l)&1)&&!reachL[l]) dR|=(mask)1<<l;
  for(int r=0;r<NN;r++) if(((Cm>>r)&1)&&reachR[r]) dC|=(mask)1<<r;
  *Rp=Rm&~dR; *Cp=Cm&~dC;
}

static void check_pair(mask P, mask Q, int is_max){
  cnt_pairs++;
  mask FP=F1(P), FQ=F2(Q);
  if(FP&FQ){ printf("internal: not compatible\n"); return; }
  /* C1 */
  for(int j=0;j<n;j++){ int col=0; for(int p=0;p<n;p++) if((P>>(p*n+j))&1) col++; if(col==n){ f1++; } }
  for(int i=0;i<n;i++){ int row=0; for(int s=0;s<n;s++) if((Q>>(i*n+s))&1) row++; if(row==n){ f1++; } }
  mask F2P=F1(FP), F2Q=F2(FQ);
  mask R=ALL&~F2Q, C=ALL&~F2P;
  if((DIAG&R)!=DIAG || (DIAG&C)!=DIAG) f2++;
  mask Rp,Cp; prune(R,C,&Rp,&Cp);
  /* C3: diagonal members conflict-free */
  for(int v=0;v<n;v++){ int dv=v*n+v; if(confR[dv]) f3++; if((conflictedC>>dv)&1) f3++; }
  /* C4 */
  if(Hm & (FP|FQ)) f4++;
  int d = NN - pc(FP) - pc(FQ);
  int M = pc(P)+pc(Q);
  /* C5 */
  int ok=1;
  if(!conflict_free(Rp,Cp)) ok=0;
  if((F1(Rp)&F2(Cp))) ok=0;
  if((Rp&DIAG)!=DIAG || (Cp&DIAG)!=DIAG) ok=0;
  int dp = NN - pc(F1(Rp)) - pc(F2(Cp));
  int unc = pc(Zm & ~cover_pts(Rp,Cp));
  if(pc(R&~Rp)+pc(C&~Cp) > pc(Hm)+unc) ok=0;
  if(pc(Rp)+pc(Cp)+dp < pc(R)+pc(C)-d) ok=0;
  if(!ok){ f5++; if(f5<5) printf("C5 fail\n"); }
  int Ssum = pc(P)+pc(F2P)-2*pc(FP) + pc(Q)+pc(F2Q)-2*pc(FQ);
  if(Ssum<0){ negS++; if(!(pc(Rp)+pc(Cp)+dp > M+d)){ f5s++; } }
  if(is_max){ maxchecked++; if(Ssum<0){ f6++; printf("C6 fail: maximiser with negative summed defect\n"); } }
}

int main(int argc,char**argv){
  const char*mode=argv[1]; long samples=atol(argv[2]); rs^=(uint64_t)atoll(argv[3])*0x2545F4914F6CDD1DULL; if(!rs) rs=1;
  int ne; unsigned long long graphs=0; int nmax=0;
  while(scanf("%d %d",&n,&ne)==2){
    memset(out,0,sizeof(out)); int indeg[16]={0};
    for(int e=0;e<ne;e++){ int a,b; if(scanf("%d %d",&a,&b)!=2) return 1; out[a]|=1u<<b; indeg[b]++; }
    int allpos=1; for(int v=0;v<n;v++) if(!indeg[v]) allpos=0;
    if(!allpos) continue;
    graphs++; if(n>nmax) nmax=n;
    NN=n*n; DIAG=0; ALL = NN==64? ~(mask)0 : (((mask)1<<NN)-1);
    for(int v=0;v<n;v++) DIAG|=(mask)1<<(v*n+v);
    for(int c=0;c<NN;c++){ int a=c/n,b=c%n; img1[c]=img2[c]=0;
      for(int i=0;i<n;i++) if(rel(a,i)) img1[c]|=(mask)1<<(i*n+b);   /* (a,b)->(i,b) */
      for(int j=0;j<n;j++) if(rel(b,j)) img2[c]|=(mask)1<<(a*n+j);   /* (a,b)->(a,j) */ }
    mask offd = ALL&~DIAG; int no=pc(offd); int offidx[64]; { int k=0; for(int c=0;c<NN;c++) if((offd>>c)&1) offidx[k++]=c; }
    if(!strcmp(mode,"exh")||!strcmp(mode,"max")){
      /* enumerate all P,Q containing diagonal */
      long np = 1L<<no;
      mask *Ps=malloc(sizeof(mask)*np), *FPs=malloc(sizeof(mask)*np), *FQs=malloc(sizeof(mask)*np);
      for(long t=0;t<np;t++){ mask S=DIAG; for(int k=0;k<no;k++) if((t>>k)&1) S|=(mask)1<<offidx[k]; Ps[t]=S; FPs[t]=F1(S); FQs[t]=F2(S); }
      int best=-1; long bp=-1,bq=-1; long ncompat=0;
      for(long a=0;a<np;a++) for(long b=0;b<np;b++) if(!(FPs[a]&FQs[b])){
        ncompat++;
        int val = pc(Ps[a])+pc(Ps[b]) + NN - pc(FPs[a]) - pc(FQs[b]);
        if(val>best){ best=val; bp=a; bq=b; }
        if(!strcmp(mode,"exh")) check_pair(Ps[a],Ps[b],0);
      }
      /* all maximisers */
      for(long a=0;a<np;a++) for(long b=0;b<np;b++) if(!(FPs[a]&FQs[b])){
        int val = pc(Ps[a])+pc(Ps[b]) + NN - pc(FPs[a]) - pc(FQs[b]);
        if(val==best) check_pair(Ps[a],Ps[b],1);
      }
      if(!strcmp(mode,"max")){
        for(long t=0;t<samples;t++){ long a=xr()%np, b=xr()%np; if(!(FPs[a]&FQs[b])) check_pair(Ps[a],Ps[b],0); }
      }
      free(Ps); free(FPs); free(FQs);
    } else {
      /* random compatible pairs grown greedily from the diagonal */
      for(long t=0;t<samples;t++){
        mask P=DIAG,Q=DIAG; int tries = xr()%(2*NN+1);
        for(int k=0;k<tries;k++){ int c=offidx[xr()%no]; if(xr()&1){ mask P2=P|((mask)1<<c); if(!(F1(P2)&F2(Q))) P=P2; } else { mask Q2=Q|((mask)1<<c); if(!(F1(P)&F2(Q2))) Q=Q2; } }
        check_pair(P,Q,0);
      }
    }
  }
  printf("mode=%s graphs(pos indeg, n<=%d)=%llu pairs=%llu C1fail=%llu C2fail=%llu C3fail=%llu C4fail=%llu C5fail=%llu negDefect=%llu strictImproveFail=%llu maximisers=%llu C6fail=%llu\n",
     mode,nmax,graphs,cnt_pairs,f1,f2,f3,f4,f5,negS,f5s,maxchecked,f6);
  return 0;
}
