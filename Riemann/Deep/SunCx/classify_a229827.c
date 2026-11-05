// classify_a229827.c -- exhaustive check of Sun's OEIS A229827 conjecture (i) for q = p^k <= 9.
// Every P: F_q -> F_q with sum_x P(x) = 0 (i.e. deg P < q-1) is reduced modulo P -> g P g^-1 + c (g in AGL(k,p) acting on the
// additive group, c a constant); for each orbit a complete search decides whether some circular order (a_1..a_q) has all
// sums P(a_i)+a_{i+1} distinct.  Build: cc -O2 -o classify_a229827 classify_a229827.c ; run: ./classify_a229827 p k
// Expected (2026-11-05):  q=2: 1 orbit, 0 failing | q=3: 3, only c-x | q=4: 6, c-x + 1 orbit (8 fns: w x + c, w^2 x + c)
//   q=5: 16, only c-x | q=7: 477, only c-x | q=8: 312, c-x + 10 orbits (20152 fns = the even permutations except x+c)
//   q=9: 11520 orbits (10.5 s, 50 MB), c-x + 5 orbits (648 odd permutations; sizes 108, 216, 216, 54, 54), max 7 distinct sums.
// Classify all P: F_q -> F_q with sum_x P(x) = 0 (i.e. deg P < q-1), q = p^k <= 9,
// up to P -> g P g^{-1} + c (g in AGL(k,p) acting on the additive group F_p^k), and decide whether some
// circular permutation (a_1..a_q) of F_q has all q sums P(a_i) + a_{i+1} distinct.
// usage: cls p k     (elements = base-p digit vectors; the problem only depends on (F_q,+))
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int p, k, q;
static int ad[9][9];
static int nG = 0;
static int G[4000][9], Gi[4000][9];
static int digit(int x, int i){ for(int j=0;j<i;j++) x/=p; return x%p; }
static int addv(int a, int b){ int r=0,m=1; for(int i=0;i<k;i++){ r += ((a%p + b%p)%p)*m; a/=p; b/=p; m*=p; } return r; }
static int P[9], A[9];
static uint32_t usedV, usedC;
static long nodes;
static int dfs(int depth, int x){
  nodes++;
  if(depth==q){ int c = ad[P[x]][A[0]]; return !((usedC>>c)&1); }
  for(int y=0;y<q;y++){
    if((usedV>>y)&1) continue;
    int c = ad[P[x]][y];
    if((usedC>>c)&1) continue;
    usedV |= 1u<<y; usedC |= 1u<<c; A[depth]=y;
    if(dfs(depth+1,y)) return 1;
    usedV &= ~(1u<<y); usedC &= ~(1u<<c);
  }
  return 0;
}
static int best;
static void dfsmax(int depth, int x, int cnt){ // max number of distinct sums over all circular orders
  if(depth==q){ int c=ad[P[x]][A[0]]; int t = cnt + !((usedC>>c)&1); if(t>best) best=t; return; }
  if(cnt + (q-depth) + 1 <= best) return;
  for(int y=0;y<q;y++){
    if((usedV>>y)&1) continue;
    int c = ad[P[x]][y]; int nw = !((usedC>>c)&1);
    usedV |= 1u<<y; uint32_t oc=usedC; usedC |= 1u<<c; A[depth]=y;
    dfsmax(depth+1,y,cnt+nw);
    usedV &= ~(1u<<y); usedC = oc;
  }
}
int main(int argc, char**argv){
  p = atoi(argv[1]); k = atoi(argv[2]); q = 1; for(int i=0;i<k;i++) q*=p;
  for(int a=0;a<q;a++) for(int b=0;b<q;b++) ad[a][b]=addv(a,b);
  // GL(k,p): all k x k matrices, keep bijective ones
  int nm = 1; for(int i=0;i<k*k;i++) nm*=p;
  static int GL[200][9]; int nGL=0;
  for(int m=0;m<nm;m++){
    int M[3][3]; int t=m; for(int i=0;i<k;i++) for(int j=0;j<k;j++){ M[i][j]=t%p; t/=p; }
    int img[9]; uint32_t seen=0; int ok=1;
    for(int x=0;x<q;x++){ int r=0,mm=1; for(int i=0;i<k;i++){ int s=0; for(int j=0;j<k;j++) s+=M[i][j]*digit(x,j); r+=(s%p)*mm; mm*=p; } img[x]=r; if((seen>>r)&1) ok=0; seen|=1u<<r; }
    if(ok){ memcpy(GL[nGL++], img, sizeof img); }
  }
  for(int g=0;g<nGL;g++) for(int v=0;v<q;v++){ for(int x=0;x<q;x++){ G[nG][x]=ad[GL[g][x]][v]; } for(int x=0;x<q;x++) Gi[nG][G[nG][x]]=x; nG++; }
  long long N=1; for(int i=0;i<q;i++) N*=q;
  uint8_t *vis = calloc((N+7)/8,1);
  long long pw[10]; pw[0]=1; for(int i=1;i<=q;i++) pw[i]=pw[i-1]*q;
  printf("q=%d p=%d k=%d |GL|=%d |AGL|=%d functions=%lld\n", q,p,k,nGL,nG,N);
  long long norb=0, nfail=0, tot=0; int digs[9]={0};
  for(long long code=0; code<N; code++){
    if(code){ int i=0; while(++digs[i]==q){ digs[i]=0; i++; } }
    if((vis[code>>3]>>(code&7))&1) continue;
    int s=0; for(int x=0;x<q;x++) s=ad[s][digs[x]];
    if(s!=0) continue;
    // new orbit
    long long osz=0;
    for(int g=0; g<nG; g++){
      int Q[9]; for(int x=0;x<q;x++) Q[x]=G[g][digs[Gi[g][x]]];
      for(int c=0;c<q;c++){ long long cd=0; for(int x=0;x<q;x++) cd += pw[x]*ad[Q[x]][c];
        if(!((vis[cd>>3]>>(cd&7))&1)){ vis[cd>>3] |= 1u<<(cd&7); osz++; } }
    }
    tot += osz; norb++;
    // excluded P = c - x ?
    int excl=1; for(int x=1;x<q;x++) if(ad[digs[x]][x]!=ad[digs[0]][0]) excl=0;
    for(int x=0;x<q;x++) P[x]=digs[x];
    A[0]=0; usedV=1; usedC=0; nodes=0;
    int ok = dfs(1,0);
    if(!ok){
      best=0; A[0]=0; usedV=1; usedC=0; dfsmax(1,0,0);
      nfail++;
      printf("%s q=%d P=[", excl?"EXCLUDED":"FAIL", q); for(int x=0;x<q;x++) printf("%d%s",P[x],x<q-1?",":""); printf("] orbit=%lld maxdistinct=%d\n", osz, best);
    }
  }
  printf("q=%d orbits=%lld functions_with_sum0=%lld (expected %lld) nonsolvable_orbits=%lld\n", q, norb, tot, N/q, nfail);
  free(vis); return 0;
}

