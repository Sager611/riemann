/* Exhaustive small-graph checks for OpenAI's Seymour second-neighbourhood paper.
   Reads oriented graphs in nauty directg -T format (n ne v1 w1 ...) on stdin.
   For every graph:
     (a) Seymour's statement: some v with |N1+(v)| <= |N2+(v)|;
     (b) if all indegrees > 0: Proposition 3.1 claims NO graph has
           |U| + |F^2 U| < 2|F U|  for every nonempty proper U;
         we count graphs where it holds (must be 0) and record the best slack;
     (c) the subset-deficit property |F^2 S \ F S| < |F S \ S| for all
         nonempty proper S (Lemma 1.2 conclusion; must never hold, since
         singletons make the graph a counterexample).
   usage: geng -q N | directg -q -o -T | seymour_graphs */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int main(void){
  int n, ne;
  unsigned long long graphs=0, posin=0, ssnc_fail=0, concave=0, deficit=0;
  int best_slack=-1000; /* max over graphs of min_U (2|FU|-|U|-|F^2U|) - 1; concavity <=> min >= 1 */
  int nmax=0;
  static unsigned F[1<<10];
  while(scanf("%d %d",&n,&ne)==2){
    unsigned out[16]={0}; int indeg[16]={0};
    for(int e=0;e<ne;e++){ int a,b; if(scanf("%d %d",&a,&b)!=2) return 1; out[a]|=1u<<b; indeg[b]++; }
    graphs++; if(n>nmax) nmax=n;
    /* (a) */
    int good=0;
    for(int v=0;v<n && !good;v++){
      unsigned N1=out[v], F2=0;
      for(int u=0;u<n;u++) if((N1>>u)&1) F2|=out[u];
      unsigned N2=F2 & ~N1 & ~(1u<<v);
      if(__builtin_popcount(N1)<=__builtin_popcount(N2)) good=1;
    }
    if(!good){ ssnc_fail++; printf("SSNC-FAIL n=%d\n",n); }
    int allpos=1; for(int v=0;v<n;v++) if(!indeg[v]) allpos=0;
    unsigned full=(1u<<n)-1;
    F[0]=0;
    for(unsigned U=1;U<=full;U++){ int lb=__builtin_ctz(U); F[U]=F[U&(U-1)]|out[lb]; }
    if(allpos && n>=1){
      posin++;
      int minsl=1000;
      for(unsigned U=1;U<full;U++){
        int sl=2*__builtin_popcount(F[U])-__builtin_popcount(U)-__builtin_popcount(F[F[U]]);
        if(sl<minsl) minsl=sl;
      }
      if(n==1) minsl=1000; /* no nonempty proper subsets; but n=1 cannot have positive indegree */
      if(minsl>=1){ concave++; printf("CONCAVITY-GRAPH n=%d ne=%d\n",n,ne); }
      if(minsl>best_slack) best_slack=minsl;
    }
    /* (c) */
    int all=1;
    for(unsigned S=1;S<full && all;S++){
      unsigned FS=F[S], F2S=F[FS];
      if(!(__builtin_popcount(F2S & ~FS) < __builtin_popcount(FS & ~S))) all=0;
    }
    if(n>=2 && all){ deficit++; printf("DEFICIT-GRAPH n=%d\n",n); }
  }
  printf("graphs=%llu (n<=%d) positive_indegree=%llu ssnc_fail=%llu concavity_graphs=%llu deficit_graphs=%llu best_min_slack=%d\n",
    graphs,nmax,posin,ssnc_fail,concave,deficit,best_slack);
  return 0;
}
