// hyp.c - minimal hyperedges for symmetric 3-colourings with no 2-coloured K_k
// input (stdin or file): n m k nroots ; roots... ; then n*n ints cls[u][v] (-1 on diagonal), symmetric
// output: "m NE" then NE lines: list of class indices of a minimal hyperedge
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef unsigned __int128 u128;
static int n, m, K, nr, roots[256], C[160][160];
static u128 bit[128];
static u128 *hs; static size_t hcap, hcnt;
static void hput(u128 x);
static size_t hh(u128 x){ uint64_t a=(uint64_t)x, b=(uint64_t)(x>>64); uint64_t h=a*0x9E3779B97F4A7C15ULL ^ (b+0x632BE59BD9B4E019ULL)*0xC2B2AE3D27D4EB4FULL; h^=h>>29; return (size_t)h; }
static void hgrow(void){ size_t oc=hcap; u128 *o=hs; hcap*=2; hs=calloc(hcap,sizeof(u128)); hcnt=0; for(size_t i=0;i<oc;i++) if(o[i]) hput(o[i]); free(o);}
static void hput(u128 x){ size_t h=hh(x)&(hcap-1); while(hs[h]){ if(hs[h]==x) return; h=(h+1)&(hcap-1);} hs[h]=x; hcnt++; if(hcnt*2>hcap) hgrow(); }
static int sel[16], root;
static void rec(int depth, int start, u128 mask){
  if(depth==K-1){ hput(mask); return; }
  for(int x=start;x<n;x++){ if(x==root) continue;
    u128 nm = mask | bit[C[root][x]];
    for(int i=0;i<depth;i++) nm |= bit[C[sel[i]][x]];
    sel[depth]=x; rec(depth+1,x+1,nm);
  }
}
static int pc(u128 x){ return __builtin_popcountll((uint64_t)x)+__builtin_popcountll((uint64_t)(x>>64)); }
static int ctz(u128 x){ uint64_t a=(uint64_t)x; return a? __builtin_ctzll(a) : 64+__builtin_ctzll((uint64_t)(x>>64)); }
static int cmpm(const void*a,const void*b){ u128 x=*(u128*)a,y=*(u128*)b; int px=pc(x),py=pc(y); if(px!=py) return px-py; return x<y?-1:x>y; }
int main(int argc,char**argv){
  FILE*f= argc>1? fopen(argv[1],"r") : stdin; if(!f){perror("open");return 1;}
  if(fscanf(f,"%d %d %d %d",&n,&m,&K,&nr)!=4) return 1;
  if(n>160||m>128){ fprintf(stderr,"too big\n"); return 1; }
  for(int i=0;i<nr;i++) if(fscanf(f,"%d",&roots[i])!=1) return 1;
  for(int a=0;a<n;a++) for(int b=0;b<n;b++) if(fscanf(f,"%d",&C[a][b])!=1){fprintf(stderr,"bad matrix\n");return 1;}
  for(int i=0;i<128;i++) bit[i]=((u128)1)<<i;
  for(int a=0;a<n;a++) for(int b=0;b<n;b++) if(a!=b && (C[a][b]!=C[b][a] || C[a][b]<0 || C[a][b]>=m)){ fprintf(stderr,"bad class %d %d\n",a,b); return 1; }
  hcap=1<<16; hs=calloc(hcap,sizeof(u128)); hcnt=0;
  for(int r=0;r<nr;r++){ root=roots[r]; rec(0,0,0); }
  size_t N=0; u128 *L=malloc(hcnt*sizeof(u128)); for(size_t i=0;i<hcap;i++) if(hs[i]) L[N++]=hs[i];
  free(hs); qsort(L,N,sizeof(u128),cmpm);
  u128 *E=malloc(N*sizeof(u128)); size_t NE=0;
  u128 **lb=calloc(128,sizeof(u128*)); size_t *nlb=calloc(128,sizeof(size_t)), *clb=calloc(128,sizeof(size_t));
  for(size_t i=0;i<N;i++){
    u128 T=L[i]; int red=0;
    for(u128 t=T; t && !red; t&=t-1){ int b=ctz(t); for(size_t j=0;j<nlb[b];j++){ if((lb[b][j]&~T)==0){ red=1; break; } } }
    if(!red){ E[NE++]=T; int b=ctz(T); if(nlb[b]==clb[b]){ clb[b]=clb[b]?2*clb[b]:64; lb[b]=realloc(lb[b],clb[b]*sizeof(u128));} lb[b][nlb[b]++]=T; }
  }
  int minsz=999,maxsz=0; for(size_t i=0;i<NE;i++){int p=pc(E[i]); if(p<minsz)minsz=p; if(p>maxsz)maxsz=p;}
  fprintf(stderr,"n=%d m=%d K=%d roots=%d distinct=%zu minimal=%zu size %d..%d\n",n,m,K,nr,N,NE,minsz,maxsz);
  printf("%d %zu\n", m, NE);
  for(size_t i=0;i<NE;i++){ for(u128 t=E[i];t;t&=t-1) printf("%d ",ctz(t)); printf("\n"); }
  return 0;
}

