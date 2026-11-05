// cay.c - hyperedges + exhaustive solver for Cayley 3-colourings with no 2-coloured K_k
// usage: cay <tablefile> <k> <mode> [maxsol]
//   tablefile: "n" then n*n ints (mul[a][b], 0 = identity)
//   mode: "edges" -> print minimal hyperedges (class masks, hex) ; "solve" -> enumerate solutions
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
static int n, K, mul[64][64], inv_[64], cls[64], m, rep[64];
static uint64_t cm[64][64];
static uint64_t *hs; static size_t hcap, hcnt;
static void hput(uint64_t x);
static void hgrow(void){ size_t oc=hcap; uint64_t *o=hs; hcap*=2; hs=calloc(hcap,8); hcnt=0; for(size_t i=0;i<oc;i++) if(o[i]) hput(o[i]); free(o);}
static void hput(uint64_t x){ size_t h=(size_t)((x*0x9E3779B97F4A7C15ULL)>>17)&(hcap-1); while(hs[h]){ if(hs[h]==x) return; h=(h+1)&(hcap-1);} hs[h]=x; hcnt++; if(hcnt*2>hcap) hgrow(); }
static int sel[16];
static void rec(int depth, int start, uint64_t mask){
  if(depth==K-1){ hput(mask); return; }
  for(int x=start;x<n;x++){
    uint64_t nm = mask | cm[0][x];
    for(int i=0;i<depth;i++) nm |= cm[sel[i]][x];
    sel[depth]=x; rec(depth+1,x+1,nm);
  }
}
static int pc(uint64_t x){ return __builtin_popcountll(x); }
static int cmpm(const void*a,const void*b){ uint64_t x=*(uint64_t*)a,y=*(uint64_t*)b; int px=pc(x),py=pc(y); if(px!=py) return px-py; return x<y?-1:x>y; }
static uint64_t *E; static size_t NE;   // minimal hyperedges
// solver
static int order[64], pos[64];
static uint64_t *byLast[64]; static int nByLast[64];
static uint64_t A[3]; static int col[64]; static long long nsol, maxsol; static long long nodes;
static void report(void){
  nsol++;
  if(nsol<=20 || nsol<=maxsol){ printf("SOL "); for(int c=1;c<=m;c++) printf("%d", col[c-1]); printf("\n"); }
}
static int used;
static int stopflag;
static void solve(int i, int ncolused){
  if(stopflag) return;
  nodes++;
  if(i==m){ report(); if(maxsol>0 && nsol>=maxsol) stopflag=1; return; }
  int v=order[i]; uint64_t b=1ULL<<v;
  int lim = ncolused<3 ? ncolused+1 : 3; // colours introduced in order
  for(int c=0;c<lim;c++){
    A[c]|=b; col[v]=c;
    int ok=1;
    for(int j=0;j<nByLast[v];j++){ uint64_t T=byLast[v][j]; if(!(T&A[0]) || !(T&A[1]) || !(T&A[2])){ ok=0; break; } }
    if(ok) solve(i+1, c==ncolused? ncolused+1 : ncolused);
    A[c]&=~b;
  }
}
int main(int argc,char**argv){
  if(argc<4){ fprintf(stderr,"usage\n"); return 1; }
  FILE*f=fopen(argv[1],"r"); if(!f){perror("open");return 1;}
  if(fscanf(f,"%d",&n)!=1) return 1;
  for(int a=0;a<n;a++) for(int b=0;b<n;b++) if(fscanf(f,"%d",&mul[a][b])!=1){fprintf(stderr,"bad table\n");return 1;}
  fclose(f);
  K=atoi(argv[2]); const char*mode=argv[3]; maxsol = argc>4? atoll(argv[4]) : 0;
  for(int a=0;a<n;a++){ inv_[a]=-1; for(int b=0;b<n;b++) if(mul[a][b]==0){ inv_[a]=b; break; } }
  for(int a=0;a<n;a++) cls[a]=-1;
  m=0; for(int a=1;a<n;a++) if(cls[a]<0){ cls[a]=m; cls[inv_[a]]=m; rep[m]=a; m++; }
  if(m>64){ fprintf(stderr,"too many classes\n"); return 1; }
  for(int a=0;a<n;a++) for(int b=0;b<n;b++) if(a!=b){ int d=mul[inv_[a]][b]; cm[a][b]=1ULL<<cls[d]; }
  hcap=1<<16; hs=calloc(hcap,8); hcnt=0;
  rec(0,1,0);
  size_t N=0; uint64_t *L=malloc(hcnt*8); for(size_t i=0;i<hcap;i++) if(hs[i]) L[N++]=hs[i];
  free(hs);
  qsort(L,N,8,cmpm);
  // minimalize: list approach indexed by lowest set bit
  E=malloc(N*8); NE=0;
  uint64_t **lb=calloc(64,sizeof(uint64_t*)); size_t *nlb=calloc(64,sizeof(size_t)), *clb=calloc(64,sizeof(size_t));
  for(size_t i=0;i<N;i++){
    uint64_t T=L[i]; int red=0;
    for(uint64_t t=T; t && !red; t&=t-1){ int b=__builtin_ctzll(t); for(size_t j=0;j<nlb[b];j++){ if((lb[b][j]&~T)==0){ red=1; break; } } }
    if(!red){ E[NE++]=T; int b=__builtin_ctzll(T); if(nlb[b]==clb[b]){ clb[b]=clb[b]?2*clb[b]:64; lb[b]=realloc(lb[b],clb[b]*8);} lb[b][nlb[b]++]=T; }
  }
  int minsz=99, maxsz=0; for(size_t i=0;i<NE;i++){ int p=pc(E[i]); if(p<minsz)minsz=p; if(p>maxsz)maxsz=p; }
  fprintf(stderr,"n=%d K=%d classes=%d distinct=%zu minimal=%zu size %d..%d\n", n,K,m,N,NE,minsz,maxsz);
  if(!strcmp(mode,"edges")){
    printf("%d %d %d %zu\n", n, K, m, NE);
    for(int a=0;a<n;a++) printf("%d ", a==0?-1:cls[a]); printf("\n");
    for(size_t i=0;i<NE;i++) printf("%llx\n",(unsigned long long)E[i]);
    return 0;
  }
  // variable order: greedy, pick class covering most weight among small hyperedges
  double w[64]={0}; for(size_t i=0;i<NE;i++){ double ww=1.0/(1ULL<<pc(E[i])); for(uint64_t t=E[i];t;t&=t-1) w[__builtin_ctzll(t)]+=ww; }
  int done[64]={0}; uint64_t chosen=0;
  for(int i=0;i<m;i++){
    // prefer variable completing most hyperedges given chosen set, tie-break by weight
    int best=-1; double bs=-1;
    for(int v=0;v<m;v++) if(!done[v]){
      uint64_t c2=chosen|(1ULL<<v); double s=0;
      for(size_t j=0;j<NE;j++){ uint64_t T=E[j]; if((T>>v)&1){ int miss=pc(T&~c2); s+= 1.0/(1ULL<<(miss*2)); } }
      s+=1e-9*w[v];
      if(s>bs){bs=s;best=v;}
    }
    order[i]=best; pos[best]=i; done[best]=1; chosen|=1ULL<<best;
  }
  for(int v=0;v<m;v++){ byLast[v]=malloc(NE*8); nByLast[v]=0; }
  for(size_t j=0;j<NE;j++){ int last=-1,lv=-1; for(uint64_t t=E[j];t;t&=t-1){ int b=__builtin_ctzll(t); if(pos[b]>last){last=pos[b];lv=b;} } byLast[lv][nByLast[lv]++]=E[j]; }
  A[0]=A[1]=A[2]=0;
  solve(0,0);
  fprintf(stderr,"solutions(colour-perm reduced)=%lld nodes=%lld\n", nsol, nodes);
  printf("COUNT %lld\n", nsol);
  return 0;
}

