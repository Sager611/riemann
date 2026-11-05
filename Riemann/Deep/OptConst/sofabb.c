// sofabb.c -- float branch-and-bound for ambidextrous sofa bounds (exploration + certificate tree).
// Usage: sofabb cfg.txt out.bits
// cfg: line1: k W T margin maxdepth ; next k lines: cnum cden snum sden (cos=cnum/cden, sin=snum/sden)
// next k lines: a1 b1 a2 b2 (domain box per corridor, rationals given as decimal doubles exactly representable)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define MAXV 96
#define MAXK 8
typedef struct { double x[MAXV], y[MAXV]; int n; } Poly;
static int K; static double Wd, T, MARGIN; static int MAXDEPTH;
static double C[MAXK], S[MAXK];
static void clip(const Poly *in, Poly *out, double a, double b, double c){
  int n=in->n; out->n=0; if(n==0) return;
  for(int i=0;i<n;i++){
    int j=(i+1)%n; double px=in->x[i],py=in->y[i],qx=in->x[j],qy=in->y[j];
    double dp=a*px+b*py-c, dq=a*qx+b*qy-c;
    if(dp<=0){ out->x[out->n]=px; out->y[out->n]=py; out->n++; }
    if((dp<0&&dq>0)||(dp>0&&dq<0)){ double t=dp/(dp-dq); out->x[out->n]=px+t*(qx-px); out->y[out->n]=py+t*(qy-py); out->n++; }
  }
  if(out->n<3) out->n=0;
}
static double parea(const Poly *p){ double s=0; for(int i=0;i<p->n;i++){int j=(i+1)%p->n; s+=p->x[i]*p->y[j]-p->x[j]*p->y[i];} return fabs(s)*0.5; }
// half-plane in rotated coords: cp*p + cq*q <= r  with p = c x + s y, q = -s x + c y
static void hp_world(int j,double cp,double cq,double r,double *A,double *B,double *Cc){
  *A=cp*C[j]-cq*S[j]; *B=cp*S[j]+cq*C[j]; *Cc=r; }
static double box[2*MAXK][2]; // per corridor: u1 interval box[2j], u2 interval box[2j+1]
static double tot;
static void rec(int j, const Poly *poly, int sign){
  if(poly->n==0) return;
  if(j==K){ tot+=sign*parea(poly); return; }
  double a1=box[2*j][0], b1=box[2*j][1], a2=box[2*j+1][0], b2=box[2*j+1][1];
  double A,B,Cc; Poly t1,t2,t3,t4;
  // P~: p <= b1+1 ; -q <= -a2 ; q <= b2+1
  hp_world(j,1,0,b1+1,&A,&B,&Cc); clip(poly,&t1,A,B,Cc);
  hp_world(j,0,-1,-a2,&A,&B,&Cc); clip(&t1,&t2,A,B,Cc);
  hp_world(j,0,1,b2+1,&A,&B,&Cc); clip(&t2,&t3,A,B,Cc);
  rec(j+1,&t3,sign);
  // Q~: -p <= -a1 ; p <= b1+1 ; q <= b2
  hp_world(j,-1,0,-a1,&A,&B,&Cc); clip(poly,&t1,A,B,Cc);
  hp_world(j,1,0,b1+1,&A,&B,&Cc); clip(&t1,&t2,A,B,Cc);
  hp_world(j,0,1,b2,&A,&B,&Cc); clip(&t2,&t3,A,B,Cc);
  rec(j+1,&t3,sign);
  // PQ: -p<=-a1, p<=b1+1, -q<=-a2, q<=b2   (subtract)
  hp_world(j,-1,0,-a1,&A,&B,&Cc); clip(poly,&t1,A,B,Cc);
  hp_world(j,1,0,b1+1,&A,&B,&Cc); clip(&t1,&t2,A,B,Cc);
  hp_world(j,0,-1,-a2,&A,&B,&Cc); clip(&t2,&t3,A,B,Cc);
  hp_world(j,0,1,b2,&A,&B,&Cc); clip(&t3,&t4,A,B,Cc);
  rec(j+1,&t4,-sign);
}
static double UB(void){
  Poly k0; k0.n=4; k0.x[0]=0;k0.y[0]=0;k0.x[1]=Wd;k0.y[1]=0;k0.x[2]=Wd;k0.y[2]=1;k0.x[3]=0;k0.y[3]=1;
  tot=0; rec(0,&k0,1); return tot; }
static unsigned char *bits; static long nbits=0, cap=0;
static void putbit(int b){ if(nbits>=cap*8){ long nc=cap?cap*2:(1<<20); bits=realloc(bits,nc); memset(bits+cap,0,nc-cap); cap=nc;} if(b) bits[nbits>>3]|=(unsigned char)(1<<(nbits&7)); nbits++; }
static unsigned char *dims; static long ndims=0, dcap=0; static int SPLITBEST=0;
static void putdim(int d){ if(ndims>=dcap){ long nc=dcap?dcap*2:(1<<20); dims=realloc(dims,nc); dcap=nc;} dims[ndims++]=(unsigned char)d; }
static long nodes=0, leaves=0; static int maxd=0; static double maxleaf=-1, maxfail=-1; static long fails=0;
static time_t t0; static double best_center=0;
static void dfs(int depth){
  nodes++; if(depth>maxd) maxd=depth;
  double ub=UB();
  if(ub<=T-MARGIN){ putbit(0); leaves++; if(ub>maxleaf) maxleaf=ub; return; }
  if(depth>=MAXDEPTH){ putbit(0); leaves++; fails++; if(ub>maxfail) maxfail=ub; return; }
  putbit(1);
  int bi=0; double bw=-1;
  if(!SPLITBEST){ for(int i=0;i<2*K;i++){ double w=box[i][1]-box[i][0]; if(w>bw){bw=w;bi=i;} } }
  else { double bestv=1e300; for(int i=0;i<2*K;i++){ double lo=box[i][0],hi=box[i][1],mid=0.5*(lo+hi); if(hi-lo<1e-12) continue;
      box[i][1]=mid; double u1=UB(); box[i][1]=hi; box[i][0]=mid; double u2=UB(); box[i][0]=lo;
      double v=(u1>u2?u1:u2)+1e-3*(u1<u2?u1:u2); if(v<bestv){bestv=v;bi=i;} }
    putdim(bi); }
  double lo=box[bi][0], hi=box[bi][1], mid=0.5*(lo+hi);
  box[bi][1]=mid; dfs(depth+1); box[bi][1]=hi;
  box[bi][0]=mid; dfs(depth+1); box[bi][0]=lo;
  if((nodes & ((1<<22)-1))==0){ fprintf(stderr,"nodes %ld leaves %ld fails %ld maxd %d t %ld\n",nodes,leaves,fails,maxd,(long)(time(NULL)-t0)); }
}
int main(int argc,char**argv){
  if(argc<3){fprintf(stderr,"usage\n");return 1;}
  FILE*f=fopen(argv[1],"r"); if(!f){perror("cfg");return 1;}
  if(fscanf(f,"%d %lf %lf %lf %d",&K,&Wd,&T,&MARGIN,&MAXDEPTH)!=5) return 2;
  for(int j=0;j<K;j++){ double cn,cd,sn,sd; if(fscanf(f,"%lf %lf %lf %lf",&cn,&cd,&sn,&sd)!=4) return 3; C[j]=cn/cd; S[j]=sn/sd; }
  for(int j=0;j<K;j++){ if(fscanf(f,"%lf %lf %lf %lf",&box[2*j][0],&box[2*j][1],&box[2*j+1][0],&box[2*j+1][1])!=4) return 4; }
  fclose(f);
  t0=time(NULL); SPLITBEST=getenv("SPLIT")!=NULL;
  double r=UB(); fprintf(stderr,"root UB %.9f\n",r);
  dfs(0);
  FILE*o=fopen(argv[2],"wb"); fwrite(&nbits,sizeof(long),1,o); fwrite(bits,1,(nbits+7)/8,o); fwrite(&ndims,sizeof(long),1,o); if(ndims) fwrite(dims,1,ndims,o); fclose(o);
  printf("nodes %ld leaves %ld fails %ld maxdepth %d maxleafUB %.12f maxfailUB %.12f time %lds\n",nodes,leaves,fails,maxd,maxleaf,maxfail,(long)(time(NULL)-t0));
  return fails?5:0;
}
