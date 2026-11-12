// fam.c (rc5audit, own code): the family {(D,chi): D<0 fundamental, X/2<|D|<=X, chi in Cl(D)^}.
// Classes = reduced forms; v_f(n) = r_f(n)/2 = #ideals of norm n in the class of f (|D|>4).
// mode 0 (det): read specs "N sigma V tp kill deriv" from stdin; b supported on admissible n (odd squarefree,
//   all primes >= 19).  For each spec prints  NB=sum|b|^2, LHS=sum_D h sum_f |S_f|^2 (= sum_D sum_chi |sum b lambda_chi|^2),
//   TRIV=sum_D |sum_f S_f|^2 (trivial chars), GEN=sum_D sum_{genus chi}|sum_f chi(f)S_f|^2 (incl. trivial).
// mode 1 (exp): export sparse v_f(n) for admissible n in [nlo,nhi) to files prefix.{trip,cls,gen}.
// mode 2 (dual): same specs on stdin; prints Q'(b) = sum over all orders |Delta| in [X/4,2X] of Phi(|Delta|/X) Q'_Delta(b)
//   via the Lemma P blocks (w_e/2) sum_C |S_C^{(Delta/e^2)}(b_e)|^2, and Q'/(X NB).
// build: cc -O2 -o fam fam.c -lm.  usage: fam X mode [nlo nhi prefix]
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

static int cmpl(const void*x,const void*y){ long a=*(const long*)x,b=*(const long*)y; return (a>b)-(a<b); }
static long gcdl(long a,long b){ if(a<0)a=-a; if(b<0)b=-b; while(b){long t=a%b;a=b;b=t;} return a; }
static int issqf(long m){ if(m<0)m=-m; for(long p=2;p*p<=m;p++){ if(m%(p*p)==0) return 0; } return 1; }
static int isfund(long D){ // D<0
  long m=-D; if(((D%4)+4)%4==1) return issqf(m);
  if(m%4) return 0; long q=D/4; long r=((q%4)+4)%4; if(r!=2 && r!=3) return 0; return issqf(q); }
static int legendre(long a,long p){ a%=p; if(a<0)a+=p; if(a==0) return 0; long r=1,e=(p-1)/2,b=a; while(e){ if(e&1) r=r*b%p; b=b*b%p; e>>=1; } return r==1?1:-1; }
// value of assigned character of prime discriminant qs (=-4,8,-8, or p* odd) at m coprime to it
static int chiq(long qs,long m){
  if(qs==-4){ long r=((m%4)+4)%4; return r==1?1:-1; }
  if(qs==8){ long r=((m%8)+8)%8; return (r==1||r==7)?1:-1; }
  if(qs==-8){ long r=((m%8)+8)%8; return (r==1||r==3)?1:-1; }
  long p=qs<0?-qs:qs; return legendre(m,p);
}
static double sstep(double t){ if(t<=0) return 0; if(t>=1) return 1; double a=exp(-1.0/t), b=exp(-1.0/(1.0-t)); return a/(a+b); }
static double phi2(double v){ if(v<=0.25||v>=4.0) return 0.0; double l=log2(v); if(l<-1) return sstep(l+2); if(l>1) return sstep(2-l); return 1.0; }

int main(int argc,char**argv){
  if(argc<3){ fprintf(stderr,"usage\n"); return 1; }
  alarm(1700);
  long X=atol(argv[1]); int mode=atoi(argv[2]);
  long nlo=0,nhi=0; const char *prefix=NULL;
  // specs
  int K=0; long SN[64]; double Ssig[64], SV[64], Stp[64]; int Skill[64], Sder[64];
  long nmax=0;
  if(mode==0||mode==2){
    while(K<64 && scanf("%ld %lf %lf %lf %d %d",&SN[K],&Ssig[K],&SV[K],&Stp[K],&Skill[K],&Sder[K])==6){ if(4*SN[K]>nmax) nmax=4*SN[K]; K++; }
  } else { nlo=atol(argv[3]); nhi=atol(argv[4]); prefix=argv[5]; nmax=nhi-1; }
  // admissible n
  char *adm=malloc(nmax+1); memset(adm,1,nmax+1); adm[0]=0;
  for(long p=2;p*p<=nmax;p++) for(long m=p*p;m<=nmax;m+=p*p) adm[m]=0;
  { long sp[7]={2,3,5,7,11,13,17}; for(int i=0;i<7;i++) for(long m=sp[i];m<=nmax;m+=sp[i]) adm[m]=0; }
  // slots: each spec -> 1 real slot (tp==0) or 2 slots (re,im)
  int slot0[64], nsl[64], NS=0; for(int k=0;k<K;k++){ slot0[k]=NS; nsl[k]=(Stp[k]!=0.0)?2:1; NS+=nsl[k]; }
  double *B=NULL; double NB[64]; memset(NB,0,sizeof NB);
  if(mode==0||mode==2){
    B=calloc((size_t)(nmax+1)*NS,sizeof(double));
    // a_V for each distinct V
    long Vmax=0; for(int k=0;k<K;k++) if((long)SV[k]>Vmax) Vmax=(long)SV[k];
    int *mu=malloc((Vmax+2)*sizeof(int)); for(long k=0;k<=Vmax+1;k++) mu[k]=1;
    { char *comp=calloc(Vmax+2,1); for(long p=2;p<=Vmax;p++) if(!comp[p]){ for(long m=p;m<=Vmax;m+=p){ if(m>p) comp[m]=1; mu[m]=-mu[m]; } for(long m=p*p;m<=Vmax;m+=p*p) mu[m]=0; } free(comp); }
    double *aV=malloc((nmax+1)*sizeof(double));
    for(int k=0;k<K;k++){
      long V=(long)SV[k]; memset(aV,0,(nmax+1)*sizeof(double));
      for(long d=1;d<=V;d++) if(mu[d]) for(long m=d;m<=nmax;m+=d) aV[m]+=mu[d];
      double sig=Ssig[k], tp=Stp[k], N=(double)SN[k], Y=(double)X, l2=pow(4.0,1.0-sig);
      for(long n=1;n<=nmax;n++){ if(!adm[n]) continue; double v=n/N;
        double w1=exp(-n/Y)*phi2(v), w2=l2*exp(-4.0*n/Y)*phi2(4*v); double ln=log((double)n);
        double pre=aV[n]*pow((double)n,-sig);
        if(tp==0.0){
          double val;
          if(Sder[k]==0) val=pre*(w1-(Skill[k]?w2:0.0));
          else val=pre*(-ln*(w1-(Skill[k]?w2:0.0)) + (Skill[k]? 2*log(2.0)*w2:0.0));
          B[(size_t)n*NS+slot0[k]]=val; NB[k]+=val*val;
        } else { // complex: pre*n^{-i tp}[w1 - 4^{-i tp} w2]
          double c1=cos(tp*ln), s1=-sin(tp*ln); double c2=cos(2*tp*log(2.0)), s2=-sin(2*tp*log(2.0));
          double re_in=w1-(Skill[k]? w2*c2:0.0), im_in=-(Skill[k]? w2*s2:0.0);
          double re=pre*(c1*re_in - s1*im_in), im=pre*(c1*im_in + s1*re_in);
          if(Sder[k]){ fprintf(stderr,"deriv with tp!=0 unsupported\n"); return 2; }
          B[(size_t)n*NS+slot0[k]]=re; B[(size_t)n*NS+slot0[k]+1]=im; NB[k]+=re*re+im*im;
        }
      }
    }
    free(aV); free(mu);
  }
  double LHS[128], TRIV[128], GEN[128]; memset(LHS,0,sizeof LHS); memset(TRIV,0,sizeof TRIV); memset(GEN,0,sizeof GEN);
  long nD=0, famsize=0;
  int maxf=20000; long *FA=malloc(maxf*sizeof(long)),*FB=malloc(maxf*sizeof(long)),*FC=malloc(maxf*sizeof(long));
  double *S=malloc((size_t)maxf*NS*sizeof(double)); signed char *chi=malloc((size_t)maxf*16);
  FILE *ftrip=NULL,*fcls=NULL,*fgen=NULL; long gclass=0;
  long *nlist=NULL; if(mode==1) nlist=malloc(4000000*sizeof(long));
  if(mode==1){ char fn[512]; sprintf(fn,"%s.trip",prefix); ftrip=fopen(fn,"wb"); sprintf(fn,"%s.cls",prefix); fcls=fopen(fn,"w"); sprintf(fn,"%s.gen",prefix); fgen=fopen(fn,"wb"); }
  if(mode==2){ // dual: Q' = sum_Delta Phi(|Delta|/X) sum_{e|f sqfree, admissible} (w_e/2) sum_C |S_C^{(Delta_e)}(b_e)|^2
    double QP[128]; memset(QP,0,sizeof QP); long nDelta=0, maxh=0;
    for(long Da=X/4; Da<=2*X; Da++){
      long D=-Da; long r4=((D%4)+4)%4; if(r4!=0 && r4!=1) continue;
      double u=(double)Da/X, Ph; if(u<=0.25||u>=2.0) continue; if(u<0.5) Ph=sstep((u-0.25)/0.25); else if(u<=1.0) Ph=1.0; else Ph=sstep(2.0-u);
      if(Ph<=0) continue; nDelta++;
      long f=1; { long Dc=D; for(long q=2;q*q<=-Dc;q++){ while(Dc%(q*q)==0){ long r=(((Dc/(q*q))%4)+4)%4; if(r!=0&&r!=1) break; Dc/=q*q; f*=q; } } }
      long fp[16]; int nfp=0; { long ff=f; for(long q=2;q<=ff;q++){ if(ff%q==0){ if(q>=19) fp[nfp++]=q; while(ff%q==0) ff/=q; } } }
      for(int mask=0; mask<(1<<nfp); mask++){
        long e=1; for(int i=0;i<nfp;i++) if(mask>>i&1) e*=fp[i];
        long De=D/(e*e), Dea=-De; long mmax=nmax/e; if(mmax<1) continue;
        int we=(De==-3)?6:((De==-4)?4:2);
        int nf=0;
        for(long a=1; 3*a*a<=Dea; a++) for(long b=-a+1;b<=a;b++){ if(((b-De)%2+2)%2) continue; long num=b*b+Dea; if(num%(4*a)) continue; long c=num/(4*a); if(c<a) continue; if(b<0 && a==c) continue; if(gcdl(gcdl(a,b),c)!=1) continue; if(nf>=maxf){fprintf(stderr,"maxf\n");return 4;} FA[nf]=a;FB[nf]=b;FC[nf]=c;nf++; }
        if(e==1 && nf>maxh) maxh=nf;
        memset(S,0,(size_t)nf*NS*sizeof(double));
        for(int fi=0;fi<nf;fi++){ long a=FA[fi],b=FB[fi],c=FC[fi]; double *Sf=S+(size_t)fi*NS;
          long ymax=(long)floor(sqrt(4.0*a*mmax/(double)Dea))+1;
          for(long y=-ymax;y<=ymax;y++){ double disc=4.0*a*mmax - (double)Dea*y*y; if(disc<0) continue; double sq=sqrt(disc);
            long x0=(long)ceil((-b*y - sq)/(2.0*a))-1, x1=(long)floor((-b*y+sq)/(2.0*a))+1;
            for(long x=x0;x<=x1;x++){ long m=a*x*x+b*x*y+c*y*y; if(m<1||m>mmax) continue; long n=e*m; if(!adm[n]) continue; if(gcdl(m,f)!=1) continue;
              const double *bn=B+(size_t)n*NS; for(int s=0;s<NS;s++) Sf[s]+=bn[s]/we; } } }
        for(int s=0;s<NS;s++){ double sumsq=0; for(int fi=0;fi<nf;fi++){ double v=S[(size_t)fi*NS+s]; sumsq+=v*v; } QP[s]+=Ph*0.5*we*sumsq; }
      }
    }
    printf("# dual X=%ld nDelta=%ld maxh(e=1 orders)=%ld\n",X,nDelta,maxh);
    for(int k=0;k<K;k++){ double q=0; for(int s=slot0[k]; s<slot0[k]+nsl[k]; s++) q+=QP[s];
      printf("%ld %.4f %.1f %.2f %d %d  NB=%.6e QP=%.6e QP/(X NB)=%.4f\n",SN[k],Ssig[k],SV[k],Stp[k],Skill[k],Sder[k],NB[k],q,q/((double)X*NB[k])); }
    return 0;
  }
  long ngencheck_bad=0;
  for(long Da=X/2+1; Da<=X; Da++){
    long D=-Da; if(!isfund(D)) continue;
    // prime discriminants
    long qs[16]; int t=0; { long m=Da; while(m%2==0) m/=2; long odd=1; long mm=m; for(long p=3;p*p<=mm;p+=2){ if(mm%p==0){ mm/=p; long ps=(p%4==1)?p:-p; qs[t++]=ps; odd*=ps; } } if(mm>1){ long ps=(mm%4==1)?mm:-mm; qs[t++]=ps; odd*=ps; }
      if(D%odd!=0){ fprintf(stderr,"fact err\n"); return 3; } long two=D/odd; if(two!=1){ if(two!=-4 && two!=8 && two!=-8){ fprintf(stderr,"2part err D=%ld two=%ld\n",D,two); return 3; } qs[t++]=two; } }
    // reduced forms
    int nf=0;
    for(long a=1; 3*a*a<=Da; a++) for(long b=-a+1;b<=a;b++){ if(((b-D)%2+2)%2) continue; long num=b*b+Da; if(num%(4*a)) continue; long c=num/(4*a); if(c<a) continue; if(b<0 && a==c) continue; if(gcdl(gcdl(a,b),c)!=1) continue; if(nf>=maxf){fprintf(stderr,"maxf\n");return 4;} FA[nf]=a;FB[nf]=b;FC[nf]=c;nf++; }
    long h=nf; nD++; famsize+=h;
    // assigned characters per form
    for(int f=0;f<nf;f++) for(int i=0;i<t;i++){ long q=qs[i]<0?-qs[i]:qs[i]; if(q==4||q==8) q=2; int got=0;
      for(long x=-3;x<=3 && !got;x++) for(long y=-3;y<=3 && !got;y++){ if(x==0&&y==0) continue; long m=FA[f]*x*x+FB[f]*x*y+FC[f]*y*y; if(gcdl(m,q)==1){ chi[f*16+i]=(signed char)chiq(qs[i],m); got=1; } }
      if(!got){ fprintf(stderr,"no m\n"); return 5; } }
    int ng=1<<(t-1); // genus chars: subsets of first t-1 prime discs
    { // check orthogonality: each nontrivial genus char sums to 0 over classes
      for(int s=1;s<ng;s++){ long sum=0; for(int f=0;f<nf;f++){ int v=1; for(int i=0;i<t-1;i++) if(s>>i&1) v*=chi[f*16+i]; sum+=v; } if(sum!=0) ngencheck_bad++; } }
    if(mode==0){
      memset(S,0,(size_t)nf*NS*sizeof(double));
      for(int f=0;f<nf;f++){ long a=FA[f],b=FB[f],c=FC[f]; double *Sf=S+(size_t)f*NS;
        long ymax=(long)floor(sqrt(4.0*a*nmax/(double)Da))+1;
        for(long y=-ymax;y<=ymax;y++){ double disc=4.0*a*nmax - (double)Da*y*y; if(disc<0) continue; double sq=sqrt(disc);
          long x0=(long)ceil((-b*y - sq)/(2.0*a))-1, x1=(long)floor((-b*y+sq)/(2.0*a))+1;
          for(long x=x0;x<=x1;x++){ long n=a*x*x+b*x*y+c*y*y; if(n<1||n>nmax||!adm[n]) continue; const double *bn=B+(size_t)n*NS; for(int s=0;s<NS;s++) Sf[s]+=0.5*bn[s]; } } }
      for(int s=0;s<NS;s++){ double sumsq=0, tot=0; for(int f=0;f<nf;f++){ double v=S[(size_t)f*NS+s]; sumsq+=v*v; tot+=v; } LHS[s]+=h*sumsq; TRIV[s]+=tot*tot;
        for(int g=0; g<ng; g++){ double acc=0; for(int f=0;f<nf;f++){ int v=1; for(int i=0;i<t-1;i++) if(g>>i&1) v*=chi[f*16+i]; acc+=v*S[(size_t)f*NS+s]; } GEN[s]+=acc*acc; } }
    } else {
      fprintf(fcls,"%ld %ld %d\n",D,h,ng);
      for(int f=0;f<nf;f++){ long a=FA[f],b=FB[f],c=FC[f]; long cnt=0;
        long ymax=(long)floor(sqrt(4.0*a*nmax/(double)Da))+1;
        for(long y=-ymax;y<=ymax;y++){ double disc=4.0*a*nmax - (double)Da*y*y; if(disc<0) continue; double sq=sqrt(disc);
          long x0=(long)ceil((-b*y - sq)/(2.0*a))-1, x1=(long)floor((-b*y+sq)/(2.0*a))+1;
          for(long x=x0;x<=x1;x++){ long n=a*x*x+b*x*y+c*y*y; if(n<nlo||n>=nhi||!adm[n]) continue; if(cnt<4000000) nlist[cnt++]=n; } }
        // merge duplicates; write (class, n, reps) with v = reps/2
        qsort(nlist,cnt,sizeof(long),cmpl);
        for(long i=0;i<cnt;){ long j=i; while(j<cnt && nlist[j]==nlist[i]) j++; int g=(int)(gclass+f); int nn=(int)nlist[i]; int rep=(int)(j-i); fwrite(&g,4,1,ftrip); fwrite(&nn,4,1,ftrip); fwrite(&rep,4,1,ftrip); i=j; }
        // genus values for this class: ng ints
        for(int g=0; g<ng; g++){ int v=1; for(int i=0;i<t-1;i++) if(g>>i&1) v*=chi[f*16+i]; signed char vv=(signed char)v; fwrite(&vv,1,1,fgen); }
      }
      gclass+=nf;
    }
  }
  if(mode==0){
    printf("# X=%ld nD=%ld famsize=%ld genus_orth_fail=%ld\n",X,nD,famsize,ngencheck_bad);
    for(int k=0;k<K;k++){ double L=0,T=0,G=0; for(int s=slot0[k]; s<slot0[k]+nsl[k]; s++){ L+=LHS[s]; T+=TRIV[s]; G+=GEN[s]; }
      printf("%ld %.4f %.1f %.2f %d %d  NB=%.6e LHS=%.6e TRIV=%.6e GEN=%.6e\n",SN[k],Ssig[k],SV[k],Stp[k],Skill[k],Sder[k],NB[k],L,T,G); }
  } else { fclose(ftrip); fclose(fcls); fclose(fgen); printf("# X=%ld nD=%ld famsize=%ld classes=%ld genus_orth_fail=%ld\n",X,nD,famsize,gclass,ngencheck_bad); }
  return 0;
}
