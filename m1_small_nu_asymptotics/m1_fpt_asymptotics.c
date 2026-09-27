/*
 * m1_fpt_asymptotics.c — first-passage moments of the reduced recovery process at high resolution (OpenMP),
 *                        Monte-Carlo checks and the FPT density, for the small-ν expansion.
 *
 * Reduced process (rescaled time s' = s/ε_eff):  dH = (1−H) ds' + ν sqrt(H(1−H)) dW,  H(0) = H0, absorbed at H*.
 * Dynkin/Siegert moments with the substitution z = y w^{1/k}, k = 2/ν²:
 *   T1(x) = ∫_x^{H*} G(y) dy,  G(y) = ∫_0^1 dw/(1 − y w^{1/k});   T2 = 2∫∫ T1/(1−z);   T3 = 3∫∫ T2/(1−z)
 * (same quadrature as m1_first_passage_theory/m1_fpt_jacobi.c, but: the y-integrals use Simpson's rule, the
 *  inner w-integral a composite Gauss–Legendre rule, both grids are parameters, and the y-loop is parallel).
 * The cubic cancellation in the third cumulant makes the skewness at small ν sensitive to the grid; the run
 * script doubles ny, nw to show convergence.
 *
 * Small-ν expansion (derived symbolically in derive_small_nu.py; a = ν²/2, z0 = 1/(1−H0), z_s = 1/(1−H*)):
 *   mean = T0 − a I1 + a² (z_s² − z0² − 3(z_s − z0) + T0) + O(a³),  T0 = ln(z_s / z0),  I1 = z_s − z0 − T0
 *   var  = 2a I1 + a² (−4(z_s² − z0²) + 14(z_s − z0) − 6 T0) + O(a³)
 *   κ3   = 6a² (z_s² − z0² − 4(z_s − z0) + 2 T0) + O(a³)
 * so sd ≈ ν sqrt(I1), mean ≈ T0 − var/2, skew ≈ (3/2) ν K3 / I1^{3/2}, skew/CV → (3/2) K3 T0 / I1².
 *
 * Usage:
 *   ./m1_fpt_asymptotics sweep=nu lo=0.01 hi=3 n=30 log=1 H0=0.064 Hs=0.408 [ny=8000 nw=4000]
 *        -> nu,H0,Hs,mean,sd,CV,skew,mean_det, and the asymptotic values mean_as,sd_as,skew_as (leading orders)
 *   ./m1_fpt_asymptotics mc=200000 nu=0.3 H0=.. Hs=.. [dt=1e-4] [hist=1 nbins=80]
 *        -> moments (and, with hist=1, the histogram of the FPT: t_center, density)
 *   ./m1_fpt_asymptotics ratio=1 H0=.. Hslo=0.15 Hshi=0.85 n=29     -> the limit ratio skew/CV vs H*
 * Threads: OpenMP (y-grid loop; MC trajectories with one xoshiro256** state per thread).
 * Build: see Makefile
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static int NY=8001, NW=4000;      /* NY odd for Simpson */

/* composite Gauss–Legendre (4-point) on NW panels of [0,1] */
static const double GX[4]={-0.8611363115940526,-0.3399810435848563,0.3399810435848563,0.8611363115940526};
static const double GW[4]={0.3478548451374538,0.6521451548625461,0.6521451548625461,0.3478548451374538};
static double inner(double y, double k, const double *f, const double *yg, int ny){
    double sum=0; int j,q; double hw=1.0/NW;
    for(j=0;j<NW;j++) for(q=0;q<4;q++){
        double w=(j+0.5*(1+GX[q]))*hw, z=y*pow(w,1.0/k), fz=1.0;
        if(f){ double p=z/yg[ny-1]*(ny-1); int i=(int)p; if(i>=ny-1)i=ny-2; double t=p-i; fz=f[i]*(1-t)+f[i+1]*t; }
        sum+=0.5*GW[q]*hw*fz/(1.0-z);
    }
    return sum;
}
/* cumulative integral from the top: T[i] = ∫_{y_i}^{Hs} g ; Simpson on even indices, half-panel 3-point rule on odd */
static void cumint(const double *g, const double *yg, int ny, double *T){
    int i; T[ny-1]=0; double h=yg[1]-yg[0];
    for(i=ny-3;i>=0;i-=2) T[i]=T[i+2]+h/3.0*(g[i]+4*g[i+1]+g[i+2]);
    for(i=ny-2;i>=1;i-=2) T[i]=T[i+1]+h/12.0*(5*g[i]+8*g[i+1]-g[i+2]);
}
static void moments(double nu, double H0, double Hs, double *mean, double *sd, double *skew, double *mean_det){
    double k=2.0/(nu*nu); int ny=NY, i;
    double *yg=malloc(sizeof(double)*ny), *T1=malloc(sizeof(double)*ny), *T2=malloc(sizeof(double)*ny), *T3=malloc(sizeof(double)*ny), *g=malloc(sizeof(double)*ny);
    for(i=0;i<ny;i++) yg[i]=Hs*i/(double)(ny-1);
#pragma omp parallel for schedule(static)
    for(i=0;i<ny;i++) g[i]=inner(yg[i],k,NULL,yg,ny);
    cumint(g,yg,ny,T1);
#pragma omp parallel for schedule(static)
    for(i=0;i<ny;i++) g[i]=2.0*inner(yg[i],k,T1,yg,ny);
    cumint(g,yg,ny,T2);
#pragma omp parallel for schedule(static)
    for(i=0;i<ny;i++) g[i]=3.0*inner(yg[i],k,T2,yg,ny);
    cumint(g,yg,ny,T3);
    double p=H0/Hs*(ny-1); int j=(int)p; if(j>=ny-1)j=ny-2; double t=p-j;
    double m1=T1[j]*(1-t)+T1[j+1]*t, m2=T2[j]*(1-t)+T2[j+1]*t, m3=T3[j]*(1-t)+T3[j+1]*t;
    double var=m2-m1*m1; *mean=m1; *sd=var>0?sqrt(var):0;
    double mu3=m3-3*m1*m2+2*m1*m1*m1; *skew=(var>0)?mu3/pow(var,1.5):0;
    *mean_det=log((1-H0)/(1-Hs));
    free(yg);free(T1);free(T2);free(T3);free(g);
}
/* leading-order asymptotics */
static void asym(double nu, double H0, double Hs, double *mean, double *sd, double *skew, double *ratio){
    double z0=1/(1-H0), zs=1/(1-Hs), T0=log(zs/z0), I1=zs-z0-T0, a=0.5*nu*nu;
    double K3=zs*zs-z0*z0-4*(zs-z0)+2*T0;
    *mean=T0-a*I1+a*a*(zs*zs-z0*z0-3*(zs-z0)+T0);
    double var=2*a*I1+a*a*(-4*(zs*zs-z0*z0)+14*(zs-z0)-6*T0); *sd=var>0?sqrt(var):0;
    *skew=1.5*nu*K3/pow(I1,1.5); *ratio=1.5*K3*T0/(I1*I1);
}

/* Monte-Carlo with per-thread RNG */
typedef struct { uint64_t s[4]; } Rng;
static inline uint64_t rotl(uint64_t x,int k){return (x<<k)|(x>>(64-k));}
static inline uint64_t next64(Rng *r){uint64_t *s=r->s; uint64_t res=rotl(s[1]*5,7)*9,t=s[1]<<17;s[2]^=s[0];s[3]^=s[1];s[1]^=s[2];s[0]^=s[3];s[2]^=t;s[3]=rotl(s[3],45);return res;}
static void seed_rng(Rng *r, uint64_t s){int i;for(i=0;i<4;i++){s+=0x9E3779B97F4A7C15ULL;uint64_t z=s;z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL;z=(z^(z>>27))*0x94D049BB133111EBULL;r->s[i]=z^(z>>31);}}
static inline double unif(Rng *r){return (next64(r)>>11)*(1.0/9007199254740992.0);}
static inline double gauss(Rng *r){double u1=unif(r),u2=unif(r);if(u1<1e-300)u1=1e-300;return sqrt(-2*log(u1))*cos(2*M_PI*u2);}
static void montecarlo(long N, double nu, double H0, double Hs, double dt, int hist, int nbins){
    double s1=0,s2=0,s3=0; long n; double *ts=hist?malloc(sizeof(double)*N):NULL;
#pragma omp parallel
    {
        Rng r; int tid=0;
#ifdef _OPENMP
        tid=omp_get_thread_num();
#endif
        seed_rng(&r,12345ULL+1000003ULL*tid);
#pragma omp for reduction(+:s1,s2,s3) schedule(static)
        for(n=0;n<N;n++){
            double h=H0,t=0;
            while(h<Hs){ double hh=fmin(fmax(h,0.0),1.0); h+=(1-h)*dt+nu*sqrt(hh*(1-hh))*sqrt(dt)*gauss(&r); if(h<0)h=-h; t+=dt; }
            s1+=t; s2+=t*t; s3+=t*t*t; if(ts) ts[n]=t;
        }
    }
    double m=s1/N, v=s2/N-m*m, mu3=s3/N-3*m*s2/N+2*m*m*m;
    printf("kind,N,mean,sd,CV,skew\nmc,%ld,%.6g,%.6g,%.6g,%.6g\n",N,m,sqrt(v),sqrt(v)/m,mu3/pow(v,1.5));
    if(ts){ double tmax=m+6*sqrt(v); long *cnt=calloc(nbins,sizeof(long)); for(n=0;n<N;n++){ int b=(int)(ts[n]/tmax*nbins); if(b>=0&&b<nbins) cnt[b]++; }
        printf("t,density\n"); int b; for(b=0;b<nbins;b++) printf("%.6g,%.6g\n",(b+0.5)*tmax/nbins,cnt[b]/(double)N/(tmax/nbins)); free(cnt); free(ts); }
}

int main(int argc,char**argv){
    double nu=0.3,H0=0.064,Hs=0.408,lo=0,hi=0,dt=1e-4,Hslo=0.15,Hshi=0.85; long mc=0; int n=0,islog=0,i,hist=0,nbins=80,ratio=0,threads=0; char sweep[8]="";
    for(i=1;i<argc;i++){
        char key[64],val[64]; if(sscanf(argv[i],"%63[^=]=%63s",key,val)!=2){fprintf(stderr,"bad arg %s\n",argv[i]);return 1;}
        if(!strcmp(key,"nu"))nu=atof(val); else if(!strcmp(key,"H0"))H0=atof(val); else if(!strcmp(key,"Hs"))Hs=atof(val);
        else if(!strcmp(key,"sweep"))strncpy(sweep,val,7); else if(!strcmp(key,"lo"))lo=atof(val); else if(!strcmp(key,"hi"))hi=atof(val);
        else if(!strcmp(key,"n"))n=atoi(val); else if(!strcmp(key,"log"))islog=atoi(val); else if(!strcmp(key,"mc"))mc=atol(val);
        else if(!strcmp(key,"dt"))dt=atof(val); else if(!strcmp(key,"ny")){NY=atoi(val); if(NY%2==0)NY++;} else if(!strcmp(key,"nw"))NW=atoi(val);
        else if(!strcmp(key,"hist"))hist=atoi(val); else if(!strcmp(key,"nbins"))nbins=atoi(val); else if(!strcmp(key,"ratio"))ratio=atoi(val);
        else if(!strcmp(key,"Hslo"))Hslo=atof(val); else if(!strcmp(key,"Hshi"))Hshi=atof(val); else if(!strcmp(key,"threads"))threads=atoi(val);
        else {fprintf(stderr,"unknown key %s\n",key);return 1;}
    }
#ifdef _OPENMP
    if(threads>0) omp_set_num_threads(threads);
    fprintf(stderr,"OpenMP: %d threads\n",omp_get_max_threads());
#endif
    if(mc>0){ montecarlo(mc,nu,H0,Hs,dt,hist,nbins); return 0; }
    if(ratio){ printf("H0,Hs,T0,I1,K3,sd_over_nu,skew_over_nu,skew_over_CV_limit\n"); int nn=n>1?n:29;
        for(i=0;i<nn;i++){ double hs=Hslo+(Hshi-Hslo)*i/(double)(nn-1), m,s,sk,r; asym(1.0,H0,hs,&m,&s,&sk,&r);
            double z0=1/(1-H0), zs=1/(1-hs), T0=log(zs/z0), I1=zs-z0-T0, K3=zs*zs-z0*z0-4*(zs-z0)+2*T0;
            printf("%.4f,%.4f,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n",H0,hs,T0,I1,K3,sqrt(I1),1.5*K3/pow(I1,1.5),r); } return 0; }
    printf("nu,H0,Hs,mean,sd,CV,skew,mean_det,mean_as,sd_as,skew_as,ratio_as\n");
    int nv=sweep[0]?n:1;
    for(i=0;i<nv;i++){
        if(sweep[0]){ double v=islog?lo*pow(hi/lo,(double)i/(nv-1)):lo+(hi-lo)*i/(double)(nv>1?nv-1:1);
            if(!strcmp(sweep,"nu"))nu=v; else if(!strcmp(sweep,"Hs"))Hs=v; else if(!strcmp(sweep,"H0"))H0=v; else {fprintf(stderr,"sweep must be nu, Hs or H0\n");return 1;} }
        double m,s,sk,md,ma,sa,ska,ra; moments(nu,H0,Hs,&m,&s,&sk,&md); asym(nu,H0,Hs,&ma,&sa,&ska,&ra);
        printf("%.6g,%.4f,%.4f,%.8g,%.8g,%.6g,%.6g,%.8g,%.8g,%.8g,%.6g,%.6g\n",nu,H0,Hs,m,s,s/m,sk,md,ma,sa,ska,ra); fflush(stdout);
    }
    return 0;
}
