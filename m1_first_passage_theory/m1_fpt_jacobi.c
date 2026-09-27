/*
 * m1_fpt_jacobi.c — first-passage moments of the reduced recovery process (Proposition 1 of the plan).
 *
 * Between spikes C ≈ C_rest, Q ≈ 1, and H recovers as   dH = (1−H) ds' + ν sqrt(H(1−H)) dW   in rescaled time
 * s' = s / ε_eff  (ε_eff = ε_H / [fP fC Q] evaluated at rest).  The next spike fires when H first reaches H* (= H_on),
 * starting from H0 (= H_min after a spike).  ISI = T_spike + ε_eff τ_FP.
 *
 * Moments by the Dynkin / Siegert formulas with scale density s(y) = y^(−k), speed density m(z) = z^(k−1)/(ν²(1−z)),
 * k = 2/ν²  (0 is an entrance boundary, H* absorbing):
 *   T1(x) = ∫_x^{H*} G(y) dy,           G(y) = ∫_0^1 dw / (1 − y w^{1/k})               (mean)
 *   T2(x) = 2 ∫_x^{H*} ∫_0^1 T1(y w^{1/k}) / (1 − y w^{1/k}) dw dy                       (second moment)
 *   T3(x) = 3 ∫_x^{H*} ∫_0^1 T2(y w^{1/k}) / (1 − y w^{1/k}) dw dy                       (third moment)
 * (the substitution z = y w^{1/k} removes the y^{±k} factors, so the quadrature is stable for any ν).
 * ν → 0 gives G = 1/(1−y) and T1 = ln((1−H0)/(1−H*)), the deterministic recovery time.
 *
 * Usage:
 *   ./m1_fpt_jacobi nu=0.3 H0=0.06 Hs=0.41                      -> one line: nu,H0,Hs,mean,sd,CV,skew,mean_det
 *   ./m1_fpt_jacobi sweep=nu lo=0.02 hi=3 n=25 log=1 H0=.. Hs=..  (or sweep=Hs / sweep=H0)
 *   ./m1_fpt_jacobi mc=20000 nu=0.3 H0=0.06 Hs=0.41 [dt=1e-4]   -> Monte-Carlo check of the 1-D SDE (Euler–Maruyama)
 * Options: ny (grid in y, default 4000), nw (grid in w, default 2000)
 * Build: gcc -O2 -o m1_fpt_jacobi m1_fpt_jacobi.c -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

static int NY=4000, NW=2000;

/* inner integral  ∫_0^1 f(y w^{1/k}) / (1 − y w^{1/k}) dw  by the midpoint rule in w (smooth integrand);
   f is given on the y-grid (linear interpolation); f==NULL means f ≡ 1 */
static double inner(double y, double k, const double *f, const double *yg, int ny){
    double sum=0; int j;
    for(j=0;j<NW;j++){
        double w=(j+0.5)/NW, z=y*pow(w,1.0/k), fz=1.0;
        if(f){ double p=z/yg[ny-1]*(ny-1); int i=(int)p; if(i>=ny-1)i=ny-2; double t=p-i; fz=f[i]*(1-t)+f[i+1]*t; }
        sum+=fz/(1.0-z);
    }
    return sum/NW;
}

/* moments of the first-passage time from H0 to Hs */
static void moments(double nu, double H0, double Hs, double *mean, double *sd, double *skew, double *mean_det){
    double k=2.0/(nu*nu); int ny=NY, i;
    double *yg=malloc(sizeof(double)*ny), *T1=malloc(sizeof(double)*ny), *T2=malloc(sizeof(double)*ny), *T3=malloc(sizeof(double)*ny), *g=malloc(sizeof(double)*ny);
    for(i=0;i<ny;i++) yg[i]=Hs*i/(double)(ny-1);
    /* T1: integrate G from Hs downward */
    for(i=0;i<ny;i++) g[i]=inner(yg[i],k,NULL,yg,ny);
    T1[ny-1]=0; for(i=ny-2;i>=0;i--) T1[i]=T1[i+1]+0.5*(g[i]+g[i+1])*(yg[i+1]-yg[i]);
    for(i=0;i<ny;i++) g[i]=2.0*inner(yg[i],k,T1,yg,ny);
    T2[ny-1]=0; for(i=ny-2;i>=0;i--) T2[i]=T2[i+1]+0.5*(g[i]+g[i+1])*(yg[i+1]-yg[i]);
    for(i=0;i<ny;i++) g[i]=3.0*inner(yg[i],k,T2,yg,ny);
    T3[ny-1]=0; for(i=ny-2;i>=0;i--) T3[i]=T3[i+1]+0.5*(g[i]+g[i+1])*(yg[i+1]-yg[i]);
    /* evaluate at H0 */
    double p=H0/Hs*(ny-1); int j=(int)p; if(j>=ny-1)j=ny-2; double t=p-j;
    double m1=T1[j]*(1-t)+T1[j+1]*t, m2=T2[j]*(1-t)+T2[j+1]*t, m3=T3[j]*(1-t)+T3[j+1]*t;
    double var=m2-m1*m1; *mean=m1; *sd=var>0?sqrt(var):0;
    double mu3=m3-3*m1*m2+2*m1*m1*m1; *skew=(var>0)?mu3/pow(var,1.5):0;
    *mean_det=log((1-H0)/(1-Hs));
    free(yg);free(T1);free(T2);free(T3);free(g);
}

/* Monte-Carlo check */
static uint64_t rs[4];
static inline uint64_t rotl(uint64_t x,int k){return (x<<k)|(x>>(64-k));}
static uint64_t next64(void){uint64_t r=rotl(rs[1]*5,7)*9,t=rs[1]<<17;rs[2]^=rs[0];rs[3]^=rs[1];rs[1]^=rs[2];rs[0]^=rs[3];rs[2]^=t;rs[3]=rotl(rs[3],45);return r;}
static void seed_rng(uint64_t s){int i;for(i=0;i<4;i++){s+=0x9E3779B97F4A7C15ULL;uint64_t z=s;z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL;z=(z^(z>>27))*0x94D049BB133111EBULL;rs[i]=z^(z>>31);}}
static double unif(void){return (next64()>>11)*(1.0/9007199254740992.0);}
static double gauss(void){double u1=unif(),u2=unif();if(u1<1e-300)u1=1e-300;return sqrt(-2*log(u1))*cos(2*M_PI*u2);}
static void montecarlo(long N, double nu, double H0, double Hs, double dt){
    double s1=0,s2=0,s3=0; long n; seed_rng(12345);
    for(n=0;n<N;n++){
        double h=H0, t=0;
        while(h<Hs){ double hh=fmin(fmax(h,0.0),1.0); h+=(1-h)*dt+nu*sqrt(hh*(1-hh))*sqrt(dt)*gauss(); if(h<0)h=-h; t+=dt; }
        s1+=t; s2+=t*t; s3+=t*t*t;
    }
    double m=s1/N, v=s2/N-m*m, mu3=s3/N-3*m*s2/N+2*m*m*m;
    printf("mc,%ld,%.6g,%.6g,%.6g,%.6g\n",N,m,sqrt(v),sqrt(v)/m,mu3/pow(v,1.5));
}

int main(int argc,char**argv){
    double nu=0.3,H0=0.06,Hs=0.41,lo=0,hi=0,dt=1e-4; long mc=0; int n=0,islog=0,i; char sweep[8]="";
    for(i=1;i<argc;i++){
        char key[64],val[64]; if(sscanf(argv[i],"%63[^=]=%63s",key,val)!=2){fprintf(stderr,"bad arg %s\n",argv[i]);return 1;}
        if(!strcmp(key,"nu"))nu=atof(val); else if(!strcmp(key,"H0"))H0=atof(val); else if(!strcmp(key,"Hs"))Hs=atof(val);
        else if(!strcmp(key,"sweep"))strncpy(sweep,val,7); else if(!strcmp(key,"lo"))lo=atof(val); else if(!strcmp(key,"hi"))hi=atof(val);
        else if(!strcmp(key,"n"))n=atoi(val); else if(!strcmp(key,"log"))islog=atoi(val); else if(!strcmp(key,"mc"))mc=atol(val);
        else if(!strcmp(key,"dt"))dt=atof(val); else if(!strcmp(key,"ny"))NY=atoi(val); else if(!strcmp(key,"nw"))NW=atoi(val);
        else {fprintf(stderr,"unknown key %s\n",key);return 1;}
    }
    if(mc>0){ printf("kind,N,mean,sd,CV,skew\n"); montecarlo(mc,nu,H0,Hs,dt); return 0; }
    printf("nu,H0,Hs,mean,sd,CV,skew,mean_det\n");
    int nv=sweep[0]?n:1;
    for(i=0;i<nv;i++){
        if(sweep[0]){ double v=islog?lo*pow(hi/lo,(double)i/(nv-1)):lo+(hi-lo)*i/(double)(nv>1?nv-1:1);
            if(!strcmp(sweep,"nu"))nu=v; else if(!strcmp(sweep,"Hs"))Hs=v; else if(!strcmp(sweep,"H0"))H0=v; else {fprintf(stderr,"sweep must be nu, Hs or H0\n");return 1;} }
        double m,s,sk,md; moments(nu,H0,Hs,&m,&s,&sk,&md);
        printf("%.6g,%.4f,%.4f,%.6g,%.6g,%.6g,%.6g,%.6g\n",nu,H0,Hs,m,s,s/m,sk,md); fflush(stdout);
    }
    return 0;
}
