/*
 * m0_joint_scan.c — joint parameter scan of the ORIGINAL AIMS reaction system M0 (no modified terms).
 *
 *   dC/dt = alpha * P/(1+P) * C/(1+C) * H - beta * C^2/(kS^2 + C^2) - mu*C
 *   dP/dt = gamma * S/(1+S) + rho * C/(kC + C) - lambda*P
 *   dH/dt = (1/epsH) * ( kH^2/(kH^2 + C^2) - H )
 *
 * Question (audit item 5/6-3): is there ANY parameter set of the original equations that gives, over a stimulus
 * range, (i) a period that falls by the data ratio (3 -> 30 µM: ~7x, train-level medians), (ii) a narrow spike
 * (duty = W50/T <= 0.15) and (iii) an amplitude that stays within a factor ~1.5?  A finite scan cannot prove
 * impossibility; it reports how close the best points come and which criteria fail together.
 *
 * Design: Latin-hypercube sample of NPTS points in log10 space over the 10 parameters (ranges below, one decade
 * either side of Table 1 unless stated), each evaluated on a log-spaced stimulus grid S in [Slo, Shi] (NS points).
 * Per point the stimulus curve is summarised:
 *   n_osc        number of grid points that oscillate
 *   T_ratio      max T(S_i)/T(S_j) over i<j with every grid point between i and j oscillating and T non-increasing
 *                (the dose -> S map is monotone but unknown, so any sub-window of the grid may be the data window)
 *   i_lo, i_hi   the window that attains T_ratio
 *   duty_max     max W50/T inside that window;  amp_ratio = amp(S_lo)/amp(S_hi);  W_ratio = W50(S_lo)/W50(S_hi)
 *   T_lo, T_hi   periods at the window ends (model units)
 *   T_1, duty_1, amp_1  features at the middle grid point (for marginal plots)
 * Threads: OpenMP over sample points (each thread integrates its own trajectories).  The integration length
 * scales with epsH (Tend = 4000 + 40 epsH, Ttrans = 0.6 Tend) because the period is ~ 2.8 epsH.
 * Numerics: RK4, dt = 0.02 (dt-convergence of T, W50 to < 0.1 % checked in the audit); half-width crossings
 * are interpolated linearly; clipped states are counted (n_clip) and reported.
 *
 * Usage:
 *   ./m0_joint_scan npts=4000 seed=1 [ns=9 Slo=0.03 Shi=30] [dt=0.02] [threads=N] > results/m0_joint_scan.csv
 *   ./m0_joint_scan curve=1 alpha=.. beta=.. ... [ns=25 Slo=.. Shi=..]      stimulus curve of one point (CSV)
 *   ./m0_joint_scan trace=1 S=1 alpha=.. ...                                  t, C, P, H of the last cycles
 * Ranges (log10): alpha [0.35,35] beta [0.026,2.6] gamma [0.42,42] rho [0.24,24] mu [1e-4,1e-2]
 *                 lambda [1e-4,10] kS [0.02,2] kC [0.008,0.8] kH [0.015,1.5] epsH [7,700]
 * Build: see Makefile (gcc -O2 -fopenmp ... -lm; without OpenMP the pragmas are ignored and it runs serially)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct { double S, alpha, beta, gamma, rho, mu, lambda, kS, kC, kH, epsH; double dt; } Par;
static const char *PNAME[10] = {"alpha","beta","gamma","rho","mu","lambda","kS","kC","kH","epsH"};
static const double PLO[10]  = {0.35, 0.026, 0.42, 0.24, 1e-4, 1e-4, 0.02, 0.008, 0.015, 7};
static const double PHI[10]  = {35.0, 2.6,   42.0, 24.0, 1e-2, 10.0, 2.0,  0.8,   1.5,   700};

static void defaults(Par *p){ p->S=1.0; p->alpha=3.50; p->beta=0.26; p->gamma=4.20; p->rho=2.40; p->mu=0.001; p->lambda=0.001; p->kS=0.20; p->kC=0.08; p->kH=0.15; p->epsH=70.0; p->dt=0.02; }
static double *pfield(Par *p, int i){ double *f[10]={&p->alpha,&p->beta,&p->gamma,&p->rho,&p->mu,&p->lambda,&p->kS,&p->kC,&p->kH,&p->epsH}; return f[i]; }
static double *field(Par *p, const char *k){ int i; for(i=0;i<10;i++) if(!strcmp(k,PNAME[i])) return pfield(p,i); if(!strcmp(k,"S")) return &p->S; if(!strcmp(k,"dt")) return &p->dt; return NULL; }

static inline void rhs(const Par *p, const double *u, double *du){
    double C=u[0],P=u[1],H=u[2],C2=C*C;
    du[0]=p->alpha*P/(1.0+P)*C/(1.0+C)*H - p->beta*C2/(p->kS*p->kS+C2) - p->mu*C;
    du[1]=p->gamma*p->S/(1.0+p->S) + p->rho*C/(p->kC+C) - p->lambda*P;
    du[2]=(p->kH*p->kH/(p->kH*p->kH+C2) - H)/p->epsH;
}
static inline int rk4(const Par *p, double *u, double dt){
    double k1[3],k2[3],k3[3],k4[3],v[3]; int i, clip=0;
    rhs(p,u,k1); for(i=0;i<3;i++) v[i]=u[i]+0.5*dt*k1[i]; rhs(p,v,k2);
    for(i=0;i<3;i++) v[i]=u[i]+0.5*dt*k2[i]; rhs(p,v,k3); for(i=0;i<3;i++) v[i]=u[i]+dt*k3[i]; rhs(p,v,k4);
    for(i=0;i<3;i++) u[i]+=dt/6.0*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
    if(u[0]<0){u[0]=0;clip=1;} if(u[1]<0){u[1]=0;clip=1;} if(u[2]<0){u[2]=0;clip=1;} if(u[2]>1){u[2]=1;clip=1;}
    return clip;
}

typedef struct { int osc; double T, Cmax, Cmin, amp, W50, duty, Hmin, Hmax, Pmean; long nclip; } Feat;
#define NPEAK 6
/* linear-interpolated crossing time of level y between samples i (above) and i-1 or i+1 */
static double cross_left(const double *y, long a, double lev, double dt){ double y0=y[a], y1=y[a+1]; return (a + (lev-y0)/(y1-y0>0?y1-y0:1e-300))*dt; }
static double cross_right(const double *y, long b, double lev, double dt){ double y0=y[b-1], y1=y[b]; return (b-1 + (y0-lev)/(y0-y1>0?y0-y1:1e-300))*dt; }

static Feat run(const Par *p, int trace){
    Feat f; memset(&f,0,sizeof f);
    double u[3]={0.10,0.05,1.0}, dt=p->dt;
    double Tend=4000.0+40.0*p->epsH, Ttrans=0.6*Tend;
    long nsteps=(long)(Tend/dt+0.5), ntr=(long)(Ttrans/dt+0.5), nrec=nsteps-ntr, n, i;
    double *Cb=malloc(sizeof(double)*nrec), *Hb=malloc(sizeof(double)*nrec), *Pb=malloc(sizeof(double)*nrec);
    for(n=0;n<nsteps;n++){ if(n>=ntr){ Cb[n-ntr]=u[0]; Pb[n-ntr]=u[1]; Hb[n-ntr]=u[2]; } f.nclip+=rk4(p,u,dt); }
    double cmax=-1,cmin=1e300,hmin=1e300,hmax=-1,ps=0;
    for(i=0;i<nrec;i++){ if(Cb[i]>cmax)cmax=Cb[i]; if(Cb[i]<cmin)cmin=Cb[i]; if(Hb[i]<hmin)hmin=Hb[i]; if(Hb[i]>hmax)hmax=Hb[i]; ps+=Pb[i]; }
    f.Cmax=cmax; f.Cmin=cmin; f.Hmin=hmin; f.Hmax=hmax; f.Pmean=ps/nrec; f.amp=cmax-cmin;
    if(cmax-cmin<1e-4 || cmax-cmin<0.05*cmax){ free(Cb);free(Hb);free(Pb); return f; }
    double mid=0.5*(cmax+cmin); long pk[NPEAK+1]; int np=0;
    for(i=nrec-2;i>0&&np<NPEAK;i--){ if(Cb[i]>mid&&Cb[i]>=Cb[i-1]&&Cb[i]>Cb[i+1]){ pk[np++]=i; long j=i; while(j>0&&Cb[j]>mid) j--; i=j; } }
    if(np<3){ free(Cb);free(Hb);free(Pb); return f; }
    f.osc=1; double Ts=0; int k; for(k=0;k<np-1;k++) Ts+=(pk[k]-pk[k+1])*dt; f.T=Ts/(np-1);
    long ip=pk[1]; double base=1e300; for(i=pk[2];i<=ip;i++) if(Cb[i]<base) base=Cb[i];
    double half=base+0.5*(Cb[ip]-base); long a=ip,b=ip; while(a>0&&Cb[a]>half) a--; while(b<nrec-1&&Cb[b]>half) b++;
    f.W50=cross_right(Cb,b,half,dt)-cross_left(Cb,a,half,dt); f.duty=f.W50/f.T; f.amp=Cb[ip]-base;
    if(trace){ long start=pk[np-1]-(long)(0.1*f.T/dt); if(start<0)start=0; printf("t,C,P,H\n"); for(i=start;i<nrec;i+=(long)fmax(1,0.1/dt)) printf("%.3f,%.6g,%.6g,%.6g\n",Ttrans+i*dt,Cb[i],Pb[i],Hb[i]); }
    free(Cb);free(Hb);free(Pb); return f;
}

/* summary of a stimulus curve */
typedef struct { int n_osc, i_lo, i_hi, monotone_dec; double T_ratio, duty_max, amp_ratio, W_ratio, T_lo, T_hi, T_mid, duty_mid, amp_mid, duty_min; long nclip; } Summ;
static Summ summarise(const Feat *F, int ns){
    Summ s; memset(&s,0,sizeof s); int i,j; s.T_ratio=1.0; s.i_lo=s.i_hi=-1; s.duty_min=1e300;
    for(i=0;i<ns;i++){ s.n_osc+=F[i].osc; s.nclip+=F[i].nclip; if(F[i].osc&&F[i].duty<s.duty_min) s.duty_min=F[i].duty; }
    if(s.duty_min>1e299) s.duty_min=0;
    for(i=0;i<ns;i++){ if(!F[i].osc) continue;
        for(j=i+1;j<ns;j++){ if(!F[j].osc || F[j].T>F[j-1].T*1.0001) break;      /* stop when non-oscillating or T increases */
            double r=F[i].T/F[j].T; if(r>s.T_ratio){ s.T_ratio=r; s.i_lo=i; s.i_hi=j; } } }
    if(s.i_lo>=0){ double dm=0; for(i=s.i_lo;i<=s.i_hi;i++) if(F[i].duty>dm) dm=F[i].duty;
        s.duty_max=dm; s.amp_ratio=F[s.i_lo].amp/F[s.i_hi].amp; s.W_ratio=F[s.i_lo].W50/F[s.i_hi].W50; s.T_lo=F[s.i_lo].T; s.T_hi=F[s.i_hi].T; }
    int m=ns/2; s.T_mid=F[m].osc?F[m].T:0; s.duty_mid=F[m].osc?F[m].duty:0; s.amp_mid=F[m].osc?F[m].amp:0;
    /* monotone decreasing over all oscillating points? */
    s.monotone_dec=1; for(i=1;i<ns;i++) if(F[i].osc&&F[i-1].osc&&F[i].T>F[i-1].T*1.0001) s.monotone_dec=0;
    return s;
}

/* xoshiro256** for the Latin hypercube */
static uint64_t rs[4];
static inline uint64_t rotl(uint64_t x,int k){return (x<<k)|(x>>(64-k));}
static uint64_t next64(void){uint64_t r=rotl(rs[1]*5,7)*9,t=rs[1]<<17;rs[2]^=rs[0];rs[3]^=rs[1];rs[1]^=rs[2];rs[0]^=rs[3];rs[2]^=t;rs[3]=rotl(rs[3],45);return r;}
static void seed_rng(uint64_t s){int i;for(i=0;i<4;i++){s+=0x9E3779B97F4A7C15ULL;uint64_t z=s;z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL;z=(z^(z>>27))*0x94D049BB133111EBULL;rs[i]=z^(z>>31);}}
static double unif(void){return (next64()>>11)*(1.0/9007199254740992.0);}

int main(int argc,char**argv){
    Par base; defaults(&base);
    long npts=0; int ns=9, curve=0, trace=0, i, threads=0; double Slo=0.03, Shi=30; uint64_t seed=1;
    for(i=1;i<argc;i++){
        char key[64],val[64]; if(sscanf(argv[i],"%63[^=]=%63s",key,val)!=2){fprintf(stderr,"bad arg %s\n",argv[i]);return 1;}
        if(!strcmp(key,"npts")) npts=atol(val); else if(!strcmp(key,"ns")) ns=atoi(val); else if(!strcmp(key,"Slo")) Slo=atof(val); else if(!strcmp(key,"Shi")) Shi=atof(val);
        else if(!strcmp(key,"seed")) seed=strtoull(val,NULL,10); else if(!strcmp(key,"curve")) curve=atoi(val); else if(!strcmp(key,"trace")) trace=atoi(val);
        else if(!strcmp(key,"threads")) threads=atoi(val);
        else { double *f=field(&base,key); if(!f){fprintf(stderr,"unknown key %s\n",key);return 1;} *f=atof(val); }
    }
#ifdef _OPENMP
    if(threads>0) omp_set_num_threads(threads);
    fprintf(stderr,"OpenMP: %d threads\n", omp_get_max_threads());
#else
    fprintf(stderr,"built without OpenMP: serial\n");
#endif
    if(trace){ run(&base,1); return 0; }
    if(curve){ printf("S,osc,T,Cmax,Cmin,amp,W50,duty,Hmin,Hmax,Pmean,nclip\n");
        for(i=0;i<ns;i++){ Par p=base; p.S=Slo*pow(Shi/Slo,(double)i/(ns-1)); Feat f=run(&p,0);
            printf("%.6g,%d,%.4f,%.5g,%.5g,%.5g,%.4f,%.4f,%.4f,%.4f,%.5g,%ld\n",p.S,f.osc,f.T,f.Cmax,f.Cmin,f.amp,f.W50,f.duty,f.Hmin,f.Hmax,f.Pmean,f.nclip); }
        return 0; }
    if(npts<=0){ fprintf(stderr,"give npts=N (scan), curve=1 or trace=1\n"); return 1; }
    /* Latin hypercube in log10 space */
    double *X=malloc(sizeof(double)*npts*10); seed_rng(seed); int d; long n;
    for(d=0;d<10;d++){ long *perm=malloc(sizeof(long)*npts); for(n=0;n<npts;n++) perm[n]=n;
        for(n=npts-1;n>0;n--){ long j=(long)(unif()*(n+1)); long t=perm[n]; perm[n]=perm[j]; perm[j]=t; }
        for(n=0;n<npts;n++){ double u=(perm[n]+unif())/npts; X[n*10+d]=pow(10.0, log10(PLO[d])+u*(log10(PHI[d])-log10(PLO[d]))); }
        free(perm); }
    Summ *R=malloc(sizeof(Summ)*npts);
    long done=0;
#pragma omp parallel for schedule(dynamic,1)
    for(n=0;n<npts;n++){
        Par p=base; int k; for(k=0;k<10;k++) *pfield(&p,k)=X[n*10+k];
        Feat *F=malloc(sizeof(Feat)*ns);
        for(k=0;k<ns;k++){ p.S=Slo*pow(Shi/Slo,(double)k/(ns-1)); F[k]=run(&p,0); }
        R[n]=summarise(F,ns); free(F);
#pragma omp atomic
        done++;
        if(done%200==0){
#pragma omp critical
            fprintf(stderr,"  %ld / %ld\n",done,npts); }
    }
    printf("idx"); for(d=0;d<10;d++) printf(",%s",PNAME[d]);
    printf(",n_osc,monotone_dec,i_lo,i_hi,S_lo,S_hi,T_ratio,T_lo,T_hi,duty_max,duty_min,amp_ratio,W_ratio,T_mid,duty_mid,amp_mid,nclip\n");
    for(n=0;n<npts;n++){ Summ *s=&R[n]; printf("%ld",n); for(d=0;d<10;d++) printf(",%.5g",X[n*10+d]);
        double slo=s->i_lo>=0?Slo*pow(Shi/Slo,(double)s->i_lo/(ns-1)):0, shi=s->i_hi>=0?Slo*pow(Shi/Slo,(double)s->i_hi/(ns-1)):0;
        printf(",%d,%d,%d,%d,%.5g,%.5g,%.5g,%.5g,%.5g,%.5g,%.5g,%.5g,%.5g,%.5g,%.5g,%.5g,%ld\n",s->n_osc,s->monotone_dec,s->i_lo,s->i_hi,slo,shi,s->T_ratio,s->T_lo,s->T_hi,s->duty_max,s->duty_min,s->amp_ratio,s->W_ratio,s->T_mid,s->duty_mid,s->amp_mid,s->nclip); }
    free(X); free(R); return 0;
}
