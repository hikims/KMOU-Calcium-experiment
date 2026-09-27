/*
 * emrich_jin_sde.c — M1 (bounded H-noise on the stereotyped reference) under reduced basal influx Jin:
 *                    the prediction for SOCE knock-outs (Emrich 2021, ORAI-TKO / STIM KO), OpenMP.
 *
 * Model = m1_isi_variability/m1_isi_sde.c (Euler–Maruyama, Jacobi noise sigH sqrt(H(1−H)) dW_H, hysteresis spike
 * detection with thresholds from a deterministic pre-run).  Three modes:
 *   sweep=Jin lo=.. hi=.. n=.. [log=0] ...     ISI statistics per value (deterministic period + stochastic CV etc.);
 *                                              the sweep values run in parallel (one RNG stream per value)
 *   cells=N Twin=T [Tpre=T0] [Spre=0] [epsH_cv=X] [seed=S] ...
 *       N cells in parallel.  Cell i has epsH_i = epsH exp(epsH_cv * xi_i) (xi ~ N(0,1), lognormal spread of the
 *       slow time scale; the rescaled noise nu = sigH sqrt(epsH) is common when nuH is given).  Protocol: rest at
 *       stimulus Spre for Tpre (noise on), then the stimulus is switched to S and spikes are counted in [0, Twin]
 *       (Emrich Fig. 1: CCh added, spikes counted over 14 min).  Output per cell: cell, epsH_i, n_spikes,
 *       t_first, Tav (mean ISI in the window), amp_mean, Hpre (H at the switch).
 *   trace=1 [tspan=X] ...                      t, C, P, H of one cell (same protocol) for a figure
 * Keys: model keys of m1_isi_sde + sigH sigP sigC nuH seed dt Tend Ttrans + Twin Tpre Spre epsH_cv cells threads
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

typedef struct {
    double S, alpha, beta, gamma, rho, mu, lambda, kS, kC, kH, epsH;
    double mech, gH, dH, qoff, kP, kR, nA, kA, Jin, rC;
    double C0, P0, H0, dt, Tend, Ttrans, tspan;
    double sigH, sigP, sigC, seed, nuH, Twin, Tpre, Spre, epsH_cv;
} Par;
static void defaults(Par *p){
    p->S=0; p->alpha=3.5; p->beta=0.26; p->gamma=4.2; p->rho=2.4; p->mu=0.001; p->lambda=0.001; p->kS=0.2; p->kC=0.08; p->kH=0.15; p->epsH=70;
    p->mech=0; p->gH=0; p->dH=1; p->qoff=1; p->kP=1; p->kR=-1; p->nA=1; p->kA=1; p->Jin=0; p->rC=0;
    p->C0=0.1; p->P0=0.05; p->H0=1; p->dt=0.01; p->Tend=60000; p->Ttrans=3000; p->tspan=600;
    p->sigH=0; p->sigP=0; p->sigC=0; p->seed=1; p->nuH=0; p->Twin=1300; p->Tpre=3000; p->Spre=0; p->epsH_cv=0;
}
#define KEY(name) if(!strcmp(k,#name)) return &p->name;
static double *field(Par *p, const char *k){
    KEY(S) KEY(alpha) KEY(beta) KEY(gamma) KEY(rho) KEY(mu) KEY(lambda) KEY(kS) KEY(kC) KEY(kH) KEY(epsH) KEY(mech) KEY(gH) KEY(dH) KEY(qoff) KEY(kP) KEY(kR) KEY(nA) KEY(kA) KEY(Jin) KEY(rC)
    KEY(C0) KEY(P0) KEY(H0) KEY(dt) KEY(Tend) KEY(Ttrans) KEY(tspan) KEY(sigH) KEY(sigP) KEY(sigC) KEY(seed) KEY(nuH) KEY(Twin) KEY(Tpre) KEY(Spre) KEY(epsH_cv)
    return NULL;
}
typedef struct { uint64_t s[4]; int has; double g2; } Rng;
static inline uint64_t rotl(uint64_t x,int k){return (x<<k)|(x>>(64-k));}
static inline uint64_t next64(Rng *r){uint64_t *s=r->s; uint64_t res=rotl(s[1]*5,7)*9,t=s[1]<<17;s[2]^=s[0];s[3]^=s[1];s[1]^=s[2];s[0]^=s[3];s[2]^=t;s[3]=rotl(s[3],45);return res;}
static void seed_rng(Rng *r, uint64_t s){int i;r->has=0;for(i=0;i<4;i++){s+=0x9E3779B97F4A7C15ULL;uint64_t z=s;z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL;z=(z^(z>>27))*0x94D049BB133111EBULL;r->s[i]=z^(z>>31);}}
static inline double unif(Rng *r){return (next64(r)>>11)*(1.0/9007199254740992.0);}
static inline double gauss(Rng *r){ if(r->has){r->has=0;return r->g2;} double u1=unif(r),u2=unif(r); if(u1<1e-300)u1=1e-300; double a=sqrt(-2*log(u1)); r->g2=a*sin(2*M_PI*u2); r->has=1; return a*cos(2*M_PI*u2); }

static inline void rhs(const Par *p, const double *u, double *du){
    double C=u[0],P=u[1],H=u[2],C2=C*C, Cn=pow(C,p->nA), act=Cn/(pow(p->kA,p->nA)+Cn);
    du[0]=p->Jin + p->alpha*P/(p->kP+P)*act*H - p->beta*C2/(p->kS*p->kS+C2) - p->mu*C;
    du[1]=p->gamma*p->S/(1.0+p->S) + p->rho*C/(p->kC+C) - p->lambda*P;
    double kh=p->kH, fP=1, fPq=1, kR=p->kR>0?p->kR:p->kH, fC=1.0+p->rC*C2/(kR*kR+C2);
    if((int)p->mech==1) kh=p->kH*(p->dH+p->gH*P)/(p->dH+P); else if((int)p->mech>=2){ double x=p->gH*P/(p->dH+P); fP=1+x; fPq=1+p->qoff*x; }
    double Q=kh*kh/(kh*kh+C2); du[2]=(fP*fC*Q*(1.0-H) - fPq*fC*(1.0-Q)*H)/p->epsH;
}
static inline void em_step(const Par *p, double *u, double dt, int stoch, Rng *r){
    double du[3]; rhs(p,u,du); double sq=sqrt(dt);
    u[0]+=du[0]*dt; u[1]+=du[1]*dt; u[2]+=du[2]*dt;
    if(stoch){ if(p->sigC>0) u[0]+=p->sigC*sqrt(fmax(u[0],0.0))*sq*gauss(r); if(p->sigP>0) u[1]+=p->sigP*fmax(u[1],0.0)*sq*gauss(r);
               if(p->sigH>0){ double h=fmin(fmax(u[2],0.0),1.0); u[2]+=p->sigH*sqrt(h*(1.0-h))*sq*gauss(r); } }
    if(u[0]<0)u[0]=-u[0]; if(u[1]<0)u[1]=-u[1]; if(u[2]<0)u[2]=-u[2]; if(u[2]>1)u[2]=2.0-u[2];
}
typedef struct { int osc; double T, base, amp, Hon, W50; double u[3]; } Det;
static Det det_run(const Par *p){
    Det d; memset(&d,0,sizeof d); double u[3]={p->C0,p->P0,p->H0}; double dt=fmin(p->dt,0.02);
    long ntr=(long)(p->Ttrans/dt), nrec=(long)(fmax(4000.0,p->Ttrans)/dt), n;
    for(n=0;n<ntr;n++) em_step(p,u,dt,0,NULL);
    double *Cb=malloc(sizeof(double)*nrec), *Hb=malloc(sizeof(double)*nrec);
    for(n=0;n<nrec;n++){ Cb[n]=u[0]; Hb[n]=u[2]; em_step(p,u,dt,0,NULL); }
    d.u[0]=u[0]; d.u[1]=u[1]; d.u[2]=u[2];
    double cmax=-1,cmin=1e300; long i; for(i=0;i<nrec;i++){ if(Cb[i]>cmax)cmax=Cb[i]; if(Cb[i]<cmin)cmin=Cb[i]; }
    d.base=cmin; d.amp=cmax-cmin;
    if(cmax-cmin<1e-4||cmax-cmin<0.05*cmax){ free(Cb);free(Hb); return d; }
    double mid=0.5*(cmax+cmin); long pk[8]; int np=0;
    for(i=nrec-2;i>0&&np<8;i--){ if(Cb[i]>mid&&Cb[i]>=Cb[i-1]&&Cb[i]>Cb[i+1]){ pk[np++]=i; long j=i; while(j>0&&Cb[j]>mid) j--; i=j; } }
    if(np<3){ free(Cb);free(Hb); return d; }
    d.osc=1; double Ts=0; int k; for(k=0;k<np-1;k++) Ts+=(pk[k]-pk[k+1])*dt; d.T=Ts/(np-1);
    long ip=pk[1]; double base=1e300; for(i=pk[2];i<=ip;i++) if(Cb[i]<base) base=Cb[i];
    d.base=base; d.amp=Cb[ip]-base; double half=base+0.5*d.amp; long a=ip; while(a>0&&Cb[a]>half) a--; long b=ip; while(b<nrec-1&&Cb[b]>half) b++;
    d.W50=(b-a)*dt; d.Hon=Hb[a]; free(Cb);free(Hb); return d;
}
typedef struct { double t, Hon, amp; } Spk;

/* stationary statistics of one parameter set (as m1_isi_sde) */
static void stats_line(Par p, double v, Rng *r, char *buf){
    if(p.nuH>0) p.sigH=p.nuH/sqrt(p.epsH);
    Det d=det_run(&p);
    if(!d.osc){ snprintf(buf,512,"%.6g,0,0,0,0,0,0,0,0,0,%.5g,0,0,0,0,0\n",v,d.base+d.amp); return; }
    double thr_hi=d.base+0.5*d.amp, thr_lo=d.base+0.2*d.amp, u[3]={d.u[0],d.u[1],d.u[2]};
    long nsteps=(long)(p.Tend/p.dt), ntr=(long)(p.Ttrans/p.dt), k; int cap=100000, ns=0; Spk *sp=malloc(sizeof(Spk)*cap); int armed=1,inspike=0; double pkC=-1;
    for(k=0;k<nsteps;k++){ double Hprev=u[2]; em_step(&p,u,p.dt,1,r);
        if(armed&&u[0]>=thr_hi){ armed=0; inspike=1; pkC=u[0]; if(k>=ntr&&ns<cap){ sp[ns].t=k*p.dt; sp[ns].Hon=Hprev; sp[ns].amp=0; ns++; } }
        if(inspike){ if(u[0]>pkC) pkC=u[0]; if(u[0]<thr_lo){ inspike=0; armed=1; if(ns>0&&k>=ntr) sp[ns-1].amp=pkC-d.base; } } }
    int m=0,i; double *x=malloc(sizeof(double)*(ns+1)), *A=malloc(sizeof(double)*(ns+1)), *Hn=malloc(sizeof(double)*(ns+1));
    for(i=1;i<ns;i++) if(sp[i].amp>0){ x[m]=sp[i].t-sp[i-1].t; A[m]=sp[i].amp; Hn[m]=sp[i].Hon; m++; }
    double s1=0,s2=0,s3=0,am=0,a2=0,hm=0,h2=0; for(i=0;i<m;i++){ s1+=x[i]; am+=A[i]; hm+=Hn[i]; }
    double Tav=m?s1/m:0, amean=m?am/m:0, hmean=m?hm/m:0;
    for(i=0;i<m;i++){ double e=x[i]-Tav; s2+=e*e; s3+=e*e*e; a2+=(A[i]-amean)*(A[i]-amean); h2+=(Hn[i]-hmean)*(Hn[i]-hmean); }
    double sd=m>1?sqrt(s2/(m-1)):0, skew=(m>2&&s2>0)?(s3/m)/pow(s2/m,1.5):0, asd=m>1?sqrt(a2/(m-1)):0, hsd=m>1?sqrt(h2/(m-1)):0;
    double cxa=0; for(i=0;i<m;i++) cxa+=(x[i]-Tav)*(A[i]-amean); double r_amp=(s2>0&&a2>0)?cxa/sqrt(s2*a2):0;
    double mx1=0,mx2=0,cxx=0,sxx1=0,sxx2=0; for(i=0;i<m-1;i++){ mx1+=x[i]; mx2+=x[i+1]; } if(m>2){ mx1/=(m-1); mx2/=(m-1); }
    for(i=0;i<m-1;i++){ cxx+=(x[i]-mx1)*(x[i+1]-mx2); sxx1+=(x[i]-mx1)*(x[i]-mx1); sxx2+=(x[i+1]-mx2)*(x[i+1]-mx2); } double r_next=(sxx1>0&&sxx2>0)?cxx/sqrt(sxx1*sxx2):0;
    snprintf(buf,512,"%.6g,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.5g,%.5g,%.4f,%.4f,%.4f,%.4f\n",v,m,d.T,Tav,sd,Tav>0?sd/Tav:0,skew,d.Hon,hmean,hsd,d.amp,amean,amean>0?asd/amean:0,r_amp,r_next,d.W50);
    free(sp);free(x);free(A);free(Hn);
}

/* one cell under the step protocol: returns n_spikes, t_first, Tav, amp_mean, Hpre */
static void cell_run(Par p, Rng *r, int trace, double *out){
    if(p.nuH>0) p.sigH=p.nuH/sqrt(p.epsH);
    Det d=det_run(&p);                                   /* thresholds at the stimulated level */
    double thr_hi=d.base+0.5*fmax(d.amp,1e-3), thr_lo=d.base+0.2*fmax(d.amp,1e-3);
    Par q=p; q.S=p.Spre; double u[3]={p.C0,p.P0,p.H0}; long n, npre=(long)(p.Tpre/p.dt), nwin=(long)(p.Twin/p.dt), every=(long)fmax(1,0.2/p.dt);
    for(n=0;n<npre;n++) em_step(&q,u,p.dt,1,r);
    double Hpre=u[2]; int ns=0, armed=1, inspike=0; double tfirst=-1, tlast=-1, pkC=-1, asum=0; int na=0; double base=d.base;
    if(trace) printf("t,C,P,H\n");
    for(n=0;n<nwin;n++){ em_step(&p,u,p.dt,1,r); double t=n*p.dt;
        if(trace&&n%every==0) printf("%.3f,%.6g,%.6g,%.6g\n",t,u[0],u[1],u[2]);
        if(armed&&u[0]>=thr_hi){ armed=0; inspike=1; pkC=u[0]; ns++; if(tfirst<0) tfirst=t; tlast=t; }
        if(inspike){ if(u[0]>pkC) pkC=u[0]; if(u[0]<thr_lo){ inspike=0; armed=1; asum+=pkC-base; na++; } } }
    out[0]=ns; out[1]=tfirst; out[2]=(ns>1)?(tlast-tfirst)/(ns-1):0; out[3]=na?asum/na:0; out[4]=Hpre;
}

int main(int argc,char**argv){
    Par p; defaults(&p); char sweep[32]=""; double lo=0,hi=0; int n=0,islog=0,trace=0,cells=0,threads=0,i;
    for(i=1;i<argc;i++){ char key[64],val[64]; if(sscanf(argv[i],"%63[^=]=%63s",key,val)!=2){fprintf(stderr,"bad arg %s\n",argv[i]);return 1;}
        if(!strcmp(key,"sweep")) strncpy(sweep,val,31); else if(!strcmp(key,"lo")) lo=atof(val); else if(!strcmp(key,"hi")) hi=atof(val); else if(!strcmp(key,"n")) n=atoi(val);
        else if(!strcmp(key,"log")) islog=atoi(val); else if(!strcmp(key,"trace")) trace=atoi(val); else if(!strcmp(key,"cells")) cells=atoi(val); else if(!strcmp(key,"threads")) threads=atoi(val);
        else { double *f=field(&p,key); if(!f){fprintf(stderr,"unknown key %s\n",key);return 1;} *f=atof(val); } }
#ifdef _OPENMP
    if(threads>0) omp_set_num_threads(threads); fprintf(stderr,"OpenMP: %d threads\n",omp_get_max_threads());
#endif
    if(trace){ Rng r; seed_rng(&r,(uint64_t)p.seed); double out[5]; cell_run(p,&r,1,out); return 0; }
    if(cells>0){
        double *O=malloc(sizeof(double)*cells*6); int c;
#pragma omp parallel for schedule(dynamic,1)
        for(c=0;c<cells;c++){ Rng r; seed_rng(&r,(uint64_t)p.seed*1000003ULL+c); Par q=p; double xi=gauss(&r); q.epsH=p.epsH*exp(p.epsH_cv*xi); double out[5]; cell_run(q,&r,0,out);
            O[c*6]=q.epsH; O[c*6+1]=out[0]; O[c*6+2]=out[1]; O[c*6+3]=out[2]; O[c*6+4]=out[3]; O[c*6+5]=out[4]; }
        printf("cell,epsH,n_spikes,t_first,Tav,amp_mean,Hpre\n");
        for(c=0;c<cells;c++) printf("%d,%.5g,%d,%.3f,%.4f,%.5g,%.5f\n",c,O[c*6],(int)O[c*6+1],O[c*6+2],O[c*6+3],O[c*6+4],O[c*6+5]);
        free(O); return 0; }
    printf("value,n,Tdet,Tav,sd,CV,skew,Hon_det,Hon_mean,Hon_sd,amp_det,amp_mean,amp_cv,corr_isi_amp,corr_isi_next,W50_det\n");
    double *sw=sweep[0]?field(&p,sweep):NULL; if(sweep[0]&&!sw){fprintf(stderr,"unknown sweep key %s\n",sweep);return 1;}
    int nv=sweep[0]?n:1, iv;
    /* parallel over sweep values: results are collected in order */
    char **lines=malloc(sizeof(char*)*nv);
#pragma omp parallel for schedule(dynamic,1)
    for(iv=0;iv<nv;iv++){ Par q=p; double v=0; if(sw){ v=islog?lo*pow(hi/lo,(double)iv/(nv-1)):lo+(hi-lo)*iv/(double)(nv>1?nv-1:1); *field(&q,sweep)=v; }
        Rng r; seed_rng(&r,(uint64_t)p.seed+7919ULL*iv);
        lines[iv]=malloc(512); stats_line(q,v,&r,lines[iv]); }
    for(iv=0;iv<nv;iv++){ fputs(lines[iv],stdout); free(lines[iv]); } free(lines);
    return 0;
}
