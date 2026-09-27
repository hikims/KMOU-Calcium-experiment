/*
 * m1_isi_sde.c — M1: bounded noise in the receptor-recovery variable H; ISI statistics.
 *
 * Deterministic part = stereotyped reference of m0ip3h_amplitude_invariance (all switches kept, defaults = original AIMS):
 *   dC = [ Jin + alpha P/(kP+P) C^nA/(kA^nA+C^nA) H - beta C^2/(kS^2+C^2) - mu C ] dt            (+ sigC sqrt(C) dW_C)
 *   dP = [ gamma S/(1+S) + rho C/(kC+C) - lambda P ] dt                                           (+ sigP P dW_P)
 *   dH = (1/epsH) [ fP fC Q (1-H) - fPq fC (1-Q) H ] dt  +  sigH sqrt(H(1-H)) dW_H               (M1: Jacobi-type, keeps 0<=H<=1)
 *   Q = kH^2/(kH^2+C^2), fP = 1+gH P/(dH+P), fPq = 1+qoff gH P/(dH+P), fC = 1+rC C^2/(kR^2+C^2)
 * Noise switches: sigH (M1), sigP (M2), sigC (Ca-triggered comparison); any combination (M12 = sigH & sigP).
 *   nuH > 0 sets sigH = nuH/sqrt(epsH) (common rescaled intensity nu^2 = sigH^2 epsH across cells that differ in epsH).
 * Euler–Maruyama, reflection at C>=0, P>=0, 0<=H<=1.  The noise amplitudes are in model time units.
 *
 * Spike detection (hysteresis): a spike starts when C crosses up through thr_hi = base + 0.5 amp_det and the next
 * one is allowed after C has fallen below thr_lo = base + 0.2 amp_det, where base, amp_det come from a deterministic
 * pre-run (all sigmas = 0) with the same parameters.  Per spike: onset time, H at onset, peak C.
 *
 * Usage:
 *   ./m1_isi_sde [key=value ...]                          one run -> one CSV line of ISI statistics
 *   ./m1_isi_sde sweep=NAME lo=A hi=B n=N [log=1] ...     sweep NAME (e.g. sigH, epsH, S), one line per value
 *   ./m1_isi_sde isi=1 ...                                list: onset time, ISI (preceding), H_onset, amplitude
 *   ./m1_isi_sde trace=1 [tspan=X] ...                    t, C, P, H after the transient (length tspan)
 * Keys: model keys as in m0ip3h_invariance_sweep + sigH sigP sigC nuH seed dt Tend Ttrans tspan
 * Output columns (statistics):
 *   value, n, Tdet, Tav, sd, CV, skew, Hon_det, Hon_mean, Hon_sd, amp_det, amp_mean, amp_cv, corr_isi_amp, corr_isi_next, W50_det
 *   Tdet/Hon_det/amp_det/W50_det = deterministic period / H at onset / amplitude / half-width
 *   corr_isi_amp = corr(preceding ISI, amplitude);  corr_isi_next = corr(ISI_k, ISI_{k+1}) (serial correlation)
 * Build: gcc -O2 -o m1_isi_sde m1_isi_sde.c -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

typedef struct {
    double S, alpha, beta, gamma, rho, mu, lambda, kS, kC, kH, epsH;
    double mech, gH, dH, qoff, kP, kR, nA, kA, Jin, rC;
    double C0, P0, H0, dt, Tend, Ttrans, tspan;
    double sigH, sigP, sigC, seed, nuH;
} Par;

static void defaults(Par *p) {
    p->S = 0.0;  p->alpha = 3.50; p->beta = 0.26; p->gamma = 4.20; p->rho = 2.40;
    p->mu = 0.001; p->lambda = 0.001; p->kS = 0.20; p->kC = 0.08; p->kH = 0.15; p->epsH = 70.0;
    p->mech = 0; p->gH = 0.0; p->dH = 1.0; p->qoff = 1.0; p->kP = 1.0; p->kR = -1.0;
    p->nA = 1.0; p->kA = 1.0; p->Jin = 0.0; p->rC = 0.0;
    p->C0 = 0.10; p->P0 = 0.05; p->H0 = 1.0;
    p->dt = 0.01; p->Tend = 60000.0; p->Ttrans = 3000.0; p->tspan = 600.0;
    p->sigH = 0.0; p->sigP = 0.0; p->sigC = 0.0; p->seed = 1; p->nuH = 0.0;
}
#define KEY(name) if(!strcmp(k,#name)) return &p->name;
static double *field(Par *p, const char *k) {
    KEY(S) KEY(alpha) KEY(beta) KEY(gamma) KEY(rho) KEY(mu) KEY(lambda) KEY(kS) KEY(kC) KEY(kH) KEY(epsH)
    KEY(mech) KEY(gH) KEY(dH) KEY(qoff) KEY(kP) KEY(kR) KEY(nA) KEY(kA) KEY(Jin) KEY(rC)
    KEY(C0) KEY(P0) KEY(H0) KEY(dt) KEY(Tend) KEY(Ttrans) KEY(tspan) KEY(sigH) KEY(sigP) KEY(sigC) KEY(seed) KEY(nuH)
    return NULL;
}

/* xoshiro256** + Box–Muller */
static uint64_t rs[4];
static inline uint64_t rotl(uint64_t x, int k){ return (x<<k)|(x>>(64-k)); }
static uint64_t next64(void){ uint64_t r=rotl(rs[1]*5,7)*9, t=rs[1]<<17; rs[2]^=rs[0]; rs[3]^=rs[1]; rs[1]^=rs[2]; rs[0]^=rs[3]; rs[2]^=t; rs[3]=rotl(rs[3],45); return r; }
static void seed_rng(uint64_t s){ int i; for(i=0;i<4;i++){ s+=0x9E3779B97F4A7C15ULL; uint64_t z=s; z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL; z=(z^(z>>27))*0x94D049BB133111EBULL; rs[i]=z^(z>>31);} }
static double unif(void){ return (next64()>>11)*(1.0/9007199254740992.0); }
static double gauss(void){ static int has=0; static double g2; if(has){has=0;return g2;} double u1=unif(),u2=unif(); if(u1<1e-300)u1=1e-300; double r=sqrt(-2*log(u1)); g2=r*sin(2*M_PI*u2); has=1; return r*cos(2*M_PI*u2); }

static void rhs(const Par *p, const double *u, double *du) {
    double C=u[0], P=u[1], H=u[2], C2=C*C;
    double Cn=pow(C,p->nA), act=Cn/(pow(p->kA,p->nA)+Cn);
    du[0]=p->Jin + p->alpha*P/(p->kP+P)*act*H - p->beta*C2/(p->kS*p->kS+C2) - p->mu*C;
    du[1]=p->gamma*p->S/(1.0+p->S) + p->rho*C/(p->kC+C) - p->lambda*P;
    double kh=p->kH, fP=1.0, fPq=1.0, kR=p->kR>0?p->kR:p->kH;
    double fC=1.0+p->rC*C2/(kR*kR+C2);
    if((int)p->mech==1) kh=p->kH*(p->dH+p->gH*P)/(p->dH+P);
    else if((int)p->mech>=2){ double x=p->gH*P/(p->dH+P); fP=1.0+x; fPq=1.0+p->qoff*x; }
    double Q=kh*kh/(kh*kh+C2);
    du[2]=(fP*fC*Q*(1.0-H) - fPq*fC*(1.0-Q)*H)/p->epsH;
}
/* one Euler–Maruyama step; sig* = 0 gives plain Euler (used for the deterministic pre-run with a small dt) */
static void em_step(const Par *p, double *u, double dt, int stoch){
    double du[3]; rhs(p,u,du);
    double sq=sqrt(dt);
    u[0]+=du[0]*dt; u[1]+=du[1]*dt; u[2]+=du[2]*dt;
    if(stoch){
        if(p->sigC>0) u[0]+=p->sigC*sqrt(fmax(u[0],0.0))*sq*gauss();
        if(p->sigP>0) u[1]+=p->sigP*fmax(u[1],0.0)*sq*gauss();
        if(p->sigH>0){ double h=fmin(fmax(u[2],0.0),1.0); u[2]+=p->sigH*sqrt(h*(1.0-h))*sq*gauss(); }
    }
    if(u[0]<0)u[0]=-u[0]; if(u[1]<0)u[1]=-u[1]; if(u[2]<0)u[2]=-u[2]; if(u[2]>1)u[2]=2.0-u[2];
}

/* deterministic pre-run: period, baseline, amplitude, H at onset (features of the last cycles) */
typedef struct { int osc; double T, base, amp, Hon, W50; double u[3]; } Det;
static Det det_run(const Par *p){
    Det d; memset(&d,0,sizeof d);
    double u[3]={p->C0,p->P0,p->H0}; double dt=fmin(p->dt,0.02);
    long ntr=(long)(p->Ttrans/dt), nrec=(long)(fmax(4000.0,p->Ttrans)/dt), n;
    for(n=0;n<ntr;n++) em_step(p,u,dt,0);
    double *Cb=malloc(sizeof(double)*nrec), *Hb=malloc(sizeof(double)*nrec);
    for(n=0;n<nrec;n++){ Cb[n]=u[0]; Hb[n]=u[2]; em_step(p,u,dt,0); }
    d.u[0]=u[0]; d.u[1]=u[1]; d.u[2]=u[2];
    double cmax=-1,cmin=1e300; long i; for(i=0;i<nrec;i++){ if(Cb[i]>cmax)cmax=Cb[i]; if(Cb[i]<cmin)cmin=Cb[i]; }
    if(cmax-cmin<1e-4 || cmax-cmin<0.05*cmax){ free(Cb);free(Hb); return d; }
    double mid=0.5*(cmax+cmin); long pk[8]; int np=0;
    for(i=nrec-2;i>0&&np<8;i--){ if(Cb[i]>mid&&Cb[i]>=Cb[i-1]&&Cb[i]>Cb[i+1]){ pk[np++]=i; long j=i; while(j>0&&Cb[j]>mid) j--; i=j; } }
    if(np<3){ free(Cb);free(Hb); return d; }
    d.osc=1; double Ts=0; int k; for(k=0;k<np-1;k++) Ts+=(pk[k]-pk[k+1])*dt; d.T=Ts/(np-1);
    long ip=pk[1]; double base=1e300; for(i=pk[2];i<=ip;i++) if(Cb[i]<base) base=Cb[i];
    d.base=base; d.amp=Cb[ip]-base; double half=base+0.5*d.amp; long a=ip; while(a>0&&Cb[a]>half) a--; long b=ip; while(b<nrec-1&&Cb[b]>half) b++;
    d.W50=(b-a)*dt; d.Hon=Hb[a];
    free(Cb);free(Hb); return d;
}

typedef struct { double t, Hon, amp; } Spk;

int main(int argc,char**argv){
    Par p; defaults(&p);
    char sweep[32]=""; double lo=0,hi=0; int n=0,islog=0,trace=0,isi=0,i;
    for(i=1;i<argc;i++){
        char key[64],val[64]; if(sscanf(argv[i],"%63[^=]=%63s",key,val)!=2){fprintf(stderr,"bad arg %s\n",argv[i]);return 1;}
        if(!strcmp(key,"sweep")) strncpy(sweep,val,31);
        else if(!strcmp(key,"lo")) lo=atof(val); else if(!strcmp(key,"hi")) hi=atof(val);
        else if(!strcmp(key,"n")) n=atoi(val); else if(!strcmp(key,"log")) islog=atoi(val);
        else if(!strcmp(key,"trace")) trace=atoi(val); else if(!strcmp(key,"isi")) isi=atoi(val);
        else { double *f=field(&p,key); if(!f){fprintf(stderr,"unknown key %s\n",key);return 1;} *f=atof(val); }
    }
    double *sw = sweep[0] ? field(&p,sweep) : NULL; if(sweep[0]&&!sw){fprintf(stderr,"unknown sweep key %s\n",sweep);return 1;}
    int nv = sweep[0] ? n : 1;
    if(!trace && !isi) printf("value,n,Tdet,Tav,sd,CV,skew,Hon_det,Hon_mean,Hon_sd,amp_det,amp_mean,amp_cv,corr_isi_amp,corr_isi_next,W50_det\n");
    int iv; for(iv=0;iv<nv;iv++){
        double v=0; if(sw){ v = islog ? lo*pow(hi/lo,(double)iv/(nv-1)) : lo+(hi-lo)*iv/(double)(nv>1?nv-1:1); *sw=v; }
        seed_rng((uint64_t)p.seed + 7919ULL*iv);
        if(p.nuH>0) p.sigH = p.nuH/sqrt(p.epsH);          /* common rescaled noise intensity nu: sigH^2 epsH = nu^2 */
        Det d=det_run(&p);
        if(!d.osc){ if(!trace&&!isi) printf("%.6g,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0\n",v); continue; }
        double thr_hi=d.base+0.5*d.amp, thr_lo=d.base+0.2*d.amp;
        double u[3]={d.u[0],d.u[1],d.u[2]};                     /* start on the deterministic cycle */
        long nsteps=(long)(p.Tend/p.dt), ntr=(long)(p.Ttrans/p.dt), k;
        int cap=100000, ns=0; Spk *sp=malloc(sizeof(Spk)*cap);
        int armed=1; double pkC=-1; int inspike=0;
        if(trace){ printf("t,C,P,H\n"); long nsp=(long)(p.tspan/p.dt), every=(long)fmax(1,0.1/p.dt); for(k=0;k<ntr+nsp;k++){ if(k>=ntr && (k-ntr)%every==0) printf("%.3f,%.6g,%.6g,%.6g\n",(k-ntr)*p.dt,u[0],u[1],u[2]); em_step(&p,u,p.dt,1);} free(sp); return 0; }
        for(k=0;k<nsteps;k++){
            double Hprev=u[2]; em_step(&p,u,p.dt,1);
            if(armed && u[0]>=thr_hi){ armed=0; inspike=1; pkC=u[0]; if(k>=ntr && ns<cap){ sp[ns].t=k*p.dt; sp[ns].Hon=Hprev; sp[ns].amp=0; ns++; } }
            if(inspike){ if(u[0]>pkC) pkC=u[0]; if(u[0]<thr_lo){ inspike=0; armed=1; if(ns>0 && k>=ntr) sp[ns-1].amp=pkC-d.base; } }
        }
        if(isi){ printf("t,ISI,Hon,amp\n"); for(i=1;i<ns;i++) if(sp[i].amp>0) printf("%.3f,%.4f,%.5f,%.5g\n",sp[i].t,sp[i].t-sp[i-1].t,sp[i].Hon,sp[i].amp); free(sp); return 0; }
        /* statistics over ISIs (preceding intervals) with completed spikes */
        int m=0; double *x=malloc(sizeof(double)*ns), *A=malloc(sizeof(double)*ns), *Hn=malloc(sizeof(double)*ns);
        for(i=1;i<ns;i++) if(sp[i].amp>0){ x[m]=sp[i].t-sp[i-1].t; A[m]=sp[i].amp; Hn[m]=sp[i].Hon; m++; }
        double s1=0,s2=0,s3=0,am=0,a2=0,hm=0,h2=0; for(i=0;i<m;i++){ s1+=x[i]; am+=A[i]; hm+=Hn[i]; }
        double Tav=m?s1/m:0, amean=m?am/m:0, hmean=m?hm/m:0;
        for(i=0;i<m;i++){ double e=x[i]-Tav; s2+=e*e; s3+=e*e*e; a2+=(A[i]-amean)*(A[i]-amean); h2+=(Hn[i]-hmean)*(Hn[i]-hmean); }
        double sd=m>1?sqrt(s2/(m-1)):0, skew=(m>2&&s2>0)?(s3/m)/pow(s2/m,1.5):0, asd=m>1?sqrt(a2/(m-1)):0, hsd=m>1?sqrt(h2/(m-1)):0;
        double cxa=0, cxx=0, sxx1=0, sxx2=0; for(i=0;i<m;i++) cxa+=(x[i]-Tav)*(A[i]-amean);
        double r_amp = (s2>0&&a2>0)? cxa/sqrt(s2*a2) : 0;
        /* serial correlation of successive ISIs */
        double mx1=0,mx2=0; for(i=0;i<m-1;i++){ mx1+=x[i]; mx2+=x[i+1]; } if(m>2){ mx1/=(m-1); mx2/=(m-1); }
        for(i=0;i<m-1;i++){ cxx+=(x[i]-mx1)*(x[i+1]-mx2); sxx1+=(x[i]-mx1)*(x[i]-mx1); sxx2+=(x[i+1]-mx2)*(x[i+1]-mx2); }
        double r_next=(sxx1>0&&sxx2>0)?cxx/sqrt(sxx1*sxx2):0;
        printf("%.6g,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.5g,%.5g,%.4f,%.4f,%.4f,%.4f\n",
               v,m,d.T,Tav,sd,Tav>0?sd/Tav:0,skew,d.Hon,hmean,hsd,d.amp,amean,amean>0?asd/amean:0,r_amp,r_next,d.W50);
        fflush(stdout); free(sp); free(x); free(A); free(Hn);
    }
    return 0;
}
