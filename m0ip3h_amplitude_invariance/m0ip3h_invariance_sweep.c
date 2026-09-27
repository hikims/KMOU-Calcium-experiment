/*
 * m0ip3h_invariance_sweep.c — spike stereotypy: which structure keeps the spike (amplitude, width)
 *                             fixed while the stimulus changes the period?
 *
 * Base model (AIMS Math 2026, eqs. 2.6–2.8 without diffusion) + the ingredients fixed in the previous steps:
 *   dC/dt = Jin + alpha * P/(kP+P) * C^nA/(kA^nA + C^nA) * H - beta * C^2/(kS^2 + C^2) - mu*C
 *   dP/dt = gamma * S/(1+S) + rho * C/(kC + C) - lambda*P
 *   H : receptor availability, written as a two-rate gate (identical to the original (Q - H)/epsH form)
 *
 *       dH/dt = (1/epsH) * [ fP * fC * Q(C) * (1 - H)   -   fPq * fC * (1 - Q(C)) * H ]
 *                             ^ recovery (on)              ^ Ca-inactivation (off)
 *       Q(C) = kH^2/(kH^2 + C^2)
 *       fP   = 1 + gH * P/(dH + P)          IP3 accelerates recovery              (gH = 0: none)
 *       fPq  = 1 + qoff * gH * P/(dH + P)   fraction qoff of it also accelerates inactivation
 *       fC   = 1 + rC * C^2/(kR^2 + C^2)    Ca accelerates the H kinetics          (rC = 0: none)
 *              (kR < 0, default: kR = kH as in the previous steps; kR >> C_rest switches the speed-up off at rest)
 *
 *   qoff = 1 : IP3 rescales the whole H clock  = mech=2 of the previous steps (period, width, amplitude all scale)
 *   qoff = 0 : IP3 acts on recovery only       -> the spike itself is P-independent (stereotyped)
 *   kP       : half-saturation of the IP3 gain of release (kP = 1 original; P is O(10^3) at lambda = 0.001,
 *              so kP ~ 10^3 makes the release gain dose-dependent, which can compensate a shorter spike)
 *   mech=1 (threshold coupling) is kept for completeness: kH(P) = kH (dH + gH P)/(dH + P), fP = fPq = 1.
 *
 * Usage:
 *   ./m0ip3h_invariance_sweep [key=value ...]                       one run, features (one CSV line)
 *   ./m0ip3h_invariance_sweep sweep=NAME lo=A hi=B n=N [log=1] ...  sweep NAME, CSV to stdout
 *   ./m0ip3h_invariance_sweep trace=1 ...                           t, C, P, H of the last cycles
 *   ./m0ip3h_invariance_sweep recovery=1 [ntau=N pulse=X evth=Y] ...
 *        recovery function: after a free spike, a spike is evoked by an instantaneous Ca pulse (C += pulse)
 *        at N delays tau spanning the free inter-spike interval; prints
 *        tau, tau_over_T, tau_over_wait, H, evoked(0/1), A, A_over_Afree, W50, W50_over_Wfree, t_peak, t_free
 *        (evoked: A > evth*Afree; the response is measured in a window that ends before the free spike would arrive;
 *        (A = evoked peak minus the free baseline; Afree = amplitude of the free spike;
 *         wait = time from the end of the free spike to the next free spike onset)
 *
 * Keys: S alpha beta gamma rho mu lambda kS kC kH epsH  mech gH dH qoff kP  nA kA Jin rC kR  C0 P0 H0  dt Tend Ttrans
 *
 * Output columns (sweep / single run; features of the last full cycle after the transient):
 *   value, osc, Tosc, Cmax, Cmin, W50, rise50, decay50, duty, Pmean, Hmin, Hmax, amp, Hon
 *   amp = Cmax - Cmin (the quantity a non-saturating indicator reports),  Hon = H at the upward half-crossing
 *
 * Build: gcc -O2 -o m0ip3h_invariance_sweep m0ip3h_invariance_sweep.c -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct {
    double S, alpha, beta, gamma, rho, mu, lambda, kS, kC, kH, epsH;
    double mech, gH, dH, qoff, kP, kR;
    double nA, kA, Jin, rC;
    double C0, P0, H0;
    double dt, Tend, Ttrans;
} Par;

static void defaults(Par *p) {
    p->S = 0.0;  p->alpha = 3.50; p->beta = 0.26; p->gamma = 4.20; p->rho = 2.40;
    p->mu = 0.001; p->lambda = 0.001; p->kS = 0.20; p->kC = 0.08; p->kH = 0.15;
    p->epsH = 70.0;
    p->mech = 0; p->gH = 0.0; p->dH = 1.0; p->qoff = 1.0; p->kP = 1.0; p->kR = -1.0;   /* original AIMS model; kR<0: kR = kH */
    p->nA = 1.0; p->kA = 1.0; p->Jin = 0.0; p->rC = 0.0;
    p->C0 = 0.10; p->P0 = 0.05; p->H0 = 1.0;
    p->dt = 0.01; p->Tend = 6000.0; p->Ttrans = 3000.0;
}

static double *field(Par *p, const char *k) {
    if (!strcmp(k,"S")) return &p->S;
    if (!strcmp(k,"alpha")) return &p->alpha;
    if (!strcmp(k,"beta")) return &p->beta;
    if (!strcmp(k,"gamma")) return &p->gamma;
    if (!strcmp(k,"rho")) return &p->rho;
    if (!strcmp(k,"mu")) return &p->mu;
    if (!strcmp(k,"lambda")) return &p->lambda;
    if (!strcmp(k,"kS")) return &p->kS;
    if (!strcmp(k,"kC")) return &p->kC;
    if (!strcmp(k,"kH")) return &p->kH;
    if (!strcmp(k,"epsH")) return &p->epsH;
    if (!strcmp(k,"mech")) return &p->mech;
    if (!strcmp(k,"gH")) return &p->gH;
    if (!strcmp(k,"dH")) return &p->dH;
    if (!strcmp(k,"qoff")) return &p->qoff;
    if (!strcmp(k,"kP")) return &p->kP;
    if (!strcmp(k,"kR")) return &p->kR;
    if (!strcmp(k,"nA")) return &p->nA;
    if (!strcmp(k,"kA")) return &p->kA;
    if (!strcmp(k,"Jin")) return &p->Jin;
    if (!strcmp(k,"rC")) return &p->rC;
    if (!strcmp(k,"C0")) return &p->C0;
    if (!strcmp(k,"P0")) return &p->P0;
    if (!strcmp(k,"H0")) return &p->H0;
    if (!strcmp(k,"dt")) return &p->dt;
    if (!strcmp(k,"Tend")) return &p->Tend;
    if (!strcmp(k,"Ttrans")) return &p->Ttrans;
    return NULL;
}

static void rhs(const Par *p, const double *u, double *du) {
    double C = u[0], P = u[1], H = u[2];
    double C2 = C*C;
    double Cn = pow(C, p->nA), act = Cn/(pow(p->kA, p->nA) + Cn);
    du[0] = p->Jin + p->alpha * P/(p->kP+P) * act * H - p->beta * C2/(p->kS*p->kS + C2) - p->mu*C;
    du[1] = p->gamma * p->S/(1.0+p->S) + p->rho * C/(p->kC + C) - p->lambda*P;
    double kh = p->kH, fP = 1.0, fPq = 1.0;
    double kR = p->kR > 0 ? p->kR : p->kH;
    double fC = 1.0 + p->rC * C2/(kR*kR + C2);
    if ((int)p->mech == 1) kh = p->kH * (p->dH + p->gH*P)/(p->dH + P);
    else if ((int)p->mech >= 2) { double x = p->gH * P/(p->dH + P); fP = 1.0 + x; fPq = 1.0 + p->qoff * x; }
    double Q = kh*kh/(kh*kh + C2);
    du[2] = ( fP*fC*Q*(1.0 - H) - fPq*fC*(1.0 - Q)*H ) / p->epsH;
}

static void rk4(const Par *p, double *u, double dt) {
    double k1[3],k2[3],k3[3],k4[3],v[3]; int i;
    rhs(p,u,k1);
    for(i=0;i<3;i++) v[i]=u[i]+0.5*dt*k1[i]; rhs(p,v,k2);
    for(i=0;i<3;i++) v[i]=u[i]+0.5*dt*k2[i]; rhs(p,v,k3);
    for(i=0;i<3;i++) v[i]=u[i]+dt*k3[i];     rhs(p,v,k4);
    for(i=0;i<3;i++) u[i]+=dt/6.0*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
    if(u[0]<0) u[0]=0; if(u[1]<0) u[1]=0; if(u[2]<0) u[2]=0; if(u[2]>1) u[2]=1;
}

#define NPEAK 6
typedef struct { int osc; double T,Cmax,Cmin,W50,rise50,decay50,duty,Pmean,Hmin,Hmax,amp,Hon;
                 double base, tend_spike, wait, ufin[3]; } Feat;

/* half-width of a spike in array y around index ip, baseline base; returns indices a,b of the half crossings */
static double halfwidth(const double *y, long n, long ip, double base, double dt, long *ia, long *ib){
    double half=base+0.5*(y[ip]-base); long a=ip, b=ip;
    while(a>0 && y[a]>half) a--;
    while(b<n-1 && y[b]>half) b++;
    if(ia)*ia=a; if(ib)*ib=b; return (b-a)*dt;
}

/* Integrate; record C,P,H after Ttrans; extract the last NPEAK peaks and the features of one full cycle.
   Also returns (for the recovery protocol) the state at the end of the analysed spike and the free waiting time. */
static Feat run(const Par *p, int trace) {
    Feat f; memset(&f,0,sizeof f);
    double u[3]={p->C0,p->P0,p->H0};
    long nsteps=(long)(p->Tend/p->dt+0.5), ntr=(long)(p->Ttrans/p->dt+0.5);
    long nrec=nsteps-ntr; if(nrec<10){fprintf(stderr,"Tend must exceed Ttrans\n"); exit(1);}
    double *Cb=malloc(sizeof(double)*nrec), *Pb=malloc(sizeof(double)*nrec), *Hb=malloc(sizeof(double)*nrec);
    long n,i;
    for(n=0;n<nsteps;n++){
        if(n>=ntr){ Cb[n-ntr]=u[0]; Pb[n-ntr]=u[1]; Hb[n-ntr]=u[2]; }
        rk4(p,u,p->dt);
    }
    double cmax=-1,cmin=1e300,hmin=1e300,hmax=-1,psum=0;
    for(i=0;i<nrec;i++){ if(Cb[i]>cmax)cmax=Cb[i]; if(Cb[i]<cmin)cmin=Cb[i]; if(Hb[i]<hmin)hmin=Hb[i]; if(Hb[i]>hmax)hmax=Hb[i]; psum+=Pb[i]; }
    f.Cmax=cmax; f.Cmin=cmin; f.Hmin=hmin; f.Hmax=hmax; f.Pmean=psum/nrec; f.amp=cmax-cmin;
    if(cmax-cmin<1e-4 || (cmax-cmin)<0.05*cmax){ f.osc=0; free(Cb);free(Pb);free(Hb); return f; }
    double mid=0.5*(cmax+cmin);
    long pk[NPEAK+1]; int np=0;
    for(i=nrec-2;i>0 && np<NPEAK;i--){
        if(Cb[i]>mid && Cb[i]>=Cb[i-1] && Cb[i]>Cb[i+1]){
            pk[np++]=i; long j=i; while(j>0 && Cb[j]>mid) j--; i=j;
        }
    }
    if(np<3){ f.osc=0; free(Cb);free(Pb);free(Hb); return f; }
    f.osc=1;
    double Tsum=0; int k; for(k=0;k<np-1;k++) Tsum+=(pk[k]-pk[k+1])*p->dt; f.T=Tsum/(np-1);
    long ip=pk[1];                                   /* second-to-last peak: complete cycle on both sides */
    double base=1e300; for(i=pk[2];i<=ip;i++) if(Cb[i]<base) base=Cb[i];
    long a,b; f.W50=halfwidth(Cb,nrec,ip,base,p->dt,&a,&b);
    f.rise50=(ip-a)*p->dt; f.decay50=(b-ip)*p->dt; f.duty=f.W50/f.T; f.Hon=Hb[a]; f.base=base;
    /* end of this spike: first index after the peak where C is back within 10% of the amplitude above baseline */
    long e=ip; while(e<nrec-1 && Cb[e]>base+0.1*(Cb[ip]-base)) e++;
    /* next free onset: next upward half-crossing after e */
    double half=base+0.5*(Cb[ip]-base); long o=e; while(o<nrec-1 && Cb[o]<half) o++;
    f.tend_spike=e*p->dt; f.wait=(o-e)*p->dt; f.ufin[0]=Cb[e]; f.ufin[1]=Pb[e]; f.ufin[2]=Hb[e];
    if(trace){
        long start=pk[np-1]-(long)(0.1*f.T/p->dt); if(start<0)start=0;
        printf("t,C,P,H\n");
        for(i=start;i<nrec;i+=(long)fmax(1,0.02/p->dt)) printf("%.3f,%.6g,%.6g,%.6g\n",p->Ttrans+i*p->dt,Cb[i],Pb[i],Hb[i]);
    }
    free(Cb);free(Pb);free(Hb); return f;
}

/* Recovery protocol: from the state at the end of a free spike, wait tau, add 'pulse' to C, and measure the evoked spike. */
static void recovery(const Par *p, int ntau, double pulse, double evth){
    Feat f=run(p,0);
    if(!f.osc){ fprintf(stderr,"no free oscillation at these parameters\n"); return; }
    printf("tau,tau_over_T,tau_over_wait,H,evoked,A,A_over_Afree,W50,W50_over_Wfree,t_peak,t_free\n");
    int k; for(k=0;k<ntau;k++){
        double tau = f.wait * (0.02 + 0.94*(double)k/(ntau-1));   /* 2%..96% of the free waiting time */
        double u[3]={f.ufin[0],f.ufin[1],f.ufin[2]};
        long nw=(long)(tau/p->dt+0.5), n;
        for(n=0;n<nw;n++) rk4(p,u,p->dt);
        double Hpre=u[2];
        u[0]+=pulse;
        /* observation window: spike-scale, and ending before the free spike would arrive */
        double tfree=f.wait-tau, win=fmin(4.0*f.W50 + 2.0*f.rise50, 0.6*tfree + 2.0*f.rise50 + 0.5*f.W50);
        long nobs=(long)(win/p->dt)+1; double *y=malloc(sizeof(double)*nobs);
        long ipk=0; for(n=0;n<nobs;n++){ y[n]=u[0]; if(u[0]>y[ipk]) ipk=n; rk4(p,u,p->dt); }
        double A=y[ipk]-f.base, tpk=ipk*p->dt;
        int ev = (A > evth*f.amp);
        double w = ev ? halfwidth(y,nobs,ipk,f.base,p->dt,NULL,NULL) : 0.0;
        printf("%.4f,%.4f,%.4f,%.5f,%d,%.5g,%.5f,%.4f,%.5f,%.4f,%.4f\n", tau, tau/f.T, tau/f.wait, Hpre, ev, A, A/f.amp, w, w/f.W50, tpk, tfree);
        free(y);
    }
}

static void print_feat(double val, const Feat *f){
    printf("%.6g,%d,%.4f,%.5g,%.5g,%.4f,%.4f,%.4f,%.4f,%.5g,%.4f,%.4f,%.5g,%.4f\n",
        val,f->osc,f->T,f->Cmax,f->Cmin,f->W50,f->rise50,f->decay50,f->duty,f->Pmean,f->Hmin,f->Hmax,f->amp,f->Hon);
}

int main(int argc,char**argv){
    Par p; defaults(&p);
    char sweep[32]=""; double lo=0,hi=0,pulse=1.0,evth=0.5; int n=0, islog=0, trace=0, rec=0, ntau=25; int i;
    for(i=1;i<argc;i++){
        char key[64],val[64]; if(sscanf(argv[i],"%63[^=]=%63s",key,val)!=2){fprintf(stderr,"bad arg %s\n",argv[i]);return 1;}
        if(!strcmp(key,"sweep")) strncpy(sweep,val,31);
        else if(!strcmp(key,"lo")) lo=atof(val);
        else if(!strcmp(key,"hi")) hi=atof(val);
        else if(!strcmp(key,"n")) n=atoi(val);
        else if(!strcmp(key,"log")) islog=atoi(val);
        else if(!strcmp(key,"trace")) trace=atoi(val);
        else if(!strcmp(key,"recovery")) rec=atoi(val);
        else if(!strcmp(key,"ntau")) ntau=atoi(val);
        else if(!strcmp(key,"pulse")) pulse=atof(val);
        else if(!strcmp(key,"evth")) evth=atof(val);
        else { double *f=field(&p,key); if(!f){fprintf(stderr,"unknown key %s\n",key);return 1;} *f=atof(val); }
    }
    if(trace){ run(&p,1); return 0; }
    if(rec){ recovery(&p,ntau,pulse,evth); return 0; }
    printf("value,osc,Tosc,Cmax,Cmin,W50,rise50,decay50,duty,Pmean,Hmin,Hmax,amp,Hon\n");
    if(sweep[0]){
        double *f=field(&p,sweep); if(!f){fprintf(stderr,"unknown sweep key %s\n",sweep);return 1;}
        for(i=0;i<n;i++){
            double v = islog ? lo*pow(hi/lo,(double)i/(n-1)) : lo+(hi-lo)*i/(double)(n-1);
            *f=v; Feat ft=run(&p,0); print_feat(v,&ft); fflush(stdout);
        }
    } else { Feat ft=run(&p,0); print_feat(0,&ft); }
    return 0;
}
