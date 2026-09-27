/*
 * m0ip3h_waveform_sweep.c — M0 + IP3->H rate coupling + cooperative Ca activation: can the spike be narrow?
 *
 * Base model (AIMS Math 2026, eqs. 2.6–2.8 without diffusion):
 *   dC/dt = Jin + alpha * P/(1+P) * C^nA/(kA^nA + C^nA) * H - beta * C^2/(kS^2 + C^2) - mu*C
 *           (Jin = basal Ca influx/leak, default 0 = original)
 *           (nA=1, kA=1 reproduces the original AIMS activation C/(1+C); nA>=2 = cooperative CICR)
 *   dP/dt = gamma * S/(1+S) + rho * C/(kC + C) - lambda*P
 *   dH/dt = (1/epsH) * ( Q(C,P) - H )                      [mech=0: Q = kH^2/(kH^2+C^2), original]
 *
 * Ca-dependent speed-up of H kinetics (Li–Rinzel-like tau_h(C); rC=0 = original):
 *   dH/dt = (1/epsH) * (1 + rC * C^2/(kH^2+C^2)) * [IP3 factor] * (Q(C,P) - H)
 *
 * IP3 -> H coupling (mech, gH, dH; gH at its neutral value reproduces the original model exactly):
 *   mech=1 (threshold, Li–Rinzel-like): kH(P) = kH * (dH + gH*P)/(dH + P)
 *           -> IP3 raises the Ca-inactivation threshold; neutral gH=1
 *   mech=2 (rate):  dH/dt = (1 + gH*P/(dH+P)) * (Q(C) - H)/epsH
 *           -> IP3 accelerates receptor recovery/inactivation; neutral gH=0
 *   dH = half-saturation of P for the coupling (P is O(1000) at lambda=0.001, O(1) at lambda=1)
 *
 * Usage:
 *   ./m0ip3h_waveform_sweep [key=value ...]                 -> one run, prints features (one CSV line)
 *   ./m0ip3h_waveform_sweep sweep=NAME lo=A hi=B n=N [log=1] [key=value ...]
 *                                            -> sweep NAME from A to B (N points), CSV to stdout
 *   ./m0ip3h_waveform_sweep trace=1 [key=value ...]         -> print t, C, P, H of the last cycles to stdout
 *
 * Keys: S alpha beta gamma rho mu lambda kS kC kH epsH  mech gH dH  nA kA Jin rC  C0 P0 H0  dt Tend Ttrans
 *
 * Output columns (features of the last full cycle after transient):
 *   value, osc(0/1), T, Cmax, Cmin, W50, rise50, decay50, duty, Pmean, Hmin, Hmax
 *   T      = period (mean of last NPEAK-1 peak intervals)
 *   W50    = full width at half maximum of the C spike (relative to Cmin)
 *   rise50 = time from upward half-crossing to peak
 *   decay50= time from peak to downward half-crossing
 *   duty   = W50 / Tosc
 *
 * Build: gcc -O2 -o m0ip3h_stimulus_sweep m0ip3h_waveform_sweep.c -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct {
    double S, alpha, beta, gamma, rho, mu, lambda, kS, kC, kH, epsH;
    double mech, gH, dH;
    double nA, kA, Jin, rC;
    double C0, P0, H0;
    double dt, Tend, Ttrans;
} Par;

static void defaults(Par *p) {
    /* Table 1 of the AIMS paper (baseline). S = constant stimulus level. */
    p->S = 0.0;  p->alpha = 3.50; p->beta = 0.26; p->gamma = 4.20; p->rho = 2.40;
    p->mu = 0.001; p->lambda = 0.001; p->kS = 0.20; p->kC = 0.08; p->kH = 0.15;
    p->epsH = 70.0;
    p->mech = 0; p->gH = 1.0; p->dH = 1.0;   /* mech=0: original AIMS model */
    p->nA = 1.0; p->kA = 1.0; p->Jin = 0.0; /* nA=1,kA=1,Jin=0: original C equation */
    p->rC = 0.0;                            /* rC=0: original H kinetics */
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
    double Cn = pow(C, p->nA), act = Cn/(pow(p->kA, p->nA) + Cn);   /* cooperative Ca activation */
    du[0] = p->Jin + p->alpha * P/(1.0+P) * act * H - p->beta * C2/(p->kS*p->kS + C2) - p->mu*C;
    du[1] = p->gamma * p->S/(1.0+p->S) + p->rho * C/(p->kC + C) - p->lambda*P;
    double kh = p->kH, rate = 1.0/p->epsH;
    if ((int)p->mech == 1) kh = p->kH * (p->dH + p->gH*P)/(p->dH + P);          /* threshold coupling */
    else if ((int)p->mech == 2) rate = (1.0 + p->gH*P/(p->dH + P)) / p->epsH;  /* rate coupling */
    rate *= (1.0 + p->rC * C2/(p->kH*p->kH + C2));                            /* Ca speeds up H kinetics */
    du[2] = rate * (kh*kh/(kh*kh + C2) - H);
}

static void rk4(const Par *p, double *u, double dt) {
    double k1[3],k2[3],k3[3],k4[3],v[3]; int i;
    rhs(p,u,k1);
    for(i=0;i<3;i++) v[i]=u[i]+0.5*dt*k1[i]; rhs(p,v,k2);
    for(i=0;i<3;i++) v[i]=u[i]+0.5*dt*k2[i]; rhs(p,v,k3);
    for(i=0;i<3;i++) v[i]=u[i]+dt*k3[i];     rhs(p,v,k4);
    for(i=0;i<3;i++) u[i]+=dt/6.0*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
    /* keep the biologically admissible set (numerical safety only) */
    if(u[0]<0) u[0]=0; if(u[1]<0) u[1]=0; if(u[2]<0) u[2]=0; if(u[2]>1) u[2]=1;
}

#define NPEAK 6

typedef struct { int osc; double T,Cmax,Cmin,W50,rise50,decay50,duty,Pmean,Hmin,Hmax; } Feat;

/* Integrate; after Ttrans record C(t) in a ring buffer covering the last 'keep' time units,
   then extract the last NPEAK peaks and the features of the last full cycle. */
static Feat run(const Par *p, int trace) {
    Feat f; memset(&f,0,sizeof f);
    double u[3]={p->C0,p->P0,p->H0};
    long nsteps=(long)(p->Tend/p->dt+0.5);
    long ntr=(long)(p->Ttrans/p->dt+0.5);
    long nrec=nsteps-ntr; if(nrec<10){fprintf(stderr,"Tend must exceed Ttrans\n"); exit(1);}
    double *Cb=malloc(sizeof(double)*nrec), *Pb=malloc(sizeof(double)*nrec), *Hb=malloc(sizeof(double)*nrec);
    long n;
    for(n=0;n<nsteps;n++){
        if(n>=ntr){ Cb[n-ntr]=u[0]; Pb[n-ntr]=u[1]; Hb[n-ntr]=u[2]; }
        rk4(p,u,p->dt);
    }
    /* global stats over recorded window */
    double cmax=-1,cmin=1e300,hmin=1e300,hmax=-1,psum=0; long i;
    for(i=0;i<nrec;i++){ if(Cb[i]>cmax)cmax=Cb[i]; if(Cb[i]<cmin)cmin=Cb[i]; if(Hb[i]<hmin)hmin=Hb[i]; if(Hb[i]>hmax)hmax=Hb[i]; psum+=Pb[i]; }
    f.Cmax=cmax; f.Cmin=cmin; f.Hmin=hmin; f.Hmax=hmax; f.Pmean=psum/nrec;
    /* require a finite-amplitude oscillation (relative amplitude > 5%) */
    if(cmax-cmin<1e-4 || (cmax-cmin)<0.05*cmax){ f.osc=0; free(Cb);free(Pb);free(Hb); return f; }
    /* find local maxima above midline (prominence filter) */
    double mid=0.5*(cmax+cmin);
    long pk[NPEAK+1]; int np=0;
    for(i=nrec-2;i>0 && np<NPEAK;i--){
        if(Cb[i]>mid && Cb[i]>=Cb[i-1] && Cb[i]>Cb[i+1]){
            /* skip the very last peak if the cycle after it is incomplete (need a downward crossing) */
            pk[np++]=i;
            /* jump back at least a quarter of an estimated spike width to avoid double counting */
            long j=i; while(j>0 && Cb[j]>mid) j--; i=j;
        }
    }
    if(np<3){ f.osc=0; free(Cb);free(Pb);free(Hb); return f; }
    f.osc=1;
    /* period: mean interval between successive peaks (pk is in reverse time order) */
    double Tsum=0; int k; for(k=0;k<np-1;k++) Tsum+=(pk[k]-pk[k+1])*p->dt; f.T=Tsum/(np-1);
    /* features on the second-to-last peak (guaranteed a complete cycle around it) */
    long ip=pk[1];
    /* baseline = local minimum between the previous peak and this peak (inter-spike rest level) */
    double base=1e300; for(i=pk[2];i<=ip;i++) if(Cb[i]<base) base=Cb[i];
    double half=base+0.5*(Cb[ip]-base);
    long a=ip; while(a>0 && Cb[a]>half) a--;
    long b=ip; while(b<nrec-1 && Cb[b]>half) b++;
    f.W50=(b-a)*p->dt; f.rise50=(ip-a)*p->dt; f.decay50=(b-ip)*p->dt; f.duty=f.W50/f.T;
    if(trace){
        long start=pk[np-1]-(long)(0.1*f.T/p->dt); if(start<0)start=0;
        printf("t,C,P,H\n");
        for(i=start;i<nrec;i+=(long)fmax(1,0.02/p->dt)) printf("%.3f,%.6g,%.6g,%.6g\n",p->Ttrans+i*p->dt,Cb[i],Pb[i],Hb[i]);
    }
    free(Cb);free(Pb);free(Hb); return f;
}

static void print_feat(double val, const Feat *f){
    printf("%.6g,%d,%.4f,%.5g,%.5g,%.4f,%.4f,%.4f,%.4f,%.5g,%.4f,%.4f\n",
        val,f->osc,f->T,f->Cmax,f->Cmin,f->W50,f->rise50,f->decay50,f->duty,f->Pmean,f->Hmin,f->Hmax);
}

int main(int argc,char**argv){
    Par p; defaults(&p);
    char sweep[32]=""; double lo=0,hi=0; int n=0, islog=0, trace=0; int i;
    for(i=1;i<argc;i++){
        char key[64],val[64]; if(sscanf(argv[i],"%63[^=]=%63s",key,val)!=2){fprintf(stderr,"bad arg %s\n",argv[i]);return 1;}
        if(!strcmp(key,"sweep")) strncpy(sweep,val,31);
        else if(!strcmp(key,"lo")) lo=atof(val);
        else if(!strcmp(key,"hi")) hi=atof(val);
        else if(!strcmp(key,"n")) n=atoi(val);
        else if(!strcmp(key,"log")) islog=atoi(val);
        else if(!strcmp(key,"trace")) trace=atoi(val);
        else { double *f=field(&p,key); if(!f){fprintf(stderr,"unknown key %s\n",key);return 1;} *f=atof(val); }
    }
    if(trace){ run(&p,1); return 0; }
    printf("value,osc,Tosc,Cmax,Cmin,W50,rise50,decay50,duty,Pmean,Hmin,Hmax\n");
    if(sweep[0]){
        double *f=field(&p,sweep); if(!f){fprintf(stderr,"unknown sweep key %s\n",sweep);return 1;}
        for(i=0;i<n;i++){
            double v = islog ? lo*pow(hi/lo,(double)i/(n-1)) : lo+(hi-lo)*i/(double)(n-1);
            *f=v; Feat ft=run(&p,0); print_feat(v,&ft); fflush(stdout);
        }
    } else { Feat ft=run(&p,0); print_feat(0,&ft); }
    return 0;
}
