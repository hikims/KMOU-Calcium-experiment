/*
 * m1_pf_sde.c — bootstrap particle filter for the 3-variable M1 SDE (gate 2: synthetic recovery of (eps_H, sigma_H, beta)).
 *
 * State (C, P, H) follows the same SDE as m1_isi_variability/m1_isi_sde.c (deterministic part = stereotyped reference,
 * noise sigH sqrt(H(1-H)) dW_H; sigP, sigC also available).  Observation at times t_k (file, CSV "t,F"):
 *     F_k = b + a C(t_k) + e_k,   e_k ~ N(0, sobs^2)
 * Hidden initial state: each particle starts at a random phase of the deterministic limit cycle (C0,P0,H0 unknown).
 *
 * Modes
 *   ./m1_pf_sde simulate=1 Tobs=2000 dtobs=2 [seed=..] [key=value ...]   -> synthetic record: t,C,P,H,F  (CSV to stdout)
 *   ./m1_pf_sde obs=FILE [npart=1000] [seed=..] [key=value ...]          -> one line: loglik,ess_min,n_obs
 *   ./m1_pf_sde obs=FILE filt=1 ...                                      -> filtered means: t,F,Cf,Pf,Hf  (+ loglik on stderr)
 * Keys: model keys (S alpha beta ... epsH mech gH dH qoff kP kR nA kA Jin rC sigH sigP sigC), a b sobs, dt, npart, seed, Ttrans
 * Systematic resampling at every observation (ESS reported).  Build: gcc -O2 -o m1_pf_sde m1_pf_sde.c -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

typedef struct {
    double S, alpha, beta, gamma, rho, mu, lambda, kS, kC, kH, epsH;
    double mech, gH, dH, qoff, kP, kR, nA, kA, Jin, rC;
    double sigH, sigP, sigC, a, b, sobs, dt, Ttrans, seed, npart, Tobs, dtobs;
} Par;
static void defaults(Par *p){
    p->S=1.0; p->alpha=20; p->beta=4; p->gamma=4.2; p->rho=0.024; p->mu=0.001; p->lambda=0.001; p->kS=0.25; p->kC=0.08; p->kH=0.15; p->epsH=1000;
    p->mech=2; p->gH=30; p->dH=10000; p->qoff=0; p->kP=1; p->kR=1; p->nA=2; p->kA=0.4; p->Jin=0.02; p->rC=300;   /* stereotyped reference */
    p->sigH=0.026; p->sigP=0; p->sigC=0; p->a=0.11; p->b=0.0; p->sobs=0.05; p->dt=0.02; p->Ttrans=3000; p->seed=1; p->npart=1000; p->Tobs=2000; p->dtobs=2;
}
#define KEY(name) if(!strcmp(k,#name)) return &p->name;
static double *field(Par *p,const char*k){
    KEY(S)KEY(alpha)KEY(beta)KEY(gamma)KEY(rho)KEY(mu)KEY(lambda)KEY(kS)KEY(kC)KEY(kH)KEY(epsH)KEY(mech)KEY(gH)KEY(dH)KEY(qoff)KEY(kP)KEY(kR)KEY(nA)KEY(kA)KEY(Jin)KEY(rC)
    KEY(sigH)KEY(sigP)KEY(sigC)KEY(a)KEY(b)KEY(sobs)KEY(dt)KEY(Ttrans)KEY(seed)KEY(npart)KEY(Tobs)KEY(dtobs)
    return NULL;
}
static uint64_t rs[4];
static inline uint64_t rotl(uint64_t x,int k){return (x<<k)|(x>>(64-k));}
static uint64_t next64(void){uint64_t r=rotl(rs[1]*5,7)*9,t=rs[1]<<17;rs[2]^=rs[0];rs[3]^=rs[1];rs[1]^=rs[2];rs[0]^=rs[3];rs[2]^=t;rs[3]=rotl(rs[3],45);return r;}
static void seed_rng(uint64_t s){int i;for(i=0;i<4;i++){s+=0x9E3779B97F4A7C15ULL;uint64_t z=s;z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL;z=(z^(z>>27))*0x94D049BB133111EBULL;rs[i]=z^(z>>31);}}
static double unif(void){return (next64()>>11)*(1.0/9007199254740992.0);}
static double gauss(void){double u1=unif(),u2=unif();if(u1<1e-300)u1=1e-300;return sqrt(-2*log(u1))*cos(2*M_PI*u2);}

static void rhs(const Par*p,const double*u,double*du){
    double C=u[0],P=u[1],H=u[2],C2=C*C, Cn=pow(C,p->nA), act=Cn/(pow(p->kA,p->nA)+Cn);
    du[0]=p->Jin+p->alpha*P/(p->kP+P)*act*H-p->beta*C2/(p->kS*p->kS+C2)-p->mu*C;
    du[1]=p->gamma*p->S/(1+p->S)+p->rho*C/(p->kC+C)-p->lambda*P;
    double kh=p->kH,fP=1,fPq=1,kR=p->kR>0?p->kR:p->kH, fC=1+p->rC*C2/(kR*kR+C2);
    if((int)p->mech==1) kh=p->kH*(p->dH+p->gH*P)/(p->dH+P); else if((int)p->mech>=2){double x=p->gH*P/(p->dH+P);fP=1+x;fPq=1+p->qoff*x;}
    double Q=kh*kh/(kh*kh+C2); du[2]=(fP*fC*Q*(1-H)-fPq*fC*(1-Q)*H)/p->epsH;
}
static void em_step(const Par*p,double*u,double dt,int stoch){
    double du[3]; rhs(p,u,du); double sq=sqrt(dt);
    u[0]+=du[0]*dt; u[1]+=du[1]*dt; u[2]+=du[2]*dt;
    if(stoch){ if(p->sigC>0)u[0]+=p->sigC*sqrt(fmax(u[0],0))*sq*gauss(); if(p->sigP>0)u[1]+=p->sigP*fmax(u[1],0)*sq*gauss();
               if(p->sigH>0){double h=fmin(fmax(u[2],0),1); u[2]+=p->sigH*sqrt(h*(1-h))*sq*gauss();} }
    if(u[0]<0)u[0]=-u[0]; if(u[1]<0)u[1]=-u[1]; if(u[2]<0)u[2]=-u[2]; if(u[2]>1)u[2]=2-u[2];
}
/* deterministic limit cycle sampled at NCYC points (for the initial-state prior) */
#define NCYC 400
static int cycle_states(const Par*p,double cyc[NCYC][3]){
    double u[3]={0.1,0.05,1.0}; long n,ntr=(long)(p->Ttrans/p->dt); for(n=0;n<ntr;n++) em_step(p,u,p->dt,0);
    /* find the period from upward crossings of C at a mid level */
    long nrec=(long)(4000.0/p->dt); double *Cb=malloc(sizeof(double)*nrec); double (*U)[3]=malloc(sizeof(double)*3*nrec);
    for(n=0;n<nrec;n++){Cb[n]=u[0];U[n][0]=u[0];U[n][1]=u[1];U[n][2]=u[2];em_step(p,u,p->dt,0);}
    double cmax=-1,cmin=1e300; for(n=0;n<nrec;n++){if(Cb[n]>cmax)cmax=Cb[n];if(Cb[n]<cmin)cmin=Cb[n];}
    double mid=0.5*(cmax+cmin); long c1=-1,c2=-1;
    for(n=nrec-2;n>0;n--){ if(Cb[n]<mid&&Cb[n+1]>=mid){ if(c2<0)c2=n; else {c1=n;break;} } }
    if(c1<0||cmax-cmin<0.05*cmax){ free(Cb);free(U); return 0; }
    int i; for(i=0;i<NCYC;i++){ long j=c1+(long)((c2-c1)*(double)i/NCYC); cyc[i][0]=U[j][0];cyc[i][1]=U[j][1];cyc[i][2]=U[j][2]; }
    free(Cb);free(U); return 1;
}

int main(int argc,char**argv){
    Par p; defaults(&p); char obs[256]=""; int sim=0,filt=0,i;
    for(i=1;i<argc;i++){ char key[64],val[200]; if(sscanf(argv[i],"%63[^=]=%199s",key,val)!=2){fprintf(stderr,"bad arg %s\n",argv[i]);return 1;}
        if(!strcmp(key,"obs"))strncpy(obs,val,255); else if(!strcmp(key,"simulate"))sim=atoi(val); else if(!strcmp(key,"filt"))filt=atoi(val);
        else {double*f=field(&p,key); if(!f){fprintf(stderr,"unknown key %s\n",key);return 1;} *f=atof(val);} }
    seed_rng((uint64_t)p.seed);
    static double cyc[NCYC][3]; if(!cycle_states(&p,cyc)){ if(sim){fprintf(stderr,"no oscillation\n");return 1;} printf("%.6g,%.3f,%d\n",-1e300,0.0,0); return 0; }
    if(sim){
        double u[3]; int k0=(int)(unif()*NCYC); u[0]=cyc[k0][0];u[1]=cyc[k0][1];u[2]=cyc[k0][2];
        long nper=(long)(p.dtobs/p.dt+0.5), nobs=(long)(p.Tobs/p.dtobs); printf("t,C,P,H,F\n");
        for(long k=0;k<=nobs;k++){ printf("%.3f,%.6g,%.6g,%.6g,%.6g\n",k*p.dtobs,u[0],u[1],u[2],p.b+p.a*u[0]+p.sobs*gauss()); for(long n=0;n<nper;n++) em_step(&p,u,p.dt,1); }
        return 0;
    }
    /* read observations */
    FILE*fp=fopen(obs,"r"); if(!fp){fprintf(stderr,"cannot open %s\n",obs);return 1;}
    double *to=malloc(sizeof(double)*200000),*Fo=malloc(sizeof(double)*200000); long nobs=0; char line[1024];
    if(!fgets(line,1024,fp)){return 1;}
    while(fgets(line,1024,fp)){ double t,F; char *q=line; int col=0; double vals[8]={0}; char*tok=strtok(q,",\n"); while(tok&&col<8){vals[col++]=atof(tok);tok=strtok(NULL,",\n");}
        t=vals[0]; F=vals[col-1]; to[nobs]=t; Fo[nobs]=F; nobs++; }
    fclose(fp);
    int N=(int)p.npart; double (*X)[3]=malloc(sizeof(double)*3*N),(*Y)[3]=malloc(sizeof(double)*3*N); double *w=malloc(sizeof(double)*N),*cw=malloc(sizeof(double)*N);
    for(i=0;i<N;i++){int k0=(int)(unif()*NCYC); X[i][0]=cyc[k0][0];X[i][1]=cyc[k0][1];X[i][2]=cyc[k0][2];}
    double loglik=0, essmin=1e300; if(filt) printf("t,F,Cf,Pf,Hf\n");
    for(long k=0;k<nobs;k++){
        if(k>0){ double span=to[k]-to[k-1]; long nst=(long)(span/p.dt+0.5); double dts=span/nst; for(i=0;i<N;i++) for(long n=0;n<nst;n++) em_step(&p,X[i],dts,1); }
        double mx=-1e300; for(i=0;i<N;i++){ double r=(Fo[k]-p.b-p.a*X[i][0])/p.sobs; w[i]=-0.5*r*r; if(w[i]>mx)mx=w[i]; }
        double sw=0,sw2=0; for(i=0;i<N;i++){ w[i]=exp(w[i]-mx); sw+=w[i]; sw2+=w[i]*w[i]; }
        loglik+=mx+log(sw/N)-0.5*log(2*M_PI)-log(p.sobs); double ess=sw*sw/sw2; if(ess<essmin)essmin=ess;
        if(filt){ double cm=0,pm=0,hm=0; for(i=0;i<N;i++){cm+=w[i]*X[i][0];pm+=w[i]*X[i][1];hm+=w[i]*X[i][2];} printf("%.3f,%.6g,%.6g,%.6g,%.6g\n",to[k],Fo[k],cm/sw,pm/sw,hm/sw); }
        /* systematic resampling */
        double c=0; for(i=0;i<N;i++){c+=w[i]/sw; cw[i]=c;} double u0=unif()/N; int j=0;
        for(i=0;i<N;i++){ double uu=u0+(double)i/N; while(j<N-1&&cw[j]<uu) j++; Y[i][0]=X[j][0];Y[i][1]=X[j][1];Y[i][2]=X[j][2]; }
        memcpy(X,Y,sizeof(double)*3*N);
    }
    if(filt) fprintf(stderr,"loglik %.6g ess_min %.1f\n",loglik,essmin); else printf("%.6g,%.3f,%ld\n",loglik,essmin,nobs);
    return 0;
}
