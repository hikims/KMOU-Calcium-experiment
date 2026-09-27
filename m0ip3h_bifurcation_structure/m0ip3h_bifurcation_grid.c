/*
 * m0ip3h_bifurcation_grid.c — oscillation region, fixed points, Jacobian eigenvalues and cycle features of the
 *                             modified deterministic model on parameter grids (OpenMP over grid cells).
 *
 * Model (all switches of m0ip3h_amplitude_invariance; defaults = original AIMS system):
 *   dC/dt = Jin + alpha P/(kP+P) C^nA/(kA^nA+C^nA) H − beta C²/(kS²+C²) − mu C
 *   dP/dt = gamma S/(1+S) + rho C/(kC+C) − lambda P
 *   dH/dt = (1/epsH) [ fP fC Q (1−H) − fPq fC (1−Q) H ],   Q = kH²/(kH²+C²), fP = 1+gH P/(dH+P),
 *           fPq = 1+qoff gH P/(dH+P), fC = 1+rC C²/(kR²+C²)   (mech=1: kH(P) threshold form, fP = fPq = 1)
 * Per parameter point:
 *   cycle features from RK4 (osc, T, amp, W50, duty, Cmin, Cmax, Hmin, Hmax; half-width interpolated; nclip)
 *   fixed points: P*(C), H*(C) are explicit, so g(C) = dC/dt|_(C,P*,H*) is scanned on a log grid in C
 *                 (sign changes bisected); for each fixed point the Jacobian (central differences) and its
 *                 eigenvalues (cubic, one real root by bisection + deflation): nfp, ReL_max (largest real part
 *                 over all fixed points), ImL (imaginary part of that eigenvalue), n_unstable_fp, C_fp of it
 *   class: 0 steady (no cycle, all FPs stable); 1 cycle + unstable FP (Hopf-born); 2 cycle + stable FP
 *          (coexistence); 3 no cycle but unstable FP (slow transient or numerical); +10 if more than one FP
 *   hopf_ratio = T · ImL / (2π): 1 near a Hopf point, ≫ 1 deep in the relaxation regime
 * Usage:
 *   ./m0ip3h_bifurcation_grid grid=1 px=S xlo=0.01 xhi=100 nx=41 logx=1 py=epsH ylo=10 yhi=10000 ny=41 logy=1 [model keys] [threads=N]
 *   ./m0ip3h_bifurcation_grid line=1 px=S xlo=0.003 xhi=300 nx=121 logx=1 [model keys]
 *   ./m0ip3h_bifurcation_grid point=1 [model keys]                    one point, also lists every fixed point
 * Model keys: S alpha beta gamma rho mu lambda kS kC kH epsH mech gH dH qoff kP kR nA kA Jin rC C0 P0 H0 dt Tend Ttrans
 * Output columns: x,y,osc,T,amp,W50,duty,Cmin,Cmax,Hmin,Hmax,nfp,C_fp,ReL,ImL,n_unstable,class,hopf_ratio,nclip
 * Build: see Makefile (OpenMP if available)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct { double S, alpha, beta, gamma, rho, mu, lambda, kS, kC, kH, epsH; double mech, gH, dH, qoff, kP, kR, nA, kA, Jin, rC; double C0, P0, H0, dt, Tend, Ttrans; } Par;
static void defaults(Par *p){ p->S=0; p->alpha=3.5; p->beta=0.26; p->gamma=4.2; p->rho=2.4; p->mu=0.001; p->lambda=0.001; p->kS=0.2; p->kC=0.08; p->kH=0.15; p->epsH=70;
    p->mech=0; p->gH=0; p->dH=1; p->qoff=1; p->kP=1; p->kR=-1; p->nA=1; p->kA=1; p->Jin=0; p->rC=0; p->C0=0.1; p->P0=0.05; p->H0=1; p->dt=0.02; p->Tend=15000; p->Ttrans=9000; }
#define KEY(name) if(!strcmp(k,#name)) return &p->name;
static double *field(Par *p, const char *k){ KEY(S) KEY(alpha) KEY(beta) KEY(gamma) KEY(rho) KEY(mu) KEY(lambda) KEY(kS) KEY(kC) KEY(kH) KEY(epsH) KEY(mech) KEY(gH) KEY(dH) KEY(qoff) KEY(kP) KEY(kR) KEY(nA) KEY(kA) KEY(Jin) KEY(rC) KEY(C0) KEY(P0) KEY(H0) KEY(dt) KEY(Tend) KEY(Ttrans) return NULL; }

static inline void rhs(const Par *p, const double *u, double *du){
    double C=u[0],P=u[1],H=u[2],C2=C*C, Cn=pow(C,p->nA), act=Cn/(pow(p->kA,p->nA)+Cn);
    du[0]=p->Jin + p->alpha*P/(p->kP+P)*act*H - p->beta*C2/(p->kS*p->kS+C2) - p->mu*C;
    du[1]=p->gamma*p->S/(1.0+p->S) + p->rho*C/(p->kC+C) - p->lambda*P;
    double kh=p->kH, fP=1, fPq=1, kR=p->kR>0?p->kR:p->kH, fC=1.0+p->rC*C2/(kR*kR+C2);
    if((int)p->mech==1) kh=p->kH*(p->dH+p->gH*P)/(p->dH+P); else if((int)p->mech>=2){ double x=p->gH*P/(p->dH+P); fP=1+x; fPq=1+p->qoff*x; }
    double Q=kh*kh/(kh*kh+C2);
    du[2]=(fP*fC*Q*(1.0-H) - fPq*fC*(1.0-Q)*H)/p->epsH;
}
static inline int rk4(const Par *p, double *u, double dt){
    double k1[3],k2[3],k3[3],k4[3],v[3]; int i,clip=0;
    rhs(p,u,k1); for(i=0;i<3;i++) v[i]=u[i]+0.5*dt*k1[i]; rhs(p,v,k2); for(i=0;i<3;i++) v[i]=u[i]+0.5*dt*k2[i]; rhs(p,v,k3);
    for(i=0;i<3;i++) v[i]=u[i]+dt*k3[i]; rhs(p,v,k4); for(i=0;i<3;i++) u[i]+=dt/6.0*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
    if(u[0]<0){u[0]=0;clip=1;} if(u[1]<0){u[1]=0;clip=1;} if(u[2]<0){u[2]=0;clip=1;} if(u[2]>1){u[2]=1;clip=1;} return clip;
}

/* ---- fixed points ---------------------------------------------------------------------------------------- */
static void steady_PH(const Par *p, double C, double *P, double *H){
    *P=(p->gamma*p->S/(1.0+p->S) + p->rho*C/(p->kC+C))/p->lambda;
    double C2=C*C, kh=p->kH, fP=1, fPq=1, kR=p->kR>0?p->kR:p->kH, fC=1.0+p->rC*C2/(kR*kR+C2);
    if((int)p->mech==1) kh=p->kH*(p->dH+p->gH*(*P))/(p->dH+(*P)); else if((int)p->mech>=2){ double x=p->gH*(*P)/(p->dH+(*P)); fP=1+x; fPq=1+p->qoff*x; }
    double Q=kh*kh/(kh*kh+C2), A=fP*fC*Q, B=fPq*fC*(1.0-Q); *H=A/(A+B);
}
static double gC(const Par *p, double C){ double u[3],du[3]; u[0]=C; steady_PH(p,C,&u[1],&u[2]); rhs(p,u,du); return du[0]; }
#define MAXFP 8
static int fixed_points(const Par *p, double *Cfp){
    int nscan=3000, i, n=0; double lo=log(1e-7), hi=log(1e3);
    double xprev=exp(lo), gprev=gC(p,xprev);
    for(i=1;i<=nscan;i++){ double x=exp(lo+(hi-lo)*i/(double)nscan), g=gC(p,x);
        if((gprev<0&&g>=0)||(gprev>0&&g<=0)){ double a=xprev,b=x,ga=gprev; int it; for(it=0;it<80;it++){ double m=0.5*(a+b), gm=gC(p,m); if((ga<0&&gm<0)||(ga>0&&gm>0)){a=m;ga=gm;} else b=m; }
            if(n<MAXFP) Cfp[n++]=0.5*(a+b); }
        xprev=x; gprev=g; }
    return n;
}
static void jacobian(const Par *p, const double *u, double J[3][3]){
    int i,j; double h[3]={fmax(1e-6*fabs(u[0]),1e-8), fmax(1e-6*fabs(u[1]),1e-8), fmax(1e-6*fabs(u[2]),1e-8)};
    for(j=0;j<3;j++){ double up[3],um[3],fp[3],fm[3]; for(i=0;i<3;i++){up[i]=u[i];um[i]=u[i];} up[j]+=h[j]; um[j]-=h[j]; if(um[j]<0&&j!=2) um[j]=0; rhs(p,up,fp); rhs(p,um,fm); for(i=0;i<3;i++) J[i][j]=(fp[i]-fm[i])/(up[j]-um[j]); }
}
static inline double cubic(double l,double a,double b,double c){ return ((l+a)*l+b)*l+c; }
/* eigenvalues of a real 3x3: returns max real part and the imaginary part of that eigenvalue */
static void eig3(double J[3][3], double *re, double *im){
    double a=-(J[0][0]+J[1][1]+J[2][2]);
    double b=J[0][0]*J[1][1]-J[0][1]*J[1][0]+J[0][0]*J[2][2]-J[0][2]*J[2][0]+J[1][1]*J[2][2]-J[1][2]*J[2][1];
    double c=-(J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])+J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]));
    /* real root of l^3 + a l^2 + b l + c by bisection on a bracket, then Newton polish */
    double R=1.0+fmax(fabs(a),fmax(fabs(b),fabs(c))); double lo=-R,hi=R; int it;
    for(it=0;it<200;it++){ double m=0.5*(lo+hi); if(cubic(lo,a,b,c)*cubic(m,a,b,c)<=0) hi=m; else lo=m; }
    double r=0.5*(lo+hi); for(it=0;it<5;it++){ double d=(3*r+2*a)*r+b; if(fabs(d)>1e-300) r-=cubic(r,a,b,c)/d; }
    /* deflate: l^2 + (a+r) l + (b + r(a+r)) */
    double q1=a+r, q0=b+r*q1, disc=q1*q1-4*q0; double re1,im1,re2,im2;
    if(disc>=0){ re1=0.5*(-q1+sqrt(disc)); re2=0.5*(-q1-sqrt(disc)); im1=im2=0; } else { re1=re2=-0.5*q1; im1=0.5*sqrt(-disc); im2=-im1; }
    double best=r, bim=0; if(re1>best){best=re1;bim=fabs(im1);} if(re2>best){best=re2;bim=fabs(im2);} *re=best; *im=bim;
}

/* ---- cycle features -------------------------------------------------------------------------------------- */
typedef struct { int osc; double T,amp,W50,duty,Cmin,Cmax,Hmin,Hmax; long nclip; } Feat;
#define NPEAK 6
static Feat run(const Par *p){
    Feat f; memset(&f,0,sizeof f); double u[3]={p->C0,p->P0,p->H0}, dt=p->dt;
    long nsteps=(long)(p->Tend/dt+0.5), ntr=(long)(p->Ttrans/dt+0.5), nrec=nsteps-ntr, n, i;
    double *Cb=malloc(sizeof(double)*nrec), *Hb=malloc(sizeof(double)*nrec);
    for(n=0;n<nsteps;n++){ if(n>=ntr){ Cb[n-ntr]=u[0]; Hb[n-ntr]=u[2]; } f.nclip+=rk4(p,u,dt); }
    double cmax=-1,cmin=1e300,hmin=1e300,hmax=-1; for(i=0;i<nrec;i++){ if(Cb[i]>cmax)cmax=Cb[i]; if(Cb[i]<cmin)cmin=Cb[i]; if(Hb[i]<hmin)hmin=Hb[i]; if(Hb[i]>hmax)hmax=Hb[i]; }
    f.Cmax=cmax; f.Cmin=cmin; f.Hmin=hmin; f.Hmax=hmax; f.amp=cmax-cmin;
    if(cmax-cmin<1e-4||cmax-cmin<0.05*cmax){ free(Cb);free(Hb); return f; }
    double mid=0.5*(cmax+cmin); long pk[NPEAK+1]; int np=0;
    for(i=nrec-2;i>0&&np<NPEAK;i--){ if(Cb[i]>mid&&Cb[i]>=Cb[i-1]&&Cb[i]>Cb[i+1]){ pk[np++]=i; long j=i; while(j>0&&Cb[j]>mid) j--; i=j; } }
    if(np<3){ free(Cb);free(Hb); return f; }
    f.osc=1; double Ts=0; int k; for(k=0;k<np-1;k++) Ts+=(pk[k]-pk[k+1])*dt; f.T=Ts/(np-1);
    long ip=pk[1]; double base=1e300; for(i=pk[2];i<=ip;i++) if(Cb[i]<base) base=Cb[i];
    double half=base+0.5*(Cb[ip]-base); long a=ip,b=ip; while(a>0&&Cb[a]>half) a--; while(b<nrec-1&&Cb[b]>half) b++;
    double ta=(a+(half-Cb[a])/fmax(Cb[a+1]-Cb[a],1e-300))*dt, tb=(b-1+(Cb[b-1]-half)/fmax(Cb[b-1]-Cb[b],1e-300))*dt;
    f.W50=tb-ta; f.duty=f.W50/f.T; f.amp=Cb[ip]-base;
    free(Cb);free(Hb); return f;
}

typedef struct { Feat f; int nfp, nunst, cls; double Cfp, ReL, ImL, hopf_ratio; } Res;
static Res evaluate(const Par *p0){
    Par pp=*p0; if(pp.Tend<25.0*pp.epsH){ pp.Tend=25.0*pp.epsH; pp.Ttrans=0.6*pp.Tend; }   /* slow cells: period ~ 0.1-0.5 epsH, keep >= 10 cycles */
    const Par *p=&pp; Res r; memset(&r,0,sizeof r); r.f=run(p);
    double Cfp[MAXFP]; int n=fixed_points(p,Cfp), i; r.nfp=n; r.ReL=-1e300; r.ImL=0;
    for(i=0;i<n;i++){ double u[3]; u[0]=Cfp[i]; steady_PH(p,Cfp[i],&u[1],&u[2]); double J[3][3]; jacobian(p,u,J); double re,im; eig3(J,&re,&im);
        if(re>0) r.nunst++; if(re>r.ReL){ r.ReL=re; r.ImL=im; r.Cfp=Cfp[i]; } }
    if(n==0){ r.ReL=0; }
    if(r.f.osc) r.cls=(r.nunst>0)?1:2; else r.cls=(r.nunst>0)?3:0;
    if(n>1) r.cls+=10;
    r.hopf_ratio=(r.f.osc&&r.ImL>0)?r.f.T*r.ImL/(2*M_PI):0;
    return r;
}
static void print_res(double x, double y, const Res *r){
    printf("%.6g,%.6g,%d,%.4f,%.5g,%.4f,%.5f,%.5g,%.5g,%.5f,%.5f,%d,%.5g,%.6g,%.6g,%d,%d,%.4f,%ld\n",x,y,r->f.osc,r->f.T,r->f.amp,r->f.W50,r->f.duty,r->f.Cmin,r->f.Cmax,r->f.Hmin,r->f.Hmax,r->nfp,r->Cfp,r->ReL,r->ImL,r->nunst,r->cls,r->hopf_ratio,r->f.nclip);
}

int main(int argc,char**argv){
    Par base; defaults(&base); char px[32]="S", py[32]="epsH"; double xlo=0.01,xhi=100,ylo=10,yhi=1e4; int nx=41,ny=41,logx=1,logy=1,grid=0,line=0,point=0,threads=0,i;
    for(i=1;i<argc;i++){ char key[64],val[64]; if(sscanf(argv[i],"%63[^=]=%63s",key,val)!=2){fprintf(stderr,"bad arg %s\n",argv[i]);return 1;}
        if(!strcmp(key,"px")) strncpy(px,val,31); else if(!strcmp(key,"py")) strncpy(py,val,31);
        else if(!strcmp(key,"xlo")) xlo=atof(val); else if(!strcmp(key,"xhi")) xhi=atof(val); else if(!strcmp(key,"ylo")) ylo=atof(val); else if(!strcmp(key,"yhi")) yhi=atof(val);
        else if(!strcmp(key,"nx")) nx=atoi(val); else if(!strcmp(key,"ny")) ny=atoi(val); else if(!strcmp(key,"logx")) logx=atoi(val); else if(!strcmp(key,"logy")) logy=atoi(val);
        else if(!strcmp(key,"grid")) grid=atoi(val); else if(!strcmp(key,"line")) line=atoi(val); else if(!strcmp(key,"point")) point=atoi(val); else if(!strcmp(key,"threads")) threads=atoi(val);
        else { double *f=field(&base,key); if(!f){fprintf(stderr,"unknown key %s\n",key);return 1;} *f=atof(val); } }
#ifdef _OPENMP
    if(threads>0) omp_set_num_threads(threads); fprintf(stderr,"OpenMP: %d threads\n",omp_get_max_threads());
#endif
    printf("x,y,osc,T,amp,W50,duty,Cmin,Cmax,Hmin,Hmax,nfp,C_fp,ReL,ImL,n_unstable,class,hopf_ratio,nclip\n");
    if(point){ Res r=evaluate(&base); print_res(0,0,&r); double Cfp[MAXFP]; int n=fixed_points(&base,Cfp); fprintf(stderr,"fixed points (C*, P*, H*, max Re lambda, Im):\n");
        for(i=0;i<n;i++){ double u[3]; u[0]=Cfp[i]; steady_PH(&base,Cfp[i],&u[1],&u[2]); double J[3][3]; jacobian(&base,u,J); double re,im; eig3(J,&re,&im); fprintf(stderr,"  %.5g %.5g %.5g %.5g %.5g\n",u[0],u[1],u[2],re,im); } return 0; }
    if(!field(&base,px)||(grid&&!field(&base,py))){ fprintf(stderr,"unknown px/py\n"); return 1; }
    int NX=nx, NYY=grid?ny:1, total=NX*NYY, c;
    Res *R=malloc(sizeof(Res)*total); double *X=malloc(sizeof(double)*total), *Y=malloc(sizeof(double)*total);
    for(c=0;c<total;c++){ int ix=c%NX, iy=c/NX; X[c]=logx?xlo*pow(xhi/xlo,(double)ix/(NX>1?NX-1:1)):xlo+(xhi-xlo)*ix/(double)(NX>1?NX-1:1); Y[c]=grid?(logy?ylo*pow(yhi/ylo,(double)iy/(NYY>1?NYY-1:1)):ylo+(yhi-ylo)*iy/(double)(NYY>1?NYY-1:1)):0; }
    long done=0;
#pragma omp parallel for schedule(dynamic,1)
    for(c=0;c<total;c++){ Par p=base; *field(&p,px)=X[c]; if(grid) *field(&p,py)=Y[c]; R[c]=evaluate(&p);
#pragma omp atomic
        done++;
        if(done%100==0){
#pragma omp critical
            fprintf(stderr,"  %ld / %d\n",done,total); } }
    for(c=0;c<total;c++) print_res(X[c],Y[c],&R[c]);
    free(R);free(X);free(Y); return 0;
}
