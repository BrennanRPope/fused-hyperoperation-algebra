/* DRAFT. Mass-aware copy of tbw.c (masses 1, 1, m3; middle body gets -2v/m3, the Li-Jing-Liao convention).
   (copied from ../three-body/tb.c, plus syzygy letters). Planar equal-mass three-body problem (G = m = 1), Dormand-Prince 5(4) adaptive.
   Built as a shared library and called from checks.py with ctypes.
   State X = (x1,y1,x2,y2,x3,y3, u1,v1,u2,v2,u3,v3). */
#include <math.h>
#include <string.h>

static double MS[3] = {1,1,1};
void set_masses(double a, double b, double c){ MS[0]=a; MS[1]=b; MS[2]=c; }
static void deriv(const double *X, double *D) {
    double ax[3] = {0,0,0}, ay[3] = {0,0,0};
    for (int i = 0; i < 3; i++) for (int j = i+1; j < 3; j++) {
        double dx = X[2*j]-X[2*i], dy = X[2*j+1]-X[2*i+1];
        double r2 = dx*dx+dy*dy, r3 = r2*sqrt(r2);
        double fx = dx/r3, fy = dy/r3;
        ax[i] += MS[j]*fx; ay[i] += MS[j]*fy; ax[j] -= MS[i]*fx; ay[j] -= MS[i]*fy;
    }
    for (int k = 0; k < 6; k++) D[k] = X[6+k];
    for (int i = 0; i < 3; i++) { D[6+2*i] = ax[i]; D[7+2*i] = ay[i]; }
}

static const double c2=1./5, c3=3./10, c4=4./5, c5=8./9;
static const double a21=1./5, a31=3./40, a32=9./40, a41=44./45, a42=-56./15, a43=32./9,
 a51=19372./6561, a52=-25360./2187, a53=64448./6561, a54=-212./729,
 a61=9017./3168, a62=-355./33, a63=46732./5247, a64=49./176, a65=-5103./18656,
 a71=35./384, a73=500./1113, a74=125./192, a75=-2187./6784, a76=11./84;
static const double e1=71./57600, e3=-71./16695, e4=71./1920, e5=-17253./339200, e6=22./525, e7=-1./40;

/* one adaptive step attempt; returns error norm, writes Xn */
static double step(const double *X, double h, double *Xn, double tol) {
    double k1[12],k2[12],k3[12],k4[12],k5[12],k6[12],k7[12],Y[12];
    int n = 12;
    deriv(X,k1);
    for(int i=0;i<n;i++) Y[i]=X[i]+h*a21*k1[i]; deriv(Y,k2);
    for(int i=0;i<n;i++) Y[i]=X[i]+h*(a31*k1[i]+a32*k2[i]); deriv(Y,k3);
    for(int i=0;i<n;i++) Y[i]=X[i]+h*(a41*k1[i]+a42*k2[i]+a43*k3[i]); deriv(Y,k4);
    for(int i=0;i<n;i++) Y[i]=X[i]+h*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]); deriv(Y,k5);
    for(int i=0;i<n;i++) Y[i]=X[i]+h*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]); deriv(Y,k6);
    for(int i=0;i<n;i++) Xn[i]=X[i]+h*(a71*k1[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]);
    deriv(Xn,k7);
    double err=0;
    for(int i=0;i<n;i++){
        double e = h*(e1*k1[i]+e3*k3[i]+e4*k4[i]+e5*k5[i]+e6*k6[i]+e7*k7[i]);
        double sc = tol*(1.0+fmax(fabs(X[i]),fabs(Xn[i])));
        err += (e/sc)*(e/sc);
    }
    return sqrt(err/n);
}

/* integrate X in place from 0 to T. returns number of steps (negative on failure) */
int integrate(double *X, double T, double tol) {
    double t=0, h=1e-3, Xn[12]; int ns=0;
    if (T <= 0) return 0;
    while (t < T) {
        if (t+h > T) h = T-t;
        double err = step(X,h,Xn,tol);
        if (err <= 1.0) { t += h; memcpy(X,Xn,sizeof Xn); ns++; }
        double fac = 0.9*pow(err > 1e-16 ? err : 1e-16, -0.2);
        if (fac > 5) fac = 5; if (fac < 0.2) fac = 0.2;
        h *= fac;
        if (h < 1e-14 || ns > 50000000) return -1;
    }
    return ns;
}


static void suvakov_ic(double p1, double p2, double *X) {
    double I[12] = {-1,0, 1,0, 0,0, p1,p2, p1,p2, -2*p1/MS[2],-2*p2/MS[2]};
    memcpy(X,I,sizeof I);
}
static double area(const double *X){ return (X[2]-X[0])*(X[5]-X[1]) - (X[3]-X[1])*(X[4]-X[0]); }
/* middle body (1..3) of a collinear configuration */
static int middle(const double *X){
    double dx=0,dy=0,best=-1;
    for(int i=0;i<3;i++)for(int j=i+1;j<3;j++){double ex=X[2*j]-X[2*i],ey=X[2*j+1]-X[2*i+1],d=ex*ex+ey*ey; if(d>best){best=d;dx=ex;dy=ey;}}
    double s[3]; for(int k=0;k<3;k++) s[k]=X[2*k]*dx+X[2*k+1]*dy;
    int lo=0,hi=0; for(int k=1;k<3;k++){ if(s[k]<s[lo])lo=k; if(s[k]>s[hi])hi=k; }
    return 6-lo-hi-3+1; /* the remaining index, 1-based */
}
/* Integrate from the slice IC; record syzygies with middle body != 3 as letters
   +1=a (body1 middle, area rising), -1=A, +2=b (body2), -2=B. Also returns times. Stops at maxev letters or Tmax.
   Returns number of letters, or -(n+1) on failure. */
int letters(double p1, double p2, double Tmax, double tol, int maxev, int *L, double *tL){
    double X[12], Xn[12]; suvakov_ic(p1,p2,X);
    double t=0,h=1e-3; int ns=0, n=0;
    double a0=area(X);
    while (t < Tmax && n < maxev) {
        if (t+h > Tmax) h = Tmax-t;
        double err = step(X,h,Xn,tol);
        if (err <= 1.0) {
            double a1=area(Xn);
            if ((a0<0 && a1>=0) || (a0>0 && a1<=0)) {
                double f = a0/(a0-a1), Y[12];
                for(int i=0;i<12;i++) Y[i]=X[i]+f*(Xn[i]-X[i]);
                int m = middle(Y);
                if (m != 3) { L[n] = (a1>a0 ? 1 : -1) * m; tL[n] = t+f*h; n++; }
            }
            a0=a1; t+=h; memcpy(X,Xn,sizeof Xn); ns++;
        }
        double fac = 0.9*pow(err > 1e-16 ? err : 1e-16, -0.2);
        if (fac > 5) fac = 5; if (fac < 0.2) fac = 0.2;
        h *= fac;
        if (h < 1e-14 || ns > 50000000) return -(n+1);
    }
    return n;
}
void letters_grid(const double *P1,const double *P2,int nq,double Tmax,double tol,int maxev,int *L,int *cnt){
    #pragma omp parallel for schedule(dynamic,8)
    for(int q=0;q<nq;q++){
        double tL[512]; if (maxev > 512) maxev = 512;
        int n = letters(P1[q],P2[q],Tmax,tol,maxev,L+(long)q*maxev,tL);
        cnt[q]=n;
    }
}

/* For each start, the closest return |X(t)-X0| for t in [T0[q], T1[q]] (checked at every accepted step,
   then refined by one cubic-free bisection on the step). Also returns the time of the closest return. */
void retdist_grid(const double *P1,const double *P2,const double *T0,const double *T1,int nq,double tol,double *dmin,double *tmin){
    #pragma omp parallel for schedule(dynamic,4)
    for(int q=0;q<nq;q++){
        double X0[12],X[12],Xn[12]; suvakov_ic(P1[q],P2[q],X0); memcpy(X,X0,sizeof X);
        double t=0,h=1e-3,best=1e9,bt=0; int ns=0;
        while(t<T1[q]){
            if(t+h>T1[q]) h=T1[q]-t;
            double err=step(X,h,Xn,tol);
            if(err<=1.0){ t+=h; memcpy(X,Xn,sizeof Xn); ns++;
                if(t>=T0[q]){ double dd=0; for(int i=0;i<12;i++){double e=X[i]-X0[i]; dd+=e*e;} dd=sqrt(dd); if(dd<best){best=dd;bt=t;} }
            }
            double fac=0.9*pow(err>1e-16?err:1e-16,-0.2); if(fac>5)fac=5; if(fac<0.2)fac=0.2; h*=fac;
            if(h<1e-14||ns>50000000){best=1e9;break;}
        }
        dmin[q]=best; tmin[q]=bt;
    }
}
