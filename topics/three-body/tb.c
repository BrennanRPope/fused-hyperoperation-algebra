/* DRAFT. Planar equal-mass three-body problem (G = m = 1), Dormand-Prince 5(4) adaptive.
   Built as a shared library and called from checks.py with ctypes.
   State X = (x1,y1,x2,y2,x3,y3, u1,v1,u2,v2,u3,v3). */
#include <math.h>
#include <string.h>

static void deriv(const double *X, double *D) {
    double ax[3] = {0,0,0}, ay[3] = {0,0,0};
    for (int i = 0; i < 3; i++) for (int j = i+1; j < 3; j++) {
        double dx = X[2*j]-X[2*i], dy = X[2*j+1]-X[2*i+1];
        double r2 = dx*dx+dy*dy, r3 = r2*sqrt(r2);
        double fx = dx/r3, fy = dy/r3;
        ax[i] += fx; ay[i] += fy; ax[j] -= fx; ay[j] -= fy;
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

/* Integrate and record states at given increasing times (nt of them) into out[nt*12]. */
int sample(double *X, const double *times, int nt, double *out, double tol) {
    double t = 0; int tot = 0;
    for (int k = 0; k < nt; k++) {
        int r = integrate(X, times[k]-t, tol);
        if (r < 0) return -1; tot += r;
        t = times[k]; memcpy(out+12*k, X, 12*sizeof(double));
    }
    return tot;
}

static void suvakov_ic(double p1, double p2, double *X) {
    double I[12] = {-1,0, 1,0, 0,0, p1,p2, p1,p2, -2*p1,-2*p2};
    memcpy(X,I,sizeof I);
}

/* Return-distance landscape over the Suvakov slice. For each (p1,p2): integrate to Tmax,
   track at each accepted step after leaving (dist > leave) the minimum of
   abs distance |X(t)-X0| and rotation-free distance min_alpha |X(t) - R_alpha X0|. */
void landscape(const double *P1, const double *P2, int n, double Tmax, double tol, double leave,
               double *dabs, double *tabs, double *drot, double *trot, double *arot) {
    #pragma omp parallel for schedule(dynamic,4)
    for (int q = 0; q < n; q++) {
        double X0[12], X[12], Xn[12];
        suvakov_ic(P1[q],P2[q],X0); memcpy(X,X0,sizeof X);
        double n0 = 0; for (int i=0;i<12;i++) n0 += X0[i]*X0[i];
        double t=0,h=1e-3; int left=0, ns=0;
        double bA=1e9,bAt=0,bR=1e9,bRt=0,bRa=0;
        while (t < Tmax) {
            if (t+h > Tmax) h = Tmax-t;
            double err = step(X,h,Xn,tol);
            if (err <= 1.0) {
                t += h; memcpy(X,Xn,sizeof Xn); ns++;
                double da=0, re=0, im=0, nx=0;
                for (int i=0;i<12;i++){ double d=X[i]-X0[i]; da+=d*d; nx+=X[i]*X[i]; }
                for (int i=0;i<6;i++){ double a=X0[2*i],b=X0[2*i+1],c=X[2*i],d=X[2*i+1];
                    re += a*c+b*d; im += a*d-b*c; }
                da = sqrt(da);
                double dr = sqrt(fmax(0, nx + n0 - 2*sqrt(re*re+im*im)));
                if (!left && da > leave) left = 1;
                if (left) {
                    if (da < bA) { bA=da; bAt=t; }
                    if (dr < bR) { bR=dr; bRt=t; bRa=atan2(im,re); }
                }
            }
            double fac = 0.9*pow(err > 1e-16 ? err : 1e-16, -0.2);
            if (fac > 5) fac = 5; if (fac < 0.2) fac = 0.2;
            h *= fac;
            if (h < 1e-13 || ns > 20000000) break;
        }
        dabs[q]=bA; tabs[q]=bAt; drot[q]=bR; trot[q]=bRt; arot[q]=bRa;
    }
}
