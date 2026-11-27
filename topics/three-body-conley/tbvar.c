/* DRAFT. Planar equal-mass three-body problem (G = m = 1) with variational equations.
   Dormand-Prince 5(4), step size controlled on the state only, so the 12x12 Jacobian is
   integrated along exactly the same steps as the orbit.
   State X = (x1,y1,x2,y2,x3,y3, u1,v1,u2,v2,u3,v3); Phi is 12x12 row-major. */
#include <math.h>
#include <string.h>
#define N 12
#define NV (N + N*N)

static void accel(const double *X, double *a, double *Hs /* 6x6 d(accel)/d(pos), may be NULL */) {
    for (int k = 0; k < 6; k++) a[k] = 0;
    if (Hs) for (int k = 0; k < 36; k++) Hs[k] = 0;
    for (int i = 0; i < 3; i++) for (int j = i+1; j < 3; j++) {
        double dx = X[2*j]-X[2*i], dy = X[2*j+1]-X[2*i+1];
        double r2 = dx*dx+dy*dy, r = sqrt(r2), r3 = r2*r, r5 = r3*r2;
        double fx = dx/r3, fy = dy/r3;
        a[2*i] += fx; a[2*i+1] += fy; a[2*j] -= fx; a[2*j+1] -= fy;
        if (Hs) {
            /* d f / d (dx,dy) where f = d/r^3 : B = I/r^3 - 3 d d^T / r^5 */
            double B[2][2] = {{1/r3 - 3*dx*dx/r5, -3*dx*dy/r5}, {-3*dx*dy/r5, 1/r3 - 3*dy*dy/r5}};
            for (int p = 0; p < 2; p++) for (int q = 0; q < 2; q++) {
                double b = B[p][q];
                /* a_i += f(d), d = r_j - r_i ; a_j -= f(d) */
                Hs[(2*i+p)*6 + 2*j+q] += b;  Hs[(2*i+p)*6 + 2*i+q] -= b;
                Hs[(2*j+p)*6 + 2*j+q] -= b;  Hs[(2*j+p)*6 + 2*i+q] += b;
            }
        }
    }
}

static void deriv(const double *Y, double *D, int var) {
    double a[6], Hs[36];
    accel(Y, a, var ? Hs : 0);
    for (int k = 0; k < 6; k++) { D[k] = Y[6+k]; D[6+k] = a[k]; }
    if (!var) return;
    const double *P = Y + N; double *DP = D + N;
    /* dPhi/dt = A Phi, A = [[0, I],[Hs, 0]] */
    for (int r = 0; r < 6; r++) for (int c = 0; c < N; c++) DP[r*N+c] = P[(6+r)*N+c];
    for (int r = 0; r < 6; r++) for (int c = 0; c < N; c++) {
        double s = 0; for (int k = 0; k < 6; k++) s += Hs[r*6+k]*P[k*N+c];
        DP[(6+r)*N+c] = s;
    }
}

static const double c2=1./5, c3=3./10, c4=4./5, c5=8./9;
static const double a21=1./5, a31=3./40, a32=9./40, a41=44./45, a42=-56./15, a43=32./9,
 a51=19372./6561, a52=-25360./2187, a53=64448./6561, a54=-212./729,
 a61=9017./3168, a62=-355./33, a63=46732./5247, a64=49./176, a65=-5103./18656,
 a71=35./384, a73=500./1113, a74=125./192, a75=-2187./6784, a76=11./84;
static const double e1=71./57600, e3=-71./16695, e4=71./1920, e5=-17253./339200, e6=22./525, e7=-1./40;

static double step(const double *X, double h, double *Xn, double tol, int var) {
    static double k1[NV],k2[NV],k3[NV],k4[NV],k5[NV],k6[NV],k7[NV],Y[NV];
    int n = var ? NV : N;
    deriv(X,k1,var);
    for(int i=0;i<n;i++) Y[i]=X[i]+h*a21*k1[i]; deriv(Y,k2,var);
    for(int i=0;i<n;i++) Y[i]=X[i]+h*(a31*k1[i]+a32*k2[i]); deriv(Y,k3,var);
    for(int i=0;i<n;i++) Y[i]=X[i]+h*(a41*k1[i]+a42*k2[i]+a43*k3[i]); deriv(Y,k4,var);
    for(int i=0;i<n;i++) Y[i]=X[i]+h*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]); deriv(Y,k5,var);
    for(int i=0;i<n;i++) Y[i]=X[i]+h*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]); deriv(Y,k6,var);
    for(int i=0;i<n;i++) Xn[i]=X[i]+h*(a71*k1[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]);
    deriv(Xn,k7,var);
    double err=0;
    for(int i=0;i<N;i++){ /* error on the state only */
        double e = h*(e1*k1[i]+e3*k3[i]+e4*k4[i]+e5*k5[i]+e6*k6[i]+e7*k7[i]);
        double sc = tol*(1+fabs(X[i]));
        err += (e/sc)*(e/sc);
    }
    return sqrt(err/N);
}

static int run(double *X, double T, double tol, int var) {
    static double Xn[NV];
    double t=0, h=1e-3; int ns=0, n = var ? NV : N;
    if (T <= 0) return 0;
    while (t < T) {
        if (t+h > T) h = T-t;
        double err = step(X,h,Xn,tol,var);
        if (err <= 1.0) { t += h; memcpy(X,Xn,n*sizeof(double)); ns++; }
        double fac = 0.9*pow(err > 1e-16 ? err : 1e-16, -0.2);
        if (fac > 5) fac = 5; if (fac < 0.2) fac = 0.2;
        h *= fac;
        if (h < 1e-15 || ns > 50000000) return -1;
    }
    return ns;
}

/* state only */
int integrate(double *X, double T, double tol) { return run(X, T, tol, 0); }
/* state + Jacobian: XP has 12 + 144 entries; caller sets XP[12..] (e.g. identity) */
int integrate_var(double *XP, double T, double tol) { return run(XP, T, tol, 1); }

/* ---- section hits: line-up with body 3 between 1 and 2, crossed with dA/dt < 0 ---- */
static double areaA(const double *X) {
    double ax = X[2]-X[0], ay = X[3]-X[1], bx = X[4]-X[0], by = X[5]-X[1];
    return ax*by - ay*bx;
}
static int hitok(const double *X) {
    double ax = X[2]-X[0], ay = X[3]-X[1], bx = X[4]-X[0], by = X[5]-X[1];
    double t = (bx*ax + by*ay)/(ax*ax + ay*ay);
    double dA = (X[8]-X[6])*by - (X[9]-X[7])*bx + ax*(X[11]-X[7]) - ay*(X[10]-X[6]);
    return (t > 0 && t < 1 && dA < 0);
}
/* Flow X forward to the next hit (at least tmin later). Returns time taken, or -1. */
double next_hit(double *X, double tol, double tmax) {
    double t = 0, h = 1e-3, Xn[N], Xs[N];
    double A0 = areaA(X);
    int ns = 0;
    while (t < tmax) {
        double err = step(X, h, Xn, tol, 0);
        if (err <= 1.0) {
            double A1 = areaA(Xn);
            if (t > 1e-9 && A0 > 0 && A1 <= 0) {
                /* bisect on the step length from X with single DP5 steps (h' <= accepted h) */
                double lo = 0, hi = h;
                for (int it = 0; it < 80 && hi - lo > 1e-16; it++) {
                    double mid = 0.5*(lo+hi);
                    step(X, mid, Xs, tol, 0);
                    if (areaA(Xs) > 0) lo = mid; else hi = mid;
                }
                step(X, hi, Xs, tol, 0);
                if (hitok(Xs)) { memcpy(X, Xs, sizeof Xs); return t + hi; }
            }
            t += h; memcpy(X, Xn, sizeof Xn); A0 = A1; ns++;
        }
        double fac = 0.9*pow(err > 1e-16 ? err : 1e-16, -0.2);
        if (fac > 5) fac = 5; if (fac < 0.2) fac = 0.2;
        h *= fac;
        if (h < 1e-15 || ns > 50000000) return -1;
    }
    return -1;
}
