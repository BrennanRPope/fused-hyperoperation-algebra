/* DRAFT. Hyperspherical bookkeeping along an orbit. Jacobi vectors rho = (r2-r1)/sqrt2, lam = sqrt(2/3)(r3-(r1+r2)/2);
   hyperradius R = sqrt(|rho|^2+|lam|^2); shape-sphere point via the Hopf map. Accumulates over [0,T]:
   out[0] = tau = int dt / R^2, out[1] = shape-sphere arclength, out[2] = number of maxima of R,
   out[3] = min R, out[4] = max R, out[5] = int R^2 dt / T (mean moment of inertia), out[6] = steps. */
#include "tbw.c"
static void jac(const double *X, double *s, double *R2){
    double r1 = (X[2]-X[0])/sqrt(2.), r2 = (X[3]-X[1])/sqrt(2.);
    double l1 = sqrt(2./3)*(X[4]-0.5*(X[0]+X[2])), l2 = sqrt(2./3)*(X[5]-0.5*(X[1]+X[3]));
    double a = r1*r1+r2*r2, b = l1*l1+l2*l2; *R2 = a+b;
    s[0] = (a-b)/(a+b); s[1] = 2*(r1*l1+r2*l2)/(a+b); s[2] = 2*(r1*l2-r2*l1)/(a+b);
}
static double dR2(const double *X){ /* d/dt of R^2 / 2 = sum m r.v about the centre, centre is fixed at 0 */
    double s=0; for(int i=0;i<3;i++) s += X[2*i]*X[6+2*i]+X[2*i+1]*X[7+2*i]; return s; }
int hyper_stats(double p1, double p2, double T, double tol, double hmax, double *out){
    double X[12], Xn[12]; suvakov_ic(p1,p2,X);
    double t=0, h=1e-3, s0[3], s1[3], R20, R21; int ns=0, nmax=0;
    jac(X,s0,&R20); double tau=0, arc=0, I=0, rmin=sqrt(R20), rmax=rmin, d0=dR2(X);
    while (t < T) {
        if (h > hmax) h = hmax; if (t+h > T) h = T-t;
        double err = step(X,h,Xn,tol);
        if (err <= 1.0) {
            jac(Xn,s1,&R21); double d1 = dR2(Xn);
            tau += 0.5*h*(1/R20+1/R21); I += 0.5*h*(R20+R21);
            arc += sqrt((s1[0]-s0[0])*(s1[0]-s0[0])+(s1[1]-s0[1])*(s1[1]-s0[1])+(s1[2]-s0[2])*(s1[2]-s0[2]));
            if (d0 > 0 && d1 <= 0) nmax++;
            double r = sqrt(R21); if (r<rmin) rmin=r; if (r>rmax) rmax=r;
            t += h; memcpy(X,Xn,sizeof Xn); memcpy(s0,s1,sizeof s1); R20=R21; d0=d1; ns++;
        }
        double fac = 0.9*pow(err > 1e-16 ? err : 1e-16, -0.2);
        if (fac > 5) fac = 5; if (fac < 0.2) fac = 0.2;
        h *= fac;
        if (h < 1e-14 || ns > 50000000) return -1;
    }
    out[0]=tau; out[1]=arc; out[2]=nmax; out[3]=rmin; out[4]=rmax; out[5]=I/T; out[6]=ns; return ns;
}
/* positions at n+1 equally spaced times on [0,T]; out has 6*(n+1) doubles */
int sample(double p1, double p2, double T, int n, double tol, double *out){
    double X[12]; suvakov_ic(p1,p2,X); memcpy(out,X,6*sizeof(double));
    for (int k=1;k<=n;k++){ if (integrate(X,T/n,tol) < 0) return -1; memcpy(out+6*k,X,6*sizeof(double)); }
    return 0;
}
