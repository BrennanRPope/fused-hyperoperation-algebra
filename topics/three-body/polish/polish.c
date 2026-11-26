/* DRAFT. Newton polish of three-body periodic orbits in the Suvakov slice, using the half-period symmetry.
   Build: gcc -O3 -o polish polish.c -lquadmath -lm
   Input (stdin): one line per orbit: id p1 p2 T (decimal). Output: one result block per orbit. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <quadmath.h>
#define R double
#define FN(x) x##_d
#define SQRT sqrt
#define FABS fabs
#define LOG log
#define POW pow
#include "taylor.h"
#undef R
#undef FN
#undef SQRT
#undef FABS
#undef LOG
#undef POW
#undef NMAX
#define R __float128
#define FN(x) x##_q
#define SQRT sqrtq
#define FABS fabsq
#define LOG logq
#define POW powq
#include "taylor.h"

static void qs(char *b, __float128 x) { quadmath_snprintf(b, 64, "%.34Qe", x); }
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

int main(int argc, char **argv) {
    double qeps1 = argc > 1 ? atof(argv[1]) : 1e-33;
    char id[64], s1[80], s2[80], s3[80];
    while (scanf("%63s %79s %79s %79s", id, s1, s2, s3) == 4) {
        double t0 = now(); double hist[64]; int nit;
        /* 0: how far the starting orbit is from the half-period condition (quad) */
        __float128 z0[3] = {strtoflt128(s1,0), strtoflt128(s2,0), strtoflt128(s3,0)/2}, G0[3];
        halfres_q(z0, 1e-33Q, G0);
        printf("ORBIT %s\nstart_halfres %.3e %.3e %.3e\n", id, (double)G0[0], (double)G0[1], (double)G0[2]);
        /* 1: double-precision polish */
        double zd[3] = {atof(s1), atof(s2), atof(s3)/2};
        int rc = newton_d(zd, 1e-16, 1e-7, 8, 1e-13, hist, &nit);
        printf("double rc %d it %d", rc, nit); for (int k=0;k<nit;k++) printf(" [%.1e %.1e]", hist[2*k], hist[2*k+1]);
        printf("\ndouble_z %.17e %.17e %.17e\n", zd[0], zd[1], 2*zd[2]);
        double Xd[12]; suvakov_d(zd[0], zd[1], Xd); double X0d[12]; for(int i=0;i<12;i++) X0d[i]=Xd[i];
        integrate_d(Xd, 2*zd[2], 1e-16, 0); double cd=0; for(int i=0;i<12;i++) cd=fmax(cd,fabs(Xd[i]-X0d[i]));
        printf("double_close %.3e time %.2f\n", cd, now()-t0);
        /* 2: quad polish, two tolerances */
        __float128 za[3], zb[3];
        for (int pass = 0; pass < 2; pass++) {
            __float128 *z = pass ? zb : za; double t1 = now();
            z[0] = zd[0]; z[1] = zd[1]; z[2] = zd[2];
            __float128 eps = pass ? (__float128)1e-30 : (__float128)qeps1;
            rc = newton_q(z, eps, 1e-17Q, 8, 1e-31Q, hist, &nit);
            printf("quad%d rc %d it %d", pass, rc, nit); for (int k=0;k<nit;k++) printf(" [%.1e %.1e]", hist[2*k], hist[2*k+1]);
            printf(" time %.2f\n", now()-t1);
        }
        char b[3][64]; for (int i=0;i<3;i++) qs(b[i], i<2 ? za[i] : 2*za[i]);
        printf("quad_z %s %s %s\n", b[0], b[1], b[2]);
        double dd = 0; for (int i=0;i<3;i++) { double r=(double)fabsq((za[i]-zb[i])/(fabsq(za[i])>1e-3Q?za[i]:1)); if (r>dd) dd=r; }
        printf("quad_tol_diff %.3e\n", dd);
        /* 3: full-period closure in quad, energy drift, closest approach */
        __float128 X[12], X0[12], rmin; suvakov_q(za[0], za[1], X); for(int i=0;i<12;i++) X0[i]=X[i];
        long ns = integrate_q(X, 2*za[2], (__float128)qeps1, &rmin);
        __float128 c=0; for(int i=0;i<12;i++){ __float128 d=fabsq(X[i]-X0[i]); if(d>c) c=d; }
        printf("quad_close %.3e steps %ld rmin %.5f total_time %.2f\nEND\n", (double)c, ns, (double)rmin, now()-t0);
        fflush(stdout);
    }
    return 0;
}
