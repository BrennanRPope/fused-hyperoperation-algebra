/* DRAFT. Taylor-series integrator for the planar equal-mass three-body problem (G = m = 1).
   Included twice by polish.c: once with R = double, once with R = __float128.
   State X = (x1,y1,x2,y2,x3,y3, u1,v1,u2,v2,u3,v3), as in ../tb.c. */
#define NMAX 64
static void FN(coeffs)(const R *X, int N, R q[12][NMAX+2]) {
    static const int P[3][2] = {{0,1},{0,2},{1,2}};
    R dx[3][NMAX+1], dy[3][NMAX+1], s[3][NMAX+1], f[3][NMAX+1];
    for (int i = 0; i < 12; i++) q[i][0] = X[i];
    for (int n = 0; n < N; n++) {
        R a[6]; for (int k = 0; k < 6; k++) a[k] = 0;
        for (int p = 0; p < 3; p++) {
            int i = P[p][0], j = P[p][1];
            dx[p][n] = q[2*j][n] - q[2*i][n];
            dy[p][n] = q[2*j+1][n] - q[2*i+1][n];
            R sn = 0; for (int k = 0; k <= n; k++) sn += dx[p][k]*dx[p][n-k] + dy[p][k]*dy[p][n-k];
            s[p][n] = sn;
            if (n == 0) f[p][0] = 1 / (s[p][0] * SQRT(s[p][0]));
            else {
                R acc = 0;
                for (int k = 0; k < n; k++) acc += ((R)(-1.5)*(n-k) - k) * s[p][n-k] * f[p][k];
                f[p][n] = acc / (n * s[p][0]);
            }
            R gx = 0, gy = 0;
            for (int k = 0; k <= n; k++) { gx += dx[p][k]*f[p][n-k]; gy += dy[p][k]*f[p][n-k]; }
            a[2*i] += gx; a[2*i+1] += gy; a[2*j] -= gx; a[2*j+1] -= gy;
        }
        for (int k = 0; k < 6; k++) { q[k][n+1] = q[6+k][n] / (n+1); q[6+k][n+1] = a[k] / (n+1); }
    }
}
static R FN(norm_c)(R q[12][NMAX+2], int n) { R m = 0; for (int i = 0; i < 12; i++) { R v = FABS(q[i][n]); if (v > m) m = v; } return m; }

/* integrate X in place from 0 to T with tolerance eps; returns steps (negative on failure).
   rmin (if not NULL) receives the smallest pair distance seen at step nodes. */
static long FN(integrate)(R *X, R T, R eps, R *rmin) {
    int N = (int)(-LOG(eps) / 2) + 2; if (N > NMAX) N = NMAX;
    R q[12][NMAX+2]; R t = 0; long ns = 0;
    if (rmin) *rmin = 1e9;
    while (t < T) {
        FN(coeffs)(X, N, q);
        R x0 = FN(norm_c)(q, 0); if (x0 < 1) x0 = 1;
        R r1 = POW(x0 / FN(norm_c)(q, N-1), (R)1/(N-1)), r2 = POW(x0 / FN(norm_c)(q, N), (R)1/N);
        R h = (r1 < r2 ? r1 : r2) * (R)0.1353352832366127 * (R)0.9;  /* e^-2, safety 0.9 */
        if (t + h > T) h = T - t;
        for (int i = 0; i < 12; i++) { R v = q[i][N]; for (int n = N-1; n >= 0; n--) v = v*h + q[i][n]; X[i] = v; }
        t += h; ns++;
        if (rmin) { for (int a = 0; a < 3; a++) for (int b = a+1; b < 3; b++) {
            R d = SQRT((X[2*a]-X[2*b])*(X[2*a]-X[2*b]) + (X[2*a+1]-X[2*b+1])*(X[2*a+1]-X[2*b+1])); if (d < *rmin) *rmin = d; } }
        if (!(h > 0) || ns > 20000000) return -1;
    }
    return ns;
}
static void FN(suvakov)(R p1, R p2, R *X) {
    R v[12] = {-1,0,1,0,0,0,p1,p2,p1,p2,-2*p1,-2*p2}; for (int i = 0; i < 12; i++) X[i] = v[i];
}
/* half-period symmetry residual: body 3 at the centre, outer velocities equal along the line. */
static int FN(halfres)(const R *z, R eps, R *G) {
    R X[12]; FN(suvakov)(z[0], z[1], X);
    if (FN(integrate)(X, z[2], eps, 0) < 0) return -1;
    G[0] = X[4]; G[1] = X[5];
    G[2] = (X[6]-X[8])*(X[0]-X[2]) + (X[7]-X[9])*(X[1]-X[3]);
    return 0;
}
static void FN(solve3)(R J[3][3], R *b, R *x) {
    R A[3][4]; for (int i=0;i<3;i++){ for(int j=0;j<3;j++) A[i][j]=J[i][j]; A[i][3]=b[i]; }
    for (int c=0;c<3;c++){ int p=c; for(int r=c+1;r<3;r++) if (FABS(A[r][c])>FABS(A[p][c])) p=r;
        for(int k=0;k<4;k++){ R t=A[c][k]; A[c][k]=A[p][k]; A[p][k]=t; }
        for(int r=0;r<3;r++) if(r!=c){ R f=A[r][c]/A[c][c]; for(int k=c;k<4;k++) A[r][k]-=f*A[c][k]; } }
    for (int i=0;i<3;i++) x[i]=A[i][3]/A[i][i];
}
/* Newton on the half-period residual. z = (p1, p2, T/2). Forward-difference Jacobian with step dh. */
static int FN(newton)(R *z, R eps, R dh, int maxit, R stop, double *hist, int *nit) {
    int k;
    for (k = 0; k < maxit; k++) {
        R G[3], J[3][3], dz[3];
        if (FN(halfres)(z, eps, G)) return -1;
        for (int j = 0; j < 3; j++) {
            R zz[3] = {z[0], z[1], z[2]}, Gj[3]; R d = dh * (FABS(z[j]) > 1 ? FABS(z[j]) : 1); zz[j] += d;
            if (FN(halfres)(zz, eps, Gj)) return -1;
            for (int i = 0; i < 3; i++) J[i][j] = (Gj[i] - G[i]) / d;
        }
        FN(solve3)(J, G, dz);
        R m = 0; for (int i = 0; i < 3; i++) { z[i] -= dz[i]; if (FABS(dz[i]) > m) m = FABS(dz[i]); }
        R g = FABS(G[0]); if (FABS(G[1])>g) g=FABS(G[1]); if (FABS(G[2])>g) g=FABS(G[2]);
        hist[2*k] = (double)g; hist[2*k+1] = (double)m;
        if (m < stop) { k++; break; }
    }
    *nit = k; return 0;
}
