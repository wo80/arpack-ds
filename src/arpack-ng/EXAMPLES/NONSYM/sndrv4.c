/* EXAMPLES\NONSYM\sndrv4.f -- translated by f2c (version 20230428). */

#include <stdlib.h>
#include "arpack_internal.h"

struct
{
    float rho;
} convct_;

#define convct_1 convct_

/* Table of constant values */

static a_int c__9 = 9;
static a_int c__1 = 1;
static a_int c__256 = 256;
static a_int c__3 = 3;
static a_int c__6 = 6;
static a_int c__25 = 25;
static a_int c_n6 = -6;
static a_int c__4 = 4;

int main()
{
    /* System generated locals */
    a_int i__1;
    float r__1;

    /* Local variables */
    a_bool select[25];
    a_int iparam[11];
    a_int ipntr[14];
    a_bool rvec, first;
    a_int j, n, ido, ncv, nev, ierr, info;
    a_int mode, nconv, ishfts, lworkl, maxitr;
    char *bmat, *which;
    float h, s, s1, s2, s3, tol, sigmai;
    float sigmar;

    /*     Simple program to illustrate the idea of reverse communication */
    /*     in shift-invert mode for a generalized nonsymmetric eigenvalue */
    /*     problem. */

    /*     We implement example four of ex-nonsym.doc in DOCUMENTS directory */

    /* \Example-4 */
    /*     ... Suppose we want to solve A*x = lambda*B*x in inverse mode, */
    /*         where A and B are derived from the finite element discretization */
    /*         of the 1-dimensional convection-diffusion operator */
    /*                           (d^2u / dx^2) + rho*(du/dx) */
    /*         on the interval [0,1] with zero Dirichlet boundary condition */
    /*         using linear elements. */

    /*     ... The shift sigma is a real number. */

    /*     ... OP = inv[A-SIGMA*M]*M  and  B = M. */

    /*     ... Use mode 3 of SNAUPD. */

    /* \BeginLib */

    /* \Routines called: */
    /*     snaupd  ARPACK reverse communication interface routine. */
    /*     sneupd  ARPACK routine that returns Ritz values and (optionally) */
    /*             Ritz vectors. */
    /*     sgttrf  LAPACK tridiagonal factorization routine. */
    /*     sgttrs  LAPACK tridiagonal linear system solve routine. */
    /*     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully. */
    /*     saxpy   Level 1 BLAS that computes y <- alpha*x+y. */
    /*     scopy   Level 1 BLAS that copies one vector to another. */
    /*     sdot    Level 1 BLAS that computes the dot product of two vectors. */
    /*     snrm2   Level 1 BLAS that computes the norm of a vector. */
    /*     av      Matrix vector multiplication routine that computes A*x. */
    /*     mv      Matrix vector multiplication routine that computes M*x. */

    /* \Author */
    /*     Richard Lehoucq */
    /*     Danny Sorensen */
    /*     Chao Yang */
    /*     Dept. of Computational & */
    /*     Applied Mathematics */
    /*     Rice University */
    /*     Houston, Texas */

    /* \SCCS Information: @(#) */
    /* FILE: ndrv4.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2 */

    /* \Remarks */
    /*     1. None */

    /* \EndLib */
    /* ----------------------------------------------------------------------- */

    /* --------------------------- */
    /* Define leading dimensions   */
    /* for all arrays.             */
    /* MAXN:   Maximum dimension   */
    /*         of the A allowed.   */
    /* MAXNEV: Maximum NEV allowed */
    /* MAXNCV: Maximum NCV allowed */
    /* --------------------------- */

    /* -------------------------------------------------- */
    /* The number N is the dimension of the matrix.  A    */
    /* generalized eigenvalue problem is solved (BMAT =   */
    /* 'G').  NEV is the number of eigenvalues (closest   */
    /* to SIGMAR) to be approximated.  Since the          */
    /* shift-invert mode is used,  WHICH is set to 'LM'.  */
    /* The user can modify NEV, NCV, SIGMAR to solve      */
    /* problems of different sizes, and to get different  */
    /* parts of the spectrum.  However, The following     */
    /* conditions must be satisfied:                      */
    /*                     N <= MAXN,                     */
    /*                   NEV <= MAXNEV,                   */
    /*               NEV + 2 <= NCV <= MAXNCV             */
    /* -------------------------------------------------- */

    n = 100;
    nev = 4;
    ncv = 10;
    if (n > 256)
    {
        printf(" ERROR with _NDRV4: N is greater than MAXN \n");
        return ierr;
    }
    else if (nev > 10)
    {
        printf(" ERROR with _NDRV4: NEV is greater than MAXNEV \n");
        return ierr;
    }
    else if (ncv > 25)
    {
        printf(" ERROR with _NDRV4: NCV is greater than MAXNCV \n");
        return ierr;
    }
    bmat = "G";
    which = "LM";
    sigmar = 1.f;
    sigmai = 0.f;

    /* ------------------------------------------------ */
    /* Construct C = A - SIGMA*M in real arithmetic,    */
    /* and factor C in real arithmetic (using LAPACK    */
    /* subroutine sgttrf). The matrix A is chosen to be */
    /* the tridiagonal matrix derived from the standard */
    /* central difference discretization of the 1-d     */
    /* convection-diffusion operator u" + rho*u' on the */
    /* interval [0, 1] with zero Dirichlet boundary     */
    /* condition.  The matrix M is the mass matrix      */
    /* formed by using piecewise linear elements on     */
    /* [0,1].                                           */
    /* ------------------------------------------------ */

    a_int* ipiv = (a_int*)malloc(sizeof(a_int) * 256);

    float* dd = (float*)malloc(sizeof(float) * 256);
    float* dl = (float*)malloc(sizeof(float) * 256);
    float* du = (float*)malloc(sizeof(float) * 256);
    float* du2 = (float*)malloc(sizeof(float) * 256);

    convct_1.rho = 10.f;
    h = 1.f / (float)(n + 1);
    s = convct_1.rho / 2.f;

    s1 = -1.f / h - s - sigmar * h / 6.f;
    s2 = 2.f / h - sigmar * 4.f * h / 6.f;
    s3 = -1.f / h + s - sigmar * h / 6.f;

    i__1 = n - 1;
    for (j = 1; j <= i__1; ++j)
    {
        dl[j - 1] = s1;
        dd[j - 1] = s2;
        du[j - 1] = s3;
        /* L10: */
    }
    dd[n - 1] = s2;

    sgttrf_(&n, dl, dd, du, du2, ipiv, &ierr);
    if (ierr != 0)
    {
        printf(" \n");
        printf(" ERROR with _gttrf in _NDRV4.\n");
        printf(" \n");
        return ierr;
    }

    /* --------------------------------------------------- */
    /* The work array WORKL is used in SNAUPD as           */
    /* workspace.  Its dimension LWORKL is set as          */
    /* illustrated below.  The parameter TOL determines    */
    /* the stopping criterion. If TOL<=0, machine          */
    /* precision is used.  The variable IDO is used for    */
    /* reverse communication, and is initially set to 0.   */
    /* Setting INFO=0 indicates that a random vector is    */
    /* generated in SNAUPD to start the Arnoldi iteration. */
    /* --------------------------------------------------- */

    lworkl = ncv * ncv * 3 + ncv * 6;
    tol = 0.f;
    ido = 0;
    info = 0;

    float* d = (float*)malloc(sizeof(float) * 25 * 3);
    float* v = (float*)malloc(sizeof(float) * 256 * 25);
    float* ax = (float*)malloc(sizeof(float) * 256);
    float* mx = (float*)malloc(sizeof(float) * 256);
    float* resid = (float*)malloc(sizeof(float) * 256);
    float* workd = (float*)malloc(sizeof(float) * 768);
    float* workl = (float*)malloc(sizeof(float) * 2025);
    float* workev = (float*)malloc(sizeof(float) * 75);

    /* ------------------------------------------------- */
    /* This program uses exact shifts with respect to    */
    /* the current Hessenberg matrix (IPARAM(1) = 1).    */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 3 of SNAUPD is used     */
    /* (IPARAM(7) = 3).  All these options can be        */
    /* changed by the user. For details, see the         */
    /* documentation in SNAUPD.                          */
    /* ------------------------------------------------- */

    ishfts = 1;
    maxitr = 300;
    mode = 3;

    iparam[0] = ishfts;
    iparam[2] = maxitr;
    iparam[6] = mode;

    /* ---------------------------------------- */
    /* M A I N   L O O P(Reverse communication) */
    /* ---------------------------------------- */

L20:

    /* ------------------------------------------- */
    /* Repeatedly call the routine SNAUPD and take */
    /* actions indicated by parameter IDO until    */
    /* either convergence is indicated or maxitr   */
    /* has been exceeded.                          */
    /* ------------------------------------------- */

    snaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, &info);

    if (ido == -1)
    {

        /* ----------------------------------------- */
        /* Perform  y <--- OP*x = inv[A-SIGMA*M]*M*x */
        /* to force starting vector into the range   */
        /* of OP.   The user should supply his/her   */
        /* own matrix vector multiplication routine  */
        /* and a linear system solver.  The matrix   */
        /* vector multiplication routine should take */
        /* workd(ipntr(1)) as the input. The final   */
        /* result should be returned to              */
        /* workd(ipntr(2)).                          */
        /* ----------------------------------------- */

        mv_(&n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
        sgttrs_("N", &n, &c__1, dl, dd, du, du2, ipiv, &workd[ipntr[1] - 1], &n, &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" ERROR with _gttrs in _NDRV4.\n");
            printf(" \n");
            return ierr;
        }

        /* --------------------------------------- */
        /* L O O P   B A C K to call SNAUPD again. */
        /* --------------------------------------- */

        goto L20;
    }
    else if (ido == 1)
    {

        /* --------------------------------------- */
        /* Perform y <-- OP*x = inv[A-sigma*M]*M*x */
        /* M*x has been saved in workd(ipntr(3)).  */
        /* The user only need the linear system    */
        /* solver here that takes workd(ipntr(3))  */
        /* as input, and returns the result to     */
        /* workd(ipntr(2)).                        */
        /* --------------------------------------- */

        scopy_(&n, &workd[ipntr[2] - 1], &c__1, &workd[ipntr[1] - 1], &c__1);
        sgttrs_("N", &n, &c__1, dl, dd, du, du2, ipiv, &workd[ipntr[1] - 1], &n, &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" ERROR with _gttrs in _NDRV4.\n");
            printf(" \n");
            return ierr;
        }

        /* --------------------------------------- */
        /* L O O P   B A C K to call SNAUPD again. */
        /* --------------------------------------- */

        goto L20;
    }
    else if (ido == 2)
    {

        /* ------------------------------------------- */
        /*          Perform  y <--- M*x                */
        /* Need matrix vector multiplication routine   */
        /* here that takes workd(ipntr(1)) as input    */
        /* and returns the result to workd(ipntr(2)).  */
        /* ------------------------------------------- */

        mv_(&n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

        /* --------------------------------------- */
        /* L O O P   B A C K to call SNAUPD again. */
        /* --------------------------------------- */

        goto L20;
    }

    /* --------------------------------------- */
    /* Either we have convergence, or there is */
    /* an error.                               */
    /* --------------------------------------- */

    if (info < 0)
    {

        /* ------------------------ */
        /* Error message, check the */
        /* documentation in SNAUPD. */
        /* ------------------------ */

        printf(" \n");
        printf(" Error with _naupd info = %d", info);
        printf(" Check the documentation in _naupd.\n");
        printf(" \n");
    }
    else
    {

        /* ----------------------------------------- */
        /* No fatal errors occurred.                 */
        /* Post-Process using SNEUPD.                */
        /*                                           */
        /* Computed eigenvalues may be extracted.    */
        /*                                           */
        /* Eigenvectors may also be computed now if  */
        /* desired.  (indicated by rvec = .true.)    */
        /* ----------------------------------------- */

        rvec = TRUE_;
        sneupd_(&rvec, "A", select, d, &d[25], v, &c__256, &sigmar, &sigmai, workev, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, &ierr);

        /* --------------------------------------------- */
        /* The real part of the eigenvalue is returned   */
        /* in the first column of the two dimensional    */
        /* array D, and the IMAGINARY part is returned   */
        /* in the second column of D.  The corresponding */
        /* eigenvectors are returned in the first NEV    */
        /* columns of the two dimensional array V if     */
        /* requested.  Otherwise, an orthogonal basis    */
        /* for the invariant subspace corresponding to   */
        /* the eigenvalues in D is returned in V.        */
        /* --------------------------------------------- */

        if (ierr != 0)
        {

            /* ---------------------------------- */
            /* Error condition:                   */
            /* Check the documentation of SNEUPD. */
            /* ---------------------------------- */

            printf(" \n");
            printf(" Error with _neupd info = %d", ierr);
            printf(" Check the documentation of _neupd. \n");
            printf(" \n");
        }
        else
        {

            first = TRUE_;
            nconv = iparam[4];
            i__1 = nconv;
            for (j = 1; j <= i__1; ++j)
            {

                /* ------------------------- */
                /* Compute the residual norm */
                /*                           */
                /*   ||  A*x - lambda*x ||   */
                /*                           */
                /* for the NCONV accurately  */
                /* computed eigenvalues and  */
                /* eigenvectors.  (iparam(5) */
                /* indicates how many are    */
                /* accurate to the requested */
                /* tolerance)                */
                /* ------------------------- */

                if (d[j + 24] == 0.f)
                {

                    /* ------------------ */
                    /* Ritz value is real */
                    /* ------------------ */

                    av_(&n, &v[(j << 8) - 256], ax);
                    mv_(&n, &v[(j << 8) - 256], mx);
                    r__1 = -d[j - 1];
                    saxpy_(&n, &r__1, mx, &c__1, ax, &c__1);
                    d[j + 49] = snrm2_(&n, ax, &c__1);
                    d[j + 49] /= (r__1 = d[j - 1], dabs(r__1));
                }
                else if (first)
                {

                    /* ---------------------- */
                    /* Ritz value is complex. */
                    /* Residual of one Ritz   */
                    /* value of the conjugate */
                    /* pair is computed.      */
                    /* ---------------------- */

                    av_(&n, &v[(j << 8) - 256], ax);
                    mv_(&n, &v[(j << 8) - 256], mx);
                    r__1 = -d[j - 1];
                    saxpy_(&n, &r__1, mx, &c__1, ax, &c__1);
                    mv_(&n, &v[(j + 1 << 8) - 256], mx);
                    saxpy_(&n, &d[j + 24], mx, &c__1, ax, &c__1);
                    d[j + 49] = snrm2_(&n, ax, &c__1);
                    av_(&n, &v[(j + 1 << 8) - 256], ax);
                    mv_(&n, &v[(j + 1 << 8) - 256], mx);
                    r__1 = -d[j - 1];
                    saxpy_(&n, &r__1, mx, &c__1, ax, &c__1);
                    mv_(&n, &v[(j << 8) - 256], mx);
                    r__1 = -d[j + 24];
                    saxpy_(&n, &r__1, mx, &c__1, ax, &c__1);
                    r__1 = snrm2_(&n, ax, &c__1);
                    d[j + 49] = slapy2_(&d[j + 49], &r__1);
                    d[j + 49] /= slapy2_(&d[j - 1], &d[j + 24]);
                    d[j + 50] = d[j + 49];
                    first = FALSE_;
                }
                else
                {
                    first = TRUE_;
                }

                /* L30: */
            }

            /* --------------------------- */
            /* Display computed residuals. */
            /* --------------------------- */

            smout_(nconv, 3, d, 25, -6, "Ritz values (Real,Imag) and relative residuals");
        }

        /* ----------------------------------------- */
        /* Print additional convergence information. */
        /* ----------------------------------------- */

        if (info == 1)
        {
            printf(" \n");
            printf(" Maximum number of iterations reached.\n");
            printf(" \n");
        }
        else if (info == 3)
        {
            printf(" \n");
            printf(" No shifts could be applied during implicit\n");
            printf(" Arnoldi update try increasing NCV.\n");
            printf(" \n");
        }

        printf(" \n");
        printf(" _NDRV4 \n");
        printf(" ====== \n");
        printf(" \n");
        printf(" Size of the matrix is %d", n);
        printf(" The number of Ritz values requested is %d", nev);
        printf(" The number of Arnoldi vectors generated (NCV) is %d", ncv);
        printf(" What portion of the spectrum: %s", which);
        printf(" The number of converged Ritz values is %d", nconv);
        printf(" The number of Implicit Arnoldi update iterations taken is %d", iparam[2]);
        printf(" The number of OP*x is %d", iparam[8]);
        printf(" The convergence criterion is %e", tol);
        printf(" \n");
    }

    /* ------------------------- */
    /* Done with program sndrv4. */
    /* ------------------------- */

    free(ipiv);
    free(d);
    free(v);
    free(ax);
    free(dd);
    free(dl);
    free(du);
    free(mx);
    free(du2);
    free(resid);
    free(workd);
    free(workl);
    free(workev);

    return 0;
}

/* ========================================================================== */

/*     matrix vector multiplication subroutine */

int mv_(a_int *n, float *v, float *w)
{
    /* System generated locals */
    a_int i__1;

    /* Local variables */
    float h;
    a_int j;

    /*     Compute the matrix vector multiplication y<---M*x */
    /*     where M is mass matrix formed by using piecewise linear elements */
    /*     on [0,1]. */

    /* Parameter adjustments */
    --w;
    --v;

    w[1] = (v[1] * 4.f + v[2] * 1.f) / 6.f;
    i__1 = *n - 1;
    for (j = 2; j <= i__1; ++j)
    {
        w[j] = (v[j - 1] * 1.f + v[j] * 4.f + v[j + 1] * 1.f) / 6.f;
        /* L10: */
    }
    w[*n] = (v[*n - 1] * 1.f + v[*n] * 4.f) / 6.f;

    h = 1.f / (float)(*n + 1);
    sscal_(n, &h, &w[1], &c__1);
    return 0;
} /* mv_ */

/* ------------------------------------------------------------------ */
int av_(a_int *n, float *v, float *w)
{
    /* System generated locals */
    a_int i__1;

    /* Local variables */
    float h;
    a_int j;
    float s, dd, dl, du;

    /*     Compute the matrix vector multiplication y<---A*x */
    /*     where A is obtained from the finite element discretization of the */
    /*     1-dimensional convection diffusion operator */
    /*                     d^u/dx^2 + rho*(du/dx) */
    /*     on the interval [0,1] with zero Dirichlet boundary condition */
    /*     using linear elements. */
    /*     This routine is only used in residual calculation. */

    /* Parameter adjustments */
    --w;
    --v;

    h = 1.f / (float)(*n + 1);
    s = convct_1.rho / 2.f;
    dd = 2.f / h;
    dl = -1.f / h - s;
    du = -1.f / h + s;

    w[1] = dd * v[1] + du * v[2];
    i__1 = *n - 1;
    for (j = 2; j <= i__1; ++j)
    {
        w[j] = dl * v[j - 1] + dd * v[j] + du * v[j + 1];
        /* L10: */
    }
    w[*n] = dl * v[*n - 1] + dd * v[*n];
    return 0;
} /* av_ */
