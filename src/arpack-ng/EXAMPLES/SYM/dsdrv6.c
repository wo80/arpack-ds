/* EXAMPLES\SYM\dsdrv6.f -- translated by f2c (version 20230428). */

#include <stdlib.h>
#include "arpack_internal.h"

static a_int i_one = 1;
static a_int c__256 = 256;

void mv_(const a_int n, double *v, double *w);
void av_(const a_int n, double *v, double *w);

/**
 * \BeginDoc
 *
 *     Program to illustrate the idea of reverse communication
 *     in Cayley mode for a generalized symmetric eigenvalue
 *     problem.  The following program uses the two LAPACK subroutines
 *     dgttrf.f and dgttrs.f to factor and solve a tridiagonal system of
 *     equations.
 *
 *     We implement example six of ex-sym.doc in DOCUMENTS directory
 *
 * \Example-6
 *     ... Suppose we want to solve A*x = lambda*M*x in inverse mode,
 *         where A and M are obtained by the finite element of the
 *         1-dimensional discrete Laplacian
 *                             d^2u / dx^2
 *         on the interval [0,1] with zero Dirichlet boundary condition
 *         using piecewise linear elements.
 *
 *     ... OP = (inv[A-sigma*M])*(A+sigma*M)  and  B = M.
 *
 *     ... Use mode 5 of DSAUPD.
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * Routines called:
 *     dsaupd  ARPACK reverse communication interface routine.
 *     dseupd  ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     dgttrf  LAPACK tridiagonal factorization routine.
 *     dgttrs  LAPACK tridiagonal solve routine.
 *     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
 *     dcopy   Level 1 BLAS that copies one vector to another.
 *     dscal   Level 1 BLAS that scales a vector by a scalar.
 *     dnrm2   Level 1 BLAS that computes the norm of a vector.
 *     av      Matrix vector multiplication routine that computes A*x.
 *     mv      Matrix vector multiplication routine that computes M*x.
 *
 * \EndLib
 */
int main()
{
    /* System generated locals */
    a_int i__1;
    double d__1;

    /* Local variables */
    a_bool select[25];
    a_int iparam[11];
    a_int ipntr[11];
    a_bool rvec;
    a_int j, n, ido, ncv, nev, info, ierr = 0;
    a_int mode, nconv, ishfts, lworkl, maxitr;
    char *bmat, *which;
    double h, r1, r2, tol, sigma;

    /* Define maximum dimensions for all arrays. */

    const int MAXN   = 256; /* Maximum dimension of the A allowed. */
    const int MAXNEV =  10; /* Maximum NEV allowed */
    const int MAXNCV =  25; /* Maximum NCV allowed */

    /* ------------------------------------------------ */
    /* The number N is the dimension of the matrix. A   */
    /* generalized eigenvalue problem is solved (BMAT = */
    /* 'G'.) NEV is the number of eigenvalues to be     */
    /* approximated.  Since the Cayley mode is used,    */
    /* WHICH is set to 'LM'.  The user can modify NEV,  */
    /* NCV, SIGMA to solve problems of different sizes, */
    /* and to get different parts of the spectrum.      */
    /* However, The following conditions must be        */
    /* satisfied:                                       */
    /*                 N <= MAXN,                       */
    /*               NEV <= MAXNEV,                     */
    /*           NEV + 1 <= NCV <= MAXNCV               */
    /* ------------------------------------------------ */

    n = 100;
    nev = 4;
    ncv = 20;
    if (n > 256)
    {
        printf(" ERROR with _SDRV6: N is greater than MAXN \n");
        return ierr;
    }
    else if (nev > 10)
    {
        printf(" ERROR with _SDRV6: NEV is greater than MAXNEV \n");
        return ierr;
    }
    else if (ncv > 25)
    {
        printf(" ERROR with _SDRV6: NCV is greater than MAXNCV \n");
        return ierr;
    }
    bmat = "G";
    which = "LM";
    sigma = 150.;

    /* ------------------------------------------------ */
    /* The work array WORKL is used in DSAUPD as        */
    /* workspace.  Its dimension LWORKL is set as       */
    /* illustrated below.  The parameter TOL determines */
    /* the stopping criterion.  If TOL<=0, machine      */
    /* precision is used.  The variable IDO is used for */
    /* reverse communication and is initially set to 0. */
    /* Setting INFO=0 indicates that a random vector is */
    /* generated in DSAUPD to start the Arnoldi         */
    /* iteration.                                       */
    /* ------------------------------------------------ */

    lworkl = ncv * (ncv + 8);
    tol = 0.;
    ido = 0;
    info = 0;
    nconv = 0;

    a_int* ipiv = (a_int*)malloc(sizeof(a_int) * 256);
    double* d = (double*)malloc(sizeof(double) * 25 * 2);
    double* v = (double*)malloc(sizeof(double) * 256 * 25);
    double* ad = (double*)malloc(sizeof(double) * 256);
    double* ax = (double*)malloc(sizeof(double) * 256);
    double* mx = (double*)malloc(sizeof(double) * 256);
    double* adl = (double*)malloc(sizeof(double) * 256);
    double* adu = (double*)malloc(sizeof(double) * 256);
    double* adu2 = (double*)malloc(sizeof(double) * 256);
    double* temp = (double*)malloc(sizeof(double) * 256);
    double* resid = (double*)malloc(sizeof(double) * 256);
    double* workd = (double*)malloc(sizeof(double) * 768);
    double* workl = (double*)malloc(sizeof(double) * 825);

    /* ------------------------------------------------- */
    /* This program uses exact shifts with respect to    */
    /* the current Hessenberg matrix (IPARAM(1) = 1).    */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 5 specified in the      */
    /* documentation of DSAUPD is used (IPARAM(7) = 5).  */
    /* All these options may be changed by the user. For */
    /* details, see the documentation in DSAUPD.         */
    /* ------------------------------------------------- */

    ishfts = 1;
    maxitr = 300;
    mode = 5;

    iparam[0] = ishfts;
    iparam[2] = maxitr;
    iparam[6] = mode;

    /* ---------------------------------------------------- */
    /* Call LAPACK routine to factor (A-sigma*M).  The      */
    /* stiffness matrix A is the 1-d discrete Laplacian.    */
    /* The mass matrix M is the associated mass matrix      */
    /* arising from using piecewise linear finite elements  */
    /* on the interval [0, 1].                              */
    /* ---------------------------------------------------- */

    h = 1. / (double)(n + 1);
    r1 = h * .66666666666666663;
    r2 = h * .16666666666666666;
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
        ad[j - 1] = 2. / h - sigma * r1;
        adl[j - 1] = -1. / h - sigma * r2;
    }
    dcopy_(&n, adl, &i_one, adu, &i_one);
    dgttrf_(&n, adl, ad, adu, adu2, ipiv, &ierr);
    if (ierr != 0)
    {
        printf(" \n");
        printf(" Error with _gttrf in _SDRV6.\n");
        printf(" \n");
        return ierr;
    }

    /* ----------------------------------------- */
    /* M A I N   L O O P (Reverse communication) */
    /* ----------------------------------------- */

L10:

    /* ------------------------------------------- */
    /* Repeatedly call the routine DSAUPD and take */
    /* actions indicated by parameter IDO until    */
    /* either convergence is indicated or maxitr   */
    /* has been exceeded.                          */
    /* ------------------------------------------- */

    dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, &info);

    if (ido == -1)
    {

        /* ----------------------------------------------------- */
        /* Perform  y <--- OP*x = (inv[A-SIGMA*M])*(A+SIGMA*M)*x */
        /* to force starting vector into the range of OP.  The   */
        /* user should provide his/her matrix vector (A*x, M*x)  */
        /* multiplication routines and a linear system solver    */
        /* here.  The matrix vector multiplication routine takes */
        /* workd(ipntr(1)) as the input vector.  The final       */
        /* result is returned to workd(ipntr(2)).                */
        /* ----------------------------------------------------- */

        av_(n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
        mv_(n, &workd[ipntr[0] - 1], temp);
        daxpy_(&n, &sigma, temp, &i_one, &workd[ipntr[1] - 1], &i_one);

        dgttrs_("N", &n, &i_one, adl, ad, adu, adu2, ipiv, &workd[ipntr[1] - 1], &n, &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" Error with _gttrs in _SDRV6.\n");
            printf(" \n");
            return ierr;
        }

        /* --------------------------------------- */
        /* L O O P   B A C K to call DSAUPD again. */
        /* --------------------------------------- */

        goto L10;
    }
    else if (ido == 1)
    {

        /* -------------------------------------------------- */
        /* Perform y <-- OP*x = inv[A-SIGMA*M]*(A+SIGMA*M)*x. */
        /* M*x has been saved in workd(ipntr(3)).  The user   */
        /* need the matrix vector multiplication (A*x)        */
        /* routine and a linear system solver here.  The      */
        /* matrix vector multiplication routine takes         */
        /* workd(ipntr(1)) as the input, and the result is    */
        /* combined with workd(ipntr(3)) to form the input    */
        /* for the linear system solver.  The final result is */
        /* returned to workd(ipntr(2)).                       */
        /* -------------------------------------------------- */

        av_(n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
        daxpy_(&n, &sigma, &workd[ipntr[2] - 1], &i_one, &workd[ipntr[1] - 1], &i_one);
        dgttrs_("N", &n, &i_one, adl, ad, adu, adu2, ipiv, &workd[ipntr[1] - 1], &n, &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" Error with _gttrs in _SDRV6. \n");
            printf(" \n");
            return ierr;
        }

        /* --------------------------------------- */
        /* L O O P   B A C K to call DSAUPD again. */
        /* --------------------------------------- */

        goto L10;
    }
    else if (ido == 2)
    {

        /* ------------------------------------------ */
        /*             Perform  y <--- M*x.           */
        /* Need matrix vector multiplication routine  */
        /* here that takes workd(ipntr(1)) as input   */
        /* and returns the result to workd(ipntr(2)). */
        /* ------------------------------------------ */

        mv_(n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

        /* --------------------------------------- */
        /* L O O P   B A C K to call DSAUPD again. */
        /* --------------------------------------- */

        goto L10;
    }

    /* --------------------------------------- */
    /* Either we have convergence, or there is */
    /* an error.                               */
    /* --------------------------------------- */

    if (info < 0)
    {

        /* ------------------------ */
        /* Error message, check the */
        /* documentation in DSAUPD  */
        /* ------------------------ */

        printf(" \n");
        printf(" Error with _saupd info = %d\n", info);
        printf(" Check the documentation of _saupd. \n");
        printf(" \n");
    }
    else
    {

        /* ----------------------------------------- */
        /* No fatal errors occurred.                 */
        /* Post-Process using DSEUPD.                */
        /*                                           */
        /* Computed eigenvalues may be extracted.    */
        /*                                           */
        /* Eigenvectors may also be computed now if  */
        /* desired.  (indicated by rvec = .true.)    */
        /* ----------------------------------------- */

        rvec = TRUE_;

        dseupd_(&rvec, "A", select, d, v, &c__256, &sigma, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, &ierr);

        /* -------------------------------------------- */
        /* Eigenvalues are returned in the first column */
        /* of the two dimensional array D and the       */
        /* corresponding eigenvectors are returned in   */
        /* the first NEV columns of the two dimensional */
        /* array V if requested.  Otherwise, an         */
        /* orthogonal basis for the invariant subspace  */
        /* corresponding to the eigenvalues in D is     */
        /* returned in V.                               */
        /* -------------------------------------------- */

        if (ierr != 0)
        {

            /* ---------------------------------- */
            /* Error condition:                   */
            /* Check the documentation of DSEUPD. */
            /* ---------------------------------- */

            printf(" \n");
            printf(" Error with _seupd info = %d\n", ierr);
            printf(" Check the documentation of _seupd \n");
            printf(" \n");
        }
        else
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

            nconv = iparam[4];
            i__1 = nconv;
            for (j = 1; j <= i__1; ++j)
            {
                av_(n, &v[(j << 8) - 256], ax);
                mv_(n, &v[(j << 8) - 256], mx);
                d__1 = -d[j - 1];
                daxpy_(&n, &d__1, mx, &i_one, ax, &i_one);
                d[j + 24] = dnrm2_(&n, ax, &i_one);
                d[j + 24] /= (d__1 = d[j - 1], abs(d__1));
            }

            dmout_(nconv, 2, d, 25, -6, "Ritz values and relative residuals");
        }

        /* ---------------------------------------- */
        /* Print additional convergence information */
        /* ---------------------------------------- */

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
        printf(" _SDRV6 \n");
        printf(" ====== \n");
        printf(" \n");
        printf(" Size of the matrix is %d\n", n);
        printf(" The number of Ritz values requested is %d\n", nev);
        printf(" The number of Arnoldi vectors generated (NCV) is %d\n", ncv);
        printf(" What portion of the spectrum: %s\n", which);
        printf(" The number of converged Ritz values is %d\n", nconv);
        printf(" The number of Implicit Arnoldi update iterations taken is %d\n", iparam[2]);
        printf(" The number of OP*x is %d\n", iparam[8]);
        printf(" The convergence criterion is %e\n", tol);
        printf(" \n");
    }

    /* ------------------------- */
    /* Done with program dsdrv6. */
    /* ------------------------- */

    free(ipiv);
    free(d);
    free(v);
    free(ad);
    free(ax);
    free(mx);
    free(adl);
    free(adu);
    free(adu2);
    free(temp);
    free(resid);
    free(workd);
    free(workl);

    return nconv < nev ? EXIT_FAILURE : EXIT_SUCCESS;
}

/* ------------------------------------------------------------------------ */
/*     Matrix vector subroutine */
/*     where the matrix used is the 1 dimensional mass matrix */
/*     arising from using the piecewise linear finite element */
/*     on the interval [0,1]. */

void mv_(const a_int n, double *v, double *w)
{
    /* System generated locals */
    a_int i__1;

    /* Local variables */
    double h;
    a_int j;

    /* Parameter adjustments */
    --w;
    --v;

    w[1] = v[1] * 4. + v[2];
    i__1 = n - 1;
    for (j = 2; j <= i__1; ++j)
    {
        w[j] = v[j - 1] + v[j] * 4. + v[j + 1];
    }
    j = n;
    w[j] = v[j - 1] + v[j] * 4.;

    /*     Scale the vector w by h. */

    h = 1. / ((double)(n + 1) * 6.);
    dscal_(&n, &h, &w[1], &i_one);
} /* mv_ */

/* ------------------------------------------------------------------------ */
/*     Matrix vector subroutine */
/*     where the matrix is the stiffness matrix obtained from the */
/*     finite element discretization of the 1-dimensional discrete Laplacian */
/*     on the interval [0,1] with zero Dirichlet boundary condition */
/*     using piecewise linear elements. */

void av_(const a_int n, double *v, double *w)
{
    /* System generated locals */
    a_int i__1;
    double d__1;

    /* Local variables */
    double h;
    a_int j;

    /* Parameter adjustments */
    --w;
    --v;

    w[1] = v[1] * 2. - v[2];
    i__1 = n - 1;
    for (j = 2; j <= i__1; ++j)
    {
        w[j] = -v[j - 1] + v[j] * 2. - v[j + 1];
    }
    j = n;
    w[j] = -v[j - 1] + v[j] * 2.;

    /*     Scale the vector w by (1/h). */

    h = 1. / (double)(n + 1);
    d__1 = 1. / h;
    dscal_(&n, &d__1, &w[1], &i_one);
} /* av_ */
