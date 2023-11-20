/* EXAMPLES\SYM\dsdrv3.f -- translated by f2c (version 20230428). */

#include <stdlib.h>
#include "arpack_internal.h"

/* Table of constant values */

static a_int c__9 = 9;
static a_int c__1 = 1;
static a_int c__256 = 256;
static a_int c__3 = 3;
static a_int c__6 = 6;
static a_int c__2 = 2;
static a_int c__25 = 25;
static a_int c_n6 = -6;
static a_int c__5 = 5;

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
    a_int j, n, ido, ncv, nev, ierr, info;
    a_int mode, nconv, ishfts, lworkl, maxitr;
    char *bmat, *which;
    double h, r1, r2, tol, sigma;

    /*     Program to illustrate the idea of reverse communication in */
    /*     inverse mode for a generalized symmetric eigenvalue problem. */
    /*     The following program uses the two LAPACK subroutines dgttrf .f */
    /*     and dgttrs .f to factor and solve a tridiagonal system of equations. */

    /*     We implement example three of ex-sym.doc in DOCUMENTS directory */

    /* \Example-3 */
    /*     ... Suppose we want to solve A*x = lambda*M*x in inverse mode, */
    /*         where A and M are obtained by the finite element of the */
    /*         1-dimensional discrete Laplacian */
    /*                             d^2u / dx^2 */
    /*         on the interval [0,1] with zero Dirichlet boundary condition */
    /*         using piecewise linear elements. */

    /*     ... OP = inv[M]*A  and  B = M. */

    /*     ... Use mode 2 of DSAUPD . */

    /* \BeginLib */

    /* \Routines called: */
    /*     dsaupd   ARPACK reverse communication interface routine. */
    /*     dseupd   ARPACK routine that returns Ritz values and (optionally) */
    /*             Ritz vectors. */
    /*     dgttrf   LAPACK tridiagonal factorization routine. */
    /*     dgttrs   LAPACK tridiagonal solve routine. */
    /*     daxpy    Level 1 BLAS that computes y <- alpha*x+y. */
    /*     dscal    Level 1 BLAS that scales a vector by a scalar. */
    /*     dcopy    Level 1 BLAS that copies one vector to another. */
    /*     dnrm2    Level 1 BLAS that computes the norm of a vector. */
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
    /* FILE: sdrv3.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2 */

    /* \Remarks */
    /*     1. None */

    /* \EndLib */
    /* -------------------------------------------------------------------------- */

    /* --------------------------- */
    /* Define leading dimensions   */
    /* for all arrays.             */
    /* MAXN:   Maximum dimension   */
    /*         of the A allowed.   */
    /* MAXNEV: Maximum NEV allowed */
    /* MAXNCV: Maximum NCV allowed */
    /* --------------------------- */

    /* -------------------------------------------------- */
    /* The number N is the dimension of the matrix. A     */
    /* generalized eigenvalue problem is solved (BMAT =   */
    /* 'G'.) NEV is the number of eigenvalues to be       */
    /* approximated.  The user can modify NEV, NCV, WHICH */
    /* to solve problems of different sizes, and to get   */
    /* different parts of the spectrum.  However, The     */
    /* following conditions must be satisfied:            */
    /*                     N <= MAXN,                     */
    /*                   NEV <= MAXNEV,                   */
    /*               NEV + 1 <= NCV <= MAXNCV             */
    /* -------------------------------------------------- */

    n = 100;
    nev = 4;
    ncv = 10;
    if (n > 256)
    {
        printf(" ERROR with _SDRV3: N is greater than MAXN \n");
        return ierr;
    }
    else if (nev > 10)
    {
        printf(" ERROR with _SDRV3: NEV is greater than MAXNEV \n");
        return ierr;
    }
    else if (ncv > 25)
    {
        printf(" ERROR with _SDRV3: NCV is greater than MAXNCV \n");
        return ierr;
    }
    bmat = "G";
    which = "LM";

    /* ------------------------------------------------ */
    /* The work array WORKL is used in DSAUPD  as        */
    /* workspace.  Its dimension LWORKL is set as       */
    /* illustrated below.  The parameter TOL determines */
    /* the stopping criterion.  If TOL<=0, machine      */
    /* precision is used.  The variable IDO is used for */
    /* reverse communication and is initially set to 0. */
    /* Setting INFO=0 indicates that a random vector is */
    /* generated in DSAUPD  to start the Arnoldi         */
    /* iteration.                                       */
    /* ------------------------------------------------ */

    lworkl = ncv * (ncv + 8);
    tol = 0.;
    ido = 0;
    info = 0;

    a_int* ipiv = (a_int*)malloc(sizeof(a_int) * 256);
    double* d = (double*)malloc(sizeof(double) * 25 * 2);
    double* v = (double*)malloc(sizeof(double) * 256 * 25);
    double* ad = (double*)malloc(sizeof(double) * 256);
    double* ax = (double*)malloc(sizeof(double) * 256);
    double* mx = (double*)malloc(sizeof(double) * 256);
    double* adl = (double*)malloc(sizeof(double) * 256);
    double* adu = (double*)malloc(sizeof(double) * 256);
    double* adu2 = (double*)malloc(sizeof(double) * 256);
    double* resid = (double*)malloc(sizeof(double) * 256);
    double* workd = (double*)malloc(sizeof(double) * 768);
    double* workl = (double*)malloc(sizeof(double) * 825);

    /* ------------------------------------------------- */
    /* This program uses exact shifts with respect to    */
    /* the current Hessenberg matrix (IPARAM(1) = 1).    */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 2 of DSAUPD  is used     */
    /* (IPARAM(7) = 2).  All these options may be        */
    /* changed by the user. For details, see the         */
    /* documentation in DSAUPD .                          */
    /* ------------------------------------------------- */

    ishfts = 1;
    maxitr = 300;
    mode = 2;

    iparam[0] = ishfts;
    iparam[2] = maxitr;
    iparam[6] = mode;

    /* ---------------------------------------------- */
    /* Call LAPACK routine to factor the mass matrix. */
    /* The mass matrix is the tridiagonal matrix      */
    /* arising from using piecewise linear finite     */
    /* elements on the interval [0, 1].               */
    /* ---------------------------------------------- */

    h = 1. / (double)(n + 1);

    r1 = h * .66666666666666663;
    r2 = h * .16666666666666666;
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
        ad[j - 1] = r1;
        adl[j - 1] = r2;
        /* L20: */
    }
    dcopy_(&n, adl, &c__1, adu, &c__1);
    dgttrf_(&n, adl, ad, adu, adu2, ipiv, &ierr);
    if (ierr != 0)
    {
        printf(" \n");
        printf(" Error with _gttrf in _SDRV3. \n");
        printf(" \n");
        return ierr;
    }

    /* ----------------------------------------- */
    /* M A I N   L O O P (Reverse communication) */
    /* ----------------------------------------- */

L10:

    /* ------------------------------------------- */
    /* Repeatedly call the routine DSAUPD  and take */
    /* actions indicated by parameter IDO until    */
    /* either convergence is indicated or maxitr   */
    /* has been exceeded.                          */
    /* ------------------------------------------- */

    dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, &info);

    if (ido == -1 || ido == 1)
    {

        /* ------------------------------------ */
        /* Perform  y <--- OP*x = inv[M]*A*x.   */
        /* The user should supply his/her own   */
        /* matrix vector multiplication (A*x)   */
        /* routine and a linear system solver   */
        /* here.  The matrix vector             */
        /* multiplication routine takes         */
        /* workd(ipntr(1)) as the input vector. */
        /* The final result is returned to      */
        /* workd(ipntr(2)). The result of A*x   */
        /* overwrites workd(ipntr(1)).          */
        /* ------------------------------------ */

        av_(&n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
        dcopy_(&n, &workd[ipntr[1] - 1], &c__1, &workd[ipntr[0] - 1], &c__1);
        dgttrs_("Notranspose", &n, &c__1, adl, ad, adu, adu2, ipiv, &workd[ipntr[1] - 1], &n, &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" Error with _gttrs in _SDRV3.\n");
            printf(" \n");
            return ierr;
        }

        /* --------------------------------------- */
        /* L O O P   B A C K to call DSAUPD  again. */
        /* --------------------------------------- */

        goto L10;
    }
    else if (ido == 2)
    {

        /* --------------------------------------- */
        /*         Perform  y <--- M*x.            */
        /* Need the matrix vector multiplication   */
        /* routine here that takes workd(ipntr(1)) */
        /* as the input and returns the result to  */
        /* workd(ipntr(2)).                        */
        /* --------------------------------------- */

        mv_(&n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

        /* --------------------------------------- */
        /* L O O P   B A C K to call DSAUPD  again. */
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
        /* documentation in DSAUPD   */
        /* ------------------------ */

        printf(" \n");
        printf(" Error with _saupd info = %d", info);
        printf(" Check the documentation of _saupd \n");
        printf(" \n");
    }
    else
    {

        /* ----------------------------------------- */
        /* No fatal errors occurred.                 */
        /* Post-Process using DSEUPD .                */
        /*                                           */
        /* Computed eigenvalues may be extracted.    */
        /*                                           */
        /* Eigenvectors may also be computed now if  */
        /* desired.  (indicated by rvec = .true.)    */
        /* ----------------------------------------- */

        rvec = TRUE_;

        dseupd_(&rvec, "All", select, d, v, &c__256, &sigma, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, &ierr);

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
            /* Check the documentation of DSEUPD . */
            /* ---------------------------------- */

            printf(" \n");
            printf(" Error with _seupd info = %d", ierr);
            printf(" Check the documentation of _seupd\n");
            printf(" \n");
        }
        else
        {

            nconv = iparam[4];
            i__1 = nconv;
            for (j = 1; j <= i__1; ++j)
            {

                /* ------------------------- */
                /* Compute the residual norm */
                /*                           */
                /*  ||  A*x - lambda*M*x ||  */
                /*                           */
                /* for the NCONV accurately  */
                /* computed eigenvalues and  */
                /* eigenvectors.  (iparam(5) */
                /* indicates how many are    */
                /* accurate to the requested */
                /* tolerance)                */
                /* ------------------------- */

                av_(&n, &v[(j << 8) - 256], ax);
                mv_(&n, &v[(j << 8) - 256], mx);
                d__1 = -d[j - 1];
                daxpy_(&n, &d__1, mx, &c__1, ax, &c__1);
                d[j + 24] = dnrm2_(&n, ax, &c__1);
                d[j + 24] /= (d__1 = d[j - 1], abs(d__1));

                /* L30: */
            }

            /* --------------------------- */
            /* Display computed residuals. */
            /* --------------------------- */

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
        printf(" _SDRV3 \n");
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
    /* Done with program dsdrv3 . */
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
    free(resid);
    free(workd);
    free(workl);

    return 0;
}

/* ------------------------------------------------------------------------ */
/*     Matrix vector subroutine */
/*     where the matrix is the 1 dimensional mass matrix */
/*     on the interval [0,1]. */

int mv_(a_int *n, double *v, double *w)
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
    i__1 = *n - 1;
    for (j = 2; j <= i__1; ++j)
    {
        w[j] = v[j - 1] + v[j] * 4. + v[j + 1];
        /* L100: */
    }
    j = *n;
    w[j] = v[j - 1] + v[j] * 4.;

    /*     Scale the vector w by h. */

    h = 1. / ((double)(*n + 1) * 6.);
    dscal_(n, &h, &w[1], &c__1);
    return 0;
} /* mv_ */

/* -------------------------------------------------------------------- */
/*     matrix vector subroutine */

/*     The matrix used is the stiffness matrix obtained from the finite */
/*     element discretization of the 1-dimensional discrete Laplacian */
/*     on the interval [0,1] with zero Dirichlet boundary condition using */
/*     piecewise linear elements. */

int av_(a_int *n, double *v, double *w)
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
    i__1 = *n - 1;
    for (j = 2; j <= i__1; ++j)
    {
        w[j] = -v[j - 1] + v[j] * 2. - v[j + 1];
        /* L100: */
    }
    j = *n;
    w[j] = -v[j - 1] + v[j] * 2.;

    /*     Scale the vector w by (1 / h). */

    h = 1. / (double)(*n + 1);
    d__1 = 1. / h;
    dscal_(n, &d__1, &w[1], &c__1);
    return 0;
} /* av_ */
