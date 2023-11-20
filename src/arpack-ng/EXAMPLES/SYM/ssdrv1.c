/* EXAMPLES\SYM\ssdrv1.f -- translated by f2c (version 20230428). */

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
static a_int c__4 = 4;
static float c_b138 = -1.f;

int main()
{
    /* System generated locals */
    a_int i__1;
    float r__1;

    /* Local variables */
    a_bool select[25];
    a_int iparam[11];
    a_int ipntr[11];
    a_bool rvec;
    a_int j, n, nx, ido, ncv, nev, ierr = 0;
    a_int info, mode, nconv, ishfts, lworkl, maxitr;
    char *bmat, *which;
    float tol, sigma;

    /*     Simple program to illustrate the idea of reverse communication */
    /*     in regular mode for a standard symmetric eigenvalue problem. */

    /*     We implement example one of ex-sym.doc in SRC directory */

    /* \Example-1 */
    /*     ... Suppose we want to solve A*x = lambda*x in regular mode, */
    /*         where A is derived from the central difference discretization */
    /*         of the 2-dimensional Laplacian on the unit square [0,1]x[0,1] */
    /*         with zero Dirichlet boundary condition. */

    /*     ... OP = A  and  B = I. */

    /*     ... Assume "call av (n,x,y)" computes y = A*x. */

    /*     ... Use mode 1 of SSAUPD. */

    /* \BeginLib */

    /* \Routines called: */
    /*     ssaupd  ARPACK reverse communication interface routine. */
    /*     sseupd  ARPACK routine that returns Ritz values and (optionally) */
    /*             Ritz vectors. */
    /*     snrm2   Level 1 BLAS that computes the norm of a vector. */
    /*     saxpy   Level 1 BLAS that computes y <- alpha*x+y. */
    /*     av      Matrix vector multiplication routine that computes A*x. */
    /*     tv      Matrix vector multiplication routine that computes T*x, */
    /*             where T is a tridiagonal matrix.  It is used in routine */
    /*             av. */

    /* \Author */
    /*     Richard Lehoucq */
    /*     Danny Sorensen */
    /*     Chao Yang */
    /*     Dept. of Computational & */
    /*     Applied Mathematics */
    /*     Rice University */
    /*     Houston, Texas */

    /* \SCCS Information: @(#) */
    /* FILE: sdrv1.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2 */

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
    /* The number NX is the number of interior points     */
    /* in the discretization of the 2-dimensional         */
    /* Laplacian on the unit square with zero Dirichlet   */
    /* boundary condition.  The number N(=NX*NX) is the   */
    /* dimension of the matrix.  A standard eigenvalue    */
    /* problem is solved (BMAT = 'I'). NEV is the number  */
    /* of eigenvalues to be approximated.  The user can   */
    /* modify NEV, NCV, WHICH to solve problems of        */
    /* different sizes, and to get different parts of the */
    /* spectrum.  However, The following conditions must  */
    /* be satisfied:                                      */
    /*                   N <= MAXN,                       */
    /*                 NEV <= MAXNEV,                     */
    /*             NEV + 1 <= NCV <= MAXNCV               */
    /* -------------------------------------------------- */

    nx = 10;
    n = nx * nx;
    nev = 4;
    ncv = 10;
    if (n > 256)
    {
        printf(" ERROR with _SDRV1: N is greater than MAXN \n");
        return ierr;
    }
    else if (nev > 10)
    {
        printf(" ERROR with _SDRV1: NEV is greater than MAXNEV \n");
        return ierr;
    }
    else if (ncv > 25)
    {
        printf(" ERROR with _SDRV1: NCV is greater than MAXNCV \n");
        return ierr;
    }
    bmat = "I";
    which = "SM";

    /* ------------------------------------------------ */
    /* The work array WORKL is used in SSAUPD as        */
    /* workspace.  Its dimension LWORKL is set as       */
    /* illustrated below.  The parameter TOL determines */
    /* the stopping criterion.  If TOL<=0, machine      */
    /* precision is used.  The variable IDO is used for */
    /* reverse communication and is initially set to 0. */
    /* Setting INFO=0 indicates that a random vector is */
    /* generated in SSAUPD to start the Arnoldi         */
    /* iteration.                                       */
    /* ------------------------------------------------ */

    lworkl = ncv * (ncv + 8);
    tol = 0.f;
    info = 0;
    ido = 0;

    float* d = (float*)malloc(sizeof(float) * 25 * 2);
    float* v = (float*)malloc(sizeof(float) * 256 * 25);
    float* ax = (float*)malloc(sizeof(float) * 256);
    float* resid = (float*)malloc(sizeof(float) * 256);
    float* workd = (float*)malloc(sizeof(float) * 768);
    float* workl = (float*)malloc(sizeof(float) * 825);

    /* ------------------------------------------------- */
    /* This program uses exact shifts with respect to    */
    /* the current Hessenberg matrix (IPARAM(1) = 1).    */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 1 of SSAUPD is used     */
    /* (IPARAM(7) = 1).  All these options may be        */
    /* changed by the user. For details, see the         */
    /* documentation in SSAUPD.                          */
    /* ------------------------------------------------- */

    ishfts = 1;
    maxitr = 300;
    mode = 1;

    iparam[0] = ishfts;
    iparam[2] = maxitr;
    iparam[6] = mode;

    /* ----------------------------------------- */
    /* M A I N   L O O P (Reverse communication) */
    /* ----------------------------------------- */

L10:

    /* ------------------------------------------- */
    /* Repeatedly call the routine SSAUPD and take */
    /* actions indicated by parameter IDO until    */
    /* either convergence is indicated or maxitr   */
    /* has been exceeded.                          */
    /* ------------------------------------------- */

    ssaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, &info);

    if (ido == -1 || ido == 1)
    {

        /* ------------------------------------ */
        /* Perform matrix vector multiplication */
        /*              y <--- OP*x             */
        /* The user should supply his/her own   */
        /* matrix vector multiplication routine */
        /* here that takes workd(ipntr(1)) as   */
        /* the input, and return the result to  */
        /* workd(ipntr(2)).                     */
        /* ------------------------------------ */

        av_(&nx, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

        /* --------------------------------------- */
        /* L O O P   B A C K to call SSAUPD again. */
        /* --------------------------------------- */

        goto L10;
    }

    /* -------------------------------------- */
    /* Either we have convergence or there is */
    /* an error.                              */
    /* -------------------------------------- */

    if (info < 0)
    {

        /* ------------------------ */
        /* Error message. Check the */
        /* documentation in SSAUPD. */
        /* ------------------------ */

        printf(" \n");
        printf(" Error with _saupd info = %d", info);
        printf(" Check documentation in _saupd \n");
        printf(" \n");
    }
    else
    {

        /* ----------------------------------------- */
        /* No fatal errors occurred.                 */
        /* Post-Process using SSEUPD.                */
        /*                                           */
        /* Computed eigenvalues may be extracted.    */
        /*                                           */
        /* Eigenvectors may also be computed now if  */
        /* desired.  (indicated by rvec = .true.)    */
        /* ----------------------------------------- */

        rvec = TRUE_;

        sseupd_(&rvec, "All", select, d, v, &c__256, &sigma, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, &ierr);
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
            /* Check the documentation of SSEUPD. */
            /* ---------------------------------- */

            printf(" \n");
            printf(" Error with _seupd info = %d", ierr);
            printf(" Check the documentation of _seupd. \n");
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
                /*   ||  A*x - lambda*x ||   */
                /*                           */
                /* for the NCONV accurately  */
                /* computed eigenvalues and  */
                /* eigenvectors.  (iparam(5) */
                /* indicates how many are    */
                /* accurate to the requested */
                /* tolerance)                */
                /* ------------------------- */

                av_(&nx, &v[(j << 8) - 256], ax);
                r__1 = -d[j - 1];
                saxpy_(&n, &r__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
                d[j + 24] = snrm2_(&n, ax, &c__1);
                d[j + 24] /= (r__1 = d[j - 1], dabs(r__1));

                /* L20: */
            }

            /* ----------------------------- */
            /* Display computed residuals    */
            /* ----------------------------- */

            smout_(nconv, 2, d, 25, -6, "Ritz values and relative residuals");
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
        printf(" _SDRV1 \n");
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
    /* Done with program ssdrv1. */
    /* ------------------------- */

    free(d);
    free(v);
    free(ax);
    free(resid);
    free(workd);
    free(workl);

    return 0;
}

/* ------------------------------------------------------------------ */
/*     matrix vector subroutine */

/*     The matrix used is the 2 dimensional discrete Laplacian on unit */
/*     square with zero Dirichlet boundary condition. */

/*     Computes w <--- OP*v, where OP is the nx*nx by nx*nx block */
/*     tridiagonal matrix */

/*                  | T -I          | */
/*                  |-I  T -I       | */
/*             OP = |   -I  T       | */
/*                  |        ...  -I| */
/*                  |           -I T| */

/*     The subroutine TV is called to computed y<---T*x. */

int av_(a_int *nx, float *v, float *w)
{
    /* System generated locals */
    a_int i__1;
    float r__1;

    /* Local variables */
    a_int j;
    float h2;
    a_int n2, lo;

    /* Parameter adjustments */
    --w;
    --v;

    tv_(nx, &v[1], &w[1]);
    saxpy_(nx, &c_b138, &v[*nx + 1], &c__1, &w[1], &c__1);

    i__1 = *nx - 1;
    for (j = 2; j <= i__1; ++j)
    {
        lo = (j - 1) * *nx;
        tv_(nx, &v[lo + 1], &w[lo + 1]);
        saxpy_(nx, &c_b138, &v[lo - *nx + 1], &c__1, &w[lo + 1], &c__1);
        saxpy_(nx, &c_b138, &v[lo + *nx + 1], &c__1, &w[lo + 1], &c__1);
        /* L10: */
    }

    lo = (*nx - 1) * *nx;
    tv_(nx, &v[lo + 1], &w[lo + 1]);
    saxpy_(nx, &c_b138, &v[lo - *nx + 1], &c__1, &w[lo + 1], &c__1);

    /*     Scale the vector w by (1/h^2), where h is the mesh size */

    n2 = *nx * *nx;
    h2 = 1.f / (float)((*nx + 1) * (*nx + 1));
    r__1 = 1.f / h2;
    sscal_(&n2, &r__1, &w[1], &c__1);
    return 0;
} /* av_ */

/* ------------------------------------------------------------------- */
int tv_(a_int *nx, float *x, float *y)
{
    /* System generated locals */
    a_int i__1;

    /* Local variables */
    a_int j;
    float dd, dl, du;

    /*     Compute the matrix vector multiplication y<---T*x */
    /*     where T is a nx by nx tridiagonal matrix with DD on the */
    /*     diagonal, DL on the subdiagonal, and DU on the superdiagonal. */

    /* Parameter adjustments */
    --y;
    --x;

    dd = 4.f;
    dl = -1.f;
    du = -1.f;

    y[1] = dd * x[1] + du * x[2];
    i__1 = *nx - 1;
    for (j = 2; j <= i__1; ++j)
    {
        y[j] = dl * x[j - 1] + dd * x[j] + du * x[j + 1];
        /* L10: */
    }
    y[*nx] = dl * x[*nx - 1] + dd * x[*nx];
    return 0;
} /* tv_ */
