/* EXAMPLES\SIMPLE\dssimp.f -- translated by f2c (version 20230428). */

#include <stdlib.h>
#include "arpack_internal.h"

static a_int i_one = 1;
static a_int c__256 = 256;

static double minus_one = -1.;

void av_(const a_int nx, double *v, double *w);
void tv_(const a_int nx, double *x, double *y);

/**
 * \BeginDoc
 *
 *     This example program is intended to illustrate the
 *     simplest case of using ARPACK in considerable detail.
 *     This code may be used to understand basic usage of ARPACK
 *     and as a template for creating an interface to ARPACK.
 *
 *     This code shows how to use ARPACK to find a few eigenvalues
 *     (lambda) and corresponding eigenvectors (x) for the standard
 *     eigenvalue problem:
 *
 *                        A*x = lambda*x
 *
 *     where A is an n by n real symmetric matrix.
 *
 *     The main points illustrated here are
 *
 *        1) How to declare sufficient memory to find NEV
 *           eigenvalues of largest magnitude.  Other options
 *           are available.
 *
 *        2) Illustration of the reverse communication interface
 *           needed to utilize the top level ARPACK routine DSAUPD
 *           that computes the quantities needed to construct
 *           the desired eigenvalues and eigenvectors(if requested).
 *
 *        3) How to extract the desired eigenvalues and eigenvectors
 *           using the ARPACK routine DSEUPD.
 *
 *     The only thing that must be supplied in order to use this
 *     routine on your problem is to change the array dimensions
 *     appropriately, to specify WHICH eigenvalues you want to compute
 *     and to supply a matrix-vector product
 *
 *                         w <-  Av
 *
 *     in place of the call to AV( ) below.
 *
 *     Once usage of this routine is understood, you may wish to explore
 *     the other available options to improve convergence, to solve generalized
 *     problems, etc.  Look at the file ex-sym.doc in DOCUMENTS directory.
 *     This codes implements
 *
 * \Example-1
 *     ... Suppose we want to solve A*x = lambda*x in regular mode,
 *         where A is derived from the central difference discretization
 *         of the 2-dimensional Laplacian on the unit square with
 *         zero Dirichlet boundary condition.
 *     ... OP = A  and  B = I.
 *     ... Assume "call av (n,x,y)" computes y = A*x
 *     ... Use mode 1 of DSAUPD.
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * Routines called:
 *     dsaupd  ARPACK reverse communication interface routine.
 *     dseupd  ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     dnrm2   Level 1 BLAS that computes the norm of a vector.
 *     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
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
    a_int j, n, nx, ido, ncv, nev, ierr = 0;
    a_int info, mode1, nconv, ishfts, lworkl, maxitr;
    char *bmat, *which;
    double tol, sigma;

    /* ---------------------------------------------------- */
    /* Storage Declarations:                                */
    /*                                                      */
    /* The maximum dimensions for all arrays are            */
    /* set here to accommodate a problem size of            */
    /* N .le. MAXN                                          */
    /*                                                      */
    /* NEV is the number of eigenvalues requested.          */
    /*     See specifications for ARPACK usage below.       */
    /*                                                      */
    /* NCV is the largest number of basis vectors that will */
    /*     be used in the Implicitly Restarted Arnoldi      */
    /*     Process.  Work per major iteration is            */
    /*     proportional to N*NCV*NCV.                       */
    /*                                                      */
    /* You must set:                                        */
    /*                                                      */
    /* MAXN:   Maximum dimension of the A allowed.          */
    /* MAXNEV: Maximum NEV allowed.                         */
    /* MAXNCV: Maximum NCV allowed.                         */
    /* ---------------------------------------------------- */

    /* ----------------------------------------------- */
    /* The following include statement and assignments */
    /* initiate trace output from the internal         */
    /* actions of ARPACK.  See debug.doc in the        */
    /* DOCUMENTS directory for usage.  Initially, the  */
    /* most useful information will be a breakdown of  */
    /* time spent in the various stages of computation */
    /* given by setting msaupd = 1.                    */
    /* ----------------------------------------------- */

    debug_1.ndigit = -3;
    debug_1.logfil = 6;
    debug_1.msgets = 0;
    debug_1.msaitr = 0;
    debug_1.msapps = 0;
    debug_1.msaupd = 1;
    debug_1.msaup2 = 0;
    debug_1.mseigt = 0;
    debug_1.mseupd = 0;

    /* ----------------------------------------------- */
    /* The following sets dimensions for this problem. */
    /* ----------------------------------------------- */

    nx = 10;
    n = nx * nx;

    /* --------------------------------------------- */
    /*                                               */
    /* Specifications for ARPACK usage are set       */
    /* below:                                        */
    /*                                               */
    /*    1) NEV = 4  asks for 4 eigenvalues to be   */
    /*       computed.                               */
    /*                                               */
    /*    2) NCV = 20 sets the length of the Arnoldi */
    /*       factorization                           */
    /*                                               */
    /*    3) This is a standard problem              */
    /*         (indicated by bmat  = 'I')            */
    /*                                               */
    /*    4) Ask for the NEV eigenvalues of          */
    /*       largest magnitude                       */
    /*         (indicated by which = 'LM')           */
    /*       See documentation in DSAUPD for the     */
    /*       other options SM, LA, SA, LI, SI.       */
    /*                                               */
    /* Note: NEV and NCV must satisfy the following  */
    /* conditions:                                   */
    /*              NEV <= MAXNEV                    */
    /*          NEV + 1 <= NCV <= MAXNCV             */
    /* --------------------------------------------- */

    nev = 4;
    ncv = 20;
    bmat = "I";
    which = "LM";

    if (n > 256)
    {
        printf(" ERROR with _SSIMP: N is greater than MAXN \n");
        return ierr;
    }
    else if (nev > 10)
    {
        printf(" ERROR with _SSIMP: NEV is greater than MAXNEV \n");
        return ierr;
    }
    else if (ncv > 25)
    {
        printf(" ERROR with _SSIMP: NCV is greater than MAXNCV \n");
        return ierr;
    }

    /* --------------------------------------------------- */
    /*                                                     */
    /* Specification of stopping rules and initial         */
    /* conditions before calling DSAUPD                    */
    /*                                                     */
    /* TOL  determines the stopping criterion.             */
    /*                                                     */
    /*      Expect                                         */
    /*           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) */
    /*               computed   true                       */
    /*                                                     */
    /*      If TOL .le. 0,  then TOL <- macheps            */
    /*           (machine precision) is used.              */
    /*                                                     */
    /* IDO  is the REVERSE COMMUNICATION parameter         */
    /*      used to specify actions to be taken on return  */
    /*      from DSAUPD. (See usage below.)                */
    /*                                                     */
    /*      It MUST initially be set to 0 before the first */
    /*      call to DSAUPD.                                */
    /*                                                     */
    /* INFO on entry specifies starting vector information */
    /*      and on return indicates error codes            */
    /*                                                     */
    /*      Initially, setting INFO=0 indicates that a     */
    /*      random starting vector is requested to         */
    /*      start the ARNOLDI iteration.  Setting INFO to  */
    /*      a nonzero value on the initial call is used    */
    /*      if you want to specify your own starting       */
    /*      vector (This vector must be placed in RESID.)  */
    /*                                                     */
    /* The work array WORKL is used in DSAUPD as           */
    /* workspace.  Its dimension LWORKL is set as          */
    /* illustrated below.                                  */
    /*                                                     */
    /* --------------------------------------------------- */

    lworkl = ncv * (ncv + 8);
    tol = 0.;
    info = 0;
    nconv = 0;
    ido = 0;

    double* d = (double*)malloc(sizeof(double) * 25 * 2);
    double* v = (double*)malloc(sizeof(double) * 256 * 25);
    double* ax = (double*)malloc(sizeof(double) * 256);
    double* resid = (double*)malloc(sizeof(double) * 256);
    double* workd = (double*)malloc(sizeof(double) * 768);
    double* workl = (double*)malloc(sizeof(double) * 825);

    /* ------------------------------------------------- */
    /* Specification of Algorithm Mode:                  */
    /*                                                   */
    /* This program uses the exact shift strategy        */
    /* (indicated by setting PARAM(1) = 1).              */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 1 of DSAUPD is used     */
    /* (IPARAM(7) = 1). All these options can be changed */
    /* by the user. For details see the documentation in */
    /* DSAUPD.                                           */
    /* ------------------------------------------------- */

    ishfts = 1;
    maxitr = 300;
    mode1 = 1;

    iparam[0] = ishfts;

    iparam[2] = maxitr;

    iparam[6] = mode1;

    /* ---------------------------------------------- */
    /* M A I N   L O O P (Reverse communication loop) */
    /* ---------------------------------------------- */

L10:

    /* ------------------------------------------- */
    /* Repeatedly call the routine DSAUPD and take */
    /* actions indicated by parameter IDO until    */
    /* either convergence is indicated or maxitr   */
    /* has been exceeded.                          */
    /* ------------------------------------------- */

    dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, &info);

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

        av_(nx, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

        /* --------------------------------------- */
        /* L O O P   B A C K to call DSAUPD again. */
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
        /* documentation in DSAUPD. */
        /* ------------------------ */

        printf(" \n");
        printf(" Error with _saupd info = %d\n", info);
        printf(" Check documentation in _saupd \n");
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
        /* Eigenvectors may be also computed now if  */
        /* desired.  (indicated by rvec = .true.)    */
        /*                                           */
        /* The routine DSEUPD now called to do this  */
        /* post processing (Other modes may require  */
        /* more complicated post processing than     */
        /* mode1.)                                   */
        /*                                           */
        /* ----------------------------------------- */

        rvec = TRUE_;

        dseupd_(&rvec, "A", select, d, v, &c__256, &sigma, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, &ierr);

        /* -------------------------------------------- */
        /* Eigenvalues are returned in the first column */
        /* of the two dimensional array D and the       */
        /* corresponding eigenvectors are returned in   */
        /* the first NCONV (=IPARAM(5)) columns of the  */
        /* two dimensional array V if requested.        */
        /* Otherwise, an orthogonal basis for the       */
        /* invariant subspace corresponding to the      */
        /* eigenvalues in D is returned in V.           */
        /* -------------------------------------------- */

        if (ierr != 0)
        {

            /* ---------------------------------- */
            /* Error condition:                   */
            /* Check the documentation of DSEUPD. */
            /* ---------------------------------- */

            printf(" \n");
            printf(" Error with _seupd info = %d\n", ierr);
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

                av_(nx, &v[(j << 8) - 256], ax);
                d__1 = -d[j - 1];
                daxpy_(&n, &d__1, &v[(j << 8) - 256], &i_one, ax, &i_one);
                d[j + 24] = dnrm2_(&n, ax, &i_one);
                d[j + 24] /= (d__1 = d[j - 1], abs(d__1));
            }

            /* --------------------------- */
            /* Display computed residuals. */
            /* --------------------------- */

            dmout_(nconv, 2, d, 25, -6, "Ritz values and relative residuals");
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
        printf(" _SSIMP \n");
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
    /* Done with program dssimp. */
    /* ------------------------- */

    free(d);
    free(v);
    free(ax);
    free(resid);
    free(workd);
    free(workl);

    return nconv < nev ? EXIT_FAILURE : EXIT_SUCCESS;
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

void av_(const a_int nx, double *v, double *w)
{
    /* System generated locals */
    a_int i__1;
    double d__1;

    /* Local variables */
    a_int j;
    double h2;
    a_int n2, lo;

    /* Parameter adjustments */
    --w;
    --v;

    tv_(nx, &v[1], &w[1]);
    daxpy_(&nx, &minus_one, &v[nx + 1], &i_one, &w[1], &i_one);

    i__1 = nx - 1;
    for (j = 2; j <= i__1; ++j)
    {
        lo = (j - 1) * nx;
        tv_(nx, &v[lo + 1], &w[lo + 1]);
        daxpy_(&nx, &minus_one, &v[lo - nx + 1], &i_one, &w[lo + 1], &i_one);
        daxpy_(&nx, &minus_one, &v[lo + nx + 1], &i_one, &w[lo + 1], &i_one);
    }

    lo = (nx - 1) * nx;
    tv_(nx, &v[lo + 1], &w[lo + 1]);
    daxpy_(&nx, &minus_one, &v[lo - nx + 1], &i_one, &w[lo + 1], &i_one);

    /*     Scale the vector w by (1/h^2), where h is the mesh size */

    n2 = nx * nx;
    h2 = 1. / (double)((nx + 1) * (nx + 1));
    d__1 = 1. / h2;
    dscal_(&n2, &d__1, &w[1], &i_one);
} /* av_ */

/* ------------------------------------------------------------------- */
void tv_(const a_int nx, double *x, double *y)
{
    /* System generated locals */
    a_int i__1;

    /* Local variables */
    a_int j;
    double dd, dl, du;

    /*     Compute the matrix vector multiplication y<---T*x */
    /*     where T is a nx by nx tridiagonal matrix with DD on the */
    /*     diagonal, DL on the subdiagonal, and DU on the superdiagonal. */

    /* Parameter adjustments */
    --y;
    --x;

    dd = 4.;
    dl = -1.;
    du = -1.;

    y[1] = dd * x[1] + du * x[2];
    i__1 = nx - 1;
    for (j = 2; j <= i__1; ++j)
    {
        y[j] = dl * x[j - 1] + dd * x[j] + du * x[j + 1];
    }
    y[nx] = dl * x[nx - 1] + dd * x[nx];
} /* tv_ */
