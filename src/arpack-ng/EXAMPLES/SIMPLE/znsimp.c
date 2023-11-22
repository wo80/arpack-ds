/* EXAMPLES\SIMPLE\znsimp.f -- translated by f2c (version 20230428). */

#include <stdlib.h>
#include "arpack_internal.h"

/* Table of constant values */

static a_int c__9 = 9;
static a_int c__1 = 1;
static a_int c__256 = 256;
static a_int c__3 = 3;
static a_int c__6 = 6;
static a_int c__30 = 30;
static a_int c_n6 = -6;
static a_int c__5 = 5;
static a_dcomplex c_b137 = {1., 0.};
static a_dcomplex c_b151 = {4., 0.};

int main()
{
    /* System generated locals */
    a_int i__1, i__2;
    a_dcomplex z__1;

    /* Local variables */
    a_bool select[30];
    a_int iparam[11];
    a_int ipntr[14];
    a_bool rvec;
    a_dcomplex sigma;
    a_int j, n, nx, ido, ncv, nev, ierr = 0;
    a_int info, mode1, nconv, ishfts, lworkl, maxitr;
    char *bmat, *which;
    double tol;

    /*     This example program is intended to illustrate the */
    /*     simplest case of using ARPACK in considerable detail. */
    /*     This code may be used to understand basic usage of ARPACK */
    /*     and as a template for creating an interface to ARPACK. */

    /*     This code shows how to use ARPACK to find a few eigenvalues */
    /*     (lambda) and corresponding eigenvectors (x) for the standard */
    /*     eigenvalue problem: */

    /*                        A*x = lambda*x */

    /*     where A is a general n by n complex matrix. */

    /*     The main points illustrated here are */

    /*        1) How to declare sufficient memory to find NEV */
    /*           eigenvalues of largest magnitude.  Other options */
    /*           are available. */

    /*        2) Illustration of the reverse communication interface */
    /*           needed to utilize the top level ARPACK routine ZNAUPD */
    /*           that computes the quantities needed to construct */
    /*           the desired eigenvalues and eigenvectors(if requested). */

    /*        3) How to extract the desired eigenvalues and eigenvectors */
    /*           using the ARPACK routine ZNEUPD . */

    /*     The only thing that must be supplied in order to use this */
    /*     routine on your problem is to change the array dimensions */
    /*     appropriately, to specify WHICH eigenvalues you want to compute */
    /*     and to supply a matrix-vector product */

    /*                         w <-  Av */

    /*     in place of the call to AV( )  below. */

    /*     Once usage of this routine is understood, you may wish to explore */
    /*     the other available options to improve convergence, to solve generalized */
    /*     problems, etc.  Look at the file ex-complex.doc in DOCUMENTS directory. */
    /*     This codes implements */

    /* \Example-1 */
    /*     ... Suppose we want to solve A*x = lambda*x in regular mode, */
    /*     ... OP = A  and  B = I. */
    /*     ... Assume "call av (nx,x,y)" computes y = A*x */
    /*     ... Use mode 1 of ZNAUPD . */

    /* \BeginLib */

    /* \Routines called */
    /*     znaupd   ARPACK reverse communication interface routine. */
    /*     zneupd   ARPACK routine that returns Ritz values and (optionally) */
    /*             Ritz vectors. */
    /*     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully. */
    /*     dznrm2   Level 1 BLAS that computes the norm of a complex vector. */
    /*     zaxpy    Level 1 BLAS that computes y <- alpha*x+y. */
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
    /* FILE: nsimp.F   SID: 2.4   DATE OF SID: 10/20/00   RELEASE: 2 */

    /* \Remarks */
    /*     1. None */

    /* \EndLib */
    /* --------------------------------------------------------------------------- */

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
    /* given by setting mcaupd = 1                     */
    /* ----------------------------------------------- */

    debug_1.ndigit = -3;
    debug_1.logfil = 6;
    debug_1.mcaitr = 0;
    debug_1.mcapps = 0;
    debug_1.mcaupd = 1;
    debug_1.mcaup2 = 0;
    debug_1.mceigh = 0;
    debug_1.mceupd = 0;

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
    /*       See documentation in ZNAUPD  for the     */
    /*       other options SM, LR, SR, LI, SI.       */
    /*                                               */
    /* Note: NEV and NCV must satisfy the following  */
    /* conditions:                                   */
    /*              NEV <= MAXNEV                    */
    /*          NEV + 2 <= NCV <= MAXNCV             */
    /*                                               */
    /* --------------------------------------------- */

    nev = 4;
    ncv = 20;
    bmat = "I";
    which = "LM";

    if (n > 256)
    {
        printf(" ERROR with _NSIMP: N is greater than MAXN \n");
        return ierr;
    }
    else if (nev > 12)
    {
        printf(" ERROR with _NSIMP: NEV is greater than MAXNEV \n");
        return ierr;
    }
    else if (ncv > 30)
    {
        printf(" ERROR with _NSIMP: NCV is greater than MAXNCV \n");
        return ierr;
    }

    /* --------------------------------------------------- */
    /*                                                     */
    /* Specification of stopping rules and initial         */
    /* conditions before calling ZNAUPD                     */
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
    /*      from ZNAUPD . (see usage below)                 */
    /*                                                     */
    /*      It MUST initially be set to 0 before the first */
    /*      call to ZNAUPD .                                */
    /*                                                     */
    /* INFO on entry specifies starting vector information */
    /*      and on return indicates error codes            */
    /*                                                     */
    /*      Initially, setting INFO=0 indicates that a     */
    /*      random starting vector is requested to         */
    /*      start the ARNOLDI iteration.  Setting INFO to  */
    /*      a nonzero value on the initial call is used    */
    /*      if you want to specify your own starting       */
    /*      vector (This vector must be placed in RESID).  */
    /*                                                     */
    /* The work array WORKL is used in ZNAUPD  as           */
    /* workspace.  Its dimension LWORKL is set as          */
    /* illustrated below.                                  */
    /*                                                     */
    /* --------------------------------------------------- */

    lworkl = ncv * ncv * 3 + ncv * 5;
    tol = 0.f;
    ido = 0;
    info = 0;
    nconv = 0;

    a_dcomplex* d = (a_dcomplex*)malloc(sizeof(a_dcomplex) * 30);
    a_dcomplex* v = (a_dcomplex*)malloc(sizeof(a_dcomplex) * 256 * 30);
    a_dcomplex* ax = (a_dcomplex*)malloc(sizeof(a_dcomplex) * 256);
    a_dcomplex* resid = (a_dcomplex*)malloc(sizeof(a_dcomplex) * 256);
    a_dcomplex* workd = (a_dcomplex*)malloc(sizeof(a_dcomplex) * 768);
    a_dcomplex* workl = (a_dcomplex*)malloc(sizeof(a_dcomplex) * 2850);
    a_dcomplex* workev = (a_dcomplex*)malloc(sizeof(a_dcomplex) * 60);
    double* rd = (double*)malloc(sizeof(double) * 30 * 3);
    double* rwork = (double*)malloc(sizeof(double) * 30);

    /* ------------------------------------------------- */
    /* Specification of Algorithm Mode:                  */
    /*                                                   */
    /* This program uses the exact shift strategy        */
    /* (indicated by setting IPARAM(1) = 1).             */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 1 of ZNAUPD  is used     */
    /* (IPARAM(7) = 1). All these options can be changed */
    /* by the user. For details see the documentation in */
    /* ZNAUPD .                                           */
    /* ------------------------------------------------- */

    ishfts = 1;
    maxitr = 300;
    mode1 = 1;

    iparam[0] = ishfts;

    iparam[2] = maxitr;

    iparam[6] = mode1;

    /* ---------------------------------------------- */
    /* M A I N   L O O P (Reverse Communication Loop) */
    /* ---------------------------------------------- */

L10:

    /* ------------------------------------------- */
    /* Repeatedly call the routine ZNAUPD  and take */
    /* actions indicated by parameter IDO until    */
    /* either convergence is indicated or maxitr   */
    /* has been exceeded.                          */
    /* ------------------------------------------- */
    znaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, rwork, &info);

    if (ido == -1 || ido == 1)
    {

        /* ----------------------------------------- */
        /* Perform matrix vector multiplication      */
        /*                                           */
        /*                y <--- A*x                 */
        /*                                           */
        /* The user should supply his/her own        */
        /* matrix vector multiplication routine here */
        /* that takes workd(ipntr(1)) as the input   */
        /* vector x , and returns the resulting      */
        /* matrix-vector product y = A*x in the      */
        /* array workd(ipntr(2)).                    */
        /* ----------------------------------------- */

        av_(&nx, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

        /* --------------------------------------- */
        /* L O O P   B A C K to call ZNAUPD  again. */
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
        /* Error message, check the */
        /* documentation in ZNAUPD   */
        /* ------------------------ */

        printf(" \n");
        printf(" Error with _naupd info = %d\n", info);
        printf(" Check the documentation of _naupd\n");
        printf(" \n");
    }
    else
    {

        /* ----------------------------------------- */
        /* No fatal errors occurred.                 */
        /* Post-Process using ZNEUPD .                */
        /*                                           */
        /* Computed eigenvalues may be extracted.    */
        /*                                           */
        /* Eigenvectors may be also computed now if  */
        /* desired.  (indicated by rvec = .true.)    */
        /*                                           */
        /* The routine ZNEUPD  now called to do this  */
        /* post processing (Other modes may require  */
        /* more complicated post processing than     */
        /* mode1.)                                   */
        /*                                           */
        /* ----------------------------------------- */

        rvec = TRUE_;

        zneupd_(&rvec, "A", select, d, v, &c__256, &sigma, workev, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, rwork, &ierr);

        /* --------------------------------------------- */
        /* Eigenvalues are returned in the one           */
        /* dimensional array D and the corresponding     */
        /* eigenvectors are returned in the first        */
        /* NCONV (=IPARAM(5)) columns of the two         */
        /* dimensional array V if requested.  Otherwise, */
        /* an orthogonal basis for the invariant         */
        /* subspace corresponding to the eigenvalues in  */
        /* D is returned in V.                           */
        /* --------------------------------------------- */

        if (ierr != 0)
        {

            /* ---------------------------------- */
            /* Error condition:                   */
            /* Check the documentation of ZNEUPD . */
            /* ---------------------------------- */

            printf(" \n");
            printf(" Error with _neupd info = %d\n", ierr);
            printf(" Check the documentation of _neupd. \n");
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
                i__2 = j - 1;
                z__1.r = -d[i__2].r, z__1.i = -d[i__2].i;
                zaxpy_(&n, &z__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
                i__2 = j - 1;
                rd[j - 1] = d[i__2].r;
                rd[j + 29] = d[j - 1].i;
                rd[j + 59] = dznrm2_(&n, ax, &c__1);
                rd[j + 59] /= dlapy2_(&rd[j - 1], &rd[j + 29]);
            }

            /* --------------------------- */
            /* Display computed residuals. */
            /* --------------------------- */

            dmout_(nconv, 3, rd, 30, -6, "Ritz values (Real, Imag) and relative residuals");
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
        printf("_NSIMP \n");
        printf("====== \n");
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
    /* Done with program znsimp . */
    /* ------------------------- */

    free(d);
    free(v);
    free(ax);
    free(resid);
    free(workd);
    free(workl);
    free(workev);
    free(rd);
    free(rwork);

    return nconv < nev ? EXIT_FAILURE : EXIT_SUCCESS;
}

/* ========================================================================== */

/*     matrix vector subroutine */

/*     The matrix used is the convection-diffusion operator */
/*     discretized using centered difference. */

int av_(a_int *nx, a_dcomplex *v, a_dcomplex *w)
{
    /* System generated locals */
    a_int i__1;
    a_dcomplex z__1, z__2;

    /* Local variables */
    a_int j;
    a_dcomplex h2;
    a_int lo;

    /*     Computes w <--- OP*v, where OP is the nx*nx by nx*nx block */
    /*     tridiagonal matrix */

    /*                  | T -I          | */
    /*                  |-I  T -I       | */
    /*             OP = |   -I  T       | */
    /*                  |        ...  -I| */
    /*                  |           -I T| */

    /*     derived from the standard central difference discretization */
    /*     of the 2-dimensional convection-diffusion operator */
    /*                  (Laplacian u) + rho*(du/dx) */
    /*     on the unit squqre with zero boundary condition. */

    /*     The subroutine TV is called to computed y<---T*x. */

    /* Parameter adjustments */
    --w;
    --v;

    i__1 = (*nx + 1) * (*nx + 1);
    z__2.r = (double)i__1, z__2.i = 0.;
    ar_z_div(&z__1, &c_b137, &z__2);
    h2.r = z__1.r, h2.i = z__1.i;

    tv_(nx, &v[1], &w[1]);
    z__2.r = -1., z__2.i = -0.;
    ar_z_div(&z__1, &z__2, &h2);
    zaxpy_(nx, &z__1, &v[*nx + 1], &c__1, &w[1], &c__1);

    i__1 = *nx - 1;
    for (j = 2; j <= i__1; ++j)
    {
        lo = (j - 1) * *nx;
        tv_(nx, &v[lo + 1], &w[lo + 1]);
        z__2.r = -1., z__2.i = -0.;
        ar_z_div(&z__1, &z__2, &h2);
        zaxpy_(nx, &z__1, &v[lo - *nx + 1], &c__1, &w[lo + 1], &c__1);
        z__2.r = -1., z__2.i = -0.;
        ar_z_div(&z__1, &z__2, &h2);
        zaxpy_(nx, &z__1, &v[lo + *nx + 1], &c__1, &w[lo + 1], &c__1);
    }

    lo = (*nx - 1) * *nx;
    tv_(nx, &v[lo + 1], &w[lo + 1]);
    z__2.r = -1., z__2.i = -0.;
    ar_z_div(&z__1, &z__2, &h2);
    zaxpy_(nx, &z__1, &v[lo - *nx + 1], &c__1, &w[lo + 1], &c__1);

    return 0;
} /* av_ */

/* ========================================================================= */
int tv_(a_int *nx, a_dcomplex *x, a_dcomplex *y)
{
    /* System generated locals */
    a_int i__1, i__2, i__3, i__4, i__5;
    a_dcomplex z__1, z__2, z__3, z__4, z__5;

    /* Local variables */
    a_dcomplex h;
    a_int j;
    a_dcomplex h2, dd, dl, du;

    /*     Compute the matrix vector multiplication y<---T*x */
    /*     where T is a nx by nx tridiagonal matrix with DD on the */
    /*     diagonal, DL on the subdiagonal, and DU on the superdiagonal */

    /* Parameter adjustments */
    --y;
    --x;

    i__1 = *nx + 1;
    z__2.r = (double)i__1, z__2.i = 0.;
    ar_z_div(&z__1, &c_b137, &z__2);
    h.r = z__1.r, h.i = z__1.i;
    z__1.r = h.r * h.r - h.i * h.i, z__1.i = h.r * h.i + h.i * h.r;
    h2.r = z__1.r, h2.i = z__1.i;
    ar_z_div(&z__1, &c_b151, &h2);
    dd.r = z__1.r, dd.i = z__1.i;
    z__3.r = -1., z__3.i = -0.;
    ar_z_div(&z__2, &z__3, &h2);
    z__5.r = 50., z__5.i = 0.;
    ar_z_div(&z__4, &z__5, &h);
    z__1.r = z__2.r - z__4.r, z__1.i = z__2.i - z__4.i;
    dl.r = z__1.r, dl.i = z__1.i;
    z__3.r = -1., z__3.i = -0.;
    ar_z_div(&z__2, &z__3, &h2);
    z__5.r = 50., z__5.i = 0.;
    ar_z_div(&z__4, &z__5, &h);
    z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
    du.r = z__1.r, du.i = z__1.i;

    z__2.r = dd.r * x[1].r - dd.i * x[1].i, z__2.i = dd.r * x[1].i + dd.i * x[1].r;
    z__3.r = du.r * x[2].r - du.i * x[2].i, z__3.i = du.r * x[2].i + du.i * x[2].r;
    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
    y[1].r = z__1.r, y[1].i = z__1.i;
    i__1 = *nx - 1;
    for (j = 2; j <= i__1; ++j)
    {
        i__2 = j;
        i__3 = j - 1;
        z__3.r = dl.r * x[i__3].r - dl.i * x[i__3].i, z__3.i = dl.r * x[i__3].i + dl.i * x[i__3].r;
        i__4 = j;
        z__4.r = dd.r * x[i__4].r - dd.i * x[i__4].i, z__4.i = dd.r * x[i__4].i + dd.i * x[i__4].r;
        z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
        i__5 = j + 1;
        z__5.r = du.r * x[i__5].r - du.i * x[i__5].i, z__5.i = du.r * x[i__5].i + du.i * x[i__5].r;
        z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
        y[i__2].r = z__1.r, y[i__2].i = z__1.i;
    }
    i__1 = *nx;
    i__2 = *nx - 1;
    z__2.r = dl.r * x[i__2].r - dl.i * x[i__2].i, z__2.i = dl.r * x[i__2].i + dl.i * x[i__2].r;
    i__3 = *nx;
    z__3.r = dd.r * x[i__3].r - dd.i * x[i__3].i, z__3.i = dd.r * x[i__3].i + dd.i * x[i__3].r;
    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
    y[i__1].r = z__1.r, y[i__1].i = z__1.i;
    return 0;
} /* tv_ */
