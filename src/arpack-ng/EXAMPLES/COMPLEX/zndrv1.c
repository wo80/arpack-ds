/* EXAMPLES\COMPLEX\zndrv1.f -- translated by f2c (version 20230428). */

#include <stdlib.h>
#include "arpack_internal.h"

static a_int i_one = 1;
static a_int c__256 = 256;

static a_dcomplex one = {1., 0.};
static a_dcomplex four = {4., 0.};

void av_(const a_int nx, a_dcomplex *v, a_dcomplex *w);
void tv_(const a_int nx, a_dcomplex *x, a_dcomplex *y);

/**
 * \BeginDoc
 *
 *     Example program to illustrate the idea of reverse communication
 *     for a standard complex nonsymmetric eigenvalue problem.
 *
 *     We implement example one of ex-complex.doc in DOCUMENTS directory
 *
 * \Example-1
 *     ... Suppose we want to solve A*x = lambda*x in regular mode,
 *         where A is obtained from the standard central difference
 *         discretization of the convection-diffusion operator
 *                 (Laplacian u) + rho*(du / dx)
 *         on the unit squre [0,1]x[0,1] with zero Dirichlet boundary
 *         condition.
 *
 *     ... OP = A  and  B = I.
 *
 *     ... Assume "call av (nx,x,y)" computes y = A*x
 *
 *     ... Use mode 1 of ZNAUPD .
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * Routines called
 *     znaupd   ARPACK reverse communication interface routine.
 *     zneupd   ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     dznrm2   Level 1 BLAS that computes the norm of a complex vector.
 *     zaxpy    Level 1 BLAS that computes y <- alpha*x+y.
 *     av      Matrix vector multiplication routine that computes A*x.
 *     tv      Matrix vector multiplication routine that computes T*x,
 *             where T is a tridiagonal matrix.  It is used in routine
 *             av.
 *
 * \EndLib
 */
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
    a_int info, mode, nconv, ishfts, lworkl, maxitr;
    char *bmat, *which;
    double tol;

    /* Define maximum dimensions for all arrays. */

    const int MAXN   = 256; /* Maximum dimension of the A allowed. */
    const int MAXNEV =  12; /* Maximum NEV allowed */
    const int MAXNCV =  30; /* Maximum NCV allowed */

    /* ------------------------------------------------ */
    /* The number NX is the number of interior points   */
    /* in the discretization of the 2-dimensional       */
    /* convection-diffusion operator on the unit        */
    /* square with zero Dirichlet boundary condition.   */
    /* The number N(=NX*NX) is the dimension of the     */
    /* matrix.  A standard eigenvalue problem is        */
    /* solved (BMAT = 'I').  NEV is the number of       */
    /* eigenvalues to be approximated.  The user can    */
    /* modify NX, NEV, NCV, WHICH to solve problems of  */
    /* different sizes, and to get different parts of   */
    /* the spectrum.  However, The following            */
    /* conditions must be satisfied:                    */
    /*                   N <= MAXN                      */
    /*                 NEV <= MAXNEV                    */
    /*           NEV + 2 <= NCV <= MAXNCV               */
    /* ------------------------------------------------ */

    nx = 10;
    n = nx * nx;
    nev = 4;
    ncv = 20;
    if (n > 256)
    {
        printf(" ERROR with _NDRV1: N is greater than MAXN \n");
        return ierr;
    }
    else if (nev > 12)
    {
        printf(" ERROR with _NDRV1: NEV is greater than MAXNEV \n");
        return ierr;
    }
    else if (ncv > 30)
    {
        printf(" ERROR with _NDRV1: NCV is greater than MAXNCV \n");
        return ierr;
    }
    bmat = "I";
    which = "LM";

    /* ------------------------------------------------- */
    /* The work array WORKL is used in ZNAUPD  as         */
    /* workspace.  Its dimension LWORKL is set as        */
    /* illustrated below.  The parameter TOL determines  */
    /* the stopping criterion. If TOL<=0, machine        */
    /* precision is used.  The variable IDO is used for  */
    /* reverse communication, and is initially set to 0. */
    /* Setting INFO=0 indicates that a random vector is  */
    /* generated to start the ARNOLDI iteration.         */
    /* ------------------------------------------------- */

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
    a_dcomplex* workev = (a_dcomplex*)malloc(sizeof(a_dcomplex) * 90);
    double* rd = (double*)malloc(sizeof(double) * 30 * 3);
    double* rwork = (double*)malloc(sizeof(double) * 30);

    /* ------------------------------------------------- */
    /* This program uses exact shift with respect to     */
    /* the current Hessenberg matrix (IPARAM(1) = 1).    */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 1 of ZNAUPD  is used     */
    /* (IPARAM(7) = 1). All these options can be changed */
    /* by the user. For details see the documentation in */
    /* ZNAUPD .                                           */
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
        /*                y <--- OP*x                */
        /* The user should supply his/her own        */
        /* matrix vector multiplication routine here */
        /* that takes workd(ipntr(1)) as the input   */
        /* vector, and return the matrix vector      */
        /* product to workd(ipntr(2)).               */
        /* ----------------------------------------- */

        av_(nx, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

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
        /* Eigenvectors may also be computed now if  */
        /* desired.  (indicated by rvec = .true.)    */
        /* ----------------------------------------- */

        rvec = TRUE_;

        zneupd_(&rvec, "A", select, d, v, &c__256, &sigma, workev, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, rwork, &ierr);

        /* -------------------------------------------- */
        /* Eigenvalues are returned in the one          */
        /* dimensional array D.  The corresponding      */
        /* eigenvectors are returned in the first NCONV */
        /* (=IPARAM(5)) columns of the two dimensional  */
        /* array V if requested.  Otherwise, an         */
        /* orthogonal basis for the invariant subspace  */
        /* corresponding to the eigenvalues in D is     */
        /* returned in V.                               */
        /* -------------------------------------------- */

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

                av_(nx, &v[(j << 8) - 256], ax);
                i__2 = j - 1;
                z__1.r = -d[i__2].r, z__1.i = -d[i__2].i;
                zaxpy_(&n, &z__1, &v[(j << 8) - 256], &i_one, ax, &i_one);
                i__2 = j - 1;
                rd[j - 1] = d[i__2].r;
                rd[j + 29] = d[j - 1].i;
                rd[j + 59] = dznrm2_(&n, ax, &i_one);
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
        printf("_NDRV1\n");
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
    /* Done with program zndrv1 . */
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

/** Matrix vector subroutine.
 *
 * The matrix used is the convection-diffusion operator
 * discretized using centered difference.
 * 
 * Computes w <--- OP*v, where OP is the nx*nx by nx*nx block
 * tridiagonal matrix
 *
 *              | T -I          |
 *              |-I  T -I       |
 *         OP = |   -I  T       |
 *              |        ...  -I|
 *              |           -I T|
 *
 * derived from the standard central difference  discretization
 * of the 2-dimensional convection-diffusion operator
 *              (Laplacian u) + rho*(du/dx)
 * on the unit squqre with zero boundary condition.
 *
 * The subroutine TV is called to computed y<---T*x.
 */
void av_(const a_int nx, a_dcomplex *v, a_dcomplex *w)
{
    /* System generated locals */
    a_int i__1;
    a_dcomplex z__1, z__2;

    /* Local variables */
    a_int j;
    a_dcomplex h2;
    a_int lo;

    /* Parameter adjustments */
    --w;
    --v;

    i__1 = (nx + 1) * (nx + 1);
    z__2.r = (double)i__1, z__2.i = 0.;
    ar_z_div(&z__1, &one, &z__2);
    h2.r = z__1.r, h2.i = z__1.i;

    tv_(nx, &v[1], &w[1]);
    z__2.r = -1., z__2.i = -0.;
    ar_z_div(&z__1, &z__2, &h2);
    zaxpy_(&nx, &z__1, &v[nx + 1], &i_one, &w[1], &i_one);

    i__1 = nx - 1;
    for (j = 2; j <= i__1; ++j)
    {
        lo = (j - 1) * nx;
        tv_(nx, &v[lo + 1], &w[lo + 1]);
        z__2.r = -1., z__2.i = -0.;
        ar_z_div(&z__1, &z__2, &h2);
        zaxpy_(&nx, &z__1, &v[lo - nx + 1], &i_one, &w[lo + 1], &i_one);
        z__2.r = -1., z__2.i = -0.;
        ar_z_div(&z__1, &z__2, &h2);
        zaxpy_(&nx, &z__1, &v[lo + nx + 1], &i_one, &w[lo + 1], &i_one);
    }

    lo = (nx - 1) * nx;
    tv_(nx, &v[lo + 1], &w[lo + 1]);
    z__2.r = -1., z__2.i = -0.;
    ar_z_div(&z__1, &z__2, &h2);
    zaxpy_(&nx, &z__1, &v[lo - nx + 1], &i_one, &w[lo + 1], &i_one);
} /* av_ */

/**
 * Compute the matrix vector multiplication y<---T*x
 * where T is a nx by nx tridiagonal matrix with DD on the
 * diagonal, DL on the subdiagonal, and DU on the superdiagonal
 */
void tv_(const a_int nx, a_dcomplex *x, a_dcomplex *y)
{
    /* System generated locals */
    a_int i__1, i__2, i__3, i__4, i__5;
    a_dcomplex z__1, z__2, z__3, z__4, z__5;

    /* Local variables */
    a_dcomplex h;
    a_int j;
    a_dcomplex h2, dd, dl, du;

    /* Parameter adjustments */
    --y;
    --x;

    i__1 = nx + 1;
    z__2.r = (double)i__1, z__2.i = 0.;
    ar_z_div(&z__1, &one, &z__2);
    h.r = z__1.r, h.i = z__1.i;
    z__1.r = h.r * h.r - h.i * h.i, z__1.i = h.r * h.i + h.i * h.r;
    h2.r = z__1.r, h2.i = z__1.i;
    ar_z_div(&z__1, &four, &h2);
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
    i__1 = nx - 1;
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
    i__1 = nx;
    i__2 = nx - 1;
    z__2.r = dl.r * x[i__2].r - dl.i * x[i__2].i, z__2.i = dl.r * x[i__2].i + dl.i * x[i__2].r;
    i__3 = nx;
    z__3.r = dd.r * x[i__3].r - dd.i * x[i__3].i, z__3.i = dd.r * x[i__3].i + dd.i * x[i__3].r;
    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
    y[i__1].r = z__1.r, y[i__1].i = z__1.i;
} /* tv_ */
