/* EXAMPLES\NONSYM\dndrv3.f -- translated by f2c (version 20230428). */

#include <stdlib.h>
#include "arpack_internal.h"

static a_int i_one = 1;
static a_int c__256 = 256;

void av_(const a_int n, double *v, double *w);
void mv_(const a_int n, double *v, double *w);

/**
 * \BeginDoc
 *
 *     Simple program to illustrate the idea of reverse communication
 *     in inverse mode for a generalized nonsymmetric eigenvalue problem.
 *
 *     We implement example three of ex-nonsym.doc in DOCUMENTS directory
 *
 * \Example-3
 *     ... Suppose we want to solve A*x = lambda*B*x in inverse mode,
 *         where A and B are derived from the finite element discretization
 *         of the 1-dimensional convection-diffusion operator
 *                           (d^2u / dx^2) + rho*(du/dx)
 *         on the interval [0,1] with zero Dirichlet boundary condition
 *         using linear elements.
 *
 *     ... So OP = inv[M]*A  and  B = M.
 *
 *     ... Use mode 2 of DNAUPD.
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * Routines called:
 *     dnaupd  ARPACK reverse communication interface routine.
 *     dneupd  ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     dpttrf  LAPACK symmetric positive definite tridiagonal factorization
 *             routine.
 *     dpttrs  LAPACK symmetric positive definite tridiagonal solve routine.
 *     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
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
    a_int ipntr[14];
    a_bool rvec, first;
    a_int j, n, ido, ncv, nev, info, ierr = 0;
    a_int mode, nconv, ishfts, lworkl, maxitr;
    char *bmat, *which;
    double h, tol, sigmai, sigmar;

    /* Define maximum dimensions for all arrays. */

    const int MAXN   = 256; /* Maximum dimension of the A allowed. */
    const int MAXNEV =  10; /* Maximum NEV allowed */
    const int MAXNCV =  25; /* Maximum NCV allowed */

    /* -------------------------------------------------- */
    /* The number N is the dimension of the matrix.  A    */
    /* generalized eigenvalue problem is solved (BMAT =   */
    /* 'G').  NEV is the number of eigenvalues to be      */
    /* approximated.  The user can modify NEV, NCV, WHICH */
    /* to solve problems of different sizes, and to get   */
    /* different parts of the spectrum.  However, The     */
    /* following conditions must be satisfied:            */
    /*                    N <= MAXN,                      */
    /*                  NEV <= MAXNEV,                    */
    /*              NEV + 2 <= NCV <= MAXNCV              */
    /* -------------------------------------------------- */

    n = 100;
    nev = 4;
    ncv = 20;
    if (n > 256)
    {
        printf(" ERROR with _NDRV3: N is greater than MAXN \n");
        return ierr;
    }
    else if (nev > 10)
    {
        printf(" ERROR with _NDRV3: NEV is greater than MAXNEV \n");
        return ierr;
    }
    else if (ncv > 25)
    {
        printf(" ERROR with _NDRV3: NCV is greater than MAXNCV \n");
        return ierr;
    }
    bmat = "G";
    which = "LM";

    /* ---------------------------------------------- */
    /* M is the mass matrix formed by using piecewise */
    /* linear elements on [0,1].                      */
    /* ---------------------------------------------- */

    double* md = (double*)malloc(sizeof(double) * 256);
    double* me = (double*)malloc(sizeof(double) * 255);

    h = 1. / (double)(n + 1);
    i__1 = n - 1;
    for (j = 1; j <= i__1; ++j)
    {
        md[j - 1] = h * 4.;
        me[j - 1] = h * 1.;
    }
    md[n - 1] = h * 4.;

    dpttrf_(&n, md, me, &ierr);
    if (ierr != 0)
    {
        printf(" \n");
        printf(" ERROR with _pttrf. \n");
        printf(" \n");
        return ierr;
    }

    /* --------------------------------------------------- */
    /* The work array WORKL is used in DNAUPD as           */
    /* workspace.  Its dimension LWORKL is set as          */
    /* illustrated below.  The parameter TOL determines    */
    /* the stopping criterion. If TOL<=0, machine          */
    /* precision is used.  The variable IDO is used for    */
    /* reverse communication, and is initially set to 0.   */
    /* Setting INFO=0 indicates that a random vector is    */
    /* generated in DNAUPD to start the Arnoldi iteration. */
    /* --------------------------------------------------- */

    lworkl = ncv * ncv * 3 + ncv * 6;
    tol = 0.f;
    ido = 0;
    info = 0;
    nconv = 0;

    double* d = (double*)malloc(sizeof(double) * 25 * 3);
    double* v = (double*)malloc(sizeof(double) * 256 * 25);
    double* ax = (double*)malloc(sizeof(double) * 256);
    double* mx = (double*)malloc(sizeof(double) * 256);
    double* resid = (double*)malloc(sizeof(double) * 256);
    double* workd = (double*)malloc(sizeof(double) * 768);
    double* workl = (double*)malloc(sizeof(double) * 2025);
    double* workev = (double*)malloc(sizeof(double) * 75);

    /* ------------------------------------------------- */
    /* This program uses exact shifts with respect to    */
    /* the current Hessenberg matrix (IPARAM(1) = 1).    */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 2 of DNAUPD is used     */
    /* (IPARAM(7) = 2).  All these options can be        */
    /* changed by the user. For details, see the         */
    /* documentation in DNAUPD.                          */
    /* ------------------------------------------------- */

    ishfts = 1;
    maxitr = 300;
    mode = 2;

    iparam[0] = ishfts;
    iparam[2] = maxitr;
    iparam[6] = mode;

    /* ----------------------------------------- */
    /* M A I N   L O O P (Reverse communication) */
    /* ----------------------------------------- */

L10:

    /* ------------------------------------------- */
    /* Repeatedly call the routine DNAUPD and take */
    /* actions indicated by parameter IDO until    */
    /* either convergence is indicated or maxitr   */
    /* has been exceeded.                          */
    /* ------------------------------------------- */

    dnaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, &info);

    if (ido == -1 || ido == 1)
    {

        /* -------------------------------------- */
        /* Perform  y <--- OP*x = inv[M]*A*x      */
        /* The user should supply his/her own     */
        /* matrix vector routine and a linear     */
        /* system solver.  The matrix-vector      */
        /* subroutine should take workd(ipntr(1)) */
        /* as input, and the final result should  */
        /* be returned to workd(ipntr(2)).        */
        /* -------------------------------------- */

        av_(n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
        dpttrs_(&n, &i_one, md, me, &workd[ipntr[1] - 1], &n, &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" ERROR with _pttrs. \n");
            printf(" \n");
            return ierr;
        }

        /* --------------------------------------- */
        /* L O O P   B A C K to call DNAUPD again. */
        /* --------------------------------------- */

        goto L10;
    }
    else if (ido == 2)
    {

        /* ----------------------------------- */
        /*        Perform  y <--- M*x          */
        /* The matrix vector multiplication    */
        /* routine should take workd(ipntr(1)) */
        /* as input and return the result to   */
        /* workd(ipntr(2)).                    */
        /* ----------------------------------- */

        mv_(n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

        /* --------------------------------------- */
        /* L O O P   B A C K to call DNAUPD again. */
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
        /* Error message. Check the */
        /* documentation in DNAUPD. */
        /* ------------------------ */

        printf(" \n");
        printf(" Error with _naupd info = %d\n", info);
        printf(" Check the documentation of _naupd.\n");
        printf(" \n");
    }
    else
    {

        /* ----------------------------------------- */
        /* No fatal errors occurred.                 */
        /* Post-Process using DNEUPD.                */
        /*                                           */
        /* Computed eigenvalues may be extracted.    */
        /*                                           */
        /* Eigenvectors may also be computed now if  */
        /* desired.  (indicated by rvec = .true.)    */
        /* ----------------------------------------- */

        rvec = TRUE_;
        dneupd_(&rvec, "A", select, d, &d[25], v, &c__256, &sigmar, &sigmai, workev, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, &ierr);

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
            /* Check the documentation of DNEUPD. */
            /* ---------------------------------- */

            printf(" \n");
            printf(" Error with _neupd info = %d\n", ierr);
            printf(" Check the documentation of _neupd\n");
            printf(" \n");
        }
        else
        {

            first = TRUE_;
            nconv = iparam[4];
            i__1 = iparam[4];
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

                if (d[j + 24] == 0.)
                {

                    /* ------------------ */
                    /* Ritz value is real */
                    /* ------------------ */

                    av_(n, &v[(j << 8) - 256], ax);
                    mv_(n, &v[(j << 8) - 256], mx);
                    d__1 = -d[j - 1];
                    daxpy_(&n, &d__1, mx, &i_one, ax, &i_one);
                    d[j + 49] = dnrm2_(&n, ax, &i_one);
                    d[j + 49] /= (d__1 = d[j - 1], abs(d__1));
                }
                else if (first)
                {

                    /* ---------------------- */
                    /* Ritz value is complex  */
                    /* Residual of one Ritz   */
                    /* value of the conjugate */
                    /* pair is computed.      */
                    /* ---------------------- */

                    av_(n, &v[(j << 8) - 256], ax);
                    mv_(n, &v[(j << 8) - 256], mx);
                    d__1 = -d[j - 1];
                    daxpy_(&n, &d__1, mx, &i_one, ax, &i_one);
                    mv_(n, &v[(j + 1 << 8) - 256], mx);
                    daxpy_(&n, &d[j + 24], mx, &i_one, ax, &i_one);
                    /* Computing 2nd power */
                    d__1 = dnrm2_(&n, ax, &i_one);
                    d[j + 49] = d__1 * d__1;
                    av_(n, &v[(j + 1 << 8) - 256], ax);
                    mv_(n, &v[(j + 1 << 8) - 256], mx);
                    d__1 = -d[j - 1];
                    daxpy_(&n, &d__1, mx, &i_one, ax, &i_one);
                    mv_(n, &v[(j << 8) - 256], mx);
                    d__1 = -d[j + 24];
                    daxpy_(&n, &d__1, mx, &i_one, ax, &i_one);
                    d__1 = dnrm2_(&n, ax, &i_one);
                    d[j + 49] = dlapy2_(&d[j + 49], &d__1);
                    d[j + 49] /= dlapy2_(&d[j - 1], &d[j + 24]);
                    d[j + 50] = d[j + 49];
                    first = FALSE_;
                }
                else
                {
                    first = TRUE_;
                }
            }

            /* --------------------------- */
            /* Display computed residuals. */
            /* --------------------------- */

            dmout_(nconv, 3, d, 25, -6, "Ritz values (Real,Imag) and relative residuals");
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
        printf(" _NDRV3 \n");
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
    /* Done with program dndrv3. */
    /* ------------------------- */

    free(d);
    free(v);
    free(ax);
    free(md);
    free(me);
    free(mx);
    free(resid);
    free(workd);
    free(workl);
    free(workev);

    return nconv < nev ? EXIT_FAILURE : EXIT_SUCCESS;
}

/* ========================================================================== */

/*     matrix vector multiplication subroutine */

void av_(const a_int n, double *v, double *w)
{
    /* System generated locals */
    a_int i__1;

    /* Local variables */
    double h;
    a_int j;
    double s, dd, dl, du;

    /*     Compute the matrix vector multiplication y<---A*x */
    /*     where A is stiffness matrix obtained from the finite element */
    /*     discretization of the 1-dimensional convection diffusion operator */
    /*                           d^2u/dx^2 + rho*(du/dx) */
    /*     on the interval [0,1] with zero Dirichlet boundary condition using */
    /*     linear elements. */

    /* Parameter adjustments */
    --w;
    --v;

    h = 1. / (double)(n + 1);
    s = 5.;
    dd = 2. / h;
    dl = -1. / h - s;
    du = -1. / h + s;

    w[1] = dd * v[1] + du * v[2];
    i__1 = n - 1;
    for (j = 2; j <= i__1; ++j)
    {
        w[j] = dl * v[j - 1] + dd * v[j] + du * v[j + 1];
    }
    w[n] = dl * v[n - 1] + dd * v[n];
} /* av_ */

/* ------------------------------------------------------------------------ */
void mv_(const a_int n, double *v, double *w)
{
    /* System generated locals */
    a_int i__1;

    /* Local variables */
    double h;
    a_int j;

    /*     Compute the matrix vector multiplication y<---M*x */
    /*     where M is the mass matrix formed by using piecewise linear */
    /*     elements on [0,1]. */

    /* Parameter adjustments */
    --w;
    --v;

    w[1] = v[1] * 4. + v[2] * 1.;
    i__1 = n - 1;
    for (j = 2; j <= i__1; ++j)
    {
        w[j] = v[j - 1] * 1. + v[j] * 4. + v[j + 1] * 1.;
    }
    w[n] = v[n - 1] * 1. + v[n] * 4.;

    h = 1. / (double)(n + 1);
    dscal_(&n, &h, &w[1], &i_one);
} /* mv_ */
