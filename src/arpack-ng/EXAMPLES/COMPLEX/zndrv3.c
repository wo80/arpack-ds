/* EXAMPLES\COMPLEX\zndrv3.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

/* Table of constant values */

static a_dcomplex c_b2 = {1., 0.};
static a_int c__9 = 9;
static a_int c__1 = 1;
static a_int c__256 = 256;
static a_int c__3 = 3;
static a_int c__6 = 6;
static a_int c__25 = 25;
static a_int c_n6 = -6;
static a_int c__5 = 5;
static a_dcomplex c_b163 = {2., 0.};
static a_dcomplex c_b164 = {10., 0.};

/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    a_int i__1, i__2;
    a_dcomplex z__1, z__2;

    /* Local variables */
    a_dcomplex d[25], h;
    a_int j, n;
    a_dcomplex v[6400] /* was [256][25] */, dd[256], dl[256];
    double rd[75] /* was [25][3] */;
    a_dcomplex ax[256], du[256];
    a_dcomplex mx[256], du2[256];
    a_int ido, ncv, nev;
    double tol;
    char* bmat;
    a_int mode, info;
    a_bool rvec;
    a_int ierr, ipiv[256];
    a_dcomplex sigma;
    char* which;
    a_dcomplex resid[256];
    a_int nconv;
    a_dcomplex workd[768];
    a_int ipntr[14];
    a_dcomplex workl[2000];
    double rwork[256];
    a_int iparam[11];
    a_bool select[25];
    a_int ishfts;
    a_int maxitr;
    a_int lworkl;
    a_dcomplex workev[50];

    /*     Simple program to illustrate the idea of reverse communication */
    /*     in inverse mode for a generalized complex nonsymmetric eigenvalue */
    /*     problem. */

    /*     We implement example three of ex-complex.doc in DOCUMENTS directory */

    /* \Example-3 */
    /*     ... Suppose we want to solve A*x = lambda*B*x in regular mode, */
    /*         where A and B are derived from the finite element discretization */
    /*         of the 1-dimensional convection-diffusion operator */
    /*                   (d^2u/dx^2) + rho*(du/dx) */
    /*         on the interval [0,1] with zero boundary condition using */
    /*         piecewise linear elements. */

    /*     ... OP = inv[M]*A  and  B = M. */

    /*     ... Use mode 2 of ZNAUPD . */

    /* \BeginLib */

    /* \Routines called: */
    /*     znaupd   ARPACK reverse communication interface routine. */
    /*     zneupd   ARPACK routine that returns Ritz values and (optionally) */
    /*             Ritz vectors. */
    /*     zgttrf   LAPACK tridiagonal factorization routine. */
    /*     zgttrs   LAPACK tridiagonal solve routine. */
    /*     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully. */
    /*     zaxpy    Level 1 BLAS that computes y <- alpha*x+y. */
    /*     dznrm2   Level 1 BLAS that computes the norm of a vector. */
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
    /* FILE: ndrv3.F   SID: 2.4   DATE OF SID: 10/18/00   RELEASE: 2 */

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
    sigma.r = 0., sigma.i = 0.;

    /* --------------------------------------------------- */
    /* The matrix M is chosen to be the symmetric tri-     */
    /* diagonal matrix with 4 on the diagonal and 1 on the */
    /* off diagonals. It is factored by LAPACK subroutine  */
    /* zgttrf .                                             */
    /* --------------------------------------------------- */

    i__1 = n + 1;
    z__2.r = (double)i__1, z__2.i = 0.;
    ar_z_div(&z__1, &c_b2, &z__2);
    h.r = z__1.r, h.i = z__1.i;
    i__1 = n - 1;
    for (j = 1; j <= i__1; ++j)
    {
        i__2 = j - 1;
        z__1.r = h.r * 1. - h.i * 0., z__1.i = h.i * 1. + h.r * 0.;
        dl[i__2].r = z__1.r, dl[i__2].i = z__1.i;
        i__2 = j - 1;
        z__1.r = h.r * 4. - h.i * 0., z__1.i = h.r * 0. + h.i * 4.;
        dd[i__2].r = z__1.r, dd[i__2].i = z__1.i;
        i__2 = j - 1;
        z__1.r = h.r * 1. - h.i * 0., z__1.i = h.i * 1. + h.r * 0.;
        du[i__2].r = z__1.r, du[i__2].i = z__1.i;
        /* L20: */
    }
    i__1 = n - 1;
    z__1.r = h.r * 4. - h.i * 0., z__1.i = h.r * 0. + h.i * 4.;
    dd[i__1].r = z__1.r, dd[i__1].i = z__1.i;

    zgttrf_(&n, dl, dd, du, du2, ipiv, &ierr);
    if (ierr != 0)
    {
        printf(" \n");
        printf(" ERROR with _gttrf. \n");
        printf(" \n");
        return ierr;
    }

    /* --------------------------------------------------- */
    /* The work array WORKL is used in ZNAUPD  as           */
    /* workspace.  Its dimension LWORKL is set as          */
    /* illustrated below.  The parameter TOL determines    */
    /* the stopping criterion. If TOL<=0, machine          */
    /* precision is used.  The variable IDO is used for    */
    /* reverse communication, and is initially set to 0.   */
    /* Setting INFO=0 indicates that a random vector is    */
    /* generated in ZNAUPD  to start the Arnoldi iteration. */
    /* --------------------------------------------------- */

    /* Computing 2nd power */
    i__1 = ncv;
    lworkl = i__1 * i__1 * 3 + ncv * 5;
    tol = 0.f;
    ido = 0;
    info = 0;

    /* ------------------------------------------------- */
    /* This program uses exact shifts with respect to    */
    /* the current Hessenberg matrix (IPARAM(1) = 1).    */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 2 of ZNAUPD  is used     */
    /* (IPARAM(7) = 2).  All these options can be        */
    /* changed by the user. For details, see the         */
    /* documentation in ZNAUPD .                          */
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
    /* Repeatedly call the routine ZNAUPD  and take */
    /* actions indicated by parameter IDO until    */
    /* either convergence is indicated or maxitr   */
    /* has been exceeded.                          */
    /* ------------------------------------------- */

    znaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, rwork, &info);

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

        av_(&n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
        zgttrs_("N", &n, &c__1, dl, dd, du, du2, ipiv, &workd[ipntr[1] - 1], &n, &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" ERROR with _gttrs. \n");
            printf(" \n");
            return ierr;
        }

        /* --------------------------------------- */
        /* L O O P   B A C K to call ZNAUPD  again. */
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

        mv_(&n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

        /* --------------------------------------- */
        /* L O O P   B A C K to call ZNAUPD  again. */
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
        /* documentation in ZNAUPD . */
        /* ------------------------ */

        printf(" \n");
        printf(" Error with _naupd info = %d", info);
        printf(" Check the documentation of _naupd.\n");
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
            printf(" Error with _neupd info = %d", ierr);
            printf(" Check the documentation of _neupd\n");
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
                i__2 = j - 1;
                z__1.r = -d[i__2].r, z__1.i = -d[i__2].i;
                zaxpy_(&n, &z__1, mx, &c__1, ax, &c__1);
                i__2 = j - 1;
                rd[j - 1] = d[i__2].r;
                rd[j + 24] = d[j - 1].i;
                rd[j + 49] = dznrm2_(&n, ax, &c__1);
                rd[j + 49] /= dlapy2_(&rd[j - 1], &rd[j + 24]);
                /* L80: */
            }

            /* --------------------------- */
            /* Display computed residuals. */
            /* --------------------------- */

            dmout_(&nconv, &c__3, rd, &c__25, &c_n6,"Ritz values (Real, Imag) and relative residuals");
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
        printf("_NDRV3 \n");
        printf("====== \n");
        printf(" \n");
        printf(" Size of the matrix is %d", n);
        printf(" The number of Ritz values requested is %d", nev);
        printf(" The number of Arnoldi vectors generated \n");
        printf(" (NCV) is %d", ncv);
        printf(" What portion of the spectrum: %s", which);
        printf(" The number of converged Ritz values is %d", nconv);
        printf(" The number of Implicit Arnoldi update iterations taken is %d", iparam[2]);
        printf(" The number of OP*x is %d", iparam[8]);
        printf(" The convergence criterion is %e", tol);
        printf(" \n");
    }

    return 0;
} /* MAIN__ */

/* ========================================================================== */

/*     matrix vector multiplication subroutine */

int av_(a_int *n, a_dcomplex *v, a_dcomplex *w)
{
    /* System generated locals */
    a_int i__1, i__2, i__3, i__4, i__5;
    a_dcomplex z__1, z__2, z__3, z__4, z__5;

    /* Builtin functions */
    void ar_z_div(a_dcomplex *, a_dcomplex *, a_dcomplex *);

    /* Local variables */
    a_dcomplex h;
    a_int j;
    a_dcomplex s, dd, dl, du;

    /*     Compute the matrix vector multiplication y<---A*x */
    /*     where A is the stiffness matrix formed by using piecewise linear */
    /*     elements on [0,1]. */

    /* Parameter adjustments */
    --w;
    --v;

    i__1 = *n + 1;
    z__2.r = (double)i__1, z__2.i = 0.;
    ar_z_div(&z__1, &c_b2, &z__2);
    h.r = z__1.r, h.i = z__1.i;
    ar_z_div(&z__1, &c_b164, &c_b163);
    s.r = z__1.r, s.i = z__1.i;
    ar_z_div(&z__1, &c_b163, &h);
    dd.r = z__1.r, dd.i = z__1.i;
    z__3.r = -1., z__3.i = -0.;
    ar_z_div(&z__2, &z__3, &h);
    z__1.r = z__2.r - s.r, z__1.i = z__2.i - s.i;
    dl.r = z__1.r, dl.i = z__1.i;
    z__3.r = -1., z__3.i = -0.;
    ar_z_div(&z__2, &z__3, &h);
    z__1.r = z__2.r + s.r, z__1.i = z__2.i + s.i;
    du.r = z__1.r, du.i = z__1.i;

    z__2.r = dd.r * v[1].r - dd.i * v[1].i, z__2.i = dd.r * v[1].i + dd.i * v[1].r;
    z__3.r = du.r * v[2].r - du.i * v[2].i, z__3.i = du.r * v[2].i + du.i * v[2].r;
    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
    w[1].r = z__1.r, w[1].i = z__1.i;
    i__1 = *n - 1;
    for (j = 2; j <= i__1; ++j)
    {
        i__2 = j;
        i__3 = j - 1;
        z__3.r = dl.r * v[i__3].r - dl.i * v[i__3].i, z__3.i = dl.r * v[i__3].i + dl.i * v[i__3].r;
        i__4 = j;
        z__4.r = dd.r * v[i__4].r - dd.i * v[i__4].i, z__4.i = dd.r * v[i__4].i + dd.i * v[i__4].r;
        z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
        i__5 = j + 1;
        z__5.r = du.r * v[i__5].r - du.i * v[i__5].i, z__5.i = du.r * v[i__5].i + du.i * v[i__5].r;
        z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
        w[i__2].r = z__1.r, w[i__2].i = z__1.i;
        /* L10: */
    }
    i__1 = *n;
    i__2 = *n - 1;
    z__2.r = dl.r * v[i__2].r - dl.i * v[i__2].i, z__2.i = dl.r * v[i__2].i + dl.i * v[i__2].r;
    i__3 = *n;
    z__3.r = dd.r * v[i__3].r - dd.i * v[i__3].i, z__3.i = dd.r * v[i__3].i + dd.i * v[i__3].r;
    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
    w[i__1].r = z__1.r, w[i__1].i = z__1.i;
    return 0;
} /* av_ */

/* ------------------------------------------------------------------------ */
int mv_(a_int *n, a_dcomplex *v, a_dcomplex *w)
{
    /* System generated locals */
    a_int i__1, i__2, i__3, i__4, i__5;
    a_dcomplex z__1, z__2, z__3, z__4, z__5;

    /* Builtin functions */
    void ar_z_div(a_dcomplex *, a_dcomplex *, a_dcomplex *);

    /* Local variables */
    a_dcomplex h;
    a_int j;

    /*     Compute the matrix vector multiplication y<---M*x */
    /*     where M is the mass matrix formed by using piecewise linear elements */
    /*     on [0,1]. */

    /* Parameter adjustments */
    --w;
    --v;

    z__2.r = v[1].r * 4. - v[1].i * 0., z__2.i = v[1].i * 4. + v[1].r * 0.;
    z__3.r = v[2].r * 1. - v[2].i * 0., z__3.i = v[2].i * 1. + v[2].r * 0.;
    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
    w[1].r = z__1.r, w[1].i = z__1.i;
    i__1 = *n - 1;
    for (j = 2; j <= i__1; ++j)
    {
        i__2 = j;
        i__3 = j - 1;
        z__3.r = v[i__3].r * 1. - v[i__3].i * 0., z__3.i = v[i__3].i * 1. + v[i__3].r * 0.;
        i__4 = j;
        z__4.r = v[i__4].r * 4. - v[i__4].i * 0., z__4.i = v[i__4].i * 4. + v[i__4].r * 0.;
        z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
        i__5 = j + 1;
        z__5.r = v[i__5].r * 1. - v[i__5].i * 0., z__5.i = v[i__5].i * 1. + v[i__5].r * 0.;
        z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
        w[i__2].r = z__1.r, w[i__2].i = z__1.i;
        /* L10: */
    }
    i__1 = *n;
    i__2 = *n - 1;
    z__2.r = v[i__2].r * 1. - v[i__2].i * 0., z__2.i = v[i__2].i * 1. + v[i__2].r * 0.;
    i__3 = *n;
    z__3.r = v[i__3].r * 4. - v[i__3].i * 0., z__3.i = v[i__3].i * 4. + v[i__3].r * 0.;
    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
    w[i__1].r = z__1.r, w[i__1].i = z__1.i;

    i__1 = *n + 1;
    z__2.r = (double)i__1, z__2.i = 0.;
    ar_z_div(&z__1, &c_b2, &z__2);
    h.r = z__1.r, h.i = z__1.i;
    zscal_(n, &h, &w[1], &c__1);
    return 0;
} /* mv_ */

/* Main program alias */ int zndrv3_()
{
    MAIN__();
    return 0;
}
