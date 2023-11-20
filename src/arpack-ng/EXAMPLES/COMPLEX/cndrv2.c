/* EXAMPLES\COMPLEX\cndrv2.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

Extern struct
{
    a_fcomplex rho;
} convct_;

#define convct_1 convct_

/* Table of constant values */

static a_fcomplex c_b1 = {1.f, 0.f};
static a_fcomplex c_b3 = {2.f, 0.f};
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
    a_int i__1, i__2;
    a_fcomplex q__1, q__2, q__3, q__4;

    /* Local variables */
    a_fcomplex d[25], h;
    a_int j, n;
    a_fcomplex s, v[6400] /* was [256][25] */, h2, s1, s2, s3, dd[256], dl[256];
    float rd[75] /* was [25][3] */;
    a_fcomplex ax[256], du[256], du2[256];
    a_int ido, ncv, nev;
    float tol;
    char* bmat;
    a_int mode, info;
    a_bool rvec;
    a_int ierr, ipiv[256];
    a_fcomplex sigma;
    char* which;
    a_fcomplex resid[256];
    a_int nconv;
    a_fcomplex workd[768];
    a_int ipntr[14];
    a_fcomplex workl[2000];
    float rwork[256];
    a_int iparam[11];
    a_bool select[25];
    a_int ishfts, maxitr;
    a_int lworkl;
    a_fcomplex workev[50];

    /*     Simple program to illustrate the idea of reverse communication */
    /*     in shift-invert mode for a standard complex nonsymmetric eigenvalue */
    /*     problem. */

    /*     We implement example two of ex-complex.doc in DOCUMENTS directory */

    /* \Example-2 */
    /*     ... Suppose we want to solve A*x = lambda*x in shift-invert mode, */
    /*         where A is derived from the central difference discretization */
    /*         of the 1-dimensional convection-diffusion operator */
    /*                   (d^2u/dx^2) + rho*(du/dx) */
    /*         on the interval [0,1] with zero Dirichlet boundary condition. */
    /*     ... The shift sigma is a complex number. */

    /*     ... OP = inv[A-sigma*I] and  B = I. */

    /*     ... Use mode 3 of CNAUPD. */

    /* \BeginLib */

    /* \Routines called: */
    /*     cnaupd  ARPACK reverse communication interface routine. */
    /*     cneupd  ARPACK routine that returns Ritz values and (optionally) */
    /*             Ritz vectors. */
    /*     cgttrf  LAPACK tridiagonal factorization routine. */
    /*     cgttrs  LAPACK tridiagonal solve routine. */
    /*     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully. */
    /*     caxpy   Level 1 BLAS that computes y <- alpha*x+y. */
    /*     ccopy   Level 1 BLAS that copies one vector to another. */
    /*     scnrm2  Level 1 BLAS that computes the norm of a vector. */
    /*     av      Matrix vector multiplication routine that computes A*x. */

    /* \Author */
    /*     Richard Lehoucq */
    /*     Danny Sorensen */
    /*     Chao Yang */
    /*     Dept. of Computational & */
    /*     Applied Mathematics */
    /*     Rice University */
    /*     Houston, Texas */

    /* \SCCS Information: @(#) */
    /* FILE: ndrv2.F   SID: 2.6   DATE OF SID: 10/18/00   RELEASE: 2 */

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

    /* ------------------------------------------------ */
    /* The number N is the dimension of the matrix.  A  */
    /* standard eigenvalue problem is solved (BMAT =    */
    /* 'I').  NEV is the number of eigenvalues (closest */
    /* to the shift SIGMA) to be approximated.  Since   */
    /* the shift-invert mode is used, WHICH is set to   */
    /* 'LM'.  The user can modify NEV, NCV, SIGMA to    */
    /* solve problems of different sizes, and to get    */
    /* different parts of the spectrum.  However, The   */
    /* following conditions must be satisfied:          */
    /*                 N <= MAXN,                       */
    /*               NEV <= MAXNEV,                     */
    /*           NEV + 2 <= NCV <= MAXNCV               */
    /* ------------------------------------------------ */

    n = 100;
    nev = 4;
    ncv = 20;
    if (n > 256)
    {
        printf(" ERROR with _NDRV2: N is greater than MAXN \n");
        return ierr;
    }
    else if (nev > 10)
    {
        printf(" ERROR with _NDRV2: NEV is greater than MAXNEV \n");
        return ierr;
    }
    else if (ncv > 25)
    {
        printf(" ERROR with _NDRV2: NCV is greater than MAXNCV \n");
        return ierr;
    }
    bmat = "I";
    which = "LM";
    sigma.r = 0.f, sigma.i = 0.f;

    /* -------------------------------------------------- */
    /* Construct C = A - SIGMA*I, factor C in complex     */
    /* arithmetic (using LAPACK subroutine cgttrf). The   */
    /* matrix A is chosen to be the tridiagonal matrix    */
    /* derived from standard central difference of the    */
    /* 1-d convection diffusion operator - u``+ rho*u` on */
    /* the interval [0, 1] with zero Dirichlet boundary   */
    /* condition.                                         */
    /* -------------------------------------------------- */

    convct_1.rho.r = 10.f, convct_1.rho.i = 0.f;
    i__1 = n + 1;
    q__2.r = (float)i__1, q__2.i = 0.f;
    ar_c_div(&q__1, &c_b1, &q__2);
    h.r = q__1.r, h.i = q__1.i;
    q__1.r = h.r * h.r - h.i * h.i, q__1.i = h.r * h.i + h.i * h.r;
    h2.r = q__1.r, h2.i = q__1.i;
    ar_c_div(&q__1, &convct_1.rho, &c_b3);
    s.r = q__1.r, s.i = q__1.i;

    q__3.r = -1.f, q__3.i = -0.f;
    ar_c_div(&q__2, &q__3, &h2);
    ar_c_div(&q__4, &s, &h);
    q__1.r = q__2.r - q__4.r, q__1.i = q__2.i - q__4.i;
    s1.r = q__1.r, s1.i = q__1.i;
    ar_c_div(&q__2, &c_b3, &h2);
    q__1.r = q__2.r - sigma.r, q__1.i = q__2.i - sigma.i;
    s2.r = q__1.r, s2.i = q__1.i;
    q__3.r = -1.f, q__3.i = -0.f;
    ar_c_div(&q__2, &q__3, &h2);
    ar_c_div(&q__4, &s, &h);
    q__1.r = q__2.r + q__4.r, q__1.i = q__2.i + q__4.i;
    s3.r = q__1.r, s3.i = q__1.i;

    i__1 = n - 1;
    for (j = 1; j <= i__1; ++j)
    {
        i__2 = j - 1;
        dl[i__2].r = s1.r, dl[i__2].i = s1.i;
        i__2 = j - 1;
        dd[i__2].r = s2.r, dd[i__2].i = s2.i;
        i__2 = j - 1;
        du[i__2].r = s3.r, du[i__2].i = s3.i;
        /* L10: */
    }
    i__1 = n - 1;
    dd[i__1].r = s2.r, dd[i__1].i = s2.i;

    cgttrf_(&n, dl, dd, du, du2, ipiv, &ierr);
    if (ierr != 0)
    {
        printf(" \n");
        printf(" ERROR with _gttrf in _NDRV2.\n");
        printf(" \n");
        return ierr;
    }

    /* --------------------------------------------------- */
    /* The work array WORKL is used in CNAUPD as           */
    /* workspace.  Its dimension LWORKL is set as          */
    /* illustrated below.  The parameter TOL determines    */
    /* the stopping criterion. If TOL<=0, machine          */
    /* precision is used.  The variable IDO is used for    */
    /* reverse communication, and is initially set to 0.   */
    /* Setting INFO=0 indicates that a random vector is    */
    /* generated in CNAUPD to start the Arnoldi iteration. */
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
    /* iterations allowed. Mode 3 of CNAUPD is used      */
    /* (IPARAM(7) = 3).  All these options can be        */
    /* changed by the user. For details see the          */
    /* documentation in CNAUPD.                          */
    /* ------------------------------------------------- */

    ishfts = 1;
    maxitr = 300;
    mode = 3;

    iparam[0] = ishfts;
    iparam[2] = maxitr;
    iparam[6] = mode;

    /* ----------------------------------------- */
    /* M A I N   L O O P (Reverse communication) */
    /* ----------------------------------------- */

L20:

    /* ------------------------------------------- */
    /* Repeatedly call the routine CNAUPD and take */
    /* actions indicated by parameter IDO until    */
    /* either convergence is indicated or maxitr   */
    /* has been exceeded.                          */
    /* ------------------------------------------- */

    cnaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, rwork, &info);

    if (ido == -1 || ido == 1)
    {

        /* ----------------------------------------- */
        /* Perform  y <--- OP*x = inv[A-SIGMA*I]*x   */
        /* The user should supply his/her own linear */
        /* system solver here that takes             */
        /* workd(ipntr(1)) as the input, and returns */
        /* the result to workd(ipntr(2)).            */
        /* ----------------------------------------- */

        ccopy_(&n, &workd[ipntr[0] - 1], &c__1, &workd[ipntr[1] - 1], &c__1);

        cgttrs_("N", &n, &c__1, dl, dd, du, du2, ipiv, &workd[ipntr[1] - 1], &n, &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" ERROR with _gttrs in _NDRV2.\n");
            printf(" \n");
            return ierr;
        }

        /* --------------------------------------- */
        /* L O O P   B A C K to call CNAUPD again. */
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
        /* documentation in CNAUPD  */
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
        /* Post-Process using CNEUPD.                */
        /*                                           */
        /* Computed eigenvalues may be extracted.    */
        /*                                           */
        /* Eigenvectors may also be computed now if  */
        /* desired.  (indicated by rvec = .true.)    */
        /* ----------------------------------------- */

        rvec = TRUE_;

        cneupd_(&rvec, "A", select, d, v, &c__256, &sigma, workev, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, rwork, &ierr);

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
            /* Check the documentation of CNEUPD. */
            /* ---------------------------------- */

            printf(" \n");
            printf(" Error with _neupd info = %d", ierr);
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

                av_(&n, &v[(j << 8) - 256], ax);
                i__2 = j - 1;
                q__1.r = -d[i__2].r, q__1.i = -d[i__2].i;
                caxpy_(&n, &q__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
                i__2 = j - 1;
                rd[j - 1] = d[i__2].r;
                rd[j + 24] = d[j - 1].i;
                rd[j + 49] = scnrm2_(&n, ax, &c__1);
                rd[j + 49] /= slapy2_(&rd[j - 1], &rd[j + 24]);
                /* L60: */
            }

            /* --------------------------- */
            /* Display computed residuals. */
            /* --------------------------- */

            smout_(&nconv, &c__3, rd, &c__25, &c_n6,"Ritz values (Real, Imag) and relative residuals");
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
        printf("_NDRV2 \n");
        printf("====== \n");
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
    /* Done with program cndrv2. */
    /* ------------------------- */

    return 0;
}

/* ------------------------------------------------------------------- */

/*     matrix vector multiplication subroutine */

int av_(a_int *n, a_fcomplex *v, a_fcomplex *w)
{
    /* System generated locals */
    a_int i__1, i__2, i__3, i__4, i__5;
    a_fcomplex q__1, q__2, q__3, q__4, q__5;

    /* Builtin functions */
    void ar_c_div(a_fcomplex *, a_fcomplex *, a_fcomplex *);

    /* Local variables */
    a_fcomplex h;
    a_int j;
    a_fcomplex s, h2, dd, dl, du;

    /* Parameter adjustments */
    --w;
    --v;

    i__1 = *n + 1;
    q__2.r = (float)i__1, q__2.i = 0.f;
    ar_c_div(&q__1, &c_b1, &q__2);
    h.r = q__1.r, h.i = q__1.i;
    q__1.r = h.r * h.r - h.i * h.i, q__1.i = h.r * h.i + h.i * h.r;
    h2.r = q__1.r, h2.i = q__1.i;
    ar_c_div(&q__1, &convct_1.rho, &c_b3);
    s.r = q__1.r, s.i = q__1.i;
    ar_c_div(&q__1, &c_b3, &h2);
    dd.r = q__1.r, dd.i = q__1.i;
    q__3.r = -1.f, q__3.i = -0.f;
    ar_c_div(&q__2, &q__3, &h2);
    ar_c_div(&q__4, &s, &h);
    q__1.r = q__2.r - q__4.r, q__1.i = q__2.i - q__4.i;
    dl.r = q__1.r, dl.i = q__1.i;
    q__3.r = -1.f, q__3.i = -0.f;
    ar_c_div(&q__2, &q__3, &h2);
    ar_c_div(&q__4, &s, &h);
    q__1.r = q__2.r + q__4.r, q__1.i = q__2.i + q__4.i;
    du.r = q__1.r, du.i = q__1.i;

    q__2.r = dd.r * v[1].r - dd.i * v[1].i, q__2.i = dd.r * v[1].i + dd.i * v[1].r;
    q__3.r = du.r * v[2].r - du.i * v[2].i, q__3.i = du.r * v[2].i + du.i * v[2].r;
    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
    w[1].r = q__1.r, w[1].i = q__1.i;
    i__1 = *n - 1;
    for (j = 2; j <= i__1; ++j)
    {
        i__2 = j;
        i__3 = j - 1;
        q__3.r = dl.r * v[i__3].r - dl.i * v[i__3].i, q__3.i = dl.r * v[i__3].i + dl.i * v[i__3].r;
        i__4 = j;
        q__4.r = dd.r * v[i__4].r - dd.i * v[i__4].i, q__4.i = dd.r * v[i__4].i + dd.i * v[i__4].r;
        q__2.r = q__3.r + q__4.r, q__2.i = q__3.i + q__4.i;
        i__5 = j + 1;
        q__5.r = du.r * v[i__5].r - du.i * v[i__5].i, q__5.i = du.r * v[i__5].i + du.i * v[i__5].r;
        q__1.r = q__2.r + q__5.r, q__1.i = q__2.i + q__5.i;
        w[i__2].r = q__1.r, w[i__2].i = q__1.i;
        /* L10: */
    }
    i__1 = *n;
    i__2 = *n - 1;
    q__2.r = dl.r * v[i__2].r - dl.i * v[i__2].i, q__2.i = dl.r * v[i__2].i + dl.i * v[i__2].r;
    i__3 = *n;
    q__3.r = dd.r * v[i__3].r - dd.i * v[i__3].i, q__3.i = dd.r * v[i__3].i + dd.i * v[i__3].r;
    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
    w[i__1].r = q__1.r, w[i__1].i = q__1.i;
    return 0;
} /* av_ */
