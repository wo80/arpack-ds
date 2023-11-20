/* EXAMPLES\BAND\znbdr3.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

/* Table of constant values */

static a_dcomplex c_b1 = {1., 0.};
static a_dcomplex c_b2 = {0., 0.};
static a_dcomplex c_b3 = {2., 0.};
static a_int c__9 = 9;
static a_int c__1 = 1;
static a_int c__50 = 50;
static a_int c__1000 = 1000;
static a_int c__3 = 3;
static a_int c__5 = 5;
static a_int c__6 = 6;
static a_int c_n6 = -6;

/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    a_int i__1, i__2;
    a_dcomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    a_int s_wsle(cilist *), do_lio(a_int *, a_int *, char *, ftnlen), e_wsle(void);

    /* Local variables */
    a_dcomplex a[50000] /* was [50][1000] */, d[50], h;
    a_int j;
    a_dcomplex m[50000] /* was [50][1000] */;
    a_int n;
    a_dcomplex v[50000] /* was [1000][50] */;
    double rd[150] /* was [50][3] */;
    a_int kl;
    a_dcomplex ax[1000];
    a_int ku;
    a_dcomplex mx[1000], fac[50000] /* was [50][1000] */;
    a_int ncv, nev;
    a_dcomplex rho;
    double tol;
    char bmat[1];
    a_int mode, info;
    a_bool rvec;
    a_int isub, isup, idiag;
    a_dcomplex sigma;
    char which[2];
    a_dcomplex resid[1000];
    a_int nconv;
    a_dcomplex workd[3000];
    a_int iwork[1000];
    a_dcomplex workl[7750];
    double rwork[1000];
    a_int iparam[11];
    a_bool select[50];
    a_int maxitr, lworkl;
    a_dcomplex workev[100];

    /* Fortran I/O blocks */
    static cilist io___4 = {0, 6, 0, 0, 0};
    static cilist io___5 = {0, 6, 0, 0, 0};
    static cilist io___6 = {0, 6, 0, 0, 0};
    static cilist io___38 = {0, 6, 0, 0, 0};
    static cilist io___39 = {0, 6, 0, 0, 0};
    static cilist io___40 = {0, 6, 0, 0, 0};
    static cilist io___41 = {0, 6, 0, 0, 0};
    static cilist io___42 = {0, 6, 0, 0, 0};
    static cilist io___43 = {0, 6, 0, 0, 0};
    static cilist io___44 = {0, 6, 0, 0, 0};
    static cilist io___45 = {0, 6, 0, 0, 0};
    static cilist io___46 = {0, 6, 0, 0, 0};
    static cilist io___47 = {0, 6, 0, 0, 0};
    static cilist io___48 = {0, 6, 0, 0, 0};
    static cilist io___49 = {0, 6, 0, 0, 0};
    static cilist io___50 = {0, 6, 0, 0, 0};
    static cilist io___54 = {0, 6, 0, 0, 0};
    static cilist io___55 = {0, 6, 0, 0, 0};
    static cilist io___56 = {0, 6, 0, 0, 0};
    static cilist io___57 = {0, 6, 0, 0, 0};

    /*     ... Construct matrices A and M in LAPACK-style band form. */
    /*         Matrices A and M are derived from the finite */
    /*         element discretization of the 1-dimensional */
    /*         convection-diffusion operator */
    /*                         (d^2u/dx^2) + rho*(du/dx) */
    /*         on the interval [0,1] with zero boundary condition using */
    /*         piecewise linear elements. */

    /*     ... Call ZNBAND  to find eigenvalues LAMBDA such that */
    /*                    A*x = M*x*LAMBDA. */

    /*     ... Eigenvalues with largest real parts are sought. */

    /*     ... Use mode 2 of ZNAUPD . */

    /* \BeginLib */

    /* \Routines called: */
    /*     znband   ARPACK banded eigenproblem solver. */
    /*     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully. */
    /*     zlaset   LAPACK routine to initialize a matrix to zero. */
    /*     zaxpy    Level 1 BLAS that computes y <- alpha*x+y. */
    /*     dznrm2   Level 1 BLAS that computes the norm of a vector. */
    /*     zgbmv    Level 2 BLAS that computes the band matrix vector product. */

    /* \Author */
    /*     Richard Lehoucq */
    /*     Danny Sorensen */
    /*     Chao Yang */
    /*     Dept. of Computational & */
    /*     Applied Mathematics */
    /*     Rice University */
    /*     Houston, Texas */

    /* \SCCS Information: @(#) */
    /* FILE: nbdr3.F   SID: 2.4   DATE OF SID: 10/20/00   RELEASE: 2 */

    /* \Remarks */
    /*     1. None */

    /* \EndLib */

    /* ------------------------------------------------------------------------- */

    /* ----------------------------------- */
    /* Define leading dimensions for all   */
    /* arrays.                             */
    /* MAXN   - Maximum size of the matrix */
    /* MAXNEV - Maximum number of          */
    /*          eigenvalues to be computed */
    /* MAXNCV - Maximum number of Arnoldi  */
    /*          vectors stored             */
    /* MAXBDW - Maximum bandwidth          */
    /* ----------------------------------- */

    /* ----------------------------------------------- */
    /* The number NX is the number of interior points  */
    /* in the discretization of the 2-dimensional      */
    /* convection-diffusion operator on the unit       */
    /* square with zero Dirichlet boundary condition.  */
    /* The number N(=NX*NX) is the dimension of the    */
    /* matrix.  A generalized eigenvalue problem is    */
    /* solved (BMAT = 'G').  NEV is the number of      */
    /* eigenvalues to be approximated. The user can    */
    /* modify NX, NEV, NCV and WHICH to solve problems */
    /* of different sizes, and to get different parts  */
    /* the spectrum.  However, The following           */
    /* conditions must be satisfied:                   */
    /*                   N <= MAXN                     */
    /*                 NEV <= MAXNEV                   */
    /*           NEV + 2 <= NCV <= MAXNCV              */
    /* ----------------------------------------------- */

    n = 100;
    nev = 4;
    ncv = 10;
    if (n > 1000)
    {
        s_wsle(&io___4);
        do_lio(&c__9, &c__1, " ERROR with _NBDR3: N is greater than MAXN ", (ftnlen)43);
        e_wsle();
        goto L9000;
    }
    else if (nev > 25)
    {
        s_wsle(&io___5);
        do_lio(&c__9, &c__1, " ERROR with _NBDR3: NEV is greater than MAXNEV ", (ftnlen)47);
        e_wsle();
        goto L9000;
    }
    else if (ncv > 50)
    {
        s_wsle(&io___6);
        do_lio(&c__9, &c__1, " ERROR with _NBDR3: NCV is greater than MAXNCV ", (ftnlen)47);
        e_wsle();
        goto L9000;
    }
    *bmat = 'G';
    strcpy(which, "LM");
    sigma.r = 0., sigma.i = 0.;

    /* -------------------------------------------------- */
    /* The work array WORKL is used in ZNAUPD  as          */
    /* workspace.  Its dimension LWORKL has to be set as  */
    /* illustrated below.  The parameter TOL determines   */
    /* the stopping criterion. If TOL<=0, machine machine */
    /* precision is used.  Setting INFO=0 indicates that  */
    /* using a randomly generated vector to start the     */
    /* the ARNOLDI process.                               */
    /* -------------------------------------------------- */

    /* Computing 2nd power */
    i__1 = ncv;
    lworkl = i__1 * i__1 * 3 + ncv * 5;
    info = 0;
    tol = 0.f;

    /* ------------------------------------------------- */
    /* IPARAm(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 2 of ZNAUPD  is used     */
    /* (IPARAm(7) = 2). All these options can be changed */
    /* by the user. For details, see the documentation   */
    /* in znband .                                        */
    /* ------------------------------------------------- */

    maxitr = 300;
    mode = 2;

    iparam[2] = maxitr;
    iparam[6] = mode;

    /* ------------------------------------------ */
    /* Construct matrices A and M in LAPACK-style */
    /* banded form.                               */
    /* ------------------------------------------ */

    /* ------------------------------------------- */
    /* Zero out the workspace for banded matrices. */
    /* ------------------------------------------- */

    zlaset_("A", &c__50, &n, &c_b2, &c_b2, a, &c__50);
    zlaset_("A", &c__50, &n, &c_b2, &c_b2, m, &c__50);
    zlaset_("A", &c__50, &n, &c_b2, &c_b2, fac, &c__50);

    /* ----------------------------------- */
    /* KU, KL are number of superdiagonals */
    /* and subdiagonals within the band of */
    /* matrices A and M.                   */
    /* ----------------------------------- */

    kl = 1;
    ku = 1;

    /* ------------- */
    /* Main diagonal */
    /* ------------- */

    i__1 = n + 1;
    z__2.r = (double)i__1, z__2.i = 0.;
    ar_z_div(&z__1, &c_b1, &z__2);
    h.r = z__1.r, h.i = z__1.i;

    idiag = kl + ku + 1;
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
        i__2 = idiag + j * 50 - 51;
        ar_z_div(&z__1, &c_b3, &h);
        a[i__2].r = z__1.r, a[i__2].i = z__1.i;
        i__2 = idiag + j * 50 - 51;
        z__1.r = h.r * 4. - h.i * 0., z__1.i = h.r * 0. + h.i * 4.;
        m[i__2].r = z__1.r, m[i__2].i = z__1.i;
        /* L30: */
    }

    /* ----------------------------------- */
    /* First subdiagonal and superdiagonal */
    /* ----------------------------------- */

    isup = kl + ku;
    isub = kl + ku + 2;
    rho.r = 10., rho.i = 0.;
    i__1 = n - 1;
    for (j = 1; j <= i__1; ++j)
    {
        i__2 = isup + (j + 1) * 50 - 51;
        z__3.r = -1., z__3.i = -0.;
        ar_z_div(&z__2, &z__3, &h);
        ar_z_div(&z__4, &rho, &c_b3);
        z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
        a[i__2].r = z__1.r, a[i__2].i = z__1.i;
        i__2 = isub + j * 50 - 51;
        z__3.r = -1., z__3.i = -0.;
        ar_z_div(&z__2, &z__3, &h);
        ar_z_div(&z__4, &rho, &c_b3);
        z__1.r = z__2.r - z__4.r, z__1.i = z__2.i - z__4.i;
        a[i__2].r = z__1.r, a[i__2].i = z__1.i;
        i__2 = isup + (j + 1) * 50 - 51;
        z__1.r = h.r * 1. - h.i * 0., z__1.i = h.i * 1. + h.r * 0.;
        m[i__2].r = z__1.r, m[i__2].i = z__1.i;
        i__2 = isub + j * 50 - 51;
        z__1.r = h.r * 1. - h.i * 0., z__1.i = h.i * 1. + h.r * 0.;
        m[i__2].r = z__1.r, m[i__2].i = z__1.i;
        /* L40: */
    }

    /* --------------------------------------------- */
    /* Call ARPACK banded solver to find eigenvalues */
    /* and eigenvectors. Eigenvalues are returned in */
    /* the one dimensional array D.  Eigenvectors    */
    /* are returned in the first NCONV (=IPARAM(5))  */
    /* columns of V.                                 */
    /* --------------------------------------------- */

    rvec = TRUE_;
    znband_(&rvec, "A", select, d, v, &c__1000, &sigma, workev, &n, a, m, &c__50, fac, &kl, &ku, which, bmat, &nev, &tol, resid, &ncv, v, &c__1000, iparam, workd, workl, &lworkl, rwork, iwork, &info);

    if (info == 0)
    {

        nconv = iparam[4];

        /* --------------------------------- */
        /* Print out convergence information */
        /* --------------------------------- */

        s_wsle(&io___38);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
        s_wsle(&io___39);
        do_lio(&c__9, &c__1, "_NBDR3 ", (ftnlen)7);
        e_wsle();
        s_wsle(&io___40);
        do_lio(&c__9, &c__1, "====== ", (ftnlen)7);
        e_wsle();
        s_wsle(&io___41);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
        s_wsle(&io___42);
        do_lio(&c__9, &c__1, " The size of the matrix is ", (ftnlen)27);
        do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___43);
        do_lio(&c__9, &c__1, " Number of eigenvalue requested is ", (ftnlen)35);
        do_lio(&c__3, &c__1, (char *)&nev, (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___44);
        do_lio(&c__9, &c__1, " The number of Arnoldi vectors generated", (ftnlen)40);
        do_lio(&c__9, &c__1, " (NCV) is ", (ftnlen)10);
        do_lio(&c__3, &c__1, (char *)&ncv, (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___45);
        do_lio(&c__9, &c__1, " The number of converged Ritz values is ", (ftnlen)40);
        do_lio(&c__3, &c__1, (char *)&nconv, (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___46);
        do_lio(&c__9, &c__1, " What portion of the spectrum ", (ftnlen)30);
        do_lio(&c__9, &c__1, which, (ftnlen)2);
        e_wsle();
        s_wsle(&io___47);
        do_lio(&c__9, &c__1, " The number of Implicit Arnoldi ", (ftnlen)32);
        do_lio(&c__9, &c__1, " update taken is ", (ftnlen)17);
        do_lio(&c__3, &c__1, (char *)&iparam[2], (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___48);
        do_lio(&c__9, &c__1, " The number of OP*x is ", (ftnlen)23);
        do_lio(&c__3, &c__1, (char *)&iparam[8], (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___49);
        do_lio(&c__9, &c__1, " The convergence tolerance is ", (ftnlen)30);
        do_lio(&c__5, &c__1, (char *)&tol, (ftnlen)sizeof(double));
        e_wsle();
        s_wsle(&io___50);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();

        i__1 = nconv;
        for (j = 1; j <= i__1; ++j)
        {

            /* -------------------------- */
            /* Compute the residual norm. */
            /*    ||  A*x - lambda*x ||   */
            /* -------------------------- */

            zgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b1, &a[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b2, ax, &c__1);
            zgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b1, &m[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b2, mx, &c__1);
            i__2 = j - 1;
            z__1.r = -d[i__2].r, z__1.i = -d[i__2].i;
            zaxpy_(&n, &z__1, mx, &c__1, ax, &c__1);
            i__2 = j - 1;
            rd[j - 1] = d[i__2].r;
            rd[j + 49] = d[j - 1].i;
            rd[j + 99] = dznrm2_(&n, ax, &c__1);
            rd[j + 99] /= dlapy2_(&rd[j - 1], &rd[j + 49]);
            /* L50: */
        }
        dmout_(&c__6, &nconv, &c__3, rd, &c__50, &c_n6,"Ritz values (Real,Imag) and relative residuals");
    }
    else
    {

        /* ----------------------------------- */
        /* Either convergence failed, or there */
        /* is error.  Check the documentation  */
        /* for znband .                         */
        /* ----------------------------------- */

        s_wsle(&io___54);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
        s_wsle(&io___55);
        do_lio(&c__9, &c__1, " Error with _band, info= ", (ftnlen)25);
        do_lio(&c__3, &c__1, (char *)&info, (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___56);
        do_lio(&c__9, &c__1, " Check the documentation of _band ", (ftnlen)34);
        e_wsle();
        s_wsle(&io___57);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
    }

L9000:
    return 0;
} /* MAIN__ */

/* Main program alias */ int znbdr3_()
{
    MAIN__();
    return 0;
}
