/* EXAMPLES\BAND\dsbdr5.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

/* Table of constant values */

static a_int c__9 = 9;
static a_int c__1 = 1;
static a_int c__50 = 50;
static double c_b15 = 0.;
static a_int c__1000 = 1000;
static a_int c__3 = 3;
static a_int c__5 = 5;
static double c_b97 = 1.;
static a_int c__6 = 6;
static a_int c__2 = 2;
static a_int c_n6 = -6;

/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    a_int i__1;
    double d__1;

    /* Builtin functions */
    a_int s_wsle(cilist *), do_lio(a_int *, a_int *, char *, ftnlen), e_wsle(void);
    int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    double a[50000] /* was [50][1000] */, d[100] /* was [50][2] */, h;
    a_int j;
    double m[50000] /* was [50][1000] */;
    a_int n;
    double v[50000] /* was [1000][50] */, r1, r2;
    a_int kl;
    double ax[1000];
    a_int ku;
    double mx[1000];
    a_int ido, ncv, nev;
    double tol, rfac[50000] /* was [50][1000] */;
    char bmat[1];
    a_int mode, info;
    a_bool rvec;
    a_int isub, isup;
    extern double dnrm2_(a_int *, double *, a_int *);
    a_int idiag;
    extern int dgbmv_(char *, a_int *, a_int *, a_int *, a_int *, double *, double *, a_int *, double *, a_int *, double *, double *, a_int *, ftnlen);
    double sigma;
    char which[2];
    double resid[1000];
    a_int nconv;
    extern int daxpy_(a_int *, double *, double *, a_int *, double *, a_int *);
    double workd[3000];
    extern int dmout_(a_int *, a_int *, a_int *, double *, a_int *, a_int *, char *, ftnlen);
    a_int iwork[1000];
    double workl[2900];
    extern int dsband_(a_bool *, char *, a_bool *, double *, double *, a_int *, double *, a_int *, double *, double *, a_int *, double *, a_int *, a_int *, char *, char *, a_int *, double *, double *, a_int *, double *, a_int *, a_int *, double *, double *, a_int *, a_int *, a_int *, ftnlen, ftnlen, ftnlen);
    a_int iparam[11];
    extern int dlaset_(char *, a_int *, a_int *, double *, double *, double *, a_int *, ftnlen);
    a_bool select[50];
    a_int maxitr, lworkl;

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
    static cilist io___53 = {0, 6, 0, 0, 0};
    static cilist io___54 = {0, 6, 0, 0, 0};
    static cilist io___55 = {0, 6, 0, 0, 0};
    static cilist io___56 = {0, 6, 0, 0, 0};

    /*     ... Construct the matrix A in LAPACK-style band form. */
    /*         The matrix A is the 1-dimensional discrete Laplacian on [0,1] */
    /*         with zero Dirichlet boundary condition, KG is the mass */
    /*         formed by using piecewise linear elements on [0,1]. */

    /*     ... Call DSBAND  with Buckling mode to find eigenvalues LAMBDA */
    /*         such that */
    /*                          A*x = M*x*LAMBDA. */

    /*     ... Use mode 4 of DSAUPD . */

    /* \BeginLib */

    /* \Routines called: */
    /*     dsband   ARPACK banded eigenproblem solver. */
    /*     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully. */
    /*     dlaset   LAPACK routine to initialize a matrix to zero. */
    /*     daxpy    Level 1 BLAS that computes y <- alpha*x+y. */
    /*     dnrm2    Level 1 BLAS that computes the norm of a vector. */
    /*     dgbmv    Level 2 BLAS that computes the band matrix vector product */

    /* \Author */
    /*     Richard Lehoucq */
    /*     Danny Sorensen */
    /*     Chao Yang */
    /*     Dept. of Computational & */
    /*     Applied Mathematics */
    /*     Rice University */
    /*     Houston, Texas */

    /* \SCCS Information: @(#) */
    /* FILE: sbdr5.F   SID: 2.5   DATE OF SID: 08/26/96   RELEASE: 2 */

    /* \Remarks */
    /*     1. None */

    /* \EndLib */

    /* ---------------------------------------------------------------------- */

    /*     %-------------------------------------% */
    /*     | Define leading dimensions for all   | */
    /*     | arrays.                             | */
    /*     | MAXN   - Maximum size of the matrix | */
    /*     | MAXNEV - Maximum number of          | */
    /*     |          eigenvalues to be computed | */
    /*     | MAXNCV - Maximum number of Arnoldi  | */
    /*     |          vectors stored             | */
    /*     | MAXBDW - Maximum bandwidth          | */
    /*     %-------------------------------------% */

    /*     %--------------% */
    /*     | Local Arrays | */
    /*     %--------------% */

    /*     %---------------% */
    /*     | Local Scalars | */
    /*     %---------------% */

    /*     %------------% */
    /*     | Parameters | */
    /*     %------------% */

    /*     %-----------------------------% */
    /*     | BLAS & LAPACK routines used | */
    /*     %-----------------------------% */

    /*     %--------------------% */
    /*     | Intrinsic function | */
    /*     %--------------------% */

    /*     %-----------------------% */
    /*     | Executable Statements | */
    /*     %-----------------------% */

    /*     %--------------------------------------------------% */
    /*     | The number N is the dimension of the matrix.  A  | */
    /*     | generalized eigenvalue problem is solved         | */
    /*     | (BMAT = 'G').  NEV is the number of eigenvalues  | */
    /*     | to be approximated.  Since the Buckling mode is  | */
    /*     | is used, WHICH is set to 'LM'.  The user can     | */
    /*     | modify N, NEV, NCV and SIGMA to solve problems   | */
    /*     | of different sizes, and to get different parts   | */
    /*     | the spectrum.  However, the following conditions | */
    /*     | must be satisfied:                               | */
    /*     |                   N <= MAXN                      | */
    /*     |                 NEV <= MAXNEV                    | */
    /*     |           NEV + 1 <= NCV <= MAXNCV               | */
    /*     %--------------------------------------------------% */

    n = 100;
    nev = 4;
    ncv = 10;
    if (n > 1000)
    {
        s_wsle(&io___4);
        do_lio(&c__9, &c__1, " ERROR with _SBDR5: N is greater than MAXN ", (ftnlen)43);
        e_wsle();
        goto L9000;
    }
    else if (nev > 25)
    {
        s_wsle(&io___5);
        do_lio(&c__9, &c__1, " ERROR with _SBDR5: NEV is greater than MAXNEV ", (ftnlen)47);
        e_wsle();
        goto L9000;
    }
    else if (ncv > 50)
    {
        s_wsle(&io___6);
        do_lio(&c__9, &c__1, " ERROR with _SBDR5: NCV is greater than MAXNCV ", (ftnlen)47);
        e_wsle();
        goto L9000;
    }
    *(unsigned char *)bmat = 'G';
    s_copy(which, "LM", (ftnlen)2, (ftnlen)2);
    sigma = 1.f;

    /*     %-----------------------------------------------------% */
    /*     | The work array WORKL is used in DSAUPD  as           | */
    /*     | workspace.  Its dimension LWORKL is set as          | */
    /*     | illustrated below.  The parameter TOL determines    | */
    /*     | the stopping criterion. If TOL<=0, machine          | */
    /*     | precision is used.  The variable IDO is used for    | */
    /*     | reverse communication, and is initially set to 0.   | */
    /*     | Setting INFO=0 indicates that a random vector is    | */
    /*     | generated in DSAUPD  to start the Arnoldi iteration. | */
    /*     %-----------------------------------------------------% */

    /* Computing 2nd power */
    i__1 = ncv;
    lworkl = i__1 * i__1 + (ncv << 3);
    tol = 0.;
    ido = 0;
    info = 0;

    /*     %---------------------------------------------------% */
    /*     | IPARAM(3) specifies the maximum number of Arnoldi | */
    /*     | iterations allowed.  Mode 4 of DSAUPD  is used     | */
    /*     | (IPARAM(7) = 4). All these options can be changed | */
    /*     | by the user. For details see the documentation in | */
    /*     | DSBAND .                                           | */
    /*     %---------------------------------------------------% */

    maxitr = 300;
    mode = 4;

    iparam[2] = maxitr;
    iparam[6] = mode;

    /*     %----------------------------------------% */
    /*     | Construct the matrix A in LAPACK-style | */
    /*     | banded form.                           | */
    /*     %----------------------------------------% */

    /*     %---------------------------------------------% */
    /*     | Zero out the workspace for banded matrices. | */
    /*     %---------------------------------------------% */

    dlaset_("A", &c__50, &n, &c_b15, &c_b15, a, &c__50, (ftnlen)1);
    dlaset_("A", &c__50, &n, &c_b15, &c_b15, m, &c__50, (ftnlen)1);
    dlaset_("A", &c__50, &n, &c_b15, &c_b15, rfac, &c__50, (ftnlen)1);

    /*     %-------------------------------------% */
    /*     | KU, KL are number of superdiagonals | */
    /*     | and subdiagonals within the band of | */
    /*     | matrices A and M.                   | */
    /*     %-------------------------------------% */

    kl = 1;
    ku = 1;

    /*     %---------------% */
    /*     | Main diagonal | */
    /*     %---------------% */

    h = 1. / (double)(n + 1);
    r1 = .66666666666666663;
    idiag = kl + ku + 1;
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
        a[idiag + j * 50 - 51] = 2. / h;
        m[idiag + j * 50 - 51] = r1 * h;
        /* L30: */
    }

    /*     %-------------------------------------% */
    /*     | First subdiagonal and superdiagonal | */
    /*     %-------------------------------------% */

    r2 = .16666666666666666;
    isup = kl + ku;
    isub = kl + ku + 2;
    i__1 = n - 1;
    for (j = 1; j <= i__1; ++j)
    {
        a[isup + (j + 1) * 50 - 51] = -1. / h;
        a[isub + j * 50 - 51] = -1. / h;
        m[isup + (j + 1) * 50 - 51] = r2 * h;
        m[isub + j * 50 - 51] = r2 * h;
        /* L60: */
    }

    /*     %-------------------------------------% */
    /*     | Call DSBAND  to find eigenvalues and | */
    /*     | eigenvectors.  Eigenvalues are      | */
    /*     | returned in the first column of D.  | */
    /*     | Eigenvectors are returned in the    | */
    /*     | first NCONV (=IPARAM(5)) columns of | */
    /*     | V.                                  | */
    /*     %-------------------------------------% */

    rvec = TRUE_;
    dsband_(&rvec, "A", select, d, v, &c__1000, &sigma, &n, a, m, &c__50, rfac, &kl, &ku, which, bmat, &nev, &tol, resid, &ncv, v, &c__1000, iparam, workd, workl, &lworkl, iwork, &info, (ftnlen)1, (ftnlen)2, (ftnlen)1);

    if (info == 0)
    {

        nconv = iparam[4];

        /*        %-----------------------------------% */
        /*        | Print out convergence information | */
        /*        %-----------------------------------% */

        s_wsle(&io___38);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
        s_wsle(&io___39);
        do_lio(&c__9, &c__1, " _SBDR5 ", (ftnlen)8);
        e_wsle();
        s_wsle(&io___40);
        do_lio(&c__9, &c__1, " ====== ", (ftnlen)8);
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
        do_lio(&c__9, &c__1, " The number of Lanczos vectors generated", (ftnlen)40);
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
        do_lio(&c__9, &c__1, " The number of Implicit Arnoldi", (ftnlen)31);
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

        /*        %----------------------------% */
        /*        | Compute the residual norm. | */
        /*        |    ||  A*x - lambda*x ||   | */
        /*        %----------------------------% */

        i__1 = nconv;
        for (j = 1; j <= i__1; ++j)
        {
            dgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b97, &a[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b15, ax, &c__1, (ftnlen)11);
            dgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b97, &m[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b15, mx, &c__1, (ftnlen)11);
            d__1 = -d[j - 1];
            daxpy_(&n, &d__1, mx, &c__1, ax, &c__1);
            d[j + 49] = dnrm2_(&n, ax, &c__1);
            d[j + 49] /= (d__1 = d[j - 1], abs(d__1));

            /* L90: */
        }
        dmout_(&c__6, &nconv, &c__2, d, &c__50, &c_n6,
               "Ritz values and re"
               "lative residuals",
               (ftnlen)34);
    }
    else
    {

        /*        %-------------------------------------% */
        /*        | Either convergence failed, or there | */
        /*        | is error.  Check the documentation  | */
        /*        | for DSBAND .                         | */
        /*        %-------------------------------------% */

        s_wsle(&io___53);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
        s_wsle(&io___54);
        do_lio(&c__9, &c__1, " Error with _sband, info= ", (ftnlen)26);
        do_lio(&c__3, &c__1, (char *)&info, (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___55);
        do_lio(&c__9, &c__1, " Check the documentation of _sband ", (ftnlen)35);
        e_wsle();
        s_wsle(&io___56);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
    }

L9000:
    return 0;
} /* MAIN__ */

/* Main program alias */ int dsbdr5_()
{
    MAIN__();
    return 0;
}
