/* D:\projects\Fortran\arpack-ng-3.9.1-patched\EXAMPLES\BAND\ssbdr5.f -- translated by f2c (version 20230428).
   You must link the resulting object file with libf2c:
    on Microsoft Windows system, link with libf2c.lib;
    on Linux or Unix systems, link with .../path/to/libf2c.a -lm
    or, if you install libf2c.a in a standard place, with -lf2c -lm
    -- in that order, at the end of the command line, as in
        cc *.o -lf2c -lm
    Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

        http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__50 = 50;
static real c_b15 = 0.f;
static integer c__1000 = 1000;
static integer c__3 = 3;
static integer c__4 = 4;
static real c_b97 = 1.f;
static integer c__6 = 6;
static integer c__2 = 2;
static integer c_n6 = -6;

/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), e_wsle(void);
    int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    real a[50000] /* was [50][1000] */, d__[100] /* was [50][2] */, h__;
    integer j;
    real m[50000] /* was [50][1000] */;
    integer n;
    real v[50000] /* was [1000][50] */, r1, r2;
    integer kl;
    real ax[1000];
    integer ku;
    real mx[1000];
    integer ido, ncv, nev;
    real tol, rfac[50000] /* was [50][1000] */;
    char bmat[1];
    integer mode, info;
    logical rvec;
    integer isub, isup;
    extern doublereal snrm2_(integer *, real *, integer *);
    integer idiag;
    real sigma;
    char which[2];
    real resid[1000];
    extern int sgbmv_(char *, integer *, integer *, integer *, integer *, real *, real *, integer *, real *, integer *, real *, real *, integer *, ftnlen);
    integer nconv;
    real workd[3000];
    integer iwork[1000];
    real workl[2900];
    extern int saxpy_(integer *, real *, real *, integer *, real *, integer *), smout_(integer *, integer *, integer *, real *, integer *, integer *, char *, ftnlen);
    integer iparam[11];
    extern int ssband_(logical *, char *, logical *, real *, real *, integer *, real *, integer *, real *, real *, integer *, real *, integer *, integer *, char *, char *, integer *, real *, real *, integer *, real *, integer *, integer *, real *, real *, integer *, integer *, integer *, ftnlen, ftnlen, ftnlen);
    logical select[50];
    extern int slaset_(char *, integer *, integer *, real *, real *, real *, integer *, ftnlen);
    integer maxitr, lworkl;

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

    /*     ... Call SSBAND with Buckling mode to find eigenvalues LAMBDA */
    /*         such that */
    /*                          A*x = M*x*LAMBDA. */

    /*     ... Use mode 4 of SSAUPD. */

    /* \BeginLib */

    /* \Routines called: */
    /*     ssband  ARPACK banded eigenproblem solver. */
    /*     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully. */
    /*     slaset  LAPACK routine to initialize a matrix to zero. */
    /*     saxpy   Level 1 BLAS that computes y <- alpha*x+y. */
    /*     snrm2   Level 1 BLAS that computes the norm of a vector. */
    /*     sgbmv   Level 2 BLAS that computes the band matrix vector product */

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
    /*     | The work array WORKL is used in SSAUPD as           | */
    /*     | workspace.  Its dimension LWORKL is set as          | */
    /*     | illustrated below.  The parameter TOL determines    | */
    /*     | the stopping criterion. If TOL<=0, machine          | */
    /*     | precision is used.  The variable IDO is used for    | */
    /*     | reverse communication, and is initially set to 0.   | */
    /*     | Setting INFO=0 indicates that a random vector is    | */
    /*     | generated in SSAUPD to start the Arnoldi iteration. | */
    /*     %-----------------------------------------------------% */

    /* Computing 2nd power */
    i__1 = ncv;
    lworkl = i__1 * i__1 + (ncv << 3);
    tol = 0.f;
    ido = 0;
    info = 0;

    /*     %---------------------------------------------------% */
    /*     | IPARAM(3) specifies the maximum number of Arnoldi | */
    /*     | iterations allowed.  Mode 4 of SSAUPD is used     | */
    /*     | (IPARAM(7) = 4). All these options can be changed | */
    /*     | by the user. For details see the documentation in | */
    /*     | SSBAND.                                           | */
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

    slaset_("A", &c__50, &n, &c_b15, &c_b15, a, &c__50, (ftnlen)1);
    slaset_("A", &c__50, &n, &c_b15, &c_b15, m, &c__50, (ftnlen)1);
    slaset_("A", &c__50, &n, &c_b15, &c_b15, rfac, &c__50, (ftnlen)1);

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

    h__ = 1.f / (real)(n + 1);
    r1 = .66666666666666663f;
    idiag = kl + ku + 1;
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
        a[idiag + j * 50 - 51] = 2.f / h__;
        m[idiag + j * 50 - 51] = r1 * h__;
        /* L30: */
    }

    /*     %-------------------------------------% */
    /*     | First subdiagonal and superdiagonal | */
    /*     %-------------------------------------% */

    r2 = .16666666666666666f;
    isup = kl + ku;
    isub = kl + ku + 2;
    i__1 = n - 1;
    for (j = 1; j <= i__1; ++j)
    {
        a[isup + (j + 1) * 50 - 51] = -1.f / h__;
        a[isub + j * 50 - 51] = -1.f / h__;
        m[isup + (j + 1) * 50 - 51] = r2 * h__;
        m[isub + j * 50 - 51] = r2 * h__;
        /* L60: */
    }

    /*     %-------------------------------------% */
    /*     | Call SSBAND to find eigenvalues and | */
    /*     | eigenvectors.  Eigenvalues are      | */
    /*     | returned in the first column of D.  | */
    /*     | Eigenvectors are returned in the    | */
    /*     | first NCONV (=IPARAM(5)) columns of | */
    /*     | V.                                  | */
    /*     %-------------------------------------% */

    rvec = TRUE_;
    ssband_(&rvec, "A", select, d__, v, &c__1000, &sigma, &n, a, m, &c__50, rfac, &kl, &ku, which, bmat, &nev, &tol, resid, &ncv, v, &c__1000, iparam, workd, workl, &lworkl, iwork, &info, (ftnlen)1, (ftnlen)2, (ftnlen)1);

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
        do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___43);
        do_lio(&c__9, &c__1, " Number of eigenvalue requested is ", (ftnlen)35);
        do_lio(&c__3, &c__1, (char *)&nev, (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___44);
        do_lio(&c__9, &c__1, " The number of Lanczos vectors generated", (ftnlen)40);
        do_lio(&c__9, &c__1, " (NCV) is ", (ftnlen)10);
        do_lio(&c__3, &c__1, (char *)&ncv, (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___45);
        do_lio(&c__9, &c__1, " The number of converged Ritz values is ", (ftnlen)40);
        do_lio(&c__3, &c__1, (char *)&nconv, (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___46);
        do_lio(&c__9, &c__1, " What portion of the spectrum ", (ftnlen)30);
        do_lio(&c__9, &c__1, which, (ftnlen)2);
        e_wsle();
        s_wsle(&io___47);
        do_lio(&c__9, &c__1, " The number of Implicit Arnoldi", (ftnlen)31);
        do_lio(&c__9, &c__1, " update taken is ", (ftnlen)17);
        do_lio(&c__3, &c__1, (char *)&iparam[2], (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___48);
        do_lio(&c__9, &c__1, " The number of OP*x is ", (ftnlen)23);
        do_lio(&c__3, &c__1, (char *)&iparam[8], (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___49);
        do_lio(&c__9, &c__1, " The convergence tolerance is ", (ftnlen)30);
        do_lio(&c__4, &c__1, (char *)&tol, (ftnlen)sizeof(real));
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
            sgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b97, &a[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b15, ax, &c__1, (ftnlen)11);
            sgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b97, &m[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b15, mx, &c__1, (ftnlen)11);
            r__1 = -d__[j - 1];
            saxpy_(&n, &r__1, mx, &c__1, ax, &c__1);
            d__[j + 49] = snrm2_(&n, ax, &c__1);
            d__[j + 49] /= (r__1 = d__[j - 1], dabs(r__1));

            /* L90: */
        }
        smout_(&c__6, &nconv, &c__2, d__, &c__50, &c_n6,
               "Ritz values and re"
               "lative residuals",
               (ftnlen)34);
    }
    else
    {

        /*        %-------------------------------------% */
        /*        | Either convergence failed, or there | */
        /*        | is error.  Check the documentation  | */
        /*        | for SSBAND.                         | */
        /*        %-------------------------------------% */

        s_wsle(&io___53);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
        s_wsle(&io___54);
        do_lio(&c__9, &c__1, " Error with _sband, info= ", (ftnlen)26);
        do_lio(&c__3, &c__1, (char *)&info, (ftnlen)sizeof(integer));
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

/* Main program alias */ int ssbdr5_()
{
    MAIN__();
    return 0;
}
