/* D:\projects\Fortran\arpack-ng-3.9.1-patched\EXAMPLES\BAND\dnbdr6.f -- translated by f2c (version 20230428).
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
static doublereal c_b15 = 0.;
static integer c__1000 = 1000;
static integer c__3 = 3;
static integer c__5 = 5;
static doublereal c_b101 = 1.;
static integer c__6 = 6;
static integer c_n6 = -6;

/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), e_wsle(void);
    int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    doublereal a[50000] /* was [50][1000] */, d__[150] /* was [50][3] */, h__;
    integer i__, j;
    doublereal m[50000] /* was [50][1000] */;
    integer n;
    doublereal v[50000] /* was [1000][50] */;
    integer kl;
    doublereal ax[1000];
    integer lo, ku;
    doublereal mx[1000];
    integer nx, ido, ncv, nev;
    doublereal rho, tol;
    doublecomplex cfac[50000] /* was [50][1000] */;
    doublereal rfac[50000] /* was [50][1000] */;
    char bmat[1];
    integer mode, info;
    logical rvec;
    integer isub, isup;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    integer idiag;
    extern int dgbmv_(char *, integer *, integer *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, ftnlen);
    char which[2];
    doublereal resid[1000];
    integer nconv;
    doublecomplex workc[1000];
    doublereal workd[3000];
    logical first;
    extern int daxpy_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *), dmout_(integer *, integer *, integer *, doublereal *, integer *, integer *, char *, ftnlen);
    integer iwork[1000];
    doublereal workl[7800];
    extern doublereal dlapy2_(doublereal *, doublereal *);
    extern int dnband_(logical *, char *, logical *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, doublereal *, integer *, doublereal *, doublecomplex *, integer *, integer *, char *, char *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, integer *, doublereal *, doublereal *, integer *, doublecomplex *, integer *, integer *, ftnlen, ftnlen, ftnlen);
    integer iparam[11];
    doublereal sigmai;
    extern int dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    logical select[50];
    doublereal sigmar;
    integer maxitr, lworkl;
    doublereal workev[150];

    /* Fortran I/O blocks */
    static cilist io___5 = {0, 6, 0, 0, 0};
    static cilist io___6 = {0, 6, 0, 0, 0};
    static cilist io___7 = {0, 6, 0, 0, 0};
    static cilist io___44 = {0, 6, 0, 0, 0};
    static cilist io___45 = {0, 6, 0, 0, 0};
    static cilist io___46 = {0, 6, 0, 0, 0};
    static cilist io___47 = {0, 6, 0, 0, 0};
    static cilist io___48 = {0, 6, 0, 0, 0};
    static cilist io___49 = {0, 6, 0, 0, 0};
    static cilist io___50 = {0, 6, 0, 0, 0};
    static cilist io___51 = {0, 6, 0, 0, 0};
    static cilist io___52 = {0, 6, 0, 0, 0};
    static cilist io___53 = {0, 6, 0, 0, 0};
    static cilist io___54 = {0, 6, 0, 0, 0};
    static cilist io___55 = {0, 6, 0, 0, 0};
    static cilist io___56 = {0, 6, 0, 0, 0};
    static cilist io___60 = {0, 6, 0, 0, 0};
    static cilist io___61 = {0, 6, 0, 0, 0};
    static cilist io___62 = {0, 6, 0, 0, 0};
    static cilist io___63 = {0, 6, 0, 0, 0};

    /*     ... Construct matrices A and M in LAPACK-style band form. */
    /*         The matrix A is a block tridiagonal matrix.  Each */
    /*         diagonal block is a tridiagonal matrix with */
    /*         4 on the diagonal, 1-rho*h/2 on the subdiagonal and */
    /*         1+rho*h/2 on the superdiagonal.  Each subdiagonal block */
    /*         of A is an identity matrix.  The matrix M is the */
    /*         tridiagonal matrix with 4 on the diagonal and 1 on the */
    /*         subdiagonal and superdiagonal. */

    /*     ... Define COMPLEX shift SIGMA=(SIGMAR,SIGMAI), SIGMAI .ne. zero. */

    /*     ... Call dnband  to find eigenvalues LAMBDA closest to SIGMA */
    /*         such that */
    /*                 A*x = LAMBDA*M*x. */

    /*     ... Use mode 4 of DNAUPD . */

    /* \BeginLib */

    /* \Routines called: */
    /*     dnband   ARPACK banded eigenproblem solver. */
    /*     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully. */
    /*     dlaset   LAPACK routine to initialize a matrix to zero. */
    /*     daxpy    Level 1 BLAS that computes y <- alpha*x+y. */
    /*     dnrm2    Level 1 BLAS that computes the norm of a vector. */
    /*     dgbmv    Level 2 BLAS that computes the band matrix vector product. */

    /* \Author */
    /*     Danny Sorensen */
    /*     Richard Lehoucq */
    /*     Chao Yang */
    /*     Dept. of Computational & */
    /*     Applied Mathematics */
    /*     Rice University */
    /*     Houston, Texas */

    /* \SCCS Information: @(#) */
    /* FILE: nbdr6.F   SID: 2.5   DATE OF SID: 08/26/96   RELEASE: 2 */

    /* \Remarks */
    /*     1. None */

    /* \EndLib */

    /* --------------------------------------------------------------------- */

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

    /*     %--------------------% */
    /*     | Intrinsic function | */
    /*     %--------------------% */

    /*     %-----------------------------% */
    /*     | BLAS & LAPACK routines used | */
    /*     %-----------------------------% */

    /*     %-----------------------% */
    /*     | Executable Statements | */
    /*     %-----------------------% */

    /*     %-----------------------------------------------------% */
    /*     | The number NX is the size of each diagonal block of | */
    /*     | A.  The number N(=NX*NX) is the dimension of the    | */
    /*     | matrix.  The number N(=NX*NX) is the dimension of   | */
    /*     | the matrix.  A generalized eigenvalue problem is    | */
    /*     | solved (BMAT = 'G').  NEV numbers of eigenvalues    | */
    /*     | closest to the COMPLEX shift (SIGMAR,SIGMAI)        | */
    /*     | (WHICH='LM') and their corresponding eigenvectors   | */
    /*     | are computed. The user can modify NX, NEV, NCV,     | */
    /*     | WHICH to solve problems of different sizes, and     | */
    /*     | to get different parts the spectrum. However, the   | */
    /*     | following rules must be satisfied:                  | */
    /*     |                   N <= MAXN                         | */
    /*     |                 NEV <= MAXNEV                       | */
    /*     |           NEV + 2 <= NCV <= MAXNCV                  | */
    /*     %-----------------------------------------------------% */

    nx = 10;
    n = nx * nx;
    nev = 4;
    ncv = 10;
    if (n > 1000)
    {
        s_wsle(&io___5);
        do_lio(&c__9, &c__1, " ERROR with _NBDR6: N is greater than MAXN ", (ftnlen)43);
        e_wsle();
        goto L9000;
    }
    else if (nev > 25)
    {
        s_wsle(&io___6);
        do_lio(&c__9, &c__1, " ERROR with _NBDR6: NEV is greater than MAXNEV ", (ftnlen)47);
        e_wsle();
        goto L9000;
    }
    else if (ncv > 50)
    {
        s_wsle(&io___7);
        do_lio(&c__9, &c__1, " ERROR with _NBDR6: NCV is greater than MAXNCV ", (ftnlen)47);
        e_wsle();
        goto L9000;
    }
    *(unsigned char *)bmat = 'G';
    s_copy(which, "LM", (ftnlen)2, (ftnlen)2);
    sigmar = .4;
    sigmai = .6;

    /*     %-----------------------------------------------------% */
    /*     | The work array WORKL is used in DNAUPD  as           | */
    /*     | workspace.  Its dimension LWORKL is set as          | */
    /*     | illustrated below.  The parameter TOL determines    | */
    /*     | the stopping criterion. If TOL<=0, machine          | */
    /*     | precision is used.  The variable IDO is used for    | */
    /*     | reverse communication, and is initially set to 0.   | */
    /*     | Setting INFO=0 indicates that a random vector is    | */
    /*     | generated in DNAUPD  to start the Arnoldi iteration. | */
    /*     %-----------------------------------------------------% */

    /* Computing 2nd power */
    i__1 = ncv;
    lworkl = i__1 * i__1 * 3 + ncv * 6;
    tol = 0.;
    ido = 0;
    info = 0;

    /*     %---------------------------------------------------% */
    /*     | IPARAM(3) specifies the maximum number of Arnoldi | */
    /*     | iterations allowed.  Mode 4 of DNAUPD  is used     | */
    /*     | (IPARAm(7) = 4). All these options can be changed | */
    /*     | by the user. For details, see the documentation   | */
    /*     | in dnband .                                        | */
    /*     %---------------------------------------------------% */

    maxitr = 300;
    mode = 4;

    iparam[2] = maxitr;
    iparam[6] = mode;

    /*     %--------------------------------------------% */
    /*     | Construct matrices A and M in LAPACK-style | */
    /*     | banded form.                               | */
    /*     %--------------------------------------------% */

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

    kl = nx;
    ku = nx;

    /*     %---------------% */
    /*     | Main diagonal | */
    /*     %---------------% */

    idiag = kl + ku + 1;
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
        a[idiag + j * 50 - 51] = 4.;
        m[idiag + j * 50 - 51] = 4.;
        /* L30: */
    }

    /*     %-------------------------------------% */
    /*     | First subdiagonal and superdiagonal | */
    /*     %-------------------------------------% */

    isup = kl + ku;
    isub = kl + ku + 2;
    h__ = 1. / (doublereal)(nx + 1);
    rho = 100.;
    i__1 = nx;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
        lo = (i__ - 1) * nx;
        i__2 = lo + nx - 1;
        for (j = lo + 1; j <= i__2; ++j)
        {
            a[isub + (j + 1) * 50 - 51] = h__ * rho / 2. - 1.;
            a[isup + j * 50 - 51] = -1. - h__ * rho / 2.;
            /* L40: */
        }
        /* L50: */
    }

    i__1 = n - 1;
    for (j = 1; j <= i__1; ++j)
    {
        m[isub + (j + 1) * 50 - 51] = 1.;
        m[isup + j * 50 - 51] = 1.;
        /* L60: */
    }

    /*     %------------------------------------% */
    /*     | KL-th subdiagonal and KU-th super- | */
    /*     | diagonal.                          | */
    /*     %------------------------------------% */

    isup = kl + 1;
    isub = (kl << 1) + ku + 1;
    i__1 = nx - 1;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
        lo = (i__ - 1) * nx;
        i__2 = lo + nx;
        for (j = lo + 1; j <= i__2; ++j)
        {
            a[isup + (nx + j) * 50 - 51] = -1.;
            a[isub + j * 50 - 51] = -1.;
            /* L70: */
        }
        /* L80: */
    }

    /*     %------------------------------------------------% */
    /*     | Call ARPACK banded solver to find eigenvalues  | */
    /*     | and eigenvectors. The real parts of the        | */
    /*     | eigenvalues are returned in the first column   | */
    /*     | of D, the imaginary parts are returned in the  | */
    /*     | second column of D.  Eigenvectors are returned | */
    /*     | in the first NCONV (=IPARAM(5)) columns of V.  | */
    /*     %------------------------------------------------% */

    rvec = TRUE_;
    dnband_(&rvec, "A", select, d__, &d__[50], v, &c__1000, &sigmar, &sigmai, workev, &n, a, m, &c__50, rfac, cfac, &ku, &kl, which, bmat, &nev, &tol, resid, &ncv, v, &c__1000, iparam, workd, workl, &lworkl, workc, iwork, &info, (ftnlen)1, (ftnlen)2, (ftnlen)1);

    if (info == 0)
    {

        /*        %-----------------------------------% */
        /*        | Print out convergence information | */
        /*        %-----------------------------------% */

        nconv = iparam[4];

        s_wsle(&io___44);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
        s_wsle(&io___45);
        do_lio(&c__9, &c__1, " _NBDR6 ", (ftnlen)8);
        e_wsle();
        s_wsle(&io___46);
        do_lio(&c__9, &c__1, " ====== ", (ftnlen)8);
        e_wsle();
        s_wsle(&io___47);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
        s_wsle(&io___48);
        do_lio(&c__9, &c__1, " The size of the matrix is ", (ftnlen)27);
        do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___49);
        do_lio(&c__9, &c__1, " Number of eigenvalue requested is ", (ftnlen)35);
        do_lio(&c__3, &c__1, (char *)&nev, (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___50);
        do_lio(&c__9, &c__1, " The number of Arnoldi vectors generated", (ftnlen)40);
        do_lio(&c__9, &c__1, " (NCV) is ", (ftnlen)10);
        do_lio(&c__3, &c__1, (char *)&ncv, (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___51);
        do_lio(&c__9, &c__1, " The number of converged Ritz values is ", (ftnlen)40);
        do_lio(&c__3, &c__1, (char *)&nconv, (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___52);
        do_lio(&c__9, &c__1, " What portion of the spectrum ", (ftnlen)30);
        do_lio(&c__9, &c__1, which, (ftnlen)2);
        e_wsle();
        s_wsle(&io___53);
        do_lio(&c__9, &c__1, " The number of Implicit Arnoldi ", (ftnlen)32);
        do_lio(&c__9, &c__1, " update taken is ", (ftnlen)17);
        do_lio(&c__3, &c__1, (char *)&iparam[2], (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___54);
        do_lio(&c__9, &c__1, " The number of OP*x is ", (ftnlen)23);
        do_lio(&c__3, &c__1, (char *)&iparam[8], (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___55);
        do_lio(&c__9, &c__1, " The convergence tolerance is ", (ftnlen)30);
        do_lio(&c__5, &c__1, (char *)&tol, (ftnlen)sizeof(doublereal));
        e_wsle();
        s_wsle(&io___56);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();

        /*        %----------------------------% */
        /*        | Compute the residual norm. | */
        /*        |    ||  A*x - lambda*x ||   | */
        /*        %----------------------------% */

        first = TRUE_;
        i__1 = nconv;
        for (j = 1; j <= i__1; ++j)
        {

            if (d__[j + 49] == 0.)
            {

                /*              %--------------------% */
                /*              | Ritz value is real | */
                /*              %--------------------% */

                dgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b101, &a[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b15, ax, &c__1, (ftnlen)11);
                dgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b101, &m[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b15, mx, &c__1, (ftnlen)11);
                d__1 = -d__[j - 1];
                daxpy_(&n, &d__1, mx, &c__1, ax, &c__1);
                d__[j + 99] = dnrm2_(&n, ax, &c__1);
                d__[j + 99] /= (d__1 = d__[j - 1], abs(d__1));
            }
            else if (first)
            {

                /*              %------------------------% */
                /*              | Ritz value is complex  | */
                /*              | Residual of one Ritz   | */
                /*              | value of the conjugate | */
                /*              | pair is computed.      | */
                /*              %------------------------% */

                dgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b101, &a[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b15, ax, &c__1, (ftnlen)11);
                dgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b101, &m[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b15, mx, &c__1, (ftnlen)11);
                d__1 = -d__[j - 1];
                daxpy_(&n, &d__1, mx, &c__1, ax, &c__1);
                dgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b101, &m[kl], &c__50, &v[(j + 1) * 1000 - 1000], &c__1, &c_b15, mx, &c__1, (ftnlen)11);
                daxpy_(&n, &d__[j + 49], mx, &c__1, ax, &c__1);
                d__[j + 99] = dnrm2_(&n, ax, &c__1);
                dgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b101, &a[kl], &c__50, &v[(j + 1) * 1000 - 1000], &c__1, &c_b15, ax, &c__1, (ftnlen)11);
                dgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b101, &m[kl], &c__50, &v[(j + 1) * 1000 - 1000], &c__1, &c_b15, mx, &c__1, (ftnlen)11);
                d__1 = -d__[j - 1];
                daxpy_(&n, &d__1, mx, &c__1, ax, &c__1);
                dgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b101, &m[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b15, mx, &c__1, (ftnlen)11);
                d__1 = -d__[j + 49];
                daxpy_(&n, &d__1, mx, &c__1, ax, &c__1);
                d__1 = dnrm2_(&n, ax, &c__1);
                d__[j + 99] = dlapy2_(&d__[j + 99], &d__1);
                d__[j + 99] /= dlapy2_(&d__[j - 1], &d__[j + 49]);
                d__[j + 100] = d__[j + 99];
                first = FALSE_;
            }
            else
            {
                first = TRUE_;
            }

            /* L90: */
        }
        dmout_(&c__6, &nconv, &c__3, d__, &c__50, &c_n6,
               "Ritz values (Real,"
               "Imag) and relative residuals",
               (ftnlen)46);
    }
    else
    {

        /*        %-------------------------------------% */
        /*        | Either convergence failed, or there | */
        /*        | is error.  Check the documentation  | */
        /*        | for dnband .                         | */
        /*        %-------------------------------------% */

        s_wsle(&io___60);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
        s_wsle(&io___61);
        do_lio(&c__9, &c__1, " Error with _nband, info= ", (ftnlen)26);
        do_lio(&c__3, &c__1, (char *)&info, (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___62);
        do_lio(&c__9, &c__1, " Check the documentation of _nband ", (ftnlen)35);
        e_wsle();
        s_wsle(&io___63);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
    }

L9000:
    return 0;
} /* MAIN__ */

/* Main program alias */ int dnbdr6_()
{
    MAIN__();
    return 0;
}
