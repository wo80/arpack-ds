/* EXAMPLES\BAND\cnbdr1.f -- translated by f2c (version 20230428). */

#include "f2c.h"

/* Table of constant values */

static complex c_b1 = {1.f, 0.f};
static complex c_b2 = {0.f, 0.f};
static complex c_b3 = {2.f, 0.f};
static integer c__9 = 9;
static integer c__1 = 1;
static integer c__50 = 50;
static complex c_b26 = {4.f, 0.f};
static integer c__1000 = 1000;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c__6 = 6;
static integer c_n6 = -6;

/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    complex q__1, q__2, q__3, q__4, q__5;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), e_wsle(void);
    int s_copy(char *, char *, ftnlen, ftnlen);
    void c_div(complex *, complex *, complex *);
    double r_imag(complex *);

    /* Local variables */
    complex a[50000] /* was [50][1000] */, d__[50], h__;
    integer i__, j;
    complex m[50000] /* was [50][1000] */;
    integer n;
    complex v[50000] /* was [1000][50] */, h2;
    real rd[150] /* was [50][3] */;
    integer kl;
    complex ax[1000];
    integer lo, ku, nx;
    complex fac[50000] /* was [50][1000] */;
    integer ncv, nev;
    complex rho;
    real tol;
    char bmat[1];
    integer mode, info;
    logical rvec;
    integer isub, isup, idiag;
    extern int cgbmv_(char *, integer *, integer *, integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, complex *, integer *, ftnlen);
    complex sigma;
    char which[2];
    complex resid[1000];
    integer nconv;
    extern int caxpy_(integer *, complex *, complex *, integer *, complex *, integer *);
    complex workd[3000];
    integer iwork[1000];
    complex workl[7750];
    real rwork[1000];
    extern int smout_(integer *, integer *, integer *, real *, integer *, integer *, char *, ftnlen);
    extern doublereal scnrm2_(integer *, complex *, integer *), slapy2_(real *, real *);
    extern int cnband_(logical *, char *, logical *, complex *, complex *, integer *, complex *, complex *, integer *, complex *, complex *, integer *, complex *, integer *, integer *, char *, char *, integer *, real *, complex *, integer *, complex *, integer *, integer *, complex *, complex *, integer *, real *, integer *, integer *, ftnlen, ftnlen, ftnlen);
    integer iparam[11];
    extern int claset_(char *, integer *, integer *, complex *, complex *, complex *, integer *, ftnlen);
    logical select[50];
    integer maxitr, lworkl;
    complex workev[100];

    /* Fortran I/O blocks */
    static cilist io___5 = {0, 6, 0, 0, 0};
    static cilist io___6 = {0, 6, 0, 0, 0};
    static cilist io___7 = {0, 6, 0, 0, 0};
    static cilist io___42 = {0, 6, 0, 0, 0};
    static cilist io___43 = {0, 6, 0, 0, 0};
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
    static cilist io___57 = {0, 6, 0, 0, 0};
    static cilist io___58 = {0, 6, 0, 0, 0};
    static cilist io___59 = {0, 6, 0, 0, 0};
    static cilist io___60 = {0, 6, 0, 0, 0};

    /*     ... Construct the matrix A in LAPACK-style band form. */
    /*         The matrix A is derived from the discretization of */
    /*         the 2-d convection-diffusion operator */

    /*              -Laplacian(u) + rho*partial(u)/partial(x). */

    /*         on the unit square with zero Dirichlet boundary condition */
    /*         using standard central difference. */

    /*     ... Call CNBAND to find eigenvalues LAMBDA such that */
    /*                          A*x = x*LAMBDA. */

    /*     ... Use mode 1 of CNAUPD. */

    /* \BeginLib */

    /*     cnband  ARPACK banded eigenproblem solver. */
    /*     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully. */
    /*     claset  LAPACK routine to initialize a matrix to zero. */
    /*     caxpy   Level 1 BLAS that computes y <- alpha*x+y. */
    /*     scnrm2  Level 1 BLAS that computes the norm of a vector. */
    /*     cgbmv   Level 2 BLAS that computes the band matrix vector product */

    /* \Author */
    /*     Richard Lehoucq */
    /*     Danny Sorensen */
    /*     Chao Yang */
    /*     Dept. of Computational & */
    /*     Applied Mathematics */
    /*     Rice University */
    /*     Houston, Texas */

    /* \SCCS Information: @(#) */
    /* FILE: nbdr1.F   SID: 2.3   DATE OF SID: 08/26/96   RELEASE: 2 */

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

    /*     %-----------------------% */
    /*     | Executable Statements | */
    /*     %-----------------------% */

    /*     %-------------------------------------------------% */
    /*     | The number NX is the number of interior points  | */
    /*     | in the discretization of the 2-dimensional      | */
    /*     | convection-diffusion operator on the unit       | */
    /*     | square with zero Dirichlet boundary condition.  | */
    /*     | The number N(=NX*NX) is the dimension of the    | */
    /*     | matrix.  A standard eigenvalue problem is       | */
    /*     | solved (BMAT = 'I').  NEV is the number of      | */
    /*     | eigenvalues to be approximated. The user can    | */
    /*     | modify NX, NEV, NCV and WHICH to solve problems | */
    /*     | of different sizes, and to get different parts  | */
    /*     | the spectrum.  However, the following           | */
    /*     | conditions must be satisfied:                   | */
    /*     |                   N <= MAXN                     | */
    /*     |                 NEV <= MAXNEV                   | */
    /*     |           NEV + 2 <= NCV <= MAXNCV              | */
    /*     %-------------------------------------------------% */

    nx = 10;
    n = nx * nx;
    nev = 4;
    ncv = 10;
    if (n > 1000)
    {
        s_wsle(&io___5);
        do_lio(&c__9, &c__1, " ERROR with _NBDR1: N is greater than MAXN ", (ftnlen)43);
        e_wsle();
        goto L9000;
    }
    else if (nev > 25)
    {
        s_wsle(&io___6);
        do_lio(&c__9, &c__1, " ERROR with _NBDR1: NEV is greater than MAXNEV ", (ftnlen)47);
        e_wsle();
        goto L9000;
    }
    else if (ncv > 50)
    {
        s_wsle(&io___7);
        do_lio(&c__9, &c__1, " ERROR with _NBDR1: NCV is greater than MAXNCV ", (ftnlen)47);
        e_wsle();
        goto L9000;
    }
    *(unsigned char *)bmat = 'I';
    s_copy(which, "LM", (ftnlen)2, (ftnlen)2);

    /*     %-----------------------------------------------------% */
    /*     | The work array WORKL is used in CNAUPD as           | */
    /*     | workspace.  Its dimension LWORKL is set as          | */
    /*     | illustrated below.  The parameter TOL determines    | */
    /*     | the stopping criterion. If TOL<=0, machine          | */
    /*     | precision is used.  Setting INFO=0 indicates that a | */
    /*     | random vector is generated in CNAUPD to start the   | */
    /*     | Arnoldi iteration.                                  | */
    /*     %-----------------------------------------------------% */

    /* Computing 2nd power */
    i__1 = ncv;
    lworkl = i__1 * i__1 * 3 + ncv * 5;
    tol = 0.f;
    info = 0;

    /*     %---------------------------------------------------% */
    /*     | IPARAM(3) specifies the maximum number of Arnoldi | */
    /*     | iterations allowed.  Mode 1 of CNAUPD is used     | */
    /*     | (IPARAM(7) = 1). All these options can be changed | */
    /*     | by the user. For details, see the documentation   | */
    /*     | in cnband.                                        | */
    /*     %---------------------------------------------------% */

    maxitr = 300;
    mode = 1;

    iparam[2] = maxitr;
    iparam[6] = mode;

    /*     %----------------------------------------% */
    /*     | Construct the matrix A in LAPACK-style | */
    /*     | banded form.                           | */
    /*     %----------------------------------------% */

    /*     %---------------------------------------------% */
    /*     | Zero out the workspace for banded matrices. | */
    /*     %---------------------------------------------% */

    claset_("A", &c__50, &n, &c_b2, &c_b2, a, &c__50, (ftnlen)1);
    claset_("A", &c__50, &n, &c_b2, &c_b2, m, &c__50, (ftnlen)1);
    claset_("A", &c__50, &n, &c_b2, &c_b2, fac, &c__50, (ftnlen)1);

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

    i__1 = nx + 1;
    q__2.r = (real)i__1, q__2.i = 0.f;
    c_div(&q__1, &c_b1, &q__2);
    h__.r = q__1.r, h__.i = q__1.i;
    q__1.r = h__.r * h__.r - h__.i * h__.i, q__1.i = h__.r * h__.i + h__.i * h__.r;
    h2.r = q__1.r, h2.i = q__1.i;

    idiag = kl + ku + 1;
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
        i__2 = idiag + j * 50 - 51;
        c_div(&q__1, &c_b26, &h2);
        a[i__2].r = q__1.r, a[i__2].i = q__1.i;
        /* L30: */
    }

    /*     %-------------------------------------% */
    /*     | First subdiagonal and superdiagonal | */
    /*     %-------------------------------------% */

    rho.r = 100.f, rho.i = 0.f;
    isup = kl + ku;
    isub = kl + ku + 2;
    i__1 = nx;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
        lo = (i__ - 1) * nx;
        i__2 = lo + nx - 1;
        for (j = lo + 1; j <= i__2; ++j)
        {
            i__3 = isup + (j + 1) * 50 - 51;
            q__3.r = -1.f, q__3.i = -0.f;
            c_div(&q__2, &q__3, &h2);
            c_div(&q__5, &rho, &c_b3);
            c_div(&q__4, &q__5, &h__);
            q__1.r = q__2.r + q__4.r, q__1.i = q__2.i + q__4.i;
            a[i__3].r = q__1.r, a[i__3].i = q__1.i;
            i__3 = isub + j * 50 - 51;
            q__3.r = -1.f, q__3.i = -0.f;
            c_div(&q__2, &q__3, &h2);
            c_div(&q__5, &rho, &c_b3);
            c_div(&q__4, &q__5, &h__);
            q__1.r = q__2.r - q__4.r, q__1.i = q__2.i - q__4.i;
            a[i__3].r = q__1.r, a[i__3].i = q__1.i;
            /* L40: */
        }
        /* L50: */
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
            i__3 = isup + (nx + j) * 50 - 51;
            q__2.r = -1.f, q__2.i = -0.f;
            c_div(&q__1, &q__2, &h2);
            a[i__3].r = q__1.r, a[i__3].i = q__1.i;
            i__3 = isub + j * 50 - 51;
            q__2.r = -1.f, q__2.i = -0.f;
            c_div(&q__1, &q__2, &h2);
            a[i__3].r = q__1.r, a[i__3].i = q__1.i;
            /* L70: */
        }
        /* L80: */
    }

    /*     %-----------------------------------------------% */
    /*     | Call ARPACK banded solver to find eigenvalues | */
    /*     | and eigenvectors. Eigenvalues are returned in | */
    /*     | the one dimensional array D.  Eigenvectors    | */
    /*     | are returned in the first NCONV (=IPARAM(5))  | */
    /*     | columns of V.                                 | */
    /*     %-----------------------------------------------% */

    rvec = TRUE_;
    cnband_(&rvec, "A", select, d__, v, &c__1000, &sigma, workev, &n, a, m, &c__50, fac, &kl, &ku, which, bmat, &nev, &tol, resid, &ncv, v, &c__1000, iparam, workd, workl, &lworkl, rwork, iwork, &info, (ftnlen)1, (ftnlen)2, (ftnlen)1);

    if (info == 0)
    {

        nconv = iparam[4];

        /*        %-----------------------------------% */
        /*        | Print out convergence information | */
        /*        %-----------------------------------% */

        s_wsle(&io___42);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
        s_wsle(&io___43);
        do_lio(&c__9, &c__1, "_NBDR1 ", (ftnlen)7);
        e_wsle();
        s_wsle(&io___44);
        do_lio(&c__9, &c__1, "====== ", (ftnlen)7);
        e_wsle();
        s_wsle(&io___45);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
        s_wsle(&io___46);
        do_lio(&c__9, &c__1, " The size of the matrix is ", (ftnlen)27);
        do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___47);
        do_lio(&c__9, &c__1, " Number of eigenvalue requested is ", (ftnlen)35);
        do_lio(&c__3, &c__1, (char *)&nev, (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___48);
        do_lio(&c__9, &c__1, " The number of Arnoldi vectors generated", (ftnlen)40);
        do_lio(&c__9, &c__1, " (NCV) is ", (ftnlen)10);
        do_lio(&c__3, &c__1, (char *)&ncv, (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___49);
        do_lio(&c__9, &c__1, " The number of converged Ritz values is ", (ftnlen)40);
        do_lio(&c__3, &c__1, (char *)&nconv, (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___50);
        do_lio(&c__9, &c__1, " What portion of the spectrum ", (ftnlen)30);
        do_lio(&c__9, &c__1, which, (ftnlen)2);
        e_wsle();
        s_wsle(&io___51);
        do_lio(&c__9, &c__1, " The number of Implicit Arnoldi ", (ftnlen)32);
        do_lio(&c__9, &c__1, " update taken is ", (ftnlen)17);
        do_lio(&c__3, &c__1, (char *)&iparam[2], (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___52);
        do_lio(&c__9, &c__1, " The number of OP*x is ", (ftnlen)23);
        do_lio(&c__3, &c__1, (char *)&iparam[8], (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___53);
        do_lio(&c__9, &c__1, " The convergence tolerance is ", (ftnlen)30);
        do_lio(&c__4, &c__1, (char *)&tol, (ftnlen)sizeof(real));
        e_wsle();
        s_wsle(&io___54);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();

        /*        %----------------------------% */
        /*        | Compute the residual norm. | */
        /*        |    ||  A*x - lambda*x ||   | */
        /*        %----------------------------% */

        i__1 = nconv;
        for (j = 1; j <= i__1; ++j)
        {

            /*           %---------------------------% */
            /*           | Compute the residual norm | */
            /*           |   ||  A*x - lambda*x ||   | */
            /*           %---------------------------% */

            cgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b1, &a[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b2, ax, &c__1, (ftnlen)11);
            i__2 = j - 1;
            q__1.r = -d__[i__2].r, q__1.i = -d__[i__2].i;
            caxpy_(&n, &q__1, &v[j * 1000 - 1000], &c__1, ax, &c__1);
            i__2 = j - 1;
            rd[j - 1] = d__[i__2].r;
            rd[j + 49] = r_imag(&d__[j - 1]);
            rd[j + 99] = scnrm2_(&n, ax, &c__1);
            rd[j + 99] /= slapy2_(&rd[j - 1], &rd[j + 49]);
            /* L90: */
        }
        smout_(&c__6, &nconv, &c__3, rd, &c__50, &c_n6,
               "Ritz values (Real,I"
               "mag) and relative residuals",
               (ftnlen)46);
    }
    else
    {

        /*        %-------------------------------------% */
        /*        | Either convergence failed, or there | */
        /*        | is error.  Check the documentation  | */
        /*        | for cnband.                         | */
        /*        %-------------------------------------% */

        s_wsle(&io___57);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
        s_wsle(&io___58);
        do_lio(&c__9, &c__1, " Error with _nband, info= ", (ftnlen)26);
        do_lio(&c__3, &c__1, (char *)&info, (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___59);
        do_lio(&c__9, &c__1, " Check the documentation of _nband ", (ftnlen)35);
        e_wsle();
        s_wsle(&io___60);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
    }

L9000:
    return 0;
} /* MAIN__ */

/* Main program alias */ int cnbdr1_()
{
    MAIN__();
    return 0;
}
