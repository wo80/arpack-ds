/* EXAMPLES\BAND\snbdr2.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

/* Table of constant values */

static a_int c__9 = 9;
static a_int c__1 = 1;
static a_int c__50 = 50;
static float c_b15 = 0.f;
static a_int c__1000 = 1000;
static a_int c__3 = 3;
static a_int c__4 = 4;
static float c_b100 = 1.f;
static a_int c__6 = 6;
static a_int c_n6 = -6;

/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    a_int i__1, i__2;
    float r__1;

    /* Builtin functions */
    a_int s_wsle(cilist *), do_lio(a_int *, a_int *, char *, ftnlen), e_wsle(void);
    int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    float a[50000] /* was [50][1000] */, d[150] /* was [50][3] */, h;
    a_int i, j;
    float m[50000] /* was [50][1000] */;
    a_int n;
    float v[50000] /* was [1000][50] */, h2;
    a_int kl;
    float ax[1000];
    a_int lo, ku, nx, ido, ncv, nev;
    float rho, tol;
    a_fcomplex cfac[50000] /* was [50][1000] */;
    float rfac[50000] /* was [50][1000] */;
    char bmat[1];
    a_int mode, info;
    a_bool rvec;
    a_int isub, isup;
    extern double snrm2_(a_int *, float *, a_int *);
    a_int idiag;
    char which[2];
    float resid[1000];
    extern int sgbmv_(char *, a_int *, a_int *, a_int *, a_int *, float *, float *, a_int *, float *, a_int *, float *, float *, a_int *, ftnlen);
    a_int nconv;
    a_fcomplex workc[1000];
    float workd[3000];
    a_bool first;
    a_int iwork[1000];
    float workl[7800];
    extern int saxpy_(a_int *, float *, float *, a_int *, float *, a_int *), smout_(a_int *, a_int *, a_int *, float *, a_int *, a_int *, char *, ftnlen);
    extern double slapy2_(float *, float *);
    extern int snband_(a_bool *, char *, a_bool *, float *, float *, float *, a_int *, float *, float *, float *, a_int *, float *, float *, a_int *, float *, a_fcomplex *, a_int *, a_int *, char *, char *, a_int *, float *, float *, a_int *, float *, a_int *, a_int *, float *, float *, a_int *, a_fcomplex *, a_int *, a_int *, ftnlen, ftnlen, ftnlen);
    a_int iparam[11];
    float sigmai;
    a_bool select[50];
    float sigmar;
    extern int slaset_(char *, a_int *, a_int *, float *, float *, float *, a_int *, ftnlen);
    a_int maxitr, lworkl;
    float workev[150];

    /* Fortran I/O blocks */
    static cilist io___5 = {0, 6, 0, 0, 0};
    static cilist io___6 = {0, 6, 0, 0, 0};
    static cilist io___7 = {0, 6, 0, 0, 0};
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
    static cilist io___57 = {0, 6, 0, 0, 0};
    static cilist io___60 = {0, 6, 0, 0, 0};
    static cilist io___61 = {0, 6, 0, 0, 0};
    static cilist io___62 = {0, 6, 0, 0, 0};
    static cilist io___63 = {0, 6, 0, 0, 0};

    /*     ... Construct matrices A in LAPACK-style band form. */
    /*         The matrix A is derived from the discretization of */
    /*         the 2-d convection-diffusion operator */

    /*               -Laplacian(u) + rho*partial(u)/partial(x). */

    /*         on the unit square with zero Dirichlet boundary condition */
    /*         using standard central difference. */

    /*     ... Define the shift SIGMA = (SIGMAR, SIGMAI). */

    /*     ... Call SNBAND to find eigenvalues LAMBDA closest to SIGMA */
    /*         such that */
    /*                       A*x = LAMBDA*x. */

    /*     ... Use mode 3 of SNAUPD. */

    /* \BeginLib */

    /* \Routines called: */
    /*     snband  ARPACK banded eigenproblem solver. */
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
    /* FILE: nbdr2.F   SID: 2.5   DATE OF SID: 08/26/96   RELEASE: 2 */

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

    /*     %-----------------------------% */
    /*     | BLAS & LAPACK routines used | */
    /*     %-----------------------------% */

    /*     %--------------------% */
    /*     | Intrinsic function | */
    /*     %--------------------% */

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
    /*     | eigenvalues (closest to (SIGMAR,SIGMAI)) to be  | */
    /*     | approximated. Since the shift-invert moded is   | */
    /*     | used, WHICH is set to 'LM'. The user can modify | */
    /*     | NX, NEV, NCV, SIGMAR, SIGMAI to solve problems  | */
    /*     | of different sizes, and to get different parts  | */
    /*     | the spectrum. However, The following conditions | */
    /*     | must be satisfied:                              | */
    /*     |                   N <= MAXN                     | */
    /*     |                 NEV <= MAXNEV                   | */
    /*     |           NEV + 2 <= NCV <= MAXNCV              | */
    /*     %-------------------------------------------------% */

    nx = 10;
    n = nx * nx;
    nev = 4;
    ncv = 20;
    if (n > 1000)
    {
        s_wsle(&io___5);
        do_lio(&c__9, &c__1, " ERROR with _NBDR2: N is greater than MAXN ", (ftnlen)43);
        e_wsle();
        goto L9000;
    }
    else if (nev > 25)
    {
        s_wsle(&io___6);
        do_lio(&c__9, &c__1, " ERROR with _NBDR2: NEV is greater than MAXNEV ", (ftnlen)47);
        e_wsle();
        goto L9000;
    }
    else if (ncv > 50)
    {
        s_wsle(&io___7);
        do_lio(&c__9, &c__1, " ERROR with _NBDR2: NCV is greater than MAXNCV ", (ftnlen)47);
        e_wsle();
        goto L9000;
    }
    *(unsigned char *)bmat = 'I';
    s_copy(which, "LM", (ftnlen)2, (ftnlen)2);
    sigmar = 1e4f;
    sigmai = 0.f;

    /*     %-----------------------------------------------------% */
    /*     | The work array WORKL is used in SNAUPD as           | */
    /*     | workspace.  Its dimension LWORKL is set as          | */
    /*     | illustrated below.  The parameter TOL determines    | */
    /*     | the stopping criterion. If TOL<=0, machine          | */
    /*     | precision is used.  The variable IDO is used for    | */
    /*     | reverse communication, and is initially set to 0.   | */
    /*     | Setting INFO=0 indicates that a random vector is    | */
    /*     | generated in SNAUPD to start the Arnoldi iteration. | */
    /*     %-----------------------------------------------------% */

    /* Computing 2nd power */
    i__1 = ncv;
    lworkl = i__1 * i__1 * 3 + ncv * 6;
    tol = 0.f;
    ido = 0;
    info = 0;

    /*     %---------------------------------------------------% */
    /*     | IPARAM(3) specifies the maximum number of Arnoldi | */
    /*     | iterations allowed.  Mode 3 of SNAUPD is used     | */
    /*     | (IPARAM(7) = 3). All these options can be changed | */
    /*     | by the user. For details, see the documentation   | */
    /*     | in SNBAND.                                        | */
    /*     %---------------------------------------------------% */

    maxitr = 300;
    mode = 3;

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
    /*     | matrices A.                   | */
    /*     %-------------------------------------% */

    kl = nx;
    ku = nx;

    /*     %---------------% */
    /*     | Main diagonal | */
    /*     %---------------% */

    h = 1.f / (float)(nx + 1);
    h2 = h * h;

    idiag = kl + ku + 1;
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
        a[idiag + j * 50 - 51] = 4.f / h2;
        /* L30: */
    }

    /*     %-------------------------------------% */
    /*     | First subdiagonal and superdiagonal | */
    /*     %-------------------------------------% */

    isup = kl + ku;
    isub = kl + ku + 2;
    rho = 10.f;
    i__1 = nx;
    for (i = 1; i <= i__1; ++i)
    {
        lo = (i - 1) * nx;
        i__2 = lo + nx - 1;
        for (j = lo + 1; j <= i__2; ++j)
        {
            a[isub + (j + 1) * 50 - 51] = -1.f / h2 + rho / 2.f / h;
            a[isup + j * 50 - 51] = -1.f / h2 - rho / 2.f / h;
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
    for (i = 1; i <= i__1; ++i)
    {
        lo = (i - 1) * nx;
        i__2 = lo + nx;
        for (j = lo + 1; j <= i__2; ++j)
        {
            a[isup + (nx + j) * 50 - 51] = -1.f / h2;
            a[isub + j * 50 - 51] = -1.f / h2;
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
    snband_(&rvec, "A", select, d, &d[50], v, &c__1000, &sigmar, &sigmai, workev, &n, a, m, &c__50, rfac, cfac, &kl, &ku, which, bmat, &nev, &tol, resid, &ncv, v, &c__1000, iparam, workd, workl, &lworkl, workc, iwork, &info, (ftnlen)1, (ftnlen)2, (ftnlen)1);

    if (info == 0)
    {

        /*        %-----------------------------------% */
        /*        | Print out convergence information | */
        /*        %-----------------------------------% */

        nconv = iparam[4];

        s_wsle(&io___45);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
        s_wsle(&io___46);
        do_lio(&c__9, &c__1, " _NBDR2 ", (ftnlen)8);
        e_wsle();
        s_wsle(&io___47);
        do_lio(&c__9, &c__1, " ====== ", (ftnlen)8);
        e_wsle();
        s_wsle(&io___48);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
        s_wsle(&io___49);
        do_lio(&c__9, &c__1, " The size of the matrix is ", (ftnlen)27);
        do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___50);
        do_lio(&c__9, &c__1, " Number of eigenvalue requested is ", (ftnlen)35);
        do_lio(&c__3, &c__1, (char *)&nev, (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___51);
        do_lio(&c__9, &c__1, " The number of Arnoldi vectors generated", (ftnlen)40);
        do_lio(&c__9, &c__1, " (NCV) is ", (ftnlen)10);
        do_lio(&c__3, &c__1, (char *)&ncv, (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___52);
        do_lio(&c__9, &c__1, " The number of converged Ritz values is ", (ftnlen)40);
        do_lio(&c__3, &c__1, (char *)&nconv, (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___53);
        do_lio(&c__9, &c__1, " What portion of the spectrum ", (ftnlen)30);
        do_lio(&c__9, &c__1, which, (ftnlen)2);
        e_wsle();
        s_wsle(&io___54);
        do_lio(&c__9, &c__1, " The number of Implicit Arnoldi ", (ftnlen)32);
        do_lio(&c__9, &c__1, " update taken is ", (ftnlen)17);
        do_lio(&c__3, &c__1, (char *)&iparam[2], (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___55);
        do_lio(&c__9, &c__1, " The number of OP*x is ", (ftnlen)23);
        do_lio(&c__3, &c__1, (char *)&iparam[8], (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___56);
        do_lio(&c__9, &c__1, " The convergence tolerance is ", (ftnlen)30);
        do_lio(&c__4, &c__1, (char *)&tol, (ftnlen)sizeof(float));
        e_wsle();
        s_wsle(&io___57);
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

            if (d[j + 49] == 0.f)
            {

                /*              %--------------------% */
                /*              | Ritz value is real | */
                /*              %--------------------% */

                sgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b100, &a[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b15, ax, &c__1, (ftnlen)11);
                r__1 = -d[j - 1];
                saxpy_(&n, &r__1, &v[j * 1000 - 1000], &c__1, ax, &c__1);
                d[j + 99] = snrm2_(&n, ax, &c__1);
                d[j + 99] /= (r__1 = d[j - 1], dabs(r__1));
            }
            else if (first)
            {

                /*              %------------------------% */
                /*              | Ritz value is complex  | */
                /*              | Residual of one Ritz   | */
                /*              | value of the conjugate | */
                /*              | pair is computed.      | */
                /*              %------------------------% */

                sgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b100, &a[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b15, ax, &c__1, (ftnlen)11);
                r__1 = -d[j - 1];
                saxpy_(&n, &r__1, &v[j * 1000 - 1000], &c__1, ax, &c__1);
                saxpy_(&n, &d[j + 49], &v[(j + 1) * 1000 - 1000], &c__1, ax, &c__1);
                d[j + 99] = snrm2_(&n, ax, &c__1);
                sgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b100, &a[kl], &c__50, &v[(j + 1) * 1000 - 1000], &c__1, &c_b15, ax, &c__1, (ftnlen)11);
                r__1 = -d[j - 1];
                saxpy_(&n, &r__1, &v[(j + 1) * 1000 - 1000], &c__1, ax, &c__1);
                r__1 = -d[j + 49];
                saxpy_(&n, &r__1, &v[j * 1000 - 1000], &c__1, ax, &c__1);
                r__1 = snrm2_(&n, ax, &c__1);
                d[j + 99] = slapy2_(&d[j + 99], &r__1);
                d[j + 99] /= slapy2_(&d[j - 1], &d[j + 49]);
                d[j + 100] = d[j + 99];
                first = FALSE_;
            }
            else
            {
                first = TRUE_;
            }

            /* L90: */
        }
        smout_(&c__6, &nconv, &c__3, d, &c__50, &c_n6,
               "Ritz values (Real,"
               "Imag) and relative residuals",
               (ftnlen)46);
    }
    else
    {

        /*        %-------------------------------------% */
        /*        | Either convergence failed, or there | */
        /*        | is error.  Check the documentation  | */
        /*        | for SNBAND.                         | */
        /*        %-------------------------------------% */

        s_wsle(&io___60);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
        s_wsle(&io___61);
        do_lio(&c__9, &c__1, " Error with _nband, info= ", (ftnlen)26);
        do_lio(&c__3, &c__1, (char *)&info, (ftnlen)sizeof(a_int));
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

/* Main program alias */ int snbdr2_()
{
    MAIN__();
    return 0;
}
