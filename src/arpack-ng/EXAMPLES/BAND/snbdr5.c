/* EXAMPLES\BAND\snbdr5.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

/* Table of constant values */

static a_int c__9 = 9;
static a_int c__1 = 1;
static a_int c__50 = 50;
static float c_b15 = 0.f;
static a_int c__1000 = 1000;
static a_int c__3 = 3;
static a_int c__4 = 4;
static float c_b101 = 1.f;
static a_int c__6 = 6;
static a_int c_n6 = -6;

int main()
{
    /* System generated locals */
    a_int i__1, i__2;
    float r__1;

    /* Local variables */
    float a[50000] /* was [50][1000] */, d[150] /* was [50][3] */, h;
    a_int i, j;
    float m[50000] /* was [50][1000] */;
    a_int n;
    float v[50000] /* was [1000][50] */;
    a_int kl;
    float ax[1000];
    a_int lo, ku, nx, ido, ncv, nev;
    float rho, tol;
    a_fcomplex cfac[50000] /* was [50][1000] */;
    float rfac[50000] /* was [50][1000] */;
    char* bmat;
    a_int mode, info;
    a_bool rvec;
    a_int isub, isup;
    a_int idiag;
    char* which;
    float resid[1000];
    a_int nconv;
    a_fcomplex workc[1000];
    float workd[3000];
    a_bool first;
    a_int iwork[1000];
    float workl[7800];
    a_int iparam[11];
    float sigmai;
    a_bool select[50];
    float sigmar;
    a_int maxitr, lworkl;
    float workev[150];

    /*     ... Construct matrices A and M in LAPACK-style band form. */
    /*         The matrix A is a block tridiagonal matrix.  Each */
    /*         diagonal block is a tridiagonal matrix with */
    /*         4 on the diagonal, 1-rho*h/2 on the subdiagonal and */
    /*         1+rho*h/2 on the superdiagonal.  Each off-diagonal block */
    /*         of A is an identity matrices. */

    /*     ... Define COMPLEX shift SIGMA = (SIGMAR,SIGMAI), SIGMAI .ne. 0. */

    /*     ... Call SNBAND to find eigenvalues LAMBDA closest to SIGMA */
    /*         such that */
    /*                       A*x = LAMBDA*x. */

    /*     ... Use mode 4 of SNAUPD. */

    /* \BeginLib */

    /* \Routines called: */
    /*     snband  ARPACK banded eigenproblem solver. */
    /*     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully. */
    /*     slaset  LAPACK routine to initialize a matrix to zero. */
    /*     saxpy   Level 1 BLAS that computes y <- alpha*x+y. */
    /*     snrm2   Level 1 BLAS that computes the norm of a vector. */
    /*     sgbmv   Level 2 BLAS that computes the band matrix vector product. */

    /* \Author */
    /*     Richard Lehoucq */
    /*     Danny Sorensen */
    /*     Chao Yang */
    /*     Dept. of Computational & */
    /*     Applied Mathematics */
    /*     Rice University */
    /*     Houston, Texas */

    /* \SCCS Information: @(#) */
    /* FILE: nbdr5.F   SID: 2.5   DATE OF SID: 08/26/96   RELEASE: 2 */

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

    /* ------------------------------------------------ */
    /* The number NX is the size of each block diagonal */
    /* of A. The number N(=NX*NX) is the dimension of   */
    /* the matrix.  A standard eigenvalue problem is    */
    /* solved (BMAT = 'I').  NEV is the number of       */
    /* eigenvalues (closest to (SIGMAR,SIGMAI)) to be   */
    /* approximated. Since the shift-invert moded is    */
    /* used, WHICH is set to 'LM'. The user can modify  */
    /* NX, NEV, NCV, SIGMAR, SIGMAI to solve problems   */
    /* of different sizes, and to get different parts   */
    /* the spectrum. However, The following conditions  */
    /* must be satisfied:                               */
    /*                   N <= MAXN                      */
    /*                 NEV <= MAXNEV                    */
    /*           NEV + 2 <= NCV <= MAXNCV               */
    /* ------------------------------------------------ */

    nx = 10;
    n = nx * nx;
    nev = 4;
    ncv = 10;
    if (n > 1000)
    {
        printf(" ERROR with _NBDR5: N is greater than MAXN \n");
        return -1;
    }
    else if (nev > 25)
    {
        printf(" ERROR with _NBDR5: NEV is greater than MAXNEV \n");
        return -1;
    }
    else if (ncv > 50)
    {
        printf(" ERROR with _NBDR5: NCV is greater than MAXNCV \n");
        return -1;
    }
    bmat = "I";
    which = "LM";
    sigmar = .4f;
    sigmai = .6f;

    /* --------------------------------------------------- */
    /* The work array WORKL is used in SNAUPD as           */
    /* workspace.  Its dimension LWORKL is set as          */
    /* illustrated below.  The parameter TOL determines    */
    /* the stopping criterion. If TOL<=0, machine          */
    /* precision is used.  The variable IDO is used for    */
    /* reverse communication, and is initially set to 0.   */
    /* Setting INFO=0 indicates that a random vector is    */
    /* generated in SNAUPD to start the Arnoldi iteration. */
    /* --------------------------------------------------- */

    /* Computing 2nd power */
    i__1 = ncv;
    lworkl = i__1 * i__1 * 3 + ncv * 6;
    tol = 0.f;
    ido = 0;
    info = 0;

    /* ------------------------------------------------- */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 4 of SNAUPD is used     */
    /* (IPARAM(7) = 4). All these options can be changed */
    /* by the user. For details, see the documentation   */
    /* in SNBAND.                                        */
    /* ------------------------------------------------- */

    maxitr = 300;
    mode = 4;

    iparam[2] = maxitr;
    iparam[6] = mode;

    /* ------------------------------------------ */
    /* Construct matrices A and M in LAPACK-style */
    /* banded form.                               */
    /* ------------------------------------------ */

    /* ------------------------------------------- */
    /* Zero out the workspace for banded matrices. */
    /* ------------------------------------------- */

    slaset_("A", &c__50, &n, &c_b15, &c_b15, a, &c__50);
    slaset_("A", &c__50, &n, &c_b15, &c_b15, m, &c__50);
    slaset_("A", &c__50, &n, &c_b15, &c_b15, rfac, &c__50);

    /* ----------------------------------- */
    /* KU, KL are number of superdiagonals */
    /* and subdiagonals within the band of */
    /* matrices A and M.                   */
    /* ----------------------------------- */

    kl = nx;
    ku = nx;

    /* ------------- */
    /* Main diagonal */
    /* ------------- */

    idiag = kl + ku + 1;
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
        a[idiag + j * 50 - 51] = 4.f;
        m[idiag + j * 50 - 51] = 4.f;
        /* L30: */
    }

    /* ----------------------------------- */
    /* First subdiagonal and superdiagonal */
    /* ----------------------------------- */

    isup = kl + ku;
    isub = kl + kl + 2;
    h = 1.f / (float)(nx + 1);
    rho = 100.f;
    i__1 = nx;
    for (i = 1; i <= i__1; ++i)
    {
        lo = (i - 1) * nx;
        i__2 = lo + nx - 1;
        for (j = lo + 1; j <= i__2; ++j)
        {
            a[isup + (j + 1) * 50 - 51] = h * rho / 2.f - 1.f;
            a[isub + j * 50 - 51] = -1.f - h * rho / 2.f;
            /* L40: */
        }
        /* L50: */
    }

    i__1 = n - 1;
    for (j = 1; j <= i__1; ++j)
    {
        m[isup + (j + 1) * 50 - 51] = 1.f;
        m[isub + j * 50 - 51] = 1.f;
        /* L60: */
    }

    /* ---------------------------------- */
    /* KL-th subdiagonal and KU-th super- */
    /* diagonal.                          */
    /* ---------------------------------- */

    isup = kl + 1;
    isub = (kl << 1) + ku + 1;
    i__1 = nx - 1;
    for (i = 1; i <= i__1; ++i)
    {
        lo = (i - 1) * nx;
        i__2 = lo + nx;
        for (j = lo + 1; j <= i__2; ++j)
        {
            a[isup + (nx + j) * 50 - 51] = -1.f;
            a[isub + j * 50 - 51] = -1.f;
            /* L70: */
        }
        /* L80: */
    }

    /* ---------------------------------------------- */
    /* Call ARPACK banded solver to find eigenvalues  */
    /* and eigenvectors. The real parts of the        */
    /* eigenvalues are returned in the first column   */
    /* of D, the imaginary parts are returned in the  */
    /* second column of D.  Eigenvectors are returned */
    /* in the first NCONV (=IPARAM(5)) columns of V.  */
    /* ---------------------------------------------- */

    rvec = TRUE_;
    snband_(&rvec, "A", select, d, &d[50], v, &c__1000, &sigmar, &sigmai, workev, &n, a, m, &c__50, rfac, cfac, &ku, &kl, which, bmat, &nev, &tol, resid, &ncv, v, &c__1000, iparam, workd, workl, &lworkl, workc, iwork, &info);

    if (info == 0)
    {

        /* --------------------------------- */
        /* Print out convergence information */
        /* --------------------------------- */

        nconv = iparam[4];

        printf(" \n");
        printf(" _NBDR5 \n");
        printf(" ====== \n");
        printf(" \n");
        printf(" The size of the matrix is %d", n);
        printf(" Number of eigenvalue requested is %d", nev);
        printf(" The number of Arnoldi vectors generated (NCV) is %d", ncv);
        printf(" The number of converged Ritz values is %d", nconv);
        printf(" What portion of the spectrum %s", which);
        printf(" The number of Implicit Arnoldi update iterations taken is %d", iparam[2]);
        printf(" The number of OP*x is %d", iparam[8]);
        printf(" The convergence tolerance is %e", tol);
        printf(" \n");

        /* -------------------------- */
        /* Compute the residual norm. */
        /*    ||  A*x - lambda*x ||   */
        /* -------------------------- */

        first = TRUE_;
        i__1 = nconv;
        for (j = 1; j <= i__1; ++j)
        {

            if (d[j + 49] == 0.f)
            {

                /* ------------------ */
                /* Ritz value is real */
                /* ------------------ */

                sgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b101, &a[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b15, ax, &c__1);
                r__1 = -d[j - 1];
                saxpy_(&n, &r__1, &v[j * 1000 - 1000], &c__1, ax, &c__1);
                d[j + 99] = snrm2_(&n, ax, &c__1);
                d[j + 99] /= (r__1 = d[j - 1], dabs(r__1));
            }
            else if (first)
            {

                /* ---------------------- */
                /* Ritz value is complex  */
                /* Residual of one Ritz   */
                /* value of the conjugate */
                /* pair is computed.      */
                /* ---------------------- */

                sgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b101, &a[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b15, ax, &c__1);
                r__1 = -d[j - 1];
                saxpy_(&n, &r__1, &v[j * 1000 - 1000], &c__1, ax, &c__1);
                saxpy_(&n, &d[j + 49], &v[(j + 1) * 1000 - 1000], &c__1, ax, &c__1);
                d[j + 99] = snrm2_(&n, ax, &c__1);
                sgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b101, &a[kl], &c__50, &v[(j + 1) * 1000 - 1000], &c__1, &c_b15, ax, &c__1);
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
        smout_(&nconv, &c__3, d, &c__50, &c_n6,"Ritz values (Real,Imag) and relative residuals");
    }
    else
    {

        /* ----------------------------------- */
        /* Either convergence failed, or there */
        /* is error.  Check the documentation  */
        /* for SNBAND.                         */
        /* ----------------------------------- */

        printf(" \n");
        printf(" Error with _nband info= %d", info);
        printf(" Check the documentation of _nband \n");
        printf(" \n");
    }

L9000:
    return 0;
}
