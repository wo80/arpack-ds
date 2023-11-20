/* EXAMPLES\BAND\dsbdr6.f -- translated by f2c (version 20230428). */

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

int main()
{
    /* System generated locals */
    a_int i__1;
    double d__1;

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
    char* bmat;
    a_int mode, info;
    a_bool rvec;
    a_int isub, isup;
    a_int idiag;
    double sigma;
    char* which;
    double resid[1000];
    a_int nconv;
    double workd[3000];
    a_int iwork[1000];
    double workl[2900];
    a_int iparam[11];
    a_bool select[50];
    a_int maxitr, lworkl;

    /*     ... Construct the matrix A in LAPACK-style band form. */
    /*         The matrix A is the 1-dimensional discrete Laplacian on [0,1] */
    /*         with zero Dirichlet boundary condition, M is the mass */
    /*         formed by using piecewise linear elements on [0,1]. */

    /*     ... Call DSBAND  with Cayley mode to find eigenvalues LAMBDA such that */
    /*                          A*x = LAMBDA*M*x. */

    /*     ... Use mode 5 of DSAUPD . */

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
    /* FILE: sbdr6.F   SID: 2.5   DATE OF SID: 08/26/96   RELEASE: 2 */

    /* \Remarks */
    /*     1. None */

    /* \EndLib */

    /* ---------------------------------------------------------------------- */

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
    /* The number N is the dimension of the matrix.  A  */
    /* generalized eigenvalue problem is solved         */
    /* (BMAT = 'G').  NEV is the number of eigenvalues  */
    /* to be approximated.  Since the Cayley mode is    */
    /* used, WHICH is set to 'LM'.  The user can        */
    /* modify N, NEV, NCV and SIGMA to solve problems   */
    /* of different sizes, and to get different parts   */
    /* the spectrum.  However, the following conditions */
    /* must be satisfied:                               */
    /*                   N <= MAXN                      */
    /*                 NEV <= MAXNEV                    */
    /*           NEV + 1 <= NCV <= MAXNCV               */
    /* ------------------------------------------------ */

    n = 100;
    nev = 4;
    ncv = 10;
    if (n > 1000)
    {
        printf(" ERROR with _SBDR6: N is greater than MAXN \n");
        return ierr;
    }
    else if (nev > 25)
    {
        printf(" ERROR with _SBDR6: NEV is greater than MAXNEV \n");
        return ierr;
    }
    else if (ncv > 50)
    {
        printf(" ERROR with _SBDR6: NCV is greater than MAXNCV \n");
        return ierr;
    }
    bmat = "G";
    which = "LM";
    sigma = 150.f;

    /* --------------------------------------------------- */
    /* The work array WORKL is used in DSAUPD  as           */
    /* workspace.  Its dimension LWORKL is set as          */
    /* illustrated below.  The parameter TOL determines    */
    /* the stopping criterion. If TOL<=0, machine          */
    /* precision is used.  The variable IDO is used for    */
    /* reverse communication, and is initially set to 0.   */
    /* Setting INFO=0 indicates that a random vector is    */
    /* generated in DSAUPD  to start the Arnoldi iteration. */
    /* --------------------------------------------------- */

    /* Computing 2nd power */
    i__1 = ncv;
    lworkl = i__1 * i__1 + (ncv << 3);
    tol = 0.;
    ido = 0;
    info = 0;

    /* ------------------------------------------------- */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 5 of DSAUPD  is used     */
    /* (IPARAM(7) = 5). All these options can be changed */
    /* by the user. For details, see the documentation   */
    /* in SBAND.                                         */
    /* ------------------------------------------------- */

    maxitr = 300;
    mode = 5;

    iparam[2] = maxitr;
    iparam[6] = mode;

    /* -------------------------------------- */
    /* Construct the matrix A in LAPACK-style */
    /* banded form.                           */
    /* -------------------------------------- */

    /* ------------------------------------------- */
    /* Zero out the workspace for banded matrices. */
    /* ------------------------------------------- */

    dlaset_("A", &c__50, &n, &c_b15, &c_b15, a, &c__50);
    dlaset_("A", &c__50, &n, &c_b15, &c_b15, m, &c__50);
    dlaset_("A", &c__50, &n, &c_b15, &c_b15, rfac, &c__50);

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

    /* ----------------------------------- */
    /* First subdiagonal and superdiagonal */
    /* ----------------------------------- */

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

    /* ----------------------------------- */
    /* Call DSBAND  to find eigenvalues and */
    /* eigenvectors.  Eigenvalues are      */
    /* returned in the first column of D.  */
    /* Eigenvectors are returned in the    */
    /* first NCONV (=IPARAM(5)) columns of */
    /* V.                                  */
    /* ----------------------------------- */

    rvec = TRUE_;
    dsband_(&rvec, "A", select, d, v, &c__1000, &sigma, &n, a, m, &c__50, rfac, &kl, &ku, which, bmat, &nev, &tol, resid, &ncv, v, &c__1000, iparam, workd, workl, &lworkl, iwork, &info);

    if (info == 0)
    {

        nconv = iparam[4];

        /* --------------------------------- */
        /* Print out convergence information */
        /* --------------------------------- */

        printf(" \n");
        printf(" _SBDR6 \n");
        printf(" ====== \n");
        printf(" \n");
        printf(" The size of the matrix is %d", n);
        printf(" Number of eigenvalue requested is %d", nev);
        printf(" The number of Lanczos vectors generated (NCV) is %d", ncv);
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

        i__1 = nconv;
        for (j = 1; j <= i__1; ++j)
        {
            dgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b97, &a[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b15, ax, &c__1);
            dgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b97, &m[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b15, mx, &c__1);
            d__1 = -d[j - 1];
            daxpy_(&n, &d__1, mx, &c__1, ax, &c__1);
            d[j + 49] = dnrm2_(&n, ax, &c__1);
            d[j + 49] /= (d__1 = d[j - 1], abs(d__1));

            /* L90: */
        }
        dmout_(&nconv, &c__2, d, &c__50, &c_n6,"Ritz values and relative residuals");
    }
    else
    {

        /* ----------------------------------- */
        /* Either convergence failed, or there */
        /* is error.  Check the documentation  */
        /* for DSBAND .                         */
        /* ----------------------------------- */

        printf(" \n");
        printf(" Error with _band info= %d", info);
        printf(" Check the documentation of _band \n");
        printf(" \n");
    }

L9000:
    return 0;
}
