/* EXAMPLES\BAND\ssbdr1.f -- translated by f2c (version 20230428). */

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
static a_int c__2 = 2;
static a_int c_n6 = -6;

int main()
{
    /* System generated locals */
    a_int i__1, i__2;
    float r__1;

    /* Local variables */
    float a[50000] /* was [50][1000] */, d[100] /* was [50][2] */;
    a_int i, j;
    float m[50000] /* was [50][1000] */;
    a_int n;
    float v[50000] /* was [1000][50] */, h2;
    a_int kl;
    float ax[1000];
    a_int lo, ku, nx, ido, ncv, nev;
    float tol, rfac[50000] /* was [50][1000] */;
    char* bmat;
    a_int mode, info;
    a_bool rvec;
    a_int isub, isup;
    a_int idiag;
    float sigma;
    char* which;
    float resid[1000];
    a_int nconv;
    float workd[3000];
    a_int iwork[1000];
    float workl[2900];
    a_int iparam[11];
    a_bool select[50];
    a_int maxitr, lworkl;

    /*     ... Construct the matrix A in LAPACK-style band form. */
    /*         The matrix A is derived from the discretization of */
    /*         the 2-dimensional Laplacian on the unit square with */
    /*         zero Dirichlet boundary condition using standard */
    /*         central difference. */

    /*     ... Call SSBAND to find eigenvalues LAMBDA such that */
    /*                          A*x = x*LAMBDA. */

    /*     ... Use mode 1 of SSAUPD. */

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
    /* FILE: sbdr1.F   SID: 2.5   DATE OF SID: 08/26/96   RELEASE: 2 */

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

    /* ----------------------------------------------- */
    /* The number NX is the number of interior points  */
    /* in the discretization of the 2-dimensional      */
    /* Laplacian operator on the unit square with zero */
    /* Dirichlet boundary condition. The number        */
    /* N(=NX*NX) is the dimension of the matrix.  A    */
    /* standard eigenvalue problem is solved           */
    /* (BMAT = 'I').  NEV is the number of eigenvalues */
    /* to be approximated. The user can modify NX,NEV, */
    /* NCV and WHICH to solve problems of different    */
    /* sizes, and to get different parts the spectrum. */
    /* However, the following conditions must be       */
    /* satisfied:                                      */
    /*                   N <= MAXN                     */
    /*                 NEV <= MAXNEV                   */
    /*           NEV + 1 <= NCV <= MAXNCV              */
    /* ----------------------------------------------- */

    nx = 10;
    n = nx * nx;
    nev = 4;
    ncv = 10;
    if (n > 1000)
    {
        printf(" ERROR with _SBDR1: N is greater than MAXN \n");
        return ierr;
    }
    else if (nev > 25)
    {
        printf(" ERROR with _SBDR1: NEV is greater than MAXNEV \n");
        return ierr;
    }
    else if (ncv > 50)
    {
        printf(" ERROR with _SBDR1: NCV is greater than MAXNCV \n");
        return ierr;
    }
    bmat = "I";
    which = "LM";

    /* --------------------------------------------------- */
    /* The work array WORKL is used in SSAUPD as           */
    /* workspace.  Its dimension LWORKL is set as          */
    /* illustrated below.  The parameter TOL determines    */
    /* the stopping criterion. If TOL<=0, machine          */
    /* precision is used.  The variable IDO is used for    */
    /* reverse communication, and is initially set to 0.   */
    /* Setting INFO=0 indicates that a random vector is    */
    /* generated in SSAUPD to start the Arnoldi iteration. */
    /* --------------------------------------------------- */

    /* Computing 2nd power */
    i__1 = ncv;
    lworkl = i__1 * i__1 + (ncv << 3);
    tol = 0.f;
    ido = 0;
    info = 0;

    /* ------------------------------------------------- */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 1 of SSAUPD is used     */
    /* (IPARAM(7) = 1). All these options can be changed */
    /* by the user. For details see the documentation in */
    /* SSBAND.                                           */
    /* ------------------------------------------------- */

    maxitr = 300;
    mode = 1;

    iparam[2] = maxitr;
    iparam[6] = mode;

    /* -------------------------------------- */
    /* Construct the matrix A in LAPACK-style */
    /* banded form.                           */
    /* -------------------------------------- */

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

    h2 = 1.f / ((nx + 1) * (nx + 1));
    idiag = kl + ku + 1;
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
        a[idiag + j * 50 - 51] = 4.f / h2;
        /* L30: */
    }

    /* ----------------------------------- */
    /* First subdiagonal and superdiagonal */
    /* ----------------------------------- */

    isup = kl + ku;
    isub = kl + ku + 2;
    i__1 = nx;
    for (i = 1; i <= i__1; ++i)
    {
        lo = (i - 1) * nx;
        i__2 = lo + nx - 1;
        for (j = lo + 1; j <= i__2; ++j)
        {
            a[isup + (j + 1) * 50 - 51] = -1.f / h2;
            a[isub + j * 50 - 51] = -1.f / h2;
            /* L40: */
        }
        /* L50: */
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
            a[isup + (nx + j) * 50 - 51] = -1.f / h2;
            a[isub + j * 50 - 51] = -1.f / h2;
            /* L70: */
        }
        /* L80: */
    }

    /* ----------------------------------- */
    /* Call SSBAND to find eigenvalues and */
    /* eigenvectors.  Eigenvalues are      */
    /* returned in the first column of D.  */
    /* Eigenvectors are returned in the    */
    /* first NCONV (=IPARAM(5)) columns of */
    /* V.                                  */
    /* ----------------------------------- */

    rvec = TRUE_;
    ssband_(&rvec, "A", select, d, v, &c__1000, &sigma, &n, a, m, &c__50, rfac, &kl, &ku, which, bmat, &nev, &tol, resid, &ncv, v, &c__1000, iparam, workd, workl, &lworkl, iwork, &info);

    if (info == 0)
    {

        nconv = iparam[4];

        /* --------------------------------- */
        /* Print out convergence information */
        /* --------------------------------- */

        printf(" \n");
        printf(" _SBDR1 \n");
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
            sgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b100, &a[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b15, ax, &c__1);
            r__1 = -d[j - 1];
            saxpy_(&n, &r__1, &v[j * 1000 - 1000], &c__1, ax, &c__1);
            d[j + 49] = snrm2_(&n, ax, &c__1);
            d[j + 49] /= (r__1 = d[j - 1], dabs(r__1));

            /* L90: */
        }
        smout_(&nconv, &c__2, d, &c__50, &c_n6,"Ritz values and relative residuals");
    }
    else
    {

        /* ----------------------------------- */
        /* Either convergence failed, or there */
        /* is error.  Check the documentation  */
        /* for SSBAND.                         */
        /* ----------------------------------- */

        printf(" \n");
        printf(" Error with _sband info= %d", info);
        printf(" Check the documentation of _sband \n");
        printf(" \n");
    }

L9000:
    return 0;
}
