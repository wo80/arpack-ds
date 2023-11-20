/* EXAMPLES\BAND\ssbdr5.f -- translated by f2c (version 20230428). */

#include <stdlib.h>
#include "arpack_internal.h"

/* Table of constant values */

static a_int c__9 = 9;
static a_int c__1 = 1;
static a_int c__50 = 50;
static float c_b15 = 0.f;
static a_int c__1000 = 1000;
static a_int c__3 = 3;
static a_int c__4 = 4;
static float c_b97 = 1.f;
static a_int c__6 = 6;
static a_int c__2 = 2;
static a_int c_n6 = -6;

int main()
{
    /* System generated locals */
    a_int i__1;
    float r__1;

    /* Local variables */
    a_bool select[50];
    a_int iparam[11];
    a_bool rvec;
    a_int j, n, kl, ku, ido, ncv, nev;
    a_int info, isub, isup, mode, idiag, nconv, lworkl;
    a_int maxitr;
    char *bmat, *which;
    float h, r1, r2, tol, sigma;

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
    /* to be approximated.  Since the Buckling mode is  */
    /* is used, WHICH is set to 'LM'.  The user can     */
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
        printf(" ERROR with _SBDR5: N is greater than MAXN \n");
        return -1;
    }
    else if (nev > 25)
    {
        printf(" ERROR with _SBDR5: NEV is greater than MAXNEV \n");
        return -1;
    }
    else if (ncv > 50)
    {
        printf(" ERROR with _SBDR5: NCV is greater than MAXNCV \n");
        return -1;
    }
    bmat = "G";
    which = "LM";
    sigma = 1.f;

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

    lworkl = ncv * ncv + (ncv << 3);
    tol = 0.f;
    ido = 0;
    info = 0;

    a_int* iwork = (a_int*)malloc(sizeof(a_int) * 1000);
    float* a = (float*)malloc(sizeof(float) * 50 * 1000);
    float* d = (float*)malloc(sizeof(float) * 50 * 2);
    float* m = (float*)malloc(sizeof(float) * 50 * 1000);
    float* v = (float*)malloc(sizeof(float) * 1000 * 50);
    float* ax = (float*)malloc(sizeof(float) * 1000);
    float* mx = (float*)malloc(sizeof(float) * 1000);
    float* rfac = (float*)malloc(sizeof(float) * 50 * 1000);
    float* resid = (float*)malloc(sizeof(float) * 1000);
    float* workd = (float*)malloc(sizeof(float) * 3000);
    float* workl = (float*)malloc(sizeof(float) * 2900);

    /* ------------------------------------------------- */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 4 of SSAUPD is used     */
    /* (IPARAM(7) = 4). All these options can be changed */
    /* by the user. For details see the documentation in */
    /* SSBAND.                                           */
    /* ------------------------------------------------- */

    maxitr = 300;
    mode = 4;

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

    kl = 1;
    ku = 1;

    /* ------------- */
    /* Main diagonal */
    /* ------------- */

    h = 1.f / (float)(n + 1);
    r1 = .66666666666666663f;
    idiag = kl + ku + 1;
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
        a[idiag + j * 50 - 51] = 2.f / h;
        m[idiag + j * 50 - 51] = r1 * h;
        /* L30: */
    }

    /* ----------------------------------- */
    /* First subdiagonal and superdiagonal */
    /* ----------------------------------- */

    r2 = .16666666666666666f;
    isup = kl + ku;
    isub = kl + ku + 2;
    i__1 = n - 1;
    for (j = 1; j <= i__1; ++j)
    {
        a[isup + (j + 1) * 50 - 51] = -1.f / h;
        a[isub + j * 50 - 51] = -1.f / h;
        m[isup + (j + 1) * 50 - 51] = r2 * h;
        m[isub + j * 50 - 51] = r2 * h;
        /* L60: */
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
        printf(" _SBDR5 \n");
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
            sgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b97, &a[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b15, ax, &c__1);
            sgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b97, &m[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b15, mx, &c__1);
            r__1 = -d[j - 1];
            saxpy_(&n, &r__1, mx, &c__1, ax, &c__1);
            d[j + 49] = snrm2_(&n, ax, &c__1);
            d[j + 49] /= (r__1 = d[j - 1], dabs(r__1));

            /* L90: */
        }
        smout_(nconv, 2, d, 50, -6, "Ritz values and relative residuals");
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

    free(iwork);
    free(a);
    free(d);
    free(m);
    free(v);
    free(ax);
    free(mx);
    free(rfac);
    free(resid);
    free(workd);
    free(workl);

    return 0;
}
