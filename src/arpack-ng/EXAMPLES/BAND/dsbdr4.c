/* EXAMPLES\BAND\dsbdr4.f -- translated by f2c (version 20230428). */

#include <stdlib.h>
#include "arpack_internal.h"

static a_int c__1 = 1;
static a_int c__50 = 50;
static a_int c__1000 = 1000;

static double zero = 0.;
static double one = 1.;

/**
 * \BeginDoc
 *
 * Construct the matrix A in LAPACK-style band form.
 * The matrix A is the 1-dimensional discrete Laplacian on [0,1]
 * with zero Dirichlet boundary condition, M is the mass
 * formed by using piecewise linear elements on [0,1].
 *
 * Call DSBAND  with shift-invert mode to find eigenvalues LAMBDA
 * closest to SIGMA such that
 *                  A*x = LAMBDA*M*x.
 *
 * Use mode 3 of DSAUPD .
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * Routines called:
 *     dsband   ARPACK banded eigenproblem solver.
 *     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     dlaset   LAPACK routine to initialize a matrix to zero.
 *     daxpy    Level 1 BLAS that computes y <- alpha*x+y.
 *     dnrm2    Level 1 BLAS that computes the norm of a vector.
 *     dgbmv    Level 2 BLAS that computes the band matrix vector product
 *
 * \EndLib */
int main()
{
    /* System generated locals */
    a_int i__1;
    double d__1;

    /* Local variables */
    a_bool select[50];
    a_int iparam[11];
    a_bool rvec;
    a_int j, n, kl, ku, ido, ncv, nev;
    a_int info, isub, isup, mode, idiag, nconv, lworkl;
    a_int maxitr;
    char *bmat, *which;
    double h, r1, r2, tol, sigma;

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
    /* (closest to the shift SIGMA) to be approximated. */
    /* Since the shift and invert mode is used, WHICH   */
    /* is set to 'LM'.  The user can modify N, NEV, NCV */
    /* and SIGMA to solve problems of different sizes,  */
    /* and to get different parts the spectrum.         */
    /* However, the following conditions must be        */
    /* satisfied:                                       */
    /*                   N <= MAXN                      */
    /*                 NEV <= MAXNEV                    */
    /*           NEV + 1 <= NCV <= MAXNCV               */
    /* ------------------------------------------------ */

    n = 100;
    nev = 4;
    ncv = 10;
    if (n > 1000)
    {
        printf(" ERROR with _SBDR4: N is greater than MAXN \n");
        return -1;
    }
    else if (nev > 25)
    {
        printf(" ERROR with _SBDR4: NEV is greater than MAXNEV \n");
        return -1;
    }
    else if (ncv > 50)
    {
        printf(" ERROR with _SBDR4: NCV is greater than MAXNCV \n");
        return -1;
    }
    bmat = "G";
    which = "LM";
    sigma = 0.;

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

    lworkl = ncv * ncv + (ncv << 3);
    tol = 0.;
    ido = 0;
    info = 0;
    nconv = 0;

    a_int* iwork = (a_int*)malloc(sizeof(a_int) * 1000);
    double* a = (double*)malloc(sizeof(double) * 50 * 1000);
    double* d = (double*)malloc(sizeof(double) * 50 * 2);
    double* m = (double*)malloc(sizeof(double) * 50 * 1000);
    double* v = (double*)malloc(sizeof(double) * 1000 * 50);
    double* ax = (double*)malloc(sizeof(double) * 1000);
    double* mx = (double*)malloc(sizeof(double) * 1000);
    double* rfac = (double*)malloc(sizeof(double) * 50 * 1000);
    double* resid = (double*)malloc(sizeof(double) * 1000);
    double* workd = (double*)malloc(sizeof(double) * 3000);
    double* workl = (double*)malloc(sizeof(double) * 2900);

    /* ------------------------------------------------- */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 3 of DSAUPD  is used     */
    /* (IPARAM(7) = 3). All these options can be changed */
    /* by the user. For details, see the documentation   */
    /* in DSBAND .                                        */
    /* ------------------------------------------------- */

    maxitr = 300;
    mode = 3;

    iparam[2] = maxitr;
    iparam[6] = mode;

    /* -------------------------------------- */
    /* Construct the matrix A in LAPACK-style */
    /* banded form.                           */
    /* -------------------------------------- */

    /* ------------------------------------------- */
    /* Zero out the workspace for banded matrices. */
    /* ------------------------------------------- */

    dlaset_("A", &c__50, &n, &zero, &zero, a, &c__50);
    dlaset_("A", &c__50, &n, &zero, &zero, m, &c__50);
    dlaset_("A", &c__50, &n, &zero, &zero, rfac, &c__50);

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
        printf(" _SBDR4 \n");
        printf(" ====== \n");
        printf(" \n");
        printf(" The size of the matrix is %d\n", n);
        printf(" Number of eigenvalue requested is %d\n", nev);
        printf(" The number of Lanczos vectors generated (NCV) is %d\n", ncv);
        printf(" The number of converged Ritz values is %d\n", nconv);
        printf(" What portion of the spectrum %s\n", which);
        printf(" The number of Implicit Arnoldi update iterations taken is %d\n", iparam[2]);
        printf(" The number of OP*x is %d\n", iparam[8]);
        printf(" The convergence tolerance is %e\n", tol);
        printf(" \n");

        /* -------------------------- */
        /* Compute the residual norm. */
        /*    ||  A*x - lambda*x ||   */
        /* -------------------------- */

        i__1 = nconv;
        for (j = 1; j <= i__1; ++j)
        {
            dgbmv_("N", &n, &n, &kl, &ku, &one, &a[kl], &c__50, &v[j * 1000 - 1000], &c__1, &zero, ax, &c__1);
            dgbmv_("N", &n, &n, &kl, &ku, &one, &m[kl], &c__50, &v[j * 1000 - 1000], &c__1, &zero, mx, &c__1);
            d__1 = -d[j - 1];
            daxpy_(&n, &d__1, mx, &c__1, ax, &c__1);
            d[j + 49] = dnrm2_(&n, ax, &c__1);
            d[j + 49] /= (d__1 = d[j - 1], abs(d__1));
        }
        dmout_(nconv, 2, d, 50, -6, "Ritz values and relative residuals");
    }
    else
    {

        /* ----------------------------------- */
        /* Either convergence failed, or there */
        /* is error.  Check the documentation  */
        /* for DSBAND .                         */
        /* ----------------------------------- */

        printf(" \n");
        printf(" Error with _sband info= %d\n", info);
        printf(" Check the documentation of _sband \n");
        printf(" \n");
    }

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

    return nconv < nev ? EXIT_FAILURE : EXIT_SUCCESS;
}
