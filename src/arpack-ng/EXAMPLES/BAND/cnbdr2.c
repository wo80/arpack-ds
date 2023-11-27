/* EXAMPLES\BAND\cnbdr2.f -- translated by f2c (version 20230428). */

#include <stdlib.h>
#include "arpack_internal.h"

static a_int c__1 = 1;
static a_int c__50 = 50;
static a_int c__1000 = 1000;

static a_fcomplex zero = {0.f, 0.f};
static a_fcomplex one = {1.f, 0.f};
static a_fcomplex two = {2.f, 0.f};
static a_fcomplex four = {4.f, 0.f};

/**
 * \BeginDoc
 *
 * Construct the matrix A in LAPACK-style band form.
 * The matrix A is derived from the discretization of
 * the 2-d convection-diffusion operator
 *
 *      -Laplacian(u) + rho*partial(u)/partial(x).
 *
 * on the unit square with zero Dirichlet boundary condition
 * using standard central difference.
 *
 * Call CNBAND to find eigenvalues LAMBDA such that
 *                  A*x = x*LAMBDA.
 *
 * Use mode 3 of CNAUPD.
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * Routines called:
 *     cnband  ARPACK banded eigenproblem solver.
 *     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     claset  LAPACK routine to initialize a matrix to zero.
 *     caxpy   Level 1 BLAS that computes y <- alpha*x+y.
 *     scnrm2  Level 1 BLAS that computes the norm of a vector.
 *     cgbmv   Level 2 BLAS that computes the band matrix vector product
 *
 * \EndLib */
int main()
{
    /* System generated locals */
    a_int i__1, i__2, i__3;
    a_fcomplex q__1, q__2, q__3, q__4, q__5;

    /* Local variables */
    a_bool select[50];
    a_int iparam[11];
    a_bool rvec;
    a_fcomplex h, h2, rho, sigma;
    a_int i, j, n, kl, ku, lo, ncv;
    a_int nev, nxi, info, isub, isup, mode, idiag;
    a_int nconv, lworkl, maxitr;
    char *bmat, *which;
    float tol;

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
    /* convection-diffusion operator on the unit       */
    /* square with zero Dirichlet boundary condition.  */
    /* The number N(=NX*NX) is the dimension of the    */
    /* matrix.  A standard eigenvalue problem is       */
    /* solved (BMAT = 'I').  NEV is the number of      */
    /* eigenvalues (closest to SIGMA) to be            */
    /* approximated.  Since the shift and invert mode  */
    /* is used, WHICH is set to 'LM'.  The user can    */
    /* modify NX, NEV and NCV to solve problems of     */
    /* different sizes, and to get different parts the */
    /* spectrum.  However, the following conditions    */
    /* must be satisfied:                              */
    /*                   N <= MAXN                     */
    /*                 NEV <= MAXNEV                   */
    /*           NEV + 2 <= NCV <= MAXNCV              */
    /* ----------------------------------------------- */

    nxi = 10;
    n = nxi * nxi;
    nev = 4;
    ncv = 10;
    if (n > 1000)
    {
        printf(" ERROR with _NBDR2: N is greater than MAXN \n");
        return -1;
    }
    else if (nev > 25)
    {
        printf(" ERROR with _NBDR2: NEV is greater than MAXNEV \n");
        return -1;
    }
    else if (ncv > 50)
    {
        printf(" ERROR with _NBDR2: NCV is greater than MAXNCV \n");
        return -1;
    }
    bmat = "I";
    which = "LM";
    sigma.r = 0.f, sigma.i = 0.f;

    /* --------------------------------------------------- */
    /* The work array WORKL is used in CNAUPD as           */
    /* workspace.  Its dimension LWORKL is set as          */
    /* illustrated below.  The parameter TOL determines    */
    /* the stopping criterion. If TOL<=0, machine          */
    /* precision is used.  Setting INFO=0 indicates that a */
    /* random vector is generated in CNAUPD to start the   */
    /* Arnoldi iteration.                                  */
    /* --------------------------------------------------- */

    lworkl = ncv * ncv * 3 + ncv * 5;
    tol = 0.f;
    info = 0;
    nconv = 0;

    a_fcomplex* a = (a_fcomplex*)malloc(sizeof(a_fcomplex) * 50 * 1000);
    a_fcomplex* d = (a_fcomplex*)malloc(sizeof(a_fcomplex) * 50);
    a_fcomplex* m = (a_fcomplex*)malloc(sizeof(a_fcomplex) * 50 * 1000);
    a_fcomplex* v = (a_fcomplex*)malloc(sizeof(a_fcomplex) * 1000 * 50);
    a_fcomplex* ax = (a_fcomplex*)malloc(sizeof(a_fcomplex) * 1000);
    a_fcomplex* fac = (a_fcomplex*)malloc(sizeof(a_fcomplex) * 50 * 1000);
    a_fcomplex* resid = (a_fcomplex*)malloc(sizeof(a_fcomplex) * 1000);
    a_fcomplex* workd = (a_fcomplex*)malloc(sizeof(a_fcomplex) * 3000);
    a_fcomplex* workl = (a_fcomplex*)malloc(sizeof(a_fcomplex) * 7750);
    a_fcomplex* workev = (a_fcomplex*)malloc(sizeof(a_fcomplex) * 100);
    a_int* iwork = (a_int*)malloc(sizeof(a_int) * 1000);
    float* rd = (float*)malloc(sizeof(float) * 50 * 3);
    float* rwork = (float*)malloc(sizeof(float) * 1000);

    /* ------------------------------------------------- */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 3 of CNAUPD is used     */
    /* (IPARAM(7) = 3). All these options can be changed */
    /* by the user. For details, see the documentation   */
    /* in cnband.                                        */
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

    claset_("A", &c__50, &n, &zero, &zero, a, &c__50);
    claset_("A", &c__50, &n, &zero, &zero, m, &c__50);
    claset_("A", &c__50, &n, &zero, &zero, fac, &c__50);

    /* ----------------------------------- */
    /* KU, KL are number of superdiagonals */
    /* and subdiagonals within the band of */
    /* matrices A and M.                   */
    /* ----------------------------------- */

    kl = nxi;
    ku = nxi;

    /* ------------- */
    /* Main diagonal */
    /* ------------- */

    i__1 = nxi + 1;
    q__2.r = (float)i__1, q__2.i = 0.f;
    ar_c_div(&q__1, &one, &q__2);
    h.r = q__1.r, h.i = q__1.i;
    q__1.r = h.r * h.r - h.i * h.i, q__1.i = h.r * h.i + h.i * h.r;
    h2.r = q__1.r, h2.i = q__1.i;

    idiag = kl + ku + 1;
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
        i__2 = idiag + j * 50 - 51;
        ar_c_div(&q__1, &four, &h2);
        a[i__2].r = q__1.r, a[i__2].i = q__1.i;
        m[i__2].r = 1.f, m[i__2].i = 0.f;
    }

    /* ----------------------------------- */
    /* First subdiagonal and superdiagonal */
    /* ----------------------------------- */

    rho.r = 100.f, rho.i = 0.f;
    isup = kl + ku;
    isub = kl + ku + 2;
    i__1 = nxi;
    for (i = 1; i <= i__1; ++i)
    {
        lo = (i - 1) * nxi;
        i__2 = lo + nxi - 1;
        for (j = lo + 1; j <= i__2; ++j)
        {
            i__3 = isup + (j + 1) * 50 - 51;
            q__3.r = -1.f, q__3.i = -0.f;
            ar_c_div(&q__2, &q__3, &h2);
            ar_c_div(&q__5, &rho, &two);
            ar_c_div(&q__4, &q__5, &h);
            q__1.r = q__2.r + q__4.r, q__1.i = q__2.i + q__4.i;
            a[i__3].r = q__1.r, a[i__3].i = q__1.i;
            i__3 = isub + j * 50 - 51;
            q__3.r = -1.f, q__3.i = -0.f;
            ar_c_div(&q__2, &q__3, &h2);
            ar_c_div(&q__5, &rho, &two);
            ar_c_div(&q__4, &q__5, &h);
            q__1.r = q__2.r - q__4.r, q__1.i = q__2.i - q__4.i;
            a[i__3].r = q__1.r, a[i__3].i = q__1.i;
        }
    }

    /* ---------------------------------- */
    /* KL-th subdiagonal and KU-th super- */
    /* diagonal.                          */
    /* ---------------------------------- */

    isup = kl + 1;
    isub = (kl << 1) + ku + 1;
    i__1 = nxi - 1;
    for (i = 1; i <= i__1; ++i)
    {
        lo = (i - 1) * nxi;
        i__2 = lo + nxi;
        for (j = lo + 1; j <= i__2; ++j)
        {
            i__3 = isup + (nxi + j) * 50 - 51;
            q__2.r = -1.f, q__2.i = -0.f;
            ar_c_div(&q__1, &q__2, &h2);
            a[i__3].r = q__1.r, a[i__3].i = q__1.i;
            i__3 = isub + j * 50 - 51;
            q__2.r = -1.f, q__2.i = -0.f;
            ar_c_div(&q__1, &q__2, &h2);
            a[i__3].r = q__1.r, a[i__3].i = q__1.i;
        }
    }

    /* --------------------------------------------- */
    /* Call ARPACK banded solver to find eigenvalues */
    /* and eigenvectors. Eigenvalues are returned in */
    /* the one dimensional array D.  Eigenvectors    */
    /* are returned in the first NCONV (=IPARAM(5))  */
    /* columns of V.                                 */
    /* --------------------------------------------- */

    rvec = TRUE_;
    cnband_(&rvec, "A", select, d, v, &c__1000, &sigma, workev, &n, a, m, &c__50, fac, &kl, &ku, which, bmat, &nev, &tol, resid, &ncv, v, &c__1000, iparam, workd, workl, &lworkl, rwork, iwork, &info);

    if (info == 0)
    {

        nconv = iparam[4];

        /* --------------------------------- */
        /* Print out convergence information */
        /* --------------------------------- */

        printf(" \n");
        printf("_NBDR2 \n");
        printf("====== \n");
        printf(" \n");
        printf(" The size of the matrix is %d\n", n);
        printf(" Number of eigenvalue requested is %d\n", nev);
        printf(" The number of Arnoldi vectors generated (NCV) is %d\n", ncv);
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

            /* ------------------------- */
            /* Compute the residual norm */
            /*   ||  A*x - lambda*x ||   */
            /* ------------------------- */

            cgbmv_("N", &n, &n, &kl, &ku, &one, &a[kl], &c__50, &v[j * 1000 - 1000], &c__1, &zero, ax, &c__1);
            i__2 = j - 1;
            q__1.r = -d[i__2].r, q__1.i = -d[i__2].i;
            caxpy_(&n, &q__1, &v[j * 1000 - 1000], &c__1, ax, &c__1);
            i__2 = j - 1;
            rd[j - 1] = d[i__2].r;
            rd[j + 49] = d[j - 1].i;
            rd[j + 99] = scnrm2_(&n, ax, &c__1);
            rd[j + 99] /= slapy2_(&rd[j - 1], &rd[j + 49]);
        }
        smout_(nconv, 3, rd, 50, -6, "Ritz values (Real,Imag) and relative residuals");
    }
    else
    {

        /* ----------------------------------- */
        /* Either convergence failed, or there */
        /* is error.  Check the documentation  */
        /* for cnband.                         */
        /* ----------------------------------- */

        printf(" \n");
        printf(" Error with _nband info= %d\n", info);
        printf(" Check the documentation of _nband \n");
        printf(" \n");
    }

    free(a);
    free(d);
    free(m);
    free(v);
    free(ax);
    free(fac);
    free(resid);
    free(workd);
    free(workl);
    free(workev);
    free(iwork);
    free(rd);
    free(rwork);

    return nconv < nev ? EXIT_FAILURE : EXIT_SUCCESS;
}
