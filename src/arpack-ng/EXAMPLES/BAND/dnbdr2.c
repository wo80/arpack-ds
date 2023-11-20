/* EXAMPLES\BAND\dnbdr2.f -- translated by f2c (version 20230428). */

#include <stdlib.h>
#include "arpack_internal.h"

/* Table of constant values */

static a_int c__9 = 9;
static a_int c__1 = 1;
static a_int c__50 = 50;
static double c_b15 = 0.;
static a_int c__1000 = 1000;
static a_int c__3 = 3;
static a_int c__5 = 5;
static double c_b100 = 1.;
static a_int c__6 = 6;
static a_int c_n6 = -6;

int main()
{
    /* System generated locals */
    a_int i__1, i__2;
    double d__1;

    /* Local variables */
    a_bool select[50];
    a_int iparam[11];
    a_bool rvec, first;
    a_int i, j, n, kl, ku, lo, nx;
    a_int ido, ncv, nev, info, isub, isup, mode;
    a_int idiag, nconv, lworkl, maxitr;
    char *bmat, *which;
    double h, h2, rho, tol, sigmai, sigmar;

    /*     ... Construct matrices A in LAPACK-style band form. */
    /*         The matrix A is derived from the discretization of */
    /*         the 2-d convection-diffusion operator */

    /*               -Laplacian(u) + rho*partial(u)/partial(x). */

    /*         on the unit square with zero Dirichlet boundary condition */
    /*         using standard central difference. */

    /*     ... Define the shift SIGMA = (SIGMAR, SIGMAI). */

    /*     ... Call DNBAND  to find eigenvalues LAMBDA closest to SIGMA */
    /*         such that */
    /*                       A*x = LAMBDA*x. */

    /*     ... Use mode 3 of DNAUPD . */

    /* \BeginLib */

    /* \Routines called: */
    /*     dnband   ARPACK banded eigenproblem solver. */
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
    /* FILE: nbdr2.F   SID: 2.5   DATE OF SID: 08/26/96   RELEASE: 2 */

    /* \Remarks */
    /*     1. None */

    /* \EndLib */

    /* --------------------------------------------------------------------- */

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
    /* eigenvalues (closest to (SIGMAR,SIGMAI)) to be  */
    /* approximated. Since the shift-invert moded is   */
    /* used, WHICH is set to 'LM'. The user can modify */
    /* NX, NEV, NCV, SIGMAR, SIGMAI to solve problems  */
    /* of different sizes, and to get different parts  */
    /* the spectrum. However, The following conditions */
    /* must be satisfied:                              */
    /*                   N <= MAXN                     */
    /*                 NEV <= MAXNEV                   */
    /*           NEV + 2 <= NCV <= MAXNCV              */
    /* ----------------------------------------------- */

    nx = 10;
    n = nx * nx;
    nev = 4;
    ncv = 20;
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
    sigmar = 1e4;
    sigmai = 0.;

    /* --------------------------------------------------- */
    /* The work array WORKL is used in DNAUPD  as           */
    /* workspace.  Its dimension LWORKL is set as          */
    /* illustrated below.  The parameter TOL determines    */
    /* the stopping criterion. If TOL<=0, machine          */
    /* precision is used.  The variable IDO is used for    */
    /* reverse communication, and is initially set to 0.   */
    /* Setting INFO=0 indicates that a random vector is    */
    /* generated in DNAUPD  to start the Arnoldi iteration. */
    /* --------------------------------------------------- */

    lworkl = ncv * ncv * 3 + ncv * 6;
    tol = 0.;
    ido = 0;
    info = 0;

    a_dcomplex* cfac = (a_dcomplex*)malloc(sizeof(a_dcomplex) * 50 * 1000);
    a_dcomplex* workc = (a_dcomplex*)malloc(sizeof(a_dcomplex) * 1000);
    a_int* iwork = (a_int*)malloc(sizeof(a_int) * 1000);
    double* a = (double*)malloc(sizeof(double) * 50 * 1000);
    double* d = (double*)malloc(sizeof(double) * 50 * 3);
    double* m = (double*)malloc(sizeof(double) * 50 * 1000);
    double* v = (double*)malloc(sizeof(double) * 1000 * 50);
    double* ax = (double*)malloc(sizeof(double) * 1000);
    double* rfac = (double*)malloc(sizeof(double) * 50 * 1000);
    double* resid = (double*)malloc(sizeof(double) * 1000);
    double* workd = (double*)malloc(sizeof(double) * 3000);
    double* workl = (double*)malloc(sizeof(double) * 7800);
    double* workev = (double*)malloc(sizeof(double) * 150);

    /* ------------------------------------------------- */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 3 of DNAUPD  is used     */
    /* (IPARAM(7) = 3). All these options can be changed */
    /* by the user. For details, see the documentation   */
    /* in DNBAND .                                        */
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

    dlaset_("A", &c__50, &n, &c_b15, &c_b15, a, &c__50);
    dlaset_("A", &c__50, &n, &c_b15, &c_b15, m, &c__50);
    dlaset_("A", &c__50, &n, &c_b15, &c_b15, rfac, &c__50);

    /* ----------------------------------- */
    /* KU, KL are number of superdiagonals */
    /* and subdiagonals within the band of */
    /* matrices A.                   */
    /* ----------------------------------- */

    kl = nx;
    ku = nx;

    /* ------------- */
    /* Main diagonal */
    /* ------------- */

    h = 1. / (double)(nx + 1);
    h2 = h * h;

    idiag = kl + ku + 1;
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
        a[idiag + j * 50 - 51] = 4. / h2;
        /* L30: */
    }

    /* ----------------------------------- */
    /* First subdiagonal and superdiagonal */
    /* ----------------------------------- */

    isup = kl + ku;
    isub = kl + ku + 2;
    rho = 10.;
    i__1 = nx;
    for (i = 1; i <= i__1; ++i)
    {
        lo = (i - 1) * nx;
        i__2 = lo + nx - 1;
        for (j = lo + 1; j <= i__2; ++j)
        {
            a[isub + (j + 1) * 50 - 51] = -1. / h2 + rho / 2. / h;
            a[isup + j * 50 - 51] = -1. / h2 - rho / 2. / h;
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
            a[isup + (nx + j) * 50 - 51] = -1. / h2;
            a[isub + j * 50 - 51] = -1. / h2;
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
    dnband_(&rvec, "A", select, d, &d[50], v, &c__1000, &sigmar, &sigmai, workev, &n, a, m, &c__50, rfac, cfac, &kl, &ku, which, bmat, &nev, &tol, resid, &ncv, v, &c__1000, iparam, workd, workl, &lworkl, workc, iwork, &info);

    if (info == 0)
    {

        /* --------------------------------- */
        /* Print out convergence information */
        /* --------------------------------- */

        nconv = iparam[4];

        printf(" \n");
        printf(" _NBDR2 \n");
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

            if (d[j + 49] == 0.)
            {

                /* ------------------ */
                /* Ritz value is real */
                /* ------------------ */

                dgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b100, &a[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b15, ax, &c__1);
                d__1 = -d[j - 1];
                daxpy_(&n, &d__1, &v[j * 1000 - 1000], &c__1, ax, &c__1);
                d[j + 99] = dnrm2_(&n, ax, &c__1);
                d[j + 99] /= (d__1 = d[j - 1], abs(d__1));
            }
            else if (first)
            {

                /* ---------------------- */
                /* Ritz value is complex  */
                /* Residual of one Ritz   */
                /* value of the conjugate */
                /* pair is computed.      */
                /* ---------------------- */

                dgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b100, &a[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b15, ax, &c__1);
                d__1 = -d[j - 1];
                daxpy_(&n, &d__1, &v[j * 1000 - 1000], &c__1, ax, &c__1);
                daxpy_(&n, &d[j + 49], &v[(j + 1) * 1000 - 1000], &c__1, ax, &c__1);
                d[j + 99] = dnrm2_(&n, ax, &c__1);
                dgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b100, &a[kl], &c__50, &v[(j + 1) * 1000 - 1000], &c__1, &c_b15, ax, &c__1);
                d__1 = -d[j - 1];
                daxpy_(&n, &d__1, &v[(j + 1) * 1000 - 1000], &c__1, ax, &c__1);
                d__1 = -d[j + 49];
                daxpy_(&n, &d__1, &v[j * 1000 - 1000], &c__1, ax, &c__1);
                d__1 = dnrm2_(&n, ax, &c__1);
                d[j + 99] = dlapy2_(&d[j + 99], &d__1);
                d[j + 99] /= dlapy2_(&d[j - 1], &d[j + 49]);
                d[j + 100] = d[j + 99];
                first = FALSE_;
            }
            else
            {
                first = TRUE_;
            }

            /* L90: */
        }
        dmout_(nconv, 3, d, 50, -6, "Ritz values (Real,Imag) and relative residuals");
    }
    else
    {

        /* ----------------------------------- */
        /* Either convergence failed, or there */
        /* is error.  Check the documentation  */
        /* for DNBAND .                         */
        /* ----------------------------------- */

        printf(" \n");
        printf(" Error with _nband info= %d", info);
        printf(" Check the documentation of _nband \n");
        printf(" \n");
    }

L9000:

    free(cfac);
    free(workc);
    free(iwork);
    free(a);
    free(d);
    free(m);
    free(v);
    free(ax);
    free(rfac);
    free(resid);
    free(workd);
    free(workl);
    free(workev);

    return 0;
}
