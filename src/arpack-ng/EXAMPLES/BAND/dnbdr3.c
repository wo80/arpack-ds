/* EXAMPLES\BAND\dnbdr3.f -- translated by f2c (version 20230428). */

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
static a_int c_n6 = -6;

int main()
{
    /* System generated locals */
    a_int i__1;
    double d__1;

    /* Local variables */
    double a[50000] /* was [50][1000] */, d[150] /* was [50][3] */, h;
    a_int j;
    double m[50000] /* was [50][1000] */;
    a_int n;
    double v[50000] /* was [1000][50] */;
    a_int kl;
    double ax[1000];
    a_int ku;
    double mx[1000];
    a_int ido, ncv, nev;
    double rho, tol;
    a_dcomplex cfac[50000] /* was [50][1000] */;
    double rfac[50000] /* was [50][1000] */;
    char* bmat;
    a_int mode, info;
    a_bool rvec;
    a_int isub, isup;
    a_int idiag;
    char* which;
    double resid[1000];
    a_int nconv;
    a_dcomplex workc[1000];
    double workd[3000];
    a_bool first;
    a_int iwork[1000];
    double workl[7800];
    a_int iparam[11];
    double sigmai;
    a_bool select[50];
    double sigmar;
    a_int maxitr, lworkl;
    double workev[150];

    /*     ... Construct matrices A and M in LAPACK-style band form. */
    /*         The matrix A and M are derived from the finite element */
    /*         discretization of the 1-dimensional convection-diffusion operator */
    /*                         (d^2u/dx^2) + rho*(du/dx) */
    /*         on the interval [0,1] with zero boundary condition, */
    /*     ... Call DNBAND  to find eigenvalues LAMBDA such that */
    /*                    A*x = LAMBDA*M*x. */

    /*     ... Eigenvalues with largest real parts are sought. */

    /*     ... Use mode 2 of DNAUPD . */

    /* \BeginLib */

    /* \Routines called: */
    /*     dnband   ARPACK banded eigenproblem solver. */
    /*     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully. */
    /*     dlaset   LAPACK routine to initialize a matrix to zero. */
    /*     daxpy    Level 1 BLAS that computes y <- alpha*x+y. */
    /*     dnrm2    Level 1 BLAS that computes the norm of a vector. */
    /*     dgbmv    Level 2 BLAS that computes the band matrix vector product. */

    /* \Author */
    /*     Richard Lehoucq */
    /*     Danny Sorensen */
    /*     Chao Yang */
    /*     Dept. of Computational & */
    /*     Applied Mathematics */
    /*     Rice University */
    /*     Houston, Texas */

    /* \SCCS Information: @(#) */
    /* FILE: nbdr3.F   SID: 2.5   DATE OF SID: 08/26/96   RELEASE: 2 */

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

    /* ----------------------------------------------- */
    /* The number N is the dimension of the matrix.  A */
    /* generalized eigenvalue problem is solved        */
    /* (BMAT = 'G').  NEV is the number of eigenvalues */
    /* to be approximated. The user can modify N, NEV, */
    /* NCV and WHICH to solve problems of different    */
    /* sizes, and to get different parts the spectrum. */
    /* However, the following conditions must be       */
    /* satisfied:                                      */
    /*                   N <= MAXN                     */
    /*                 NEV <= MAXNEV                   */
    /*           NEV + 2 <= NCV <= MAXNCV              */
    /* ----------------------------------------------- */

    n = 100;
    nev = 4;
    ncv = 10;
    if (n > 1000)
    {
        printf(" ERROR with _NBDR3: N is greater than MAXN \n");
        return ierr;
    }
    else if (nev > 25)
    {
        printf(" ERROR with _NBDR3: NEV is greater than MAXNEV \n");
        return ierr;
    }
    else if (ncv > 50)
    {
        printf(" ERROR with _NBDR3: NCV is greater than MAXNCV \n");
        return ierr;
    }
    bmat = "G";
    which = "LM";

    /* -------------------------------------------------- */
    /* The work array WORKL is used in DNAUPD  as          */
    /* workspace.  Its dimension LWORKL has to be set as  */
    /* illustrated below.  The parameter TOL determines   */
    /* the stopping criterion. If TOL<=0, machine machine */
    /* precision is used.  The number IDO is used for     */
    /* reverse communication and has to be set to 0 at    */
    /* the beginning.  Setting INFO=0 indicates that we   */
    /* using a randomly generated vector to start the     */
    /* the ARNOLDI process.                               */
    /* -------------------------------------------------- */

    /* Computing 2nd power */
    i__1 = ncv;
    lworkl = i__1 * i__1 * 3 + ncv * 6;
    info = 0;
    tol = 0.;
    ido = 0;

    /* ------------------------------------------------- */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 2 of DNAUPD  is used     */
    /* (IPARAM(7) = 2). All these options can be changed */
    /* by the user. For details, see the documentation   */
    /* in DNBAND .                                        */
    /* ------------------------------------------------- */

    mode = 2;
    maxitr = 300;

    iparam[2] = maxitr;
    iparam[6] = mode;

    /* ------------------------------------------ */
    /* Construct matrices A and M in LAPACK-style */
    /* banded form.                               */
    /* ------------------------------------------ */

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

    idiag = kl + ku + 1;
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
        a[idiag + j * 50 - 51] = 2. / h;
        m[idiag + j * 50 - 51] = h * 4.;
        /* L30: */
    }

    /* ----------------------------------- */
    /* First subdiagonal and superdiagonal */
    /* ----------------------------------- */

    isup = kl + ku;
    isub = kl + ku + 2;
    rho = 10.;
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
        a[isup + (j + 1) * 50 - 51] = -1. / h + rho / 2.;
        a[isub + j * 50 - 51] = -1. / h - rho / 2.;
        m[isup + (j + 1) * 50 - 51] = h * 1.;
        m[isub + j * 50 - 51] = h * 1.;
        /* L50: */
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
        printf(" _NBDR3 \n");
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

                dgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b97, &a[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b15, ax, &c__1);
                dgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b97, &m[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b15, mx, &c__1);
                d__1 = -d[j - 1];
                daxpy_(&n, &d__1, mx, &c__1, ax, &c__1);
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

                dgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b97, &a[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b15, ax, &c__1);
                dgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b97, &m[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b15, mx, &c__1);
                d__1 = -d[j - 1];
                daxpy_(&n, &d__1, mx, &c__1, ax, &c__1);
                dgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b97, &m[kl], &c__50, &v[(j + 1) * 1000 - 1000], &c__1, &c_b15, mx, &c__1);
                daxpy_(&n, &d[j + 49], mx, &c__1, ax, &c__1);
                d[j + 99] = dnrm2_(&n, ax, &c__1);
                dgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b97, &a[kl], &c__50, &v[(j + 1) * 1000 - 1000], &c__1, &c_b15, ax, &c__1);
                dgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b97, &m[kl], &c__50, &v[(j + 1) * 1000 - 1000], &c__1, &c_b15, mx, &c__1);
                d__1 = -d[j - 1];
                daxpy_(&n, &d__1, mx, &c__1, ax, &c__1);
                dgbmv_("Notranspose", &n, &n, &kl, &ku, &c_b97, &m[kl], &c__50, &v[j * 1000 - 1000], &c__1, &c_b15, mx, &c__1);
                d__1 = -d[j + 49];
                daxpy_(&n, &d__1, mx, &c__1, ax, &c__1);
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
        dmout_(&nconv, &c__3, d, &c__50, &c_n6,"Ritz values (Real,Imag) and relative residuals");
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
    return 0;
}
