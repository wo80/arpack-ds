/* EXAMPLES\SYM\dsdrv1.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

/* Table of constant values */

static a_int c__9 = 9;
static a_int c__1 = 1;
static a_int c__256 = 256;
static a_int c__3 = 3;
static a_int c__6 = 6;
static a_int c__2 = 2;
static a_int c__25 = 25;
static a_int c_n6 = -6;
static a_int c__5 = 5;
static double c_b138 = -1.;

/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    a_int i__1;
    double d__1;

    /* Builtin functions */
    a_int s_wsle(cilist *), do_lio(a_int *, a_int *, char *, ftnlen), e_wsle(void);

    /* Local variables */
    double d[50] /* was [25][2] */;
    a_int j, n;
    double v[6400] /* was [256][25] */;
    double ax[256];
    a_int nx, ido, ncv, nev;
    double tol;
    char bmat[1];
    a_int mode, info;
    a_bool rvec;
    a_int ierr;
    double sigma;
    char which[2];
    double resid[256];
    a_int nconv;
    double workd[768];
    a_int ipntr[11];
    double workl[825];
    a_int iparam[11];
    a_bool select[25];
    a_int ishfts, maxitr, lworkl;

    /* Fortran I/O blocks */
    static cilist io___5 = {0, 6, 0, 0, 0};
    static cilist io___6 = {0, 6, 0, 0, 0};
    static cilist io___7 = {0, 6, 0, 0, 0};
    static cilist io___23 = {0, 6, 0, 0, 0};
    static cilist io___24 = {0, 6, 0, 0, 0};
    static cilist io___25 = {0, 6, 0, 0, 0};
    static cilist io___26 = {0, 6, 0, 0, 0};
    static cilist io___32 = {0, 6, 0, 0, 0};
    static cilist io___33 = {0, 6, 0, 0, 0};
    static cilist io___34 = {0, 6, 0, 0, 0};
    static cilist io___35 = {0, 6, 0, 0, 0};
    static cilist io___39 = {0, 6, 0, 0, 0};
    static cilist io___40 = {0, 6, 0, 0, 0};
    static cilist io___41 = {0, 6, 0, 0, 0};
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
    static cilist io___55 = {0, 6, 0, 0, 0};
    static cilist io___56 = {0, 6, 0, 0, 0};
    static cilist io___57 = {0, 6, 0, 0, 0};

    /*     Simple program to illustrate the idea of reverse communication */
    /*     in regular mode for a standard symmetric eigenvalue problem. */

    /*     We implement example one of ex-sym.doc in SRC directory */

    /* \Example-1 */
    /*     ... Suppose we want to solve A*x = lambda*x in regular mode, */
    /*         where A is derived from the central difference discretization */
    /*         of the 2-dimensional Laplacian on the unit square [0,1]x[0,1] */
    /*         with zero Dirichlet boundary condition. */

    /*     ... OP = A  and  B = I. */

    /*     ... Assume "call av (n,x,y)" computes y = A*x. */

    /*     ... Use mode 1 of DSAUPD. */

    /* \BeginLib */

    /* \Routines called: */
    /*     dsaupd  ARPACK reverse communication interface routine. */
    /*     dseupd  ARPACK routine that returns Ritz values and (optionally) */
    /*             Ritz vectors. */
    /*     dnrm2   Level 1 BLAS that computes the norm of a vector. */
    /*     daxpy   Level 1 BLAS that computes y <- alpha*x+y. */
    /*     av      Matrix vector multiplication routine that computes A*x. */
    /*     tv      Matrix vector multiplication routine that computes T*x, */
    /*             where T is a tridiagonal matrix.  It is used in routine */
    /*             av. */

    /* \Author */
    /*     Richard Lehoucq */
    /*     Danny Sorensen */
    /*     Chao Yang */
    /*     Dept. of Computational & */
    /*     Applied Mathematics */
    /*     Rice University */
    /*     Houston, Texas */

    /* \SCCS Information: @(#) */
    /* FILE: sdrv1.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2 */

    /* \Remarks */
    /*     1. None */

    /* \EndLib */

    /* ----------------------------------------------------------------------- */

    /* --------------------------- */
    /* Define leading dimensions   */
    /* for all arrays.             */
    /* MAXN:   Maximum dimension   */
    /*         of the A allowed.   */
    /* MAXNEV: Maximum NEV allowed */
    /* MAXNCV: Maximum NCV allowed */
    /* --------------------------- */

    /* -------------------------------------------------- */
    /* The number NX is the number of interior points     */
    /* in the discretization of the 2-dimensional         */
    /* Laplacian on the unit square with zero Dirichlet   */
    /* boundary condition.  The number N(=NX*NX) is the   */
    /* dimension of the matrix.  A standard eigenvalue    */
    /* problem is solved (BMAT = 'I'). NEV is the number  */
    /* of eigenvalues to be approximated.  The user can   */
    /* modify NEV, NCV, WHICH to solve problems of        */
    /* different sizes, and to get different parts of the */
    /* spectrum.  However, The following conditions must  */
    /* be satisfied:                                      */
    /*                   N <= MAXN,                       */
    /*                 NEV <= MAXNEV,                     */
    /*             NEV + 1 <= NCV <= MAXNCV               */
    /* -------------------------------------------------- */

    nx = 10;
    n = nx * nx;
    nev = 4;
    ncv = 10;
    if (n > 256)
    {
        s_wsle(&io___5);
        do_lio(&c__9, &c__1, " ERROR with _SDRV1: N is greater than MAXN ", (ftnlen)43);
        e_wsle();
        goto L9000;
    }
    else if (nev > 10)
    {
        s_wsle(&io___6);
        do_lio(&c__9, &c__1, " ERROR with _SDRV1: NEV is greater than MAXNEV ", (ftnlen)47);
        e_wsle();
        goto L9000;
    }
    else if (ncv > 25)
    {
        s_wsle(&io___7);
        do_lio(&c__9, &c__1, " ERROR with _SDRV1: NCV is greater than MAXNCV ", (ftnlen)47);
        e_wsle();
        goto L9000;
    }
    *bmat = 'I';
    strcpy(which, "SM");

    /* ------------------------------------------------ */
    /* The work array WORKL is used in DSAUPD as        */
    /* workspace.  Its dimension LWORKL is set as       */
    /* illustrated below.  The parameter TOL determines */
    /* the stopping criterion.  If TOL<=0, machine      */
    /* precision is used.  The variable IDO is used for */
    /* reverse communication and is initially set to 0. */
    /* Setting INFO=0 indicates that a random vector is */
    /* generated in DSAUPD to start the Arnoldi         */
    /* iteration.                                       */
    /* ------------------------------------------------ */

    lworkl = ncv * (ncv + 8);
    tol = 0.;
    info = 0;
    ido = 0;

    /* ------------------------------------------------- */
    /* This program uses exact shifts with respect to    */
    /* the current Hessenberg matrix (IPARAM(1) = 1).    */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 1 of DSAUPD is used     */
    /* (IPARAM(7) = 1).  All these options may be        */
    /* changed by the user. For details, see the         */
    /* documentation in DSAUPD.                          */
    /* ------------------------------------------------- */

    ishfts = 1;
    maxitr = 300;
    mode = 1;

    iparam[0] = ishfts;
    iparam[2] = maxitr;
    iparam[6] = mode;

    /* ----------------------------------------- */
    /* M A I N   L O O P (Reverse communication) */
    /* ----------------------------------------- */

L10:

    /* ------------------------------------------- */
    /* Repeatedly call the routine DSAUPD and take */
    /* actions indicated by parameter IDO until    */
    /* either convergence is indicated or maxitr   */
    /* has been exceeded.                          */
    /* ------------------------------------------- */

    dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, &info);

    if (ido == -1 || ido == 1)
    {

        /* ------------------------------------ */
        /* Perform matrix vector multiplication */
        /*              y <--- OP*x             */
        /* The user should supply his/her own   */
        /* matrix vector multiplication routine */
        /* here that takes workd(ipntr(1)) as   */
        /* the input, and return the result to  */
        /* workd(ipntr(2)).                     */
        /* ------------------------------------ */

        av_(&nx, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

        /* --------------------------------------- */
        /* L O O P   B A C K to call DSAUPD again. */
        /* --------------------------------------- */

        goto L10;
    }

    /* -------------------------------------- */
    /* Either we have convergence or there is */
    /* an error.                              */
    /* -------------------------------------- */

    if (info < 0)
    {

        /* ------------------------ */
        /* Error message. Check the */
        /* documentation in DSAUPD. */
        /* ------------------------ */

        s_wsle(&io___23);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
        s_wsle(&io___24);
        do_lio(&c__9, &c__1, " Error with _saupd, info = ", (ftnlen)27);
        do_lio(&c__3, &c__1, (char *)&info, (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___25);
        do_lio(&c__9, &c__1, " Check documentation in _saupd ", (ftnlen)31);
        e_wsle();
        s_wsle(&io___26);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
    }
    else
    {

        /* ----------------------------------------- */
        /* No fatal errors occurred.                 */
        /* Post-Process using DSEUPD.                */
        /*                                           */
        /* Computed eigenvalues may be extracted.    */
        /*                                           */
        /* Eigenvectors may also be computed now if  */
        /* desired.  (indicated by rvec = .true.)    */
        /* ----------------------------------------- */

        rvec = TRUE_;

        dseupd_(&rvec, "All", select, d, v, &c__256, &sigma, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, &ierr);
        /* -------------------------------------------- */
        /* Eigenvalues are returned in the first column */
        /* of the two dimensional array D and the       */
        /* corresponding eigenvectors are returned in   */
        /* the first NEV columns of the two dimensional */
        /* array V if requested.  Otherwise, an         */
        /* orthogonal basis for the invariant subspace  */
        /* corresponding to the eigenvalues in D is     */
        /* returned in V.                               */
        /* -------------------------------------------- */

        if (ierr != 0)
        {

            /* ---------------------------------- */
            /* Error condition:                   */
            /* Check the documentation of DSEUPD. */
            /* ---------------------------------- */

            s_wsle(&io___32);
            do_lio(&c__9, &c__1, " ", (ftnlen)1);
            e_wsle();
            s_wsle(&io___33);
            do_lio(&c__9, &c__1, " Error with _seupd, info = ", (ftnlen)27);
            do_lio(&c__3, &c__1, (char *)&ierr, (ftnlen)sizeof(a_int));
            e_wsle();
            s_wsle(&io___34);
            do_lio(&c__9, &c__1, " Check the documentation of _seupd. ", (ftnlen)36);
            e_wsle();
            s_wsle(&io___35);
            do_lio(&c__9, &c__1, " ", (ftnlen)1);
            e_wsle();
        }
        else
        {

            nconv = iparam[4];
            i__1 = nconv;
            for (j = 1; j <= i__1; ++j)
            {

                /* ------------------------- */
                /* Compute the residual norm */
                /*                           */
                /*   ||  A*x - lambda*x ||   */
                /*                           */
                /* for the NCONV accurately  */
                /* computed eigenvalues and  */
                /* eigenvectors.  (iparam(5) */
                /* indicates how many are    */
                /* accurate to the requested */
                /* tolerance)                */
                /* ------------------------- */

                av_(&nx, &v[(j << 8) - 256], ax);
                d__1 = -d[j - 1];
                daxpy_(&n, &d__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
                d[j + 24] = dnrm2_(&n, ax, &c__1);
                d[j + 24] /= (d__1 = d[j - 1], abs(d__1));

                /* L20: */
            }

            /* ----------------------------- */
            /* Display computed residuals    */
            /* ----------------------------- */

            dmout_(&c__6, &nconv, &c__2, d, &c__25, &c_n6,"Ritz values and relative residuals");
        }

        /* ---------------------------------------- */
        /* Print additional convergence information */
        /* ---------------------------------------- */

        if (info == 1)
        {
            s_wsle(&io___39);
            do_lio(&c__9, &c__1, " ", (ftnlen)1);
            e_wsle();
            s_wsle(&io___40);
            do_lio(&c__9, &c__1, " Maximum number of iterations reached.", (ftnlen)38);
            e_wsle();
            s_wsle(&io___41);
            do_lio(&c__9, &c__1, " ", (ftnlen)1);
            e_wsle();
        }
        else if (info == 3)
        {
            s_wsle(&io___42);
            do_lio(&c__9, &c__1, " ", (ftnlen)1);
            e_wsle();
            s_wsle(&io___43);
            do_lio(&c__9, &c__1, " No shifts could be applied during implicit", (ftnlen)43);
            do_lio(&c__9, &c__1, " Arnoldi update, try increasing NCV.", (ftnlen)36);
            e_wsle();
            s_wsle(&io___44);
            do_lio(&c__9, &c__1, " ", (ftnlen)1);
            e_wsle();
        }

        s_wsle(&io___45);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
        s_wsle(&io___46);
        do_lio(&c__9, &c__1, " _SDRV1 ", (ftnlen)8);
        e_wsle();
        s_wsle(&io___47);
        do_lio(&c__9, &c__1, " ====== ", (ftnlen)8);
        e_wsle();
        s_wsle(&io___48);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
        s_wsle(&io___49);
        do_lio(&c__9, &c__1, " Size of the matrix is ", (ftnlen)23);
        do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___50);
        do_lio(&c__9, &c__1, " The number of Ritz values requested is ", (ftnlen)40);
        do_lio(&c__3, &c__1, (char *)&nev, (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___51);
        do_lio(&c__9, &c__1, " The number of Arnoldi vectors generated", (ftnlen)40);
        do_lio(&c__9, &c__1, " (NCV) is ", (ftnlen)10);
        do_lio(&c__3, &c__1, (char *)&ncv, (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___52);
        do_lio(&c__9, &c__1, " What portion of the spectrum: ", (ftnlen)31);
        do_lio(&c__9, &c__1, which, (ftnlen)2);
        e_wsle();
        s_wsle(&io___53);
        do_lio(&c__9, &c__1, " The number of converged Ritz values is ", (ftnlen)40);
        do_lio(&c__3, &c__1, (char *)&nconv, (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___54);
        do_lio(&c__9, &c__1, " The number of Implicit Arnoldi update", (ftnlen)38);
        do_lio(&c__9, &c__1, " iterations taken is ", (ftnlen)21);
        do_lio(&c__3, &c__1, (char *)&iparam[2], (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___55);
        do_lio(&c__9, &c__1, " The number of OP*x is ", (ftnlen)23);
        do_lio(&c__3, &c__1, (char *)&iparam[8], (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___56);
        do_lio(&c__9, &c__1, " The convergence criterion is ", (ftnlen)30);
        do_lio(&c__5, &c__1, (char *)&tol, (ftnlen)sizeof(double));
        e_wsle();
        s_wsle(&io___57);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
    }

    /* ------------------------- */
    /* Done with program dsdrv1. */
    /* ------------------------- */

L9000:

    return 0;
} /* MAIN__ */

/* ------------------------------------------------------------------ */
/*     matrix vector subroutine */

/*     The matrix used is the 2 dimensional discrete Laplacian on unit */
/*     square with zero Dirichlet boundary condition. */

/*     Computes w <--- OP*v, where OP is the nx*nx by nx*nx block */
/*     tridiagonal matrix */

/*                  | T -I          | */
/*                  |-I  T -I       | */
/*             OP = |   -I  T       | */
/*                  |        ...  -I| */
/*                  |           -I T| */

/*     The subroutine TV is called to computed y<---T*x. */

int av_(a_int *nx, double *v, double *w)
{
    /* System generated locals */
    a_int i__1;
    double d__1;

    /* Local variables */
    a_int j;
    double h2;
    a_int n2, lo;

    /* Parameter adjustments */
    --w;
    --v;

    tv_(nx, &v[1], &w[1]);
    daxpy_(nx, &c_b138, &v[*nx + 1], &c__1, &w[1], &c__1);

    i__1 = *nx - 1;
    for (j = 2; j <= i__1; ++j)
    {
        lo = (j - 1) * *nx;
        tv_(nx, &v[lo + 1], &w[lo + 1]);
        daxpy_(nx, &c_b138, &v[lo - *nx + 1], &c__1, &w[lo + 1], &c__1);
        daxpy_(nx, &c_b138, &v[lo + *nx + 1], &c__1, &w[lo + 1], &c__1);
        /* L10: */
    }

    lo = (*nx - 1) * *nx;
    tv_(nx, &v[lo + 1], &w[lo + 1]);
    daxpy_(nx, &c_b138, &v[lo - *nx + 1], &c__1, &w[lo + 1], &c__1);

    /*     Scale the vector w by (1/h^2), where h is the mesh size */

    n2 = *nx * *nx;
    h2 = 1. / (double)((*nx + 1) * (*nx + 1));
    d__1 = 1. / h2;
    dscal_(&n2, &d__1, &w[1], &c__1);
    return 0;
} /* av_ */

/* ------------------------------------------------------------------- */
int tv_(a_int *nx, double *x, double *y)
{
    /* System generated locals */
    a_int i__1;

    /* Local variables */
    a_int j;
    double dd, dl, du;

    /*     Compute the matrix vector multiplication y<---T*x */
    /*     where T is a nx by nx tridiagonal matrix with DD on the */
    /*     diagonal, DL on the subdiagonal, and DU on the superdiagonal. */

    /* Parameter adjustments */
    --y;
    --x;

    dd = 4.;
    dl = -1.;
    du = -1.;

    y[1] = dd * x[1] + du * x[2];
    i__1 = *nx - 1;
    for (j = 2; j <= i__1; ++j)
    {
        y[j] = dl * x[j - 1] + dd * x[j] + du * x[j + 1];
        /* L10: */
    }
    y[*nx] = dl * x[*nx - 1] + dd * x[*nx];
    return 0;
} /* tv_ */

/* Main program alias */ int dsdrv1_()
{
    MAIN__();
    return 0;
}
