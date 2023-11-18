/* EXAMPLES\SYM\dsdrv3.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__256 = 256;
static integer c__3 = 3;
static integer c__6 = 6;
static integer c__2 = 2;
static integer c__25 = 25;
static integer c_n6 = -6;
static integer c__5 = 5;

/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), e_wsle(void);
    int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    doublereal d__[50] /* was [25][2] */, h__;
    integer j, n;
    doublereal v[6400] /* was [256][25] */, r1, r2, ad[256];
    extern int av_(integer *, doublereal *, doublereal *);
    doublereal ax[256];
    extern int mv_(integer *, doublereal *, doublereal *);
    doublereal mx[256], adl[256], adu[256];
    integer ido, ncv, nev;
    doublereal tol, adu2[256];
    char bmat[1];
    integer mode, info;
    logical rvec;
    integer ierr, ipiv[256];
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    doublereal sigma;
    char which[2];
    doublereal resid[256];
    extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    integer nconv;
    extern int daxpy_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *);
    doublereal workd[768];
    extern int dmout_(integer *, integer *, integer *, doublereal *, integer *, integer *, char *, ftnlen);
    integer ipntr[11];
    doublereal workl[825];
    integer iparam[11];
    logical select[25];
    extern int dsaupd_(integer *, char *, integer *, char *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, integer *, integer *, doublereal *, doublereal *, integer *, integer *, ftnlen, ftnlen), dseupd_(logical *, char *, logical *, doublereal *, doublereal *, integer *, doublereal *, char *, integer *, char *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, integer *, integer *, doublereal *, doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen),
        dgttrf_(integer *, doublereal *, doublereal *, doublereal *, doublereal *, integer *, integer *);
    integer ishfts, maxitr;
    extern int dgttrs_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, integer *, integer *, ftnlen);
    integer lworkl;

    /* Fortran I/O blocks */
    static cilist io___4 = {0, 6, 0, 0, 0};
    static cilist io___5 = {0, 6, 0, 0, 0};
    static cilist io___6 = {0, 6, 0, 0, 0};
    static cilist io___27 = {0, 6, 0, 0, 0};
    static cilist io___28 = {0, 6, 0, 0, 0};
    static cilist io___29 = {0, 6, 0, 0, 0};
    static cilist io___35 = {0, 6, 0, 0, 0};
    static cilist io___36 = {0, 6, 0, 0, 0};
    static cilist io___37 = {0, 6, 0, 0, 0};
    static cilist io___38 = {0, 6, 0, 0, 0};
    static cilist io___39 = {0, 6, 0, 0, 0};
    static cilist io___40 = {0, 6, 0, 0, 0};
    static cilist io___41 = {0, 6, 0, 0, 0};
    static cilist io___46 = {0, 6, 0, 0, 0};
    static cilist io___47 = {0, 6, 0, 0, 0};
    static cilist io___48 = {0, 6, 0, 0, 0};
    static cilist io___49 = {0, 6, 0, 0, 0};
    static cilist io___53 = {0, 6, 0, 0, 0};
    static cilist io___54 = {0, 6, 0, 0, 0};
    static cilist io___55 = {0, 6, 0, 0, 0};
    static cilist io___56 = {0, 6, 0, 0, 0};
    static cilist io___57 = {0, 6, 0, 0, 0};
    static cilist io___58 = {0, 6, 0, 0, 0};
    static cilist io___59 = {0, 6, 0, 0, 0};
    static cilist io___60 = {0, 6, 0, 0, 0};
    static cilist io___61 = {0, 6, 0, 0, 0};
    static cilist io___62 = {0, 6, 0, 0, 0};
    static cilist io___63 = {0, 6, 0, 0, 0};
    static cilist io___64 = {0, 6, 0, 0, 0};
    static cilist io___65 = {0, 6, 0, 0, 0};
    static cilist io___66 = {0, 6, 0, 0, 0};
    static cilist io___67 = {0, 6, 0, 0, 0};
    static cilist io___68 = {0, 6, 0, 0, 0};
    static cilist io___69 = {0, 6, 0, 0, 0};
    static cilist io___70 = {0, 6, 0, 0, 0};
    static cilist io___71 = {0, 6, 0, 0, 0};

    /*     Program to illustrate the idea of reverse communication in */
    /*     inverse mode for a generalized symmetric eigenvalue problem. */
    /*     The following program uses the two LAPACK subroutines dgttrf .f */
    /*     and dgttrs .f to factor and solve a tridiagonal system of equations. */

    /*     We implement example three of ex-sym.doc in DOCUMENTS directory */

    /* \Example-3 */
    /*     ... Suppose we want to solve A*x = lambda*M*x in inverse mode, */
    /*         where A and M are obtained by the finite element of the */
    /*         1-dimensional discrete Laplacian */
    /*                             d^2u / dx^2 */
    /*         on the interval [0,1] with zero Dirichlet boundary condition */
    /*         using piecewise linear elements. */

    /*     ... OP = inv[M]*A  and  B = M. */

    /*     ... Use mode 2 of DSAUPD . */

    /* \BeginLib */

    /* \Routines called: */
    /*     dsaupd   ARPACK reverse communication interface routine. */
    /*     dseupd   ARPACK routine that returns Ritz values and (optionally) */
    /*             Ritz vectors. */
    /*     dgttrf   LAPACK tridiagonal factorization routine. */
    /*     dgttrs   LAPACK tridiagonal solve routine. */
    /*     daxpy    Level 1 BLAS that computes y <- alpha*x+y. */
    /*     dscal    Level 1 BLAS that scales a vector by a scalar. */
    /*     dcopy    Level 1 BLAS that copies one vector to another. */
    /*     dnrm2    Level 1 BLAS that computes the norm of a vector. */
    /*     av      Matrix vector multiplication routine that computes A*x. */
    /*     mv      Matrix vector multiplication routine that computes M*x. */

    /* \Author */
    /*     Richard Lehoucq */
    /*     Danny Sorensen */
    /*     Chao Yang */
    /*     Dept. of Computational & */
    /*     Applied Mathematics */
    /*     Rice University */
    /*     Houston, Texas */

    /* \SCCS Information: @(#) */
    /* FILE: sdrv3.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2 */

    /* \Remarks */
    /*     1. None */

    /* \EndLib */
    /* -------------------------------------------------------------------------- */

    /*     %-----------------------------% */
    /*     | Define leading dimensions   | */
    /*     | for all arrays.             | */
    /*     | MAXN:   Maximum dimension   | */
    /*     |         of the A allowed.   | */
    /*     | MAXNEV: Maximum NEV allowed | */
    /*     | MAXNCV: Maximum NCV allowed | */
    /*     %-----------------------------% */

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
    /*     | Executable statements | */
    /*     %-----------------------% */

    /*     %----------------------------------------------------% */
    /*     | The number N is the dimension of the matrix. A     | */
    /*     | generalized eigenvalue problem is solved (BMAT =   | */
    /*     | 'G'.) NEV is the number of eigenvalues to be       | */
    /*     | approximated.  The user can modify NEV, NCV, WHICH | */
    /*     | to solve problems of different sizes, and to get   | */
    /*     | different parts of the spectrum.  However, The     | */
    /*     | following conditions must be satisfied:            | */
    /*     |                     N <= MAXN,                     | */
    /*     |                   NEV <= MAXNEV,                   | */
    /*     |               NEV + 1 <= NCV <= MAXNCV             | */
    /*     %----------------------------------------------------% */

    n = 100;
    nev = 4;
    ncv = 10;
    if (n > 256)
    {
        s_wsle(&io___4);
        do_lio(&c__9, &c__1, " ERROR with _SDRV3: N is greater than MAXN ", (ftnlen)43);
        e_wsle();
        goto L9000;
    }
    else if (nev > 10)
    {
        s_wsle(&io___5);
        do_lio(&c__9, &c__1, " ERROR with _SDRV3: NEV is greater than MAXNEV ", (ftnlen)47);
        e_wsle();
        goto L9000;
    }
    else if (ncv > 25)
    {
        s_wsle(&io___6);
        do_lio(&c__9, &c__1, " ERROR with _SDRV3: NCV is greater than MAXNCV ", (ftnlen)47);
        e_wsle();
        goto L9000;
    }
    *(unsigned char *)bmat = 'G';
    s_copy(which, "LM", (ftnlen)2, (ftnlen)2);

    /*     %--------------------------------------------------% */
    /*     | The work array WORKL is used in DSAUPD  as        | */
    /*     | workspace.  Its dimension LWORKL is set as       | */
    /*     | illustrated below.  The parameter TOL determines | */
    /*     | the stopping criterion.  If TOL<=0, machine      | */
    /*     | precision is used.  The variable IDO is used for | */
    /*     | reverse communication and is initially set to 0. | */
    /*     | Setting INFO=0 indicates that a random vector is | */
    /*     | generated in DSAUPD  to start the Arnoldi         | */
    /*     | iteration.                                       | */
    /*     %--------------------------------------------------% */

    lworkl = ncv * (ncv + 8);
    tol = 0.;
    ido = 0;
    info = 0;

    /*     %---------------------------------------------------% */
    /*     | This program uses exact shifts with respect to    | */
    /*     | the current Hessenberg matrix (IPARAM(1) = 1).    | */
    /*     | IPARAM(3) specifies the maximum number of Arnoldi | */
    /*     | iterations allowed.  Mode 2 of DSAUPD  is used     | */
    /*     | (IPARAM(7) = 2).  All these options may be        | */
    /*     | changed by the user. For details, see the         | */
    /*     | documentation in DSAUPD .                          | */
    /*     %---------------------------------------------------% */

    ishfts = 1;
    maxitr = 300;
    mode = 2;

    iparam[0] = ishfts;
    iparam[2] = maxitr;
    iparam[6] = mode;

    /*     %------------------------------------------------% */
    /*     | Call LAPACK routine to factor the mass matrix. | */
    /*     | The mass matrix is the tridiagonal matrix      | */
    /*     | arising from using piecewise linear finite     | */
    /*     | elements on the interval [0, 1].               | */
    /*     %------------------------------------------------% */

    h__ = 1. / (doublereal)(n + 1);

    r1 = h__ * .66666666666666663;
    r2 = h__ * .16666666666666666;
    i__1 = n;
    for (j = 1; j <= i__1; ++j)
    {
        ad[j - 1] = r1;
        adl[j - 1] = r2;
        /* L20: */
    }
    dcopy_(&n, adl, &c__1, adu, &c__1);
    dgttrf_(&n, adl, ad, adu, adu2, ipiv, &ierr);
    if (ierr != 0)
    {
        s_wsle(&io___27);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
        s_wsle(&io___28);
        do_lio(&c__9, &c__1, " Error with _gttrf in _SDRV3. ", (ftnlen)30);
        e_wsle();
        s_wsle(&io___29);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
        goto L9000;
    }

    /*     %-------------------------------------------% */
    /*     | M A I N   L O O P (Reverse communication) | */
    /*     %-------------------------------------------% */

L10:

    /*        %---------------------------------------------% */
    /*        | Repeatedly call the routine DSAUPD  and take | */
    /*        | actions indicated by parameter IDO until    | */
    /*        | either convergence is indicated or maxitr   | */
    /*        | has been exceeded.                          | */
    /*        %---------------------------------------------% */

    dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, &info, (ftnlen)1, (ftnlen)2);

    if (ido == -1 || ido == 1)
    {

        /*           %--------------------------------------% */
        /*           | Perform  y <--- OP*x = inv[M]*A*x.   | */
        /*           | The user should supply his/her own   | */
        /*           | matrix vector multiplication (A*x)   | */
        /*           | routine and a linear system solver   | */
        /*           | here.  The matrix vector             | */
        /*           | multiplication routine takes         | */
        /*           | workd(ipntr(1)) as the input vector. | */
        /*           | The final result is returned to      | */
        /*           | workd(ipntr(2)). The result of A*x   | */
        /*           | overwrites workd(ipntr(1)).          | */
        /*           %--------------------------------------% */

        av_(&n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
        dcopy_(&n, &workd[ipntr[1] - 1], &c__1, &workd[ipntr[0] - 1], &c__1);
        dgttrs_("Notranspose", &n, &c__1, adl, ad, adu, adu2, ipiv, &workd[ipntr[1] - 1], &n, &ierr, (ftnlen)11);
        if (ierr != 0)
        {
            s_wsle(&io___35);
            do_lio(&c__9, &c__1, " ", (ftnlen)1);
            e_wsle();
            s_wsle(&io___36);
            do_lio(&c__9, &c__1, " Error with _gttrs in _SDRV3.", (ftnlen)29);
            e_wsle();
            s_wsle(&io___37);
            do_lio(&c__9, &c__1, " ", (ftnlen)1);
            e_wsle();
            goto L9000;
        }

        /*           %-----------------------------------------% */
        /*           | L O O P   B A C K to call DSAUPD  again. | */
        /*           %-----------------------------------------% */

        goto L10;
    }
    else if (ido == 2)
    {

        /*           %-----------------------------------------% */
        /*           |         Perform  y <--- M*x.            | */
        /*           | Need the matrix vector multiplication   | */
        /*           | routine here that takes workd(ipntr(1)) | */
        /*           | as the input and returns the result to  | */
        /*           | workd(ipntr(2)).                        | */
        /*           %-----------------------------------------% */

        mv_(&n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

        /*           %-----------------------------------------% */
        /*           | L O O P   B A C K to call DSAUPD  again. | */
        /*           %-----------------------------------------% */

        goto L10;
    }

    /*     %-----------------------------------------% */
    /*     | Either we have convergence, or there is | */
    /*     | an error.                               | */
    /*     %-----------------------------------------% */

    if (info < 0)
    {

        /*        %--------------------------% */
        /*        | Error message, check the | */
        /*        | documentation in DSAUPD   | */
        /*        %--------------------------% */

        s_wsle(&io___38);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
        s_wsle(&io___39);
        do_lio(&c__9, &c__1, " Error with _saupd, info = ", (ftnlen)27);
        do_lio(&c__3, &c__1, (char *)&info, (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___40);
        do_lio(&c__9, &c__1, " Check the documentation of _saupd ", (ftnlen)35);
        e_wsle();
        s_wsle(&io___41);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
    }
    else
    {

        /*        %-------------------------------------------% */
        /*        | No fatal errors occurred.                 | */
        /*        | Post-Process using DSEUPD .                | */
        /*        |                                           | */
        /*        | Computed eigenvalues may be extracted.    | */
        /*        |                                           | */
        /*        | Eigenvectors may also be computed now if  | */
        /*        | desired.  (indicated by rvec = .true.)    | */
        /*        %-------------------------------------------% */

        rvec = TRUE_;

        dseupd_(&rvec, "All", select, d__, v, &c__256, &sigma, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, &ierr, (ftnlen)3, (ftnlen)1, (ftnlen)2);

        /*        %----------------------------------------------% */
        /*        | Eigenvalues are returned in the first column | */
        /*        | of the two dimensional array D and the       | */
        /*        | corresponding eigenvectors are returned in   | */
        /*        | the first NEV columns of the two dimensional | */
        /*        | array V if requested.  Otherwise, an         | */
        /*        | orthogonal basis for the invariant subspace  | */
        /*        | corresponding to the eigenvalues in D is     | */
        /*        | returned in V.                               | */
        /*        %----------------------------------------------% */

        if (ierr != 0)
        {

            /*           %------------------------------------% */
            /*           | Error condition:                   | */
            /*           | Check the documentation of DSEUPD . | */
            /*           %------------------------------------% */

            s_wsle(&io___46);
            do_lio(&c__9, &c__1, " ", (ftnlen)1);
            e_wsle();
            s_wsle(&io___47);
            do_lio(&c__9, &c__1, " Error with _seupd, info = ", (ftnlen)27);
            do_lio(&c__3, &c__1, (char *)&ierr, (ftnlen)sizeof(integer));
            e_wsle();
            s_wsle(&io___48);
            do_lio(&c__9, &c__1, " Check the documentation of _seupd", (ftnlen)34);
            e_wsle();
            s_wsle(&io___49);
            do_lio(&c__9, &c__1, " ", (ftnlen)1);
            e_wsle();
        }
        else
        {

            nconv = iparam[4];
            i__1 = nconv;
            for (j = 1; j <= i__1; ++j)
            {

                /*              %---------------------------% */
                /*              | Compute the residual norm | */
                /*              |                           | */
                /*              |  ||  A*x - lambda*M*x ||  | */
                /*              |                           | */
                /*              | for the NCONV accurately  | */
                /*              | computed eigenvalues and  | */
                /*              | eigenvectors.  (iparam(5) | */
                /*              | indicates how many are    | */
                /*              | accurate to the requested | */
                /*              | tolerance)                | */
                /*              %---------------------------% */

                av_(&n, &v[(j << 8) - 256], ax);
                mv_(&n, &v[(j << 8) - 256], mx);
                d__1 = -d__[j - 1];
                daxpy_(&n, &d__1, mx, &c__1, ax, &c__1);
                d__[j + 24] = dnrm2_(&n, ax, &c__1);
                d__[j + 24] /= (d__1 = d__[j - 1], abs(d__1));

                /* L30: */
            }

            /*           %-----------------------------% */
            /*           | Display computed residuals. | */
            /*           %-----------------------------% */

            dmout_(&c__6, &nconv, &c__2, d__, &c__25, &c_n6,
                   "Ritz values an"
                   "d relative residuals",
                   (ftnlen)34);
        }

        /*        %------------------------------------------% */
        /*        | Print additional convergence information | */
        /*        %------------------------------------------% */

        if (info == 1)
        {
            s_wsle(&io___53);
            do_lio(&c__9, &c__1, " ", (ftnlen)1);
            e_wsle();
            s_wsle(&io___54);
            do_lio(&c__9, &c__1, " Maximum number of iterations reached.", (ftnlen)38);
            e_wsle();
            s_wsle(&io___55);
            do_lio(&c__9, &c__1, " ", (ftnlen)1);
            e_wsle();
        }
        else if (info == 3)
        {
            s_wsle(&io___56);
            do_lio(&c__9, &c__1, " ", (ftnlen)1);
            e_wsle();
            s_wsle(&io___57);
            do_lio(&c__9, &c__1, " No shifts could be applied during implicit", (ftnlen)43);
            do_lio(&c__9, &c__1, " Arnoldi update, try increasing NCV.", (ftnlen)36);
            e_wsle();
            s_wsle(&io___58);
            do_lio(&c__9, &c__1, " ", (ftnlen)1);
            e_wsle();
        }

        s_wsle(&io___59);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
        s_wsle(&io___60);
        do_lio(&c__9, &c__1, " _SDRV3 ", (ftnlen)8);
        e_wsle();
        s_wsle(&io___61);
        do_lio(&c__9, &c__1, " ====== ", (ftnlen)8);
        e_wsle();
        s_wsle(&io___62);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
        s_wsle(&io___63);
        do_lio(&c__9, &c__1, " Size of the matrix is ", (ftnlen)23);
        do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___64);
        do_lio(&c__9, &c__1, " The number of Ritz values requested is ", (ftnlen)40);
        do_lio(&c__3, &c__1, (char *)&nev, (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___65);
        do_lio(&c__9, &c__1, " The number of Arnoldi vectors generated", (ftnlen)40);
        do_lio(&c__9, &c__1, " (NCV) is ", (ftnlen)10);
        do_lio(&c__3, &c__1, (char *)&ncv, (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___66);
        do_lio(&c__9, &c__1, " What portion of the spectrum: ", (ftnlen)31);
        do_lio(&c__9, &c__1, which, (ftnlen)2);
        e_wsle();
        s_wsle(&io___67);
        do_lio(&c__9, &c__1, " The number of converged Ritz values is ", (ftnlen)40);
        do_lio(&c__3, &c__1, (char *)&nconv, (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___68);
        do_lio(&c__9, &c__1, " The number of Implicit Arnoldi update", (ftnlen)38);
        do_lio(&c__9, &c__1, " iterations taken is ", (ftnlen)21);
        do_lio(&c__3, &c__1, (char *)&iparam[2], (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___69);
        do_lio(&c__9, &c__1, " The number of OP*x is ", (ftnlen)23);
        do_lio(&c__3, &c__1, (char *)&iparam[8], (ftnlen)sizeof(integer));
        e_wsle();
        s_wsle(&io___70);
        do_lio(&c__9, &c__1, " The convergence criterion is ", (ftnlen)30);
        do_lio(&c__5, &c__1, (char *)&tol, (ftnlen)sizeof(doublereal));
        e_wsle();
        s_wsle(&io___71);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
    }

    /*     %---------------------------% */
    /*     | Done with program dsdrv3 . | */
    /*     %---------------------------% */

L9000:

    return 0;
} /* MAIN__ */

/* ------------------------------------------------------------------------ */
/*     Matrix vector subroutine */
/*     where the matrix is the 1 dimensional mass matrix */
/*     on the interval [0,1]. */

int mv_(integer *n, doublereal *v, doublereal *w)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    doublereal h__;
    integer j;
    extern int dscal_(integer *, doublereal *, doublereal *, integer *);

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    w[1] = v[1] * 4. + v[2];
    i__1 = *n - 1;
    for (j = 2; j <= i__1; ++j)
    {
        w[j] = v[j - 1] + v[j] * 4. + v[j + 1];
        /* L100: */
    }
    j = *n;
    w[j] = v[j - 1] + v[j] * 4.;

    /*     Scale the vector w by h. */

    h__ = 1. / ((doublereal)(*n + 1) * 6.);
    dscal_(n, &h__, &w[1], &c__1);
    return 0;
} /* mv_ */

/* -------------------------------------------------------------------- */
/*     matrix vector subroutine */

/*     The matrix used is the stiffness matrix obtained from the finite */
/*     element discretization of the 1-dimensional discrete Laplacian */
/*     on the interval [0,1] with zero Dirichlet boundary condition using */
/*     piecewise linear elements. */

int av_(integer *n, doublereal *v, doublereal *w)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    doublereal h__;
    integer j;
    extern int dscal_(integer *, doublereal *, doublereal *, integer *);

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    w[1] = v[1] * 2. - v[2];
    i__1 = *n - 1;
    for (j = 2; j <= i__1; ++j)
    {
        w[j] = -v[j - 1] + v[j] * 2. - v[j + 1];
        /* L100: */
    }
    j = *n;
    w[j] = -v[j - 1] + v[j] * 2.;

    /*     Scale the vector w by (1 / h). */

    h__ = 1. / (doublereal)(*n + 1);
    d__1 = 1. / h__;
    dscal_(n, &d__1, &w[1], &c__1);
    return 0;
} /* av_ */

/* Main program alias */ int dsdrv3_()
{
    MAIN__();
    return 0;
}
