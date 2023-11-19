/* EXAMPLES\NONSYM\sndrv1.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

/* Table of constant values */

static a_int c__9 = 9;
static a_int c__1 = 1;
static a_int c__256 = 256;
static a_int c__3 = 3;
static a_int c__6 = 6;
static a_int c__30 = 30;
static a_int c_n6 = -6;
static a_int c__4 = 4;

/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    a_int i__1;
    float r__1;

    /* Builtin functions */
    a_int s_wsle(cilist *), do_lio(a_int *, a_int *, char *, ftnlen), e_wsle(void);
    int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    float d[90] /* was [30][3] */;
    a_int j, n;
    float v[7680] /* was [256][30] */;
    float ax[256];
    a_int nx, ido, ncv, nev;
    float tol;
    char bmat[1];
    a_int mode, info;
    a_bool rvec;
    a_int ierr;
    char which[2];
    float resid[256];
    a_int nconv;
    float workd[768];
    a_bool first;
    a_int ipntr[14];
    float workl[2880];
    a_int iparam[11];
    float sigmai;
    a_bool select[30];
    float sigmar;
    a_int ishfts, maxitr, lworkl;
    float workev[90];

    /* Fortran I/O blocks */
    static cilist io___5 = {0, 6, 0, 0, 0};
    static cilist io___6 = {0, 6, 0, 0, 0};
    static cilist io___7 = {0, 6, 0, 0, 0};
    static cilist io___23 = {0, 6, 0, 0, 0};
    static cilist io___24 = {0, 6, 0, 0, 0};
    static cilist io___25 = {0, 6, 0, 0, 0};
    static cilist io___26 = {0, 6, 0, 0, 0};
    static cilist io___34 = {0, 6, 0, 0, 0};
    static cilist io___35 = {0, 6, 0, 0, 0};
    static cilist io___36 = {0, 6, 0, 0, 0};
    static cilist io___37 = {0, 6, 0, 0, 0};
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
    static cilist io___58 = {0, 6, 0, 0, 0};
    static cilist io___59 = {0, 6, 0, 0, 0};
    static cilist io___60 = {0, 6, 0, 0, 0};

    /*     Example program to illustrate the idea of reverse communication */
    /*     for a standard nonsymmetric eigenvalue problem. */

    /*     We implement example one of ex-nonsym.doc in DOCUMENTS directory */

    /* \Example-1 */
    /*     ... Suppose we want to solve A*x = lambda*x in regular mode, */
    /*         where A is obtained from the standard central difference */
    /*         discretization of the convection-diffusion operator */
    /*                 (Laplacian u) + rho*(du / dx) */
    /*         on the unit square [0,1]x[0,1] with zero Dirichlet boundary */
    /*         condition. */

    /*     ... OP = A  and  B = I. */

    /*     ... Assume "call av (nx,x,y)" computes y = A*x.c */

    /*     ... Use mode 1 of SNAUPD. */

    /* \BeginLib */

    /* \Routines called: */
    /*     snaupd  ARPACK reverse communication interface routine. */
    /*     sneupd  ARPACK routine that returns Ritz values and (optionally) */
    /*             Ritz vectors. */
    /*     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully. */
    /*     saxpy   Level 1 BLAS that computes y <- alpha*x+y. */
    /*     snrm2   Level 1 BLAS that computes the norm of a vector. */
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
    /* FILE: ndrv1.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2 */

    /* \Remarks */
    /*     1. None */

    /* \EndLib */
    /* --------------------------------------------------------------------------- */

    /*     %-----------------------------% */
    /*     | Define maximum dimensions   | */
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
    /*     | Executable Statements | */
    /*     %-----------------------% */

    /*     %--------------------------------------------------% */
    /*     | The number NX is the number of interior points   | */
    /*     | in the discretization of the 2-dimensional       | */
    /*     | convection-diffusion operator on the unit        | */
    /*     | square with zero Dirichlet boundary condition.   | */
    /*     | The number N(=NX*NX) is the dimension of the     | */
    /*     | matrix.  A standard eigenvalue problem is        | */
    /*     | solved (BMAT = 'I').  NEV is the number of       | */
    /*     | eigenvalues to be approximated.  The user can    | */
    /*     | modify NX, NEV, NCV, WHICH to solve problems of  | */
    /*     | different sizes, and to get different parts of   | */
    /*     | the spectrum.  However, The following            | */
    /*     | conditions must be satisfied:                    | */
    /*     |                   N <= MAXN                      | */
    /*     |                 NEV <= MAXNEV                    | */
    /*     |           NEV + 2 <= NCV <= MAXNCV               | */
    /*     %--------------------------------------------------% */

    nx = 10;
    n = nx * nx;
    nev = 4;
    ncv = 20;
    if (n > 256)
    {
        s_wsle(&io___5);
        do_lio(&c__9, &c__1, " ERROR with _NDRV1: N is greater than MAXN ", (ftnlen)43);
        e_wsle();
        goto L9000;
    }
    else if (nev > 12)
    {
        s_wsle(&io___6);
        do_lio(&c__9, &c__1, " ERROR with _NDRV1: NEV is greater than MAXNEV ", (ftnlen)47);
        e_wsle();
        goto L9000;
    }
    else if (ncv > 30)
    {
        s_wsle(&io___7);
        do_lio(&c__9, &c__1, " ERROR with _NDRV1: NCV is greater than MAXNCV ", (ftnlen)47);
        e_wsle();
        goto L9000;
    }
    *(unsigned char *)bmat = 'I';
    s_copy(which, "SM", (ftnlen)2, (ftnlen)2);

    /*     %-----------------------------------------------------% */
    /*     | The work array WORKL is used in SNAUPD as           | */
    /*     | workspace.  Its dimension LWORKL is set as          | */
    /*     | illustrated below.  The parameter TOL determines    | */
    /*     | the stopping criterion. If TOL<=0, machine          | */
    /*     | precision is used.  The variable IDO is used for    | */
    /*     | reverse communication, and is initially set to 0.   | */
    /*     | Setting INFO=0 indicates that a random vector is    | */
    /*     | generated in SNAUPD to start the Arnoldi iteration. | */
    /*     %-----------------------------------------------------% */

    /* Computing 2nd power */
    i__1 = ncv;
    lworkl = i__1 * i__1 * 3 + ncv * 6;
    tol = 0.f;
    ido = 0;
    info = 0;

    /*     %---------------------------------------------------% */
    /*     | This program uses exact shifts with respect to    | */
    /*     | the current Hessenberg matrix (IPARAM(1) = 1).    | */
    /*     | IPARAM(3) specifies the maximum number of Arnoldi | */
    /*     | iterations allowed.  Mode 1 of SNAUPD is used     | */
    /*     | (IPARAM(7) = 1). All these options can be changed | */
    /*     | by the user. For details see the documentation in | */
    /*     | SNAUPD.                                           | */
    /*     %---------------------------------------------------% */

    ishfts = 1;
    maxitr = 300;
    mode = 1;

    iparam[0] = ishfts;
    iparam[2] = maxitr;
    iparam[6] = mode;

    /*     %-------------------------------------------% */
    /*     | M A I N   L O O P (Reverse communication) | */
    /*     %-------------------------------------------% */

L10:

    /*        %---------------------------------------------% */
    /*        | Repeatedly call the routine SNAUPD and take | */
    /*        | actions indicated by parameter IDO until    | */
    /*        | either convergence is indicated or maxitr   | */
    /*        | has been exceeded.                          | */
    /*        %---------------------------------------------% */

    snaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, &info);

    if (ido == -1 || ido == 1)
    {

        /*           %-------------------------------------------% */
        /*           | Perform matrix vector multiplication      | */
        /*           |                y <--- OP*x                | */
        /*           | The user should supply his/her own        | */
        /*           | matrix vector multiplication routine here | */
        /*           | that takes workd(ipntr(1)) as the input   | */
        /*           | vector, and return the matrix vector      | */
        /*           | product to workd(ipntr(2)).               | */
        /*           %-------------------------------------------% */

        av_(&nx, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

        /*           %-----------------------------------------% */
        /*           | L O O P   B A C K to call SNAUPD again. | */
        /*           %-----------------------------------------% */

        goto L10;
    }

    /*     %----------------------------------------% */
    /*     | Either we have convergence or there is | */
    /*     | an error.                              | */
    /*     %----------------------------------------% */

    if (info < 0)
    {

        /*        %--------------------------% */
        /*        | Error message, check the | */
        /*        | documentation in SNAUPD. | */
        /*        %--------------------------% */

        s_wsle(&io___23);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
        s_wsle(&io___24);
        do_lio(&c__9, &c__1, " Error with _naupd, info = ", (ftnlen)27);
        do_lio(&c__3, &c__1, (char *)&info, (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___25);
        do_lio(&c__9, &c__1, " Check the documentation of _naupd", (ftnlen)34);
        e_wsle();
        s_wsle(&io___26);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
    }
    else
    {

        /*        %-------------------------------------------% */
        /*        | No fatal errors occurred.                 | */
        /*        | Post-Process using SNEUPD.                | */
        /*        |                                           | */
        /*        | Computed eigenvalues may be extracted.    | */
        /*        |                                           | */
        /*        | Eigenvectors may also be computed now if  | */
        /*        | desired.  (indicated by rvec = .true.)    | */
        /*        %-------------------------------------------% */

        rvec = TRUE_;

        sneupd_(&rvec, "A", select, d, &d[30], v, &c__256, &sigmar, &sigmai, workev, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, &ierr);

        /*        %-----------------------------------------------% */
        /*        | The real part of the eigenvalue is returned   | */
        /*        | in the first column of the two dimensional    | */
        /*        | array D, and the imaginary part is returned   | */
        /*        | in the second column of D.  The corresponding | */
        /*        | eigenvectors are returned in the first NEV    | */
        /*        | columns of the two dimensional array V if     | */
        /*        | requested.  Otherwise, an orthogonal basis    | */
        /*        | for the invariant subspace corresponding to   | */
        /*        | the eigenvalues in D is returned in V.        | */
        /*        %-----------------------------------------------% */

        if (ierr != 0)
        {

            /*           %------------------------------------% */
            /*           | Error condition:                   | */
            /*           | Check the documentation of SNEUPD. | */
            /*           %------------------------------------% */

            s_wsle(&io___34);
            do_lio(&c__9, &c__1, " ", (ftnlen)1);
            e_wsle();
            s_wsle(&io___35);
            do_lio(&c__9, &c__1, " Error with _neupd, info = ", (ftnlen)27);
            do_lio(&c__3, &c__1, (char *)&ierr, (ftnlen)sizeof(a_int));
            e_wsle();
            s_wsle(&io___36);
            do_lio(&c__9, &c__1, " Check the documentation of _neupd. ", (ftnlen)36);
            e_wsle();
            s_wsle(&io___37);
            do_lio(&c__9, &c__1, " ", (ftnlen)1);
            e_wsle();
        }
        else
        {

            first = TRUE_;
            nconv = iparam[4];
            i__1 = nconv;
            for (j = 1; j <= i__1; ++j)
            {

                /*               %---------------------------% */
                /*               | Compute the residual norm | */
                /*               |                           | */
                /*               |   ||  A*x - lambda*x ||   | */
                /*               |                           | */
                /*               | for the NCONV accurately  | */
                /*               | computed eigenvalues and  | */
                /*               | eigenvectors.  (iparam(5) | */
                /*               | indicates how many are    | */
                /*               | accurate to the requested | */
                /*               | tolerance)                | */
                /*               %---------------------------% */

                if (d[j + 29] == 0.f)
                {

                    /*                  %--------------------% */
                    /*                  | Ritz value is real | */
                    /*                  %--------------------% */

                    av_(&nx, &v[(j << 8) - 256], ax);
                    r__1 = -d[j - 1];
                    saxpy_(&n, &r__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
                    d[j + 59] = snrm2_(&n, ax, &c__1);
                    d[j + 59] /= (r__1 = d[j - 1], dabs(r__1));
                }
                else if (first)
                {

                    /*                  %------------------------% */
                    /*                  | Ritz value is complex. | */
                    /*                  | Residual of one Ritz   | */
                    /*                  | value of the conjugate | */
                    /*                  | pair is computed.      | */
                    /*                  %------------------------% */

                    av_(&nx, &v[(j << 8) - 256], ax);
                    r__1 = -d[j - 1];
                    saxpy_(&n, &r__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
                    saxpy_(&n, &d[j + 29], &v[(j + 1 << 8) - 256], &c__1, ax, &c__1);
                    d[j + 59] = snrm2_(&n, ax, &c__1);
                    av_(&nx, &v[(j + 1 << 8) - 256], ax);
                    r__1 = -d[j + 29];
                    saxpy_(&n, &r__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
                    r__1 = -d[j - 1];
                    saxpy_(&n, &r__1, &v[(j + 1 << 8) - 256], &c__1, ax, &c__1);
                    r__1 = snrm2_(&n, ax, &c__1);
                    d[j + 59] = slapy2_(&d[j + 59], &r__1);
                    d[j + 59] /= slapy2_(&d[j - 1], &d[j + 29]);
                    d[j + 60] = d[j + 59];
                    first = FALSE_;
                }
                else
                {
                    first = TRUE_;
                }

                /* L20: */
            }

            /*            %-----------------------------% */
            /*            | Display computed residuals. | */
            /*            %-----------------------------% */

            smout_(&c__6, &nconv, &c__3, d, &c__30, &c_n6,"Ritz values (Real,Imag) and relative residuals");
        }

        /*        %-------------------------------------------% */
        /*        | Print additional convergence information. | */
        /*        %-------------------------------------------% */

        if (info == 1)
        {
            s_wsle(&io___42);
            do_lio(&c__9, &c__1, " ", (ftnlen)1);
            e_wsle();
            s_wsle(&io___43);
            do_lio(&c__9, &c__1, " Maximum number of iterations reached.", (ftnlen)38);
            e_wsle();
            s_wsle(&io___44);
            do_lio(&c__9, &c__1, " ", (ftnlen)1);
            e_wsle();
        }
        else if (info == 3)
        {
            s_wsle(&io___45);
            do_lio(&c__9, &c__1, " ", (ftnlen)1);
            e_wsle();
            s_wsle(&io___46);
            do_lio(&c__9, &c__1, " No shifts could be applied during implicit", (ftnlen)43);
            do_lio(&c__9, &c__1, " Arnoldi update, try increasing NCV.", (ftnlen)36);
            e_wsle();
            s_wsle(&io___47);
            do_lio(&c__9, &c__1, " ", (ftnlen)1);
            e_wsle();
        }

        s_wsle(&io___48);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
        s_wsle(&io___49);
        do_lio(&c__9, &c__1, " _NDRV1 ", (ftnlen)8);
        e_wsle();
        s_wsle(&io___50);
        do_lio(&c__9, &c__1, " ====== ", (ftnlen)8);
        e_wsle();
        s_wsle(&io___51);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
        s_wsle(&io___52);
        do_lio(&c__9, &c__1, " Size of the matrix is ", (ftnlen)23);
        do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___53);
        do_lio(&c__9, &c__1, " The number of Ritz values requested is ", (ftnlen)40);
        do_lio(&c__3, &c__1, (char *)&nev, (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___54);
        do_lio(&c__9, &c__1, " The number of Arnoldi vectors generated", (ftnlen)40);
        do_lio(&c__9, &c__1, " (NCV) is ", (ftnlen)10);
        do_lio(&c__3, &c__1, (char *)&ncv, (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___55);
        do_lio(&c__9, &c__1, " What portion of the spectrum: ", (ftnlen)31);
        do_lio(&c__9, &c__1, which, (ftnlen)2);
        e_wsle();
        s_wsle(&io___56);
        do_lio(&c__9, &c__1, " The number of converged Ritz values is ", (ftnlen)40);
        do_lio(&c__3, &c__1, (char *)&nconv, (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___57);
        do_lio(&c__9, &c__1, " The number of Implicit Arnoldi update", (ftnlen)38);
        do_lio(&c__9, &c__1, " iterations taken is ", (ftnlen)21);
        do_lio(&c__3, &c__1, (char *)&iparam[2], (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___58);
        do_lio(&c__9, &c__1, " The number of OP*x is ", (ftnlen)23);
        do_lio(&c__3, &c__1, (char *)&iparam[8], (ftnlen)sizeof(a_int));
        e_wsle();
        s_wsle(&io___59);
        do_lio(&c__9, &c__1, " The convergence criterion is ", (ftnlen)30);
        do_lio(&c__4, &c__1, (char *)&tol, (ftnlen)sizeof(float));
        e_wsle();
        s_wsle(&io___60);
        do_lio(&c__9, &c__1, " ", (ftnlen)1);
        e_wsle();
    }

    /*     %---------------------------% */
    /*     | Done with program sndrv1. | */
    /*     %---------------------------% */

L9000:

    return 0;
} /* MAIN__ */

/* ========================================================================== */

/*     matrix vector subroutine */

/*     The matrix used is the 2 dimensional convection-diffusion */
/*     operator discretized using central difference. */

int av_(a_int *nx, float *v, float *w)
{
    /* System generated locals */
    a_int i__1;
    float r__1;

    /* Local variables */
    a_int j;
    float h2;
    a_int lo;

    /*     Computes w <--- OP*v, where OP is the nx*nx by nx*nx block */
    /*     tridiagonal matrix */

    /*                  | T -I          | */
    /*                  |-I  T -I       | */
    /*             OP = |   -I  T       | */
    /*                  |        ...  -I| */
    /*                  |           -I T| */

    /*     derived from the standard central difference discretization */
    /*     of the 2 dimensional convection-diffusion operator */
    /*     (Laplacian u) + rho*(du/dx) on a unit square with zero boundary */
    /*     condition. */

    /*     When rho*h/2 <= 1, the discrete convection-diffusion operator */
    /*     has real eigenvalues.  When rho*h/2 > 1, it has COMPLEX */
    /*     eigenvalues. */

    /*     The subroutine TV is called to compute y<---T*x. */

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    h2 = 1.f / (float)((*nx + 1) * (*nx + 1));

    tv_(nx, &v[1], &w[1]);
    r__1 = -1.f / h2;
    saxpy_(nx, &r__1, &v[*nx + 1], &c__1, &w[1], &c__1);

    i__1 = *nx - 1;
    for (j = 2; j <= i__1; ++j)
    {
        lo = (j - 1) * *nx;
        tv_(nx, &v[lo + 1], &w[lo + 1]);
        r__1 = -1.f / h2;
        saxpy_(nx, &r__1, &v[lo - *nx + 1], &c__1, &w[lo + 1], &c__1);
        r__1 = -1.f / h2;
        saxpy_(nx, &r__1, &v[lo + *nx + 1], &c__1, &w[lo + 1], &c__1);
        /* L10: */
    }

    lo = (*nx - 1) * *nx;
    tv_(nx, &v[lo + 1], &w[lo + 1]);
    r__1 = -1.f / h2;
    saxpy_(nx, &r__1, &v[lo - *nx + 1], &c__1, &w[lo + 1], &c__1);

    return 0;
} /* av_ */

/* ========================================================================= */
int tv_(a_int *nx, float *x, float *y)
{
    /* System generated locals */
    a_int i__1;

    /* Local variables */
    float h;
    a_int j;
    float h2, dd, dl, du;

    /*     Compute the matrix vector multiplication y<---T*x */
    /*     where T is a nx by nx tridiagonal matrix with DD on the */
    /*     diagonal, DL on the subdiagonal, and DU on the superdiagonal. */

    /*     When rho*h/2 <= 1, the discrete convection-diffusion operator */
    /*     has real eigenvalues.  When rho*h/2 > 1, it has COMPLEX */
    /*     eigenvalues. */

    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    h = 1.f / (float)(*nx + 1);
    h2 = h * h;
    dd = 4.f / h2;
    dl = -1.f / h2 - 0.f / h;
    du = -1.f / h2 + 0.f / h;

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

/* Main program alias */ int sndrv1_()
{
    MAIN__();
    return 0;
}
