/* D:\projects\Fortran\arpack-ng-3.9.1-patched\EXAMPLES\NONSYM\dndrv1.f -- translated by f2c (version 20230428).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__256 = 256;
static integer c__3 = 3;
static integer c__6 = 6;
static integer c__30 = 30;
static integer c_n6 = -6;
static integer c__5 = 5;

/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    doublereal d__[90]	/* was [30][3] */;
    integer j, n;
    doublereal v[7680]	/* was [256][30] */;
    extern /* Subroutine */ int av_(integer *, doublereal *, doublereal *);
    doublereal ax[256];
    integer nx, ido, ncv, nev;
    doublereal tol;
    char bmat[1];
    integer mode, info;
    logical rvec;
    integer ierr;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    char which[2];
    doublereal resid[256];
    integer nconv;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    doublereal workd[768];
    logical first;
    extern /* Subroutine */ int dmout_(integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, char *, ftnlen);
    integer ipntr[14];
    doublereal workl[2880];
    extern doublereal dlapy2_(doublereal *, doublereal *);
    integer iparam[11];
    doublereal sigmai;
    extern /* Subroutine */ int dnaupd_(integer *, char *, integer *, char *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen);
    logical select[30];
    extern /* Subroutine */ int dneupd_(logical *, char *, logical *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, char *, integer *, char *, integer *,
	     doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    doublereal sigmar;
    integer ishfts, maxitr, lworkl;
    doublereal workev[90];

    /* Fortran I/O blocks */
    static cilist io___5 = { 0, 6, 0, 0, 0 };
    static cilist io___6 = { 0, 6, 0, 0, 0 };
    static cilist io___7 = { 0, 6, 0, 0, 0 };
    static cilist io___23 = { 0, 6, 0, 0, 0 };
    static cilist io___24 = { 0, 6, 0, 0, 0 };
    static cilist io___25 = { 0, 6, 0, 0, 0 };
    static cilist io___26 = { 0, 6, 0, 0, 0 };
    static cilist io___34 = { 0, 6, 0, 0, 0 };
    static cilist io___35 = { 0, 6, 0, 0, 0 };
    static cilist io___36 = { 0, 6, 0, 0, 0 };
    static cilist io___37 = { 0, 6, 0, 0, 0 };
    static cilist io___42 = { 0, 6, 0, 0, 0 };
    static cilist io___43 = { 0, 6, 0, 0, 0 };
    static cilist io___44 = { 0, 6, 0, 0, 0 };
    static cilist io___45 = { 0, 6, 0, 0, 0 };
    static cilist io___46 = { 0, 6, 0, 0, 0 };
    static cilist io___47 = { 0, 6, 0, 0, 0 };
    static cilist io___48 = { 0, 6, 0, 0, 0 };
    static cilist io___49 = { 0, 6, 0, 0, 0 };
    static cilist io___50 = { 0, 6, 0, 0, 0 };
    static cilist io___51 = { 0, 6, 0, 0, 0 };
    static cilist io___52 = { 0, 6, 0, 0, 0 };
    static cilist io___53 = { 0, 6, 0, 0, 0 };
    static cilist io___54 = { 0, 6, 0, 0, 0 };
    static cilist io___55 = { 0, 6, 0, 0, 0 };
    static cilist io___56 = { 0, 6, 0, 0, 0 };
    static cilist io___57 = { 0, 6, 0, 0, 0 };
    static cilist io___58 = { 0, 6, 0, 0, 0 };
    static cilist io___59 = { 0, 6, 0, 0, 0 };
    static cilist io___60 = { 0, 6, 0, 0, 0 };




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

/*     ... Use mode 1 of DNAUPD. */

/* \BeginLib */

/* \Routines called: */
/*     dnaupd  ARPACK reverse communication interface routine. */
/*     dneupd  ARPACK routine that returns Ritz values and (optionally) */
/*             Ritz vectors. */
/*     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully. */
/*     daxpy   Level 1 BLAS that computes y <- alpha*x+y. */
/*     dnrm2   Level 1 BLAS that computes the norm of a vector. */
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
    if (n > 256) {
	s_wsle(&io___5);
	do_lio(&c__9, &c__1, " ERROR with _NDRV1: N is greater than MAXN ", (
		ftnlen)43);
	e_wsle();
	goto L9000;
    } else if (nev > 12) {
	s_wsle(&io___6);
	do_lio(&c__9, &c__1, " ERROR with _NDRV1: NEV is greater than MAXNEV "
		, (ftnlen)47);
	e_wsle();
	goto L9000;
    } else if (ncv > 30) {
	s_wsle(&io___7);
	do_lio(&c__9, &c__1, " ERROR with _NDRV1: NCV is greater than MAXNCV "
		, (ftnlen)47);
	e_wsle();
	goto L9000;
    }
    *(unsigned char *)bmat = 'I';
    s_copy(which, "SM", (ftnlen)2, (ftnlen)2);

/*     %-----------------------------------------------------% */
/*     | The work array WORKL is used in DNAUPD as           | */
/*     | workspace.  Its dimension LWORKL is set as          | */
/*     | illustrated below.  The parameter TOL determines    | */
/*     | the stopping criterion. If TOL<=0, machine          | */
/*     | precision is used.  The variable IDO is used for    | */
/*     | reverse communication, and is initially set to 0.   | */
/*     | Setting INFO=0 indicates that a random vector is    | */
/*     | generated in DNAUPD to start the Arnoldi iteration. | */
/*     %-----------------------------------------------------% */

/* Computing 2nd power */
    i__1 = ncv;
    lworkl = i__1 * i__1 * 3 + ncv * 6;
    tol = 0.;
    ido = 0;
    info = 0;

/*     %---------------------------------------------------% */
/*     | This program uses exact shifts with respect to    | */
/*     | the current Hessenberg matrix (IPARAM(1) = 1).    | */
/*     | IPARAM(3) specifies the maximum number of Arnoldi | */
/*     | iterations allowed.  Mode 1 of DNAUPD is used     | */
/*     | (IPARAM(7) = 1). All these options can be changed | */
/*     | by the user. For details see the documentation in | */
/*     | DNAUPD.                                           | */
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
/*        | Repeatedly call the routine DNAUPD and take | */
/*        | actions indicated by parameter IDO until    | */
/*        | either convergence is indicated or maxitr   | */
/*        | has been exceeded.                          | */
/*        %---------------------------------------------% */

    dnaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, 
	    iparam, ipntr, workd, workl, &lworkl, &info, (ftnlen)1, (ftnlen)2)
	    ;

    if (ido == -1 || ido == 1) {

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
/*           | L O O P   B A C K to call DNAUPD again. | */
/*           %-----------------------------------------% */

	goto L10;

    }

/*     %----------------------------------------% */
/*     | Either we have convergence or there is | */
/*     | an error.                              | */
/*     %----------------------------------------% */

    if (info < 0) {

/*        %--------------------------% */
/*        | Error message, check the | */
/*        | documentation in DNAUPD. | */
/*        %--------------------------% */

	s_wsle(&io___23);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___24);
	do_lio(&c__9, &c__1, " Error with _naupd, info = ", (ftnlen)27);
	do_lio(&c__3, &c__1, (char *)&info, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___25);
	do_lio(&c__9, &c__1, " Check the documentation of _naupd", (ftnlen)34)
		;
	e_wsle();
	s_wsle(&io___26);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();

    } else {

/*        %-------------------------------------------% */
/*        | No fatal errors occurred.                 | */
/*        | Post-Process using DNEUPD.                | */
/*        |                                           | */
/*        | Computed eigenvalues may be extracted.    | */
/*        |                                           | */
/*        | Eigenvectors may also be computed now if  | */
/*        | desired.  (indicated by rvec = .true.)    | */
/*        %-------------------------------------------% */

	rvec = TRUE_;

	dneupd_(&rvec, "A", select, d__, &d__[30], v, &c__256, &sigmar, &
		sigmai, workev, bmat, &n, which, &nev, &tol, resid, &ncv, v, &
		c__256, iparam, ipntr, workd, workl, &lworkl, &ierr, (ftnlen)
		1, (ftnlen)1, (ftnlen)2);

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

	if (ierr != 0) {

/*           %------------------------------------% */
/*           | Error condition:                   | */
/*           | Check the documentation of DNEUPD. | */
/*           %------------------------------------% */

	    s_wsle(&io___34);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___35);
	    do_lio(&c__9, &c__1, " Error with _neupd, info = ", (ftnlen)27);
	    do_lio(&c__3, &c__1, (char *)&ierr, (ftnlen)sizeof(integer));
	    e_wsle();
	    s_wsle(&io___36);
	    do_lio(&c__9, &c__1, " Check the documentation of _neupd. ", (
		    ftnlen)36);
	    e_wsle();
	    s_wsle(&io___37);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();

	} else {

	    first = TRUE_;
	    nconv = iparam[4];
	    i__1 = nconv;
	    for (j = 1; j <= i__1; ++j) {

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

		if (d__[j + 29] == 0.) {

/*                  %--------------------% */
/*                  | Ritz value is real | */
/*                  %--------------------% */

		    av_(&nx, &v[(j << 8) - 256], ax);
		    d__1 = -d__[j - 1];
		    daxpy_(&n, &d__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
		    d__[j + 59] = dnrm2_(&n, ax, &c__1);
		    d__[j + 59] /= (d__1 = d__[j - 1], abs(d__1));

		} else if (first) {

/*                  %------------------------% */
/*                  | Ritz value is complex. | */
/*                  | Residual of one Ritz   | */
/*                  | value of the conjugate | */
/*                  | pair is computed.      | */
/*                  %------------------------% */

		    av_(&nx, &v[(j << 8) - 256], ax);
		    d__1 = -d__[j - 1];
		    daxpy_(&n, &d__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
		    daxpy_(&n, &d__[j + 29], &v[(j + 1 << 8) - 256], &c__1, 
			    ax, &c__1);
		    d__[j + 59] = dnrm2_(&n, ax, &c__1);
		    av_(&nx, &v[(j + 1 << 8) - 256], ax);
		    d__1 = -d__[j + 29];
		    daxpy_(&n, &d__1, &v[(j << 8) - 256], &c__1, ax, &c__1);
		    d__1 = -d__[j - 1];
		    daxpy_(&n, &d__1, &v[(j + 1 << 8) - 256], &c__1, ax, &
			    c__1);
		    d__1 = dnrm2_(&n, ax, &c__1);
		    d__[j + 59] = dlapy2_(&d__[j + 59], &d__1);
		    d__[j + 59] /= dlapy2_(&d__[j - 1], &d__[j + 29]);
		    d__[j + 60] = d__[j + 59];
		    first = FALSE_;
		} else {
		    first = TRUE_;
		}

/* L20: */
	    }

/*            %-----------------------------% */
/*            | Display computed residuals. | */
/*            %-----------------------------% */

	    dmout_(&c__6, &nconv, &c__3, d__, &c__30, &c_n6, "Ritz values (R"
		    "eal,Imag) and relative residuals", (ftnlen)46);
	}

/*        %-------------------------------------------% */
/*        | Print additional convergence information. | */
/*        %-------------------------------------------% */

	if (info == 1) {
	    s_wsle(&io___42);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___43);
	    do_lio(&c__9, &c__1, " Maximum number of iterations reached.", (
		    ftnlen)38);
	    e_wsle();
	    s_wsle(&io___44);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	} else if (info == 3) {
	    s_wsle(&io___45);
	    do_lio(&c__9, &c__1, " ", (ftnlen)1);
	    e_wsle();
	    s_wsle(&io___46);
	    do_lio(&c__9, &c__1, " No shifts could be applied during implicit"
		    , (ftnlen)43);
	    do_lio(&c__9, &c__1, " Arnoldi update, try increasing NCV.", (
		    ftnlen)36);
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
	do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___53);
	do_lio(&c__9, &c__1, " The number of Ritz values requested is ", (
		ftnlen)40);
	do_lio(&c__3, &c__1, (char *)&nev, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___54);
	do_lio(&c__9, &c__1, " The number of Arnoldi vectors generated", (
		ftnlen)40);
	do_lio(&c__9, &c__1, " (NCV) is ", (ftnlen)10);
	do_lio(&c__3, &c__1, (char *)&ncv, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___55);
	do_lio(&c__9, &c__1, " What portion of the spectrum: ", (ftnlen)31);
	do_lio(&c__9, &c__1, which, (ftnlen)2);
	e_wsle();
	s_wsle(&io___56);
	do_lio(&c__9, &c__1, " The number of converged Ritz values is ", (
		ftnlen)40);
	do_lio(&c__3, &c__1, (char *)&nconv, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___57);
	do_lio(&c__9, &c__1, " The number of Implicit Arnoldi update", (
		ftnlen)38);
	do_lio(&c__9, &c__1, " iterations taken is ", (ftnlen)21);
	do_lio(&c__3, &c__1, (char *)&iparam[2], (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___58);
	do_lio(&c__9, &c__1, " The number of OP*x is ", (ftnlen)23);
	do_lio(&c__3, &c__1, (char *)&iparam[8], (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___59);
	do_lio(&c__9, &c__1, " The convergence criterion is ", (ftnlen)30);
	do_lio(&c__5, &c__1, (char *)&tol, (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___60);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();

    }

/*     %---------------------------% */
/*     | Done with program dndrv1. | */
/*     %---------------------------% */

L9000:

    return 0;
} /* MAIN__ */


/* ========================================================================== */

/*     matrix vector subroutine */

/*     The matrix used is the 2 dimensional convection-diffusion */
/*     operator discretized using central difference. */

/* Subroutine */ int av_(integer *nx, doublereal *v, doublereal *w)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    integer j;
    doublereal h2;
    integer lo;
    extern /* Subroutine */ int tv_(integer *, doublereal *, doublereal *), 
	    daxpy_(integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *);


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
    h2 = 1. / (doublereal) ((*nx + 1) * (*nx + 1));

    tv_(nx, &v[1], &w[1]);
    d__1 = -1. / h2;
    daxpy_(nx, &d__1, &v[*nx + 1], &c__1, &w[1], &c__1);

    i__1 = *nx - 1;
    for (j = 2; j <= i__1; ++j) {
	lo = (j - 1) * *nx;
	tv_(nx, &v[lo + 1], &w[lo + 1]);
	d__1 = -1. / h2;
	daxpy_(nx, &d__1, &v[lo - *nx + 1], &c__1, &w[lo + 1], &c__1);
	d__1 = -1. / h2;
	daxpy_(nx, &d__1, &v[lo + *nx + 1], &c__1, &w[lo + 1], &c__1);
/* L10: */
    }

    lo = (*nx - 1) * *nx;
    tv_(nx, &v[lo + 1], &w[lo + 1]);
    d__1 = -1. / h2;
    daxpy_(nx, &d__1, &v[lo - *nx + 1], &c__1, &w[lo + 1], &c__1);

    return 0;
} /* av_ */

/* ========================================================================= */
/* Subroutine */ int tv_(integer *nx, doublereal *x, doublereal *y)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    doublereal h__;
    integer j;
    doublereal h2, dd, dl, du;




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
    h__ = 1. / (doublereal) (*nx + 1);
    h2 = h__ * h__;
    dd = 4. / h2;
    dl = -1. / h2 - 0. / h__;
    du = -1. / h2 + 0. / h__;

    y[1] = dd * x[1] + du * x[2];
    i__1 = *nx - 1;
    for (j = 2; j <= i__1; ++j) {
	y[j] = dl * x[j - 1] + dd * x[j] + du * x[j + 1];
/* L10: */
    }
    y[*nx] = dl * x[*nx - 1] + dd * x[*nx];
    return 0;
} /* tv_ */

/* Main program alias */ int dndrv1_ () { MAIN__ (); return 0; }
