/* EXAMPLES\NONSYM\dndrv6.f -- translated by f2c (version 20230428). */

#include <stdlib.h>
#include "arpack_internal.h"

static a_int i_one = 1;
static a_int c__256 = 256;

void mv_(const a_int n, double *v, double *w);
void av_(const a_int n, double *v, double *w);

/**
 * \BeginDoc
 *
 *     Simple program to illustrate the idea of reverse communication
 *     in shift-invert mode for a generalized nonsymmetric eigenvalue problem.
 *
 *     We implement example six of ex-nonsym.doc in DOCUMENTS directory
 *
 * \Example-6
 *
 *     ... Suppose we want to solve A*x = lambda*B*x in shift-invert mode
 *         The matrix A is the tridiagonal matrix with 2 on the diagonal,
 *         -2 on the subdiagonal and 3 on the superdiagonal.  The matrix M
 *         is the tridiagonal matrix with 4 on the diagonal and 1 on the
 *         off-diagonals.
 *     ... The shift sigma is a complex number (sigmar, sigmai).
 *     ... OP = Imaginary_Part{inv[A-(SIGMAR,SIGMAI)*M]*M  and  B = M.
 *     ... Use mode 4 of DNAUPD.
 *
 * \EndDoc
 *
 * \BeginLib
 *
 * Routines called:
 *     dnaupd  ARPACK reverse communication interface routine.
 *     dneupd  ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     zgttrf  LAPACK complex matrix factorization routine.
 *     zgttrs  LAPACK complex linear system solve routine.
 *     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
 *     ddot    Level 1 BLAS that computes the dot product of two vectors.
 *     dnrm2   Level 1 BLAS that computes the norm of a vector.
 *     av      Matrix vector subroutine that computes A*x.
 *     mv      Matrix vector subroutine that computes M*x.
 *
 * \EndLib
 */
int main()
{
    /* System generated locals */
    a_int i__1, i__2, i__3;
    double d__1, d__2;
    a_dcomplex z__1;

    /* Local variables */
    a_bool select[25];
    a_int iparam[11];
    a_int ipntr[14];
    a_bool rvec, first;
    a_dcomplex c1, c2, c3;
    a_int j, n, ido, ncv, nev, info, ierr = 0;
    a_int mode, nconv, ishfts, lworkl, maxitr;
    char *bmat, *which;
    double tol, deni, denr, numi, numr, sigmai, sigmar;

    /* Define maximum dimensions for all arrays. */

    const int MAXN   = 256; /* Maximum dimension of the A allowed. */
    const int MAXNEV =  10; /* Maximum NEV allowed */
    const int MAXNCV =  25; /* Maximum NCV allowed */

    /* -------------------------------------------------- */
    /* The number N is the dimension of the matrix.  A    */
    /* generalized eigenvalue problem is solved (BMAT =   */
    /* 'G').  NEV is the number of eigenvalues (closest   */
    /* to the shift (SIGMAR,SIGMAI)) to be approximated.  */
    /* Since the shift-invert mode is used, WHICH is set  */
    /* to 'LM'.  The user can modify NEV, NCV, SIGMA to   */
    /* solve problems of different sizes, and to get      */
    /* different parts of the spectrum. However, The      */
    /* following conditions must be satisfied:            */
    /*                     N <= MAXN,                     */
    /*                   NEV <= MAXNEV,                   */
    /*               NEV + 2 <= NCV <= MAXNCV             */
    /* -------------------------------------------------- */

    n = 100;
    nev = 4;
    ncv = 20;
    if (n > 256)
    {
        printf(" ERROR with _NDRV6: N is greater than MAXN \n");
        return ierr;
    }
    else if (nev > 10)
    {
        printf(" ERROR with _NDRV6: NEV is greater than MAXNEV \n");
        return ierr;
    }
    else if (ncv > 25)
    {
        printf(" ERROR with _NDRV6: NCV is greater than MAXNCV \n");
        return ierr;
    }
    bmat = "G";
    which = "LM";
    sigmar = .4;
    sigmai = .6;

    /* -------------------------------------------------- */
    /* Construct C = A - (SIGMAR,SIGMAI)*M in complex     */
    /* arithmetic, and factor C in complex arithmetic     */
    /* (using LAPACK subroutine zgttrf). The matrix A is  */
    /* chosen to be the tridiagonal matrix with -2 on the */
    /* subdiagonal, 2 on the diagonal and 3 on the        */
    /* superdiagonal. The matrix M is chosen to be the    */
    /* symmetric tridiagonal matrix with 4 on the         */
    /* diagonal and 1 on the off-diagonals.               */
    /* -------------------------------------------------- */

    a_int* ipiv = (a_int*)malloc(sizeof(a_int) * 256);

    a_dcomplex* cdd = (a_dcomplex*)malloc(sizeof(a_dcomplex) * 256);
    a_dcomplex* cdl = (a_dcomplex*)malloc(sizeof(a_dcomplex) * 256);
    a_dcomplex* cdu = (a_dcomplex*)malloc(sizeof(a_dcomplex) * 256);
    a_dcomplex* cdu2 = (a_dcomplex*)malloc(sizeof(a_dcomplex) * 256);

    d__1 = -2. - sigmar;
    d__2 = -sigmai;
    z__1.r = d__1, z__1.i = d__2;
    c1.r = z__1.r, c1.i = z__1.i;
    d__1 = 2. - sigmar * 4.;
    d__2 = sigmai * -4.;
    z__1.r = d__1, z__1.i = d__2;
    c2.r = z__1.r, c2.i = z__1.i;
    d__1 = 3. - sigmar;
    d__2 = -sigmai;
    z__1.r = d__1, z__1.i = d__2;
    c3.r = z__1.r, c3.i = z__1.i;

    i__1 = n - 1;
    for (j = 1; j <= i__1; ++j)
    {
        i__2 = j - 1;
        cdl[i__2].r = c1.r, cdl[i__2].i = c1.i;
        i__2 = j - 1;
        cdd[i__2].r = c2.r, cdd[i__2].i = c2.i;
        i__2 = j - 1;
        cdu[i__2].r = c3.r, cdu[i__2].i = c3.i;
    }
    i__1 = n - 1;
    cdd[i__1].r = c2.r, cdd[i__1].i = c2.i;

    zgttrf_(&n, cdl, cdd, cdu, cdu2, ipiv, &ierr);
    if (ierr != 0)
    {
        printf(" \n");
        printf(" ERROR with _gttrf in _NDRV6.\n");
        printf(" \n");
        return ierr;
    }

    /* --------------------------------------------------- */
    /* The work array WORKL is used in DNAUPD as           */
    /* workspace.  Its dimension LWORKL is set as          */
    /* illustrated below.  The parameter TOL determines    */
    /* the stopping criterion. If TOL<=0, machine          */
    /* precision is used.  The variable IDO is used for    */
    /* reverse communication, and is initially set to 0.   */
    /* Setting INFO=0 indicates that a random vector is    */
    /* generated in DNAUPD to start the Arnoldi iteration. */
    /* --------------------------------------------------- */

    lworkl = ncv * ncv * 3 + ncv * 6;
    tol = 0.;
    ido = 0;
    info = 0;
    nconv = 0;

    a_dcomplex* ctemp = (a_dcomplex*)malloc(sizeof(a_dcomplex) * 256);
    double* d = (double*)malloc(sizeof(double) * 25 * 3);
    double* v = (double*)malloc(sizeof(double) * 256 * 25);
    double* ax = (double*)malloc(sizeof(double) * 256);
    double* mx = (double*)malloc(sizeof(double) * 256);
    double* resid = (double*)malloc(sizeof(double) * 256);
    double* workd = (double*)malloc(sizeof(double) * 768);
    double* workl = (double*)malloc(sizeof(double) * 2025);
    double* workev = (double*)malloc(sizeof(double) * 75);

    /* ------------------------------------------------- */
    /* This program uses exact shift with respect to     */
    /* the current Hessenberg matrix (IPARAM(1) = 1).    */
    /* IPARAM(3) specifies the maximum number of Arnoldi */
    /* iterations allowed.  Mode 3 of DNAUPD is used     */
    /* (IPARAM(7) = 3).  All these options can be        */
    /* changed by the user. For details, see the         */
    /* documentation in DNAUPD.                          */
    /* ------------------------------------------------- */

    ishfts = 1;
    maxitr = 300;
    mode = 3;

    iparam[0] = ishfts;
    iparam[2] = maxitr;
    iparam[6] = mode;

    /* ---------------------------------------- */
    /* M A I N   L O O P(Reverse communication) */
    /* ---------------------------------------- */

L20:

    /* ------------------------------------------- */
    /* Repeatedly call the routine DNAUPD and take */
    /* actions indicated by parameter IDO until    */
    /* either convergence is indicated or maxitr   */
    /* has been exceeded.                          */
    /* ------------------------------------------- */

    dnaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, &info);

    if (ido == -1)
    {

        /* ---------------------------------------------------------- */
        /*                           Perform                          */
        /* y <--- OP*x = Imaginary_Part{inv[A-(SIGMAR,SIGMAI)*M]*M*x} */
        /* to force starting vector into the range of OP. The user    */
        /* should supply his/her own matrix vector multiplication     */
        /* routine and a complex linear system solver.  The matrix    */
        /* vector multiplication routine should take workd(ipntr(1))  */
        /* as the input. The final result (a real vector) should be   */
        /* returned to workd(ipntr(2)).                               */
        /* ---------------------------------------------------------- */

        mv_(n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
        i__1 = n;
        for (j = 1; j <= i__1; ++j)
        {
            i__2 = j - 1;
            i__3 = ipntr[1] + j - 2;
            z__1.r = workd[i__3], z__1.i = 0.;
            ctemp[i__2].r = z__1.r, ctemp[i__2].i = z__1.i;
        }

        zgttrs_("N", &n, &i_one, cdl, cdd, cdu, cdu2, ipiv, ctemp, &c__256, &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" ERROR with _gttrs in _NDRV6.\n");
            printf(" \n");
            return ierr;
        }
        i__1 = n;
        for (j = 1; j <= i__1; ++j)
        {
            workd[ipntr[1] + j - 2] = ctemp[j - 1].i;
        }

        /* --------------------------------------- */
        /* L O O P   B A C K to call DNAUPD again. */
        /* --------------------------------------- */

        goto L20;
    }
    else if (ido == 1)
    {

        /* ---------------------------------------------------------- */
        /*                          Perform                           */
        /* y <--- OP*x = Imaginary_Part{inv[A-(SIGMAR,SIGMAI)*M]*M*x} */
        /* M*x has been saved in workd(ipntr(3)). The user only need  */
        /* the complex linear system solver here that takes           */
        /* complex[workd(ipntr(3))] as input, and returns the result  */
        /* to workd(ipntr(2)).                                        */
        /* ---------------------------------------------------------- */

        i__1 = n;
        for (j = 1; j <= i__1; ++j)
        {
            i__2 = j - 1;
            i__3 = ipntr[2] + j - 2;
            z__1.r = workd[i__3], z__1.i = 0.;
            ctemp[i__2].r = z__1.r, ctemp[i__2].i = z__1.i;
        }
        zgttrs_("N", &n, &i_one, cdl, cdd, cdu, cdu2, ipiv, ctemp, &c__256, &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" ERROR with _gttrs in _NDRV6.\n");
            printf(" \n");
            return ierr;
        }
        i__1 = n;
        for (j = 1; j <= i__1; ++j)
        {
            workd[ipntr[1] + j - 2] = ctemp[j - 1].i;
        }

        /* --------------------------------------- */
        /* L O O P   B A C K to call DNAUPD again. */
        /* --------------------------------------- */

        goto L20;
    }
    else if (ido == 2)
    {

        /* ------------------------------------------- */
        /*          Perform  y <--- M*x                */
        /* Need matrix vector multiplication routine   */
        /* here that takes workd(ipntr(1)) as input    */
        /* and returns the result to workd(ipntr(2)).  */
        /* ------------------------------------------- */

        mv_(n, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);

        /* --------------------------------------- */
        /* L O O P   B A C K to call DNAUPD again. */
        /* --------------------------------------- */

        goto L20;
    }

    /* --------------------------------------- */
    /* Either we have convergence, or there is */
    /* an error.                               */
    /* --------------------------------------- */

    if (info < 0)
    {

        /* ------------------------ */
        /* Error message, check the */
        /* documentation in DNAUPD. */
        /* ------------------------ */

        printf(" \n");
        printf(" Error with _naupd info = %d\n", info);
        printf(" Check the documentation of _naupd.\n");
        printf(" \n");
    }
    else
    {

        /* ----------------------------------------- */
        /* No fatal errors occurred.                 */
        /* Post-Process using DNEUPD.                */
        /*                                           */
        /* Computed eigenvalues may be extracted.    */
        /*                                           */
        /* Eigenvectors may also be computed now if  */
        /* desired.  (indicated by rvec = .true.)    */
        /* ----------------------------------------- */

        rvec = TRUE_;
        dneupd_(&rvec, "A", select, d, &d[25], v, &c__256, &sigmar, &sigmai, workev, bmat, &n, which, &nev, &tol, resid, &ncv, v, &c__256, iparam, ipntr, workd, workl, &lworkl, &ierr);

        /* --------------------------------------------- */
        /* The real part of the eigenvalue is returned   */
        /* in the first column of the two dimensional    */
        /* array D, and the IMAGINARY part is returned   */
        /* in the second column of D.  The corresponding */
        /* eigenvectors are returned in the first NEV    */
        /* columns of the two dimensional array V if     */
        /* requested.  Otherwise, an orthogonal basis    */
        /* for the invariant subspace corresponding to   */
        /* the eigenvalues in D is returned in V.        */
        /* --------------------------------------------- */

        if (ierr != 0)
        {

            /* ---------------------------------- */
            /* Error condition:                   */
            /* Check the documentation of DNEUPD. */
            /* ---------------------------------- */

            printf(" \n");
            printf(" Error with _neupd info = %d\n", ierr);
            printf(" Check the documentation of _neupd. \n");
            printf(" \n");
        }
        else
        {

            first = TRUE_;
            nconv = iparam[4];
            i__1 = nconv;
            for (j = 1; j <= i__1; ++j)
            {

                /* ----------------------------------- */
                /* Use Rayleigh Quotient to recover    */
                /* eigenvalues of the original problem.*/
                /* ----------------------------------- */

                if (d[j + 24] == 0.)
                {

                    /* -------------------------- */
                    /* Eigenvalue is real.        */
                    /* Compute d = x'(Ax)/x'(Mx). */
                    /* -------------------------- */

                    av_(n, &v[(j << 8) - 256], ax);
                    numr = ddot_(&n, &v[(j << 8) - 256], &i_one, ax, &i_one);
                    mv_(n, &v[(j << 8) - 256], ax);
                    denr = ddot_(&n, &v[(j << 8) - 256], &i_one, ax, &i_one);
                    d[j - 1] = numr / denr;
                }
                else if (first)
                {

                    /* ---------------------- */
                    /* Eigenvalue is complex. */
                    /* Compute the first one  */
                    /* of the conjugate pair. */
                    /* ---------------------- */

                    /* -------------- */
                    /* Compute x'(Ax) */
                    /* -------------- */

                    av_(n, &v[(j << 8) - 256], ax);
                    numr = ddot_(&n, &v[(j << 8) - 256], &i_one, ax, &i_one);
                    numi = ddot_(&n, &v[(j + 1 << 8) - 256], &i_one, ax, &i_one);
                    av_(n, &v[(j + 1 << 8) - 256], ax);
                    numr += ddot_(&n, &v[(j + 1 << 8) - 256], &i_one, ax, &i_one);
                    numi = -numi + ddot_(&n, &v[(j << 8) - 256], &i_one, ax, &i_one);

                    /* -------------- */
                    /* Compute x'(Mx) */
                    /* -------------- */

                    mv_(n, &v[(j << 8) - 256], ax);
                    denr = ddot_(&n, &v[(j << 8) - 256], &i_one, ax, &i_one);
                    deni = ddot_(&n, &v[(j + 1 << 8) - 256], &i_one, ax, &i_one);
                    mv_(n, &v[(j + 1 << 8) - 256], ax);
                    denr += ddot_(&n, &v[(j + 1 << 8) - 256], &i_one, ax, &i_one);
                    deni = -deni + ddot_(&n, &v[(j << 8) - 256], &i_one, ax, &i_one);

                    /* -------------- */
                    /* d=x'(Ax)/x'(Mx)*/
                    /* -------------- */

                    d[j - 1] = (numr * denr + numi * deni) / dlapy2_(&denr, &deni);
                    d[j + 24] = (numi * denr - numr * deni) / dlapy2_(&denr, &deni);
                    first = FALSE_;
                }
                else
                {

                    /* ---------------------------- */
                    /* Get the second eigenvalue of */
                    /* the conjugate pair by taking */
                    /* the conjugate of the last    */
                    /* eigenvalue computed.         */
                    /* ---------------------------- */

                    d[j - 1] = d[j - 2];
                    d[j + 24] = -d[j + 23];
                    first = TRUE_;
                }
            }

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

            first = TRUE_;
            nconv = iparam[4];
            i__1 = nconv;
            for (j = 1; j <= i__1; ++j)
            {

                if (d[j + 24] == 0.)
                {

                    /* ------------------ */
                    /* Ritz value is real */
                    /* ------------------ */

                    av_(n, &v[(j << 8) - 256], ax);
                    mv_(n, &v[(j << 8) - 256], mx);
                    d__1 = -d[j - 1];
                    daxpy_(&n, &d__1, mx, &i_one, ax, &i_one);
                    d[j + 49] = dnrm2_(&n, ax, &i_one);
                    d[j + 49] /= (d__1 = d[j - 1], abs(d__1));
                }
                else if (first)
                {

                    /* ---------------------- */
                    /* Ritz value is complex  */
                    /* Residual of one Ritz   */
                    /* value of the conjugate */
                    /* pair is computed.      */
                    /* ---------------------- */

                    av_(n, &v[(j << 8) - 256], ax);
                    mv_(n, &v[(j << 8) - 256], mx);
                    d__1 = -d[j - 1];
                    daxpy_(&n, &d__1, mx, &i_one, ax, &i_one);
                    mv_(n, &v[(j + 1 << 8) - 256], mx);
                    daxpy_(&n, &d[j + 24], mx, &i_one, ax, &i_one);
                    d[j + 49] = dnrm2_(&n, ax, &i_one);
                    av_(n, &v[(j + 1 << 8) - 256], ax);
                    mv_(n, &v[(j + 1 << 8) - 256], mx);
                    d__1 = -d[j - 1];
                    daxpy_(&n, &d__1, mx, &i_one, ax, &i_one);
                    mv_(n, &v[(j << 8) - 256], mx);
                    d__1 = -d[j + 24];
                    daxpy_(&n, &d__1, mx, &i_one, ax, &i_one);
                    d__1 = dnrm2_(&n, ax, &i_one);
                    d[j + 49] = dlapy2_(&d[j + 49], &d__1);
                    d[j + 49] /= dlapy2_(&d[j - 1], &d[j + 24]);
                    d[j + 50] = d[j + 49];
                    first = FALSE_;
                }
                else
                {
                    first = TRUE_;
                }
            }

            /* --------------------------- */
            /* Display computed residuals. */
            /* --------------------------- */

            dmout_(nconv, 3, d, 25, -6, "Ritz values (Real,Imag) and relative residuals");
        }

        /* ----------------------------------------- */
        /* Print additional convergence information. */
        /* ----------------------------------------- */

        if (info == 1)
        {
            printf(" \n");
            printf(" Maximum number of iterations reached.\n");
            printf(" \n");
        }
        else if (info == 3)
        {
            printf(" \n");
            printf(" No shifts could be applied during implicit\n");
            printf(" Arnoldi update try increasing NCV.\n");
            printf(" \n");
        }

        printf(" \n");
        printf(" _NDRV6 \n");
        printf(" ====== \n");
        printf(" \n");
        printf(" Size of the matrix is %d\n", n);
        printf(" The number of Ritz values requested is %d\n", nev);
        printf(" The number of Arnoldi vectors generated (NCV) is %d\n", ncv);
        printf(" What portion of the spectrum: %s\n", which);
        printf(" The number of converged Ritz values is %d\n", nconv);
        printf(" The number of Implicit Arnoldi update iterations taken is %d\n", iparam[2]);
        printf(" The number of OP*x is %d\n", iparam[8]);
        printf(" The convergence criterion is %e\n", tol);
        printf(" \n");
    }

    /* ------------------------- */
    /* Done with program dndrv6. */
    /* ------------------------- */

    free(cdd);
    free(cdl);
    free(cdu);
    free(cdu2);
    free(ctemp);
    free(ipiv);
    free(d);
    free(v);
    free(ax);
    free(mx);
    free(resid);
    free(workd);
    free(workl);
    free(workev);

    return nconv < nev ? EXIT_FAILURE : EXIT_SUCCESS;
}

/* ========================================================================== */

/*     matrix vector multiplication subroutine */

void mv_(const a_int n, double *v, double *w)
{
    /* System generated locals */
    a_int i__1;

    /* Local variables */
    a_int j;

    /*     Compute the matrix vector multiplication y<---M*x */
    /*     where M is a n by n symmetric tridiagonal matrix with 4 on the */
    /*     diagonal, 1 on the subdiagonal and superdiagonal. */

    /* Parameter adjustments */
    --w;
    --v;

    w[1] = v[1] * 4. + v[2] * 1.;
    i__1 = n - 1;
    for (j = 2; j <= i__1; ++j)
    {
        w[j] = v[j - 1] * 1. + v[j] * 4. + v[j + 1] * 1.;
    }
    w[n] = v[n - 1] * 1. + v[n] * 4.;
} /* mv_ */

/* ------------------------------------------------------------------ */
void av_(const a_int n, double *v, double *w)
{
    /* System generated locals */
    a_int i__1;

    /* Local variables */
    a_int j;

    /*     Compute the matrix vector multiplication y<---A*x */
    /*     where M is a n by n symmetric tridiagonal matrix with 2 on the */
    /*     diagonal, -2 on the subdiagonal and 3 on the superdiagonal. */

    /* Parameter adjustments */
    --w;
    --v;

    w[1] = v[1] * 2. + v[2] * 3.;
    i__1 = n - 1;
    for (j = 2; j <= i__1; ++j)
    {
        w[j] = v[j - 1] * -2. + v[j] * 2. + v[j + 1] * 3.;
    }
    w[n] = v[n - 1] * -2. + v[n] * 2.;
} /* av_ */
