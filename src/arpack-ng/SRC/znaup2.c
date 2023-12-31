/* SRC\znaup2.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

/* Table of constant values */

static a_int i_one = 1;
static a_int i_zero = 0;
static a_int i_three = 3;
static a_bool b_true = TRUE_;
static a_int i_two = 2;
/**
 * \BeginDoc
 *
 * \Name: znaup2
 *
 * \Description:
 *  Intermediate level interface called by znaupd .
 *
 * \Usage:
 *  call znaup2
 *     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,
 *       ISHIFT, MXITER, V, LDV, H, LDH, RITZ, BOUNDS,
 *       Q, LDQ, WORKL, IPNTR, WORKD, RWORK, INFO )
 *
 * \Arguments
 *
 *  IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in znaupd .
 *  MODE, ISHIFT, MXITER: see the definition of IPARAM in znaupd .
 *
 *  NP      Integer.  (INPUT/OUTPUT)
 *          Contains the number of implicit shifts to apply during
 *          each Arnoldi iteration.
 *          If ISHIFT=1, NP is adjusted dynamically at each iteration
 *          to accelerate convergence and prevent stagnation.
 *          This is also roughly equal to the number of matrix-vector
 *          products (involving the operator OP) per Arnoldi iteration.
 *          The logic for adjusting is contained within the current
 *          subroutine.
 *          If ISHIFT=0, NP is the number of shifts the user needs
 *          to provide via reverse communication. 0 < NP < NCV-NEV.
 *          NP may be less than NCV-NEV since a leading block of the current
 *          upper Hessenberg matrix has split off and contains "unwanted"
 *          Ritz values.
 *          Upon termination of the IRA iteration, NP contains the number
 *          of "converged" wanted Ritz values.
 *
 *  IUPD    Integer.  (INPUT)
 *          IUPD .EQ. 0: use explicit restart instead implicit update.
 *          IUPD .NE. 0: use implicit update.
 *
 *  V       Complex*16  N by (NEV+NP) array.  (INPUT/OUTPUT)
 *          The Arnoldi basis vectors are returned in the first NEV
 *          columns of V.
 *
 *  LDV     Integer.  (INPUT)
 *          Leading dimension of V exactly as declared in the calling
 *          program.
 *
 *  H       Complex*16  (NEV+NP) by (NEV+NP) array.  (OUTPUT)
 *          H is used to store the generated upper Hessenberg matrix
 *
 *  LDH     Integer.  (INPUT)
 *          Leading dimension of H exactly as declared in the calling
 *          program.
 *
 *  RITZ    Complex*16  array of length NEV+NP.  (OUTPUT)
 *          RITZ(1:NEV)  contains the computed Ritz values of OP.
 *
 *  BOUNDS  Complex*16  array of length NEV+NP.  (OUTPUT)
 *          BOUNDS(1:NEV) contain the error bounds corresponding to
 *          the computed Ritz values.
 *
 *  Q       Complex*16  (NEV+NP) by (NEV+NP) array.  (WORKSPACE)
 *          Private (replicated) work array used to accumulate the
 *          rotation in the shift application step.
 *
 *  LDQ     Integer.  (INPUT)
 *          Leading dimension of Q exactly as declared in the calling
 *          program.
 *
 *  WORKL   Complex*16  work array of length at least
 *          (NEV+NP)**2 + 3*(NEV+NP).  (WORKSPACE)
 *          Private (replicated) array on each PE or array allocated on
 *          the front end.  It is used in shifts calculation, shifts
 *          application and convergence checking.
 *
 *  IPNTR   Integer array of length 3.  (OUTPUT)
 *          Pointer to mark the starting locations in the WORKD for
 *          vectors used by the Arnoldi iteration.
 *          -------------------------------------------------------------
 *          IPNTR(1): pointer to the current operand vector X.
 *          IPNTR(2): pointer to the current result vector Y.
 *          IPNTR(3): pointer to the vector B * X when used in the
 *                    shift-and-invert mode.  X is the current operand.
 *          -------------------------------------------------------------
 *
 *  WORKD   Complex*16  work array of length 3*N.  (WORKSPACE)
 *          Distributed array to be used in the basic Arnoldi iteration
 *          for reverse communication.  The user should not use WORKD
 *          as temporary workspace during the iteration !!!!!!!!!!
 *          See Data Distribution Note in ZNAUPD .
 *
 *  RWORK   Double precision    work array of length  NEV+NP ( WORKSPACE)
 *          Private (replicated) array on each PE or array allocated on
 *          the front end.
 *
 *  INFO    Integer.  (INPUT/OUTPUT)
 *          If INFO .EQ. 0, a randomly initial residual vector is used.
 *          If INFO .NE. 0, RESID contains the initial residual vector,
 *                          possibly from a previous run.
 *          Error flag on output.
 *          =     0: Normal return.
 *          =     1: Maximum number of iterations taken.
 *                   All possible eigenvalues of OP has been found.
 *                   NP returns the number of converged Ritz values.
 *          =     2: No shifts could be applied.
 *          =    -8: Error return from LAPACK eigenvalue calculation;
 *                   This should never happen.
 *          =    -9: Starting vector is zero.
 *          = -9999: Could not build an Arnoldi factorization.
 *                   Size that was built in returned in NP.
 *
 * \EndDoc
 *
 * -----------------------------------------------------------------------
 *
 * \BeginLib
 *
 * \Local variables:
 *     xxxxxx  Complex*16
 *
 * \References:
 *  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
 *     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
 *     pp 357-385.
 *  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
 *     Restarted Arnoldi Iteration", Rice University Technical Report
 *     TR95-13, Department of Computational and Applied Mathematics.
 *
 * \Routines called:
 *     zgetv0   ARPACK initial vector generation routine.
 *     znaitr   ARPACK Arnoldi factorization routine.
 *     znapps   ARPACK application of implicit shifts routine.
 *     zneigh   ARPACK compute Ritz values and error bounds routine.
 *     zngets   ARPACK reorder Ritz values and error bounds routine.
 *     zsortc   ARPACK sorting routine.
 *     ivout   ARPACK utility routine that prints integers.
 *     arscnd  ARPACK utility routine for timing.
 *     zmout    ARPACK utility routine that prints matrices
 *     zvout    ARPACK utility routine that prints vectors.
 *     dvout    ARPACK utility routine that prints vectors.
 *     dlamch   LAPACK routine that determines machine constants.
 *     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     zcopy    Level 1 BLAS that copies one vector to another .
 *     zdotc    Level 1 BLAS that computes the scalar product of two vectors.
 *     zswap    Level 1 BLAS that swaps two vectors.
 *     dznrm2   Level 1 BLAS that computes the norm of a vector.
 *
 * \Author
 *     Danny Sorensen               Phuong Vu
 *     Richard Lehoucq              CRPC / Rice Universitya
 *     Chao Yang                    Houston, Texas
 *     Dept. of Computational &
 *     Applied Mathematics
 *     Rice University
 *     Houston, Texas
 *
 * \SCCS Information: @(#)
 * FILE: naup2.F   SID: 2.6   DATE OF SID: 06/01/00   RELEASE: 2
 *
 * \Remarks
 *     1. None
 *
 * \EndLib
 */
int znaup2_(a_int *ido, const char *bmat, a_int *n, const char *which, a_int *nev, a_int *np,
     double *tol, a_dcomplex *resid, a_int *mode, a_int *iupd, a_int *ishift, a_int *mxiter,
     a_dcomplex *v, a_int *ldv, a_dcomplex *h, a_int *ldh, a_dcomplex *ritz, a_dcomplex *bounds,
     a_dcomplex *q, a_int *ldq, a_dcomplex *workl, a_int *ipntr, a_dcomplex *workd,
     double *rwork, a_int *info)
{
    /* System generated locals */
    a_int h_dim1, h_offset, q_dim1, q_offset, v_dim1, v_offset, i__1, i__2, i__3;
    double d__1, d__2, d__3, d__4;
    a_dcomplex z__1;

    /* Local variables */
    a_int i, j;
    static float t0, t1, t2, t3;
    a_int kp[3];
    static a_int np0, nev0;
    static double eps23;
    a_int ierr;
    static a_int iter;
    static a_bool getv0, cnorm;
    static a_int nconv;
    double rtemp;
    static a_bool initv;
    static double rnorm;
    static a_int nevbef;
    static a_bool update, ushift;
    static a_int kplusp, msglvl;
    a_int nptemp;
    char wprime[3];
    a_dcomplex cmpnorm;

    /* Parameter adjustments */
    --workd;
    --resid;
    --rwork;
    --workl;
    --bounds;
    --ritz;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h -= h_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --ipntr;

    if (*ido == 0)
    {

#ifndef NO_TIMER
        arscnd_(&t0);
#endif

        msglvl = debug_1.mcaup2;

        nev0 = *nev;
        np0 = *np;

        /* ----------------------------------- */
        /* kplusp is the bound on the largest  */
        /*        Lanczos factorization built. */
        /* nconv is the current number of      */
        /*        "converged" eigenvalues.     */
        /* iter is the counter on the current  */
        /*      iteration step.                */
        /* ----------------------------------- */

        kplusp = *nev + *np;
        nconv = 0;
        iter = 0;

        /* ------------------------------- */
        /* Get machine dependent constant. */
        /* ------------------------------- */

        eps23 = dlamch_("E");
        eps23 = pow(eps23, TWO_THIRDS);

        /* ------------------------------------- */
        /* Set flags for computing the first NEV */
        /* steps of the Arnoldi factorization.   */
        /* ------------------------------------- */

        getv0 = TRUE_;
        update = FALSE_;
        ushift = FALSE_;
        cnorm = FALSE_;

        if (*info != 0)
        {

            /* ------------------------------------------ */
            /* User provides the initial residual vector. */
            /* ------------------------------------------ */

            initv = TRUE_;
            *info = 0;
        }
        else
        {
            initv = FALSE_;
        }
    }

    /* ------------------------------------------- */
    /* Get a possibly random starting vector and   */
    /* force it into the range of the operator OP. */
    /* ------------------------------------------- */

    /* L10: */

    if (getv0)
    {
        zgetv0_(ido, bmat, &i_one, &initv, n, &i_one, &v[v_offset], ldv, &resid[1], &rnorm, &ipntr[1], &workd[1], info);

        if (*ido != 99)
        {
            goto L9000;
        }

        if (rnorm == 0.)
        {

            /* --------------------------------------- */
            /* The initial vector is zero. Error exit. */
            /* --------------------------------------- */

            *info = -9;
            goto L1100;
        }
        getv0 = FALSE_;
        *ido = 0;
    }

    /* --------------------------------- */
    /* Back from reverse communication : */
    /* continue with update step         */
    /* --------------------------------- */

    if (update)
    {
        goto L20;
    }

    /* ----------------------------------------- */
    /* Back from computing user specified shifts */
    /* ----------------------------------------- */

    if (ushift)
    {
        goto L50;
    }

    /* ----------------------------------- */
    /* Back from computing residual norm   */
    /* at the end of the current iteration */
    /* ----------------------------------- */

    if (cnorm)
    {
        goto L100;
    }

    /* -------------------------------------------------------- */
    /* Compute the first NEV steps of the Arnoldi factorization */
    /* -------------------------------------------------------- */

    znaitr_(ido, bmat, n, &i_zero, nev, mode, &resid[1], &rnorm, &v[v_offset], ldv, &h[h_offset], ldh, &ipntr[1], &workd[1], info);

    if (*ido != 99)
    {
        goto L9000;
    }

    if (*info > 0)
    {
        *np = *info;
        *mxiter = iter;
        *info = -9999;
        goto L1200;
    }

    /* ------------------------------------------------------------ */
    /*                                                              */
    /*           M A I N  ARNOLDI  I T E R A T I O N  L O O P       */
    /*           Each iteration implicitly restarts the Arnoldi     */
    /*           factorization in place.                            */
    /*                                                              */
    /* ------------------------------------------------------------ */

L1000:

    ++iter;

#ifndef NO_TRACE
    if (msglvl > 0)
    {
        ivout_(1, &iter, debug_1.ndigit, "_naup2: **** Start of major iteration number ****");
    }
#endif

    /* --------------------------------------------------------- */
    /* Compute NP additional steps of the Arnoldi factorization. */
    /* Adjust NP since NEV might have been updated by last call  */
    /* to the shift application routine znapps .                  */
    /* --------------------------------------------------------- */

    *np = kplusp - *nev;

#ifndef NO_TRACE
    if (msglvl > 1)
    {
        ivout_(1, nev, debug_1.ndigit, "_naup2: The length of the current Arnoldi factorization");
        ivout_(1, np, debug_1.ndigit, "_naup2: Extend the Arnoldi factorization by");
    }
#endif

    /* --------------------------------------------------------- */
    /* Compute NP additional steps of the Arnoldi factorization. */
    /* --------------------------------------------------------- */

    *ido = 0;
L20:
    update = TRUE_;

    znaitr_(ido, bmat, n, nev, np, mode, &resid[1], &rnorm, &v[v_offset], ldv, &h[h_offset], ldh, &ipntr[1], &workd[1], info);

    if (*ido != 99)
    {
        goto L9000;
    }

    if (*info > 0)
    {
        *np = *info;
        *mxiter = iter;
        *info = -9999;
        goto L1200;
    }
    update = FALSE_;

#ifndef NO_TRACE
    if (msglvl > 1)
    {
        dvout_(1, &rnorm, debug_1.ndigit, "_naup2: Corresponding B-norm of the residual");
    }
#endif

    /* ------------------------------------------------------ */
    /* Compute the eigenvalues and corresponding error bounds */
    /* of the current upper Hessenberg matrix.                */
    /* ------------------------------------------------------ */

    zneigh_(&rnorm, &kplusp, &h[h_offset], ldh, &ritz[1], &bounds[1], &q[q_offset], ldq, &workl[1], &rwork[1], &ierr);

    if (ierr != 0)
    {
        *info = -8;
        goto L1200;
    }

    /* ------------------------------------------------- */
    /* Select the wanted Ritz values and their bounds    */
    /* to be used in the convergence test.               */
    /* The wanted part of the spectrum and corresponding */
    /* error bounds are in the last NEV loc. of RITZ,    */
    /* and BOUNDS respectively.                          */
    /* ------------------------------------------------- */

    *nev = nev0;
    *np = np0;

    /* ------------------------------------------------ */
    /* Make a copy of Ritz values and the corresponding */
    /* Ritz estimates obtained from zneigh .             */
    /* ------------------------------------------------ */

    /* Computing 2nd power */
    i__1 = kplusp;
    zcopy_(&kplusp, &ritz[1], &i_one, &workl[i__1 * i__1 + 1], &i_one);
    /* Computing 2nd power */
    i__1 = kplusp;
    zcopy_(&kplusp, &bounds[1], &i_one, &workl[i__1 * i__1 + kplusp + 1], &i_one);

    /* ------------------------------------------------- */
    /* Select the wanted Ritz values and their bounds    */
    /* to be used in the convergence test.               */
    /* The wanted part of the spectrum and corresponding */
    /* bounds are in the last NEV loc. of RITZ           */
    /* BOUNDS respectively.                              */
    /* ------------------------------------------------- */

    zngets_(ishift, which, nev, np, &ritz[1], &bounds[1]);

    /* ---------------------------------------------------------- */
    /* Convergence test: currently we use the following criteria. */
    /* The relative accuracy of a Ritz value is considered        */
    /* acceptable if:                                             */
    /*                                                            */
    /* error_bounds(i) .le. tol*max(eps23, magnitude_of_ritz(i)). */
    /*                                                            */
    /* ---------------------------------------------------------- */

    nconv = 0;

    i__1 = *nev;
    for (i = 1; i <= i__1; ++i)
    {
        /* Computing MAX */
        i__2 = *np + i;
        d__3 = ritz[i__2].r;
        d__4 = ritz[*np + i].i;
        d__1 = eps23, d__2 = dlapy2_(&d__3, &d__4);
        rtemp = max(d__1, d__2);
        i__2 = *np + i;
        d__1 = bounds[i__2].r;
        d__2 = bounds[*np + i].i;
        if (dlapy2_(&d__1, &d__2) <= *tol * rtemp)
        {
            ++nconv;
        }
        /* L25: */
    }

#ifndef NO_TRACE
    if (msglvl > 2)
    {
        kp[0] = *nev;
        kp[1] = *np;
        kp[2] = nconv;
        ivout_(3, kp, debug_1.ndigit, "_naup2: NEV, NP, NCONV are");
        zvout_(kplusp, &ritz[1], debug_1.ndigit, "_naup2: The eigenvalues of H");
        zvout_(kplusp, &bounds[1], debug_1.ndigit, "_naup2: Ritz estimates of the current NCV Ritz values");
    }
#endif

    /* ------------------------------------------------------- */
    /* Count the number of unwanted Ritz values that have zero */
    /* Ritz estimates. If any Ritz estimates are equal to zero */
    /* then a leading block of H of order equal to at least    */
    /* the number of Ritz values with zero Ritz estimates has  */
    /* split off. None of these Ritz values may be removed by  */
    /* shifting. Decrease NP the number of shifts to apply. If */
    /* no shifts may be applied, then prepare to exit          */
    /* ------------------------------------------------------- */

    nptemp = *np;
    i__1 = nptemp;
    for (j = 1; j <= i__1; ++j)
    {
        i__2 = j;
        if (bounds[i__2].r == 0. && bounds[i__2].i == 0.)
        {
            --(*np);
            ++(*nev);
        }
        /* L30: */
    }

    if (nconv >= nev0 || iter > *mxiter || *np == 0)
    {

#ifndef NO_TRACE
        if (msglvl > 4)
        {
            /* Computing 2nd power */
            i__1 = kplusp;
            zvout_(kplusp, &workl[i__1 * i__1 + 1], debug_1.ndigit, "_naup2: Eigenvalues computed by _neigh:");
            /* Computing 2nd power */
            i__1 = kplusp;
            zvout_(kplusp, &workl[i__1 * i__1 + kplusp + 1], debug_1.ndigit, "_naup2: Ritz estimates computed by _neigh:");
        }
#endif

        /* ---------------------------------------------- */
        /* Prepare to exit. Put the converged Ritz values */
        /* and corresponding bounds in RITZ(1:NCONV) and  */
        /* BOUNDS(1:NCONV) respectively. Then sort. Be    */
        /* careful when NCONV > NP                        */
        /* ---------------------------------------------- */

        /* ---------------------------------------- */
        /*  Use h( 3,1 ) as storage to communicate  */
        /*  rnorm to zneupd  if needed               */
        /* ---------------------------------------- */
        i__1 = h_dim1 + 3;
        z__1.r = rnorm, z__1.i = 0.;
        h[i__1].r = z__1.r, h[i__1].i = z__1.i;

        /* -------------------------------------------- */
        /* Sort Ritz values so that converged Ritz      */
        /* values appear within the first NEV locations */
        /* of ritz and bounds, and the most desired one */
        /* appears at the front.                        */
        /* -------------------------------------------- */

        if (strcmp(which, "LM") == 0)
        {
            strcpy(wprime, "SM");
        }
        if (strcmp(which, "SM") == 0)
        {
            strcpy(wprime, "LM");
        }
        if (strcmp(which, "LR") == 0)
        {
            strcpy(wprime, "SR");
        }
        if (strcmp(which, "SR") == 0)
        {
            strcpy(wprime, "LR");
        }
        if (strcmp(which, "LI") == 0)
        {
            strcpy(wprime, "SI");
        }
        if (strcmp(which, "SI") == 0)
        {
            strcpy(wprime, "LI");
        }

        zsortc_(wprime, &b_true, &kplusp, &ritz[1], &bounds[1]);

        /* ------------------------------------------------ */
        /* Scale the Ritz estimate of each Ritz value       */
        /* by 1 / max(eps23, magnitude of the Ritz value).  */
        /* ------------------------------------------------ */

        i__1 = nev0;
        for (j = 1; j <= i__1; ++j)
        {
            /* Computing MAX */
            i__2 = j;
            d__3 = ritz[i__2].r;
            d__4 = ritz[j].i;
            d__1 = eps23, d__2 = dlapy2_(&d__3, &d__4);
            rtemp = max(d__1, d__2);
            i__2 = j;
            i__3 = j;
            z__1.r = bounds[i__3].r / rtemp, z__1.i = bounds[i__3].i / rtemp;
            bounds[i__2].r = z__1.r, bounds[i__2].i = z__1.i;
            /* L35: */
        }

        /* ------------------------------------------------- */
        /* Sort the Ritz values according to the scaled Ritz */
        /* estimates.  This will push all the converged ones */
        /* towards the front of ritz, bounds (in the case    */
        /* when NCONV < NEV.)                                */
        /* ------------------------------------------------- */

        strcpy(wprime, "LM");
        zsortc_(wprime, &b_true, &nev0, &bounds[1], &ritz[1]);

        /* -------------------------------------------- */
        /* Scale the Ritz estimate back to its original */
        /* value.                                       */
        /* -------------------------------------------- */

        i__1 = nev0;
        for (j = 1; j <= i__1; ++j)
        {
            /* Computing MAX */
            i__2 = j;
            d__3 = ritz[i__2].r;
            d__4 = ritz[j].i;
            d__1 = eps23, d__2 = dlapy2_(&d__3, &d__4);
            rtemp = max(d__1, d__2);
            i__2 = j;
            i__3 = j;
            z__1.r = rtemp * bounds[i__3].r, z__1.i = rtemp * bounds[i__3].i;
            bounds[i__2].r = z__1.r, bounds[i__2].i = z__1.i;
            /* L40: */
        }

        /* --------------------------------------------- */
        /* Sort the converged Ritz values again so that  */
        /* the "threshold" value appears at the front of */
        /* ritz and bound.                               */
        /* --------------------------------------------- */

        zsortc_(which, &b_true, &nconv, &ritz[1], &bounds[1]);

#ifndef NO_TRACE
        if (msglvl > 1)
        {
            zvout_(kplusp, &ritz[1], debug_1.ndigit, "_naup2: Sorted eigenvalues");
            zvout_(kplusp, &bounds[1], debug_1.ndigit, "_naup2: Sorted ritz estimates.");
        }
#endif

        /* ---------------------------------- */
        /* Max iterations have been exceeded. */
        /* ---------------------------------- */

        if (iter > *mxiter && nconv < nev0)
        {
            *info = 1;
        }

        /* ------------------- */
        /* No shifts to apply. */
        /* ------------------- */

        if (*np == 0 && nconv < nev0)
        {
            *info = 2;
        }

        *np = nconv;
        goto L1100;
    }
    else if (nconv < nev0 && *ishift == 1)
    {

        /* ----------------------------------------------- */
        /* Do not have all the requested eigenvalues yet.  */
        /* To prevent possible stagnation, adjust the size */
        /* of NEV.                                         */
        /* ----------------------------------------------- */

        nevbef = *nev;
        /* Computing MIN */
        i__1 = nconv, i__2 = *np / 2;
        *nev += min(i__1, i__2);
        if (*nev == 1 && kplusp >= 6)
        {
            *nev = kplusp / 2;
        }
        else if (*nev == 1 && kplusp > 3)
        {
            *nev = 2;
        }
        *np = kplusp - *nev;

        /* ------------------------------------- */
        /* If the size of NEV was just increased */
        /* resort the eigenvalues.               */
        /* ------------------------------------- */

        if (nevbef < *nev)
        {
            zngets_(ishift, which, nev, np, &ritz[1], &bounds[1]);
        }
    }

#ifndef NO_TRACE
    if (msglvl > 0)
    {
        ivout_(1, &nconv, debug_1.ndigit, "_naup2: no. of \"converged\" Ritz values at this iter.");
        if (msglvl > 1)
        {
            kp[0] = *nev;
            kp[1] = *np;
            ivout_(2, kp, debug_1.ndigit, "_naup2: NEV and NP are");
            zvout_(*nev, &ritz[*np + 1], debug_1.ndigit, "_naup2: \"wanted\" Ritz values ");
            zvout_(*nev, &bounds[*np + 1], debug_1.ndigit, "_naup2: Ritz estimates of the \"wanted\" values ");
        }
    }
#endif

    if (*ishift == 0)
    {

        /* ----------------------------------------------------- */
        /* User specified shifts: pop back out to get the shifts */
        /* and return them in the first 2*NP locations of WORKL. */
        /* ----------------------------------------------------- */

        ushift = TRUE_;
        *ido = 3;
        goto L9000;
    }
L50:
    ushift = FALSE_;

    if (*ishift != 1)
    {

        /* -------------------------------- */
        /* Move the NP shifts from WORKL to */
        /* RITZ, to free up WORKL           */
        /* for non-exact shift case.        */
        /* -------------------------------- */

        zcopy_(np, &workl[1], &i_one, &ritz[1], &i_one);
    }

#ifndef NO_TRACE
    if (msglvl > 2)
    {
        ivout_(1, np, debug_1.ndigit, "_naup2: The number of shifts to apply ");
        zvout_(*np, &ritz[1], debug_1.ndigit, "_naup2: values of the shifts");
        if (*ishift == 1)
        {
            zvout_(*np, &bounds[1], debug_1.ndigit, "_naup2: Ritz estimates of the shifts");
        }
    }
#endif

    /* ------------------------------------------------------- */
    /* Apply the NP implicit shifts by QR bulge chasing.       */
    /* Each shift is applied to the whole upper Hessenberg     */
    /* matrix H.                                               */
    /* The first 2*N locations of WORKD are used as workspace. */
    /* ------------------------------------------------------- */

    znapps_(n, nev, np, &ritz[1], &v[v_offset], ldv, &h[h_offset], ldh, &resid[1], &q[q_offset], ldq, &workl[1], &workd[1]);

    /* ------------------------------------------- */
    /* Compute the B-norm of the updated residual. */
    /* Keep B*RESID in WORKD(1:N) to be used in    */
    /* the first step of the next call to znaitr .  */
    /* ------------------------------------------- */

    cnorm = TRUE_;
#ifndef NO_TIMER
    arscnd_(&t2);
#endif

    if (*bmat == 'G')
    {
        ++timing_1.nbx;
        zcopy_(n, &resid[1], &i_one, &workd[*n + 1], &i_one);
        ipntr[1] = *n + 1;
        ipntr[2] = 1;
        *ido = 2;

        /* -------------------------------- */
        /* Exit in order to compute B*RESID */
        /* -------------------------------- */

        goto L9000;
    }
    else if (*bmat == 'I')
    {
        zcopy_(n, &resid[1], &i_one, &workd[1], &i_one);
    }

L100:

    /* -------------------------------- */
    /* Back from reverse communication; */
    /* WORKD(1:N) := B*RESID            */
    /* -------------------------------- */

#ifndef NO_TIMER
    if (*bmat == 'G')
    {
        arscnd_(&t3);
        timing_1.tmvbx += t3 - t2;
    }
#endif

    if (*bmat == 'G')
    {
        zzdotc_(&z__1, n, &resid[1], &i_one, &workd[1], &i_one);
        cmpnorm.r = z__1.r, cmpnorm.i = z__1.i;
        d__1 = cmpnorm.r;
        d__2 = cmpnorm.i;
        rnorm = sqrt(dlapy2_(&d__1, &d__2));
    }
    else if (*bmat == 'I')
    {
        rnorm = dznrm2_(n, &resid[1], &i_one);
    }
    cnorm = FALSE_;

#ifndef NO_TRACE
    if (msglvl > 2)
    {
        dvout_(1, &rnorm, debug_1.ndigit, "_naup2: B-norm of residual for compressed factorization");
        zmout_(*nev, *nev, &h[h_offset], *ldh, debug_1.ndigit, "_naup2: Compressed upper Hessenberg matrix H");
    }
#endif

    goto L1000;

    /* ------------------------------------------------------------- */
    /*                                                               */
    /*  E N D     O F     M A I N     I T E R A T I O N     L O O P  */
    /*                                                               */
    /* ------------------------------------------------------------- */

L1100:

    *mxiter = iter;
    *nev = nconv;

L1200:
    *ido = 99;

    /* ---------- */
    /* Error Exit */
    /* ---------- */

#ifndef NO_TIMER
    arscnd_(&t1);
    timing_1.tcaup2 = t1 - t0;
#endif

L9000:

    /* ------------- */
    /* End of znaup2  */
    /* ------------- */

    return 0;
} /* znaup2_ */
