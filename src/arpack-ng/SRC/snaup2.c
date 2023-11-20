/* SRC\snaup2.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

/* Table of constant values */

static a_int i_one = 1;
static a_int i_zero = 0;
static a_int i_four = 4;
static a_bool b_true = TRUE_;
static a_int i_two = 2;
/**
 * \BeginDoc
 *
 * \Name: snaup2
 *
 * \Description:
 *  Intermediate level interface called by snaupd.
 *
 * \Usage:
 *  call snaup2
 *     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,
 *       ISHIFT, MXITER, V, LDV, H, LDH, RITZR, RITZI, BOUNDS,
 *       Q, LDQ, WORKL, IPNTR, WORKD, INFO )
 *
 * \Arguments
 *
 *  IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in snaupd.
 *  MODE, ISHIFT, MXITER: see the definition of IPARAM in snaupd.
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
 *          NP may be less than NCV-NEV for two reasons. The first, is
 *          to keep complex conjugate pairs of "wanted" Ritz values
 *          together. The second, is that a leading block of the current
 *          upper Hessenberg matrix has split off and contains "unwanted"
 *          Ritz values.
 *          Upon termination of the IRA iteration, NP contains the number
 *          of "converged" wanted Ritz values.
 *
 *  IUPD    Integer.  (INPUT)
 *          IUPD .EQ. 0: use explicit restart instead implicit update.
 *          IUPD .NE. 0: use implicit update.
 *
 *  V       Real  N by (NEV+NP) array.  (INPUT/OUTPUT)
 *          The Arnoldi basis vectors are returned in the first NEV
 *          columns of V.
 *
 *  LDV     Integer.  (INPUT)
 *          Leading dimension of V exactly as declared in the calling
 *          program.
 *
 *  H       Real  (NEV+NP) by (NEV+NP) array.  (OUTPUT)
 *          H is used to store the generated upper Hessenberg matrix
 *
 *  LDH     Integer.  (INPUT)
 *          Leading dimension of H exactly as declared in the calling
 *          program.
 *
 *  RITZR,  Real  arrays of length NEV+NP.  (OUTPUT)
 *  RITZI   RITZR(1:NEV) (resp. RITZI(1:NEV)) contains the real (resp.
 *          imaginary) part of the computed Ritz values of OP.
 *
 *  BOUNDS  Real  array of length NEV+NP.  (OUTPUT)
 *          BOUNDS(1:NEV) contain the error bounds corresponding to
 *          the computed Ritz values.
 *
 *  Q       Real  (NEV+NP) by (NEV+NP) array.  (WORKSPACE)
 *          Private (replicated) work array used to accumulate the
 *          rotation in the shift application step.
 *
 *  LDQ     Integer.  (INPUT)
 *          Leading dimension of Q exactly as declared in the calling
 *          program.
 *
 *  WORKL   Real  work array of length at least
 *          (NEV+NP)**2 + 3*(NEV+NP).  (INPUT/WORKSPACE)
 *          Private (replicated) array on each PE or array allocated on
 *          the front end.  It is used in shifts calculation, shifts
 *          application and convergence checking.
 *
 *          On exit, the last 3*(NEV+NP) locations of WORKL contain
 *          the Ritz values (real,imaginary) and associated Ritz
 *          estimates of the current Hessenberg matrix.  They are
 *          listed in the same order as returned from sneigh.
 *
 *          If ISHIFT .EQ. O and IDO .EQ. 3, the first 2*NP locations
 *          of WORKL are used in reverse communication to hold the user
 *          supplied shifts.
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
 *  WORKD   Real  work array of length 3*N.  (WORKSPACE)
 *          Distributed array to be used in the basic Arnoldi iteration
 *          for reverse communication.  The user should not use WORKD
 *          as temporary workspace during the iteration !!!!!!!!!!
 *          See Data Distribution Note in DNAUPD.
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
 *     xxxxxx  real
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
 *     sgetv0  ARPACK initial vector generation routine.
 *     snaitr  ARPACK Arnoldi factorization routine.
 *     snapps  ARPACK application of implicit shifts routine.
 *     snconv  ARPACK convergence of Ritz values routine.
 *     sneigh  ARPACK compute Ritz values and error bounds routine.
 *     sngets  ARPACK reorder Ritz values and error bounds routine.
 *     ssortc  ARPACK sorting routine.
 *     ivout   ARPACK utility routine that prints integers.
 *     arscnd  ARPACK utility routine for timing.
 *     smout   ARPACK utility routine that prints matrices
 *     svout   ARPACK utility routine that prints vectors.
 *     slamch  LAPACK routine that determines machine constants.
 *     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     scopy   Level 1 BLAS that copies one vector to another .
 *     sdot    Level 1 BLAS that computes the scalar product of two vectors.
 *     snrm2   Level 1 BLAS that computes the norm of a vector.
 *     sswap   Level 1 BLAS that swaps two vectors.
 *
 * \Author
 *     Danny Sorensen               Phuong Vu
 *     Richard Lehoucq              CRPC / Rice University
 *     Dept. of Computational &     Houston, Texas
 *     Applied Mathematics
 *     Rice University
 *     Houston, Texas
 *
 * \SCCS Information: @(#)
 * FILE: naup2.F   SID: 2.8   DATE OF SID: 10/17/00   RELEASE: 2
 *
 * \Remarks
 *     1. None
 *
 * \EndLib
 */
int snaup2_(a_int *ido, const char *bmat, a_int *n, const char *which, a_int *nev, a_int *np,
     float *tol, float *resid, a_int *mode, a_int *iupd, a_int *ishift, a_int *mxiter,
     float *v, a_int *ldv, float *h, a_int *ldh, float *ritzr, float *ritzi, float *bounds,
     float *q, a_int *ldq, float *workl, a_int *ipntr, float *workd, a_int *info)
{
    /* System generated locals */
    a_int h_dim1, h_offset, q_dim1, q_offset, v_dim1, v_offset, i__1, i__2;
    float r__1, r__2;
    double d__1;

    /* Local variables */
    a_int j;
    static float t0, t1, t2, t3;
    a_int kp[4];
    static a_int np0, nev0;
    static float eps23;
    a_int ierr;
    static a_int iter;
    float temp;
    static a_bool getv0;
    static a_bool cnorm;
    static a_int nconv;
    static a_bool initv;
    static float rnorm;
    static a_int nevbef;
    static a_bool update;
    char wprime[2];
    static a_bool ushift;
    static a_int kplusp, msglvl;
    a_int nptemp;
    static a_int numcnv;

    /* Parameter adjustments */
    --workd;
    --resid;
    --workl;
    --bounds;
    --ritzi;
    --ritzr;
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

        arscnd_(&t0);

        msglvl = debug_1.mnaup2;

        /* ----------------------------------- */
        /* Get the machine dependent constant. */
        /* ----------------------------------- */

        eps23 = slamch_("E");
        eps23 = pow((double)eps23, TWO_THIRDS);

        nev0 = *nev;
        np0 = *np;

        /* ----------------------------------- */
        /* kplusp is the bound on the largest  */
        /*        Lanczos factorization built. */
        /* nconv is the current number of      */
        /*        "converged" eigenvlues.      */
        /* iter is the counter on the current  */
        /*      iteration step.                */
        /* ----------------------------------- */

        kplusp = *nev + *np;
        nconv = 0;
        iter = 0;

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
        sgetv0_(ido, bmat, &i_one, &initv, n, &i_one, &v[v_offset], ldv, &resid[1], &rnorm, &ipntr[1], &workd[1], info);

        if (*ido != 99)
        {
            goto L9000;
        }

        if (rnorm == 0.f)
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

    snaitr_(ido, bmat, n, &i_zero, nev, mode, &resid[1], &rnorm, &v[v_offset], ldv, &h[h_offset], ldh, &ipntr[1], &workd[1], info);

    /* ------------------------------------------------- */
    /* ido .ne. 99 implies use of reverse communication  */
    /* to compute operations involving OP and possibly B */
    /* ------------------------------------------------- */

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

    if (msglvl > 0)
    {
        ivout_(1, &iter, debug_1.ndigit, "_naup2: **** Start of major iteration number ****");
    }

    /* --------------------------------------------------------- */
    /* Compute NP additional steps of the Arnoldi factorization. */
    /* Adjust NP since NEV might have been updated by last call  */
    /* to the shift application routine snapps.                  */
    /* --------------------------------------------------------- */

    *np = kplusp - *nev;

    if (msglvl > 1)
    {
        ivout_(1, nev, debug_1.ndigit, "_naup2: The length of the current Arnoldi factorization");
        ivout_(1, np, debug_1.ndigit, "_naup2: Extend the Arnoldi factorization by");
    }

    /* --------------------------------------------------------- */
    /* Compute NP additional steps of the Arnoldi factorization. */
    /* --------------------------------------------------------- */

    *ido = 0;
L20:
    update = TRUE_;

    snaitr_(ido, bmat, n, nev, np, mode, &resid[1], &rnorm, &v[v_offset], ldv, &h[h_offset], ldh, &ipntr[1], &workd[1], info);

    /* ------------------------------------------------- */
    /* ido .ne. 99 implies use of reverse communication  */
    /* to compute operations involving OP and possibly B */
    /* ------------------------------------------------- */

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

    if (msglvl > 1)
    {
        svout_(1, &rnorm, debug_1.ndigit, "_naup2: Corresponding B-norm of the residual");
    }

    /* ------------------------------------------------------ */
    /* Compute the eigenvalues and corresponding error bounds */
    /* of the current upper Hessenberg matrix.                */
    /* ------------------------------------------------------ */

    sneigh_(&rnorm, &kplusp, &h[h_offset], ldh, &ritzr[1], &ritzi[1], &bounds[1], &q[q_offset], ldq, &workl[1], &ierr);

    if (ierr != 0)
    {
        *info = -8;
        goto L1200;
    }

    /* -------------------------------------------------- */
    /* Make a copy of eigenvalues and corresponding error */
    /* bounds obtained from sneigh.                       */
    /* -------------------------------------------------- */

    /* Computing 2nd power */
    i__1 = kplusp;
    scopy_(&kplusp, &ritzr[1], &i_one, &workl[i__1 * i__1 + 1], &i_one);
    /* Computing 2nd power */
    i__1 = kplusp;
    scopy_(&kplusp, &ritzi[1], &i_one, &workl[i__1 * i__1 + kplusp + 1], &i_one);
    /* Computing 2nd power */
    i__1 = kplusp;
    scopy_(&kplusp, &bounds[1], &i_one, &workl[i__1 * i__1 + (kplusp << 1) + 1], &i_one);

    /* ------------------------------------------------- */
    /* Select the wanted Ritz values and their bounds    */
    /* to be used in the convergence test.               */
    /* The wanted part of the spectrum and corresponding */
    /* error bounds are in the last NEV loc. of RITZR,   */
    /* RITZI and BOUNDS respectively. The variables NEV  */
    /* and NP may be updated if the NEV-th wanted Ritz   */
    /* value has a non zero imaginary part. In this case */
    /* NEV is increased by one and NP decreased by one.  */
    /* NOTE: The last two arguments of sngets are no     */
    /* longer used as of version 2.1.                    */
    /* ------------------------------------------------- */

    *nev = nev0;
    *np = np0;
    numcnv = *nev;
    sngets_(ishift, which, nev, np, &ritzr[1], &ritzi[1], &bounds[1], &workl[1], &workl[*np + 1]);
    if (*nev == nev0 + 1)
    {
        numcnv = nev0 + 1;
    }

    /* ----------------- */
    /* Convergence test. */
    /* ----------------- */

    scopy_(nev, &bounds[*np + 1], &i_one, &workl[(*np << 1) + 1], &i_one);
    snconv_(nev, &ritzr[*np + 1], &ritzi[*np + 1], &workl[(*np << 1) + 1], tol, &nconv);

    if (msglvl > 2)
    {
        kp[0] = *nev;
        kp[1] = *np;
        kp[2] = numcnv;
        kp[3] = nconv;
        ivout_(4, kp, debug_1.ndigit, "_naup2: NEV, NP, NUMCNV, NCONV are");
        svout_(kplusp, &ritzr[1], debug_1.ndigit, "_naup2: Real part of the eigenvalues of H");
        svout_(kplusp, &ritzi[1], debug_1.ndigit, "_naup2: Imaginary part of the eigenvalues of H");
        svout_(kplusp, &bounds[1], debug_1.ndigit, "_naup2: Ritz estimates of the current NCV Ritz values");
    }

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
        if (bounds[j] == 0.f)
        {
            --(*np);
            ++(*nev);
        }
        /* L30: */
    }

    if (nconv >= numcnv || iter > *mxiter || *np == 0)
    {

        if (msglvl > 4)
        {
            /* Computing 2nd power */
            i__1 = kplusp;
            svout_(kplusp, &workl[i__1 * i__1 + 1], debug_1.ndigit, "_naup2: Real part of the eig computed by _neigh:");
            /* Computing 2nd power */
            i__1 = kplusp;
            svout_(kplusp, &workl[i__1 * i__1 + kplusp + 1], debug_1.ndigit, "_naup2: Imag part of the eig computed by _neigh:");
            /* Computing 2nd power */
            i__1 = kplusp;
            svout_(kplusp, &workl[i__1 * i__1 + (kplusp << 1) + 1], debug_1.ndigit, "_naup2: Ritz eistmates computed by _neigh:");
        }

        /* ---------------------------------------------- */
        /* Prepare to exit. Put the converged Ritz values */
        /* and corresponding bounds in RITZ(1:NCONV) and  */
        /* BOUNDS(1:NCONV) respectively. Then sort. Be    */
        /* careful when NCONV > NP                        */
        /* ---------------------------------------------- */

        /* ---------------------------------------- */
        /*  Use h( 3,1 ) as storage to communicate  */
        /*  rnorm to _neupd if needed               */
        /* ---------------------------------------- */
        h[h_dim1 + 3] = rnorm;

        /* -------------------------------------------- */
        /* To be consistent with sngets, we first do a  */
        /* pre-processing sort in order to keep complex */
        /* conjugate pairs together.  This is similar   */
        /* to the pre-processing sort used in sngets    */
        /* except that the sort is done in the opposite */
        /* order.                                       */
        /* -------------------------------------------- */

        if (strcmp(which, "LM") == 0)
        {
            strcpy(wprime, "SR");
        }
        if (strcmp(which, "SM") == 0)
        {
            strcpy(wprime, "LR");
        }
        if (strcmp(which, "LR") == 0)
        {
            strcpy(wprime, "SM");
        }
        if (strcmp(which, "SR") == 0)
        {
            strcpy(wprime, "LM");
        }
        if (strcmp(which, "LI") == 0)
        {
            strcpy(wprime, "SM");
        }
        if (strcmp(which, "SI") == 0)
        {
            strcpy(wprime, "LM");
        }

        ssortc_(wprime, &b_true, &kplusp, &ritzr[1], &ritzi[1], &bounds[1]);

        /* -------------------------------------------- */
        /* Now sort Ritz values so that converged Ritz  */
        /* values appear within the first NEV locations */
        /* of ritzr, ritzi and bounds, and the most     */
        /* desired one appears at the front.            */
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

        ssortc_(wprime, &b_true, &kplusp, &ritzr[1], &ritzi[1], &bounds[1]);

        /* ------------------------------------------------ */
        /* Scale the Ritz estimate of each Ritz value       */
        /* by 1 / max(eps23,magnitude of the Ritz value).   */
        /* ------------------------------------------------ */

        i__1 = numcnv;
        for (j = 1; j <= i__1; ++j)
        {
            /* Computing MAX */
            r__1 = eps23, r__2 = slapy2_(&ritzr[j], &ritzi[j]);
            temp = dmax(r__1, r__2);
            bounds[j] /= temp;
            /* L35: */
        }

        /* -------------------------------------------------- */
        /* Sort the Ritz values according to the scaled Ritz  */
        /* estimates.  This will push all the converged ones  */
        /* towards the front of ritzr, ritzi, bounds          */
        /* (in the case when NCONV < NEV.)                    */
        /* -------------------------------------------------- */

        strcpy(wprime, "LR");
        ssortc_(wprime, &b_true, &numcnv, &bounds[1], &ritzr[1], &ritzi[1]);

        /* -------------------------------------------- */
        /* Scale the Ritz estimate back to its original */
        /* value.                                       */
        /* -------------------------------------------- */

        i__1 = numcnv;
        for (j = 1; j <= i__1; ++j)
        {
            /* Computing MAX */
            r__1 = eps23, r__2 = slapy2_(&ritzr[j], &ritzi[j]);
            temp = dmax(r__1, r__2);
            bounds[j] *= temp;
            /* L40: */
        }

        /* ---------------------------------------------- */
        /* Sort the converged Ritz values again so that   */
        /* the "threshold" value appears at the front of  */
        /* ritzr, ritzi and bound.                        */
        /* ---------------------------------------------- */

        ssortc_(which, &b_true, &nconv, &ritzr[1], &ritzi[1], &bounds[1]);

        if (msglvl > 1)
        {
            svout_(kplusp, &ritzr[1], debug_1.ndigit, "_naup2: Sorted float part of the eigenvalues");
            svout_(kplusp, &ritzi[1], debug_1.ndigit, "_naup2: Sorted imaginary part of the eigenvalues");
            svout_(kplusp, &bounds[1], debug_1.ndigit, "_naup2: Sorted ritz estimates.");
        }

        /* ---------------------------------- */
        /* Max iterations have been exceeded. */
        /* ---------------------------------- */

        if (iter > *mxiter && nconv < numcnv)
        {
            *info = 1;
        }

        /* ------------------- */
        /* No shifts to apply. */
        /* ------------------- */

        if (*np == 0 && nconv < numcnv)
        {
            *info = 2;
        }

        *np = nconv;
        goto L1100;
    }
    else if (nconv < numcnv && *ishift == 1)
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
        /*           %---- Scipy fix ------------------------------------------------ */
        /*           | We must keep nev below this value, as otherwise we can get */
        /*           | np == 0 (note that sngets below can bump nev by 1). If np == 0, */
        /*           | the next call to `snaitr` will write out-of-bounds. */
        /*           | */
        if (*nev > kplusp - 2)
        {
            *nev = kplusp - 2;
        }
        /*           | */
        /*           %---- Scipy fix end -------------------------------------------- */

        *np = kplusp - *nev;

        /* ------------------------------------- */
        /* If the size of NEV was just increased */
        /* resort the eigenvalues.               */
        /* ------------------------------------- */

        if (nevbef < *nev)
        {
            sngets_(ishift, which, nev, np, &ritzr[1], &ritzi[1], &bounds[1], &workl[1], &workl[*np + 1]);
        }
    }

    if (msglvl > 0)
    {
        ivout_(1, &nconv, debug_1.ndigit, "_naup2: no. of \"converged\" Ritz values at this iter.");
        if (msglvl > 1)
        {
            kp[0] = *nev;
            kp[1] = *np;
            ivout_(2, kp, debug_1.ndigit, "_naup2: NEV and NP are");
            svout_(*nev, &ritzr[*np + 1], debug_1.ndigit, "_naup2: \"wanted\" Ritz values -- float part");
            svout_(*nev, &ritzi[*np + 1], debug_1.ndigit, "_naup2: \"wanted\" Ritz values -- imag part");
            svout_(*nev, &bounds[*np + 1], debug_1.ndigit, "_naup2: Ritz estimates of the \"wanted\" values ");
        }
    }

    if (*ishift == 0)
    {

        /* ----------------------------------------------------- */
        /* User specified shifts: reverse communication to       */
        /* compute the shifts. They are returned in the first    */
        /* 2*NP locations of WORKL.                              */
        /* ----------------------------------------------------- */

        ushift = TRUE_;
        *ido = 3;
        goto L9000;
    }

L50:

    /* ---------------------------------- */
    /* Back from reverse communication;   */
    /* User specified shifts are returned */
    /* in WORKL(1:2*NP)                   */
    /* ---------------------------------- */

    ushift = FALSE_;

    if (*ishift == 0)
    {

        /* -------------------------------- */
        /* Move the NP shifts from WORKL to */
        /* RITZR, RITZI to free up WORKL    */
        /* for non-exact shift case.        */
        /* -------------------------------- */

        scopy_(np, &workl[1], &i_one, &ritzr[1], &i_one);
        scopy_(np, &workl[*np + 1], &i_one, &ritzi[1], &i_one);
    }

    if (msglvl > 2)
    {
        ivout_(1, np, debug_1.ndigit, "_naup2: The number of shifts to apply ");
        svout_(*np, &ritzr[1], debug_1.ndigit, "_naup2: Real part of the shifts");
        svout_(*np, &ritzi[1], debug_1.ndigit, "_naup2: Imaginary part of the shifts");
        if (*ishift == 1)
        {
            svout_(*np, &bounds[1], debug_1.ndigit, "_naup2: Ritz estimates of the shifts");
        }
    }

    /* ------------------------------------------------------- */
    /* Apply the NP implicit shifts by QR bulge chasing.       */
    /* Each shift is applied to the whole upper Hessenberg     */
    /* matrix H.                                               */
    /* The first 2*N locations of WORKD are used as workspace. */
    /* ------------------------------------------------------- */

    snapps_(n, nev, np, &ritzr[1], &ritzi[1], &v[v_offset], ldv, &h[h_offset], ldh, &resid[1], &q[q_offset], ldq, &workl[1], &workd[1]);

    /* ------------------------------------------- */
    /* Compute the B-norm of the updated residual. */
    /* Keep B*RESID in WORKD(1:N) to be used in    */
    /* the first step of the next call to snaitr.  */
    /* ------------------------------------------- */

    cnorm = TRUE_;
    arscnd_(&t2);
    if (*bmat == 'G')
    {
        ++timing_1.nbx;
        scopy_(n, &resid[1], &i_one, &workd[*n + 1], &i_one);
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
        scopy_(n, &resid[1], &i_one, &workd[1], &i_one);
    }

L100:

    /* -------------------------------- */
    /* Back from reverse communication; */
    /* WORKD(1:N) := B*RESID            */
    /* -------------------------------- */

    if (*bmat == 'G')
    {
        arscnd_(&t3);
        timing_1.tmvbx += t3 - t2;
    }

    if (*bmat == 'G')
    {
        rnorm = sdot_(n, &resid[1], &i_one, &workd[1], &i_one);
        rnorm = sqrt((dabs(rnorm)));
    }
    else if (*bmat == 'I')
    {
        rnorm = snrm2_(n, &resid[1], &i_one);
    }
    cnorm = FALSE_;

    if (msglvl > 2)
    {
        svout_(1, &rnorm, debug_1.ndigit, "_naup2: B-norm of residual for compressed factorization");
        smout_(*nev, *nev, &h[h_offset], ldh, debug_1.ndigit, "_naup2: Compressed upper Hessenberg matrix H");
    }

    goto L1000;

    /* ------------------------------------------------------------- */
    /*                                                               */
    /*  E N D     O F     M A I N     I T E R A T I O N     L O O P  */
    /*                                                               */
    /* ------------------------------------------------------------- */

L1100:

    *mxiter = iter;
    *nev = numcnv;

L1200:
    *ido = 99;

    /* ---------- */
    /* Error Exit */
    /* ---------- */

    arscnd_(&t1);
    timing_1.tnaup2 = t1 - t0;

L9000:

    /* ------------- */
    /* End of snaup2 */
    /* ------------- */

    return 0;
} /* snaup2_ */
