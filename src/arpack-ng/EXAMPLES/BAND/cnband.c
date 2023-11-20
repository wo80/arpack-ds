/* EXAMPLES\BAND\cnband.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

/* Table of constant values */

static a_fcomplex c_b1 = {1.f, 0.f};
static a_fcomplex c_b2 = {0.f, 0.f};
static a_int c__9 = 9;
static a_int c__1 = 1;
static a_int c__3 = 3;
/**
 * \BeginDoc
 *
 * \Name: cnband
 *
 * \Description:
 *  This subroutine returns the converged approximations to eigenvalues
 *  of A*z = lambda*B*z and (optionally):
 *
 *      (1) The corresponding approximate eigenvectors;
 *
 *      (2) An orthonormal basis for the associated approximate
 *          invariant subspace;
 *
 *      (3) Both.
 *
 *  Matrices A and B are stored in LAPACK-style banded form.
 *
 *  There is negligible additional cost to obtain eigenvectors.  An orthonormal
 *  basis is always computed.  There is an additional storage cost of n*nev
 *  if both are requested (in this case a separate array Z must be supplied).
 *
 *  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
 *  are commonly called Ritz values and Ritz vectors respectively.  They are
 *  referred to as such in the comments that follow.  The computed orthonormal
 *  basis for the invariant subspace corresponding to these Ritz values is
 *  referred to as a Schur basis.
 *
 *  cnband can be called with one of the following modes:
 *
 *  Mode 1:  A*z = lambda*z.
 *           ===> OP = A  and  B = I.
 *
 *  Mode 2:  A*z = lambda*M*z, M symmetric positive definite
 *           ===> OP = inv[M]*A  and  B = M.
 *
 *  Mode 3:  A*z = lambda*M*z, M symmetric semi-definite
 *           ===> OP = inv[A - sigma*M]*M   and  B = M.
 *           ===> shift-and-invert mode.
 *
 *  Choice of different modes can be specified in IPARAM(7) defined below.
 *
 * \Usage
 *   call cnband
 *      ( RVEC, HOWMNY, SELECT, D , Z, LDZ, SIGMA, WORKEV, N, AB,
 *        MB, LDA, FAC, KL, KU, WHICH, BMAT, NEV, TOL, RESID, NCV,
 *        V, LDV, IPARAM, WORKD, WORKL, LWORKL, RWORK, IWORK, INFO )
 *
 * \Arguments
 *  RVEC    LOGICAL  (INPUT)
 *          Specifies whether a basis for the invariant subspace corresponding
 *          to the converged Ritz value approximations for the eigenproblem
 *          A*z = lambda*B*z is computed.
 *
 *             RVEC = .FALSE.     Compute Ritz values only.
 *
 *             RVEC = .TRUE.      Compute Ritz vectors or Schur vectors.
 *                                See Remarks below.
 *
 *  HOWMNY  Character*1  (INPUT)
 *          Specifies the form of the invariant subspace to be computed
 *          corresponding to the converged Ritz values.
 *          = 'A': Compute NEV Ritz vectors;
 *          = 'P': Compute NEV Schur vectors;
 *          = 'S': compute some of the Ritz vectors, specified
 *                 by the logical array SELECT.
 *
 *  SELECT  Logical array of dimension NCV.  (INPUT)
 *          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
 *          computed. To select the real Ritz vector corresponding to a
 *          Ritz value D(j), SELECT(j) must be set to .TRUE..
 *          If HOWMNY = 'A' or 'P', SELECT need not be initialized
 *          but it is used as internal workspace.
 *
 *  D       Complex  array of dimension NEV+1.  (OUTPUT)
 *          On exit, D contains the  Ritz  approximations
 *          to the eigenvalues lambda for A*z = lambda*B*z.
 *
 *  Z       Complex  N by NEV array     (OUTPUT)
 *          On exit, if RVEC = .TRUE. and HOWMNY = 'A', then the columns of
 *          Z represents approximate eigenvectors (Ritz vectors) corresponding
 *          to the NCONV=IPARAM(5) Ritz values for eigensystem
 *          A*z = lambda*B*z.
 *
 *          If RVEC = .FALSE. or HOWMNY = 'P', then Z is NOT REFERENCED.
 *
 *          NOTE: If if RVEC = .TRUE. and a Schur basis is not required,
 *          the array Z may be set equal to first NEV columns of the
 *          array V.
 *
 *  LDZ     Integer.  (INPUT)
 *          The leading dimension of the array Z.  If Ritz vectors are
 *          desired, then  LDZ .ge.  max( 1, N ) is required.
 *          In any case,  LDZ .ge. 1 is required.
 *
 *  SIGMA   Complex   (INPUT)
 *          If IPARAM(7) = 3 then SIGMA represents the shift.
 *          Not referenced if IPARAM(7) = 1 or 2.
 *
 *  WORKEV  Complex  work array of dimension NCV.  (WORKSPACE)
 *
 *  N       Integer.  (INPUT)
 *          Dimension of the eigenproblem.
 *
 *  AB      Complex  array of dimension LDA by N. (INPUT)
 *          The matrix A in band storage, in rows KL+1 to
 *          2*KL+KU+1; rows 1 to KL of the array need not be set.
 *          The j-th column of A is stored in the j-th column of the
 *          array AB as follows:
 *          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
 *
 *  MB      Complex  array of dimension LDA by N. (INPUT)
 *          The matrix M in band storage, in rows KL+1 to
 *          2*KL+KU+1; rows 1 to KL of the array need not be set.
 *          The j-th column of M is stored in the j-th column of the
 *          array MB as follows:
 *          MB(kl+ku+1+i-j,j) = M(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
 *          Not referenced if IPARAM(7)=1.
 *
 *  LDA     Integer. (INPUT)
 *          Leading dimension of AB, MB, FAC.
 *
 *  FAC     Complex  array of LDA by N. (WORKSPACE/OUTPUT)
 *          FAC is used to store the LU factors of MB when mode 2
 *          is invoked.  It is used to store the LU factors of
 *          (A-sigma*M) when mode 3 is invoked.
 *          It is not referenced when IPARAM(7)=1.
 *
 *  KL      Integer. (INPUT)
 *          Max(number of subdiagonals of A, number of subdiagonals of M)
 *
 *  KU      Integer. (OUTPUT)
 *          Max(number of superdiagonals of A, number of superdiagonals of M)
 *
 *  WHICH   Character*2.  (INPUT)
 *          When mode 1,2 are used, WHICH can be set to any one of
 *          the following.
 *
 *            'LM' -> want the NEV eigenvalues of largest magnitude.
 *            'SM' -> want the NEV eigenvalues of smallest magnitude.
 *            'LR' -> want the NEV eigenvalues of largest real part.
 *            'SR' -> want the NEV eigenvalues of smallest real part.
 *            'LI' -> want the NEV eigenvalues of largest imaginary part.
 *            'SI' -> want the NEV eigenvalues of smallest imaginary part.
 *
 *          When mode 3 is used, WHICH should be set to 'LM' only.
 *
 *  BMAT    Character*1.  (INPUT)
 *          BMAT specifies the type of the matrix B that defines the
 *          semi-inner product for the operator OP.
 *          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
 *          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*M*x
 *  NEV     Integer. (INPUT)
 *          Number of eigenvalues of to be computed.
 *
 *  TOL     Real  scalar.  (INPUT)
 *          Stopping criteria: the relative accuracy of the Ritz value
 *          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I))
 *          where ABS(RITZ(I)) is the magnitude when RITZ(I) is complex.
 *          DEFAULT = slamch('EPS')  (machine precision as computed
 *                    by the LAPACK auxiliary subroutine slamch).
 *
 *  RESID   Complex  array of length N.  (INPUT/OUTPUT)
 *          On INPUT:
 *          If INFO .EQ. 0, a random initial residual vector is used.
 *          If INFO .NE. 0, RESID contains the initial residual vector,
 *                          possibly from a previous run.
 *          On OUTPUT:
 *          RESID contains the final residual vector.
 *
 *  NCV     Integer.  (INPUT)
 *          Number of columns of the matrix V. NCV must satisfy the two
 *          inequalities 2 <= NCV-NEV and NCV <= N.
 *          This will indicate how many Arnoldi vectors are generated
 *          at each iteration.  After the startup phase in which NEV
 *          Arnoldi vectors are generated, the algorithm generates
 *          approximately NCV-NEV Arnoldi vectors at each subsequent update
 *          iteration. Most of the cost in generating each Arnoldi vector is
 *          in the matrix-vector operation OP*x.
 *
 *  V       Complex  array N by NCV.  (OUTPUT)
 *          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
 *                       contain approximate Schur vectors that span the
 *                       desired invariant subspace.
 *
 *          NOTE: If the array Z has been set equal to first NEV+1 columns
 *          of the array V and RVEC=.TRUE. and HOWMNY= 'A', then
 *          the first NCONV=IPARAM(5) columns of V will contain Ritz vectors
 *          of the eigensystem A*z = lambda*B*z.
 *
 *  LDV     Integer.  (INPUT)
 *          Leading dimension of V exactly as declared in the calling
 *          program.  LDV must be great than or equal to N.
 *
 *  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
 *          IPARAM(1) = ISHIFT:
 *          The shifts selected at each iteration are used to restart
 *          the Arnoldi iteration in an implicit fashion.
 *          It is set to 1 in this subroutine.  The user do not need
 *          to set this parameter.
 *           ----------------------------------------------------------
 *          ISHIFT = 1: exact shift with respect to the current
 *                      Hessenberg matrix H.  This is equivalent to
 *                      restarting the iteration from the beginning
 *                      after updating the starting vector with a linear
 *                      combination of Ritz vectors associated with the
 *                      "wanted" eigenvalues.
 *          -------------------------------------------------------------
 *
 *          IPARAM(2) = Not referenced.
 *
 *          IPARAM(3) = MXITER
 *          On INPUT:  max number of Arnoldi update iterations allowed.
 *          On OUTPUT: actual number of Arnoldi update iterations taken.
 *
 *          IPARAM(4) = NB: blocksize to be used in the recurrence.
 *          The code currently works only for NB = 1.
 *
 *          IPARAM(5) = NCONV: number of "converged" eigenvalues.
 *
 *          IPARAM(6) = IUPD
 *          Not referenced. Implicit restarting is ALWAYS used.
 *
 *          IPARAM(7) = MODE
 *          On INPUT determines what type of eigenproblem is being solved.
 *          Must be 1,2 or 3; See under \Description of cnband for the
 *          three modes available.
 *
 * WORKD    Complex  work array of length at least 3*n. (WORKSPACE)
 *
 * WORKL    Complex  work array of length LWORKL. (WORKSPACE)
 *
 * LWORKL   Integer.  (INPUT)
 *          LWORKL must be at least 3*NCV**2 + 5*NCV.
 *
 * RWORK    Real  array of length N (WORKSPACE)
 *          Workspace used in cnaupd.
 *
 * IWORK    Integer array of dimension at least N. (WORKSPACE)
 *          Used to mode 2,3. Store the pivot information in the
 *          factorization of M or (A-SIGMA*M).
 *
 * INFO     Integer.  (INPUT/OUTPUT)
 *          Error flag on output.
 *          =  0: Normal exit.
 *          = -1: N must be positive.
 *          = -2: NEV must be positive.
 *          = -3: NCV-NEV >= 2 and less than or equal to N.
 *          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
 *          = -6: BMAT must be one of 'I' or 'G'.
 *          = -7: Length of private work WORKL array is not sufficient.
 *          = -8: Error return from LAPACK eigenvalue calculation.
 *                This should never happened.
 *          = -10: IPARAM(7) must be 1,2,3.
 *          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
 *          = -12: HOWMNY = 'S' not yet implemented
 *          = -13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.
 *          = -14: CNAUPD did not find any eigenvalues to sufficient
 *                 accuracy.
 *
 * \EndDoc
 *
 * ------------------------------------------------------------------------
 *
 * \BeginLib
 *
 * \Routines called
 *     cnaupd  ARPACK reverse communication interface routine.
 *     cneupd  ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     cgbtrf  LAPACK band matrix factorization routine.
 *     cgbtrs  LAPACK band linear system solve routine.
 *     clacpy  LAPACK matrix copy routine.
 *     ccopy   Level 1 BLAS that copies one vector to another.
 *     scnrm2  Level 1 BLAS that computes the norm of a vector.
 *     cgbmv   Level 2 BLAS that computes the band matrix vector product.
 *
 * \References:
 *  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
 *     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
 *     pp 357-385.
 *  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
 *     Restarted Arnoldi Iteration", Ph.D thesis, TR95-13, Rice Univ,
 *     May 1995.
 *
 * \Author
 *     Richard Lehoucq
 *     Danny Sorensen
 *     Chao Yang
 *     Dept. of Computational &
 *     Applied Mathematics
 *     Rice University
 *     Houston, Texas
 *
 * \SCCS Information: @(#)
 * FILE: nband.F   SID: 2.3   DATE OF SID: 10/17/00   RELEASE: 2
 *
 * \EndLib
 */
int cnband_(a_bool *rvec, char *howmny, a_bool *select, a_fcomplex *d, a_fcomplex *z, a_int *ldz, a_fcomplex *sigma, a_fcomplex *workev, a_int *n, a_fcomplex *ab, a_fcomplex *mb, a_int *lda, a_fcomplex *fac, a_int *kl, a_int *ku, char *which, char *bmat, a_int *nev, float *tol, a_fcomplex *resid, a_int *ncv, a_fcomplex *v, a_int *ldv, a_int *iparam, a_fcomplex *workd, a_fcomplex *workl, a_int *lworkl, float *rwork, a_int *iwork, a_int *info, ftnlen howmny_len, ftnlen which_len, ftnlen bmat_len)
{
    /* System generated locals */
    a_int v_dim1, v_offset, z_dim1, z_offset, ab_dim1, ab_offset, mb_dim1, mb_offset, fac_dim1, fac_offset, i__1, i__2, i__3, i__4, i__5;
    a_fcomplex q__1, q__2;

    /* Local variables */
    a_int i, j, ido, imid, mode, ibot, ierr, itop;
    a_int ipntr[14];

    /* Parameter adjustments */
    --select;
    --d;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z -= z_offset;
    --workev;
    fac_dim1 = *lda;
    fac_offset = 1 + fac_dim1;
    fac -= fac_offset;
    mb_dim1 = *lda;
    mb_offset = 1 + mb_dim1;
    mb -= mb_offset;
    ab_dim1 = *lda;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    --resid;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --iparam;
    --workd;
    --workl;
    --rwork;
    --iwork;

    mode = iparam[7];

    /* ---------------------- */
    /* Initialize the reverse */
    /* communication flag.    */
    /* ---------------------- */

    ido = 0;

    /* -------------- */
    /* Exact shift is */
    /* used.          */
    /* -------------- */

    iparam[1] = 1;

    /* --------------------------------- */
    /* Both matrices A and M are stored  */
    /* between rows itop and ibot.  Imid */
    /* is the index of the row that      */
    /* stores the diagonal elements.     */
    /* --------------------------------- */

    itop = *kl + 1;
    imid = *kl + *ku + 1;
    ibot = (*kl << 1) + *ku + 1;

    if (mode == 2)
    {

        /* --------------------------------------------- */
        /* Copy M to fac and Call LAPACK routine cgbtrf  */
        /* to factor M.                                  */
        /* --------------------------------------------- */

        clacpy_("A", &ibot, n, &mb[mb_offset], lda, &fac[fac_offset], lda);
        cgbtrf_(n, n, kl, ku, &fac[fac_offset], lda, &iwork[1], &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf("_band:  error in _gbtrf\n");
            printf(" \n");
            return ierr;
        }
    }
    else if (mode == 3)
    {

        if (*bmat == 'I')
        {

            /* ----------------------- */
            /* Construct (A - sigma*I) */
            /* ----------------------- */

            clacpy_("A", &ibot, n, &ab[ab_offset], lda, &fac[fac_offset], lda);
            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                i__2 = imid + j * fac_dim1;
                i__3 = imid + j * ab_dim1;
                q__1.r = ab[i__3].r - sigma->r, q__1.i = ab[i__3].i - sigma->i;
                fac[i__2].r = q__1.r, fac[i__2].i = q__1.i;
                /* L10: */
            }
        }
        else
        {

            /* ------------------------- */
            /* Construct (A - sigma*M)   */
            /* ------------------------- */

            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                i__2 = ibot;
                for (i = itop; i <= i__2; ++i)
                {
                    i__3 = i + j * fac_dim1;
                    i__4 = i + j * ab_dim1;
                    i__5 = i + j * mb_dim1;
                    q__2.r = sigma->r * mb[i__5].r - sigma->i * mb[i__5].i, q__2.i = sigma->r * mb[i__5].i + sigma->i * mb[i__5].r;
                    q__1.r = ab[i__4].r - q__2.r, q__1.i = ab[i__4].i - q__2.i;
                    fac[i__3].r = q__1.r, fac[i__3].i = q__1.i;
                    /* L20: */
                }
                /* L30: */
            }
        }

        /* ---------------------- */
        /* Factor (A - sigma*M)   */
        /* ---------------------- */

        cgbtrf_(n, n, kl, ku, &fac[fac_offset], lda, &iwork[1], &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf("_band: error in _gbtrf.\n");
            printf(" \n");
            return ierr;
        }
    }

    /* ------------------------------------------ */
    /*  M A I N   L O O P (reverse communication) */
    /* ------------------------------------------ */

L40:

    cnaupd_(&ido, bmat, n, which, nev, tol, &resid[1], ncv, &v[v_offset], ldv, &iparam[1], ipntr, &workd[1], &workl[1], lworkl, &rwork[1], info);

    if (ido == -1)
    {

        if (mode == 1)
        {

            /* -------------------------- */
            /* Perform  y <--- OP*x = A*x */
            /* -------------------------- */

            cgbmv_("Notranspose", n, n, kl, ku, &c_b1, &ab[itop + ab_dim1], lda, &workd[ipntr[0]], &c__1, &c_b2, &workd[ipntr[1]], &c__1);
        }
        else if (mode == 2)
        {

            /* --------------------------------- */
            /* Perform  y <--- OP*x = inv[M]*A*x */
            /* --------------------------------- */

            cgbmv_("Notranspose", n, n, kl, ku, &c_b1, &ab[itop + ab_dim1], lda, &workd[ipntr[0]], &c__1, &c_b2, &workd[ipntr[1]], &c__1);

            cgbtrs_("Notranspose", n, kl, ku, &c__1, &fac[fac_offset], lda, &iwork[1], &workd[ipntr[1]], n, &ierr);
            if (ierr != 0)
            {
                printf(" \n");
                printf("_band: error in sbgtrs.\n");
                printf(" \n");
                return ierr;
            }
        }
        else if (mode == 3)
        {

            /* --------------------------------------- */
            /* Perform y <-- OP*x                      */
            /*           |           = inv[A-SIGMA*M]*M* x */
            /*           | to force the starting vector into the   | */
            /*           | range of OP.                            | */
            /* --------------------------------------- */

            cgbmv_("Notranspose", n, n, kl, ku, &c_b1, &mb[itop + mb_dim1], lda, &workd[ipntr[0]], &c__1, &c_b2, &workd[ipntr[1]], &c__1);

            cgbtrs_("Notranspose", n, kl, ku, &c__1, &fac[fac_offset], lda, &iwork[1], &workd[ipntr[1]], n, &ierr);
            if (ierr != 0)
            {
                printf(" \n");
                printf("_band: error in _gbtrs.\n");
                printf(" \n");
                return ierr;
            }
        }
    }
    else if (ido == 1)
    {

        if (mode == 1)
        {

            /* -------------------------- */
            /* Perform  y <--- OP*x = A*x */
            /* -------------------------- */

            cgbmv_("Notranspose", n, n, kl, ku, &c_b1, &ab[itop + ab_dim1], lda, &workd[ipntr[0]], &c__1, &c_b2, &workd[ipntr[1]], &c__1);
        }
        else if (mode == 2)
        {

            /* --------------------------------- */
            /* Perform  y <--- OP*x = inv[M]*A*x */
            /* --------------------------------- */

            cgbmv_("Notranspose", n, n, kl, ku, &c_b1, &ab[itop + ab_dim1], lda, &workd[ipntr[0]], &c__1, &c_b2, &workd[ipntr[1]], &c__1);

            cgbtrs_("Notranspose", n, kl, ku, &c__1, &fac[fac_offset], lda, &iwork[1], &workd[ipntr[1]], ldv, &ierr);
            if (ierr != 0)
            {
                printf(" \n");
                printf("_band: error in sbgtrs.\n");
                printf(" \n");
                return ierr;
            }
        }
        else if (mode == 3)
        {

            if (*bmat == 'I')
            {

                /* -------------------------------- */
                /* Perform  y <-- inv(A-sigma*I)*x. */
                /* -------------------------------- */

                ccopy_(n, &workd[ipntr[0]], &c__1, &workd[ipntr[1]], &c__1);
                cgbtrs_("Notranspose", n, kl, ku, &c__1, &fac[fac_offset], lda, &iwork[1], &workd[ipntr[1]], n, &ierr);
                if (ierr != 0)
                {
                    printf(" \n");
                    printf("_band: error in _gbtrs.\n");
                    printf(" \n");
                    return ierr;
                }
            }
            else
            {

                /* ------------------------------------ */
                /* Perform  y <-- inv(A-sigma*M)*(M*x). */
                /* (M*x) has been computed and stored   */
                /* in workd(ipntr(3)).                  */
                /* ------------------------------------ */

                ccopy_(n, &workd[ipntr[2]], &c__1, &workd[ipntr[1]], &c__1);
                cgbtrs_("Notranspose", n, kl, ku, &c__1, &fac[fac_offset], lda, &iwork[1], &workd[ipntr[1]], n, &ierr);
                if (ierr != 0)
                {
                    printf(" \n");
                    printf("_band: error in _gbtrs.\n");
                    printf(" \n");
                    return ierr;
                }
            }
        }
    }
    else if (ido == 2)
    {

        /* ------------------ */
        /* Perform y <-- M*x  */
        /* ------------------ */

        cgbmv_("Notranspose", n, n, kl, ku, &c_b1, &mb[itop + mb_dim1], lda, &workd[ipntr[0]], &c__1, &c_b2, &workd[ipntr[1]], &c__1);
    }
    else
    {

        /* ----------------------------------------- */
        /*   Either we have convergence, or there is */
        /*   error.                                  */
        /* ----------------------------------------- */

        if (*info != 0)
        {

            /* ------------------------ */
            /* Error message, check the */
            /* documentation in dnaupd  */
            /* ------------------------ */

            printf(" \n");
            printf(" Error with _naupd info = %d", (*info));
            printf(" Check the documentation of _naupd \n");
            printf(" \n");
        }
        else
        {

            cneupd_(rvec, howmny, &select[1], &d[1], &z[z_offset], ldz, sigma, &workev[1], bmat, n, which, nev, tol, &resid[1], ncv, &v[v_offset], ldv, &iparam[1], ipntr, &workd[1], &workl[1], lworkl, &rwork[1], info);

            if (*info != 0)
            {

                /* ---------------------------------- */
                /* Check the documentation of cneupd. */
                /* ---------------------------------- */

                printf(" \n");
                printf(" Error with _neupd = %d", (*info));
                printf(" Check the documentation of _neupd \n");
                printf(" \n");
            }
        }

        return ierr;
    }

    /* -------------------------------------- */
    /* L O O P  B A C K to call cnaupd again. */
    /* -------------------------------------- */

    goto L40;

    return 0;
} /* cnband_ */
