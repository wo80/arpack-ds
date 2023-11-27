/* EXAMPLES\BAND\ssband.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

/* Table of constant values */

static a_int i_one = 1;
static float one = 1.f;
static float zero = 0.f;

/**
 * \BeginDoc
 *
 * \Name: ssband
 *
 * \Description:
 *
 *  This subroutine returns the converged approximations to eigenvalues
 *  of A*z = lambda*B*z and (optionally):
 *
 *      (1) The corresponding approximate eigenvectors;
 *
 *      (2) An orthonormal (Lanczos) basis for the associated approximate
 *          invariant subspace;
 *
 *      (3) Both.
 *
 *  Matrices A and B are stored in LAPACK-style band form.
 *
 *  There is negligible additional cost to obtain eigenvectors.  An orthonormal
 *  (Lanczos) basis is always computed.  There is an additional storage cost
 *  of n*nev if both are requested (in this case a separate array Z must be
 *  supplied).
 *
 *  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
 *  are called Ritz values and Ritz vectors respectively.  They are referred
 *  to as such in the comments that follow.  The computed orthonormal basis
 *  for the invariant subspace corresponding to these Ritz values is referred
 *  to as a Lanczos basis.
 *
 *  ssband can be called with one of the following modes:
 *
 *  Mode 1:  A*x = lambda*x, A symmetric
 *           ===> OP = A  and  B = I.
 *
 *  Mode 2:  A*x = lambda*M*x, A symmetric, M symmetric positive definite
 *           ===> OP = inv[M]*A  and  B = M.
 *           ===> (If M can be factored see remark 3 in SSAUPD)
 *
 *  Mode 3:  K*x = lambda*M*x, K symmetric, M symmetric positive semi-definite
 *           ===> OP = (inv[K - sigma*M])*M  and  B = M.
 *           ===> Shift-and-Invert mode
 *
 *  Mode 4:  K*x = lambda*KG*x, K symmetric positive semi-definite,
 *           KG symmetric indefinite
 *           ===> OP = (inv[K - sigma*KG])*K  and  B = K.
 *           ===> Buckling mode
 *
 *  Mode 5:  A*x = lambda*M*x, A symmetric, M symmetric positive semi-definite
 *           ===> OP = inv[A - sigma*M]*[A + sigma*M]  and  B = M.
 *           ===> Cayley transformed mode
 *
 *  The choice of mode must be specified in IPARAM(7) defined below.
 *
 * \Usage
 *   call ssband
 *      ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, N, AB, MB, LDA,
 *        RFAC, KL, KU, WHICH, BMAT, NEV, TOL, RESID, NCV, V,
 *        LDV, IPARAM, WORKD, WORKL, LWORKL, IWORK, INFO )
 *
 * \Arguments
 *
 *  RVEC    Logical (INPUT)
 *          Specifies whether Ritz vectors corresponding to the Ritz value
 *          approximations to the eigenproblem A*z = lambda*B*z are computed.
 *
 *             RVEC = .FALSE.     Compute Ritz values only.
 *
 *             RVEC = .TRUE.      Compute the associated Ritz vectors.
 *
 *  HOWMNY  Character*1  (INPUT)
 *          Specifies how many Ritz vectors are wanted and the form of Z
 *          the matrix of Ritz vectors. See remark 1 below.
 *          = 'A': compute all Ritz vectors;
 *          = 'S': compute some of the Ritz vectors, specified
 *                 by the logical array SELECT.
 *
 *  SELECT  Logical array of dimension NCV.  (INPUT)
 *          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
 *          computed. To select the Ritz vector corresponding to a
 *          Ritz value D(j), SELECT(j) must be set to .TRUE..
 *          If HOWMNY = 'A' , SELECT is not referenced.
 *
 *  D       Real array of dimension NEV.  (OUTPUT)
 *          On exit, D contains the Ritz value approximations to the
 *          eigenvalues of A*z = lambda*B*z. The values are returned
 *          in ascending order. If IPARAM(7) = 3,4,5 then D represents
 *          the Ritz values of OP computed by ssaupd transformed to
 *          those of the original eigensystem A*z = lambda*B*z. If
 *          IPARAM(7) = 1,2 then the Ritz values of OP are the same
 *          as the those of A*z = lambda*B*z.
 *
 *  Z       Real N by NEV array if HOWMNY = 'A'.  (OUTPUT)
 *          On exit, Z contains the B-orthonormal Ritz vectors of the
 *          eigensystem A*z = lambda*B*z corresponding to the Ritz
 *          value approximations.
 *
 *          If  RVEC = .FALSE. then Z is not referenced.
 *          NOTE: The array Z may be set equal to first NEV columns of the
 *          Lanczos basis array V computed by SSAUPD.
 *
 *  LDZ     Integer.  (INPUT)
 *          The leading dimension of the array Z.  If Ritz vectors are
 *          desired, then  LDZ .ge.  max( 1, N ).  In any case,  LDZ .ge. 1.
 *
 *  SIGMA   Real  (INPUT)
 *          If IPARAM(7) = 3,4,5 represents the shift. Not referenced if
 *          IPARAM(7) = 1 or 2.
 *
 *  N       Integer.  (INPUT)
 *          Dimension of the eigenproblem.
 *
 *  AB      Real array of dimension LDA by N. (INPUT)
 *          The matrix A in band storage, in rows KL+1 to
 *          2*KL+KU+1; rows 1 to KL of the array need not be set.
 *          The j-th column of A is stored in the j-th column of the
 *          array AB as follows:
 *          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
 *
 *  MB      Real array of dimension LDA by N. (INPUT)
 *          The matrix M in band storage, in rows KL+1 to
 *          2*KL+KU+1; rows 1 to KL of the array need not be set.
 *          The j-th column of M is stored in the j-th column of the
 *          array AB as follows:
 *          MB(kl+ku+1+i-j,j) = M(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
 *          Not referenced if IPARAM(7) = 1
 *
 *  LDA     Integer. (INPUT)
 *          Leading dimension of AB, MB, RFAC.
 *
 *  RFAC    Real array of LDA by N. (WORKSPACE/OUTPUT)
 *          RFAC is used to store the LU factors of MB when IPARAM(7) = 2
 *          is invoked.  It is used to store the LU factors of
 *          (A-sigma*M) when IPARAM(7) = 3,4,5 is invoked.
 *          It is not referenced when IPARAM(7) = 1.
 *
 *  KL      Integer. (INPUT)
 *          Max(number of subdiagonals of A, number of subdiagonals of M)
 *
 *  KU      Integer. (OUTPUT)
 *          Max(number of superdiagonals of A, number of superdiagonals of M)
 *
 *  WHICH   Character*2.  (INPUT)
 *          When IPARAM(7)= 1 or 2,  WHICH can be set to any one of
 *          the following.
 *
 *            'LM' -> want the NEV eigenvalues of largest magnitude.
 *            'SM' -> want the NEV eigenvalues of smallest magnitude.
 *            'LA' -> want the NEV eigenvalues of largest REAL part.
 *            'SA' -> want the NEV eigenvalues of smallest REAL part.
 *            'BE' -> Compute NEV eigenvalues, half from each end of the
 *                    spectrum.  When NEV is odd, compute one more from
 *                    the high end than from the low end.
 *
 *          When IPARAM(7) = 3, 4, or 5,  WHICH should be set to 'LM' only.
 *
 *  BMAT    Character*1.  (INPUT)
 *          BMAT specifies the type of the matrix B that defines the
 *          semi-inner product for the operator OP.
 *          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
 *          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*M*x
 *  NEV     Integer. (INPUT)
 *          Number of eigenvalues of OP to be computed.
 *
 *  TOL     Real scalar.  (INPUT)
 *          Stopping criterion: the relative accuracy of the Ritz value
 *          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I)).
 *          If TOL .LE. 0. is passed a default is set:
 *          DEFAULT = SLAMCH('EPS')  (machine precision as computed
 *                    by the LAPACK auxiliary subroutine SLAMCH).
 *
 *  RESID   Real array of length N.  (INPUT/OUTPUT)
 *          On INPUT:
 *          If INFO .EQ. 0, a random initial residual vector is used.
 *          If INFO .NE. 0, RESID contains the initial residual vector,
 *                          possibly from a previous run.
 *          On OUTPUT:
 *          RESID contains the final residual vector.
 *
 *  NCV     Integer.  (INPUT)
 *          Number of columns of the matrix V (less than or equal to N).
 *          Represents the dimension of the Lanczos basis constructed
 *          by ssaupd for OP.
 *
 *  V       Real array N by NCV.  (OUTPUT)
 *          Upon INPUT: the NCV columns of V contain the Lanczos basis
 *                      vectors as constructed by ssaupd for OP.
 *          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
 *                       represent the Ritz vectors that span the desired
 *                       invariant subspace.
 *          NOTE: The array Z may be set equal to first NEV columns of the
 *          Lanczos basis vector array V computed by ssaupd. In this case
 *          if RVEC=.TRUE., the first NCONV=IPARAM(5) columns of V contain
 *          the desired Ritz vectors.
 *
 *  LDV     Integer.  (INPUT)
 *          Leading dimension of V exactly as declared in the calling
 *          program.
 *
 *  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
 *          IPARAM(1) = ISHIFT:
 *          The shifts selected at each iteration are used to restart
 *          the Arnoldi iteration in an implicit fashion.
 *          It is set to 1 in this subroutine.  The user do not need
 *          to set this parameter.
 *          ------------------------------------------------------------
 *          ISHIFT = 1: exact shifts with respect to the reduced
 *                      tridiagonal matrix T.  This is equivalent to
 *                      restarting the iteration with a starting vector
 *                      that is a linear combination of Ritz vectors
 *                      associated with the "wanted" Ritz values.
 *          -------------------------------------------------------------
 *
 *          IPARAM(2) = No longer referenced.
 *
 *          IPARAM(3) = MXITER
 *          On INPUT:  max number of Arnoldi update iterations allowed.
 *          On OUTPUT: actual number of Arnoldi update iterations taken.
 *
 *          IPARAM(4) = NB: blocksize to be used in the recurrence.
 *          The code currently works only for NB = 1.
 *
 *          IPARAM(5) = NCONV: number of "converged" eigenvalues.
 *          This represents the number of Ritz values that satisfy
 *          the convergence criterion.
 *
 *          IPARAM(6) = IUPD
 *          No longer referenced. Implicit restarting is ALWAYS used.
 *
 *          IPARAM(7) = MODE
 *          On INPUT determines what type of eigenproblem is being solved.
 *          Must be 1,2,3,4,5; See under \Description of ssband for the
 *          five modes available.
 *
 *          IPARAM(8) = NP
 *          Not referenced.
 *
 *          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
 *          OUTPUT: NUMOP  = total number of OP*x operations,
 *                  NUMOPB = total number of B*x operations if BMAT='G',
 *                  NUMREO = total number of steps of re-orthogonalization.
 *
 * WORKD    Real work array of length at least 3*n. (WORKSPACE)
 *
 * WORKL    Real work array of length LWORKL.  (WORKSPACE)
 *
 * LWORKL   Integer.  (INPUT)
 *          LWORKL must be at least NCV**2 + 8*NCV.
 *
 * IWORK    Integer array of dimension at least N. (WORKSPACE)
 *          Used when IPARAM(7)=2,3,4,5 to store the pivot information in the
 *          factorization of M or (A-SIGMA*M).
 *
 * INFO     Integer.  (INPUT/OUTPUT)
 *          Error flag on output.
 *          =  0: Normal exit.
 *          =  1: Maximum number of iterations taken.
 *                All possible eigenvalues of OP has been found. IPARAM(5)
 *                returns the number of wanted converged Ritz values.
 *          =  3: No shifts could be applied during a cycle of the
 *                Implicitly restarted Arnoldi iteration. One possibility
 *                is to increase the size of NCV relative to NEV.
 *                See remark 4 in SSAUPD.
 *
 *          = -1: N must be positive.
 *          = -2: NEV must be positive.
 *          = -3: NCV-NEV >= 2 and less than or equal to N.
 *          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
 *          = -6: BMAT must be one of 'I' or 'G'.
 *          = -7: Length of private work WORKL array is not sufficient.
 *          = -8: Error return from trid. eigenvalue calculation;
 *                Informational error from LAPACK routine ssteqr.
 *          = -9: Starting vector is zero.
 *          = -10: IPARAM(7) must be 1,2,3,4,5.
 *          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
 *          = -12: NEV and WHICH = 'BE' are incompatible.
 *          = -13: HOWMNY must be one of 'A' or 'P'
 *          = -14: SSAUPD did not find any eigenvalues to sufficient
 *                 accuracy.
 *          = -9999: Could not build an Arnoldi factorization.
 *                   IPARAM(5) returns the size of the current
 *                   Arnoldi factorization.
 *
 * \EndDoc
 *
 * ------------------------------------------------------------------------
 *
 * \BeginLib
 *
 * \References:
 *  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
 *     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
 *     pp 357-385.
 *
 *  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
 *     Restarted Arnoldi Iteration", Ph.D thesis, TR95-13, Rice Univ,
 *     May 1995.
 *
 * Routines called:
 *     ssaupd  ARPACK reverse communication interface routine.
 *     sseupd  ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     sgbtrf  LAPACK band matrix factorization routine.
 *     sgbtrs  LAPACK band linear system solve routine.
 *     slacpy  LAPACK matrix copy routine.
 *     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     scopy   Level 1 BLAS that copies one vector to another.
 *     sdot    Level 1 BLAS that computes the dot product of two vectors.
 *     snrm2   Level 1 BLAS that computes the norm of a vector.
 *     sgbmv   Level 2 BLAS that computes the band matrix vector product.
 *
 * \Remarks
 *  1. The converged Ritz values are always returned in increasing
 *     (algebraic) order.
 *
 *  2. Currently only HOWMNY = 'A' is implemented. It is included at this
 *     stage for the user who wants to incorporate it.
 *
 * \Author
 *     Danny Sorensen
 *     Richard Lehoucq
 *     Chao Yang
 *     Dept. of Computational &
 *     Applied Mathematics
 *     Rice University
 *     Houston, Texas
 *
 * \SCCS Information: @(#)
 * FILE: sband.F   SID: 2.3   DATE OF SID: 10/17/00   RELEASE: 2
 *
 * \EndLib
 */
int ssband_(a_bool *rvec, char *howmny, a_bool *select, float *d, float *z, a_int *ldz, float *sigma, a_int *n, float *ab, float *mb, a_int *lda, float *rfac, a_int *kl, a_int *ku, char *which, char *bmat, a_int *nev, float *tol, float *resid, a_int *ncv, float *v, a_int *ldv, a_int *iparam, float *workd, float *workl, a_int *lworkl, a_int *iwork, a_int *info)
{
    /* System generated locals */
    a_int v_dim1, v_offset, z_dim1, z_offset, ab_dim1, ab_offset, mb_dim1, mb_offset, rfac_dim1, rfac_offset, i__1, i__2;

    /* Local variables */
    a_int i, j, ido, imid, ibot, itop, type, ierr = 0;
    a_int ipntr[14];

    /* -------------------------------------------------------------- */
    /* Set type of the problem to be solved. Check consistency        */
    /* between BMAT and IPARAM(7).                                    */
    /* type = 1 --> Solving standard problem in regular mode.         */
    /* type = 2 --> Solving standard problem in shift-invert mode.    */
    /* type = 3 --> Solving generalized problem in regular mode.      */
    /* type = 4 --> Solving generalized problem in shift-invert mode. */
    /* type = 5 --> Solving generalized problem in Buckling mode.     */
    /* type = 6 --> Solving generalized problem in Cayley mode.       */
    /* -------------------------------------------------------------- */

    /* Parameter adjustments */
    --select;
    --d;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z -= z_offset;
    rfac_dim1 = *lda;
    rfac_offset = 1 + rfac_dim1;
    rfac -= rfac_offset;
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
    --iwork;

    if (iparam[7] == 1)
    {
        type = 1;
    }
    else if (iparam[7] == 3 && *bmat == 'I')
    {
        type = 2;
    }
    else if (iparam[7] == 2)
    {
        type = 3;
    }
    else if (iparam[7] == 3 && *bmat == 'G')
    {
        type = 4;
    }
    else if (iparam[7] == 4)
    {
        type = 5;
    }
    else if (iparam[7] == 5)
    {
        type = 6;
    }
    else
    {
        printf(" \n");
        printf("BMAT is inconsistent with IPARAM(7).\n");
        printf(" \n");
        return ierr;
    }

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

    if (type == 2 || type == 6 && *bmat == 'I')
    {

        /* -------------------------------- */
        /* Solving a standard eigenvalue    */
        /* problem in shift-invert or       */
        /* Cayley mode. Factor (A-sigma*I). */
        /* -------------------------------- */

        slacpy_("A", &ibot, n, &ab[ab_offset], lda, &rfac[rfac_offset], lda);
        i__1 = *n;
        for (j = 1; j <= i__1; ++j)
        {
            rfac[imid + j * rfac_dim1] = ab[imid + j * ab_dim1] - *sigma;
        }
        sgbtrf_(n, n, kl, ku, &rfac[rfac_offset], lda, &iwork[1], &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf(" _SBAND: Error with _gbtrf. \n");
            printf(" \n");
            return ierr;
        }
    }
    else if (type == 3)
    {

        /* -------------------------------------------- */
        /* Solving generalized eigenvalue problem in    */
        /* regular mode. Copy M to rfac and Call LAPACK */
        /* routine sgbtrf to factor M.                  */
        /* -------------------------------------------- */

        slacpy_("A", &ibot, n, &mb[mb_offset], lda, &rfac[rfac_offset], lda);
        sgbtrf_(n, n, kl, ku, &rfac[rfac_offset], lda, &iwork[1], &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf("_SBAND:  Error with _gbtrf.\n");
            printf(" \n");
            return ierr;
        }
    }
    else if (type == 4 || type == 5 || type == 6 && *bmat == 'G')
    {

        /* ----------------------------------------- */
        /* Solving generalized eigenvalue problem in */
        /* shift-invert, Buckling, or Cayley mode.   */
        /* ----------------------------------------- */

        /* ----------------------------------- */
        /* Construct and factor (A - sigma*M). */
        /* ----------------------------------- */

        i__1 = *n;
        for (j = 1; j <= i__1; ++j)
        {
            i__2 = ibot;
            for (i = itop; i <= i__2; ++i)
            {
                rfac[i + j * rfac_dim1] = ab[i + j * ab_dim1] - *sigma * mb[i + j * mb_dim1];
            }
        }

        sgbtrf_(n, n, kl, ku, &rfac[rfac_offset], lda, &iwork[1], &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf("_SBAND: Error with _gbtrf.\n");
            printf(" \n");
            return ierr;
        }
    }

    /* ------------------------------------------ */
    /*  M A I N   L O O P (reverse communication) */
    /* ------------------------------------------ */

L90:

    ssaupd_(&ido, bmat, n, which, nev, tol, &resid[1], ncv, &v[v_offset], ldv, &iparam[1], ipntr, &workd[1], &workl[1], lworkl, info);

    if (ido == -1)
    {

        if (type == 1)
        {

            /* -------------------------- */
            /* Perform  y <--- OP*x = A*x */
            /* -------------------------- */

            sgbmv_("N", n, n, kl, ku, &one, &ab[itop + ab_dim1], lda, &workd[ipntr[0]], &i_one, &zero, &workd[ipntr[1]], &i_one);
        }
        else if (type == 2)
        {

            /* -------------------------------- */
            /*             Perform              */
            /* y <--- OP*x = inv[A-sigma*I]*x   */
            /* to force the starting vector     */
            /* into the range of OP.            */
            /* -------------------------------- */

            scopy_(n, &workd[ipntr[0]], &i_one, &workd[ipntr[1]], &i_one);
            sgbtrs_("N", n, kl, ku, &i_one, &rfac[rfac_offset], lda, &iwork[1], &workd[ipntr[1]], n, &ierr);
            if (ierr != 0)
            {
                printf(" \n");
                printf(" _SBAND: Error with _bgtrs. \n");
                printf(" \n");
                return ierr;
            }
        }
        else if (type == 3)
        {

            /* --------------------------------- */
            /* Perform  y <--- OP*x = inv[M]*A*x */
            /* to force the starting vector into */
            /* the range of OP.                  */
            /* --------------------------------- */

            sgbmv_("N", n, n, kl, ku, &one, &ab[itop + ab_dim1], lda, &workd[ipntr[0]], &i_one, &zero, &workd[ipntr[1]], &i_one);
            scopy_(n, &workd[ipntr[1]], &i_one, &workd[ipntr[0]], &i_one);
            sgbtrs_("N", n, kl, ku, &i_one, &rfac[rfac_offset], lda, &iwork[1], &workd[ipntr[1]], n, &ierr);
            if (ierr != 0)
            {
                printf(" \n");
                printf("_SBAND: Error with sbgtrs.\n");
                printf(" \n");
                return ierr;
            }
        }
        else if (type == 4)
        {

            /* --------------------------------------- */
            /* Perform y <-- OP*x                      */
            /*           = inv[A-SIGMA*M]*M            */
            /* to force the starting vector into the   */
            /* range of OP.                            */
            /* --------------------------------------- */

            sgbmv_("N", n, n, kl, ku, &one, &mb[itop + mb_dim1], lda, &workd[ipntr[0]], &i_one, &zero, &workd[ipntr[1]], &i_one);
            sgbtrs_("N", n, kl, ku, &i_one, &rfac[rfac_offset], lda, &iwork[1], &workd[ipntr[1]], n, &ierr);
            if (ierr != 0)
            {
                printf(" \n");
                printf("_SBAND: Error with _gbtrs.\n");
                printf(" \n");
                return ierr;
            }
        }
        else if (type == 5)
        {

            /* ------------------------------------- */
            /* Perform y <-- OP*x                    */
            /*    = inv[A-SIGMA*M]*A                 */
            /* to force the starting vector into the */
            /* range of OP.                          */
            /* ------------------------------------- */

            sgbmv_("N", n, n, kl, ku, &one, &ab[itop + ab_dim1], lda, &workd[ipntr[0]], &i_one, &zero, &workd[ipntr[1]], &i_one);
            sgbtrs_("N", n, kl, ku, &i_one, &rfac[rfac_offset], lda, &iwork[1], &workd[ipntr[1]], n, &ierr);

            if (ierr != 0)
            {
                printf(" \n");
                printf(" _SBAND: Error with _gbtrs. \n");
                printf(" \n");
                return ierr;
            }
        }
        else if (type == 6)
        {

            /* ------------------------------------- */
            /* Perform y <-- OP*x                    */
            /* = (inv[A-SIGMA*M])*(A+SIGMA*M)*x      */
            /* to force the starting vector into the */
            /* range of OP.                          */
            /* ------------------------------------- */

            if (*bmat == 'G')
            {
                sgbmv_("N", n, n, kl, ku, &one, &ab[itop + ab_dim1], lda, &workd[ipntr[0]], &i_one, &zero, &workd[ipntr[1]], &i_one);
                sgbmv_("N", n, n, kl, ku, sigma, &mb[itop + mb_dim1], lda, &workd[ipntr[0]], &i_one, &one, &workd[ipntr[1]], &i_one);
            }
            else
            {
                scopy_(n, &workd[ipntr[0]], &i_one, &workd[ipntr[1]], &i_one);
                sgbmv_("N", n, n, kl, ku, &one, &ab[itop + ab_dim1], lda, &workd[ipntr[0]], &i_one, sigma, &workd[ipntr[1]], &i_one);
            }

            sgbtrs_("N", n, kl, ku, &i_one, &rfac[rfac_offset], lda, &iwork[1], &workd[ipntr[1]], n, &ierr);

            if (ierr != 0)
            {
                printf(" \n");
                printf("_SBAND: Error with _gbtrs.\n");
                printf(" \n");
                return ierr;
            }
        }
    }
    else if (ido == 1)
    {

        if (type == 1)
        {

            /* -------------------------- */
            /* Perform  y <--- OP*x = A*x */
            /* -------------------------- */

            sgbmv_("N", n, n, kl, ku, &one, &ab[itop + ab_dim1], lda, &workd[ipntr[0]], &i_one, &zero, &workd[ipntr[1]], &i_one);
        }
        else if (type == 2)
        {

            /* -------------------------------- */
            /*             Perform              */
            /* y <--- OP*x = inv[A-sigma*I]*x.  */
            /* -------------------------------- */

            scopy_(n, &workd[ipntr[0]], &i_one, &workd[ipntr[1]], &i_one);
            sgbtrs_("N", n, kl, ku, &i_one, &rfac[rfac_offset], lda, &iwork[1], &workd[ipntr[1]], n, &ierr);
            if (ierr != 0)
            {
                printf(" \n");
                printf("_SBAND: Error with _gbtrs.\n");
                printf(" \n");
                return ierr;
            }
        }
        else if (type == 3)
        {

            /* --------------------------------- */
            /* Perform  y <--- OP*x = inv[M]*A*x */
            /* --------------------------------- */

            sgbmv_("N", n, n, kl, ku, &one, &ab[itop + ab_dim1], lda, &workd[ipntr[0]], &i_one, &zero, &workd[ipntr[1]], &i_one);
            scopy_(n, &workd[ipntr[1]], &i_one, &workd[ipntr[0]], &i_one);
            sgbtrs_("N", n, kl, ku, &i_one, &rfac[rfac_offset], lda, &iwork[1], &workd[ipntr[1]], n, &ierr);
            if (ierr != 0)
            {
                printf(" \n");
                printf("_SBAND: error with _bgtrs.\n");
                printf(" \n");
                return ierr;
            }
        }
        else if (type == 4)
        {

            /* ----------------------------------- */
            /* Perform y <-- inv(A-sigma*M)*(M*x). */
            /* (M*x) has been computed and stored  */
            /* in workd(ipntr(3)).                 */
            /* ----------------------------------- */

            scopy_(n, &workd[ipntr[2]], &i_one, &workd[ipntr[1]], &i_one);
            sgbtrs_("N", n, kl, ku, &i_one, &rfac[rfac_offset], lda, &iwork[1], &workd[ipntr[1]], n, &ierr);
            if (ierr != 0)
            {
                printf(" \n");
                printf("_SBAND: Error with _gbtrs.\n");
                printf(" \n");
                return ierr;
            }
        }
        else if (type == 5)
        {

            /* ----------------------------- */
            /* Perform y <-- OP*x            */
            /*    = inv[A-SIGMA*M]*A*x       */
            /* B*x = A*x has been computed   */
            /* and saved in workd(ipntr(3)). */
            /* ----------------------------- */

            scopy_(n, &workd[ipntr[2]], &i_one, &workd[ipntr[1]], &i_one);
            sgbtrs_("N", n, kl, ku, &i_one, &rfac[rfac_offset], lda, &iwork[1], &workd[ipntr[1]], n, &ierr);
            if (ierr != 0)
            {
                printf(" \n");
                printf(" _SBAND: Error with _gbtrs. \n");
                printf(" \n");
                return ierr;
            }
        }
        else if (type == 6)
        {

            /* ------------------------------- */
            /* Perform y <-- OP*x              */
            /* = inv[A-SIGMA*M]*(A+SIGMA*M)*x. */
            /* (M*x) has been saved in         */
            /* workd(ipntr(3)).                */
            /* ------------------------------- */

            if (*bmat == 'G')
            {
                sgbmv_("N", n, n, kl, ku, &one, &ab[itop + ab_dim1], lda, &workd[ipntr[0]], &i_one, &zero, &workd[ipntr[1]], &i_one);
                saxpy_(n, sigma, &workd[ipntr[2]], &i_one, &workd[ipntr[1]], &i_one);
            }
            else
            {
                scopy_(n, &workd[ipntr[0]], &i_one, &workd[ipntr[1]], &i_one);
                sgbmv_("N", n, n, kl, ku, &one, &ab[itop + ab_dim1], lda, &workd[ipntr[0]], &i_one, sigma, &workd[ipntr[1]], &i_one);
            }
            sgbtrs_("N", n, kl, ku, &i_one, &rfac[rfac_offset], lda, &iwork[1], &workd[ipntr[1]], n, &ierr);
        }
    }
    else if (ido == 2)
    {

        /* -------------------------------- */
        /*        Perform y <-- B*x         */
        /* Note when Buckling mode is used, */
        /* B = A, otherwise B=M.            */
        /* -------------------------------- */

        if (type == 5)
        {

            /* ------------------- */
            /* Buckling Mode, B=A. */
            /* ------------------- */

            sgbmv_("N", n, n, kl, ku, &one, &ab[itop + ab_dim1], lda, &workd[ipntr[0]], &i_one, &zero, &workd[ipntr[1]], &i_one);
        }
        else
        {
            sgbmv_("N", n, n, kl, ku, &one, &mb[itop + mb_dim1], lda, &workd[ipntr[0]], &i_one, &zero, &workd[ipntr[1]], &i_one);
        }
    }
    else
    {

        /* --------------------------------------- */
        /* Either we have convergence, or there is */
        /* error.                                  */
        /* --------------------------------------- */

        if (*info < 0)
        {

            /* ------------------------ */
            /* Error message, check the */
            /* documentation in SSAUPD  */
            /* ------------------------ */

            printf(" \n");
            printf(" Error with _saupd info = %d", (*info));
            printf(" Check the documentation of _saupd \n");
            printf(" \n");
            return ierr;
        }
        else
        {

            if (*info == 1)
            {
                printf(" \n");
                printf(" Maximum number of iterations reached.\n");
                printf(" \n");
            }
            else if (*info == 3)
            {
                printf(" \n");
                printf(" No shifts could be applied during implicit\n");
                printf(" Arnoldi update try increasing NCV.\n");
                printf(" \n");
            }

            if (iparam[5] > 0)
            {

                sseupd_(rvec, "A", &select[1], &d[1], &z[z_offset], ldz, sigma, bmat, n, which, nev, tol, &resid[1], ncv, &v[v_offset], ldv, &iparam[1], ipntr, &workd[1], &workl[1], lworkl, info);

                if (*info != 0)
                {

                    /* ---------------------------------- */
                    /* Check the documentation of sneupd. */
                    /* ---------------------------------- */

                    printf(" \n");
                    printf(" Error with _neupd = %d", (*info));
                    printf(" Check the documentation of _neupd \n");
                    printf(" \n");
                    return ierr;
                }
            }
        }

        return ierr;
    }

    /* -------------------------------------- */
    /* L O O P  B A C K to call SSAUPD again. */
    /* -------------------------------------- */

    goto L90;

    return 0;
} /* ssband_ */
