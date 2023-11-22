/* EXAMPLES\BAND\snband.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

/* Table of constant values */

static a_int c__9 = 9;
static a_int c__1 = 1;
static float c_b83 = 1.f;
static float c_b85 = 0.f;
static a_int c__3 = 3;
/**
 * \BeginDoc
 *
 * \Name: snband
 *
 * \Description:
 *
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
 *  The approximate eigenvalues and vectors are commonly called Ritz
 *  values and Ritz vectors respectively.  They are referred to as such
 *  in the comments that follow.  The computed orthonormal basis for the
 *  invariant subspace corresponding to these Ritz values is referred to as a
 *  Schur basis.
 *
 *  snband can be called with one of the following modes:
 *
 *  Mode 1:  A*z = lambda*z.
 *           ===> OP = A  and  B = I.
 *
 *  Mode 2:  A*z = lambda*M*z, M symmetric positive definite
 *           ===> OP = inv[M]*A  and  B = M.
 *
 *  Mode 3:  A*z = lambda*M*z, M symmetric semi-definite
 *           ===> OP = Real_Part{ inv[A - sigma*M]*M }  and  B = M.
 *           ===> shift-and-invert mode (in real arithmetic)
 *           If OP*z = amu*z, then
 *           amu = 1/2 * [ 1/(lambda-sigma) + 1/(lambda-conjg(sigma)) ].
 *           Note: If sigma is real, i.e. imaginary part of sigma is zero;
 *                 Real_Part{ inv[A - sigma*M]*M } == inv[A - sigma*M]*M
 *                 amu == 1/(lambda-sigma).
 *
 *  Mode 4:  A*z = lambda*M*z, M symmetric semi-definite
 *           ===> OP = Imaginary_Part{ inv[A - sigma*M]*M }  and  B = M.
 *           ===> shift-and-invert mode (in real arithmetic)
 *           If OP*z = amu*z, then
 *           amu = 1/2i * [ 1/(lambda-sigma) - 1/(lambda-conjg(sigma)) ].
 *
 *  The choice of mode must be specified in IPARAM(7) defined below.
 *
 * \Usage
 *   call snband
 *      ( RVEC, HOWMNY, SELECT, DR, DI, Z, LDZ, SIGMAR, SIGMAI,
 *        WORKEV, V, N, AB, MB, LDA, RFAC, CFAC, KL, KU, WHICH,
 *        BMAT, NEV, TOL, RESID, NCV, V, LDV, IPARAM, WORKD,
 *        WORKL, LWORKL, WORKC, IWORK, INFO )
 *
 * \Arguments
 *
 *  RVEC    LOGICAL  (INPUT)
 *          Specifies whether a basis for the invariant subspace corresponding
 *          to the converged Ritz value approximations for the eigenproblem
 *          A*z = lambda*B*z is computed.
 *
 *             RVEC = .FALSE.     Compute Ritz values only.
 *
 *             RVEC = .TRUE.      Compute the Ritz vectors or Schur vectors.
 *                                See Remarks below.
 *
 *  HOWMNY  Character*1  (INPUT)
 *          Specifies the form of the basis for the invariant subspace
 *          corresponding to the converged Ritz values that is to be computed.
 *
 *          = 'A': Compute NEV Ritz vectors;
 *          = 'P': Compute NEV Schur vectors;
 *          = 'S': compute some of the Ritz vectors, specified
 *                 by the logical array SELECT.
 *
 *  SELECT  Logical array of dimension NCV.  (INPUT)
 *          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
 *          computed. To select the Ritz vector corresponding to a
 *          Ritz value (DR(j), DI(j)), SELECT(j) must be set to .TRUE..
 *          If HOWMNY = 'A' or 'P', SELECT is used as internal workspace.
 *
 *  DR      Real array of dimension NEV+1.  (OUTPUT)
 *          On exit, DR contains the real part of the Ritz value approximations
 *          to the eigenvalues of A*z = lambda*B*z.
 *
 *  DI      Real array of dimension NEV+1.  (OUTPUT)
 *          On exit, DI contains the imaginary part of the Ritz value
 *          approximations to the eigenvalues of A*z = lambda*B*z associated
 *          with DR.
 *
 *          NOTE: When Ritz values are complex, they will come in complex
 *                conjugate pairs.  If eigenvectors are requested, the
 *                corresponding Ritz vectors will also come in conjugate
 *                pairs and the real and imaginary parts of these are
 *                represented in two consecutive columns of the array Z
 *                (see below).
 *
 *  Z       Real N by NEV+1 array if RVEC = .TRUE. and HOWMNY = 'A'. (OUTPUT)
 *          On exit,
 *          if RVEC = .TRUE. and HOWMNY = 'A', then the columns of
 *          Z represent approximate eigenvectors (Ritz vectors) corresponding
 *          to the NCONV=IPARAM(5) Ritz values for eigensystem
 *          A*z = lambda*B*z computed by SNAUPD.
 *
 *          The complex Ritz vector associated with the Ritz value
 *          with positive imaginary part is stored in two consecutive
 *          columns.  The first column holds the real part of the Ritz
 *          vector and the second column holds the imaginary part.  The
 *          Ritz vector associated with the Ritz value with negative
 *          imaginary part is simply the complex conjugate of the Ritz vector
 *          associated with the positive imaginary part.
 *
 *          If  RVEC = .FALSE. or HOWMNY = 'P', then Z is not referenced.
 *
 *          NOTE: If if RVEC = .TRUE. and a Schur basis is not required,
 *          the array Z may be set equal to first NEV+1 columns of the Arnoldi
 *          basis array V computed by SNAUPD.  In this case the Arnoldi basis
 *          will be destroyed and overwritten with the eigenvector basis.
 *
 *  LDZ     Integer.  (INPUT)
 *          The leading dimension of the array Z.  If Ritz vectors are
 *          desired, then  LDZ >= max( 1, N ).  In any case,  LDZ >= 1.
 *
 *  SIGMAR  Real  (INPUT)
 *          If IPARAM(7) = 3 or 4, represents the real part of the shift.
 *          Not referenced if IPARAM(7) = 1 or 2.
 *
 *  SIGMAI  Real  (INPUT)
 *          If IPARAM(7) = 3 or 4, represents the imaginary part of the
 *          shift.
 *          Not referenced if IPARAM(7) = 1 or 2.
 *
 *  WORKEV  Real work array of dimension 3*NCV.  (WORKSPACE)
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
 *          Leading dimension of AB, MB, RFAC and CFAC.
 *
 *  RFAC    Real array of LDA by N. (WORKSPACE/OUTPUT)
 *          RFAC is used to store the LU factors of MB when IPARAM(7) = 2
 *          is invoked.  It is used to store the LU factors of
 *          (A-sigma*M) when IPARAM(7) = 3 is invoked with a real shift.
 *          It is not referenced when IPARAM(7) = 1 or 4.
 *
 *  CFAC    Complex array of LDA by N. (WORKSPACE/OUTPUT)
 *          CFAC is used to store (A-SIGMA*M) and its LU factors
 *          when IPARAM(7) = 3 or 4 are used with a complex shift SIGMA.
 *          On exit, it contains the LU factors of (A-SIGMA*M).
 *          It is not referenced when IPARAM(7) = 1 or 2.
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
 *            'LR' -> want the NEV eigenvalues of largest real part.
 *            'SR' -> want the NEV eigenvalues of smallest real part.
 *            'LI' -> want the NEV eigenvalues of largest imaginary part.
 *            'SI' -> want the NEV eigenvalues of smallest imaginary part.
 *
 *          When IPARAM(7) = 3 or 4, WHICH should be set to 'LM' only.
 *
 *  BMAT    Character*1.  (INPUT)
 *          BMAT specifies the type of the matrix B that defines the
 *          semi-inner product for the operator OP.
 *          BMAT = 'I' -> standard eigenvalue problem A*z = lambda*z
 *          BMAT = 'G' -> generalized eigenvalue problem A*z = lambda*M*z
 *  NEV     Integer. (INPUT)
 *          Number of eigenvalues to be computed.
 *
 *  TOL     Real scalar.  (INPUT)
 *          Stopping criteria: the relative accuracy of the Ritz value
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
 *          Represents the dimension of the Arnoldi basis constructed
 *          by snaupd for OP.
 *
 *  V       Real array N by NCV+1.  (OUTPUT)
 *          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
 *                       represent approximate Schur vectors that span the
 *                       desired invariant subspace.
 *          NOTE: The array Z may be set equal to first NEV+1 columns of the
 *          Arnoldi basis vector array V computed by SNAUPD. In this case
 *          if RVEC = .TRUE. and HOWMNY='A', then the first NCONV=IPARAM(5)
 *          are the desired Ritz vectors.
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
 *           ----------------------------------------------------------
 *          ISHIFT = 1: exact shift with respect to the current
 *                      Hessenberg matrix H.  This is equivalent to
 *                      restarting the iteration from the beginning
 *                      after updating the starting vector with a linear
 *                      combination of Ritz vectors associated with the
 *                      "wanted" eigenvalues.
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
 *
 *          IPARAM(6) = IUPD
 *          Not referenced. Implicit restarting is ALWAYS used.
 *
 *          IPARAM(7) = IPARAM(7):
 *          On INPUT determines what type of eigenproblem is being solved.
 *          Must be 1,2,3,4; See under \Description of snband for the
 *          four modes available.
 *
 *          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
 *          OUTPUT: NUMOP  = total number of OP*z operations,
 *                  NUMOPB = total number of B*z operations if BMAT='G',
 *                  NUMREO = total number of steps of re-orthogonalization.
 *
 * WORKD    Real work array of length at least 3*n. (WORKSPACE)
 *
 * WORKL    Real work array of length LWORKL. (WORKSPACE)
 *
 * LWORKL   Integer.  (INPUT)
 *          LWORKL must be at least 3*NCV**2 + 6*NCV.
 *
 * WORKC    Complex array of length N. (WORKSPACE)
 *          Workspace used when IPARAM(7) = 3 or 4 for storing a temporary
 *          complex vector.
 *
 * IWORK    Integer array of dimension at least N. (WORKSPACE)
 *          Used when IPARAM(7)=2,3,4 to store the pivot information in the
 *          factorization of M or (A-SIGMA*M).
 *
 * INFO     Integer.  (INPUT/OUTPUT)
 *          Error flag on output.
 *          =  0: Normal exit.
 *          =  1: The Schur form computed by LAPACK routine slahqr
 *                could not be reordered by LAPACK routine strsen.
 *                Re-enter subroutine SNEUPD with IPARAM(5)=NCV and
 *                increase the size of the arrays DR and DI to have
 *                dimension at least NCV and allocate at least NCV
 *                columns for Z. NOTE: Not necessary if Z and V share
 *                the same space. Please notify the authors.
 *
 *          = -1: N must be positive.
 *          = -2: NEV must be positive.
 *          = -3: NCV-NEV >= 2 and less than or equal to N.
 *          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
 *          = -6: BMAT must be one of 'I' or 'G'.
 *          = -7: Length of private work WORKL array is not sufficient.
 *          = -8: Error return from calculation of a real Schur form.
 *                Informational error from LAPACK routine slahqr.
 *          = -9: Error return from calculation of eigenvectors.
 *                Informational error from LAPACK routine strevc.
 *          = -10: IPARAM(7) must be 1,2,3,4.
 *          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
 *          = -12: HOWMNY = 'S' not yet implemented
 *          = -13: HOWMNY must be one of 'A' or 'P'
 *          = -14: SNAUPD did not find any eigenvalues to sufficient
 *                 accuracy.
 *          = -15: Overflow occurs when we try to transform the Ritz
 *                 values returned from SNAUPD to those of the original
 *                 problem using Rayleigh Quotient.
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
 * \Routines called:
 *     snaupd  ARPACK reverse communication interface routine.
 *     sneupd  ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     sgbtrf  LAPACK band matrix factorization routine.
 *     sgbtrs  LAPACK band linear system solve routine.
 *     cgbtrf  LAPACK complex band matrix factorization routine.
 *     cgbtrs  LAPACK complex linear system solve routine.
 *     slacpy  LAPACK matrix copy routine.
 *     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     slamch  LAPACK routine to compute the underflow threshold.
 *     scopy   Level 1 BLAS that copies one vector to another.
 *     sdot    Level 1 BLAS that computes the dot product of two vectors.
 *     snrm2   Level 1 BLAS that computes the norm of a vector.
 *     sgbmv   Level 2 BLAS that computes the band matrix vector product.
 *
 * \Remarks
 *
 *  1. Currently only HOWMNY = 'A' and 'P' are implemented.
 *
 *     Let X' denote the transpose of X.
 *
 *  2. Schur vectors are an orthogonal representation for the basis of
 *     Ritz vectors. Thus, their numerical properties are often superior.
 *     If RVEC = .TRUE. then the relationship
 *             A * V(:,1:IPARAM(5)) = V(:,1:IPARAM(5)) * T, and
 *     V(:,1:IPARAM(5))' * V(:,1:IPARAM(5)) = I are approximately satisfied.
 *     Here T is the leading submatrix of order IPARAM(5) of the real
 *     upper quasi-triangular matrix stored workl(ipntr(12)). That is,
 *     T is block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;
 *     each 2-by-2 diagonal block has its diagonal elements equal and its
 *     off-diagonal elements of opposite sign.  Corresponding to each 2-by-2
 *     diagonal block is a complex conjugate pair of Ritz values. The real
 *     Ritz values are stored on the diagonal of T.
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
 * FILE: nband.F   SID: 2.3   DATE OF SID: 10/17/00   RELEASE: 2
 *
 * \EndLib
 */
int snband_(a_bool *rvec, char *howmny, a_bool *select, float *dr, float *di, float *z, a_int *ldz, float *sigmar, float *sigmai, float *workev, a_int *n, float *ab, float *mb, a_int *lda, float *rfac, a_fcomplex *cfac, a_int *kl, a_int *ku, char *which, char *bmat, a_int *nev, float *tol, float *resid, a_int *ncv, float *v, a_int *ldv, a_int *iparam, float *workd, float *workl, a_int *lworkl, a_fcomplex *workc, a_int *iwork, a_int *info)
{
    /* System generated locals */
    a_int v_dim1, v_offset, z_dim1, z_offset, ab_dim1, ab_offset, mb_dim1, mb_offset, rfac_dim1, rfac_offset, cfac_dim1, cfac_offset, i__1, i__2, i__3, i__4, i__5;
    float r__1, r__2, r__3;
    a_fcomplex q__1, q__2;

    /* Local variables */
    a_int i, j, ido;
    float deni;
    a_int imid;
    float denr;
    a_int ibot, ierr = 0;
    a_int itop, type;
    float numr;
    float dmdul;
    a_bool first;
    a_int ipntr[14];
    float safmin;

    /* ------------------------------ */
    /* safmin = safe minimum is such  */
    /* that 1/sfmin does not overflow */
    /* ------------------------------ */

    /* Parameter adjustments */
    --select;
    --dr;
    --di;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z -= z_offset;
    --workev;
    cfac_dim1 = *lda;
    cfac_offset = 1 + cfac_dim1;
    cfac -= cfac_offset;
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
    --workc;
    --iwork;

    safmin = slamch_("S");

    /* -------------------------------------------------------------- */
    /* Set type of the problem to be solved. Check consistency        */
    /* between BMAT and IPARAM(7).                                    */
    /* type = 1 --> Solving standard problem in regular mode.         */
    /* type = 2 --> Solving standard problem in shift-invert mode.    */
    /* type = 3 --> Solving generalized problem in regular mode.      */
    /* type = 4 --> Solving generalized problem in shift-invert mode. */
    /* type = 5 --> Solving standard problem in shift-invert mode     */
    /*              using iparam(7) = 4 in SNAUPD.                    */
    /* type = 6 --> Solving generalized problem in shift-invert mode. */
    /*              using iparam(7) = 4 in SNAUPD.                    */
    /* -------------------------------------------------------------- */

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
    else if (iparam[7] == 4 && *bmat == 'I')
    {
        type = 5;
    }
    else if (iparam[7] == 4 && *bmat == 'G')
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

    /* -------------------------------- */
    /* When type = 5,6 are used, sigmai */
    /* must be nonzero.                 */
    /* -------------------------------- */

    if (type == 5 || type == 6)
    {
        if (*sigmai == 0.f)
        {
            printf(" \n");
            printf("_NBAND: sigmai must be nonzero when type 5 or 6                   is used. \n");
            printf(" \n");
            return ierr;
        }
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

    if (type == 2 || type == 5)
    {

        /* ----------------------------- */
        /* Solving a standard eigenvalue */
        /* problem in shift-invert mode. */
        /* Factor (A-sigma*I).           */
        /* ----------------------------- */

        if (*sigmai == 0.f)
        {

            /* --------------------------------- */
            /* Construct (A-sigmar*I) and factor */
            /* in real arithmetic.               */
            /* --------------------------------- */

            slacpy_("A", &ibot, n, &ab[ab_offset], lda, &rfac[rfac_offset], lda);
            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                rfac[imid + j * rfac_dim1] = ab[imid + j * ab_dim1] - *sigmar;
            }
            sgbtrf_(n, n, kl, ku, &rfac[rfac_offset], lda, &iwork[1], &ierr);
            if (ierr != 0)
            {
                printf(" \n");
                printf(" _NBAND: Error with _gbtrf. \n");
                printf(" \n");
                return ierr;
            }
        }
        else
        {

            /* --------------------------------- */
            /* Construct (A-sigmar*I) and factor */
            /* in COMPLEX arithmetic.            */
            /* --------------------------------- */

            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                i__2 = ibot;
                for (i = itop; i <= i__2; ++i)
                {
                    i__3 = i + j * cfac_dim1;
                    i__4 = i + j * ab_dim1;
                    q__1.r = ab[i__4], q__1.i = 0.f;
                    cfac[i__3].r = q__1.r, cfac[i__3].i = q__1.i;
                }
            }

            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                i__2 = imid + j * cfac_dim1;
                i__3 = imid + j * cfac_dim1;
                q__2.r = *sigmar, q__2.i = *sigmai;
                q__1.r = cfac[i__3].r - q__2.r, q__1.i = cfac[i__3].i - q__2.i;
                cfac[i__2].r = q__1.r, cfac[i__2].i = q__1.i;
            }

            cgbtrf_(n, n, kl, ku, &cfac[cfac_offset], lda, &iwork[1], &ierr);
            if (ierr != 0)
            {
                printf(" \n");
                printf(" _NBAND: Error with _gbtrf. \n");
                printf(" \n");
                return ierr;
            }
        }
    }
    else if (type == 3)
    {

        /* --------------------------------------------- */
        /* Solving generalized eigenvalue problem in     */
        /* regular mode. Copy M to rfac, and call LAPACK */
        /* routine sgbtrf to factor M.                   */
        /* --------------------------------------------- */

        slacpy_("A", &ibot, n, &mb[mb_offset], lda, &rfac[rfac_offset], lda);
        sgbtrf_(n, n, kl, ku, &rfac[rfac_offset], lda, &iwork[1], &ierr);
        if (ierr != 0)
        {
            printf(" \n");
            printf("_NBAND:  Error with _gbtrf.\n");
            printf(" \n");
            return ierr;
        }
    }
    else if (type == 4 || type == 6)
    {

        /* ----------------------------------------- */
        /* Solving generalized eigenvalue problem in */
        /* shift-invert mode.                        */
        /* ----------------------------------------- */

        if (*sigmai == 0.f)
        {

            /* ------------------------------------------ */
            /* Construct (A - sigma*M) and factor in real */
            /* arithmetic.                                */
            /* ------------------------------------------ */

            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                i__2 = ibot;
                for (i = itop; i <= i__2; ++i)
                {
                    rfac[i + j * rfac_dim1] = ab[i + j * ab_dim1] - *sigmar * mb[i + j * mb_dim1];
                }
            }

            sgbtrf_(n, n, kl, ku, &rfac[rfac_offset], lda, &iwork[1], &ierr);
            if (ierr != 0)
            {
                printf(" \n");
                printf("_NBAND: Error with _gbtrf.\n");
                printf(" \n");
                return ierr;
            }
        }
        else
        {

            /* --------------------------------------------- */
            /* Construct (A - sigma*M) and factor in complex */
            /* arithmetic.                                   */
            /* --------------------------------------------- */

            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                i__2 = ibot;
                for (i = itop; i <= i__2; ++i)
                {
                    i__3 = i + j * cfac_dim1;
                    r__1 = ab[i + j * ab_dim1] - *sigmar * mb[i + j * mb_dim1];
                    r__2 = -(*sigmai) * mb[i + j * mb_dim1];
                    q__1.r = r__1, q__1.i = r__2;
                    cfac[i__3].r = q__1.r, cfac[i__3].i = q__1.i;
                }
            }

            cgbtrf_(n, n, kl, ku, &cfac[cfac_offset], lda, &iwork[1], &ierr);
            if (ierr != 0)
            {
                printf(" \n");
                printf("_NBAND: Error with _gbtrf.\n");
                printf(" \n");
                return ierr;
            }
        }
    }

    /* ------------------------------------------ */
    /*  M A I N   L O O P (reverse communication) */
    /* ------------------------------------------ */

L90:

    snaupd_(&ido, bmat, n, which, nev, tol, &resid[1], ncv, &v[v_offset], ldv, &iparam[1], ipntr, &workd[1], &workl[1], lworkl, info);

    if (ido == -1)
    {

        if (type == 1)
        {

            /* -------------------------- */
            /* Perform  y <--- OP*x = A*x */
            /* -------------------------- */

            sgbmv_("N", n, n, kl, ku, &c_b83, &ab[itop + ab_dim1], lda, &workd[ipntr[0]], &c__1, &c_b85, &workd[ipntr[1]], &c__1);
        }
        else if (type == 2)
        {

            if (*sigmai == 0.f)
            {

                /* -------------------------------- */
                /* Shift is real.  Perform          */
                /* y <--- OP*x = inv[A-sigmar*I]*x  */
                /* to force the starting vector     */
                /* into the range of OP.            */
                /* -------------------------------- */

                scopy_(n, &workd[ipntr[0]], &c__1, &workd[ipntr[1]], &c__1);
                sgbtrs_("N", n, kl, ku, &c__1, &rfac[rfac_offset], lda, &iwork[1], &workd[ipntr[1]], n, &ierr);
                if (ierr != 0)
                {
                    printf(" \n");
                    printf(" _NBAND: Error with _bgtrs. \n");
                    printf(" \n");
                    return ierr;
                }
            }
            else
            {

                /* ------------------------------------------ */
                /* Shift is COMPLEX. Perform                  */
                /* y <--- OP*x = Real_Part{inv[A-sigma*I]*x}  */
                /* to force the starting vector into the      */
                /* range of OP.                               */
                /* ------------------------------------------ */

                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    i__2 = j;
                    i__3 = ipntr[0] + j - 1;
                    q__1.r = workd[i__3], q__1.i = 0.f;
                    workc[i__2].r = q__1.r, workc[i__2].i = q__1.i;
                }

                cgbtrs_("N", n, kl, ku, &c__1, &cfac[cfac_offset], lda, &iwork[1], &workc[1], n, &ierr);
                if (ierr != 0)
                {
                    printf(" \n");
                    printf(" _NBAND: Error with _gbtrs. \n");
                    printf(" \n");
                    return ierr;
                }

                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    i__2 = j;
                    workd[ipntr[1] + j - 1] = workc[i__2].r;
                }
            }
        }
        else if (type == 3)
        {

            /* --------------------------------- */
            /* Perform  y <--- OP*x = inv[M]*A*x */
            /* to force the starting vector into */
            /* the range of OP.                  */
            /* --------------------------------- */

            sgbmv_("N", n, n, kl, ku, &c_b83, &ab[itop + ab_dim1], lda, &workd[ipntr[0]], &c__1, &c_b85, &workd[ipntr[1]], &c__1);

            sgbtrs_("N", n, kl, ku, &c__1, &rfac[rfac_offset], lda, &iwork[1], &workd[ipntr[1]], n, &ierr);
            if (ierr != 0)
            {
                printf(" \n");
                printf("_NBAND: Error with _bgtrs.\n");
                printf(" \n");
                return ierr;
            }
        }
        else if (type == 4)
        {

            /* --------------------------------------- */
            /* Perform y <-- OP*x                      */
            /*         = Real_part{inv[A-SIGMA*M]*M}*x */
            /* to force the starting vector into the   */
            /* range of OP.                            */
            /* --------------------------------------- */

            sgbmv_("N", n, n, kl, ku, &c_b83, &mb[itop + mb_dim1], lda, &workd[ipntr[0]], &c__1, &c_b85, &workd[ipntr[1]], &c__1);

            if (*sigmai == 0.f)
            {

                /* ------------------- */
                /* Shift is real, stay */
                /* in real arithmetic. */
                /* ------------------- */

                sgbtrs_("N", n, kl, ku, &c__1, &rfac[rfac_offset], lda, &iwork[1], &workd[ipntr[1]], n, &ierr);
                if (ierr != 0)
                {
                    printf(" \n");
                    printf("_NBAND: Error with _gbtrs.\n");
                    printf(" \n");
                    return ierr;
                }
            }
            else
            {

                /* ------------------------ */
                /* Goto complex arithmetic. */
                /* ------------------------ */

                i__1 = *n;
                for (i = 1; i <= i__1; ++i)
                {
                    i__2 = i;
                    i__3 = ipntr[1] + i - 1;
                    q__1.r = workd[i__3], q__1.i = 0.f;
                    workc[i__2].r = q__1.r, workc[i__2].i = q__1.i;
                }

                cgbtrs_("N", n, kl, ku, &c__1, &cfac[cfac_offset], lda, &iwork[1], &workc[1], n, &ierr);
                if (ierr != 0)
                {
                    printf(" \n");
                    printf("_NBAND: Error with _gbtrs.\n");
                    printf(" \n");
                    return ierr;
                }

                i__1 = *n;
                for (i = 1; i <= i__1; ++i)
                {
                    i__2 = i;
                    workd[ipntr[1] + i - 1] = workc[i__2].r;
                }
            }
        }
        else if (type == 5)
        {

            /* ------------------------------------- */
            /* Perform y <-- OP*x                    */
            /*    = Imaginary_part{inv[A-SIGMA*I]}*x */
            /* to force the starting vector into the */
            /* range of OP.                          */
            /* ------------------------------------- */

            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                i__2 = j;
                i__3 = ipntr[0] + j - 1;
                q__1.r = workd[i__3], q__1.i = 0.f;
                workc[i__2].r = q__1.r, workc[i__2].i = q__1.i;
            }

            cgbtrs_("N", n, kl, ku, &c__1, &cfac[cfac_offset], lda, &iwork[1], &workc[1], n, &ierr);
            if (ierr != 0)
            {
                printf(" \n");
                printf(" _NBAND: Error with _gbtrs. \n");
                printf(" \n");
                return ierr;
            }

            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                workd[ipntr[1] + j - 1] = workc[j].i;
            }
        }
        else if (type == 6)
        {

            /* -------------------------------------- */
            /* Perform y <-- OP*x                     */
            /*       Imaginary_part{inv[A-SIGMA*M]*M} */
            /* to force the starting vector into the  */
            /* range of OP.                           */
            /* -------------------------------------- */

            sgbmv_("N", n, n, kl, ku, &c_b83, &mb[itop + mb_dim1], lda, &workd[ipntr[0]], &c__1, &c_b85, &workd[ipntr[1]], &c__1);

            i__1 = *n;
            for (i = 1; i <= i__1; ++i)
            {
                i__2 = i;
                i__3 = ipntr[1] + i - 1;
                q__1.r = workd[i__3], q__1.i = 0.f;
                workc[i__2].r = q__1.r, workc[i__2].i = q__1.i;
            }

            cgbtrs_("N", n, kl, ku, &c__1, &cfac[cfac_offset], lda, &iwork[1], &workc[1], n, &ierr);
            if (ierr != 0)
            {
                printf(" \n");
                printf("_NBAND: Error with _gbtrs.\n");
                printf(" \n");
                return ierr;
            }

            i__1 = *n;
            for (i = 1; i <= i__1; ++i)
            {
                workd[ipntr[1] + i - 1] = workc[i].i;
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

            sgbmv_("N", n, n, kl, ku, &c_b83, &ab[itop + ab_dim1], lda, &workd[ipntr[0]], &c__1, &c_b85, &workd[ipntr[1]], &c__1);
        }
        else if (type == 2)
        {

            if (*sigmai == 0.f)
            {

                /* -------------------------------- */
                /* Shift is real.  Perform          */
                /* y <--- OP*x = inv[A-sigmar*I]*x. */
                /* -------------------------------- */

                scopy_(n, &workd[ipntr[0]], &c__1, &workd[ipntr[1]], &c__1);
                sgbtrs_("N", n, kl, ku, &c__1, &rfac[rfac_offset], lda, &iwork[1], &workd[ipntr[1]], n, &ierr);
            }
            else
            {

                /* ---------------------------------------- */
                /* Shift is COMPLEX. Perform                */
                /* y <-- OP*x = Real_Part{inv[A-sigma*I]*x} */
                /* in COMPLEX arithmetic.                   */
                /* ---------------------------------------- */

                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    i__2 = j;
                    i__3 = ipntr[0] + j - 1;
                    q__1.r = workd[i__3], q__1.i = 0.f;
                    workc[i__2].r = q__1.r, workc[i__2].i = q__1.i;
                }

                cgbtrs_("N", n, kl, ku, &c__1, &cfac[cfac_offset], lda, &iwork[1], &workc[1], n, &ierr);
                if (ierr != 0)
                {
                    printf(" \n");
                    printf("_NBAND: Error with _gbtrs.\n");
                    printf(" \n");
                    return ierr;
                }

                i__1 = *n;
                for (j = 1; j <= i__1; ++j)
                {
                    i__2 = j;
                    workd[ipntr[1] + j - 1] = workc[i__2].r;
                }
            }
        }
        else if (type == 3)
        {

            /* --------------------------------- */
            /* Perform  y <--- OP*x = inv[M]*A*x */
            /* --------------------------------- */

            sgbmv_("N", n, n, kl, ku, &c_b83, &ab[itop + ab_dim1], lda, &workd[ipntr[0]], &c__1, &c_b85, &workd[ipntr[1]], &c__1);

            sgbtrs_("N", n, kl, ku, &c__1, &rfac[rfac_offset], lda, &iwork[1], &workd[ipntr[1]], n, &ierr);
            if (ierr != 0)
            {
                printf(" \n");
                printf("_NBAND: Error with _bgtrs.\n");
                printf(" \n");
                return ierr;
            }
        }
        else if (type == 4)
        {

            /* ------------------------------------ */
            /* Perform  y <-- inv(A-sigma*M)*(M*x). */
            /* (M*x) has been computed and stored   */
            /* in workd(ipntr(3)).                  */
            /* ------------------------------------ */

            if (*sigmai == 0.f)
            {

                /* ---------------------- */
                /* Shift is real, stay in */
                /* real arithmetic.       */
                /* ---------------------- */

                scopy_(n, &workd[ipntr[2]], &c__1, &workd[ipntr[1]], &c__1);
                sgbtrs_("N", n, kl, ku, &c__1, &rfac[rfac_offset], lda, &iwork[1], &workd[ipntr[1]], n, &ierr);
                if (ierr != 0)
                {
                    printf(" \n");
                    printf("_NBAND: Error with _gbtrs.\n");
                    printf(" \n");
                    return ierr;
                }
            }
            else
            {

                /* ------------------------- */
                /* Go to COMPLEX arithmetic. */
                /* ------------------------- */

                i__1 = *n;
                for (i = 1; i <= i__1; ++i)
                {
                    i__2 = i;
                    i__3 = ipntr[2] + i - 1;
                    q__1.r = workd[i__3], q__1.i = 0.f;
                    workc[i__2].r = q__1.r, workc[i__2].i = q__1.i;
                }

                cgbtrs_("N", n, kl, ku, &c__1, &cfac[cfac_offset], lda, &iwork[1], &workc[1], n, &ierr);
                if (ierr != 0)
                {
                    printf(" \n");
                    printf("_NBAND: Error in _gbtrs.\n");
                    printf(" \n");
                    return ierr;
                }

                i__1 = *n;
                for (i = 1; i <= i__1; ++i)
                {
                    i__2 = i;
                    workd[ipntr[1] + i - 1] = workc[i__2].r;
                }
            }
        }
        else if (type == 5)
        {

            /* ------------------------------------- */
            /* Perform y <-- OP*x                    */
            /*    = Imaginary_part{inv[A-SIGMA*I]*x} */
            /* ------------------------------------- */

            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                i__2 = j;
                i__3 = ipntr[0] + j - 1;
                q__1.r = workd[i__3], q__1.i = 0.f;
                workc[i__2].r = q__1.r, workc[i__2].i = q__1.i;
            }

            cgbtrs_("N", n, kl, ku, &c__1, &cfac[cfac_offset], lda, &iwork[1], &workc[1], n, &ierr);
            if (ierr != 0)
            {
                printf(" \n");
                printf(" _NBAND: Error with _gbtrs. \n");
                printf(" \n");
                return ierr;
            }

            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                workd[ipntr[1] + j - 1] = workc[j].i;
            }
        }
        else if (type == 6)
        {

            /* --------------------------------------- */
            /* Perform y <-- OP*x                      */
            /*   = Imaginary_part{inv[A-SIGMA*M]*M}*x. */
            /* --------------------------------------- */

            i__1 = *n;
            for (i = 1; i <= i__1; ++i)
            {
                i__2 = i;
                i__3 = ipntr[2] + i - 1;
                q__1.r = workd[i__3], q__1.i = 0.f;
                workc[i__2].r = q__1.r, workc[i__2].i = q__1.i;
            }

            cgbtrs_("N", n, kl, ku, &c__1, &cfac[cfac_offset], lda, &iwork[1], &workc[1], n, &ierr);
            if (ierr != 0)
            {
                printf(" \n");
                printf("_NBAND: Error with _gbtrs.\n");
                printf(" \n");
                return ierr;
            }

            i__1 = *n;
            for (i = 1; i <= i__1; ++i)
            {
                workd[ipntr[1] + i - 1] = workc[i].i;
            }
        }
    }
    else if (ido == 2)
    {

        /* ------------------ */
        /* Perform y <-- M*x  */
        /* Not used when      */
        /* type = 1,2.        */
        /* ------------------ */

        sgbmv_("N", n, n, kl, ku, &c_b83, &mb[itop + mb_dim1], lda, &workd[ipntr[0]], &c__1, &c_b85, &workd[ipntr[1]], &c__1);
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
            /* documentation in SNAUPD  */
            /* ------------------------ */

            printf(" \n");
            printf(" Error with _naupd info = %d", (*info));
            printf(" Check the documentation of _naupd \n");
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

                sneupd_(rvec, "A", &select[1], &dr[1], &di[1], &z[z_offset], ldz, sigmar, sigmai, &workev[1], bmat, n, which, nev, tol, &resid[1], ncv, &v[v_offset], ldv, &iparam[1], ipntr, &workd[1], &workl[1], lworkl, info);

                if (*info != 0)
                {

                    /* ---------------------------------- */
                    /* Check the documentation of SNEUPD. */
                    /* ---------------------------------- */

                    printf(" \n");
                    printf(" Error with _neupd = %d", (*info));
                    printf(" Check the documentation of _neupd \n");
                    printf(" \n");
                    return ierr;
                }
                else if (*sigmai != 0.f)
                {

                    if (type == 4 || type == 6)
                    {

                        first = TRUE_;
                        i__1 = iparam[5];
                        for (j = 1; j <= i__1; ++j)
                        {

                            /* -------------------------------- */
                            /* Use Rayleigh Quotient to recover */
                            /* eigenvalues of the original      */
                            /* generalized eigenvalue problem.  */
                            /* -------------------------------- */

                            if (di[j] == 0.f)
                            {

                                /* ------------------------------------ */
                                /* Eigenvalue is real. Compute          */
                                /* d = (x'*inv[A-sigma*M]*M*x) / (x'*x) */
                                /* ------------------------------------ */

                                sgbmv_("Nontranspose", n, n, kl, ku, &c_b83, &mb[itop + mb_dim1], lda, &z[j * z_dim1 + 1], &c__1, &c_b85, &workd[1], &c__1);
                                i__2 = *n;
                                for (i = 1; i <= i__2; ++i)
                                {
                                    i__3 = i;
                                    i__4 = i;
                                    q__1.r = workd[i__4], q__1.i = 0.f;
                                    workc[i__3].r = q__1.r, workc[i__3].i = q__1.i;
                                }
                                cgbtrs_("N", n, kl, ku, &c__1, &cfac[cfac_offset], lda, &iwork[1], &workc[1], n, info);
                                i__2 = *n;
                                for (i = 1; i <= i__2; ++i)
                                {
                                    i__3 = i;
                                    workd[i] = workc[i__3].r;
                                    workd[i + *n] = workc[i].i;
                                }
                                denr = sdot_(n, &z[j * z_dim1 + 1], &c__1, &workd[1], &c__1);
                                deni = sdot_(n, &z[j * z_dim1 + 1], &c__1, &workd[*n + 1], &c__1);
                                /* Computing 2nd power */
                                r__1 = snrm2_(n, &z[j * z_dim1 + 1], &c__1);
                                numr = r__1 * r__1;
                                /* Computing 2nd power */
                                r__1 = slapy2_(&denr, &deni);
                                dmdul = r__1 * r__1;
                                if (dmdul >= safmin)
                                {
                                    dr[j] = *sigmar + numr * denr / dmdul;
                                }
                                else
                                {

                                    /* ------------------- */
                                    /* dmdul is too small. */
                                    /* Exit to avoid       */
                                    /* overflow.           */
                                    /* ------------------- */

                                    *info = -15;
                                    return ierr;
                                }
                            }
                            else if (first)
                            {

                                /* ---------------------- */
                                /* Eigenvalue is complex. */
                                /* Compute the first one  */
                                /* of the conjugate pair. */
                                /* ---------------------- */

                                /* ----------- */
                                /* Compute M*x */
                                /* ----------- */

                                sgbmv_("Nontranspose", n, n, kl, ku, &c_b83, &mb[itop + mb_dim1], lda, &z[j * z_dim1 + 1], &c__1, &c_b85, &workd[1], &c__1);
                                sgbmv_("Nontranspose", n, n, kl, ku, &c_b83, &mb[itop + mb_dim1], lda, &z[(j + 1) * z_dim1 + 1], &c__1, &c_b85, &workd[*n + 1], &c__1);
                                i__2 = *n;
                                for (i = 1; i <= i__2; ++i)
                                {
                                    i__3 = i;
                                    i__4 = i;
                                    i__5 = i + *n;
                                    q__1.r = workd[i__4], q__1.i = workd[i__5];
                                    workc[i__3].r = q__1.r, workc[i__3].i = q__1.i;
                                }

                                /* -------------------------- */
                                /* Compute inv(A-sigma*M)*M*x */
                                /* -------------------------- */

                                cgbtrs_("N", n, kl, ku, &c__1, &cfac[cfac_offset], lda, &iwork[1], &workc[1], n, info);

                                /* ----------------------------- */
                                /* Compute x'*inv(A-sigma*M)*M*x */
                                /* ----------------------------- */

                                i__2 = *n;
                                for (i = 1; i <= i__2; ++i)
                                {
                                    i__3 = i;
                                    workd[i] = workc[i__3].r;
                                    workd[i + *n] = workc[i].i;
                                }
                                denr = sdot_(n, &z[j * z_dim1 + 1], &c__1, &workd[1], &c__1);
                                denr += sdot_(n, &z[(j + 1) * z_dim1 + 1], &c__1, &workd[*n + 1], &c__1);
                                deni = sdot_(n, &z[j * z_dim1 + 1], &c__1, &workd[*n + 1], &c__1);
                                deni -= sdot_(n, &z[(j + 1) * z_dim1 + 1], &c__1, &workd[1], &c__1);

                                /* -------------- */
                                /* Compute (x'*x) */
                                /* -------------- */

                                r__2 = snrm2_(n, &z[j * z_dim1 + 1], &c__1);
                                r__3 = snrm2_(n, &z[(j + 1) * z_dim1 + 1], &c__1);
                                /* Computing 2nd power */
                                r__1 = slapy2_(&r__2, &r__3);
                                numr = r__1 * r__1;

                                /* -------------------------------------- */
                                /* Compute (x'x) / (x'*inv(A-sigma*M)*Mx) */
                                /* -------------------------------------- */

                                /* Computing 2nd power */
                                r__1 = slapy2_(&denr, &deni);
                                dmdul = r__1 * r__1;
                                if (dmdul >= safmin)
                                {
                                    dr[j] = *sigmar + numr * denr / dmdul;
                                    di[j] = *sigmai - numr * deni / dmdul;
                                    first = FALSE_;
                                }
                                else
                                {

                                    /* ------------------- */
                                    /* dmdul is too small. */
                                    /* Exit to avoid       */
                                    /* overflow.           */
                                    /* ------------------- */

                                    *info = -15;
                                    return ierr;
                                }
                            }
                            else
                            {

                                /* ------------------------- */
                                /* Get the second eigenvalue */
                                /* of the conjugate pair by  */
                                /* taking the conjugate of   */
                                /* previous one.             */
                                /* ------------------------- */

                                dr[j] = dr[j - 1];
                                di[j] = -di[j - 1];
                                first = TRUE_;
                            }
                        }
                    }
                    else if (type == 2 || type == 5)
                    {

                        first = TRUE_;
                        i__1 = iparam[5];
                        for (j = 1; j <= i__1; ++j)
                        {

                            /* -------------------------------- */
                            /* Use Rayleigh Quotient to recover */
                            /* eigenvalues of the original      */
                            /* standard eigenvalue problem.     */
                            /* -------------------------------- */

                            if (di[j] == 0.f)
                            {

                                /* ----------------------------------- */
                                /* Eigenvalue is real. Compute         */
                                /* d = (x'*inv[A-sigma*I]*x) / (x'*x). */
                                /* ----------------------------------- */

                                i__2 = *n;
                                for (i = 1; i <= i__2; ++i)
                                {
                                    i__3 = i;
                                    i__4 = i + j * z_dim1;
                                    q__1.r = z[i__4], q__1.i = 0.f;
                                    workc[i__3].r = q__1.r, workc[i__3].i = q__1.i;
                                }
                                cgbtrs_("N", n, kl, ku, &c__1, &cfac[cfac_offset], lda, &iwork[1], &workc[1], n, info);
                                i__2 = *n;
                                for (i = 1; i <= i__2; ++i)
                                {
                                    i__3 = i;
                                    workd[i] = workc[i__3].r;
                                    workd[i + *n] = workc[i].i;
                                }
                                denr = sdot_(n, &z[j * z_dim1 + 1], &c__1, &workd[1], &c__1);
                                deni = sdot_(n, &z[j * z_dim1 + 1], &c__1, &workd[*n + 1], &c__1);
                                /* Computing 2nd power */
                                r__1 = snrm2_(n, &z[j * z_dim1 + 1], &c__1);
                                numr = r__1 * r__1;
                                /* Computing 2nd power */
                                r__1 = slapy2_(&denr, &deni);
                                dmdul = r__1 * r__1;
                                if (dmdul >= safmin)
                                {
                                    dr[j] = *sigmar + numr * denr / dmdul;
                                }
                                else
                                {

                                    /* ------------------- */
                                    /* dmdul is too small. */
                                    /* Exit to avoid       */
                                    /* overflow.           */
                                    /* ------------------- */

                                    *info = -15;
                                    return ierr;
                                }
                            }
                            else if (first)
                            {

                                /* ---------------------- */
                                /* Eigenvalue is complex. */
                                /* Compute the first one  */
                                /* of the conjugate pair. */
                                /* ---------------------- */

                                i__2 = *n;
                                for (i = 1; i <= i__2; ++i)
                                {
                                    i__3 = i;
                                    i__4 = i + j * z_dim1;
                                    i__5 = i + (j + 1) * z_dim1;
                                    q__1.r = z[i__4], q__1.i = z[i__5];
                                    workc[i__3].r = q__1.r, workc[i__3].i = q__1.i;
                                }

                                /* ------------------------- */
                                /* Compute inv[A-sigma*I]*x. */
                                /* ------------------------- */

                                cgbtrs_("N", n, kl, ku, &c__1, &cfac[cfac_offset], lda, &iwork[1], &workc[1], n, info);

                                /* --------------------------- */
                                /* Compute x'*inv(A-sigma*I)*x */
                                /* --------------------------- */

                                i__2 = *n;
                                for (i = 1; i <= i__2; ++i)
                                {
                                    i__3 = i;
                                    workd[i] = workc[i__3].r;
                                    workd[i + *n] = workc[i].i;
                                }
                                denr = sdot_(n, &z[j * z_dim1 + 1], &c__1, &workd[1], &c__1);
                                denr += sdot_(n, &z[(j + 1) * z_dim1 + 1], &c__1, &workd[*n + 1], &c__1);
                                deni = sdot_(n, &z[j * z_dim1 + 1], &c__1, &workd[*n + 1], &c__1);
                                deni -= sdot_(n, &z[(j + 1) * z_dim1 + 1], &c__1, &workd[1], &c__1);

                                /* -------------- */
                                /* Compute (x'*x) */
                                /* -------------- */

                                r__2 = snrm2_(n, &z[j * z_dim1 + 1], &c__1);
                                r__3 = snrm2_(n, &z[(j + 1) * z_dim1 + 1], &c__1);
                                /* Computing 2nd power */
                                r__1 = slapy2_(&r__2, &r__3);
                                numr = r__1 * r__1;

                                /* -------------------------------------- */
                                /* Compute (x'x) / (x'*inv(A-sigma*I)*x). */
                                /* -------------------------------------- */

                                /* Computing 2nd power */
                                r__1 = slapy2_(&denr, &deni);
                                dmdul = r__1 * r__1;
                                if (dmdul >= safmin)
                                {
                                    dr[j] = *sigmar + numr * denr / dmdul;
                                    di[j] = *sigmai - numr * deni / dmdul;
                                    first = FALSE_;
                                }
                                else
                                {

                                    /* ------------------- */
                                    /* dmdul is too small. */
                                    /* Exit to avoid       */
                                    /* overflow.           */
                                    /* ------------------- */

                                    *info = -15;
                                    return ierr;
                                }
                            }
                            else
                            {

                                /* ------------------------- */
                                /* Get the second eigenvalue */
                                /* of the conjugate pair by  */
                                /* taking the conjugate of   */
                                /* previous one.             */
                                /* ------------------------- */

                                dr[j] = dr[j - 1];
                                di[j] = -di[j - 1];
                                first = TRUE_;
                            }
                        }
                    }
                }
            }
        }

        return ierr;
    }

    /* -------------------------------------- */
    /* L O O P  B A C K to call SNAUPD again. */
    /* -------------------------------------- */

    goto L90;

    return 0;
} /* snband_ */
