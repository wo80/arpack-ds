/* EXAMPLES\BAND\dnband.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

static a_int i_one = 1;
static double one = 1.;
static double zero = 0.;

/**
 * \BeginDoc
 *
 * \Name: dnband
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
 *  dnband can be called with one of the following modes:
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
 *   call dnband
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
 *  DR      Double precision array of dimension NEV+1.  (OUTPUT)
 *          On exit, DR contains the real part of the Ritz value approximations
 *          to the eigenvalues of A*z = lambda*B*z.
 *
 *  DI      Double precision array of dimension NEV+1.  (OUTPUT)
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
 *          A*z = lambda*B*z computed by DNAUPD.
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
 *          basis array V computed by DNAUPD.  In this case the Arnoldi basis
 *          will be destroyed and overwritten with the eigenvector basis.
 *
 *  LDZ     Integer.  (INPUT)
 *          The leading dimension of the array Z.  If Ritz vectors are
 *          desired, then  LDZ >= max( 1, N ).  In any case,  LDZ >= 1.
 *
 *  SIGMAR  Double precision  (INPUT)
 *          If IPARAM(7) = 3 or 4, represents the real part of the shift.
 *          Not referenced if IPARAM(7) = 1 or 2.
 *
 *  SIGMAI  Double precision  (INPUT)
 *          If IPARAM(7) = 3 or 4, represents the imaginary part of the
 *          shift.
 *          Not referenced if IPARAM(7) = 1 or 2.
 *
 *  WORKEV  Double precision work array of dimension 3*NCV.  (WORKSPACE)
 *
 *  N       Integer.  (INPUT)
 *          Dimension of the eigenproblem.
 *
 *  AB      Double precision array of dimension LDA by N. (INPUT)
 *          The matrix A in band storage, in rows KL+1 to
 *          2*KL+KU+1; rows 1 to KL of the array need not be set.
 *          The j-th column of A is stored in the j-th column of the
 *          array AB as follows:
 *          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
 *
 *  MB      Double precision array of dimension LDA by N. (INPUT)
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
 *  RFAC    Double precision array of LDA by N. (WORKSPACE/OUTPUT)
 *          RFAC is used to store the LU factors of MB when IPARAM(7) = 2
 *          is invoked.  It is used to store the LU factors of
 *          (A-sigma*M) when IPARAM(7) = 3 is invoked with a real shift.
 *          It is not referenced when IPARAM(7) = 1 or 4.
 *
 *  CFAC    Complex*16 array of LDA by N. (WORKSPACE/OUTPUT)
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
 *  TOL     Double precision scalar.  (INPUT)
 *          Stopping criteria: the relative accuracy of the Ritz value
 *          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I)).
 *          If TOL .LE. 0. is passed a default is set:
 *          DEFAULT = DLAMCH('EPS')  (machine precision as computed
 *                    by the LAPACK auxiliary subroutine DLAMCH).
 *
 *  RESID   Double precision array of length N.  (INPUT/OUTPUT)
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
 *          by dnaupd for OP.
 *
 *  V       Double precision array N by NCV+1.  (OUTPUT)
 *          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
 *                       represent approximate Schur vectors that span the
 *                       desired invariant subspace.
 *          NOTE: The array Z may be set equal to first NEV+1 columns of the
 *          Arnoldi basis vector array V computed by DNAUPD. In this case
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
 *          Must be 1,2,3,4; See under \Description of dnband for the
 *          four modes available.
 *
 *          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
 *          OUTPUT: NUMOP  = total number of OP*z operations,
 *                  NUMOPB = total number of B*z operations if BMAT='G',
 *                  NUMREO = total number of steps of re-orthogonalization.
 *
 * WORKD    Double precision work array of length at least 3*n. (WORKSPACE)
 *
 * WORKL    Double precision work array of length LWORKL. (WORKSPACE)
 *
 * LWORKL   Integer.  (INPUT)
 *          LWORKL must be at least 3*NCV**2 + 6*NCV.
 *
 * WORKC    Complex*16 array of length N. (WORKSPACE)
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
 *          =  1: The Schur form computed by LAPACK routine dlahqr
 *                could not be reordered by LAPACK routine dtrsen.
 *                Re-enter subroutine DNEUPD with IPARAM(5)=NCV and
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
 *                Informational error from LAPACK routine dlahqr.
 *          = -9: Error return from calculation of eigenvectors.
 *                Informational error from LAPACK routine dtrevc.
 *          = -10: IPARAM(7) must be 1,2,3,4.
 *          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
 *          = -12: HOWMNY = 'S' not yet implemented
 *          = -13: HOWMNY must be one of 'A' or 'P'
 *          = -14: DNAUPD did not find any eigenvalues to sufficient
 *                 accuracy.
 *          = -15: Overflow occurs when we try to transform the Ritz
 *                 values returned from DNAUPD to those of the original
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
 * Routines called:
 *     dnaupd  ARPACK reverse communication interface routine.
 *     dneupd  ARPACK routine that returns Ritz values and (optionally)
 *             Ritz vectors.
 *     dgbtrf  LAPACK band matrix factorization routine.
 *     dgbtrs  LAPACK band linear system solve routine.
 *     zgbtrf  LAPACK complex band matrix factorization routine.
 *     zgbtrs  LAPACK complex linear system solve routine.
 *     dlacpy  LAPACK matrix copy routine.
 *     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
 *     dlamch  LAPACK routine to compute the underflow threshold.
 *     dcopy   Level 1 BLAS that copies one vector to another.
 *     ddot    Level 1 BLAS that computes the dot product of two vectors.
 *     dnrm2   Level 1 BLAS that computes the norm of a vector.
 *     dgbmv   Level 2 BLAS that computes the band matrix vector product.
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
int dnband_(a_bool *rvec, char *howmny, a_bool *select, double *dr, double *di, double *z, a_int *ldz, double *sigmar, double *sigmai, double *workev, a_int *n, double *ab, double *mb, a_int *lda, double *rfac, a_dcomplex *cfac, a_int *kl, a_int *ku, char *which, char *bmat, a_int *nev, double *tol, double *resid, a_int *ncv, double *v, a_int *ldv, a_int *iparam, double *workd, double *workl, a_int *lworkl, a_dcomplex *workc, a_int *iwork, a_int *info)
{
    /* System generated locals */
    a_int v_dim1, v_offset, z_dim1, z_offset, ab_dim1, ab_offset, mb_dim1, mb_offset, rfac_dim1, rfac_offset, cfac_dim1, cfac_offset, i__1, i__2, i__3, i__4, i__5;
    double d__1, d__2, d__3;
    a_dcomplex z__1, z__2;

    /* Local variables */
    a_int i, j, ido;
    double deni;
    a_int imid;
    double denr;
    a_int ibot, itop, type, ierr = 0;
    double numr;
    double dmdul;
    a_bool first;
    a_int ipntr[14];
    double safmin;

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

    safmin = dlamch_("S");

    /* -------------------------------------------------------------- */
    /* Set type of the problem to be solved. Check consistency        */
    /* between BMAT and IPARAM(7).                                    */
    /* type = 1 --> Solving standard problem in regular mode.         */
    /* type = 2 --> Solving standard problem in shift-invert mode.    */
    /* type = 3 --> Solving generalized problem in regular mode.      */
    /* type = 4 --> Solving generalized problem in shift-invert mode. */
    /* type = 5 --> Solving standard problem in shift-invert mode     */
    /*              using iparam(7) = 4 in DNAUPD.                    */
    /* type = 6 --> Solving generalized problem in shift-invert mode. */
    /*              using iparam(7) = 4 in DNAUPD.                    */
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
        if (*sigmai == 0.)
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

        if (*sigmai == 0.)
        {

            /* --------------------------------- */
            /* Construct (A-sigmar*I) and factor */
            /* in real arithmetic.               */
            /* --------------------------------- */

            dlacpy_("A", &ibot, n, &ab[ab_offset], lda, &rfac[rfac_offset], lda);
            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                rfac[imid + j * rfac_dim1] = ab[imid + j * ab_dim1] - *sigmar;
            }
            dgbtrf_(n, n, kl, ku, &rfac[rfac_offset], lda, &iwork[1], &ierr);
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
                    z__1.r = ab[i__4], z__1.i = 0.;
                    cfac[i__3].r = z__1.r, cfac[i__3].i = z__1.i;
                }
            }

            i__1 = *n;
            for (j = 1; j <= i__1; ++j)
            {
                i__2 = imid + j * cfac_dim1;
                i__3 = imid + j * cfac_dim1;
                z__2.r = *sigmar, z__2.i = *sigmai;
                z__1.r = cfac[i__3].r - z__2.r, z__1.i = cfac[i__3].i - z__2.i;
                cfac[i__2].r = z__1.r, cfac[i__2].i = z__1.i;
            }

            zgbtrf_(n, n, kl, ku, &cfac[cfac_offset], lda, &iwork[1], &ierr);
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
        /* routine dgbtrf to factor M.                   */
        /* --------------------------------------------- */

        dlacpy_("A", &ibot, n, &mb[mb_offset], lda, &rfac[rfac_offset], lda);
        dgbtrf_(n, n, kl, ku, &rfac[rfac_offset], lda, &iwork[1], &ierr);
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

        if (*sigmai == 0.)
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

            dgbtrf_(n, n, kl, ku, &rfac[rfac_offset], lda, &iwork[1], &ierr);
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
                    d__1 = ab[i + j * ab_dim1] - *sigmar * mb[i + j * mb_dim1];
                    d__2 = -(*sigmai) * mb[i + j * mb_dim1];
                    z__1.r = d__1, z__1.i = d__2;
                    cfac[i__3].r = z__1.r, cfac[i__3].i = z__1.i;
                }
            }

            zgbtrf_(n, n, kl, ku, &cfac[cfac_offset], lda, &iwork[1], &ierr);
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

    dnaupd_(&ido, bmat, n, which, nev, tol, &resid[1], ncv, &v[v_offset], ldv, &iparam[1], ipntr, &workd[1], &workl[1], lworkl, info);

    if (ido == -1)
    {

        if (type == 1)
        {

            /* -------------------------- */
            /* Perform  y <--- OP*x = A*x */
            /* -------------------------- */

            dgbmv_("N", n, n, kl, ku, &one, &ab[itop + ab_dim1], lda, &workd[ipntr[0]], &i_one, &zero, &workd[ipntr[1]], &i_one);
        }
        else if (type == 2)
        {

            if (*sigmai == 0.)
            {

                /* -------------------------------- */
                /* Shift is real.  Perform          */
                /* y <--- OP*x = inv[A-sigmar*I]*x  */
                /* to force the starting vector     */
                /* into the range of OP.            */
                /* -------------------------------- */

                dcopy_(n, &workd[ipntr[0]], &i_one, &workd[ipntr[1]], &i_one);
                dgbtrs_("N", n, kl, ku, &i_one, &rfac[rfac_offset], lda, &iwork[1], &workd[ipntr[1]], n, &ierr);
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
                    z__1.r = workd[i__3], z__1.i = 0.;
                    workc[i__2].r = z__1.r, workc[i__2].i = z__1.i;
                }

                zgbtrs_("N", n, kl, ku, &i_one, &cfac[cfac_offset], lda, &iwork[1], &workc[1], n, &ierr);
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

            dgbmv_("N", n, n, kl, ku, &one, &ab[itop + ab_dim1], lda, &workd[ipntr[0]], &i_one, &zero, &workd[ipntr[1]], &i_one);

            dgbtrs_("N", n, kl, ku, &i_one, &rfac[rfac_offset], lda, &iwork[1], &workd[ipntr[1]], n, &ierr);
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

            dgbmv_("N", n, n, kl, ku, &one, &mb[itop + mb_dim1], lda, &workd[ipntr[0]], &i_one, &zero, &workd[ipntr[1]], &i_one);

            if (*sigmai == 0.)
            {

                /* ------------------- */
                /* Shift is real, stay */
                /* in real arithmetic. */
                /* ------------------- */

                dgbtrs_("N", n, kl, ku, &i_one, &rfac[rfac_offset], lda, &iwork[1], &workd[ipntr[1]], n, &ierr);
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
                    z__1.r = workd[i__3], z__1.i = 0.;
                    workc[i__2].r = z__1.r, workc[i__2].i = z__1.i;
                }

                zgbtrs_("N", n, kl, ku, &i_one, &cfac[cfac_offset], lda, &iwork[1], &workc[1], n, &ierr);
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
                z__1.r = workd[i__3], z__1.i = 0.;
                workc[i__2].r = z__1.r, workc[i__2].i = z__1.i;
            }

            zgbtrs_("N", n, kl, ku, &i_one, &cfac[cfac_offset], lda, &iwork[1], &workc[1], n, &ierr);
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

            dgbmv_("N", n, n, kl, ku, &one, &mb[itop + mb_dim1], lda, &workd[ipntr[0]], &i_one, &zero, &workd[ipntr[1]], &i_one);

            i__1 = *n;
            for (i = 1; i <= i__1; ++i)
            {
                i__2 = i;
                i__3 = ipntr[1] + i - 1;
                z__1.r = workd[i__3], z__1.i = 0.;
                workc[i__2].r = z__1.r, workc[i__2].i = z__1.i;
            }

            zgbtrs_("N", n, kl, ku, &i_one, &cfac[cfac_offset], lda, &iwork[1], &workc[1], n, &ierr);
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

            dgbmv_("N", n, n, kl, ku, &one, &ab[itop + ab_dim1], lda, &workd[ipntr[0]], &i_one, &zero, &workd[ipntr[1]], &i_one);
        }
        else if (type == 2)
        {

            if (*sigmai == 0.)
            {

                /* -------------------------------- */
                /* Shift is real.  Perform          */
                /* y <--- OP*x = inv[A-sigmar*I]*x. */
                /* -------------------------------- */

                dcopy_(n, &workd[ipntr[0]], &i_one, &workd[ipntr[1]], &i_one);
                dgbtrs_("N", n, kl, ku, &i_one, &rfac[rfac_offset], lda, &iwork[1], &workd[ipntr[1]], n, &ierr);
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
                    z__1.r = workd[i__3], z__1.i = 0.;
                    workc[i__2].r = z__1.r, workc[i__2].i = z__1.i;
                }

                zgbtrs_("N", n, kl, ku, &i_one, &cfac[cfac_offset], lda, &iwork[1], &workc[1], n, &ierr);
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

            dgbmv_("N", n, n, kl, ku, &one, &ab[itop + ab_dim1], lda, &workd[ipntr[0]], &i_one, &zero, &workd[ipntr[1]], &i_one);

            dgbtrs_("N", n, kl, ku, &i_one, &rfac[rfac_offset], lda, &iwork[1], &workd[ipntr[1]], n, &ierr);
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

            if (*sigmai == 0.)
            {

                /* ---------------------- */
                /* Shift is real, stay in */
                /* real arithmetic.       */
                /* ---------------------- */

                dcopy_(n, &workd[ipntr[2]], &i_one, &workd[ipntr[1]], &i_one);
                dgbtrs_("N", n, kl, ku, &i_one, &rfac[rfac_offset], lda, &iwork[1], &workd[ipntr[1]], n, &ierr);
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
                    z__1.r = workd[i__3], z__1.i = 0.;
                    workc[i__2].r = z__1.r, workc[i__2].i = z__1.i;
                }

                zgbtrs_("N", n, kl, ku, &i_one, &cfac[cfac_offset], lda, &iwork[1], &workc[1], n, &ierr);
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
                z__1.r = workd[i__3], z__1.i = 0.;
                workc[i__2].r = z__1.r, workc[i__2].i = z__1.i;
            }

            zgbtrs_("N", n, kl, ku, &i_one, &cfac[cfac_offset], lda, &iwork[1], &workc[1], n, &ierr);
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
                z__1.r = workd[i__3], z__1.i = 0.;
                workc[i__2].r = z__1.r, workc[i__2].i = z__1.i;
            }

            zgbtrs_("N", n, kl, ku, &i_one, &cfac[cfac_offset], lda, &iwork[1], &workc[1], n, &ierr);
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

        dgbmv_("N", n, n, kl, ku, &one, &mb[itop + mb_dim1], lda, &workd[ipntr[0]], &i_one, &zero, &workd[ipntr[1]], &i_one);
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
            /* documentation in DNAUPD  */
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

                dneupd_(rvec, "A", &select[1], &dr[1], &di[1], &z[z_offset], ldz, sigmar, sigmai, &workev[1], bmat, n, which, nev, tol, &resid[1], ncv, &v[v_offset], ldv, &iparam[1], ipntr, &workd[1], &workl[1], lworkl, info);

                if (*info != 0)
                {

                    /* ---------------------------------- */
                    /* Check the documentation of DNEUPD. */
                    /* ---------------------------------- */

                    printf(" \n");
                    printf(" Error with _neupd = %d", (*info));
                    printf(" Check the documentation of _neupd \n");
                    printf(" \n");
                    return ierr;
                }
                else if (*sigmai != 0.)
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

                            if (di[j] == 0.)
                            {

                                /* ------------------------------------ */
                                /* Eigenvalue is real. Compute          */
                                /* d = (x'*inv[A-sigma*M]*M*x) / (x'*x) */
                                /* ------------------------------------ */

                                dgbmv_("Nontranspose", n, n, kl, ku, &one, &mb[itop + mb_dim1], lda, &z[j * z_dim1 + 1], &i_one, &zero, &workd[1], &i_one);
                                i__2 = *n;
                                for (i = 1; i <= i__2; ++i)
                                {
                                    i__3 = i;
                                    i__4 = i;
                                    z__1.r = workd[i__4], z__1.i = 0.;
                                    workc[i__3].r = z__1.r, workc[i__3].i = z__1.i;
                                }
                                zgbtrs_("N", n, kl, ku, &i_one, &cfac[cfac_offset], lda, &iwork[1], &workc[1], n, info);
                                i__2 = *n;
                                for (i = 1; i <= i__2; ++i)
                                {
                                    i__3 = i;
                                    workd[i] = workc[i__3].r;
                                    workd[i + *n] = workc[i].i;
                                }
                                denr = ddot_(n, &z[j * z_dim1 + 1], &i_one, &workd[1], &i_one);
                                deni = ddot_(n, &z[j * z_dim1 + 1], &i_one, &workd[*n + 1], &i_one);
                                /* Computing 2nd power */
                                d__1 = dnrm2_(n, &z[j * z_dim1 + 1], &i_one);
                                numr = d__1 * d__1;
                                /* Computing 2nd power */
                                d__1 = dlapy2_(&denr, &deni);
                                dmdul = d__1 * d__1;
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

                                dgbmv_("Nontranspose", n, n, kl, ku, &one, &mb[itop + mb_dim1], lda, &z[j * z_dim1 + 1], &i_one, &zero, &workd[1], &i_one);
                                dgbmv_("Nontranspose", n, n, kl, ku, &one, &mb[itop + mb_dim1], lda, &z[(j + 1) * z_dim1 + 1], &i_one, &zero, &workd[*n + 1], &i_one);
                                i__2 = *n;
                                for (i = 1; i <= i__2; ++i)
                                {
                                    i__3 = i;
                                    i__4 = i;
                                    i__5 = i + *n;
                                    z__1.r = workd[i__4], z__1.i = workd[i__5];
                                    workc[i__3].r = z__1.r, workc[i__3].i = z__1.i;
                                }

                                /* -------------------------- */
                                /* Compute inv(A-sigma*M)*M*x */
                                /* -------------------------- */

                                zgbtrs_("N", n, kl, ku, &i_one, &cfac[cfac_offset], lda, &iwork[1], &workc[1], n, info);

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
                                denr = ddot_(n, &z[j * z_dim1 + 1], &i_one, &workd[1], &i_one);
                                denr += ddot_(n, &z[(j + 1) * z_dim1 + 1], &i_one, &workd[*n + 1], &i_one);
                                deni = ddot_(n, &z[j * z_dim1 + 1], &i_one, &workd[*n + 1], &i_one);
                                deni -= ddot_(n, &z[(j + 1) * z_dim1 + 1], &i_one, &workd[1], &i_one);

                                /* -------------- */
                                /* Compute (x'*x) */
                                /* -------------- */

                                d__2 = dnrm2_(n, &z[j * z_dim1 + 1], &i_one);
                                d__3 = dnrm2_(n, &z[(j + 1) * z_dim1 + 1], &i_one);
                                /* Computing 2nd power */
                                d__1 = dlapy2_(&d__2, &d__3);
                                numr = d__1 * d__1;

                                /* -------------------------------------- */
                                /* Compute (x'x) / (x'*inv(A-sigma*M)*Mx) */
                                /* -------------------------------------- */

                                /* Computing 2nd power */
                                d__1 = dlapy2_(&denr, &deni);
                                dmdul = d__1 * d__1;
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

                            if (di[j] == 0.)
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
                                    z__1.r = z[i__4], z__1.i = 0.;
                                    workc[i__3].r = z__1.r, workc[i__3].i = z__1.i;
                                }
                                zgbtrs_("N", n, kl, ku, &i_one, &cfac[cfac_offset], lda, &iwork[1], &workc[1], n, info);
                                i__2 = *n;
                                for (i = 1; i <= i__2; ++i)
                                {
                                    i__3 = i;
                                    workd[i] = workc[i__3].r;
                                    workd[i + *n] = workc[i].i;
                                }
                                denr = ddot_(n, &z[j * z_dim1 + 1], &i_one, &workd[1], &i_one);
                                deni = ddot_(n, &z[j * z_dim1 + 1], &i_one, &workd[*n + 1], &i_one);
                                /* Computing 2nd power */
                                d__1 = dnrm2_(n, &z[j * z_dim1 + 1], &i_one);
                                numr = d__1 * d__1;
                                /* Computing 2nd power */
                                d__1 = dlapy2_(&denr, &deni);
                                dmdul = d__1 * d__1;
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
                                    z__1.r = z[i__4], z__1.i = z[i__5];
                                    workc[i__3].r = z__1.r, workc[i__3].i = z__1.i;
                                }

                                /* ------------------------- */
                                /* Compute inv[A-sigma*I]*x. */
                                /* ------------------------- */

                                zgbtrs_("N", n, kl, ku, &i_one, &cfac[cfac_offset], lda, &iwork[1], &workc[1], n, info);

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
                                denr = ddot_(n, &z[j * z_dim1 + 1], &i_one, &workd[1], &i_one);
                                denr += ddot_(n, &z[(j + 1) * z_dim1 + 1], &i_one, &workd[*n + 1], &i_one);
                                deni = ddot_(n, &z[j * z_dim1 + 1], &i_one, &workd[*n + 1], &i_one);
                                deni -= ddot_(n, &z[(j + 1) * z_dim1 + 1], &i_one, &workd[1], &i_one);

                                /* -------------- */
                                /* Compute (x'*x) */
                                /* -------------- */

                                d__2 = dnrm2_(n, &z[j * z_dim1 + 1], &i_one);
                                d__3 = dnrm2_(n, &z[(j + 1) * z_dim1 + 1], &i_one);
                                /* Computing 2nd power */
                                d__1 = dlapy2_(&d__2, &d__3);
                                numr = d__1 * d__1;

                                /* -------------------------------------- */
                                /* Compute (x'x) / (x'*inv(A-sigma*I)*x). */
                                /* -------------------------------------- */

                                /* Computing 2nd power */
                                d__1 = dlapy2_(&denr, &deni);
                                dmdul = d__1 * d__1;
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
    /* L O O P  B A C K to call DNAUPD again. */
    /* -------------------------------------- */

    goto L90;

    return 0;
} /* dnband_ */
