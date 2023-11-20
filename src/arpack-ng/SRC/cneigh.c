/* SRC\cneigh.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

/* Table of constant values */

static a_fcomplex c_one = {1.f, 0.f};
static a_fcomplex c_zero = {0.f, 0.f};
static a_bool b_true = TRUE_;
static a_int i_one = 1;
/**
 * \BeginDoc
 *
 * \Name: cneigh
 *
 * \Description:
 *  Compute the eigenvalues of the current upper Hessenberg matrix
 *  and the corresponding Ritz estimates given the current residual norm.
 *
 * \Usage:
 *  call cneigh
 *     ( RNORM, N, H, LDH, RITZ, BOUNDS, Q, LDQ, WORKL, RWORK, IERR )
 *
 * \Arguments
 *  RNORM   Real scalar.  (INPUT)
 *          Residual norm corresponding to the current upper Hessenberg
 *          matrix H.
 *
 *  N       Integer.  (INPUT)
 *          Size of the matrix H.
 *
 *  H       Complex N by N array.  (INPUT)
 *          H contains the current upper Hessenberg matrix.
 *
 *  LDH     Integer.  (INPUT)
 *          Leading dimension of H exactly as declared in the calling
 *          program.
 *
 *  RITZ    Complex array of length N.  (OUTPUT)
 *          On output, RITZ(1:N) contains the eigenvalues of H.
 *
 *  BOUNDS  Complex array of length N.  (OUTPUT)
 *          On output, BOUNDS contains the Ritz estimates associated with
 *          the eigenvalues held in RITZ.  This is equal to RNORM
 *          times the last components of the eigenvectors corresponding
 *          to the eigenvalues in RITZ.
 *
 *  Q       Complex N by N array.  (WORKSPACE)
 *          Workspace needed to store the eigenvectors of H.
 *
 *  LDQ     Integer.  (INPUT)
 *          Leading dimension of Q exactly as declared in the calling
 *          program.
 *
 *  WORKL   Complex work array of length N**2 + 3*N.  (WORKSPACE)
 *          Private (replicated) array on each PE or array allocated on
 *          the front end.  This is needed to keep the full Schur form
 *          of H and also in the calculation of the eigenvectors of H.
 *
 *  RWORK   Real  work array of length N (WORKSPACE)
 *          Private (replicated) array on each PE or array allocated on
 *          the front end.
 *
 *  IERR    Integer.  (OUTPUT)
 *          Error exit flag from clahqr or ctrevc.
 *
 * \EndDoc
 *
 * -----------------------------------------------------------------------
 *
 * \BeginLib
 *
 * \Local variables:
 *     xxxxxx  Complex
 *
 * \Routines called:
 *     ivout   ARPACK utility routine that prints integers.
 *     arscnd  ARPACK utility routine for timing.
 *     cmout   ARPACK utility routine that prints matrices
 *     cvout   ARPACK utility routine that prints vectors.
 *     svout   ARPACK utility routine that prints vectors.
 *     clacpy  LAPACK matrix copy routine.
 *     clahqr  LAPACK routine to compute the Schur form of an
 *             upper Hessenberg matrix.
 *     claset  LAPACK matrix initialization routine.
 *     ctrevc  LAPACK routine to compute the eigenvectors of a matrix
 *             in upper triangular form
 *     ccopy   Level 1 BLAS that copies one vector to another.
 *     csscal  Level 1 BLAS that scales a complex vector by a real number.
 *     scnrm2  Level 1 BLAS that computes the norm of a vector.
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
 * FILE: neigh.F   SID: 2.2   DATE OF SID: 4/20/96   RELEASE: 2
 *
 *
 * \EndLib
 */
int cneigh_(float *rnorm, a_int *n, a_fcomplex *h, a_int *ldh, a_fcomplex *ritz, a_fcomplex *bounds, a_fcomplex *q, a_int *ldq, a_fcomplex *workl, float *rwork, a_int *ierr)
{
    /* System generated locals */
    a_int h_dim1, h_offset, q_dim1, q_offset, i__1;
    float r__1;

    /* Local variables */
    a_int j;
    static float t0, t1;
    a_fcomplex vl[1];
    float temp;
    a_bool select[1];
    a_int msglvl;

    /* ----------------------------- */
    /* Initialize timing statistics  */
    /* & message level for debugging */
    /* ----------------------------- */

    /* Parameter adjustments */
    --rwork;
    --workl;
    --bounds;
    --ritz;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h -= h_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;

    arscnd_(&t0);
    msglvl = debug_1.mceigh;

    if (msglvl > 2)
    {
        cmout_(n, n, &h[h_offset], ldh, debug_1.ndigit, "_neigh: Entering upper Hessenberg matrix H ");
    }

    /* -------------------------------------------------------- */
    /* 1. Compute the eigenvalues, the last components of the   */
    /*    corresponding Schur vectors and the full Schur form T */
    /*    of the current upper Hessenberg matrix H.             */
    /*    clahqr returns the full Schur form of H               */
    /*    in WORKL(1:N**2), and the Schur vectors in q.         */
    /* -------------------------------------------------------- */

    clacpy_("A", n, n, &h[h_offset], ldh, &workl[1], n);
    claset_("A", n, n, &c_zero, &c_one, &q[q_offset], ldq);
    clahqr_(&b_true, &b_true, n, &i_one, n, &workl[1], ldh, &ritz[1], &i_one, n, &q[q_offset], ldq, ierr);
    if (*ierr != 0)
    {
        goto L9000;
    }

    ccopy_(n, &q[*n - 1 + q_dim1], ldq, &bounds[1], &i_one);
    if (msglvl > 1)
    {
        cvout_(n, &bounds[1], debug_1.ndigit, "_neigh: last row of the Schur matrix for H");
    }

    /* -------------------------------------------------------- */
    /* 2. Compute the eigenvectors of the full Schur form T and */
    /*    apply the Schur vectors to get the corresponding      */
    /*    eigenvectors.                                         */
    /* -------------------------------------------------------- */

    ctrevc_("R", "B", select, n, &workl[1], n, vl, n, &q[q_offset], ldq, n, n, &workl[*n * *n + 1], &rwork[1], ierr);

    if (*ierr != 0)
    {
        goto L9000;
    }

    /* ---------------------------------------------- */
    /* Scale the returning eigenvectors so that their */
    /* Euclidean norms are all one. LAPACK subroutine */
    /* ctrevc returns each eigenvector normalized so  */
    /* that the element of largest magnitude has      */
    /* magnitude 1; here the magnitude of a complex   */
    /* number (x,y) is taken to be |x| + |y|.         */
    /* ---------------------------------------------- */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j)
    {
        temp = scnrm2_(n, &q[j * q_dim1 + 1], &i_one);
        r__1 = 1.f / temp;
        csscal_(n, &r__1, &q[j * q_dim1 + 1], &i_one);
        /* L10: */
    }

    if (msglvl > 1)
    {
        ccopy_(n, &q[*n + q_dim1], ldq, &workl[1], &i_one);
        cvout_(n, &workl[1], debug_1.ndigit, "_neigh: Last row of the eigenvector matrix for H");
    }

    /* -------------------------- */
    /* Compute the Ritz estimates */
    /* -------------------------- */

    ccopy_(n, &q[*n + q_dim1], n, &bounds[1], &i_one);
    csscal_(n, rnorm, &bounds[1], &i_one);

    if (msglvl > 2)
    {
        cvout_(n, &ritz[1], debug_1.ndigit, "_neigh: The eigenvalues of H");
        cvout_(n, &bounds[1], debug_1.ndigit, "_neigh: Ritz estimates for the eigenvalues of H");
    }

    arscnd_(&t1);
    timing_1.tceigh += t1 - t0;

L9000:
    return 0;

    /* ------------- */
    /* End of cneigh */
    /* ------------- */

} /* cneigh_ */
