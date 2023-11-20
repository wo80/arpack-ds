/* SRC\dseigt.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

/* Table of constant values */

static a_int i_one = 1;

/**
 * \BeginDoc
 *
 * \Name: dseigt
 *
 * \Description:
 *  Compute the eigenvalues of the current symmetric tridiagonal matrix
 *  and the corresponding error bounds given the current residual norm.
 *
 * \Usage:
 *  call dseigt
 *     ( RNORM, N, H, LDH, EIG, BOUNDS, WORKL, IERR )
 *
 * \Arguments
 *  RNORM   Double precision scalar.  (INPUT)
 *          RNORM contains the residual norm corresponding to the current
 *          symmetric tridiagonal matrix H.
 *
 *  N       Integer.  (INPUT)
 *          Size of the symmetric tridiagonal matrix H.
 *
 *  H       Double precision N by 2 array.  (INPUT)
 *          H contains the symmetric tridiagonal matrix with the
 *          subdiagonal in the first column starting at H(2,1) and the
 *          main diagonal in second column.
 *
 *  LDH     Integer.  (INPUT)
 *          Leading dimension of H exactly as declared in the calling
 *          program.
 *
 *  EIG     Double precision array of length N.  (OUTPUT)
 *          On output, EIG contains the N eigenvalues of H possibly
 *          unsorted.  The BOUNDS arrays are returned in the
 *          same sorted order as EIG.
 *
 *  BOUNDS  Double precision array of length N.  (OUTPUT)
 *          On output, BOUNDS contains the error estimates corresponding
 *          to the eigenvalues EIG.  This is equal to RNORM times the
 *          last components of the eigenvectors corresponding to the
 *          eigenvalues in EIG.
 *
 *  WORKL   Double precision work array of length 3*N.  (WORKSPACE)
 *          Private (replicated) array on each PE or array allocated on
 *          the front end.
 *
 *  IERR    Integer.  (OUTPUT)
 *          Error exit flag from dstqrb.
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
 * \Routines called:
 *     dstqrb  ARPACK routine that computes the eigenvalues and the
 *             last components of the eigenvectors of a symmetric
 *             and tridiagonal matrix.
 *     arscnd  ARPACK utility routine for timing.
 *     dvout   ARPACK utility routine that prints vectors.
 *     dcopy   Level 1 BLAS that copies one vector to another.
 *
 * \Author
 *     Danny Sorensen               Phuong Vu
 *     Richard Lehoucq              CRPC / Rice University
 *     Dept. of Computational &     Houston, Texas
 *     Applied Mathematics
 *     Rice University
 *     Houston, Texas
 *
 * \Revision history:
 *     xx/xx/92: Version ' 2.4'
 *
 * \SCCS Information: @(#)
 * FILE: seigt.F   SID: 2.4   DATE OF SID: 8/27/96   RELEASE: 2
 *
 *
 * \EndLib
 */
int dseigt_(double *rnorm, a_int *n, double *h, a_int *ldh, double *eig, double *bounds, double *workl, a_int *ierr)
{
    /* System generated locals */
    a_int h_dim1, h_offset, i__1;
    double d__1;

    /* Local variables */
    a_int k;
    static float t0, t1;
    a_int msglvl;

    /* ----------------------------- */
    /* Initialize timing statistics  */
    /* & message level for debugging */
    /* ----------------------------- */

    /* Parameter adjustments */
    --workl;
    --bounds;
    --eig;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h -= h_offset;

#ifndef NO_TIMER
    arscnd_(&t0);
#endif

    msglvl = debug_1.mseigt;

#ifndef NO_TRACE
    if (msglvl > 0)
    {
        dvout_(n, &h[(h_dim1 << 1) + 1], debug_1.ndigit, "_seigt: main diagonal of matrix H");
        if (*n > 1)
        {
            i__1 = *n - 1;
            dvout_(i__1, &h[h_dim1 + 2], debug_1.ndigit, "_seigt: sub diagonal of matrix H");
        }
    }
#endif

    dcopy_(n, &h[(h_dim1 << 1) + 1], &i_one, &eig[1], &i_one);
    i__1 = *n - 1;
    dcopy_(&i__1, &h[h_dim1 + 2], &i_one, &workl[1], &i_one);
    dstqrb_(n, &eig[1], &workl[1], &bounds[1], &workl[*n + 1], ierr);
    if (*ierr != 0)
    {
        goto L9000;
    }
#ifndef NO_TRACE
    if (msglvl > 1)
    {
        dvout_(n, &bounds[1], debug_1.ndigit, "_seigt: last row of the eigenvector matrix for H");
    }
#endif

    /* --------------------------------------------- */
    /* Finally determine the error bounds associated */
    /* with the n Ritz values of H.                  */
    /* --------------------------------------------- */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k)
    {
        bounds[k] = *rnorm * (d__1 = bounds[k], abs(d__1));
        /* L30: */
    }

#ifndef NO_TIMER
    arscnd_(&t1);
    timing_1.tseigt += t1 - t0;
#endif

L9000:
    return 0;

    /* ------------- */
    /* End of dseigt */
    /* ------------- */

} /* dseigt_ */
