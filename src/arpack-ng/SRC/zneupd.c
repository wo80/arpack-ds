/* SRC\zneupd.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

/* Table of constant values */

static a_dcomplex z_one = {1., 0.};
static a_dcomplex z_zero = {0., 0.};
static a_int i_one = 1;
static a_bool b_true = TRUE_;
/**
 * \BeginDoc
 *
 * \Name: zneupd
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
 *  There is negligible additional cost to obtain eigenvectors.  An orthonormal
 *  basis is always computed.  There is an additional storage cost of n*nev
 *  if both are requested (in this case a separate array Z must be supplied).
 *
 *  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
 *  are derived from approximate eigenvalues and eigenvectors of
 *  of the linear operator OP prescribed by the MODE selection in the
 *  call to ZNAUPD.  ZNAUPD must be called before this routine is called.
 *  These approximate eigenvalues and vectors are commonly called Ritz
 *  values and Ritz vectors respectively.  They are referred to as such
 *  in the comments that follow.   The computed orthonormal basis for the
 *  invariant subspace corresponding to these Ritz values is referred to as a
 *  Schur basis.
 *
 *  The definition of OP as well as other terms and the relation of computed
 *  Ritz values and vectors of OP with respect to the given problem
 *  A*z = lambda*B*z may be found in the header of ZNAUPD.  For a brief
 *  description, see definitions of IPARAM(7), MODE and WHICH in the
 *  documentation of ZNAUPD.
 *
 * \Usage:
 *  call zneupd
 *     ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, WORKEV, BMAT,
 *       N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD,
 *       WORKL, LWORKL, RWORK, INFO )
 *
 * \Arguments:
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
 *          computed. To select the  Ritz vector corresponding to a
 *          Ritz value D(j), SELECT(j) must be set to .TRUE..
 *          If HOWMNY = 'A' or 'P', SELECT need not be initialized
 *          but it is used as internal workspace.
 *
 *  D       Complex*16 array of dimension NEV+1.  (OUTPUT)
 *          On exit, D contains the  Ritz  approximations
 *          to the eigenvalues lambda for A*z = lambda*B*z.
 *
 *  Z       Complex*16 N by NEV array    (OUTPUT)
 *          On exit, if RVEC = .TRUE. and HOWMNY = 'A', then the columns of
 *          Z represents approximate eigenvectors (Ritz vectors) corresponding
 *          to the NCONV=IPARAM(5) Ritz values for eigensystem
 *          A*z = lambda*B*z.
 *
 *          If RVEC = .FALSE. or HOWMNY = 'P', then Z is NOT REFERENCED.
 *
 *          NOTE: If if RVEC = .TRUE. and a Schur basis is not required,
 *          the array Z may be set equal to first NEV+1 columns of the Arnoldi
 *          basis array V computed by ZNAUPD.  In this case the Arnoldi basis
 *          will be destroyed and overwritten with the eigenvector basis.
 *
 *  LDZ     Integer.  (INPUT)
 *          The leading dimension of the array Z.  If Ritz vectors are
 *          desired, then  LDZ .ge.  max( 1, N ) is required.
 *          In any case,  LDZ .ge. 1 is required.
 *
 *  SIGMA   Complex*16  (INPUT)
 *          If IPARAM(7) = 3 then SIGMA represents the shift.
 *          Not referenced if IPARAM(7) = 1 or 2.
 *
 *  WORKEV  Complex*16 work array of dimension 2*NCV.  (WORKSPACE)
 *
 *  **** The remaining arguments MUST be the same as for the   ****
 *  **** call to ZNAUPD that was just completed.               ****
 *
 *  NOTE: The remaining arguments
 *
 *           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR,
 *           WORKD, WORKL, LWORKL, RWORK, INFO
 *
 *         must be passed directly to ZNEUPD following the last call
 *         to ZNAUPD.  These arguments MUST NOT BE MODIFIED between
 *         the the last call to ZNAUPD and the call to ZNEUPD.
 *
 *  Three of these parameters (V, WORKL and INFO) are also output parameters:
 *
 *  V       Complex*16 N by NCV array.  (INPUT/OUTPUT)
 *
 *          Upon INPUT: the NCV columns of V contain the Arnoldi basis
 *                      vectors for OP as constructed by ZNAUPD .
 *
 *          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
 *                       contain approximate Schur vectors that span the
 *                       desired invariant subspace.
 *
 *          NOTE: If the array Z has been set equal to first NEV+1 columns
 *          of the array V and RVEC=.TRUE. and HOWMNY= 'A', then the
 *          Arnoldi basis held by V has been overwritten by the desired
 *          Ritz vectors.  If a separate array Z has been passed then
 *          the first NCONV=IPARAM(5) columns of V will contain approximate
 *          Schur vectors that span the desired invariant subspace.
 *
 *  WORKL   Double precision work array of length LWORKL.  (OUTPUT/WORKSPACE)
 *          WORKL(1:ncv*ncv+2*ncv) contains information obtained in
 *          znaupd.  They are not changed by zneupd.
 *          WORKL(ncv*ncv+2*ncv+1:3*ncv*ncv+4*ncv) holds the
 *          untransformed Ritz values, the untransformed error estimates of
 *          the Ritz values, the upper triangular matrix for H, and the
 *          associated matrix representation of the invariant subspace for H.
 *
 *          Note: IPNTR(9:13) contains the pointer into WORKL for addresses
 *          of the above information computed by zneupd.
 *          -------------------------------------------------------------
 *          IPNTR(9):  pointer to the NCV RITZ values of the
 *                     original system.
 *          IPNTR(10): Not used
 *          IPNTR(11): pointer to the NCV corresponding error estimates.
 *          IPNTR(12): pointer to the NCV by NCV upper triangular
 *                     Schur matrix for H.
 *          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
 *                     of the upper Hessenberg matrix H. Only referenced by
 *                     zneupd if RVEC = .TRUE. See Remark 2 below.
 *          -------------------------------------------------------------
 *
 *  INFO    Integer.  (OUTPUT)
 *          Error flag on output.
 *          =  0: Normal exit.
 *
 *          =  1: The Schur form computed by LAPACK routine csheqr
 *                could not be reordered by LAPACK routine ztrsen.
 *                Re-enter subroutine zneupd with IPARAM(5)=NCV and
 *                increase the size of the array D to have
 *                dimension at least dimension NCV and allocate at least NCV
 *                columns for Z. NOTE: Not necessary if Z and V share
 *                the same space. Please notify the authors if this error
 *                occurs.
 *
 *          = -1: N must be positive.
 *          = -2: NEV must be positive.
 *          = -3: NCV-NEV >= 1 and less than or equal to N.
 *          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
 *          = -6: BMAT must be one of 'I' or 'G'.
 *          = -7: Length of private work WORKL array is not sufficient.
 *          = -8: Error return from LAPACK eigenvalue calculation.
 *                This should never happened.
 *          = -9: Error return from calculation of eigenvectors.
 *                Informational error from LAPACK routine ztrevc.
 *          = -10: IPARAM(7) must be 1,2,3
 *          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
 *          = -12: HOWMNY = 'S' not yet implemented
 *          = -13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.
 *          = -14: ZNAUPD did not find any eigenvalues to sufficient
 *                 accuracy.
 *          = -15: ZNEUPD got a different count of the number of converged
 *                 Ritz values than ZNAUPD got.  This indicates the user
 *                 probably made an error in passing data from ZNAUPD to
 *                 ZNEUPD or that the data was modified before entering
 *                 ZNEUPD
 *
 * \BeginLib
 *
 * \References:
 *  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
 *     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
 *     pp 357-385.
 *  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
 *     Restarted Arnoldi Iteration", Rice University Technical Report
 *     TR95-13, Department of Computational and Applied Mathematics.
 *  3. B. Nour-Omid, B. N. Parlett, T. Ericsson and P. S. Jensen,
 *     "How to Implement the Spectral Transformation", Math Comp.,
 *     Vol. 48, No. 178, April, 1987 pp. 664-673.
 *
 * \Routines called:
 *     ivout   ARPACK utility routine that prints integers.
 *     zmout   ARPACK utility routine that prints matrices
 *     zvout   ARPACK utility routine that prints vectors.
 *     zgeqr2  LAPACK routine that computes the QR factorization of
 *             a matrix.
 *     zlacpy  LAPACK matrix copy routine.
 *     zlahqr  LAPACK routine that computes the Schur form of a
 *             upper Hessenberg matrix.
 *     zlaset  LAPACK matrix initialization routine.
 *     ztrevc  LAPACK routine to compute the eigenvectors of a matrix
 *             in upper triangular form.
 *     ztrsen  LAPACK routine that re-orders the Schur form.
 *     zunm2r  LAPACK routine that applies an orthogonal matrix in
 *             factored form.
 *     dlamch  LAPACK routine that determines machine constants.
 *     ztrmm   Level 3 BLAS matrix times an upper triangular matrix.
 *     zgeru   Level 2 BLAS rank one update to a matrix.
 *     zcopy   Level 1 BLAS that copies one vector to another .
 *     zscal   Level 1 BLAS that scales a vector.
 *     zdscal  Level 1 BLAS that scales a complex vector by a real number.
 *     dznrm2  Level 1 BLAS that computes the norm of a complex vector.
 *
 * \Remarks
 *
 *  1. Currently only HOWMNY = 'A' and 'P' are implemented.
 *
 *  2. Schur vectors are an orthogonal representation for the basis of
 *     Ritz vectors. Thus, their numerical properties are often superior.
 *     If RVEC = .true. then the relationship
 *             A * V(:,1:IPARAM(5)) = V(:,1:IPARAM(5)) * T, and
 *       transpose( V(:,1:IPARAM(5)) ) * V(:,1:IPARAM(5)) = I
 *     are approximately satisfied.
 *     Here T is the leading submatrix of order IPARAM(5) of the
 *     upper triangular matrix stored workl(ipntr(12)).
 *
 * \Authors
 *     Danny Sorensen               Phuong Vu
 *     Richard Lehoucq              CRPC / Rice University
 *     Chao Yang                    Houston, Texas
 *     Dept. of Computational &
 *     Applied Mathematics
 *     Rice University
 *     Houston, Texas
 *
 * \SCCS Information: @(#)
 * FILE: neupd.F   SID: 2.8   DATE OF SID: 07/21/02   RELEASE: 2
 *
 * \EndLib
 */
int zneupd_(a_bool *rvec, const char *howmny, a_bool *select, a_dcomplex *d, a_dcomplex *z,
     a_int *ldz, a_dcomplex *sigma, a_dcomplex *workev, const char *bmat, a_int *n, const char *which,
     a_int *nev, double *tol, a_dcomplex *resid, a_int *ncv, a_dcomplex *v, a_int *ldv,
     a_int *iparam, a_int *ipntr, a_dcomplex *workd, a_dcomplex *workl, a_int *lworkl,
     double *rwork, a_int *info)
{
    /* System generated locals */
    a_int v_dim1, v_offset, z_dim1, z_offset, i__1, i__2;
    double d__1, d__2, d__3, d__4;
    a_dcomplex z__1, z__2;

    /* Local variables */
    a_int j, k, ih, jj, iq, np;
    a_dcomplex vl[1];
    a_int wr, ibd, ldh, ldq;
    double sep;
    a_int irz, mode;
    double eps23;
    a_int ierr;
    a_dcomplex temp;
    a_int iwev;
    char type[7];
    a_int ritz, iheig, ihbds;
    double conds;
    a_bool reord;
    a_int nconv;
    double rtemp;
    a_dcomplex rnorm;
    a_int nconv2;
    a_int bounds, invsub, iuptri, msglvl, outncv, numcnv, ishift;

    /* ---------------------- */
    /* Set default parameters */
    /* ---------------------- */

    /* Parameter adjustments */
    --workd;
    --resid;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z -= z_offset;
    --d;
    --rwork;
    --workev;
    --select;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --iparam;
    --ipntr;
    --workl;

    msglvl = debug_1.mceupd;
    mode = iparam[7];
    nconv = iparam[5];
    *info = 0;

    /* ------------------------------- */
    /* Get machine dependent constant. */
    /* ------------------------------- */

    eps23 = dlamch_("E");
    eps23 = pow(eps23, TWO_THIRDS);

    /* ----------------------------- */
    /* Quick return                  */
    /* Check for incompatible input  */
    /* ----------------------------- */

    ierr = 0;

    if (nconv <= 0)
    {
        ierr = -14;
    }
    else if (*n <= 0)
    {
        ierr = -1;
    }
    else if (*nev <= 0)
    {
        ierr = -2;
    }
    else if (*ncv <= *nev + 1 || *ncv > *n)
    {
        ierr = -3;
    }
    else if (strcmp(which, "LM") != 0 && strcmp(which, "SM") != 0 && strcmp(which, "LR") != 0 &&
             strcmp(which, "SR") != 0 && strcmp(which, "LI") != 0 && strcmp(which, "SI") != 0)
    {
        ierr = -5;
    }
    else if (*bmat != 'I' && *bmat != 'G')
    {
        ierr = -6;
    }
    else /* if(complicated condition) */
    {
        /* Computing 2nd power */
        i__1 = *ncv;
        if (*lworkl < i__1 * i__1 * 3 + (*ncv << 2))
        {
            ierr = -7;
        }
        else if (*howmny != 'A' && *howmny != 'P' && *howmny != 'S' && *rvec)
        {
            ierr = -13;
        }
        else if (*howmny == 'S')
        {
            ierr = -12;
        }
    }

    if (mode == 1 || mode == 2)
    {
        strcpy(type, "REGULR");
    }
    else if (mode == 3)
    {
        strcpy(type, "SHIFTI");
    }
    else
    {
        ierr = -10;
    }
    if (mode == 1 && *bmat == 'G')
    {
        ierr = -11;
    }

    /* ---------- */
    /* Error Exit */
    /* ---------- */

    if (ierr != 0)
    {
        *info = ierr;
        goto L9000;
    }

    /* ------------------------------------------------------ */
    /* Pointer into WORKL for address of H, RITZ, WORKEV, Q   */
    /* etc... and the remaining workspace.                    */
    /* Also update pointer to be used on output.              */
    /* Memory is laid out as follows:                         */
    /* workl(1:ncv*ncv) := generated Hessenberg matrix        */
    /* workl(ncv*ncv+1:ncv*ncv+ncv) := ritz values            */
    /* workl(ncv*ncv+ncv+1:ncv*ncv+2*ncv) := error bounds     */
    /* ------------------------------------------------------ */

    /* --------------------------------------------------------- */
    /* The following is used and set by ZNEUPD.                 */
    /* workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := The untransformed */
    /*                                      Ritz values.         */
    /* workl(ncv*ncv+3*ncv+1:ncv*ncv+4*ncv) := The untransformed */
    /*                                      error bounds of      */
    /*                                      the Ritz values      */
    /* workl(ncv*ncv+4*ncv+1:2*ncv*ncv+4*ncv) := Holds the upper */
    /*                                      triangular matrix    */
    /*                                      for H.               */
    /* workl(2*ncv*ncv+4*ncv+1: 3*ncv*ncv+4*ncv) := Holds the    */
    /*                                      associated matrix    */
    /*                                      representation of    */
    /*                                      the invariant        */
    /*                                      subspace for H.      */
    /* GRAND total of NCV * ( 3 * NCV + 4 ) locations.           */
    /* --------------------------------------------------------- */

    ih = ipntr[5];
    ritz = ipntr[6];
    iq = ipntr[7];
    bounds = ipntr[8];
    ldh = *ncv;
    ldq = *ncv;
    iheig = bounds + ldh;
    ihbds = iheig + ldh;
    iuptri = ihbds + ldh;
    invsub = iuptri + ldh * *ncv;
    ipntr[9] = iheig;
    ipntr[11] = ihbds;
    ipntr[12] = iuptri;
    ipntr[13] = invsub;
    wr = 1;
    iwev = wr + *ncv;

    /* --------------------------------------- */
    /* irz points to the Ritz values computed  */
    /*     by _neigh before exiting _naup2.    */
    /* ibd points to the Ritz estimates        */
    /*     computed by _neigh before exiting   */
    /*     _naup2.                             */
    /* --------------------------------------- */

    irz = ipntr[14] + *ncv * *ncv;
    ibd = irz + *ncv;

    /* ---------------------------------- */
    /* RNORM is B-norm of the RESID(1:N). */
    /* ---------------------------------- */

    i__1 = ih + 2;
    rnorm.r = workl[i__1].r, rnorm.i = workl[i__1].i;
    i__1 = ih + 2;
    workl[i__1].r = 0., workl[i__1].i = 0.;

#ifndef NO_TRACE
    if (msglvl > 2)
    {
        zvout_(*ncv, &workl[irz], debug_1.ndigit, "_neupd: Ritz values passed in from _NAUPD.");
        zvout_(*ncv, &workl[ibd], debug_1.ndigit, "_neupd: Ritz estimates passed in from _NAUPD.");
    }
#endif

    if (*rvec)
    {

        reord = FALSE_;

        /* ------------------------------------------------- */
        /* Use the temporary bounds array to store indices   */
        /* These will be used to mark the select array later */
        /* ------------------------------------------------- */

        i__1 = *ncv;
        for (j = 1; j <= i__1; ++j)
        {
            i__2 = bounds + j - 1;
            workl[i__2].r = (double)j, workl[i__2].i = 0.;
            select[j] = FALSE_;
            /* L10: */
        }

        /* ----------------------------------- */
        /* Select the wanted Ritz values.      */
        /* Sort the Ritz values so that the    */
        /* wanted ones appear at the tailing   */
        /* NEV positions of workl(irr) and     */
        /* workl(iri).  Move the corresponding */
        /* error estimates in workl(ibd)       */
        /* accordingly.                        */
        /* ----------------------------------- */

        np = *ncv - *nev;
        ishift = 0;
        zngets_(&ishift, which, nev, &np, &workl[irz], &workl[bounds]);

#ifndef NO_TRACE
        if (msglvl > 2)
        {
            zvout_(*ncv, &workl[irz], debug_1.ndigit, "_neupd: Ritz values after calling _NGETS.");
            zvout_(*ncv, &workl[bounds], debug_1.ndigit, "_neupd: Ritz value indices after calling _NGETS.");
        }
#endif

        /* --------------------------------------------------- */
        /* Record indices of the converged wanted Ritz values  */
        /* Mark the select array for possible reordering       */
        /* --------------------------------------------------- */

        numcnv = 0;
        i__1 = *ncv;
        for (j = 1; j <= i__1; ++j)
        {
            /* Computing MAX */
            i__2 = irz + *ncv - j;
            d__3 = workl[i__2].r;
            d__4 = workl[irz + *ncv - j].i;
            d__1 = eps23, d__2 = dlapy2_(&d__3, &d__4);
            rtemp = max(d__1, d__2);
            i__2 = bounds + *ncv - j;
            jj = (a_int)workl[i__2].r;
            i__2 = ibd + jj - 1;
            d__1 = workl[i__2].r;
            d__2 = workl[ibd + jj - 1].i;
            if (numcnv < nconv && dlapy2_(&d__1, &d__2) <= *tol * rtemp)
            {
                select[jj] = TRUE_;
                ++numcnv;
                if (jj > nconv)
                {
                    reord = TRUE_;
                }
            }
            /* L11: */
        }

        /* --------------------------------------------------------- */
        /* Check the count (numcnv) of converged Ritz values with    */
        /* the number (nconv) reported by dnaupd.  If these two      */
        /* are different then there has probably been an error       */
        /* caused by incorrect passing of the dnaupd data.           */
        /* --------------------------------------------------------- */

#ifndef NO_TRACE
        if (msglvl > 2)
        {
            ivout_(1, &numcnv, debug_1.ndigit, "_neupd: Number of specified eigenvalues");
            ivout_(1, &nconv, debug_1.ndigit, "_neupd: Number of \"converged\" eigenvalues");
        }
#endif

        if (numcnv != nconv)
        {
            *info = -15;
            goto L9000;
        }

        /* ----------------------------------------------------- */
        /* Call LAPACK routine zlahqr to compute the Schur form */
        /* of the upper Hessenberg matrix returned by ZNAUPD.   */
        /* Make a copy of the upper Hessenberg matrix.           */
        /* Initialize the Schur vector matrix Q to the identity. */
        /* ----------------------------------------------------- */

        i__1 = ldh * *ncv;
        zcopy_(&i__1, &workl[ih], &i_one, &workl[iuptri], &i_one);
        zlaset_("A", ncv, ncv, &z_zero, &z_one, &workl[invsub], &ldq);
        zlahqr_(&b_true, &b_true, ncv, &i_one, ncv, &workl[iuptri], &ldh, &workl[iheig], &i_one, ncv, &workl[invsub], &ldq, &ierr);
        zcopy_(ncv, &workl[invsub + *ncv - 1], &ldq, &workl[ihbds], &i_one);

        if (ierr != 0)
        {
            *info = -8;
            goto L9000;
        }

#ifndef NO_TRACE
        if (msglvl > 1)
        {
            zvout_(*ncv, &workl[iheig], debug_1.ndigit, "_neupd: Eigenvalues of H");
            zvout_(*ncv, &workl[ihbds], debug_1.ndigit, "_neupd: Last row of the Schur vector matrix");
            if (msglvl > 3)
            {
                zmout_(*ncv, *ncv, &workl[iuptri], ldh, debug_1.ndigit, "_neupd: The upper triangular matrix ");
            }
        }
#endif

        if (reord)
        {

            /* --------------------------------------------- */
            /* Reorder the computed upper triangular matrix. */
            /* --------------------------------------------- */

            ztrsen_("N", "V", &select[1], ncv, &workl[iuptri], &ldh, &workl[invsub], &ldq, &workl[iheig], &nconv2, &conds, &sep, &workev[1], ncv, &ierr);

            if (nconv2 < nconv)
            {
                nconv = nconv2;
            }
            if (ierr == 1)
            {
                *info = 1;
                goto L9000;
            }

#ifndef NO_TRACE
            if (msglvl > 2)
            {
                zvout_(*ncv, &workl[iheig], debug_1.ndigit, "_neupd: Eigenvalues of H--reordered");
                if (msglvl > 3)
                {
                    zmout_(*ncv, *ncv, &workl[iuptri], ldq, debug_1.ndigit, "_neupd: Triangular matrix after re-ordering");
                }
            }
#endif

        }

        /* ------------------------------------------- */
        /* Copy the last row of the Schur basis matrix */
        /* to workl(ihbds).  This vector will be used  */
        /* to compute the Ritz estimates of converged  */
        /* Ritz values.                                */
        /* ------------------------------------------- */

        zcopy_(ncv, &workl[invsub + *ncv - 1], &ldq, &workl[ihbds], &i_one);

        /* ------------------------------------------ */
        /* Place the computed eigenvalues of H into D */
        /* if a spectral transformation was not used. */
        /* ------------------------------------------ */

        if (strcmp(type, "REGULR") == 0)
        {
            zcopy_(&nconv, &workl[iheig], &i_one, &d[1], &i_one);
        }

        /* -------------------------------------------------------- */
        /* Compute the QR factorization of the matrix representing  */
        /* the wanted invariant subspace located in the first NCONV */
        /* columns of workl(invsub,ldq).                            */
        /* -------------------------------------------------------- */

        zgeqr2_(ncv, &nconv, &workl[invsub], &ldq, &workev[1], &workev[*ncv + 1], &ierr);

        /* ------------------------------------------------------ */
        /* * Postmultiply V by Q using zunm2r.                    */
        /* * Copy the first NCONV columns of VQ into Z.           */
        /* * Postmultiply Z by R.                                 */
        /* The N by NCONV matrix Z is now a matrix representation */
        /* of the approximate invariant subspace associated with  */
        /* the Ritz values in workl(iheig). The first NCONV       */
        /* columns of V are now approximate Schur vectors         */
        /* associated with the upper triangular matrix of order   */
        /* NCONV in workl(iuptri).                                */
        /* ------------------------------------------------------ */

        zunm2r_("R", "N", n, ncv, &nconv, &workl[invsub], &ldq, &workev[1], &v[v_offset], ldv, &workd[*n + 1], &ierr);
        zlacpy_("A", n, &nconv, &v[v_offset], ldv, &z[z_offset], ldz);

        i__1 = nconv;
        for (j = 1; j <= i__1; ++j)
        {

            /* ------------------------------------------------- */
            /* Perform both a column and row scaling if the      */
            /* diagonal element of workl(invsub,ldq) is negative */
            /* I'm lazy and don't take advantage of the upper    */
            /* triangular form of workl(iuptri,ldq).             */
            /* Note that since Q is orthogonal, R is a diagonal  */
            /* matrix consisting of plus or minus ones.          */
            /* ------------------------------------------------- */

            i__2 = invsub + (j - 1) * ldq + j - 1;
            if (workl[i__2].r < 0.)
            {
                z__1.r = -1., z__1.i = -0.;
                zscal_(&nconv, &z__1, &workl[iuptri + j - 1], &ldq);
                z__1.r = -1., z__1.i = -0.;
                zscal_(&nconv, &z__1, &workl[iuptri + (j - 1) * ldq], &i_one);
            }

            /* L20: */
        }

        if (*howmny == 'A')
        {

            /* ------------------------------------------ */
            /* Compute the NCONV wanted eigenvectors of T */
            /* located in workl(iuptri,ldq).              */
            /* ------------------------------------------ */

            i__1 = *ncv;
            for (j = 1; j <= i__1; ++j)
            {
                if (j <= nconv)
                {
                    select[j] = TRUE_;
                }
                else
                {
                    select[j] = FALSE_;
                }
                /* L30: */
            }

            ztrevc_("R", "S", &select[1], ncv, &workl[iuptri], &ldq, vl, &i_one, &workl[invsub], &ldq, ncv, &outncv, &workev[1], &rwork[1], &ierr);

            if (ierr != 0)
            {
                *info = -9;
                goto L9000;
            }

            /* ---------------------------------------------- */
            /* Scale the returning eigenvectors so that their */
            /* Euclidean norms are all one. LAPACK subroutine */
            /* ztrevc returns each eigenvector normalized so  */
            /* that the element of largest magnitude has      */
            /* magnitude 1.                                   */
            /* ---------------------------------------------- */

            i__1 = nconv;
            for (j = 1; j <= i__1; ++j)
            {
                rtemp = dznrm2_(ncv, &workl[invsub + (j - 1) * ldq], &i_one);
                rtemp = 1. / rtemp;
                zdscal_(ncv, &rtemp, &workl[invsub + (j - 1) * ldq], &i_one);

                /* ---------------------------------------- */
                /* Ritz estimates can be obtained by taking */
                /* the inner product of the last row of the */
                /* Schur basis of H with eigenvectors of T. */
                /* Note that the eigenvector matrix of T is */
                /* upper triangular, thus the length of the */
                /* inner product can be set to j.           */
                /* ---------------------------------------- */

                i__2 = j;
                zzdotc_(&z__1, &j, &workl[ihbds], &i_one, &workl[invsub + (j - 1) * ldq], &i_one);
                workev[i__2].r = z__1.r, workev[i__2].i = z__1.i;
                /* L40: */
            }

#ifndef NO_TRACE
            if (msglvl > 2)
            {
                zcopy_(&nconv, &workl[invsub + *ncv - 1], &ldq, &workl[ihbds], &i_one);
                zvout_(nconv, &workl[ihbds], debug_1.ndigit, "_neupd: Last row of the eigenvector matrix for T");
                if (msglvl > 3)
                {
                    zmout_(*ncv, *ncv, &workl[invsub], ldq, debug_1.ndigit, "_neupd: The eigenvector matrix for T");
                }
            }
#endif

            /* ------------------------------------- */
            /* Copy Ritz estimates into workl(ihbds) */
            /* ------------------------------------- */

            zcopy_(&nconv, &workev[1], &i_one, &workl[ihbds], &i_one);

            /* -------------------------------------------- */
            /* The eigenvector matrix Q of T is triangular. */
            /* Form Z*Q.                                    */
            /* -------------------------------------------- */

            ztrmm_("R", "U", "N", "N", n, &nconv, &z_one, &workl[invsub], &ldq, &z[z_offset], ldz);
        }
    }
    else
    {

        /* ------------------------------------------------ */
        /* An approximate invariant subspace is not needed. */
        /* Place the Ritz values computed ZNAUPD into D.    */
        /* ------------------------------------------------ */

        zcopy_(&nconv, &workl[ritz], &i_one, &d[1], &i_one);
        zcopy_(&nconv, &workl[ritz], &i_one, &workl[iheig], &i_one);
        zcopy_(&nconv, &workl[bounds], &i_one, &workl[ihbds], &i_one);
    }

    /* ---------------------------------------------- */
    /* Transform the Ritz values and possibly vectors */
    /* and corresponding error bounds of OP to those  */
    /* of A*x = lambda*B*x.                           */
    /* ---------------------------------------------- */

    if (strcmp(type, "REGULR") == 0)
    {

        if (*rvec)
        {
            zscal_(ncv, &rnorm, &workl[ihbds], &i_one);
        }
    }
    else
    {

        /* ------------------------------------- */
        /*   A spectral transformation was used. */
        /* * Determine the Ritz estimates of the */
        /*   Ritz values in the original system. */
        /* ------------------------------------- */

        if (*rvec)
        {
            zscal_(ncv, &rnorm, &workl[ihbds], &i_one);
        }

        i__1 = *ncv;
        for (k = 1; k <= i__1; ++k)
        {
            i__2 = iheig + k - 1;
            temp.r = workl[i__2].r, temp.i = workl[i__2].i;
            i__2 = ihbds + k - 1;
            ar_z_div(&z__2, &workl[ihbds + k - 1], &temp);
            ar_z_div(&z__1, &z__2, &temp);
            workl[i__2].r = z__1.r, workl[i__2].i = z__1.i;
            /* L50: */
        }
    }

    /* --------------------------------------------------------- */
    /* *  Transform the Ritz values back to the original system. */
    /*    For TYPE = 'SHIFTI' the transformation is              */
    /*             lambda = 1/theta + sigma                      */
    /* NOTES:                                                    */
    /* *The Ritz vectors are not affected by the transformation. */
    /* --------------------------------------------------------- */

    if (strcmp(type, "SHIFTI") == 0)
    {
        i__1 = nconv;
        for (k = 1; k <= i__1; ++k)
        {
            i__2 = k;
            ar_z_div(&z__2, &z_one, &workl[iheig + k - 1]);
            z__1.r = z__2.r + sigma->r, z__1.i = z__2.i + sigma->i;
            d[i__2].r = z__1.r, d[i__2].i = z__1.i;
            /* L60: */
        }
    }

    if (strcmp(type, "REGULR") != 0 && msglvl > 1)
    {
        zvout_(nconv, &d[1], debug_1.ndigit, "_neupd: Untransformed Ritz values.");
        zvout_(nconv, &workl[ihbds], debug_1.ndigit, "_neupd: Ritz estimates of the untransformed Ritz values.");
    }
    else if (msglvl > 1)
    {
        zvout_(nconv, &d[1], debug_1.ndigit, "_neupd: Converged Ritz values.");
        zvout_(nconv, &workl[ihbds], debug_1.ndigit, "_neupd: Associated Ritz estimates.");
    }

    /* ----------------------------------------------- */
    /* Eigenvector Purification step. Formally perform */
    /* one of inverse subspace iteration. Only used    */
    /* for MODE = 3. See reference 3.                  */
    /* ----------------------------------------------- */

    if (*rvec && *howmny == 'A' && strcmp(type, "SHIFTI") == 0)
    {

        /* ---------------------------------------------- */
        /* Purify the computed Ritz vectors by adding a   */
        /* little bit of the residual vector:             */
        /*                      T                         */
        /*          resid(:)*( e    s ) / theta           */
        /*                      NCV                       */
        /* where H s = s theta.                           */
        /* ---------------------------------------------- */

        i__1 = nconv;
        for (j = 1; j <= i__1; ++j)
        {
            i__2 = iheig + j - 1;
            if (workl[i__2].r != 0. || workl[i__2].i != 0.)
            {
                i__2 = j;
                ar_z_div(&z__1, &workl[invsub + (j - 1) * ldq + *ncv - 1], &workl[iheig + j - 1]);
                workev[i__2].r = z__1.r, workev[i__2].i = z__1.i;
            }
            /* L100: */
        }
        /* ------------------------------------- */
        /* Perform a rank one update to Z and    */
        /* purify all the Ritz vectors together. */
        /* ------------------------------------- */

        zgeru_(n, &nconv, &z_one, &resid[1], &i_one, &workev[1], &i_one, &z[z_offset], ldz);
    }

L9000:

    return 0;

    /* ------------- */
    /* End of zneupd*/
    /* ------------- */

} /* zneupd_ */
