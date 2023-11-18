/* SRC\dneigh.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

/* Table of constant values */

static a_bool b_true = TRUE_;
static a_int i_one = 1;
static double d_one = 1.;
static double d_zero = 0.;

/* ----------------------------------------------------------------------- */
/* \BeginDoc */

/* \Name: dneigh */

/* \Description: */
/*  Compute the eigenvalues of the current upper Hessenberg matrix */
/*  and the corresponding Ritz estimates given the current residual norm. */

/* \Usage: */
/*  call dneigh */
/*     ( RNORM, N, H, LDH, RITZR, RITZI, BOUNDS, Q, LDQ, WORKL, IERR ) */

/* \Arguments */
/*  RNORM   Double precision scalar.  (INPUT) */
/*          Residual norm corresponding to the current upper Hessenberg */
/*          matrix H. */

/*  N       Integer.  (INPUT) */
/*          Size of the matrix H. */

/*  H       Double precision N by N array.  (INPUT) */
/*          H contains the current upper Hessenberg matrix. */

/*  LDH     Integer.  (INPUT) */
/*          Leading dimension of H exactly as declared in the calling */
/*          program. */

/*  RITZR,  Double precision arrays of length N.  (OUTPUT) */
/*  RITZI   On output, RITZR(1:N) (resp. RITZI(1:N)) contains the real */
/*          (respectively imaginary) parts of the eigenvalues of H. */

/*  BOUNDS  Double precision array of length N.  (OUTPUT) */
/*          On output, BOUNDS contains the Ritz estimates associated with */
/*          the eigenvalues RITZR and RITZI.  This is equal to RNORM */
/*          times the last components of the eigenvectors corresponding */
/*          to the eigenvalues in RITZR and RITZI. */

/*  Q       Double precision N by N array.  (WORKSPACE) */
/*          Workspace needed to store the eigenvectors of H. */

/*  LDQ     Integer.  (INPUT) */
/*          Leading dimension of Q exactly as declared in the calling */
/*          program. */

/*  WORKL   Double precision work array of length N**2 + 3*N.  (WORKSPACE) */
/*          Private (replicated) array on each PE or array allocated on */
/*          the front end.  This is needed to keep the full Schur form */
/*          of H and also in the calculation of the eigenvectors of H. */

/*  IERR    Integer.  (OUTPUT) */
/*          Error exit flag from dlahqr or dtrevc. */

/* \EndDoc */

/* ----------------------------------------------------------------------- */

/* \BeginLib */

/* \Local variables: */
/*     xxxxxx  real */

/* \Routines called: */
/*     dlahqr  LAPACK routine to compute the real Schur form of an */
/*             upper Hessenberg matrix and last row of the Schur vectors. */
/*     arscnd  ARPACK utility routine for timing. */
/*     dmout   ARPACK utility routine that prints matrices */
/*     dvout   ARPACK utility routine that prints vectors. */
/*     dlacpy  LAPACK matrix copy routine. */
/*     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully. */
/*     dtrevc  LAPACK routine to compute the eigenvectors of a matrix */
/*             in upper quasi-triangular form */
/*     dgemv   Level 2 BLAS routine for matrix vector multiplication. */
/*     dcopy   Level 1 BLAS that copies one vector to another . */
/*     dnrm2   Level 1 BLAS that computes the norm of a vector. */
/*     dscal   Level 1 BLAS that scales a vector. */

/* \Author */
/*     Danny Sorensen               Phuong Vu */
/*     Richard Lehoucq              CRPC / Rice University */
/*     Dept. of Computational &     Houston, Texas */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \Revision history: */
/*     xx/xx/92: Version ' 2.1' */

/* \SCCS Information: @(#) */
/* FILE: neigh.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2 */

/* \Remarks */
/*     None */

/* \EndLib */

/* ----------------------------------------------------------------------- */

int dneigh_(double *rnorm, a_int *n, double *h__, a_int *ldh, double *ritzr, double *ritzi, double *bounds, double *q, a_int *ldq, double *workl, a_int *ierr)
{
    /* System generated locals */
    a_int h_dim1, h_offset, q_dim1, q_offset, i__1;
    double d__1, d__2;

    /* Local variables */
    a_int i__, j;
    static float t0, t1;
    double vl[1], temp;
    extern double dnrm2_(a_int *, double *, a_int *);
    extern int dscal_(a_int *, double *, double *, a_int *);
    a_int iconj;
    extern int dgemv_(char *, a_int *, a_int *, double *, double *, a_int *, double *, a_int *, double *, double *, a_int *, ftnlen), dmout_(a_int *, a_int *, a_int *, double *, a_int *, a_int *, char *, ftnlen), dvout_(a_int *, a_int *, double *, a_int *, char *, ftnlen);
    extern double dlapy2_(double *, double *);
    extern int arscnd_(float *), dlahqr_(a_bool *, a_bool *, a_int *, a_int *, a_int *, double *, a_int *, double *, double *, a_int *, a_int *, double *, a_int *, a_int *);
    a_bool select[1];
    a_int msglvl;
    extern int dlacpy_(char *, a_int *, a_int *, double *, a_int *, double *, a_int *, ftnlen), dtrevc_(char *, char *, a_bool *, a_int *, double *, a_int *, double *, a_int *, double *, a_int *, a_int *, a_int *, double *, a_int *, ftnlen, ftnlen);

    /*     %----------------------------------------------------% */
    /*     | Include files for debugging and timing information | */
    /*     %----------------------------------------------------% */

    /* \SCCS Information: @(#) */
    /* FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 */

    /*     %---------------------------------% */
    /*     | See debug.doc for documentation | */
    /*     %---------------------------------% */

    /*     %------------------% */
    /*     | Scalar Arguments | */
    /*     %------------------% */

    /*     %--------------------------------% */
    /*     | See stat.doc for documentation | */
    /*     %--------------------------------% */

    /* \SCCS Information: @(#) */
    /* FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2 */

    /*     %-----------------% */
    /*     | Array Arguments | */
    /*     %-----------------% */

    /*     %------------% */
    /*     | Parameters | */
    /*     %------------% */

    /*     %------------------------% */
    /*     | Local Scalars & Arrays | */
    /*     %------------------------% */

    /*     %----------------------% */
    /*     | External Subroutines | */
    /*     %----------------------% */

    /*     %--------------------% */
    /*     | External Functions | */
    /*     %--------------------% */

    /*     %---------------------% */
    /*     | Intrinsic Functions | */
    /*     %---------------------% */

    /*     %-----------------------% */
    /*     | Executable Statements | */
    /*     %-----------------------% */

    /*     %-------------------------------% */
    /*     | Initialize timing statistics  | */
    /*     | & message level for debugging | */
    /*     %-------------------------------% */

    /* Parameter adjustments */
    --workl;
    --bounds;
    --ritzi;
    --ritzr;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;

    /* Function Body */
    arscnd_(&t0);
    msglvl = debug_1.mneigh;

    if (msglvl > 2)
    {
        dmout_(&debug_1.logfil, n, n, &h__[h_offset], ldh, &debug_1.ndigit, "_neigh: Entering upper Hessenberg matrix H ", (ftnlen)43);
    }

    /*     %-----------------------------------------------------------% */
    /*     | 1. Compute the eigenvalues, the last components of the    | */
    /*     |    corresponding Schur vectors and the full Schur form T  | */
    /*     |    of the current upper Hessenberg matrix H.              | */
    /*     | dlahqr returns the full Schur form of H in WORKL(1:N**2)  | */
    /*     | and the last components of the Schur vectors in BOUNDS.   | */
    /*     %-----------------------------------------------------------% */

    dlacpy_("All", n, n, &h__[h_offset], ldh, &workl[1], n, (ftnlen)3);
    i__1 = *n - 1;
    for (j = 1; j <= i__1; ++j)
    {
        bounds[j] = 0.;
        /* L5: */
    }
    bounds[*n] = 1.;
    dlahqr_(&b_true, &b_true, n, &i_one, n, &workl[1], n, &ritzr[1], &ritzi[1], &i_one, &i_one, &bounds[1], &i_one, ierr);
    if (*ierr != 0)
    {
        goto L9000;
    }

    if (msglvl > 1)
    {
        dvout_(&debug_1.logfil, n, &bounds[1], &debug_1.ndigit,
               "_neigh: las"
               "t row of the Schur matrix for H",
               (ftnlen)42);
    }

    /*     %-----------------------------------------------------------% */
    /*     | 2. Compute the eigenvectors of the full Schur form T and  | */
    /*     |    apply the last components of the Schur vectors to get  | */
    /*     |    the last components of the corresponding eigenvectors. | */
    /*     | Remember that if the i-th and (i+1)-st eigenvalues are    | */
    /*     | complex conjugate pairs, then the real & imaginary part   | */
    /*     | of the eigenvector components are split across adjacent   | */
    /*     | columns of Q.                                             | */
    /*     %-----------------------------------------------------------% */

    dtrevc_("R", "A", select, n, &workl[1], n, vl, n, &q[q_offset], ldq, n, n, &workl[*n * *n + 1], ierr, (ftnlen)1, (ftnlen)1);

    if (*ierr != 0)
    {
        goto L9000;
    }

    /*     %------------------------------------------------% */
    /*     | Scale the returning eigenvectors so that their | */
    /*     | euclidean norms are all one. LAPACK subroutine | */
    /*     | dtrevc returns each eigenvector normalized so  | */
    /*     | that the element of largest magnitude has      | */
    /*     | magnitude 1; here the magnitude of a complex   | */
    /*     | number (x,y) is taken to be |x| + |y|.         | */
    /*     %------------------------------------------------% */

    iconj = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
        if ((d__1 = ritzi[i__], abs(d__1)) <= 0.)
        {

            /*           %----------------------% */
            /*           | Real eigenvalue case | */
            /*           %----------------------% */

            temp = dnrm2_(n, &q[i__ * q_dim1 + 1], &i_one);
            d__1 = 1. / temp;
            dscal_(n, &d__1, &q[i__ * q_dim1 + 1], &i_one);
        }
        else
        {

            /*           %-------------------------------------------% */
            /*           | Complex conjugate pair case. Note that    | */
            /*           | since the real and imaginary part of      | */
            /*           | the eigenvector are stored in consecutive | */
            /*           | columns, we further normalize by the      | */
            /*           | square root of two.                       | */
            /*           %-------------------------------------------% */

            if (iconj == 0)
            {
                d__1 = dnrm2_(n, &q[i__ * q_dim1 + 1], &i_one);
                d__2 = dnrm2_(n, &q[(i__ + 1) * q_dim1 + 1], &i_one);
                temp = dlapy2_(&d__1, &d__2);
                d__1 = 1. / temp;
                dscal_(n, &d__1, &q[i__ * q_dim1 + 1], &i_one);
                d__1 = 1. / temp;
                dscal_(n, &d__1, &q[(i__ + 1) * q_dim1 + 1], &i_one);
                iconj = 1;
            }
            else
            {
                iconj = 0;
            }
        }
        /* L10: */
    }

    dgemv_("T", n, n, &d_one, &q[q_offset], ldq, &bounds[1], &i_one, &d_zero, &workl[1], &i_one, (ftnlen)1);

    if (msglvl > 1)
    {
        dvout_(&debug_1.logfil, n, &workl[1], &debug_1.ndigit,
               "_neigh: Last"
               " row of the eigenvector matrix for H",
               (ftnlen)48);
    }

    /*     %----------------------------% */
    /*     | Compute the Ritz estimates | */
    /*     %----------------------------% */

    iconj = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
        if ((d__1 = ritzi[i__], abs(d__1)) <= 0.)
        {

            /*           %----------------------% */
            /*           | Real eigenvalue case | */
            /*           %----------------------% */

            bounds[i__] = *rnorm * (d__1 = workl[i__], abs(d__1));
        }
        else
        {

            /*           %-------------------------------------------% */
            /*           | Complex conjugate pair case. Note that    | */
            /*           | since the real and imaginary part of      | */
            /*           | the eigenvector are stored in consecutive | */
            /*           | columns, we need to take the magnitude    | */
            /*           | of the last components of the two vectors | */
            /*           %-------------------------------------------% */

            if (iconj == 0)
            {
                bounds[i__] = *rnorm * dlapy2_(&workl[i__], &workl[i__ + 1]);
                bounds[i__ + 1] = bounds[i__];
                iconj = 1;
            }
            else
            {
                iconj = 0;
            }
        }
        /* L20: */
    }

    if (msglvl > 2)
    {
        dvout_(&debug_1.logfil, n, &ritzr[1], &debug_1.ndigit,
               "_neigh: Real"
               " part of the eigenvalues of H",
               (ftnlen)41);
        dvout_(&debug_1.logfil, n, &ritzi[1], &debug_1.ndigit,
               "_neigh: Imag"
               "inary part of the eigenvalues of H",
               (ftnlen)46);
        dvout_(&debug_1.logfil, n, &bounds[1], &debug_1.ndigit,
               "_neigh: Rit"
               "z estimates for the eigenvalues of H",
               (ftnlen)47);
    }

    arscnd_(&t1);
    timing_1.tneigh += t1 - t0;

L9000:
    return 0;

    /*     %---------------% */
    /*     | End of dneigh | */
    /*     %---------------% */

} /* dneigh_ */
