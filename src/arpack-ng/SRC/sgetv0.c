/* SRC\sgetv0.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

/* Table of constant values */

static a_int i_one = 1;
static float s_one = 1.f;
static float s_zero = 0.f;
static float s_n1 = -1.f;

/* ----------------------------------------------------------------------- */
/* \BeginDoc */

/* \Name: sgetv0 */

/* \Description: */
/*  Generate a random initial residual vector for the Arnoldi process. */
/*  Force the residual vector to be in the range of the operator OP. */

/* \Usage: */
/*  call sgetv0 */
/*     ( IDO, BMAT, ITRY, INITV, N, J, V, LDV, RESID, RNORM, */
/*       IPNTR, WORKD, IERR ) */

/* \Arguments */
/*  IDO     Integer.  (INPUT/OUTPUT) */
/*          Reverse communication flag.  IDO must be zero on the first */
/*          call to sgetv0. */
/*          ------------------------------------------------------------- */
/*          IDO =  0: first call to the reverse communication interface */
/*          IDO = -1: compute  Y = OP * X  where */
/*                    IPNTR(1) is the pointer into WORKD for X, */
/*                    IPNTR(2) is the pointer into WORKD for Y. */
/*                    This is for the initialization phase to force the */
/*                    starting vector into the range of OP. */
/*          IDO =  2: compute  Y = B * X  where */
/*                    IPNTR(1) is the pointer into WORKD for X, */
/*                    IPNTR(2) is the pointer into WORKD for Y. */
/*          IDO = 99: done */
/*          ------------------------------------------------------------- */

/*  BMAT    Character*1.  (INPUT) */
/*          BMAT specifies the type of the matrix B in the (generalized) */
/*          eigenvalue problem A*x = lambda*B*x. */
/*          B = 'I' -> standard eigenvalue problem A*x = lambda*x */
/*          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x */

/*  ITRY    Integer.  (INPUT) */
/*          ITRY counts the number of times that sgetv0 is called. */
/*          It should be set to 1 on the initial call to sgetv0. */

/*  INITV   Logical variable.  (INPUT) */
/*          .TRUE.  => the initial residual vector is given in RESID. */
/*          .FALSE. => generate a random initial residual vector. */

/*  N       Integer.  (INPUT) */
/*          Dimension of the problem. */

/*  J       Integer.  (INPUT) */
/*          Index of the residual vector to be generated, with respect to */
/*          the Arnoldi process.  J > 1 in case of a "restart". */

/*  V       Real N by J array.  (INPUT) */
/*          The first J-1 columns of V contain the current Arnoldi basis */
/*          if this is a "restart". */

/*  LDV     Integer.  (INPUT) */
/*          Leading dimension of V exactly as declared in the calling */
/*          program. */

/*  RESID   Real array of length N.  (INPUT/OUTPUT) */
/*          Initial residual vector to be generated.  If RESID is */
/*          provided, force RESID into the range of the operator OP. */

/*  RNORM   Real scalar.  (OUTPUT) */
/*          B-norm of the generated residual. */

/*  IPNTR   Integer array of length 3.  (OUTPUT) */

/*  WORKD   Real work array of length 2*N.  (REVERSE COMMUNICATION). */
/*          On exit, WORK(1:N) = B*RESID to be used in SSAITR. */

/*  IERR    Integer.  (OUTPUT) */
/*          =  0: Normal exit. */
/*          = -1: Cannot generate a nontrivial restarted residual vector */
/*                in the range of the operator OP. */

/* \EndDoc */

/* ----------------------------------------------------------------------- */

/* \BeginLib */

/* \Local variables: */
/*     xxxxxx  real */

/* \References: */
/*  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in */
/*     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992), */
/*     pp 357-385. */
/*  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly */
/*     Restarted Arnoldi Iteration", Rice University Technical Report */
/*     TR95-13, Department of Computational and Applied Mathematics. */

/* \Routines called: */
/*     arscnd  ARPACK utility routine for timing. */
/*     svout   ARPACK utility routine for vector output. */
/*     slarnv  LAPACK routine for generating a random vector. */
/*     sgemv   Level 2 BLAS routine for matrix vector multiplication. */
/*     scopy   Level 1 BLAS that copies one vector to another. */
/*     sdot    Level 1 BLAS that computes the scalar product of two vectors. */
/*     snrm2   Level 1 BLAS that computes the norm of a vector. */

/* \Author */
/*     Danny Sorensen               Phuong Vu */
/*     Richard Lehoucq              CRPC / Rice University */
/*     Dept. of Computational &     Houston, Texas */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \SCCS Information: @(#) */
/* FILE: getv0.F   SID: 2.7   DATE OF SID: 04/07/99   RELEASE: 2 */

/* \EndLib */

/* ----------------------------------------------------------------------- */

int sgetv0_(a_int *ido, const char *bmat, a_int *itry, a_bool *initv, a_int *n, a_int *j,
     float *v, a_int *ldv, float *resid, float *rnorm, a_int *ipntr, float *workd, a_int *ierr)
{
    /* Initialized data */

    static a_bool inits = TRUE_;

    /* System generated locals */
    a_int v_dim1, v_offset, i__1;

    /* Builtin functions */

    /* Local variables */
    static float t0, t1, t2, t3;
    a_int jj;
    static a_int iter;
    static a_bool orth;
    static a_int iseed[4];
    a_int idist;
    static a_bool first;
    static float rnorm0;
    static a_int msglvl;

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

    /*     %-----------------% */
    /*     | Data Statements | */
    /*     %-----------------% */

    /* Parameter adjustments */
    --workd;
    --resid;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --ipntr;

    /* Function Body */

    /*     %-----------------------% */
    /*     | Executable Statements | */
    /*     %-----------------------% */

    /*     %-----------------------------------% */
    /*     | Initialize the seed of the LAPACK | */
    /*     | random number generator           | */
    /*     %-----------------------------------% */

    if (inits)
    {
        iseed[0] = 1;
        iseed[1] = 3;
        iseed[2] = 5;
        iseed[3] = 7;
        inits = FALSE_;
    }

    if (*ido == 0)
    {

        /*        %-------------------------------% */
        /*        | Initialize timing statistics  | */
        /*        | & message level for debugging | */
        /*        %-------------------------------% */

        arscnd_(&t0);
        msglvl = debug_1.mgetv0;

        *ierr = 0;
        iter = 0;
        first = FALSE_;
        orth = FALSE_;

        /*        %-----------------------------------------------------% */
        /*        | Possibly generate a random starting vector in RESID | */
        /*        | Use a LAPACK random number generator used by the    | */
        /*        | matrix generation routines.                         | */
        /*        |    idist = 1: uniform (0,1)  distribution;          | */
        /*        |    idist = 2: uniform (-1,1) distribution;          | */
        /*        |    idist = 3: normal  (0,1)  distribution;          | */
        /*        %-----------------------------------------------------% */

        if (!(*initv))
        {
            idist = 2;
            slarnv_(&idist, iseed, n, &resid[1]);
        }

        /*        %----------------------------------------------------------% */
        /*        | Force the starting vector into the range of OP to handle | */
        /*        | the generalized problem when B is possibly (singular).   | */
        /*        %----------------------------------------------------------% */

        arscnd_(&t2);
        if (*itry == 1)
        {
            ++timing_1.nopx;
            ipntr[1] = 1;
            ipntr[2] = *n + 1;
            scopy_(n, &resid[1], &i_one, &workd[1], &i_one);
            *ido = -1;
            goto L9000;
        }
        else if (*itry > 1 && *bmat == 'G')
        {
            scopy_(n, &resid[1], &i_one, &workd[*n + 1], &i_one);
        }
    }

    /*     %-----------------------------------------% */
    /*     | Back from computing OP*(initial-vector) | */
    /*     %-----------------------------------------% */

    if (first)
    {
        goto L20;
    }

    /*     %-----------------------------------------------% */
    /*     | Back from computing OP*(orthogonalized-vector) | */
    /*     %-----------------------------------------------% */

    if (orth)
    {
        goto L40;
    }

    if (*bmat == 'G')
    {
        arscnd_(&t3);
        timing_1.tmvopx += t3 - t2;
    }

    /*     %------------------------------------------------------% */
    /*     | Starting vector is now in the range of OP; r = OP*r; | */
    /*     | Compute B-norm of starting vector.                   | */
    /*     %------------------------------------------------------% */

    arscnd_(&t2);
    first = TRUE_;
    if (*itry == 1)
    {
        scopy_(n, &workd[*n + 1], &i_one, &resid[1], &i_one);
    }
    if (*bmat == 'G')
    {
        ++timing_1.nbx;
        ipntr[1] = *n + 1;
        ipntr[2] = 1;
        *ido = 2;
        goto L9000;
    }
    else if (*bmat == 'I')
    {
        scopy_(n, &resid[1], &i_one, &workd[1], &i_one);
    }

L20:

    if (*bmat == 'G')
    {
        arscnd_(&t3);
        timing_1.tmvbx += t3 - t2;
    }

    first = FALSE_;
    if (*bmat == 'G')
    {
        rnorm0 = sdot_(n, &resid[1], &i_one, &workd[1], &i_one);
        rnorm0 = sqrt((dabs(rnorm0)));
    }
    else if (*bmat == 'I')
    {
        rnorm0 = snrm2_(n, &resid[1], &i_one);
    }
    *rnorm = rnorm0;

    /*     %---------------------------------------------% */
    /*     | Exit if this is the very first Arnoldi step | */
    /*     %---------------------------------------------% */

    if (*j == 1)
    {
        goto L50;
    }

    /*     %---------------------------------------------------------------- */
    /*     | Otherwise need to B-orthogonalize the starting vector against | */
    /*     | the current Arnoldi basis using Gram-Schmidt with iter. ref.  | */
    /*     | This is the case where an invariant subspace is encountered   | */
    /*     | in the middle of the Arnoldi factorization.                   | */
    /*     |                                                               | */
    /*     |       s = V^{T}*B*r;   r = r - V*s;                           | */
    /*     |                                                               | */
    /*     | Stopping criteria used for iter. ref. is discussed in         | */
    /*     | Parlett's book, page 107 and in Gragg & Reichel TOMS paper.   | */
    /*     %---------------------------------------------------------------% */

    orth = TRUE_;
L30:

    i__1 = *j - 1;
    sgemv_("T", n, &i__1, &s_one, &v[v_offset], ldv, &workd[1], &i_one, &s_zero, &workd[*n + 1], &i_one);
    i__1 = *j - 1;
    sgemv_("N", n, &i__1, &s_n1, &v[v_offset], ldv, &workd[*n + 1], &i_one, &s_one, &resid[1], &i_one);

    /*     %----------------------------------------------------------% */
    /*     | Compute the B-norm of the orthogonalized starting vector | */
    /*     %----------------------------------------------------------% */

    arscnd_(&t2);
    if (*bmat == 'G')
    {
        ++timing_1.nbx;
        scopy_(n, &resid[1], &i_one, &workd[*n + 1], &i_one);
        ipntr[1] = *n + 1;
        ipntr[2] = 1;
        *ido = 2;
        goto L9000;
    }
    else if (*bmat == 'I')
    {
        scopy_(n, &resid[1], &i_one, &workd[1], &i_one);
    }

L40:

    if (*bmat == 'G')
    {
        arscnd_(&t3);
        timing_1.tmvbx += t3 - t2;
    }

    if (*bmat == 'G')
    {
        *rnorm = sdot_(n, &resid[1], &i_one, &workd[1], &i_one);
        *rnorm = sqrt((dabs(*rnorm)));
    }
    else if (*bmat == 'I')
    {
        *rnorm = snrm2_(n, &resid[1], &i_one);
    }

    /*     %--------------------------------------% */
    /*     | Check for further orthogonalization. | */
    /*     %--------------------------------------% */

    if (msglvl > 2)
    {
        svout_(1, &rnorm0, debug_1.ndigit, "_getv0: re-orthonalization ; rnorm0 is");
        svout_(1, rnorm, debug_1.ndigit, "_getv0: re-orthonalization ; rnorm is");
    }

    if (*rnorm > rnorm0 * .717f)
    {
        goto L50;
    }

    ++iter;
    if (iter <= 5)
    {

        /*        %-----------------------------------% */
        /*        | Perform iterative refinement step | */
        /*        %-----------------------------------% */

        rnorm0 = *rnorm;
        goto L30;
    }
    else
    {

        /*        %------------------------------------% */
        /*        | Iterative refinement step "failed" | */
        /*        %------------------------------------% */

        i__1 = *n;
        for (jj = 1; jj <= i__1; ++jj)
        {
            resid[jj] = 0.f;
            /* L45: */
        }
        *rnorm = 0.f;
        *ierr = -1;
    }

L50:

    if (msglvl > 0)
    {
        svout_(1, rnorm, debug_1.ndigit, "_getv0: B-norm of initial / restarted starting vector");
    }
    if (msglvl > 3)
    {
        svout_(n, &resid[1], debug_1.ndigit, "_getv0: initial / restarted starting vector");
    }
    *ido = 99;

    arscnd_(&t1);
    timing_1.tgetv0 += t1 - t0;

L9000:
    return 0;

    /*     %---------------% */
    /*     | End of sgetv0 | */
    /*     %---------------% */

} /* sgetv0_ */
