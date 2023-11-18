/* SRC\dsconv.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

/* Table of constant values */

static double TWO_THIRDS = .66666666666666663;

/* ----------------------------------------------------------------------- */
/* \BeginDoc */

/* \Name: dsconv */

/* \Description: */
/*  Convergence testing for the symmetric Arnoldi eigenvalue routine. */

/* \Usage: */
/*  call dsconv */
/*     ( N, RITZ, BOUNDS, TOL, NCONV ) */

/* \Arguments */
/*  N       Integer.  (INPUT) */
/*          Number of Ritz values to check for convergence. */

/*  RITZ    Double precision array of length N.  (INPUT) */
/*          The Ritz values to be checked for convergence. */

/*  BOUNDS  Double precision array of length N.  (INPUT) */
/*          Ritz estimates associated with the Ritz values in RITZ. */

/*  TOL     Double precision scalar.  (INPUT) */
/*          Desired relative accuracy for a Ritz value to be considered */
/*          "converged". */

/*  NCONV   Integer scalar.  (OUTPUT) */
/*          Number of "converged" Ritz values. */

/* \EndDoc */

/* ----------------------------------------------------------------------- */

/* \BeginLib */

/* \Routines called: */
/*     arscnd  ARPACK utility routine for timing. */
/*     dlamch  LAPACK routine that determines machine constants. */

/* \Author */
/*     Danny Sorensen               Phuong Vu */
/*     Richard Lehoucq              CRPC / Rice University */
/*     Dept. of Computational &     Houston, Texas */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \SCCS Information: @(#) */
/* FILE: sconv.F   SID: 2.4   DATE OF SID: 4/19/96   RELEASE: 2 */

/* \Remarks */
/*     1. Starting with version 2.4, this routine no longer uses the */
/*        Parlett strategy using the gap conditions. */

/* \EndLib */

/* ----------------------------------------------------------------------- */

int dsconv_(a_int *n, double *ritz, double *bounds, double *tol, a_int *nconv)
{
    /* System generated locals */
    a_int i__1;
    double d__1, d__2, d__3;

    /* Builtin functions */
    double pow_dd(double *, double *);

    /* Local variables */
    a_int i__;
    static float t0, t1;
    double eps23, temp;
    extern double dlamch_(char *, ftnlen);
    extern int arscnd_(float *);

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

    /*     %---------------% */
    /*     | Local Scalars | */
    /*     %---------------% */

    /*     %-------------------% */
    /*     | External routines | */
    /*     %-------------------% */

    /*     %---------------------% */
    /*     | Intrinsic Functions | */
    /*     %---------------------% */

    /*     %-----------------------% */
    /*     | Executable Statements | */
    /*     %-----------------------% */

    /* Parameter adjustments */
    --bounds;
    --ritz;

    /* Function Body */
    arscnd_(&t0);

    eps23 = dlamch_("Epsilon-Machine", (ftnlen)15);
    eps23 = pow_dd(&eps23, &TWO_THIRDS);

    *nconv = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__)
    {

        /*        %-----------------------------------------------------% */
        /*        | The i-th Ritz value is considered "converged"       | */
        /*        | when: bounds(i) .le. TOL*max(eps23, abs(ritz(i)))   | */
        /*        %-----------------------------------------------------% */

        /* Computing MAX */
        d__2 = eps23, d__3 = (d__1 = ritz[i__], abs(d__1));
        temp = max(d__2, d__3);
        if (bounds[i__] <= *tol * temp)
        {
            ++(*nconv);
        }

        /* L10: */
    }

    arscnd_(&t1);
    timing_1.tsconv += t1 - t0;

    return 0;

    /*     %---------------% */
    /*     | End of dsconv | */
    /*     %---------------% */

} /* dsconv_ */
