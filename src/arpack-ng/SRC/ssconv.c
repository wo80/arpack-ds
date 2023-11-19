/* SRC\ssconv.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

/* Table of constant values */

static double TWO_THIRDS = .66666666666666663;

/* ----------------------------------------------------------------------- */
/* \BeginDoc */

/* \Name: ssconv */

/* \Description: */
/*  Convergence testing for the symmetric Arnoldi eigenvalue routine. */

/* \Usage: */
/*  call ssconv */
/*     ( N, RITZ, BOUNDS, TOL, NCONV ) */

/* \Arguments */
/*  N       Integer.  (INPUT) */
/*          Number of Ritz values to check for convergence. */

/*  RITZ    Real array of length N.  (INPUT) */
/*          The Ritz values to be checked for convergence. */

/*  BOUNDS  Real array of length N.  (INPUT) */
/*          Ritz estimates associated with the Ritz values in RITZ. */

/*  TOL     Real scalar.  (INPUT) */
/*          Desired relative accuracy for a Ritz value to be considered */
/*          "converged". */

/*  NCONV   Integer scalar.  (OUTPUT) */
/*          Number of "converged" Ritz values. */

/* \EndDoc */

/* ----------------------------------------------------------------------- */

/* \BeginLib */

/* \Routines called: */
/*     arscnd  ARPACK utility routine for timing. */
/*     slamch  LAPACK routine that determines machine constants. */

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

int ssconv_(a_int *n, float *ritz, float *bounds, float *tol, a_int *nconv)
{
    /* System generated locals */
    a_int i__1;
    float r__1, r__2, r__3;
    double d__1;

    /* Builtin functions */
    double pow_dd(double *, double *);

    /* Local variables */
    a_int i;
    static float t0, t1;
    float eps23, temp;

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

    eps23 = slamch_("Epsilon-Machine", (ftnlen)15);
    d__1 = (double)eps23;
    eps23 = pow_dd(&d__1, &TWO_THIRDS);

    *nconv = 0;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {

        /*        %-----------------------------------------------------% */
        /*        | The i-th Ritz value is considered "converged"       | */
        /*        | when: bounds(i) .le. TOL*max(eps23, abs(ritz(i)))   | */
        /*        %-----------------------------------------------------% */

        /* Computing MAX */
        r__2 = eps23, r__3 = (r__1 = ritz[i], dabs(r__1));
        temp = dmax(r__2, r__3);
        if (bounds[i] <= *tol * temp)
        {
            ++(*nconv);
        }

        /* L10: */
    }

    arscnd_(&t1);
    timing_1.tsconv += t1 - t0;

    return 0;

    /*     %---------------% */
    /*     | End of ssconv | */
    /*     %---------------% */

} /* ssconv_ */
