/* SRC\cngets.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

/* Table of constant values */

static a_bool b_true = TRUE_;
static a_int i_one = 1;

/* \BeginDoc */

/* \Name: cngets */

/* \Description: */
/*  Given the eigenvalues of the upper Hessenberg matrix H, */
/*  computes the NP shifts AMU that are zeros of the polynomial of */
/*  degree NP which filters out components of the unwanted eigenvectors */
/*  corresponding to the AMU's based on some given criteria. */

/*  NOTE: call this even in the case of user specified shifts in order */
/*  to sort the eigenvalues, and error bounds of H for later use. */

/* \Usage: */
/*  call cngets */
/*      ( ISHIFT, WHICH, KEV, NP, RITZ, BOUNDS ) */

/* \Arguments */
/*  ISHIFT  Integer.  (INPUT) */
/*          Method for selecting the implicit shifts at each iteration. */
/*          ISHIFT = 0: user specified shifts */
/*          ISHIFT = 1: exact shift with respect to the matrix H. */

/*  WHICH   Character*2.  (INPUT) */
/*          Shift selection criteria. */
/*          'LM' -> want the KEV eigenvalues of largest magnitude. */
/*          'SM' -> want the KEV eigenvalues of smallest magnitude. */
/*          'LR' -> want the KEV eigenvalues of largest REAL part. */
/*          'SR' -> want the KEV eigenvalues of smallest REAL part. */
/*          'LI' -> want the KEV eigenvalues of largest imaginary part. */
/*          'SI' -> want the KEV eigenvalues of smallest imaginary part. */

/*  KEV     Integer.  (INPUT) */
/*          The number of desired eigenvalues. */

/*  NP      Integer.  (INPUT) */
/*          The number of shifts to compute. */

/*  RITZ    Complex array of length KEV+NP.  (INPUT/OUTPUT) */
/*          On INPUT, RITZ contains the the eigenvalues of H. */
/*          On OUTPUT, RITZ are sorted so that the unwanted */
/*          eigenvalues are in the first NP locations and the wanted */
/*          portion is in the last KEV locations.  When exact shifts are */
/*          selected, the unwanted part corresponds to the shifts to */
/*          be applied. Also, if ISHIFT .eq. 1, the unwanted eigenvalues */
/*          are further sorted so that the ones with largest Ritz values */
/*          are first. */

/*  BOUNDS  Complex array of length KEV+NP.  (INPUT/OUTPUT) */
/*          Error bounds corresponding to the ordering in RITZ. */

/* \EndDoc */

/* ----------------------------------------------------------------------- */

/* \BeginLib */

/* \Local variables: */
/*     xxxxxx  Complex */

/* \Routines called: */
/*     csortc  ARPACK sorting routine. */
/*     ivout   ARPACK utility routine that prints integers. */
/*     arscnd  ARPACK utility routine for timing. */
/*     cvout   ARPACK utility routine that prints vectors. */

/* \Author */
/*     Danny Sorensen               Phuong Vu */
/*     Richard Lehoucq              CRPC / Rice University */
/*     Dept. of Computational &     Houston, Texas */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \SCCS Information: @(#) */
/* FILE: ngets.F   SID: 2.2   DATE OF SID: 4/20/96   RELEASE: 2 */

/* \Remarks */
/*     1. This routine does not keep complex conjugate pairs of */
/*        eigenvalues together. */

/* \EndLib */

/* ----------------------------------------------------------------------- */

int cngets_(a_int *ishift, char *which, a_int *kev, a_int *np, a_fcomplex *ritz, a_fcomplex *bounds, ftnlen which_len)
{
    /* System generated locals */
    a_int i__1;

    /* Local variables */
    static float t0, t1;
    extern int cvout_(a_int *, a_int *, a_fcomplex *, a_int *, char *, ftnlen), ivout_(a_int *, a_int *, a_int *, a_int *, char *, ftnlen), arscnd_(float *);
    a_int msglvl;
    extern int csortc_(char *, a_bool *, a_int *, a_fcomplex *, a_fcomplex *, ftnlen);

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

    /*     %---------------% */
    /*     | Local Scalars | */
    /*     %---------------% */

    /*     %----------------------% */
    /*     | External Subroutines | */
    /*     %----------------------% */

    /*     %-----------------------% */
    /*     | Executable Statements | */
    /*     %-----------------------% */

    /*     %-------------------------------% */
    /*     | Initialize timing statistics  | */
    /*     | & message level for debugging | */
    /*     %-------------------------------% */

    /* Parameter adjustments */
    --bounds;
    --ritz;

    /* Function Body */
    arscnd_(&t0);
    msglvl = debug_1.mcgets;

    i__1 = *kev + *np;
    csortc_(which, &b_true, &i__1, &ritz[1], &bounds[1], (ftnlen)2);

    if (*ishift == 1)
    {

        /*        %-------------------------------------------------------% */
        /*        | Sort the unwanted Ritz values used as shifts so that  | */
        /*        | the ones with largest Ritz estimates are first        | */
        /*        | This will tend to minimize the effects of the         | */
        /*        | forward instability of the iteration when the shifts  | */
        /*        | are applied in subroutine cnapps.                     | */
        /*        | Be careful and use 'SM' since we want to sort BOUNDS! | */
        /*        %-------------------------------------------------------% */

        csortc_("SM", &b_true, np, &bounds[1], &ritz[1], (ftnlen)2);
    }

    arscnd_(&t1);
    timing_1.tcgets += t1 - t0;

    if (msglvl > 0)
    {
        ivout_(&debug_1.logfil, &i_one, kev, &debug_1.ndigit, "_ngets: KEV is", (ftnlen)14);
        ivout_(&debug_1.logfil, &i_one, np, &debug_1.ndigit, "_ngets: NP is", (ftnlen)13);
        i__1 = *kev + *np;
        cvout_(&debug_1.logfil, &i__1, &ritz[1], &debug_1.ndigit,
               "_ngets: E"
               "igenvalues of current H matrix ",
               (ftnlen)40);
        i__1 = *kev + *np;
        cvout_(&debug_1.logfil, &i__1, &bounds[1], &debug_1.ndigit,
               "_ngets:"
               " Ritz estimates of the current KEV+NP Ritz values",
               (ftnlen)56);
    }

    return 0;

    /*     %---------------% */
    /*     | End of cngets | */
    /*     %---------------% */

} /* cngets_ */
