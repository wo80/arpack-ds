/* SRC\sstatn.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

/* Common Block Declarations */

Extern struct
{
    integer nopx, nbx, nrorth, nitref, nrstrt;
    real tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv, tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv, tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv, tmvopx, tmvbx, tgetv0, titref, trvec;
} timing_;

#define timing_1 timing_

/*     %---------------------------------------------% */
/*     | Initialize statistic and timing information | */
/*     | for nonsymmetric Arnoldi code.              | */
/*     %---------------------------------------------% */

/* \Author */
/*     Danny Sorensen               Phuong Vu */
/*     Richard Lehoucq              CRPC / Rice University */
/*     Dept. of Computational &     Houston, Texas */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \SCCS Information: @(#) */
/* FILE: statn.F   SID: 2.4   DATE OF SID: 4/20/96   RELEASE: 2 */

int sstatn_(void)
{

    /*     %--------------------------------% */
    /*     | See stat.doc for documentation | */
    /*     %--------------------------------% */

    /*     %-----------------------% */
    /*     | Executable Statements | */
    /*     %-----------------------% */

    /*     %--------------------------------% */
    /*     | See stat.doc for documentation | */
    /*     %--------------------------------% */

    /* \SCCS Information: @(#) */
    /* FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2 */

    timing_1.nopx = 0;
    timing_1.nbx = 0;
    timing_1.nrorth = 0;
    timing_1.nitref = 0;
    timing_1.nrstrt = 0;

    timing_1.tnaupd = 0.f;
    timing_1.tnaup2 = 0.f;
    timing_1.tnaitr = 0.f;
    timing_1.tneigh = 0.f;
    timing_1.tngets = 0.f;
    timing_1.tnapps = 0.f;
    timing_1.tnconv = 0.f;
    timing_1.titref = 0.f;
    timing_1.tgetv0 = 0.f;
    timing_1.trvec = 0.f;

    /*     %----------------------------------------------------% */
    /*     | User time including reverse communication overhead | */
    /*     %----------------------------------------------------% */

    timing_1.tmvopx = 0.f;
    timing_1.tmvbx = 0.f;

    return 0;

    /*     %---------------% */
    /*     | End of sstatn | */
    /*     %---------------% */

} /* sstatn_ */
