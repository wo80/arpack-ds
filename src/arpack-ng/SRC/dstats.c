/* SRC\dstats.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

/* ------------------------------------------- */
/* Initialize statistic and timing information */
/* for symmetric Arnoldi code.                 */
/* ------------------------------------------- */
int dstats_(void)
{

    timing_1.nopx = 0;
    timing_1.nbx = 0;
    timing_1.nrorth = 0;
    timing_1.nitref = 0;
    timing_1.nrstrt = 0;
    timing_1.tsaupd = 0.f;
    timing_1.tsaup2 = 0.f;
    timing_1.tsaitr = 0.f;
    timing_1.tseigt = 0.f;
    timing_1.tsgets = 0.f;
    timing_1.tsapps = 0.f;
    timing_1.tsconv = 0.f;
    timing_1.titref = 0.f;
    timing_1.tgetv0 = 0.f;
    timing_1.trvec = 0.f;
    /* -------------------------------------------------- */
    /* User time including reverse communication overhead */
    /* -------------------------------------------------- */
    timing_1.tmvopx = 0.f;
    timing_1.tmvbx = 0.f;
    return 0;

    /*     End of dstats */

} /* dstats_ */
