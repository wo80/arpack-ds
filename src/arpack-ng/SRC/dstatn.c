/* SRC\dstatn.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

/* ------------------------------------------- */
/* Initialize statistic and timing information */
/* for nonsymmetric Arnoldi code.              */
/* ------------------------------------------- */
int dstatn_(void)
{

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

    /* -------------------------------------------------- */
    /* User time including reverse communication overhead */
    /* -------------------------------------------------- */

    timing_1.tmvopx = 0.f;
    timing_1.tmvbx = 0.f;

    return 0;

    /* ------------- */
    /* End of dstatn */
    /* ------------- */

} /* dstatn_ */
