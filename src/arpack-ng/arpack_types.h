#pragma once

/* f2c.h */

typedef struct { float  r, i; } a_fcomplex;
typedef struct { double r, i; } a_dcomplex;
typedef int a_int;
typedef int a_bool;

#define TRUE_ (1)
#define FALSE_ (0)

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (double)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (double)min(a,b)
#define dmax(a,b) (double)max(a,b)

/* Common blocks */

extern struct {
    int logfil, ndigit, mgetv0, msaupd, msaup2, msaitr, mseigt, msapps,
        msgets, mseupd, mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets,
        mneupd, mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd;
} debug_;

#define debug_1 debug_

extern struct {
    int nopx, nbx, nrorth, nitref, nrstrt;
    float tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv, tnaupd,
        tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv, tcaupd, tcaup2,
        tcaitr, tceigh, tcgets, tcapps, tcconv, tmvopx, tmvbx, tgetv0,
        titref, trvec;
} timing_;

#define timing_1 timing_
