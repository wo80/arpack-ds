
struct {
    int logfil, ndigit, mgetv0, msaupd, msaup2, msaitr, mseigt, msapps,
        msgets, mseupd, mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets,
        mneupd, mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd;
} debug_;

struct {
    int nopx, nbx, nrorth, nitref, nrstrt;
    float tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv, tnaupd,
        tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv, tcaupd, tcaup2,
        tcaitr, tceigh, tcgets, tcapps, tcconv, tmvopx, tmvbx, tgetv0,
        titref, trvec;
} timing_;
