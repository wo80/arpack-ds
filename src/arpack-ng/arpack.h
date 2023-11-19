#pragma once

#include "arpack_types.h"

#ifdef __cplusplus
extern "C"
{
#endif

    /* BEGIN: public interface */

    int cnaupd_(a_int* ido, char* bmat, a_int* n, char* which, a_int* nev, float* tol,
        a_fcomplex* resid, a_int* ncv, a_fcomplex* v, a_int* ldv, a_int* iparam, a_int* ipntr,
        a_fcomplex* workd, a_fcomplex* workl, a_int* lworkl, float* rwork, a_int* info);

    int cneupd_(a_bool* rvec, char* howmny, a_bool* select, a_fcomplex* d, a_fcomplex* z, a_int* ldz,
        a_fcomplex* sigma, a_fcomplex* workev, char* bmat, a_int* n, char* which, a_int* nev,
        float* tol, a_fcomplex* resid, a_int* ncv, a_fcomplex* v, a_int* ldv, a_int* iparam,
        a_int* ipntr, a_fcomplex* workd, a_fcomplex* workl, a_int* lworkl, float* rwork,
        a_int* info);

    int dnaupd_(a_int* ido, char* bmat, a_int* n, char* which, a_int* nev, double* tol,
        double* resid, a_int* ncv, double* v, a_int* ldv, a_int* iparam, a_int* ipntr,
        double* workd, double* workl, a_int* lworkl, a_int* info);

    int dneupd_(a_bool* rvec, char* howmny, a_bool* select, double* dr, double* di, double* z,
        a_int* ldz, double* sigmar, double* sigmai, double* workev, char* bmat, a_int* n,
        char* which, a_int* nev, double* tol, double* resid, a_int* ncv, double* v,
        a_int* ldv, a_int* iparam, a_int* ipntr, double* workd, double* workl, a_int* lworkl,
        a_int* info);

    int dsaupd_(a_int* ido, char* bmat, a_int* n, char* which, a_int* nev, double* tol,
        double* resid, a_int* ncv, double* v, a_int* ldv, a_int* iparam, a_int* ipntr,
        double* workd, double* workl, a_int* lworkl, a_int* info);

    int dseupd_(a_bool* rvec, char* howmny, a_bool* select, double* d, double* z, a_int* ldz,
        double* sigma, char* bmat, a_int* n, char* which, a_int* nev, double* tol,
        double* resid, a_int* ncv, double* v, a_int* ldv, a_int* iparam, a_int* ipntr,
        double* workd, double* workl, a_int* lworkl, a_int* info);

    int snaupd_(a_int* ido, char* bmat, a_int* n, char* which, a_int* nev, float* tol,
        float* resid, a_int* ncv, float* v, a_int* ldv, a_int* iparam, a_int* ipntr,
        float* workd, float* workl, a_int* lworkl, a_int* info);

    int sneupd_(a_bool* rvec, char* howmny, a_bool* select, float* dr, float* di, float* z,
        a_int* ldz, float* sigmar, float* sigmai, float* workev, char* bmat, a_int* n,
        char* which, a_int* nev, float* tol, float* resid, a_int* ncv, float* v,
        a_int* ldv, a_int* iparam, a_int* ipntr, float* workd, float* workl, a_int* lworkl,
        a_int* info);

    int ssaupd_(a_int* ido, char* bmat, a_int* n, char* which, a_int* nev, float* tol,
        float* resid, a_int* ncv, float* v, a_int* ldv, a_int* iparam, a_int* ipntr,
        float* workd, float* workl, a_int* lworkl, a_int* info);

    int sseupd_(a_bool* rvec, char* howmny, a_bool* select, float* d, float* z, a_int* ldz,
        float* sigma, char* bmat, a_int* n, char* which, a_int* nev, float* tol,
        float* resid, a_int* ncv, float* v, a_int* ldv, a_int* iparam, a_int* ipntr,
        float* workd, float* workl, a_int* lworkl, a_int* info);

    int znaupd_(a_int* ido, char* bmat, a_int* n, char* which, a_int* nev, double* tol,
        a_dcomplex* resid, a_int* ncv, a_dcomplex* v, a_int* ldv, a_int* iparam, a_int* ipntr,
        a_dcomplex* workd, a_dcomplex* workl, a_int* lworkl, double* rwork, a_int* info);

    int zneupd_(a_bool* rvec, char* howmny, a_bool* select, a_dcomplex* d, a_dcomplex* z, a_int* ldz,
        a_dcomplex* sigma, a_dcomplex* workev, char* bmat, a_int* n, char* which, a_int* nev,
        double* tol, a_dcomplex* resid, a_int* ncv, a_dcomplex* v, a_int* ldv, a_int* iparam,
        a_int* ipntr, a_dcomplex* workd, a_dcomplex* workl, a_int* lworkl, double* rwork,
        a_int* info);

    /* END: public interface */

#ifdef __cplusplus
}
#endif
