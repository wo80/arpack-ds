#pragma once

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "arpack.h"
#include "blas.h"
#include "lapack.h"

#ifdef __cplusplus
extern "C"
{
#endif

    /* BEGIN: private interface */

    int arscnd_(float*);

    int ivout_(a_int, a_int*, a_int, char*);

    /* a_fcomplex */

    int cgetv0_(a_int*, char*, a_int*, a_bool*, a_int*, a_int*, a_fcomplex*, a_int*, a_fcomplex*, float*, a_int*, a_fcomplex*, a_int*);
    int cnaitr_(a_int*, char*, a_int*, a_int*, a_int*, a_int*, a_fcomplex*, float*, a_fcomplex*, a_int*, a_fcomplex*, a_int*, a_int*, a_fcomplex*, a_int*);
    int cnapps_(a_int*, a_int*, a_int*, a_fcomplex*, a_fcomplex*, a_int*, a_fcomplex*, a_int*, a_fcomplex*, a_fcomplex*, a_int*, a_fcomplex*, a_fcomplex*);
    int cnaup2_(a_int*, char*, a_int*, char*, a_int*, a_int*, float*, a_fcomplex*, a_int*, a_int*, a_int*, a_int*, a_fcomplex*, a_int*, a_fcomplex*, a_int*, a_fcomplex*, a_fcomplex*, a_fcomplex*, a_int*, a_fcomplex*, a_int*, a_fcomplex*, float*, a_int*);
    int cneigh_(float*, a_int*, a_fcomplex*, a_int*, a_fcomplex*, a_fcomplex*, a_fcomplex*, a_int*, a_fcomplex*, float*, a_int*);
    int cngets_(a_int*, char*, a_int*, a_int*, a_fcomplex*, a_fcomplex*);
    int csortc_(char*, a_bool*, a_int*, a_fcomplex*, a_fcomplex*);
    int cstatn_(void);

    int cmout_(a_int, a_int, a_fcomplex*, a_int, a_int, char*);
    int cvout_(a_int, a_fcomplex*, a_int, char*);

    /* double */

    int dgetv0_(a_int*, char*, a_int*, a_bool*, a_int*, a_int*, double*, a_int*, double*, double*, a_int*, double*, a_int*);
    int dnaitr_(a_int*, char*, a_int*, a_int*, a_int*, a_int*, double*, double*, double*, a_int*, double*, a_int*, a_int*, double*, a_int*);
    int dnapps_(a_int*, a_int*, a_int*, double*, double*, double*, a_int*, double*, a_int*, double*, double*, a_int*, double*, double*);
    int dnaup2_(a_int*, char*, a_int*, char*, a_int*, a_int*, double*, double*, a_int*, a_int*, a_int*, a_int*, double*, a_int*, double*, a_int*, double*, double*, double*, double*, a_int*, double*, a_int*, double*, a_int*);
    int dnconv_(a_int*, double*, double*, double*, double*, a_int*);
    int dneigh_(double*, a_int*, double*, a_int*, double*, double*, double*, double*, a_int*, double*, a_int*);
    int dngets_(a_int*, char*, a_int*, a_int*, double*, double*, double*, double*, double*);
    int dsaitr_(a_int*, char*, a_int*, a_int*, a_int*, a_int*, double*, double*, double*, a_int*, double*, a_int*, a_int*, double*, a_int*);
    int dsapps_(a_int*, a_int*, a_int*, double*, double*, a_int*, double*, a_int*, double*, double*, a_int*, double*);
    int dsaup2_(a_int*, char*, a_int*, char*, a_int*, a_int*, double*, double*, a_int*, a_int*, a_int*, a_int*, double*, a_int*, double*, a_int*, double*, double*, double*, a_int*, double*, a_int*, double*, a_int*);
    int dsconv_(a_int*, double*, double*, double*, a_int*);
    int dseigt_(double*, a_int*, double*, a_int*, double*, double*, double*, a_int*);
    int dsesrt_(char*, a_bool*, a_int*, double*, a_int*, double*, a_int*);
    int dsgets_(a_int*, char*, a_int*, a_int*, double*, double*, double*);
    int dsortc_(char*, a_bool*, a_int*, double*, double*, double*);
    int dsortr_(char*, a_bool*, a_int*, double*, double*);
    int dstatn_(void);
    int dstats_(void);
    int dstqrb_(a_int*, double*, double*, double*, double*, a_int*);

    int dmout_(a_int, a_int, double*, a_int, a_int, char*);
    int dvout_(a_int, double*, a_int, char*);

    /* single */

    int sgetv0_(a_int*, char*, a_int*, a_bool*, a_int*, a_int*, float*, a_int*, float*, float*, a_int*, float*, a_int*);
    int snaitr_(a_int*, char*, a_int*, a_int*, a_int*, a_int*, float*, float*, float*, a_int*, float*, a_int*, a_int*, float*, a_int*);
    int snapps_(a_int*, a_int*, a_int*, float*, float*, float*, a_int*, float*, a_int*, float*, float*, a_int*, float*, float*);
    int snaup2_(a_int*, char*, a_int*, char*, a_int*, a_int*, float*, float*, a_int*, a_int*, a_int*, a_int*, float*, a_int*, float*, a_int*, float*, float*, float*, float*, a_int*, float*, a_int*, float*, a_int*);
    int snconv_(a_int*, float*, float*, float*, float*, a_int*);
    int sneigh_(float*, a_int*, float*, a_int*, float*, float*, float*, float*, a_int*, float*, a_int*);
    int sngets_(a_int*, char*, a_int*, a_int*, float*, float*, float*, float*, float*);
    int ssaitr_(a_int*, char*, a_int*, a_int*, a_int*, a_int*, float*, float*, float*, a_int*, float*, a_int*, a_int*, float*, a_int*);
    int ssapps_(a_int*, a_int*, a_int*, float*, float*, a_int*, float*, a_int*, float*, float*, a_int*, float*);
    int ssaup2_(a_int*, char*, a_int*, char*, a_int*, a_int*, float*, float*, a_int*, a_int*, a_int*, a_int*, float*, a_int*, float*, a_int*, float*, float*, float*, a_int*, float*, a_int*, float*, a_int*);
    int ssconv_(a_int*, float*, float*, float*, a_int*);
    int sseigt_(float*, a_int*, float*, a_int*, float*, float*, float*, a_int*);
    int ssesrt_(char*, a_bool*, a_int*, float*, a_int*, float*, a_int*);
    int ssgets_(a_int*, char*, a_int*, a_int*, float*, float*, float*);
    int ssortc_(char*, a_bool*, a_int*, float*, float*, float*);
    int ssortr_(char*, a_bool*, a_int*, float*, float*);
    int sstatn_(void);
    int sstats_(void);
    int sstqrb_(a_int*, float*, float*, float*, float*, a_int*);

    int smout_(a_int, a_int, float*, a_int, a_int, char*);
    int svout_(a_int, float*, a_int, char*);

    /* a_dcomplex */

    int zgetv0_(a_int*, char*, a_int*, a_bool*, a_int*, a_int*, a_dcomplex*, a_int*, a_dcomplex*, double*, a_int*, a_dcomplex*, a_int*);
    int znaitr_(a_int*, char*, a_int*, a_int*, a_int*, a_int*, a_dcomplex*, double*, a_dcomplex*, a_int*, a_dcomplex*, a_int*, a_int*, a_dcomplex*, a_int*);
    int znapps_(a_int*, a_int*, a_int*, a_dcomplex*, a_dcomplex*, a_int*, a_dcomplex*, a_int*, a_dcomplex*, a_dcomplex*, a_int*, a_dcomplex*, a_dcomplex*);
    int znaup2_(a_int*, char*, a_int*, char*, a_int*, a_int*, double*, a_dcomplex*, a_int*, a_int*, a_int*, a_int*, a_dcomplex*, a_int*, a_dcomplex*, a_int*, a_dcomplex*, a_dcomplex*, a_dcomplex*, a_int*, a_dcomplex*, a_int*, a_dcomplex*, double*, a_int*);
    int zneigh_(double*, a_int*, a_dcomplex*, a_int*, a_dcomplex*, a_dcomplex*, a_dcomplex*, a_int*, a_dcomplex*, double*, a_int*);
    int zngets_(a_int*, char*, a_int*, a_int*, a_dcomplex*, a_dcomplex*);
    int zsortc_(char*, a_bool*, a_int*, a_dcomplex*, a_dcomplex*);
    int zstatn_(void);

    int zmout_(a_int, a_int, a_dcomplex*, a_int, a_int, char*);
    int zvout_(a_int, a_dcomplex*, a_int, char*);

    /* END: private interface */

#ifdef __cplusplus
}
#endif
