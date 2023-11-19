/* SRC\ccdotc.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

/* Complex */ VOID ccdotc_(complex *ret_val, integer *n, complex *zx, integer *incx, complex *zy, integer *incy)
{
    /* System generated locals */
    a_int i__1, i__2;
    a_fcomplex q__1, q__2, q__3;

    /* Builtin functions */

    /* Local variables */
    a_int i, ix, iy;
    a_fcomplex ztemp;

    /*     forms the dot product of a vector. */
    /*     jack dongarra, 3/11/78. */
    /*     modified 12/3/93, array(1) declarations changed to array(*) */

    /* Parameter adjustments */
    --zy;
    --zx;

    /* Function Body */
    ztemp.r = 0.f, ztemp.i = 0.f;
    ret_val->r = 0.f, ret_val->i = 0.f;
    if (*n <= 0)
    {
        return;
    }
    if (*incx == 1 && *incy == 1)
    {
        goto L20;
    }

    /*        code for unequal increments or equal increments */
    /*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0)
    {
        ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0)
    {
        iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
        ar_r_cnjg(&q__3, &zx[ix]);
        i__2 = iy;
        q__2.r = q__3.r * zy[i__2].r - q__3.i * zy[i__2].i, q__2.i = q__3.r * zy[i__2].i + q__3.i * zy[i__2].r;
        q__1.r = ztemp.r + q__2.r, q__1.i = ztemp.i + q__2.i;
        ztemp.r = q__1.r, ztemp.i = q__1.i;
        ix += *incx;
        iy += *incy;
        /* L10: */
    }
    ret_val->r = ztemp.r, ret_val->i = ztemp.i;
    return;

    /*        code for both increments equal to 1 */

L20:
    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
        ar_r_cnjg(&q__3, &zx[i]);
        i__2 = i;
        q__2.r = q__3.r * zy[i__2].r - q__3.i * zy[i__2].i, q__2.i = q__3.r * zy[i__2].i + q__3.i * zy[i__2].r;
        q__1.r = ztemp.r + q__2.r, q__1.i = ztemp.i + q__2.i;
        ztemp.r = q__1.r, ztemp.i = q__1.i;
        /* L30: */
    }
    ret_val->r = ztemp.r, ret_val->i = ztemp.i;
    return;
} /* ccdotc_ */
