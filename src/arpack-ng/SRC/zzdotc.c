/* SRC\zzdotc.f -- translated by f2c (version 20230428). */

#include "arpack_internal.h"

void zzdotc_(a_dcomplex* ret_val, a_int* n, a_dcomplex* zx, a_int* incx, a_dcomplex* zy, a_int* incy)
{
    /* System generated locals */
    a_int i__1, i__2;
    a_dcomplex z__1, z__2, z__3;

    /* Local variables */
    a_int i, ix, iy;
    a_dcomplex ztemp;

    /*     forms the dot product of a vector. */
    /*     jack dongarra, 3/11/78. */
    /*     modified 12/3/93, array(1) declarations changed to array(*) */

    /* Parameter adjustments */
    --zy;
    --zx;

    ztemp.r = 0., ztemp.i = 0.;
    ret_val->r = 0., ret_val->i = 0.;
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
        ar_d_cnjg(&z__3, &zx[ix]);
        i__2 = iy;
        z__2.r = z__3.r * zy[i__2].r - z__3.i * zy[i__2].i, z__2.i = z__3.r * zy[i__2].i + z__3.i * zy[i__2].r;
        z__1.r = ztemp.r + z__2.r, z__1.i = ztemp.i + z__2.i;
        ztemp.r = z__1.r, ztemp.i = z__1.i;
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
        ar_d_cnjg(&z__3, &zx[i]);
        i__2 = i;
        z__2.r = z__3.r * zy[i__2].r - z__3.i * zy[i__2].i, z__2.i = z__3.r * zy[i__2].i + z__3.i * zy[i__2].r;
        z__1.r = ztemp.r + z__2.r, z__1.i = ztemp.i + z__2.i;
        ztemp.r = z__1.r, ztemp.i = z__1.i;
        /* L30: */
    }
    ret_val->r = ztemp.r, ret_val->i = ztemp.i;
    return;
} /* zzdotc_ */
