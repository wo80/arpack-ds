/* D:\projects\Fortran\arpack-ng-3.9.1-patched\SRC\ccdotc.f -- translated by f2c (version 20230428).
   You must link the resulting object file with libf2c:
    on Microsoft Windows system, link with libf2c.lib;
    on Linux or Unix systems, link with .../path/to/libf2c.a -lm
    or, if you install libf2c.a in a standard place, with -lf2c -lm
    -- in that order, at the end of the command line, as in
        cc *.o -lf2c -lm
    Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

        http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Complex */ VOID ccdotc_(complex *ret_val, integer *n, complex *zx, integer *incx, complex *zy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2;
    complex q__1, q__2, q__3;

    /* Builtin functions */
    void r_cnjg(complex *, complex *);

    /* Local variables */
    integer i__, ix, iy;
    complex ztemp;

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
    for (i__ = 1; i__ <= i__1; ++i__)
    {
        r_cnjg(&q__3, &zx[ix]);
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
    for (i__ = 1; i__ <= i__1; ++i__)
    {
        r_cnjg(&q__3, &zx[i__]);
        i__2 = i__;
        q__2.r = q__3.r * zy[i__2].r - q__3.i * zy[i__2].i, q__2.i = q__3.r * zy[i__2].i + q__3.i * zy[i__2].r;
        q__1.r = ztemp.r + q__2.r, q__1.i = ztemp.i + q__2.i;
        ztemp.r = q__1.r, ztemp.i = q__1.i;
        /* L30: */
    }
    ret_val->r = ztemp.r, ret_val->i = ztemp.i;
    return;
} /* ccdotc_ */
