#include "arpack_types.h"

double ar_d_sign(const double* a, const double* b)
{
    double x = fabs(*a);
    return *b >= 0 ? x : -x;
}

double ar_r_sign(const float* a, const float* b)
{
    double x = (*a >= 0 ? *a : -*a);
    return *b >= 0 ? x : -x;
}

void ar_r_cnjg(a_fcomplex* r, const a_fcomplex* z)
{
    float zi = z->i;
    r->r = z->r;
    r->i = -zi;
}

void ar_d_cnjg(a_dcomplex* r, const a_dcomplex* z)
{
    double zi = z->i;
    r->r = z->r;
    r->i = -zi;
}

void ar_c_div(a_fcomplex* c, const a_fcomplex* a, const a_fcomplex* b)
{
    double ratio, den;
    double abr, abi, cr;

    if ((abr = b->r) < 0.)
        abr = -abr;
    if ((abi = b->i) < 0.)
        abi = -abi;
    if (abr <= abi)
    {
        if (abi == 0)
        {
            float af, bf;
            af = bf = abr;
            if (a->i != 0 || a->r != 0)
                af = 1.;
            c->i = c->r = af / bf;
            return;
        }
        ratio = (double)b->r / b->i;
        den = b->i * (1 + ratio * ratio);
        cr = (a->r * ratio + a->i) / den;
        c->i = (a->i * ratio - a->r) / den;
    }

    else
    {
        ratio = (double)b->i / b->r;
        den = b->r * (1 + ratio * ratio);
        cr = (a->r + a->i * ratio) / den;
        c->i = (a->i - a->r * ratio) / den;
    }
    c->r = cr;
}

void ar_z_div(a_dcomplex* c, const a_dcomplex* a, const a_dcomplex* b)
{
    double ratio, den;
    double abr, abi, cr;

    if ((abr = b->r) < 0.)
        abr = -abr;
    if ((abi = b->i) < 0.)
        abi = -abi;
    if (abr <= abi)
    {
        if (abi == 0) {
            if (a->i != 0 || a->r != 0)
                abi = 1.;
            c->i = c->r = abi / abr;
            return;
        }
        ratio = b->r / b->i;
        den = b->i * (1 + ratio * ratio);
        cr = (a->r * ratio + a->i) / den;
        c->i = (a->i * ratio - a->r) / den;
    }

    else
    {
        ratio = b->i / b->r;
        den = b->r * (1 + ratio * ratio);
        cr = (a->r + a->i * ratio) / den;
        c->i = (a->i - a->r * ratio) / den;
    }
    c->r = cr;
}
