
#include "matrices.h"

/*****************************************************************************
 *
 *
 *
 *
 *
 *****************************************************************************/

float *VEC_crossp(float *v1, float *v2, float *result)
{
    
    float tmp[3];
    
    tmp[0] = v1[1]*v2[2]-v2[1]*v1[2];
    tmp[1] = v1[2]*v2[0]-v2[2]*v1[0];
    tmp[2] = v1[0]*v2[1]-v2[0]*v1[1];
    
    result[0] = tmp[0];
    result[1] = tmp[1];
    result[2] = tmp[2];
    
    return result;
    
}

/*****************************************************************************
 *
 *
 *
 *
 *
 *****************************************************************************/

float *VEC_make_unit(float *v) 
{
    float norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if (norm != 0) {
        v[0] /= norm;
        v[1] /= norm;
        v[2] /= norm;
    }
    
    return v;
}

/*****************************************************************************
 *
 *
 *
 *
 *
 *****************************************************************************/

SYMMAT33 *MAT33_to_SYMMAT33(MAT33 *src, SYMMAT33 *dest)
{
    dest->data[0] = src->data[0];
    dest->data[1] = src->data[4];
    dest->data[2] = src->data[8];
    dest->data[3] = src->data[1];
    dest->data[4] = src->data[2];
    dest->data[5] = src->data[5];
    
    return dest;
    
}

/*****************************************************************************
 *
 *
 *
 *
 *
 *****************************************************************************/

MAT33 *MAT33_mult(MAT33 *A, MAT33 *B, MAT33 *C)
{
    
    MAT33 tmp;
    float *t = tmp.data;
    
    float *a = A->data;
    float *b = B->data;
    float *c = C->data;
    
    t[0] = a[0]*b[0] + a[1]*b[3] + a[2]*b[6];
    t[1] = a[0]*b[1] + a[1]*b[4] + a[2]*b[7];
    t[2] = a[0]*b[2] + a[1]*b[5] + a[2]*b[8];
    
    t[3] = a[3]*b[0] + a[4]*b[3] + a[5]*b[6];
    t[4] = a[3]*b[1] + a[4]*b[4] + a[5]*b[7];
    t[5] = a[3]*b[2] + a[4]*b[5] + a[5]*b[8];
    
    t[6] = a[6]*b[0] + a[7]*b[3] + a[8]*b[6];
    t[7] = a[6]*b[1] + a[7]*b[4] + a[8]*b[7];
    t[8] = a[6]*b[2] + a[7]*b[5] + a[8]*b[8];
    
    c[0] = t[0] ; c[1] = t[1]; c[2] = t[2];
    c[3] = t[3] ; c[4] = t[4]; c[5] = t[5];
    c[6] = t[6] ; c[7] = t[7]; c[8] = t[8];
    
    return C;
    
}

/*****************************************************************************
 *
 *
 *
 *
 *
 *****************************************************************************/

SYMMAT33 *SYMMAT3_mult(SYMMAT33 *A, SYMMAT33 *B, SYMMAT33 *C)
{
    
    return C;
    
}

/*****************************************************************************
 *
 *
 *
 *
 *
 *****************************************************************************/

SYMMAT33 *SYMMAT33_make(float a00, float a01, float a02,
                        float a10, float a11, float a12,
                        float a20, float a21, float a22)
{
    
    SYMMAT33 *m = malloc(sizeof(SYMMAT33));
    /*
     a00 a01 a02
     a01 a11 a12
     a02 a12 a22
     */
    
    m->data[0] = a00; 
    m->data[1] = a11;
    m->data[2] = a22;
    m->data[3] = a01;
    m->data[4] = a02;
    m->data[5] = a12;
    
    return m;
    
}

/*****************************************************************************
 *
 *
 *
 *
 *
 *****************************************************************************/

MAT33 *MAT33_make(float a00, float a01, float a02,
                  float a10, float a11, float a12,
                  float a20, float a21, float a22)
{
    
    MAT33 *m = malloc(sizeof(MAT33));
    
    MAT33_assign(a00, a01, a02,
                 a10, a11, a12,
                 a20, a21, a22, m);
    
    return m;
    
}

/*****************************************************************************
 *
 *
 *
 *
 *
 *****************************************************************************/

void MAT33_assign(float a00, float a01, float a02,
                  float a10, float a11, float a12,
                  float a20, float a21, float a22, MAT33 *dest)
{
    
    dest->data[0] = a00; dest->data[1] = a01; dest->data[2] = a02;
    dest->data[3] = a10; dest->data[4] = a11; dest->data[5] = a12;
    dest->data[6] = a20; dest->data[7] = a21; dest->data[8] = a22;
    
}