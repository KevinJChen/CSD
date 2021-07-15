//
//  SH2RH.c
//  csd
//
//  Created by Kevin Chen.
//  Copyright © 2021 Kevin Chen. All rights reserved.
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "matrix.h"
#include "legendre.h"

#include "SH2RH.h"

/*
 
 Calculate the rotational harmonic decomposition up to harmonic order "lmax"
 
 */

float* SH2RH(float* SH, size_t size_sh)
{
    
    int lmax = 8;
    
    // el = [0]
    float *el = malloc(sizeof(float));
    el[0] = 0;
    
    // az = [0]
    float azs[1] = {0};
    matrix az = assignMat(sizeof(azs)/sizeof(azs[0]), 1, azs);
    
    float* D_SH = gen_delta(el, 1, az, lmax);
    
    int nonzero = 0;
    // count number of nonzero elements in D_SH
    for (int i = 1; i < D_SH[0]+1; i++)
    {
        if (D_SH[i] != 0)
        {
            nonzero++;
        }
    }
    
    // allocate memory
    float* zD_SH = malloc(nonzero*sizeof(float));
    float* zSH = malloc(nonzero*sizeof(float));
    float* RH = malloc(nonzero*sizeof(float));
    
    // find nonzero elements in D_SH
    int counter = 0;
    for (int i = 1; i < D_SH[0]+1; i++)
    {
        if (D_SH[i] != 0)
        {
            zD_SH[counter] = D_SH[i];
            zSH[counter] = SH[i-1];
            counter++;
        }
    }
    
    // right array division
    for (int i = 0; i < nonzero; i++)
    {
        RH[i] = zSH[i] / zD_SH[i];
    }
    
    return RH;
}

/*
 
 generate the SH coefficients for a delta function pointing along [el az] up to 'lmax'
 
 */

float* gen_delta(float el[], size_t size_e, matrix az, int lmax)
{
    int lD_SH = 0;
    // allocate memory
    matrix *SH = malloc(lmax/2+1*sizeof(float));
    
    int counter = 0;
    for (int i = 0; i < lmax+1; i=i+2)
    {
        matrix temp = eval_SH(i, el, size_e, az);
        lD_SH = lD_SH + temp.row*temp.col;
        SH[counter] = temp;
        counter++;
    }
    
    float* D_SH = malloc(lD_SH*sizeof(float));
    
    counter = 1;
    for (int i = 0; i < lmax/2+1; i++)
    {
        for (int j = 0; j < SH[i].row*SH[i].col; j++)
        {
            D_SH[counter] = SH[i].data[j];
            counter++;
        }
    }
    
    // first element = size of float array
    D_SH[0] = counter-1;
    
    return D_SH;
}



/*
 
 evaluates the lth order spherical harmonic coefficient at positions [ el az ]
 
 */

matrix eval_SH(int l, float el[], size_t size_e, matrix az)
{
    
    // memory allocation
    float * azs1 = malloc(l * az.row * sizeof(float));
    float * azs2 = malloc(l * az.row * sizeof(float));
    float * o = malloc(az.row * sizeof(float));
    float * ss = malloc(az.row * (2*l+1) * sizeof(float));
    
    // az*(l:-1:1)
    int counter = 0;
    for (int i = 0; i < az.row; i++)
    {
        for (int j = l; j > 0; j--)
        {
            azs1[counter] = j * az.data[i];
            counter++;
        }
    }
    
    // az*(1:l)
    counter = 0;
    for (int i = 0 ; i < az.row; i++)
    {
        for (int j = 1; j < l + 1; j++)
        {
            azs2[counter] = j * az.data[i];
            counter++;
        }
    }
                
    // ones(size(az,1),1)
    for (int i = 0; i < az.row; i++)
    {
        o[i] = 1;
    }
                
    // sqrt(2* sin(az * l:-1:1))
    matrix paz = scalarMat(sqrt(2), sinMat(assignMat(az.row, l, azs1)));
    // s
    matrix ones = assignMat(az.row, 1, o);
    // sqrt(2) * cos(az*(1:l))
    matrix aaz = scalarMat(sqrt(2), cosMat(assignMat(az.row, l, azs2)));
    
    matrix s = assignMat(az.row, 1, o);
    
    if (l != 0)
    {
        counter = 0;
        for (int i = 0; i < az.row; i++)
        {
            for (int j = 0; j < paz.col; j++)
            {
                ss[counter] = paz.data[i*paz.col+j];
                counter++;
            }
            for (int j = 0; j < ones.col; j++)
            {
                ss[counter] = ones.data[i*ones.col+j];
                counter++;
                
            }
            for (int j = 0; j < aaz.col; j++)
            {
                ss[counter] = aaz.data[i*aaz.col+j];
                counter++;
            }
        }
        
        // s = [paz ones aaz]
        s = assignMat(az.row, 2*l+1, ss);
    }
    
    // evaluate legendre polynomials
    matrix ev = eval_ALP(l, el, size_e);
    // spherical harmonics
    matrix SH = eleMultiplyMat(ev, transposeMat(s));
    
    return SH;
}


/*
 
 evaluates the Associated Legendre Polynommial at elevations 'el' for harmonic order l
 
 */
matrix eval_ALP(int l, float el[], size_t size)
{
    matrix els = assignMat(size, 1, el);
    
    // cos (el)
    float * cosEls = cosEl(el, size);
    
    // cosEls to double
    double * dEls = malloc(size * sizeof(double));
    for (int i = 0; i < size; i++)
    {
        dEls[i] = (double) cosEls[i];
    }
    
    // legendre(l, els)
    matrix leg = legendre(l, dEls, size);
    
    
    for (int m = 0; m <= l; m++)
    {
        // sqrt[ {2*l+1 * (l-m)!} / {4pi * (l+m)!} ]
        float mult = sqrt( ((2*l+1)*factorial(l-m)) / ((4*M_PI) * (factorial(l+m))) );
        // printf("mult %f\n", mult);
        for (int i = 0; i < leg.col; i++)
        {
            // printf("index %d\n", m+1*leg.col + i);
            leg.data[m*leg.col + i] = mult*leg.data[m*leg.col + i];
        }
    }
    
    if (l != 0)
    {
        
        float * leg1data = malloc((leg.row*2-1)*leg.col * sizeof(float));
        
        
        // flipped first section
        for (int i = 0; i < leg.row-1; i++)
        {
            for (int j = 0; j < leg.col; j++)
            {
                leg1data[i*leg.col + j] = leg.data[(leg.row-i-1) * leg.col + j];
            }
        }
        
        // normal second section
        for (int i = leg.row-1; i < leg.row*2-1; i++)
        {
            for (int j = 0; j < leg.col; j++)
            {
                leg1data[i*leg.col+j] = leg.data[(i-(leg.row-1)) * leg.col + j];
            }
        }
        
        matrix leg1 = assignMat(leg.row*2-1, leg.col, leg1data);

        return leg1;
    }
    
    return leg;
}


matrix amp2SH(float* S, size_t size)
{
    // float* SH = malloc(sizeof(float)*size);
    
    matrix mS = assignMat(size, 1, S);
    matrix mI = mpinverseMat(mS);
    
    matrix SH = multiplyMat(mI, mS);
    
    return SH;
}

float * dir3002SH(float * dir300, size_t size)
{
    
    int lmax = 8;
    
    /* copy of dir300 */
    float * cdir300 = malloc(sizeof(float) * size);
    for (int i = 0; i < size; i++)
    {
        cdir300[i] = dir300[i];
    }
    
    /* square everything */
    for (int i = 0; i < size; i++)
    {
        float temp = dir300[i] * dir300[i];
        dir300[i] = temp;
    }
    
    
    float * n = malloc(sizeof(float) * size);
    
    /* sum of rows */
    int counter = 0;
    for (int i = 0; i < size; i=i+3)
    {
    
        n[counter] = dir300[i] + dir300[i+1] + dir300[i+2];
        counter++;
    }
    
    
    int k_size = 0;
    /* square root */
    for (int i = 0; i < size/3; i++)
    {
        n[i] = sqrt(n[i]);
        if (n[i] != 0)
        {
            k_size++;
        }
    }
    if (k_size != size/3)
    {
        printf("zero elements found");
    }
    
    counter = 0;
    for (int i = size/3; i < size/3*2; i++)
    {
        n[i] = n[counter];
        n[i+size/3] = n[counter];
    }
    
    /* cdir300 ./ n (element-wise division) */
    float *P = malloc(sizeof(float) * size);
    for (int i = 0; i < size; i++)
    {
        P[i] = cdir300[i] / n[i];
    }
    
    /* cartesian to spherical */
    for (int i = 0; i < size; i=i+3)
    {
        float x = P[i];
        float y = P[i+1];
        float z = P[i+2];

        P[i] = atan2(sqrt(x*x + y*y), z);
        P[i+1] = atan2(y, x);
        P[i+2] = sqrt(x*x + y*y + z*z);
    }
    
    float *el = malloc(sizeof(float) * size/3);
    for (int i = 0; i < size; i=i+3)
    {
        el[i] = P[i];
    }
    float *daz = malloc(sizeof(float) * size/3);
    for (int i = 1; i < size; i=i+3)
    {
        daz[i] = P[i];
    }
    matrix az = assignMat(size/3, 1, daz);
    
    
    for (int i = 0; i < lmax; i=i+2)
    {
        matrix add_temp = eval_SH(i, el, size/3, az);
    }
    matrix add_temp = eval_SH(0, el, size/3, az);
    printMat(add_temp);

    return P;


}


/*
 
 HELPER FUNCTIONS
 
 */

/* get cosine of each element in float array */
float * cosEl(float el[], size_t size)
{
    float * els = malloc(size * sizeof(float));
    for (int i = 0; i < size; i++)
    {
        els[i] = cos(el[i]);
    }
    return els;
}

/* factorial */
float factorial(int num)
{
    if (num == 0 || num == 1)
    {
        return 1;
    }
    else if (num < 0)
    {
        printf("ERROR: Can not take factorial of negative number");
    }
    else
    {
        return num * factorial(num-1);
    }
    return 0;
}
