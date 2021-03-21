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

int main()
{
    int lmax = 8;
    float els[5] = {1, 2, 3, 4, 5};
    float azs[5] = {0};
    matrix az = assignMat(sizeof(azs)/sizeof(azs[0]), 1, azs);
    
    
    float *dyntest = malloc(sizeof(float));
    dyntest[0] = 0;
//    float* test = SH2RH(dyntest, 1);
    
    
    gen_delta(dyntest, 1, az, 8);
    

    printf("----\n");
    matrix test = eval_SH(2, dyntest, 1, az);
    printMat(test);
    
    
    return 0;
}

/*
 
 Calculate the rotational harmonic decomposition up to harmonic order "lmax"
 
 */

float* SH2RH(float* SH, size_t size_sh)
{
    
    int lmax = 8;
    
    float *one = malloc(sizeof(float));
    one[0] = 0;
    float azs[1] = {0};
    matrix az = assignMat(sizeof(azs)/sizeof(azs[0]), 1, azs);
    
    matrix* D_SH = gen_delta(one, 0, az, lmax);
    
    for (int i = 0; i < 5; i++)
    {
        printf("i: %d\n", i);
        printMat(D_SH[i]);
    }
    
    // allocate memory
    float* RH = malloc(5*sizeof(float));
    
    return RH;
}

/*
 
 generate the SH coefficients for a delta function pointing along [el az] up to 'lmax'
 
 */

matrix* gen_delta(float el[], size_t size_e, matrix az, int lmax)
{
    // allocate memory
    matrix *SH = malloc(lmax/2+1*sizeof(float));
    
    int counter = 0;
    for (int i = 0; i < lmax+1; i=i+2)
    {
        matrix temp = eval_SH(i, el, size_e, az);
        printf("i: %d\n", i);
        printMat(temp);
        SH[counter] = temp;
        counter ++;
    }
    return SH;
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
