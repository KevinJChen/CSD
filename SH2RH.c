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
    float els[5] = {1, 2, 3, 4, 5};
    float azs[5] = {1, 2, 3, 4, 5};
    matrix az = assignMat(sizeof(azs)/sizeof(azs[0]), 1, azs);
    matrix test = eval_SH(3, els, sizeof(els)/sizeof(els[0]), az);

    return 0;
}


/*
 
 evaluates the lth order spherical harmonic coefficient at positions [ el az ]
 
 */

matrix eval_SH(int l, float el[], size_t size_e, matrix az)
{
    
    // memory allocation
    float * azs1 = malloc(l * az.row * sizeof(float));
    float * azs2 = malloc(l * az.row * sizeof(float));
    float * ones = malloc(az.row * sizeof(float));
    
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
        for (int j = 0; j < l; j++)
        {
            azs2[counter] = j * az.data[i];
            counter++;
        }
    }
                
    // ones(size(az,1),1)
    for (int i = 0; i < az.row; i++)
    {
        ones[i] = 1;
    }
                
    matrix m = assignMat(1, l, azs1);
    printMat(m);
    printf("\n");
    matrix s = multiplyMat(az, m);
    printMat(s);
    printf("\n");
    matrix m1 = assignMat(1, l, azs2);
    matrix s1 = multiplyMat(az, m1);
    printMat(s1);

    
    if (l != 0)
    {
        
    }
    
    
    return s;
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
