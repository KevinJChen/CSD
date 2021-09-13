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

//float* SH2RH(float* SH, size_t size_sh)
//{
//
//    int lmax = 8;
//
//    // el = [0]
//    float *el = malloc(sizeof(float));
//    el[0] = 0;
//
//    // az = [0]
//    float azs[1] = {0};
//    matrix az = assignMat(sizeof(azs)/sizeof(azs[0]), 1, azs);
//
//    float* D_SH = gen_delta(el, 1, az, lmax);
//
//    int nonzero = 0;
//    // count number of nonzero elements in D_SH
//    for (int i = 1; i < D_SH[0]+1; i++)
//    {
//        if (D_SH[i] != 0)
//        {
//            nonzero++;
//        }
//    }
//
//    // allocate memory
//    float* zD_SH = malloc(nonzero*sizeof(float));
//    float* zSH = malloc(nonzero*sizeof(float));
//    float* RH = malloc(nonzero*sizeof(float));
//
//    // find nonzero elements in D_SH
//    int counter = 0;
//    for (int i = 1; i < D_SH[0]+1; i++)
//    {
//        if (D_SH[i] != 0)
//        {
//            zD_SH[counter] = D_SH[i];
//            zSH[counter] = SH[i-1];
//            counter++;
//        }
//    }
//
//    // right array division
//    for (int i = 0; i < nonzero; i++)
//    {
//        RH[i] = zSH[i] / zD_SH[i];
//    }
//
//    return RH;
//}

/*

 generate the SH coefficients for a delta function pointing along [el az] up to 'lmax'

 */

float* gen_delta(float* el, size_t size_e, matrix* az, int lmax)
{
    
    return el;
//    int lD_SH = 0;
//    float ** SH = malloc((lmax/2+1)*(lmax+1)sizeof(float));
//
//    int counter = 0;
//    matrix* temp = malloc(sizeof(float));
//    for (int i = 0; i < lmax+1; i=i+2)
//    {
//        eval_SH(i, el, size_e, az, temp);
//        float* tempdata = malloc(sizeof(float)*i);
//
//        for (int i = 0; i < temp->row)
//
//        lD_SH = lD_SH + temp->row*temp->col;
//        SH[counter] = temp;
//        counter++;
//    }
//
//
//
//
//
//
//
//
//
//    int lD_SH = 0;
//    // allocate memory
//    matrix *SH = malloc(lmax/2+1*sizeof(float));
//
//    int counter = 0;
//    matrix* temp = malloc(sizeof(matrix));
//    for (int i = 0; i < lmax+1; i=i+2)
//    {
//        eval_SH(i, el, size_e, az, temp);
//        lD_SH = lD_SH + temp->row*temp->col;
//        SH[counter] = temp;
//        counter++;
//    }
//
//    float* D_SH = malloc(lD_SH*sizeof(float));
//
//    counter = 1;
//    for (int i = 0; i < lmax/2+1; i++)
//    {
//        for (int j = 0; j < SH[i].row*SH[i].col; j++)
//        {
//            D_SH[counter] = SH[i]->data[j];
//            counter++;
//        }
//    }
//
//    // first element = size of float array
//    D_SH[0] = counter-1;
//
//    return D_SH;
}



/*

 evaluates the lth order spherical harmonic coefficient at positions [ el az ]

 */

matrix* eval_SH(int l, float* el, size_t size_e, matrix* az, matrix* SH)
{

    // memory allocation
    float * azs1 = malloc(l * az->row * sizeof(float));
    float * azs2 = malloc(l * az->row * sizeof(float));
    float * o = malloc(az->row * sizeof(float));
    float * ss = malloc(az->row * (2*l+1) * sizeof(float));

    // az*(l:-1:1)
    int counter = 0;
    for (int i = 0; i < az->row; i++)
    {
        for (int j = l; j > 0; j--)
        {
            azs1[counter] = j * az->data[i];
            counter++;
        }
    }

    // az*(1:l)
    counter = 0;
    for (int i = 0 ; i < az->row; i++)
    {
        for (int j = 1; j < l + 1; j++)
        {
            azs2[counter] = j * az->data[i];
            counter++;
        }
    }

    // ones(size(az,1),1)
    for (int i = 0; i < az->row; i++)
    {
        o[i] = 1;
    }

    // sqrt(2* sin(az * l:-1:1))
    matrix* paz = malloc(sizeof(matrix));
    assignMat(az->row, l, azs1, paz);
    matrix* scapaz = malloc(sizeof(matrix));
    sinMat(paz, scapaz);
    scalarMat(sqrt(2), scapaz);
        
    // s
    matrix* ones = malloc(sizeof(matrix));
    assignMat(az->row, 1, o, ones);
    
    // sqrt(2) * cos(az*(1:l))
    matrix* aaz = malloc(sizeof(matrix));
    assignMat(az->row, l, azs2, aaz);
    matrix* scaaaz = malloc(sizeof(matrix));
    cosMat(aaz, scaaaz);
    scalarMat(sqrt(2), scaaaz);
    
    matrix* s = malloc(sizeof(matrix));
    assignMat(az->row, 1, o, s);
    //printMat(s);
    
    if (l != 0)
    {
        counter = 0;
        for (int i = 0; i < az->row; i++)
        {
            for (int j = 0; j < scapaz->col; j++)
            {
                ss[counter] = scapaz->data[i*scapaz->col+j];
                //printf("p-ss[%d]: %f\n", counter, ss[counter]);
                counter++;
            }
            for (int j = 0; j < ones->col; j++)
            {
                ss[counter] = ones->data[i*ones->col+j];
                //printf("o-ss[%d]: %f\n", counter, ss[counter]);
                counter++;

            }
            for (int j = 0; j < scaaaz->col; j++)
            {
                ss[counter] = scaaaz->data[i*scaaaz->col+j];
                //printf("c-ss[%d]: %f\n", counter, ss[counter]);
                counter++;
            }
        }
        
        // s = [paz ones aaz]
        assignMat(az->row, 2*l+1, ss, s);
//        printf("l: %d\n", l);
//        printf("scapaz->col: %d\n", scapaz->col);
//        printf("ones->col: %d\n", ones->col);
//        printf("scaaaz->col: %d\n", scaaaz->col);
        //printf("counter: %d\n", counter);
//        for (int i = 0; i < az->row * (2*l+1); i++)
//        {
//            printf("ss[%d] = %f\n", i, ss[i]);
//        }
             
    }
    
    // evaluate legendre polynomials
    matrix* ev = malloc(sizeof(matrix));
    eval_ALP(l, el, size_e, ev);
    // spherical harmonics
    matrix* trans = malloc(sizeof(matrix));
    transposeMat(s, trans);
    
    eleMultiplyMat(ev, trans, SH);
    
    freeMat(aaz);
    freeMat(paz);
    freeMat(ev);
    freeMat(trans);
    freeMat(scapaz);
    freeMat(scaaaz);

    return SH;
}


/*

 evaluates the Associated Legendre Polynommial at elevations 'el' for harmonic order l

 */
matrix* eval_ALP(int l, float* el, size_t size, matrix* leg)
{
    matrix* els = malloc(sizeof(matrix));
    assignMat(size, 1, el, els);

    // cos (el)
    float* cosEls = malloc(sizeof(float) * size);
    cosEls = cosEl(el, size);

    // cosEls to double
    double * dEls = malloc(size * sizeof(double));
    for (int i = 0; i < size; i++)
    {
        dEls[i] = (double) cosEls[i];
    }

    // legendre(l, els)
    legendre(l, dEls, size, leg);

    for (int m = 0; m <= l; m++)
    {
        // sqrt[ {2*l+1 * (l-m)!} / {4pi * (l+m)!} ]
        float mult = sqrt( ((2*l+1)*factorial(l-m)) / ((4*M_PI) * (factorial(l+m))) );
        // printf("mult %f\n", mult);
        for (int i = 0; i < leg->col; i++)
        {
            // printf("index %d\n", m+1*leg.col + i);
            leg->data[m*leg->col + i] = mult*leg->data[m*leg->col + i];
        }
    }

    
    if (l != 0)
    {
        float * leg1data = malloc((leg->row*2-1)*leg->col * sizeof(float));


        // flipped first section
        for (int i = 0; i < leg->row-1; i++)
        {
            for (int j = 0; j < leg->col; j++)
            {
                leg1data[i*leg->col + j] = leg->data[(leg->row-i-1) * leg->col + j];
            }
        }

        // normal second section
        for (int i = leg->row-1; i < leg->row*2-1; i++)
        {
            for (int j = 0; j < leg->col; j++)
            {
                leg1data[i*leg->col+j] = leg->data[(i-(leg->row-1)) * leg->col + j];
            }
        }

        assignMat(leg->row*2-1, leg->col, leg1data, leg);

        free(leg1data);
    }
    
    freeMat(els);
    // free(cosEls);
    free(dEls);

    return leg;
}


matrix* amp2SH(float* S, size_t size, matrix* SH)
{
    // float* SH = malloc(sizeof(float)*size);

    matrix* mS = malloc(sizeof(matrix));
    assignMat(size, 1, S, mS);
    matrix* mI = malloc(sizeof(matrix));
    mpinverseMat(mS, mI);
    multiplyMat(mI, mS, SH);
    
    freeMat(mS);
    freeMat(mI);
    
    return SH;
}

matrix* dir3002SH(float * dir300, size_t size, matrix* SH)
{

    /* higher res SH */
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
        printf("zero elements found\n");
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
        el[i/3] = P[i];
    }
    float *daz = malloc(sizeof(float) * size/3);
    for (int i = 1; i < size; i=i+3)
    {
        daz[i/3] = P[i];
    }
    

    matrix* az = malloc(sizeof(matrix));
    assignMat(size/3, 1, daz, az);
    
    matrix* transaz = malloc(sizeof(float));
    eval_SH(0, el, size/3, az, transaz);
    
    transposeMat(transaz, SH);
    
    free(cdir300);
    free(n);
    free(P);
    freeMat(transaz);

   //  matrix SH = transposeMat(eval_SH(0, el, size/3, az));
    matrix* add_temp = malloc(sizeof(matrix));
    matrix* trans_temp = malloc(sizeof(matrix));
    for (int i = 2; i < lmax+1; i=i+2)
    {
        
        eval_SH(i, el, size/3, az, trans_temp);
        transposeMat(trans_temp, add_temp);
        
        // matrix add_temp = transposeMat(eval_SH(i, el, size/3, az));
        conjoinMat(SH, add_temp, SH);
    }
    
    freeMat(add_temp);
    freeMat(trans_temp);
    free(el);
    free(daz);
    freeMat(az);

    return SH;
}

matrix* bvec2SH(float* bvec, size_t size, matrix* SH)
{

    int lmax = 8;

    /* copy of bvec */
    float *cbvec = malloc(sizeof(float) * size);
    for (int i = 0; i < size; i++)
    {
        cbvec[i] = bvec[i];
    }

    /* square everything */
    for (int i = 0; i < size; i++)
    {
        float temp = bvec[i] * bvec[i];
        bvec[i] = temp;
    }

    float * n = malloc(sizeof(float) * size);

    /* sum of rows */
    for (int i = 0; i < size/3; i++)
    {
        n[i] = bvec[i] + bvec[i+size/3] + bvec[i+size/3*2];
        // printf("n[%d]: %f\n", i, n[i]);
    }

    int k_size = 0;
    /* square root */
    for (int i = 0; i < size/3; i++)
    {
        float temp = sqrt(n[i]);
        n[i] = temp;
        n[i+size/3] = temp;
        n[i+size/3*2] = temp;
        if (n[i] != 0)
        {
            k_size++;
        }
    }
//

    /* */


    /* cbvec ./ n (element-wise division) */
    float *P = malloc(sizeof(float) * k_size*3);
    int counter = 0;
    for (int i = 0; i < size; i++)
    {
        if (cbvec[i] == 0 && n[i] == 0)
        {
            continue;
        }
        P[counter] = cbvec[i] / n[i];
        counter++;
    }

    /* cartesian to spherical */
    for (int i = 0; i < k_size; i++)
    {
        float x = P[i];
        float y = P[i+k_size];
        float z = P[i+k_size*2];

        P[i] = atan2(sqrt(x*x + y*y), z);
        P[i+k_size] = atan2(y, x);
        P[i+k_size*2] = sqrt(x*x + y*y + z*z);
    }

    float *el = malloc(sizeof(float) * k_size);
    for (int i = 0; i < k_size; i++)
    {
        el[i] = P[i];
    }
    float *daz = malloc(sizeof(float) * k_size);
    for (int i = k_size; i < k_size*2; i++)
    {
        daz[i-k_size] = P[i];
    }


    //matrix az = assignMat(k_size, 1, daz);
    matrix* az = malloc(sizeof(matrix));
    assignMat(k_size, 1, daz, az);


    matrix* eval_az = malloc(sizeof(matrix));
    eval_SH(0, el, k_size, az, eval_az);

    transposeMat(eval_az, SH);

    // matrix SH = transposeMat(eval_SH(0, el, k_size, az));
    matrix* add_temp = malloc(sizeof(matrix));
    matrix* eval_temp = malloc(sizeof(matrix));
    for (int i = 2; i < lmax+1; i=i+2)
    {
        // matrix add_temp = transposeMat(eval_SH(i, el, k_size, az));
        eval_SH(i, el, k_size, az, eval_temp);
        // printMat(eval_temp);
        transposeMat(eval_temp, add_temp);
        conjoinMat(SH, add_temp, SH);
    }
    
    freeMat(eval_az);
    free(cbvec);
    free(n);
    free(P);
    free(daz);
    free(el);
    freeMat(az);
    freeMat(add_temp);
    freeMat(eval_temp);
    
    return SH;

}


/*
 
 HELPER FUNCTIONS
 
 */

/* get cosine of each element in float array */
float* cosEl(float *el, size_t size)
{
    for (int i = 0; i < size; i++)
    {
        el[i] = cos(el[i]);
    }
    return el;
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
