//
//  csd.c
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
#include "qr_solve.h"

#include "csd.h"

matrix csdeconv(float* R_RH, matrix DW_SH, matrix HR_SH, matrix S, float lambda, float tau)
{
    // lambda default is 1
    // tau default is 0.1
    
    // appropriate lmax value
    int lmax = 8;
    
    int s = 0;
    for (int i = 0; i < lmax+1; i=i+2)
    {
        s=s+2*i+1;
    }
    float* m = malloc(s*sizeof(float));
    int counter = 0;
    // buliding spherical convolution matrix
    for (int i = 0; i < lmax+1; i=i+2)
    {
        for (int j = 0; j < 2*i+1; j++)
        {
            m[counter] = R_RH[i/2];
            counter++;
        }
    }
    
    if (DW_SH.col != s)
    {
        printf("SH.col is NOT equal to R_RH");
    }
    
    float *dfconv = malloc(DW_SH.row*DW_SH.col*sizeof(float));
    counter = 0;
    for (int i = 0; i < DW_SH.row; i++)
    {
        for (int j = 0; j < DW_SH.col; j++)
        {
            dfconv[counter] = m[j];
            counter++;
        }
    }
    matrix rfconv = assignMat(DW_SH.row, DW_SH.col, dfconv);
    matrix fconv = eleMultiplyMat(DW_SH, rfconv);
    
    
    // generate initial FOD estimate
    // truncated at lmax = 4
    double* dF_SH1 = malloc(DW_SH.row * nSH_for_lmax(4)*sizeof(double));
    double* dS1 = malloc(S.row*sizeof(double));

    counter = 0;
    for (int i = 0; i < DW_SH.row; i++)
    {
        for (int j = 0; j < nSH_for_lmax(4); j++)
        {
            dF_SH1[counter] = (double) fconv.data[i * DW_SH.col + j];
            counter++;
        }
    }
    //printf("DW_SH.row: %d\n", DW_SH.row);
    //printf("S.row: %d\n", S.row);
    for (int i = 0; i < S.row; i++)
    {
        dS1[i] = (double) S.data[i];
    }

    double*dF_SH2 = malloc(nSH_for_lmax(4) * sizeof(double));
    dF_SH2 = qr_solve(DW_SH.row, nSH_for_lmax(4), dF_SH1, dS1);

    float *dataF_SH = malloc(HR_SH.col * sizeof(float));
    for (int i = 0; i < nSH_for_lmax(4); i++)
    {
        dataF_SH[i] = (float) dF_SH2[i];
        // printf("dataF_SH[%d]: %f\n", i, dataF_SH[i]);
    }
    
//    dataF_SH[0] = 0.5651;
//    dataF_SH[1] = 0.5454;
//    dataF_SH[2] = -0.0021;
//    dataF_SH[3] = -0.6328;
//    dataF_SH[4] =  -0.0010;
//    dataF_SH[5] =  0.5453;
//    dataF_SH[6] = -0.0219;
//    dataF_SH[7] = -0.0000;
//    dataF_SH[8] = -0.4668;
//    dataF_SH[9] = 0.0036;
//    dataF_SH[10] =  0.6466;
//    dataF_SH[11] = 0.0043;
//    dataF_SH[12] = -0.4768;
//    dataF_SH[13] = 0.0075;
//    dataF_SH[14] = -0.0042;

    for (int i = nSH_for_lmax(4); i < HR_SH.col; i++)
    {
        dataF_SH[i] = 0;
    }
    matrix F_SH = assignMat(HR_SH.col, 1, dataF_SH);

    // set threshold on FOD amplitude used to identify 'negative' values:
    matrix mt = scalarMat(tau, meanMat(multiplyMat(HR_SH, F_SH)));
    float threshold = mt.data[0];

    // scale lambda to account for differences in the number of
    // DW directions and number of mapped directions
    lambda = lambda * fconv.row * R_RH[0] / HR_SH.row;

    /*


        main iteration loop


     */

    float* ndfconv = malloc(fconv.row*HR_SH.col*sizeof(float));
    for (int i = 0; i < fconv.row; i++)
    {
        for (int j = 0; j < HR_SH.col; j++)
        {
            if (j >= fconv.col)
            {
                ndfconv[i*HR_SH.col+j] = 0;
            }
            else
            {
                ndfconv[i*HR_SH.col+j] = fconv.data[i*fconv.col+j];
            }
        }
    }
    matrix nfconv = assignMat(fconv.row, HR_SH.col, ndfconv);

    float* dk = malloc(sizeof(float));
    matrix k = assignMat(1, 1, dk);

    for (int i = 0; i < 50; i++)
    {
        matrix A = multiplyMat(HR_SH, F_SH);
        matrix k2 = findMat(A, threshold);

//        printf("k.row: %d\n", k.row);
//        printf("k.col: %d\n", k.col);
//        printf("k2.row: %d\n", k2.row);
//        printf("k2.col: %d\n", k2.col);
//        printf("-------\n");

        if ((k2.row + nfconv.row) < HR_SH.col)
        {
            printf("ERROR: Too few negative directions identified - failed to converge\n");
            return F_SH;
        }

        if (k2.row == k.row && k2.col == k.col)
        {
            // return F_SH;
            bool same = true;
            for (int i = 0; i < k.row*k.col; i++)
            {
                if (k.data[i] != k2.data[i])
                {
                    same = false;
                }
            }
            if (same == true)
            {
                // printf("same was true\n");
                return F_SH;
            }
        }
        k = assignMat(k2.row, k2.col, k2.data);
        //printMat(k);
        //printf("k.row: %d\n", k.row);
        //printf("k.col: %d\n", k.col);

        int counter = 0;
        int counter2 = 0;
        float* dHRSH = malloc(k.row*HR_SH.col*sizeof(float));
        for (int i = 0; i < HR_SH.row; i++)
        {
            if (k.data[counter] == i)
            {
                for (int j = 0; j < HR_SH.col; j++)
                {
                    dHRSH[counter2] = HR_SH.data[i*HR_SH.col+j];
                    counter2++;
                }
                counter++;
            }
        }
        matrix second = scalarMat(lambda, assignMat(k.row, HR_SH.col, dHRSH));

        float* dM = malloc((second.row+nfconv.row)*second.col*sizeof(float));
        double* ddM = malloc((second.row+nfconv.row)*second.col*sizeof(double));

        int dMcounter = 0;
        int fcounter = 0;
        int secondcounter = 0;
        for (int i = 0; i < second.row+nfconv.row; i++)
        {
            for (int j = 0; j< second.col; j++)
            {
                if (i < nfconv.row)
                {
                    dM[dMcounter] = nfconv.data[fcounter];
                    ddM[dMcounter] = nfconv.data[fcounter];
                    fcounter++;
                }
                else
                {
                    dM[dMcounter] = second.data[secondcounter];
                    ddM[dMcounter] = second.data[secondcounter];
                    secondcounter++;
                }
                dMcounter++;
            }
        }

        matrix M = assignMat((second.row+nfconv.row), second.col, dM);

        // float* dfF_SH = malloc(M.row*M.col*sizeof(float));

        double* dS2 = malloc((S.row+k.row)*sizeof(double));

        for (int i = 0; i < S.row+k.row; i++)
        {
            if (i < S.row)
            {
                dS2[i] = dS1[i];
            } else
            {
                dS2[i] = 0;
            }
            // printf("dS2[%d]: %f\n", i, dS2[i]);
        }
        // printf("S.row: %d\n", S.row);
        // printf("k.row: %d\n", k.row);

        double*dF_SH2 = malloc(M.col*M.row*sizeof(double));

        dF_SH2 = qr_solve(M.row, M.col, ddM, dS1);

        float* ddfF_SH = malloc(M.col*sizeof(float));
        for (int i = 0; i < M.col; i++)
        {
            ddfF_SH[i] = dF_SH2[i];
            // printf("ddfF_SH[%d]: %f\n", i, ddfF_SH[i]);
        }
        F_SH = assignMat(M.col, 1, ddfF_SH);
        // printf("F_SH.row: %d\n", F_SH.row);
        // printf("F_SH.col: %d\n", F_SH.col);
    }

     printf("maximum number of iterations exceeded - FAILED TO CONVERGE\n");
    // matrix F_SH;
    return F_SH;
    
}

/* returns number of even SH coefficients for 'lmax'*/

int nSH_for_lmax(int lmax)
{
    return (lmax+1)*(lmax+2)/2;
}

