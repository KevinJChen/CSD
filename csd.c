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
    
    float *dF_SH = malloc(DW_SH.row*nSH_for_lmax(4)*sizeof(float));
    double* dF_SH1 = malloc(DW_SH.row * nSH_for_lmax(4) * sizeof(double));
    double* dS1 = malloc(DW_SH.row*sizeof(double));
    
    counter = 0;
    for (int i = 0; i < DW_SH.row; i++)
    {
        for (int j = 0; j < nSH_for_lmax(4); j++)
        {
            dF_SH[counter] = fconv.data[counter];
            dF_SH1[counter] = fconv.data[counter];
            counter++;
        }
    }
    
    double*dF_SH2 = malloc(nSH_for_lmax(4) * sizeof(double));
    dF_SH2 = qr_solve(DW_SH.row, nSH_for_lmax(4), dF_SH1, dS1);
    float *dataF_SH = malloc(HR_SH.col * sizeof(float));
    for (int i = 0; i < nSH_for_lmax(4); i++)
    {
        dataF_SH[i] = dF_SH2[i];
    }
    for (int i = nSH_for_lmax(4); i < HR_SH.col; i++)
    {
        dataF_SH[i] = 0;
    }
    matrix F_SH = assignMat(HR_SH.col, 1, dataF_SH);
    
    
    // set threshold on FOD amplitude used to identify 'negative' values:
    matrix threshold = scalarMat(tau, meanMat(multiplyMat(F_SH, HR_SH)));
    
    // scale lambda to account for differences in the number of
    // DW directions and number of mapped directions
    lambda = lambda * fconv.row * R_RH[0] / HR_SH.col;
    
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
                printf("HERERERASDLKF\n");
                printf("%d\n", fconv.col);
                // ndfconv[i*HR_SH.col+j] = 0;
            }
            else
            {
                ndfconv[i*HR_SH.col+j] = fconv.data[i*fconv.col+j];
            }
        }
    }
    matrix nfconv = assignMat(fconv.row, HR_SH.col, ndfconv);
    
    matrix k;
    for (int i = 0; i < 50; i++)
    {
        matrix A = multiplyMat(HR_SH, F_SH);
        matrix k2 = findMat(A, threshold);
        
        if ((k2.row + nfconv.row) < HR_SH.col)
        {
            printf("ERROR: Too few negative directions identified - failed to converge\n");
            matrix error;
            return error;
        }
        
        if (k2.row == k.row && k2.col == k.col)
        {
            printf("ERROR: \n");
            matrix error;
            return error;
        }
        k = k2;
        
        int counter = 0;
        int counter2 = 0;
        float* dHRSH = malloc(HR_SH.col*k.row*sizeof(float));
        for (int i = 0; i < HR_SH.row; i++)
        {
            if (k.data[counter] == i + 1)
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
        
        float* dM = malloc((second.row+fconv.row)*second.col*sizeof(float));
        double* ddM = malloc((second.row+fconv.row)*second.col*sizeof(double));
        
        
        int dMcounter = 0;
        int fcounter = 0;
        int secondcounter = 0;
        for (int i = 0; i < second.row+fconv.row; i++)
        {
            for (int j = 0; j< second.col; j++)
            {
                if (i < fconv.row)
                {
                    dM[dMcounter] = fconv.data[fcounter];
                    ddM[dMcounter] = fconv.data[fcounter];
                    fcounter++;
                }
                else
                {
                    dM[dMcounter] = second.data[secondcounter];
                    ddM[dMcounter] = second.data[secondcounter];
                    secondcounter++;
                }
            }
        }
        
        matrix M = assignMat((second.row+fconv.row), second.col, dM);
        
        dF_SH1 = malloc(M.row*M.col*sizeof(float));
        
        double* dS2 = malloc(S.row*k.row*sizeof(float));
        
        for (int i = 0; i < S.row*k.row; i++)
        {
            if (i < S.row)
            {
                dS1[i] = S.data[i];
            }
            else
            {
                dS1[i] = 0;
            }
        }

        double*dF_SH2 = malloc(M.col*M.row*sizeof(double));
        dF_SH2 = qr_solve(M.row, M.col, ddM, dS1);

        float* ddfF_SH = malloc(M.col*M.row*sizeof(float));
        for (int i = 0; i < M.col*M.row; i++)
        {
            ddfF_SH[i] = dF_SH2[i];
        }
        matrix finalF_SH = assignMat(M.row, M.col, ddfF_SH);
        
        return finalF_SH;

    }
    
    printf("maximum number of iterations exceeded - FAILED TO CONVERGE\n");
    matrix error;
    return error;
    
    
}

/* returns number of even SH coefficients for 'lmax'*/

int nSH_for_lmax(int lmax)
{
    return (lmax+1)*(lmax+2)/2;
}

