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

matrix* csdeconv(float* R_RH, matrix* DW_SH, matrix* HR_SH, matrix* S, float lambda, float tau, matrix* F_SH)
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
    int counter = 0;
    /*
     
     
        building spherical convolution matrix
     
     
     */
    /* --- */ float* m = malloc(s*sizeof(float));
    for (int i = 0; i < lmax+1; i=i+2)
    {
        for (int j = 0; j < 2*i+1; j++)
        {
            m[counter] = R_RH[i/2];
            counter++;
        }
    }

    if (DW_SH->col != s)
    {
        printf("SH.col is NOT equal to R_RH\n");
    }

    /* --- */ float *dfconv = malloc(DW_SH->row*DW_SH->col*sizeof(float));
    counter = 0;
    for (int i = 0; i < DW_SH->row; i++)
    {
        for (int j = 0; j < DW_SH->col; j++)
        {
            dfconv[counter] = m[j];
            counter++;
        }
    }
    
    matrix* rfconv = malloc(sizeof(matrix));
    assignMat(DW_SH->row, DW_SH->col, dfconv, rfconv);
    
    freeMat(rfconv);
    free(dfconv);
    free(m);
    
//    /* FCONV matrix = spherical convultion */
//    matrix* fconv = malloc(sizeof(matrix));
//    eleMultiplyMat(DW_SH, rfconv, fconv);
//
//    free(dfconv);
//    free(m);
//    freeMat(rfconv);
//
//    /*
//
//
//        generate initial FOD estimate
//        truncated at lmax=4
//
//
//     */
//    /* --- */ double* dF_SH1 = malloc(DW_SH->row * nSH_for_lmax(4)*sizeof(double));
//    /* --- */ double* dS1 = malloc(S->row*sizeof(double));
//
//    counter = 0;
//    for (int i = 0; i < DW_SH->row; i++)
//    {
//        for (int j = 0; j < nSH_for_lmax(4); j++)
//        {
//            dF_SH1[counter] = (double) fconv->data[i * DW_SH->col + j];
//            counter++;
//        }
//    }
////    //printf("DW_SH.row: %d\n", DW_SH.row);
////    //printf("S.row: %d\n", S.row);
//    for (int i = 0; i < S->row; i++)
//    {
//        dS1[i] = (double) S->data[i];
//    }
//
//    /* --- */ double*dF_SH2 = malloc(nSH_for_lmax(4) * sizeof(double));
//    dF_SH2 = qr_solve(DW_SH->row, nSH_for_lmax(4), dF_SH1, dS1);
//
//    /* --- */ float *dataF_SH = malloc(HR_SH->col * sizeof(float));
//    for (int i = 0; i < nSH_for_lmax(4); i++)
//    {
//        dataF_SH[i] = (float) dF_SH2[i];
//        // printf("dataF_SH[%d]: %f\n", i, dataF_SH[i]);
//    }
//
//    for (int i = nSH_for_lmax(4); i < HR_SH->col; i++)
//    {
//        dataF_SH[i] = 0;
//    }
//    /* --- */ assignMat(HR_SH->col, 1, dataF_SH, F_SH);
//
//    free(dF_SH1);
//    // free(dS1);
//    free(dF_SH2);
//    free(dataF_SH);
//
//    // set threshold on FOD amplitude used to identify 'negative' values:
//    matrix* multmt = malloc(sizeof(matrix));
//    multiplyMat(HR_SH, F_SH, multmt);
//
//    matrix* mt = malloc(sizeof(matrix));
//    meanMat(multmt, mt);
//    scalarMat(tau, mt);
//    // matrix mt = scalarMat(tau, meanMat(multiplyMat(HR_SH, F_SH)));
//    float threshold = mt->data[0];
//    freeMat(multmt);
//    freeMat(mt);
//
//
//    // scale lambda to account for differences in the number of
//    // DW directions and number of mapped directions
//    lambda = lambda * fconv->row * R_RH[0] / HR_SH->row;
//
//    /*
//
//
//        main iteration loop
//
//
//     */
//
//    /* --- */ float* ndfconv = malloc(fconv->row*HR_SH->col*sizeof(float));
//    for (int i = 0; i < fconv->row; i++)
//    {
//        for (int j = 0; j < HR_SH->col; j++)
//        {
//            if (j >= fconv->col)
//            {
//                ndfconv[i*HR_SH->col+j] = 0;
//            }
//            else
//            {
//                ndfconv[i*HR_SH->col+j] = fconv->data[i*fconv->col+j];
//            }
//        }
//    }
//
//    matrix* nfconv = malloc(sizeof(matrix));
//    assignMat(fconv->row, HR_SH->col, ndfconv, nfconv);
////    /* --- */ matrix nfconv = assignMat(fconv.row, HR_SH.col, ndfconv);
//    free(ndfconv);
//    freeMat(fconv);
//
//
//    // /* --- */ float* dk = malloc(sizeof(float));
//    /* --- */ matrix* k = malloc(sizeof(matrix));
//    /* --- */ matrix* A = malloc(sizeof(matrix));
//    /* --- */ matrix* k2 = malloc(sizeof(matrix));
//////    free(dk);
//
//    for (int i = 0; i < 50; i++)
//    {
//
//
////        /* --- */ matrix A = multiplyMat(HR_SH, F_SH);
//        multiplyMat(HR_SH, F_SH, A);
////        /* --- */ matrix k2 = findMat(A, threshold);
//        findMat(A, threshold, k2);
//
//        if ((k2->row + nfconv->row) < HR_SH->col)
//        {
//            printf("ERROR: Too few negative directions identified - failed to converge\n");
//            freeMat(k);
//            freeMat(A);
//            freeMat(k2);
//            return F_SH;
//        }
//
//        if (k2->row == k->row && k2->col == k->col)
//        {
//            // return F_SH;
//            bool same = true;
//            for (int i = 0; i < k->row*k->col; i++)
//            {
//                if (k->data[i] != k2->data[i])
//                {
//                    same = false;
//                }
//            }
//            if (same == true)
//            {
//                // printf("same was true\n");
//                freeMat(k);
//                freeMat(A);
//                freeMat(k2);
//                return F_SH;
//            }
//        }
//        assignMat(k2->row, k2->col, k2->data, k);
//        //printMat(k);
//        //printf("k.row: %d\n", k.row);
//        //printf("k.col: %d\n", k.col);
//
//        int counter = 0;
//        int counter2 = 0;
//        /* --- */ float* dHRSH = malloc(k->row*HR_SH->col*sizeof(float));
//        for (int i = 0; i < k->row; i++)
//        {
//            if (k->data[counter] == i)
//            {
//                for (int j = 0; j < HR_SH->col; j++)
//                {
//                    dHRSH[counter2] = HR_SH->data[i*HR_SH->col+j];
//                    counter2++;
//                }
//                counter++;
//            }
//        }
////        /* --- */ matrix second = scalarMat(lambda, assignMat(k.row, HR_SH.col, dHRSH));
//
//        matrix* second = malloc(sizeof(matrix));
//        assignMat(k->row, HR_SH->col, dHRSH, second);
//        scalarMat(lambda, second);
//        free(dHRSH);
//
//        /* --- */ float* dM = malloc((second->row+nfconv->row)*second->col*sizeof(float));
//        /* --- */ double* ddM = malloc((second->row+nfconv->row)*second->col*sizeof(double));
//
//        int dMcounter = 0;
//        int fcounter = 0;
//        int secondcounter = 0;
//        for (int i = 0; i < second->row+nfconv->row; i++)
//        {
//            for (int j = 0; j< second->col; j++)
//            {
//                if (i < nfconv->row)
//                {
//                    dM[dMcounter] = nfconv->data[fcounter];
//                    ddM[dMcounter] = nfconv->data[fcounter];
//                    fcounter++;
//                }
//                else
//                {
//                    dM[dMcounter] = second->data[secondcounter];
//                    ddM[dMcounter] = second->data[secondcounter];
//                    secondcounter++;
//                }
//                dMcounter++;
//            }
//        }
//
////        /* --- */ matrix M = assignMat((second.row+nfconv.row), second.col, dM);
//
//        matrix* M = malloc(sizeof(matrix));
//        assignMat((second->row+nfconv->row), second->col, dM, M);
//
//        float* dfF_SH = malloc(M->row*M->col*sizeof(float));
//
//
//        /* --- */ double* dS2 = malloc((S->row+k->row)*sizeof(double));
//
//        for (int i = 0; i < S->row+k->row; i++)
//        {
//            if (i < S->row)
//            {
//                dS2[i] = dS1[i];
//            } else
//            {
//                dS2[i] = 0;
//            }
//            // printf("dS2[%d]: %f\n", i, dS2[i]);
//        }
////        // printf("S.row: %d\n", S.row);
////        // printf("k.row: %d\n", k.row);
//
//        /* --- */ double*dF_SH2 = malloc(M->col*M->row*sizeof(double));
//
//        dF_SH2 = qr_solve(M->row, M->col, ddM, dS1);
//
//        /* --- */ float* ddfF_SH = malloc(M->col*sizeof(float));
//        for (int i = 0; i < M->col; i++)
//        {
//            ddfF_SH[i] = dF_SH2[i];
//            // printf("ddfF_SH[%d]: %f\n", i, ddfF_SH[i]);
//        }
//        assignMat(M->col, 1, ddfF_SH, F_SH);
//        free(dM);
//        free(ddM);
//        free(dfF_SH);
//        free(dS2);
//        free(dF_SH2);
//        // printf("F_SH.row: %d\n", F_SH.row);
//        // printf("F_SH.col: %d\n", F_SH.col);
//    }
////
////    printf("maximum number of iterations exceeded - FAILED TO CONVERGE\n");

    return F_SH;
    
}

matrix* SH2amp(matrix* csd, float* V, matrix* amp)
{
    
    matrix* temp_csd = malloc(sizeof(matrix));
    float* dtemp_csd = malloc(sizeof(float)*csd->row*csd->col);
    /* takes n x csd->row */
    int counter = 0;
    for (int i = 0; i < csd->row; i++)
    {
        for (int j = 0; j < csd->col; j++)
        {
            dtemp_csd[counter] = temp_csd->data[i*csd->col + j];
            counter++;
        }
    }
    assignMat(csd->row, csd->col, dtemp_csd, temp_csd);
    
    multiplyMat(temp_csd, csd, amp);
    
    return amp;
}

/* returns number of even SH coefficients for 'lmax'*/

int nSH_for_lmax(int lmax)
{
    return (lmax+1)*(lmax+2)/2;
}

