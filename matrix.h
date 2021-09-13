//
//  matrix.h
//  csd
//
//  Created by Kevin Chen.
//  Copyright © 2020 Kevin Chen. All rights reserved.
//

#ifndef matrix_h
#define matrix_h

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct matrix
{
    int row;
    int col;
    float* data;
} matrix;

/* creation */
matrix* assignMat(int m, int n, float* data1, matrix* mat1);
//
//
///* operations */
matrix* transposeMat(matrix* mat1, matrix* trans);
matrix* multiplyMat(matrix* mat1, matrix* mat2, matrix* prod);
float detMat(matrix* m);
matrix* scalarMat(float s, matrix* mat1);
matrix* inverseMat(matrix* mat1, matrix* inv);
matrix* mpinverseMat(matrix* mat1, matrix* mpinv);
matrix* sinMat(matrix* mat1, matrix* sinM);
matrix* cosMat(matrix* mat1, matrix* cosM);
matrix* eleMultiplyMat(matrix* mat1, matrix* mat2, matrix* prod);
matrix* rightDivideMat(matrix* mat1, matrix* mat2, matrix* quot);
matrix* meanMat(matrix* mat1, matrix* mean);
matrix* findMat(matrix* mat1, float threshold, matrix* found);
matrix* conjoinMat(matrix* mat1, matrix* mat2, matrix* conjoined);
void freeMat(matrix* mat1);
//
///* helper */
float multProd (int len, float * row, float * col);
matrix* minorM(matrix* mat1, matrix* matco);
matrix* cofactor(matrix* mat1, matrix* matmin);
matrix* adjugate(matrix* mat1, matrix* admat);
//
/* print */
void printMat(matrix* mat1);


#endif /* matrix_h */
