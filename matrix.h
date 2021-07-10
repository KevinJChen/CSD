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
matrix makeMat(int m, int n);
matrix assignMat(int m, int n, float* data1);


/* operations */
matrix transposeMat(matrix mat1);
matrix multiplyMat(matrix mat1, matrix mat2);
float detMat(matrix m);
matrix scalarMat(float s, matrix mat1);
matrix inverseMat(matrix mat1);
matrix mpinverseMat(matrix mat1);
matrix sinMat(matrix mat1);
matrix cosMat(matrix mat1);
matrix eleMultiplyMat(matrix mat1, matrix mat2);
matrix rightDivideMat(matrix mat1, matrix mat2);
matrix meanMat(matrix mat1);
matrix findMat(matrix mat1, matrix mat2);

/* helper */
float multProd (int len, float * row, float * col);
matrix minorM(matrix mat1);
matrix cofactor(matrix mat1);
matrix adjugate(matrix mat1);

/* print */
void printMat(matrix mat1);


#endif /* matrix_h */
