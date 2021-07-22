//
//  matrix.c
//  csd
//
//  Created by Kevin Chen.
//  Copyright © 2020 Kevin Chen. All rights reserved.
//

#include "matrix.h"


/*
 
    CREATION
 
 */

/* allocate memory for matrix */
matrix makeMat(int m, int n)
{
    matrix matrix1;
    matrix1.row = m;
    matrix1.col = n;
    matrix1.data = malloc(m*n * sizeof(float));
    return matrix1;
}

/* assign values for matrix data */
matrix assignMat(int m, int n, float* data1)
{
    matrix mat1 = makeMat(m, n);
    mat1.data = data1;
    return mat1;
}

/*
 
    OPERATIONS
 
 */

/* returns the transpose of the matrix */
matrix transposeMat(matrix mat1)
{
    matrix trans = makeMat(mat1.col, mat1.row);
    for (int i = 0; i < mat1.col; i++)
    {
        for (int j = 0; j < mat1.row; j++)
        {
            trans.data[i*mat1.row+j] = mat1.data[j*mat1.col+i];
        }
    }
    return trans;
}

/* returns the product of the two matrices */
matrix multiplyMat(matrix mat1, matrix mat2)
{
    matrix temp1 = makeMat(mat1.row, mat2.col);
    for (int i = 0; i < temp1.row; i++)
    {
        // create array for the row
        float *row1 = malloc(mat1.col*sizeof(float));
        for (int i2 = 0; i2 < mat1.col; i2++) { row1[i2] = mat1.data[i*mat1.col+i2]; }
        
        for (int j = 0; j < temp1.col; j++)
        {
            // create array for the col
            float *col1 = malloc(mat2.row*sizeof(float));
            for (int j2 = 0; j2 < mat2.row; j2++) { col1[j2] = mat2.data[j2*mat2.col+j]; }
            
            
//            printf("\n------\n");
//            printf("row\n");
//            for (int i = 0; i < mat1.col; i++)
//            {
//                printf("%f ", row1[i]);
//            }
//            printf("\ncol\n");
//            for (int i = 0;i < mat2.row; i++)
//            {
//                printf("%f ", col1[i]);
//            }
//            printf("\n------\n");
            temp1.data[i*temp1.col+j] = multProd(mat1.col, row1, col1);
//            printf("i: %d\n", i);
//            printf("temp1.row: %d\n", temp1.row);
//            printf("j: %d\n", j);
//            printf("location: %d\n", i*temp1.col+j);
//            printf("multProd: %f\n", multProd(mat1.col, row1, col1));
//            printf("data: %f\n", temp1.data[i*temp1.row+j]);
//
//            printf("temp1\n");
//            printMat(temp1);
            
        }
    }
    return temp1;
}

/* recursive function for determinant of matrix (with Cramer's Rule) */
float detMat(matrix m)
{
    if (m.row != m.col)
    {
        // printf("ERROR: Can not compute determinant of non-square matrix");
        return -1;
    }
    else if (m.row == 1)
    {
        return m.data[0];
    }
    else if (m.row == 2)
    {
        return m.data[0]*m.data[3] - m.data[1]*m.data[2];
    }
    else
    {
        float det = 0;
        for (int i = 0; i < m.row; i++)
        {
            float * argdet = malloc((m.row-1)*(m.col-1)*sizeof(float));
            int pos = 0;
            for (int j = 0; j < m.row*m.col; j++)
            {
                if (j%m.row != i && m.row <= j)
                {
                    argdet[pos] = m.data[j];
                    pos++;
                }
            }
            
            matrix mdet = assignMat(m.row-1, m.col-1, argdet);
            
            /* checkerboard +/- */
            if (i%2 != 0) { det += -1*m.data[i]*detMat(mdet); }
            else { det += m.data[i]*detMat(mdet); }
        }
        
        return det;
    }
}


/* multiples matrix by the given scalar */
matrix scalarMat(float s, matrix mat1)
{
    matrix matsca = makeMat(mat1.row, mat1.col);
    for (int i = 0; i < mat1.row*mat1.col; i++)
    {
        matsca.data[i] = mat1.data[i] * s;
    }
    return matsca;
}

/* returns inverse matrix */
matrix inverseMat(matrix mat1)
{
    printf("mat1.col: %d\n", mat1.col);
    printf("mat1.row: %d\n", mat1.row);
    if (mat1.col != mat1.row)
    {
        printf("ERROR: Can not compute inverse of non-square matrix");
    }
    return scalarMat(1/detMat(mat1), adjugate(mat1));
    
}

/* returns Moore-Penrose inverse matrix */
matrix mpinverseMat(matrix mat1)
{

    // (X^t * X)^-1 * X^t
    matrix m1 = multiplyMat(transposeMat(mat1), mat1);

    if (detMat(m1) == 0)
    {
        printf("Moore-Penrose Inverse of Matrix Does Not Exist\n");
        return m1;
    }

    matrix m2 = inverseMat(m1);

    matrix m3 = multiplyMat(m2, transposeMat(mat1));
    

    return m3;
}

/* takes sin of each element in matrix */
matrix sinMat(matrix mat1)
{
    matrix mat2 = makeMat(mat1.row, mat1.col);
    for (int i = 0; i < mat1.row; i++)
    {
        for (int j = 0; j < mat1.col; j++)
        {
            mat2.data[i*mat1.col+j] = sin(mat1.data[i*mat1.col+j]);
        }
    }
    return mat2;
}

/* takes cos of each element in matrix*/
matrix cosMat(matrix mat1)
{
    matrix mat2 = makeMat(mat1.row, mat1.col);
    for (int i = 0; i < mat1.row; i++)
    {
        for (int j = 0; j < mat1.col; j++)
        {
            mat2.data[i*mat1.col+j] = cos(mat1.data[i*mat1.col+j]);
        }
    }
    return mat2;
}

/* element-wise matrix multiplication */
matrix eleMultiplyMat(matrix mat1, matrix mat2)
{
    // check to see if matrices are same size
    if (mat1.row != mat2.row && mat1.col != mat2.col)
    {
        printf("INCOMPATIBLE MATRICE DIMENSIONS (element-wise)");
    }
    
    // allocate memory
    float *s = malloc(mat1.row*mat1.col*sizeof(float));
    
    for (int i = 0; i < mat1.row*mat1.col; i++)
    {
        s[i] = mat1.data[i] * mat2.data[i];
    }
    
    matrix emat = assignMat(mat1.row, mat1.col, s);
    
    return emat;
}

/* right array matrix division (A / B) */
matrix rightDivideMat(matrix mat1, matrix mat2)
{
    
    // check to see if matrices are same size
    if (mat1.row != mat2.row && mat1.col != mat2.col)
    {
        printf("INCOMPATIBLE MATRICE DIMENSIONS (right array)");
    }
    
    // allocate memory
    float *s = malloc(mat1.row*mat1.col*sizeof(float));
    
    for (int i = 0; i < mat1.row*mat1.col; i++)
    {
        s[i] = mat1.data[i] / mat2.data[i];
    }
    
    matrix rmat = assignMat(mat1.row, mat1.col, s);
    
    return rmat;
    
}

/* return mean of colums */
matrix meanMat(matrix mat1)
{
    
    float * dmat = malloc(mat1.col*sizeof(float));
    for (int i = 0; i < mat1.col; i++)
    {
        float sum = 0;
        for (int j = 0; j < mat1.row; j++)
        {
            sum += mat1.data[i+j*mat1.col];
        }
        dmat[i] = sum/mat1.row;
    }
    matrix mat2 = assignMat(1, mat1.col, dmat);
    return mat2;
}

/* find where mat2 > mat1 */
matrix findMat(matrix mat1, matrix mat2)
{
    float* dm = malloc(mat1.row*mat1.col*sizeof(float));
    int counter = 0;
    
    for (int i = 0; i < mat1.col; i++)
    {
        for (int j = 0; j < mat1.row; j++)
        {
            if (mat2.data[j*mat1.col+i] > mat1.data[j*mat1.col+i])
            {
                dm[counter] = i*mat1.col+j + 1;
                counter++;
            }
        }
    }
    
    matrix m = assignMat(counter, 1, dm);
    
    return m;
}

/* conjoin matrices (side by side) */
matrix conjoinMat(matrix mat1, matrix mat2)
{
    if (mat1.row != mat2.row)
    {
        printf("ERROR: The matrices do not have the same number of rows");
    }
    float* dmat3 = malloc(sizeof(float) * mat1.row * (mat1.col+mat2.col));
    
    int counter = 0;
    for (int i = 0; i < mat1.row; i++)
    {
        for (int j = 0; j < mat1.col; j++)
        {
            dmat3[counter] = mat1.data[i*mat1.col+j];
            counter++;
        }
        for (int k = 0; k < mat2.col; k++)
        {
            dmat3[counter] = mat2.data[i*mat2.col+k];
            counter++;
        }
    }
    matrix mat3 = assignMat(mat1.row, mat1.col+mat2.col, dmat3);
    return mat3;
}


/*
    
    HELPER
 
 */
                         
/* dot product */
float multProd (int len, float * row, float * col)
{
    float prod = 0.0;
    for (int i = 0; i < len; i++)
    {
        prod += row[i] * col[i];
    }
    return prod;
}

/* returns minor matrix */
matrix minorM(matrix mat1)
{
    matrix matco = makeMat(mat1.row, mat1.col);
    
    for (int i = 0; i < matco.row*matco.col; i++)
    {
        float * argdet = malloc((mat1.row-1)*(mat1.col-1)*sizeof(float));
        int pos = 0;
        for (int j = 0; j < matco.row*matco.col; j++)
        {
        if (j%mat1.row != i%mat1.row && j/mat1.row != i/mat1.row)
          {
              argdet[pos] = mat1.data[j];
              pos++;
          }
        }
        matrix mdet = assignMat(mat1.row-1, mat1.col-1, argdet);
        
        matco.data[i] = detMat(mdet);
    }

    return matco;
}

/* returns cofactor matrix */
matrix cofactor(matrix mat1)
{
    matrix matmin = minorM(mat1);
    int mult1 = -1;
    int mult2 = 1;
    for (int i = 0; i < mat1.row; i++)
    {
        for (int j = 0; j < mat1.col; j++)
        {
            if (i%2==0) {mult1=-1; mult2=1;}
            else {mult1=1; mult2=-1;}
            
            if (j%2!=0)
            {
                matmin.data[i*mat1.col+j] = matmin.data[i*mat1.col+j] * mult1;
            }
            else
            {
                matmin.data[i*mat1.col+j] = matmin.data[i*mat1.col+j] * mult2;
            }
        }
    }
    
    return matmin;
}

/* returns adjugate (adjoint matrix) */
matrix adjugate(matrix mat1) { return transposeMat(cofactor(mat1)); }


/*
 
    PRINT
 
 */

/* prints matrix (for debugging) */
void printMat(matrix mat1)
{
    for (int i = 0; i < mat1.row; i++)
    {
        for (int j = 0; j < mat1.col; j++)
        {
            printf("%f, ", mat1.data[i*mat1.col+j]);
        }
        printf("\n");
    }
}



