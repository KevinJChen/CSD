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


/* assign values for matrix data */
matrix* assignMat(int m, int n, float* data1, matrix* mat1)
{
    //matrix* mat1 = malloc(sizeof(matrix));
    mat1->row = m;
    mat1->col = n;
    mat1->data = malloc(sizeof(float)*m*n);
    for (int i = 0; i < m*n; i++)
    {
        mat1->data[i] = data1[i];
    }
    return mat1;
}

/*
 
    OPERATIONS
 
 */

/* returns the transpose of the matrix */
matrix* transposeMat(matrix* mat1, matrix* trans)
{
    trans->row = mat1->col;
    trans->col = mat1->row;
    trans->data = malloc(sizeof(float)*mat1->row*mat1->col);
    for (int i = 0; i < mat1->col; i++)
    {
        for (int j = 0; j < mat1->row; j++)
        {
            trans->data[i*mat1->row+j] = mat1->data[j*mat1->col+i];
        }
    }
    return trans;
}

/* returns the product of the two matrices */
matrix* multiplyMat(matrix* mat1, matrix* mat2, matrix* prod)
{
    prod->row = mat1->row;
    prod->col = mat2->col;
    prod->data = malloc(sizeof(float)*mat1->row*mat2->col);
    // matrix temp1 = makeMat(mat1.row, mat2.col);
    for (int i = 0; i < prod->row; i++)
    {
        // create array for the row
        float *row1 = malloc(mat1->col*sizeof(float));
        for (int i2 = 0; i2 < mat1->col; i2++) { row1[i2] = mat1->data[i*mat1->col+i2]; }

        for (int j = 0; j < prod->col; j++)
        {
            // create array for the col
            float *col1 = malloc(mat2->row*sizeof(float));
            for (int j2 = 0; j2 < mat2->row; j2++) { col1[j2] = mat2->data[j2*mat2->col+j]; }
            prod->data[i*prod->col+j] = multProd(mat1->col, row1, col1);
            free(col1);
        }
        free(row1);
    }
    return prod;
}

/* recursive function for determinant of matrix (with Cramer's Rule) */
float detMat(matrix* m)
{
    if (m->row != m->col)
    {
        // printf("ERROR: Can not compute determinant of non-square matrix");
        return -1;
    }
    else if (m->row == 1)
    {
        return m->data[0];
    }
    else if (m->row == 2)
    {
        return m->data[0]*m->data[3] - m->data[1]*m->data[2];
    }
    else
    {
        float det = 0;
        for (int i = 0; i < m->row; i++)
        {
            float * argdet = malloc((m->row-1)*(m->col-1)*sizeof(float));
            int pos = 0;
            for (int j = 0; j < m->row*m->col; j++)
            {
                if (j%m->row != i && m->row <= j)
                {
                    argdet[pos] = m->data[j];
                    pos++;
                }
            }

            matrix* mdet = malloc(sizeof(matrix));
            assignMat(m->row-1, m->col-1, argdet, mdet);

            /* checkerboard +/- */
            if (i%2 != 0) { det += -1*m->data[i]*detMat(mdet); }
            else { det += m->data[i]*detMat(mdet); }
            free(argdet);
            free(mdet);
            
        }

        return det;
    }
}


/* multiples matrix by the given scalar */
matrix* scalarMat(float s, matrix* mat1)
{
    for (int i = 0; i < mat1->row*mat1->col; i++)
    {
        mat1->data[i] = mat1->data[i] * s;
    }
    return mat1;
}

/* returns inverse matrix */
matrix* inverseMat(matrix* mat1, matrix* inv)
{
    //printf("mat1.col: %d\n", mat1.col);
    //printf("mat1.row: %d\n", mat1.row);
    if (mat1->col != mat1->row)
    {
        printf("ERROR: Can not compute inverse of non-square matrix");
    }
    
    return scalarMat(1/detMat(mat1), adjugate(mat1, inv));

}

/* returns Moore-Penrose inverse matrix */
matrix* mpinverseMat(matrix* mat1, matrix* mpinv)
{
    
    // (X^t * X)^-1 * X^t
    matrix* trans = malloc(sizeof(matrix));
    transposeMat(mat1, trans);
    matrix* prod = malloc(sizeof(matrix));
    multiplyMat(trans, mat1, prod);
    if (detMat(prod) == 0)
    {
        printf("Moore-Penrose Inverse of Matrix Does Not Exist\n");
        return mpinv;
    }
    matrix* inv = malloc(sizeof(matrix));
    inverseMat(prod, inv);
    //multiplyMat(inv, trans, mpinv);
    
    
    //matrix m2 = inverseMat(m1);
    // matrix m3 = multiplyMat(m2, transposeMat(mat1));
    freeMat(trans);
    freeMat(prod);
    freeMat(inv);
    return mpinv;
}

/* takes sin of each element in matrix */
matrix* sinMat(matrix* mat1, matrix* sinM)
{
    //matrix mat2 = makeMat(mat1.row, mat1.col);
    sinM->row = mat1->row;
    sinM->col = mat1->col;
    sinM->data = malloc(sizeof(float) * mat1->row * mat1->col);
    for (int i = 0; i < mat1->row; i++)
    {
        for (int j = 0; j < mat1->col; j++)
        {
            sinM->data[i*mat1->col+j] = sin(mat1->data[i*mat1->col+j]);
        }
    }
    return sinM;
}

/* takes cos of each element in matrix*/
matrix* cosMat(matrix* mat1, matrix* cosM)
{
    // matrix mat2 = makeMat(mat1.row, mat1.col);
    cosM->row = mat1->row;
    cosM->col = mat1->col;
    cosM->data = malloc(sizeof(float) * mat1->row * mat1->col);
    for (int i = 0; i < mat1->row; i++)
    {
        for (int j = 0; j < mat1->col; j++)
        {
            cosM->data[i*mat1->col+j] = cos(mat1->data[i*mat1->col+j]);
        }
    }
    return cosM;
}

/* element-wise matrix multiplication */
matrix* eleMultiplyMat(matrix* mat1, matrix* mat2, matrix* prod)
{
    // check to see if matrices are same size
    if (mat1->row != mat2->row && mat1->col != mat2->col)
    {
        printf("INCOMPATIBLE MATRICE DIMENSIONS (element-wise)");
    }

    // allocate memory
    float *s = malloc(mat1->row*mat1->col*sizeof(float));

    for (int i = 0; i < mat1->row*mat1->col; i++)
    {
        s[i] = mat1->data[i] * mat2->data[i];
    }
    
    assignMat(mat1->row, mat1->col, s, prod);
    free(s);
    return prod;
}

/* right array matrix division (A / B) */
matrix* rightDivideMat(matrix* mat1, matrix* mat2, matrix* quot)
{

    // check to see if matrices are same size
    if (mat1->row != mat2->row && mat1->col != mat2->col)
    {
        printf("INCOMPATIBLE MATRICE DIMENSIONS (right array)");
    }

    // allocate memory
    float *s = malloc(mat1->row*mat1->col*sizeof(float));

    for (int i = 0; i < mat1->row*mat1->col; i++)
    {
        s[i] = mat1->data[i] / mat2->data[i];
    }

    assignMat(mat1->row, mat1->col, s, quot);
    free(s);

    return quot;

}

/* return mean of colums */
matrix* meanMat(matrix* mat1, matrix* mean)
{

    float * dmat = malloc(mat1->col*sizeof(float));
    for (int i = 0; i < mat1->col; i++)
    {
        float sum = 0;
        for (int j = 0; j < mat1->row; j++)
        {
            sum += mat1->data[i+j*mat1->col];
        }
        dmat[i] = sum/mat1->row;
    }
    assignMat(1, mat1->col, dmat, mean);
    free(dmat);
    return mean;
}

/* find where mat2.data < threshold */
matrix* findMat(matrix* mat1, float threshold, matrix* found)
{
    // mat 1 should always have one column
    int counter = 0;
    for (int i = 0 ; i < mat1->row; i++)
    {
        if (mat1->data[i] < threshold)
        {
            counter++;
        }
    }

    int dcounter = 0;
    float* dm = malloc(counter*sizeof(float));
    for (int i = 0; i < mat1->row; i++)
    {
        if (mat1->data[i] < threshold)
        {
            dm[dcounter] = i;
            dcounter++;
        }
    }
    assignMat(counter, 1, dm, found);
    free(dm);
    return found;
}

/* conjoin matrices (side by side) */
matrix* conjoinMat(matrix* mat1, matrix* mat2, matrix* conjoined)
{
    if (mat1->row != mat2->row)
    {
        printf("ERROR: The matrices do not have the same number of rows (conjoining matrices)");
    }
    float* dmat3 = malloc(sizeof(float) * mat1->row * (mat1->col+mat2->col));

    int counter = 0;
    for (int i = 0; i < mat1->row; i++)
    {
        for (int j = 0; j < mat1->col; j++)
        {
            dmat3[counter] = mat1->data[i*mat1->col+j];
            counter++;
        }
        for (int k = 0; k < mat2->col; k++)
        {
            dmat3[counter] = mat2->data[i*mat2->col+k];
            counter++;
        }
    }
    assignMat(mat1->row, mat1->col+mat2->col, dmat3, conjoined);
    free(dmat3);
    return conjoined;
}

void freeMat(matrix* mat1)
{
    free(mat1);
    free(mat1->data);
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
matrix* minorM(matrix* mat1, matrix* matco)
{
    //matrix matco = makeMat(mat1.row, mat1.col);
    matco->row = mat1->row;
    matco->col = mat1->col;
    for (int i = 0; i < matco->row*matco->col; i++)
    {
        float * argdet = malloc((mat1->row-1)*(mat1->col-1)*sizeof(float));
        int pos = 0;
        for (int j = 0; j < matco->row*matco->col; j++)
        {
        if (j%mat1->row != i%mat1->row && j/mat1->row != i/mat1->row)
          {
              argdet[pos] = mat1->data[j];
              pos++;
          }
        }
        matrix* mdet = malloc(sizeof(matrix));
        assignMat(mat1->row-1, mat1->col-1, argdet, mdet);

        matco->data[i] = detMat(mdet);
        free(argdet);
        freeMat(mdet);
    }
    return matco;
}

/* returns cofactor matrix */
matrix* cofactor(matrix* mat1, matrix* matmin)
{
    //matrix matmin = minorM(mat1);
    matmin = minorM(mat1, matmin);
    int mult1 = -1;
    int mult2 = 1;
    for (int i = 0; i < mat1->row; i++)
    {
        for (int j = 0; j < mat1->col; j++)
        {
            if (i%2==0) {mult1=-1; mult2=1;}
            else {mult1=1; mult2=-1;}

            if (j%2!=0)
            {
                matmin->data[i*mat1->col+j] = matmin->data[i*mat1->col+j] * mult1;
            }
            else
            {
                matmin->data[i*mat1->col+j] = matmin->data[i*mat1->col+j] * mult2;
            }
        }
    }
    return matmin;
}

/* returns adjugate (adjoint matrix) */
matrix* adjugate(matrix* mat1, matrix* admat)
{
    matrix* temp = malloc(sizeof(temp));
    cofactor(mat1,temp);
    transposeMat(temp, admat);
    freeMat(temp);
    return admat;
}


/*

    PRINT

 */

/* prints matrix (for debugging) */
void printMat(matrix* mat1)
{
    for (int i = 0; i < mat1->row; i++)
    {
        for (int j = 0; j < mat1->col; j++)
        {
            printf("%f, ", mat1->data[i*mat1->col+j]);
        }
        printf("\n");
    }
}



