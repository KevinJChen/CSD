//
//  csd.h
//  csd
//
//  Created by Kevin Chen.
//  Copyright © 2021 Kevin Chen. All rights reserved.
//

#ifndef CSD_h
#define CSD_h

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "legendre.h"
#include "SH2RH.h"
#include "qr_solve.h"

matrix* csdeconv(float* R_RH, matrix* DW_SH, matrix* HR_SH, matrix* S, float lambda, float tau, matrix* F_SH);
matrix* SH2amp(matrix* csd, float* V, matrix* amp);
int nSH_for_lmax(int lmax);

#endif /* CSD_h */
