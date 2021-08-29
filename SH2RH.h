//
//  SH2RH.h
//  csd
//
//  Created by Kevin Chen.
//  Copyright © 2021 Kevin Chen. All rights reserved.
//

#ifndef SH2RH_h
#define SH2RH_h

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "legendre.h"

//float* SH2RH(float* SH, size_t size_sh);
//float* gen_delta(float el[], size_t size_e, matrix az, int lmax);
//matrix eval_SH(int l, float el[], size_t size_e, matrix az);
matrix* eval_ALP(int l, float *el, size_t size, matrix* ALP);
matrix* amp2SH(float* S, size_t size, matrix* SH);
//matrix dir3002SH(float* dir300, size_t size);
//matrix bvec2SH(float* bvec, size_t size);

/* helper functions */
float * cosEl(float *el, size_t size);
float factorial(int num);

#endif /* SH2RH_h */
