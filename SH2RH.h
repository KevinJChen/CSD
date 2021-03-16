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


matrix eval_SH(int l, float el[], size_t size, matrix az);
matrix eval_ALP(int l, float el[], size_t size);

/* helper functions */
float * cosEl(float el[], size_t size);
float factorial(int num);

#endif /* SH2RH_h */
