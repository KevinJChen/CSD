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

void csdeconv(float* R_RH, matrix SH, matrix S, float lambda, float tau);
int nSH_for_lmax(int lmax);


#endif /* CSD_h */
