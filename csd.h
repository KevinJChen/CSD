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

void csdeconv(float* R_RH, float lambda, float tau);


#endif /* CSD_h */
