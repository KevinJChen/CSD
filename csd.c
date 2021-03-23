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

#include "csd.h"

int main()
{

    float* R_RH = malloc(5*sizeof(float));
    for (int i = 0; i < 5; i++)
    {
        R_RH[i] = i+1;
    }
    printf("Success!\n");
    return 0;
}

void csdeconv(float* R_RH, float lambda, float tau)
{
    
}


