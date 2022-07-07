/* Maximum-likelihood estimation of relatedness coefficient for polyploids */

#pragma once
#include "vcfpop.h"

void MLRelatednessAssign(int sumploidy, double *f, int *a, double *c, int IBS);