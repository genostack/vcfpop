/* Method-of-moment estimation of relatedness coefficient for polyploids */

#pragma once
#include "vcfpop.h"

double mp(double base, int index);
void MOMRelatednessAssign(int p, int refmode, double *e, double *f, int *xx);
void MOMRelatednessAssign1(int refmode, double *e, double *f, int *xx);
void MOMRelatednessAssign2(int refmode, double *e, double *f, int *xx);
void MOMRelatednessAssign3(int refmode, double *e, double *f, int *xx);
void MOMRelatednessAssign4(int refmode, double *e, double *f, int *xx);
void MOMRelatednessAssign5(int refmode, double *e, double *f, int *xx);
void MOMRelatednessAssign6(int refmode, double *e, double *f, int *xx);
void MOMRelatednessAssign7(int refmode, double *e, double *f, int *xx);
