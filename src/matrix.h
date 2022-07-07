/* Matrix Functions */

#pragma once
#include "vcfpop.h"

/* Maximum diagonal element */
TARGET double MaxDiag(double *a, uint m, uint n);

/* Matrix multiplication */
TARGET int MatrixMul(double *l, uint lr, uint lc, double *r, uint rr, uint rc, double *res);

/* Matrix inverstion */
TARGET int MatrixInv(double *M, uint m);

/* SVD decomposition */
TARGET void S(double fg[2], double cs[2]);

/* SVD decomposition */
TARGET void D(double *a, double *b, uint m, uint n, uint k, double *c);

/* SVD decomposition */
TARGET void P(double *a, double *e, double *s, double *v, uint m, uint n);

/* SVD decomposition */
TARGET int MatrixSVD(double *a, uint m, uint n, double *u, double *v, double eps = EPSILON_SVD);

/* L2 norm */
TARGET double MatrixNorm(double *a, uint m, uint n);

/* Condition number */
TARGET double MatrixCond(double *a, uint m, uint n);

/* Solve Ax = B */
TARGET bool SolveEquation(double *A, double *B, double *x, uint n);