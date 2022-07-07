/* AVX512 Instruction Set Functions */

#pragma once
#include "vcfpop.h"

TARGET512 uint64 GetMinIdx512(double *A, int64 n, double &val);

TARGET512 void GetMinMaxVal512(double *A, int64 n, double &minv, double &maxv);

TARGET512 double GetMaxVal512(double *A, int64 n);

TARGET512 double GetMinVal512(double *A, int64 n);

TARGET512 void SetVal512(uint *a, ushort *b, int64 n);

TARGET512 double LogProd512(double *A, int64 n);

TARGET512 double LogProd512(double *A, int64 n, int64 sep);

TARGET512 double LogProdDiv512(double *A, double *B, int64 n, int64 sep);

TARGET512 uint64 CountNonZero512(byte *A, int64 n);

TARGET512 double Sum512(double *A, int64 n);

TARGET512 double Sum512(byte *A, int64 n);

TARGET512 double Sum512(double *A, int64 n, int64 sep);

TARGET512 void Sum512(double *A, double **B, int64 k, int64 n);

TARGET512 double Prod512(double *A, int64 n);

TARGET512 double Prod512(double *A, int64 n, int64 sep);

TARGET512 double SumSquare512(double *A, int64 n);

TARGET512 double SumSquare512(byte *A, int64 n);

TARGET512 void SumSumSquare512(double *A, int64 n, double &sum, double &sumsq);

TARGET512 double SumProd512(double *A, double *B, int64 sep, int64 n);

TARGET512 double SumProd512(double *A, double *B, int64 n);

TARGET512 void Add512(double *A, double *B, int64 n);

TARGET512 void Add512(int *A, int *B, int64 n);

TARGET512 void Add512(double *A, double B, int64 n);

TARGET512 void Mul512(double *C, double *A, double *B, int64 n);

TARGET512 void Mul512(double *C, double *A, double B, int64 n);

TARGET512 void Mul512(double *A, double B, int64 n);

TARGET512 void AddProd512(double *C, double *A, double *B, int64 n);

TARGET512 void AddProd512(double *C, double *A, double B, int64 n);

TARGET512 void Unify512(double *A, int64 n);

TARGET512 char *StrNextIdx512(char *A, char val, int64 rep, int64 n);

TARGET512 uint64 CountChar512(char *A, char val, int64 n);