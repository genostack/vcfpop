/* SSE Instruction Set Functions */

#pragma once
#include "vcfpop.h"

TARGETSSE uint64 GetMinIdxSSE(double *A, int64 n, double &val);

TARGETSSE void GetMinMaxValSSE(double *A, int64 n, double &minv, double &maxv);

TARGETSSE double GetMaxValSSE(double *A, int64 n);

TARGETSSE double GetMinValSSE(double *A, int64 n);

TARGETSSE void SetValSSE(uint *a, ushort *b, int64 n);

TARGETSSE double LogProdSSE(double *A, int64 n);

TARGETSSE double LogProdSSE(double *A, int64 n, int64 sep);

TARGETSSE double LogProdDivSSE(double *A, double *B, int64 n, int64 sep);

TARGETSSE uint64 CountNonZeroSSE(byte *A, int64 n);

TARGETSSE double SumSSE(double *A, int64 n);

TARGETSSE double SumSSE(byte *A, int64 n);

TARGETSSE double SumSSE(double *A, int64 n, int64 sep);

TARGETSSE void SumSSE(double *A, double **B, int64 k, int64 n);

TARGETSSE double ProdSSE(double *A, int64 n);

TARGETSSE double ProdSSE(double *A, int64 n, int64 sep);

TARGETSSE double SumSquareSSE(double *A, int64 n);

TARGETSSE double SumSquareSSE(byte *A, int64 n);

TARGETSSE void SumSumSquareSSE(double *A, int64 n, double &sum, double &sumsq);

TARGETSSE double SumProdSSE(double *A, double *B, int64 sep, int64 n);

TARGETSSE double SumProdSSE(double *A, double *B, int64 n);

TARGETSSE void AddSSE(double *A, double *B, int64 n);

TARGETSSE void AddSSE(int *A, int *B, int64 n);

TARGETSSE void AddSSE(double *A, double B, int64 n);

TARGETSSE void MulSSE(double *C, double *A, double *B, int64 n);

TARGETSSE void MulSSE(double *C, double *A, double B, int64 n);

TARGETSSE void MulSSE(double *A, double B, int64 n);

TARGETSSE void AddProdSSE(double *C, double *A, double *B, int64 n);

TARGETSSE void AddProdSSE(double *C, double *A, double B, int64 n);

TARGETSSE void UnifySSE(double *A, int64 n);

TARGETSSE char *StrNextIdxSSE(char *A, char val, int64 rep, int64 n);

TARGETSSE uint64 CountCharSSE(char *A, char val, int64 n);