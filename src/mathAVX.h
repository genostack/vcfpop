/* AVX Instruction Set Functions */

#pragma once
#include "vcfpop.h"

TARGETAVX uint64 GetMinIdxAVX(double *A, int64 n, double &val);

TARGETAVX void GetMinMaxValAVX(double *A, int64 n, double &minv, double &maxv);

TARGETAVX double GetMaxValAVX(double *A, int64 n);

TARGETAVX double GetMinValAVX(double *A, int64 n);

TARGETAVX void SetValAVX(uint *a, ushort *b, int64 n);

TARGETAVX double LogProdAVX(double *A, int64 n);

TARGETAVX double LogProdAVX(double *A, int64 n, int64 sep);

TARGETAVX double LogProdDivAVX(double *A, double *B, int64 n, int64 sep);

TARGETAVX uint64 CountNonZeroAVX(byte *A, int64 n);

TARGETAVX double SumAVX(double *A, int64 n);

TARGETAVX double SumAVX(byte *A, int64 n);

TARGETAVX double SumAVX(double *A, int64 n, int64 sep);

TARGETAVX void SumAVX(double *A, double **B, int64 k, int64 n);

TARGETAVX double ProdAVX(double *A, int64 n);

TARGETAVX double ProdAVX(double *A, int64 n, int64 sep);

TARGETAVX double SumSquareAVX(double *A, int64 n);

TARGETAVX double SumSquareAVX(byte *A, int64 n);

TARGETAVX void SumSumSquareAVX(double *A, int64 n, double &sum, double &sumsq);

TARGETAVX double SumProdAVX(double *A, double *B, int64 sep, int64 n);

TARGETAVX double SumProdAVX(double *A, double *B, int64 n);

TARGETAVX void AddAVX(double *A, double *B, int64 n);

TARGETAVX void AddAVX(int *A, int *B, int64 n);

TARGETAVX void AddAVX(double *A, double B, int64 n);

TARGETAVX void MulAVX(double *C, double *A, double *B, int64 n);

TARGETAVX void MulAVX(double *C, double *A, double B, int64 n);

TARGETAVX void MulAVX(double *A, double B, int64 n);

TARGETAVX void AddProdAVX(double *C, double *A, double *B, int64 n);

TARGETAVX void AddProdAVX(double *C, double *A, double B, int64 n);

TARGETAVX void UnifyAVX(double *A, int64 n);

TARGETAVX char *StrNextIdxAVX(char *A, char val, int64 rep, int64 n);

TARGETAVX uint64 CountCharAVX(char *A, char val, int64 n);