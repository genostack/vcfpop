/* Math  Functions */

#include "vcfpop.h"

/* Is an erroneous real number */
TARGET bool IsError(double x)
{
	return isnan(x) || isinf(x);
}

/* Is a normal real number */
TARGET bool IsNormal(double x)
{
	return !IsError(x);
}

/* Find the index of the mimumum element */
TARGETMMX int64 GetMinIdx(double *A, int64 n, double &val)
{
	if (SIMD_TYPE >= 2 && n >= 8)
	{
		if (SIMD_TYPE >= 4) return GetMinIdx512(A, n, val);
		else if (SIMD_TYPE >= 3) return GetMinIdxAVX(A, n, val);
		else if (SIMD_TYPE >= 2) return GetMinIdxSSE(A, n, val);
	}
	int64 i = 0;
	val = 1e300;
	int64 idx = -1;
	for (; i < n; ++i, ++A)
	{
		if (*A > val) continue;
		val = *A;
		idx = i;
	}
	return idx;
}

/* Find maximum and minimum element of A */
TARGETMMX void GetMinMaxVal(double *A, int64 n, double &minv, double &maxv)
{
	if (SIMD_TYPE >= 2 && n >= 8)
	{
		if (SIMD_TYPE >= 4) return GetMinMaxVal512(A, n, minv, maxv);
		else if (SIMD_TYPE >= 3) return GetMinMaxValAVX(A, n, minv, maxv);
		else if (SIMD_TYPE >= 2) return GetMinMaxValSSE(A, n, minv, maxv);
	}
	int64 i = 0;
	minv = 1e300;
	maxv = -1e300;
	for (; i < n; ++i, ++A)
	{
		if (*A < minv) minv = *A;
		if (*A > maxv) maxv = *A;
	}
}

/* Find minimum element of A */
TARGETMMX double GetMinVal(double *A, int64 n)
{
	if (SIMD_TYPE >= 2 && n >= 8)
	{
		if (SIMD_TYPE >= 4) return GetMinVal512(A, n);
		else if (SIMD_TYPE >= 3) return GetMinValAVX(A, n);
		else if (SIMD_TYPE >= 2) return GetMinValSSE(A, n);
	}
	int64 i = 0;
	double val = 1e300;
	for (; i < n; ++i, ++A)
	{
		if (*A > val) continue;
		val = *A;
	}
	return val;
}

/* Find Maximum element of A */
TARGETMMX double GetMaxVal(double *A, int64 n)
{
	if (SIMD_TYPE >= 2 && n >= 8)
	{
		if (SIMD_TYPE >= 4) return GetMaxVal512(A, n);
		else if (SIMD_TYPE >= 3) return GetMaxValAVX(A, n);
		else if (SIMD_TYPE >= 2) return GetMaxValSSE(A, n);
	}
	int64 i = 0;
	double val = -1e300;
	for (; i < n; ++i, ++A)
	{
		if (*A < val) continue;
		val = *A;
	}
	return val;
}

/* A[i] = B[i] */
TARGETMMX void SetVal(uint *A, ushort *B, int64 n)
{
	if (SIMD_TYPE >= 4 && n >= 32) return SetVal512(A, B, n);
	if (SIMD_TYPE >= 3 && n >= 16) return SetValAVX(A, B, n);
	if (SIMD_TYPE >= 2 && n >= 8) return SetValSSE(A, B, n);
	int64 i = 0;
	for (; i < n; ++i)
		*A++ = *B++;
}

/* log(prod(A(1 to n))) */
TARGETMMX double LogProd(double *A, int64 n)
{
	int64 i = 0;
	double li = 0, li2 = 1;
	if (SIMD_TYPE >= 2 && n >= 8)
	{
		if (SIMD_TYPE >= 4) return LogProd512(A, n);
		else if (SIMD_TYPE >= 3) return LogProdAVX(A, n);
		else if (SIMD_TYPE >= 2) return LogProdSSE(A, n);
	}
	for (; i < n; ++i, ++A)
	{
		if (*A < DOUBLE_UNDERFLOW || *A > DOUBLE_OVERFLOW)
			li += log(*A);
		else
		{
			li2 *= *A;
			if (li2 < DOUBLE_UNDERFLOW || li2 > DOUBLE_OVERFLOW)
			{
				li += log(li2);
				li2 = 1.0;
			}
		}
	}
	li += log(li2);
	li2 = 1;
	return li;
}

/* log(prod(A(1, 1+sep, ..., n))) */
TARGETMMX double LogProd(double *A, int64 n, int64 sep)
{
	if (SIMD_TYPE >= 2 && n >= 8)
	{
		if (SIMD_TYPE >= 4) return LogProd512(A, n, sep);
		else if (SIMD_TYPE >= 3) return LogProdAVX(A, n, sep);
		else if (SIMD_TYPE >= 2) return LogProdSSE(A, n, sep);
	}
	int64 i = 0;
	double li = 0, li2 = 1;
	for (; i < n; ++i, A += sep)
	{
		if (*A < DOUBLE_UNDERFLOW || *A > DOUBLE_OVERFLOW)
			li += log(*A);
		else
		{
			li2 *= *A;
			if (li2 < DOUBLE_UNDERFLOW || li2 > DOUBLE_OVERFLOW)
			{
				li += log(li2);
				li2 = 1;
			}
		}
	}
	li += log(li2);
	li2 = 1;
	return li;
}

/* log(prod(A/B (1, 1+sep, ..., n) )) */
TARGETMMX double LogProdDiv(double *A, double *B, int64 n, int64 sep)
{
	if (SIMD_TYPE >= 2 && n >= 8)
	{
		if (SIMD_TYPE >= 4) return LogProdDiv512(A, B, n, sep);
		else if (SIMD_TYPE >= 3) return LogProdDivAVX(A, B, n, sep);
		else if (SIMD_TYPE >= 2) return LogProdDivSSE(A, B, n, sep);
	}
	int64 i = 0;
	double li = 0, li2 = 1;
	for (; i < n; ++i, A += sep, B += sep)
	{
		double val = *A / *B;
		if (val < DOUBLE_UNDERFLOW || val > DOUBLE_OVERFLOW)
			li += log(val);
		else
		{
			li2 *= val;
			if (li2 < DOUBLE_UNDERFLOW || li2 > DOUBLE_OVERFLOW)
			{
				li += log(li2);
				li2 = 1.0;
			}
		}
	}
	li += log(li2);
	li2 = 1;
	return li;
}

/* Count non-zero elements */
TARGETMMX int64 CountNonZero(byte *A, int64 n)
{
	// for ploidy
	if (SIMD_TYPE >= 4 && n >= 64) return CountNonZero512(A, n);
	else if (SIMD_TYPE >= 3 && n >= 32) return CountNonZeroAVX(A, n);
	else if (SIMD_TYPE >= 2 && n >= 16) return CountNonZeroSSE(A, n);
	int64 re = 0;
	int64 i = 0;
	for (; i < n; ++i, ++A)
		if (*A) re++;
	return re;
}

/* Sum of A */
TARGETMMX double Sum(double *A, int64 n)
{
	if (SIMD_TYPE >= 2 && n >= 8)
	{
		if (SIMD_TYPE >= 4) return Sum512(A, n);
		else if (SIMD_TYPE >= 3) return SumAVX(A, n);
		else if (SIMD_TYPE >= 2) return SumSSE(A, n);
	}
	int64 i = 0;
	double re = 0;
	for (; i < n; ++i)
		re += *A++;
	return re;
}

/* Sum of A */
TARGETMMX double Sum(byte *A, int64 n)
{
	// for ploidy level
	if (SIMD_TYPE >= 4 && n >= 64) return Sum512(A, n);
	if (SIMD_TYPE >= 3 && n >= 32) return SumAVX(A, n);
	if (SIMD_TYPE >= 2 && n >= 16) return SumSSE(A, n);
	uint64 re = 0;
	int64 i = 0;
	for (; i < n; ++i)
		re += *A++;
	return (double)re;
}

/* re += A[i += sep] */
TARGETMMX double Sum(double *A, int64 n, int64 sep)
{
	if (SIMD_TYPE >= 2 && n >= 8)
	{
		if (SIMD_TYPE >= 4) return Sum512(A, n, sep);
		else if (SIMD_TYPE >= 3) return SumAVX(A, n, sep);
		else if (SIMD_TYPE >= 2) return SumSSE(A, n, sep);
	}
	int64 i = 0;
	double re = 0;
	for (; i < n; ++i, A += sep)
		re += *A;
	return re;
}

/* A[i] = B[0][i] + ... + B[k][i] */
TARGETMMX void Sum(double *A, double **B, int64 k, int64 n)
{
	if (SIMD_TYPE >= 2 && n >= 8)
	{
		if (SIMD_TYPE >= 4) return Sum512(A, B, k, n);
		else if (SIMD_TYPE >= 3) return SumAVX(A, B, k, n);
		else if (SIMD_TYPE >= 2) return SumSSE(A, B, k, n);
	}
	int64 i = 0;
	for (; i < n; ++i, A++)
	{
		*A = B[0][i];
		for (int64 j = 1; j < k; ++j)
			*A += B[j][i];
	}
}

/* Product of A */
TARGETMMX double Prod(double *A, int64 n)
{
	if (SIMD_TYPE >= 2 && n >= 8)
	{
		if (SIMD_TYPE >= 4) return Prod512(A, n);
		else if (SIMD_TYPE >= 3) return ProdAVX(A, n);
		else if (SIMD_TYPE >= 2) return ProdSSE(A, n);
	}
	int64 i = 0;
	double re = 1;
	for (; i < n; ++i)
		re *= *A++;
	return re;
}

/* re *= A[i += sep] */
TARGETMMX double Prod(double *A, int64 n, int64 sep)
{
	if (SIMD_TYPE >= 2 && n >= 8)
	{
		if (SIMD_TYPE >= 4) return Prod512(A, n, sep);
		else if (SIMD_TYPE >= 3) return ProdAVX(A, n, sep);
		else if (SIMD_TYPE >= 2) return ProdSSE(A, n, sep);
	}
	int64 i = 0;
	double re = 1;
	for (; i < n; ++i, A += sep)
		re *= *A;
	return re;
}

/* Sum of squared A */
TARGETMMX double SumSquare(double *A, int64 n)
{
	if (SIMD_TYPE >= 2 && n >= 8)
	{
		if (SIMD_TYPE >= 4) return SumSquare512(A, n);
		else if (SIMD_TYPE >= 3) return SumSquareAVX(A, n);
		else if (SIMD_TYPE >= 2) return SumSquareSSE(A, n);
	}
	int64 i = 0;
	double re = 0;
	for (; i < n; ++i, ++A)
		re += *A * *A;
	return re;
}

/* Sum of squared A */
TARGETMMX double SumSquare(byte *A, int64 n)
{
	// for ploidy
	if (SIMD_TYPE >= 4 && n >= 64) return SumSquare512(A, n);
	if (SIMD_TYPE >= 3 && n >= 32) return SumSquareAVX(A, n);
	if (SIMD_TYPE >= 2 && n >= 16) return SumSquareSSE(A, n);
	int64 i = 0;
	uint64 re = 0;
	for (; i < n; ++i, ++A)
		re += *A * *A;
	return (double)re;
}

/* Sum of A and Sum of squared A */
TARGETMMX void SumSumSquare(double *A, int64 n, double &sum, double &sumsq)
{
	if (SIMD_TYPE >= 2 && n >= 8)
	{
		if (SIMD_TYPE >= 4) return SumSumSquare512(A, n, sum, sumsq);
		else if (SIMD_TYPE >= 3) return SumSumSquareAVX(A, n, sum, sumsq);
		else if (SIMD_TYPE >= 2) return SumSumSquareSSE(A, n, sum, sumsq);
	}
	int64 i = 0;
	sum = sumsq = 0;
	for (; i < n; ++i, ++A)
	{
		sum += *A;
		sumsq += *A * *A;
	}
}

/* re = Sum(A[i++] * B[j += sep]) */
TARGETMMX double SumProd(double *A, double *B, int64 sep, int64 n)
{

	if (SIMD_TYPE >= 2 && n >= 8)
	{
		if (SIMD_TYPE >= 4) return SumProd512(A, B, sep, n);
		else if (SIMD_TYPE >= 3) return SumProdAVX(A, B, sep, n);
		else if (SIMD_TYPE >= 2) return SumProdSSE(A, B, sep, n);
	}
	int64 i = 0;
	double re = 0;
	for (; i < n; ++i, A++, B += sep)
		re += *A * *B;
	return re;
}

/* re = Sum(A[i] * B[i]) */
TARGETMMX double SumProd(double *A, double *B, int64 n)
{
	if (SIMD_TYPE >= 2 && n >= 8)
	{
		if (SIMD_TYPE >= 4) return SumProd512(A, B, n);
		else if (SIMD_TYPE >= 3) return SumProdAVX(A, B, n);
		else if (SIMD_TYPE >= 2) return SumProdSSE(A, B, n);
	}
	int64 i = 0;
	double re = 0;
	for (; i < n; ++i, ++A, ++B)
		re += *A * *B;
	return re;
}

/* Add B into A, A[i] += B[i] */
TARGETMMX void Add(double *A, double *B, int64 n)
{
	if (SIMD_TYPE >= 2 && n >= 8)
	{
		if (SIMD_TYPE >= 4) return Add512(A, B, n);
		else if (SIMD_TYPE >= 3) return AddAVX(A, B, n);
		else if (SIMD_TYPE >= 2) return AddSSE(A, B, n);
	}
	int64 i = 0;
	for (; i < n; ++i, A++, B++)
		*A += *B;
}

/* Add B into A, A[i] += B[i] */
TARGETMMX void Add(int *A, int *B, int64 n)
{
	if (SIMD_TYPE >= 2 && n >= 16)
	{
		if (SIMD_TYPE >= 4) return Add512(A, B, n);
		else if (SIMD_TYPE >= 3) return AddAVX(A, B, n);
		else if (SIMD_TYPE >= 2) return AddSSE(A, B, n);
	}
	int64 i = 0;
	for (; i < n; ++i, A++, B++)
		*A += *B;
}

/* Add B into A, A[i] += B */
TARGETMMX void Add(double *A, double B, int64 n)
{
	if (SIMD_TYPE >= 2 && n >= 8)
	{
		if (SIMD_TYPE >= 4) return Add512(A, B, n);
		else if (SIMD_TYPE >= 3) return AddAVX(A, B, n);
		else if (SIMD_TYPE >= 2) return AddSSE(A, B, n);
	}
	for (int64 i = 0; i < n; ++i, A++)
		*A += B;
}

/* C[i] = A[i] * B[i] */
TARGETMMX void Mul(double *C, double *A, double *B, int64 n)
{
	if (SIMD_TYPE >= 2 && n >= 8)
	{
		if (SIMD_TYPE >= 4) return Mul512(C, A, B, n);
		else if (SIMD_TYPE >= 3) return MulAVX(C, A, B, n);
		else if (SIMD_TYPE >= 2) return MulSSE(C, A, B, n);
	}
	int64 i = 0;
	for (; i < n; ++i)
		*C++ = *A++ * *B++;
}

/* C[i] = A[i] * B */
TARGETMMX void Mul(double *C, double *A, double B, int64 n)
{
	if (SIMD_TYPE >= 2 && n >= 8)
	{
		if (SIMD_TYPE >= 4) return Mul512(C, A, B, n);
		else if (SIMD_TYPE >= 3) return MulAVX(C, A, B, n);
		else if (SIMD_TYPE >= 2) return MulSSE(C, A, B, n);
	}
	int64 i = 0;
	for (; i < n; ++i)
		*C++ = *A++ * B;
}

/* A[i] *= B */
TARGETMMX void Mul(double *A, double B, int64 n)
{
	if (SIMD_TYPE >= 2 && n >= 8)
	{
		if (SIMD_TYPE >= 4) return Mul512(A, B, n);
		else if (SIMD_TYPE >= 3) return MulAVX(A, B, n);
		else if (SIMD_TYPE >= 2) return MulSSE(A, B, n);
	}
	int64 i = 0;
	for (; i < n; ++i)
		*A++ *= B;
}

/* C[i] += A[i] * B[i] */
TARGETMMX void AddProd(double *C, double *A, double *B, int64 n)
{
	if (SIMD_TYPE >= 2 && n >= 8)
	{
		if (SIMD_TYPE >= 4) return AddProd512(C, A, B, n);
		else if (SIMD_TYPE >= 3) return AddProdAVX(C, A, B, n);
		else if (SIMD_TYPE >= 2) return AddProdSSE(C, A, B, n);
	}
	int64 i = 0;
	for (; i < n; ++i)
		*C++ += *A++ * *B++;
}

/* C[i] += A[i] * B */
TARGETMMX void AddProd(double *C, double *A, double B, int64 n)
{
	if (SIMD_TYPE >= 2 && n >= 8)
	{
		if (SIMD_TYPE >= 4) return AddProd512(C, A, B, n);
		else if (SIMD_TYPE >= 3) return AddProdAVX(C, A, B, n);
		else if (SIMD_TYPE >= 2) return AddProdSSE(C, A, B, n);
	}
	int64 i = 0;
	for (; i < n; ++i)
		*C++ += *A++ * B;
}

/* Set the sum of A to one */
TARGETMMX void Unify(double *A, int64 n)
{
	//A[i] = A[i] * invs + MIN_FREQ * invs
	if (SIMD_TYPE >= 2 && n >= 8)
	{
		if (SIMD_TYPE >= 4) return Unify512(A, n);
		else if (SIMD_TYPE >= 3) return UnifyAVX(A, n);
		else if (SIMD_TYPE >= 2) return UnifySSE(A, n);
	}
	int64 i = 0;
	double invsum = 1.0 / (Sum(A, n) + n * MIN_FREQ);
	for (; i < n; ++i, ++A)
		*A = (*A + MIN_FREQ) * invsum;
}

/* Find next position of val in string A*/
TARGETMMX char *StrNextIdx(char *A, char val, int64 rep, int64 n)
{
	if (!n) return NULL;
	if (SIMD_TYPE >= 4 && n >= 64) return StrNextIdx512(A, val, rep, n);
	else if (SIMD_TYPE >= 3 && n >= 32) return StrNextIdxAVX(A, val, rep, n);
	else if (SIMD_TYPE >= 2 && n >= 16) return StrNextIdxSSE(A, val, rep, n);
	A++; n--;
	for (int64 i = 0; i < n; ++i, A++)
		if (*A == val && !--rep)
			return A;
	return NULL;
}

/* Count val in string A */
TARGETMMX int64 CountChar(char *A, char val, int64 n)
{
	if (SIMD_TYPE >= 4 && n >= 64) return CountChar512(A, val, n);
	else if (SIMD_TYPE >= 3 && n >= 32) return CountCharAVX(A, val, n);
	else if (SIMD_TYPE >= 2 && n >= 16) return CountCharSSE(A, val, n);
	int64 re = 0;
	int64 i = 0;
	for (; i < n; ++i, A++)
		if (*A == val) re++;
	return re;
}

/* A[i] = B */
TARGET void SetVal(double *A, double B, int64 n, int64 sep)
{
	for (int64 i = 0; i < n; ++i, A += sep)
		*A = B;
}

/* re += A[i] * B[j] * alen[i * k + j] */
TARGET double SumProdSMM(ushort *alen, double *A, double *B, int64 k)
{
	alen += k;
	double re = 0;
	for (int64 i = 0; i < k; ++i)
		for (int64 j = 0; j < k; ++j)
			re += A[i] * B[j] * alen[i * k + j];
	return re;
}

/* re += A[i] * alen[i * k] */
TARGET double SumProdSMM(ushort *alen, double *A, ushort j, int64 k)
{
	alen += k;
	double re = 0;
	alen += j;
	for (int64 i = 0; i < k; ++i, ++A, alen += k)
		re += *A * *alen;
	return re;
}

/* Natural logarithm with bounds */
TARGET double MyLog(double val)
{
	if (val < 1e-300)
		return -6.90775527898214000E+02;
	else if (val > 1e300)
		return +6.90775527898214000E+02;
	return log(val);
}

/* Charge a value to fast sum log */
TARGET void OpenLog(double &v1, double &v2)
{
	v1 = 0;
	v2 = 1;
}

/* Charge a value to fast sum log */
TARGET void OpenLog(double *v1, double *v2, int64 n)
{
	SetZero(v1, n);
	SetVal(v2, 1.0, n);
}

/* Charge a value to fast sum log */
TARGET void ChargeLog(double &v1, double &v2, double val)
{
	if (val < DOUBLE_UNDERFLOW || val > DOUBLE_OVERFLOW)
		v1 += log(val);
	else
	{
		v2 *= val;
		if (v2 < DOUBLE_UNDERFLOW || v2 > DOUBLE_OVERFLOW)
		{
			v1 += log(v2);
			v2 = 1;
		}
	}
}

/* Charge a value to fast sum log */
TARGET void CloseLog(double &v1, double &v2)
{
	v1 += log(v2);
	v2 = 1;
}

/* Charge a value to fast sum log */
TARGET void ChargeLog(double *v1, double *v2, double *val, int64 n, int64 sep)
{
	for (int64 i = 0; i < n; ++i, val += sep)
		ChargeLog(v1[i], v2[i], *val);
}

/* Finalize fast sum log */
TARGET void CloseLog(double *v1, double *v2, int64 n)
{
	for (int64 i = 0; i < n; ++i)
	{
		v1[i] += log(v2[i]);
		v2[i] = 1;
	}
}

/* Count number of non-zero elements */
TARGET int64 CountNonZero(double *A, int64 n)
{
	//allele freq
	int64 count = 0;
	for (int64 i = 0; i < n; ++i)
		if (A[i]) count++;
	return count;
}

/* Set the sum of A to one */
TARGET void Unify(double *A, int64 m, int64 n)
{
	for (int64 i = 0; i < m; ++i)
	{
		Unify(A, n);
		A += n;
	}
}

/* Calculate SSWP in AMOVA */
TARGET double SSP(double *p, int64 k, int64 nhap, bool isiam, ushort *alen2)
{
	//freq array, without missing alleles
	if (!k || !nhap) return 0;
	if (isiam)
	{
		double s1 = (double)nhap, s2 = SumSquare(p, k) * nhap * nhap;
		if (!s1) return 0;
		return (s1 * s1 - s2) / (double)(2 * s1);
	}
	else
	{
		alen2 += k;
		double re = 0;
		for (int64 i = 0; i < k; ++i)
			for (int64 j = 0; j < i; ++j)
				re += p[i] * p[j] * alen2[i * k + j];
		return re * nhap;
	}
}

/* Calculate SSTOT in AMOVA */
TARGET double SSC(double *a, int64 k, bool isiam, ushort *alen2)
{
	//allele count array, without missing alleles
	if (!k) return 0;
	double s1 = 0, s2 = 0;
	if (isiam)
	{
		SumSumSquare(a, k, s1, s2);
		if (!s1) return 0;
		return (s1 * s1 - s2) / (double)(2 * s1);
	}
	else
	{
		alen2 += k;
		double re = 0, nt = 0;
		for (int64 i = 0; i < k; ++i)
		{
			nt += a[i];
			for (int64 j = 0; j < i; ++j)
				re += a[i] * a[j] * alen2[i * k + j];
		}
		return nt > 0 ? re / nt : 0;
	}
}

