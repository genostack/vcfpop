/* AVX512 Instruction Set Functions */

#ifndef __SSE3__
#define __SSE3__
#define __SSE4_1__
#define __SSE4_2__
#define __POPCNT__
#define __LZCNT__
#endif
#ifndef __AVX__
#define __AVX__
#define __AVX2__
#define __FMA__
#define __AVX512F__
#define __AVX512BW__
#endif
#include "vcfpop.h"

TARGET512 uint64 GetMinIdx512(double *A, int64 n, double &val)
{
	int64 i = 0;
	val = 1e300;
	uint64 idx = (uint64)-1;
	__m512d mval = _mm512_set1_pd(val), a;
	__m512i midx = _mm512_set1_epi8((char)0xFF);
	__m512i nidx = _mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0);
	__m512i msep = _mm512_set1_epi64(8);
	for (int64 l1 = n - 8; i <= l1; i += 8, A += 8)
	{
		a = _mm512_loadu_pd(A);
		midx = _mm512_mask_and_epi64(midx, _mm512_cmp_pd_mask(mval, a, 14), nidx, nidx);
		mval = _mm512_min_pd(mval, a);
		nidx = _mm512_add_epi64(nidx, msep);
	}
	for (int64 j = 0; j < 8; ++j)
	{
		if (m512d_f64(mval, j) > val) continue;
		val = m512d_f64(mval, j);
		idx = m512i_u64(midx, j);
	}
	for (; i < n; ++i, ++A)
	{
		if (*A > val) continue;
		val = *A;
		idx = i;
	}
	return idx;
}

TARGET512 void GetMinMaxVal512(double *A, int64 n, double &minv, double &maxv)
{
	int64 i = 0;
	minv = 1e300;
	maxv = -1e300;
	__m512d mval1 = _mm512_set1_pd(minv), mval2 = _mm512_set1_pd(maxv), a;
	for (int64 l1 = n - 8; i <= l1; i += 8, A += 8)
	{
		a = _mm512_loadu_pd(A);
		mval1 = _mm512_min_pd(mval1, a);
		mval2 = _mm512_max_pd(mval2, a);
	}
	for (int64 j = 0; j < 8; ++j)
	{
		if (m512d_f64(mval1, j) < minv) minv = m512d_f64(mval1, j);
		if (m512d_f64(mval2, j) > maxv) maxv = m512d_f64(mval2, j);
	}
	for (; i < n; ++i, ++A)
	{
		if (*A < minv) minv = *A;
		if (*A > maxv) maxv = *A;
	}
}

TARGET512 double GetMaxVal512(double *A, int64 n)
{
	int64 i = 0;
	double val = -1e300;
	__m512d mval = _mm512_set1_pd(val);
	for (int64 l1 = n - 8; i <= l1; i += 8, A += 8)
		mval = _mm512_max_pd(mval, _mm512_loadu_pd(A));
	for (int64 j = 0; j < 8; ++j)
	{
		if (m512d_f64(mval, j) < val) continue;
		val = m512d_f64(mval, j);
	}
	for (; i < n; ++i, ++A)
	{
		if (*A < val) continue;
		val = *A;
	}
	return val;
}

TARGET512 double GetMinVal512(double *A, int64 n)
{
	int64 i = 0;
	double val = 1e300;
	__m512d mval = _mm512_set1_pd(val);
	for (int64 l1 = n - 8; i <= l1; i += 8, A += 8)
		mval = _mm512_min_pd(mval, _mm512_loadu_pd(A));
	for (int64 j = 0; j < 8; ++j)
	{
		if (m512d_f64(mval, j) > val) continue;
		val = m512d_f64(mval, j);
	}
	for (; i < n; ++i, ++A)
	{
		if (*A > val) continue;
		val = *A;
	}
	return val;
}

TARGET512 void SetVal512(uint *A, ushort *B, int64 n)
{
	int64 i = 0;
	__m512i t0;
	for (int64 l1 = n - 32; i <= l1; i += 32, A += 32, B += 32)
	{
		_mm512_storeu_si512((__m512i*)A, _mm512_cvtepu16_epi32(_mm512_castsi512_si256(t0 = _mm512_loadu_si512((__m512i*)(B)))));
		_mm512_storeu_si512((__m512i*)(A + 16), _mm512_cvtepu16_epi32(_mm512_extracti64x4_epi64(t0, 1)));
	}
	for (; i < n; ++i)
		*A++ = *B++;
}

TARGET512 double LogProd512(double *A, int64 n)
{
	int64 i = 0;
	double li = 0, li2 = 1;
	__m512d pd = _mm512_set1_pd(1.0), duner = _mm512_set1_pd(DOUBLE_UNDERFLOW), dover = _mm512_set1_pd(DOUBLE_OVERFLOW);
	__mmask8 flag;
	for (int64 l1 = n - 8; i <= l1; i += 8, A += 8)
	{
		pd = _mm512_mul_pd(pd, _mm512_loadu_pd(A));
		flag = _mm512_cmp_pd_mask(pd, duner, 1) | _mm512_cmp_pd_mask(dover, pd, 1);
		if (flag)
		{
			if (flag & 1) { li += log(m512d_f64(pd, 0)); m512d_f64(pd, 0) = 1.0; }
			if (flag & 2) { li += log(m512d_f64(pd, 1)); m512d_f64(pd, 1) = 1.0; }
			if (flag & 4) { li += log(m512d_f64(pd, 2)); m512d_f64(pd, 2) = 1.0; }
			if (flag & 8) { li += log(m512d_f64(pd, 3)); m512d_f64(pd, 3) = 1.0; }
			if (flag & 16) { li += log(m512d_f64(pd, 4)); m512d_f64(pd, 4) = 1.0; }
			if (flag & 32) { li += log(m512d_f64(pd, 5)); m512d_f64(pd, 5) = 1.0; }
			if (flag & 64) { li += log(m512d_f64(pd, 6)); m512d_f64(pd, 6) = 1.0; }
			if (flag & 128) { li += log(m512d_f64(pd, 7)); m512d_f64(pd, 7) = 1.0; }
		}
	}
	for (int64 j = 0; j < 8; ++j)
	{
		li2 *= m512d_f64(pd, j);
		if (li2 < DOUBLE_UNDERFLOW || li2 > DOUBLE_OVERFLOW)
		{
			li += log(li2);
			li2 = 1.0;
		}
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

TARGET512 double LogProd512(double *A, int64 n, int64 sep)
{
	int64 i = 0;
	double li = 0, li2 = 1;
	__m512d pd = _mm512_set1_pd(1.0), duner = _mm512_set1_pd(DOUBLE_UNDERFLOW), dover = _mm512_set1_pd(DOUBLE_OVERFLOW);
	__m512i vindex = _mm512_set_epi64(sep * 7, sep * 6, sep * 5, sep * 4, sep * 3, sep * 2, sep, 0);
	__mmask8 flag;
	for (int64 l1 = n - 8, a8 = sep << 3; i <= l1; i += 8, A += a8)
	{
		pd = _mm512_mul_pd(pd, _mm512_i64gather_pd(vindex, A, 8));
		flag = _mm512_cmp_pd_mask(pd, duner, 1) | _mm512_cmp_pd_mask(dover, pd, 1);
		if (flag)
		{
			if (flag & 1) { li += log(m512d_f64(pd, 0)); m512d_f64(pd, 0) = 1.0; }
			if (flag & 2) { li += log(m512d_f64(pd, 1)); m512d_f64(pd, 1) = 1.0; }
			if (flag & 4) { li += log(m512d_f64(pd, 2)); m512d_f64(pd, 2) = 1.0; }
			if (flag & 8) { li += log(m512d_f64(pd, 3)); m512d_f64(pd, 3) = 1.0; }
			if (flag & 16) { li += log(m512d_f64(pd, 4)); m512d_f64(pd, 4) = 1.0; }
			if (flag & 32) { li += log(m512d_f64(pd, 5)); m512d_f64(pd, 5) = 1.0; }
			if (flag & 64) { li += log(m512d_f64(pd, 6)); m512d_f64(pd, 6) = 1.0; }
			if (flag & 128) { li += log(m512d_f64(pd, 7)); m512d_f64(pd, 7) = 1.0; }
		}
	}
	for (int64 j = 0; j < 8; ++j)
	{
		li2 *= m512d_f64(pd, j);
		if (li2 < DOUBLE_UNDERFLOW || li2 > DOUBLE_OVERFLOW)
		{
			li += log(li2);
			li2 = 1.0;
		}
	}
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

TARGET512 double LogProdDiv512(double *A, double *B, int64 n, int64 sep)
{
	int64 i = 0;
	double li = 0, li2 = 1;
	__m512d pd = _mm512_set1_pd(1.0), duner = _mm512_set1_pd(DOUBLE_UNDERFLOW), dover = _mm512_set1_pd(DOUBLE_OVERFLOW);
	__m512i vindex = _mm512_set_epi64(sep * 7, sep * 6, sep * 5, sep * 4, sep * 3, sep * 2, sep, 0);
	__mmask8 flag;
	for (int64 l1 = n - 8, s8 = sep << 3; i <= l1; i += 8, A += s8, B += s8)
	{
		pd = _mm512_mul_pd(pd, _mm512_div_pd(_mm512_i64gather_pd(vindex, A, 8), _mm512_i64gather_pd(vindex, B, 8)));
		flag = _mm512_cmp_pd_mask(pd, duner, 1) | _mm512_cmp_pd_mask(dover, pd, 1);
		if (flag)
		{
			if (flag & 1) { li += log(m512d_f64(pd, 0)); m512d_f64(pd, 0) = 1.0; }
			if (flag & 2) { li += log(m512d_f64(pd, 1)); m512d_f64(pd, 1) = 1.0; }
			if (flag & 4) { li += log(m512d_f64(pd, 2)); m512d_f64(pd, 2) = 1.0; }
			if (flag & 8) { li += log(m512d_f64(pd, 3)); m512d_f64(pd, 3) = 1.0; }
			if (flag & 16) { li += log(m512d_f64(pd, 4)); m512d_f64(pd, 4) = 1.0; }
			if (flag & 32) { li += log(m512d_f64(pd, 5)); m512d_f64(pd, 5) = 1.0; }
			if (flag & 64) { li += log(m512d_f64(pd, 6)); m512d_f64(pd, 6) = 1.0; }
			if (flag & 128) { li += log(m512d_f64(pd, 7)); m512d_f64(pd, 7) = 1.0; }
		}
	}
	for (int64 j = 0; j < 8; ++j)
	{
		li2 *= m512d_f64(pd, j);
		if (li2 < DOUBLE_UNDERFLOW || li2 > DOUBLE_OVERFLOW)
		{
			li += log(li2);
			li2 = 1.0;
		}
	}
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

TARGET512 uint64 CountNonZero512(byte *A, int64 n)
{
	uint64 re = 0;
	int64 i = 0;
	__m512i z = _mm512_setzero_si512();
	for (int64 l1 = n - 64; i <= l1; i += 64, A += 64)
		re += _mm_popcnt_u64(_mm512_cmpgt_epu8_mask(_mm512_loadu_si512(A), z));
	for (; i < n; ++i, ++A)
		if (*A) re++;
	return re;
}

TARGET512 double Sum512(double *A, int64 n)
{
	int64 i = 0;
	double re = 0;
	__m512d s = _mm512_setzero_pd();
	for (int64 l1 = n - 8; i <= l1; i += 8, A += 8)
		s = _mm512_add_pd(s, _mm512_loadu_pd(A));
	re = _mm512_reduce_add_pd(s);
	for (; i < n; ++i)
		re += *A++;
	return re;
}

TARGET512 double Sum512(byte *A, int64 n)
{
	uint64 re = 0;
	int64 i = 0;
	__m512i s = _mm512_setzero_si512(), z = _mm512_setzero_si512();
	for (int64 l1 = n - 64; i <= l1; i += 64, A += 64)
		s = _mm512_add_epi64(s, _mm512_sad_epu8(_mm512_loadu_si512((__m512i*)A), z));
	re += _mm512_reduce_add_epi64(s);
	for (; i < n; ++i)
		re += *A++;
	return (double)re;
}

TARGET512 double Sum512(double *A, int64 n, int64 sep)
{
	int64 i = 0;
	double re = 0;
	__m512d s = _mm512_setzero_pd();
	__m512i vindex = _mm512_set_epi64(sep * 7, sep * 6, sep * 5, sep * 4, sep * 3, sep * 2, sep, 0);
	for (int64 l1 = n - 8, a8 = sep << 3; i <= l1; i += 8, A += a8)
		s = _mm512_add_pd(s, _mm512_i64gather_pd(vindex, A, 8));
	re = _mm512_reduce_add_pd(s);
	for (; i < n; ++i, A += sep)
		re += *A;
	return re;
}

TARGET512 void Sum512(double *A, double **B, int64 k, int64 n)
{
	int64 i = 0;
	for (int64 l1 = n - 8; i <= l1; i += 8, A += 8)
	{
		__m512d a = _mm512_loadu_pd(B[0] + i);
		for (int64 j = 1; j < k; ++j)
			a = _mm512_add_pd(a, _mm512_loadu_pd(B[j] + i));
		_mm512_storeu_pd(A, a);
	}
	for (; i < n; ++i, A++)
	{
		*A = B[0][i];
		for (int64 j = 1; j < k; ++j)
			*A += B[j][i];
	}
}

TARGET512 double Prod512(double *A, int64 n)
{
	int64 i = 0;
	double re = 1;
	__m512d pd = _mm512_set1_pd(1.0);
	for (int64 l1 = n - 8; i <= l1; i += 8, A += 8)
		pd = _mm512_mul_pd(pd, _mm512_loadu_pd(A));
	re = _mm512_reduce_mul_pd(pd);
	for (; i < n; ++i)
		re *= *A++;
	return re;
}

TARGET512 double Prod512(double *A, int64 n, int64 sep)
{
	int64 i = 0;
	double re = 1;
	__m512d pd = _mm512_set1_pd(1.0);
	__m512i vindex = _mm512_set_epi64(sep * 7, sep * 6, sep * 5, sep * 4, sep * 3, sep * 2, sep, 0);
	for (int64 l1 = n - 8, a8 = sep << 3; i <= l1; i += 8, A += a8)
		pd = _mm512_mul_pd(pd, _mm512_i64gather_pd(vindex, A, 8));
	re = _mm512_reduce_mul_pd(pd);
	for (; i < n; ++i, A += sep)
		re *= *A;
	return re;
}

TARGET512 double SumSquare512(double *A, int64 n)
{
	int64 i = 0;
	double re = 0;
	__m512d s = _mm512_setzero_pd();
	for (int64 l1 = n - 8; i <= l1; i += 8, A += 8)
	{
		__m512d a = _mm512_loadu_pd(A);
		s = _mm512_fmadd_pd(a, a, s);
	}
	re = _mm512_reduce_add_pd(s);
	for (; i < n; ++i, ++A)
		re += *A * *A;
	return re;
}

TARGET512 double SumSquare512(byte *A, int64 n)
{
	int64 i = 0;
	uint64 re = 0;
	__m512i a; 
	uint t0;
	for (int64 l1 = n - 64; i <= l1; i += 64, A += 64)
	{
		a = _mm512_loadu_si512((__m512i*)A);
		t0 = _mm512_reduce_add_epi32(_mm512_maddubs_epi16(a, a));
		re += (t0 & 0xFFFF) + (t0 >> 16);
	}
	for (; i < n; ++i, ++A)
		re += *A * *A;
	return (double)re;
}

TARGET512 void SumSumSquare512(double *A, int64 n, double &sum, double &sumsq)
{
	int64 i = 0;
	sum = sumsq = 0;
	__m512d s = _mm512_setzero_pd(), sq = _mm512_setzero_pd();
	for (int64 l1 = n - 8; i <= l1; i += 8, A += 8)
	{
		__m512d a = _mm512_loadu_pd(A);
		s = _mm512_add_pd(a, s);
		sq = _mm512_fmadd_pd(a, a, sq);
	}
	sum = _mm512_reduce_add_pd(s);
	sumsq = _mm512_reduce_add_pd(sq);
	for (; i < n; ++i, ++A)
	{
		sum += *A;
		sumsq += *A * *A;
	}
}

TARGET512 double SumProd512(double *A, double *B, int64 sep, int64 n)
{
	int64 i = 0;
	double re = 0;
	__m512d s = _mm512_setzero_pd();
	__m512i vindex = _mm512_set_epi64(sep * 7, sep * 6, sep * 5, sep * 4, sep * 3, sep * 2, sep, 0);
	for (int64 l1 = n - 8, b8 = sep << 3; i <= l1; i += 8, A += 8, B += b8)
		s = _mm512_fmadd_pd(_mm512_loadu_pd(A), _mm512_i64gather_pd(vindex, B, 8), s);
	re = _mm512_reduce_add_pd(s);
	for (; i < n; ++i, A++, B += sep)
		re += *A * *B;
	return re;
}

TARGET512 double SumProd512(double *A, double *B, int64 n)
{
	int64 i = 0;
	double re = 0;
	__m512d s = _mm512_setzero_pd();
	for (int64 l1 = n - 8; i <= l1; i += 8, A += 8, B += 8)
		s = _mm512_fmadd_pd(_mm512_loadu_pd(A), _mm512_loadu_pd(B), s);
	re = _mm512_reduce_add_pd(s);
	for (; i < n; ++i, ++A, ++B)
		re += *A * *B;
	return re;
}

TARGET512 void Add512(double *A, double *B, int64 n)
{
	int64 i = 0;
	for (int64 l1 = n - 8; i <= l1; i += 8, A += 8, B += 8)
		_mm512_storeu_pd(A, _mm512_add_pd(_mm512_loadu_pd(A), _mm512_loadu_pd(B)));
	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGET512 void Add512(int *A, int *B, int64 n)
{
	int64 i = 0;
	for (int64 l1 = n - 16; i <= l1; i += 16, A += 16, B += 16)
		_mm512_storeu_si512((__m512i*)A, _mm512_add_epi32(_mm512_loadu_si512((__m512i*)A), _mm512_loadu_si512((__m512i*)B)));
	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGET512 void Add512(double *A, double B, int64 n)
{
	int64 i = 0;
	__m512d b = _mm512_set1_pd(B);
	for (int64 l1 = n - 8; i <= l1; i += 8, A += 8)
		_mm512_storeu_pd(A, _mm512_add_pd(_mm512_loadu_pd(A), b));
	for (; i < n; ++i, A++)
		*A += B;
}

TARGET512 void Mul512(double *C, double *A, double *B, int64 n)
{
	int64 i = 0;
	for (int64 l1 = n - 8; i <= l1; i += 8, A += 8, B += 8, C += 8)
		_mm512_storeu_pd(C, _mm512_mul_pd(_mm512_loadu_pd(A), _mm512_loadu_pd(B)));
	for (; i < n; ++i)
		*C++ = *A++ * *B++;
}

TARGET512 void Mul512(double *C, double *A, double B, int64 n)
{
	int64 i = 0;
	__m512d b = _mm512_set1_pd(B);
	for (int64 l1 = n - 8; i <= l1; i += 8, A += 8, C += 8)
		_mm512_storeu_pd(C, _mm512_mul_pd(_mm512_loadu_pd(A), b));
	for (; i < n; ++i)
		*C++ = *A++ * B;
}

TARGET512 void Mul512(double *A, double B, int64 n)
{
	int64 i = 0;
	__m512d b = _mm512_set1_pd(B);
	for (int64 l1 = n - 8; i <= l1; i += 8, A += 8)
		_mm512_storeu_pd(A, _mm512_mul_pd(_mm512_loadu_pd(A), b));
	for (; i < n; ++i)
		*A++ *= B;
}

TARGET512 void AddProd512(double *C, double *A, double *B, int64 n)
{
	int64 i = 0;
	for (int64 l1 = n - 8; i <= l1; i += 8, A += 8, B += 8, C += 8)
		_mm512_storeu_pd(C, _mm512_fmadd_pd(_mm512_loadu_pd(A), _mm512_loadu_pd(B), _mm512_loadu_pd(C)));
	for (; i < n; ++i)
		*C++ += *A++ * *B++;
}

TARGET512 void AddProd512(double *C, double *A, double B, int64 n)
{
	int64 i = 0;
	__m512d b = _mm512_set1_pd(B);
	for (int64 l1 = n - 8; i <= l1; i += 8, A += 8, C += 8)
		_mm512_storeu_pd(C, _mm512_fmadd_pd(_mm512_loadu_pd(A), b, _mm512_loadu_pd(C)));
	for (; i < n; ++i)
		*C++ += *A++ * B;
}

TARGET512 void Unify512(double *A, int64 n)
{
	int64 i = 0;
	double invsum = 1.0 / (Sum512(A, n) + n * MIN_FREQ);
	__m512d minv = _mm512_set1_pd(invsum * MIN_FREQ), invs = _mm512_set1_pd(invsum);
	for (int64 l1 = n - 8; i <= l1; i += 8, A += 8)
		_mm512_storeu_pd(A, _mm512_fmadd_pd(_mm512_loadu_pd(A), invs, minv));
	for (; i < n; ++i, ++A)
		*A = (*A + MIN_FREQ) * invsum;
}

TARGET512 char *StrNextIdx512(char *A, char val, int64 rep, int64 n)
{
	A++; n--;
	int64 i = 0;
	__m512i v = _mm512_set1_epi8(val);
	for (int64 l1 = n - 64; i <= l1; i += 64, A += 64)
	{
		__m512i a = _mm512_loadu_si512(A);
		__mmask64 mask = _mm512_cmpeq_epu8_mask(a, v);
		int64 count = (int64)_mm_popcnt_u64(mask);
		if (rep > count)
		{
			rep -= count;
			continue;
		}
		else for (;; A++, mask >>= 1)
			if ((mask & 1) && !--rep)
				return A;
	}
	for (; i < n; ++i, A++)
		if (*A == val && !--rep)
			return A;
	return NULL;
}

TARGET512 uint64 CountChar512(char *A, char val, int64 n)
{
	uint64 re = 0;
	int64 i = 0;
	__m512i v = _mm512_set1_epi8(val);
	for (int64 l1 = n - 64; i <= l1; i += 64, A += 64)
		re += _mm_popcnt_u64(_mm512_cmpeq_epu8_mask(_mm512_loadu_si512(A), v));
	for (; i < n; ++i, A++)
		if (*A == val) re++;
	return re;
}