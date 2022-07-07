/* AVX Instruction Set Functions */

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
#endif
#include "vcfpop.h"

TARGETAVX uint64 GetMinIdxAVX(double *A, int64 n, double &val)
{
	int64 i = 0;
	val = 1e300;
	uint64 idx = (uint64)-1;
	__m256d mval = _mm256_set1_pd(val), a;
	__m256i midx = _mm256_set1_epi8((char)0xFF), nidx = _mm256_set_epi64x(3, 2, 1, 0), msep = _mm256_set1_epi64x(4), ma;
	for (int64 l1 = n - 4; i <= l1; i += 4, A += 4)
	{
		a = _mm256_loadu_pd(A);
		ma = _mm256_castpd_si256(_mm256_cmp_pd(mval, a, 1));
		mval = _mm256_min_pd(mval, a);
		midx = _mm256_or_si256(_mm256_and_si256(ma, midx), _mm256_andnot_si256(ma, nidx));
		nidx = _mm256_add_epi64(nidx, msep);
	}
	for (int64 j = 0; j < 4; ++j)
	{
		if (m256d_f64(mval, j) > val) continue;
		val = m256d_f64(mval, j);
		idx = m256i_u64(midx, j);
	}
	for (; i < n; ++i, ++A)
	{
		if (*A > val) continue;
		val = *A;
		idx = i;
	}
	return idx;
}

TARGETAVX void GetMinMaxValAVX(double *A, int64 n, double &minv, double &maxv)
{
	int64 i = 0;
	minv = 1e300;
	maxv = -1e300;
	__m256d mval1 = _mm256_set1_pd(minv), mval2 = _mm256_set1_pd(maxv), a;
	for (int64 l1 = n - 4; i <= l1; i += 4, A += 4)
	{
		a = _mm256_loadu_pd(A);
		mval1 = _mm256_min_pd(mval1, a);
		mval2 = _mm256_max_pd(mval2, a);
	}
	for (int64 j = 0; j < 4; ++j)
	{
		if (m256d_f64(mval1, j) < minv)  minv = m256d_f64(mval1, j);
		if (m256d_f64(mval2, j) > maxv)  maxv = m256d_f64(mval2, j);
	}
	for (; i < n; ++i, ++A)
	{
		if (*A < minv) minv = *A;
		if (*A > maxv) maxv = *A;
	}
}

TARGETAVX double GetMaxValAVX(double *A, int64 n)
{
	int64 i = 0;
	double val = -1e300;
	__m256d mval = _mm256_set1_pd(val);
	for (int64 l1 = n - 4; i <= l1; i += 4, A += 4)
		mval = _mm256_max_pd(mval, _mm256_loadu_pd(A));
	for (int64 j = 0; j < 4; ++j)
	{
		if (m256d_f64(mval, j) < val) continue;
		val = m256d_f64(mval, j);
	}
	for (; i < n; ++i, ++A)
	{
		if (*A < val) continue;
		val = *A;
	}
	return val;
}

TARGETAVX double GetMinValAVX(double *A, int64 n)
{
	int64 i = 0;
	double val = 1e300;
	__m256d mval = _mm256_set1_pd(val);
	for (int64 l1 = n - 4; i <= l1; i += 4, A += 4)
		mval = _mm256_min_pd(mval, _mm256_loadu_pd(A));
	for (int64 j = 0; j < 4; ++j)
	{
		if (m256d_f64(mval, j) > val) continue;
		val = m256d_f64(mval, j);
	}
	for (; i < n; ++i, ++A)
	{
		if (*A > val) continue;
		val = *A;
	}
	return val;
}

TARGETAVX void SetValAVX(uint *A, ushort *B, int64 n)
{
	int64 i = 0;
	__m256i t0;
	for (int64 l1 = n - 16; i <= l1; i += 16, A += 16, B += 16)
	{
		_mm256_storeu_si256((__m256i*)(A), _mm256_cvtepu16_epi32(_mm256_castsi256_si128(t0 = _mm256_loadu_si256((__m256i*)(B)))));
		_mm256_storeu_si256((__m256i*)(A + 8), _mm256_cvtepu16_epi32(_mm256_extracti128_si256(t0, 1)));
	}
	for (; i < n; ++i)
		*A++ = *B++;
}

TARGETAVX double LogProdAVX(double *A, int64 n)
{
	int64 i = 0;
	double li = 0, li2 = 1;
	__m256d pd = _mm256_set1_pd(1.0), duner = _mm256_set1_pd(DOUBLE_UNDERFLOW), dover = _mm256_set1_pd(DOUBLE_OVERFLOW);
	__m256i flag;
	for (int64 l1 = n - 4; i <= l1; i += 4, A += 4)
	{
		pd = _mm256_mul_pd(pd, _mm256_loadu_pd(A));
		flag = _mm256_castpd_si256(_mm256_or_pd(_mm256_cmp_pd(pd, duner, _CMP_LT_OQ), _mm256_cmp_pd(dover, pd, _CMP_LT_OQ)));
		if (m256i_u64(flag, 0)) { li += log(m256d_f64(pd, 0)); m256d_f64(pd, 0) = 1.0; }
		if (m256i_u64(flag, 1)) { li += log(m256d_f64(pd, 1)); m256d_f64(pd, 1) = 1.0; }
		if (m256i_u64(flag, 2)) { li += log(m256d_f64(pd, 2)); m256d_f64(pd, 2) = 1.0; }
		if (m256i_u64(flag, 3)) { li += log(m256d_f64(pd, 3)); m256d_f64(pd, 3) = 1.0; }
	}
	for (int64 j = 0; j < 4; ++j)
	{
		li2 *= m256d_f64(pd, j);
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

TARGETAVX double LogProdAVX(double *A, int64 n, int64 sep)
{
	int64 i = 0;
	double li = 0, li2 = 1;
	__m256d pd = _mm256_set1_pd(1.0), duner = _mm256_set1_pd(DOUBLE_UNDERFLOW), dover = _mm256_set1_pd(DOUBLE_OVERFLOW);
	__m256i vindex = _mm256_set_epi64x(3 * sep, 2 * sep, sep, 0), flag;
	for (int64 l1 = n - 4, a4 = sep << 2; i <= l1; i += 4, A += a4)
	{
		pd = _mm256_mul_pd(pd, _mm256_i64gather_pd(A, vindex, 8));
		flag = _mm256_castpd_si256(_mm256_or_pd(_mm256_cmp_pd(pd, duner, _CMP_LT_OQ), _mm256_cmp_pd(dover, pd, _CMP_LT_OQ)));
		if (m256i_u64(flag, 0)) { li += log(m256d_f64(pd, 0)); m256d_f64(pd, 0) = 1.0; }
		if (m256i_u64(flag, 1)) { li += log(m256d_f64(pd, 1)); m256d_f64(pd, 1) = 1.0; }
		if (m256i_u64(flag, 2)) { li += log(m256d_f64(pd, 2)); m256d_f64(pd, 2) = 1.0; }
		if (m256i_u64(flag, 3)) { li += log(m256d_f64(pd, 3)); m256d_f64(pd, 3) = 1.0; }
	}
	for (int64 j = 0; j < 4; ++j)
	{
		li2 *= m256d_f64(pd, j);
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

TARGETAVX double LogProdDivAVX(double *A, double *B, int64 n, int64 sep)
{
	int64 i = 0;
	double li = 0, li2 = 1;
	__m256d pd = _mm256_set1_pd(1.0), duner = _mm256_set1_pd(DOUBLE_UNDERFLOW), dover = _mm256_set1_pd(DOUBLE_OVERFLOW);
	__m256i vindex = _mm256_set_epi64x(3 * sep, 2 * sep, sep, 0), flag;
	for (int64 l1 = n - 4, s4 = sep << 2; i <= l1; i += 4, A += s4, B += s4)
	{
		pd = _mm256_mul_pd(pd, _mm256_div_pd(_mm256_i64gather_pd(A, vindex, 8), _mm256_i64gather_pd(B, vindex, 8)));
		flag = _mm256_castpd_si256(_mm256_or_pd(_mm256_cmp_pd(pd, duner, _CMP_LT_OQ), _mm256_cmp_pd(dover, pd, _CMP_LT_OQ)));
		if (m256i_u64(flag, 0)) { li += log(m256d_f64(pd, 0)); m256d_f64(pd, 0) = 1.0; }
		if (m256i_u64(flag, 1)) { li += log(m256d_f64(pd, 1)); m256d_f64(pd, 1) = 1.0; }
		if (m256i_u64(flag, 2)) { li += log(m256d_f64(pd, 2)); m256d_f64(pd, 2) = 1.0; }
		if (m256i_u64(flag, 3)) { li += log(m256d_f64(pd, 3)); m256d_f64(pd, 3) = 1.0; }
	}
	for (int64 j = 0; j < 4; ++j)
	{
		li2 *= m256d_f64(pd, j);
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

TARGETAVX uint64 CountNonZeroAVX(byte *A, int64 n)
{
	uint64 re = 0;
	int64 i = 0;
	__m256i a = _mm256_setzero_si256(), z = _mm256_setzero_si256();
	for (int64 l1 = n - 32; i <= l1; i += 32, A += 32)
		a = _mm256_add_epi64(a, _mm256_sad_epu8(z, _mm256_sub_epi8(z, _mm256_cmpeq_epi8(_mm256_loadu_si256((__m256i*)A), z))));
	re = i - (m256i_u64(a, 0) + m256i_u64(a, 1) + m256i_u64(a, 2) + m256i_u64(a, 3));
	for (; i < n; ++i, ++A)
		if (*A) re++;
	return re;
}

TARGETAVX double SumAVX(double *A, int64 n)
{
	int64 i = 0;
	double re = 0;
	__m256d s = _mm256_setzero_pd();
	for (int64 l1 = n - 4; i <= l1; i += 4, A += 4)
		s = _mm256_add_pd(s, _mm256_loadu_pd(A));
	s = _mm256_hadd_pd(s, s);
	re = m256d_f64(s, 0) + m256d_f64(s, 2);
	for (; i < n; ++i)
		re += *A++;
	return re;
}

TARGETAVX double SumAVX(byte *A, int64 n)
{
	uint64 re = 0;
	int64 i = 0;
	__m256i s = _mm256_setzero_si256(), z = _mm256_setzero_si256();
	for (int64 l1 = n - 32; i <= l1; i += 32, A += 32)
		s = _mm256_add_epi64(s, _mm256_sad_epu8(_mm256_loadu_si256((__m256i*)A), z));
	re += m256i_u64(s, 0) + m256i_u64(s, 1) + m256i_u64(s, 2) + m256i_u64(s, 3);
	for (; i < n; ++i)
		re += *A++;
	return (double)re;
}

TARGETAVX double SumAVX(double *A, int64 n, int64 sep)
{
	int64 i = 0;
	double re = 0;
	__m256d s = _mm256_setzero_pd();
	__m256i vindex = _mm256_set_epi64x(3 * sep, 2 * sep, sep, 0);
	for (int64 l1 = n - 4, a4 = sep << 2; i <= l1; i += 4, A += a4)
		s = _mm256_add_pd(s, _mm256_i64gather_pd(A, vindex, 8));
	s = _mm256_hadd_pd(s, s);
	re = m256d_f64(s, 0) + m256d_f64(s, 2);
	for (; i < n; ++i, A += sep)
		re += *A;
	return re;
}

TARGETAVX void SumAVX(double *A, double **B, int64 k, int64 n)
{
	int64 i = 0;
	for (int64 l1 = n - 4; i <= l1; i += 4, A += 4)
	{
		__m256d a = _mm256_loadu_pd(B[0] + i);
		for (int64 j = 1; j < k; ++j)
			a = _mm256_add_pd(a, _mm256_loadu_pd(B[j] + i));
		_mm256_storeu_pd(A, a);
	}
	for (; i < n; ++i, A++)
	{
		*A = B[0][i];
		for (int64 j = 1; j < k; ++j)
			*A += B[j][i];
	}
}

TARGETAVX double ProdAVX(double *A, int64 n)
{
	int64 i = 0;
	double re = 1;
	__m256d pd = _mm256_set1_pd(1.0);
	for (int64 l1 = n - 4; i <= l1; i += 4, A += 4)
		pd = _mm256_mul_pd(pd, _mm256_loadu_pd(A));
	re = m256d_f64(pd, 0) * m256d_f64(pd, 1) * m256d_f64(pd, 2) * m256d_f64(pd, 3);
	for (; i < n; ++i)
		re *= *A++;
	return re;
}

TARGETAVX double ProdAVX(double *A, int64 n, int64 sep)
{
	int64 i = 0;
	double re = 1;
	__m256d pd = _mm256_set1_pd(1.0);
	__m256i vindex = _mm256_set_epi64x(3 * sep, 2 * sep, sep, 0);
	for (int64 l1 = n - 4, a4 = sep << 2; i <= l1; i += 4, A += a4)
		pd = _mm256_mul_pd(pd, _mm256_i64gather_pd(A, vindex, 8));
	re = m256d_f64(pd, 0) * m256d_f64(pd, 1) * m256d_f64(pd, 2) * m256d_f64(pd, 3);
	for (; i < n; ++i, A += sep)
		re *= *A;
	return re;
}

TARGETAVX double SumSquareAVX(double *A, int64 n)
{
	int64 i = 0;
	double re = 0;
	__m256d s = _mm256_setzero_pd();
	for (int64 l1 = n - 4; i <= l1; i += 4, A += 4)
	{
		__m256d a = _mm256_loadu_pd(A);
		s = _mm256_fmadd_pd(a, a, s);
	}
	s = _mm256_hadd_pd(s, s);
	re = m256d_f64(s, 0) + m256d_f64(s, 2);
	for (; i < n; ++i, ++A)
		re += *A * *A;
	return re;
}

TARGETAVX double SumSquareAVX(byte *A, int64 n)
{
	int64 i = 0;
	uint64 re = 0;
	__m256i a, s = _mm256_setzero_si256();
	for (int64 l1 = n - 32; i <= l1; i += 32, A += 32)
	{
		a = _mm256_loadu_si256((__m256i*)A);
		a = _mm256_maddubs_epi16(a, a);
		s = _mm256_add_epi32(s, _mm256_add_epi32(_mm256_cvtepu16_epi32(_mm256_castsi256_si128(a)), _mm256_cvtepu16_epi32(_mm256_extractf128_si256(a, 1))));
	}
	s = _mm256_hadd_epi32(s, s);
	s = _mm256_hadd_epi32(s, s);
	re += m256i_u32(s, 0) + m256i_u32(s, 4);
	for (; i < n; ++i, ++A)
		re += *A * *A;
	return (double)re;
}

TARGETAVX void SumSumSquareAVX(double *A, int64 n, double &sum, double &sumsq)
{
	int64 i = 0;
	sum = sumsq = 0;
	__m256d s = _mm256_setzero_pd(), sq = _mm256_setzero_pd();
	for (int64 l1 = n - 4; i <= l1; i += 4, A += 4)
	{
		__m256d a = _mm256_loadu_pd(A);
		s = _mm256_add_pd(a, s);
		sq = _mm256_fmadd_pd(a, a, sq);
	}
	s = _mm256_hadd_pd(s, s);
	sum = m256d_f64(s, 0) + m256d_f64(s, 2);
	sq = _mm256_hadd_pd(sq, sq);
	sumsq = m256d_f64(sq, 0) + m256d_f64(sq, 2);
	for (; i < n; ++i, ++A)
	{
		sum += *A;
		sumsq += *A * *A;
	}
}

TARGETAVX double SumProdAVX(double *A, double *B, int64 sep, int64 n)
{
	int64 i = 0;
	double re = 0;
	__m256d s = _mm256_setzero_pd();
	__m256i vindex = _mm256_set_epi64x(3 * sep, 2 * sep, sep, 0);
	for (int64 l1 = n - 4, b4 = sep << 2; i <= l1; i += 4, A += 4, B += b4)
		s = _mm256_fmadd_pd(_mm256_loadu_pd(A), _mm256_i64gather_pd(B, vindex, 8), s);
	s = _mm256_hadd_pd(s, s);
	re = m256d_f64(s, 0) + m256d_f64(s, 2);
	for (; i < n; ++i, A++, B += sep)
		re += *A * *B;
	return re;
}

TARGETAVX double SumProdAVX(double *A, double *B, int64 n)
{
	int64 i = 0;
	double re = 0;
	__m256d s = _mm256_setzero_pd();
	for (int64 l1 = n - 4; i <= l1; i += 4, A += 4, B += 4)
		s = _mm256_fmadd_pd(_mm256_loadu_pd(A), _mm256_loadu_pd(B), s);
	s = _mm256_hadd_pd(s, s);
	re = m256d_f64(s, 0) + m256d_f64(s, 2);
	for (; i < n; ++i, ++A, ++B)
		re += *A * *B;
	return re;
}

TARGETAVX void AddAVX(double *A, double *B, int64 n)
{
	int64 i = 0;
	for (int64 l1 = n - 4; i <= l1; i += 4, A += 4, B += 4)
		_mm256_storeu_pd(A, _mm256_add_pd(_mm256_loadu_pd(A), _mm256_loadu_pd(B)));
	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETAVX void AddAVX(int *A, int *B, int64 n)
{
	int64 i = 0;
	for (int64 l1 = n - 8; i <= l1; i += 8, A += 8, B += 8)
		_mm256_storeu_si256((__m256i*)A, _mm256_add_epi32(_mm256_loadu_si256((__m256i*)A), _mm256_loadu_si256((__m256i*)B)));
	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETAVX void AddAVX(double *A, double B, int64 n)
{
	int64 i = 0;
	__m256d b = _mm256_set1_pd(B);
	for (int64 l1 = n - 4; i <= l1; i += 4, A += 4)
		_mm256_storeu_pd(A, _mm256_add_pd(_mm256_loadu_pd(A), b));
	for (; i < n; ++i, A++)
		*A += B;
}

TARGETAVX void MulAVX(double *C, double *A, double *B, int64 n)
{
	int64 i = 0;
	for (int64 l1 = n - 4; i <= l1; i += 4, A += 4, B += 4, C += 4)
		_mm256_storeu_pd(C, _mm256_mul_pd(_mm256_loadu_pd(A), _mm256_loadu_pd(B)));
	for (; i < n; ++i)
		*C++ = *A++ * *B++;
}

TARGETAVX void MulAVX(double *C, double *A, double B, int64 n)
{
	int64 i = 0;
	__m256d b = _mm256_set1_pd(B);
	for (int64 l1 = n - 4; i <= l1; i += 4, A += 4, C += 4)
		_mm256_storeu_pd(C, _mm256_mul_pd(_mm256_loadu_pd(A), b));
	for (; i < n; ++i)
		*C++ = *A++ * B;
}

TARGETAVX void MulAVX(double *A, double B, int64 n)
{
	int64 i = 0;
	__m256d b = _mm256_set1_pd(B);
	for (int64 l1 = n - 4; i <= l1; i += 4, A += 4)
		_mm256_storeu_pd(A, _mm256_mul_pd(_mm256_loadu_pd(A), b));
	for (; i < n; ++i)
		*A++ *= B;
}

TARGETAVX void AddProdAVX(double *C, double *A, double *B, int64 n)
{
	int64 i = 0;
	for (int64 l1 = n - 4; i <= l1; i += 4, A += 4, B += 4, C += 4)
		_mm256_storeu_pd(C, _mm256_fmadd_pd(_mm256_loadu_pd(A), _mm256_loadu_pd(B), _mm256_loadu_pd(C)));
	for (; i < n; ++i)
		*C++ += *A++ * *B++;
}

TARGETAVX void AddProdAVX(double *C, double *A, double B, int64 n)
{
	int64 i = 0;
	__m256d b = _mm256_set1_pd(B);
	for (int64 l1 = n - 4; i <= l1; i += 4, A += 4, C += 4)
		_mm256_storeu_pd(C, _mm256_fmadd_pd(_mm256_loadu_pd(A), b, _mm256_loadu_pd(C)));
	for (; i < n; ++i)
		*C++ += *A++ * B;
}

TARGETAVX void UnifyAVX(double *A, int64 n)
{
	int64 i = 0;
	double invsum = 1.0 / (SumAVX(A, n) + n * MIN_FREQ);
	__m256d minv = _mm256_set1_pd(invsum * MIN_FREQ), invs = _mm256_set1_pd(invsum);
	for (int64 l1 = n - 4; i <= l1; i += 4, A += 4)
		_mm256_storeu_pd(A, _mm256_fmadd_pd(_mm256_loadu_pd(A), invs, minv));
	for (; i < n; ++i, ++A)
		*A = (*A + MIN_FREQ) * invsum;
}

TARGETAVX char *StrNextIdxAVX(char *A, char val, int64 rep, int64 n)
{
	A++; n--;
	int64 i = 0;
	__m256i r = _mm256_setzero_si256(), v = _mm256_set1_epi8(val), o = _mm256_setzero_si256();
	for (int64 l1 = n - 32; i <= l1; )
	{
		__m256i a = _mm256_loadu_si256((__m256i*)A);
		r = _mm256_cmpeq_epi8(a, v);
		if (_mm256_testz_si256(r, r)) { A += 32; i += 32; continue; }
		r = _mm256_sad_epu8(o, _mm256_sub_epi8(o, r)); //4 int 64, each 8 bytes
		for (int64 j = 0; j < 4; ++j, A += 8, i += 8)
			if (rep > (int64)m256i_u64(r, j))
				rep -= m256i_u64(r, j);
			else for (;; A++)
				if (*A == val && !--rep) return A;
	}
	for (; i < n; ++i, A++)
		if (*A == val && !--rep)
			return A;
	return NULL;
}

TARGETAVX uint64 CountCharAVX(char *A, char val, int64 n)
{
	uint64 re = 0;
	int64 i = 0;
	__m256i r = _mm256_setzero_si256(), v = _mm256_set1_epi8(val), o = _mm256_setzero_si256();
	for (int64 l1 = n - 32; i <= l1; i += 32, A += 32)
		r = _mm256_add_epi64(r, _mm256_sad_epu8(o, _mm256_sub_epi8(o, _mm256_cmpeq_epi8(_mm256_loadu_si256((__m256i*)A), v))));
	re = m256i_u64(r, 0) + m256i_u64(r, 1) + m256i_u64(r, 2) + m256i_u64(r, 3);
	for (; i < n; ++i, A++)
		if (*A == val) re++;
	return re;
}