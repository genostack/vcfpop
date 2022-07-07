/* SSE Instruction Set Functions */

#ifndef __SSE3__
#define __SSE3__
#define __SSE4_1__
#define __SSE4_2__
#define __POPCNT__
#define __LZCNT__
#endif
#include "vcfpop.h"

TARGETSSE uint64 GetMinIdxSSE(double *A, int64 n, double &val)
{
	int64 i = 0;
	val = 1e300;
	uint64 idx = (uint64)-1;
	__m128d mval = _mm_set1_pd(val), a;
	__m128i midx = _mm_set1_epi8((char)0xFF), nidx = _mm_set_epi64x(1, 0), msep = _mm_set1_epi64x(2), ma;
	for (int64 l1 = n - 2; i <= l1; i += 2, A += 2)
	{
		a = _mm_loadu_pd(A);
		ma = _mm_castpd_si128(_mm_cmplt_pd(mval, a));
		mval = _mm_min_pd(mval, a);
		midx = _mm_or_si128(_mm_and_si128(ma, midx), _mm_andnot_si128(ma, nidx));
		nidx = _mm_add_epi64(nidx, msep);
	}
	for (int64 j = 0; j < 2; ++j)
	{
		if (m128d_f64(mval, j) > val) continue;
		val = m128d_f64(mval, j);
		idx = m128i_u64(midx, j);
	}
	for (; i < n; ++i, ++A)
	{
		if (*A > val) continue;
		val = *A;
		idx = i;
	}
	return idx;
}

TARGETSSE void GetMinMaxValSSE(double *A, int64 n, double &minv, double &maxv)
{
	int64 i = 0;
	minv = 1e300;
	maxv = -1e300;
	__m128d mval1 = _mm_set1_pd(minv), mval2 = _mm_set1_pd(maxv), a;
	for (int64 l1 = n - 2; i <= l1; i += 2, A += 2)
	{
		a = _mm_loadu_pd(A);
		mval1 = _mm_min_pd(mval1, a);
		mval2 = _mm_max_pd(mval2, a);
	}
	for (int64 j = 0; j < 2; ++j)
	{
		if (m128d_f64(mval1, j) < minv) minv = m128d_f64(mval1, j);
		if (m128d_f64(mval2, j) > maxv) maxv = m128d_f64(mval2, j);
	}
	for (; i < n; ++i, ++A)
	{
		if (*A < minv) minv = *A;
		if (*A > maxv) maxv = *A;
	}
}

TARGETSSE double GetMaxValSSE(double *A, int64 n)
{
	int64 i = 0;
	double val = -1e300;
	__m128d mval = _mm_set1_pd(val);
	for (int64 l1 = n - 2; i <= l1; i += 2, A += 2)
		mval = _mm_max_pd(mval, _mm_loadu_pd(A));
	for (int64 j = 0; j < 2; ++j)
	{
		if (m128d_f64(mval, j) < val) continue;
		val = m128d_f64(mval, j);
	}
	for (; i < n; ++i, ++A)
	{
		if (*A < val) continue;
		val = *A;
	}
	return val;
}

TARGETSSE double GetMinValSSE(double *A, int64 n)
{
	int64 i = 0;
	double val = 1e300;
	__m128d mval = _mm_set1_pd(val);
	for (int64 l1 = n - 2; i <= l1; i += 2, A += 2)
		mval = _mm_min_pd(mval, _mm_loadu_pd(A));
	for (int64 j = 0; j < 2; ++j)
	{
		if (m128d_f64(mval, j) > val) continue;
		val = m128d_f64(mval, j);
	}
	for (; i < n; ++i, ++A)
	{
		if (*A > val) continue;
		val = *A;
	}
	return val;
}

TARGETSSE void SetValSSE(uint *A, ushort *B, int64 n)
{
	int64 i = 0;
	__m128i t0;
	for (int64 l1 = n - 8; i <= l1; i += 8, A += 8, B += 8)
	{
		_mm_storeu_si128((__m128i*)A, _mm_cvtepu16_epi32(t0 = _mm_loadu_si128((__m128i*)B)));
		_mm_storeu_si128((__m128i*)(A + 4), _mm_cvtepu16_epi32(_mm_srli_si128(t0, 8)));
	}
	for (; i < n; ++i)
		*A++ = *B++;
}

TARGETSSE double LogProdSSE(double *A, int64 n)
{
	int64 i = 0;
	double li = 0, li2 = 1;
	__m128d pd = _mm_set1_pd(1.0), duner = _mm_set1_pd(DOUBLE_UNDERFLOW), dover = _mm_set1_pd(DOUBLE_OVERFLOW);
	__m128i flag;
	for (int64 l1 = n - 2; i <= l1; i += 2, A += 2)
	{
		pd = _mm_mul_pd(pd, _mm_loadu_pd(A));
		flag = _mm_castpd_si128(_mm_or_pd(_mm_cmplt_pd(pd, duner), _mm_cmplt_pd(dover, pd)));
		if (m128i_u64(flag, 0)) { li += log(m128d_f64(pd, 0)); m128d_f64(pd, 0) = 1.0; }
		if (m128i_u64(flag, 1)) { li += log(m128d_f64(pd, 1)); m128d_f64(pd, 1) = 1.0; }
	}
	for (int64 j = 0; j < 2; ++j)
	{
		li2 *= m128d_f64(pd, j);
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

TARGETSSE double LogProdSSE(double *A, int64 n, int64 sep)
{
	int64 i = 0;
	double li = 0, li2 = 1;
	__m128d pd = _mm_set1_pd(1.0), duner = _mm_set1_pd(DOUBLE_UNDERFLOW), dover = _mm_set1_pd(DOUBLE_OVERFLOW);
	__m128i flag;
	for (int64 l1 = n - 2, a2 = sep << 1; i <= l1; i += 2, A += a2)
	{
		pd = _mm_mul_pd(pd, _mm_set_pd(A[sep], A[0]));
		flag = _mm_castpd_si128(_mm_or_pd(_mm_cmplt_pd(pd, duner), _mm_cmplt_pd(dover, pd)));
		if (m128i_u64(flag, 0)) { li += log(m128d_f64(pd, 0)); m128d_f64(pd, 0) = 1.0; }
		if (m128i_u64(flag, 1)) { li += log(m128d_f64(pd, 1)); m128d_f64(pd, 1) = 1.0; }
	}
	for (int64 j = 0; j < 2; ++j)
	{
		li2 *= m128d_f64(pd, j);
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

TARGETSSE double LogProdDivSSE(double *A, double *B, int64 n, int64 sep)
{
	int64 i = 0;
	double li = 0, li2 = 1;
	__m128d pd = _mm_set1_pd(1.0), duner = _mm_set1_pd(DOUBLE_UNDERFLOW), dover = _mm_set1_pd(DOUBLE_OVERFLOW);
	__m128i flag;
	for (int64 l1 = n - 2, s2 = sep << 1; i <= l1; i += 2, A += s2, B += s2)
	{
		pd = _mm_mul_pd(pd, _mm_div_pd(_mm_set_pd(A[sep], A[0]), _mm_set_pd(B[sep], B[0])));
		flag = _mm_castpd_si128(_mm_or_pd(_mm_cmplt_pd(pd, duner), _mm_cmplt_pd(dover, pd)));
		if (m128i_u64(flag, 0)) { li += log(m128d_f64(pd, 0)); m128d_f64(pd, 0) = 1.0; }
		if (m128i_u64(flag, 1)) { li += log(m128d_f64(pd, 1)); m128d_f64(pd, 1) = 1.0; }
	}
	for (int64 j = 0; j < 2; ++j)
	{
		li2 *= m128d_f64(pd, j);
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

TARGETSSE uint64 CountNonZeroSSE(byte *A, int64 n)
{
	uint64 re = 0;
	int64 i = 0;
	__m128i a = _mm_setzero_si128(), z = _mm_setzero_si128();
	for (int64 l1 = n - 16; i <= l1; i += 16, A += 16)
		a = _mm_add_epi64(a, _mm_sad_epu8(z, _mm_sub_epi8(z, _mm_cmpeq_epi8(_mm_loadu_si128((__m128i*)A), z))));
	re = i - (m128i_u64(a, 0) + m128i_u64(a, 1));
	for (; i < n; ++i, ++A)
		if (*A) re++;
	return re;
}

TARGETSSE double SumSSE(double *A, int64 n)
{
	int64 i = 0;
	double re = 0;
	__m128d s = _mm_setzero_pd();
	for (int64 l1 = n - 2; i <= l1; i += 2, A += 2)
		s = _mm_add_pd(s, _mm_loadu_pd(A));
	s = _mm_hadd_pd(s, s);
	re = m128d_f64(s, 0);
	for (; i < n; ++i)
		re += *A++;
	return re;
}

TARGETSSE double SumSSE(byte *A, int64 n)
{
	uint64 re = 0;
	int64 i = 0;
	__m128i s = _mm_setzero_si128(), z = _mm_setzero_si128();
	for (int64 l1 = n - 16; i <= l1; i += 16, A += 16)
		s = _mm_add_epi64(s, _mm_sad_epu8(_mm_loadu_si128((__m128i*)A), z));
	re += m128i_u64(s, 0) + m128i_u64(s, 1);
	for (; i < n; ++i)
		re += *A++;
	return (double)re;
}

TARGETSSE double SumSSE(double *A, int64 n, int64 sep)
{
	int64 i = 0;
	double re = 0;
	__m128d s = _mm_setzero_pd();
	for (int64 l1 = n - 2, a2 = sep << 1; i <= l1; i += 2, A += a2)
		s = _mm_add_pd(s, _mm_set_pd(A[sep], A[0]));
	s = _mm_hadd_pd(s, s);
	re = m128d_f64(s, 0);
	for (; i < n; ++i, A += sep)
		re += *A;
	return re;
}

TARGETSSE void SumSSE(double *A, double **B, int64 k, int64 n)
{
	int64 i = 0;
	for (int64 l1 = n - 2; i <= l1; i += 2, A += 2)
	{
		__m128d a = _mm_loadu_pd(B[0] + i);
		for (int64 j = 1; j < k; ++j)
			a = _mm_add_pd(a, _mm_loadu_pd(B[j] + i));
		_mm_storeu_pd(A, a);
	}
	for (; i < n; ++i, A++)
	{
		*A = B[0][i];
		for (int64 j = 1; j < k; ++j)
			*A += B[j][i];
	}
}

TARGETSSE double ProdSSE(double *A, int64 n)
{
	int64 i = 0;
	double re = 1;
	__m128d pd = _mm_set1_pd(1.0);
	for (int64 l1 = n - 2; i <= l1; i += 2, A += 2)
		pd = _mm_mul_pd(pd, _mm_loadu_pd(A));
	re = m128d_f64(pd, 0) * m128d_f64(pd, 1);
	for (; i < n; ++i)
		re *= *A++;
	return re;
}

TARGETSSE double ProdSSE(double *A, int64 n, int64 sep)
{
	int64 i = 0;
	double re = 1;
	__m128d pd = _mm_set1_pd(1.0);
	for (int64 l1 = n - 2, a2 = sep << 1; i <= l1; i += 2, A += a2)
		pd = _mm_mul_pd(pd, _mm_set_pd(A[sep], A[0]));
	re = m128d_f64(pd, 0) * m128d_f64(pd, 1);
	for (; i < n; ++i, A += sep)
		re *= *A;
	return re;
}

TARGETSSE double SumSquareSSE(double *A, int64 n)
{
	int64 i = 0;
	double re = 0;
	__m128d s = _mm_setzero_pd();
	for (int64 l1 = n - 2; i <= l1; i += 2, A += 2)
	{
		__m128d a = _mm_loadu_pd(A);
		s = _mm_add_pd(s, _mm_mul_pd(a, a));
	}
	s = _mm_hadd_pd(s, s);
	re = m128d_f64(s, 0);
	for (; i < n; ++i, ++A)
		re += *A * *A;
	return re;
}

TARGETSSE double SumSquareSSE(byte *A, int64 n)
{
	int64 i = 0;
	uint64 re = 0;
	__m128i a, s = _mm_setzero_si128();
	for (int64 l1 = n - 16; i <= l1; i += 16, A += 16)
	{
		a = _mm_loadu_si128((__m128i*)A);
		a = _mm_maddubs_epi16(a, a);
		s = _mm_add_epi32(s, _mm_add_epi32(_mm_cvtepu8_epi16(a), _mm_cvtepu8_epi16(_mm_srli_si128(a, 8))));
	}
	s = _mm_hadd_epi32(s, s);
	s = _mm_hadd_epi32(s, s);
	re += m128i_u32(s, 0);
	for (; i < n; ++i, ++A)
		re += *A * *A;
	return (double)re;
}

TARGETSSE void SumSumSquareSSE(double *A, int64 n, double &sum, double &sumsq)
{
	int64 i = 0;
	sum = sumsq = 0;
	__m128d s = _mm_setzero_pd(), sq = _mm_setzero_pd();
	for (int64 l1 = n - 2; i <= l1; i += 2, A += 2)
	{
		__m128d a = _mm_loadu_pd(A);
		s = _mm_add_pd(a, s);
		sq = _mm_add_pd(sq, _mm_mul_pd(a, a));
	}
	s = _mm_hadd_pd(s, s);
	sum = m128d_f64(s, 0);
	sq = _mm_hadd_pd(sq, sq);
	sumsq = m128d_f64(sq, 0);
	for (; i < n; ++i, ++A)
	{
		sum += *A;
		sumsq += *A * *A;
	}
}

TARGETSSE double SumProdSSE(double *A, double *B, int64 sep, int64 n)
{
	int64 i = 0;
	double re = 0;
	__m128d s = _mm_setzero_pd();
	for (int64 l1 = n - 2, b2 = sep << 1; i <= l1; i += 2, A += 2, B += b2)
		s = _mm_add_pd(s, _mm_mul_pd(_mm_loadu_pd(A), _mm_set_pd(B[sep], B[0])));
	s = _mm_hadd_pd(s, s);
	re = m128d_f64(s, 0);
	for (; i < n; ++i, A++, B += sep)
		re += *A * *B;
	return re;
}

TARGETSSE double SumProdSSE(double *A, double *B, int64 n)
{
	int64 i = 0;
	double re = 0;
	__m128d s = _mm_setzero_pd();
	for (int64 l1 = n - 2; i <= l1; i += 2, A += 2, B += 2)
		s = _mm_add_pd(s,
			_mm_mul_pd(_mm_loadu_pd(A), _mm_loadu_pd(B)));
	s = _mm_hadd_pd(s, s);
	re = m128d_f64(s, 0);
	for (; i < n; ++i, ++A, ++B)
		re += *A * *B;
	return re;
}

TARGETSSE void AddSSE(double *A, double *B, int64 n)
{
	int64 i = 0;
	for (int64 l1 = n - 2; i <= l1; i += 2, A += 2, B += 2)
		_mm_storeu_pd(A, _mm_add_pd(_mm_loadu_pd(A), _mm_loadu_pd(B)));
	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETSSE void AddSSE(int *A, int *B, int64 n)
{
	int64 i = 0;
	for (int64 l1 = n - 4; i <= l1; i += 4, A += 4, B += 4)
		_mm_storeu_si128((__m128i*)A, _mm_add_epi32(_mm_loadu_si128((__m128i*)A), _mm_loadu_si128((__m128i*)B)));
	for (; i < n; ++i, A++, B++)
		*A += *B;
}

TARGETSSE void AddSSE(double *A, double B, int64 n)
{
	int64 i = 0;
	__m128d b = _mm_set1_pd(B);
	for (int64 l1 = n - 2; i <= l1; i += 2, A += 2)
		_mm_storeu_pd(A, _mm_add_pd(_mm_loadu_pd(A), b));
	for (; i < n; ++i, A++)
		*A += B;
}

TARGETSSE void MulSSE(double *C, double *A, double *B, int64 n)
{
	int64 i = 0;
	for (int64 l1 = n - 2; i <= l1; i += 2, A += 2, B += 2, C += 2)
		_mm_storeu_pd(C, _mm_mul_pd(_mm_loadu_pd(A), _mm_loadu_pd(B)));
	for (; i < n; ++i)
		*C++ = *A++ * *B++;
}

TARGETSSE void MulSSE(double *C, double *A, double B, int64 n)
{
	int64 i = 0;
	__m128d b = _mm_set1_pd(B);
	for (int64 l1 = n - 2; i <= l1; i += 2, A += 2, C += 2)
		_mm_storeu_pd(C, _mm_mul_pd(_mm_loadu_pd(A), b));
	for (; i < n; ++i)
		*C++ = *A++ * B;
}

TARGETSSE void MulSSE(double *A, double B, int64 n)
{
	int64 i = 0;
	__m128d b = _mm_set1_pd(B);
	for (int64 l1 = n - 2; i <= l1; i += 2, A += 2)
		_mm_storeu_pd(A, _mm_mul_pd(_mm_loadu_pd(A), b));
	for (; i < n; ++i)
		*A++ *= B;
}

TARGETSSE void AddProdSSE(double *C, double *A, double *B, int64 n)
{
	int64 i = 0;
	for (int64 l1 = n - 2; i <= l1; i += 2, A += 2, B += 2, C += 2)
		_mm_storeu_pd(C, _mm_add_pd(_mm_loadu_pd(C), _mm_mul_pd(_mm_loadu_pd(A), _mm_loadu_pd(B))));
	for (; i < n; ++i)
		*C++ += *A++ * *B++;
}

TARGETSSE void AddProdSSE(double *C, double *A, double B, int64 n)
{
	int64 i = 0;
	__m128d b = _mm_set1_pd(B);
	for (int64 l1 = n - 2; i <= l1; i += 2, A += 2, C += 2)
		_mm_storeu_pd(C, _mm_add_pd(_mm_loadu_pd(C), _mm_mul_pd(_mm_loadu_pd(A), b)));
	for (; i < n; ++i)
		*C++ += *A++ * B;
}

TARGETSSE void UnifySSE(double *A, int64 n)
{
	int64 i = 0;
	double invsum = 1.0 / (SumSSE(A, n) + n * MIN_FREQ);
	__m128d minv = _mm_set1_pd(MIN_FREQ), invs = _mm_set1_pd(invsum);
	for (int64 l1 = n - 2; i <= l1; i += 2, A += 2)
		_mm_storeu_pd(A, _mm_mul_pd(invs, _mm_add_pd(_mm_loadu_pd(A), minv)));
	for (; i < n; ++i, ++A)
		*A = (*A + MIN_FREQ) * invsum;
}

TARGETSSE char *StrNextIdxSSE(char *A, char val, int64 rep, int64 n)
{
	A++; n--;
	int64 i = 0;
	__m128i r = _mm_setzero_si128(), v = _mm_set1_epi8(val), o = _mm_setzero_si128();
	for (int64 l1 = n - 16; i <= l1; )
	{
		__m128i a = _mm_loadu_si128((__m128i*)A);
		r = _mm_cmpeq_epi8(a, v);
		if (_mm_testz_si128(r, r)) { A += 16; i += 16; continue; }
		r = _mm_sad_epu8(o, _mm_sub_epi8(o, r)); //4 int 64, each 8 bytes
		for (int64 j = 0; j < 2; ++j, A += 8, i += 8)
			if (rep > (int64)m128i_u64(r, j))
				rep -= m128i_u64(r, j);
			else for (;; A++)
				if (*A == val && !--rep) return A;
	}
	for (; i < n; ++i, A++)
		if (*A == val && !--rep)
			return A;
	return NULL;
}

TARGETSSE uint64 CountCharSSE(char *A, char val, int64 n)
{
	uint64 re = 0;
	int64 i = 0;
	__m128i r = _mm_setzero_si128(), v = _mm_set1_epi8(val), o = _mm_setzero_si128();
	for (int64 l1 = n - 16; i <= l1; i += 16, A += 16)
		r = _mm_add_epi64(r, _mm_sad_epu8(o, _mm_sub_epi8(o, _mm_cmpeq_epi8(_mm_loadu_si128((__m128i*)A), v))));
	re = m128i_u64(r, 0) + m128i_u64(r, 1);
	for (; i < n; ++i, A++)
		if (*A == val) re++;
	return re;
}