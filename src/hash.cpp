/* Hash Functions */

#include "vcfpop.h"

#ifdef HASH64

/* 64-bit hash to reduce hash collision */

/* Get crc32 value for a genotype */
TARGET HASH HashGenotype(ushort *allele, int ploidy)
{
	uint crc1 = _mm_crc32_u64(0xB385FDA9ull, ploidy), crc2 = _mm_crc32_u64(0x7FED5EA6ull, ploidy);
	switch (ploidy)
	{
	case 1:
		crc1 = _mm_crc32_u64(crc1, (uint64) * (ushort*)allele);
		crc2 = _mm_crc32_u64(crc2, (uint64) * (ushort*)allele);
		break;
	case 2:
		crc1 = _mm_crc32_u64(crc1, (uint64) * (uint*)allele);
		crc2 = _mm_crc32_u64(crc2, (uint64) * (uint*)allele);
		break;
	case 3:
		crc1 = _mm_crc32_u64(crc1, (uint64) * (uint64*)allele & 0xFFFFFFFFFFFFull);
		crc2 = _mm_crc32_u64(crc2, (uint64) * (uint64*)allele & 0xFFFFFFFFFFFFull);
		break;
	case 4:
		crc1 = _mm_crc32_u64(crc1, (uint64) * (uint64*)allele);
		crc2 = _mm_crc32_u64(crc2, (uint64) * (uint64*)allele);
		break;
	case 5:
		crc1 = _mm_crc32_u64(_mm_crc32_u64(crc1, (uint64) * (uint64*)allele),
			(uint64) * (ushort*)(allele + 4));
		crc2 = _mm_crc32_u64(_mm_crc32_u64(crc2, (uint64) * (uint64*)allele),
			(uint64) * (ushort*)(allele + 4));
		break;
	case 6:
		crc1 = _mm_crc32_u64(_mm_crc32_u64(crc1, (uint64) * (uint64*)allele),
			(uint64) * (uint*)(allele + 4));
		crc2 = _mm_crc32_u64(_mm_crc32_u64(crc2, (uint64) * (uint64*)allele),
			(uint64) * (uint*)(allele + 4));
		break;
	case 7:
		crc1 = _mm_crc32_u64(_mm_crc32_u64(crc1, (uint64) * (uint64*)allele),
			(uint64) * (uint64*)(allele + 4) & 0xFFFFFFFFFFFFull);
		crc2 = _mm_crc32_u64(_mm_crc32_u64(crc2, (uint64) * (uint64*)allele),
			(uint64) * (uint64*)(allele + 4) & 0xFFFFFFFFFFFFull);
		break;
	case 8:
		crc1 = _mm_crc32_u64(_mm_crc32_u64(crc1, (uint64) * (uint64*)allele),
			(uint64) * (uint64*)(allele + 4));
		crc2 = _mm_crc32_u64(_mm_crc32_u64(crc2, (uint64) * (uint64*)allele),
			(uint64) * (uint64*)(allele + 4));
		break;
	case 9:
		crc1 = _mm_crc32_u64(_mm_crc32_u64(_mm_crc32_u64(crc1, (uint64) * (uint64*)allele),
			(uint64) * (uint64*)(allele + 4)),
			(uint64) * (ushort*)(allele + 8));
		crc2 = _mm_crc32_u64(_mm_crc32_u64(_mm_crc32_u64(crc2, (uint64) * (uint64*)allele),
			(uint64) * (uint64*)(allele + 4)),
			(uint64) * (ushort*)(allele + 8));
		break;
	case 10:
		crc1 = _mm_crc32_u64(_mm_crc32_u64(_mm_crc32_u64(crc1, (uint64) * (uint64*)allele),
			(uint64) * (uint64*)(allele + 4)),
			(uint64) * (uint*)(allele + 8));
		crc2 = _mm_crc32_u64(_mm_crc32_u64(_mm_crc32_u64(crc2, (uint64) * (uint64*)allele),
			(uint64) * (uint64*)(allele + 4)),
			(uint64) * (uint*)(allele + 8));
		break;
	default: break;
	}
	return ((uint64)crc1 << 32) | (uint64)crc2;
}

/* Get crc32 value for a string */
TARGET HASH HashString(char *str, int len)
{
	uint crc1 = 0xB385FDA9, crc2 = 0x7FED5EA6;
	if (len == -1) len = (int)strlen(str);
	for (; len >= 8; len -= 8, str += 8)
	{
		crc1 = _mm_crc32_u64(crc1, *(uint64*)str);
		crc2 = _mm_crc32_u64(crc2, *(uint64*)str);
	}

	switch (len)
	{
	case 1:
		crc1 = _mm_crc32_u64(crc1, (uint64) * (byte*)str);
		crc2 = _mm_crc32_u64(crc2, (uint64) * (byte*)str);
		break;
	case 2:
		crc1 = _mm_crc32_u64(crc1, (uint64) * (ushort*)str);
		crc2 = _mm_crc32_u64(crc2, (uint64) * (ushort*)str);
		break;
	case 3:
		crc1 = _mm_crc32_u64(crc1, ((uint64)(*(ushort*)str) << 8) | (uint64) * (byte*)(str + 2));
		crc2 = _mm_crc32_u64(crc2, ((uint64)(*(ushort*)str) << 8) | (uint64) * (byte*)(str + 2));
		break;
	case 4:
		crc1 = _mm_crc32_u64(crc1, (uint64) * (uint*)str);
		crc2 = _mm_crc32_u64(crc2, (uint64) * (uint*)str);
		break;
	case 5:
		crc1 = _mm_crc32_u64(crc1, ((uint64)(*(uint*)str) << 8) | (uint64) * (byte*)(str + 4));
		crc2 = _mm_crc32_u64(crc2, ((uint64)(*(uint*)str) << 8) | (uint64) * (byte*)(str + 4));
		break;
	case 6:
		crc1 = _mm_crc32_u64(crc1, ((uint64)(*(uint*)str) << 16) | (uint64) * (ushort*)(str + 4));
		crc2 = _mm_crc32_u64(crc2, ((uint64)(*(uint*)str) << 16) | (uint64) * (ushort*)(str + 4));
		break;
	case 7:
		crc1 = _mm_crc32_u64(crc1, ((uint64)(*(uint*)str) << 24) | ((uint64)(*(ushort*)(str + 4)) << 8) | (uint64)(*(byte*)(str + 6)));
		crc2 = _mm_crc32_u64(crc2, ((uint64)(*(uint*)str) << 24) | ((uint64)(*(ushort*)(str + 4)) << 8) | (uint64)(*(byte*)(str + 6)));
		break;
	default:
		break;
	}
	return ((uint64)crc1 << 32) | (uint64)crc2;
}

/* Get crc32 value for unsigned long integer */
TARGET HASH HashULong(uint64 val)
{
	return ((uint64)_mm_crc32_u64(0xB385FDA9, val) << 32) | (uint64)_mm_crc32_u64(0x7FED5EA6, val);
}

/* Get crc32 value for two unsigned integer */
TARGET HASH HashUInt(uint val, uint val2)
{
	return HashULong(((uint64)val << 32) | (uint64)val2);
}

#else

/* 32-bit hash */

/* Get crc32 value for a genotype */
TARGET HASH HashGenotype(ushort *allele, int ploidy)
{
	HASH crc1 = _mm_crc32_u64(0xB385FDA9ull, ploidy);
	switch (ploidy)
	{
	case 1:
		crc1 = _mm_crc32_u64(crc1, (uint64)*(ushort*)allele);
		break;
	case 2:
		crc1 = _mm_crc32_u64(crc1, (uint64)*(uint*)allele);
		break;
	case 3:
		crc1 = _mm_crc32_u64(crc1, (uint64)*(uint64*)allele & 0xFFFFFFFFFFFFull);
		break;
	case 4:
		crc1 = _mm_crc32_u64(crc1, (uint64)*(uint64*)allele);
		break;
	case 5:
		crc1 = _mm_crc32_u64(_mm_crc32_u64(crc1, (uint64)*(uint64*)allele),
			(uint64)*(ushort*)(allele + 4));
		break;
	case 6:
		crc1 = _mm_crc32_u64(_mm_crc32_u64(crc1, (uint64)*(uint64*)allele),
			(uint64)*(uint*)(allele + 4));
		break;
	case 7:
		crc1 = _mm_crc32_u64(_mm_crc32_u64(crc1, (uint64)*(uint64*)allele),
			(uint64)*(uint64*)(allele + 4) & 0xFFFFFFFFFFFFull);
		break;
	case 8:
		crc1 = _mm_crc32_u64(_mm_crc32_u64(crc1, (uint64)*(uint64*)allele),
			(uint64)*(uint64*)(allele + 4));
		break;
	case 9:
		crc1 = _mm_crc32_u64(_mm_crc32_u64(_mm_crc32_u64(crc1, (uint64)*(uint64*)allele), 
			(uint64)*(uint64*)(allele + 4)),
			(uint64)*(ushort*)(allele + 8));
		break;
	case 10:
		crc1 = _mm_crc32_u64(_mm_crc32_u64(_mm_crc32_u64(crc1, (uint64) * (uint64*)allele),
			(uint64)*(uint64*)(allele + 4)),
			(uint64)*(uint*)(allele + 8));
		break;
	default: break;
	}
	return (uint)crc1;
}

/* Get crc32 value for a string */
TARGET HASH HashString(char *str, int len)
{
	HASH crc1 = 0xB385FDA9;
	if (len == -1) len = (int)strlen(str);

	for (; len >= 8; len -= 8, str += 8)
		crc1 = _mm_crc32_u64(crc1, *(uint64*)str);

	switch (len)
	{
	case 1:
		crc1 = _mm_crc32_u64(crc1, (uint64)*(byte*)str);
		break;
	case 2:
		crc1 = _mm_crc32_u64(crc1, (uint64)*(ushort*)str);
		break;
	case 3:
		crc1 = _mm_crc32_u64(crc1, ((uint64)(*(ushort*)str) << 8) | (uint64)*(byte*)(str + 2));
		break;
	case 4:
		crc1 = _mm_crc32_u64(crc1, (uint64)*(uint*)str);
		break;
	case 5:
		crc1 = _mm_crc32_u64(crc1, ((uint64)(*(uint*)str) << 8) | (uint64)*(byte*)(str + 4));
		break;
	case 6:
		crc1 = _mm_crc32_u64(crc1, ((uint64)(*(uint*)str) << 16) | (uint64)*(ushort*)(str + 4));
		break;
	case 7:
		crc1 = _mm_crc32_u64(crc1, ((uint64)(*(uint*)str) << 24) | ((uint64)(*(ushort*)(str + 4)) << 8) | (uint64)(*(byte*)(str + 6)));
		break;
	default:
		break;
	}
	return (uint)crc1;
}

/* Get crc32 value for unsigned long integer */
TARGET HASH HashULong(uint64 val)
{
	return _mm_crc32_u64(0xB385FDA9, val);
}

/* Get crc32 value for two unsigned integer */
TARGET HASH HashUInt(uint val, uint val2)
{
	return HashULong(((uint64)val << 32) | (uint64)val2);
}

#endif

/* Initialize crypt table for Huang 2015 maximum-likelihood polyploid relatedness estimator */
TARGET void InitCryptTable()
{
	if (cryptTable != NULL)  return;
	cryptTable = new uint[0x500];

	uint dwHih, dwLow, seed = 0x00100001, index1 = 0, index2 = 0, i;
	for (index1 = 0; index1 < 0x100; ++index1)
	{
		for (index2 = index1, i = 0; i < 5; ++i, index2 += 0x100)
		{
			seed = (seed * 125 + 3) % 0x2AAAAB;
			dwHih = (seed & 0xFFFF) << 16;
			seed = (seed * 125 + 3) % 0x2AAAAB;
			dwLow = (seed & 0xFFFF);
			cryptTable[index2] = (dwHih | dwLow);
		}
	}
}

/* Get IBS model of a genotype */
TARGET uint GetSingleIBS(int *x, byte ploidy)
{
	int n = 0, t[8], tt;
	for (int i = 0; i < ploidy; )
	{
		int j = i + 1;
		for (; j < ploidy && x[j] == x[i]; ++j);
		t[n++] = j - i;
		i = j;
	}
	for (int i = 0; i < n; ++i)
		for (int j = i + 1; j < n; ++j)
			if (t[j] < t[i])
				Swap(t[i], t[j]);
	tt = 0;
	for (int i = 0; i < n; ++i)
		tt = tt * 10 + t[i];
	return tt;
}

/* Get hash of IBS mode of a pair of genotypes */
TARGET uint HashString32(char *s1, char *s2, char *s3)
{
	uint dwSeed1 = 0x7FED7FED, dwSeed2 = 0xEEEEEEEE;
	uint *cryptTable2 = (uint*)cryptTable;
	byte b1;
	for (int i = 0; s1[i]; ++i)
	{
		b1 = s1[i];
		dwSeed1 = cryptTable2[b1 + 0x100] ^ (dwSeed1 + dwSeed2);
		dwSeed2 = b1 + dwSeed1 + dwSeed2 + (dwSeed2 << 5) + 3;
	}
	for (int i = 0; s2[i]; ++i)
	{
		b1 = s2[i];
		dwSeed1 = cryptTable2[b1 + 0x100] ^ (dwSeed1 + dwSeed2);
		dwSeed2 = b1 + dwSeed1 + dwSeed2 + (dwSeed2 << 5) + 3;
	}
	for (int i = 0; s3[i]; ++i)
	{
		b1 = s3[i];
		dwSeed1 = cryptTable2[b1 + 0x100] ^ (dwSeed1 + dwSeed2);
		dwSeed2 = b1 + dwSeed1 + dwSeed2 + (dwSeed2 << 5) + 3;
	}
	return dwSeed1;
}

/* Get IBS mode of a pair of genotypes */
TARGET uint GetHuang2015Hash(int *x, int *y, int p)
{
	//xibs + yibs + crosssame + crossfactors
	int crosssame = 0;
	int crfa[8], crfb[8], crfan = 0, crfbn = 0;
	for (int i = 0; i < p; ++i)
	{
		int c1 = 0, c2 = 0;
		for (int j = 0; j < p; ++j)
		{
			if (x[i] == y[j])
			{
				crosssame++;
				c1++;
			}
			if (x[j] == y[i])
				c2++;
		}
		if (c1) crfa[crfan++] = c1;
		if (c2) crfb[crfbn++] = c2;
	}

	for (int i = 0; i < crfan; ++i)
		for (int j = i + 1; j < crfan; ++j)
			if (crfa[i] > crfa[j])
				Swap(crfa[i], crfa[j]);

	for (int i = 0; i < crfbn; ++i)
		for (int j = i + 1; j < crfbn; ++j)
			if (crfb[i] > crfb[j])
				Swap(crfb[i], crfb[j]);

	char s1[20], s2[20], s3[20];
	s2[0] = '\0';
	s3[0] = '\0';

	char *ss1 = s1;
	AppendInt(ss1, GetSingleIBS(x, p)); *ss1++ = ',';
	AppendInt(ss1, GetSingleIBS(y, p)); *ss1++ = ',';
	AppendInt(ss1, crosssame); *ss1++ = '\0';

	for (int i = 0; i < crfan; ++i)
	{
		s2[i * 2] = (char)('0' + crfa[i]);
		s2[i * 2 + 1] = ',';
	}
	if (crfan) s2[crfan * 2 - 1] = 0;
	for (int i = 0; i < crfbn; ++i)
	{
		s3[i * 2] = (char)('0' + crfb[i]);
		s3[i * 2 + 1] = ',';
	}
	if (crfbn) s3[crfbn * 2 - 1] = 0;
	return HashString32(s1, s2, s3);
}

/* Get hash of a pair of genotype id */
TARGET HASH HashDyadGenotypeIndex(HASH ha)
{
#ifdef HASH64
	return (HASH)_mm_crc32_u64(0xB385FDA97FED5EA6ull, ha);
#else
	return (HASH)_mm_crc32_u64(0xB385FDA9, ha);
#endif
}

/* Get hash of a haplotype */
TARGET void HashHaplotype(IND* ti, int64 st, int64 ed, HASH* hash, int &ploidy)
{
	//allele should be sorted
#ifdef HASH64
	HASH Seed1[N_MAX_PLOIDY] = { 0xB385FDA97FED5EA6ull, 0xB385FDA97FED5EA6ull, 0xB385FDA97FED5EA6ull, 0xB385FDA97FED5EA6ull, 0xB385FDA97FED5EA6ull, 0xB385FDA97FED5EA6ull, 0xB385FDA97FED5EA6ull, 0xB385FDA97FED5EA6ull, 0xB385FDA97FED5EA6ull, 0xB385FDA97FED5EA6ull };
	uint *Seed2 = (uint*)Seed1;
#else
	HASH Seed1[N_MAX_PLOIDY] = { 0xB385FDA9, 0xB385FDA9, 0xB385FDA9, 0xB385FDA9, 0xB385FDA9, 0xB385FDA9, 0xB385FDA9, 0xB385FDA9, 0xB385FDA9, 0xB385FDA9 };
#endif
	uint64 alleles[N_MAX_PLOIDY] = { 0 };

	for (int64 l = st, hc = 0; l <= ed; ++l)
	{
		//load four alleles and hash
		GENOTYPE &gt = ti->GetGenotype(GetLocId(l), GetLoc(l).GetGtab());
		if (hc == 0) ploidy = gt.Ploidy();

		if (gt.Nalleles() == 0)
		{
			SetFF(hash, ploidy);
			return;
		}

		ushort *als = gt.GetAlleleArray();
		for (int vi = 0; vi < ploidy; ++vi)
			alleles[vi] = alleles[vi] << 16 | als[vi];

		if (++hc % 4 == 0)
		{
			for (int vi = 0; vi < ploidy; ++vi)
			{
#ifdef HASH64
				Seed2[(vi << 1)]     = (uint)_mm_crc32_u64(Seed2[(vi << 1)    ], alleles[vi]);
				Seed2[(vi << 1) + 1] = (uint)_mm_crc32_u64(Seed2[(vi << 1 + 1)], alleles[vi]);
#else
				Seed1[vi] = (HASH)_mm_crc32_u64(Seed1[vi], alleles[vi]);
#endif
			}
			SetZero(Seed1, N_MAX_PLOIDY);
		}
	}

	for (int vi = 0; vi < ploidy; ++vi)
	{
#ifdef HASH64
		Seed2[(vi << 1)    ] = (uint)_mm_crc32_u64(Seed2[(vi << 1)    ], alleles[vi]);
		Seed2[(vi << 1) + 1] = (uint)_mm_crc32_u64(Seed2[(vi << 1 + 1)], alleles[vi]);
#else
		Seed1[vi] = (HASH)_mm_crc32_u64(Seed1[vi], alleles[vi]);
#endif
	}

	SetVal(hash, Seed1, ploidy);
}