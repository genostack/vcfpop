/* Misc */

#include "vcfpop.h"
#pragma pack(push, 1)

/* Allocate genotype table in virtual memory */
TARGET void VAllocGenotype(int64 newlen)
{
	if (newlen > 0x280000000) //10GiB
		BIG_FILE = true;
	else
		BIG_FILE = false;

	//Len at 4kb unit
	VIRTUAL_MEMORY = true;
	static byte *addr = (byte*)V_BASE_GENOTYPE;

	//first alloc
	if (addr == (byte*)V_BASE_GENOTYPE)
	{
		V_ALLOC_CSIZE_GENOTYPE = 0;
		V_ALLOC_SIZE_GENOTYPE = 16;
		V_ALLOC_LEN_GENOTYPE = new int64[V_ALLOC_SIZE_GENOTYPE];
	}

	newlen -= addr - (byte*)V_BASE_GENOTYPE;
	while (newlen)
	{
		int64 clen = newlen > 0x20000000 ? 0x20000000 : newlen;
#ifdef _WIN64
		if (!VirtualAlloc(addr, clen, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE))
#else		
		if (mmap(addr, clen, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON, -1, 0) == (void*)-1)
#endif	
			Exit("\nError: cannot alloc %lld Mib memory.", ((int64)addr - V_BASE_GENOTYPE) / 1024 / 1024);
		
		if (V_ALLOC_CSIZE_GENOTYPE == V_ALLOC_SIZE_GENOTYPE)
		{
			V_ALLOC_SIZE_GENOTYPE <<= 1;
			int64 *tlen = new int64[V_ALLOC_SIZE_GENOTYPE];
			memcpy(tlen, V_ALLOC_LEN_GENOTYPE, V_ALLOC_CSIZE_GENOTYPE * sizeof(int64));
			delete[] V_ALLOC_LEN_GENOTYPE; V_ALLOC_LEN_GENOTYPE = tlen;
		}

		V_ALLOC_LEN_GENOTYPE[V_ALLOC_CSIZE_GENOTYPE++] = clen;
		addr += clen;
		newlen -= clen;
	}

}

/* Unallocate genotype table in virtual memory */
TARGET void VUnAllocGenotype()
{
	byte *addr = (byte*)V_BASE_GENOTYPE;
	for (int i = 0; i < V_ALLOC_CSIZE_GENOTYPE; ++i)
	{
#ifdef _WIN64
		if (!VirtualFree(addr, 0, MEM_RELEASE))
#else
		if (munmap(addr, V_ALLOC_LEN_GENOTYPE[i]))
#endif
			Exit("\nError: cannot unalloc %lld Mib memory at %llx.", V_ALLOC_LEN_GENOTYPE[i] / 1024 / 1024, addr);
		addr += V_ALLOC_LEN_GENOTYPE[i];
	}
	delete[] V_ALLOC_LEN_GENOTYPE;
	V_ALLOC_CSIZE_GENOTYPE = 0;
	V_ALLOC_SIZE_GENOTYPE = 0;
	VIRTUAL_MEMORY = false;
}

/* Allocate allele depth table in virtual memory */
TARGET void VAllocAlleleDepth(int64 newlen)
{
	if (!ploidyinfer) return;
	//Len at 4kb unit
	VIRTUAL_MEMORY = true;
	static byte *addr = (byte*)V_BASE_ALLELEDEPTH;
	//Alloc At 0x8000000000
	if (addr == (byte*)V_BASE_ALLELEDEPTH)
	{
		V_ALLOC_CSIZE_ALLELEDEPTH = 0;
		V_ALLOC_SIZE_ALLELEDEPTH = 16;
		V_ALLOC_LEN_ALLELEDEPTH = new int64[V_ALLOC_SIZE_ALLELEDEPTH];
	}

	newlen -= addr - (byte*)V_BASE_ALLELEDEPTH;
	while (newlen)
	{
		int64 clen = newlen > 0x20000000 ? 0x20000000 : newlen;
#ifdef _WIN64
		if (!VirtualAlloc(addr, clen, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE))
#else		
		if (mmap(addr, clen, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON, -1, 0) == (void*)-1)
#endif	
			Exit("\nError: cannot alloc %lld Mib memory.", ((int64)addr - V_BASE_ALLELEDEPTH) / 1024 / 1024);
		if (V_ALLOC_CSIZE_ALLELEDEPTH == V_ALLOC_SIZE_ALLELEDEPTH)
		{
			V_ALLOC_SIZE_ALLELEDEPTH <<= 1;
			int64 *tlen = new int64[V_ALLOC_SIZE_ALLELEDEPTH];
			memcpy(tlen, V_ALLOC_LEN_ALLELEDEPTH, V_ALLOC_CSIZE_ALLELEDEPTH * sizeof(int64));
			delete[] V_ALLOC_LEN_ALLELEDEPTH; V_ALLOC_LEN_ALLELEDEPTH = tlen;
		}
		V_ALLOC_LEN_ALLELEDEPTH[V_ALLOC_CSIZE_ALLELEDEPTH++] = clen;
		addr += clen;
		newlen -= clen;
	}
}

/* Unallocate allele depth table in virtual memory */
TARGET void VUnAllocAlleleDepth()
{
	if (!ploidyinfer) return;
	byte *addr = (byte*)V_BASE_ALLELEDEPTH;
	for (int i = 0; i < V_ALLOC_CSIZE_ALLELEDEPTH; ++i)
	{
#ifdef _WIN64
		if (!VirtualFree(addr, 0, MEM_RELEASE))
#else
		if (munmap(addr, V_ALLOC_LEN_ALLELEDEPTH[i]))
#endif
			Exit("\nError: cannot unalloc %lld Mib memory at %llx.", V_ALLOC_LEN_ALLELEDEPTH[i] / 1024 / 1024, addr);
		addr += V_ALLOC_LEN_ALLELEDEPTH[i];
	}
	delete[] V_ALLOC_LEN_ALLELEDEPTH;
	V_ALLOC_CSIZE_ALLELEDEPTH = 0;
	V_ALLOC_SIZE_ALLELEDEPTH = 0;
	VIRTUAL_MEMORY = false;
}

/* Count number of alleles from frequency array*/
TARGET int CountK(double *fre, int k2)
{
	int re = 0;
	for (int i = 0; i < k2; ++i)
		if (fre[i] > MIN_FREQ) re++;
	return re;
}

/* Set a bit in a variable */
TARGET void SetBit(byte &b, int pos, bool val)
{
	if (val)
		b |= (1 << pos);
	else
		b &= ~(1 << pos);
}

/* Get a bit in a variable */
TARGET bool GetBit(byte b, int pos)
{
	return b & (1 << pos);
}

/* Get number of different elements */
TARGET int GetNalleles(ushort *alleles, int ploidy)
{
	int nalleles = 0;
	for (int i = 0; i < ploidy; ++i)
	{
		bool isnew = true;
		ushort asi = alleles[i];
		for (int j = i - 1; j >= 0; --j)
		{
			if (asi == alleles[j])
			{
				isnew = false;
				break;
			}
		}
		if (isnew) nalleles++;
	}
	return nalleles;
}

/* Get current directory */
TARGET void GetCurDir(char *path)
{
	int dirlen = PATH_LEN - 1;
#ifdef _WIN64
	GetCurrentDirectoryA(dirlen, path);
	dirlen = (int)strlen(path);
	if (path[dirlen - 1] != '\\')
	{
		path[dirlen] = '\\';
		path[dirlen + 1] = 0;
	}
	char *t = ReplaceStr(path, "\\\\", "\\");
	strcpy(path, t);
	delete[] t;
#else
	getcwd(path, dirlen);
	strcat(path, "/");
	char *t = ReplaceStr(path, "\\", "/");
	strcpy(path, t);
	delete[] t;
	t = ReplaceStr(path, "//", "/");
	strcpy(path, t);
	delete[] t;
#endif
}

/* Clear temp files */
TARGET void ClearTempFiles(char *path)
{
	char filename[FILE_NAME_LEN];
#ifdef _WIN64
	_finddata_t file;
	char pattern[FILE_NAME_LEN];
	sprintf(pattern, "%s%s", path, "vcfpop_*.tmp");
	int64 handle = _findfirst(pattern, &file);
	if (handle != -1) do
	{
		sprintf(filename, "%s%s", path, file.name);
		remove(filename);
	} while (_findnext(handle, &file) != -1);
	_findclose(handle);
#else
	DIR* dir = opendir(path);
	if (!dir)
		Exit("\nError: cannot open temp dir %s. \n", path);
	while (true)
	{
		dirent* p = readdir(dir);
		if (!p) break;
		if (strstr(p->d_name, "vcfpop_*.tmp"))
			remove(filename);
	}
	closedir(dir);
#endif
}

/* Pause console */
TARGET void Pause(void)
{
#ifdef _WIN64
	system("pause");
#else
	printf("Press ENTER to continue...\n");
	getchar();
#endif
}

/* Exit program with a message */
TARGET void Exit(const char *fmt, ...)
{
	static bool first = true;
	if (first)
	{
		first = false;
		va_list ap;
		int n = 0;
		va_start(ap, fmt);
		n = vprintf(fmt, ap);
		va_end(ap);
		Pause();
		exit(0);
	}
	else
		SLEEP(200000);
}

#ifndef _CPOINT
TARGET CPOINT::CPOINT(int de)
{
	diff = 0;
	dim = de;
	SetVal(image, 0.0, 19);
}

TARGET bool CPOINT::operator >(CPOINT &a)
{
	return li > a.li;
}

TARGET bool CPOINT::operator >=(CPOINT &a)
{
	return li >= a.li;
}

TARGET bool CPOINT::operator <(CPOINT &a)
{
	return li < a.li;
}

TARGET bool CPOINT::operator <=(CPOINT &a)
{
	return li <= a.li;
}

TARGET bool CPOINT::operator ==(CPOINT &a)
{
	return li == a.li;
}

TARGET bool CPOINT::operator !=(CPOINT &a)
{
	return li != a.li;
}

TARGET CPOINT CPOINT::operator =(const CPOINT &a)
{
	SetVal(this, (CPOINT*)&a, 1);
	return *this;
}

TARGET CPOINT CPOINT::operator +(const CPOINT &b)
{
	CPOINT re(*this);
	for (int i = 0; i < dim; ++i)
		re.image[i] += b.image[i];
	return re;
}

TARGET CPOINT CPOINT::operator -(const CPOINT &b)
{
	CPOINT re(*this);
	for (int i = 0; i < dim; ++i)
		re.image[i] -= b.image[i];
	return re;
}

TARGET CPOINT CPOINT::operator *(const double b)
{
	CPOINT re(*this);
	for (int i = 0; i < dim; ++i)
		re.image[i] *= b;
	return re;
}

TARGET CPOINT CPOINT::operator /(const double b)
{
	CPOINT re(*this);
	for (int i = 0; i < dim; ++i)
		re.image[i] /= b;
	return re;
}

TARGET CPOINT &CPOINT::operator +=(CPOINT &a)
{
	for (int i = 0; i < dim; ++i)
		image[i] += a.image[i];
	return *this;
}

TARGET CPOINT &CPOINT::operator -=(CPOINT &a)
{
	for (int i = 0; i < dim; ++i)
		image[i] -= a.image[i];
	return *this;
}

TARGET CPOINT &CPOINT::operator *=(double a)
{
	for (int i = 0; i < dim; ++i)
		image[i] *= a;
	return *this;
}

TARGET CPOINT &CPOINT::operator /=(double a)
{
	for (int i = 0; i < dim; ++i)
		image[i] /= a;
	return *this;
}

/* Distance in real space */
TARGET double CPOINT::DistanceReal(CPOINT &a)
{
	double s = 0;
	for (int i = 0; i < dim; ++i)
		s += (real[i] - a.real[i]) * (real[i] - a.real[i]);
	return sqrt(s);
}

/* Distance in image space */
TARGET double CPOINT::DistanceImage(CPOINT &a)
{
	double s = 0;
	for (int i = 0; i < dim; ++i)
		s += (image[i] - a.image[i]) * (image[i] - a.image[i]);
	return sqrt(s);
}

/* Calculate real points */
TARGET void CPOINT::Image2Real()
{
	if (confine && !diff && dim % 2 == 0)
	{
		if (dim == 1)
		{
			real[0] = 1 / (1 + exp(-image[0]));
			real[1] = 1 - real[0];
		}
		else if (dim == 2)
		{
			double p1 = 1 / (1 + exp(-image[0]));
			double q1 = 1 / (1 + exp(-image[1]));
			double p0 = 1 - p1;
			double q0 = 1 - q1;
			real[0] = p1 * q1;
			real[1] = p0 * q1 + p1 * q0;
			real[2] = p0 * q0;
		}
		else if (dim == 4)
		{
			double p2 = 1 / (1 + exp(-image[0]));
			double p1 = (1 - p2) / (1 + exp(-image[1]));
			double q2 = 1 / (1 + exp(-image[2]));
			double q1 = (1 - q2) / (1 + exp(-image[3]));
			double p0 = 1 - p2 - p1;
			double q0 = 1 - q2 - q1;
			real[0] = p2 * q2;
			real[1] = p2 * q1 + p1 * q2;
			real[2] = p2 * q0 + p0 * q2 + p1 * q1;
			real[3] = p0 * q1 + p1 * q0;
			real[4] = p0 * q0;
		}
		else if (dim == 6)
		{
			double p3 = 1 / (1 + exp(-image[0]));
			double p2 = (1 - p3) / (1 + exp(-image[1]));
			double p1 = (1 - p3 - p2) / (1 + exp(-image[2]));
			double q3 = 1 / (1 + exp(-image[3]));
			double q2 = (1 - q3) / (1 + exp(-image[4]));
			double q1 = (1 - q3 - q2) / (1 + exp(-image[5]));
			double p0 = 1 - p3 - p2 - p1;
			double q0 = 1 - q3 - q2 - q1;
			real[0] = p3 * q3;
			real[1] = p3 * q2 + p2 * q3;
			real[2] = p3 * q1 + p2 * q2 * p1 * q3;
			real[3] = p3 * q0 + p2 * q1 + p1 * q2 + p0 * q3;
			real[4] = p2 * q0 + p1 * q1 + p0 * q2;
			real[5] = p1 * q0 + p0 * q1;
			real[6] = p0 * q0;
		}
		else if (dim == 8)
		{
			double p4 = 1 / (1 + exp(-image[0]));
			double p3 = (1 - p4) / (1 + exp(-image[1]));
			double p2 = (1 - p4 - p3) / (1 + exp(-image[2]));
			double p1 = (1 - p4 - p3 - p2) / (1 + exp(-image[3]));
			double q4 = 1 / (1 + exp(-image[4]));
			double q3 = (1 - q4) / (1 + exp(-image[5]));
			double q2 = (1 - q4 - q3) / (1 + exp(-image[6]));
			double q1 = (1 - q4 - q3 - q2) / (1 + exp(-image[7]));
			double p0 = 1 - p4 - p3 - p2 - p1;
			double q0 = 1 - q4 - q3 - q2 - q1;
			real[0] = p4 * q4;
			real[1] = p4 * q3 + p3 * q4;
			real[2] = p4 * q2 + p3 * q3 * p2 * q4;
			real[3] = p4 * q1 + p3 * q2 + p2 * q3 + p1 * q4;
			real[4] = p4 * q0 + p3 * q1 + p2 * q2 + p1 * q3 + p0 * q4;
			real[5] = p3 * q0 + p2 * q1 + p1 * q2 + p0 * q3;
			real[6] = p2 * q0 + p1 * q1 + p0 * q2;
			real[7] = p1 * q0 + p0 * q1;
			real[8] = p0 * q0;
		}
	}
	else
	{
		real[dim] = 1;
		for (int i = 0; i < dim; ++i)
		{
			real[i] = real[dim] / (1 + exp(-image[i]));
			real[dim] -= real[i];
		}
	}
}

/* Calculate real points */
TARGET void CPOINT::Image2RealSelfing()
{
	real[0] = 1.0 / (1 + exp(image[0]));
	if (real[0] <     MIN_FREQ) real[0] =     MIN_FREQ;
	if (real[0] > 1 - MIN_FREQ) real[0] = 1 - MIN_FREQ;
}

/* Break iteration if distance between points are smaller than eps */
TARGET bool /*static*/ CPOINT::IsBreak(CPOINT xx[9], double eps)
{
	return xx[0].DistanceReal(xx[xx[0].dim]) < eps;
}

/* Sort points */
TARGET void CPOINT::Order(CPOINT xx[9])
{
	int dim = xx[0].dim;
	for (int i = 0; i <= dim; ++i)
		for (int j = i + 1; j <= dim; ++j)
			if (xx[i] < xx[j])
				Swap(xx[i], xx[j]);
}

/* Down-Hill Simplex algorithm */
TARGET CPOINT CPOINT::DownHillSimplex(int dim, int diff, bool confine, double sep, int nrep, double (*Likelihood)(CPOINT&, void**), void **Param)
{
	CPOINT xx[20] = { 0 };

	for (int i = 0; i <= dim; ++i)
	{
		xx[i].dim = dim;
		xx[i].diff = diff;
		xx[i].confine = confine;
		if (i > 0) xx[i].image[i - 1] = 0.1;
		xx[i].li = Likelihood(xx[i], Param);
	}

	double likestop = LIKELIHOOD_TERM;
	for (int kk = 0; kk < nrep; ++kk)
	{
		//Order
		for (int searchcount = 0; ; ++searchcount)
		{
			CPOINT::Order(xx);
			if (searchcount >= MAX_ITER_DOWNHILL || 
				CPOINT::IsBreak(xx, LIKELIHOOD_TERM) && abs(xx[0].real[0] - xx[1].real[0]) < LIKELIHOOD_TERM * 1e-3)
				break;

			//Reflect
			CPOINT x0 = xx[0];
			for (int i = 1; i < dim; ++i)
				x0 += xx[i];
			x0 /= dim;
			
			CPOINT xr = x0 + (x0 - xx[dim]);
			xr.li = Likelihood(xr, Param); 

			//Expansion
			//best
			if (xr > xx[0])
			{
				CPOINT xe = x0 + (xr - x0) * 2; 
				xe.li = Likelihood(xe, Param); 
				SetVal(xx + 1, xx, dim);
				xx[0] = xe > xr ? xe : xr;
				continue;
			}

			//better than second worst
			if (xr > xx[dim - 1])
			{
				xx[dim] = xr;
				continue;
			}

			//worse than second worst
			//Contraction
			CPOINT xc = x0 + (xx[dim] - x0) * 0.5;
			xc.li = Likelihood(xc, Param); 
			if (xc > xx[1])
			{
				xx[dim] = xc;
				continue;
			}

			//Reduction
			for (int i = 1; i <= dim; ++i)
			{
				xx[i] = (xx[0] + xx[i]) * 0.5;
				xx[i].li = Likelihood(xx[i], Param);
			}
			continue;
		}

		//step 2
		for (int i = 1; i <= dim; ++i)
		{
			xx[i] = xx[0];
			xx[i].image[i - 1] *= (1 - sep);
			xx[i].li = Likelihood(xx[i], Param);
		}

		//Order
		sep /= 2;
		likestop /= 2;
	}

	CPOINT::Order(xx);
	return xx[0];
}
#endif

#ifndef _MEMORY
/* Initialize */
TARGET MEMORY::MEMORY()
{
	cblock = 0;
	block_size = 2;
	blocks = new MEMBLOCK[block_size];
	SetZero(blocks, block_size);
	blocks[0].size = 65536;
	blocks[0].bucket = new bool[blocks[0].size];
	InitLock(lock);
}

/* Uninitialize */
TARGET MEMORY::~MEMORY()
{
	if (cblock != 0xFFFFFFFF)
		for (int i = 0; i <= cblock; ++i)
			delete[] blocks[i].bucket;
	if (blocks) delete[] blocks;
	blocks = NULL;
	cblock = 0xFFFFFFFF;
}

// Clear entries
TARGET void MEMORY::ClearMemory()
{
	for (int i = 0; i <= cblock; ++i)
		blocks[i].used = 0;
	cblock = 0;
}

// Expand number of blocks
TARGET void MEMORY::Expand()
{
	MEMBLOCK* nblock = new MEMBLOCK[(uint64)block_size << 1];
	SetVal(nblock, blocks, cblock);
	SetZero(nblock + cblock, (block_size << 1) - cblock);
	delete[] blocks; blocks = nblock;
	block_size <<= 1;
}

// Allocate a small piece of memory
TARGET byte* MEMORY::Alloc(int size)
{
	//if (islock) 
	Lock(lock);
	while (size + blocks[cblock].used > blocks[cblock].size)
	{
		if (++cblock == block_size) Expand();
		if (!blocks[cblock].bucket)
		{
			blocks[cblock].size = Max(Min(blocks[cblock - 1].size << 1, BIG_FILE ? 256 * 1024 * 1024 : 8 * 1024 * 1024), size);
			blocks[cblock].bucket = new bool[blocks[cblock].size];
		}
	}
	byte* addr = (byte*)(blocks[cblock].bucket + blocks[cblock].used);
	blocks[cblock].used += size;
	//if (islock) 
	UnLock(lock);
	return addr;
}
#endif

/* Initialize double reduction rates */
TARGET void InitAlpha()
{
	memset(ALPHA, 0, sizeof(ALPHA));

	//rcs
	for (int i = 1; i <= N_DRE_MODELT; ++i)
		for (int j = 1; j <= N_MAX_PLOIDY; ++j)
			ALPHA[i][j][0] = 1;

	//prcs
	ALPHA[2][4][0] = 6.0 / 7;		ALPHA[2][4][1] = 1.0 / 7;
	ALPHA[2][6][0] = 8.0 / 11;		ALPHA[2][6][1] = 3.0 / 11;
	ALPHA[2][8][0] = 40.0 / 65;		ALPHA[2][8][1] = 24.0 / 65;		ALPHA[2][8][2] = 1.0 / 65;
	ALPHA[2][10][0] = 168.0 / 323;	ALPHA[2][10][1] = 140.0 / 323;	ALPHA[2][10][2] = 15.0 / 323;

	//ces
	ALPHA[3][4][0] = 5.0 / 6;		ALPHA[3][4][1] = 1.0 / 6;
	ALPHA[3][6][0] = 7.0 / 10;		ALPHA[3][6][1] = 3.0 / 10;
	ALPHA[3][8][0] = 83.0 / 140;	ALPHA[3][8][1] = 54.0 / 140;	ALPHA[3][8][2] = 3.0 / 140;
	ALPHA[3][10][0] = 127.0 / 252;	ALPHA[3][10][1] = 110.0 / 252;	ALPHA[3][10][2] = 15.0 / 252;

	for (int E = 0; E <= 100; ++E)
	{
		double e = (double)E / 100;
		for (int p = 4; p <= N_MAX_PLOIDY; p += 2)
		{
			int D = p / 4;
			double apie[4] = { 0 };

			for (int i = 0; i <= D; ++i)
				apie[i] = BINOMIAL[p / 2][i] * BINOMIAL[p / 2 - i][p / 2 - 2 * i] * pow(2.0, p / 2 - 2 * i) / BINOMIAL[p][p / 2];

			for (int dd = 0; dd <= D; ++dd)
			{
				ALPHA[4 + E][p][dd] = 0;
				for (int j = dd; j <= D; ++j)
					for (int i = j; i <= D; ++i)
						ALPHA[4 + E][p][dd] += apie[i] * BINOMIAL[i][j] * pow(2.0, (double)-i) * BINOMIAL[j][dd] * pow(e, dd) * pow((double)1 - e, (double)j - dd);
			}
		}
	}
}

/* Calculate task in multiple threads */
TARGET void RunThreads(void (*Func) (int), void (*GuardFunc1) (int), void (*GuardFunc2) (int),
	int64 ntot, int64 nct, const char *info, int nthreads, bool isfirst, int nprogress)
{
	if (isfirst)
	{
		printf("%s", info);
		fflush(stdout);
		PROGRESS_VALUE = PROGRESS_NOUTPUTED = PROGRESS_NOUTPUTED2 = 0;
		PROGRESS_TOTAL = ntot;
	}
	progress1 = progress2 = 0;
	SetZero(state_lock, NBUF);
	for (int i = 0; i < NBUF; ++i)
		state_lock[i] = (int64)i << 2;

	PROGRESS_CEND = PROGRESS_VALUE + nct;
	PROGRESS_CSTART = PROGRESS_VALUE;
	VLA_NEW(hThread, thread, nthreads + 2);
	VLA_NEW(fThread, bool, nthreads + 2);

	for (int i = 0; i < nthreads + 2; ++i)
	{
		fThread[i] = true;
		if (i == nthreads && GuardFunc1)
			hThread[i] = thread(GuardFunc1, i);
		else if (i == nthreads + 1 && GuardFunc2)
			hThread[i] = thread(GuardFunc2, i);
		else if (i < nthreads && Func)
			hThread[i] = thread(Func, i);
		else
			fThread[i] = false;
	}

	while (PROGRESS_VALUE != PROGRESS_CEND)
	{
		SLEEP(SLEEP_TIME_TINY);
		while (PROGRESS_VALUE * (double)nprogress / PROGRESS_TOTAL > PROGRESS_NOUTPUTED)
		{
			if ((PROGRESS_NOUTPUTED + PROGRESS_NOUTPUTED2) % 10 == 0)
				printf("%lld", ((PROGRESS_NOUTPUTED + PROGRESS_NOUTPUTED2) / 10) % 10);
			else if ((PROGRESS_NOUTPUTED + PROGRESS_NOUTPUTED2) % 5 == 0)
				printf("o");
			else
				printf(".");
			PROGRESS_NOUTPUTED++;
			fflush(stdout);
		}
	}

	if (PROGRESS_CEND == PROGRESS_TOTAL)
	{
		while (PROGRESS_NOUTPUTED < nprogress)
		{
			if ((PROGRESS_NOUTPUTED + PROGRESS_NOUTPUTED2) % 10 == 0)
				printf("%lld", ((PROGRESS_NOUTPUTED + PROGRESS_NOUTPUTED2) / 10) % 10);
			else if ((PROGRESS_NOUTPUTED + PROGRESS_NOUTPUTED2) % 5 == 0)
				printf("o");
			else
				printf(".");
			PROGRESS_NOUTPUTED++;
			fflush(stdout);
		}
	}

	for (int i = 0; i < nthreads + 2; ++i)
		if (fThread[i])
			hThread[i].join();

	VLA_DELETE(hThread);
	VLA_DELETE(fThread);
}

#pragma pack(pop)