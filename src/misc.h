/* Misc */

#pragma once
#include "vcfpop.h"
#pragma pack(push, 1)

/* Initialize a lock */
TARGET inline void InitLock(LOCK& x)
{
#ifdef LOCKFREE

#else 
	#ifdef _WIN64
		InitializeCriticalSection(&x);
	#else
		x = PTHREAD_MUTEX_INITIALIZER;
	#endif
#endif
}

/* Lock a lock */
TARGET inline void Lock(LOCK& x)
{
#ifdef LOCKFREE
	while (x.test_and_set());
#else 
	#ifdef _WIN64
		EnterCriticalSection(&x);
	#else
		pthread_mutex_lock(&x);
	#endif
#endif
}

/* UnLock a lock */
TARGET inline void UnLock(LOCK& x)
{
#ifdef LOCKFREE
	x.clear();
#else 
	#ifdef _WIN64
		LeaveCriticalSection(&x);
	#else
		pthread_mutex_unlock(&x);
	#endif
#endif
}

/* Allocate genotype table in virtual memory */
TARGET void VAllocGenotype(int64 extralen);

/* Unallocate genotype table in virtual memory */
TARGET void VUnAllocGenotype();

/* Allocate allele depth table in virtual memory */
TARGET void VAllocAlleleDepth(int64 extralen);

/* Unallocate allele depth table in virtual memory */
TARGET void VUnAllocAlleleDepth();

/* Count number of alleles from frequency array*/
TARGET int CountK(double *fre, int k2);

/* Set a bit in a variable */
TARGET void SetBit(byte &b, int pos, bool val);

/* Get a bit in a variable */
TARGET bool GetBit(byte b, int pos);

/* Get number of unique elements */
TARGET int GetNalleles(ushort *alleles, int ploidy);

/* Get current directory */
TARGET void GetCurDir(char *path);

/* Try to delete an array */
template<typename T>
TARGET void TryDelete(T *&pointer)
{
	if (pointer) delete[] pointer;
	pointer = NULL;
}

/* Clear temp files */
TARGET void ClearTempFiles(char *path);

/* Pause console */
TARGET void Pause(void);

/* Exit program with a message */
TARGET void Exit(const char *fmt, ...);

/* Atomic add val to ref */
template<typename T>
TARGET inline void AtomicAdd4(volatile T &ref, T val)
{
#if defined(__clang__) || defined(__GNUC__)
	__sync_fetch_and_add(&ref, (int)val);
#else
	InterlockedAdd((volatile int*)&ref, (int)val);
#endif
}

/* Atomic add val to ref */
template<typename T>
TARGET inline void AtomicAdd4(volatile T *ref, T *val, int len)
{
	for (int i = 0; i < len; ++i)
#if defined(__clang__) || defined(__GNUC__)
		__sync_fetch_and_add(&ref[i], (int)val[i]);
#else
		InterlockedAdd((volatile int*)&ref[i], (int)val[i]);
#endif
}

/* Atomic add val to ref */
template<typename T>
TARGET inline void AtomicAdd8(volatile T &ref, T val)
{
#if defined(__clang__) || defined(__GNUC__)
	__sync_fetch_and_add(&ref, (int64)val);
#else
	InterlockedAdd64((volatile int64*)&ref, (int64)val);
#endif
}

/* Atomic add val to ref */
template<typename T>
TARGET inline void AtomicAdd8(volatile T *ref, T *val, int len)
{
	for (int i = 0; i < len; ++i)
#if defined(__clang__) || defined(__GNUC__)
		__sync_fetch_and_add(&ref[i], (int64)val[i]);
#else
		InterlockedAdd64((volatile int64*)&ref[i], (int64)val[i]);
#endif
}

/* Atomic set min */
TARGET inline void AtomicMin1(volatile byte &ref, byte val)
{
#if defined(__clang__) || defined(__GNUC__)
	for (byte ov = ref, nv = Min(ov, val);
		!__sync_bool_compare_and_swap((char*)&ref, (char)ov, (char)nv);
		ov = ref, nv = Max(ov, val));
#else
	for (uint ov = *(volatile uint*)&ref,
		nv = Min(ov & 0xFF, (uint)val) | (ov & 0xFFFFFF00);
		ov != InterlockedCompareExchange((uint*)&ref, nv, ov);
		ov = *(uint*)&ref,
		nv = Max(ov & 0xFF, (uint)val) | (ov & 0xFFFFFF00));
#endif
}

/* Atomic set max */
TARGET inline void AtomicMax1(volatile byte &ref, byte val)
{
#if defined(__clang__) || defined(__GNUC__)
	for (byte ov = ref, nv = Max(ov, val);
		!__sync_bool_compare_and_swap((char*)&ref, (char)ov, (char)nv);
		      ov = ref, nv = Max(ov, val));
#else
	for (uint ov = *(volatile uint*)&ref,
		nv = Max(ov & 0xFF, (uint)val) | (ov & 0xFFFFFF00);
		ov != InterlockedCompareExchange((uint*)&ref, nv, ov);
		ov = *(uint*)&ref,
		nv = Max(ov & 0xFF, (uint)val) | (ov & 0xFFFFFF00));
#endif
}

/* Atomic add val to ref */
TARGET inline void AtomicAddD(volatile double &ref, double val)
{
	for (double ov = ref, nv = ov + val;
#if defined(__clang__) || defined(__GNUC__)
		!__sync_bool_compare_and_swap((uint64*)&ref, *(uint64*)&ov, *(uint64*)&nv);
#else
		*(uint64*) & ov != InterlockedCompareExchange((uint64*)&ref, *(uint64*)&nv, *(uint64*)&ov);
#endif
		        ov = ref, nv = ov + val);
}

/* Atomic add val to ref */
TARGET inline void AtomicAddD(volatile double *ref, double *val, int len)
{
	for (int i = 0; i < len; ++i)
		for (double ov = ref[i], nv = ov + val[i];
#if defined(__clang__) || defined(__GNUC__)
			!__sync_bool_compare_and_swap((uint64*)&ref[i], *(uint64*)&ov, *(uint64*)&nv);
#else
			*(uint64*)&ov != InterlockedCompareExchange((uint64*)&ref[i], *(uint64*)&nv, *(uint64*)&ov);
#endif
			        ov = ref[i], nv = ov + val[i]);
}

/* Point in a simplex to optimize by Down-Hill Simplex algorithm */
class CPOINT
{
public:
	double image[8];				//Coordinates in converted image space
	double real[9];					//Coordinates in real space
	double li;						//Logarithm of likelihood
	int dim;						//Dimensions
	int diff;						//Different ploidy levels
	bool confine;					//Subject to an additioanl constraint in Anderson 2007 relatedness estimator

	TARGET CPOINT(int de = 4);

	TARGET bool operator >(CPOINT &a);

	TARGET bool operator >=(CPOINT &a);

	TARGET bool operator <(CPOINT &a);

	TARGET bool operator <=(CPOINT &a);

	TARGET bool operator ==(CPOINT &a);

	TARGET bool operator !=(CPOINT &a);

	TARGET CPOINT operator =(const CPOINT &a);

	TARGET CPOINT operator +(const CPOINT &b);

	TARGET CPOINT operator -(const CPOINT &b);

	TARGET CPOINT operator *(const double b);

	TARGET CPOINT operator /(const double b);

	TARGET CPOINT &operator +=(CPOINT &a);

	TARGET CPOINT &operator -=(CPOINT &a);

	TARGET CPOINT &operator *=(double a);

	TARGET CPOINT &operator /=(double a);

	/* Distance in real space */
	TARGET double DistanceReal(CPOINT &a);

	/* Distance in image space */
	TARGET double DistanceImage(CPOINT &a);

	/* Calculate real points */
	TARGET void Image2Real();

	/* Calculate real points */
	TARGET void Image2RealSelfing();

	/* Break iteration if distance between points are smaller than eps */
	TARGET bool static IsBreak(CPOINT xx[9], double eps = 1e-7);

	/* Sort points */
	TARGET void static Order(CPOINT xx[9]);

	/* Down-Hill Simplex algorithm */
	TARGET CPOINT static DownHillSimplex(int dim, int diff, bool confine, double sep, int nrep, double (*Likelihood)(CPOINT&, void**), void **Param);
};

/* Block of local memory management class */
struct MEMBLOCK
{
	bool* bucket;							//Data
	int size;								//Bucket size in bytes
	int used;								//Number of used bytes
};

/* Local memory management class for locus and genotype */
class MEMORY
{
public:
	MEMBLOCK* blocks;						//Blocks array
	int block_size;							//Blocks size
	int cblock;								//Current block index
	LOCK lock; 								//Lock during expansion and allocation

	/* Initialize */
	TARGET MEMORY();

	/* Uninitialize */
	TARGET ~MEMORY();

	// Clear entries
	TARGET void ClearMemory();

	// Expand number of blocks
	TARGET void Expand();

	// Allocate a small piece of memory
	template <typename T>
	TARGET T* Alloc(T*& addr, int size, bool islock = true)
	{
		if (islock) Lock(lock);
		size *= sizeof(T);
		while (size + blocks[cblock].used > blocks[cblock].size)
		{
			if (++cblock == block_size) Expand();
			if (!blocks[cblock].bucket)
			{
				blocks[cblock].size = Max(Min(blocks[cblock - 1].size << 1, BIG_FILE ? 256 * 1024 * 1024 : 8 * 1024 * 1024), size);
				blocks[cblock].bucket = new bool[blocks[cblock].size];
			}
		}
		addr = (T*)(blocks[cblock].bucket + blocks[cblock].used);
		blocks[cblock].used += size;
		if (islock) UnLock(lock);
		return addr;
	}

	// Allocate a small piece of memory
	TARGET byte* Alloc(int size);

	// Free a small piece of memory
	template <typename T>
	TARGET void Free(T* addr, int size, bool islock = true)
	{
		if (islock) Lock(lock);
		if (addr + size == (T*)(blocks[cblock].bucket + blocks[cblock].used))
			blocks[cblock].used -= size * sizeof(T);
		if (islock) UnLock(lock);
	}

	// ReAllocate a small piece of memory
	template <typename T>
	TARGET void ReAlloc(T*& addr, int size, T* oaddr, int osize)
	{
		Lock(lock);
		Free(oaddr, osize, false);
		Alloc(addr, size, false);
		if (oaddr != addr)
			SetVal(addr, oaddr, osize);
		UnLock(lock);
	}

	// Move from new mem to odd mem, free new mem, and increase old mem
	template <typename T>
	TARGET void Move(T*& oaddr, int osize, T* naddr, int nsize)
	{
		Lock(lock);

		// oaddr and naddr are adjacent and naddr is the newest
		if (oaddr + osize == naddr &&
			naddr + nsize == (T*)(blocks[cblock].bucket + blocks[cblock].used))
		{
			SetVal(oaddr, naddr, nsize); //move mem
			blocks[cblock].used -= osize * sizeof(T);
		}
		else
			oaddr = naddr;

		UnLock(lock);
	}
};

/* Fast list class */
template <typename T>
class LIST
{
public:
	T *bucket;							//Bucket
	int bucket_size;					//Size of bucket
	int size;							//Number of entries
	MEMORY *memory;						//Memory management class

	/* Deep copy a list */
	TARGET LIST<T> &operator=(LIST<T> &ref)
	{
		bucket_size = ref.bucket_size;
		size = ref.size;
		memory = ref.memory;

		if (memory)
			bucket = ref.bucket;
		else
			bucket = (T*)malloc(bucket_size * sizeof(T));

		SetVal(bucket, ref.bucket, bucket_size);
		return *this;
	}

	/* Set Zero */
	/*TARGET LIST()
	{
		SetZero(this, 1);
	}
	*/

	/* Initialize */
	TARGET LIST(MEMORY *_memory = NULL)
	{
		memory = _memory;
		bucket_size = 0x8;
		size = 0;
		bucket = memory ? (T*)memory->Alloc(bucket_size * sizeof(T)) : (T*)malloc(bucket_size * sizeof(T));
		SetFF(bucket, bucket_size);
	}

	/* Clear entries but do not shrink bucket */
	TARGET void Clear()
	{
		SetFF(bucket, bucket_size);
		size = 0;
	}

	/* Uninitialize and free memory */
	TARGET ~LIST()
	{
		if (memory && bucket)
			memory->Free(bucket, bucket_size);

		if (!memory && bucket) 
			free(bucket);

		bucket = NULL;
		size = 0;
		bucket_size = 0;
	}

	/* Expand list */
	TARGET void Expand()
	{
		int nsize = BIG_FILE ? bucket_size * 2 : Min(bucket_size * 2, bucket_size + 1024);
		if (memory)
			memory->ReAlloc(bucket, nsize, bucket, bucket_size);
		else
			bucket = (T*)realloc(bucket, nsize * sizeof(T));
		SetFF(bucket + bucket_size, nsize - bucket_size);
		bucket_size = nsize;
	}

	/* Access by index */
	TARGET T &operator[](uint64 i)
	{
		if (i >= size || i < 0)
			Exit("\nError: LIST class access an unexist element.\n");
		return bucket[i];
	}

	/* Add an entry to the end */
	TARGET void Push(T &val)
	{
		if (size >= bucket_size) Expand();
		bucket[size++] = val;
	}

	/* Return and delete the entry at the end */
	TARGET T Pop()
	{
		return bucket[--size];
	}

	/* Remove an entry by index */
	TARGET void Erase(int i)
	{
		if (i >= bucket_size || i < 0) Exit("\nError: FastVector class erase an unexist element.\n");
		SetVal(bucket + i, bucket + i + 1, size - i - 1);
		SetFF(bucket + size - 1, 1);
		size--;
	}

	/* Alloc memory */
	TARGET void SetSize(int _nsize)
	{
		int nsize = Max(bucket_size, _nsize + 1024);
		if (memory)
			memory->ReAlloc(bucket, nsize, bucket, bucket_size);
		else
			bucket = (T*)realloc(bucket, nsize * sizeof(T));
		SetFF(bucket + bucket_size, nsize - bucket_size);
		bucket_size = nsize;
	}
};

/* Hash table entry: key-val pair */
template <typename T, typename T2>
struct TABLE_ENTRY
{
	T key;								//Hash key
	T2 val;								//Table value
};

/* Fast hash table class */
template <typename T, typename T2>
class TABLE
{
public:
	TABLE_ENTRY<T,T2> *bucket;		//Key-val pair bucket
	uint *index;					//The list saves the position in bucket of ith entry
	uint mask;						//mask + 1 is the bucket size
	int size;						//Number of entries
	MEMORY *memory;					//Memory management class

	/* Set Zero */
	TARGET TABLE()
	{
		SetZero(this, 1);
	}

	/* Initialize a table */
	TARGET TABLE(bool haslist, MEMORY *_memory, int _size = 0)
	{
		memory = _memory;
		if (_size <= 4)
			mask = 3;
		else
			mask = (1 << CeilLog2(_size + (_size >> 1))) - 1;

		size = 0;
		byte *buf = memory ? 
			memory->Alloc((sizeof(TABLE_ENTRY<T, T2>) + (haslist ? sizeof(uint) : 0)) * (mask + 1)) : 
			(byte*)malloc((sizeof(TABLE_ENTRY<T, T2>) + (haslist ? sizeof(uint) : 0)) * (mask + 1));

		bucket = (TABLE_ENTRY<T, T2>*)buf;
		index = (uint*)(haslist ? bucket + (mask + 1) : NULL);
		Clear();
	}

	/* Clear entries but do not release memory */
	TARGET void Clear()
	{
		SetFF((byte*)bucket, (sizeof(TABLE_ENTRY<T, T2>) + (index ? sizeof(uint) : 0)) * (mask + 1));
		size = 0;
	}

	/* Free memory */
	TARGET ~TABLE()
	{
		if (!memory && bucket) 
			free(bucket);

		if (memory && bucket)
			memory->Free((byte*)bucket, (sizeof(TABLE_ENTRY<T, T2>) + (index ? sizeof(uint) : 0)) * (mask + 1));

		bucket = NULL;
		index = NULL;
		size = 0;
		mask = 0;
	}

	/* Expand table */
	TARGET void Expand()
	{
		uint mask2 = ((mask + 1) << 1) - 1;
		int olen = (sizeof(TABLE_ENTRY<T, T2>) + (index ? sizeof(uint) : 0)) * (mask + 1);
		int nlen = (sizeof(TABLE_ENTRY<T, T2>) + (index ? sizeof(uint) : 0)) * (mask2 + 1);

		byte*obuf = (byte*)bucket;
		byte *nbuf = memory ? memory->Alloc(nlen) : (byte*)malloc(nlen);

		SetFF(nbuf, nlen);
		TABLE_ENTRY<T, T2> *nbucket = (TABLE_ENTRY<T, T2>*)nbuf;
		uint *nindex = (uint*)(index ? nbucket + (mask2 + 1) : 0);

		//move to new bucket
		if (index)
		{
			//get an entry from old table
			for (int i = 0; i < size; ++i)
			{
				//find an empty emtry in the new table
				uint st = bucket[index[i]].key & mask2;
				for (uint j = 0; j <= mask2; ++j)
					if (nbucket[(j + st) & mask2].key == (T)-1)
					{
						//insert
						nbucket[nindex[i] = ((j + st) & mask2)] = bucket[index[i]];
						break;
					}
			}
		}
		else for (uint i = 0; i <= mask; ++i)
		{
			//get an entry from old table
			T key = bucket[i].key;
			if (key == (T)-1) continue;

			//find an empty slot in the new table
			uint st = key & mask2;
			for (uint j = 0; j <= mask2; ++j)
				if (nbucket[(j + st) & mask2].key == (T)-1)
				{
					//insert
					nbucket[(j + st) & mask2] = bucket[i];
					break;
				}
		}
		if (memory)
		{
			memory->Move(obuf, olen, nbuf, nlen);
			bucket = (TABLE_ENTRY<T, T2>*)obuf;
		}
		else
		{
			free(bucket);
			bucket = nbucket;
		}
		index = index ? (uint*)(bucket + (mask2 + 1)) : NULL;
		mask = mask2;
	}

	/* Get array [k + k * k] with k allele sizes and k * k squared distance between alleles for SMM model */
	TARGET void GetLength(T2 *&alen, MEMORY &locus_mem)
	{
		//For non-vcf input, allele identifier is the size
		//Return an array [size] + [size * size]
		//Table is the allele length table, key:value = allele:length

		locus_mem.Alloc(alen, size + size * size);
		T2* alen2 = alen + size;
		for (int i = 0; i < size; ++i)
		{
			alen[i] = bucket[index[i]].key;
			alen2[i * size + i] = 0;
			for (int j = 0; j < i; ++j)
				alen2[i * size + j] = alen2[j * size + i] = (alen[i] - alen[j]) * (alen[i] - alen[j]);
		}
	}

	/* Return value if entry exists, otherwise return null */
	TARGET T2 Try(T key)
	{
		uint st = key & mask;
		for (uint i = 0; i <= mask; ++i)
		{
			uint idx = (i + st) & mask;
			T k = bucket[idx].key;
			if (k == (T)-1) return (T2)NULL;
			if (k == key) return bucket[idx].val;
		}
		return (T2)NULL;
	}

	/* Use table as list where value is id, insert if entry doesn't exists, and return id */
	TARGET T2 PushIndex(T key)
	{
	restart:
		uint st = key & mask;
		for (uint i = 0; i <= mask; ++i)
		{
			uint idx = (i + st) & mask;
			T k = bucket[idx].key;
			if (k == (T)-1)
			{
				if (mask != 3 && (uint)(size + (size >> 1)) > mask)
				{
					Expand();
					goto restart;
				}
				bucket[idx].key = key;
				bucket[idx].val = (T2)size;
				if (index) index[size] = idx;
				return (T2)size++;
			}
			if (k == key) return bucket[idx].val;
		}

		if (mask == 3 && size == 4)
		{
			Expand();
			goto restart;
		}

		return (T2)-1;
	}

	/* Has an entry with the key */
	TARGET bool ContainsKey(T key)
	{
		uint st = key & mask;
		for (uint i = 0; i <= mask; ++i)
		{
			uint idx = (i + st) & mask;
			T k = bucket[idx].key;
			if (k == (T)-1) return false;
			if (k == key) return true;
		}
		return false;
	}

	/* Access entry by index */
	TARGET T2& operator()(uint64 i)
	{
		if (i >= size || i < 0)
			Exit("\nError: index exceed limit in interal LIST class.\n");
		return bucket[index[i]].val;
	}

	/* Access entry by key, single thread, not need to lock */
	TARGET T2& operator[](T key)
	{
	restart:
		//search from hash % (mask + 1)
		uint st = key & mask;

		for (uint i = 0; i <= mask; ++i)
		{
			uint idx = (i + st) & mask;
			T k = bucket[idx].key;
			if (k == (T)-1)
			{
				//do not exist this key
				//insert a new key
				//return ref of val
				if (mask != 3 && (uint)(size + (size >> 1)) > mask)
				{
					Expand();
					goto restart;
				}
				bucket[idx].key = key;
				if (index) index[size] = idx;
				size++;
				return bucket[idx].val;
			}
			if (k == key)
				return bucket[idx].val;
		}

		if (mask == 3 && size == 4)
		{
			Expand();
			goto restart;
		}

		Exit("\nError: key not found.\n");
		return bucket[0].val;
	}
};

/* Initialize double reduction rates */
TARGET void InitAlpha();

/* Calculate task in multiple threads */
TARGET void RunThreads(void (*Func) (int), void (*GuardFunc1) (int), void (*GuardFunc2) (int),
	int64 ntot, int64 nct, const char *info, int nthreads, bool isfirst, int nprogress = g_progress_val);
#pragma pack(pop)