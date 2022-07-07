#pragma once

#define WINDOWS_IGNORE_PACKING_MISMATCH
#define _CRT_SECURE_NO_WARNINGS

#include <immintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <new>
#include <algorithm>
#include <utility>
#include <float.h>
#include <math.h>
#include <time.h>
#include <atomic> 
#include <shared_mutex>
#include <mutex> 
#include <thread>
#include <cmath>
#include <ctime>
#include <cstdint>
#include <omp.h>

#include "Eigen/Eigen"
#include "Eigen/LU"
#include "Eigen/QR"
#include "Eigen/SVD"
#include "Spectra/SymEigsSolver.h"

using namespace Eigen;
using namespace Spectra;
using std::mutex;
using std::shared_mutex;
using std::map;
using std::thread;
using std::atomic;
using std::atomic_flag;

#define EIGEN_DONT_ALIGN_STATICALLY
#define Z_LARGE64
#define ZLIB_WINAPI
#define LOCKFREE_

typedef unsigned char byte;
typedef unsigned int uint;
typedef unsigned short ushort;
typedef unsigned long long uint64;
typedef long long int64;

#ifdef OVER4GLOC
	typedef int64 LOCN;
#else
	typedef uint LOCN;
#endif

#ifdef __APPLE__
	#include <mach-o/dyld.h>
	#define fseeko64 fseeko
	#define ftello64 ftello
	#define off64_t long long
#endif

#ifdef _WIN64
	#pragma warning(disable:26451)
	#pragma warning(disable:26495)
	#pragma warning(disable:4706)
	#pragma warning(disable:4996)
	#pragma warning(disable:4366)
	#pragma warning(disable:4310)
	#pragma warning(disable:4244)
	#pragma warning(disable:4503)

	#define TARGET
	#define TARGETMMX
	#define TARGETSSE
	#define TARGETAVX
	#define TARGET512

	#define fseeko64 _fseeki64
	#define ftello64 _ftelli64 

	//#include <windows.h>
	//#include <psapi.h>
	//#pragma comment (lib,"psapi.lib")
	#include "zlib/zlib.h"  
	#pragma comment (lib,"zlibstat.lib")
	//#include <process.h>
	//#include <direct.h>
	#include <io.h>

	#ifdef LOCKFREE
		typedef atomic_flag LOCK;
	#else
		typedef CRITICAL_SECTION LOCK;
	#endif

	#define SLEEP(x) Sleep(x)
#else
	#define VLA
	#define TARGET                              __attribute__((__target__("sse", "sse2", "sse3", "sse4", "popcnt", "lzcnt")))
	#define TARGETMMX __attribute__((noinline)) __attribute__((__target__("sse", "sse2", "sse3", "sse4", "popcnt", "lzcnt")))
	#define TARGETSSE __attribute__((noinline)) __attribute__((__target__("sse", "sse2", "sse3", "sse4", "popcnt", "lzcnt")))
	#define TARGETAVX __attribute__((noinline)) __attribute__((__target__("sse", "sse2", "sse3", "sse4", "popcnt", "lzcnt", "avx", "avx2", "fma")))
	#define TARGET512 __attribute__((noinline)) __attribute__((__target__("sse", "sse2", "sse3", "sse4", "popcnt", "lzcnt", "avx", "avx2", "fma", "avx512f", "avx512bw")))

	//#include <pthread.h>
	#include <sys/types.h>
	#include <sys/mman.h>
	#include <unistd.h>
	#include <stdarg.h>
	#include <dirent.h>
	#include <zlib.h>
	#pragma comment (lib,"libz.a")

	#ifdef LOCKFREE
		typedef atomic_flag LOCK;
	#else
		typedef pthread_mutex_t LOCK;
	#endif

	#define SLEEP(x) usleep((x)*1000)
#endif


#if defined(__clang__) || defined(__GNUC__)
	#define VLA
	#define m512d_f64(x,y) (x)[y]
	#define m256d_f64(x,y) (x)[y]
	#define m128d_f64(x,y) (x)[y]

	#define m512i_u64(x,y) ((uint64*)&(x))[y]
	#define m256i_u64(x,y) ((uint64*)&(x))[y]
	#define m128i_u64(x,y) ((uint64*)&(x))[y]

	#define m512i_u32(x,y) ((uint*)&(x))[y]
	#define m256i_u32(x,y) ((uint*)&(x))[y]
	#define m128i_u32(x,y) ((uint*)&(x))[y]

	#define m512i_u16(x,y) ((ushort*)&(x))[y]
	#define m256i_u16(x,y) ((ushort*)&(x))[y]
	#define m128i_u16(x,y) ((ushort*)&(x))[y]
#else
	#define m512d_f64(x,y) (x).m512d_f64[y]
	#define m256d_f64(x,y) (x).m256d_f64[y]
	#define m128d_f64(x,y) (x).m128d_f64[y]

	#define m512i_u64(x,y) (x).m512i_u64[y]
	#define m256i_u64(x,y) (x).m256i_u64[y]
	#define m128i_u64(x,y) (x).m128i_u64[y]

	#define m512i_u32(x,y) (x).m512i_u32[y]
	#define m256i_u32(x,y) (x).m256i_u32[y]
	#define m128i_u32(x,y) (x).m128i_u32[y]

	#define m512i_u16(x,y) (x).m512i_u16[y]
	#define m256i_u16(x,y) (x).m256i_u16[y]
	#define m128i_u16(x,y) (x).m128i_u16[y]
#endif

#ifdef VLA
	#define VLA_NEW(name,type,size) type name##VLA[size]; type* name = (type*)name##VLA
	#define VLA_DELETE(name)
#else
	#define VLA_NEW(name,type,size) type* name = new type[size]
	#define VLA_DELETE(name) delete[] name
#endif

#define NEW(mem,type,...)  new(mem.Alloc(sizeof(type))) type(mem,__VA_ARGS__)
#define NEWW(mem,type,...) new(mem.Alloc(sizeof(type))) type(__VA_ARGS__)

#define GetLoc(x) (slocus[(x)])
		//(useslocus ? slocus[l] : locus[l])
#define GetLocPos(x) (locus_pos[(x)])
		//(useslocus ? locus_pos[l] : locus[l].pos)
#define GetLocId(x) (locus_id[(x)])
		//(useslocus ? locus_id[l] : locus[l].id)

#ifdef HASH64
	typedef uint64 HASH;
#else
	typedef uint HASH;
#endif

#define THREADH(x) \
		TARGET static void x(int id);\
		TARGET static void x##In(void)

#define THREAD(x) \
		TARGET void x(int id)\
		{\
			threadid = (uint)(uint64)id;\
			x##In();\
		}\
		TARGET void x##In(void)

/* Class and Struct names*/
class BCFHEADER;
class VCF;
class MEMORY;
struct GENOTYPE;
class LOCUS;
class VESSEL;
struct POP;
class IND;
template <typename T> class LIST;
template <typename T, typename T2> struct TABLE_ENTRY;
template <typename T, typename T2> class TABLE;
struct Huang2015ENTRY;

#include "dre.h"
#include "ml.h"
#include "ml2.h"
#include "mlbin.h"
#include "mom.h"
#include "mom2.h"
#include "global.h"
#include "hash.h"
#include "math2.h"
#include "mathSSE.h"
#include "mathAVX.h"
#include "math512.h"
#include "misc.h"
#include "file.h"
#include "string2.h"
#include "statistics.h"
#include "matrix.h"
#include "parameters.h"
#include "function.h"
#include "menu.h"


TARGET int main(int _argc, char **_argv);