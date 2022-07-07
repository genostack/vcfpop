/*
 vcfpop v1.0
 -- Perform population genetics analyses based on NGS data for haploids, diploids and polyploids.

 Author      Huang Kang
 Affiliation Northwest University
 Email       huangkang@nwu.edu.cn
 Update      2022/03/31
 */

#include "vcfpop.h"

 /* Main function */
TARGET int main(int _argc, char **_argv)
{
 	GetCurDir(EXEDIR);

	//title
	printf("\nvcfpop v %s   %s\n", VERSION, DATE);
	printf("    Perform population genetics analyses for haploids, diploids, polyploids and ansioploids based on NGS data.\n");
	printf("    Huang Kang, Ph.D., Associate Prof.\n");
	printf("    Northwest University\n");

	//parse parameters
	if (PrintHelp(_argc, _argv)) return 0;

	//load parameters
	SetParameters(_argc, _argv, false);

	ClearTempFiles(g_tmpdir_val);

	//show directories
	printf("    Working directory:%s\n", CURDIR);
	printf("    Temp directory:%s\n\n", g_tmpdir_val);

	threadid = 9999999;
	if (g_benchmark_val == 1) SimdBenchmark();

	clock_t start, end;
	start = clock();

	Calculate();

	end = clock();
	double endtime = (double)(end - start) / CLOCKS_PER_SEC;
	printf("\n\nCalculations finished in %0.3f seconds.\n", endtime);

	ReleaseParameters();
	//_ASSERTE(_CrtCheckMemory());
	return 0;
}