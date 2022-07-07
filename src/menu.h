/* Menu */

#pragma once
#include "vcfpop.h"

/* Print help information */
bool PrintHelp(uint _argc, char **_argv);

/* Run benchmark for SIMD instruction set */
TARGET void SimdBenchmark();
