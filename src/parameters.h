/* Parameter Functions */

#pragma once
#include "vcfpop.h"

/* Initialize parameters */
TARGET void SetDefaultParameters();

/* Set parameters */
TARGET void SetParameters(int _argc, char **_argv, bool isinparfile);

/* Release parameter memory */
TARGET void ReleaseParameters();

/* Write parameters to result file */
TARGET void WriteParameters(FILE *f1, char *type, char *prefix);