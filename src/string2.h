/* String Functions */

#pragma once
#include "vcfpop.h"

/*Read bcf genotype string to an array */
TARGET void ReadBCFGenoString(char *gt, ushort *alleles, bool &phase, int &v, int asize, int vlen, uint64 l, char *name);

/*Read vcf genotype string to an array */
TARGET void ReadVCFGenoString(ushort *alleles, char *genostr, int ploidy, uint64 l, char *name);

/* Append val to str */
TARGET void AppendString(char *&str, const char *_val);

/* Read bcf typed int from string */
TARGET uint ReadTypedInt(char *&str);

/* write bcf typed int into string */
TARGET void AppendTypedInt(char *&str, uint val);

/* Fast print an integer to a string */
TARGET void AppendInt(char *&str, int val);

/* Print with indent */
TARGET void printi(const char *a);

/* Trim double quotes in a string */
TARGET char *TrimQuote(char *a);

/* Parse range '1-2' */
TARGET void ParseTwoNumber(char *names, uint &v1, uint &v2);

/* Parse range '1-2' */
TARGET void ParseTwoNumber(char *names, int &v1, int &v2);

/* Compare two lines, and termined if A ends with linebreak comma | \0 or specified deliminator */
TARGET int LineCmpAterm(const char *a, const char *b, char termin);

/* Compare two lines, and termined if A or B end with linebreak comma | \0 or specified deliminator */
TARGET int LineCmpABterm(const char *a, const char *b, char termin);

/* Compare two lines, and termined if A ends with linebreak  |  \0 */
TARGET int LineCmp(const char *a, const char *b);

/* Compare the lower case of two lines, and termined if A ends with linebreak or \0 */
TARGET int LwrLineCmp(const char *a, const char *b);

/* Compare the lower case of two parameters */
TARGET int LwrParCmp(const char *a, const char *b);

/* Compare the lower case of two string */
TARGET int LwrStrCmp(const char *a, const char *b);

/* Read all text of a file, allocate memory */
TARGET char *ReadAllText(char *file);

/* Replace substrings, allocate memory */
TARGET char *ReplaceStr(char *text, const char *a, const char *b);

/* Replace characters, not allocate memory */
TARGET void ReplaceChar(char *text, char a, char b);

/* Count char val in a string */
TARGET int CountChar(char *text, char val);

/* Count any elements of val in a string */
TARGET int CountChars(char *text, const char *val);

/* Count any elements of val in a string */
TARGET int CountChars(char *text, const char *val, int64 len);

/* Skip rep lines */
TARGET char *LineNextIdx(char *text, const char *_val, int64 rep);

/* Skip rep vals */
TARGET char *StrNextIdx(char *text, const char *_val, int64 rep);

/* Skip rep vals */
TARGET char *StrNextIdx(char *text, char val, int64 rep);

/* Skip rep \0s */
TARGET char *StrNextIdx0(char *text, int64 rep);

/* Replace sep with \0 and return pointers (new allocated) to each string */
TARGET char **SplitStr(char *text, char sep, int &count);

/* Replace sep with \0 and return pointers (new allocated) to each string */
TARGET char **SplitStr(char *text, char sep, int64 &count);

/* Convert to lower case characters */
TARGET void StrLwr(char *a);

/* Read an allele from spagedi files*/
TARGET int ReadIntegerSpagedi(char *&d, int maxdigit);

/* Fast read an integer from a string */
TARGET int64 ReadLong(char *&d);

/* Parse an integer in binary format */
TARGET uint ReadBinInteger(char *&d, int len);

/* Parse a integer from string and move pointer */
TARGET int ReadInteger(char *&d);

/* Parse a integer from string and do not move pointer */
TARGET int ReadIntegerKeep(char *d);

/* Parse a real number from string and move pointer */
TARGET double ReadDouble(char *&d);

/* Parse a real number from string and do not move pointer */
TARGET double ReadDoubleKeep(char *d);

/* Read real range parameter */
TARGET void GetRangeParDouble(char *gpar, bool &parid, double &minv, double &maxv, double min, double max);

/* Read 32bit integer range parameter */
TARGET void GetRangeParInteger(char *gpar, bool &parid, uint &minv, uint &maxv, uint min, uint max);

/* Read 32bit integer range parameter */
TARGET void GetRangeParInteger(char *gpar, bool &parid, int &minv, int &maxv, int min, int max);

/* Read 64bit integer range parameter */
TARGET void GetRangeParLong(char *gpar, bool &parid, uint64 &minv, uint64 &maxv, uint64 min, uint64 max);

/* Read 64bit integer range parameter */
TARGET void GetRangeParLong(char *gpar, bool &parid, int64 &minv, int64 &maxv, int64 min, int64 max);

/* Read real parameter */
TARGET void GetParDouble(char *gpar, bool &parid, double &val, double min, double max);

/* Read 32bit integer parameter */
TARGET void GetParInteger(char *gpar, bool &parid, uint &val, uint min, uint max);

/* Read 32bit integer parameter */
TARGET void GetParInteger(char *gpar, bool &parid, int &val, int min, int max);

/* Read 64bit integer parameter */
TARGET void GetParLong(char *gpar, bool &parid, uint64 &val, uint64 min, uint64 max);

/* Read 64bit integer parameter */
TARGET void GetParLong(char *gpar, bool &parid, int64 &val, int64 min, int64 max);

/* Read boolean parameter */
TARGET void GetParBool(char *gpar, bool &parid);

/* Read string parameter */
TARGET void GetParString(char *gpar, const char *ref, bool &parid, int &val);

/* Read multiple section parameter */
TARGET void GetParStringMultiSel(char *gpar, const char *ref, bool &parid, byte *val);

/* Print a real number to a file */
TARGET void WriteReal(FILE *fout, double val);

/* Print a real number to a string */
TARGET void WriteReal(char *&fout, double val);

/* Append a real number to a string */
TARGET void AppendReal(char *sout, double val);