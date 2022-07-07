/* String Functions */

#include "vcfpop.h"

/*Read bcf genotype string to an array */
TARGET void ReadBCFGenoString(char *gt, ushort *alleles, bool &phase, int &v, int asize, int vlen, uint64 l, char *name)
{
	for (; v < vlen; ++v)
	{
		switch (asize)
		{
		case 1:
			alleles[v] = *(byte*)gt++;
			if (alleles[v] == 0x80 || alleles[v] == 0x81) 
				return;
			break;
		case 2:
			alleles[v] = *(ushort*)gt; gt += 2;
			if (alleles[v] == 0x8000 || alleles[v] == 0x8001) 
				return;
			break;
		}
		if (v && (alleles[v] & 0x1) == 0) 
			phase = false;

		alleles[v] = (alleles[v] >> 1) - 1;

		if (name && alleles[v] >= locus[l].k)
			Exit("\nError: allele index exceed number of alleles (%d), in individual %s, variant %s\n", locus[l].k, name, locus[l].GetName());

		if (alleles[v] == 0xFFFF)
		{
			SetFF(alleles, N_MAX_PLOIDY);
			return;
		}
	}
}

/*Read vcf genotype string to an array */
TARGET void ReadVCFGenoString(ushort *alleles, char *genostr, int ploidy, uint64 l, char *name)
{
	for (int i = 0; i < ploidy; ++i)
	{
		int tallele = -1;
		if (*genostr != '.')
		{
			tallele = ReadInteger(genostr);
			if (name && tallele >= locus[l].k)
				Exit("\nError: allele index exceed number of alleles (%d), in individual %s, variant %s\n", locus[l].k, name, locus[l].GetName());
			genostr++;
		}

		if (tallele == -1)
		{
			SetFF(alleles, ploidy);
			break;
		}
		alleles[i] = (ushort)tallele;
	}
}

/* Append val to str */
TARGET void AppendString(char *&str, const char *_val)
{
	char *val = (char*)_val;
	while (*val) *str++ = *val++;
}

/* Read bcf typed int from string */
TARGET uint ReadTypedInt(char *&str)
{
	byte type = *(byte*)str++;
	uint re = 0;
	switch (type & 0xF)
	{
	case 1: re = *(byte*)str; str += 1; break;
	case 2: re = *(ushort*)str; str += 2; break;
	case 3: re = *(uint*)str; str += 4; break;
	default: Exit("\nError: can not read typed integer.\n");
	}
	return re;
}

/* write bcf typed int into string */
TARGET void AppendTypedInt(char *&str, uint val)
{
	if (val >> 8 == 0)
	{
		*str++ = 0x11;
		*str++ = (byte)val;
	}
	else if (val >> 16 == 0)
	{
		*str++ = 0x12;
		*(ushort*)str = (ushort)val;
		str += 2;
	}
	else
	{
		*str++ = 0x13;
		*(uint*)str = val;
		str += 4;
	}
}

/* Fast print an integer to a string */
TARGET void AppendInt(char *&str, int val)
{
	if (val < 0)
	{
		*str++ = '-';
		val = -val;
	}
	int mask = 1000000000;
	while (mask > val) mask /= 10;
	if (!mask)
	{
		*str++ = '0';
		return;
	}
	while (mask)
	{
		*str++ = (char)('0' + (val / mask));
		val %= mask;
		mask /= 10;
	}
}

/* Print with indent */
TARGET void printi(const char *a)
{
	int64 len = strlen(a);
	char abuf[121];
	while (len > 0)
	{
		printf("    ");
		int64 plen = Min(len, 116ll);
		memcpy(abuf, a, plen);
		abuf[plen] = '\0';
		printf((const char*)abuf);
		len -= 116;
		a += 116;
	}
}

/* Trim double quotes in a string */
TARGET char *TrimQuote(char *a)
{
	char *re = a;
	if (re[0] == '\'' || re[0] == '\"')
		re++;

	uint64 len = strlen(re);
	if (re[len - 1] == '\'' || re[len - 1] == '\"')
		re[len - 1] = 0;

	return re;
}

/* Parse range '1-2' */
TARGET void ParseTwoNumber(char *names, uint &v1, uint &v2)
{
	v1 = v2 = 0xFFFFFFFF;
	if (names[0] == '#')
		sscanf(names, "#%d-#%d", &v1, &v2);
}

/* Parse range '1-2' */
TARGET void ParseTwoNumber(char *names, int &v1, int &v2)
{
	v1 = v2 = -1;
	if (names[0] == '#')
		sscanf(names, "#%d-#%d", &v1, &v2);
}

/* Compare two lines, and termined if A ends with linebreak comma | \0 or specified deliminator */
TARGET int LineCmpAterm(const char *a, const char *b, char termin)
{
	for (; ; a++, b++)
	{
		if (*a == termin || *a == '\r' || *a == '\n' || *a == '|' || *a == ',' || !*a)
			return 0;
		if (*a != *b)
			return 1;
	}
}

/* Compare two lines, and termined if A or B end with linebreak comma | \0 or specified deliminator */
TARGET int LineCmpABterm(const char *a, const char *b, char termin)
{
	for (; ; a++, b++)
	{
		if ((*a == termin || *a == '\r' || *a == '\n' || *a == '|' || *a == ',' || !*a) &&
			(*b == termin || *b == '\r' || *b == '\n' || *b == '|' || *b == ',' || !*b))
			return 0;
		if (*a != *b)
			return 1;
	}
}

/* Compare two lines, and termined if A ends with linebreak  |  \0 */
TARGET int LineCmp(const char *a, const char *b)
{
	for (; ; a++, b++)
	{
		if (*a == '\r' || *a == '\n' || *a == '|' || *a == ',' || !*a)
			return 0;
		if (*a != *b)
			return 1;
	}
}

/* Compare the lower case of two lines, and termined if A ends with linebreak or \0 */
TARGET int LwrLineCmp(const char *a, const char *b)
{
	for (; ; a++, b++)
	{
		if (*a == '\r' || *a == '\n' || !*a)
			return 0;
		if (*a + ((*a >= 65 && *a <= 90) ? 32 : 0) !=
			*b + ((*b >= 65 && *b <= 90) ? 32 : 0))
			return 1;
	}
}

/* Compare the lower case of two parameters */
TARGET int LwrParCmp(const char *a, const char *b)
{
	for (; ; a++, b++)
	{
		if ((*a == '\r' || *a == '\n' || !*a || *a == '|' || *a == ',') ||
			(*b == '\r' || *b == '\n' || !*b || *b == '|' || *b == ','))
			return !((*a == '\r' || *a == '\n' || !*a || *a == '|' || *a == ',') &&
				(*b == '\r' || *b == '\n' || !*b || *b == '|' || *b == ','));
		if (*a + ((*a >= 65 && *a <= 90) ? 32 : 0) != *b + ((*b >= 65 && *b <= 90) ? 32 : 0))
			return 1;
	}
}

/* Compare the lower case of two string */
TARGET int LwrStrCmp(const char *a, const char *b)
{
	for (; ; a++, b++)
	{
		if (!*a || !*b)
			return !(!*a && !*b);
		if (*a + ((*a >= 65 && *a <= 90) ? 32 : 0) !=
			*b + ((*b >= 65 && *b <= 90) ? 32 : 0))
			return 1;
	}
}

/* Read all text of a file, allocate memory */
TARGET char *ReadAllText(char *file)
{
	FILE *f1 = fopen(file, "rb");
	if (!f1) Exit("\nError: Cannot open file %s.\n", file);
	fseeko64(f1, 0ll, SEEK_END);
	uint64 filelen = ftello64(f1);
	char *re = new char[filelen + 1];
	fseeko64(f1, 0ll, SEEK_SET);
	fread(re, 1, filelen, f1);
	fclose(f1);
	re[filelen] = '\0';
	return re;
}

/* Replace substrings, allocate memory */
TARGET char *ReplaceStr(char *text, const char *a, const char *b)
{
	char *re = new char[strlen(text) * 2 + 1];
	uint64 p1 = 0, p2 = 0;
	uint64 alen = strlen(a), blen = strlen(b);
	while (text[p1])
	{
		if (text[p1] == a[0] && !memcmp(a, text + p1, alen))
		{
			p1 += alen;
			for (uint64 p3 = 0; p3 < blen; ++p3)
				re[p2++] = b[p3];
		}
		else
			re[p2++] = text[p1++];
	}
	re[p2++] = '\0';
	return re;
}

/* Replace characters, not allocate memory */
TARGET void ReplaceChar(char *text, char a, char b)
{
	for (;;)
	{
		if (*text == a) *text = b;
		if (*++text == '\0') return;
	}
}

/* Count char val in a string */
TARGET int CountChar(char *text, char val)
{
	uint64 count = 0;
	while (*text)
		if (*text++ == val)
			count++;
	return count;
}

/* Count any elements of val in a string */
TARGET int CountChars(char *text, const char *val)
{
	uint64 count = 0;
	for (; *text; text++)
		for (uint64 i = 0; val[i]; ++i)
			if (*text == val[i] && ++count) break;
	return count;
}

/* Count any elements of val in a string */
TARGET int CountChars(char *text, const char *val, int64 len)
{
	uint64 count = 0;
	for (; len; --len, text++)
		for (uint64 i = 0; val[i]; ++i)
			if (*text == val[i] && ++count) break;
	return count;
}

/* Skip rep lines */
TARGET char *LineNextIdx(char *text, const char *_val, int64 rep)
{
	char *val = (char*)_val;
	for (; ; text++)
	{
		if (*text == *val && !LwrLineCmp(val, text) && !--rep) return text;
		if (*text == '\0' || *text == '\n') return NULL;
	}
}

/* Skip rep vals */
TARGET char *StrNextIdx(char *text, const char *_val, int64 rep)
{
	char *val = (char*)_val;
	for (; ; text++)
	{
		if (*text == *val && !LwrLineCmp(val, text) && !--rep) return text;
		if (!*text) return NULL;
	}
}

/* Skip rep vals */
TARGET char *StrNextIdx(char *text, char val, int64 rep)
{
	for (; ; text++)
	{
		if (*text == val && !--rep) return text;
		if (!*text) return NULL;
	}
}

/* Skip rep \0s */
TARGET char *StrNextIdx0(char *text, int64 rep)
{
	for (; ; text++)
		if (*text == '\0' && !--rep) return text;
}

/* Replace sep with \0 and return pointers (new allocated) to each string */
TARGET char **SplitStr(char *text, char sep, int &count)
{
	count = (uint)CountChar(text, sep) + 1;
	char **re = new char*[count];
	SetZero(re, count);
	count = 0;
	bool inquote = false;
	while (*text)
	{
		while (!inquote && *text == sep) text++;
		if (*text == '\"')
		{
			inquote = true;
			memmove(text, text + 1, strlen(text));
		}

		if (*text) re[count++] = text;
		else break;
		while (*text && (inquote || *text != sep))
		{
			if (*text == '\"')
			{
				inquote = !inquote;
				memmove(text, text + 1, strlen(text));
				text--;
			}
			text++;
		}
		if (*text == sep) *text++ = '\0';
		else break;
	}
	return re;
}

/* Replace sep with \0 and return pointers (new allocated) to each string */
TARGET char **SplitStr(char *text, char sep, int64 &count)
{
	count = CountChar(text, sep) + 1;
	char **re = new char*[count];
	SetZero(re, count);
	count = 0;
	bool inquote = false;
	while (*text)
	{
		while (!inquote && *text == sep) text++;
		if (*text == '\"')
		{
			inquote = true;
			memmove(text, text + 1, strlen(text));
		}

		if (*text) re[count++] = text;
		else break;
		while (*text && (inquote || *text != sep))
		{
			if (*text == '\"')
			{
				inquote = !inquote;
				memmove(text, text + 1, strlen(text));
				text--;
			}
			text++;
		}
		if (*text == sep) *text++ = '\0';
		else break;
	}
	return re;
}

/* Convert to lower case characters */
TARGET void StrLwr(char *a)
{
	do
	{
		if (*a >= 65 && *a <= 90)
			*a += 32;
	} while (*a++);
}

/* Read an allele from spagedi files*/
TARGET int ReadIntegerSpagedi(char *&d, int maxdigit)
{
	while ((*d < '0' || *d > '9') && *d != '-') d++;

	int sign = 0, re = 0;
	if (*d == '-')
	{
		d++;
		sign = -1;
	}
	else if (*d >= '0' && *d <= '9')
	{
		sign = 1;
	}
	else
	{
		d++;
		return -1;
	}

	int md = 0;
	while (*d >= '0' && *d <= '9')
	{
		re = re * 10 + (*d++ - '0');
		if (++md == maxdigit) break;
	}

	return re * sign;
}

/* Fast read an integer from a string */
TARGET int64 ReadLong(char *&d)
{
	while ((*d < '0' || *d > '9') && *d != '-') d++;

	int64 sign = 0, re = 0;
	if (*d == '-')
	{
		d++;
		sign = -1;
	}
	else if (*d >= '0' && *d <= '9')
	{
		sign = 1;
	}
	else
	{
		d++;
		return -1;
	}

	while (*d >= '0' && *d <= '9')
		re = re * 10ll + (*d++ - '0');
	return re * sign;
}

/* Parse an integer in binary format */
TARGET uint ReadBinInteger(char *&d, int len)
{
	uint re = 0;
	switch (len)
	{
	case 1: re = *(byte*)d; break;
	case 2: re = *(ushort*)d; break;
	case 3: re = *(uint*)d & 0xFFFFFF; break;
	case 4: re = *(uint*)d; break;
	}
	d += len;
	return re;
}

/* Parse a integer from string and move pointer */
TARGET int ReadInteger(char *&d)
{
	while ((*d < '0' || *d > '9') && *d != '-' && *d != '?') d++;
	
	if (*d == '?')
	{
		while (*d++ == '?');
		return -1;
	}

	int sign = 0, re = 0;
	if (*d == '-')
	{
		d++;
		sign = -1;
	}
	else if (*d >= '0' && *d <= '9')
	{
		sign = 1;
	}
	else
	{
		d++;
		return -1;
	}

	while (*d >= '0' && *d <= '9')
		re = re * 10 + (*d++ - '0');
	return re * sign;
}

/* Parse a integer from string and do not move pointer */
TARGET int ReadIntegerKeep(char *d)
{
	while ((*d < '0' || *d > '9') && *d != '-') d++;

	int sign = 0, re = 0;
	if (*d == '-')
	{
		d++;
		sign = -1;
	}
	else if (*d >= '0' && *d <= '9')
	{
		sign = 1;
	}
	else
	{
		d++;
		return -1;
	}

	while (*d >= '0' && *d <= '9')
		re = re * 10 + (*d++ - '0');
	return re * sign;
}

/* Parse a real number from string and move pointer */
TARGET double ReadDouble(char *&d)
{
	int sign = 0;
	double re = 0, mask = 1;
	if (*d == '-')
	{
		d++;
		sign = -1;
	}
	else if (*d >= '0' && *d <= '9')
		sign = 1;
	else
	{
		d++;
		return -1;
	}

	while (*d >= '0' && *d <= '9')
		re = re * 10 + (*d++ - '0');
	if (*d == '.')
	{
		d++;
		while (*d >= '0' && *d <= '9')
		{
			mask *= 0.1;
			re = re + (*d++ - '0') * mask;
		}
	}
	return re * sign;
}

/* Parse a real number from string and do not move pointer */
TARGET double ReadDoubleKeep(char *d)
{
	int sign = 0;
	double re = 0, mask = 1;
	if (*d == '-')
	{
		d++;
		sign = -1;
	}
	else if (*d >= '0' && *d <= '9')
	{
		sign = 1;
	}
	else
	{
		d++;
		return -1;
	}

	while (*d >= '0' && *d <= '9')
		re = re * 10 + (*d++ - '0');
	if (*d == '.')
	{
		d++;
		while (*d >= '0' && *d <= '9')
		{
			mask *= 0.1;
			re = re + (*d++ - '0') * mask;
		}
	}
	return re * sign;
}

/* Read real range parameter */
TARGET void GetRangeParDouble(char *gpar, bool &parid, double &minv, double &maxv, double min, double max)
{
	if (parid) Exit("\nError: parameter %s has been assigned twice.\n", gpar);
	parid = true;
	int s1 = 0, e1 = 0, s2 = 0, e2 = 0;
	char b = '\0';
	while (gpar[s1] && gpar[s1] != '[') s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	e1 = s1 + 1;
	while (gpar[e1] && gpar[e1] != ',') e1++;
	if (!gpar[e1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	b = gpar[e1];
	gpar[e1] = '\0';
	sscanf(gpar + s1, "%lf", &minv);
	gpar[e1] = b;

	s2 = e1 + 1;
	e2 = s2 + 1;
	while (gpar[e2] && gpar[e2] != ']') e2++;
	if (!gpar[e2]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	b = gpar[e2];
	gpar[e2] = '\0';
	sscanf(gpar + s2, "%lf", &maxv);
	gpar[e2] = b;
	if (minv > maxv || minv < min || maxv > max)
		Exit("\nError: parameter %s out of range, whose value should between %lf and %lf.\n", gpar, min, max);
}

/* Read 32bit integer range parameter */
TARGET void GetRangeParInteger(char *gpar, bool &parid, uint &minv, uint &maxv, uint min, uint max)
{
	if (parid) Exit("\nError: parameter %s has been assigned twice.\n", gpar);
	parid = true;
	int s1 = 0, e1 = 0, s2 = 0, e2 = 0;
	char b = '\0';
	while (gpar[s1] && gpar[s1] != '[') s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	e1 = s1 + 1;
	while (gpar[e1] && gpar[e1] != ',') e1++;
	if (!gpar[e1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	b = gpar[e1];
	sscanf(gpar + s1, "%u", &minv);
	gpar[e1] = b;

	s2 = e1 + 1;
	e2 = s2 + 1;
	while (gpar[e2] && gpar[e2] != ']') e2++;
	if (!gpar[e2]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	b = gpar[e2];
	sscanf(gpar + s2, "%u", &maxv);
	gpar[e2] = b;
	if (minv > maxv || minv < min || maxv > max)
		Exit("\nError: parameter %s out of range, whose value should between %u and %u.\n", gpar, min, max);
}

/* Read 32bit integer range parameter */
TARGET void GetRangeParInteger(char *gpar, bool &parid, int &minv, int &maxv, int min, int max)
{
	if (parid) Exit("\nError: parameter %s has been assigned twice.\n", gpar);
	parid = true;
	int s1 = 0, e1 = 0, s2 = 0, e2 = 0;
	char b = '\0';
	while (gpar[s1] && gpar[s1] != '[') s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	e1 = s1 + 1;
	while (gpar[e1] && gpar[e1] != ',') e1++;
	if (!gpar[e1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	b = gpar[e1];
	sscanf(gpar + s1, "%d", &minv);
	gpar[e1] = b;

	s2 = e1 + 1;
	e2 = s2 + 1;
	while (gpar[e2] && gpar[e2] != ']') e2++;
	if (!gpar[e2]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	b = gpar[e2];
	sscanf(gpar + s2, "%d", &maxv);
	gpar[e2] = b;
	if (minv > maxv || minv < min || maxv > max)
		Exit("\nError: parameter %s out of range, whose value should between %d and %d.\n", gpar, min, max);
}

/* Read 64bit integer range parameter */
TARGET void GetRangeParLong(char *gpar, bool &parid, uint64 &minv, uint64 &maxv, uint64 min, uint64 max)
{
	if (parid) Exit("\nError: parameter %s has been assigned twice.\n", gpar);
	parid = true;
	int s1 = 0, e1 = 0, s2 = 0, e2 = 0;
	char b = '\0';
	while (gpar[s1] && gpar[s1] != '[') s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	e1 = s1 + 1;
	while (gpar[e1] && gpar[e1] != ',') e1++;
	if (!gpar[e1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	b = gpar[e1];
	sscanf(gpar + s1, "%llu", &minv);
	gpar[e1] = b;

	s2 = e1 + 1;
	e2 = s2 + 1;
	while (gpar[e2] && gpar[e2] != ']') e2++;
	if (!gpar[e2]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	b = gpar[e2];
	sscanf(gpar + s2, "%llu", &maxv);
	gpar[e2] = b;
	if (minv > maxv || minv < min || maxv > max)
		Exit("\nError: parameter %s out of range, whose value should between %llu and %llu.\n", gpar, min, max);
}

/* Read 64bit integer range parameter */
TARGET void GetRangeParLong(char *gpar, bool &parid, int64 &minv, int64 &maxv, int64 min, int64 max)
{
	if (parid) Exit("\nError: parameter %s has been assigned twice.\n", gpar);
	parid = true;
	int s1 = 0, e1 = 0, s2 = 0, e2 = 0;
	char b = '\0';
	while (gpar[s1] && gpar[s1] != '[') s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	e1 = s1 + 1;
	while (gpar[e1] && gpar[e1] != ',') e1++;
	if (!gpar[e1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	b = gpar[e1];
	sscanf(gpar + s1, "%lld", &minv);
	gpar[e1] = b;

	s2 = e1 + 1;
	e2 = s2 + 1;
	while (gpar[e2] && gpar[e2] != ']') e2++;
	if (!gpar[e2]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	b = gpar[e2];
	sscanf(gpar + s2, "%lld", &maxv);
	gpar[e2] = b;
	if (minv > maxv || minv < min || maxv > max)
		Exit("\nError: parameter %s out of range, whose value should between %lld and %lld.\n", gpar, min, max);
}

/* Read real parameter */
TARGET void GetParDouble(char *gpar, bool &parid, double &val, double min, double max)
{
	if (parid) Exit("\nError: parameter %s has been assigned twice.\n", gpar);
	parid = true;
	int s1 = 0;
	while (gpar[s1] && gpar[s1] != '=') s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);

	sscanf(gpar + s1, "%lf", &val);
	if (val < min || val > max)
		Exit("\nError: parameter %s out of range, whose value should between %lf and %lf.\n", gpar, min, max);
}

/* Read 32bit integer parameter */
TARGET void GetParInteger(char *gpar, bool &parid, uint &val, uint min, uint max)
{
	if (parid) Exit("\nError: parameter %s has been assigned twice.\n", gpar);
	parid = true;
	int s1 = 0;
	while (gpar[s1] && gpar[s1] != '=') s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);

	sscanf(gpar + s1, "%u", &val);
	if (val < min || val > max)
		Exit("\nError: parameter %s out of range, whose value should between %u and %u.\n", gpar, min, max);
}

/* Read 32bit integer parameter */
TARGET void GetParInteger(char *gpar, bool &parid, int &val, int min, int max)
{
	if (parid) Exit("\nError: parameter %s has been assigned twice.\n", gpar);
	parid = true;
	int s1 = 0;
	while (gpar[s1] && gpar[s1] != '=') s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);

	sscanf(gpar + s1, "%d", &val);
	if (val < min || val > max)
		Exit("\nError: parameter %s out of range, whose value should between %d and %d.\n", gpar, min, max);
}

/* Read 64bit integer parameter */
TARGET void GetParLong(char *gpar, bool &parid, uint64 &val, uint64 min, uint64 max)
{
	if (parid) Exit("\nError: parameter %s has been assigned twice.\n", gpar);
	parid = true;
	int s1 = 0;
	while (gpar[s1] && gpar[s1] != '=') s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	sscanf(gpar + s1, "%llu", &val);
	if (val < min || val > max)
		Exit("\nError: parameter %s out of range, whose value should between %llu and %llu.\n", gpar, min, max);
}

/* Read 64bit integer parameter */
TARGET void GetParLong(char *gpar, bool &parid, int64 &val, int64 min, int64 max)
{
	if (parid) Exit("\nError: parameter %s has been assigned twice.\n", gpar);
	parid = true;
	int s1 = 0;
	while (gpar[s1] && gpar[s1] != '=') s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	s1++;
	if (!gpar[s1]) Exit("\nError: cannot parse parameter %s, check format.\n", gpar);
	sscanf(gpar + s1, "%lld", &val);
	if (val < min || val > max)
		Exit("\nError: parameter %s out of range, whose value should between %lld and %lld.\n", gpar, min, max);
}

/* Read boolean parameter */
TARGET void GetParBool(char *gpar, bool &parid)
{
	if (parid)
		Exit("\nError: parameter %s has been assigned twice.\n", gpar);
	parid = true;
}

/* Read string parameter */
TARGET void GetParString(char *gpar, const char *ref, bool &parid, int &val)
{
	if (parid) Exit("\nError: parameter %s has been assigned twice.\n", gpar);
	parid = true;

	char *gparb = gpar;
	while (*gpar && *gpar != '=') gpar++;

	if (*gpar == '=')
		gpar++;
	else
		Exit("\nError: cannot parse parameter %s, check format.\n", gpar);

	val = -1;
	int refp = 0, refi = 0;
	do
	{
		refi++;

		if (!LwrParCmp(ref + refp, gpar))
		{
			val = refi;
			break;
		}

		while (ref[refp] && ref[refp] != '|') refp++;
		if (ref[refp] == '|') refp++;
	} while (ref[refp]);
	if (val == -1)
		Exit("\nError: Unrecognized parameter: %s\n, ", gparb);
}

/* Read multiple section parameters */
TARGET void GetParStringMultiSel(char *gpar, const char *ref, bool &parid, byte *val)
{
	if (parid) Exit("\nError: parameter %s has been assigned twice.\n", gpar);
	parid = true;

	char *gparb = gpar;
	memset(val, 0, N_MAX_OPTION);
	while (*gpar && *gpar != '=') gpar++;

	if (*gpar == '=')
		gpar++;
	else
		Exit("\nError: cannot parse parameter %s, check format.\n", gpar);

	while (*gpar)
	{
		int refp = 0, refi = 0;
		val[refi] = 0;
		do
		{
			refi++;

			if (!LwrParCmp(ref + refp, gpar))
			{
				val[refi] = 1;
				break;
			}

			while (ref[refp] && ref[refp] != '|') refp++;
			if (ref[refp] == '|') refp++;
		} while (ref[refp]);

		if (val[refi] == 0)
		{
			char *commapos = StrNextIdx(gpar, ",", 1);
			if ((uint64)commapos > 0) *commapos = 0;
			Exit("\nError: Unrecognized option: %s in parameter: \n%s\n.", gpar, gparb);
		}
		while (*gpar && *gpar != ',') gpar++;
		if (*gpar == ',') gpar++;
		else break;
	}
}

/* Print a real number to a file */
TARGET void WriteReal(FILE *fout, double val)
{
	if (IsError(val)) fprintf(fout, "nan");
	else fprintf(fout, g_decimal_str, val);
}

/* Print a real number to a string */
TARGET void WriteReal(char *&fout, double val)
{
	if (IsError(val)) sprintf(fout, "nan");
	else sprintf(fout, g_decimal_str, val);
	fout += strlen(fout);
}

/* Append a real number to a string */
TARGET void AppendReal(char *sout, double val)
{
	if (IsError(val)) sprintf(sout, "nan");
	else sprintf(sout, g_decimal_str, val);
}