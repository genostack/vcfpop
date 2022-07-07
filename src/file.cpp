/* File IO */

#include "vcfpop.h"

/* Get file size */
TARGET int64 GetFileLen(FILE *file)
{
	int64 pos = ftello64(file);
	fseeko64(file, 0, SEEK_END);
	int64 re = ftello64(file);
	fseeko64(file, pos, SEEK_SET);
	return re;
}

/* Read a double from file */
TARGET double FGetDouble(FILE *file)
{
	double buf;
	if (g_format_val < 0)
		gzread((gzFile)file, &buf, 8);
	else
		fread(&buf, 1, 8, file);
	return buf;
}

/* Read a float from file */
TARGET float FGetFloat(FILE *file)
{
	float buf;
	if (g_format_val < 0)
		gzread((gzFile)file, &buf, 4);
	else
		fread(&buf, 1, 4, file);
	return buf;
}

/* Read an int64 from file */
TARGET int64 FGetInt64(FILE *file)
{
	int64 buf;
	if (g_format_val < 0)
		gzread((gzFile)file, &buf, 8);
	else
		fread(&buf, 1, 8, file);
	return buf;
}

/* Read an uint64 from file */
TARGET uint64 FGetUint64(FILE *file)
{
	uint64 buf;
	if (g_format_val < 0)
		gzread((gzFile)file, &buf, 8);
	else
		fread(&buf, 1, 8, file);
	return buf;
}

/* Read an int from file */
TARGET int FGetInt(FILE *file)
{
	int buf;
	if (g_format_val < 0)
		gzread((gzFile)file, &buf, 4);
	else
		fread(&buf, 1, 4, file);
	return buf;
}

/* Read an uint from file */
TARGET uint FGetUint(FILE *file)
{
	uint buf;
	if (g_format_val < 0)
		gzread((gzFile)file, &buf, 4);
	else
		fread(&buf, 1, 4, file);
	return buf;
}

/* Read a short from file */
TARGET short FGetShort(FILE *file)
{
	short buf;
	if (g_format_val < 0)
		gzread((gzFile)file, &buf, 2);
	else
		fread(&buf, 1, 2, file);
	return buf;
}

/* Read an ushort from file */
TARGET ushort FGetUshort(FILE *file)
{
	ushort buf;
	if (g_format_val < 0)
		gzread((gzFile)file, &buf, 2);
	else
		fread(&buf, 1, 2, file);
	return buf;
}

/* Read a byte from file */
TARGET byte FGetByte(FILE *file)
{
	byte buf;
	if (g_format_val < 0)
		gzread((gzFile)file, &buf, 1);
	else
		fread(&buf, 1, 1, file);
	return buf;
}

/* Read a char from file */
TARGET char FGetChar(FILE *file)
{
	char buf;
	if (g_format_val < 0)
		gzread((gzFile)file, &buf, 1);
	else
		fread(&buf, 1, 1, file);
	return buf;
}

/* Read a typed int from bcf file */
TARGET uint FGetTypedInt(FILE *file)
{
	byte type = FGetByte(file);
	switch (type & 0xF)
	{
	case 1: return FGetByte(file);
	case 2: return FGetUshort(file);
	case 3: return FGetUint(file);
	}
	return 0;
}

/* Read data from file */
TARGET int64 FGets(char *&buf, int64 &buflen, FILE *file)
{
	for (int64 nread = 0, pos = FTell(file);;)
	{
		buf[0] = '\0';
		if (g_format_val < 0)
			gzgets((gzFile)file, buf, (int)buflen);
		else
			fgets(buf, (int)buflen, file);

		nread = FTell(file) - pos;
		if (nread >= buflen - 1)
		{
			delete[] buf;
			buflen <<= 1;
			buf = new char[buflen];
			FSeek(file, pos, SEEK_SET);
		}
		else
			return nread;
	} 
}

/* Read data from file and extended insufficient buffer if necessary */
TARGET int64 FGets(char *&buf, int64 &buflen, FILE *file, char *&base)
{
	for (int64 nread = 0, reslen = buflen - (buf - base), pos = FTell(file);;)
	{
		buf[0] = '\0';
		if (g_format_val < 0)
			gzgets((gzFile)file, buf, (int)reslen);
		else
			fgets(buf, (int)reslen, file);
		nread = FTell(file) - pos;
		if (nread >= reslen - 1)
		{
			buflen <<= 1;
			char *base2 = new char[buflen];
			memcpy(base2, base, buf - base);
			reslen = buflen - (buf - base);
			delete[] base;
			buf = base2 + (buf - base);
			base = base2;
			FSeek(file, pos, SEEK_SET);
		}
		else
			return nread;
	} 
	
}

TARGET int64 FGetsSmall(char *buf, int64 buflen, FILE *file)
{
	int64 pos = FTell(file);
	buf[0] = '\0';
	if (g_format_val < 0)
		gzgets((gzFile)file, buf, (int)buflen);
	else
		fgets(buf, (int)buflen, file);

	return FTell(file) - pos;
}

/* Is end of file */
TARGET bool FEof(FILE *file)
{
	if (g_format_val < 0)
		return gzeof((gzFile)file);
	else
		return feof(file);
}

/* Get compressed data offset */
TARGET int64 FOffset(FILE *file)
{
	if (g_format_val < 0)
		return gzoffset64((gzFile)file);
	else
		return ftello64(file);
}

/* Get uncompressed data offset */
TARGET int64 FTell(FILE *file)
{
	if (g_format_val < 0)
		return gztell64((gzFile)file);
	else
		return ftello64(file);
}

/* Read uncompressed data */
TARGET int64 FRead(void *buf, uint64 element, uint64 count, FILE *file)
{
	if (g_format_val < 0)
		return gzread((gzFile)file, buf, (uint)(element * count));
	else
		return fread(buf, element, count, file);
}

/* Change file offset for uncompressed data */
TARGET int64 FSeek(FILE *file, int64 offset, int origin)
{
	if (g_format_val < 0)
		return gzseek64((gzFile)file, offset, origin);
	else
		return fseeko64(file, offset, origin);
}

/* Open an uncompressed or compressed file */
TARGET FILE *FOpen(const char *file, const char *spec)
{
	if (g_format_val < 0)
	{
		gzFile f = gzopen(file, spec);
		gzbuffer(f, LINE_BUFFER);
		return (FILE*)f;
	}
	else
		return fopen(file, spec);
}

/* Close an uncompressed or compressed file */
TARGET int FClose(FILE *file)
{
	if (g_format_val < 0)
		return gzclose((gzFile)file);
	else
		return fclose(file);
}

/* Format file name and open file */
TARGET FILE *FOpen(char *buf, const char *type, const char *fmt ...)
{
	va_list ap;
	int n = 0;
	va_start(ap, fmt);
	n = vsprintf(buf, fmt, ap);
	va_end(ap);
	FILE *f1 = fopen(buf, type);
	if (!f1) Exit("\nError: Cannot open file %s.\n", buf);
	return f1;
}

/* Open result file and write parameters */
TARGET void OpenResFile(const char *_spec, const char *title)
{
	char *spec = (char*)_spec;
	time_t tt1;
	time(&tt1);
	FRES_TIME = localtime(&tt1);
	FRES = FOpen(FRES_NAME, "wb", "%s.%s.txt", g_output_val, spec + 1);
	setvbuf(FRES, FRES_BUF, _IOFBF, LINE_BUFFER);
	const char *prefix = title[0] == '#' ? "#" : "";

	fprintf(FRES, "%s%s calculated by vcfpop v%s%s", prefix, title, VERSION, g_linebreak_val);
	fprintf(FRES, "%s  Input: %s%s", prefix, g_input_val, g_linebreak_val);
	fprintf(FRES, "%s  Output: %s%s", prefix, g_output_val, g_linebreak_val);
	fprintf(FRES, "%s  Time: %04d-%02d-%02d %02d:%02d:%02d%s", prefix, FRES_TIME->tm_year + 1900, FRES_TIME->tm_mon + 1, FRES_TIME->tm_mday, FRES_TIME->tm_hour, FRES_TIME->tm_min, FRES_TIME->tm_sec, g_linebreak_val);
	WriteParameters(FRES, spec, (char*)prefix);
}

/* Close result file and write parameters */
TARGET void CloseResFile()
{
	fclose(FRES);
}

/* Open temp files */
TARGET void OpenTempFiles(uint n, const char *spec, byte flag[])
{
	TEMP_NAMES = new char*[n];
	TEMP_FILES = new FILE*[n];
	TEMP_BUFS = new char*[n];

	for (uint i = 0; i < n; ++i)
	{
		if (flag && !flag[i]) continue;
		TEMP_BUFS[i] = new char[LINE_BUFFER];
		TEMP_NAMES[i] = new char[FILE_NAME_LEN];
		TEMP_FILES[i] = FOpen(TEMP_NAMES[i], "wb+", "%svcfpop_%04d_%02d_%02d_%02d_%02d_%02d%s%d%s",
			g_tmpdir_val, FRES_TIME->tm_year + 1900, FRES_TIME->tm_mon + 1, FRES_TIME->tm_mday, FRES_TIME->tm_hour, FRES_TIME->tm_min, FRES_TIME->tm_sec, spec, i + 1, ".tmp");
		setvbuf(TEMP_FILES[i], TEMP_BUFS[i], _IOFBF, LINE_BUFFER);
	}
}

/* Merge and close temp files */
TARGET void JoinTempFiles(uint n, byte flag[])
{
	int64 buflen = LINE_BUFFER;
	char *buf = new char[buflen];
	for (uint i = 0; i < n; ++i)
	{
		if (flag && !flag[i]) continue;
		fseeko64(TEMP_FILES[i], 0ll, SEEK_SET);
		while (!feof(TEMP_FILES[i]))
		{
			int64 t = fread(buf, 1, buflen, TEMP_FILES[i]);
			fwrite(buf, 1, t, FRES);
		}
		fclose(TEMP_FILES[i]);
		remove(TEMP_NAMES[i]);
		delete[] TEMP_NAMES[i];
		delete[] TEMP_BUFS[i];
	}
	delete[] TEMP_BUFS;
	delete[] TEMP_NAMES;
	delete[] TEMP_FILES;
	delete[] buf;
}