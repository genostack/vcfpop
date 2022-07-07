/* Functions */

#include "vcfpop.h"

#pragma pack(push, 1)


/* Lock-free thread competiation */
#define GUARD_BEGIN     \
						while (state_lock[ii % NBUF] >> 1 != ii * 2 + 1) \
							SLEEP(SLEEP_TIME_TINY);

#define GUARD_END		\
						int64 test_val = ii * 4 + 4, next_val = (ii + NBUF) * 4; \
						state_lock[ii % NBUF]++; \
						state_lock[ii % NBUF].compare_exchange_strong(test_val, next_val);

#define THREAD_BEGIN    \
						while (ii >= progress1 + NBUF)\
							SLEEP(SLEEP_TIME_TINY);\
						int64 avail_val = ii << 2, calc_val = avail_val + 1;\
						if (state_lock[ii % NBUF].compare_exchange_strong(avail_val, calc_val)){

#define THREAD_END(x)    \
						state_lock[ii % NBUF] += x;}

#define THREAD_BEGIN2    \
						while (ii >= progress1 + NBUF || ii >= progress2 + NBUF)\
							SLEEP(SLEEP_TIME_TINY);\
						int64 avail_val = ii << 2, calc_val = avail_val + 1;\
						if (state_lock[ii % NBUF].compare_exchange_strong(avail_val, calc_val)){

#define sfprintf(f,...) {sprintf(f,__VA_ARGS__); f += strlen(f);}

/* DEBUG */
TARGET void CheckGenotypeId()
{
#ifdef _DEBUG
	for (int64 l = 0; l < nloc; ++l)
	{
		int ngeno = (useslocus ? slocus[l] : locus[l]).ngeno;
		GENO_ITERATOR iter(0u, l, true);//OK
		for (int i = 0; i < nind; ++i)
		{
			int gid1 = iter.Read();
			int gid2 = ainds[i]->GetGenotypeId(l);
			if (gid1 != gid2 || gid1 >= ngeno || gid1 < 0)
				Exit("\nError, genotype id mismatch.");
		}
	}
#endif
}

#define extern

/* Thread-specific variables */

extern _thread MEMORY *amova_memory;				//memory class for amova vessels
extern _thread double *Anderson2007_Coef;			//Anderson 2007 relatedness estimator coefficients
extern _thread double *Huang2015_Coef;				//Huang 2015 relatedness estimator coefficients
extern _thread int threadid;						//Thread index
extern _thread bool *spa_valid;
extern _thread int spa_tn;

/* Functions */
extern TABLE<int, Huang2015ENTRY> *Huang2015_maps;	//Huang2015_maps[ploidylevel][hash] is a entry saves the ibs modex index and genotype pair pattern

/* Hash */
extern uint *cryptTable;							//crypt table for calculate hash for genotypes

/* Functions */
extern int NBUF;									//CALC_THREAD_BUFFER * g_nthread_val;
extern MEMORY *individual_memory;					//Individual memory class
extern MEMORY *locus_memory;						//Locus memory class
extern MEMORY *nvcf_memory;							//Locus memory for first round counting ngeno
extern TABLE<HASH, uint> *nvcf_gfid;				//Hash table for counting ngeno 
extern MEMORY *conversion_memory;					//Memory class for conversion_string
extern MEMORY *conversion_memory2;					//Memory class for genotype_string and convert_buf
extern LIST<char*> *conversion_string;				//Genotype string for each genotype for file conversion
extern MEMORY *gd_memory;							//Genetic distance memory class

extern TABLE<HASH, INDGD> *gdtab;					//Hash table saves the genetic distance between genotypes
extern shared_mutex* gdlock;						//gdtab lock2
extern double* amova_gd;							//Genetic distance matrix used in AMOVA

extern byte *genotype_bucket;						//Genotype index data bucket
extern uint64 genotype_size;						//Size of bucket in bytes
extern uint64 genotype_coffset;						//Size of used bucket memory
extern OFFSET *genotype_offset;						//Genotype index data at each locus

extern ushort missing_array[N_MAX_PLOIDY];			//Allele array of the missing genotypes
extern GENOTYPE missing_genotype[N_MAX_PLOIDY + 1];	//Missing genotype at different ploidy level
extern HASH missing_hash[N_MAX_PLOIDY + 1];			//Hash of missing genotype

/* Allelic depth, TEST*/
extern byte *alleledepth_bucket;					//Allele depth data bucket
extern uint64 alleledepth_size;						//Size of bucket in bytes
extern uint64 alleledepth_coffset;					//Size of used bucket memory
extern OFFSET *alleledepth_offset;					//Allele depth at each locus

extern LOCN *allele_freq_offset;					//Allele frequency array at locus l
extern int maxK;									//Max number of alleles
extern int64 KT;									//Total number of alleles

extern LOCN *genotype_count_offset;					//Genotype count array at locus l
extern int maxG;									//Max number of genotypes at all loci
extern int64 GT;									//Total number of genotypes

extern char **load_buf;								//Circle buffer for loading vcf/bcf files, NBUF
extern int *load_buf_size;							//Size of each circle buffer for loading, NBUF
extern char *vcf_header;							//Vcf header row

extern LIST<POP> pop;								//Original input population
extern LIST<LIST<POP>> reg;							//Original input regions
extern LOCUS *locus;								//Locus information
extern SLOCUS *slocus;								//SLocus information
extern bool useslocus;								//Use small locus
extern uint64 *locus_pos;							//Locus pos for small locus used in haplotype extraction
extern LOCN *locus_id;								//Locus id for small locus used in haplotype extraction
extern TABLE<HASH, GENOTYPE*> emptygftab;			//Empty genotype table

extern IND **ainds;									//Individuals
extern POP **apops;									//Rearranged populations
extern POP ***aregs;								//Rearranged regions

extern int npop;									//Number of populations
extern int lreg;									//Level of regions
extern int nreg[N_MAX_REG];							//Number of regions in each level
extern int nregt;									//Total number of regions
extern int nregt2;									//Total number of region pairs across levels
extern POP *total_pop;								//Total population

extern LOCN nloc;									//Number of loci
extern LOCN nfilter;								//Number of filtered loci
extern int nind;									//Number of individuals
extern atomic<int> nfilterind;						//Number of filtered individuals
extern bool reassigned;								//Is ploidy assigned, to distiguish diversity filter and diversity estimation
extern IND **rinds;							//Rearranged individuals according to population source

extern POP *cpop;									//Current population

extern byte maxploidy;								//Max ploidy in all genotypes
extern byte minploidy;								//Min ploidy in all genotypes

extern int64 progress1, progress2;					//Progress value for read / write in a circle buffer
extern atomic<int64> *state_lock;					//State of circle buffer

extern atomic<int64> haplotype_contig;				//1st locus in the current contig 
extern OFFSET *haplotype_offset;					//Offset for extract dummy locus, nl elements
extern byte *haplotype_bucket;						//Bucket of extract dummy genotype id
extern uint64 haplotype_size;						//Bucket size
extern LIST<QSLPAR> qslstack;						//Sort loci in haplotype extraction
extern LIST<HAPLO_DUMMY_LOCUS> haplotype_locus;		//Dummy locus information in haplotype extraction
extern LIST<HAPLO_DUMMY_LOCUS> haplotype_slocus;	//Dummy locus information in haplotype extraction
extern SLOCUS *haplotype_nslocus;					//Extracted locus
extern LOCUS *haplotype_nlocus;						//Extracted locus

extern FILE *convert_file;							//Convert file pointer
extern char **convert_buf;							//Circle buffer for file conversion, NBUF
extern int64 convert_linesize;						//Max length of each line in converted file

extern DIVERSITY* diversity_buf;					//Circle buffer for diversity estimation, NBUF
extern DIVSUM diversity_sum;						//Diversity sum
extern int diversity_stage;							//Diversity level, 3 total, 2 pop, 1 reg

extern FST *fst_buf[6];								//Read/Write buffer for each type of fst
extern int fst_type;								//1 Among regions, 2 among pops, 3 among pops/regs in region, 4 between regions, 5 between pops 

extern GDIST *gdist_buf;							//Circle buffer for genetic distance estimation, NBUF
extern int gdist_type;								//1 between inds, 2 between pops, 3 + between regions

extern AMOVA* amova_buf;							//Read/Write buffer, nthread

extern RELATEDNESS* relatedness_buf;				//Circle buffer for relatedness estimation, NBUF

extern KINSHIP* kinship_buf;						//Circle buffer for kinship estimation, NBUF

extern int gdindex[N_GD_ESTIMATOR + 1];				//Index of ith used estimator

extern double *pcoa_matrix;							//Genetic distance array of genetic distance to perform PCoA

extern SRUNINFO *structure_par;						//Genetic distance for pcoa and hierarchical clustering
extern int structure_totalruns;						//Total number of runs in Bayesian clustering

extern double *spa_x;								//Spatial pattern, TEST, don't use, n * (dim + 1)
extern int spa_n;									//Spatial pattern, TEST, don't use, n
extern int spa_np;									//Spatial pattern, TEST, don't use, n



#ifndef _BCFHEADER
	TARGET BCFHEADER::~BCFHEADER()
	{
		if (contig_name == NULL) return;
		for (int64 i = 0; i < contig_size; ++i)
			delete[] contig_name[i];
		delete[] contig_name;
		contig_name = NULL;
	}

	TARGET BCFHEADER::BCFHEADER(char *header)
	{
		int64 line = 0;
		contig_size = 0;
		contig_name = NULL;
		filter_passidx = format_gtid = format_gqid = format_dpid = format_adid = 0xFFFFFFFF;
		char *bak = header;

		int filteridx = 0, fmtidx = 0, contigidx = 0;
		while (*header)
		{
			line++;
			char *hbak = header;
			while (*header == '#') header++;

			if (LwrLineCmp("filter=<", header) == 0)
			{
				header = LineNextIdx(header, "ID=", 1);
				if (header == NULL)
				{
					header = StrNextIdx(header, '\n', 1); if (header) *header = '\0';
					Exit("\nError: bcf header format error, at line %d:\n%s\n", line, hbak);
				}
				header += 3;
				if (LwrLineCmp("pass", header) == 0)
				{
					header = LineNextIdx(header, "IDX=", 1) + 4;
					if (header == (char*)4)
						filter_passidx = ReadInteger(header);
					else
						filter_passidx = filteridx;
				}
				filteridx++;
			}
			else if (LwrLineCmp("contig=<", header) == 0)
			{
				header = LineNextIdx(header, "IDX=", 1) + 4;
				int64 idx = header == (char*)4 ? contigidx : ReadInteger(header);
				contig_size = Max(idx + 1, contig_size);
				contigidx++;
			}
			else if (LwrLineCmp("format=<", header) == 0)
			{
				header = LineNextIdx(header, "ID=", 1);
				if (header == NULL)
				{
					header = StrNextIdx(header, '\n', 1); if (header) *header = '\0';
					Exit("\nError: bcf header format error, at line %d:\n%s\n", line, hbak);
				}
				header += 3;

				if ((LwrLineCmp("gt", header) && LwrLineCmp("dp", header) && LwrLineCmp("gq", header) && LwrLineCmp("ad", header)) == 0)
				{
					char ftype = header[1];
					header = LineNextIdx(header, "IDX=", 1) + 4;
					switch (ftype)
					{
					case 'T':
					case 't': format_gtid = header == (char*)4 ? fmtidx : ReadInteger(header); break;
					case 'Q':
					case 'q': format_gqid = header == (char*)4 ? fmtidx : ReadInteger(header); break;
					case 'P':
					case 'p': format_dpid = header == (char*)4 ? fmtidx : ReadInteger(header); break;
					case 'D':
					case 'd': format_adid = header == (char*)4 ? fmtidx : ReadInteger(header); break;
					}
				}
				fmtidx++;
			}
			header = StrNextIdx(header, '\n', 1) + 1;
		}

		if (format_gtid == -1)
			Exit("\nError: bcf header format error, no GT format field found.", line);

		if (contig_size == 0)
			Exit("\nError: bcf header format error, no contig/chromosome found.", line);

		contig_name = new char*[contig_size];

		line = 0; contigidx = 0; header = bak;
		while (*header)
		{
			line++;
			while (*header == '#') header++;

			if (LwrLineCmp("contig=<", header) == 0)
			{
				header = LineNextIdx(header, "ID=", 1) + 3;
				int64 tlen = StrNextIdx(header, ',', 1) - header;
				char *tstr = new char[tlen + 1];
				memcpy(tstr, header, tlen);
				tstr[tlen] = '\0';

				header = LineNextIdx(header, "IDX=", 1) + 4;
				int64 idx = 0;
				if (header == (char*)4)
					idx = contigidx;
				else
					idx = ReadInteger(header);
				contig_name[idx] = tstr;
				contigidx++;
			}
			header = StrNextIdx(header, '\n', 1) + 1;
		}
	}
#endif

#ifndef _GENOTYPE
	/* Do nothing */
	TARGET GENOTYPE::GENOTYPE()
	{

	}

	/* Create genotype from alleles and hash */
	TARGET GENOTYPE::GENOTYPE(ushort *&gatab, ushort *alleles, int ploidy)
	{
		//memory == NULL for reconstruct: index alleles for non-vcf input
		//alleles should not be sorted for phased genoytpe
		
		patternid = 0;
		
		//zero ploidy, not set
		if (ploidy == 0)
		{
			SetAlleleArray(0);
			return;
		}

		if (alleles[0] == 0xFFFF)
		{
			//missing but have ploidy info
			patternid = N_PATTERN_END + ploidy;
			SetAlleleArray((ushort*)-1);
			return;
		}

		byte acount[N_MAX_PLOIDY] = { 0 };
		ushort akeys[N_MAX_PLOIDY] = { 0 };
		int nalleles = 0;

		//counting number of alleles
		for (int i = 0; i < ploidy; ++i)
		{
			bool isnew = true;
			ushort asi = alleles[i];
			for (int j = i - 1; j >= 0; --j)
				if (asi == alleles[j])
				{
					isnew = false;
					break;
				}
			if (!isnew) continue;

			//asi is a new allele, count dosage
			nalleles++;

			int ac = 1;
			akeys[nalleles - 1] = asi;
			for (int j = i + 1; j < ploidy; ++j)
				if (asi == alleles[j]) ac++;
			acount[nalleles - 1] = ac;
		}

		//sorting alleles in descending order of dosage
		for (int i = 0; i < nalleles; ++i)
			for (int j = i + 1; j < nalleles; ++j)
				if (acount[i] < acount[j] ||
					acount[i] == acount[j] && akeys[i] < akeys[j])
				{
					Swap(acount[i], acount[j]);
					Swap(akeys[i], akeys[j]);
				}

		ushort *als;
		if (gatab)
		{
			als = gatab;
			gatab += ploidy + nalleles;
			SetAlleleArray(als);
		}
		else //reconstruct, replace in indexing alleles, do not allloc new memory
			als = GetAlleleArray();

		SetVal(als, alleles, ploidy);
		SetVal(als + ploidy, akeys, nalleles);

		//allele pattern
		uint64 pattern = 0;
		for (int i = 0; i < nalleles; ++i)
			pattern = (pattern << 4) | acount[i];

		switch (pattern)
		{
		default: patternid = 0; //missing but have ploidy info
		case 0x1ull: patternid = 1; break;
		case 0x2ull: patternid = 2; break;
		case 0x11ull: patternid = 3; break;
		case 0x3ull: patternid = 4; break;
		case 0x21ull: patternid = 5; break;
		case 0x111ull: patternid = 6; break;
		case 0x4ull: patternid = 7; break;
		case 0x31ull: patternid = 8; break;
		case 0x22ull: patternid = 9; break;
		case 0x211ull: patternid = 10; break;
		case 0x1111ull: patternid = 11; break;
		case 0x5ull: patternid = 12; break;
		case 0x41ull: patternid = 13; break;
		case 0x32ull: patternid = 14; break;
		case 0x311ull: patternid = 15; break;
		case 0x221ull: patternid = 16; break;
		case 0x2111ull: patternid = 17; break;
		case 0x11111ull: patternid = 18; break;
		case 0x6ull: patternid = 19; break;
		case 0x51ull: patternid = 20; break;
		case 0x42ull: patternid = 21; break;
		case 0x411ull: patternid = 22; break;
		case 0x33ull: patternid = 23; break;
		case 0x321ull: patternid = 24; break;
		case 0x3111ull: patternid = 25; break;
		case 0x222ull: patternid = 26; break;
		case 0x2211ull: patternid = 27; break;
		case 0x21111ull: patternid = 28; break;
		case 0x111111ull: patternid = 29; break;
		case 0x7ull: patternid = 30; break;
		case 0x61ull: patternid = 31; break;
		case 0x52ull: patternid = 32; break;
		case 0x511ull: patternid = 33; break;
		case 0x43ull: patternid = 34; break;
		case 0x421ull: patternid = 35; break;
		case 0x4111ull: patternid = 36; break;
		case 0x331ull: patternid = 37; break;
		case 0x322ull: patternid = 38; break;
		case 0x3211ull: patternid = 39; break;
		case 0x31111ull: patternid = 40; break;
		case 0x2221ull: patternid = 41; break;
		case 0x22111ull: patternid = 42; break;
		case 0x211111ull: patternid = 43; break;
		case 0x1111111ull: patternid = 44; break;
		case 0x8ull: patternid = 45; break;
		case 0x71ull: patternid = 46; break;
		case 0x62ull: patternid = 47; break;
		case 0x611ull: patternid = 48; break;
		case 0x53ull: patternid = 49; break;
		case 0x521ull: patternid = 50; break;
		case 0x5111ull: patternid = 51; break;
		case 0x44ull: patternid = 52; break;
		case 0x431ull: patternid = 53; break;
		case 0x422ull: patternid = 54; break;
		case 0x4211ull: patternid = 55; break;
		case 0x41111ull: patternid = 56; break;
		case 0x332ull: patternid = 57; break;
		case 0x3311ull: patternid = 58; break;
		case 0x3221ull: patternid = 59; break;
		case 0x32111ull: patternid = 60; break;
		case 0x311111ull: patternid = 61; break;
		case 0x2222ull: patternid = 62; break;
		case 0x22211ull: patternid = 63; break;
		case 0x221111ull: patternid = 64; break;
		case 0x2111111ull: patternid = 65; break;
		case 0x11111111ull: patternid = 66; break;
		case 0x9ull: patternid = 67; break;
		case 0x81ull: patternid = 68; break;
		case 0x72ull: patternid = 69; break;
		case 0x711ull: patternid = 70; break;
		case 0x63ull: patternid = 71; break;
		case 0x621ull: patternid = 72; break;
		case 0x6111ull: patternid = 73; break;
		case 0x54ull: patternid = 74; break;
		case 0x531ull: patternid = 75; break;
		case 0x522ull: patternid = 76; break;
		case 0x5211ull: patternid = 77; break;
		case 0x51111ull: patternid = 78; break;
		case 0x441ull: patternid = 79; break;
		case 0x432ull: patternid = 80; break;
		case 0x4311ull: patternid = 81; break;
		case 0x4221ull: patternid = 82; break;
		case 0x42111ull: patternid = 83; break;
		case 0x411111ull: patternid = 84; break;
		case 0x333ull: patternid = 85; break;
		case 0x3321ull: patternid = 86; break;
		case 0x33111ull: patternid = 87; break;
		case 0x3222ull: patternid = 88; break;
		case 0x32211ull: patternid = 89; break;
		case 0x321111ull: patternid = 90; break;
		case 0x3111111ull: patternid = 91; break;
		case 0x22221ull: patternid = 92; break;
		case 0x222111ull: patternid = 93; break;
		case 0x2211111ull: patternid = 94; break;
		case 0x21111111ull: patternid = 95; break;
		case 0x111111111ull: patternid = 96; break;
		case 0xAull: patternid = 97; break;
		case 0x91ull: patternid = 98; break;
		case 0x82ull: patternid = 99; break;
		case 0x811ull: patternid = 100; break;
		case 0x73ull: patternid = 101; break;
		case 0x721ull: patternid = 102; break;
		case 0x7111ull: patternid = 103; break;
		case 0x64ull: patternid = 104; break;
		case 0x631ull: patternid = 105; break;
		case 0x622ull: patternid = 106; break;
		case 0x6211ull: patternid = 107; break;
		case 0x61111ull: patternid = 108; break;
		case 0x55ull: patternid = 109; break;
		case 0x541ull: patternid = 110; break;
		case 0x532ull: patternid = 111; break;
		case 0x5311ull: patternid = 112; break;
		case 0x5221ull: patternid = 113; break;
		case 0x52111ull: patternid = 114; break;
		case 0x511111ull: patternid = 115; break;
		case 0x442ull: patternid = 116; break;
		case 0x4411ull: patternid = 117; break;
		case 0x433ull: patternid = 118; break;
		case 0x4321ull: patternid = 119; break;
		case 0x43111ull: patternid = 120; break;
		case 0x4222ull: patternid = 121; break;
		case 0x42211ull: patternid = 122; break;
		case 0x421111ull: patternid = 123; break;
		case 0x4111111ull: patternid = 124; break;
		case 0x3331ull: patternid = 125; break;
		case 0x3322ull: patternid = 126; break;
		case 0x33211ull: patternid = 127; break;
		case 0x331111ull: patternid = 128; break;
		case 0x32221ull: patternid = 129; break;
		case 0x322111ull: patternid = 130; break;
		case 0x3211111ull: patternid = 131; break;
		case 0x31111111ull: patternid = 132; break;
		case 0x22222ull: patternid = 133; break;
		case 0x222211ull: patternid = 134; break;
		case 0x2221111ull: patternid = 135; break;
		case 0x22111111ull: patternid = 136; break;
		case 0x211111111ull: patternid = 137; break;
		case 0x1111111111ull: patternid = 138; break;
		}
	}

	/* Copy from a reference genotype */
	TARGET GENOTYPE::GENOTYPE(ushort *&gatab, GENOTYPE &ref)
	{
		patternid = ref.patternid;
		int nalleles = Nalleles();

		if (nalleles)
		{
			int ploidy = Ploidy();
			ushort *als = gatab;
			gatab += ploidy + nalleles;
			SetVal(als, ref.GetAlleleArray(), ploidy + nalleles);
			SetAlleleArray(als);
		}
		else
			SetAlleleArray((ushort*)-1);
	}

	/* Get allele copy at ith haplotype */
	TARGET ushort GENOTYPE::GetAlleleCopy(int i)
	{
		switch (offset)
		{
		case 0:
			return NULL;
		case 0xFFFFFF:
			return 0xFFFF;
		default:
			return *((ushort*)((byte*)this + offset) + i);
		}
	}

	/* Get allele array */
	TARGET ushort *GENOTYPE::GetAlleleArray()
	{
		switch (offset)
		{
		case 0:
			return NULL;
		case 0xFFFFFF:
			return missing_array;
		default:
			return (ushort*)((byte*)this + offset);
		}
	}

	/* Set allele array */
	TARGET void GENOTYPE::SetAlleleArray(ushort *alleles)
	{
		if (alleles == NULL) //not set
			offset = 0;
		else if (alleles == (ushort*)-1)//empty
			offset = 0xFFFFFF;
		else
			offset = (byte*)alleles - (byte*)this;
	}

	/* Get pattern code */
	TARGET uint64 GENOTYPE::GetPattern()
	{
		switch (patternid)
		{
		default: 	return 0ull;
		case	1: 	return 0x1ull;
		case	2: 	return 0x2ull;
		case	3: 	return 0x11ull;
		case	4: 	return 0x3ull;
		case	5: 	return 0x21ull;
		case	6: 	return 0x111ull;
		case	7: 	return 0x4ull;
		case	8: 	return 0x31ull;
		case	9: 	return 0x22ull;
		case	10: 	return 0x211ull;
		case	11: 	return 0x1111ull;
		case	12: 	return 0x5ull;
		case	13: 	return 0x41ull;
		case	14: 	return 0x32ull;
		case	15: 	return 0x311ull;
		case	16: 	return 0x221ull;
		case	17: 	return 0x2111ull;
		case	18: 	return 0x11111ull;
		case	19: 	return 0x6ull;
		case	20: 	return 0x51ull;
		case	21: 	return 0x42ull;
		case	22: 	return 0x411ull;
		case	23: 	return 0x33ull;
		case	24: 	return 0x321ull;
		case	25: 	return 0x3111ull;
		case	26: 	return 0x222ull;
		case	27: 	return 0x2211ull;
		case	28: 	return 0x21111ull;
		case	29: 	return 0x111111ull;
		case	30: 	return 0x7ull;
		case	31: 	return 0x61ull;
		case	32: 	return 0x52ull;
		case	33: 	return 0x511ull;
		case	34: 	return 0x43ull;
		case	35: 	return 0x421ull;
		case	36: 	return 0x4111ull;
		case	37: 	return 0x331ull;
		case	38: 	return 0x322ull;
		case	39: 	return 0x3211ull;
		case	40: 	return 0x31111ull;
		case	41: 	return 0x2221ull;
		case	42: 	return 0x22111ull;
		case	43: 	return 0x211111ull;
		case	44: 	return 0x1111111ull;
		case	45: 	return 0x8ull;
		case	46: 	return 0x71ull;
		case	47: 	return 0x62ull;
		case	48: 	return 0x611ull;
		case	49: 	return 0x53ull;
		case	50: 	return 0x521ull;
		case	51: 	return 0x5111ull;
		case	52: 	return 0x44ull;
		case	53: 	return 0x431ull;
		case	54: 	return 0x422ull;
		case	55: 	return 0x4211ull;
		case	56: 	return 0x41111ull;
		case	57: 	return 0x332ull;
		case	58: 	return 0x3311ull;
		case	59: 	return 0x3221ull;
		case	60: 	return 0x32111ull;
		case	61: 	return 0x311111ull;
		case	62: 	return 0x2222ull;
		case	63: 	return 0x22211ull;
		case	64: 	return 0x221111ull;
		case	65: 	return 0x2111111ull;
		case	66: 	return 0x11111111ull;
		case	67: 	return 0x9ull;
		case	68: 	return 0x81ull;
		case	69: 	return 0x72ull;
		case	70: 	return 0x711ull;
		case	71: 	return 0x63ull;
		case	72: 	return 0x621ull;
		case	73: 	return 0x6111ull;
		case	74: 	return 0x54ull;
		case	75: 	return 0x531ull;
		case	76: 	return 0x522ull;
		case	77: 	return 0x5211ull;
		case	78: 	return 0x51111ull;
		case	79: 	return 0x441ull;
		case	80: 	return 0x432ull;
		case	81: 	return 0x4311ull;
		case	82: 	return 0x4221ull;
		case	83: 	return 0x42111ull;
		case	84: 	return 0x411111ull;
		case	85: 	return 0x333ull;
		case	86: 	return 0x3321ull;
		case	87: 	return 0x33111ull;
		case	88: 	return 0x3222ull;
		case	89: 	return 0x32211ull;
		case	90: 	return 0x321111ull;
		case	91: 	return 0x3111111ull;
		case	92: 	return 0x22221ull;
		case	93: 	return 0x222111ull;
		case	94: 	return 0x2211111ull;
		case	95: 	return 0x21111111ull;
		case	96: 	return 0x111111111ull;
		case	97: 	return 0xAull;
		case	98: 	return 0x91ull;
		case	99: 	return 0x82ull;
		case	100: 	return 0x811ull;
		case	101: 	return 0x73ull;
		case	102: 	return 0x721ull;
		case	103: 	return 0x7111ull;
		case	104: 	return 0x64ull;
		case	105: 	return 0x631ull;
		case	106: 	return 0x622ull;
		case	107: 	return 0x6211ull;
		case	108: 	return 0x61111ull;
		case	109: 	return 0x55ull;
		case	110: 	return 0x541ull;
		case	111: 	return 0x532ull;
		case	112: 	return 0x5311ull;
		case	113: 	return 0x5221ull;
		case	114: 	return 0x52111ull;
		case	115: 	return 0x511111ull;
		case	116: 	return 0x442ull;
		case	117: 	return 0x4411ull;
		case	118: 	return 0x433ull;
		case	119: 	return 0x4321ull;
		case	120: 	return 0x43111ull;
		case	121: 	return 0x4222ull;
		case	122: 	return 0x42211ull;
		case	123: 	return 0x421111ull;
		case	124: 	return 0x4111111ull;
		case	125: 	return 0x3331ull;
		case	126: 	return 0x3322ull;
		case	127: 	return 0x33211ull;
		case	128: 	return 0x331111ull;
		case	129: 	return 0x32221ull;
		case	130: 	return 0x322111ull;
		case	131: 	return 0x3211111ull;
		case	132: 	return 0x31111111ull;
		case	133: 	return 0x22222ull;
		case	134: 	return 0x222211ull;
		case	135: 	return 0x2221111ull;
		case	136: 	return 0x22111111ull;
		case	137: 	return 0x211111111ull;
		case	138: 	return 0x1111111111ull;
		}
	}

	/* Number of alleles */
	TARGET int GENOTYPE::Nalleles()
	{
		/*
		switch (pattern)
		{
			default:	return 0;
			case 0x1ull:	return 1;
			case 0x2ull:	return 1;
			case 0x11ull:	return 2;
			case 0x3ull:	return 1;
			case 0x21ull:	return 2;
			case 0x111ull:	return 3;
			case 0x4ull:	return 1;
			case 0x31ull:	return 2;
			case 0x22ull:	return 2;
			case 0x211ull:	return 3;
			case 0x1111ull:	return 4;
			case 0x5ull:	return 1;
			case 0x41ull:	return 2;
			case 0x32ull:	return 2;
			case 0x311ull:	return 3;
			case 0x221ull:	return 3;
			case 0x2111ull:	return 4;
			case 0x11111ull:	return 5;
			case 0x6ull:	return 1;
			case 0x51ull:	return 2;
			case 0x42ull:	return 2;
			case 0x411ull:	return 3;
			case 0x33ull:	return 2;
			case 0x321ull:	return 3;
			case 0x3111ull:	return 4;
			case 0x222ull:	return 3;
			case 0x2211ull:	return 4;
			case 0x21111ull:	return 5;
			case 0x111111ull:	return 6;
			case 0x7ull:	return 1;
			case 0x61ull:	return 2;
			case 0x52ull:	return 2;
			case 0x511ull:	return 3;
			case 0x43ull:	return 2;
			case 0x421ull:	return 3;
			case 0x4111ull:	return 4;
			case 0x331ull:	return 3;
			case 0x322ull:	return 3;
			case 0x3211ull:	return 4;
			case 0x31111ull:	return 5;
			case 0x2221ull:	return 4;
			case 0x22111ull:	return 5;
			case 0x211111ull:	return 6;
			case 0x1111111ull:	return 7;
			case 0x8ull:	return 1;
			case 0x71ull:	return 2;
			case 0x62ull:	return 2;
			case 0x611ull:	return 3;
			case 0x53ull:	return 2;
			case 0x521ull:	return 3;
			case 0x5111ull:	return 4;
			case 0x44ull:	return 2;
			case 0x431ull:	return 3;
			case 0x422ull:	return 3;
			case 0x4211ull:	return 4;
			case 0x41111ull:	return 5;
			case 0x332ull:	return 3;
			case 0x3311ull:	return 4;
			case 0x3221ull:	return 4;
			case 0x32111ull:	return 5;
			case 0x311111ull:	return 6;
			case 0x2222ull:	return 4;
			case 0x22211ull:	return 5;
			case 0x221111ull:	return 6;
			case 0x2111111ull:	return 7;
			case 0x11111111ull:	return 8;
			case 0x9ull:	return 1;
			case 0x81ull:	return 2;
			case 0x72ull:	return 2;
			case 0x711ull:	return 3;
			case 0x63ull:	return 2;
			case 0x621ull:	return 3;
			case 0x6111ull:	return 4;
			case 0x54ull:	return 2;
			case 0x531ull:	return 3;
			case 0x522ull:	return 3;
			case 0x5211ull:	return 4;
			case 0x51111ull:	return 5;
			case 0x441ull:	return 3;
			case 0x432ull:	return 3;
			case 0x4311ull:	return 4;
			case 0x4221ull:	return 4;
			case 0x42111ull:	return 5;
			case 0x411111ull:	return 6;
			case 0x333ull:	return 3;
			case 0x3321ull:	return 4;
			case 0x33111ull:	return 5;
			case 0x3222ull:	return 4;
			case 0x32211ull:	return 5;
			case 0x321111ull:	return 6;
			case 0x3111111ull:	return 7;
			case 0x22221ull:	return 5;
			case 0x222111ull:	return 6;
			case 0x2211111ull:	return 7;
			case 0x21111111ull:	return 8;
			case 0x111111111ull:	return 9;
			case 0xAull:	return 1;
			case 0x91ull:	return 2;
			case 0x82ull:	return 2;
			case 0x811ull:	return 3;
			case 0x73ull:	return 2;
			case 0x721ull:	return 3;
			case 0x7111ull:	return 4;
			case 0x64ull:	return 2;
			case 0x631ull:	return 3;
			case 0x622ull:	return 3;
			case 0x6211ull:	return 4;
			case 0x61111ull:	return 5;
			case 0x55ull:	return 2;
			case 0x541ull:	return 3;
			case 0x532ull:	return 3;
			case 0x5311ull:	return 4;
			case 0x5221ull:	return 4;
			case 0x52111ull:	return 5;
			case 0x511111ull:	return 6;
			case 0x442ull:	return 3;
			case 0x4411ull:	return 4;
			case 0x433ull:	return 3;
			case 0x4321ull:	return 4;
			case 0x43111ull:	return 5;
			case 0x4222ull:	return 4;
			case 0x42211ull:	return 5;
			case 0x421111ull:	return 6;
			case 0x4111111ull:	return 7;
			case 0x3331ull:	return 4;
			case 0x3322ull:	return 4;
			case 0x33211ull:	return 5;
			case 0x331111ull:	return 6;
			case 0x32221ull:	return 5;
			case 0x322111ull:	return 6;
			case 0x3211111ull:	return 7;
			case 0x31111111ull:	return 8;
			case 0x22222ull:	return 5;
			case 0x222211ull:	return 6;
			case 0x2221111ull:	return 7;
			case 0x22111111ull:	return 8;
			case 0x211111111ull:	return 9;
			case 0x1111111111ull:	return 10;
		}
		*/
		switch (patternid)
		{
			default: 	return 0;
			case	1: 	return 1;
			case	2: 	return 1;
			case	3: 	return 2;
			case	4: 	return 1;
			case	5: 	return 2;
			case	6: 	return 3;
			case	7: 	return 1;
			case	8: 	return 2;
			case	9: 	return 2;
			case	10: 	return 3;
			case	11: 	return 4;
			case	12: 	return 1;
			case	13: 	return 2;
			case	14: 	return 2;
			case	15: 	return 3;
			case	16: 	return 3;
			case	17: 	return 4;
			case	18: 	return 5;
			case	19: 	return 1;
			case	20: 	return 2;
			case	21: 	return 2;
			case	22: 	return 3;
			case	23: 	return 2;
			case	24: 	return 3;
			case	25: 	return 4;
			case	26: 	return 3;
			case	27: 	return 4;
			case	28: 	return 5;
			case	29: 	return 6;
			case	30: 	return 1;
			case	31: 	return 2;
			case	32: 	return 2;
			case	33: 	return 3;
			case	34: 	return 2;
			case	35: 	return 3;
			case	36: 	return 4;
			case	37: 	return 3;
			case	38: 	return 3;
			case	39: 	return 4;
			case	40: 	return 5;
			case	41: 	return 4;
			case	42: 	return 5;
			case	43: 	return 6;
			case	44: 	return 7;
			case	45: 	return 1;
			case	46: 	return 2;
			case	47: 	return 2;
			case	48: 	return 3;
			case	49: 	return 2;
			case	50: 	return 3;
			case	51: 	return 4;
			case	52: 	return 2;
			case	53: 	return 3;
			case	54: 	return 3;
			case	55: 	return 4;
			case	56: 	return 5;
			case	57: 	return 3;
			case	58: 	return 4;
			case	59: 	return 4;
			case	60: 	return 5;
			case	61: 	return 6;
			case	62: 	return 4;
			case	63: 	return 5;
			case	64: 	return 6;
			case	65: 	return 7;
			case	66: 	return 8;
			case	67: 	return 1;
			case	68: 	return 2;
			case	69: 	return 2;
			case	70: 	return 3;
			case	71: 	return 2;
			case	72: 	return 3;
			case	73: 	return 4;
			case	74: 	return 2;
			case	75: 	return 3;
			case	76: 	return 3;
			case	77: 	return 4;
			case	78: 	return 5;
			case	79: 	return 3;
			case	80: 	return 3;
			case	81: 	return 4;
			case	82: 	return 4;
			case	83: 	return 5;
			case	84: 	return 6;
			case	85: 	return 3;
			case	86: 	return 4;
			case	87: 	return 5;
			case	88: 	return 4;
			case	89: 	return 5;
			case	90: 	return 6;
			case	91: 	return 7;
			case	92: 	return 5;
			case	93: 	return 6;
			case	94: 	return 7;
			case	95: 	return 8;
			case	96: 	return 9;
			case	97: 	return 1;
			case	98: 	return 2;
			case	99: 	return 2;
			case	100: 	return 3;
			case	101: 	return 2;
			case	102: 	return 3;
			case	103: 	return 4;
			case	104: 	return 2;
			case	105: 	return 3;
			case	106: 	return 3;
			case	107: 	return 4;
			case	108: 	return 5;
			case	109: 	return 2;
			case	110: 	return 3;
			case	111: 	return 3;
			case	112: 	return 4;
			case	113: 	return 4;
			case	114: 	return 5;
			case	115: 	return 6;
			case	116: 	return 3;
			case	117: 	return 4;
			case	118: 	return 3;
			case	119: 	return 4;
			case	120: 	return 5;
			case	121: 	return 4;
			case	122: 	return 5;
			case	123: 	return 6;
			case	124: 	return 7;
			case	125: 	return 4;
			case	126: 	return 4;
			case	127: 	return 5;
			case	128: 	return 6;
			case	129: 	return 5;
			case	130: 	return 6;
			case	131: 	return 7;
			case	132: 	return 8;
			case	133: 	return 5;
			case	134: 	return 6;
			case	135: 	return 7;
			case	136: 	return 8;
			case	137: 	return 9;
			case	138: 	return 10;
		}
	}

	/* Ploidy level */
	TARGET int GENOTYPE::Ploidy()
	{
		/*
		switch (pattern)
		{
		default: return 0;
		case 0x1ull: return 1;
		case 0x2ull: 
		case 0x11ull: return 2;
		case 0x3ull: 
		case 0x21ull: 
		case 0x111ull: return 3;
		case 0x4ull: 
		case 0x31ull: 
		case 0x22ull: 
		case 0x211ull: 
		case 0x1111ull: return 4;
		case 0x5ull: 
		case 0x41ull: 
		case 0x32ull: 
		case 0x311ull: 
		case 0x221ull: 
		case 0x2111ull: 
		case 0x11111ull: return 5;
		case 0x6ull: 
		case 0x51ull: 
		case 0x42ull: 
		case 0x411ull: 
		case 0x33ull: 
		case 0x321ull: 
		case 0x3111ull: 
		case 0x222ull:
		case 0x2211ull: 
		case 0x21111ull: 
		case 0x111111ull: return 6;
		case 0x7ull: 
		case 0x61ull: 
		case 0x52ull: 
		case 0x511ull: 
		case 0x43ull: 
		case 0x421ull: 
		case 0x4111ull: 
		case 0x331ull: 
		case 0x322ull:
		case 0x3211ull: 
		case 0x31111ull: 
		case 0x2221ull: 
		case 0x22111ull: 
		case 0x211111ull: 
		case 0x1111111ull: return 7;
		case 0x8ull: 
		case 0x71ull: 
		case 0x62ull: 
		case 0x611ull:
		case 0x53ull: 
		case 0x521ull: 
		case 0x5111ull: 
		case 0x44ull: 
		case 0x431ull: 
		case 0x422ull: 
		case 0x4211ull: 
		case 0x41111ull: 
		case 0x332ull: 
		case 0x3311ull: 
		case 0x3221ull: 
		case 0x32111ull: 
		case 0x311111ull: 
		case 0x2222ull: 
		case 0x22211ull: 
		case 0x221111ull: 
		case 0x2111111ull: 
		case 0x11111111ull: return 8;
		case 0x9ull: 
		case 0x81ull: 
		case 0x72ull: 
		case 0x711ull: 
		case 0x63ull: 
		case 0x621ull: 
		case 0x6111ull: 
		case 0x54ull: 
		case 0x531ull: 
		case 0x522ull: 
		case 0x5211ull: 
		case 0x51111ull: 
		case 0x441ull: 
		case 0x432ull: 
		case 0x4311ull: 
		case 0x4221ull: 
		case 0x42111ull: 
		case 0x411111ull: 
		case 0x333ull: 
		case 0x3321ull: 
		case 0x33111ull: 
		case 0x3222ull: 
		case 0x32211ull: 
		case 0x321111ull: 
		case 0x3111111ull: 
		case 0x22221ull: 
		case 0x222111ull: 
		case 0x2211111ull: 
		case 0x21111111ull: 
		case 0x111111111ull: return 9;
		case 0xAull: 
		case 0x91ull: 
		case 0x82ull: 
		case 0x811ull: 
		case 0x73ull: 
		case 0x721ull: 
		case 0x7111ull: 
		case 0x64ull: 
		case 0x631ull: 
		case 0x622ull: 
		case 0x6211ull: 
		case 0x61111ull: 
		case 0x55ull: 
		case 0x541ull: 
		case 0x532ull: 
		case 0x5311ull: 
		case 0x5221ull: 
		case 0x52111ull: 
		case 0x511111ull: 
		case 0x442ull: 
		case 0x4411ull: 
		case 0x433ull: 
		case 0x4321ull: 
		case 0x43111ull: 
		case 0x4222ull: 
		case 0x42211ull: 
		case 0x421111ull: 
		case 0x4111111ull: 
		case 0x3331ull: 
		case 0x3322ull: 
		case 0x33211ull: 
		case 0x331111ull: 
		case 0x32221ull: 
		case 0x322111ull: 
		case 0x3211111ull: 
		case 0x31111111ull:
		case 0x22222ull: 
		case 0x222211ull: 
		case 0x2221111ull: 
		case 0x22111111ull: 
		case 0x211111111ull: 
		case 0x1111111111ull: return 10;
		}
		*/
		switch (patternid)
		{
		default: return patternid >= N_PATTERN_END ? patternid - N_PATTERN_END : 0;
			case 1: return 1;
			case 2:
			case 3: return 2;
			case 4:
			case 5:
			case 6: return 3;
			case 7:
			case 8:
			case 9:
			case 10:
			case 11: return 4;
			case 12:
			case 13:
			case 14:
			case 15:
			case 16:
			case 17:
			case 18: return 5;
			case 19:
			case 20:
			case 21:
			case 22:
			case 23:
			case 24:
			case 25:
			case 26:
			case 27:
			case 28:
			case 29: return 6;
			case 30:
			case 31:
			case 32:
			case 33:
			case 34:
			case 35:
			case 36:
			case 37:
			case 38:
			case 39:
			case 40:
			case 41:
			case 42:
			case 43:
			case 44: return 7;
			case 45:
			case 46:
			case 47:
			case 48:
			case 49:
			case 50:
			case 51:
			case 52:
			case 53:
			case 54:
			case 55:
			case 56:
			case 57:
			case 58:
			case 59:
			case 60:
			case 61:
			case 62:
			case 63:
			case 64:
			case 65:
			case 66: return 8;
			case 67:
			case 68:
			case 69:
			case 70:
			case 71:
			case 72:
			case 73:
			case 74:
			case 75:
			case 76:
			case 77:
			case 78:
			case 79:
			case 80:
			case 81:
			case 82:
			case 83:
			case 84:
			case 85:
			case 86:
			case 87:
			case 88:
			case 89:
			case 90:
			case 91:
			case 92:
			case 93:
			case 94:
			case 95:
			case 96: return 9;
			case 97:
			case 98:
			case 99:
			case 100:
			case 101:
			case 102:
			case 103:
			case 104:
			case 105:
			case 106:
			case 107:
			case 108:
			case 109:
			case 110:
			case 111:
			case 112:
			case 113:
			case 114:
			case 115:
			case 116:
			case 117:
			case 118:
			case 119:
			case 120:
			case 121:
			case 122:
			case 123:
			case 124:
			case 125:
			case 126:
			case 127:
			case 128:
			case 129:
			case 130:
			case 131:
			case 132:
			case 133:
			case 134:
			case 135:
			case 136:
			case 137:
			case 138: return 10;
		}
	}

	/* Crc32 hash */
	TARGET HASH GENOTYPE::Hash()
	{
		return HashGenotype(GetAlleleArray(), Ploidy());
	}

	/* Heterozygosity in this genotype */
	TARGET double GENOTYPE::HIndex()
	{
		/*
		switch (pattern)
		{
			default: return NA;
			case 0x1ull: return 0.000000000000000;
			case 0x2ull: return 0.000000000000000;
			case 0x11ull: return 1.000000000000000;
			case 0x3ull: return 0.000000000000000;
			case 0x21ull: return 0.666666666666667;
			case 0x111ull: return 1.000000000000000;
			case 0x4ull: return 0.000000000000000;
			case 0x31ull: return 0.500000000000000;
			case 0x22ull: return 0.666666666666667;
			case 0x211ull: return 0.833333333333333;
			case 0x1111ull: return 1.000000000000000;
			case 0x5ull: return 0.000000000000000;
			case 0x41ull: return 0.400000000000000;
			case 0x32ull: return 0.600000000000000;
			case 0x311ull: return 0.700000000000000;
			case 0x221ull: return 0.800000000000000;
			case 0x2111ull: return 0.900000000000000;
			case 0x11111ull: return 1.000000000000000;
			case 0x6ull: return 0.000000000000000;
			case 0x51ull: return 0.333333333333333;
			case 0x42ull: return 0.533333333333333;
			case 0x411ull: return 0.600000000000000;
			case 0x33ull: return 0.600000000000000;
			case 0x321ull: return 0.733333333333333;
			case 0x3111ull: return 0.800000000000000;
			case 0x222ull: return 0.800000000000000;
			case 0x2211ull: return 0.866666666666667;
			case 0x21111ull: return 0.933333333333333;
			case 0x111111ull: return 1.000000000000000;
			case 0x7ull: return 0.000000000000000;
			case 0x61ull: return 0.285714285714286;
			case 0x52ull: return 0.476190476190476;
			case 0x511ull: return 0.523809523809524;
			case 0x43ull: return 0.571428571428572;
			case 0x421ull: return 0.666666666666667;
			case 0x4111ull: return 0.714285714285714;
			case 0x331ull: return 0.714285714285714;
			case 0x322ull: return 0.761904761904762;
			case 0x3211ull: return 0.809523809523809;
			case 0x31111ull: return 0.857142857142857;
			case 0x2221ull: return 0.857142857142857;
			case 0x22111ull: return 0.904761904761905;
			case 0x211111ull: return 0.952380952380952;
			case 0x1111111ull: return 1.000000000000000;
			case 0x8ull: return 0.000000000000000;
			case 0x71ull: return 0.250000000000000;
			case 0x62ull: return 0.428571428571429;
			case 0x611ull: return 0.464285714285714;
			case 0x53ull: return 0.535714285714286;
			case 0x521ull: return 0.607142857142857;
			case 0x5111ull: return 0.642857142857143;
			case 0x44ull: return 0.571428571428571;
			case 0x431ull: return 0.678571428571429;
			case 0x422ull: return 0.714285714285714;
			case 0x4211ull: return 0.750000000000000;
			case 0x41111ull: return 0.785714285714286;
			case 0x332ull: return 0.750000000000000;
			case 0x3311ull: return 0.785714285714286;
			case 0x3221ull: return 0.821428571428571;
			case 0x32111ull: return 0.857142857142857;
			case 0x311111ull: return 0.892857142857143;
			case 0x2222ull: return 0.857142857142857;
			case 0x22211ull: return 0.892857142857143;
			case 0x221111ull: return 0.928571428571429;
			case 0x2111111ull: return 0.964285714285714;
			case 0x11111111ull: return 1.000000000000000;
			case 0x9ull: return 0.000000000000000;
			case 0x81ull: return 0.222222222222222;
			case 0x72ull: return 0.388888888888889;
			case 0x711ull: return 0.416666666666667;
			case 0x63ull: return 0.500000000000000;
			case 0x621ull: return 0.555555555555556;
			case 0x6111ull: return 0.583333333333333;
			case 0x54ull: return 0.555555555555556;
			case 0x531ull: return 0.638888888888889;
			case 0x522ull: return 0.666666666666667;
			case 0x5211ull: return 0.694444444444444;
			case 0x51111ull: return 0.722222222222222;
			case 0x441ull: return 0.666666666666667;
			case 0x432ull: return 0.722222222222222;
			case 0x4311ull: return 0.750000000000000;
			case 0x4221ull: return 0.777777777777778;
			case 0x42111ull: return 0.805555555555556;
			case 0x411111ull: return 0.833333333333333;
			case 0x333ull: return 0.750000000000000;
			case 0x3321ull: return 0.805555555555556;
			case 0x33111ull: return 0.833333333333333;
			case 0x3222ull: return 0.833333333333333;
			case 0x32211ull: return 0.861111111111111;
			case 0x321111ull: return 0.888888888888889;
			case 0x3111111ull: return 0.916666666666667;
			case 0x22221ull: return 0.888888888888889;
			case 0x222111ull: return 0.916666666666667;
			case 0x2211111ull: return 0.944444444444444;
			case 0x21111111ull: return 0.972222222222222;
			case 0x111111111ull: return 1.000000000000000;
			case 0xAull: return 0.000000000000000;
			case 0x91ull: return 0.200000000000000;
			case 0x82ull: return 0.355555555555555;
			case 0x811ull: return 0.377777777777778;
			case 0x73ull: return 0.466666666666667;
			case 0x721ull: return 0.511111111111111;
			case 0x7111ull: return 0.533333333333333;
			case 0x64ull: return 0.533333333333333;
			case 0x631ull: return 0.600000000000000;
			case 0x622ull: return 0.622222222222222;
			case 0x6211ull: return 0.644444444444444;
			case 0x61111ull: return 0.666666666666667;
			case 0x55ull: return 0.555555555555556;
			case 0x541ull: return 0.644444444444444;
			case 0x532ull: return 0.688888888888889;
			case 0x5311ull: return 0.711111111111111;
			case 0x5221ull: return 0.733333333333333;
			case 0x52111ull: return 0.755555555555555;
			case 0x511111ull: return 0.777777777777778;
			case 0x442ull: return 0.711111111111111;
			case 0x4411ull: return 0.733333333333333;
			case 0x433ull: return 0.733333333333333;
			case 0x4321ull: return 0.777777777777778;
			case 0x43111ull: return 0.800000000000000;
			case 0x4222ull: return 0.800000000000000;
			case 0x42211ull: return 0.822222222222222;
			case 0x421111ull: return 0.844444444444444;
			case 0x4111111ull: return 0.866666666666667;
			case 0x3331ull: return 0.800000000000000;
			case 0x3322ull: return 0.822222222222222;
			case 0x33211ull: return 0.844444444444444;
			case 0x331111ull: return 0.866666666666667;
			case 0x32221ull: return 0.866666666666667;
			case 0x322111ull: return 0.888888888888889;
			case 0x3211111ull: return 0.911111111111111;
			case 0x31111111ull: return 0.933333333333333;
			case 0x22222ull: return 0.888888888888889;
			case 0x222211ull: return 0.911111111111111;
			case 0x2221111ull: return 0.933333333333333;
			case 0x22111111ull: return 0.955555555555555;
			case 0x211111111ull: return 0.977777777777778;
			case 0x1111111111ull: return 1.000000000000000;
		}
		*/
		switch (patternid)
		{
			default: 	return NA;
			case	1: 	return 0.000000000000000;
			case	2: 	return 0.000000000000000;
			case	3: 	return 1.000000000000000;
			case	4: 	return 0.000000000000000;
			case	5: 	return 0.666666666666667;
			case	6: 	return 1.000000000000000;
			case	7: 	return 0.000000000000000;
			case	8: 	return 0.500000000000000;
			case	9: 	return 0.666666666666667;
			case	10: 	return 0.833333333333333;
			case	11: 	return 1.000000000000000;
			case	12: 	return 0.000000000000000;
			case	13: 	return 0.400000000000000;
			case	14: 	return 0.600000000000000;
			case	15: 	return 0.700000000000000;
			case	16: 	return 0.800000000000000;
			case	17: 	return 0.900000000000000;
			case	18: 	return 1.000000000000000;
			case	19: 	return 0.000000000000000;
			case	20: 	return 0.333333333333333;
			case	21: 	return 0.533333333333333;
			case	22: 	return 0.600000000000000;
			case	23: 	return 0.600000000000000;
			case	24: 	return 0.733333333333333;
			case	25: 	return 0.800000000000000;
			case	26: 	return 0.800000000000000;
			case	27: 	return 0.866666666666667;
			case	28: 	return 0.933333333333333;
			case	29: 	return 1.000000000000000;
			case	30: 	return 0.000000000000000;
			case	31: 	return 0.285714285714286;
			case	32: 	return 0.476190476190476;
			case	33: 	return 0.523809523809524;
			case	34: 	return 0.571428571428572;
			case	35: 	return 0.666666666666667;
			case	36: 	return 0.714285714285714;
			case	37: 	return 0.714285714285714;
			case	38: 	return 0.761904761904762;
			case	39: 	return 0.809523809523809;
			case	40: 	return 0.857142857142857;
			case	41: 	return 0.857142857142857;
			case	42: 	return 0.904761904761905;
			case	43: 	return 0.952380952380952;
			case	44: 	return 1.000000000000000;
			case	45: 	return 0.000000000000000;
			case	46: 	return 0.250000000000000;
			case	47: 	return 0.428571428571429;
			case	48: 	return 0.464285714285714;
			case	49: 	return 0.535714285714286;
			case	50: 	return 0.607142857142857;
			case	51: 	return 0.642857142857143;
			case	52: 	return 0.571428571428571;
			case	53: 	return 0.678571428571429;
			case	54: 	return 0.714285714285714;
			case	55: 	return 0.750000000000000;
			case	56: 	return 0.785714285714286;
			case	57: 	return 0.750000000000000;
			case	58: 	return 0.785714285714286;
			case	59: 	return 0.821428571428571;
			case	60: 	return 0.857142857142857;
			case	61: 	return 0.892857142857143;
			case	62: 	return 0.857142857142857;
			case	63: 	return 0.892857142857143;
			case	64: 	return 0.928571428571429;
			case	65: 	return 0.964285714285714;
			case	66: 	return 1.000000000000000;
			case	67: 	return 0.000000000000000;
			case	68: 	return 0.222222222222222;
			case	69: 	return 0.388888888888889;
			case	70: 	return 0.416666666666667;
			case	71: 	return 0.500000000000000;
			case	72: 	return 0.555555555555556;
			case	73: 	return 0.583333333333333;
			case	74: 	return 0.555555555555556;
			case	75: 	return 0.638888888888889;
			case	76: 	return 0.666666666666667;
			case	77: 	return 0.694444444444444;
			case	78: 	return 0.722222222222222;
			case	79: 	return 0.666666666666667;
			case	80: 	return 0.722222222222222;
			case	81: 	return 0.750000000000000;
			case	82: 	return 0.777777777777778;
			case	83: 	return 0.805555555555556;
			case	84: 	return 0.833333333333333;
			case	85: 	return 0.750000000000000;
			case	86: 	return 0.805555555555556;
			case	87: 	return 0.833333333333333;
			case	88: 	return 0.833333333333333;
			case	89: 	return 0.861111111111111;
			case	90: 	return 0.888888888888889;
			case	91: 	return 0.916666666666667;
			case	92: 	return 0.888888888888889;
			case	93: 	return 0.916666666666667;
			case	94: 	return 0.944444444444444;
			case	95: 	return 0.972222222222222;
			case	96: 	return 1.000000000000000;
			case	97: 	return 0.000000000000000;
			case	98: 	return 0.200000000000000;
			case	99: 	return 0.355555555555555;
			case	100: 	return 0.377777777777778;
			case	101: 	return 0.466666666666667;
			case	102: 	return 0.511111111111111;
			case	103: 	return 0.533333333333333;
			case	104: 	return 0.533333333333333;
			case	105: 	return 0.600000000000000;
			case	106: 	return 0.622222222222222;
			case	107: 	return 0.644444444444444;
			case	108: 	return 0.666666666666667;
			case	109: 	return 0.555555555555556;
			case	110: 	return 0.644444444444444;
			case	111: 	return 0.688888888888889;
			case	112: 	return 0.711111111111111;
			case	113: 	return 0.733333333333333;
			case	114: 	return 0.755555555555555;
			case	115: 	return 0.777777777777778;
			case	116: 	return 0.711111111111111;
			case	117: 	return 0.733333333333333;
			case	118: 	return 0.733333333333333;
			case	119: 	return 0.777777777777778;
			case	120: 	return 0.800000000000000;
			case	121: 	return 0.800000000000000;
			case	122: 	return 0.822222222222222;
			case	123: 	return 0.844444444444444;
			case	124: 	return 0.866666666666667;
			case	125: 	return 0.800000000000000;
			case	126: 	return 0.822222222222222;
			case	127: 	return 0.844444444444444;
			case	128: 	return 0.866666666666667;
			case	129: 	return 0.866666666666667;
			case	130: 	return 0.888888888888889;
			case	131: 	return 0.911111111111111;
			case	132: 	return 0.933333333333333;
			case	133: 	return 0.888888888888889;
			case	134: 	return 0.911111111111111;
			case	135: 	return 0.933333333333333;
			case	136: 	return 0.955555555555555;
			case	137: 	return 0.977777777777778;
			case	138: 	return 1.000000000000000;
		}
	}

	/* SS within genotype under IAM model */
	TARGET double GENOTYPE::SS_IAM()
	{
		switch (patternid)
		{
		default: return NA;
		case 2: return 0.0000000000000000;
		case 3: return 0.5000000000000000;
		case 4: return 0.0000000000000000;
		case 5: return 0.6666666666666670;
		case 6: return 1.0000000000000000;
		case 7: return 0.0000000000000000;
		case 8: return 0.7500000000000000;
		case 9: return 1.0000000000000000;
		case 10: return 1.2500000000000000;
		case 11: return 1.5000000000000000;
		case 12: return 0.0000000000000000;
		case 13: return 0.8000000000000000;
		case 14: return 1.2000000000000000;
		case 15: return 1.4000000000000000;
		case 16: return 1.6000000000000000;
		case 17: return 1.8000000000000000;
		case 18: return 2.0000000000000000;
		case 19: return 0.0000000000000000;
		case 20: return 0.8333333333333320;
		case 21: return 1.3333333333333300;
		case 22: return 1.5000000000000000;
		case 23: return 1.5000000000000000;
		case 24: return 1.8333333333333300;
		case 25: return 2.0000000000000000;
		case 26: return 2.0000000000000000;
		case 27: return 2.1666666666666700;
		case 28: return 2.3333333333333300;
		case 29: return 2.5000000000000000;
		case 30: return 0.0000000000000000;
		case 31: return 0.8571428571428580;
		case 32: return 1.4285714285714300;
		case 33: return 1.5714285714285700;
		case 34: return 1.7142857142857200;
		case 35: return 2.0000000000000000;
		case 36: return 2.1428571428571400;
		case 37: return 2.1428571428571400;
		case 38: return 2.2857142857142900;
		case 39: return 2.4285714285714300;
		case 40: return 2.5714285714285700;
		case 41: return 2.5714285714285700;
		case 42: return 2.7142857142857100;
		case 43: return 2.8571428571428600;
		case 44: return 3.0000000000000000;
		case 45: return 0.0000000000000000;
		case 46: return 0.8750000000000000;
		case 47: return 1.5000000000000000;
		case 48: return 1.6250000000000000;
		case 49: return 1.8750000000000000;
		case 50: return 2.1250000000000000;
		case 51: return 2.2500000000000000;
		case 52: return 2.0000000000000000;
		case 53: return 2.3750000000000000;
		case 54: return 2.5000000000000000;
		case 55: return 2.6250000000000000;
		case 56: return 2.7500000000000000;
		case 57: return 2.6250000000000000;
		case 58: return 2.7500000000000000;
		case 59: return 2.8750000000000000;
		case 60: return 3.0000000000000000;
		case 61: return 3.1250000000000000;
		case 62: return 3.0000000000000000;
		case 63: return 3.1250000000000000;
		case 64: return 3.2500000000000000;
		case 65: return 3.3750000000000000;
		case 66: return 3.5000000000000000;
		case 67: return 0.0000000000000000;
		case 68: return 0.8888888888888880;
		case 69: return 1.5555555555555600;
		case 70: return 1.6666666666666700;
		case 71: return 2.0000000000000000;
		case 72: return 2.2222222222222200;
		case 73: return 2.3333333333333300;
		case 74: return 2.2222222222222200;
		case 75: return 2.5555555555555600;
		case 76: return 2.6666666666666700;
		case 77: return 2.7777777777777800;
		case 78: return 2.8888888888888900;
		case 79: return 2.6666666666666700;
		case 80: return 2.8888888888888900;
		case 81: return 3.0000000000000000;
		case 82: return 3.1111111111111100;
		case 83: return 3.2222222222222200;
		case 84: return 3.3333333333333300;
		case 85: return 3.0000000000000000;
		case 86: return 3.2222222222222200;
		case 87: return 3.3333333333333300;
		case 88: return 3.3333333333333300;
		case 89: return 3.4444444444444400;
		case 90: return 3.5555555555555600;
		case 91: return 3.6666666666666700;
		case 92: return 3.5555555555555600;
		case 93: return 3.6666666666666700;
		case 94: return 3.7777777777777800;
		case 95: return 3.8888888888888900;
		case 96: return 4.0000000000000000;
		case 97: return 0.0000000000000000;
		case 98: return 0.9000000000000000;
		case 99: return 1.6000000000000000;
		case 100: return 1.7000000000000000;
		case 101: return 2.1000000000000000;
		case 102: return 2.3000000000000000;
		case 103: return 2.4000000000000000;
		case 104: return 2.4000000000000000;
		case 105: return 2.7000000000000000;
		case 106: return 2.8000000000000000;
		case 107: return 2.9000000000000000;
		case 108: return 3.0000000000000000;
		case 109: return 2.5000000000000000;
		case 110: return 2.9000000000000000;
		case 111: return 3.1000000000000000;
		case 112: return 3.2000000000000000;
		case 113: return 3.3000000000000000;
		case 114: return 3.4000000000000000;
		case 115: return 3.5000000000000000;
		case 116: return 3.2000000000000000;
		case 117: return 3.3000000000000000;
		case 118: return 3.3000000000000000;
		case 119: return 3.5000000000000000;
		case 120: return 3.6000000000000000;
		case 121: return 3.6000000000000000;
		case 122: return 3.7000000000000000;
		case 123: return 3.8000000000000000;
		case 124: return 3.9000000000000000;
		case 125: return 3.6000000000000000;
		case 126: return 3.7000000000000000;
		case 127: return 3.8000000000000000;
		case 128: return 3.9000000000000000;
		case 129: return 3.9000000000000000;
		case 130: return 4.0000000000000000;
		case 131: return 4.1000000000000000;
		case 132: return 4.2000000000000000;
		case 133: return 4.0000000000000000;
		case 134: return 4.1000000000000000;
		case 135: return 4.2000000000000000;
		case 136: return 4.3000000000000000;
		case 137: return 4.4000000000000000;
		case 138: return 4.5000000000000000;
		}
		//return HIndex() * ((Ploidy() - 1) * 0.5);
	}

	/* SS within genotype under SMM model */
	TARGET double GENOTYPE::SS_SMM(ushort *alen2, int k2)
	{
		double smm = 0;
		ushort *als = GetAlleleArray();
		int ploidy = Ploidy();
		for (int i1 = 0; i1 < ploidy; ++i1)
			for (int i2 = 0; i2 < i1; ++i2)
				smm += alen2[als[i1] * k2 + als[i2]];
		return smm / ploidy;
	}

	/* The multinomial coefficient for HWE/RCS genotype frequency */
	TARGET double GENOTYPE::HWECoef()
	{
		/*
		switch (pattern)
		{
		default: return NA;
		case 0x1ull: return 1;
		case 0x2ull: return 1;
		case 0x11ull: return 2;
		case 0x3ull: return 1;
		case 0x21ull: return 3;
		case 0x111ull: return 6;
		case 0x4ull: return 1;
		case 0x31ull: return 4;
		case 0x22ull: return 6;
		case 0x211ull: return 12;
		case 0x1111ull: return 24;
		case 0x5ull: return 1;
		case 0x41ull: return 5;
		case 0x32ull: return 10;
		case 0x311ull: return 20;
		case 0x221ull: return 30;
		case 0x2111ull: return 60;
		case 0x11111ull: return 120;
		case 0x6ull: return 1;
		case 0x51ull: return 6;
		case 0x42ull: return 15;
		case 0x411ull: return 30;
		case 0x33ull: return 20;
		case 0x321ull: return 60;
		case 0x3111ull: return 120;
		case 0x222ull: return 90;
		case 0x2211ull: return 180;
		case 0x21111ull: return 360;
		case 0x111111ull: return 720;
		case 0x7ull: return 1;
		case 0x61ull: return 7;
		case 0x52ull: return 21;
		case 0x511ull: return 42;
		case 0x43ull: return 35;
		case 0x421ull: return 105;
		case 0x4111ull: return 210;
		case 0x331ull: return 140;
		case 0x322ull: return 210;
		case 0x3211ull: return 420;
		case 0x31111ull: return 840;
		case 0x2221ull: return 630;
		case 0x22111ull: return 1260;
		case 0x211111ull: return 2520;
		case 0x1111111ull: return 5040;
		case 0x8ull: return 1;
		case 0x71ull: return 8;
		case 0x62ull: return 28;
		case 0x611ull: return 56;
		case 0x53ull: return 56;
		case 0x521ull: return 168;
		case 0x5111ull: return 336;
		case 0x44ull: return 70;
		case 0x431ull: return 280;
		case 0x422ull: return 420;
		case 0x4211ull: return 840;
		case 0x41111ull: return 1680;
		case 0x332ull: return 560;
		case 0x3311ull: return 1120;
		case 0x3221ull: return 1680;
		case 0x32111ull: return 3360;
		case 0x311111ull: return 6720;
		case 0x2222ull: return 2520;
		case 0x22211ull: return 5040;
		case 0x221111ull: return 10080;
		case 0x2111111ull: return 20160;
		case 0x11111111ull: return 40320;
		case 0x9ull: return 1;
		case 0x81ull: return 9;
		case 0x72ull: return 36;
		case 0x711ull: return 72;
		case 0x63ull: return 84;
		case 0x621ull: return 252;
		case 0x6111ull: return 504;
		case 0x54ull: return 126;
		case 0x531ull: return 504;
		case 0x522ull: return 756;
		case 0x5211ull: return 1512;
		case 0x51111ull: return 3024;
		case 0x441ull: return 630;
		case 0x432ull: return 1260;
		case 0x4311ull: return 2520;
		case 0x4221ull: return 3780;
		case 0x42111ull: return 7560;
		case 0x411111ull: return 15120;
		case 0x333ull: return 1680;
		case 0x3321ull: return 5040;
		case 0x33111ull: return 10080;
		case 0x3222ull: return 7560;
		case 0x32211ull: return 15120;
		case 0x321111ull: return 30240;
		case 0x3111111ull: return 60480;
		case 0x22221ull: return 22680;
		case 0x222111ull: return 45360;
		case 0x2211111ull: return 90720;
		case 0x21111111ull: return 181440;
		case 0x111111111ull: return 362880;
		case 0xAull: return 1;
		case 0x91ull: return 10;
		case 0x82ull: return 45;
		case 0x811ull: return 90;
		case 0x73ull: return 120;
		case 0x721ull: return 360;
		case 0x7111ull: return 720;
		case 0x64ull: return 210;
		case 0x631ull: return 840;
		case 0x622ull: return 1260;
		case 0x6211ull: return 2520;
		case 0x61111ull: return 5040;
		case 0x55ull: return 252;
		case 0x541ull: return 1260;
		case 0x532ull: return 2520;
		case 0x5311ull: return 5040;
		case 0x5221ull: return 7560;
		case 0x52111ull: return 15120;
		case 0x511111ull: return 30240;
		case 0x442ull: return 3150;
		case 0x4411ull: return 6300;
		case 0x433ull: return 4200;
		case 0x4321ull: return 12600;
		case 0x43111ull: return 25200;
		case 0x4222ull: return 18900;
		case 0x42211ull: return 37800;
		case 0x421111ull: return 75600;
		case 0x4111111ull: return 151200;
		case 0x3331ull: return 16800;
		case 0x3322ull: return 25200;
		case 0x33211ull: return 50400;
		case 0x331111ull: return 100800;
		case 0x32221ull: return 75600;
		case 0x322111ull: return 151200;
		case 0x3211111ull: return 302400;
		case 0x31111111ull: return 604800;
		case 0x22222ull: return 113400;
		case 0x222211ull: return 226800;
		case 0x2221111ull: return 453600;
		case 0x22111111ull: return 907200;
		case 0x211111111ull: return 1814400;
		case 0x1111111111ull: return 3628800;
		}
		*/
		switch (patternid)
		{
			default: 	return NA;
			case	1: 	return 1;
			case	2: 	return 1;
			case	3: 	return 2;
			case	4: 	return 1;
			case	5: 	return 3;
			case	6: 	return 6;
			case	7: 	return 1;
			case	8: 	return 4;
			case	9: 	return 6;
			case	10: 	return 12;
			case	11: 	return 24;
			case	12: 	return 1;
			case	13: 	return 5;
			case	14: 	return 10;
			case	15: 	return 20;
			case	16: 	return 30;
			case	17: 	return 60;
			case	18: 	return 120;
			case	19: 	return 1;
			case	20: 	return 6;
			case	21: 	return 15;
			case	22: 	return 30;
			case	23: 	return 20;
			case	24: 	return 60;
			case	25: 	return 120;
			case	26: 	return 90;
			case	27: 	return 180;
			case	28: 	return 360;
			case	29: 	return 720;
			case	30: 	return 1;
			case	31: 	return 7;
			case	32: 	return 21;
			case	33: 	return 42;
			case	34: 	return 35;
			case	35: 	return 105;
			case	36: 	return 210;
			case	37: 	return 140;
			case	38: 	return 210;
			case	39: 	return 420;
			case	40: 	return 840;
			case	41: 	return 630;
			case	42: 	return 1260;
			case	43: 	return 2520;
			case	44: 	return 5040;
			case	45: 	return 1;
			case	46: 	return 8;
			case	47: 	return 28;
			case	48: 	return 56;
			case	49: 	return 56;
			case	50: 	return 168;
			case	51: 	return 336;
			case	52: 	return 70;
			case	53: 	return 280;
			case	54: 	return 420;
			case	55: 	return 840;
			case	56: 	return 1680;
			case	57: 	return 560;
			case	58: 	return 1120;
			case	59: 	return 1680;
			case	60: 	return 3360;
			case	61: 	return 6720;
			case	62: 	return 2520;
			case	63: 	return 5040;
			case	64: 	return 10080;
			case	65: 	return 20160;
			case	66: 	return 40320;
			case	67: 	return 1;
			case	68: 	return 9;
			case	69: 	return 36;
			case	70: 	return 72;
			case	71: 	return 84;
			case	72: 	return 252;
			case	73: 	return 504;
			case	74: 	return 126;
			case	75: 	return 504;
			case	76: 	return 756;
			case	77: 	return 1512;
			case	78: 	return 3024;
			case	79: 	return 630;
			case	80: 	return 1260;
			case	81: 	return 2520;
			case	82: 	return 3780;
			case	83: 	return 7560;
			case	84: 	return 15120;
			case	85: 	return 1680;
			case	86: 	return 5040;
			case	87: 	return 10080;
			case	88: 	return 7560;
			case	89: 	return 15120;
			case	90: 	return 30240;
			case	91: 	return 60480;
			case	92: 	return 22680;
			case	93: 	return 45360;
			case	94: 	return 90720;
			case	95: 	return 181440;
			case	96: 	return 362880;
			case	97: 	return 1;
			case	98: 	return 10;
			case	99: 	return 45;
			case	100: 	return 90;
			case	101: 	return 120;
			case	102: 	return 360;
			case	103: 	return 720;
			case	104: 	return 210;
			case	105: 	return 840;
			case	106: 	return 1260;
			case	107: 	return 2520;
			case	108: 	return 5040;
			case	109: 	return 252;
			case	110: 	return 1260;
			case	111: 	return 2520;
			case	112: 	return 5040;
			case	113: 	return 7560;
			case	114: 	return 15120;
			case	115: 	return 30240;
			case	116: 	return 3150;
			case	117: 	return 6300;
			case	118: 	return 4200;
			case	119: 	return 12600;
			case	120: 	return 25200;
			case	121: 	return 18900;
			case	122: 	return 37800;
			case	123: 	return 75600;
			case	124: 	return 151200;
			case	125: 	return 16800;
			case	126: 	return 25200;
			case	127: 	return 50400;
			case	128: 	return 100800;
			case	129: 	return 75600;
			case	130: 	return 151200;
			case	131: 	return 302400;
			case	132: 	return 604800;
			case	133: 	return 113400;
			case	134: 	return 226800;
			case	135: 	return 453600;
			case	136: 	return 907200;
			case	137: 	return 1814400;
			case	138: 	return 3628800;
		}
	}

	/* Genotypic frequency for zygotes under inbreeding */
	TARGET double GENOTYPE::GFZ(int *allele_count, int sum, double f)
	{
		// allele_count[i] / sum is the allele frequency

		int ploidy = Ploidy(), nalleles = Nalleles();
		if (nalleles == 0) return 1;

		double re = HWECoef(), re2 = 1, F = (1.0 / f - 1);
		if (F < 1e-20) F = 1e-20;  if (F > 1e20) F = 1e20;
		double FP = F / sum;
		uint64 ap = GetPattern();
		ushort *als = GetAlleleArray();

		for (int a = nalleles - 1; a >= 0; --a)
		{
			for (int j = 0; j < (int)(ap & 0xF); ++j)
				re *= allele_count[als[ploidy + a]] * FP + j;
			ap >>= 4;
		}
		for (int j = 0; j < ploidy; ++j)
			re2 *= F + j;
		return re / re2;
	}

	/* Genotypic frequency for zygotes under specific double-reduction model */
	TARGET double GENOTYPE::GFZ(int DR_MODE, double *f)
	{
		int ploidy = Ploidy(), nalleles = Nalleles();
		if (nalleles == 0) return 1;

		ushort *a = GetAlleleArray() + ploidy;
		if (DR_MODE == 0 || ploidy <= 2 || ploidy % 2 == 1 || ploidy > N_MAX_PLOIDY)
		{
			double re = HWECoef();
			uint64 ap = GetPattern();
			for (int ai = nalleles - 1; ai >= 0; --ai)
			{
				re *= pow(f[a[ai]], (int)(ap & 0xF));
				ap >>= 4;
			}
			return re;
		}

		double *alpha = &ALPHA[DR_MODE][ploidy][0];
		switch (patternid)
		{
			case	7: 		return GFZ4_iiii(alpha[1], f[a[0]]);
			case	8: 		return GFZ4_iiij(alpha[1], f[a[0]], f[a[1]]);
			case	9: 		return GFZ4_iijj(alpha[1], f[a[0]], f[a[1]]);
			case	10: 	return GFZ4_iijk(alpha[1], f[a[0]], f[a[1]], f[a[2]]);
			case	11: 	return GFZ4_ijkl(alpha[1], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);

			case	19: 	return GFZ6_iiiiii(alpha[1], f[a[0]]);
			case	20: 	return GFZ6_iiiiij(alpha[1], f[a[0]], f[a[1]]);
			case	21: 	return GFZ6_iiiijj(alpha[1], f[a[0]], f[a[1]]);
			case	22: 	return GFZ6_iiiijk(alpha[1], f[a[0]], f[a[1]], f[a[2]]);
			case	23: 	return GFZ6_iiijjj(alpha[1], f[a[0]], f[a[1]]);
			case	24: 	return GFZ6_iiijjk(alpha[1], f[a[0]], f[a[1]], f[a[2]]);
			case	25: 	return GFZ6_iiijkl(alpha[1], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	26: 	return GFZ6_iijjkk(alpha[1], f[a[0]], f[a[1]], f[a[2]]);
			case	27: 	return GFZ6_iijjkl(alpha[1], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	28: 	return GFZ6_iijklm(alpha[1], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]]);
			case	29: 	return GFZ6_ijklmn(alpha[1], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]]);

			case	45: 	return GFZ8_iiiiiiii(alpha[1], alpha[2], f[a[0]]);
			case	46: 	return GFZ8_iiiiiiij(alpha[1], alpha[2], f[a[0]], f[a[1]]);
			case	47: 	return GFZ8_iiiiiijj(alpha[1], alpha[2], f[a[0]], f[a[1]]);
			case	48: 	return GFZ8_iiiiiijk(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]]);
			case	49: 	return GFZ8_iiiiijjj(alpha[1], alpha[2], f[a[0]], f[a[1]]);
			case	50: 	return GFZ8_iiiiijjk(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]]);
			case	51: 	return GFZ8_iiiiijkl(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	52: 	return GFZ8_iiiijjjj(alpha[1], alpha[2], f[a[0]], f[a[1]]);
			case	53: 	return GFZ8_iiiijjjk(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]]);
			case	54: 	return GFZ8_iiiijjkk(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]]);
			case	55: 	return GFZ8_iiiijjkl(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	56: 	return GFZ8_iiiijklm(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]]);
			case	57: 	return GFZ8_iiijjjkk(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]]);
			case	58: 	return GFZ8_iiijjjkl(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	59: 	return GFZ8_iiijjkkl(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	60: 	return GFZ8_iiijjklm(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]]);
			case	61: 	return GFZ8_iiijklmn(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]]);
			case	62: 	return GFZ8_iijjkkll(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	63: 	return GFZ8_iijjkklm(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]]);
			case	64: 	return GFZ8_iijjklmn(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]]);
			case	65: 	return GFZ8_iijklmno(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]], f[a[6]]);
			case	66: 	return GFZ8_ijklmnop(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]], f[a[6]], f[a[7]]);

			case	97: 	return GFZ10_iiiiiiiiii(alpha[1], alpha[2], f[a[0]]);
			case	98: 	return GFZ10_iiiiiiiiij(alpha[1], alpha[2], f[a[0]], f[a[1]]);
			case	99: 	return GFZ10_iiiiiiiijj(alpha[1], alpha[2], f[a[0]], f[a[1]]);
			case	100: 	return GFZ10_iiiiiiiijk(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]]);
			case	101: 	return GFZ10_iiiiiiijjj(alpha[1], alpha[2], f[a[0]], f[a[1]]);
			case	102: 	return GFZ10_iiiiiiijjk(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]]);
			case	103: 	return GFZ10_iiiiiiijkl(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	104: 	return GFZ10_iiiiiijjjj(alpha[1], alpha[2], f[a[0]], f[a[1]]);
			case	105: 	return GFZ10_iiiiiijjjk(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]]);
			case	106: 	return GFZ10_iiiiiijjkk(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]]);
			case	107: 	return GFZ10_iiiiiijjkl(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	108: 	return GFZ10_iiiiiijklm(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]]);
			case	109: 	return GFZ10_iiiiijjjjj(alpha[1], alpha[2], f[a[0]], f[a[1]]);
			case	110: 	return GFZ10_iiiiijjjjk(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]]);
			case	111: 	return GFZ10_iiiiijjjkk(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]]);
			case	112: 	return GFZ10_iiiiijjjkl(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	113: 	return GFZ10_iiiiijjkkl(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	114: 	return GFZ10_iiiiijjklm(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]]);
			case	115: 	return GFZ10_iiiiijklmn(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]]);
			case	116: 	return GFZ10_iiiijjjjkk(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]]);
			case	117: 	return GFZ10_iiiijjjjkl(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	118: 	return GFZ10_iiiijjjkkk(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]]);
			case	119: 	return GFZ10_iiiijjjkkl(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	120: 	return GFZ10_iiiijjjklm(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]]);
			case	121: 	return GFZ10_iiiijjkkll(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	122: 	return GFZ10_iiiijjkklm(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]]);
			case	123: 	return GFZ10_iiiijjklmn(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]]);
			case	124: 	return GFZ10_iiiijklmno(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]], f[a[6]]);
			case	125: 	return GFZ10_iiijjjkkkl(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	126: 	return GFZ10_iiijjjkkll(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]]);
			case	127: 	return GFZ10_iiijjjkklm(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]]);
			case	128: 	return GFZ10_iiijjjklmn(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]]);
			case	129: 	return GFZ10_iiijjkkllm(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]]);
			case	130: 	return GFZ10_iiijjkklmn(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]]);
			case	131: 	return GFZ10_iiijjklmno(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]], f[a[6]]);
			case	132: 	return GFZ10_iiijklmnop(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]], f[a[6]], f[a[7]]);
			case	133: 	return GFZ10_iijjkkllmm(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]]);
			case	134: 	return GFZ10_iijjkkllmn(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]]);
			case	135: 	return GFZ10_iijjkklmno(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]], f[a[6]]);
			case	136: 	return GFZ10_iijjklmnop(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]], f[a[6]], f[a[7]]);
			case	137: 	return GFZ10_iijklmnopq(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]], f[a[6]], f[a[7]], f[a[8]]);
			case	138: 	return GFZ10_ijklmnopqr(alpha[1], alpha[2], f[a[0]], f[a[1]], f[a[2]], f[a[3]], f[a[4]], f[a[5]], f[a[6]], f[a[7]], f[a[8]], f[a[9]]);

			default	:	return -1;
		}
	}

	/* Number of copies of target allele */
	TARGET int GENOTYPE::GetAlleleCount(int a)
	{
		int ploidy = Ploidy(), nalleles = Nalleles();
		if (nalleles == 0)
			return 0;

		ushort *als = GetAlleleArray();
		int count = 0;
		for (int i = 0; i < ploidy; ++i)
			if (als[i] == a)
				count++;
		return count;
	}

	/* Frequency of target allele in this genotype */
	TARGET double GENOTYPE::GetFreq(int a)
	{
		int ploidy = Ploidy(), nalleles = Nalleles();
		if (nalleles == 0)
			return NA;

		int count = 0;
		ushort *als = GetAlleleArray();
		for (int i = 0; i < ploidy; ++i)
			if (als[i] == a)
				count++;
		return count / (double)ploidy;
	}

	/* Frequencies of all alleles in this genotype */
	TARGET void GENOTYPE::GetFreq(double *p, int k2)
	{
		int ploidy = Ploidy(), nalleles = Nalleles();
		SetZero(p, k2);
		if (nalleles == 0) return;

		double f = 1.0 / ploidy;
		ushort *als = GetAlleleArray();
		for (int i = 0; i < ploidy; ++i)
			p[als[i]] += f;
	}

	/* Obtain Genepop genotype string */
	TARGET char *GENOTYPE::GetGenepopStr()
	{
		int ploidy = Ploidy(), nalleles = Nalleles();
		if (ploidy != 2) Exit("\nError: Genepop supports diploids only.\n");

		char *str = NULL;
		conversion_memory2[threadid].Alloc(str, 8);

		ushort *als = GetAlleleArray();
		if (nalleles == 0)
			sprintf(str, " 000000");
		else
			sprintf(str, " %03d%03d", als[0] + 1, als[1] + 1);

		return str;
	}

	/* Obtain Arlequin genotype string */
	TARGET char *GENOTYPE::GetArlequinStr()
	{
		int ploidy = Ploidy(), nalleles = Nalleles();
		if (ploidy != 2) Exit("\nError: Arlequin supports diploids only.\n");
		char *str = NULL;
		conversion_memory2[threadid].Alloc(str, 10);

		ushort *als = GetAlleleArray();
		if (nalleles == 0)
		{
			sprintf(str, " ?");
			sprintf(str + 5, " ?");
		}
		else
		{
			sprintf(str, " %d", als[0] + 1);
			if (ploidy > 1)
				sprintf(str + 5, " %d", als[1] + 1);
			else
				sprintf(str + 5, " ?");
		}
		return str;
	}

	/* Obtain Structure genotype string */
	TARGET char *GENOTYPE::GetStructureStr()
	{
		int ploidy = Ploidy(), nalleles = Nalleles();
		char *str = NULL;
		conversion_memory2[threadid].Alloc(str, ploidy * 5);
		char *re2 = str;

		ushort *als = GetAlleleArray();
		if (nalleles == 0)
			for (byte i = 0; i < ploidy; ++i)
			{
				sprintf(re2, " -9");
				re2 += 5;
			}
		else
			for (byte i = 0; i < ploidy; ++i)
			{
				sprintf(re2, " %d", als[i] + 1);
				re2 += 5;
			}
		return str;
	}

	/* Obtain Spagedi genotype string */
	TARGET char *GENOTYPE::GetSpagediStr()
	{
		int ploidy = Ploidy(), nalleles = Nalleles();
		char *str = NULL;
		conversion_memory2[threadid].Alloc(str, ploidy * 3 + 3);

		ushort *als = GetAlleleArray();
		if (nalleles == 0)
			sprintf(str, "\t0");
		else
		{
			char *re2 = str;
			sprintf(re2, "\t");
			while (*re2) re2++;
			for (byte i = 0; i < ploidy; ++i)
			{
				sprintf(re2, "%03d", als[i] + 1);
				while (*re2) re2++;
			}
		}
		return str;
	}

	/* Obtain PolyRelatedness genotype string */
	TARGET char *GENOTYPE::GetPolyrelatednessStr()
	{
		int ploidy = Ploidy(), nalleles = Nalleles();
		char *str = NULL;
		conversion_memory2[threadid].Alloc(str, ploidy * 3 + 3);

		ushort *als = GetAlleleArray();
		char *re2 = str;
		sprintf(re2, "\t");  while (*re2) re2++;
		for (byte i = 0; i < ploidy; ++i)
		{
			if (nalleles == 0) sprintf(re2, "000");
			else sprintf(re2, "%03d", als[i] + 1);
			while (*re2) re2++;
		}
		return str;
	}

	/* Obtain PolyGene genotype string */
	TARGET char *GENOTYPE::GetPolygeneStr()
	{
		int ploidy = Ploidy(), nalleles = Nalleles();
		char *str = NULL;
		conversion_memory2[threadid].Alloc(str, ploidy * 4 + 3);

		ushort *als = GetAlleleArray();
		if (nalleles == 0)
			sprintf(str, "\t");
		else
		{
			char *re2 = str;
			sprintf(re2, "\t");
			while (*re2) re2++;
			for (byte i = 0; i < ploidy; ++i)
			{
				if (i != 0) sprintf(re2, ",%d", als[i] + 1);
				else sprintf(re2, "%d", als[i] + 1);
				while (*re2) re2++;
			}
		}
		return str;
	}

	/* Obtain Cervus genotype string */
	TARGET char *GENOTYPE::GetCervusStr()
	{
		int ploidy = Ploidy(), nalleles = Nalleles();
		if (ploidy != 2) Exit("\nError: Cervus supports diploids only.\n");
		char *str = NULL;
		conversion_memory2[threadid].Alloc(str, 9);

		ushort *als = GetAlleleArray();
		if (nalleles == 0)
			sprintf(str, ",,");
		else
			sprintf(str, ",%d,%d", als[0] + 1, als[1] + 1);
		return str;
	}
#endif

#ifndef _SLOCUS
	TARGET SLOCUS::SLOCUS()
	{
		SetZero(this, 1);
	};

	/* Convert from LOCUS */
	TARGET SLOCUS::SLOCUS(MEMORY &memory, LOCUS &ref)
	{
		SetVal(this, (SLOCUS*)&ref, 1);

		int gsize = ngeno;
		int asize = flag_alen ? (k + k * k) : 0;
		int ssize = ref.GetEnd() - ref.GetChrom();
		int gasize = ref.GetGenoAlleleSize();//genotype allele array

		byte *bucket = NULL;
		memory.Alloc(bucket, gsize * sizeof(GENOTYPE) + asize * sizeof(ushort) + ssize + gasize * sizeof(ushort));
		bits1 = (uint64)bucket;

		//copy genotype array
		SetVal(GetGtab(), ref.GetGtab(), gsize);

		//copy allele length array
		if (flag_alen)
			SetVal(GetAlenArray(), ref.GetAlenArray(), asize);

		//copy chrom ...
		SetVal(GetChrom(), ref.GetChrom(), ssize);

		//copy genotype allele array
		SetVal(GetGenoAlleleArray(), ref.GetGenoAlleleArray(), gasize);

		//set genotype allele offset
		GENOTYPE *gtab = GetGtab();
		ushort *gatab = GetGenoAlleleArray();
		for (uint gi = 0; gi < ngeno; ++gi)
		{
			GENOTYPE &gt = gtab[gi];
			if (gt.Nalleles())
			{
				gt.SetAlleleArray(gatab);
				gatab += gt.Nalleles() + gt.Ploidy();
			}
			else  //missing genotype
				gt.SetAlleleArray((ushort*)-1);
		}
	}

	/* Create unphase locus */
	TARGET SLOCUS::SLOCUS(MEMORY &memory, SLOCUS &ref, TABLE<HASH, uint> &gitab, ushort *gtmap)
	{
		SetVal(this, &ref, 1);
		ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
		gitab.Clear();
		GENOTYPE *refgtab = ref.GetGtab();

		int gasize = 0, gsize = 0;
		for (uint i = 0; i < ref.ngeno; ++i)
		{
			GENOTYPE &gt = refgtab[i];
			int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();
			HASH ha;
			if (nalleles)
			{
				ushort *als = gt.GetAlleleArray();
				SetVal(alleles, als, ploidy);
				Sort(alleles, ploidy);
				ha = HashGenotype(alleles, ploidy);
			}
			else
				ha = missing_hash[ploidy];

			if (!gitab.ContainsKey(ha))
			{
				gitab[ha] = gtmap[i] = gitab.size;
				gasize += ploidy + nalleles;
				gsize++;
			}
		}

		int asize = flag_alen ? (k + k * k) : 0;
		int ssize = ref.GetEnd() - ref.GetChrom();
		int tsize = gsize * sizeof(GENOTYPE) + asize * sizeof(ushort) + ssize + gasize * sizeof(ushort);

		byte *bucket = NULL;
		memory.Alloc(bucket, tsize);
		bits1 = (uint64)bucket;
		ngeno = gitab.size;

		//copy genotype array
		SetVal(GetGtab(), ref.GetGtab(), gsize);

		//copy allele length array
		if (flag_alen)
			SetVal(GetAlenArray(), ref.GetAlenArray(), asize);

		//copy chrom ...
		SetVal(GetChrom(), ref.GetChrom(), ssize);

		//copy genotype allele array
		SetVal(GetGenoAlleleArray(), ref.GetGenoAlleleArray(), gasize);

		//construct gtab
		GENOTYPE *gtab = GetGtab();
		ushort *gatab = GetGenoAlleleArray();
		for (uint i = 0, ng = 0; i < ref.ngeno; ++i)
		{
			GENOTYPE &gt = refgtab[i];
			int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();

			HASH ha;
			if (nalleles)
			{
				ushort *als = gt.GetAlleleArray();
				SetVal(alleles, als, ploidy);
				Sort(alleles, ploidy);
				ha = HashGenotype(alleles, ploidy);
			}
			else
			{
				ha = missing_hash[ploidy];
				SetFF(alleles, ploidy);
			}

			if (gtmap[i] == ng)
				new(&gtab[ng++]) GENOTYPE(gatab, alleles, ploidy);
		}
	}

	/* Create locus for haplotype extraction and Chi-square test */
	TARGET SLOCUS::SLOCUS(MEMORY &memory, SLOCUS &ref, int64 _id, int _ngeno, int _gasize, TABLE<HASH, TEMP_GENOTYPE> &temptab)
	{
		//do not need copy alen
		SetZero(this, 1);

		ngeno = _ngeno;
		flag_pass = true;
		flag_alen = false;

		int gsize = ngeno = _ngeno;
		int asize = 0;
		int ssize = (int)strlen(ref.GetChrom()) + 1 + 3 + CeilLog10(_id + 1) + 1 + 1;
		int gasize = _gasize;
		int tsize = gsize * sizeof(GENOTYPE) + asize * sizeof(ushort) + ssize + gasize * sizeof(ushort);

		byte *bucket = memory.Alloc(tsize);
		bits1 = (uint64)bucket;

		//no alen
		
		//chrom
		SetZero(GetChrom(), ssize);
		strcpy(GetChrom(), ref.GetChrom());
		sprintf(GetName(), "Loc%llu", _id + 1);
		
		//no allele identifier!!!
		GENOTYPE *gtab = GetGtab();
		ushort *gatab = (ushort*)(bucket + gsize * sizeof(GENOTYPE) + asize * sizeof(ushort) + ssize);
		for (uint gi = 0; gi < ngeno; ++gi)
		{
			TEMP_GENOTYPE &tgt = temptab(gi);
			new(gtab++) GENOTYPE(gatab, tgt.alleles, tgt.ploidy);
		}
	}

	/* Deep copy from SLOCUS */
	TARGET SLOCUS::SLOCUS(MEMORY &memory, SLOCUS &ref)
	{
		SetVal(this, &ref, 1);
		int tsize = 0;

		//genotype alleles array
		GENOTYPE *refgtab = ref.GetGtab();
		for (uint gi = 0; gi < ngeno; ++gi)
		{
			int nalleles = refgtab[gi].Nalleles();
			tsize += nalleles ? nalleles + refgtab[gi].Ploidy() : 0;
		}
		tsize *= sizeof(ushort);
		tsize += (uint64)ref.GetEnd() - ref.bits1;

		byte *bucket = NULL;
		memory.Alloc(bucket, tsize);
		bits1 = (uint64)bucket;
		SetVal(bucket, (byte*)ref.GetGtab(), tsize);
	}

	/* Get Genotype array */
	TARGET GENOTYPE *SLOCUS::GetGtab()
	{
		return (GENOTYPE*)bits1;
	}

	/* Get end of chrom \0 name \0 {(allele identifiers \0)[k] */
	TARGET char *SLOCUS::GetEnd()
	{
		return StrNextIdx0(GetChrom(), 2 + (ALLELE_IDENTIFIER ? k : 0)) + 1;
	}

	/* Get chrom string */
	TARGET char *SLOCUS::GetChrom()
	{
		return (char*)bits1 + ngeno * sizeof(GENOTYPE) + (flag_alen ? (k + k * k) * sizeof(ushort) : 0);
	}

	/* Get locus identifier */
	TARGET char *SLOCUS::GetName()
	{
		return StrNextIdx0(GetChrom(), 1) + 1;
	}

	/* Get allele name for vcf/bcf */
	TARGET char *SLOCUS::GetAlleleName(int a)
	{
		return ALLELE_IDENTIFIER ? StrNextIdx0(GetChrom(), 2 + a) + 1 : NULL;
	}

	/* Get alen array */
	TARGET ushort *SLOCUS::GetAlenArray()
	{
		return flag_alen ? (ushort*)(bits1 + ngeno * sizeof(GENOTYPE)) : 0;
	}

	/* Get genotype allele array */
	TARGET ushort *SLOCUS::GetGenoAlleleArray()
	{
		return (ushort*)GetEnd();
	}

	/* Get Genotype alleles array size Sum(ploidy+nalleles)*/
	TARGET uint SLOCUS::GetGenoAlleleSize()
	{
		uint gasize = 0;
		GENOTYPE *gtab = GetGtab();
		for (uint gi = 0; gi < ngeno; ++gi)
		{
			GENOTYPE gt = gtab[gi];
			if (gt.Nalleles() != 0)
				gasize += gt.Nalleles() + gt.Ploidy();
		}
		return gasize;
	}

	/* Get SMM distance */
	TARGET ushort SLOCUS::GetSMMDist(int a, int b)
	{
		return *((ushort*)(bits1 + ngeno * sizeof(GENOTYPE)) + k + a * k + b);
	}

	/*
	TARGET void SLOCUS::GetGfTab(TABLE<HASH, GENOTYPE*> &gftab)
	{
		GENOTYPE *gtab = GetGtab();
		for (uint gi = 0; gi < ngeno; ++gi)
		{
			GENOTYPE &gt = gtab[gi];
			//if (!gt.Nalleles()) continue;
			HASH hash = gt.Hash();
			if (!gftab.ContainsKey(hash))
				gftab[hash] = &gt;
		}
	}
	*/
#endif

#ifndef _LOCUS

	/* Initialize */
	TARGET LOCUS::LOCUS()
	{
		SetZero(this, 1);
	}

	/* Deep Copy Locus */
	TARGET LOCUS::LOCUS(MEMORY &memory, int64 _id, LOCUS &ref)
	{
		SetVal(this, &ref, 1);
		id = _id;

		//copy SLOCUS
		new(this) SLOCUS(memory, ref);

		GENOTYPE *refgtab = ref.GetGtab();
		int gsize = ref.ngeno; gasize = 0;
		for (uint gi = 0; gi < ngeno; ++gi)
		{
			GENOTYPE &gt = refgtab[gi];
			int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();
			if (nalleles)
				gasize += ploidy + nalleles;
		}
		int asize = flag_alen ? (k + k * k) : 0;
		int ssize = ref.GetEnd() - ref.GetChrom();
		int tsize = gsize * sizeof(GENOTYPE) + asize * sizeof(ushort) + ssize + gasize * sizeof(ushort);

		//alloc memory
		byte *bucket = NULL;
		memory.Alloc(bucket, tsize);
		bits1 = (uint64)bucket;
		ngeno = gsize;

		//set interal variables
		_chrom = (char*)bits1 + ngeno * sizeof(GENOTYPE) + (flag_alen ? (k + k * k) * sizeof(ushort) : 0);
		_alen = flag_alen ? (ushort*)(bits1 + ngeno * sizeof(GENOTYPE)) : 0;
		_genoallele = (ushort*)(StrNextIdx0(GetChrom(), 2 + (ALLELE_IDENTIFIER ? k : 0)) + 1);

		//copy genotype array
		SetVal(GetGtab(), ref.GetGtab(), gsize);

		//copy allele length array
		if (flag_alen)
			SetVal(GetAlenArray(), ref.GetAlenArray(), asize);

		//copy chrom ...
		SetVal(GetChrom(), ref.GetChrom(), ssize);

		//copy genotype allele array
		SetVal(GetGenoAlleleArray(), ref.GetGenoAlleleArray(), gasize);

		//set genotype allele offset
		GENOTYPE *gtab = GetGtab();
		ushort *gatab = GetGenoAlleleArray();
		for (uint gi = 0; gi < ngeno; ++gi)
		{
			GENOTYPE &gt = gtab[gi];
			if (gt.Nalleles())
			{
				gt.SetAlleleArray(gatab);
				gatab += gt.Nalleles() + gt.Ploidy();
			}
			else  //missing genotype
				gt.SetAlleleArray((ushort*)-1);
		}

		//construct gftab
		new(&gftab) TABLE<HASH, GENOTYPE*>(true, &memory, ngeno);
		auto bucket2 = ref.gftab.bucket;
		auto index2  = ref.gftab.index;
		for (uint gi = 0; gi < ngeno; ++gi)
			gftab[bucket2[index2[gi]].key] = &gtab[gi];
	}

	/* Create unphase locus */
	TARGET LOCUS::LOCUS(MEMORY &memory, LOCUS &ref, TABLE<HASH, uint> &gitab, ushort *gtmap)
	{
		SetVal(this, &ref, 1);

		ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
		gitab.Clear();

		GENOTYPE *refgtab = ref.GetGtab();
		int gsize = 0; gasize = 0;
		for (uint gi = 0; gi < ngeno; ++gi)
		{
			GENOTYPE &gt = refgtab[gi];
			int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();
			HASH ha;
			if (nalleles)
			{
				ushort *als = gt.GetAlleleArray();
				SetVal(alleles, als, ploidy);
				Sort(alleles, ploidy);
				ha = HashGenotype(alleles, ploidy);
			}
			else
				ha = missing_hash[ploidy];

			if (!gitab.ContainsKey(ha))
			{
				gitab[ha] = gtmap[gi] = gitab.size;
				gasize += ploidy + nalleles;
				gsize++;
			}
		}

		int asize = flag_alen ? (k + k * k) : 0;
		int ssize = ref.GetEnd() - ref.GetChrom();
		int tsize = gsize * sizeof(GENOTYPE) + asize * sizeof(ushort) + ssize + gasize * sizeof(ushort);

		//alloc memory
		byte *bucket = NULL;
		memory.Alloc(bucket, tsize);
		bits1 = (uint64)bucket;
		ngeno = gsize;

		//set internal variables
		_chrom = (char*)bits1 + ngeno * sizeof(GENOTYPE) + (flag_alen ? (k + k * k) * sizeof(ushort) : 0);
		_alen = flag_alen ? (ushort*)(bits1 + ngeno * sizeof(GENOTYPE)) : 0;
		_genoallele = (ushort*)(StrNextIdx0(GetChrom(), 2 + (ALLELE_IDENTIFIER ? k : 0)) + 1);

		//copy genotype array
		SetVal(GetGtab(), ref.GetGtab(), gsize);

		//copy allele length array
		if (flag_alen)
			SetVal(GetAlenArray(), ref.GetAlenArray(), asize);

		//copy chrom ...
		SetVal(GetChrom(), ref.GetChrom(), ssize);

		//copy genotype allele array
		SetVal(GetGenoAlleleArray(), ref.GetGenoAlleleArray(), gasize);

		//construct gftab
		new(&gftab) TABLE<HASH, GENOTYPE*>(true, &memory, ngeno);

		//construct genotypes
		GENOTYPE *gtab = GetGtab();
		ushort *gatab = GetGenoAlleleArray();
		for (uint i = 0, ng = 0; i < ref.ngeno; ++i)
		{
			GENOTYPE &gt = refgtab[i];
			int ploidy = gt.Ploidy(), nalleles = gt.Nalleles();

			HASH ha;
			if (nalleles)
			{
				ushort *als = gt.GetAlleleArray();
				SetVal(alleles, als, ploidy);
				Sort(alleles, ploidy);
				ha = HashGenotype(alleles, ploidy);
			}
			else
			{
				ha = missing_hash[ploidy];
				SetFF(alleles, ploidy);
			}

			if (gtmap[i] == ng)
			{
				new(&gtab[ng]) GENOTYPE(gatab, alleles, ploidy);
				gftab[ha] = &gtab[ng++];
			}
		}
	}

	/* Create locus for haplotype extraction and Chi-square test */
	TARGET LOCUS::LOCUS(MEMORY &memory, LOCUS &ref, int64 _id, int _ngeno, int _gasize, TABLE<HASH, TEMP_GENOTYPE> &temptab)
	{
		//do not need copy alen
		SetZero(this, 1);
		ngeno = _ngeno;
		flag_pass = true;
		flag_alen = false;

		id = _id;
		pos = ref.pos;

		int gsize = ngeno;
		int asize = 0;
		int ssize = (int)strlen(ref.GetChrom()) + 1 + 3 + CeilLog10(id + 1) + 1 + 1;
		gasize = _gasize;

		byte *bucket = memory.Alloc(ngeno * sizeof(GENOTYPE) + gasize * sizeof(ushort));
		bits1 = (uint64)bucket;
		GENOTYPE *gtab = (GENOTYPE*)bits1;

		_genoallele = (ushort*)(gtab + ngeno);
		_alen = NULL;
		 memory.Alloc(_chrom, ssize);

		//no alen

		//chrom
		SetZero(GetChrom(), ssize);
		strcpy(GetChrom(), ref.GetChrom());
		sprintf(GetName(), "Loc%lu", id + 1);

		//no allele identifier!!!

		//gftab
		new(&gftab) TABLE<HASH, GENOTYPE*>(true, &memory, ngeno);
		ushort *gatab = GetGenoAlleleArray();
		for (uint gi = 0; gi < ngeno; ++gi)
		{
			TEMP_GENOTYPE &tgt = temptab(gi);
			gftab[tgt.hash] = new(gtab++) GENOTYPE(gatab, tgt.alleles, tgt.ploidy);
		}
	}

	/* For dummy locus for collapse alleles during testing genotype distributions */
	TARGET LOCUS::LOCUS(MEMORY &memory, SLOCUS &ref)
	{
		//do not need copy alen
		SetZero(this, 1);
		ngeno = ref.ngeno;
		id = (LOCN)-1;
		flag_pass = true;
		flag_original = true;
		flag_qual = true;

		//use original genotype and genotype allele memory
		bits1 = (uint64)ref.GetGtab();
		_genoallele = ref.GetGenoAlleleArray();

		//set chrom, alen, genoallelearray
		_chrom = ref.GetChrom();
		_alen = ref.GetAlenArray();
		_genoallele = ref.GetGenoAlleleArray();
		gasize = 0;

		//create gftab
		GENOTYPE *gtab = GetGtab();
		new(&gftab) TABLE<HASH, GENOTYPE*>(true, &memory, ngeno);
		for (uint gi = 0; gi < ngeno; ++gi)
		{
			GENOTYPE &gt = gtab[gi];
			HASH hash = gt.Hash();
			if (gt.Nalleles())
				gasize += gt.Nalleles() + gt.Ploidy();
			if (!gftab.ContainsKey(hash))
				gftab[hash] = &gt;
		}
	}

	/* For non-vcf input, set locus name and id, SLOCUS */
	TARGET LOCUS::LOCUS(MEMORY &memory, char *line, int64 _id, int _ngenotype, GENOTYPE *&gtab, ushort *&gatab)
	{
		//_alen is set at IndexAlleleLength
		int _gasize = gasize;
		SetZero(this, 1);
		ngeno = _ngenotype;
		gasize = _gasize;
		id = _id;
		flag_pass = true;
		flag_original = true;
		flag_qual = true;
		flag_indel = false;

		gtab = (GENOTYPE*)memory.Alloc(ngeno * sizeof(GENOTYPE) + gasize * sizeof(ushort));
		bits1 = (uint64)gtab;
		_genoallele = gatab = (ushort*)(gtab + ngeno);

		new(&gftab) TABLE<HASH, GENOTYPE*>(true, &memory, ngeno);

		if (line)
		{
			//no chrom, set locus name
			int namelen = 0;
			while (line[namelen] &&
				line[namelen] != ' '  && line[namelen] != '\t' &&
				line[namelen] != '\r' && line[namelen] != '\n')
				namelen++;

			memory.Alloc(_chrom, 1 + namelen + 1 + 1);
			SetZero(_chrom, namelen + 3);
			SetVal(_chrom + 1, line, namelen);
		}
		else
		{
			//no chrom, no locus name, name it
			memory.Alloc(_chrom, 1 + 3 + CeilLog10(id + 1) + 1 + 1);
			SetZero(_chrom, CeilLog10(id + 1) + 5);
			sprintf(_chrom + 1, "Loc%lu", id + 1);
		}
	}

	/* For vcf input, set locus name and id, SLOCUS */
	TARGET LOCUS::LOCUS(MEMORY &memory, char *&line, uint64 _mask, int _ngenotype, GENOTYPE *&gtab, ushort *&gatab)
	{
		//do not need copy alen
		int _gasize = gasize;
		SetZero(this, 1);
		ngeno = _ngenotype;
		gasize = _gasize;
		pos = (uint64)-1;
		flag_pass = true;
		gtab = (GENOTYPE*)memory.Alloc(ngeno * sizeof(GENOTYPE) + gasize * sizeof(ushort));
		bits1 = (uint64)gtab;
		_genoallele = gatab = (ushort*)(gtab + ngeno);

		new(&gftab) TABLE<HASH, GENOTYPE*>(true, &memory, ngeno);

		//locate address of each string
		char *base = line, *pstr = NULL;
		while (*line != '\t') line++; *line++ = '\0';

		if (line[0] != '.') pos = ReadLong(line);
		while (*line != '\t') line++; *line++ = '\0';

		char *name = line;
		while (*line != '\t') line++; *line++ = '\0';

		char *ref = line;
		while (*line != '\t') line++; *line++ = '\0';

		char *alt = line;
		while (*line != '\t') line++; *line++ = '\0';

		//write chrom string
		if ((name[0] == '.' && name[1] == '\0') || name[0] == '\0')
		{
			//no loc name, name it
			int slen = (int)(strlen(base) + 1 + strlen(ref) + 1 + strlen(alt) + 1) * 2 + CeilLog10(pos) + 1;

			//alloc
			memory.Alloc(_chrom, slen);
			pstr = _chrom;
			SetZero(pstr, slen);

			//copy chrom
			strcpy(pstr, base);
			pstr += strlen(base) + 1;

			//set name
			sprintf(pstr, "%s_%lld_%s_%s", pstr, pos, ref, alt);
		}
		else
		{
			//has loc name
			int slen = (int)(strlen(base) + 1 + strlen(name) + 1 + strlen(ref) + 1 + strlen(alt) + 1);

			//alloc
			memory.Alloc(_chrom, slen);
			pstr = _chrom;
			SetZero(pstr, slen);

			//copy chrom
			strcpy(pstr, base);
			pstr += strlen(base) + 1;

			//copy name
			strcpy(pstr, name);
		}

		//relocate name
		name = pstr;

		//copy allele identifiers
		pstr += strlen(name) + 1;
		sprintf(pstr, "%s,%s", ref, alt);

		//obtain number of alleles
		k = (ushort)(CountChar(alt, ',') + 2);
		ReplaceChar(pstr, ',', '\0');
		ReplaceChar(alt, ',', '\0');

		//SNP or indel
		int maxal = 0, minal = 65536000;
		for (int i = 0; i < k; ++i)
		{
			int tl = (int)strlen(pstr);
			if (tl > maxal) maxal = tl;
			if (tl < minal) minal = tl;
			pstr += tl + 1;
		}
		flag_indel = minal != maxal;

		//QUAL to FILTER
		while (*line != '\t') line++; *line++ = '\0';

		if (info_filter)								//MASK INFO Filter
		{
			if (f_qual_b && (_mask & 0x4))				//fail in qual filter
				flag_pass = false;
			if (f_type_b)
			{
				if (f_type_val == 1 && flag_indel)		//snp
					flag_pass = false;
				else if (f_type_val == 2 && !flag_indel) //indel
					flag_pass = false;
			}
			if (f_original_b && f_original_val == 1 && (_mask & 0x2)) //fail in original filter
				flag_pass = false;
		}

		//FILTER to INFO
		while (*line != '\t') line++; *line++ = '\0';

		//INFO to FORMAT
		for (;;)
		{
			while (*line != ';' && *line != '\t') line++;
			if (*line == ';')
				*line++ = '\0';
			else
			{
				*line++ = '\0';
				break;
			}
		}

		char *end = StrNextIdx(line + 1, '\t', 1);
		*end = '\0';
		format_size = (ushort)CountChar(line, ':') + 1;
		*end = '\t';

		//find offset of each format string (from base)
		VLA_NEW(format_offset, ushort, format_size);
		for (int i = 0; line < end; ++i)
		{
			format_offset[i] = (ushort)(line - base);
			while (*line != ':' && *line != '\t') line++;
			if (*line == ':')
				*line++ = '\0';
			else
			{
				*line++ = '\0';
				break;
			}
		}

		dpid = gqid = adid = 0xFFFF;
		if (abs(g_format_val) <= 2)
		{
			gtid = GetFormatId(base, (char*)"GT", format_offset);
			if (gtid == 0xFFFF) //must have gt
				Exit("\nExit: variant %s does not has GT format_offset tag. \n", name);

			gqid = GetFormatId(base, (char*)"GQ", format_offset);
			flag_hasgq = gqid != 0xFFFF;

			dpid = GetFormatId(base, (char*)"DP", format_offset);
			flag_hasdp = dpid != 0xFFFF;

			adid = GetFormatId(base, (char*)"AD", format_offset);
			flag_hasad = adid != 0xFFFF;
		}
		VLA_DELETE(format_offset);
	}

	/* Get end of chrom \0 name \0 {(allele identifiers \0)[k] */
	TARGET char *LOCUS::GetEnd()
	{
		return StrNextIdx0(GetChrom(), 2 + (ALLELE_IDENTIFIER ? k : 0)) + 1;
	}

	/* Get chrom string */
	TARGET char *LOCUS::GetChrom()
	{
		return _chrom;
	}

	/* Get locus identifier */
	TARGET char *LOCUS::GetName()
	{
		return StrNextIdx0(GetChrom(), 1) + 1;
	}

	/* Get allele name for vcf/bcf */
	TARGET char *LOCUS::GetAlleleName(int a)
	{
		return ALLELE_IDENTIFIER ? StrNextIdx0(GetChrom(), 2 + a) + 1 : NULL;
	}

	/* Get alen array */
	TARGET ushort *LOCUS::GetAlenArray()
	{
		return _alen;
	}

	/* Get genotype allele array */
	TARGET ushort *LOCUS::GetGenoAlleleArray()
	{
		return _genoallele;
	}

	/* Get Genotype alleles array size Sum(ploidy+nalleles)*/
	TARGET uint LOCUS::GetGenoAlleleSize()
	{
		if (gasize == 0)
		{
			GENOTYPE *gtab = GetGtab();
			for (uint gi = 0; gi < ngeno; ++gi)
			{
				GENOTYPE gt = gtab[gi];
				if (gt.Nalleles() != 0)
					gasize += gt.Nalleles() + gt.Ploidy();
			}
		}
		return gasize;
	}

	/* Get index of a target format */
	TARGET ushort LOCUS::GetFormatId(char *base, char *fmt_name, ushort *fmt_offset)
	{
		for (int i = 0; i < format_size; ++i)
			if (LwrParCmp(fmt_name, base + fmt_offset[i]) == 0)
				return i;
		return 0xFFFF;
	}
#endif

#ifndef _GENOTYPE_ITERATOR

	/* Do nothing */
	TARGET GENO_ITERATOR::GENO_ITERATOR()
	{

	}

	/* Initialize reader/writer */
	TARGET GENO_ITERATOR::GENO_ITERATOR(int indid, int64 l, bool isread, byte *bucket, OFFSET *offset)
	{
		if (bucket == NULL)
		{
			bucket = genotype_bucket;
			offset = genotype_offset;
		}
		size = (uint)offset[l].size;
		mask64 = mask = (1u << size) - 1u;
		pos = (uint*)(bucket + offset[l].offset);

		nskipbits = indid * size;
		pos += nskipbits / 32;
		nskipbits %= 32;

		if (isread)
		{
			//read 32 bits into data and skip skipnbits bits
			data = (*pos++) >> nskipbits;
			nbits = 32 - nskipbits;
			nskipbits = 0;
		}
		else
		{
			data = 0;
			mask64 <<= nskipbits;
			nbits = nskipbits;
		}
	}

	/* Get id of next ind (order by indid) */
	TARGET int GENO_ITERATOR::Read()
	{
		//if data is full, read additional 32 bits into high bits
		if (nbits <= 32)
		{
			data = (data << (32 - nbits)) | ((uint64)(*pos++) << 32);
			data >>= (32 - nbits);
			nbits += 32;
		}

		//read low (size) bits and shift right
		uint gid = (uint)data & mask;
		data >>= size;
		nbits -= size;
		return (int)gid;
	}

	/* Write id of next ind to buffer */
	TARGET void GENO_ITERATOR::Write(uint gid)
	{
		data |= ((uint64)gid << nbits) & mask64;
		mask64 <<= size;
		nbits += size;

		//Write back into pos if data has more than 32 bits
		if (nbits >= 32)
		{
			if (nskipbits)
			{
				//Don't change low (skipnbits) bits
				*pos++ = (uint)data | (*pos & ((1u << nskipbits) - 1u));
				nskipbits = 0;
			}
			else
				*pos++ = (uint)data;

			data >>= 32;
			mask64 >>= 32;
			nbits -= 32;
		}
	}

	/* Write all remaining bits to buffer */
	TARGET void GENO_ITERATOR::FinishWrite()
	{
		atomic<uint> &posa = *(atomic<uint>*)pos;
		uint odata = 0, pdata = 0, qdata = 0;

		//read 32 bits
		pdata = odata = *pos;
		if (nskipbits)
		{
			do 
			{
				//save original nskipbits bits
				qdata = pdata & ((1u << nskipbits) - 1u);

				//set low nbits to 0
				pdata >>= nbits;
				pdata <<= nbits;

				//set low bits to data and write back
			} 
			while (!posa.compare_exchange_strong(odata, data | pdata | qdata));
		}
		else
		{
			do
			{
				//set low nbits to 0
				pdata >>= nbits;
				pdata <<= nbits;

				//set low bits to data and write back
			}
			while (!posa.compare_exchange_strong(odata, data | pdata));
		}

		SetZero(this, 1);
	}
#endif

#ifndef _IND
	/* Initialize */
	TARGET IND::IND()
	{
		indid = 0xFFFF;
		name = NULL;
	}

	/* Create individual for non-vcf input */
	TARGET IND::IND(char *t, bool iscount, int id, GENOTYPE **gtab, ushort **gatab, GENO_ITERATOR *iter)
	{
		indid = id;
		switch (abs(g_format_val))//genepop|spagedi|cervus|arlequin|structure|polygene|polyrelatedness
		{
		case 3: genepop(t, iscount, gtab, gatab, iter); break;
		case 4: spagedi(t, iscount, gtab, gatab, iter); break;
		case 5: cervus(t, iscount, gtab, gatab, iter); break;
		case 6: arlequin(t, iscount, gtab, gatab, iter); break;
		case 7: structure(t, iscount, gtab, gatab, iter); break;
		case 8: polygene(t, iscount, gtab, gatab, iter); break;
		case 9: polyrelatedness(t, iscount, gtab, gatab, iter); break;
		}
	}

	/* Create individual for vcf/bcf input */
	TARGET IND::IND(char *&title, int id)
	{
		//from VCF line
		name = title;

		while (*title != '\t' && *title != '\n' && *title != '\0') title++;
		*title++ = '\0';

		char *tname;
		individual_memory->Alloc(tname, (int)strlen(name) + 1);
		//tname = new char[strlen(name) + 1];
		strcpy(tname, name);
		name = tname;
		indid = id;
	}

	/* Create individual from a reference individual */
	TARGET IND::IND(IND& ref)
	{
		SetVal(this, &ref, 1);
		individual_memory->Alloc(name, (int)strlen(ref.name) + 1);
		strcpy(name, ref.name);
	}

	/* Unnitialize */
	TARGET IND::~IND()
	{

	}

	/* Set allele sequencing depth, for ad, TEST */
	TARGET /*static*/ void IND::SetAlleleDepth(int64 l, uint *depth, int K, int indid)
	{
		uint size = (uint)alleledepth_offset[l].size, mask = (1u << size) - 1u;
		uint64 offset = size * K * indid;
		uint *pos = (uint*)(alleledepth_bucket + alleledepth_offset[l].offset + (offset >> 3));
		offset &= 7;
		for (int k = 0; k < K; ++k)
		{
			*pos = (*pos & (~(mask << offset))) | (depth[k] << offset);
			offset += size;
			if (offset > 7)
			{
				pos = (uint*)((byte*)pos + 1);
				offset -= 8;
			}
		}
	}

	/* Set allele sequencing depth, for ad, TEST */
	TARGET void IND::SetAlleleDepth(int64 l, uint *depth, int K)
	{
		uint size = (uint)alleledepth_offset[l].size, mask = (1u << size) - 1u;
		uint64 offset = size * K * indid;
		uint *pos = (uint*)(alleledepth_bucket + alleledepth_offset[l].offset + (offset >> 3));
		offset &= 7;
		for (int k = 0; k < K; ++k)
		{
			*pos = (*pos & (~(mask << offset))) | (depth[k] << offset);
			offset += size;
			if (offset > 7)
			{
				pos = (uint*)((byte*)pos + 1);
				offset -= 8;
			}
		}
	}

	/* Set allele sequencing depth, for ad, TEST */
	TARGET void IND::SetAlleleDepth(int64 l, uint *depth, int K, OFFSET *_offset, byte *bucket)
	{
		uint size = (uint)_offset[l].size;
		uint mask = (1u << size) - 1u;
		uint64 offset = size * K * indid;
		uint *pos = (uint*)(bucket + _offset[l].offset + (offset >> 3));
		offset &= 7;
		for (int k = 0; k < K; ++k)
		{
			*pos = (*pos & (~(mask << offset))) | (depth[k] << offset);
			offset += size;
			if (offset > 7)
			{
				pos = (uint*)((byte*)pos + 1);
				offset -= 8;
			}
		}
	}

	/* Set allele sequencing depth, for ad, TEST */
	TARGET void IND::GetAlleleDepth(int64 l, uint *depth)
	{
		uint K = GetLoc(l).k;
		uint size = (uint)alleledepth_offset[l].size;
		uint mask = (1u << size) - 1u;
		uint64 offset = size * K * indid;
		uint *pos = (uint*)(alleledepth_bucket + alleledepth_offset[l].offset + (offset >> 3));
		offset &= 7;
		for (uint k = 0; k < K; ++k)
		{
			depth[k] = (*pos >> offset) & mask;
			offset += size;
			if (offset > 7)
			{
				pos = (uint*)((byte*)pos + 1);
				offset -= 8;
			}
		}
	}

	/* Set individual genotype with default bucket */
	TARGET void IND::SetGenotype(int64 l, uint gid)
	{
		uint size = (uint)genotype_offset[l].size;
		uint64 offset = size * indid;
		uint *pos = (uint*)(genotype_bucket + genotype_offset[l].offset + (offset >> 3));
		offset &= 7;
		*pos = (*pos & (~(((1u << size) - 1u) << offset))) | (gid << offset);
	}

	/* Set individual genotype with local bucket */
	TARGET void IND::SetGenotype(int64 l, uint gid, OFFSET *_offset, byte *bucket)
	{
		uint size = (uint)_offset[l].size;
		uint64 offset = size * indid;
		uint *pos = (uint*)(bucket + _offset[l].offset + (offset >> 3));
		offset &= 7;
		*pos = (*pos & (~(((1u << size) - 1u) << offset))) | (gid << offset);
	}

	/* Get index for a pair of genotype */
	TARGET /*static*/ void IND::GetDyadGenotypeIdx(int &id1, int &id2, int64 l)
	{
		OFFSET ot = genotype_offset[l];
		uint size = (uint)ot.size;
		byte *bucket = genotype_bucket + ot.offset;
		uint mask = (1u << size) - 1u;
		uint64 o1 = size * id1, o2 = size * id2;
		uint *pos1 = (uint*)(bucket + (o1 >> 3));
		uint *pos2 = (uint*)(bucket + (o2 >> 3));
		id1 = (*pos1 >> (o1 & 7)) & mask;
		id2 = (*pos2 >> (o2 & 7)) & mask;
	}

	/* Get individual genotype index from default table */
	TARGET int IND::GetGenotypeId(int64 l)
	{
		uint size = (uint)genotype_offset[l].size;
		uint64 offset = size * indid;
		uint *pos = (uint*)(genotype_bucket + genotype_offset[l].offset + (offset >> 3));
		return (*pos >> (offset & 7)) & ((1u << size) - 1u);
	}

	/* Get individual genotype from default table */
	TARGET GENOTYPE &IND::GetGenotype(int64 l)
	{
		return GetLoc(l).GetGtab()[GetGenotypeId(l)];
	}

	/* Get individual genotype from local table */
	TARGET GENOTYPE &IND::GetGenotype(int64 l, TABLE<HASH, GENOTYPE*> &gftab)
	{
		return *gftab(GetGenotypeId(l));
	}

	/* Get individual genotype from local table */
	TARGET GENOTYPE &IND::GetGenotype(int64 l, GENOTYPE *gtab)
	{
		return gtab[GetGenotypeId(l)];
	}

	/* Create individual from genepop */
	TARGET void IND::genepop(char *t, bool iscount, GENOTYPE**gtab, ushort **gatab, GENO_ITERATOR *iter)
	{
		//two rounds, first round obtain the number of genotypes
		name = t;
		char *genstr = NULL;

		genstr = StrNextIdx(t, ',', 1);
		char *tname = genstr - 1;
		*genstr++ = '\0';

		if (!iscount)
		{
			while (*tname == ' ') *tname-- = '\0';
			individual_memory->Alloc(tname, (int)strlen(name) + 1);
			strcpy(tname, name);
			name = tname;
		}
		else name = NULL;

		int ploidy = 2;
		for (int64 l = 0; l < nloc; ++l)
		{
			TABLE<HASH, GENOTYPE*> &gftab = locus[l].gftab;
			TABLE<HASH, uint> &gfid = nvcf_gfid[l];

			while (*genstr == ' ' || *genstr == '\t') genstr++;
			char *gend = genstr;
			while (*gend >= '0' && *gend <= '9') gend++;
			*gend = '\0';

			int64 glen = strlen(genstr);
			ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
			if (glen == 4)
			{
				char b = genstr[2];
				genstr[2] = '\0';
				alleles[0] = (ushort)ReadIntegerKeep(genstr);
				genstr[2] = b;
				alleles[1] = (ushort)ReadIntegerKeep(genstr + 2);
			}
			else if (glen == 6)
			{
				char b = genstr[3];
				genstr[3] = '\0';
				alleles[0] = (ushort)ReadIntegerKeep(genstr);
				genstr[3] = b;
				alleles[1] = (ushort)ReadIntegerKeep(genstr + 3);
			}
			else
				Exit("\nError: Format error in individual %s at locus %s.\n", name, locus[l].GetName());

			if (alleles[0] == '\0' && alleles[1] == '\0')
				SetFF(alleles, ploidy);
			Sort(alleles, ploidy);//unphase
			
			if (iscount)
			{
				HASH mha = missing_hash[ploidy];
				if (!gfid.ContainsKey(mha))
				{
					int tid = gfid.size;
					gfid[mha] = tid;
					locus[l].gasize += 0;
				}

				HASH hash = HashGenotype(alleles, ploidy);
				if (!gfid.ContainsKey(hash))
				{
					int tid = gfid.size;
					gfid[hash] = tid;
					locus[l].gasize += ploidy + GetNalleles(alleles, ploidy);
				}
			}
			else
			{
				HASH mha = missing_hash[ploidy];
				uint mid = gfid[mha];
				if (!gftab.ContainsKey(mha))
					gftab[mha] = new(gtab[l]++) GENOTYPE(gatab[l], missing_genotype[ploidy]);

				HASH hash = HashGenotype(alleles, ploidy);
				uint gid = gfid[hash];
				if (!gftab.ContainsKey(hash))
					gftab[hash] = new(gtab[l]++) GENOTYPE(gatab[l], alleles, ploidy);

				iter[l].Write(genotype_filter && f_ploidy_b && (f_ploidy_min > ploidy || ploidy > f_ploidy_max) ? mid : gid);
			}

			genstr = gend + 1;
		}
	}

	/* Create individual from cervus */
	TARGET void IND::cervus(char *t, bool iscount, GENOTYPE **gtab, ushort **gatab, GENO_ITERATOR *iter)
	{
		int extracol = genotype_extracol;
		name = t;
		char *genstr = NULL;

		genstr = StrNextIdx(t, ',', 1);
		char *tname = genstr - 1;
		*genstr++ = '\0';

		if (!iscount)
		{
			while (*tname == ' ') *tname-- = '\0';
			individual_memory->Alloc(tname, (int)strlen(name) + 1);
			strcpy(tname, name);
			name = tname;
		}
		else name = NULL;

		genstr = StrNextIdx(genstr, ',', extracol) + 1;

		int v = 2;
		for (int64 l = 0; l < nloc; ++l)
		{
			TABLE<HASH, GENOTYPE*> &gftab = locus[l].gftab;
			TABLE<HASH, uint> &gfid = nvcf_gfid[l];

			ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
			if (*genstr != ',' && *genstr != '\r' && *genstr != '\n' && *genstr != '\0')
				alleles[0] = (ushort)ReadInteger(genstr);
			genstr++;
			if (*genstr != ',' && *genstr != '\r' && *genstr != '\n' && *genstr != '\0')
				alleles[1] = (ushort)ReadInteger(genstr);
			genstr++;

			if (alleles[0] == 0xFFFF && alleles[1] == 0xFFFF)
				SetFF(alleles, v);
			Sort(alleles, v);//unphase

			if (iscount)
			{
				HASH mha = missing_hash[v];
				if (!gfid.ContainsKey(mha))
				{
					int tid = gfid.size;
					gfid[mha] = tid;
					locus[l].gasize += 0;
				}

				HASH hash = HashGenotype(alleles, v);
				if (!gfid.ContainsKey(hash))
				{
					int tid = gfid.size;
					gfid[hash] = tid;
					locus[l].gasize += v + GetNalleles(alleles, v);
				}
			}
			else
			{
				HASH mha = missing_hash[v];
				uint mid = gfid[mha];
				if (!gftab.ContainsKey(mha))
					gftab[mha] = new(gtab[l]++) GENOTYPE(gatab[l], missing_genotype[v]);

				HASH hash = HashGenotype(alleles, v);
				uint gid = gfid[hash];
				if (!gftab.ContainsKey(hash))
					gftab[hash] = new(gtab[l]++) GENOTYPE(gatab[l], alleles, v);

				iter[l].Write(genotype_filter && f_ploidy_b && (f_ploidy_min > v || v > f_ploidy_max) ? mid : gid);
			}
		}
	}

	/* Create individual from spagedi */
	TARGET void IND::spagedi(char *t, bool iscount, GENOTYPE **gtab, ushort **gatab, GENO_ITERATOR *iter)
	{
		name = t;
		char *genstr = NULL;

		genstr = StrNextIdx(t, '\t', 1);
		char *tname = genstr - 1;
		*genstr++ = '\0';

		if (!iscount)
		{
			while (*tname == ' ') *tname-- = '\0';
			individual_memory->Alloc(tname, (int)strlen(name) + 1);
			//tname = new char[strlen(name) + 1];
			strcpy(tname, name);
			name = tname;
		}
		else name = NULL;

		genstr = StrNextIdx(genstr, '\t', genotype_extracol + 1) + 1; //1 (pop) + extracols

		int ndigit = genotype_digit;
		int maxv = 0;
		for (int64 l = 0; l < nloc; ++l)
		{
			TABLE<HASH, GENOTYPE*> &gftab = locus[l].gftab;
			TABLE<HASH, uint> &gfid = nvcf_gfid[l];

			while (*genstr == '\t') genstr++;
			ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
			int v = 0;
			bool ismissing = false;
			while (*genstr != '\t' && *genstr != '\r' && *genstr != '\n' && *genstr != '\0')
			{
				alleles[v] = (ushort)ReadIntegerSpagedi(genstr, ndigit);
				if (alleles[v] == 0)
					ismissing = true;
				if (++v > N_MAX_PLOIDY)
					Exit("\nError: ploidy level greater than %d in individual %s at locus %s.\n", N_MAX_PLOIDY, name, locus[l].GetName());
			}
			maxv = Max(v, maxv);

			if (!iscount && ismissing)
			{
				v = maxploidy;
				SetFF(alleles, v);
			}
			Sort(alleles, v);//unphase

			if (iscount)
			{
				HASH mha = missing_hash[v];
				if (!gfid.ContainsKey(mha))
				{
					int tid = gfid.size;
					gfid[mha] = tid;
					locus[l].gasize += 0;
				}

				HASH hash = HashGenotype(alleles, v);
				if (!gfid.ContainsKey(hash))
				{
					int tid = gfid.size;
					gfid[hash] = tid;
					locus[l].gasize += v + GetNalleles(alleles, v);
				}
			}
			else
			{
				HASH mha = missing_hash[v];
				uint mid = gfid[mha];
				if (!gftab.ContainsKey(mha))
					gftab[mha] = new(gtab[l]++) GENOTYPE(gatab[l], missing_genotype[v]);

				HASH hash = HashGenotype(alleles, v);
				uint gid = gfid[hash];
				if (!gftab.ContainsKey(hash))
					gftab[hash] = new(gtab[l]++) GENOTYPE(gatab[l], alleles, v);

				iter[l].Write(genotype_filter && f_ploidy_b && (f_ploidy_min > v || v > f_ploidy_max) ? mid : gid);
			}
		}
		maxploidy = maxv;
	}

	/* Create individual from arlequin */
	TARGET void IND::arlequin(char *t, bool iscount, GENOTYPE **gtab, ushort **gatab, GENO_ITERATOR *iter)
	{
		ReplaceChar(t, '\t', ' ');
		while (*t == ' ') t++;
		name = t;
		char *genstr = t + strlen(t) - 1;
		while (*genstr == ' ' || *genstr == '\r' || *genstr == '\n') *genstr-- = '\0';

		int64 nline = CountChar(t, '\n') + 1;
		if (t[strlen(t) - 1] == '\n')  nline--;

		char *lines[N_MAX_PLOIDY];
		genstr = t;
		for (int i = 0; i < nline; ++i)
		{
			while (*genstr == ' ') genstr++;
			if (i == 0) genstr = StrNextIdx(genstr, ' ', 2) + 1; //name 1 geno
			lines[i] = genstr;
			genstr = StrNextIdx(genstr, '\n', 1) + 1;
		}

		genstr = StrNextIdx(t, ' ', 1);
		char *tname = genstr - 1;
		*genstr++ = '\0';

		if (!iscount)
		{
			while (*tname == ' ') *tname-- = '\0';
			individual_memory->Alloc(tname, (int)strlen(name) + 1);
			//tname = new char[strlen(name) + 1];
			strcpy(tname, name);
			name = tname;
		}
		else name = NULL;

		int v = (int)nline;
		for (int64 l = 0; l < nloc; ++l)
		{
			TABLE<HASH, GENOTYPE*> &gftab = locus[l].gftab;
			TABLE<HASH, uint> &gfid = nvcf_gfid[l];

			if (v > N_MAX_PLOIDY)
				Exit("\nError: ploidy level greater than %d in individual %s at locus %s.\n", N_MAX_PLOIDY, name, locus[l].GetName());

			ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
			bool ismissing = false;

			for (int j = 0; j < v; ++j)
			{
				alleles[j] = (ushort)ReadInteger(lines[j]);
				if (alleles[j] == 0xFFFF) 
					ismissing = true; 
			}

			if (ismissing) SetFF(alleles, v);
			Sort(alleles, v);//unphase

			if (iscount)
			{
				HASH mha = missing_hash[v];
				if (!gfid.ContainsKey(mha))
				{
					int tid = gfid.size;
					gfid[mha] = tid;
					locus[l].gasize += 0;
				}

				HASH hash = HashGenotype(alleles, v);
				if (!gfid.ContainsKey(hash))
				{
					int tid = gfid.size;
					gfid[hash] = tid;
					locus[l].gasize += v + GetNalleles(alleles, v);
				}
			}
			else
			{
				HASH mha = missing_hash[v];
				uint mid = gfid[mha];
				if (!gftab.ContainsKey(mha))
					gftab[mha] = new(gtab[l]++) GENOTYPE(gatab[l], missing_genotype[v]);

				HASH hash = HashGenotype(alleles, v);
				uint gid = gfid[hash];
				if (!gftab.ContainsKey(hash))
					gftab[hash] = new(gtab[l]++) GENOTYPE(gatab[l], alleles, v);

				iter[l].Write(genotype_filter && f_ploidy_b && (f_ploidy_min > v || v > f_ploidy_max) ? mid : gid);
			}
		}
	}

	/* Create individual from structure */
	TARGET void IND::structure(char *t, bool iscount, GENOTYPE **gtab, ushort **gatab, GENO_ITERATOR *iter)
	{
		name = t;

		int64 nline = CountChar(t, '\n') + 1;
		if (t[strlen(t) - 1] == '\n')  nline--;

		char *lines[N_MAX_PLOIDY];
		char *genstr = t;
		for (int i = 0; i < nline; ++i)
		{
			genstr = StrNextIdx(genstr, ' ', 2) + 1; //skip name, pop
			lines[i] = genstr;
			genstr = StrNextIdx(genstr, '\n', 1) + 1;
		}

		genstr = StrNextIdx(t, ' ', 1);
		char *tname = genstr - 1;
		*genstr++ = '\0';

		if (!iscount)
		{
			while (*tname == ' ') *tname-- = '\0';
			individual_memory->Alloc(tname, (int)strlen(name) + 1);
			//tname = new char[strlen(name) + 1];
			strcpy(tname, name);
			name = tname;
		}
		else name = NULL;

		int v = (int)nline;
		for (int64 l = 0; l < nloc; ++l)
		{
			TABLE<HASH, GENOTYPE*> &gftab = locus[l].gftab;
			TABLE<HASH, uint> &gfid = nvcf_gfid[l];

			if (v > N_MAX_PLOIDY)
				Exit("\nError: ploidy level greater than %d in individual %s at locus %s.\n", N_MAX_PLOIDY, name, locus[l].GetName());

			ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
			bool ismissing = false;

			for (int j = 0; j < v; ++j)
			{
				alleles[j] = (ushort)ReadInteger(lines[j]);
				if (alleles[j] == 0xFFF7)//-9 missing
					ismissing = true;
			}

			if (ismissing) SetFF(alleles, v);
			Sort(alleles, v);//unphase

			if (iscount)
			{
				HASH mha = missing_hash[v];
				if (!gfid.ContainsKey(mha))
				{
					int tid = gfid.size;
					gfid[mha] = tid;
					locus[l].gasize += 0;
				}

				HASH hash = HashGenotype(alleles, v);
				if (!gfid.ContainsKey(hash))
				{
					int tid = gfid.size;
					gfid[hash] = tid;
					locus[l].gasize += v + GetNalleles(alleles, v);
				}
			}
			else
			{
				HASH mha = missing_hash[v];
				uint mid = gfid[mha];
				if (!gftab.ContainsKey(mha))
					gftab[mha] = new(gtab[l]++) GENOTYPE(gatab[l], missing_genotype[v]);

				HASH hash = HashGenotype(alleles, v);
				uint gid = gfid[hash];
				if (!gftab.ContainsKey(hash))
					gftab[hash] = new(gtab[l]++) GENOTYPE(gatab[l], alleles, v);

				iter[l].Write(genotype_filter && f_ploidy_b && (f_ploidy_min > v || v > f_ploidy_max) ? mid : gid);
			}
		}
	}

	/* Create individual from polyrelatedness */
	TARGET void IND::polyrelatedness(char *t, bool iscount, GENOTYPE **gtab, ushort **gatab, GENO_ITERATOR *iter)
	{
		name = t;
		char *genstr = NULL;

		genstr = StrNextIdx(t, '\t', 1);
		char *tname = genstr - 1;
		*genstr++ = '\0';

		if (!iscount)
		{
			while (*tname == ' ') *tname-- = '\0';
			individual_memory->Alloc(tname, (int)strlen(name) + 1);
			//tname = new char[strlen(name) + 1];
			strcpy(tname, name);
			name = tname;
		}
		else name = NULL;

		genstr = StrNextIdx(genstr, '\t', 1) + 1;

		int ndigit = genotype_digit;
		for (int64 l = 0; l < nloc; ++l)
		{
			TABLE<HASH, GENOTYPE*> &gftab = locus[l].gftab;
			TABLE<HASH, uint> &gfid = nvcf_gfid[l];

			while (*genstr == '\t') genstr++;
			ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
			int v = 0;
			bool ismissing = false;
			while (*genstr != '\t' && *genstr != '\r' && *genstr != '\n' && *genstr != ' ' && *genstr != '\0')
			{
				alleles[v] = (ushort)ReadIntegerSpagedi(genstr, ndigit);
				if (alleles[v] == genotype_missing)
					ismissing = true;
				if (alleles[v] == genotype_ambiguous)
					Exit("\nError: do not support ambiguous alleles in individual %s at locus %s\n.", name, locus[l].GetName());
				if (++v > N_MAX_PLOIDY)
					Exit("\nError: ploidy level greater than %d in individual %s at locus %s.\n", N_MAX_PLOIDY, name, locus[l].GetName());
			}

			if (ismissing) SetFF(alleles, v);
			Sort(alleles, v);//unphase

			if (iscount)
			{
				HASH mha = missing_hash[v];
				if (!gfid.ContainsKey(mha))
				{
					int tid = gfid.size;
					gfid[mha] = tid;
					locus[l].gasize += 0;
				}

				HASH hash = HashGenotype(alleles, v);
				if (!gfid.ContainsKey(hash))
				{
					int tid = gfid.size;
					gfid[hash] = tid;
					locus[l].gasize += v + GetNalleles(alleles, v);
				}
			}
			else
			{
				HASH mha = missing_hash[v];
				uint mid = gfid[mha];
				if (!gftab.ContainsKey(mha))
					gftab[mha] = new(gtab[l]++) GENOTYPE(gatab[l], missing_genotype[v]);

				HASH hash = HashGenotype(alleles, v);
				uint gid = gfid[hash];
				if (!gftab.ContainsKey(hash))
					gftab[hash] = new(gtab[l]++) GENOTYPE(gatab[l], alleles, v);

				iter[l].Write(genotype_filter && f_ploidy_b && (f_ploidy_min > v || v > f_ploidy_max) ? mid : gid);
			}
		}
	}

	/* Create individual from polygene */
	TARGET void IND::polygene(char *t, bool iscount, GENOTYPE **gtab, ushort **gatab, GENO_ITERATOR *iter)
	{
		name = t;
		char *genstr = NULL;

		genstr = StrNextIdx(t, '\t', 1);
		char *tname = genstr - 1;
		*genstr++ = '\0';

		if (!iscount)
		{
			while (*tname == ' ') *tname-- = '\0';
			individual_memory->Alloc(tname, (int)strlen(name) + 1);
			//tname = new char[strlen(name) + 1];
			strcpy(tname, name);
			name = tname;
		}
		else name = NULL;

		genstr = StrNextIdx(genstr, '\t', 1) + 1;
		int ploidy = ReadInteger(genstr);
		if (ploidy <= 0 || ploidy > N_MAX_PLOIDY)
			Exit("\nError: ploidy level greater than %d in individual %s.\n", N_MAX_PLOIDY, name);

		for (int64 l = 0; l < nloc; ++l)
		{
			TABLE<HASH, GENOTYPE*> &gftab = locus[l].gftab;
			TABLE<HASH, uint> &gfid = nvcf_gfid[l];
			if (*genstr == '\t') genstr++;
			ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };

			int64 len = StrNextIdx(genstr, '\t', 1) ? StrNextIdx(genstr, '\t', 1) - genstr : strlen(genstr);
			int v = (int)CountChar(genstr, ',', len);

			bool ismissing = false;
			if (len == 0)
				ismissing = true;
			else if (v && v != ploidy - 1)
				Exit("\nError: allelic phenotype in polygene format is not supported.\n");

			v = 0;
			while (*genstr != '\t' && *genstr != '\r' && *genstr != '\n' && *genstr != ' ' && *genstr != '\0')
				alleles[v++] = (ushort)ReadInteger(genstr);

			v = ploidy;
			if (ismissing) SetFF(alleles, v);
			Sort(alleles, v);//unphase

			if (iscount)
			{
				HASH mha = missing_hash[v];
				if (!gfid.ContainsKey(mha))
				{
					int tid = gfid.size;
					gfid[mha] = tid;
					locus[l].gasize += 0;
				}

				HASH hash = HashGenotype(alleles, v);
				if (!gfid.ContainsKey(hash))
				{
					int tid = gfid.size;
					gfid[hash] = tid;
					locus[l].gasize += v + GetNalleles(alleles, v);
				}
			}
			else
			{
				HASH mha = missing_hash[v];
				uint mid = gfid[mha];
				if (!gftab.ContainsKey(mha))
					gftab[mha] = new(gtab[l]++) GENOTYPE(gatab[l], missing_genotype[v]);

				HASH hash = HashGenotype(alleles, v);
				uint gid = gfid[hash];
				if (!gftab.ContainsKey(hash))
					gftab[hash] = new(gtab[l]++) GENOTYPE(gatab[l], alleles, v);

				iter[l].Write(genotype_filter && f_ploidy_b && (f_ploidy_min > v || v > f_ploidy_max) ? mid : gid);
			}
		}
	}

	/* Calculate the likelihood of genotype data */
	TARGET double IND::Likelihood(POP *grp, int model, int64 loc, double e)
	{
		//Logarithm
		double re = 0, re2 = 1;
		int64 st = loc == (int64)-1 ? 0 : loc;
		int64 ed = loc == (int64)-1 ? nloc : loc + 1;
		double e2 = e;
		double e1 = 1 - e2;

		OpenLog(re, re2);
		for (int64 l = st; l < ed; ++l)
		{
			GENOTYPE &gt = GetGenotype(l);//fine
			if (gt.Nalleles() == 0) continue;
			int m = model == 4 ? GetLoc(l).pes_model : model;
			double pr = gt.GFZ(m, grp->GetFreq(l));
			ChargeLog(re, re2, e1 * pr + e2 * gt.GFZ(m, total_pop->GetFreq(l)) * (1 - pr));
		}
		CloseLog(re, re2);

		if (loc >= 0 && re == 0) re = NA;
		return re;
	}

	/* Calculate the individual kinship coefficient */
	TARGET void IND::Theta(POP *grp, double &f_ritland, double &f_loiselle, double &f_weir, double &t_ritland, double &t_loiselle, double &t_weir, int64 loc)
	{
		double sw = 0, sr = 0, sr2 = 0, sw2 = 0, sr3 = 0, sw3 = 0;
		int ploidy = -1, N = 0;
		bool varploidy = false;
		int64 st = loc == (int64)-1 ? 0 : loc;
		int64 ed = loc == (int64)-1 ? nloc : loc + 1;

		for (int64 l = st; l < ed; ++l)
		{
			GENOTYPE &gt = GetGenotype(l);//fine
			if (gt.Nalleles() == 0) continue;
			double *p = grp->GetFreq(l);
			LOCSTAT &stat = grp->loc_stat[l];
			int k2 = stat.k;
			if (k2 <= 1) continue;

			int v = gt.Ploidy();
			if (ploidy == -1) ploidy = v;
			if (ploidy != v) varploidy = true;

			int k = GetLoc(l).k;
			double tsw2 = 0, tsw3 = 1, t1 = 0, t2 = 0, t3 = 0;

			for (int i = 0; i < k; ++i)
			{
				double pr = p[i];
				if (p[i] * stat.nhaplo <= 1e-5) continue;

				double af = gt.GetFreq(i);

				t1 += af * af / pr;

				t2 += (af - pr) * (af - pr);
				tsw2 += pr * (1 - pr);

				t3 += af * af - pr * pr;
				tsw3 -= pr * pr;
			}

			//Ritland
			sr += t1 - 1;
			sw += k2 - 1;

			//Loiselle
			sr2 += t2;
			sw2 += tsw2;

			//Weir
			sr3 += t3;
			sw3 += tsw3;

			double t_ritland1 = (t1 - 1) / (k2 - 1);
			double t_loiselle1 = t2 / tsw2;
			double t_weir1 = t3 / tsw3;

			t_ritland1 = (v * t_ritland1 - 1) / (v - 1);
			t_loiselle1 = (v * t_loiselle1 - 1) / (v - 1);
			t_weir1 = (v * t_weir1 - 1) / (v - 1);

			f_ritland += t_ritland1;
			f_loiselle += t_loiselle1;
			f_weir += t_weir1;

			N++;
		}
		t_ritland = sr / sw;
		t_loiselle = sr2 / sw2;
		t_weir = sr3 / sw3;

		if (!varploidy)
		{
			f_ritland = ploidy > 1 ? (ploidy * t_ritland - 1) / (ploidy - 1) : NA;
			f_loiselle = ploidy > 1 ? (ploidy * t_loiselle - 1) / (ploidy - 1) : NA;
			f_weir = ploidy > 1 ? (ploidy * t_weir - 1) / (ploidy - 1) : NA;
		}
		else
		{
			f_ritland = N > 0 ? f_ritland / N : NA;
			f_loiselle = N > 0 ? f_loiselle / N : NA;
			f_weir = N > 0 ? f_weir / N : NA;
		}
	}

	/* Write header row for individual statistics */
	TARGET /*static*/ void IND::IndividualStatisticsHeader(FILE *fout)
	{
		fprintf(fout, "%s%s%c%c%c%c%c", g_linebreak_val, g_linebreak_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val);
		int tn = (indstat_type_val[1]) +
				CountNonZero(indstat_ref_val, N_MAX_OPTION) *
				(indstat_type_val[2] * CountNonZero(indstat_model_val, N_MAX_OPTION) +
				(indstat_type_val[3] + indstat_type_val[4]) * CountNonZero(indstat_estimator_val, N_MAX_OPTION)) - 1;

		int64 tloc = indstat_locus_val[1] + indstat_locus_val[2] * nloc;

		if (indstat_locus_val[1])
		{
			for (int rl = 0; rl < lreg; ++rl)
				fprintf(fout, "%c", g_delimiter_val);
			fprintf(fout, "%cAll loci", g_delimiter_val);
			for (int i = 0; i < tn; ++i)
				fprintf(fout, "%c", g_delimiter_val);
		}

		if (indstat_locus_val[2]) for (int64 l = 0; l < nloc; ++l)
		{
			fprintf(fout, "%c%s", g_delimiter_val, GetLoc(l).GetName());
			for (int i = 0; i < tn; ++i)
				fprintf(fout, "%c", g_delimiter_val);
		}

		fprintf(fout, "%sInd%cPop", g_linebreak_val, g_delimiter_val);
		for (int rl = 0; rl < lreg; ++rl)
			fprintf(fout, "%cRegL%d", g_delimiter_val, rl);
		fprintf(fout, "%c#typed%c#miss%cPloidy%c#Hap", g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val);

		for (int64 l = 0; l < tloc; ++l)
		{
			if (indstat_type_val[1])
				fprintf(fout, "%cH-idx", g_delimiter_val);

			for (int i = 1; i <= N_DRE_MODEL; ++i)
			{
				if (indstat_model_val[i] == 0) continue;
				if (indstat_type_val[2])
				{
					if (indstat_ref_val[1]) fprintf(fout, "%clnL_pop_%s", g_delimiter_val, DRE_MODEL[i]);
					if (indstat_ref_val[2]) fprintf(fout, "%clnL_reg_%s", g_delimiter_val, DRE_MODEL[i]);
					if (indstat_ref_val[3]) fprintf(fout, "%clnL_tot_%s", g_delimiter_val, DRE_MODEL[i]);
				}
			}

			if (indstat_type_val[3])
			{
				if (indstat_ref_val[1])
				{
					if (indstat_estimator_val[1]) fprintf(fout, "%cF_pop_RI", g_delimiter_val);
					if (indstat_estimator_val[2]) fprintf(fout, "%cF_pop_LO", g_delimiter_val);
					if (indstat_estimator_val[3]) fprintf(fout, "%cF_pop_WE", g_delimiter_val);
				}
				if (indstat_ref_val[2])
				{
					if (indstat_estimator_val[1]) fprintf(fout, "%cF_reg_RI", g_delimiter_val);
					if (indstat_estimator_val[2]) fprintf(fout, "%cF_reg_LO", g_delimiter_val);
					if (indstat_estimator_val[3]) fprintf(fout, "%cF_reg_WE", g_delimiter_val);
				}
				if (indstat_ref_val[3])
				{
					if (indstat_estimator_val[1]) fprintf(fout, "%cF_tot_RI", g_delimiter_val);
					if (indstat_estimator_val[2]) fprintf(fout, "%cF_tot_LO", g_delimiter_val);
					if (indstat_estimator_val[3]) fprintf(fout, "%cF_tot_WE", g_delimiter_val);
				}
			}

			if (indstat_type_val[4])
			{
				if (indstat_ref_val[1])
				{
					if (indstat_estimator_val[1]) fprintf(fout, "%cTheta_pop_RI", g_delimiter_val);
					if (indstat_estimator_val[2]) fprintf(fout, "%cTheta_pop_LO", g_delimiter_val);
					if (indstat_estimator_val[3]) fprintf(fout, "%cTheta_pop_WE", g_delimiter_val);
				}
				if (indstat_ref_val[2])
				{
					if (indstat_estimator_val[1]) fprintf(fout, "%cTheta_reg_RI", g_delimiter_val);
					if (indstat_estimator_val[2]) fprintf(fout, "%cTheta_reg_LO", g_delimiter_val);
					if (indstat_estimator_val[3]) fprintf(fout, "%cTheta_reg_WE", g_delimiter_val);
				}
				if (indstat_ref_val[3])
				{
					if (indstat_estimator_val[1]) fprintf(fout, "%cTheta_tot_RI", g_delimiter_val);
					if (indstat_estimator_val[2]) fprintf(fout, "%cTheta_tot_LO", g_delimiter_val);
					if (indstat_estimator_val[3]) fprintf(fout, "%cTheta_tot_WE", g_delimiter_val);
				}
			}
		}
	}

	/* Write result row for individual statistics */
	TARGET void IND::PrintIndividualStatistics(FILE *fout)
	{
		int ntype = 0, minv = 99999, maxv = 0, nhaplo = 0, nmiss = 0;
		double hidx = 0;
		double f_tot_Ritland = 0, f_pop_Ritland = 0, f_reg_Ritland[N_MAX_REG];
		double f_tot_Loiselle = 0, f_pop_Loiselle = 0, f_reg_Loiselle[N_MAX_REG];
		double f_tot_Weir = 0, f_pop_Weir = 0, f_reg_Weir[N_MAX_REG];
		double t_tot_Ritland = 0, t_pop_Ritland = 0, t_reg_Ritland[N_MAX_REG];
		double t_tot_Loiselle = 0, t_pop_Loiselle = 0, t_reg_Loiselle[N_MAX_REG];
		double t_tot_Weir = 0, t_pop_Weir = 0, t_reg_Weir[N_MAX_REG];

		POP *ttot = total_pop;
		POP *tpop = apops[popid];
		POP *treg[N_MAX_REG]; //20200505
		if (lreg >= 1) treg[0] = aregs[0][tpop->rid]; //20200505

		for (int rl = 1; rl < lreg; ++rl)
			treg[rl] = aregs[rl][treg[rl - 1]->rid];

		for (int64 l = 0; l < nloc; ++l)
		{
			GENOTYPE &gt = GetGenotype(l);//fine
			if (gt.Nalleles() == 0) { nmiss++; continue; }
			ntype++;

			int v = gt.Ploidy();
			nhaplo += v;
			if (v > maxv) maxv = v;
			if (v < minv) minv = v;
			hidx += gt.HIndex();
		}
		hidx /= ntype;

		fprintf(fout, "%s%s%c%s%c",
			g_linebreak_val, name, g_delimiter_val,
			tpop->name, g_delimiter_val);

		for (int rl = 0; rl < lreg; ++rl)
			fprintf(fout, "%s%c", treg[rl]->name, g_delimiter_val);

		fprintf(fout, "%d%c%d%c%d-%d%c%d",
			ntype, g_delimiter_val,
			nmiss, g_delimiter_val,
			minv, maxv, g_delimiter_val,
			nhaplo);

		if (indstat_locus_val[1])
		{
			if (indstat_type_val[1]) {
				fprintf(fout, "%c", g_delimiter_val); WriteReal(fout, hidx);
			}

			for (int i = 1; i <= N_DRE_MODEL; ++i)
			{
				if (indstat_model_val[i] == 0) continue;
				if (indstat_type_val[2]) {
					if (indstat_ref_val[1]) { fprintf(fout, "%c", g_delimiter_val); WriteReal(fout, Likelihood(tpop, i, -1, 0)); }
					if (indstat_ref_val[2]) for (int rl = 0; rl < lreg; ++rl) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, Likelihood(treg[rl], i, -1, 0)); }
					if (indstat_ref_val[3]) { fprintf(fout, "%c", g_delimiter_val); WriteReal(fout, Likelihood(ttot, i, -1, 0)); }
				}
			}

			if (indstat_type_val[3] || indstat_type_val[4]) {
				if (indstat_ref_val[1])
					Theta(tpop, f_pop_Ritland, f_pop_Loiselle, f_pop_Weir, t_pop_Ritland, t_pop_Loiselle, t_pop_Weir);
				if (indstat_ref_val[2]) for (int rl = 0; rl < lreg; ++rl)
					Theta(treg[rl], f_reg_Ritland[rl], f_reg_Loiselle[rl], f_reg_Weir[rl], t_reg_Ritland[rl], t_reg_Loiselle[rl], t_reg_Weir[rl]);
				if (indstat_ref_val[3])
					Theta(ttot, f_tot_Ritland, f_tot_Loiselle, f_tot_Weir, t_tot_Ritland, t_tot_Loiselle, t_tot_Weir);
			}

			if (indstat_type_val[3]) {
				if (indstat_ref_val[1]) {
					if (indstat_estimator_val[1]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_pop_Ritland); }
					if (indstat_estimator_val[2]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_pop_Loiselle); }
					if (indstat_estimator_val[3]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_pop_Weir); }
				}
				if (indstat_ref_val[2]) for (int rl = 0; rl < lreg; ++rl) {
					if (indstat_estimator_val[1]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_reg_Ritland[rl]); }
					if (indstat_estimator_val[2]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_reg_Loiselle[rl]); }
					if (indstat_estimator_val[3]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_reg_Weir[rl]); }
				}
				if (indstat_ref_val[3]) {
					if (indstat_estimator_val[1]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_tot_Ritland); }
					if (indstat_estimator_val[2]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_tot_Loiselle); }
					if (indstat_estimator_val[3]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_tot_Weir); }
				}
			}

			if (indstat_type_val[4]) {
				if (indstat_ref_val[1]) {
					if (indstat_estimator_val[1]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_pop_Ritland); }
					if (indstat_estimator_val[2]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_pop_Loiselle); }
					if (indstat_estimator_val[3]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_pop_Weir); }
				}
				if (indstat_ref_val[2]) for (int rl = 0; rl < lreg; ++rl) {
					if (indstat_estimator_val[1]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_reg_Ritland[rl]); }
					if (indstat_estimator_val[2]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_reg_Loiselle[rl]); }
					if (indstat_estimator_val[3]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_reg_Weir[rl]); }
				}
				if (indstat_ref_val[3]) {
					if (indstat_estimator_val[1]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_tot_Ritland); }
					if (indstat_estimator_val[2]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_tot_Loiselle); }
					if (indstat_estimator_val[3]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_tot_Weir); }
				}
			}
		}

		////////////////////////////////////////////////////////////////////////////////////

		if (indstat_locus_val[2]) for (int64 l = 0; l < nloc; ++l)
		{
			if (indstat_type_val[1]) {
				fprintf(fout, "%c", g_delimiter_val); WriteReal(fout, GetGenotype(l).HIndex());//fine
			}

			for (int i = 1; i <= N_DRE_MODEL; ++i)
			{
				if (indstat_model_val[i] == 0) continue;
				if (indstat_type_val[2]) {
					if (indstat_ref_val[1]) { fprintf(fout, "%c", g_delimiter_val); WriteReal(fout, Likelihood(tpop, i, l, 0)); }
					if (indstat_ref_val[2]) for (int rl = 0; rl < lreg; ++rl) { fprintf(fout, "%c", g_delimiter_val); WriteReal(fout, Likelihood(treg[rl], i, l, 0)); }
					if (indstat_ref_val[3]) { fprintf(fout, "%c", g_delimiter_val); WriteReal(fout, Likelihood(ttot, i, l, 0)); }
				}
			}

			if (indstat_type_val[3] || indstat_type_val[4]) {
				if (indstat_ref_val[1]) Theta(tpop, f_pop_Ritland, f_pop_Loiselle, f_pop_Weir, t_pop_Ritland, t_pop_Loiselle, t_pop_Weir, l);
				if (indstat_ref_val[2]) for (int rl = 0; rl < lreg; ++rl) Theta(treg[rl], f_reg_Ritland[rl], f_reg_Loiselle[rl], f_reg_Weir[rl], t_reg_Ritland[rl], t_reg_Loiselle[rl], t_reg_Weir[rl], l);
				if (indstat_ref_val[3]) Theta(ttot, f_tot_Ritland, f_tot_Loiselle, f_tot_Weir, t_tot_Ritland, t_tot_Loiselle, t_tot_Weir, l);
			}

			if (indstat_type_val[3]) {
				if (indstat_ref_val[1]) {
					if (indstat_estimator_val[1]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_pop_Ritland); }
					if (indstat_estimator_val[2]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_pop_Loiselle); }
					if (indstat_estimator_val[3]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_pop_Weir); }
				}
				if (indstat_ref_val[2]) for (int rl = 0; rl < lreg; ++rl) {
					if (indstat_estimator_val[1]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_reg_Ritland[rl]); }
					if (indstat_estimator_val[2]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_reg_Loiselle[rl]); }
					if (indstat_estimator_val[3]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_reg_Weir[rl]); }
				}
				if (indstat_ref_val[3]) {
					if (indstat_estimator_val[1]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_tot_Ritland); }
					if (indstat_estimator_val[2]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_tot_Loiselle); }
					if (indstat_estimator_val[3]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, f_tot_Weir); }
				}
			}

			if (indstat_type_val[4]) {
				if (indstat_ref_val[1]) {
					if (indstat_estimator_val[1]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_pop_Ritland); }
					if (indstat_estimator_val[2]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_pop_Loiselle); }
					if (indstat_estimator_val[3]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_pop_Weir); }
				}
				if (indstat_ref_val[2]) for (int rl = 0; rl < lreg; ++rl) {
					if (indstat_estimator_val[1]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_reg_Ritland[rl]); }
					if (indstat_estimator_val[2]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_reg_Loiselle[rl]); }
					if (indstat_estimator_val[3]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_reg_Weir[rl]); }
				}
				if (indstat_ref_val[3]) {
					if (indstat_estimator_val[1]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_tot_Ritland); }
					if (indstat_estimator_val[2]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_tot_Loiselle); }
					if (indstat_estimator_val[3]) { fprintf(fout, "%c", g_delimiter_val);  WriteReal(fout, t_tot_Weir); }
				}
			}
		}
	}

	/* Write header row for ploidy inference, TEST */
	TARGET /*static*/ void IND::PloidyInferenceHeader(FILE *fout)
	{
		const char *ploidy_name[] = { "", "Haploid", "Diploid", "Triploid", "Tetraploid", "Pentaploid", "Hexaploid", "Heptaploid", "Octoploid", "Nonaploid", "Decaploid", "Undecaploid", "Dodecaploid" };

		fprintf(fout, "%s%s%c%c%c%c%c%c%c", g_linebreak_val, g_linebreak_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val);
		for (int rl = 0; rl < lreg; ++rl)
			fprintf(fout, "%c", g_delimiter_val);
		for (int i = 1; i <= N_MAX_PLOIDY; ++i)
			if (ploidyinfer_type_val[i])
				fprintf(fout, "%s%c%c%c", ploidy_name[i], g_delimiter_val, g_delimiter_val, g_delimiter_val);
		if (ploidyinfer_histogram_val)
			fprintf(fout, "Histogram Data");

		fprintf(fout, "%sInd%cPop", g_linebreak_val, g_delimiter_val);
		for (int rl = 0; rl < lreg; ++rl)
			fprintf(fout, "%cRegL%d", g_delimiter_val, rl);
		fprintf(fout, "%c#used%c#unused%c#miss%cEstPloidy%cPosterProb", g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val);
		for (int i = 1; i <= N_MAX_PLOIDY; ++i)
			if (ploidyinfer_type_val[i])
				fprintf(fout, "%clnL%cf%cPosterProb", g_delimiter_val, g_delimiter_val, g_delimiter_val);

		if (ploidyinfer_histogram_val)
		{
			double sep = 1.0 / ploidyinfer_nbins_val;
			for (int i = 0; i < ploidyinfer_nbins_val; ++i)
			{
				fprintf(fout, "%c", g_delimiter_val);
				WriteReal(fout, i * sep);
				fprintf(fout, "~");
				WriteReal(fout, (i + 1) * sep);
			}
		}
	}

	/* Write result row for ploidy inference, TEST */
	TARGET void IND::PrintPloidyInference(FILE *fout)
	{
		POP *tpop = apops[popid];
		POP *treg[N_MAX_REG] = { 0 }; //20200505
		if (lreg >= 1) treg[0] = aregs[0][tpop->rid]; //20200505

		for (int rl = 1; rl < lreg; ++rl)
			treg[rl] = aregs[rl][treg[rl - 1]->rid];

		int estploidy = 0;
		double sumprob = 0, max_lnL = -1e300;
		double lnL[N_MAX_PLOIDY + 1] = { 0 }, f[N_MAX_PLOIDY + 1] = { 0 }, poster_prob[N_MAX_PLOIDY + 1] = { 0 };
		int64 used = 0, unused = 0, miss = 0;

		//extract ad
		map<int64, SPF> depth;
		SPF tspf;
		tspf.count = 0;
		uint adb[16];
		uint *ad2 = adb;
		uint bins[101] = { 0 };
		char filename[4096];
		sprintf(filename, "used2-%s.txt", name);
		FILE *f1 = fopen(filename, "wb");
		fprintf(f1, "dA\tdB\r\n");

		for (int64 l = 0; l < nloc; ++l)
		{
			GENOTYPE &gt = GetGenotype(l);//fine

			if (GetLoc(l).k != 2)
			{ 
				unused++; 
				continue; 
			}

			if (gt.Nalleles() == 0) 
			{ 
				miss++; 
				continue; 
			}

			GetAlleleDepth(l, ad2);
			if (ad2[0] + ad2[1] == 0 ||
				(f_dp_b && (ad2[0] + ad2[1] < f_dp_min || ad2[0] + ad2[1] > f_dp_max)))
			{
				miss++; 
					continue;
			}

			if (ad2[0] == 0 || ad2[1] == 0) { unused++; continue; }
			int64 ha = (int64)ad2[0] << 32 | ad2[1];
			if (depth.find(ha) == depth.end()) depth[ha] = tspf;
			depth[ha].count++;
			used++;
			fprintf(f1, "%d\t%d\r\n", ad2[0], ad2[1]);
		}

		fclose(f1);

		for (int i = 1; i <= N_MAX_PLOIDY; ++i)
			if (ploidyinfer_type_val[i])
			{
				PloidyInference(i, lnL[i], f[i], depth);
				if (max_lnL < lnL[i])
				{
					max_lnL = lnL[i];
					estploidy = i;
				}
			}

		for (int i = 0; i <= N_MAX_PLOIDY; ++i)
			if (ploidyinfer_type_val[i])
				sumprob += exp(lnL[i] - max_lnL);

		for (int i = 0; i <= N_MAX_PLOIDY; ++i)
			if (ploidyinfer_type_val[i])
				poster_prob[i] = exp(lnL[i] - max_lnL) / sumprob;

		fprintf(fout, "%s%s%c%s", g_linebreak_val, name, g_delimiter_val, tpop->name);
		for (int rl = 0; rl < lreg; ++rl)
			fprintf(fout, "%c%s", g_delimiter_val, treg[rl]->name);

		fprintf(fout, "%c%lld%c%lld%c%lld%c%d%c", g_delimiter_val, used, g_delimiter_val, unused, g_delimiter_val, miss, g_delimiter_val, estploidy, g_delimiter_val);
		WriteReal(fout, poster_prob[estploidy]);

		for (int i = 0; i <= N_MAX_PLOIDY; ++i)
			if (ploidyinfer_type_val[i])
			{
				fprintf(fout, "%c", g_delimiter_val);
				WriteReal(fout, lnL[i]);
				fprintf(fout, "%c", g_delimiter_val);
				WriteReal(fout, f[i]);
				fprintf(fout, "%c", g_delimiter_val);
				WriteReal(fout, poster_prob[i]);
			}

		if (ploidyinfer_histogram_val)
		{
			for (auto d = depth.begin(); d != depth.end(); ++d)
			{
				int dA = d->first >> 32, dB = d->first & 0xFFFFFFFF;
				double fr = dA / (double)(dA + dB);
				int bid = (int)(fr * ploidyinfer_nbins_val);
				bins[bid] += d->second.count;
				bins[ploidyinfer_nbins_val - 1 - bid] += d->second.count;
			}

			for (int i = 0; i < ploidyinfer_nbins_val; ++i)
			{
				fprintf(fout, "%c", g_delimiter_val);
				fprintf(fout, "%d", bins[i]);
			}
		}

	}

	/* Calculate the likelihood for ploidy inference */
	TARGET void IND::PloidyInferlnL3(map<int64, SPF> &depth, int v, double f0, double f1, double f2, double &l0, double &l1, double &l2)
	{
		double re0 = 0, re1 = 0, re2 = 0;
		for (auto d = depth.begin(); d != depth.end(); ++d)
		{
			double *val = &(d->second.val1[0]);
			double *val2 = &(d->second.val2[0]);
			int dA = d->first >> 32, dB = d->first & 0xFFFFFFFF;
			double bicoef = LogBinomial(dA + dB, dA);
			double pr0A = 0, pr1A = 0, pr2A = 0;
			double pr0B = 1, pr1B = 1, pr2B = 1;
			for (int m = 0; m <= v; ++m)
			{
				double G = 0;
				G = AvgG(v, m, f0);
				pr0A += val[m] * G;
				pr0B -= val2[m] * G;
				G = AvgG(v, m, f1);
				pr1A += val[m] * G;
				pr1B -= val2[m] * G;
				G = AvgG(v, m, f2);
				pr2A += val[m] * G;
				pr2B -= val2[m] * G;
			}
			re0 += (log(pr0A / pr0B) + bicoef) * d->second.count;
			re1 += (log(pr1A / pr1B) + bicoef) * d->second.count;
			re2 += (log(pr2A / pr2B) + bicoef) * d->second.count;
		}
		l0 = re0;
		l1 = re1;
		l2 = re2;
	}

	/* Calculate the likelihood for ploidy inference */
	TARGET double IND::PloidyInferlnL(map<int64, SPF> &depth, int v, double f)
	{
		double re = 0;
		for (auto d = depth.begin(); d != depth.end(); ++d)
		{
			double *val = &(d->second.val1[0]);
			double *val2 = &(d->second.val2[0]);
			double prA = 0, prB = 1;
			int dA = d->first >> 32, dB = d->first & 0xFFFFFFFF;
			for (int m = 0; m <= v; ++m)
			{
				double G = AvgG(v, m, f);
				prA += val[m] * G;
				prB -= val2[m] * G;
			}
			re += (log(prA / prB) + LogBinomial(dA + dB, dA)) * d->second.count;
		}
		return re;
	}

	/* Infer the ploidy level from allelic depth distribution */
	TARGET void IND::PloidyInference(int v, double &lnL, double &f, map<int64, SPF> &depth)
	{
		for (auto d = depth.begin(); d != depth.end(); ++d)
		{
			int dA = d->first >> 32, dB = d->first & 0xFFFFFFFF, dAB = dA + dB;
			double *val = &(d->second.val1[0]);
			double *val2 = &(d->second.val2[0]);
			for (int m = 0; m <= v; ++m)
			{
				val[m] = ((m == 0 && dA) || (m == v && dB)) ? 0 : exp(LOGFRAC[m][v] * dA + LOGFRAC[v - m][v] * dB);
				val2[m] = (m == 0 ? 0 : exp(LOGFRAC[m][v] * dAB)) + (m == v ? 0 : exp(LOGFRAC[v - m][v] * dAB));
			}
		}

		//Newton's method
		double dx = 1e-5, x0 = 0.1, maxy0 = -1e300, maxx0 = 0, lastx0 = -1;

		for (int m = 0; m < MAX_ITER_SVD; ++m)
		{
			double y0 = 0, y1 = 0, y2 = 0;
			PloidyInferlnL3(depth, v, x0, x0 + dx, x0 + dx + dx, y0, y1, y2);

			if (y0 > maxy0)
			{
				maxx0 = x0;
				maxy0 = y0;
			}

			double yd0 = (y1 - y0) / dx;
			double yd1 = (y2 - y1) / dx;
			double ydd0 = (yd1 - yd0) / dx;

			lastx0 = x0;
			x0 -= yd0 / ydd0;
			if (x0 <= 0) x0 = 1e-10;
			if (x0 >= 1) x0 = 0.999999999;

			if (fabs(lastx0 - x0) < 1e-6) break;
		}

		lnL = maxy0;
		f = maxx0;
	}

	/* average genotypic frequency at a diallelic locus given m copies of A */
	TARGET double IND::AvgG(int v, int m, double f)
	{
		m = m > v - m ? v - m : m;

		switch (v)
		{
		case 1: return 0.5;
		case 2:
		{
			double div = 6;
			switch (m)
			{
			case 0: return (2 + f) / div;
			case 1: return (2 - 2 * f) / div;
			default: return -1;
			}
		}
		case 3:
		{
			double div = 4;
			switch (m)
			{
			case 0: return (1 + f) / div;
			case 1: return (1 - f) / div;
			default: return -1;
			}
		}
		case 4:
		{
			double div = 30 * (1 + f) * (1 + 2 * f);
			double f2 = f * f, f3 = f2 * f;
			switch (m)
			{
			case 0: return (6 + 27 * f + 38 * f2 + 19 * f3) / div;
			case 1: return (6 + 12 * f - 2 * f2 - 16 * f3) / div;
			case 2: return (6 + 12 * f - 12 * f2 - 6 * f3) / div;
			default: return -1;
			}
		}
		case 5:
		{
			double div = 12 * (1 + f) * (1 + 2 * f);
			double f2 = f * f, f3 = f2 * f;
			switch (m)
			{
			case 0: return (2 + 10 * f + 15 * f2 + 9 * f3) / div;
			case 1: return (2 + 4 * f + f2 - 7 * f3) / div;
			case 2: return (2 + 4 * f - 4 * f2 - 2 * f3) / div;
			default: return -1;
			}
		}
		case 6:
		{
			double div = 84 * (1 + f) * (1 + 2 * f) * (1 + 3 * f) * (1 + 4 * f);
			double f2 = f * f, f3 = f2 * f, f4 = f3 * f, f5 = f4 * f;
			switch (m)
			{
			case 0: return (12 + 150 * f + 708 * f2 + 1581 * f3 + 1726 * f4 + 863 * f5) / div;
			case 1: return (12 + 108 * f + 330 * f2 + 342 * f3 - 150 * f4 - 642 * f5) / div;
			case 2: return (12 + 108 * f + 288 * f2 + 153 * f3 - 402 * f4 - 159 * f5) / div;
			case 3: return (12 + 108 * f + 288 * f2 + 48 * f3 - 332 * f4 - 124 * f5) / div;
			default: return -1;
			}
		}
		case 7:
		{
			double div = 24 * (1 + f) * (1 + 2 * f) * (1 + 3 * f) * (1 + 4 * f);
			double f2 = f * f, f3 = f2 * f, f4 = f3 * f, f5 = f4 * f;
			switch (m)
			{
			case 0: return (3 + 39 * f + 190 * f2 + 438 * f3 + 495 * f4 + 275 * f5) / div;
			case 1: return (3 + 27 * f + 86 * f2 + 96 * f3 - 13 * f4 - 199 * f5) / div;
			case 2: return (3 + 27 * f + 72 * f2 + 54 * f3 - 111 * f4 - 45 * f5) / div;
			case 3: return (3 + 27 * f + 72 * f2 + 12 * f3 - 83 * f4 - 31 * f5) / div;
			default: return -1;
			}
		}
		case 8:
		{
			double div = 90 * (1 + f) * (1 + 2 * f) * (1 + 3 * f) * (1 + 4 * f) * (1 + 5 * f) * (1 + 6 * f);
			double f2 = f * f, f3 = f2 * f, f4 = f3 * f, f5 = f4 * f, f6 = f5 * f, f7 = f6 * f;
			switch (m)
			{
			case 0: return (10 + 245 * f + 2460 * f2 + 13075 * f3 + 39692 * f4 + 69459 * f5 + 67906 * f6 + 33953 * f7) / div;
			case 1: return (10 + 200 * f + 1590 * f2 + 6340 * f3 + 12854 * f4 + 10128 * f5 - 6998 * f6 - 24124 * f7) / div;
			case 2: return (10 + 200 * f + 1530 * f2 + 5590 * f3 + 9356 * f4 + 2622 * f5 - 14192 * f6 - 5116 * f7) / div;
			case 3: return (10 + 200 * f + 1530 * f2 + 5380 * f3 + 7718 * f4 - 1704 * f5 - 9866 * f6 - 3268 * f7) / div;
			case 4: return (10 + 200 * f + 1530 * f2 + 5380 * f3 + 6920 * f4 - 2250 * f5 - 8900 * f6 - 2890 * f7) / div;
			default: return -1;
			}
		}
		case 9:
		{
			double div = 20 * (1 + f) * (1 + 2 * f) * (1 + 3 * f) * (1 + 4 * f) * (1 + 5 * f) * (1 + 6 * f);
			double f2 = f * f, f3 = f2 * f, f4 = f3 * f, f5 = f4 * f, f6 = f5 * f, f7 = f6 * f;
			switch (m)
			{
			case 0: return (2 + 50 * f + 511 * f2 + 2761 * f3 + 8518 * f4 + 15178 * f5 + 15197 * f6 + 8183 * f7) / div;
			case 1: return (2 + 40 * f + 321 * f2 + 1301 * f3 + 2722 * f4 + 2316 * f5 - 961 * f6 - 5741 * f7) / div;
			case 2: return (2 + 40 * f + 306 * f2 + 1136 * f3 + 1966 * f4 + 864 * f5 - 3154 * f6 - 1160 * f7) / div;
			case 3: return (2 + 40 * f + 306 * f2 + 1076 * f3 + 1650 * f4 - 268 * f5 - 2102 * f6 - 704 * f7) / div;
			case 4: return (2 + 40 * f + 306 * f2 + 1076 * f3 + 1384 * f4 - 450 * f5 - 1780 * f6 - 578 * f7) / div;
			default: return -1;
			}
		}
		case 10:
		{
			double div = 132 * (1 + f) * (1 + 2 * f) * (1 + 3 * f) * (1 + 4 * f) * (1 + 5 * f) * (1 + 6 * f) * (1 + 7 * f) * (1 + 8 * f);
			double f2 = f * f, f3 = f2 * f, f4 = f3 * f, f5 = f4 * f, f6 = f5 * f, f7 = f6 * f, f8 = f7 * f, f9 = f8 * f;
			switch (m)
			{
			case 0: return (12 + 486 * f + 8440 * f2 + 82229 * f3 + 493806 * f4 + 1891753 * f5 + 4631876 * f6 + 7090179 * f7 + 6500866 * f8 + 3250433 * f9) / div;
			case 1: return (12 + 420 * f + 6218 * f2 + 50626 * f3 + 246174 * f4 + 721694 * f5 + 1192990 * f6 + 781206 * f7 - 739378 * f8 - 2259962 * f9) / div;
			case 2: return (12 + 420 * f + 6108 * f2 + 47931 * f3 + 219246 * f4 + 580443 * f5 + 776508 * f6 + 99393 * f7 - 1290522 * f8 - 439539 * f9) / div;
			case 3: return (12 + 420 * f + 6108 * f2 + 47436 * f3 + 210468 * f4 + 519492 * f5 + 567156 * f6 - 266940 * f7 - 827136 * f8 - 257016 * f9) / div;
			case 4: return (12 + 420 * f + 6108 * f2 + 47436 * f3 + 207960 * f4 + 489066 * f5 + 431064 * f6 - 312348 * f7 - 669000 * f8 - 200718 * f9) / div;
			case 5: return (12 + 420 * f + 6108 * f2 + 47436 * f3 + 207960 * f4 + 476592 * f5 + 393180 * f6 - 317892 * f7 - 627420 * f8 - 186396 * f9) / div;
			default: return -1;
			}
		}
		default: return -1;
		}
	}

	/* Write header row for population assignment */
	TARGET /*static*/ void IND::AssignmentHeader(FILE *fout)
	{
		fprintf(fout, "%s%sInd%cPop", g_linebreak_val, g_linebreak_val, g_delimiter_val);
		for (int rl = 0; rl < lreg; ++rl)
			fprintf(fout, "%cRegL%d", g_delimiter_val, rl + 1);
		fprintf(fout, "%c#typed%c#miss%cPloidy%c#Hap",
			g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val);
		for (int rl = -1; rl <= lreg; ++rl)
		{
			if (popas_level_val[rl == -1 ? 1 : 2] == 0) continue;
			POP **grps = rl == -1 ? apops : aregs[rl];
			int ngrp = rl == -1 ? npop : nreg[rl];

			//level: pop reg
			for (int l2 = 1; l2 <= N_DRE_MODEL; ++l2)
			{
				if (popas_model_val[l2] == 0) continue;
				fprintf(fout, "%cassign_%s", g_delimiter_val, rl == -1 ? "pop" : "reg");
				if (rl >= 0) fprintf(fout, "L%d", rl + 1);
				fprintf(fout, "_%s", DRE_MODEL[l2]);
				for (int id = 0; id < ngrp; ++id)
					fprintf(fout, "%clnL_%s_%s", g_delimiter_val, grps[id]->name, DRE_MODEL[l2]);
			}
		}
	}

	/* Write result row for population assignment */
	TARGET void IND::PrintAssignment(FILE *fout)
	{
		int64 ntype = 0, nhaplo = 0, nmiss = 0;
		int ms = npop;
		byte minv = 100, maxv = 0;

		VLA_NEW(lnL, double, ms);
		SetZero(lnL, ms);

		for (int64 l = 0; l < nloc; ++l)
		{
			GENOTYPE &gt = GetGenotype(l);//fine
			if (gt.Nalleles() == 0) 
			{ 
				nmiss++; 
				continue; 
			}

			int v = gt.Ploidy();
			ntype++;
			nhaplo += v;

			if (v > maxv) maxv = v;
			if (v < minv) minv = v;
		}

		fprintf(fout, "%s%s%c%s",
			g_linebreak_val, name,
			g_delimiter_val, apops[popid]->name);

		POP *tr = lreg >= 0 ? aregs[0][apops[popid]->rid] : NULL;
		for (int rl = 0; rl < lreg; ++rl)
		{
			fprintf(fout, "%c%s", g_delimiter_val, tr->name);
			tr = aregs[rl + 1][tr->rid];
		}

		fprintf(fout, "%c%lld%c%lld%c%d - %d%c%lld",
			g_delimiter_val, ntype,
			g_delimiter_val, nmiss,
			g_delimiter_val, minv, maxv,
			g_delimiter_val, nhaplo);

		for (int rl = -1; rl <= lreg; ++rl)
		{
			if (popas_level_val[rl == -1 ? 1 : 2] == 0) continue;
			POP **grps = rl == -1 ? apops : aregs[rl];
			int ngrp = rl == -1 ? npop : nreg[rl];

			for (int m2 = 1; m2 <= N_DRE_MODEL; ++m2)
			{
				if (popas_model_val[m2] == 0) continue;

				SetZero(lnL, ms);

				for (int j = 0; j < ngrp; ++j)
					lnL[j] = Likelihood(grps[j], m2, -1, popas_error_val);

				int64 mid = GetMaxID(lnL, ngrp);
				fprintf(fout, "%c%s", g_delimiter_val, grps[mid]->name);
				for (int j = 0; j < ngrp; ++j)
				{
					fprintf(fout, "%c", g_delimiter_val);
					WriteReal(fout, lnL[j]);
				}
			}
		}

		VLA_DELETE(lnL);
	}
	
	/* Read and set genotype from bcf input */
	TARGET void IND::AddBCFGenotype(int64 l, char *&gtstr, char *&gqstr, char *&dpstr, char *&adstr, int vlen, int asize, int gqlen, int dplen, int adlen, uint *&depth, TABLE<HASH, uint> &gfid, GENOTYPE *&gtab, ushort *&gatab, GENO_ITERATOR &iter)
	{
		LOCUS &loc = locus[l];
		TABLE<HASH, GENOTYPE*> &gftab = loc.gftab;
		ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };

		//convert hash
		char *gtnext = gtstr + asize * vlen, *gqnext = gqstr ? gqstr + gqlen : NULL, *dpnext = dpstr ? dpstr + dplen : NULL, *adnext = adstr ? adstr + adlen : NULL;
		bool phase = true;
		int v = 0;
		ReadBCFGenoString(gtstr, alleles, phase, v, asize, vlen, l, name);

		if (haplotype && !phase)
			loc.flag_pass = false;		  //exclude unphased genotype

		if (!haplotype) Sort(alleles, v); //phaseed, do not sort alleles

		HASH mha = missing_hash[v];
		uint mid = gfid[mha];
		if (!gftab.ContainsKey(mha))
			gftab[mha] = new(gtab++) GENOTYPE(gatab, missing_genotype[v]);
		
		HASH hash = HashGenotype(alleles, v);
		uint gid = gfid[hash];
		if (!gftab.ContainsKey(hash))
			gftab[hash] = new(gtab++) GENOTYPE(gatab, alleles, v);

		if (genotype_filter)
		{
			if (gid != mid && f_dp_b && dpstr && locus[l].dpid != 0xFFFF)
			{
				uint dp = ReadBinInteger(dpstr, dplen);
				if (dp != -1 && (dp < f_dp_min || dp > f_dp_max))
					gid = mid;
			}

			if (gid != mid && f_gq_b && gqstr && locus[l].gqid != 0xFFFF)
			{
				int gq2 = ReadBinInteger(gqstr, gqlen);
				if (gq2 != -1 && (gq2 < f_gq_min || gq2 > f_gq_max))
					gid = mid;
			}

			if (gid != mid && f_ploidy_b && (v < f_ploidy_min || v > f_ploidy_max))
				gid = mid;
		}

		if (ad == 2 && gid != mid && locus[l].adid != 0xFFFF)
		{
			double *fre = pop[popid].GetFreq(l);
			for (int j = 0; j < locus[l].k; ++j)
				fre[j] += ReadBinInteger(adstr, adlen);
		}
		else if (ploidyinfer && locus[l].adid != 0xFFFF)
		{
			uint maxdepth = 0;
			for (int j = 0; j < locus[l].k; ++j)
			{
				*depth = ReadBinInteger(adstr, adlen);
				maxdepth = Max(*depth++, maxdepth);
			}
			locus[l].maxdepth = maxdepth;
		}

		gtstr = gtnext; gqstr = gqnext; dpstr = dpnext; adstr = adnext;
		iter.Write(gid);
	}

	/* Read and set alleles from vcf input */
	TARGET void IND::AddVCFGenotype(char *&line, int64 l, uint *&depth, TABLE<HASH,uint> &gfid, GENOTYPE *&gtab, ushort *&gatab, GENO_ITERATOR &iter)
	{
		LOCUS &loc = locus[l];
		TABLE<HASH, GENOTYPE*> &gftab = loc.gftab;
		int fmtc = 0, gtid = loc.gtid;
		char *linebegin = line;
		for (;;)
		{
			fmtc++;
			while (*line != ':' && *line != '\t' && *line != '\n' && *line != '\0') line++;
			if (*line == ':')
				*line++ = '\0';
			else
			{
				*line++ = '\0';
				break;
			}
		}

		if (fmtc == 0)
			Exit("\nError: Format error in individual %s at locus %s.\n", name, loc.GetName());

		//convert hash
		char *genostr = GetTagValue(linebegin, gtid);
		int v = CountChars(genostr, "/|") + 1;
		if (v > N_MAX_PLOIDY) Exit("\nError: ploidy level greater than %d in individual %s at locus %s.\n", N_MAX_PLOIDY, name, loc.GetName());

		if (haplotype && CountChar(genostr, '/') > 0)
			loc.flag_pass = false;		//exclude unphased genotype

		ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
		ReadVCFGenoString(alleles, genostr, v, l, name);
		
		if (!haplotype) Sort(alleles, v); //phaseed, do not sort alleles

		HASH mha = missing_hash[v];
		uint mid = gfid[mha];
		if (!gftab.ContainsKey(mha))
			gftab[mha] = new(gtab++) GENOTYPE(gatab, missing_genotype[v]);

		HASH hash = HashGenotype(alleles, v);
		uint gid = gfid[hash];
		if (!gftab.ContainsKey(hash))
			gftab[hash] = new(gtab++) GENOTYPE(gatab, alleles, v);
		
		if (genotype_filter)
		{
			if (gid != mid && f_dp_b && locus[l].dpid != 0xFFFF && locus[l].dpid < fmtc)
			{
				uint dp = (uint)ReadIntegerKeep(GetTagValue(linebegin, locus[l].dpid));
				if (dp != -1 && (dp < f_dp_min || dp > f_dp_max))
					gid = mid;
			}

			if (gid != mid && f_gq_b && locus[l].gqid != 0xFFFF && locus[l].gqid < fmtc)
			{
				int gq = ReadIntegerKeep(GetTagValue(linebegin, locus[l].gqid));
				if (gq != -1 && (gq < f_gq_min || gq > f_gq_max))
					gid = mid;
			}

			if (gid != mid && f_ploidy_b && (v < f_ploidy_min || v > f_ploidy_max))
				gid = mid;
		}

		if (ad == 2 && gid != mid && locus[l].adid != 0xFFFF && locus[l].adid < fmtc)
		{
			double *fre = pop[popid].GetFreq(l);
			char *adval = GetTagValue(linebegin, locus[l].adid) - 1;
			for (int j = 0; j < locus[l].k; ++j)
				fre[j] += ReadInteger(++adval);
		}
		else if (ploidyinfer && locus[l].adid != 0xFFFF)
		{
			if (locus[l].adid < fmtc)
			{
				uint maxdepth = 0;
				char *adstr = GetTagValue(linebegin, locus[l].adid) - 1;
				for (int j = 0; j < locus[l].k; ++j)
				{
					*depth = ReadInteger(++adstr);
					maxdepth = Max(*depth++, maxdepth);
				}
				locus[l].maxdepth = maxdepth;
			}
			else SetZero(depth, locus[l].k);
		}

		iter.Write(gid);
	}

	/* Get tag value */
	TARGET char *IND::GetTagValue(char *re, int tagid)
	{
		if (tagid == -1) return 0;
		for (int i = 0; i < tagid; ++i)
			while (*re++);
		return re;
	}
#endif

#ifndef _DIVSUM

	/* Initialize sum */
	TARGET void DIVSUM::Init()
	{
		SetZero(this, 1);
		NE1P = 1;
		NE2P = 1;
		NEPP = 1;
		NEID = 1;
		NESI = 1;
		InitLock(lock);
	}

	/* Add loc to the sum */
	TARGET void DIVSUM::Add(DIVERSITY& loc)
	{
		Lock(lock);

		if (IsNormal(loc.NE1P)) NE1P *= loc.NE1P;
		if (IsNormal(loc.NE2P)) NE2P *= loc.NE2P;
		if (IsNormal(loc.NEPP)) NEPP *= loc.NEPP;
		if (IsNormal(loc.NEID)) NEID *= loc.NEID;
		if (IsNormal(loc.NESI)) NESI *= loc.NESI;

		ChargeSum(loc.k, k, kc);
		ChargeSum(loc.n, n, nc);
		ChargeSum(loc.nhaplo, nhaplo, nhaploc);
		ChargeSum(loc.ptype, ptype, ptypec);
		ChargeSum(loc.pval, pval, pvalc);
		ChargeSum(loc.he, he, hec);
		ChargeSum(loc.ho, ho, hoc);
		ChargeSum(loc.pic, pic, picc);
		ChargeSum(loc.ae, ae, aec);
		ChargeSum(loc.I, I, Ic);
		ChargeSum(loc.fis, fis, fisc);
		ChargeSum(loc.minploidy, minploidy, minploidyc);
		ChargeSum(loc.maxploidy, maxploidy, maxploidyc);

		UnLock(lock);
	}

	/* Add loc to the sum */
	TARGET void DIVSUM::Write(FILE *f, const char *name)
	{
		k /= kc;
		n /= nc;
		nhaplo /= nhaploc;
		ptype /= ptypec;
		pval /= pvalc;
		he /= hec;
		ho /= hoc;
		pic /= picc;
		ae /= aec;
		I /= Ic;
		fis /= fisc;
		minploidy /= minploidyc;
		maxploidy /= maxploidyc;

		fprintf(f, "%s%s%c(%lu)%c", g_linebreak_val, name, g_delimiter_val, nloc, g_delimiter_val);
		WriteReal(f, minploidy); fprintf(f, "-");
		WriteReal(f, maxploidy); fprintf(f, "%c", g_delimiter_val);
		WriteReal(f, k); fprintf(f, "%c", g_delimiter_val);
		WriteReal(f, n); fprintf(f, "%c", g_delimiter_val);
		WriteReal(f, nhaplo); fprintf(f, "%c", g_delimiter_val);
		WriteReal(f, ho); fprintf(f, "%c", g_delimiter_val);
		WriteReal(f, he); fprintf(f, "%c", g_delimiter_val);
		WriteReal(f, pic); fprintf(f, "%c", g_delimiter_val);
		WriteReal(f, ae); fprintf(f, "%c", g_delimiter_val);
		WriteReal(f, I); fprintf(f, "%c", g_delimiter_val);
		WriteReal(f, NE1P); fprintf(f, "%c", g_delimiter_val);
		WriteReal(f, NE2P); fprintf(f, "%c", g_delimiter_val);
		WriteReal(f, NEPP); fprintf(f, "%c", g_delimiter_val);
		WriteReal(f, NEID); fprintf(f, "%c", g_delimiter_val);
		WriteReal(f, NESI); fprintf(f, "%c", g_delimiter_val);
		WriteReal(f, fis); fprintf(f, "%c", g_delimiter_val);
		fprintf(f, "-%c-%c-", g_delimiter_val, g_delimiter_val);
	}

#endif
#ifndef _DIVERSITY

	/* Initialize */
	TARGET DIVERSITY::DIVERSITY()
	{
		SetZero(this, 1);
	}

	/* Calculate diveristy indices */
	TARGET void DIVERSITY::CalcDiversity(int64 _l)
	{
		l = _l;
		minploidy = 100;
		maxploidy = 0;

		int k2 = GetLoc(l).k, ploidy = 0, ni = cpop->nind;

		varploidy = false;
		n = 0;

		if (ni == 0)
		{
			nhaplo = 0;
			minploidy = 0;
			maxploidy = 0;
			bmaf = NA;
			k = 0;
			ptype = NA;
			pval = NA;
			he = NA;
			ho = NA;
			pic = NA;
			ae = NA;
			I = NA;
			NE1P = NA;
			NE2P = NA;
			NEPP = NA;
			NEID = NA;
			NESI = NA;
			fis = NA;
			return;
		}

		double *fre = cpop->GetFreq(l);
		ho = 0; nhaplo = 0;

		ushort *gcount = cpop->GetGenoCount(l);
		LOCSTAT &stat = cpop->loc_stat[l];

		//during diversity filter
		if (!reassigned)
		{
			GENO_ITERATOR iter(cpop->ind0id, l, true);//20220316
			for (int i = 0; i < nind; ++i)
				gcount[iter.Read()]++;
		}

		GENOTYPE *gtab = GetLoc(l).GetGtab();
		int ngeno = GetLoc(l).ngeno;

		v2i = 0;
		if (ad == 0) for (int gi = 0; gi < ngeno; ++gi)
		{
			GENOTYPE &gt = gtab[gi];
			if (gt.Nalleles() == 0 || gcount[gi] == 0) continue;
			uint c = gcount[gi];
			int v = gt.Ploidy();
			nhaplo += v * c;

			if (v > maxploidy) maxploidy = v;
			if (v < minploidy) minploidy = v;
			if (ploidy == 0) ploidy = v;
			if (maxploidy != minploidy) varploidy = true;

			ho += gt.HIndex() * c * v * (v - 1) / 2;
			n += c;
			v2i += c * v * (v - 1) / 2;

			if (!reassigned)
			{
				ushort *als = gt.GetAlleleArray();
				for (int i = 0; i < v; ++i)
					fre[als[i]] += c;
			}
		}
		else for (int gi = 0; gi < ngeno; ++gi)
		{
			GENOTYPE &gt = gtab[gi];
			if (gt.Nalleles() == 0 || gcount[gi] == 0) continue;
			uint c = gcount[gi];
			int v = gt.Ploidy();

			if (v > maxploidy) maxploidy = v;
			if (v < minploidy) minploidy = v;
			if (ploidy == 0) ploidy = v;
			if (maxploidy != minploidy) varploidy = true;

			ho += gt.HIndex() * c * v * (v - 1) / 2;
			n += c;
			v2i += c * v * (v - 1) / 2;

			if (!reassigned)
			{
				ushort *als = gt.GetAlleleArray();
				for (int i = 0; i < v; ++i)
					fre[als[i]] += c;
			}
		}

		if (maxploidy == 0) minploidy = maxploidy = 0;

		ho /= v2i;
		ptype = n / (double)ni;

		if (!reassigned)
			Unify(fre, k2);
		bmaf = 1;
		he = 1;
		pic = 0;
		ae = 0;
		I = 0;
		k = stat.k = (ushort)CountK(fre, k2);
		if (ad == 0) stat.nhaplo = nhaplo;
		else nhaplo = stat.nhaplo;
		int nh = nhaplo;

		double a2 = 0, a3 = 0, a4 = 0, a5 = 0, a6 = 0;
		for (int a = 0; a < k2; ++a)
		{
			double af = fre[a];
			if (af * nh < 1e-5) continue;
			if (af > 1e-5) I += -af * log(af);
			if (bmaf > af) bmaf = af;
			a2 += af * af;
			a3 += af * af * af;
			a4 += af * af * af * af;
			a5 += af * af * af * af * af;
			a6 += af * af * af * af * af * af;

			for (int b = a + 1; b < k2; ++b)
				if (fre[b] > 0)
					pic += 2 * af * fre[b] * (1 - af * fre[b]);
		}

		he -= a2;
		ae = 1 / a2;
		fis = 1 - ho / he;
		NE1P = 1 - (1 - 4 * a2 + 2 * a2 * a2 + 4 * a3 - 3 * a4);
		NE2P = 1 - (1 - 2 * a2 * a2 - 2 * a2 + a3 + 3 * a2 * a3 - 3 * a5 + 2 * a4);
		NEPP = 1 - (1 + 4 * a4 - 4 * a5 - 3 * a6 - 8 * a2 * a2 + 8 * a2 * a3 + 2 * a3 * a3);
		NESI = 1 - (0.75 - 0.5 * a2 - 0.5 * a2 * a2 + 0.25 * a4);
		NEID = 1 - (1 + 3 * a4 - 4 * a2 * a2);

		//chi test
		g = NA;
		df = NA;
		pval = NA;

		//do perform HWE test for locus diversity filter
		if (ad == 0 && !varploidy && !haplotype)
		{
			// initial mapping
			ploidy = maxploidy;
			MEMORY tlocus_memory;
			TABLE<HASH, uint> gitab(false, &tlocus_memory);
			VLA_NEW(gtmap2, ushort, maxG);

			//original
			LOCUS *locus1 = NULL;
			if (useslocus)
				locus1 = NEW(tlocus_memory, LOCUS, slocus[l]);
			else
				locus1 = &locus[l];

			VLA_NEW(fre1, double, k2);					SetVal(fre1, fre, k2);
			VLA_NEW(gcount1, ushort, ngeno);			SetVal(gcount1, gcount, ngeno);
			VLA_NEW(fre2, double, k2);
			VLA_NEW(gcount2, ushort, ngeno);
			VLA_NEW(amap, int, k2);

			ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };


			bool flag = false;
			for (int kk = k2; kk >= 1 && !flag; kk = (kk <= k ? kk - 1 : k))
			{
				if (kk == 1)
				{
					g = NA;
					df = NA;
					pval = NA;
					break;
				}

				TABLE<HASH, GENOTYPE*> &gftab1 = locus1->gftab;
				int ngeno1 = gftab1.size;
				SetZero(gcount2, ngeno1);
				SetZero(fre2, kk);

				if (kk <= k)
				{
					// test locus1
					flag = true;
					g = 0;

					for (int gi = 0; gi < ngeno1; ++gi)
					{
						GENOTYPE &gt = *gftab1(gi);
						if (gt.Nalleles() == 0 || gcount1[gi] == 0) continue;

						double obs = gcount1[gi];
						double exp = gt.GFZ(f_model_val, fre1) * n;
						if (exp < 5)
						{
							flag = false;
							break;
						}
						g += 2 * obs * log(obs / exp);
					}
				}

				if (!flag)
				{
					// fail in test, collapse
					if (kk <= k) 
					{
						// secondary map
						// find two minor alleles
						int mid1 = GetMinID(fre1, kk);
						fre1[mid1]++;
						int mid2 = GetMinID(fre1, kk);
						fre1[mid1]--;
						if (mid1 > mid2) Swap(mid1, mid2);

						// create allele map
						for (int j = 0, p = 0; j < kk; ++j)
						{
							if (j == mid2)
							{
								amap[j] = mid1;
								fre2[mid1] = fre1[mid1] + fre1[mid2];
							}
							else
							{
								amap[j] = p;
								fre2[p++] = fre1[j];
							}
						}
					}
					else
					{
						// initial map, exclude abscent alleles
						int mid1 = GetMaxID(fre1, kk);
						for (int j = 0, p = 0; j < kk; ++j)
							amap[j] = fre1[j] < MIN_FREQ ? -1 : p++;
						for (int j = 0; j < kk; ++j)
							amap[j] = amap[j] == -1 ? amap[mid1] : amap[j];
					}

					// create dummy locus
					TABLE<HASH, TEMP_GENOTYPE> tab2(true, &tlocus_memory);
					int gasize = 0;

					// create new alleles
					for (int gi = 0; gi < ngeno1; ++gi)
					{
						GENOTYPE *gt1 = gftab1(gi);
						if (gt1->Nalleles() == 0 || gcount1[gi] == 0) continue;

						ushort *als = gt1->GetAlleleArray();
						for (int ai = 0; ai < ploidy; ++ai)
							alleles[ai] = amap[als[ai]];

						Sort(alleles, ploidy);//unphase

						for (int j = 0; j < ploidy; ++j)
							fre2[alleles[j]] += gcount1[gi];

						HASH hash2 = HashGenotype(alleles, ploidy);
						if (!tab2.ContainsKey(hash2))
						{
							TEMP_GENOTYPE &gt2 = tab2[hash2];
							SetVal(gt2.alleles, alleles, N_MAX_PLOIDY);
							gt2.hash = hash2;
							gt2.ploidy = ploidy;
							gt2.gid = tab2.size - 1;
							gcount2[gt2.gid] += gcount1[gi];
							gasize += ploidy + GetNalleles(alleles, ploidy);
						}
						else
						{
							TEMP_GENOTYPE &gt2 = tab2[hash2];
							gcount2[gt2.gid] += gcount1[gi];
						}
					}


					// create dummy locus
					int ngeno2 = tab2.size;
					LOCUS *locus2 = new(locus_memory[threadid].Alloc(sizeof(LOCUS))) 
						LOCUS(locus_memory[threadid], *locus1, 1, ngeno2, gasize, tab2);//chi2 test

					// swap
					Unify(fre2, kk);
					locus1 = locus2;
					Swap(gcount1, gcount2);
					Swap(fre1, fre2);
				}
				else
				{
					// pass test
					df = (int)(Binomial(ploidy + kk - 1, ploidy) - kk + 0.5);
					pval = ChiSquareProb(g, df);
					break;
				}
			}


			VLA_DELETE(fre1);
			VLA_DELETE(fre2);
			VLA_DELETE(gcount1);
			VLA_DELETE(gcount2);
			VLA_DELETE(amap);
			VLA_DELETE(gtmap2);
		}
	}

	/* Write header to the result file */
	TARGET /*static*/ void DIVERSITY::WriteHeader(FILE *f)
	{
		fprintf(f, "%s%sPop", g_linebreak_val, g_linebreak_val);
		fprintf(f, "%cLocus", g_delimiter_val);
		fprintf(f, "%cPloidy", g_delimiter_val);
		fprintf(f, "%ck", g_delimiter_val);
		fprintf(f, "%cn", g_delimiter_val);
		fprintf(f, "%c#Hap", g_delimiter_val);
		fprintf(f, "%cHo", g_delimiter_val);
		fprintf(f, "%cHe", g_delimiter_val);
		fprintf(f, "%cPIC", g_delimiter_val);
		fprintf(f, "%cAe", g_delimiter_val);
		fprintf(f, "%cI", g_delimiter_val);
		fprintf(f, "%cNE1P", g_delimiter_val);
		fprintf(f, "%cNE2P", g_delimiter_val);
		fprintf(f, "%cNEPP", g_delimiter_val);
		fprintf(f, "%cNEID", g_delimiter_val);
		fprintf(f, "%cNESID", g_delimiter_val);
		fprintf(f, "%cFis", g_delimiter_val);
		fprintf(f, "%cG", g_delimiter_val);
		fprintf(f, "%cd.f.", g_delimiter_val);
		fprintf(f, "%cP-val", g_delimiter_val);
	}

	/* Write diversity of a locus to the result file */
	TARGET void DIVERSITY::WriteLocus(FILE *f, const char *name)
	{
		//sum pop, sum locus
		fprintf(f, "%s%s%c%s", g_linebreak_val, name, g_delimiter_val, GetLoc(l).GetName());
		fprintf(f, "%c%d-%d", g_delimiter_val, minploidy, maxploidy);
		fprintf(f, "%c%d", g_delimiter_val, k);
		fprintf(f, "%c%d", g_delimiter_val, n);
		fprintf(f, "%c%d%c", g_delimiter_val, nhaplo, g_delimiter_val);
		WriteReal(f, ho); fprintf(f, "%c", g_delimiter_val);
		WriteReal(f, he); fprintf(f, "%c", g_delimiter_val);
		WriteReal(f, pic); fprintf(f, "%c", g_delimiter_val);
		WriteReal(f, ae); fprintf(f, "%c", g_delimiter_val);
		WriteReal(f, I); fprintf(f, "%c", g_delimiter_val);
		WriteReal(f, NE1P); fprintf(f, "%c", g_delimiter_val);
		WriteReal(f, NE2P); fprintf(f, "%c", g_delimiter_val);
		WriteReal(f, NEPP); fprintf(f, "%c", g_delimiter_val);
		WriteReal(f, NEID); fprintf(f, "%c", g_delimiter_val);
		WriteReal(f, NESI); fprintf(f, "%c", g_delimiter_val);
		WriteReal(f, fis); fprintf(f, "%c", g_delimiter_val);
		WriteReal(f, g); fprintf(f, "%c", g_delimiter_val);
		if (IsError(df)) fprintf(f, "0%c", g_delimiter_val);
		else fprintf(f, "%0.0f%c", df, g_delimiter_val);
		WriteReal(f, pval);
	}

#endif

#ifndef _POP
	/* Initialize */
	TARGET POP::POP()
	{
		SetZero(this, 1);
	}

	/* Create a pop */
	TARGET POP::POP(char *_name, char **_names, int _nind, int _regid, int _npop, int _id, bool _ispop)
	{
		ispop = _ispop;
		id = _id;
		name = _name;
		names = _names;
		nind = _nind;
		rid = _regid;
		inds = NULL;
		loc_stat = NULL;
		allelefreq = NULL;
		genocount = NULL;

		npop = _npop;
		vpop = NULL;
	}

	/* Get genotype count array */
	TARGET ushort *POP::GetGenoCount(int64 l)
	{
		return genocount + genotype_count_offset[l];
	}

	/* Get genotype count of a genotype */
	TARGET ushort POP::GetGenoCount(int64 l, int gid)
	{
		return *(genocount + genotype_count_offset[l] + gid);
	}

	/* Get allele frequencies array */
	TARGET double *POP::GetFreq(int64 l)
	{
		return allelefreq + allele_freq_offset[l];
	}

	/* Get allele frequency of an allele */
	TARGET double POP::GetFreq(int64 l, int a)
	{
		return *(allelefreq + allele_freq_offset[l] + a);
	}

	/* Uninitialize a pop */
	TARGET void POP::Uninitialize()
	{
		if (allelefreq) 
		{
			delete[] allelefreq; 
			allelefreq = NULL;
		}

		if (genocount)
		{
			delete[] genocount;
			genocount = NULL;
		}

		if (loc_stat) 
		{
			delete[] loc_stat; 
			loc_stat = NULL;
		}
		
		if (names) 
		{ 
			delete[] names;
			names = NULL;
		}
	}

	/* Uncllocate memory for locstat, allele frequency, genotype count */
	TARGET void POP::UnAllocFreq()
	{
		if (loc_stat) 
		{ 
			delete[] loc_stat; 
			loc_stat = NULL;
		}

		if (allelefreq) 
		{ 
			delete[] allelefreq; 
			allelefreq = NULL; 
		}

		if (genocount) 
		{ 
			delete[] genocount; 
			genocount = NULL; 
		}
	}

	/* Allocate memory for locstat, allele frequency, genotype count */
	TARGET void POP::AllocFreq()
	{
		//to apply diversity filter and test genotype distributions
		loc_stat = new LOCSTAT[nloc];  SetZero(loc_stat, nloc);
		allelefreq = new double[KT];   SetZero(allelefreq, KT);
		genocount = new ushort[GT];    SetZero(genocount, GT);
	}

	/* Move after filter locus */
	TARGET void POP::MoveFreq(LOCN *nafoffset, LOCN *ngcoffset)
	{
		LOCSTAT *nloc_stat = NULL;
		if (loc_stat)
		{
			new LOCSTAT[nfilter];
			SetZero(nloc_stat, nfilter);
		}

		double *nallelefreq = NULL;
		if (allelefreq)
		{
			nallelefreq = new double[KT];
			SetZero(nallelefreq, KT);
		}

		ushort *ngenocount = NULL;
		if (genocount)
		{
			ngenocount = new ushort[GT];
			SetZero(ngenocount, GT);
		}

		//move freq and count to another compat memory
		for (int64 l = 0, nl = 0; l < nloc; ++l)
		{
			if (GetLoc(l).flag_pass)
			{
				if (loc_stat)
					nloc_stat[nl] = loc_stat[l];
				if (allelefreq)
					SetVal(nallelefreq + nafoffset[nl], GetFreq(l), GetLoc(l).k);
				if (genocount)
					SetVal(ngenocount + ngcoffset[nl], GetGenoCount(l), GetLoc(l).ngeno);
				nl++;
			}
		}

		if (loc_stat)
		{
			delete[] allelefreq; 
			allelefreq = nallelefreq;
		}

		if (allelefreq)
		{
			delete[] genocount;  
			genocount = ngenocount;
		}

		if (genocount)
		{
			delete[] loc_stat;   
			loc_stat = nloc_stat;
		}
	}

	/* Clear memory for locstat, allele frequency, genotype count */
	TARGET void POP::ClearFreqGcount()
	{
		SetZero(loc_stat, nloc);
		SetZero(allelefreq, KT);
		SetZero(genocount, GT);
	}

	/* Calculate loc stat, allele frequency, genotype count */
	TARGET void POP::CalcFreqGcount()
	{
		if (nind == 0) return;

		if (ad == 0)
		{
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
			for (int64 l = 0; l < nloc; ++l)
			{
				GENOTYPE *gtab = GetLoc(l).GetGtab();
				int ngeno = GetLoc(l).ngeno;

				double *p = GetFreq(l);
				ushort *gcount = GetGenoCount(l);
				LOCSTAT &stat = loc_stat[l];

				GENO_ITERATOR iter(ind0id, l, true);//20220316
				for (int i = 0; i < nind; ++i)
				{
					int gtid = iter.Read();
					gcount[gtid]++;
				}

				int nhaplo = 0;
				for (int gi = 0; gi < ngeno; ++gi)
				{
					int gc = gcount[gi];
					if (gc == 0) continue;

					GENOTYPE &gt = gtab[gi]; 
					if (gt.Nalleles() == 0) continue;

					ushort *als = gt.GetAlleleArray();
					int v = gt.Ploidy();
					nhaplo += v * gc;

					for (int j = 0; j < v; ++j)
						p[als[j]] += gc;
				}

				int k2 = GetLoc(l).k;
				stat.k = (ushort)CountK(p, k2);
				stat.nhaplo = nhaplo;
				Unify(p, k2);

				stat.a2 = stat.a3 = stat.a4 = 0;
				for (int a = 0; a < k2; ++a)
				{
					double af = p[a];
					if (af * nhaplo <= 1e-5) continue;
					stat.a2 += af * af;
					stat.a3 += af * af * af;
					stat.a4 += af * af * af * af;
				}
				PROGRESS_VALUE += nind;
			}
		}
		else
		{
			loc_stat = new LOCSTAT[nloc];
			genocount = new ushort[GT];
			SetZero(genocount, GT);

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
			for (int64 l = 0; l < nloc; ++l)
			{
				double *p = GetFreq(l);
				ushort *gcount = GetGenoCount(l);
				LOCSTAT &stat = loc_stat[l];

				GENO_ITERATOR iter(ind0id, l, true);//20220316
				for (int i = 0; i < nind; ++i)
					gcount[iter.Read()]++; 

				int k2 = GetLoc(l).k;
				stat.k = (ushort)CountK(p, k2);
				int nhaplo = stat.nhaplo = Round(Sum(p, k2));
				Unify(p, k2);

				stat.a2 = stat.a3 = stat.a4 = 0;
				for (int a = 0; a < k2; ++a)
				{
					double af = p[a];
					if (af * nhaplo <= 1e-5) continue;
					stat.a2 += af * af;
					stat.a3 += af * af * af;
					stat.a4 += af * af * af * af;
				}
				PROGRESS_VALUE += nind;
			}
		}
	}
#endif

#ifndef _FST
	/* Fst estimator warpper */
	TARGET /*static*/ double FST::FstEstimator(POP **grps, int n, int e, double *each, double *buf)
	{
		switch (e)
		{
		case 1: return Fst_Nei1973(grps, n, each, buf);
		case 2: return Fst_Weir1984(grps, n, each);
		case 3: return Fst_Hudson1992(grps, n, each, buf);
		case 4: return Fst_Slatkin1995(grps, n, each, buf);
		case 5: return Fst_Hedrick2005(grps, n, each, buf);
		case 6: return Fst_Jost2008(grps, n, each, buf);
		case 7: return Fst_Huang2021_homo(grps, n, 2, true, each);//homoploid
		case 8: return Fst_Huang2021_aniso(grps, n, 2, true, true, each, buf);//anisoploid
		}
		return 0;
	}

	/* Estimate Fst and test differentiation for two pops */
	TARGET void FST::CalcFst(POP *a, POP *b)
	{
		POP *c[] = { a, b };
		CalcFst(c, 2);
	}

	/* Estimate Fst and test differentiation for multiple pops */
	TARGET void FST::CalcFst(POP **grps, int n)
	{
		Genotype_G = Genotype_P = Allele_G = Allele_P = 0; Genotype_DF = Allele_DF = 0;//Difference
		Nei1973 = Weir1984 = Hudson1992 = Slatkin1995 = Hedrick2005 = Jost2008 = Huang2021_homo = Huang2021_aniso = 0;

		VLA_NEW(buf, double, maxK);
		byte *estimator = fst_estimator_val;
		SetZero(&Nei1973, N_FST_ESTIMATOR);

		for (int i = 1; i <= N_FST_ESTIMATOR; ++i)
		{
			if (estimator[i] == 0) continue;
			*(&Nei1973 + i - 1) = fst_locus_val[2] ? new double[nloc] : NULL;
			*(&Nei1973T + i - 1) = FstEstimator(grps, n, i, *(&Nei1973 + i - 1), buf);
		}

		//test Fst
		if (fst_locus_val[1] == 0 && fst_locus_val[2] == 0)
			Genotype_PT = Allele_PT = NA;

		if ((fst_locus_val[1] || fst_locus_val[2]) && fst_test_val[1])
		{
			Genotype_DFT = 0; Genotype_GT = 0;
			Genotype_G = new double[nloc];
			Genotype_DF = new int[nloc];
			Genotype_P = new double[nloc];

			int64 gcsize = maxG + 2;
			VLA_NEW(obs, double, gcsize * (4 * n + 1));
			double *obs2 = obs + gcsize * n;
			double *exp = obs2 + gcsize * n;
			double *exp2 = exp + gcsize * n;
			double *colsum = exp2 + gcsize * n;
			VLA_NEW(rowsum, double, n);

			for (int64 l = 0; l < nloc; ++l)
			{
				bool sameploidy = true;
				int ploidy = -1;
				GENOTYPE *gtab = GetLoc(l).GetGtab();
				int ngeno = GetLoc(l).ngeno;

				SetZero(obs, ngeno * n);
				if (ngeno >= 2 && n >= 2)
				{
					for (int gi = 0; gi < ngeno; ++gi)
					{
						GENOTYPE &gt = gtab[gi];
						if (gt.Nalleles() == 0) continue;

						int v = gt.Ploidy();
						if (ploidy == -1) ploidy = v;
						if (ploidy != v)
						{
							sameploidy = false;
							break;
						}

						for (int p = 0; p < n; ++p)
							obs[p * ngeno + gi] = grps[p]->GetGenoCount(l, gi);
					}
					if (sameploidy)
					{
						CombineTable(obs, n, ngeno, Genotype_G[l], Genotype_DF[l], Genotype_P[l], (bool)fst_locus_val[2], obs2, exp, exp2, rowsum, colsum);
						Genotype_DFT += Genotype_DF[l];
						Genotype_GT += Genotype_G[l];
					}
				}
				if (ngeno < 2 || n < 2 || !sameploidy)
				{
					Genotype_DF[l] = 0;
					Genotype_G[l] = 0;
					Genotype_P[l] = NA;
				}
			}
			VLA_DELETE(obs);
			VLA_DELETE(rowsum);
			Genotype_PT = Genotype_DFT > 0 ? ChiSquareProb(Genotype_GT, Genotype_DFT) : NA;
		}

		if ((fst_locus_val[1] || fst_locus_val[2]) && fst_test_val[2])
		{
			Allele_DFT = 0; Allele_GT = 0;
			Allele_G = new double[nloc];
			Allele_DF = new int[nloc];
			Allele_P = new double[nloc];
			VLA_NEW(obs, double, maxK * n * 4 + maxK + n);
			double *obs2 = obs + maxK * n;
			double *exp = obs2 + maxK * n;
			double *exp2 = exp + maxK * n;
			double *colsum = exp2 + maxK * n;
			double *rowsum = colsum + maxK;

			for (int64 l = 0; l < nloc; ++l)
			{
				int m = grps[0]->loc_stat[l].k;

				if (m >= 2 && n >= 2)
				{
					SetZero(obs, m * n);
					for (int p = 0; p < n; ++p)
					{
						double *freq = grps[p]->GetFreq(l);
						int nhaplo = grps[p]->loc_stat[l].nhaplo;
						for (int i = 0; i < m; ++i)
							obs[p * m + i] = nhaplo * freq[i];
					}
					CombineTable(obs, n, m, Allele_G[l], Allele_DF[l], Allele_P[l], (bool)fst_locus_val[2], obs2, exp, exp2, rowsum, colsum);
					Allele_DFT += Allele_DF[l];
					Allele_GT += Allele_G[l];
				}
				else
				{
					Allele_DF[l] = 0; Allele_G[l] = 0;
					Allele_P[l] = NA;
				}
			}
			Allele_PT = Allele_DFT > 0 ? ChiSquareProb(Allele_GT, Allele_DFT) : NA;
			VLA_DELETE(obs);
		}

		PROGRESS_VALUE++;

		VLA_DELETE(buf);
		return;
	}

	/* Uninitialize */
	TARGET void FST::Uninitialize()
	{
		if (Nei1973)		delete[] Nei1973;		 Nei1973 = NULL;
		if (Weir1984)		delete[] Weir1984;		 Weir1984 = NULL;
		if (Hudson1992)		delete[] Hudson1992;	 Hudson1992 = NULL;
		if (Slatkin1995)	delete[] Slatkin1995;	 Slatkin1995 = NULL;
		if (Hedrick2005)	delete[] Hedrick2005;	 Hedrick2005 = NULL;
		if (Jost2008)		delete[] Jost2008;		 Jost2008 = NULL;
		if (Huang2021_homo)	delete[] Huang2021_homo; Huang2021_homo = NULL;
		if (Huang2021_aniso)delete[] Huang2021_aniso; Huang2021_aniso = NULL;

		if (Genotype_G)		delete[] Genotype_G;	Genotype_G = NULL;
		if (Genotype_DF)	delete[] Genotype_DF;	Genotype_DF = NULL;
		if (Genotype_P)		delete[] Genotype_P;	Genotype_P = NULL;
		if (Allele_G)		delete[] Allele_G;		Allele_G = NULL;
		if (Allele_DF)		delete[] Allele_DF;		Allele_DF = NULL;
		if (Allele_P)		delete[] Allele_P;		Allele_P = NULL;
	}

	/* Nei 1973 Gst estimator based on heterozgysotiy */
	TARGET /*static*/ double FST::Fst_Nei1973(POP **grps, int n, double *each, double *buf)
	{
		//weight by pop size chakraborty 1974, allow ad
		double Dst = 0, Ht = 0;
		double *pi = buf;
		for (int64 l = 0; l < nloc; ++l)
		{
			if (each) each[l] = NA;
			int k2 = GetLoc(l).k;
			if (k2 < 2) continue;

			double ht = 0, hs = 0, sw = 0;
			int s = 0;
			SetZero(pi, k2);
			for (int i = 0; i < n; ++i)
			{
				int nhaplo = grps[i]->loc_stat[l].nhaplo;
				double *p = grps[i]->GetFreq(l);
				if (nhaplo == 0) continue;
				s++;
				AddProd(pi, p, nhaplo, k2);
				sw += nhaplo;
				hs += nhaplo * SumSquare(p, k2);
			}

			if (s < 2) continue;
			Unify(pi, k2);
			hs = 1 - hs / sw;
			ht = 1 - SumSquare(pi, k2);
			Dst += ht - hs;
			Ht += ht;
			if (each) each[l] = IsError((ht - hs) / ht) ? NA : (ht - hs) / ht;
		}
		return Dst / Ht;
	}

	/* Weir 1984 Fst estimator based on variance components */
	TARGET /*static*/ double FST::Fst_Weir1984(POP **grps, int n, double *each)
	{
		//for diploid only, allow ad
		if (minploidy != 2 || maxploidy != 2)
			Exit("\nError: Weir1984 Fst estimator can only be applied for diploids\n");

		double theta_W1 = 0, theta_W2 = 0;
		for (int64 l = 0; l < nloc; ++l)
		{
			if (each) each[l] = NA;
			int k2 = GetLoc(l).k;
			if (k2 < 2) continue;

			double nt = 0, nb = 0, nc = 0;
			double r = 0, ni2 = 0;
			for (int i = 0; i < n; ++i)
			{
				int nhaplo = grps[i]->loc_stat[l].nhaplo;
				if (nhaplo == 0) continue;
				r++;
				int ni = nhaplo >> 1;
				nt += ni;
				ni2 += ni * ni;
			}
			if (r < 2 || nt == 0) continue;

			nb = nt / r;
			nc = (r * nb - ni2 / nt) / (r - 1);

			double w1 = 0, w2 = 0;
			for (int a = 0; a < k2; ++a)
			{
				double p1 = 0, p2 = 0, s2 = 0, hb = 0;
				for (int i = 0; i < n; ++i)
				{
					int nhaplo = grps[i]->loc_stat[l].nhaplo;
					if (nhaplo == 0) continue;
					double p = grps[i]->GetFreq(l, a);
					int ni = nhaplo >> 1;
					p1 += p * ni;
					p2 += p * p * ni;
				}
				p1 /= nt;
				p2 /= nt;
				s2 = (p2 - p1 * p1) * r / (r - 1);
				hb = 2 * (p1 - p2);

				if (nb <= 1 || hb <= 0) break;

				double va = nb / nc * (s2 - 1 / (nb - 1) * (p1 * (1 - p1) - (r - 1) / r * s2 - 0.25 * hb));
				double vb = nb / (nb - 1) * (p1 * (1 - p1) - (r - 1) / r * s2 - (2 * nb - 1) / (4 * nb) * hb);
				double vc = 0.5 * hb;

				if (IsError(va) || IsError(vb) || IsError(vc)) break;

				w1 += va;
				w2 += va + vb + vc;
			}
			theta_W1 += w1;
			theta_W2 += w2;
			if (each) each[l] = IsError(w1 / w2) ? NA : w1 / w2;
		}

		return theta_W1 / theta_W2;
	}

	/* Nei 1973 Fst estimator based on mean allele difference */
	TARGET /*static*/ double FST::Fst_Hudson1992(POP **grps, int n, double *each, double *buf)
	{
		//do not need to weight, allow ad
		double *Ni = buf;
		double f1 = 0, f2 = 0;
		for (int64 l = 0; l < nloc; ++l)
		{
			if (each) each[l] = NA;
			int k2 = GetLoc(l).k;
			if (k2 < 2) continue;

			int Nh = 0, Np = 0;
			double Dw = 0, Nw = 0;
			SetZero(Ni, k2);
			for (int i = 0; i < n; ++i)
			{
				int nhaplo = grps[i]->loc_stat[l].nhaplo;
				if (nhaplo == 0) continue;
				double *p = grps[i]->GetFreq(l);
				AddProd(Ni, p, nhaplo, k2);

				Nh += nhaplo;
				Np++;

				int nh = nhaplo;
				Dw += nh * (nh - 1); //2x
				Nw += nh * (nh - 1); //2x
				for (int a = 0; a < k2; ++a)
				{
					double na = nh * p[a];
					if (na > 1) Dw -= na * (na - 1); //substrate same
				}
			}

			if (Np < 2) continue;

			double hw = Dw / Nw; //mean diff
			double Db = (double)(Nh * (Nh - 1)), Nb = Db;
			for (int a = 0; a < k2; ++a)
				Db -= Ni[a] * (Ni[a] - 1);
			double hb = (Db - Dw) / (Nb - Nw);

			if (each) each[l] = IsError(1 - hw / hb) ? NA : 1 - hw / hb;
			if (IsNormal(hw) && IsNormal(hb))
			{
				f1 += hb - hw;
				f2 += hb;
			}
		}
		return f1 / f2;
	}

	/* Slatkin 1995 Fst estimator based on allele size */
	TARGET /*static*/ double FST::Fst_Slatkin1995(POP **grps, int n, double *each, double *buf)
	{
		if (abs(g_format_val) <= 2)
			Exit("\nError: Slatkin1995 Fst estimator can only be applied for non-vcf input file, and should use size as allele identifier. \n");

		//do not need to weight, allow ad
		double *nt = buf;
		double f1 = 0, f2 = 0;
		for (int64 l = 0; l < nloc; ++l)
		{
			double Sw1 = 0, Sw2 = 0, St1 = 0, St2 = 0;
			if (each) each[l] = NA;
			int k2 = GetLoc(l).k;
			if (k2 < 2) continue;

			int Nh = 0, Np = 0;
			SetZero(nt, k2);
			ushort *alen2 = GetLoc(l).GetAlenArray() + k2;

			for (int i = 0; i < n; ++i)
			{
				int nhaplo = grps[i]->loc_stat[l].nhaplo;
				if (nhaplo == 0) continue;
				double *p = grps[i]->GetFreq(l);
				AddProd(nt, p, nhaplo, k2);

				Nh += nhaplo;
				Np++;
				Sw2 += nhaplo * (nhaplo - 1) * 0.5;

				for (int ai = 0; ai < k2; ++ai)
				{
					double ni = nhaplo * p[ai];
					for (int aj = 0; aj < ai; ++aj)
					{
						double nj = nhaplo * p[aj];
						Sw1 += ni * nj * alen2[ai * k2 + aj];
					}
				}
			}

			if (Np < 2) continue;
			St2 = Nh * (Nh - 1) * 0.5;
			for (int i = 0; i < k2; ++i)
				for (int j = 0; j < i; ++j)
					St1 += nt[i] * nt[j] * alen2[i * k2 + j];

			double Sw = Sw1 / Sw2, St = St1 / St2;
			if (IsNormal(Sw) && IsNormal(St))
			{
				if (each) each[l] = IsError((St - Sw) / St) ? NA : (St - Sw) / St;
				f1 += St - Sw;
				f2 += St;
			}
		}
		return f1 / f2;
	}

	/* Hedrick 2005 G'st */
	TARGET /*static*/ double FST::Fst_Hedrick2005(POP **grps, int n, double *each, double *buf)
	{
		//weight by pop size chakraborty 1974, allow ad
		double DA = 0, DB = 0;
		double *pi = buf;
		for (int64 l = 0; l < nloc; ++l)
		{
			if (each) each[l] = NA;
			int k2 = GetLoc(l).k;
			if (k2 < 2) continue;

			double ht = 0, hs = 0, sw = 0;
			int s = 0;
			SetZero(pi, k2);
			for (int i = 0; i < n; ++i)
			{
				int nhaplo = grps[i]->loc_stat[l].nhaplo;
				if (nhaplo == 0) continue;
				double *p = grps[i]->GetFreq(l);
				s++;
				sw += nhaplo;
				AddProd(pi, p, nhaplo, k2);
				hs += nhaplo * SumSquare(p, k2);
			}
			if (s < 2) continue;
			Unify(pi, k2);
			hs = 1 - hs / sw;
			ht = 1 - SumSquare(pi, k2);
			double da = (ht - hs) * (s - 1 + hs);
			double db = ht * (s - 1) * (1 - hs);
			DA += da;
			DB += db;
			if (each) each[l] = IsError(da / db) ? NA : da / db;
		}
		return DA / DB;
	}

	/* Jost 2008 D */
	TARGET /*static*/ double FST::Fst_Jost2008(POP **grps, int n, double *each, double *buf)
	{
		//weight by pop size chakraborty 1974, allow ad
		double DA = 0, DB = 0;
		double *pi = buf;
		for (int64 l = 0; l < nloc; ++l)
		{
			if (each) each[l] = NA;
			int k2 = GetLoc(l).k;
			if (k2 < 2) continue;

			double ht = 0, hs = 0, sw = 0;
			int s = 0;
			SetZero(pi, k2);
			for (int i = 0; i < n; ++i)
			{
				int nhaplo = grps[i]->loc_stat[l].nhaplo;
				if (nhaplo == 0) continue;
				double *p = grps[i]->GetFreq(l);
				s++;
				sw += nhaplo;
				AddProd(pi, p, nhaplo, k2);
				hs += nhaplo * SumSquare(p, k2);
			}
			if (s < 2) continue;
			Unify(pi, k2);
			hs = 1 - hs / sw;
			ht = 1 - SumSquare(pi, k2);
			double da = (ht - hs) * s;
			double db = ht * (s - 1);
			DA += da;
			DB += db;
			if (each) each[l] = IsError(da / db) ? NA : da / db;
		}
		return DA / DB;
	}

	/* Huang 2021 Fst estimator based on multi-level amova */
	TARGET /*static*/ double FST::Fst_Huang2021_homo(POP **grps, int n, int layer, bool isiam, double *each)
	{
		//2 or 3 Level, homoploids, forbid ad
		if (ad && layer == 3) Exit("\nError: Huang2021 Fst estimator (-fst_estimator=Huang2021_homo) is incompatible with allelic depth (-ad) option.\n");
		double svi2dvp = 0, svi2dvt = 0, svp2dvt = 0;
		int Ni = 0, Np = n, Nh = 0;
		VLA_NEW(pop_nhaplo, int, n);
		SetZero(pop_nhaplo, n);

		for (int p = 0; p < n; ++p)
		{
			Ni += grps[p]->nind;
			IND **vind = grps[p]->inds;
			for (int i = 0; i < grps[p]->nind; ++i)
			{
				int v1 = vind[i]->vmin, v2 = vind[i]->vmax;
				if (v1 != v2) 
					Exit("\nError: Huang2021 Fst estimator (-fst_estimator=Huang2021_homo) do not support anisoploids, in individual %s.\n", vind[i]->name);
				pop_nhaplo[p] += v1;
				Nh += v1;
			}
		}

		struct HAP
		{
			double invi;
			double invp;
			int indid;
			ushort popid;
		};

		VLA_NEW(hap, HAP, Nh);
		double invt = 1.0 / Nh;
		ushort *hap_bucket = new ushort[Nh * nloc];
		SetFF(hap_bucket, Nh * nloc);

		//place alleles, xxx
		for (int p = 0, ph = 0; p < n; ++p)
		{
			IND **vind = grps[p]->inds;
			int iend = grps[p]->nind;

			for (int64 l = 0; l < nloc; ++l)
			{
				GENO_ITERATOR rt(grps[p]->ind0id, l, true);//20220316
				GENOTYPE *gtab = GetLoc(l).GetGtab();

				for (int i = 0, ph2 = ph; i < iend; ++i)
				{
					ushort *als = gtab[rt.Read()].GetAlleleArray();
					for (int j = 0, vi = vind[i]->vmin; j < vi; ++j, ++ph2)
						hap_bucket[ph2 * nloc + l] = als[j];
				}
			}


			double invp = 1.0 / pop_nhaplo[p];
			for (int i = 0; i < iend; ++i)
			{
				int vi = vind[i]->vmin;

				hap[ph].indid = i;
				hap[ph].popid = (ushort)p;
				hap[ph].invi = 1.0 / vi;
				hap[ph].invp = invp;

				for (int j = 1; j < vi; ++j)
					hap[ph + j] = hap[ph];

				ph += vi;
			}
		}

		VLA_DELETE(pop_nhaplo);

		//calc ns
		for (int p = 0; p < n; ++p)
		{
			int vp = 0;
			double psvi2 = 0;
			IND **vind2 = grps[p]->inds;
			for (int i = 0, iend = grps[p]->nind; i < iend; ++i)
			{
				int vi = vind2[i]->vmin; 
				psvi2 += vi * vi;
				vp += vi;
			}
			svi2dvt += psvi2;
			svi2dvp += psvi2 / vp;
			svp2dvt += vp * vp;
		}
		svi2dvt /= Nh;
		svp2dvt /= Nh;

		double n11 = 0;
		double n21 = 0, n22 = 0;
		double n31 = 0, n32 = 0, n33 = 0;

		switch (layer)
		{
		case 2:
			n11 = Nh - Np;
			n21 = Nh - svp2dvt;
			n22 = Nh - 1;
			break;
		case 3:
			n11 = Nh - Ni;
			n21 = Nh - svi2dvp; n22 = Nh - Np;
			n31 = Nh - svp2dvt; n32 = Nh - svi2dvt; n33 = Nh - 1;
			break;
		}

		VLA_NEW(lsswi, double, each ? nloc : 0);
		VLA_NEW(lsswp, double, each ? nloc : 0);
		VLA_NEW(lsstot, double, each ? nloc : 0);
		double sswi = 0, sswp = 0, sstot = 0;
		double vwi = 0, vai = 0, vwp = 0, vap = 0, vtot = 0;

		for (int i = 0; i < Nh; ++i)
		{
			HAP& hi = hap[i];
			ushort *ai = hap_bucket + i * nloc;
			for (int j = 0; j < i; ++j)
			{
				HAP& hj = hap[j];
				ushort *aj = hap_bucket + j * nloc;
				double gd = 0;
				for (int64 l = 0; l < nloc; ++l)
				{
					int k2 = GetLoc(l).k;
					double tgd = 0;
					if (isiam)
					{
						if (ai[l] == 0xFFFF && aj[l] == 0xFFFF)
							tgd = 1 - SumProd(grps[hap[i].popid]->GetFreq(l), grps[hap[j].popid]->GetFreq(l), k2);
						else if (ai[l] == 0xFFFF)
							tgd = 1 - grps[hap[i].popid]->GetFreq(l, aj[l]);
						else if (aj[l] == 0xFFFF)
							tgd = 1 - grps[hap[j].popid]->GetFreq(l, ai[l]);
						else if (ai[l] != aj[l])
							tgd = 1;
					}
					else
					{
						if (ai[l] == 0xFFFF && aj[l] == 0xFFFF)
							tgd = SumProdSMM(GetLoc(l).GetAlenArray(), grps[hap[i].popid]->GetFreq(l), grps[hap[j].popid]->GetFreq(l), k2);
						else if (ai[l] == 0xFFFF)
							tgd = SumProdSMM(GetLoc(l).GetAlenArray(), grps[hap[i].popid]->GetFreq(l), aj[l], k2);
						else if (aj[l] == 0xFFFF)
							tgd = SumProdSMM(GetLoc(l).GetAlenArray(), grps[hap[j].popid]->GetFreq(l), ai[l], k2);
						else if (ai[l] != aj[l])
							tgd = GetLoc(l).GetSMMDist(ai[l], aj[l]);
					}

					if (each)
					{
						if (hi.indid == hj.indid) lsswi[l] += tgd * hi.invi;
						if (hi.popid == hj.popid) lsswp[l] += tgd * hi.invp;
						lsstot[l] += tgd * invt;
					}
					gd += tgd;
				}
				if (hi.indid == hj.indid) sswi += gd * hi.invi;
				if (hi.popid == hj.popid) sswp += gd * hi.invp;
				sstot += gd * invt;
			}
		}
		delete[] hap_bucket;
		VLA_DELETE(hap);

		if (each) for (int64 l = 0; l < nloc; ++l)
		{
			switch (layer)
			{
			case 2:
				vwp = lsswp[l] / n11;
				vap = (n11 * lsstot[l] - n22 * lsswp[l]) / (n11 * n21);
				vtot = vwp + vap;
				break;
			case 3:
				vwi = lsswi[l] / n11;
				vai = (n11 * lsswp[l] - n22 * lsswi[l]) / (n11 * n21);
				vap = (n22 * n32 * lsswi[l] - n21 * n33 * lsswi[l] - n11 * n32 * lsswp[l] + n11 * n21 * lsstot[l]) / (n11 * n21 * n31);
				vtot = vwi + vai + vap;
				break;
			}
			each[l] = vap / vtot;
		}
		VLA_DELETE(lsswi);
		VLA_DELETE(lsswp);
		VLA_DELETE(lsstot);

		switch (layer)
		{
		case 2:
			vwp = sswp / n11;
			vap = (n11 * sstot - n22 * sswp) / (n11 * n21);
			vtot = vwp + vap;
			break;
		case 3:
			vwi = sswi / n11;
			vai = (n11 * sswp - n22 * sswi) / (n11 * n21);
			vap = (n22 * n32 * sswi - n21 * n33 * sswi - n11 * n32 * sswp + n11 * n21 * sstot) / (n11 * n21 * n31);
			vtot = vwi + vai + vap;
			break;
		}
		return vap / vtot;
	}

	/* Huang 2021 Fst estimator based on multi-level amova */
	TARGET /*static*/ double FST::Fst_Huang2021_aniso(POP **grps, int n, int layer, bool isiam, bool sumss, double *each, double *buf)
	{
		//2 or 3 Level, ansioploids
		if (ad && layer == 3) Exit("\nError: Huang2021 Fst estimator (-fst_estimator=Huang2021_aniso) is incompatible with allelic depth (-ad) option.\n");
		double Vap = 0, Vtot = 0;
		double *Na = buf;
		double SStot = 0, SSwp = 0, SSwi = 0;
		double N11 = 0, N21 = 0, N22 = 0, N31 = 0, N32 = 0, N33 = 0;
		for (int64 l = 0; l < nloc; ++l)
		{
			if (each) each[l] = NA;
			int k2 = GetLoc(l).k;
			ushort *alen2 = GetLoc(l).GetAlenArray();

			if (k2 < 2) continue;

			//assign allele
			int Nh = 0, Np = 0, Ni = 0;
			double n11 = 0, n21 = 0, n22 = 0, n31 = 0, n32 = 0, n33 = 0;
			double svi2dvp = 0, svi2dvt = 0, svp2dvt = 0;
			double sstot = 0, sswi = 0, sswp = 0;

			SetZero(Na, k2);
			for (int p = 0; p < n; ++p)
			{
				int vp = grps[p]->loc_stat[l].nhaplo;
				if (vp == 0) continue;

				double *fp = grps[p]->GetFreq(l);
				IND **pind = grps[p]->inds;
				double vi2 = 0;

				if (layer == 3)
				{
					ushort *gcount = grps[p]->GetGenoCount(l);//20220316
					GENOTYPE *gtab = GetLoc(l).GetGtab();
					int ngeno = GetLoc(l).ngeno;

					for (int gi = 0; gi < ngeno; ++gi)
					{
						int gc = gcount[gi];
						if (gc == 0) continue;

						GENOTYPE &gt = gtab[gi];
						if (gt.Nalleles() == 0) continue;

						sswi += gc * (isiam ? gt.SS_IAM() : gt.SS_SMM(alen2, k2));
						Ni += gc;
						vi2 += gc * gt.Ploidy() * gt.Ploidy();
					}
				}

				Np++;
				Nh += vp;
				svi2dvp += vi2 / vp;
				svi2dvt += vi2;
				svp2dvt += vp * vp;
				AddProd(Na, fp, vp, k2);
				sswp += SSP(fp, k2, vp, true, NULL);
			}

			if (Np < 2) continue;

			svi2dvt /= Nh;
			svp2dvt /= Nh;
			sstot = SSC(Na, k2, true, NULL);

			//ns
			if (layer == 3)
			{
				n11 = Nh - Ni;
				n21 = Nh - svi2dvp; n22 = Nh - Np;
				n31 = Nh - svp2dvt; n32 = Nh - svi2dvt; n33 = Nh - 1;

				SSwi += sswp; SSwp += sswp; SStot += sstot;

				N11 += n11;
				N21 += n21; N22 += n22;
				N31 += n31; N32 += n32; N33 += n33;

				double vwi = sswi / n11;
				double vai = (n11 * sswp - n22 * sswi) / (n11 * n21);
				double vap = (n22 * n32 * sswi - n21 * n33 * sswi - n11 * n32 * sswp + n11 * n21 * sstot) / (n11 * n21 * n31);

				if (IsError(vwi) || IsError(vai) || IsError(vap)) continue;

				Vap += vap;
				Vtot += vwi + vai + vap;
				if (each) each[l] = IsError(vap / (vwi + vai + vap)) ? NA : vap / (vwi + vai + vap);
			}
			else
			{
				n11 = Nh - Np; n21 = Nh - svp2dvt; n22 = Nh - 1;

				SSwp += sswp; SStot += sstot;

				N11 += n11; N21 += n21; N22 += n22;

				double vwp = sswp / n11;
				double vap = (n11 * sstot - n22 * sswp) / (n11 * n21);

				if (IsError(vwp) || IsError(vap)) continue;

				Vap += vap;
				Vtot += vwp + vap;
				if (each) each[l] = IsError(vap / (vwp + vap)) ? NA : vap / (vwp + vap);
			}

		}

		if (layer == 3 && sumss)
		{
			double Vwi = SSwi / N11;
			double Vai = (N11 * SSwp - N22 * SSwi) / (N11 * N21);
			Vap = (N22 * N32 * SSwi - N21 * N33 * SSwi - N11 * N32 * SSwp + N11 * N21 * SStot) / (N11 * N21 * N31);
			Vtot = Vwi + Vai + Vap;
		}
		else if (sumss)
		{
			double Vwp = SSwp / N11;
			Vap = (N11 * SStot - N22 * SSwp) / (N11 * N21);
			Vtot = Vwp + Vap;
		}
		return Vap / Vtot;
	}

	/* Write results file in column format */
	TARGET /*static*/ void FST::ColumnPrint(FILE *fout)
	{
		byte *estimator = fst_estimator_val;
#define FORBEGIN \
		for (int type = 1; type <= 5; ++type)\
		{\
			if (fst_level_val[type] == 0) continue;\
			if ((type == 1 || type == 4 || type == 3) && lreg == 0) continue;\
			FST *Fst = fst_buf[type];\
			int n0 = 0, n1 = 0, n2 = 0;\
			switch (type)\
			{\
			case 1: n0 = lreg; n1 = 1; n2 = 1; break;\
			case 2: n0 = 1; n1 = 1; n2 = 1; break;\
			case 3: n0 = lreg; n1 = 1; n2 = 1; break;\
			case 4: n0 = lreg; n1 = nregt2; n2 = nregt2; break;\
			case 5: n0 = 1; n1 = npop; n2 = npop; break;\
			}\
			for (int rl = (int)n0 - 1; rl >= 0; --rl)\
				{\
					if (type == 3) n2 = nreg[rl];\
					if (type == 4) n1 = n2 = nreg[rl];\
					for (int i = 0; i < n1; ++i)\
					{\
						for (int j = type <= 3 ? 0 : i + 1; j < n2; ++j)\
						{\
							FST *g = Fst + i * n1 + j;
#define FOREND }}}}

		//type 1 among reg in tot, 2 among pop in tot, 3 among pop/reg in reg, 4 between reg, 5 between pop

		fprintf(fout, "%s%s", g_linebreak_val, g_linebreak_val);

		//Line 1
		fprintf(fout, "Locus%cA", g_delimiter_val);
		FORBEGIN
			switch (type)
			{
			case 1: fprintf(fout, "%cAmong all regL%d", g_delimiter_val, rl + 1); break;
			case 2: fprintf(fout, "%cAmong all pops", g_delimiter_val); break;
			case 3:
				if (rl == 0) fprintf(fout, "%cAmong pops in %s", g_delimiter_val, aregs[rl][j]->name);
				else		 fprintf(fout, "%cAmong regsL%d in %s", g_delimiter_val, rl, aregs[rl][j]->name);
				break;
			case 4:  fprintf(fout, "%c%s", g_delimiter_val, aregs[rl][i]->name); break;
			case 5:  fprintf(fout, "%c%s", g_delimiter_val, apops[i]->name); break;
			}
		FOREND
			fprintf(fout, "%s", g_linebreak_val);

		//Line 2
		fprintf(fout, "%cB", g_delimiter_val);
		FORBEGIN
			switch (type)
			{
			case 1: fprintf(fout, "%c", g_delimiter_val); break;
			case 2: fprintf(fout, "%c", g_delimiter_val); break;
			case 3:
				if (rl == 0) fprintf(fout, "%c", g_delimiter_val);
				else		 fprintf(fout, "%c", g_delimiter_val);
				break;
			case 4:  fprintf(fout, "%c%s", g_delimiter_val, aregs[rl][j]->name); break;
			case 5:  fprintf(fout, "%c%s", g_delimiter_val, apops[j]->name); break;
			}
		FOREND
			fprintf(fout, "%s", g_linebreak_val);

		//Line 3
		if (fst_locus_val[1])
		{
			fprintf(fout, "All loci");

			for (int k = 1; k <= N_FST_ESTIMATOR; ++k)
			{
				if (estimator[k] == 0) continue;
				fprintf(fout, "%c%s", g_delimiter_val, FST_ESTIMATOR[k]);

				FORBEGIN
					fprintf(fout, "%c", g_delimiter_val);
				WriteReal(fout, *((&g->Nei1973T) + k - 1));
				FOREND
					fprintf(fout, "%s", g_linebreak_val);
			}

			if (fst_test_val[1])
			{
				fprintf(fout, "%cGenotype G", g_delimiter_val);
				FORBEGIN
					fprintf(fout, "%c", g_delimiter_val);
				WriteReal(fout, g->Genotype_GT);
				FOREND
					fprintf(fout, "%s", g_linebreak_val);

				fprintf(fout, "%cd.f.", g_delimiter_val);
				FORBEGIN
					fprintf(fout, "%c", g_delimiter_val);
				fprintf(fout, "%d", g->Genotype_DFT);
				FOREND
					fprintf(fout, "%s", g_linebreak_val);

				fprintf(fout, "%cP", g_delimiter_val);
				FORBEGIN
					fprintf(fout, "%c", g_delimiter_val);
				WriteReal(fout, g->Genotype_PT);
				FOREND
					fprintf(fout, "%s", g_linebreak_val);
			}

			if (fst_test_val[2])
			{
				fprintf(fout, "%cAllele G", g_delimiter_val);
				FORBEGIN
					fprintf(fout, "%c", g_delimiter_val);
				WriteReal(fout, g->Allele_GT);
				FOREND
					fprintf(fout, "%s", g_linebreak_val);

				fprintf(fout, "%cd.f.", g_delimiter_val);
				FORBEGIN
					fprintf(fout, "%c", g_delimiter_val);
				fprintf(fout, "%d", g->Allele_DFT);
				FOREND
					fprintf(fout, "%s", g_linebreak_val);

				fprintf(fout, "%cP", g_delimiter_val);
				FORBEGIN
					fprintf(fout, "%c", g_delimiter_val);
				WriteReal(fout, g->Allele_PT);
				FOREND
					fprintf(fout, "%s", g_linebreak_val);
			}
		}


		if (fst_locus_val[2]) for (int64 l = 0; l < nloc; ++l)
		{
			fprintf(fout, "%s", GetLoc(l).GetName());

			for (int k = 1; k <= N_FST_ESTIMATOR; ++k)
			{
				if (estimator[k] == 0) continue;
				fprintf(fout, "%c%s", g_delimiter_val, FST_ESTIMATOR[k]);

				FORBEGIN
					fprintf(fout, "%c", g_delimiter_val);
				WriteReal(fout, *(*((&g->Nei1973) + k - 1) + l));
				FOREND
					fprintf(fout, "%s", g_linebreak_val);
			}


			if (fst_test_val[1])
			{
				fprintf(fout, "%cGenotype G", g_delimiter_val);
				FORBEGIN
					fprintf(fout, "%c", g_delimiter_val);
				WriteReal(fout, g->Genotype_G[l]);
				FOREND
					fprintf(fout, "%s", g_linebreak_val);

				fprintf(fout, "%cd.f.", g_delimiter_val);
				FORBEGIN
					fprintf(fout, "%c", g_delimiter_val);
				fprintf(fout, "%d", g->Genotype_DF[l]);
				FOREND
					fprintf(fout, "%s", g_linebreak_val);

				fprintf(fout, "%cP", g_delimiter_val);
				FORBEGIN
					fprintf(fout, "%c", g_delimiter_val);
				WriteReal(fout, g->Genotype_P[l]);
				FOREND
					fprintf(fout, "%s", g_linebreak_val);
			}

			if (fst_test_val[2])
			{
				fprintf(fout, "%cAllele G", g_delimiter_val);
				FORBEGIN
					fprintf(fout, "%c", g_delimiter_val);
				WriteReal(fout, g->Allele_G[l]);
				FOREND
					fprintf(fout, "%s", g_linebreak_val);

				fprintf(fout, "%cd.f.", g_delimiter_val);
				FORBEGIN
					fprintf(fout, "%c", g_delimiter_val);
				fprintf(fout, "%d", g->Allele_DF[l]);
				FOREND
					fprintf(fout, "%s", g_linebreak_val);

				fprintf(fout, "%cP", g_delimiter_val);
				FORBEGIN
					fprintf(fout, "%c", g_delimiter_val);
				WriteReal(fout, g->Allele_P[l]);
				FOREND
					fprintf(fout, "%s", g_linebreak_val);
			}
		}
	}

	/* Write results file in matrix format */
	TARGET /*static*/ void FST::MatrixPrint(FILE *fout, FST *Fst, int n, int type)
	{
		byte *estimator = fst_estimator_val;
		for (int rl = type == 4 ? lreg - 1 : 0; rl >= 0; --rl)
		{
			int n1 = type == 4 ? nreg[rl] : n;
			for (int k = 1; k <= N_FST_ESTIMATOR; ++k)
			{
				if (k <= N_FST_ESTIMATOR && estimator[k] == 0) continue;

				fprintf(fout, "%s%s%s", g_linebreak_val, g_linebreak_val, FST_ESTIMATOR[k]);

				for (int i = 0; i < n1; ++i)
				{
					switch (type)
					{
					case 4: fprintf(fout, "%c%s", g_delimiter_val, aregs[rl][i]->name);  break;
					case 5: fprintf(fout, "%c%s", g_delimiter_val, apops[i]->name);  break;
					}
				}

				FST *g = Fst;
				for (int i = 0; i < n1; ++i)
				{
					switch (type)
					{
					case 4: fprintf(fout, "%s%s", g_linebreak_val, aregs[rl][i]->name);  break;
					case 5: fprintf(fout, "%s%s", g_linebreak_val, apops[i]->name);  break;
					}

					for (int j = 0; j < n1; ++j)
					{
						fprintf(fout, "%c", g_delimiter_val);
						WriteReal(fout, *((&g->Nei1973T) + k - 1));
						g++;
					}
				}
			}
			Fst += n1 * n1;
		}
	}
#endif

#ifndef _GDDIST
	/* Write column format header row for genetic distance estimation */
	TARGET /*static*/ void GDIST::ColumnPrintHeader()
	{
		fprintf(FRES, "%s%sA", g_linebreak_val, g_linebreak_val);
		for (int rl = gdist_type - 2; rl < lreg; ++rl)
			if (rl >= 0)
				fprintf(FRES, "%cRegL%d", g_delimiter_val, rl + 1);
			else
				fprintf(FRES, "%cPop", g_delimiter_val);

		fprintf(FRES, "%cB", g_delimiter_val);

		for (int rl = gdist_type - 2; rl < lreg; ++rl)
			if (rl >= 0)
				fprintf(FRES, "%cRegL%d", g_delimiter_val, rl + 1);
			else
				fprintf(FRES, "%cPop", g_delimiter_val);

		for (int k = 1; k <= (gdist_type == 1 ? N_GD_ESTIMATOR - 2 * N_FST_ESTIMATOR : N_GD_ESTIMATOR); ++k)
			if (gdist_estimator_val[k])
				fprintf(FRES, "%c%s", g_delimiter_val, GD_ESTIMATOR[k]);
	}

	/* Write column format result row for genetic distance estimation */
	TARGET void GDIST::ColumnPrintLine(int i, int j)
	{
		POP *tp = NULL;
		switch (gdist_type)
		{
		case 1:
			fprintf(FRES, "%s%s", g_linebreak_val, ainds[i]->name);
			tp = apops[ainds[i]->popid];
			break;
		case 2:
			fprintf(FRES, "%s%s", g_linebreak_val, apops[i]->name);
			tp = lreg >= 0 ? aregs[0][apops[i]->rid] : NULL;
			break;
		case 3:
		default:
			fprintf(FRES, "%s%s", g_linebreak_val, aregs[gdist_type - 3][i]->name);
			tp = aregs[gdist_type - 2][aregs[gdist_type - 3][i]->rid];
			break;
		}

		for (int rl = gdist_type - 2; rl < lreg; ++rl)
		{
			fprintf(FRES, "%c%s", g_delimiter_val, tp->name);
			tp = aregs[rl + 1][tp->rid];
		}

		switch (gdist_type)
		{
		case 1:
			fprintf(FRES, "%c%s", g_delimiter_val, ainds[j]->name);
			tp = apops[ainds[j]->popid];
			break;
		case 2:
			fprintf(FRES, "%c%s", g_delimiter_val, apops[j]->name);
			tp = lreg >= 0 ? aregs[0][apops[j]->rid] : NULL;
			break;
		case 3:
		default:
			fprintf(FRES, "%c%s", g_delimiter_val, aregs[gdist_type - 3][j]->name);
			tp = aregs[gdist_type - 2][aregs[gdist_type - 3][j]->rid];
			break;
		}

		for (int rl = gdist_type - 2; rl < lreg; ++rl)
		{
			fprintf(FRES, "%c%s", g_delimiter_val, tp->name);
			tp = aregs[rl + 1][tp->rid];
		}

		for (int k = 1; k <= (gdist_type == 1 ? N_GD_ESTIMATOR - 2 * N_FST_ESTIMATOR : N_GD_ESTIMATOR); ++k)
			if (gdist_estimator_val[k])
			{
				fprintf(FRES, "%c", g_delimiter_val);
				WriteReal(FRES, *((&Nei1972) + k - 1));
			}
	}

	/* Write matrix format header for genetic distance estimation */
	TARGET /*static*/ void GDIST::MatrixPrintMatrixHeader(int k, int n)
	{
		if (gdist_estimator_val[k] == 0) return;
		fprintf(TEMP_FILES[k], "%s%s%s", g_linebreak_val, g_linebreak_val, GD_ESTIMATOR[k]);
		for (int i = 0; i < n; ++i)
		{
			switch (gdist_type)
			{
			case 1: fprintf(TEMP_FILES[k], "%c%s", g_delimiter_val, ainds[i]->name); break;
			case 2: fprintf(TEMP_FILES[k], "%c%s", g_delimiter_val, apops[i]->name);  break;
			case 3:
			default:
				fprintf(TEMP_FILES[k], "%c%s", g_delimiter_val, aregs[gdist_type - 3][i]->name);  break;
			}
		}
	}

	/* Write matrix format row header for genetic distance estimation */
	TARGET /*static*/ void GDIST::MatrixPrintRowHeader(int k, int i)
	{
		if (gdist_estimator_val[k] == 0) return;
		switch (gdist_type)
		{
		case 1: fprintf(TEMP_FILES[k], "%s%s", g_linebreak_val, ainds[i]->name); break;
		case 2: fprintf(TEMP_FILES[k], "%s%s", g_linebreak_val, apops[i]->name);  break;
		case 3:
		default:
			fprintf(TEMP_FILES[k], "%s%s", g_linebreak_val, aregs[gdist_type - 3][i]->name);  break;
		}
	}

	/* Write matrix format grid for genetic distance estimation */
	TARGET void GDIST::MatrixPrintCell(int k)
	{
		if (gdist_estimator_val[k] == 0) return;
		fprintf(TEMP_FILES[k], "%c", g_delimiter_val);
		WriteReal(TEMP_FILES[k], *((&Nei1972) + k - 1));
	}

	/* Use population/region allele frequency as the missing data */
	TARGET void GDIST::GetMissingFreq(IND *a, GENOTYPE &gt, int64 l, double *pbuf)
	{
		int k2 = GetLoc(l).k;
		if (gt.Nalleles())
		{
			gt.GetFreq(pbuf, k2);
			return;
		}

		POP *tp = apops[a->popid];
		if (tp->loc_stat[l].nhaplo)
		{
			SetVal(pbuf, tp->GetFreq(l), k2);
			return;
		}

		POP *tr = lreg >= 0 ? aregs[0][tp->rid] : NULL;
		for (int rl = 0; rl < lreg; ++rl)
		{
			if (tp->loc_stat[l].nhaplo)
			{
				SetVal(pbuf, tr->GetFreq(l), k2);
				return;
			}
			tr = aregs[rl + 1][tr->rid];
		}
	}

	/* Calculate genetic distance between two individuals */
	TARGET void GDIST::CalcGD(IND *a, IND *b, double *p1, double *p2)
	{
		byte *estimator = NULL;
		switch (GDIST_METHOD)
		{
		default: break;
		case 1: estimator = gdist_estimator_val; break;
		case 2: estimator = pcoa_estimator_val; break;
		case 3: estimator = cluster_estimator_val; break;
		}

		if ((estimator[12] || estimator[20]) && abs(g_format_val) <= 2)
			Exit("\nError: Reynolds_Slatkin1995 or Slatkin_Slatkin1995 genetic distance estimator uses Stepwise mutation model (smm), which can only be applied for non-vcf input file, and should use the allele size as the identifier. \n");
		if ((estimator[6]) && abs(g_format_val) <= 2)
			Exit("\nError: Goldstein1995 genetic distance estimator uses Stepwise mutation model (smm), which can only be applied for non-vcf input file, and should use the allele size as the identifier. \n");

		if (ad) Exit("\nError: individual genetic distance (-gdist_level=ind -pcoa_level=ind -cluster_level=ind) is incompatible with allelic depth (-ad) option.\n");

		SetZero(&Nei1972, N_GD_ESTIMATOR);

		if (a == b) return;
		INDGD sumgd, tgd;
		SetZero(&sumgd, 1);
		int aid = a->indid, bid = b->indid;
		int64 eL = 0;

		for (int64 l = 0; l < nloc; ++l)
		{
			if (total_pop->loc_stat[l].nhaplo == 0) continue;

			int gid1 = aid, gid2 = bid;
			IND::GetDyadGenotypeIdx(gid1, gid2, l);
			bool cis = gid1 <= gid2;
			bool usetab = GetLoc(l).ngeno < N_MAX_GDTAB;
			HASH hash = usetab ? HashDyadGenotypeIndex(cis ? (gid1 << 16) | gid2 : (gid2 << 16) | gid1) : 0;

			if (usetab && gdtab[l].ContainsKey(hash))
			{
				gdlock[l].lock_shared();
				INDGD gd = gdtab[l][hash];
				gdlock[l].unlock_shared();

				Add((double*)&sumgd, (double*)&gd, N_INDGD);

				if (!cis)
				{
					sumgd.Jx1 += gd.Jx2 - gd.Jx1;
					sumgd.Jx2 += gd.Jx1 - gd.Jx2;
				}
				continue;
			}
			SetZero(&tgd, 1);


			GENOTYPE *gtab = GetLoc(l).GetGtab();
			GENOTYPE &gt1 = gtab[gid1], &gt2 = gtab[gid2];

			if (gt1.Nalleles() && gt2.Nalleles())
				tgd.ABtype = 1;

			if (tgd.ABtype == 1 || gdist_weightmissing_val == 1)
			{
				eL++;
				GetMissingFreq(a, gt1, l, p1);
				GetMissingFreq(b, gt2, l, p2);

				double Sx1 = 0, Sx2 = 0;
				int k2 = GetLoc(l).k;
				ushort *alen2 = GetLoc(l).GetAlenArray();

				for (int i = 0; i < k2; ++i)
				{
					if (estimator[1] || estimator[7])
					{
						//Nei1972, Nei1973
						tgd.Jx1 += p1[i] * p1[i];
						tgd.Jx2 += p2[i] * p2[i];
						tgd.Jxy += p1[i] * p2[i];
					}

					if (estimator[2])
						tgd.Cavalli1967 += sqrt(Max(DOUBLE_UNDERFLOW, p1[i] * p2[i]));

					if (estimator[3])
					{
						tgd.t1 += (p1[i] - p2[i]) * (p1[i] - p2[i]);
						tgd.t2 += p1[i] * p2[i];
					}

					if (estimator[4])
						tgd.Nei1983 += sqrt(Max(DOUBLE_UNDERFLOW, p1[i] * p2[i]));

					if (estimator[5])
						tgd.Euclidean += (p1[i] - p2[i]) * (p1[i] - p2[i]);

					if (abs(g_format_val) > 2 && estimator[6])
					{
						Sx1 += p1[i] * alen2[i];
						Sx2 += p2[i] * alen2[i];
					}

					if (estimator[8])
						tgd.Roger1972 += (p1[i] - p2[i]) * (p1[i] - p2[i]) / 2;
				}
				if (estimator[5])
					tgd.Euclidean = Max(DOUBLE_UNDERFLOW, tgd.Euclidean);

				if (estimator[6])
					tgd.Goldstein1995 = (Sx1 - Sx2) * (Sx1 - Sx2);

				if (estimator[8])
					tgd.Roger1972 = sqrt(Max(DOUBLE_UNDERFLOW, tgd.Roger1972));
			}

			Add((double*)&sumgd, (double*)&tgd, N_INDGD);

			if (usetab && !gdtab[l].ContainsKey(hash))
			{
				if (!cis) Swap(tgd.Jx1, tgd.Jx2);

				gdlock[l].lock();
				gdtab[l][hash] = tgd;
				gdlock[l].unlock();
			}
		}

		if (eL == 0)
		{
			SetZero(&Nei1972, N_GD_ESTIMATOR);
			return;
		}

		if (estimator[1] || estimator[7])
		{
			sumgd.Jxy /= eL;
			sumgd.Jx1 /= eL;
			sumgd.Jx2 /= eL;

			Nei1972 = -log(sumgd.Jxy / sqrt(sumgd.Jx1 * sumgd.Jx2));
			Nei1974 = (sumgd.Jx1 + sumgd.Jx2) / 2 - sumgd.Jxy;
		}

		if (estimator[2])
			Cavalli1967 = 2 / PI * sqrt(2 * (1 - sumgd.Cavalli1967 / eL));

		if (estimator[3])
		{
			sumgd.t2 = 2 * (eL - sumgd.t2);
			Reynolds1983 = sumgd.t2 > 0 ? sqrt(sumgd.t1 / sumgd.t2) : 0;
		}

		if (estimator[4])
			Nei1983 = 1 - sumgd.Nei1983 / eL;

		if (estimator[5])
			Euclidean = sumgd.Euclidean = sqrt(sumgd.Euclidean);

		if (estimator[6])
			Goldstein1995 = sumgd.Goldstein1995 / eL;

		if (estimator[8])
			Roger1972 = sumgd.Roger1972 / eL;
	}

	/* Calculate genetic distance between two populations/regions */
	TARGET void GDIST::CalcGD(POP *a, POP *b, double *buf)
	{
		byte *estimator = NULL;
		switch (GDIST_METHOD)
		{
		default: break;
		case 1: estimator = gdist_estimator_val; break;
		case 2: estimator = pcoa_estimator_val; break;
		case 3: estimator = cluster_estimator_val; break;
		}
		if ((estimator[12] || estimator[20]) && abs(g_format_val) <= 2)
			Exit("\nError: Reynolds_Slatkin1995 or Slatkin_Slatkin1995 genetic distance estimator uses Stepwise mutation model (smm), which can only be applied for non-vcf input file, and should use the allele size as the identifier. \n");
		if ((estimator[6]) && abs(g_format_val) <= 2)
			Exit("\nError: Goldstein1995 genetic distance estimator uses Stepwise mutation model (smm), which can only be applied for non-vcf input file, and should use the allele size as the identifier. \n");

		double Jx1 = 0, Jx2 = 0, Jxy = 0, t1 = 0, t2 = 0;
		SetZero(&Nei1972, N_GD_ESTIMATOR);

		if (a == b) return;
		double ABtype = 0;

		for (int64 l = 0; l < nloc; ++l)
		{
			LOCSTAT *f1 = a->loc_stat + l, *f2 = b->loc_stat + l;
			if (f1->nhaplo && f2->nhaplo) ABtype++;
			else continue;

			double *p1 = a->GetFreq(l), *p2 = b->GetFreq(l);
			int k2 = GetLoc(l).k;
			ushort *alen2 = GetLoc(l).GetAlenArray();

			double Roger1972t = 0, Euclideant = 0, Sx1 = 0, Sx2 = 0;
			for (int i = 0; i < k2; ++i)
			{
				if (estimator[1] || estimator[7])
				{
					//Nei1972, Nei1973
					Jx1 += p1[i] * p1[i];
					Jx2 += p2[i] * p2[i];
					Jxy += p1[i] * p2[i];
				}

				if (estimator[2])
					Cavalli1967 += sqrt(Max(DOUBLE_UNDERFLOW, p1[i] * p2[i]));

				if (estimator[3])
				{
					t1 += (p1[i] - p2[i]) * (p1[i] - p2[i]);
					t2 += p1[i] * p2[i];
				}

				if (estimator[4])
					Nei1983 += sqrt(Max(DOUBLE_UNDERFLOW, p1[i] * p2[i]));

				if (estimator[5])
					Euclideant += (p1[i] - p2[i]) * (p1[i] - p2[i]);

				if (abs(g_format_val) > 2 && estimator[6])
				{
					Sx1 += p1[i] * alen2[i];
					Sx2 += p2[i] * alen2[i];
				}

				if (estimator[8])
					Roger1972t += (p1[i] - p2[i]) * (p1[i] - p2[i]) / 2;
			}
			Roger1972 += sqrt(Max(DOUBLE_UNDERFLOW, Roger1972t));
			Euclidean += Max(DOUBLE_UNDERFLOW, Euclideant);

			if (estimator[6])
				Goldstein1995 += (Sx1 - Sx2) * (Sx1 - Sx2);
		}

		if (ABtype == 0)
		{
			SetZero(&Nei1972, N_GD_ESTIMATOR);
			return;
		}

		if (estimator[1] || estimator[7])
		{
			Jxy /= ABtype;
			Jx1 /= ABtype;
			Jx2 /= ABtype;

			Nei1972 = -log(Jxy / sqrt(Jx1 * Jx2));
			Nei1974 = (Jx1 + Jx2) / 2 - Jxy;
		}

		if (estimator[2])
			Cavalli1967 = 2 / PI * sqrt(2 * (1 - Cavalli1967 / ABtype));

		if (estimator[3])
		{
			t2 = 2 * (ABtype - t2);
			Reynolds1983 = t2 > 0 ? sqrt(t1 / t2) : 0;
		}

		if (estimator[4])
			Nei1983 = 1 - Nei1983 / ABtype;

		if (estimator[5])
			Euclidean = sqrt(Euclidean);

		if (estimator[6])
			Goldstein1995 /= ABtype;

		if (estimator[8])
			Roger1972 /= ABtype;

		POP *gs[] = { a, b };
		for (int i = 1; i <= N_FST_ESTIMATOR; ++i)
		{
			if (estimator[8 + i] || estimator[8 + N_FST_ESTIMATOR + i])
			{
				double Fst = FST::FstEstimator(gs, 2, i, NULL, buf);
				if (estimator[8 + i])				    *(double*)(&Slatkin_Nei1973 + i - 1) = Fst / (1 - Fst);
				if (estimator[8 + N_FST_ESTIMATOR + i]) *(double*)(&Reynolds_Nei1973 + i - 1) = -log(1 - Fst);
			}
		}
	}
#endif

#ifndef _VESSEL_ITERATOR
	/* Go to start */
	TARGET void VESSEL_ITERATOR::Rewind(int nlay)
	{
		SetZero(relative_id, nlay + 1);
		SetZero(universal_id, nlay + 1);
		for (int clay = nlay - 1; clay >= lay; --clay)
			trace[clay] = trace[clay + 1]->subunits[0];
	}

	/* Copy from a reference*/
	TARGET void VESSEL_ITERATOR::Copy(VESSEL_ITERATOR &ref, int nlay)
	{
		SetVal(relative_id, ref.relative_id, nlay + 1);
		SetVal(universal_id, ref.universal_id, nlay + 1);
		SetVal(trace, ref.trace, nlay + 1);
		lay = ref.lay;
	}

	/* Initialize */
	TARGET VESSEL_ITERATOR::VESSEL_ITERATOR()
	{

	}

	/* Initialize */
	TARGET VESSEL_ITERATOR::VESSEL_ITERATOR(int _lay, VESSEL &root, int nlay)
	{
		lay = _lay;
		SetZero(relative_id, nlay + 1);
		SetZero(universal_id, nlay + 1);
		trace[nlay] = &root;
		for (int clay = nlay - 1; clay >= lay; --clay)
			trace[clay] = trace[clay + 1]->subunits[0];
	}

	/* Uninitialize */
	TARGET VESSEL_ITERATOR::~VESSEL_ITERATOR()
	{

	}

	/* Go to next vessel */
	TARGET void VESSEL_ITERATOR::Next(int nlay)
	{
		//Higher level vessel have remaining vessels, use rest vessels
		if (trace[lay + 1]->nsubunits > relative_id[lay] + 1)
		{
			relative_id[lay]++;
			universal_id[lay]++;
			trace[lay] = trace[lay + 1]->subunits[relative_id[lay]];
		}
		else
		{
			//Higher level vessel do not have remaining vessels
			for (int clay = lay; clay < nlay; ++clay)
			{
				if (trace[clay + 1]->nsubunits > relative_id[clay] + 1)
				{
					universal_id[clay]++;
					relative_id[clay]++;
					trace[clay] = trace[clay + 1]->subunits[relative_id[clay]];
					for (int clay2 = clay - 1; clay2 >= lay; --clay2)
						trace[clay2] = trace[clay2 + 1]->subunits[relative_id[clay2]];
					break;
				}
				else
				{
					universal_id[clay]++;
					relative_id[clay] = 0;
					trace[lay] = NULL;
				}
			}
		}
		//return trace[lay];
	}

	/* Get haplotype index to calculate genetic distance */
	TARGET int VESSEL_ITERATOR::GetHapId()
	{
		return trace[lay]->hid;
	}

	/* Get allele to fetch calculate distance */
	TARGET ushort VESSEL_ITERATOR::GetAllele()
	{
		return trace[lay]->allele;
	}

	/* Get subpopulation in print SS */
	TARGET POP *VESSEL_ITERATOR::GetSubpop(int nlay, int tlay)
	{
		POP *tpop = total_pop;
		for (int clay = nlay - 1; clay >= tlay; --clay)
			tpop = tpop->vpop[relative_id[clay]];
		return tpop;
	}

	/* Get individual in print SS */
	TARGET IND *VESSEL_ITERATOR::GetInd(int nlay, int tlay)
	{
		POP *tpop = total_pop;
		for (int clay = nlay - 1; clay >= tlay; --clay)
			if (clay > tlay)
				tpop = tpop->vpop[relative_id[clay]];
			else
				return tpop->inds[relative_id[clay]];
		return NULL;
	}
#endif

#ifndef _VESSEL
	/* Uninitialize */
	TARGET VESSEL::~VESSEL()
	{

	}

	/* Initialize */
	TARGET VESSEL::VESSEL()
	{
		subunits = NULL;
		nhaplos = NULL;
		allelecount = NULL;
		nsubunits = nhaplo = 0;
		lay = -1;
		hid = -1;
		allele = 0xFFFF;
	}

	/* Deep copy a vessel */
	TARGET VESSEL::VESSEL(VESSEL &r)
	{
		subunits = NULL;
		nhaplos = NULL;
		allelecount = NULL;
		nhaplo = r.nhaplo;
		lay = r.lay;
		hid = r.hid;
		allele = r.allele;
		nsubunits = r.nsubunits;

		if (r.subunits != NULL)
		{
			amova_memory->Alloc(subunits, nsubunits);
			for (int i = 0; i < nsubunits; ++i)
				subunits[i] = new(amova_memory->Alloc(sizeof(VESSEL))) VESSEL(*r.subunits[i]);
		}

		if (r.nhaplos != NULL)
		{
			amova_memory->Alloc(nhaplos, nloc);
			SetVal(nhaplos, r.nhaplos, nloc);
		}

		if (r.allelecount != NULL)
		{
			amova_memory->Alloc(allelecount, KT);
			SetVal(allelecount, r.allelecount, KT);
		}
	}

	/* Create vessel from population */
	TARGET VESSEL::VESSEL(POP *s, int _lay, int &_hid, int64 loc, int method)
	{
		subunits = NULL;
		nhaplos = NULL;
		allelecount = NULL;
		lay = _lay;
		nsubunits = 0;
		nhaplo = 0;
		hid = -1;
		allele = 0xFFFF;

		if (!s->ispop)
		{
			if (loc == -1)
			{
				//region homo, ml
				nsubunits = s->npop;
				amova_memory->Alloc(subunits, nsubunits);
				nhaplo = s->nhaplotypes;

				if (method == 4)
				{
					amova_memory->Alloc(nhaplos, nloc);
					SetZero(nhaplos, nloc);
					amova_memory->Alloc(allelecount, KT);
					SetZero(allelecount, KT);
				}

				for (int i = 0; i < nsubunits; ++i)
				{
					subunits[i] = new(amova_memory->Alloc(sizeof(VESSEL))) VESSEL(s->vpop[i], lay - 1, _hid, loc, method);
					if (method == 4)
					{
						Add(nhaplos, subunits[i]->nhaplos, nloc);
						Add(allelecount, subunits[i]->allelecount, KT);
					}
				}
			}
			else
			{
				//region aniso
				for (int i = 0; i < s->npop; ++i)
					if (s->vpop[i]->loc_stat[loc].nhaplo > 0)
						nsubunits++;

				amova_memory->Alloc(subunits, nsubunits);
				nhaplo = s->loc_stat[loc].nhaplo;

				for (int i = 0, ic = 0; i < s->npop; ++i)
				{
					if (s->vpop[i]->loc_stat[loc].nhaplo > 0)
						subunits[ic++] = new(amova_memory->Alloc(sizeof(VESSEL))) VESSEL(s->vpop[i], lay - 1, _hid, loc, method);
					else
						_hid += s->vpop[i]->nhaplotypes;
				}
			}
		}
		else if (amova_cind_val == 1)
		{
			//pop, ind
			if (loc == -1)
			{
				//pop, ind, homo/ml
				nsubunits = s->nind;
				amova_memory->Alloc(subunits, nsubunits);
				nhaplo = s->nhaplotypes;
				if (method < 3) //homo, aniso
					for (int i = 0; i < nsubunits; ++i)
						subunits[i] = new(amova_memory->Alloc(sizeof(VESSEL))) VESSEL(s->inds[i], lay - 1, _hid, loc, method);
				else
				{
					//ml
					amova_memory->Alloc(nhaplos, nloc);
					SetZero(nhaplos, nloc);
					amova_memory->Alloc(allelecount, KT);
					SetZero(allelecount, KT);

					for (int i = 0; i < nsubunits; ++i)
					{
						subunits[i] = new(amova_memory->Alloc(sizeof(VESSEL))) VESSEL(s->inds[i], lay - 1, _hid, loc, method);
						Add(nhaplos, subunits[i]->nhaplos, nloc);
					}

					for (int64 l = 0; l < nloc; ++l)
					{
						int *acount = GetAlleleCount(l);
						ushort *gcount = s->GetGenoCount(l);//20220316
						GENOTYPE *gtab = GetLoc(l).GetGtab();
						int ngeno = GetLoc(l).ngeno;

						for (int gi = 0; gi < ngeno; ++gi)
						{
							int gc = gcount[gi];
							if (gc == 0) continue;

							GENOTYPE &gt = gtab[gi];
							if (gt.Nalleles() == 0) continue;

							ushort *als = gt.GetAlleleArray();
							for (int j = 0, vi = gt.Ploidy(); j < vi; ++j)
								acount[als[j]] += gc;
						}
					}
				}
			}
			else
			{
				//pop, ind, aniso
				ushort *gcount = s->GetGenoCount(loc);//20220316
				GENOTYPE *gtab = GetLoc(loc).GetGtab();
				int ngeno = GetLoc(loc).ngeno;

				for (int gi = 0; gi < ngeno; ++gi)
				{
					int gc = gcount[gi];
					if (gc == 0) continue;

					GENOTYPE &gt = gtab[gi];
					if (gt.Nalleles() == 0) continue;
					nsubunits += gc;
				}

				amova_memory->Alloc(subunits, nsubunits);
				nhaplo = s->loc_stat[loc].nhaplo;

				GENO_ITERATOR rt(s->ind0id, loc, true);
				for (int i = 0, ic = 0, iend = s->nind; i < iend; ++i)
				{
					GENOTYPE &gt = gtab[rt.Read()];
					if (gt.Nalleles())
						subunits[ic++] = new(amova_memory->Alloc(sizeof(VESSEL))) VESSEL(s->inds[i], lay - 1, _hid, loc, method);
					else
						_hid += gt.Ploidy();
				}
			}
		}
		else
		{
			//pop, no ind
			if (loc == -1)
			{
				//pop, no ind, homo/ml
				if (method < 3)
				{
					//homo
					nsubunits = nhaplo = s->nhaplotypes;
					amova_memory->Alloc(subunits, nsubunits);
					for (int i = 0, ic = 0, iend = s->nind; i < iend; ++i)
						for (int j = 0, vi = s->inds[i]->vmin; j < vi; ++j)
							subunits[ic++] = new(amova_memory->Alloc(sizeof(VESSEL))) VESSEL(lay - 1, _hid, 0xFFFF);
				}
				else
				{
					//ml
					nhaplo = s->nhaplotypes;
					nsubunits = s->nind;
					amova_memory->Alloc(subunits, nsubunits);

					amova_memory->Alloc(nhaplos, nloc);
					SetZero(nhaplos, nloc);
					amova_memory->Alloc(allelecount, KT);
					SetZero(allelecount, KT);

					for (int i = 0; i < nsubunits; ++i)
					{
						subunits[i] = new(amova_memory->Alloc(sizeof(VESSEL))) VESSEL(s->inds[i], lay - 1, _hid, loc, method);
						Add(nhaplos, subunits[i]->nhaplos, nloc);
					}

					for (int64 l = 0; l < nloc; ++l)
					{
						int* acount = GetAlleleCount(l);
						ushort* gcount = s->GetGenoCount(l);//20220316
						GENOTYPE* gtab = GetLoc(l).GetGtab();
						int ngeno = GetLoc(l).ngeno;

						for (int gi = 0; gi < ngeno; ++gi)
						{
							int gc = gcount[gi];
							if (gc == 0) continue;

							GENOTYPE& gt = gtab[gi];
							if (gt.Nalleles() == 0) continue;

							ushort* als = gt.GetAlleleArray();
							for (int j = 0, vi = gt.Ploidy(); j < vi; ++j)
								acount[als[j]] += gc;
						}
					}
				}
			}
			else
			{
				//pop, no ind, aniso
				ushort *gcount = s->GetGenoCount(loc);//20220316
				GENOTYPE *gtab = GetLoc(loc).GetGtab();
				int ngeno = GetLoc(loc).ngeno;

				for (int gi = 0; gi < ngeno; ++gi)
				{
					int gc = gcount[gi];
					if (gc == 0) continue;

					GENOTYPE &gt = gtab[gi];
					if (gt.Nalleles() == 0) continue;
					nsubunits += gc * gt.Ploidy();
				}

				amova_memory->Alloc(subunits, nsubunits);
				nhaplo = s->loc_stat[loc].nhaplo;

				GENO_ITERATOR rt(s->ind0id, loc, true);//20220316
				for (int i = 0, ic = 0, iend = s->nind; i < iend; ++i)
				{
					GENOTYPE &gt = gtab[rt.Read()];
					ushort *alleles = gt.GetAlleleArray();
					if (gt.Nalleles()) for (int j = 0, vi = gt.Ploidy(); j < vi; ++j)
						subunits[ic++] = new(amova_memory->Alloc(sizeof(VESSEL))) VESSEL(lay - 1, _hid, alleles[j]);
					else
						_hid += gt.Ploidy();
				}
			}
		}
	}

	/* Create vessel from individual */
	TARGET VESSEL::VESSEL(IND *s, int _lay, int &_hid, int64 loc, int method)
	{
		subunits = NULL;
		nhaplos = NULL;
		allelecount = NULL;
		lay = _lay;
		nhaplo = nsubunits = 0;
		hid = -1;
		allele = 0xFFFF;

		if (method < 3)
		{
			//homom aniso
			GENOTYPE &gt = s->GetGenotype(loc == -1 ? 0 : loc);
			ushort *alleles = gt.GetAlleleArray();
			nsubunits = nhaplo = gt.Ploidy();
			amova_memory->Alloc(subunits, nhaplo);
			for (int i = 0; i < nhaplo; ++i)
				subunits[i] = new(amova_memory->Alloc(sizeof(VESSEL))) VESSEL(lay - 1, _hid, loc == -1 ? 0xFFFF : alleles[i]);
		}
		else
		{
			//ml
			amova_memory->Alloc(nhaplos, nloc);
			amova_memory->Alloc(allelecount, KT);
			SetZero(allelecount, KT);

			hid = _hid++;
			nhaplo = s->vmin;
			for (int64 l = 0; l < nloc; ++l)
			{
				int *acount = GetAlleleCount(l);
				GENOTYPE &gt = s->GetGenotype(l);//fine
				nhaplos[l] = gt.Nalleles() ? gt.Ploidy() : 0;
				if (gt.Nalleles())
				{
					ushort *als = gt.GetAlleleArray();
					for (int j = 0, vi = gt.Ploidy(); j < vi; ++j)
						acount[als[j]]++;
				}
			}
		}
	}

	/* Create vessel from haplotype */
	TARGET VESSEL::VESSEL(int _lay, int &_hid, ushort _allele)
	{
		subunits = NULL;
		nhaplos = NULL;
		allelecount = NULL;
		nsubunits = 0;
		nhaplo = 1;
		lay = _lay;
		hid = _hid++;
		allele = _allele;
	}

	/* Get the allele count array at locus l*/
	TARGET int *VESSEL::GetAlleleCount(int l)
	{
		return allelecount + allele_freq_offset[l];
	}

	/* Save all vellels in level fa into an array */
	TARGET void VESSEL::GetVessels(VESSEL **vs, int &nvessels, int fa)
	{
		if (fa == lay) for (int i = 0; i < nsubunits; ++i)
			vs[nvessels++] = subunits[i];
		else for (int i = 0; i < nsubunits; ++i)
			subunits[i]->GetVessels(vs, nvessels, fa);
	}

	/* Replace with shuffled vessels */
	TARGET int VESSEL::Replace(VESSEL **vs, int &nvessels, int fa, int method)
	{
		if (fa == lay)
		{
			nhaplo = 0;
			if (method == 4)
			{
				SetZero(nhaplos, nloc);
				SetZero(allelecount, KT);
			}

			for (int i = 0; i < nsubunits; ++i)
			{
				subunits[i] = vs[nvessels++];
				nhaplo += subunits[i]->nhaplo;

				if (method == 4)
				{
					Add(nhaplos, subunits[i]->nhaplos, nloc);
					Add(allelecount, subunits[i]->allelecount, KT);
				}
			}
		}
		else
		{
			nhaplo = 0;
			if (method == 4)
			{
				SetZero(nhaplos, nloc);
				SetZero(allelecount, KT);
			}

			for (int i = 0; i < nsubunits; ++i)
			{
				nhaplo += subunits[i]->Replace(vs, nvessels, fa, method);

				if (method == 4)
				{
					Add(nhaplos, subunits[i]->nhaplos, nloc);
					Add(allelecount, subunits[i]->allelecount, KT);
				}
			}
		}
		return nhaplo;
	}

	/* Shuffle fa level vessels among fb level vessels */
	TARGET void VESSEL::Shuffle(RNG &rng, int fa, int fb, int method, VESSEL **buf)
	{
		if (lay > fb) for (int i = 0; i < nsubunits; ++i)
			subunits[i]->Shuffle(rng, fa, fb, method, buf);
		else if (fb == lay)
		{
			int nvessels = 0;
			//VLA_NEW(vs, VESSEL*, nhaplo);
			GetVessels(buf, nvessels, fa);
			rng.Permute(buf, nvessels);
			nvessels = 0;
			Replace(buf, nvessels, fa, method);
			//VLA_DELETE(vs);
		}
	}

	/* Calculate matrix C for maximum-likelihood method */
	TARGET void VESSEL::GetCML(double *C, int64 l, int *tid, double *tw, int Nh, int nlay, double **W)
	{
		InitC(C, tid, Nh, nlay);
		GetC(tw, tid, C, nlay, W, l);
	}

	/* Initialize matrix C, the coefficient matrix of S = CV */
	TARGET void VESSEL::InitC(double *C, int *tid, int Nh, int nlay)
	{
		SetZero(C, nlay * nlay);
		for (int i = 0; i < nlay; ++i)
			for (int j = 0; j <= i; ++j)
				C[i * nlay + j] = Nh;
		SetZero(tid, nlay + 1);
		return;
	}

	/* Calculate matrix C */
	TARGET void VESSEL::GetC(double *tw, int *tid, double *C, int nlay, double **W, int64 l)
	{
		if (l == -1)
		{
			if (nhaplo == 0) return;
			tw[lay] = 1.0 / nhaplo;
			for (int i = lay; i < nlay; ++i)
				C[i * nlay + lay] -= nhaplo * nhaplo * tw[i + 1];
			for (int i = 0; i < nsubunits; ++i)
				subunits[i]->GetC(tw, tid, C, nlay, W, -1);
			W[lay][tid[lay]++] = tw[lay];
		}
		else
		{
			if (nhaplos[l] == 0) { tid[lay]++; return; }
			tw[lay] = 1.0 / nhaplos[l];
			for (int i = lay; i < nlay; ++i)
				C[i * nlay + lay] -= nhaplos[l] * nhaplos[l] * tw[i + 1];

			if (amova_cind_val == 2)
			{
				if (lay > 1) for (int i = 0; i < nsubunits; ++i)
					subunits[i]->GetC(tw, tid, C, nlay, W, l);
				else for (int i = lay - 1; i < nlay; ++i)
					C[i * nlay + lay - 1] -= nhaplos[l] * tw[i + 1];
			}
			else
			{
				if (subunits != NULL) for (int i = 0; i < nsubunits; ++i)
					subunits[i]->GetC(tw, tid, C, nlay, W, l);
				else for (int i = lay - 1; i < nlay; ++i)
					C[i * nlay + lay - 1] -= nhaplos[l] * tw[i + 1];
			}
			W[lay][tid[lay]++] = tw[lay];
		}
	}

	/* Count number of vessels in each level */
	TARGET void VESSEL::CountVessels(int *count)
	{
		count[lay]++;
		for (int i = 0; i < nsubunits; ++i)
			subunits[i]->CountVessels(count);
	}

	/* Initialize W, W[lay][tid[lay]] = 1 / nhaplo */
	TARGET void VESSEL::InitW(MEMORY &mem, double **&W, int nlay)
	{
		VLA_NEW(count, int, nlay + 1);
		SetZero(count, nlay + 1);
		CountVessels(count);
		mem.Alloc(W, nlay + 1);
		for (int clay = 0; clay <= nlay; ++clay)
			mem.Alloc(W[clay], count[clay]);
		VLA_DELETE(count);
	}

	/* Calculate SS for homoploid method */
	TARGET void VESSEL::GetSSHomo(double *SS, double *gd, int Nh, double **W, int nlay, VESSEL_ITERATOR &ve1, VESSEL_ITERATOR &ve2)
	{
		SetZero(SS, nlay);
		ve1.Rewind(nlay);
		for (int i = 0; i < Nh; ++i)
		{
			double *gd1 = gd + ve1.GetHapId() * Nh;
			ve2.Copy(ve1, nlay); ve2.Next(nlay);
			for (int j = i + 1; j < Nh; ++j)
			{
				double GDt = gd1[ve2.GetHapId()];
				for (int clay = 1; clay <= nlay; ++clay)
					if (ve1.universal_id[clay] == ve2.universal_id[clay])
						SS[clay - 1] += GDt * W[clay][ve1.universal_id[clay]];
				ve2.Next(nlay);
			}
			ve1.Next(nlay);
		}
	}

	/* Calculate SS for anisoploid method */
	TARGET void VESSEL::GetSSAniso(ushort *hap_bucket, double *SS, bool isiam, int nh, int64 l, double **W, int nlay, VESSEL_ITERATOR &ve1, VESSEL_ITERATOR &ve2)
	{
		SetZero(SS, nlay);
		ve1.Rewind(nlay);
		for (int i = 0; i < nh; ++i)
		{
			ve2.Copy(ve1, nlay); ve2.Next(nlay);
			ushort *hi = hap_bucket + ve1.GetHapId() * nloc;
			for (int j = i + 1; j < nh; ++j)
			{
				ushort *hj = hap_bucket + ve2.GetHapId() * nloc;
				double tdist = isiam ? hi[l] != hj[l] : GetLoc(l).GetSMMDist(hi[l], hj[l]);
				for (int clay = 1; clay <= nlay; ++clay)
					if (ve1.universal_id[clay] == ve2.universal_id[clay])
						SS[clay - 1] += tdist * W[clay][ve1.universal_id[clay]];
				ve2.Next(nlay);
			}
			ve1.Next(nlay);
		}
	}

	/* Calculate SS for anisoploid method */
	TARGET void VESSEL::GetSSAniso(double *SS, bool isiam, int nh, int64 l, double **W, int nlay, VESSEL_ITERATOR &ve1, VESSEL_ITERATOR &ve2)
	{
		SetZero(SS, nlay);
		ve1.Rewind(nlay);
		for (int i = 0; i < nh; ++i)
		{
			ve2.Copy(ve1, nlay); ve2.Next(nlay);
			ushort a1 = ve1.GetAllele();
			for (int j = i + 1; j < nh; ++j)
			{
				ushort a2 = ve2.GetAllele();
				double tdist = isiam ? a1 != a2 : GetLoc(l).GetSMMDist(a1, a2);
				for (int clay = 1; clay <= nlay; ++clay)
					if (ve1.universal_id[clay] == ve2.universal_id[clay])
						SS[clay - 1] += tdist * W[clay][ve1.universal_id[clay]];
				ve2.Next(nlay);
			}
			ve1.Next(nlay);
		}
	}


	/* Calculate variance component matrix V */
	TARGET /*static*/ void VESSEL::GetV(double *C, double *SS, double *&V, int nlay)
	{
		Map<MatrixXd> MC(C, nlay, nlay);//symmetric
		MC = MC.inverse();
		//MatrixInv(C, nlay);
		MatrixMul(C, nlay, nlay, SS, nlay, 1, V);
	}

	/* Calculate F-statistics */
	TARGET /*static*/ void VESSEL::GetF(double *V, double *F, double *vs, int nlay)
	{
		vs[0] = V[0];
		for (int i = 1; i < nlay; ++i)
			vs[i] = vs[i - 1] + V[i];
		for (int i = 0; i < nlay; ++i)
			for (int j = i + 1; j < nlay; ++j)
				F[i * nlay + j] = 1.0 - vs[i] / vs[j];
	}

	/* Calculate F-statistics */
	TARGET /*static*/ void VESSEL::GetF(double *Fi, double *F, int nlay)
	{
		if (amova_cind_val == 2)
		{
			for (int i = 0; i < nlay; ++i)
				for (int j = i + 1; j < nlay; ++j)
					F[i * nlay + j] = 1 - (1 - Fi[j]) / (1 - Fi[i]);
		}
		else
		{
			for (int i = 1; i < nlay; ++i)
				F[0 * nlay + i] = Fi[i];
			for (int i = 1; i < nlay; ++i)
				for (int j = i + 1; j < nlay; ++j)
					F[i * nlay + j] = 1 - (1 - Fi[j]) / (1 - Fi[i]);
		}

	}
#endif

#ifndef _AMOVA

	/* Initialize */
	TARGET AMOVA::AMOVA()
	{
		Lind = amova_cind_val == 1;
		nlay = Lind + 2 + lreg;
		SSW = new double *[nlay];									nSSW = new int[nlay];
		if (Lind)
		{
			SSW[0] = new double[nind];								SetZero(SSW[0], nSSW[0] = nind);
			SSW[1] = new double[npop];								SetZero(SSW[1], nSSW[1] = npop);
			for (int clay = 2; clay < nlay; ++clay)
			{ SSW[clay] = new double[nreg[clay - 2]];				SetZero(SSW[clay], nSSW[clay] = nreg[clay - 2]); }
		}
		else
		{
			SSW[0] = new double[npop];								SetZero(SSW[0], nSSW[0] = npop);
			for (int clay = 1; clay < nlay; ++clay)
			{ SSW[clay] = new double[nreg[clay - 1]];				SetZero(SSW[clay], nSSW[clay] = nreg[clay - 1]); }
		}

		DF = new double[nlay + 1];									SetZero(DF, nlay + 1);
		V = new double[nlay];										SetZero(V, nlay);
		SS = new double[nlay];										SetZero(SS, nlay);
		F = new double[nlay * nlay]; 								SetZero(F, nlay * nlay);
		EF = new double[nlay * nlay];								SetZero(EF, nlay * nlay);
		EF2 = new double[nlay * nlay];								SetZero(EF2, nlay * nlay);
		G = new int[nlay * nlay];									SetZero(G, nlay * nlay);
		E = new int[nlay * nlay];									SetZero(E, nlay * nlay);
	}

	/* Extract dummy haplotype for homoploid method */
	TARGET void AMOVA::GetHaplotype(POP *tpop, ushort *bucket)
	{

		for (int i = 0; i < nind; ++i)
		{
			IND *tind = ainds[i];
			int vi = tind->vmin, vi2 = tind->vmax;
			if (vi != vi2)
				Exit("\nError: Homoploid amova estimator can only be used for homoploids. Error in individual %s.\n", tind->name);
		}

		for (int64 l = 0; l < nloc; ++l)
		{
			GENOTYPE *gtab = GetLoc(l).GetGtab();//20220316
			GENO_ITERATOR rt(0u, l, true);

			for (int i = 0, ph = 0; i < nind; ++i)
			{
				GENOTYPE &gt = gtab[rt.Read()]; 
				ushort *als = gt.GetAlleleArray();
				for (int j = 0, vi = gt.Ploidy(); j < vi; ++j, ++ph)
					bucket[ph * nloc + l] = als[j];
			}
		}
	}

	/* Perform AMOVA using homoploid method */
	TARGET void AMOVA::CalcAMOVA_homo()
	{
		//single thread
		MEMORY mem1;		VLA_NEW(Gmem2, MEMORY, g_nthread_val);
		int64 nperm = amova_nperm_val * (amova_test_val == 1);

		method = amova_cmethod_val;
		Lind = amova_cind_val == 1;
		nlay = Lind + 2 + lreg;

		//Dummy haplotype, homoploids
		bool isiam = amova_cmutation_val == 1;
		if (!isiam && abs(g_format_val) <= 2)
			Exit("\nError: Stepwise mutation model (smm) in AMOVA can only be applied for non-vcf input file, and should use size as allele identifier. \n");

		for (int i = 0; i < nind; ++i)
			if (ainds[i]->vmin != ainds[i]->vmax)
				Exit("\nError: Homoploid AMOVA method do not support anisoploids, in individual %s.\n", ainds[i]->name);

		int Nh = total_pop->nhaplotypes;
		VESSEL **Gpermbuf = new VESSEL*[Nh * g_nthread_val];		SetZero(Gpermbuf, Nh * g_nthread_val);
		double *dist = new double[Nh * Nh];							SetZero(dist, Nh * Nh);
#define distW(x,y) dist[(x)*Nh+(y)]

		VLA_NEW(Gss,  double, nlay * g_nthread_val);				SetZero(Gss,  nlay * g_nthread_val);
		VLA_NEW(Gv,   double, nlay * g_nthread_val);				SetZero(Gv,   nlay * g_nthread_val);
		VLA_NEW(Gvs,  double, nlay * g_nthread_val);				SetZero(Gvs,  nlay * g_nthread_val);
		VLA_NEW(Gtid, int,    (nlay + 1) * g_nthread_val);			SetZero(Gtid, (nlay + 1) * g_nthread_val);
		VLA_NEW(Gtw,  double, (nlay + 1) * g_nthread_val);			SetZero(Gtw,  (nlay + 1) * g_nthread_val);
		VLA_NEW(Gf,   double, nlay * nlay * g_nthread_val);			SetZero(Gf,   nlay * nlay * g_nthread_val);
		VLA_NEW(GC,   double, nlay * nlay * g_nthread_val);			SetZero(GC,   nlay * nlay * g_nthread_val);
		VLA_NEW(Gg,   int,    nlay * nlay * g_nthread_val);			SetZero(Gg,   nlay * nlay * g_nthread_val);
		VLA_NEW(Ge,   int,    nlay * nlay * g_nthread_val);			SetZero(Ge,   nlay * nlay * g_nthread_val);
		VLA_NEW(Gef,  double, nlay * nlay * g_nthread_val);			SetZero(Gef,  nlay * nlay * g_nthread_val);
		VLA_NEW(Gef2, double, nlay * nlay * g_nthread_val);			SetZero(Gef2, nlay * nlay * g_nthread_val);

		int tref = 0;
		amova_memory = &mem1; amova_memory->ClearMemory();
		VESSEL _vs(total_pop, nlay, tref, -1, 1);
		
		double **W;  _vs.InitW(mem1, W, nlay);
		_vs.InitC(GC, Gtid, Nh, nlay);
		_vs.GetC(Gtw, Gtid, GC, nlay, W, -1);

		for (int i = 0; i < nlay; ++i)
			DF[i] = Round(i > 0 ? (GC[i * nlay + 0] - GC[(i - 1) * nlay + 0]) : (GC[i * nlay + 0]));
		DF[nlay] = GC[(nlay - 1) * nlay + 0];

		//Calculate genetic distance
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			GENOTYPE *gtab = GetLoc(l).GetGtab();
			GENO_ITERATOR r1(0u, l, true), r1b;

			for (int i1 = 0, h1 = 0; i1 < nind; ++i1)
			{
				r1b = r1;
				POP *p1 = apops[ainds[i1]->popid];
				ushort *als1 = gtab[r1.Read()].GetAlleleArray();

				for (int j1 = 0, v1 = ainds[i1]->vmin; j1 < v1; ++j1, ++h1)
				{
					ushort a1 = als1[j1];
					GENO_ITERATOR r2 = r1b;

					for (int i2 = i1, h2 = h1 + 1; i2 < nind; ++i2)
					{
						POP *p2 = apops[ainds[i2]->popid];
						ushort *als2 = gtab[r2.Read()].GetAlleleArray();

						for (int j2 = i2 == i1 ? j1 + 1 : 0, v2 = ainds[i2]->vmin; j2 < v2; ++j2, ++h2)
						{
							ushort a2 = als2[j2];

							double tdist = 0;
							if (isiam)
							{
								if (a1 == 0xFFFF && a2 == 0xFFFF)	tdist = 1 - SumProd(p1->GetFreq(l), p2->GetFreq(l), GetLoc(l).k);
								else if (a1 == 0xFFFF)				tdist = 1 - p1->GetFreq(l, a2);
								else if (a2 == 0xFFFF)				tdist = 1 - p2->GetFreq(l, a1);
								else if (a1 != a2)					tdist = 1;
							}
							else
							{
								if (a1 == 0xFFFF && a2 == 0xFFFF)	tdist = SumProdSMM(GetLoc(l).GetAlenArray(), p1->GetFreq(l), p2->GetFreq(l), GetLoc(l).k);
								else if (a1 == 0xFFFF)				tdist = SumProdSMM(GetLoc(l).GetAlenArray(), p1->GetFreq(l), a2, GetLoc(l).k);
								else if (a2 == 0xFFFF)				tdist = SumProdSMM(GetLoc(l).GetAlenArray(), p2->GetFreq(l), a1, GetLoc(l).k);
								else if (a1 != a2)					tdist = GetLoc(l).GetSMMDist(a1, a2);
							}

							AtomicAddD(distW(h1, h2), tdist);
						}
					}
				}
			}
			PROGRESS_VALUE += 10;
		}

		//Calculate SS within each vessel
#pragma omp parallel  num_threads(g_nthread_val)
		{
			threadid = omp_get_thread_num();
			VESSEL_ITERATOR ve1(0, _vs, nlay), ve2(0, _vs, nlay);
			for (int i = 0; i < Nh; ++i)
			{
				if (i % g_nthread_val != threadid) { ve1.Next(nlay); continue; }

				int hid1 = ve1.GetHapId();
				ve2.Copy(ve1, nlay);  ve2.Next(nlay);

				for (int j = i + 1; j < Nh; ++j)
				{
					int hid2 = ve2.GetHapId();
					double tdist = distW(hid2, hid1) = distW(hid1, hid2);

					for (int clay = 1; clay <= nlay; ++clay)
						if (ve1.universal_id[clay] == ve2.universal_id[clay])
							AtomicAddD(SSW[clay - 1][ve1.universal_id[clay]], tdist * W[clay][ve1.universal_id[clay]]);

					ve2.Next(nlay);
				}
				ve1.Next(nlay);

				PROGRESS_VALUE += Nh - i - 1;
			}
		}

		//Calculate initial variance components and F-statistics
		{
			VESSEL_ITERATOR ve1(0, _vs, nlay), ve2(0, _vs, nlay);
			_vs.InitC(GC, Gtid, Nh, nlay);
			_vs.GetC(Gtw, Gtid, GC, nlay, W, -1);
			_vs.GetSSHomo(Gss, dist, Nh, W, nlay, ve1, ve2);

			VESSEL::GetV(GC, Gss, Gv, nlay);
			VESSEL::GetF(Gv, Gf, Gvs, nlay);

			SetVal(SS, Gss, nlay);
			SetVal(V, Gv, nlay);
			SetVal(F, Gf, nlay * nlay);
		}


		for (int fa = 1; fa <= nlay; ++fa)
			for (int fb = fa + 1; fb <= nlay; ++fb)
			{
				int pairid = (fa - 1) * nlay + (fb - 1);
				VLA_NEW(_vs2, VESSEL, g_nthread_val);

				for(int i = 0; i < g_nthread_val; ++i)
					new(&_vs2[i]) VESSEL(_vs);

				//Begin permutations
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(static, 1)
				for (int64 m = 0; m < nperm; ++m)
				{
					threadid = omp_get_thread_num();
					double*	ss  = Gss  + nlay * threadid;
					double *v   = Gv   + nlay * threadid;
					double *vs  = Gvs  + nlay * threadid;
					int*	tid = Gtid + (nlay + 1) * threadid;
					double *tw  = Gtw  + (nlay + 1) * threadid;
					double *f   = Gf   + nlay * nlay * threadid;
					double *C   = GC   + nlay * nlay * threadid;
					int*	g   = Gg   + nlay * nlay * threadid;
					int*	e   = Ge   + nlay * nlay * threadid;
					double *ef  = Gef  + nlay * nlay * threadid;
					double *ef2 = Gef2 + nlay * nlay * threadid;
					VESSEL **permbuf = Gpermbuf + Nh * threadid;

					RNG rng(HashULong(((uint64)g_seed_val << 32 | 'AMHO') ^ (pairid * nperm + m)));
					amova_memory = &Gmem2[threadid];	amova_memory->ClearMemory();
					VESSEL &vs2 = _vs2[threadid];
					VESSEL_ITERATOR ve1(0, vs2, nlay), ve2(0, vs2, nlay);

					vs2.Shuffle(rng, fa, fb, 1, permbuf);
					vs2.InitC(C, tid, Nh, nlay);
					vs2.GetC(tw, tid, C, nlay, W, -1);
					vs2.GetSSHomo(ss, dist, Nh, W, nlay, ve1, ve2);

					VESSEL::GetV(C, ss, v, nlay);
					VESSEL::GetF(v, f, vs, nlay);

					if      (f[pairid] > F[pairid] + 1e-7) g[pairid]++;
					else if (f[pairid] > F[pairid] - 1e-7) e[pairid]++;
					ef [pairid] += f[pairid];
					ef2[pairid] += f[pairid] * f[pairid];

					PROGRESS_VALUE += 100;
				}

				//Sum results of multiple threads
				for (int i = 0; i < g_nthread_val; ++i)
				{
					int *g = Gg + nlay * nlay * i;
					int *e = Ge + nlay * nlay * i;
					double *ef = Gef + nlay * nlay * i;
					double *ef2 = Gef2 + nlay * nlay * i;

					G[pairid]	+= g[pairid];
					E[pairid]	+= e[pairid];
					EF[pairid]  += ef[pairid];
					EF2[pairid] += ef2[pairid];
				}

				VLA_DELETE(_vs2);
			}
		
		VLA_DELETE(Gss);
		VLA_DELETE(Gv);
		VLA_DELETE(Gvs);
		VLA_DELETE(Gtid);
		VLA_DELETE(Gtw);
		VLA_DELETE(Gf);
		VLA_DELETE(GC);
		VLA_DELETE(Gg);
		VLA_DELETE(Ge);
		VLA_DELETE(Gef);
		VLA_DELETE(Gef2);
		VLA_DELETE(Gmem2);
		VLA_DELETE(Gpermbuf);
		delete[] dist;
	}

	/* Perform AMOVA using anisoploid method */
	TARGET void AMOVA::CalcAMOVA_aniso()
	{
		//Anisoploids, sum SS and C across locus
		method = amova_cmethod_val;
		Lind = amova_cind_val == 1;
		nlay = Lind + 2 + lreg;

		bool isiam = amova_cmutation_val == 1;
		if (!isiam && abs(g_format_val) <= 2)
			Exit("\nError: Stepwise mutation model (smm) in AMOVA can only be applied for non-vcf input file, and should use size as allele identifier. \n");

		//Number of pseudo permuations at each locus
		int64 M     = amova_pseudo_val * (amova_test_val == 1);
		int64 nperm = amova_nperm_val  * (amova_test_val == 1);
		bool PseudoPerm = amova_pseudo_val > 0;
		if (!PseudoPerm) M = nperm;

		int npair = BINOMIAL[nlay][2];
		double *SSL, *CL;

		//Approach 1: cache results (SS,C) for each locus
		//and each permute will randomly select one result out of M results at each locus and take their sum
#define SSLWLocus(lid,pair,mid) SSL[((lid) * npair * M + (pair) * M + (mid)) * nlay]
#define  CLWLocus(lid,pair,mid)  CL[((lid) * npair * M + (pair) * M + (mid)) * nlay * nlay]
		bool CacheLocus = PseudoPerm && (nloc * M < nperm);
		if (CacheLocus)
		{
			//[L][npair][M]
			//l * npair * M + pairid * M
			int npt = nloc * npair * M * nlay;
			SSL = new double[npt];										SetZero(SSL, npt);
			CL  = new double[npt * nlay];								SetZero(CL,  npt * nlay);
		}

		//Approach 2: cache results (SS,C) for each permute
		//each locus have M results, randomly distriubte to amova_nperm_val permutes
#define SSLWPerm(pair,perm) SSL[((pair) * nperm + (perm)) * nlay]
#define  CLWPerm(pair,perm)  CL[((pair) * nperm + (perm)) * nlay * nlay]
		bool CachePerm = !CacheLocus;
		if (CachePerm)
		{
			//[npair][nperm]
			//pairid * nperm
			int npt = npair * nperm * nlay;
			SSL = new double[npt];										SetZero(SSL, npt);
			CL  = new double[npt * nlay];								SetZero(CL,  npt * nlay);
		}

		//Allocate single thread memory
		VLA_NEW(C, double, nlay * nlay);								SetZero(C, nlay * nlay);

		//Allocate memory for each thread
		int Nht = total_pop->nhaplotypes;
		VESSEL **Gpermbuf = new VESSEL * [Nht * g_nthread_val];			SetZero(Gpermbuf, Nht * g_nthread_val);

		VLA_NEW(Gss,  double, nlay * g_nthread_val);					SetZero(Gss,  nlay * g_nthread_val);
		VLA_NEW(Gv,   double, nlay * g_nthread_val);					SetZero(Gv,   nlay * g_nthread_val);
		VLA_NEW(Gvs,  double, nlay * g_nthread_val);					SetZero(Gvs,  nlay * g_nthread_val);
		VLA_NEW(Gf,   double, nlay * nlay * g_nthread_val);				SetZero(Gf,   nlay * nlay * g_nthread_val);
		VLA_NEW(Gc,   double, nlay * nlay * g_nthread_val);				SetZero(Gc,   nlay * nlay * g_nthread_val);
		VLA_NEW(Gtid, int,    (nlay + 1) * g_nthread_val);				SetZero(Gtid, (nlay + 1) * g_nthread_val);
		VLA_NEW(Gtw,  double, (nlay + 1) * g_nthread_val);				SetZero(Gtw,  (nlay + 1) * g_nthread_val);
		VLA_NEW(Gg,   int,    nlay * nlay * g_nthread_val);				SetZero(Gg,   nlay * nlay * g_nthread_val);
		VLA_NEW(Ge,   int,    nlay * nlay * g_nthread_val);				SetZero(Ge,   nlay * nlay * g_nthread_val);
		VLA_NEW(Gef,  double, nlay * nlay * g_nthread_val);				SetZero(Gef,  nlay * nlay * g_nthread_val);
		VLA_NEW(Gef2, double, nlay * nlay * g_nthread_val);				SetZero(Gef2, nlay * nlay * g_nthread_val);
		VLA_NEW(Gmem1, MEMORY, g_nthread_val);
		double *Gmss = NULL, *Gmc = NULL;

		if (CachePerm && PseudoPerm)
		{
			Gmss = new double[nlay * M * g_nthread_val];				SetZero(Gmss, nlay * M * g_nthread_val);
			Gmc  = new double[nlay * nlay * M * g_nthread_val];			SetZero(Gmc, nlay * nlay * M * g_nthread_val);
		}

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			threadid = omp_get_thread_num();

			int nh = total_pop->loc_stat[l].nhaplo;
			MEMORY &mem1 = Gmem1[threadid];
			amova_memory = &mem1; amova_memory->ClearMemory();

			double *ss  = Gss  + nlay * threadid;
			double *v   = Gv   + nlay * threadid;
			double *vs  = Gvs  + nlay * threadid;
			double *f   = Gf   + nlay * nlay * threadid;
			double *c   = Gc   + nlay * nlay * threadid;
			int *  tid  = Gtid + (nlay + 1) * threadid;
			double *tw  = Gtw  + (nlay + 1) * threadid;
			double *mss = Gmss + nlay * M * threadid;
			double *mc  = Gmc  + nlay * nlay * M * threadid;
			VESSEL **permbuf = Gpermbuf + Nht * threadid;

			int tref = 0;
			VESSEL _vs(total_pop, nlay, tref, l, 1);
			double **W;  _vs.InitW(mem1, W, nlay);

			_vs.InitC(c, tid, nh, nlay);
			_vs.GetC(tw, tid, c, nlay, W, -1);

			for (int i = 0; i < nlay; ++i)
				AtomicAddD(DF[i], Round(i > 0 ? (c[i * nlay + 0] - c[(i - 1) * nlay + 0]) : (c[i * nlay + 0])));
			AtomicAddD(DF[nlay], c[(nlay - 1) * nlay + 0]);

			//Calculate SS within each vessel
			{
				VESSEL_ITERATOR ve1(0, _vs, nlay), ve2(0, _vs, nlay);
				GENOTYPE *gtab = GetLoc(l).GetGtab();
				GENO_ITERATOR r1(0, l, true), r1b;

				for (int i1 = 0, h1 = 0; i1 < nind; ++i1)
				{
					r1b = r1;
					POP *p1 = apops[ainds[i1]->popid];
					GENOTYPE &g1 = gtab[r1.Read()];
					ushort *als1 = g1.GetAlleleArray();

					for (int j1 = 0, v1 = g1.Ploidy(), na1 = g1.Nalleles(); j1 < v1; ++j1, ++h1)
					{
						if (na1 == 0) continue; 

						ushort a1 = als1[j1];
						GENO_ITERATOR r2 = r1b;
						ve2.Copy(ve1, nlay); ve2.Next(nlay);

						for (int i2 = i1, h2 = h1 + 1; i2 < nind; ++i2)
						{
							POP *p2 = apops[ainds[i2]->popid];
							GENOTYPE &g2 = gtab[r2.Read()];
							ushort *als2 = g2.GetAlleleArray();

							for (int j2 = i2 == i1 ? j1 + 1 : 0, v2 = g2.Ploidy(), na2 = g2.Nalleles(); j2 < v2; ++j2, ++h2)
							{
								if (na2 == 0) continue; 
								ushort a2 = als2[j2];
								double tdist = isiam ? a1 != a2 : GetLoc(l).GetSMMDist(a1, a2);

								for (int clay = 1; clay <= nlay; ++clay)
									if (ve1.universal_id[clay] == ve2.universal_id[clay])
										AtomicAddD(SSW[clay - 1][ve1.universal_id[clay]], tdist * W[clay][ve1.universal_id[clay]]);

								ve2.Next(nlay);
							}
						}
						ve1.Next(nlay);
					}
				}
			}

			//Calculate variance components and F-statistics
			{
				VESSEL &vs1 = _vs;
				VESSEL_ITERATOR ve1(0, vs1, nlay), ve2(0, vs1, nlay);

				vs1.InitC(c, tid, nh, nlay);
				vs1.GetC(tw, tid, c, nlay, W, -1);
				vs1.GetSSAniso(ss, isiam, nh, l, W, nlay, ve1, ve2);

				AtomicAddD(SS, ss, nlay);
				AtomicAddD(C, c, nlay * nlay);
			}
						
			//Calculate SS and C of each pseudo permutation at this locus
			for (int fa = 1, pairid = 0; fa <= nlay; ++fa)
			{
				for (int fb = fa + 1; fb <= nlay; ++fb, ++pairid)
				{
					VESSEL vs2(_vs);
					VESSEL_ITERATOR ve1(0, vs2, nlay), ve2(0, vs2, nlay);

					double *ssl, *cl;

					if (CacheLocus)
					{
						ssl = &SSLWLocus(l, pairid, 0);	 //SSL + (l * npair * M + pairid * M) * nlay;
						cl  =  &CLWLocus(l, pairid, 0);  // CL + (l * npair * M + pairid * M) * nlay * nlay;
					}

					if (CachePerm)
					{
						ssl = &SSLWPerm(pairid, 0);		 //SSL + (pairid * nperm) * nlay;
						cl =   &CLWPerm(pairid, 0);		 // CL + (pairid * nperm) * nlay * nlay;
					}

					for (int64 m2 = 0; m2 < M; ++m2)
					{
						RNG rng(HashULong(((uint64)g_seed_val << 32 | 'AMAN') ^ (pairid * nloc * M + l * M + m2)));

						vs2.Shuffle(rng, fa, fb, 2, permbuf);
						vs2.InitC(c, tid, nh, nlay);
						vs2.GetC(tw, tid, c, nlay, W, -1);
						vs2.GetSSAniso(ss, isiam, nh, l, W, nlay, ve1, ve2);

						//Save results for M pseudo permuations and will be add to real permuations after all loci are finished
						if (CacheLocus)
						{
							Add(ssl + nlay * m2, ss, nlay);//locus specific, do not need atomic operations
							Add(cl +  nlay * nlay * m2, c, nlay * nlay);
						}

						//Save results for M pseudo permuations and will be distributed to nperm real permuations
						if (CachePerm && PseudoPerm)
						{
							Add(mss + nlay * m2, ss, nlay);//thread specific, do not need atomic operations
							Add(mc + nlay * nlay * m2, c, nlay * nlay);
						}

						//Save results for nperm real permuations
						if (CachePerm && !PseudoPerm)
						{
							AtomicAddD(ssl + nlay * m2, ss, nlay);
							AtomicAddD(cl  + nlay * nlay * m2, c, nlay * nlay);
						}

						PROGRESS_VALUE += 50;
					}

					//Distribute M pseudo perms into nperm real perms
					if (CachePerm && PseudoPerm)
					{
						for (int64 m = 0; m < nperm; ++m)
						{
							uint data[4] = { (uint)g_seed_val, (uint)pairid, (uint)l, (uint)m};
							int64 m2 = HashString((char*)data, 4 * sizeof(uint)) % M;

							AtomicAddD(ssl + nlay * m,		  mss + nlay * m2,        nlay);
							AtomicAddD(cl  + nlay * nlay * m, mc + nlay * nlay * m2,  nlay * nlay);
						}
					}
				}
			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//Sum SS and C across loci for each permuation
		VESSEL::GetV(C, SS, V, nlay);
		VESSEL::GetF(V, F, Gvs, nlay);

		VLA_NEW(GpermSS, double, nlay * g_nthread_val);
		VLA_NEW(GpermC, double, nlay * nlay* g_nthread_val);

		for (int fa = 1, pairid = 0; fa <= nlay; ++fa)
		{
			for (int fb = fa + 1; fb <= nlay; ++fb, ++pairid)
			{
				int idx = (fa - 1) * nlay + (fb - 1);
				double * clPerm =  &CLWPerm(pairid, 0);		// CL + pairid * (nperm * nlay * nlay);
				double *sslPerm = &SSLWPerm(pairid, 0);		//SSL + pairid * (nperm * nlay);

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
				for (int64 m = 0; m < nperm; ++m)
				{
					threadid = omp_get_thread_num();

					double *v = Gv + nlay * threadid;
					double *f = Gf + nlay * nlay * threadid;
					double *vs = Gvs + nlay * threadid;
					int *g = Gg + nlay * nlay * threadid;
					int *e = Ge + nlay * nlay * threadid;
					double *ef = Gef + nlay * nlay * threadid;
					double *ef2 = Gef2 + nlay * nlay * threadid;

					if (CacheLocus)
					{
						double *permSS = GpermSS + nlay * threadid;
						double *permC = GpermC + nlay * nlay * threadid;
						SetZero(permSS, nlay);
						SetZero(permC, nlay * nlay);

						//Random sample one SS and C out of M results at each locus and sum them across loci
						for (int64 l = 0; l < nloc; ++l)
						{
							uint data[4] = { (uint)g_seed_val, (uint)pairid, (uint)l, (uint)m };
							int64 m2 = HashString((char*)data, 4 * sizeof(uint)) % M;
							double *ssl = &SSLWLocus(l, pairid, 0);//SSL + (l * npair * M + pairid * M) * nlay;
							double *cl =   &CLWLocus(l, pairid, 0);// CL + (l * npair * M + pairid * M) * nlay * nlay;
							Add(permSS, ssl + nlay * m2, nlay);    //thread specific, do not need atomic operations
							Add(permC, cl + nlay * nlay * m2, nlay * nlay);
						}

						VESSEL::GetV(permC, permSS, v, nlay);
					}

					if (CachePerm)
						VESSEL::GetV(clPerm + nlay * nlay * m, sslPerm + nlay * m, v, nlay);

					VESSEL::GetF(v, f, vs, nlay);

					if      (f[idx] > F[idx] + 1e-7)  g[idx]++;
					else if (f[idx] > F[idx] - 1e-7)  e[idx]++;

					ef[idx]  += f[idx];
					ef2[idx] += f[idx] * f[idx];

					PROGRESS_VALUE++;
				}

				//Sum results of multiple threads
				for (int i = 0; i < g_nthread_val; ++i)
				{
					int *g      = Gg + nlay * nlay * i;
					int *e      = Ge + nlay * nlay * i;
					double *ef  = Gef + nlay * nlay * i;
					double *ef2 = Gef2 + nlay * nlay * i;

					G[idx] += g[idx];
					E[idx] += e[idx];
					EF[idx] += ef[idx];
					EF2[idx] += ef2[idx];
				}
			}
		}

		VLA_DELETE(GpermSS);
		VLA_DELETE(GpermC);
		VLA_DELETE(Gss);
		VLA_DELETE(Gv);
		VLA_DELETE(Gvs);
		VLA_DELETE(Gf);
		VLA_DELETE(Gc);
		VLA_DELETE(Gtid);
		VLA_DELETE(Gtw);
		VLA_DELETE(Gg);
		VLA_DELETE(Ge);
		VLA_DELETE(Gef);
		VLA_DELETE(Gef2);
		VLA_DELETE(C);
		VLA_DELETE(Gmem1);
		VLA_DELETE(Gpermbuf);
		delete[] SSL;
		delete[] CL;
		if (CachePerm && PseudoPerm)
		{
			delete[] Gmss;
			delete[] Gmc;
		}
	}

	/* Calculate likelihood for permuated data */
	TARGET double AMOVA::Likelihood(CPOINT &xx, void **Param)
	{
		int clay = *(int*)Param[0];
		int nlay = *(int*)Param[1];
		VESSEL_ITERATOR &ve = *(VESSEL_ITERATOR*)Param[2];

		xx.Image2RealSelfing();
		ve.Rewind(nlay);
		double re = 0, re2 = 1, f = xx.real[0];
		OpenLog(re, re2);

		for (int i = 0; i < nind; ++i)
		{
			IND *ind = ainds[ve.GetHapId()];
			for (int64 l = 0; l < nloc; ++l)
			{
				GENOTYPE &gt = ind->GetGenotype(l);//fine
				if (gt.Nalleles())
				{
					double gfz = gt.GFZ(ve.trace[clay + 1]->GetAlleleCount(l), ve.trace[clay + 1]->nhaplos[l], f);
					ChargeLog(re, re2, gfz);
				}
			}
			ve.Next(nlay);
		}

		CloseLog(re, re2);
		return re;
	}

	/* Perform AMOVA using maximum-likelihood method */
	TARGET void AMOVA::CalcAMOVA_ml()
	{
		MEMORY mem1, mem2;

		method = amova_cmethod_val;
		Lind = amova_cind_val == 1;
		nlay = Lind + 2 + lreg;

		//Dummy haplotype, homoploids
		bool isiam = amova_cmutation_val == 1;
		if (!isiam && abs(g_format_val) <= 2)
			Exit("\nError: Stepwise mutation model (smm) in AMOVA can only be applied for non-vcf input file, and should use size as allele identifier. \n");

		for (int i = 0; i < nind; ++i)
			if (ainds[i]->vmin != ainds[i]->vmax)
				Exit("\nError: Likelihood AMOVA method do not support anisoploids, in individual %s.\n", ainds[i]->name);

		int tref = 0;
		amova_memory = &mem1; amova_memory->ClearMemory();
		VESSEL _vs(total_pop, nlay, tref, -1, 4);
		double **W;  _vs.InitW(mem1, W, nlay);

		int Nht = total_pop->nhaplotypes;
		VESSEL **Gpermbuf = new VESSEL*[Nht * g_nthread_val];		SetZero(Gpermbuf, Nht * g_nthread_val);

		VLA_NEW(Fi,  double,  nlay);								SetZero(Fi,   nlay);
		VLA_NEW(C,   double,  nlay * nlay);							SetZero(C,    nlay * nlay);
		VLA_NEW(Gfi, double,  nlay * g_nthread_val);				SetZero(Gfi,  nlay * g_nthread_val);
		VLA_NEW(Gf,  double,  nlay * nlay * g_nthread_val);			SetZero(Gf,   nlay * nlay * g_nthread_val);
		VLA_NEW(Gc,  double,  nlay * nlay * g_nthread_val);			SetZero(Gc,   nlay * nlay * g_nthread_val);
		VLA_NEW(Gc2, double,  nlay * nlay * g_nthread_val);			SetZero(Gc2,  nlay * nlay * g_nthread_val);
		VLA_NEW(Gtid, int,   (nlay + 1) * g_nthread_val);			SetZero(Gtid, (nlay + 1) * g_nthread_val);
		VLA_NEW(Gtw, double, (nlay + 1) * g_nthread_val);			SetZero(Gtw,  (nlay + 1) * g_nthread_val);
		VLA_NEW(Gg,   int,    nlay * nlay * g_nthread_val);			SetZero(Gg,   nlay * nlay * g_nthread_val);
		VLA_NEW(Ge,   int,    nlay * nlay * g_nthread_val);			SetZero(Ge,   nlay * nlay * g_nthread_val);
		VLA_NEW(Gef,  double, nlay * nlay * g_nthread_val);			SetZero(Gef,  nlay * nlay * g_nthread_val);
		VLA_NEW(Gef2, double, nlay * nlay * g_nthread_val);			SetZero(Gef2, nlay * nlay * g_nthread_val);

		//Calculate original F-statistics
		for (int clay = Lind; clay < nlay; ++clay)
		{
			VESSEL_ITERATOR ve(amova_cind_val == 2 ? 0 : 1, _vs, nlay);
			void *Param[] = { (void*)&clay, (void*)&nlay, (void*)&ve };
			CPOINT xx = CPOINT::DownHillSimplex(1, 0, false, 0.1, 10, AMOVA::Likelihood, Param);
			xx.Image2RealSelfing();
			Fi[clay] = xx.real[0];
		}
		VESSEL::GetF(Fi, F, nlay);

		int nthread = amova_test_val == 1 ? g_nthread_val : 1;
		int64 nperm = amova_nperm_val * (amova_test_val == 1);
		VLA_NEW(_vs2, VESSEL, g_nthread_val);

		for (int fa = 1 + Lind; fa <= nlay; ++fa)
			for (int fb = fa + 1; fb <= nlay; ++fb)
			{
				int pairid = (fa - 1) * nlay + (fb - 1);

				for (int i = 0; i < g_nthread_val; ++i)
					new(&_vs2[i]) VESSEL(_vs);

				//Permute
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(static, 1)
				for (int64 m = 0; m < nperm; ++m)
				{
					RNG rng(HashULong(((uint64)g_seed_val << 32 | 'AMLI') ^ (pairid * nperm + m)));

					threadid = omp_get_thread_num();
					double *fi  = Gfi  + nlay * threadid;
					double *f   = Gf   + nlay * nlay * threadid;
					int *g      = Gg   + nlay * nlay * threadid;
					int *e      = Ge   + nlay * nlay * threadid;
					double *ef  = Gef  + nlay * nlay * threadid;
					double *ef2 = Gef2 + nlay * nlay * threadid;
					VESSEL **permbuf = Gpermbuf + Nht * threadid;

					VESSEL &vs2 = _vs2[threadid];
					VESSEL_ITERATOR ve(amova_cind_val == 2 ? 0 : 1, vs2, nlay);

					vs2.Shuffle(rng, fa, fb, 4, permbuf);
					for (int clay = Lind; clay < nlay; ++clay)
					{
						void *Param[] = { (void*)&clay, (void*)&nlay, (void*)&ve };
						CPOINT xx = CPOINT::DownHillSimplex(1, 0, false, 0.1, 10, AMOVA::Likelihood, Param);
						xx.Image2RealSelfing();
						fi[clay] = xx.real[0];
					}

					VESSEL::GetF(fi, f, nlay);

					if      (f[pairid] > F[pairid] + 1e-7) g[pairid]++;
					else if (f[pairid] > F[pairid] - 1e-7) e[pairid]++;
					ef[pairid] +=	f[pairid];
					ef2[pairid] += f[pairid] * f[pairid];

					PROGRESS_VALUE += 500;
				}

			}
		VLA_DELETE(_vs2);

		if (Lind)
		{
			V[0] = 1 - Fi[nlay - 1];
			V[1] = nlay == 2 ? 1 - V[0] : Fi[1] * (1.0 - Fi[nlay - 1]) / (1.0 - Fi[1]);
			if (nlay > 2) V[nlay - 1] = (Fi[nlay - 1] - Fi[nlay - 2]) / (1.0 - Fi[nlay - 2]);
			for (int clay = 2; clay < nlay - 1; ++clay)
				V[clay] = (Fi[clay] - Fi[clay - 1]) * (1 - Fi[nlay - 1]) / ((1 - Fi[clay - 1]) * (1 - Fi[clay]));
		}
		else
		{
			double VV[100];
			VV[0] = 1 - Fi[nlay - 1];
			VV[1] = nlay == 2 ? 1 - VV[0] : Fi[0] * (1.0 - Fi[nlay - 1]) / (1.0 - Fi[0]);
			if (nlay > 1) VV[nlay] = (Fi[nlay - 1] - Fi[nlay - 2]) / (1.0 - Fi[nlay - 2]);
			for (int clay = 1; clay < nlay - 1; ++clay)
				VV[clay + 1] = (Fi[clay] - Fi[clay - 1]) * (1 - Fi[nlay - 1]) / ((1 - Fi[clay - 1]) * (1 - Fi[clay]));
			for (int i = 0; i < nlay; ++i)
				V[i] = VV[i + 1];
			V[0] += VV[0];
		}

		double SSTOT = 0;
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			double nh = total_pop->loc_stat[l].nhaplo;
			if (nh == 0) continue;

			double *c  = Gc   + nlay * nlay * threadid;
			double *c2 = Gc2  + nlay * nlay * threadid;
			int *tid   = Gtid + (nlay + 1) * threadid;
			double *tw = Gtw  + (nlay + 1) * threadid;

			AtomicAddD(SSTOT, (1.0 - total_pop->loc_stat[l].a2) * 0.5 * nh);
			_vs.InitC(c, tid, total_pop->loc_stat[l].nhaplo, nlay);
			_vs.GetC(tw, tid, c, nlay, W, l);
			Add(c2, c, nlay * nlay);
		}

		for (int i = 0; i < g_nthread_val; ++i)
		{
			Add(C,   Gc2  + nlay * nlay * i, nlay * nlay);
			Add(G,   Gg   + nlay * nlay * i, nlay * nlay);
			Add(E,   Ge   + nlay * nlay * i, nlay * nlay);
			Add(EF,  Gef  + nlay * nlay * i, nlay * nlay);
			Add(EF2, Gef2 + nlay * nlay * i, nlay * nlay);
		}

		for (int i = 0; i < nlay; ++i)
			DF[i] = Round(i > 0 ? (C[i * nlay + 0] - C[(i - 1) * nlay + 0]) : (C[i * nlay + 0]));
		DF[nlay] = C[(nlay - 1) * nlay + 0];

		double sstot = 0;
		for (int i = 0; i < nlay; ++i)
			sstot += C[(nlay - 1) * nlay + i] * V[i];
		double VTOT = SSTOT / sstot;
		Mul(V, VTOT, nlay);

		MatrixMul(C, nlay, nlay, V, nlay, 1, SS);

		VLA_DELETE(Fi);
		VLA_DELETE(C);
		VLA_DELETE(Gfi);
		VLA_DELETE(Gf);
		VLA_DELETE(Gc);
		VLA_DELETE(Gtid);
		VLA_DELETE(Gtw);
		VLA_DELETE(Gg);
		VLA_DELETE(Ge);
		VLA_DELETE(Gef);
		VLA_DELETE(Gef2);
		VLA_DELETE(Gpermbuf);
	}

	/* Destructor */
	TARGET AMOVA::~AMOVA()
	{
		for (int clay = 0; clay < nlay; ++clay)
			delete[] SSW[clay];
		delete[] nSSW;
		delete[] SSW;
		delete[] EF;
		delete[] EF2;
		delete[] G;
		delete[] E;
		delete[] DF;
		delete[] SS;
		delete[] V;
		delete[] F;
	}

	/* Write results */
	TARGET void AMOVA::PrintAMOVA(FILE *fout)
	{
		MEMORY mem1;
		amova_memory = &mem1;
		Lind = amova_cind_val == 1;
		double Vtot = Sum(V, nlay);
		VLA_NEW(SSB, double, nlay + 1);
		VLA_NEW(MS, double, nlay + 1);
		for (int i = 0; i < nlay; ++i)
		{
			SSB[i] = i > 0 ? SS[i] - SS[i - 1] : SS[i];
			MS[i] = SSB[i] / DF[i];
		}
		SSB[nlay] = Sum(SSB, nlay);
		MS[nlay] = SSB[nlay] / DF[nlay];

		switch (method)
		{
		case 1:
			fprintf(fout, "%s%sAMOVA Summary, method: homoploid, mutation model: %s, ind-level=%s",
				g_linebreak_val, g_linebreak_val,
				amova_cmutation_val == 1 ? "IAM" : "SMM",
				Lind == 1 ? "yes" : "no");
			break;
		case 2:
			fprintf(fout, "%s%sAMOVA Summary, method: anisoploid, mutation model: %s, ind-level=%s",
				g_linebreak_val, g_linebreak_val,
				amova_cmutation_val == 1 ? "IAM" : "SMM",
				Lind == 1 ? "yes" : "no");
			break;
		case 3:
			fprintf(fout, "%s%sAMOVA Summary, method: likelihood, mutation model: IAM, ind-level=%s",
				g_linebreak_val, g_linebreak_val,
				Lind == 1 ? "yes" : "no");
			break;
		}

		fprintf(fout, "%sSource%cd.f.%cSS%cMS%cVar%cPercentage", g_linebreak_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val);

		double inv_npermed = 1.0 / (amova_nperm_val * (amova_test_val == 1) + 1);

		VLA_NEW(level_name, char*, nlay);
		VLA_NEW(level_short, char*, nlay);

		if (Lind)
		{
			level_name[0] = (char*)"Individual";
			level_name[1] = (char*)"Population";
			level_name[nlay - 1] = (char*)"Total";
			level_short[0] = (char*)"I";
			level_short[1] = (char*)"S";
			level_short[nlay - 1] = (char*)"T";
			for (int rl = 0; rl < lreg; ++rl)
			{
				level_name[rl + 2] = new char[20];
				level_short[rl + 2] = new char[6];
				sprintf(level_name[rl + 2], "Region Level %d", rl + 1);
				sprintf(level_short[rl + 2], "C%d", rl + 1);
			}
		}
		else
		{
			level_name[0] = (char*)"Population";
			level_name[nlay - 1] = (char*)"Total";
			level_short[0] = (char*)"S";
			level_short[nlay - 1] = (char*)"T";
			for (int rl = 0; rl < lreg; ++rl)
			{
				level_name[rl + 1] = new char[20];
				level_short[rl + 1] = new char[6];
				sprintf(level_name[rl + 1], "Region Level %d", rl + 1);
				sprintf(level_short[rl + 1], "C%d", rl + 1);
			}
		}

		fprintf(fout, "%sWithin %s%c", g_linebreak_val, level_name[0], g_delimiter_val);
		fprintf(fout, "%0.0lf", DF[0]);										fprintf(fout, "%c", g_delimiter_val);
		WriteReal(fout, SSB[0]);											fprintf(fout, "%c", g_delimiter_val);
		WriteReal(fout, MS[0]);												fprintf(fout, "%c", g_delimiter_val);
		WriteReal(fout, V[0]);												fprintf(fout, "%c", g_delimiter_val);
		WriteReal(fout, V[0] * 100.0 / Vtot);								fprintf(fout, "%c", g_delimiter_val);

		for (int i = 1; i < nlay; ++i)
		{
			fprintf(fout, "%sAmong %s%c", g_linebreak_val, level_name[i - 1], g_delimiter_val);
			fprintf(fout, "%0.0lf", DF[i]);									fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, SSB[i]);										fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, MS[i]);											fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, V[i]);											fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, V[i] * 100.0 / Vtot);							fprintf(fout, "%c", g_delimiter_val);
		}
		fprintf(fout, "%sTotal%c", g_linebreak_val, g_delimiter_val);
		fprintf(fout, "%0.0lf", DF[nlay]);									fprintf(fout, "%c", g_delimiter_val);
		WriteReal(fout, SSB[nlay]);											fprintf(fout, "%c", g_delimiter_val);
		WriteReal(fout, MS[nlay]);											fprintf(fout, "%c", g_delimiter_val);
		WriteReal(fout, Vtot);												fprintf(fout, "%c", g_delimiter_val);
		WriteReal(fout, 100.0);												fprintf(fout, "%c", g_delimiter_val);

		fprintf(fout, "%s%sF-statistics", g_linebreak_val, g_linebreak_val);
		fprintf(fout, "%sStatistics%cValue%cPermute Mean%cPermute Var%cPr(rand>obs)%cPr(rand=obs)", g_linebreak_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val);

		for (int i = 0; i < nlay; ++i)
			for (int j = i + 1; j < nlay; ++j)
			{
				fprintf(fout, "%sF%s%s%c", g_linebreak_val, level_short[i], level_short[j], g_delimiter_val);
				WriteReal(fout, F[i * nlay + j]);							fprintf(fout, "%c", g_delimiter_val);
				if (method != 3 || (method == 3 && Lind != 0 && i > 0) || (method == 3 && Lind == 0))
				{
					WriteReal(fout, EF[i * nlay + j] * inv_npermed);		fprintf(fout, "%c", g_delimiter_val);
					WriteReal(fout, EF2[i * nlay + j] * inv_npermed - EF[i * nlay + j] * inv_npermed * EF[i * nlay + j] * inv_npermed);
																			fprintf(fout, "%c", g_delimiter_val);
					WriteReal(fout, G[i * nlay + j] * inv_npermed);			fprintf(fout, "%c", g_delimiter_val);
					WriteReal(fout, E[i * nlay + j] * inv_npermed);			fprintf(fout, "%c", g_delimiter_val);
				}
				else
					fprintf(fout, "-%c-%c-%c-%c", g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val);
			}

		if (method != 3 && amova_printss_val == 1)
		{
			for (int rl = lreg - 1; rl >= 0; --rl)
			{
				int tref = 0;
				VESSEL vs(total_pop, nlay, tref, -1, 1);
				VESSEL_ITERATOR ve1(2 + Lind + rl, vs, nlay);

				fprintf(fout, "%s%sSS within regions level %d %s", g_linebreak_val, g_linebreak_val, rl + 1, g_linebreak_val);
				for (int rl2 = rl; rl2 < lreg; ++rl2)
					fprintf(fout, "RegL%d%c", rl2 + 1, g_delimiter_val);
				fprintf(fout, "SS");
				for (int i = 0; i < nreg[rl]; ++i)
				{
					fprintf(fout, "%s", g_linebreak_val);
					POP *tr = ve1.GetSubpop(nlay, 2 + Lind + rl);
					for (int rl2 = rl; rl2 < lreg; ++rl2)
					{
						fprintf(fout, "%s%c", tr->name, g_delimiter_val);
						tr = aregs[rl2 + 1][tr->rid];
					}
					WriteReal(fout, SSW[Lind + rl + 1][i]);
					ve1.Next(nlay);
				}
			}


			{
				int tref = 0;
				VESSEL vs(total_pop, nlay, tref, -1, 1);
				VESSEL_ITERATOR ve1(1 + Lind, vs, nlay);

				fprintf(fout, "%s%sSS within populations%sPop", g_linebreak_val, g_linebreak_val, g_linebreak_val);
				for (int rl = 0; rl < lreg; ++rl)
					fprintf(fout, "%cRegL%d", g_delimiter_val, rl + 1);
				fprintf(fout, "%cSS", g_delimiter_val);
				for (int i = 0; i < npop; ++i)
				{
					POP *tr = ve1.GetSubpop(nlay, 1 + Lind);
					fprintf(fout, "%s%s%c", g_linebreak_val, tr->name, g_delimiter_val);
					for (int rl = 0; rl < lreg; ++rl)
					{
						tr = aregs[rl][tr->rid];
						fprintf(fout, "%s%c", tr->name, g_delimiter_val);
					}
					WriteReal(fout, SSW[Lind][i]);
					ve1.Next(nlay);
				}
			}

			if (Lind)
			{
				int tref = 0;
				VESSEL vs(total_pop, nlay, tref, -1, 1);
				VESSEL_ITERATOR ve1(1, vs, nlay);

				fprintf(fout, "%s%sSS within individuals%sInd%cPop", g_linebreak_val, g_linebreak_val, g_linebreak_val, g_delimiter_val);
				for (int rl = 0; rl < lreg; ++rl)
					fprintf(fout, "%cRegL%d", g_delimiter_val, rl + 1);
				fprintf(fout, "%cSS", g_delimiter_val);
				for (int i = 0; i < nind; ++i)
				{
					POP *tr = ve1.GetSubpop(nlay, 2);
					IND *ti = ve1.GetInd(nlay, 1);
					fprintf(fout, "%s%s%c%s%c", g_linebreak_val, ti->name, g_delimiter_val, tr->name, g_delimiter_val);
					for (int rl = 0; rl < lreg; ++rl)
					{
						tr = aregs[rl][tr->rid];
						fprintf(fout, "%s%c", tr->name, g_delimiter_val);
					}
					WriteReal(fout, SSW[0][i]);
					ve1.Next(nlay);
				}
			}
		}

		VLA_DELETE(SSB);
		VLA_DELETE(MS);
		for (int rl = 0; rl < lreg; ++rl)
		{
			delete[] level_name[rl + 1 + Lind];
			delete[] level_short[rl + 1 + Lind];
		}
		VLA_DELETE(level_name);
		VLA_DELETE(level_short);
	}
#endif

#ifndef _RELATEDNESS
	/* Write header row for relatedness estimation */
	TARGET /*static*/ void RELATEDNESS::ColumnPrintHeader()
	{
		fprintf(FRES, "%s%s%s%sA%cpop",
			g_linebreak_val, g_linebreak_val,
			cpop->name, g_linebreak_val,
			g_delimiter_val);
		for (int rl = 0; rl < lreg; ++rl)
			fprintf(FRES, "%cregL%d", g_delimiter_val, rl + 1);
		fprintf(FRES, "%cB%cpop", g_delimiter_val, g_delimiter_val);
		for (int rl = 0; rl < lreg; ++rl)
			fprintf(FRES, "%cregL%d", g_delimiter_val, rl + 1);
		fprintf(FRES, "%cAB_typed%cA_typed%cB_typed", g_delimiter_val, g_delimiter_val, g_delimiter_val);

		for (int k = 1; k <= N_RELATEDNESS_ESTIMATOR; ++k)
			if (relatedness_estimator_val[k])
				fprintf(FRES, "%c%s", g_delimiter_val, RELATEDNESS_ESTIMATOR[k]);
	}

	/* Write result row for relatedness estimation */
	TARGET void RELATEDNESS::ColumnPrintLine(int i, int j)
	{
		fprintf(FRES, "%s%s%c%s%c",
			g_linebreak_val,
			ainds[i]->name, g_delimiter_val,
			apops[ainds[i]->popid]->name, g_delimiter_val);

		POP *tr = lreg >= 0 ? aregs[0][apops[ainds[i]->popid]->rid] : NULL;
		for (int rl = 0; rl < lreg; ++rl)
		{
			fprintf(FRES, "%s%c", tr->name, g_delimiter_val);
			tr = aregs[rl + 1][tr->rid];
		}

		fprintf(FRES, "%s%c%s%c",
			ainds[j]->name, g_delimiter_val,
			apops[ainds[j]->popid]->name, g_delimiter_val);

		tr = lreg >= 0 ? aregs[0][apops[ainds[j]->popid]->rid] : NULL;
		for (int rl = 0; rl < lreg; ++rl)
		{
			fprintf(FRES, "%s%c", tr->name, g_delimiter_val);
			tr = aregs[rl + 1][tr->rid];
		}

		fprintf(FRES, "%d%c%d%c%d",
			ABtype, g_delimiter_val,
			Atype, g_delimiter_val,
			Btype);

		for (int k = 1; k <= N_RELATEDNESS_ESTIMATOR; ++k)
			if (relatedness_estimator_val[k])
			{
				fprintf(FRES, "%c", g_delimiter_val);
				WriteReal(FRES, *((&Lynch1999) + k - 1));
			}
	}

	/* Write matrix format header for relatedness estimation */
	TARGET /*static*/ void RELATEDNESS::MatrixPrintMatrixHeader(int k, int n)
	{
		if (relatedness_estimator_val[k] == 0) return;
		fprintf(TEMP_FILES[k], "%s%s%s%s%s", g_linebreak_val, g_linebreak_val, cpop->name, g_linebreak_val, RELATEDNESS_ESTIMATOR[k]);

		for (int i = 0; i < n; ++i)
			fprintf(TEMP_FILES[k], "%c%s", g_delimiter_val, cpop->inds[i]->name);
	}

	/* Write matrix format row header for relatedness estimation */
	TARGET /*static*/ void RELATEDNESS::MatrixPrintRowHeader(int k, int i)
	{
		if (relatedness_estimator_val[k] == 0) return;
		fprintf(TEMP_FILES[k], "%s%s", g_linebreak_val, cpop->inds[i]->name);
	}

	/* Write matrix format grid for relatedness estimation */
	TARGET void RELATEDNESS::MatrixPrintCell(int k)
	{
		if (relatedness_estimator_val[k] == 0) return;
		fprintf(TEMP_FILES[k], "%c", g_delimiter_val);
		WriteReal(TEMP_FILES[k], *((&Lynch1999) + k - 1));
	}

	/* Calculate relatedness coefficient */
	TARGET void RELATEDNESS::CalcRelatedness(IND *x, IND *y)
	{
		SetZero(&Lynch1999, N_RELATEDNESS_ESTIMATOR);
		ABtype = Atype = Btype = 0;
		byte *estimator = relatedness_estimator_val;

		for (int64 l = 0; l < nloc; ++l)
		{
			GENOTYPE &gt1 = x->GetGenotype(l), &gt2 = y->GetGenotype(l);//fine
			if (gt1.Nalleles()) Atype++;
			if (gt2.Nalleles()) Btype++;
			if (gt1.Nalleles() && gt2.Nalleles()) ABtype++;
		}

		if (ABtype == 0)
		{
			for (int k = 1; k <= N_RELATEDNESS_ESTIMATOR; ++k)
				*(&Lynch1999 + k - 1) = NA;
			return;
		}

		for (int k = 1; k <= N_RELATEDNESS_ESTIMATOR; ++k)
			if (estimator[k])
				*(&Lynch1999 + k - 1) = RelatednessEstimator(k, x, y);
	}

	/* Relatedness estimator Warpper */
	TARGET /*static*/ double RELATEDNESS::RelatednessEstimator(int k, IND *x, IND *y)
	{
		switch (k)
		{
		case 1: return R_Lynch1999(x, y);
		case 2: return R_Wang2002(x, y);
		case 3: return R_Thomas2010(x, y);
		case 4: return R_Li1993(x, y);
		case 5: return R_Queller1989(x, y);
		case 6: return R_Huang2016A(x, y);
		case 7: return R_Huang2016B(x, y);
		case 8: return R_Milligan2003(x, y);
		case 9: return R_Anderson2007(x, y, true);
		case 10: return R_Huang2014(x, y);
		case 11: return R_Huang2015(x, y);
		case 12: return R_Ritland1996(x, y, true, true);
		case 13: return R_Loiselle1995(x, y, true, true);
		case 14: return R_Ritland1996(x, y, false, true);
		case 15: return R_Loiselle1995(x, y, false, true);
		case 16: return R_Weir1996(x, y, true);
		}
		return NA;
	}

	/* Lynch 1999 relatedness estimator */
	TARGET /*static*/ double RELATEDNESS::R_Lynch1999(IND *x, IND *y)
	{
		if (x->vmax > 2 || x->vmin < 2)
			Exit("\nError: Lynch1999 relatedness estimators only supports diploids, in individual %s.\n", x->name);
		if (y->vmax > 2 || y->vmin < 2)
			Exit("\nError: Lynch1999 relatedness estimators only supports diploids, in individual %s.\n", y->name);

		double srx = 0, swx = 0, sry = 0, swy = 0;

		for (int64 l = 0; l < nloc; ++l)
		{
			GENOTYPE &gx = x->GetGenotype(l),  &gy = y->GetGenotype(l);//fine
			if (gx.Nalleles() == 0 || gy.Nalleles() == 0) continue;
			LOCSTAT &stat = cpop->loc_stat[l];
			if (stat.k <= 1) continue;

			double *p = cpop->GetFreq(l);
			int a = gx.GetAlleleCopy(0), b = gx.GetAlleleCopy(1), c = gy.GetAlleleCopy(0), d = gy.GetAlleleCopy(1);
			int dab = a == b, dcd = c == d, dbc = b == c, dbd = b == d, dac = a == c, dad = a == d;
			double wx = ((1 + dab) * (p[a] + p[b]) - 4 * p[a] * p[b]) / (2 * p[a] * p[b]);
			double rx = (p[a] * (dbc + dbd) + p[b] * (dac + dad) - 4 * p[a] * p[b]) / ((1 + dab) * (p[a] + p[b]) - 4 * p[a] * p[b]);
			if (IsNormal(rx) && IsNormal(wx))
			{
				srx += wx * rx;
				swx += wx;
			}
			double wy = ((1 + dcd) * (p[c] + p[d]) - 4 * p[c] * p[d]) / (2 * p[c] * p[d]);
			double ry = (p[c] * (dad + dbd) + p[d] * (dac + dbc) - 4 * p[c] * p[d]) / ((1 + dcd) * (p[c] + p[d]) - 4 * p[c] * p[d]);
			if (IsNormal(ry) && IsNormal(wy))
			{
				sry += wy * ry;
				swy += wy;
			}
		}

		return (srx / swx + sry / swy) * 0.5;
	}

	/* Wang 2002 relatedness estimator */
	TARGET /*static*/ double RELATEDNESS::R_Wang2002(IND *x, IND *y)
	{
		if (x->vmax > 2 || x->vmin < 2)
			Exit("\nError: Wang2002 relatedness estimators only supports diploids, in individual %s.\n", x->name);
		if (y->vmax > 2 || y->vmin < 2)
			Exit("\nError: Wang2002 relatedness estimators only supports diploids, in individual %s.\n", y->name);

		double sw = 0, P1 = 0, P2 = 0, P3 = 0, P4 = 0;
		double a2 = 0, a3 = 0, a4 = 0, a22 = 0;

		for (int64 l = 0; l < nloc; ++l)
		{
			GENOTYPE &gx = x->GetGenotype(l),  &gy = y->GetGenotype(l);//fine
			if (gx.Nalleles() == 0 || gy.Nalleles() == 0) continue;
			LOCSTAT &stat = cpop->loc_stat[l];
			if (stat.k <= 1) continue;

			double w = 1.0 / (2 * stat.a2 - stat.a3);
			sw += w;

			a2 += w * stat.a2;
			a3 += w * stat.a3;
			a4 += w * stat.a4;
			a22 += w * stat.a2 * stat.a2;

			int a = gx.GetAlleleCopy(0), b = gx.GetAlleleCopy(1), c = gy.GetAlleleCopy(0), d = gy.GetAlleleCopy(1);

			if ((a == c && b == d) || (a == d && b == c)) P1 += w;
			else if ((a == b || c == d) && (a == c || b == d || a == d || b == c)) P2 += w;
			else if (a == c || a == d || b == c || b == d) P3 += w;
			else P4 += w;
		}

		a2 /= sw; a3 /= sw; a4 /= sw; a22 /= sw;
		P1 /= sw; P2 /= sw; P3 /= sw; P4 /= sw;

		double b = 2 * a22 - a4;//e1
		double c = a2 - 2 * a22 + a4;//e5
		double d = 4 * (a3 - a4);//e2
		double e = 2 * (a2 - 3 * a3 + 2 * a4);//e4
		double f = 8 * (a4 - a3) + 4 * (a2 - a22);//e3
		double g = 1 - 3 * a2 + 2 * a3 - f;//e6
		double V = (1 - b) * (1 - b) * (e * e * f + d * g * g) - (1 - b) * (e * f - d * g) * (e * f - d * g) + 2 * c * d * f * (1 - b) * (g + e) + c * c * d * f * (d + f);
		double phi = ((d * f * ((e + g) * (1 - b) + c * (d + f)) * (P1 - 1)
			+ d * (1 - b) * (g * (1 - b - d) + f * (c + e)) * P3
			+ f * (1 - b) * (e * (1 - b - f) + d * (c + g)) * P2) / V);
		double delta = ((c * d * f * (e + g) * (P1 + 1 - 2 * b)
			+ ((1 - b) * (f * e * e + d * g * g) - (e * f - d * g) * (e * f - d * g)) * (P1 - b)
			+ c * (d * g - e * f) * (d * P3 - f * P2) - c * c * d * f * (P3 + P2 - d - f) - c * (1 - b) * (d * g * P3 + e * f * P2)) / V);

		return delta + phi / 2;
	}

	/* Thomas 2010 relatedness estimator */
	TARGET /*static*/ double RELATEDNESS::R_Thomas2010(IND *x, IND *y)
	{
		if (x->vmax > 2 || x->vmin < 2)
			Exit("\nError: Thomas2010 relatedness estimators only supports diploids, in individual %s.\n", x->name);
		if (y->vmax > 2 || y->vmin < 2)
			Exit("\nError: Thomas2010 relatedness estimators only supports diploids, in individual %s.\n", y->name);

		int t1, t2;
		double sr = 0, sw = 0;//, D = 0, w4 = 0;
		double P[3] = { 0 };

		for (int64 l = 0; l < nloc; ++l)
		{
			GENOTYPE &gx = x->GetGenotype(l),  &gy = y->GetGenotype(l);//fine
			if (gx.Nalleles() == 0 || gy.Nalleles() == 0) continue;
			LOCSTAT &stat = cpop->loc_stat[l];
			if (stat.k <= 1) continue;

			int va = gx.GetAlleleCopy(0), vb = gx.GetAlleleCopy(1), vc = gy.GetAlleleCopy(0), vd = gy.GetAlleleCopy(1);
			t1 = (va == vc) + (vb == vd);
			t2 = (va == vd) + (vb == vc);
			P[0] = P[1] = P[2] = 0;
			P[t1 >= t2 ? 2 - t1 : 2 - t2] = 1;

			double b = 2 * stat.a2 * stat.a2 - stat.a4;
			double c = stat.a2 - 2 * stat.a2 * stat.a2 + stat.a4;
			double d = 4 * (stat.a2 - stat.a2 * stat.a2 - stat.a3 + stat.a4);
			double e = 1 - 5 * stat.a2 + 4 * stat.a2 * stat.a2 + 4 * stat.a3 - 4 * stat.a4;
			double div = c * d + (1 - b) * e;
			double vr = (4 * b * e * e - 4 * (c * d - b * e) * (c * d - b * e) - 4 * c * d * e + d * (1 - b) * (1 - b) - d * d * (1 - b)) / (4 * div * div);

			double rr = ((P[0] - b) * (d * 0.5 + e) + (P[1] - d) * (0.5 - 0.5 * b - c)) / vr / div;
			if (IsNormal(rr) && IsNormal(vr))
			{
				sw += 1 / vr;
				sr += rr;
			}
			//double v4 = (d * c * c + b * e * e - (c * d - e * b) * (c * d - e * b)) / (div * div);
			//D += (e * (P[0] - b) - c * (P[1] - d)) / v4 / div;
			//w4 += 1 / v4;
		}
		//D /= w4;
		return sr / sw;
	}

	/* Li 1993 relatedness estimator */
	TARGET /*static*/ double RELATEDNESS::R_Li1993(IND *x, IND *y)
	{
		if (x->vmax > 2 || x->vmin < 2)
			Exit("\nError: Li1993 relatedness estimators only supports diploids, in individual %s.\n", x->name);
		if (y->vmax > 2 || y->vmin < 2)
			Exit("\nError: Li1993 relatedness estimators only supports diploids, in individual %s.\n", y->name);

		double sr = 0, sw = 0;

		for (int64 l = 0; l < nloc; ++l)
		{
			GENOTYPE &gx = x->GetGenotype(l),  &gy = y->GetGenotype(l);//fine
			if (gx.Nalleles() == 0 || gy.Nalleles() == 0) continue;
			LOCSTAT &stat = cpop->loc_stat[l];
			if (stat.k <= 1) continue;

			int va = gx.GetAlleleCopy(0), vb = gx.GetAlleleCopy(1), vc = gy.GetAlleleCopy(0), vd = gy.GetAlleleCopy(1);

			double a2 = stat.a2, a3 = stat.a3;
			double S0 = 2 * a2 - a3;
			double Sxy = 0;

			if ((va == vc && vb == vd) || (va == vd && vb == vc))								Sxy = 1;
			else if ((va == vb || vc == vd) && (va == vc || vb == vd || va == vd || vb == vc))	Sxy = 0.75;
			else if (va == vc || va == vd || vb == vc || vb == vd)								Sxy = 0.5;
			else																				Sxy = 0;

			double rr = (Sxy - S0) / (1 - S0);
			if (IsNormal(rr)) { sr += (Sxy - S0) / (1 - S0); sw++; }
		}
		return sr / sw;
	}

	/* Queller 1989 relatedness estimator */
	TARGET /*static*/ double RELATEDNESS::R_Queller1989(IND *x, IND *y)
	{
		if (x->vmax > 2 || x->vmin < 2)
			Exit("\nError: Queller1989 relatedness estimators only supports diploids, in individual %s.\n", x->name);
		if (y->vmax > 2 || y->vmin < 2)
			Exit("\nError: Queller1989 relatedness estimators only supports diploids, in individual %s.\n", y->name);

		double srx = 0, swx = 0, sry = 0, swy = 0;

		for (int64 l = 0; l < nloc; ++l)
		{
			GENOTYPE &gx = x->GetGenotype(l),  &gy = y->GetGenotype(l);//fine
			if (gx.Nalleles() == 0 || gy.Nalleles() == 0) continue;
			LOCSTAT &stat = cpop->loc_stat[l];
			if (stat.k <= 2) continue;//cannot be used for diallelic loci

			int va = gx.GetAlleleCopy(0), vb = gx.GetAlleleCopy(1), vc = gy.GetAlleleCopy(0), vd = gy.GetAlleleCopy(1);
			int Sac = va == vc, Sad = va == vd, Sbc = vb == vc, Sbd = vb == vd, Sab = va == vb, Scd = vc == vd;
			double *p = cpop->GetFreq(l);

			double rx = (0.5 * (Sac + Sad + Sbc + Sbd) - p[va] - p[vb]) / (1 + Sab - p[va] - p[vb]);
			if (IsNormal(rx)) { srx += rx; swx++; }
			double ry = (0.5 * (Sac + Sad + Sbc + Sbd) - p[vc] - p[vd]) / (1 + Scd - p[vc] - p[vd]);
			if (IsNormal(ry)) { sry += ry; swy++; }
		}
		return (srx / swx + sry / swy) * 0.5;
	}

	/* Huang 2016 relatedness estimator A */
	TARGET /*static*/ double RELATEDNESS::R_Huang2016A(IND *x, IND *y)
	{
		if (x->vmax > 2 || x->vmin < 2)
			Exit("\nError: Huang2016A relatedness estimators only supports diploids, in individual %s.\n", x->name);
		if (y->vmax > 2 || y->vmin < 2)
			Exit("\nError: Huang2016A relatedness estimators only supports diploids, in individual %s.\n", y->name);

		double c1 = 0, c2 = 0, c3 = 0, c4 = 0, c5 = 0, c6 = 0;
		double S = 0, S2 = 0, sw = 0;

		for (int64 l = 0; l < nloc; ++l)
		{
			GENOTYPE &gx = x->GetGenotype(l),  &gy = y->GetGenotype(l);//fine
			if (gx.Nalleles() == 0 || gy.Nalleles() == 0) continue;
			LOCSTAT &stat = cpop->loc_stat[l];
			if (stat.k <= 1) continue;

			int va = gx.GetAlleleCopy(0), vb = gx.GetAlleleCopy(1), vc = gy.GetAlleleCopy(0), vd = gy.GetAlleleCopy(1);
			double w = 1.0 / (2 * stat.a2 - stat.a3);
			if (IsError(w)) continue;

			sw += w;
			double s = 0;
			if ((va == vc && vb == vd) || (va == vd && vb == vc))								s = 1;
			else if ((va == vb || vc == vd) && (va == vc || vb == vd || va == vd || vb == vc))	s = 0.75;
			else if (va == vc || va == vd || vb == vc || vb == vd)								s = 0.5;
			else																				s = 0;

			S += w * s;
			S2 += w * s * s;

			c1 += w * (1 - 2 * stat.a2 + stat.a3);
			c2 += w * (0.5 * (1 - 2 * stat.a2 + stat.a3));
			c3 += w * (2 * stat.a2 - stat.a3);
			c4 += w * (0.25 * (4 - 4 * stat.a2 - 4 * stat.a2 * stat.a2 - stat.a3 + 5 * stat.a4));
			c5 += w * (0.125 * (2 + 3 * stat.a2 - 8 * stat.a2 * stat.a2 - 7 * stat.a3 + 10 * stat.a4));
			c6 += w * (stat.a2 + stat.a2 * stat.a2 + 0.25 * (stat.a3 - 5 * stat.a4));
		}
		S /= sw; S2 /= sw;
		c1 /= sw; c2 /= sw; c3 /= sw; c4 /= sw; c5 /= sw; c6 /= sw;
		double delta = (c5 * S - c2 * S2 + c2 * c6 - c3 * c5) / (c1 * c5 - c2 * c4);
		double phi = (-c4 * S + c1 * S2 - c1 * c6 + c3 * c4) / (c1 * c5 - c2 * c4);
		return phi / 2 + delta;
	}

	/* Huang 2016 relatedness estimator B */
	TARGET /*static*/ double RELATEDNESS::R_Huang2016B(IND *x, IND *y)
	{
		if (x->vmax > 2 || x->vmin < 2)
			Exit("\nError: Huang2016B relatedness estimators only supports diploids, in individual %s.\n", x->name);
		if (y->vmax > 2 || y->vmin < 2)
			Exit("\nError: Huang2016B relatedness estimators only supports diploids, in individual %s.\n", y->name);

#define RELAT(X) (((2 * c5 - c4) * ((X) - c3) - (2 * c2 - c1) * ((X)*(X) - c6)) / (c1 * c5 - c2 * c4) / 2)

		double c1 = 0, c2 = 0, c3 = 0, c4 = 0, c5 = 0, c6 = 0;
		double srx = 0, swx = 0, sry = 0, swy = 0;

		for (int64 l = 0; l < nloc; ++l)
		{
			GENOTYPE &gx = x->GetGenotype(l),  &gy = y->GetGenotype(l);//fine
			if (gx.Nalleles() == 0 || gy.Nalleles() == 0) continue;
			LOCSTAT &stat = cpop->loc_stat[l];
			if (stat.k <= 1) continue;

			int va = gx.GetAlleleCopy(0), vb = gx.GetAlleleCopy(1), vc = gy.GetAlleleCopy(0), vd = gy.GetAlleleCopy(1);
			double *p = cpop->GetFreq(l);
			double S = 0;
			if ((va == vc && vb == vd) || (va == vd && vb == vc)) S = 1;
			else if ((va == vb || vc == vd) && (va == vc || vb == vd || va == vd || vb == vc)) S = 0.75;
			else if (va == vc || va == vd || vb == vc || vb == vd) S = 0.5;
			else S = 0;

			double pi = 0, pj = 0, px = 0, X = 0, X2 = 0, w = 0;

			if (va != vb)
			{
				pi = p[va]; pj = p[vb]; px = 1 - pi - pj;
				double b = 2 * pi * pj, d = pj * pj + pi * pi, f = 2 * px * (pi + pj), c = (pi + pj) / 2, e = c, g = px;

				c1 = 1 - b - 0.75 * d - 0.5 * f; c4 = 1 - b - 0.5625 * d - 0.25 * f; c2 = c - b + 0.75 * (e - d) + 0.5 * (g - f);
				c5 = c - b + 0.5625 * (e - d) + 0.25 * (g - f); c3 = b + 0.75 * d + 0.5 * f; c6 = b + 0.5625 * d + 0.25 * f;

				double x1 = RELAT(1), x2 = RELAT(0.75), x3 = RELAT(0.5), x4 = RELAT(0);
				X = x1 * b + x2 * d + x3 * f + x4 * (1 - b - d - f);
				X2 = x1 * x1 * b + x2 * x2 * d + x3 * x3 * f + x4 * x4 * (1 - b - d - f);
				w = 1.0 / (X2 - X * X);
			}
			else
			{
				pi = p[va]; px = 1 - pi;
				double b = pi * pi, d = 2 * pi * px, c = (2 * pi * pi) / (2 * pi), e = (2 * pi * px) / (2 * pi);

				c1 = 1 - b - 0.75 * d; c4 = 1 - b - 0.5625 * d; c2 = c - b + 0.75 * (e - d);
				c5 = c - b + 0.5625 * (e - d); c3 = b + 0.75 * d; c6 = b + 0.5625 * d;

				double x1 = RELAT(1), x2 = RELAT(0.75), x3 = RELAT(0);
				X = x1 * b + x2 * d + x3 * (1 - b - d);
				X2 = x1 * x1 * b + x2 * x2 * d + x3 * x3 * (1 - b - d);
				w = 1.0 / (X2 - X * X);
			}

			double rx = w * RELAT(S);
			if (abs(X) < 1e-10 && IsNormal(w) && IsNormal(rx))
			{
				srx += rx;
				swx += w;
			}

			if (vc != vd)
			{
				pi = p[vc]; pj = p[vd]; px = 1 - pi - pj;
				double b = 2 * pi * pj, d = pj * pj + pi * pi, f = 2 * px * (pi + pj), c = (pi + pj) / 2, e = c, g = px;

				c1 = 1 - b - 0.75 * d - 0.5 * f; c4 = 1 - b - 0.5625 * d - 0.25 * f; c2 = c - b + 0.75 * (e - d) + 0.5 * (g - f);
				c5 = c - b + 0.5625 * (e - d) + 0.25 * (g - f); c3 = b + 0.75 * d + 0.5 * f; c6 = b + 0.5625 * d + 0.25 * f;

				double x1 = RELAT(1), x2 = RELAT(0.75), x3 = RELAT(0.5), x4 = RELAT(0);
				X = x1 * b + x2 * d + x3 * f + x4 * (1 - b - d - f);
				X2 = x1 * x1 * b + x2 * x2 * d + x3 * x3 * f + x4 * x4 * (1 - b - d - f);
				w = 1.0 / (X2 - X * X);
			}
			else
			{
				pi = p[vc]; px = 1 - pi;
				double b = pi * pi, d = 2 * pi * px, c = (2 * pi * pi) / (2 * pi), e = (2 * pi * px) / (2 * pi);

				c1 = 1 - b - 0.75 * d; c4 = 1 - b - 0.5625 * d; c2 = c - b + 0.75 * (e - d);
				c5 = c - b + 0.5625 * (e - d); c3 = b + 0.75 * d; c6 = b + 0.5625 * d;

				double x1 = RELAT(1), x2 = RELAT(0.75), x3 = RELAT(0);
				X = x1 * b + x2 * d + x3 * (1 - b - d);
				X2 = x1 * x1 * b + x2 * x2 * d + x3 * x3 * (1 - b - d);
				w = 1.0 / (X2 - X * X);
			}
			double ry = w * RELAT(S);
			if (abs(X) < 1e-10 && IsNormal(w) && IsNormal(ry))
			{
				sry += ry;
				swy += w;
			}
		}
		return 0.5 * (srx / swx + sry / swy);
#undef RELAT
	}

	/* Initialize Anderson 2007 relatedness estimator */
	TARGET /*static*/ void RELATEDNESS::R_AndersonInitialize(IND *x, IND *y)
	{
		if (x->vmax > 2 || x->vmin < 2)
			Exit("\nError: Anderson2007 and Milligan2003 relatedness estimators only supports diploids, in individual %s.\n", x->name);
		if (y->vmax > 2 || y->vmin < 2)
			Exit("\nError: Anderson2007 and Milligan2003 relatedness estimators only supports diploids, in individual %s.\n", y->name);

		SetFF(Anderson2007_Coef, nloc * 3);
		for (int64 l = 0; l < nloc; ++l)
		{
			GENOTYPE &gx = x->GetGenotype(l),  &gy = y->GetGenotype(l);//fine
			if (gx.Nalleles() == 0 || gy.Nalleles() == 0) continue;
			LOCSTAT &stat = cpop->loc_stat[l];
			if (stat.k <= 1) continue;

			double *p = cpop->GetFreq(l);
			int va = gx.GetAlleleCopy(0), vb = gx.GetAlleleCopy(1), vc = gy.GetAlleleCopy(0), vd = gy.GetAlleleCopy(1);
			int ibs = 0;
			double pi = 0, pj = 0, pk = 0, pl = 0;

			if (va == vb) {
				if (vc == vd) {
					if (vc == va)		{ ibs = 1; pi = p[va]; }
					else				{ ibs = 2; pi = p[va]; pj = p[vc]; }
				}
				else {
					if (vc == va)		{ ibs = 3; pi = p[va]; pj = p[vd]; }
					else if (vd == va)	{ ibs = 3; pi = p[va]; pj = p[vc]; }
					else				{ ibs = 4; pi = p[va]; pj = p[vc]; pk = p[vd]; }
				}
			}
			else if (vc == vd) {
				if (vc == va)			{ibs = 5; pi = p[vc]; pj = p[vb]; }
				else if (vc == vb)		{ ibs = 5; pi = p[vc]; pj = p[va]; }
				else					{ ibs = 6; pi = p[vc]; pj = p[va]; pk = p[vb]; }
			}
			else {
				if ((va == vc && vb == vd) || (va == vd && vb == vc)) 
										{ ibs = 7; pi = p[va]; pj = p[vb]; }
				else if (va == vc)		{ ibs = 8; pi = p[va]; pj = p[vb]; pk = p[vd]; }
				else if (vb == vd)		{ ibs = 8; pi = p[vb]; pj = p[va]; pk = p[vc]; }
				else if (va == vd)		{ ibs = 8; pi = p[va]; pj = p[vb]; pk = p[vc]; }
				else if (vb == vc)		{ ibs = 8; pi = p[vb]; pj = p[va]; pk = p[vd]; }
				else					{ ibs = 9; pi = p[va]; pj = p[vb]; pk = p[vc]; pl = p[vd]; }
			}

			switch (ibs)
			{
			case 1:  Anderson2007_Coef[l * 3 + 0] = pi * pi;	 Anderson2007_Coef[l * 3 + 1] = pi * pi * pi;			Anderson2007_Coef[l * 3 + 2] = pi * pi * pi * pi; break;
			case 2:  Anderson2007_Coef[l * 3 + 0] = 0;			 Anderson2007_Coef[l * 3 + 1] = 0;						Anderson2007_Coef[l * 3 + 2] = pi * pj * pi * pj; break;
			case 3:  Anderson2007_Coef[l * 3 + 0] = 0;			 Anderson2007_Coef[l * 3 + 1] = pi * pi * pj;			Anderson2007_Coef[l * 3 + 2] = 2 * pi * pi * pi * pj; break;
			case 4:  Anderson2007_Coef[l * 3 + 0] = 0;			 Anderson2007_Coef[l * 3 + 1] = 0;						Anderson2007_Coef[l * 3 + 2] = 2 * pi * pi * pj * pk; break;
			case 5:  Anderson2007_Coef[l * 3 + 0] = 0;			 Anderson2007_Coef[l * 3 + 1] = pi * pi * pj;			Anderson2007_Coef[l * 3 + 2] = 2 * pi * pi * pi * pj; break;
			case 6:  Anderson2007_Coef[l * 3 + 0] = 0;			 Anderson2007_Coef[l * 3 + 1] = 0;						Anderson2007_Coef[l * 3 + 2] = 2 * pi * pi * pj * pk; break;
			case 7:  Anderson2007_Coef[l * 3 + 0] = 2 * pi * pj; Anderson2007_Coef[l * 3 + 1] = pi * pj * (pi + pj);	Anderson2007_Coef[l * 3 + 2] = 4 * pi * pi * pj * pj; break;
			case 8:  Anderson2007_Coef[l * 3 + 0] = 0;			 Anderson2007_Coef[l * 3 + 1] = pi * pj * pk;			Anderson2007_Coef[l * 3 + 2] = 4 * pi * pi * pj * pk; break;
			case 9:  Anderson2007_Coef[l * 3 + 0] = 0;			 Anderson2007_Coef[l * 3 + 1] = 0;						Anderson2007_Coef[l * 3 + 2] = 4 * pi * pj * pk * pl; break;
			default: break;
			}
		}
	}

	/* Milligan 2003 relatedness estimator */
	TARGET /*static*/ double RELATEDNESS::R_Milligan2003(IND *x, IND *y)
	{
		return R_Anderson2007(x, y, false);
	}

	/* Anderson 2007 relatedness estimator */
	TARGET /*static*/ double RELATEDNESS::R_Anderson2007(IND *x, IND *y, bool confine)
	{
		R_AndersonInitialize(x, y);
		int dim = 2;
		CPOINT xx0 = CPOINT::DownHillSimplex(dim, 0, confine, 0.0001, 15, L_Anderson, NULL);
		xx0.Image2Real();
		return xx0.real[0] + xx0.real[1] / 2;
	}

	/* Calculate Anderson 2007 likelihood */
	TARGET /*static*/ double RELATEDNESS::L_Anderson(CPOINT &x, void **unusued)
	{
		x.Image2Real();
		double *S = &x.real[0];

		double re = 0, re2 = 1;
		OpenLog(re, re2);
		for (int64 l = 0; l < nloc; ++l)
		{
			if (*(uint*)(Anderson2007_Coef + l * 3) == 0xFFFFFFFF) continue;
			ChargeLog(re, re2, SumProd(Anderson2007_Coef + l * 3, S, 3));
		}
		CloseLog(re, re2);

		x.li = re < -1e100 ? -1e100 : re;
		return x.li;
	}

	/* Ritland 1996 kinship estimator, convert into relatedness */
	TARGET /*static*/ double RELATEDNESS::R_Ritland1996(IND *x, IND *y, bool iscorrect, bool mulv)
	{
		double sr = 0, sw = 0;

		for (int64 l = 0; l < nloc; ++l)
		{
			GENOTYPE &gx = x->GetGenotype(l), &gy = y->GetGenotype(l);//fine
			if (gx.Nalleles() == 0 || gy.Nalleles() == 0) continue;
			LOCSTAT &stat = cpop->loc_stat[l];
			if (stat.k <= 1) continue;

			int vx = gx.Ploidy(), vy = gy.Ploidy();
			int k = GetLoc(l).k;
			double *p = cpop->GetFreq(l);
			double tx = -1, ty = -1, txy = -1;

			for (int i = 0; i < k; ++i)
			{
				if (p[i] * stat.nhaplo <= 1e-5) continue;
				double ax = gx.GetFreq(i);
				double ay = gy.GetFreq(i);

				tx  += ax * ax / p[i];
				ty  += ay * ay / p[i];
				txy += ax * ay / p[i];
			}

			int minv = Min(vx, vy);
			double r = iscorrect ? 
				2 * minv * txy * (tx + ty) / (2 * tx * ty * (vx + vy) * stat.a2) :
				(mulv ? minv : 1) * txy;
			double w = iscorrect ? 1.0 / stat.a2 : stat.k - 1;
			if (IsNormal(r))
			{
				sr += r;
				sw += w;
			}
		}
		return sr / sw;
	}

	/* Loiselle 1995 kinship estimator, convert into relatedness */
	TARGET /*static*/ double RELATEDNESS::R_Loiselle1995(IND *x, IND *y, bool iscorrect, bool mulv)
	{
		double sr = 0, sw = 0;
		for (int64 l = 0; l < nloc; ++l)
		{
			GENOTYPE &gx = x->GetGenotype(l),  &gy = y->GetGenotype(l);//fine
			if (gx.Nalleles() == 0 || gy.Nalleles() == 0) continue;
			LOCSTAT &stat = cpop->loc_stat[l];
			if (stat.k <= 1) continue;

			int k = GetLoc(l).k;
			double *p = cpop->GetFreq(l);
			double txy = 0, tb = 0, tx = 0, ty = 0;
			int vx = gx.Ploidy(), vy = gy.Ploidy();

			for (int i = 0; i < k; ++i)
			{
				if (p[i] * stat.nhaplo <= 1e-5) continue;
				double ax = gx.GetFreq(i);
				double ay = gy.GetFreq(i);
				txy += (ax - p[i]) * (ay - p[i]);
				tx += (ax - p[i]) * (ax - p[i]);
				ty += (ay - p[i]) * (ay - p[i]);
				tb += p[i] * (1 - p[i]);
			}

			int minv = Min(vx, vy), maxv = Max(vx, vy);
			double r = iscorrect ? 
				 minv * txy * (tx + ty) / (tx * ty * stat.a2 * (vx + vy)) :
				(mulv ? minv : 1) * txy;
			double w = iscorrect ? 1.0 / stat.a2 : tb;
			if (IsNormal(r) && IsNormal(w))
			{
				sr += r;
				sw += w;
			}
		}
		return sr / sw;
	}

	/* Weir 1996 kinship estimator, convert into relatedness */
	TARGET /*static*/ double RELATEDNESS::R_Weir1996(IND *x, IND *y, bool mulv)
	{
		int N = 0;
		double J = 0, S = 0;

		for (int64 l = 0; l < nloc; ++l)
		{
			GENOTYPE &gx = x->GetGenotype(l),  &gy = y->GetGenotype(l);//fine
			if (gx.Nalleles() == 0 || gy.Nalleles() == 0) continue;
			LOCSTAT &stat = cpop->loc_stat[l];
			if (stat.k <= 1) continue;

			int vx = gx.Ploidy(), vy = gy.Ploidy();
			int k = GetLoc(l).k;
			int minv = Min(vx, vy);
			double *p = cpop->GetFreq(l);

			for (int i = 0; i < k; ++i)
			{
				if (p[i] * stat.nhaplo <= 1e-5) continue;
				S += (mulv ? minv : 1.0) * (gx.GetFreq(i) * gy.GetFreq(i) - p[i] * p[i]);
				J += p[i] * p[i];
			}
			N++;
		}
		return S / (N - J);
	}

	/* Huang 2014 relatedness estimator */
	TARGET /*static*/ double RELATEDNESS::R_Huang2014(IND *x, IND *y)
	{
		double srx = 0, swx = 0, sry = 0, swy = 0;
		if (x->vmax > 8) Exit("\nError: Huang2014 estimator do not support ploidy level > 8, in individual %s.\n", x->name);
		if (y->vmax > 8) Exit("\nError: Huang2014 estimator do not support ploidy level > 8, in individual %s.\n", y->name);

		for (int64 l = 0; l < nloc; ++l)
		{
			GENOTYPE &gx = x->GetGenotype(l),  &gy = y->GetGenotype(l);//fine
			if (gx.Nalleles() == 0 || gy.Nalleles() == 0) continue;
			if (cpop->loc_stat[l].k <= 1) continue;

			double r = 0, w = 0;
			r = HuangMoment(gx, gy, l, w);
			if (IsNormal(r))
			{
				srx += w * r;
				swx += w;
			}

			r = HuangMoment(gy, gx, l, w);
			if (IsNormal(r))
			{
				sry += w * r;
				swy += w;
			}
		}
		return (srx / swx + sry / swy) * 0.5;
	}

	/* Huang 2014 relatedness estimator : similarity index */
	TARGET /*static*/ double RELATEDNESS::S_Index(int *c, int *d, int ploidyx, int ploidyy)
	{
		int a[8] = { 0 };
		int b[8] = { 0 };
		SetVal(a, c, ploidyx);
		SetVal(b, d, ploidyy);
		int S = 0;
		for (int i = 0; i < ploidyx; ++i)
		{
			if (a[i] >= 0)
			{
				for (int j = 0; j < ploidyy; ++j)
				{
					if (a[i] == b[j])
					{
						a[i] = b[j] = 0xFFFFFFFF;
						S++;
						break;
					}
				}
			}
		}
		return S * 1.0 / Max(ploidyx, ploidyy);
	}

	/* Huang 2014 relatedness estimator : get genotype pattern for reference individual */
	TARGET /*static*/ int RELATEDNESS::GetRefMode(int *a, int ploidy)
	{
		int score = 0;
		for (int i = 0; i < ploidy; ++i)
			for (int j = i + 1; j < ploidy; ++j)
				if (a[i] == a[j])
					score++;

		int kind = 1;
		int last = 0;
		int tkind[8], tl = 0;
		int tkind2[8], tl2 = 0;

		for (int i = 1; i < ploidy; ++i)
			if (a[i] != a[i - 1])
			{
				kind++;
				tkind[tl++] = i - last;
				tkind2[tl2++] = a[last];
				last = i;
			}
		tkind[tl++] = ploidy - last;
		tkind2[tl2++] = a[last];

		for (int i = 0; i < kind; ++i)
		{
			for (int j = i + 1; j < kind; ++j)
			{
				if (tkind[i] < tkind[j])
				{
					Swap(tkind[i], tkind[j]);
					Swap(tkind2[i], tkind2[j]);
				}
			}
		}

		int tc = 0;
		for (int i = 0; i < kind; ++i)
			for (int j = 0; j < tkind[i]; ++j)
				a[tc++] = tkind2[i];

		switch (ploidy)
		{
		case 1: switch (score)
		{
		case 0: return 1;
		default: return 0;
		}
		case 2: switch (score)
		{
		case 0: return 1;
		case 1: return 2;
		default: return 0;
		}
		case 3: switch (score)
		{
		case 0: return 1;
		case 1: return 2;
		case 3: return 3;
		default: return 0;
		}
		case 4: switch (score)
		{
		case 0: return 1;
		case 1: return 2;
		case 2: return 3;
		case 3: return 4;
		case 6: return 5;
		default: return 0;
		}
		case 5: switch (score)
		{
		case 0: return 1;
		case 1: return 2;
		case 2: return 3;
		case 3: return 4;
		case 4: return 5;
		case 6: return 6;
		case 10: return 7;
		default: return 0;
		}
		case 6: switch (score)
		{
		case 0: return 1;
		case 1: return 2;
		case 2: return 3;
		case 3: return kind == 3 ? 4 : 5;
		case 4: return 6;
		case 6: return kind == 2 ? 7 : 8;
		case 7: return 9;
		case 10: return 10;
		case 15: return 11;
		default: return 0;
		}
		case 7: switch (score)
		{
		case 0: return 1;
		case 1: return 2;
		case 2: return 3;
		case 3: return kind == 4 ? 4 : 5;
		case 4: return 6;
		case 5: return 7;
		case 6: return kind == 3 ? 8 : 9;
		case 7: return 10;
		case 9: return 11;
		case 10: return 12;
		case 11: return 13;
		case 15: return 14;
		case 21: return 15;
		default: return 0;
		}
		case 8: switch (score)
		{
		case 0: return 1;
		case 1: return 2;
		case 2: return 3;
		case 3: return kind == 5 ? 4 : 6;
		case 4: return kind == 4 ? 5 : 7;
		case 5: return 8;
		case 6: return kind == 4 ? 9 : 11;
		case 7: return kind == 3 ? 10 : 12;
		case 8: return 13;
		case 9: return 14;
		case 10: return 16;
		case 11: return 17;
		case 12: return 15;
		case 13: return 18;
		case 15: return 19;
		case 16: return 20;
		case 21: return 21;
		case 28: return 22;
		default: return 0;
		}
		default: return 0;
		}
	}

	/* Huang 2014 relatedness estimator : calculate relatedness */
	TARGET /*static*/ double RELATEDNESS::HuangMoment(GENOTYPE &gx, GENOTYPE &gy, int64 l, double &weight)
	{
		weight = 0;
		int vx = gx.Ploidy(), vy = gy.Ploidy();
		int minv = Min(vx, vy), maxv = Max(vx, vy);
		int cp = maxv, cpp = maxv + 1;
		int start = cp - minv;

		int xx[8], yy[8];

		double M[64] = { 0 };
		double A[8] = { 0 };
		double E[8] = { 0 };

		double e[81] = { 0 };
		double P[9] = { 0 };
		double Delta[9] = { 0 };
		double Sol[9] = { 0 };
		double Diff[9] = { 0 };

		double E1 = 0, E2 = 0;
		double tmb[8];

		ushort *gxals = gx.GetAlleleArray(), *gyals = gy.GetAlleleArray();

		for (int j = 0; j < cp; ++j)
		{
			tmb[j] = 1.0;
			xx[j] = gxals[j];
			yy[j] = gyals[j];
		}

		int refmode = 10000 * vx + 100 * vx + GetRefMode(xx, vx);

		E[0] = S_Index(xx, yy, vx, vy);
		for (int j = 1; j < cp; ++j)
			E[j] = E[j - 1] * E[0];

		MOMRelatednessAssign(cp, refmode, e, cpop->GetFreq(l), (int*)xx);

		for (int j1 = 0; j1 < cpp; ++j1)
			for (int j2 = 0; j2 < cp; ++j2)
				e[j1 * cpp + j2] -= e[j1 * cpp + cp];

		for (int j1 = 0; j1 < cp; ++j1)
		{
			for (int j2 = 0; j2 < cp; ++j2)
				tmb[j2] *= (cp - j2) / (double)cp;

			double t1 = 0;
			for (int j2 = 0; j2 < cp; ++j2)
			{
				double t = 0;
				for (int j3 = 0; j3 < cp; ++j3)
					t += e[j3 * cpp + j2] * tmb[j3];
				M[j1 * cp + j2] += t;
				t1 += e[j2 * cpp + cp] * tmb[j2];
			}
			A[j1] += t1;
		}

		for (int j1 = 0; j1 <= cp; ++j1)
			P[j1] += e[j1 * cpp + cp];

		for (int j1 = 0; j1 < minv; ++j1)
			Diff[j1] = E[j1] - A[j1];

		for (int j1 = minv; j1 < cp; ++j1)
		{
			Diff[j1] = 0;
			for (int j2 = 0; j2 < cp; ++j2)
				M[j1 * cp + j2] = M[j2 * cp + cp - 1 - j1] = 0;
		}

		for (int j1 = minv; j1 < cp; ++j1)
			M[j1 * cp + cp - 1 - j1] = 1;

		if (!SolveEquation((double*)M, Diff, Delta, cp))
			return 0;

		for (int j = 0; j < cp; ++j)
			Delta[cp] += Delta[j] * (cp - j) / (double)cp;

		//calc variance
		for (int j1 = start; j1 < cpp; ++j1)
		{
			double s = (cp - j1) / (double)cp;

			for (int j2 = 0; j2 < minv; ++j2)
				Diff[j2] = mp(s, j2 + 1) - A[j2];
			for (int j2 = minv; j2 < cp; ++j2)
				Diff[j2] = 0;

			MatrixMul((double*)M, cp, cp, (double*)Diff, cp, 1, (double*)Sol);

			Sol[cp] = 0;
			for (int j = 0; j < cp; ++j)
				Sol[cp] += Sol[j] * (cp - j) / (double)cp;

			E1 += P[j1] * Sol[cp];
			E2 += P[j1] * Sol[cp] * Sol[cp];
		}

		weight = 1 / (E2 - E1 * E1);
		if (weight > DOUBLE_OVERFLOW) weight = DOUBLE_OVERFLOW;

		return Delta[cp];
	}

	/* Initialize Huang 2015 relatedness estimator */
	TARGET /*static*/ void RELATEDNESS::Huang2015_Initialize()
	{
		Huang2015_maps = new TABLE<int, Huang2015ENTRY>[9];

		struct ts {
			uint64 pattern;
			uint hash;
		};
		int len[] = { 0, 2, 9, 31, 109, 339, 1043, 2998, 8405 };
		ts* tdata = (ts*)&mlbin_data[0];//151Kib

		for (int p = 1, count = 0; p <= 8; ++p)
		{
			new(&Huang2015_maps[p]) TABLE<int, Huang2015ENTRY>(false, NULL, len[p]);
			TABLE<int, Huang2015ENTRY> &m = Huang2015_maps[p];
			for (int i = 1; i <= len[p]; ++i)
			{
				Huang2015ENTRY tentry = { i, tdata[count].pattern };
				m[tdata[count++].hash] = tentry;
			}
		}
	}

	/* Uninitialize Huang 2015 relatedness estimator */
	TARGET /*static*/ void RELATEDNESS::Huang2015_Uninitialize()
	{
		delete[] Huang2015_maps;
	}

	/* Huang 2015 relatedness estimator */
	TARGET /*static*/ double RELATEDNESS::R_Huang2015(IND *x, IND *y)
	{
		int vx = x->vmax, vy = y->vmax;
		int minv = vx > vy ? vy : vx, maxv = vy > vx ? vy : vx;
		if (vx > 8) 
			Exit("\nError: Huang2015 relatedness estimator do not support ploidy level > 8, in individual %s at locus %s.\n", x->name);
		if (vy > 8) 
			Exit("\nError: Huang2015 relatedness estimator do not support ploidy level > 8, in individual %s at locus %s.\n", y->name);
		if (x->vmax != x->vmax) 
			Exit("\nError: Huang2015 relatedness estimator do not support anisoploids, in individual %s.\n", x->name);
		if (y->vmax != y->vmax) 
			Exit("\nError: Huang2015 relatedness estimator do not support anisoploids, in individual %s.\n", y->name);

		for (int64 l = 0; l < nloc; ++l)
		{
			Huang2015_Coef[l * 9] = -999;

			GENOTYPE &gx = x->GetGenotype(l), &gy = y->GetGenotype(l);//fine
			if (gx.Nalleles() == 0 || gy.Nalleles() == 0) continue;
			if (cpop->loc_stat[l].k <= 1) continue; //monomorphic in current population

			int xx[8], yy[8], alleles[16];
			ushort *xals = gx.GetAlleleArray(), *yals = gy.GetAlleleArray();
			for (int i = 0; i < vx; ++i) xx[i] = xals[i];
			for (int i = 0; i < vy; ++i) yy[i] = yals[i];

			Huang2015ENTRY& entry = Huang2015_maps[maxv][GetHuang2015Hash(xx, yy, maxv)];
			if (entry.ibs == 0) continue;

			Huang2015_MatchAllele(entry.pattern, xx, yy, alleles, maxv);
			MLRelatednessAssign(maxv, cpop->GetFreq(l), alleles, Huang2015_Coef + l * 9, entry.ibs);
		}

		int dim = minv, diff = abs(vx - vy);
		bool confine = vx == vy && vx % 2 == 0;
		CPOINT xx0 = CPOINT::DownHillSimplex(dim, diff, confine, 0.0001, 15, RELATEDNESS::L_Huang2015, NULL);
		xx0.Image2Real();

		double re = 0;
		for (int i = 0; i < minv; ++i)
			re += xx0.real[i] * (minv - i) / maxv;

		return re;
	}

	/* Calculate Huang 2015 likelihood */
	TARGET /*static*/ double RELATEDNESS::L_Huang2015(CPOINT &xx, void **unusued)
	{
		xx.Image2Real();

		double re = 0, re2 = 1;
		OpenLog(re, re2);
		for (int64 l = 0; l < nloc; ++l)
		{
			if (Huang2015_Coef[l * 9] == -999) continue;
			double lt = 0;
			for (int k = 0; k <= xx.dim; ++k)
				lt += Huang2015_Coef[l * 9 + xx.diff + k] * xx.real[k] / BINOMIAL[xx.diff + k][xx.diff];
			ChargeLog(re, re2, lt);
		}
		CloseLog(re, re2);

		xx.li = re < -1e100 ? -1e100 : re;
		return xx.li;
	}

	/* Huang 2015 likelihood estimator: Match genotype-pair pattern and assign alleles */
	TARGET /*static*/ void RELATEDNESS::Huang2015_MatchAllele(int64 pattern, int *gx, int *gy, int *alleles, int p)
	{
		for (int i = 0; i < 16; ++i) 
			alleles[i] = 0;

		int n = Max((int)pattern & 0xF, (int)(pattern >> (p * 4)) & 0xF) + 1;
		int a1[16], a2[16], a22[16], a22n = 0;
		int cx[8] = { 0 }, cy[8] = { 0 };
		for (int i = p - 1; i >= 0; --i)
		{
			cy[i] = pattern & 0xF;
			pattern >>= 4;
		}
		for (int i = p - 1; i >= 0; --i)
		{
			cx[i] = pattern & 0xF;
			pattern >>= 4;
		}
		for (int i = 0; i < n; ++i)
		{
			int a = 0, b = 0;
			for (int j = 0; j < p; ++j)
			{
				if (cx[j] == i) a++;
				if (cy[j] == i) b++;
			}
			a1[i] = (a << 20) | (b << 16) | i;
		}

		for (int i = 0; i < p + p; ++i)
		{
			int ta = i >= p ? gy[i - p] : gx[i];
			//find is ta in a22
			bool flag = false;
			for (int j = 0; j < a22n; ++j)
				if (a22[j] == ta) flag = true;
			if (flag) continue;
			a22[a22n] = ta;
			int a = 0, b = 0;
			for (int j = 0; j < p; ++j)
			{
				if (gx[j] == ta) a++;
				if (gy[j] == ta) b++;
			}
			a2[a22n++] = (a << 20) | (b << 16) | ta;
		}

		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
			{
				if (a1[j] < a1[i]) Swap(a1[i], a1[j]);
				if (a2[j] < a2[i]) Swap(a2[i], a2[j]);
			}

		for (int i = 0; i < n; ++i)
			alleles[a1[i] & 0xF] = a2[i] & 0xFFFF;
	}
#endif

#ifndef _KINSHIP
	/* Write header row for kinship estimation */
	TARGET /*static*/ void KINSHIP::ColumnPrintHeader()
	{
		fprintf(FRES, "%s%s%s%sA%cpop",
			g_linebreak_val, g_linebreak_val,
			cpop->name, g_linebreak_val,
			g_delimiter_val);
		for (int rl = 0; rl < lreg; ++rl)
			fprintf(FRES, "%cregL%d", g_delimiter_val, rl + 1);
		fprintf(FRES, "%cB%cpop", g_delimiter_val, g_delimiter_val);
		for (int rl = 0; rl < lreg; ++rl)
			fprintf(FRES, "%cregL%d", g_delimiter_val, rl + 1);
		fprintf(FRES, "%cAB_typed%cA_typed%cB_typed", g_delimiter_val, g_delimiter_val, g_delimiter_val);

		for (int k = 1; k <= N_KINSHIP_ESTIMATOR; ++k)
			if (kinship_estimator_val[k])
				fprintf(FRES, "%c%s", g_delimiter_val, KINSHIP_ESTIMATOR[k]);
	}

	/* Write result row for kinship estimation */
	TARGET void KINSHIP::ColumnPrintLine(int i, int j)
	{
		fprintf(FRES, "%s%s%c%s%c",
			g_linebreak_val,
			ainds[i]->name, g_delimiter_val,
			apops[ainds[i]->popid]->name, g_delimiter_val);

		POP *tr = lreg >= 0 ? aregs[0][apops[ainds[i]->popid]->rid] : NULL;
		for (int rl = 0; rl < lreg; ++rl)
		{
			fprintf(FRES, "%s%c", tr->name, g_delimiter_val);
			tr = aregs[rl + 1][tr->rid];
		}

		fprintf(FRES, "%s%c%s%c",
			ainds[j]->name, g_delimiter_val,
			apops[ainds[j]->popid]->name, g_delimiter_val);

		tr = lreg >= 0 ? aregs[0][apops[ainds[j]->popid]->rid] : NULL;
		for (int rl = 0; rl < lreg; ++rl)
		{
			fprintf(FRES, "%s%c", tr->name, g_delimiter_val);
			tr = aregs[rl + 1][tr->rid];
		}

		fprintf(FRES, "%d%c%d%c%d",
			ABtype, g_delimiter_val,
			Atype, g_delimiter_val,
			Btype);

		for (int k = 1; k <= N_KINSHIP_ESTIMATOR; ++k)
			if (kinship_estimator_val[k])
			{
				fprintf(FRES, "%c", g_delimiter_val);
				WriteReal(FRES, *((&Ritland1996) + k - 1));
			}
	}

	/* Write matrix format header for kinship estimation */
	TARGET /*static*/ void KINSHIP::MatrixPrintMatrixHeader(int k, int n)
	{
		if (kinship_estimator_val[k] == 0) return;
		fprintf(TEMP_FILES[k], "%s%s%s%s%s", g_linebreak_val, g_linebreak_val, cpop->name, g_linebreak_val, KINSHIP_ESTIMATOR[k]);

		for (int i = 0; i < n; ++i)
			fprintf(TEMP_FILES[k], "%c%s", g_delimiter_val, cpop->inds[i]->name);
	}

	/* Write matrix format row header for kinship estimation */
	TARGET /*static*/ void KINSHIP::MatrixPrintRowHeader(int k, int i)
	{
		if (kinship_estimator_val[k] == 0) return;
		fprintf(TEMP_FILES[k], "%s%s", g_linebreak_val, cpop->inds[i]->name);
	}

	/* Write matrix format grid for kinship estimation */
	TARGET void KINSHIP::MatrixPrintCell(int k)
	{
		if (kinship_estimator_val[k] == 0) return;
		fprintf(TEMP_FILES[k], "%c", g_delimiter_val);
		WriteReal(TEMP_FILES[k], *((&Ritland1996) + k - 1));
	}

	/* Calculate relatedness coefficient */
	TARGET void KINSHIP::CalcKinship(IND *a, IND *b)
	{
		Ritland1996 = Loiselle1995 = Weir1996 = 0;

		ABtype = Atype = Btype = 0;

		byte *estimator = kinship_estimator_val;

		for (int64 l = 0; l < nloc; ++l)
		{
			GENOTYPE &gt1 = a->GetGenotype(l), &gt2 = b->GetGenotype(l);//fine
			if (gt1.Nalleles()) Atype++;
			if (gt2.Nalleles()) Btype++;
			if (gt1.Nalleles() && gt2.Nalleles()) ABtype++;
		}

		if (ABtype == 0)
		{
			Ritland1996 = Loiselle1995 = Weir1996 = NA;
			return;
		}

		if (estimator[1]) Ritland1996 = RELATEDNESS::R_Ritland1996(a, b, false, false);
		if (estimator[2]) Loiselle1995 = RELATEDNESS::R_Loiselle1995(a, b, false, false);
		if (estimator[3]) Weir1996 = RELATEDNESS::R_Weir1996(a, b, false);
	}
#endif

#ifndef _PCOA
	/* Do nothing */
	TARGET PCOA::PCOA()
	{

	}

	/* Destructor */
	TARGET PCOA::~PCOA()
	{
		if (U) { delete[] U; U = NULL; }
		if (V) { delete[] V; V = NULL; }
	}

	/* Perform PCoA */
	TARGET int PCOA::CalcPCoA(int _maxp)
	{
		// performing pcoa
		//http://www.esapubs.org/archive/ecol/E084/011/CAP_UserNotes.pdf
		maxp = _maxp = N - 1 > _maxp ? _maxp : N - 1;

		double *D1 = new double[N * N];
		double *C = new double[N * N];
		SetVal(D1, D, N * N);

		// Maximum normal distance
		double ma = -1e300, ex = 0;
		int npair = 0;
		for (int i = 0; i < N; ++i)
			for (int j = i; j < N; ++j)
			{
				double val = D1[i * N + j];
				if (IsError(val)) continue;

				if (ma < val && val < 1e300)
					ma = val;
				ex += val;
				npair++;
			}
		ex /= npair;

		// Check distances
		for (int i = 0; i < N; ++i)
			for (int j = i; j < N; ++j)
			{
				if (D1[i * N + j] < 0 || IsError(D1[i * N + j]))
					D1[i * N + j] = D1[i * N + j] = ex;
				if (D1[i * N + j] > 1e300)
					D1[i * N + j] = D1[i * N + j] = ma * 1.2;
			}

		// Total variance
		Vt = 0;
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < i; ++j)
				Vt += D1[i * N + j] * D1[i * N + j];
		Vt /= N * (N - 1);

		// Centering
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < N; ++j)
				D1[i * N + j] = -0.5 * D1[i * N + j] * D1[i * N + j];

		SetVal(C, -1.0 / N, N * N);
		for (int i = 0; i < N; ++i)
			C[i * N + i] = 1 - 1.0 / N;

		// Eigen-value decomposition
		Map<MatrixXd> Cm(C, N, N);//symmetric
		Map<MatrixXd> G(D1, N, N);//symmetric
		G = Cm * G * Cm;

		if (N > 2)
		{
			DenseSymMatProd<double> op(G);
			SymEigsSolver<DenseSymMatProd<double>> eigs(op, maxp, Min(maxp * 2, N));

			eigs.init();
			eigs.compute(SortRule::LargestAlge);

			// Retrieve results
			auto u = eigs.eigenvectors();
			auto v = eigs.eigenvalues();

			//sort by eigen values in descending order
			int nz = v.size();
			VLA_NEW(idx, int, nz);
			for (int i = 0; i < nz; ++i)
				idx[i] = i;

			for (int i = 0; i < nz; ++i)
				for (int j = i + 1; j < nz; ++j)
					if (v(idx[i], 0) < v(idx[j], 0))
						Swap(idx[i], idx[j]);

			// calculate number of axises
			maxp = 0;
			for (int i = 0; i < nz; ++i)
			{
				if (v(idx[i], 0) > 0) maxp++;
				else break;
			}
			p = maxp;

			//copy
			U = new double[N * maxp];
			V = new double[N];

			for (int i = 0; i < maxp; ++i)
			{
				V[i] = v(idx[i], 0);
				double s = sqrt(V[i]);
				for (int j = 0; j < N; ++j)
					U[j * maxp + i] = s * u(j, idx[i]);
			}
			VLA_DELETE(idx);
		}
		else if (N == 2)
		{
			double b = abs(G(0, 0));
			maxp = p = 1;
			U = new double[N * maxp];
			V = new double[N];
			V[0] = b * 2;
			U[0 * maxp + 0] = b * sqrt(2);
			U[1 * maxp + 0] = b * sqrt(2);
		}
		else
		{
			maxp = p = 0;
			U = new double[1];
			V = new double[1];
			V[0] = U[0] = 0;
		}

		delete[] D1;
		delete[] C;
		return maxp;
	}

	/* Print PCoA */
	TARGET void PCOA::PrintPCoA(FILE *fout, double *d, int n, int _est, int _type)
	{
		D = d;
		N = n;
		estimator = _est;
		type = _type;
		U = V = 0;

		CalcPCoA(pcoa_dim_val);

		fprintf(fout, "%s%s%s%sTotal variance%c", g_linebreak_val, g_linebreak_val, GD_ESTIMATOR[estimator], g_linebreak_val, g_delimiter_val);
		WriteReal(fout, Vt);

		fprintf(fout, "%sVariance", g_linebreak_val);
		for (int i = 0; i < p; ++i)
		{
			fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, V[i] / (N - 1));
		}

		switch (type)
		{
		case 1: fprintf(fout, "%sInd", g_linebreak_val); break;
		case 2: fprintf(fout, "%sPop", g_linebreak_val); break;
		case 3:
		default:
			fprintf(fout, "%sRegL%d", g_linebreak_val, type - 2); break;
		}

		for (int i = 0; i < p; ++i)
			fprintf(fout, "%cPC%d", g_delimiter_val, i + 1);

		if (type == 1)
		{
			int nn = nind;
			for (int i = 0; i < nn; ++i)
			{
				fprintf(fout, "%s%s", g_linebreak_val, ainds[i]->name);
				for (int j = 0; j < p; ++j)
				{
					fprintf(fout, "%c", g_delimiter_val);
					WriteReal(fout, U[i * maxp + j]);
				}
			}
		}

		if (type == 2)
		{
			int nn = npop;
			for (int i = 0; i < nn; ++i)
			{
				fprintf(fout, "%s%s", g_linebreak_val, apops[i]->name);
				for (int j = 0; j < p; ++j)
				{
					fprintf(fout, "%c", g_delimiter_val);
					WriteReal(fout, U[i * maxp + j]);
				}
			}
		}

		if (type >= 3)
		{
			int nn = nreg[type - 3];
			for (int i = 0; i < nn; ++i)
			{
				fprintf(fout, "%s%s", g_linebreak_val, aregs[type - 3][i]->name);
				for (int j = 0; j < p; ++j)
				{
					fprintf(fout, "%c", g_delimiter_val);
					WriteReal(fout, U[i * maxp + j]);
				}
			}
		}
	}
#endif

#ifndef _HCLUSTERING
	TARGET HCLUSTERING::HCLUSTERING(double *d, IND **obj, int n, int m, MEMORY *_memory)
	{
		memory = _memory;
		method = m;
		nori = n;
		ncur = n;
		dori = d;
		memory->Alloc(dcur, n * n);
		memory->Alloc(dnew, n * n);
		SetVal(dcur, d, n * n);
		SetVal(dcur, 1e300, n, n + 1);

		node.SetSize(n);
		for (int i = 0; i < n; ++i)
		{
			HCLUSTER *tc;
			memory->Alloc(tc, 1);
			tc->isend = true;
			tc->endname = obj[i]->name;
			tc->x = i;
			tc->y = 0;
			tc->idlen = 1;
			memory->Alloc(tc->id, 1);
			tc->id[0] = (ushort)i;
			tc->left = NULL;
			tc->right = NULL;
			node.Push(tc);
		}
		CalcClustering();
	}

	TARGET HCLUSTERING::HCLUSTERING(double *d, POP **obj, int n, int m, MEMORY *_memory)
	{
		memory = _memory;
		method = m;
		nori = n;
		ncur = n;
		dori = d;
		memory->Alloc(dcur, n * n);
		memory->Alloc(dnew, n * n);
		SetVal(dcur, d, n * n);
		SetVal(dcur, 1e300, n, n + 1);

		for (int i = 0; i < n; ++i)
		{
			HCLUSTER *tc;
			memory->Alloc(tc, 1);
			tc->isend = true;
			tc->endname = obj[i]->name;
			tc->x = i;
			tc->y = 0;
			tc->idlen = 1;
			memory->Alloc(tc->id, 1);
			tc->id[0] = (ushort)i;
			tc->left = NULL;
			tc->right = NULL;
			node.Push(tc);
		}
		CalcClustering();
	}

	TARGET HCLUSTERING::~HCLUSTERING()
	{
		node.~LIST();
	}

	TARGET void HCLUSTERING::CalcClustering()
	{
		while (ncur > 1)
		{
			HCLUSTER *tc;
			memory->Alloc(tc, 1);
			int a, b;
			tc->isend = false;
			tc->endname = NULL;
			tc->y = FindMinIdx(a, b) * 0.5;
			tc->x = (node[a]->x + node[b]->x) * 0.5;
			tc->idlen = node[a]->idlen + node[b]->idlen;
			memory->Alloc(tc->id, tc->idlen);
			SetVal(tc->id, node[a]->id, node[a]->idlen);
			SetVal(tc->id + node[a]->idlen, node[b]->id, node[b]->idlen);
			tc->left = node[a];
			tc->right = node[b];
			ReduceMatrix(a, b);

			Swap(dcur, dnew);
			ncur--;

			node[a] = tc;
			node.Erase(b);
		}
	}

	TARGET double HCLUSTERING::FindMinIdx(int &a, int &b)
	{
		double minv = 0;
		int id = GetMinIdx(dcur, ncur * ncur, minv);
		int a1 = id / ncur;
		int b1 = id % ncur;
		a = Min(a1, b1);
		b = Max(a1, b1);
		return minv;
	}

	TARGET void HCLUSTERING::ReduceMatrix(int _a, int _b)
	{
		int nnew = ncur - 1;
		int a = Min(_a, _b), b = Max(_a, _b);

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int i = 0; i < ncur; ++i)
		{
			int i2 = i > b ? i - 1 : i;
			SetVal(dnew + i2 * nnew, dcur + i * ncur, b);
			SetVal(dnew + i2 * nnew + b, dcur + i * ncur + b + 1, ncur - b - 1);
		}

		if (method == 1)
		{
			//NEAREST
			for (int c = 0; c < ncur; ++c)
			{
				if (c == a || c == b) continue;
				int cnew = c > b ? c - 1 : c;
				double dac = dcur[a * ncur + c], dbc = dcur[b * ncur + c];
				dnew[cnew * nnew + a] = dnew[a * nnew + cnew] = Min(dac, dbc);
			}
		}

		if (method == 2)
		{
			//FURTHEST
			for (int c = 0; c < ncur; ++c)
			{
				if (c == a || c == b) continue;
				int cnew = c > b ? c - 1 : c;
				double dac = dcur[a * ncur + c], dbc = dcur[b * ncur + c];
				dnew[cnew * nnew + a] = dnew[a * nnew + cnew] = Max(dac, dbc);
			}
		}

		if (method == 3)
		{
			//UPGMA
			for (int c = 0; c < ncur; ++c)
			{
				if (c == a || c == b) continue;
				int cnew = c > b ? c - 1 : c;
				double dac = dcur[a * ncur + c], dbc = dcur[b * ncur + c];
				int na = node[a]->idlen, nb = node[b]->idlen;
				dnew[cnew * nnew + a] = dnew[a * nnew + cnew] = (dac * na + dbc * nb) / (na + nb);
			}
		}

		if (method == 4)
		{
			//WPGMA
			for (int c = 0; c < ncur; ++c)
			{
				if (c == a || c == b) continue;
				int cnew = c > b ? c - 1 : c;
				double dac = dcur[a * ncur + c], dbc = dcur[b * ncur + c];
				dnew[cnew * nnew + a] = dnew[a * nnew + cnew] = (dac + dbc) * 0.5;
			}
		}

		if (method == 5)
		{
			//UPGMC
			for (int c = 0; c < ncur; ++c)
			{
				if (c == a || c == b) continue;
				int cnew = c > b ? c - 1 : c;
				double dac = dcur[a * ncur + c], dbc = dcur[b * ncur + c], dab = dcur[a * ncur + b];
				int na = node[a]->idlen, nb = node[b]->idlen;
				dnew[cnew * nnew + a] = dnew[a * nnew + cnew] = sqrt((na * dac * dac + nb * dbc * dbc) / (na + nb) - (dab * dab * na * nb) / (na + nb) / (na + nb));
			}
		}

		if (method == 6)
		{
			//WPGMC
			for (int c = 0; c < ncur; ++c)
			{
				if (c == a || c == b) continue;
				int cnew = c > b ? c - 1 : c;
				double dac = dcur[a * ncur + c], dbc = dcur[b * ncur + c], dab = dcur[a * ncur + b];
				dnew[cnew * nnew + a] = dnew[a * nnew + cnew] = sqrt(dac * dac * 0.5 + dbc * dbc * 0.5 - dab * dab * 0.25);
			}
		}

		if (method == 7)
		{
			//WARD
			for (int c = 0; c < ncur; ++c)
			{
				if (c == a || c == b) continue;
				int cnew = c > b ? c - 1 : c;
				double dac = dcur[a * ncur + c], dbc = dcur[b * ncur + c], dab = dcur[a * ncur + b];
				int na = node[a]->idlen, nb = node[b]->idlen, nc = node[c]->idlen;
				dnew[cnew * nnew + a] = dnew[a * nnew + cnew] = sqrt(((na + nc) * dac * dac + (nb + nc) * dbc * dbc - nc * dab * dab) / (na + nb + nc));
			}
		}
	}

	TARGET void HCLUSTERING::PrintClustering(FILE *fout, HCLUSTER *c, double cy)
	{
		if (c == NULL)
		{
			c = node[0];
			if (c->isend)
			{
				fprintf(fout, "%s:", c->endname);
				WriteReal(fout, cy - c->y);
			}
			else
			{
				fprintf(fout, "(");
				PrintClustering(fout, c->left, c->y);
				fprintf(fout, ",");
				PrintClustering(fout, c->right, c->y);
				fprintf(fout, ");%s", g_linebreak_val);
			}
		}
		else if (c->isend)
		{
			fprintf(fout, "%s:", c->endname);
			WriteReal(fout, cy - c->y);
		}
		else
		{
			fprintf(fout, "(");
			PrintClustering(fout, c->left, c->y);
			fprintf(fout, ",");
			PrintClustering(fout, c->right, c->y);
			fprintf(fout, "):");
			WriteReal(fout, cy - c->y);
		}
	}
#endif

#ifndef _STCLUSTER
	/* Get allele frequency array */
	TARGET double *SCLUSTER::GetFreq(int64 l)
	{
		return bucket + allele_freq_offset[l];
	}

	/* Get allele frequency */
	TARGET double SCLUSTER::GetFreq(int64 l, int allele)
	{
		return bucket[allele_freq_offset[l] + allele];
	}

	/* Set allele frequency pointer */
	TARGET void SCLUSTER::SetFreq(double *_bucket)
	{
		bucket = _bucket;
	}
#endif

#ifndef _STRUCTURE

	/* Set all bits to 0 */
	TARGET STRUCTURE::STRUCTURE()
	{
		SetZero(this, 1);
	}

	/* Write results for a run */
	TARGET void STRUCTURE::PrintStructure()
	{
		char filename[FILE_NAME_LEN];
		sprintf(filename, "%s.structure.k=%d_rep=%d_id=%d.txt", g_output_val, K, par2->rep, par2->id);
		FILE *fout = fopen(filename, "wb");

		fprintf(fout, "%s%sParameters:%s", g_linebreak_val, g_linebreak_val, g_linebreak_val);
		fprintf(fout, "Seed=%d%s", rng.seed, g_linebreak_val);
		fprintf(fout, "Model=%s,%s,%s%s",
			admix ? "ADM" : "NOADM",
			locpriori ? "LOC" : "NOLOC",
			fmodel ? "F" : "NOF", g_linebreak_val);

		/*
		fprintf(fout, "#Inds=%d%s", N, g_linebreak_val);
		fprintf(fout, "#Loci=%d%s", L, g_linebreak_val);
		fprintf(fout, "#Location=%d%s", S, g_linebreak_val);
		fprintf(fout, "#K=%d%s", K, g_linebreak_val);
		fprintf(fout, "#Burnin=%d%s", nburnin, g_linebreak_val);
		fprintf(fout, "#Reps=%d%s", nreps, g_linebreak_val);
		fprintf(fout, "#Thinning=%d%s", nthinning, g_linebreak_val);
		fprintf(fout, "#Runs=%d%s", nruns, g_linebreak_val);
		fprintf(fout, "#Admburnin=%d%s", nadmburnin, g_linebreak_val);

		fprintf(fout, "lambda="); WriteReal(fout, lambda); fprintf(fout, "%s", g_linebreak_val);
		fprintf(fout, "stdlambda="); WriteReal(fout, stdlambda); fprintf(fout, "%s", g_linebreak_val);
		fprintf(fout, "maxlambda="); WriteReal(fout, maxlambda); fprintf(fout, "%s", g_linebreak_val);
		fprintf(fout, "inferlambda=%s%s", inferlambda ? "yes" : "no", g_linebreak_val);
		fprintf(fout, "difflambda=%s%s", difflambda ? "yes" : "no", g_linebreak_val);
		fprintf(fout, "diversity=%s%s", diversity == 1 ? "yes" : "no", g_linebreak_val);

		fprintf(fout, "alpha="); WriteReal(fout, alpha); fprintf(fout, "%s", g_linebreak_val);
		fprintf(fout, "inferalpha=%s%s", inferalpha ? "yes" : "no", g_linebreak_val);
		fprintf(fout, "diffalpha=%s%s", diffalpha ? "yes" : "no", g_linebreak_val);
		fprintf(fout, "uniformalpha=%s%s", uniformalpha ? "yes" : "no", g_linebreak_val);
		fprintf(fout, "stdalpha="); WriteReal(fout, stdalpha); fprintf(fout, "%s", g_linebreak_val);
		fprintf(fout, "maxalpha="); WriteReal(fout, maxalpha); fprintf(fout, "%s", g_linebreak_val);
		fprintf(fout, "alphapriora="); WriteReal(fout, alphapriora); fprintf(fout, "%s", g_linebreak_val);
		fprintf(fout, "alphapriorb="); WriteReal(fout, alphapriorb); fprintf(fout, "%s", g_linebreak_val);
		fprintf(fout, "metrofreq=%d%s", metrofreq, g_linebreak_val);

		fprintf(fout, "r="); WriteReal(fout, r); fprintf(fout, "%s", g_linebreak_val);
		fprintf(fout, "maxr="); WriteReal(fout, maxr); fprintf(fout, "%s", g_linebreak_val);
		fprintf(fout, "epsr="); WriteReal(fout, epsr); fprintf(fout, "%s", g_linebreak_val);
		fprintf(fout, "epseta="); WriteReal(fout, epseta); fprintf(fout, "%s", g_linebreak_val);
		fprintf(fout, "epsgamma="); WriteReal(fout, epsgamma); fprintf(fout, "%s", g_linebreak_val);

		fprintf(fout, "pmeanf="); WriteReal(fout, pmeanf); fprintf(fout, "%s", g_linebreak_val);
		fprintf(fout, "pstdf="); WriteReal(fout, pstdf); fprintf(fout, "%s", g_linebreak_val);
		fprintf(fout, "stdf="); WriteReal(fout, stdf); fprintf(fout, "%s", g_linebreak_val);
		fprintf(fout, "fsame=%s%s", inferalpha ? "yes" : "no", g_linebreak_val);
		*/


		par2->MeanlnL = rout[0];
		par2->VarlnL = rout[2];
		par2->lnPD = rout[3];

		fprintf(fout, "Number of cluster=%d%s", K, g_linebreak_val);
		fprintf(fout, "Replicate=%d%s", par2->rep, g_linebreak_val);
		fprintf(fout, "Run id=%d%s", par2->id, g_linebreak_val);
		fprintf(fout, "Mean value of ln likelihood="); WriteReal(fout, par2->MeanlnL); fprintf(fout, "%s", g_linebreak_val);
		fprintf(fout, "Variance of ln likelihood="); WriteReal(fout, par2->VarlnL); fprintf(fout, "%s", g_linebreak_val);
		fprintf(fout, "Estimated Ln Prob of Data="); WriteReal(fout, par2->lnPD); fprintf(fout, "%s%s", g_linebreak_val, g_linebreak_val);


		int nf = fmodel ? (fsame ? 1 : K) : 0;
		int nl = difflambda ? K : 1;

		if (inferlambda)
		{
			if (!difflambda)
			{
				fprintf(fout, "%sMean value of lambda=", g_linebreak_val); WriteReal(fout, rout[4]);
			}
			else for (int j = 0; j < K; ++j)
			{
				fprintf(fout, "%sMean value of lambda_%d=", g_linebreak_val, j + 1); WriteReal(fout, rout[4 + j]);
			}
		}

		if (locpriori)
		{
			fprintf(fout, "%sMean value of r=", g_linebreak_val); WriteReal(fout, rout[4 + nl]);
			fprintf(fout, "%sMean value of global %s", g_linebreak_val, admix ? "alpha" : "eta");
			fprintf(fout, "%s%cCluster%s", g_linebreak_val, g_delimiter_val, g_linebreak_val);
			for (int j = 0; j < K; ++j)
				fprintf(fout, "%c%d", g_delimiter_val, j + 1);
			fprintf(fout, "%s", g_linebreak_val);
			for (int j = 0; j < K; ++j)
			{
				fprintf(fout, "%c", g_delimiter_val); WriteReal(fout, rout[4 + nl + 1 + j]);
			}
			fprintf(fout, "%sMean value of local %s for each location", g_linebreak_val, admix ? "alpha" : "gamma");
			fprintf(fout, "%s%cCluster%sPop", g_linebreak_val, g_delimiter_val, g_linebreak_val);
			for (int j = 0; j < K; ++j)
				fprintf(fout, "%c%d", g_delimiter_val, j + 1);
			for (int i = 0; i < S; ++i)
			{
				fprintf(fout, "%s%s", g_linebreak_val, apops[i]->name);
				for (int j = 0; j < K; ++j)
				{
					fprintf(fout, "%c", g_delimiter_val); WriteReal(fout, rout[4 + nl + 1 + K + i * K + j]);
				}
			}
		}
		else if (admix && inferalpha)
		{
			if (!diffalpha)
			{
				fprintf(fout, "%sMean value of alpha=", g_linebreak_val); WriteReal(fout, rout[4 + nl]);
			}
			else for (int j = 0; j < K; ++j)
			{
				fprintf(fout, "%sMean value of alpha_%d=", g_linebreak_val, j + 1); WriteReal(fout, rout[4 + nl + j]);
			}
		}

		if (fmodel)
		{
			if (fsame)
			{
				fprintf(fout, "%sMean value of Fst=", g_linebreak_val); WriteReal(fout, rout[rlen - 1]);
			}
			else
			{
				fprintf(fout, "%sMean value of Fst", g_linebreak_val);
				fprintf(fout, "%s%cCluster%s", g_linebreak_val, g_delimiter_val, g_linebreak_val);
				for (int j = 0; j < K; ++j)
					fprintf(fout, "%c%d", g_delimiter_val, j + 1);
				fprintf(fout, "%s", g_linebreak_val);
				for (int i = 0; i < nf; ++i)
				{
					fprintf(fout, "%c", g_delimiter_val); WriteReal(fout, rout[rlen - nf + i]);
				}
			}
		}

		VLA_NEW(O, double, S * K);
		SetZero(O, S * K);
		fprintf(fout, "%s%sProportion of membership of each pre-defined population in each of the %d clusters", g_linebreak_val, g_linebreak_val, K);
		fprintf(fout, "%s%cCluster%sPop", g_linebreak_val, g_delimiter_val, g_linebreak_val);
		for (int i = 0; i < N; ++i)
			Add(O + ainds[i]->popid * K, Q + i * K, K);
		Unify(O, S, K);

		for (int k = 0; k < K; ++k)
			fprintf(fout, "%c%d", g_delimiter_val, k + 1);

		for (int s = 0; s < S; ++s)
		{
			fprintf(fout, "%s%s", g_linebreak_val, apops[s]->name);
			for (int k = 0; k < K; ++k)
			{
				fprintf(fout, "%c", g_delimiter_val);
				WriteReal(fout, O[s * K + k]);
			}
		}

		fprintf(fout, "%s%sInferred ancestry of individuals", g_linebreak_val, g_linebreak_val);
		fprintf(fout, "%s%c%cCluster", g_linebreak_val, g_delimiter_val, g_delimiter_val);
		fprintf(fout, "%sInd%cPop", g_linebreak_val, g_delimiter_val);
		for (int j = 0; j < K; ++j)
			fprintf(fout, "%c%d", g_delimiter_val, j + 1);

		for (int i = 0; i < N; ++i)
		{
			fprintf(fout, "%s%s%c%s", g_linebreak_val, ainds[i]->name, g_delimiter_val, apops[ainds[i]->popid]->name);
			for (int k = 0; k < K; ++k)
			{
				fprintf(fout, "%c", g_delimiter_val);
				WriteReal(fout, Z_cumu[i * K + k]);
			}
		}
		VLA_DELETE(O);

		VLA_NEW(H, double, K);
		VLA_NEW(D, double, K * K);
		double *H2 = new double[K * L];

		SetZero(H, K);
		for (int k = 0; k < K; ++k)
		{
			double *p = clusterb[k].bucket;
			for (int64 l = 0; l < L; ++l)
			{
				int k2 = GetLoc(l).k;
				double h = 1 - SumSquare(p, k2);
				H2[k * L + l] = h;
				H[k] += h;
				p += k2;
			}
			H[k] /= L;

			for (int j = 0; j <= k; ++j)
			{
				double *p1 = clusterb[k].bucket;
				double *p2 = clusterb[j].bucket;
				double dt = 0;
				for (int64 l = 0; l < L; ++l)
				{
					int k2 = GetLoc(l).k;
					dt += 1 - SumProd(p1, p2, k2);
					p1 += k2;
					p2 += k2;
				}
				D[j * K + k] = D[k * K + j] = dt / L - (H[k] + H[j]) / 2;
			}
		}

		fprintf(fout, "%s%sAllele-frequency divergence among pops (Net nucleotide distance)", g_linebreak_val, g_linebreak_val);
		fprintf(fout, "%s%cCluster%sCluster", g_linebreak_val, g_delimiter_val, g_linebreak_val);
		for (int j = 0; j < K; ++j)
			fprintf(fout, "%c%d", g_delimiter_val, j + 1);

		for (int i = 0; i < K; ++i)
		{
			fprintf(fout, "%s%d", g_linebreak_val, i + 1);
			for (int j = 0; j < K; ++j)
			{
				fprintf(fout, "%c", g_delimiter_val);
				WriteReal(fout, D[i * K + j]);
			}
		}

		if (structure_diversity_val == 1)
		{
			fprintf(fout, "%s%sEstimated genetic diversity in each cluster", g_linebreak_val, g_linebreak_val);
			fprintf(fout, "%sLocus%cCluster%cHeterozygosity", g_linebreak_val, g_delimiter_val, g_delimiter_val);
			for (int a = 0; a < maxK; ++a)
				fprintf(fout, "%cAllele%d", g_delimiter_val, a);

			for (int k = 0; k < K; ++k)
			{
				fprintf(fout, "%s%s%c%d", g_linebreak_val, k ? "" : "Average", g_delimiter_val, k + 1);
				WriteReal(fout, H[k]);
				buf2[k] = clusterb[k].bucket;
			}

			for (int64 l = 0; l < L; ++l)
			{
				for (int k = 0; k < K; ++k)
				{
					fprintf(fout, "%s%s%c%d%c", g_linebreak_val, k ? "" : GetLoc(l).GetName(), g_delimiter_val, k + 1, g_delimiter_val);
					WriteReal(fout, H2[k * L + l]);
					for (int a = 0; a < GetLoc(l).k; ++a)
					{
						fprintf(fout, "%c", g_delimiter_val);
						WriteReal(fout, *buf2[k]++);
					}
				}
			}
		}

		VLA_DELETE(H);
		VLA_DELETE(D);
		delete[] H2;

		fclose(fout);
		Uninit();
	}

	/* Write results summary for all runs */
	TARGET /*static*/ void STRUCTURE::PrintSummary(SRUNINFO *sp, int len)
	{
		OpenResFile("-structure", "#Bayesian clustering");

		fprintf(FRES, "%s%sID%ck%crep%cMean lnL%cVar lnL%clnP(D)",
			g_linebreak_val, g_linebreak_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val);
		for (int i = 0; i < len; ++i)
		{
			fprintf(FRES, "%s%d%c%d%c%d%c",
				g_linebreak_val, sp[i].id,
				g_delimiter_val, sp[i].k,
				g_delimiter_val, sp[i].rep, g_delimiter_val);
			WriteReal(FRES, sp[i].MeanlnL);
			fprintf(FRES, "%c", g_delimiter_val);
			WriteReal(FRES, sp[i].VarlnL);
			fprintf(FRES, "%c", g_delimiter_val);
			WriteReal(FRES, sp[i].lnPD);
		}

		CloseResFile();
	}

	/* Copy parameters */
	TARGET void STRUCTURE::ReadPar(SRUNINFO *_par2)
	{
		par2 = _par2;
		new (&rng) RNG(HashULong(((uint64)g_seed_val << 32 | 'BAYE') ^ (par2->id)));

		S = npop;

		K = par2->k;
		N = nind;
		L = nloc;

		admix = structure_admix_val == 1;
		locpriori = structure_locpriori_val == 1;
		fmodel = structure_f_val == 1;

		nburnin = structure_nburnin_val;
		nreps = structure_nreps_val;
		nthinning = structure_nthinning_val;
		nruns = structure_nruns_val;
		nadmburnin = structure_nadmburnin_val;

		lambda = structure_lambda_val;
		stdlambda = structure_stdlambda_val;
		maxlambda = structure_maxlambda_val;
		inferlambda = structure_inferlambda_val == 0;
		difflambda = structure_difflambda_val == 1;
		diversity = structure_diversity_val;

		alpha = structure_alpha_val;
		inferalpha = structure_inferalpha_val == 1;
		diffalpha = structure_diffalpha_val == 1;
		uniformalpha = structure_uniformalpha_val == 1;
		stdalpha = structure_stdalpha_val;
		maxalpha = structure_maxalpha_val;
		alphapriora = structure_alphapriora_val;
		alphapriorb = structure_alphapriorb_val;
		metrofreq = structure_metrofreq_val;

		r = structure_r_val;
		maxr = structure_maxr_val;
		epsr = structure_stdr_val;
		epseta = structure_epseta_val;
		epsgamma = structure_epsgamma_val;

		pmeanf = structure_pmeanf_val;
		pstdf = structure_pstdf_val;
		stdf = structure_stdf_val;
		fsame = structure_singlef_val == 1;

		kdis = new int[maxK + 1];
		SetZero(kdis, maxK + 1);

		for (int64 l = 0; l < L; ++l)
			kdis[GetLoc(l).k]++;
	}

	/* Initialize MCMC */
	TARGET void STRUCTURE::InitFix()
	{
		nr = 0;
		r = 1;
		buf = new double[K];
		bufb = new double[K];
		buf2 = new double *[K];

		Z_cumu = new double[N * K];
		SetZero(Z_cumu, N * K);

		Base = new double[(K * 2 + 3) * KT];
		SetZero(Base, (K * 2 + 3) * KT);

		cluster = new SCLUSTER[K];
		clusterb = new SCLUSTER[K];
		for (int k = 0; k < K; ++k)
		{
			cluster[k].SetFreq(Base + k * KT);
			clusterb[k].SetFreq(Base + k * KT + K * KT);
			double *p = cluster[k].bucket;
			for (int64 l = 0; l < L; ++l)
			{
				int k2 = GetLoc(l).k;
				SetVal(p, 1.0 / k2, k2);
				p += k2;
			}
		}

		Ni = new int[K * KT];
		SetZero(Ni, K * KT);

		Mi = new double[N * K];
		SetZero(Mi, N * K);

		Q = new double[N * K];
		SetVal(Q, 1.0 / K, N * K);

		Alpha = new double[K];
		SetVal(Alpha, alpha, K);

		Lambda = new double[K];
		SetVal(Lambda, lambda, K);

		int nf = fmodel ? (fsame ? 1 : K) : 0;
		int na = admix ? (diffalpha ? K : 1) : 0;
		int nl = difflambda ? K : 1;

		if (locpriori) rlen = 4 + nl + 1 + K + K * S + nf;
		else rlen = 4 + nl + na + nf;

		rout = new double[rlen];
		SetZero(rout, rlen);
	}

	/* Initialize MCMC for admix model */
	TARGET void STRUCTURE::InitAdmix()
	{
		// Init Z, Z, Mi and Ni

		//Adm
		Z = new ushort*[N];
		ZZ = new ushort*[N];
		for (int i = 0; i < N; ++i)
			Z[i] = new ushort[ainds[i]->vt];
		SetVal(ZZ, Z, N);

		for (int64 l = 0, o = 0; l < L; ++l, o += GetLoc(l).k)
		{
			GENO_ITERATOR rt(0u, l, true);//20220316
			GENOTYPE *gtab = GetLoc(l).GetGtab();
			double *mi = Mi;

			for (int i = 0; i < N; ++i)
			{
				ushort *&z = ZZ[i];
				GENOTYPE &gt = gtab[rt.Read()];

				if (gt.Nalleles() == 0) continue;
				ushort *als = gt.GetAlleleArray();
				for (int a = 0, vi = gt.Ploidy(); a < vi; ++a)
				{
					ushort k2 = (ushort)rng.Next(K);
					*z++ = k2;
					mi[k2]++;
					Ni[k2 * KT + o + als[a]]++;
				}
				mi += K;
			}
		}

		/*
		double *mi = Mi;
		for (int i = 0; i < N; ++i)
		{
			ushort *z = Z[i];
			for (int64 l = 0, o = 0; l < L; ++l, o += GetLoc(l).k)
			{
				GENOTYPE &gt = ainds[i]->GetGenotype(l);
				if (gt.Nalleles() == 0) continue;

				ushort *als = gt.GetAlleleArray();
				for (int a = 0, vi = gt.Ploidy(); a < vi; ++a)
				{
					ushort k = (ushort)rng.Next(K);
					*z++ = k;
					mi[k]++;
					Ni[k * KT + o + als[a]]++;
				}
			}
			mi += K;
		}
		*/
	}

	/* Initialize MCMC for locpriori model */
	TARGET void STRUCTURE::InitLocPriori()
	{
		if (locpriori)
		{
			SumAlpha = new double[S];
			SetZero(SumAlpha, S);

			AlphaLocal = new double[S * K];
			SetVal(AlphaLocal, alpha, S * K);

			if (!admix)
			{
				Eta = new double[K];
				SetVal(Eta, 1.0 / K, K);

				Gamma = new double[S * K];
				SetVal(Gamma, 1.0 / K, S * K);

				Di = new double[S * K];
				SetZero(Di, S * K);
			}
		}
	}

	/* Initialize MCMC for F model */
	TARGET void STRUCTURE::InitFmodel()
	{
		if (fmodel)
		{
			f = new double[K];
			F = new double[K];

			PA.SetFreq(Base + 2 * K * KT);
			PA1.SetFreq(Base + 2 * K * KT + KT);
			PA2.SetFreq(Base + 2 * K * KT + 2 * KT);

			for (int k = 0; k < K; ++k)
			{
				F[k] = pmeanf;
				f[k] = (1 - F[k]) / F[k];
			}

			for (int k = 0; k < K; ++k)
				SetVal(cluster[k].bucket, Lambda[k], KT);

			double *p = PA.bucket;
			SetVal(p, lambda, KT);
			for (int64 l = 0; l < L; ++l)
			{
				int k2 = GetLoc(l).k;
				double *p2 = total_pop->GetFreq(l);
				double nhaplo = total_pop->loc_stat[l].nhaplo;
				for (int a = 0; a < k2; ++a)
					p[a] = (lambda + p2[a] * nhaplo) / (k2 * lambda + nhaplo);
				p += k2;
			}
		}
	}

	/* Update allele frequency for all clusters */
	TARGET void STRUCTURE::UpdateP()
	{
		//set lambda, dirichlet distribution priori parameter
		if (fmodel)
		{
			if (fsame)
			{
				Mul(cluster[0].bucket, PA.bucket, f[0], KT);
				for (int k = 1; k < K; ++k)
					SetVal(cluster[k].bucket, cluster[0].bucket, KT);
			}
			else for (int k = 0; k < K; ++k)
				Mul(cluster[k].bucket, PA.bucket, f[k], KT);
		}
		else for (int k = 0; k < K; ++k)
			SetVal(cluster[k].bucket, Lambda[k], KT);

		int *ni = Ni;
		for (int k = 0; k < K; ++k)
		{
			double *p = cluster[k].bucket;
			for (int64 l = 0; l < L; ++l)
			{
				int k2 = GetLoc(l).k;
				rng.Dirichlet(p, p, ni, k2);
				p += k2;
				ni += k2;
			}
		}
	}

	/* Update a priori ancetral proportion for non-admix model */
	TARGET void STRUCTURE::UpdateQNoAdmix()
	{
		SetZero(Q, N * K);
		double *q = Q;
		for (int i = 0; i < N; ++i)
		{
			double *p = Base;
			OpenLog(buf, bufb, K);

			//add priori probability
			if (locpriori)
				for (int k = 0; k < K; ++k)
					ChargeLog(buf[k], bufb[k], Gamma[ainds[i]->popid * K + k]);

			for (int64 l = 0; l < L; ++l, p += GetLoc(l).k)
			{
				GENOTYPE &gt = ainds[i]->GetGenotype(l);//fine
				if (gt.Nalleles() == 0) continue;

				ushort *als = gt.GetAlleleArray();
				for (int a = 0, vi = gt.Ploidy(); a < vi; ++a)
					ChargeLog(buf, bufb, p + als[a], K, KT);
			}
			CloseLog(buf, bufb, K);

			int k = rng.PolyLog(buf, K);
			q[k] = 1;
			Z[i][0] = (ushort)k;
			q += K;
		}
		binaryq = true;
	}

	/* Update a priori ancetral proportion for admix model */
	TARGET void STRUCTURE::UpdateQAdmix()
	{
		//used by noadmix (in admburnin) and admix model
		if (locpriori)
			for (int i = 0; i < N; ++i)
				rng.Dirichlet(Q + i * K, AlphaLocal + ainds[i]->popid * K, Mi + i * K, K);
		else
			for (int i = 0; i < N; ++i)
				rng.Dirichlet(Q + i * K, Alpha, Mi + i * K, K);
		binaryq = false;
	}

	/* Update a priori ancetral proportion by Metropolis-Hastings for admix model*/
	TARGET void STRUCTURE::UpdateQMetro()
	{
		double *q = Q;
		for (int i = 0; i < N; ++i)
		{
			if (locpriori)
				rng.Dirichlet(buf, AlphaLocal + ainds[i]->popid * K, K);
			else
				rng.Dirichlet(buf, Alpha, K);

			double *p = Base;
			double dlnL = 0, dlnL2 = 1;

			OpenLog(dlnL, dlnL2);
			for (int64 l = 0; l < L; ++l, p += GetLoc(l).k)
			{
				GENOTYPE &gt = ainds[i]->GetGenotype(l);//fine
				if (gt.Nalleles() == 0) continue;

				ushort *als = gt.GetAlleleArray();
				for (int a = 0, vi = gt.Ploidy(); a < vi; ++a)
				{
					int al = als[a];
					ChargeLog(dlnL, dlnL2,
						SumProd(buf, p + al, KT, K) / SumProd(q, p + al, KT, K));
				}
			}
			CloseLog(dlnL, dlnL2);

			if (dlnL >= 0 || rng.Uniform() < exp(dlnL))
				SetVal(q, buf, K);
			q += K;
		}
		binaryq = false;
	}

	/* Update a priori ancetral proportion */
	TARGET void STRUCTURE::UpdateQ()
	{
		//draw gene proportion for each individual
		if (!admix)
		{
			if (m >= nadmburnin || locpriori)
				UpdateQNoAdmix();
			else
				UpdateQAdmix();

			if (locpriori)
			{
				SetZero(Di, S * K);
				double *q = Q;
				for (int i = 0; i < N; ++i)
				{
					Add(Di + ainds[i]->popid * K, q, K);
					q += K;
				}
			}
		}
		else
		{
			if (metrofreq > 0 && m % metrofreq == 0)
				UpdateQMetro(); 
			else
				UpdateQAdmix(); 
		}
	}

	/* Update locpriori parameters */
	TARGET void STRUCTURE::UpdateLocPriori()
	{
		if (!locpriori) return;
		if (admix)
		{
			double rm = rng.Uniform(r - epsr, r + epsr);
			if (rm > 0 && rm < maxr)
			{
				double dlnL = 0, d = rm - r, dt1 = rm * MyLog(rm) - r * MyLog(r);
				for (int k = 0; k < K; ++k)
				{
					dlnL += S * (Alpha[k] * dt1 - LogGamma1(rm * Alpha[k]) + LogGamma1(r * Alpha[k]));
					dlnL += Alpha[k] * d * LogProd(AlphaLocal + k, S, K) - Sum(AlphaLocal + k, S, K) * d;
				}
				if (dlnL >= 0 || rng.Uniform() < exp(dlnL))
					r = rm;
			}
		}
		else
		{
			double rm = rng.Uniform(r - epsr, r + epsr);
			if (rm > 0 && rm < maxr)
			{
				double dlnL = LogGamma1(rm) - LogGamma1(r);
				for (int k = 0; k < K; ++k)
					dlnL += LogGamma1(r * Eta[k]) - LogGamma1(rm * Eta[k]);
				dlnL *= S;

				for (int k = 0; k < K; ++k)
					dlnL += (rm - r) * Eta[k] * LogProd(Gamma + k, S, K);

				if (dlnL >= 0 || rng.Uniform() < exp(dlnL))
					r = rm;
			}

			//update eta
			if (K >= 2)
			{
				int i1 = rng.Next(K), j1 = rng.NextAvoid(K, i1);
				double delta = rng.Uniform(epseta);
				double e1 = Eta[i1] + delta;
				double e2 = Eta[j1] - delta;
				if (e1 < 1 && e2 > 0)
				{
					double dlnL = S * (LogGamma1(r * Eta[i1]) + LogGamma1(r * Eta[j1]) - LogGamma1(r * e1) - LogGamma1(r * e2));
					dlnL += r * delta * LogProdDiv(Gamma + i1, Gamma + j1, S, K);

					if (dlnL >= 0 || rng.Uniform() < exp(dlnL))
					{
						Eta[i1] += delta;
						Eta[j1] -= delta;
					}
				}
			}

			//update gamma
			if (K >= 2)
			{
				for (int s = 0; s < S; ++s)
				{
					int i1 = rng.Next(K), j1 = rng.NextAvoid(K, i1);
					double delta = rng.Uniform(epsgamma);

					if (Gamma[s * K + i1] + delta < 1 && Gamma[s * K + j1] - delta > 0)
					{
						double dlnL = (r * Eta[i1] - 1.0 + Di[s * K + i1]) * MyLog(1 + delta / Gamma[s * K + i1]) +
							(r * Eta[j1] - 1.0 + Di[s * K + j1]) * MyLog(1 - delta / Gamma[s * K + j1]);

						if (dlnL >= 0 || rng.Uniform() < exp(dlnL))
						{
							Gamma[s * K + i1] += delta;
							Gamma[s * K + j1] -= delta;
						}
					}
				}
			}

			if (m % 10 == 0)
			{
				Unify(Eta, K);
				Unify(Gamma, S, K);
			}
		}
	}

	/* Update ancestral proportion for each allele or for each individual */
	TARGET void STRUCTURE::UpdateZ()
	{
		//Update Z, Mi, Ni
		SetZero(Mi, N * K);
		SetZero(Ni, K * KT);
		if (binaryq)
		{
			//Count number of allele copies in individual i and cluster k
			for (int i = 0; i < N; ++i)
			{
				int k = Z[i][0];
				SetVal(Z[i], (ushort)k, ainds[i]->vt);
				Mi[i * K + k] = ainds[i]->vt;
			}

			//Count number of allele k2 in cluster k
			for (int64 l = 0, o = 0; l < L; ++l, o += GetLoc(l).k)
			{
				GENO_ITERATOR rt(0u, l, true);//20220316
				GENOTYPE *gtab = GetLoc(l).GetGtab();

				for (int i = 0; i < N; ++i)
				{
					GENOTYPE &gt = gtab[rt.Read()];
					if (gt.Nalleles() == 0) continue;

					ushort *als = gt.GetAlleleArray();
					int *ni = Ni + Z[i][0] * KT + o;
					for (int a = 0, vi = gt.Ploidy(); a < vi; ++a)
						ni[als[a]]++;
				}
			}
		}
		else
		{
			SetVal(ZZ, Z, N);
			for (int64 l = 0, o = 0; l < L; ++l, o += GetLoc(l).k)
			{
				GENO_ITERATOR rt(0u, l, true);//20220316
				GENOTYPE *gtab = GetLoc(l).GetGtab();

				for (int i = 0; i < N; ++i)
				{
					GENOTYPE &gt = gtab[rt.Read()];
					if (gt.Nalleles() == 0) continue;

					double *q = Q + i * K;
					ushort *&z = ZZ[i];
					ushort *als = gt.GetAlleleArray();
					double *mi = Mi + i * K;

					for (int a = 0, vi = gt.Ploidy(); a < vi; ++a)
					{
						if (a == 0 || als[a] != als[a - 1])//save time
							for (int k = 0; k < K; ++k)
								buf[k] = q[k] * cluster[k].GetFreq(l, als[a]);

						//draw cluster for each allele copy
						ushort k2 = (ushort)rng.Poly(buf, K);
						*z++ = k2;
						mi[k2]++;
						Ni[k2 * KT + o + als[a]]++;
					}
				}
			}
		}
	}

	/* Update Dirichlet parameter alpha (to draw admixture proportion Q) in the admix model */
	TARGET void STRUCTURE::UpdateAlpha()
	{
		if (locpriori && admix)
		{
			//UpdateAlphaLocPrior
			//update AlphaGlobal
			for (int k = 0; k < K; ++k)
			{
				double am = rng.Normal(Alpha[k], stdalpha);
				if (am > 0 && am < maxalpha)
				{
					double d = am - Alpha[k];
					double dlnL = S * (d * r * MyLog(r) + LogGamma1(r * Alpha[k]) - LogGamma1(r * am));

					dlnL += r * d * LogProd(AlphaLocal + k, S, K);

					if (dlnL >= 0 || rng.Uniform() < exp(dlnL))
						Alpha[k] = am;
				}
			}

			//update SumAlpha
			if (m % 10 == 0) for (int s = 0; s < S; ++s)
			{
				SumAlpha[s] = 0;
				for (int k = 0; k < K; ++k)
					SumAlpha[s] += AlphaLocal[s * K + k];
			}

			//update AlphaLocal
			for (int s = 0; s < S; ++s)
			{
				for (int k = 0; k < K; ++k)
				{
					double ao = AlphaLocal[s * K + k];
					double am = rng.Normal(ao, stdalpha);
					if (am > 0 && am < maxalpha)
					{
						double d = am - ao;
						double dlnL = (r * Alpha[k] - 1) * MyLog(am / ao) - d * r +
							apops[s]->nind * (LogGamma1(SumAlpha[s] + d) - LogGamma1(am)
								- LogGamma1(SumAlpha[s]) + LogGamma1(ao));

						IND **ti = apops[s]->inds;
						double li = 0, li2 = 1;
						OpenLog(li, li2);
						for (int id = 0; id < apops[s]->nind; ++id)
							ChargeLog(li, li2, Q[ti[id]->indid * K + k]);
						CloseLog(li, li2);
						dlnL += d * li;

						if (dlnL >= 0 || rng.Uniform() < exp(dlnL))
						{
							AlphaLocal[s * K + k] = am;
							SumAlpha[s] += d;
						}
					}
				}
			}
		}
		else if (inferalpha)
		{
			if (!admix && m >= nadmburnin) return;
			//for admix and not locpriori

			if (diffalpha)
			{
				double sumalpha = 0;
				for (int k = 0; k < K; ++k)
					sumalpha += Alpha[k];

				for (int k = 0; k < K; ++k)
				{
					double am = rng.Normal(Alpha[k], stdalpha);
					if (am <= 0 || (am >= maxalpha && uniformalpha)) continue;
					double ao = Alpha[k];
					double d = am - ao;
					double dlnL = uniformalpha ? 0 : ((alphapriora - 1) * MyLog(am / ao) + (ao - am) / alphapriorb);

					dlnL += (am - ao) * LogProd(Q + k, N, K) + (LogGamma1(sumalpha + d) - LogGamma1(sumalpha) - LogGamma1(am) + LogGamma1(ao)) * N;
					if (dlnL >= 0 || rng.Uniform() < exp(dlnL))
						Alpha[k] = am;
				}
			}
			else
			{
				double am = rng.Normal(Alpha[0], stdalpha);
				if (am <= 0 || (am >= maxalpha && uniformalpha)) return;
				double ao = Alpha[0];
				double dlnL = uniformalpha ? 0 : ((alphapriora - 1) * MyLog(am / ao) + (ao - am) / alphapriorb);

				dlnL += (am - ao) * LogProd(Q, N * K) + (LogGamma1(K * am) - LogGamma1(K * ao) - K * LogGamma1(am) + K * LogGamma1(ao)) * N;
				if (dlnL >= 0 || rng.Uniform() < exp(dlnL))
					for (int k = 0; k < K; ++k)
						Alpha[k] = am;
			}
		}
	}

	/* Update Dirichlet parameter lambda (to draw allele frequency) */
	TARGET void STRUCTURE::UpdateLambda()
	{
		if (!inferlambda) return;
		if (difflambda) for (int k = 0; k < K; ++k)
		{
			double lo = Lambda[k];
			double lm = rng.Normal(lo, stdlambda);
			if (lm <= 0 || lm >= maxlambda) continue;
			double dlnL = 0;
			double gd = LogGamma1(lo) - LogGamma1(lm);

			for (int i = 2; i <= maxK; ++i)
				dlnL += kdis[i] * (LogGamma1(i * lm) - LogGamma1(i * lo) + i * gd);

			dlnL += (lm - lo) * LogProd(cluster[k].bucket, KT);

			if (dlnL >= 0 || rng.Uniform() < exp(dlnL))
				Lambda[k] = lm;
		}
		else
		{
			int np = fmodel ? 1 : K;
			double lo = Lambda[0];
			double lm = rng.Normal(lo, stdlambda);
			if (lm <= 0 || lm >= maxlambda) return;
			double dlnL = 0;
			double gd = LogGamma1(lo) - LogGamma1(lm);

			for (int i = 2; i <= maxK; ++i)
				dlnL += kdis[i] * np * (LogGamma1(i * lm) - LogGamma1(i * lo) + i * gd);

			dlnL += (lm - lo) * (fmodel ? LogProd(PA.bucket, KT) : LogProd(Base, K * KT));

			if (dlnL >= 0 || rng.Uniform() < exp(dlnL))
				SetVal(Lambda, lm, K);
		}
	}

	/* Update allele frequency of ancestral population for the F model */
	TARGET void STRUCTURE::UpdatePA()
	{
		if (!fmodel) return;
		double lam = Lambda[0];//single lambda in F model?

		if (rng.Uniform() < 0.5)
		{
			//Independent update PA
			SetVal(PA1.bucket, lam, KT);

			for (int k = 0; k < K; ++k)
				AddProd(PA1.bucket, cluster[k].bucket, f[k], KT);

			double *pa = PA.bucket, *t1 = PA1.bucket, *t2 = PA2.bucket;

			//unable to optimize
			for (int k = 0; k < K; ++k)
				buf2[k] = cluster[k].bucket;

			for (int64 l = 0; l < L; ++l)
			{
				int k2 = GetLoc(l).k;
				rng.Dirichlet(t2, t1, k2);

				double dlnL = 0;
				for (int a = 0; a < k2; ++a)
					dlnL += (t1[a] - lam) * (MyLog(pa[a] / t2[a]));

				for (int k = 0; k < K; ++k)
				{
					double fr = f[k];
					for (int a = 0; a < k2; ++a)
						dlnL += LogGamma1(fr * pa[a]) - LogGamma1(fr * t2[a])
						+ fr * (t2[a] - pa[a]) * MyLog(*buf2[k]++);
				}

				if (dlnL >= 0 || rng.Uniform() < exp(dlnL))
					SetVal(pa, t2, k2);

				pa += k2; t1 += k2; t2 += k2;
			}
		}
		else
		{
			double *pa = PA.bucket, *t1 = PA1.bucket, *t2 = PA2.bucket;
			double *p = Base;

			for (int64 l = 0; l < L; ++l)
			{
				int k2 = GetLoc(l).k;
				int m1 = rng.Next(k2), n1 = rng.NextAvoid(k2, m1);

				double delta = rng.Uniform(pow(N, -0.5));
				double pm0 = pa[m1], pm1 = pm0 + delta, pn0 = pa[n1], pn1 = pn0 - delta;
				if (pm1 >= 1 || pn1 <= 0)
				{
					pa += k2; t1 += k2; t2 += k2;
					continue;
				}
				double dlnL = (lam - 1) * MyLog(pm1 * pn1 / pm0 / pn0);

				for (int k = 0; k < K; ++k)
				{
					double fr = f[k];
					dlnL += LogGamma1(fr * pm0) + LogGamma1(fr * pn0)
						- LogGamma1(fr * pm1) - LogGamma1(fr * pn1) +
						(fr * delta) * MyLog(p[k * KT + m1] / p[k * KT + n1]);
					p += k2;
				}

				if (dlnL >= 0 || rng.Uniform() < exp(dlnL))
				{
					pa[m1] = pm1;
					pa[n1] = pn1;
				}

				pa += k2; t1 += k2; t2 += k2;
			}
		}
	}

	/* Update population-specific Fst for the F model */
	TARGET void STRUCTURE::UpdateF()
	{
		//F model
		if (!fmodel) return;

		// Update F/Fst
		for (int k = 0; k < K; ++k)
		{
			double Fo = F[k];
			double Fm = rng.Normal(Fo, stdf);
			if (Fm < 0 || Fm > 1)
			{
				if (fsame) break;
				else continue;
			}
			double fo = f[k], fm = (1 - Fm) / Fm;
			double dlnL = pmeanf * (Fo - Fm) / (pstdf * pstdf) +
				(pmeanf * pmeanf / (pstdf * pstdf) - 1) * MyLog(Fm / Fo) +
				(fsame ? K : 1) * nloc * (LogGamma1(fm) - LogGamma1(fo));//eL add in 2018.09

			for (int k2 = k; k2 < K; ++k2)
			{
				double *p = cluster[k2].bucket;
				double *pa = PA.bucket;//2018.09
				for (int a = 0; a < KT; ++a)
				{
					double pa2 = *pa++;
					dlnL += pa2 * (fm - fo) * MyLog(*p++) + LogGamma1(fo * pa2) - LogGamma1(fm * pa2);
				}
				if (!fsame) break;
			}

			if (dlnL >= 0 || rng.Uniform() < exp(dlnL))
			{
				if (fsame)
				{
					SetVal(F, Fm, K);
					SetVal(f, fm, K);
				}
				else
				{
					F[k] = Fm;
					f[k] = fm;
				}
			}

			if (fsame) break;
		}
	}

	/* Finalize records */
	TARGET void STRUCTURE::Arrange()
	{
		Unify(Z_cumu, N, K);
		Mul(rout, 1.0 / nr, rlen);
		rout[2] = rout[1] - rout[0] * rout[0];
		rout[3] = rout[0] - rout[2] / 2;
		Mul(Base + K * KT, 1.0 / nr, K * KT);
	}

	/* Record updated MCMC parameters */
	TARGET void STRUCTURE::Record()
	{
		if (m >= nburnin && (m - nburnin) % nthinning == 0)
		{
			nr++;
			Add(Z_cumu, Mi, N * K);
			double li = 0, li2 = 1;

			double *p = Base;
			li = 0; li2 = 1;

			OpenLog(li, li2);
			for (int64 l = 0; l < L; ++l)
			{
				GENO_ITERATOR rt(0u, l, true);
				GENOTYPE *gtab = GetLoc(l).GetGtab();
				double *q = Q;

				for (int i = 0; i < N; ++i, q += K)
				{
					GENOTYPE &gt = gtab[rt.Read()];
					if (gt.Nalleles() == 0) continue;
					
					ushort *als = gt.GetAlleleArray();
					for (int a = 0, vi = gt.Ploidy(); a < vi; ++a)
						ChargeLog(li, li2, SumProd(q, p + als[a], KT, K));
				}
				p += GetLoc(l).k;
			}
			CloseLog(li, li2);

			//add to result
			rout[0] += li;
			rout[1] += li * li;

			int nf = fmodel ? (fsame ? 1 : K) : 0;
			int na = admix ? (diffalpha ? K : 1) : 0;
			int nl = difflambda ? K : 1;

			Add(rout + 4, Lambda, nl);

			if (locpriori)
			{
				rout[4 + nl] += r;
				Add(rout + 4 + nl + 1, admix ? Alpha : Eta, K);
				Add(rout + 4 + nl + 1 + K, admix ? AlphaLocal : Gamma, S * K);
			}
			else if (admix)
				Add(rout + 4 + nl, Alpha, na);
			Add(rout + rlen - nf, F, nf);
			Add(Base + K * KT, Base, K * KT);
		}
		PROGRESS_VALUE++;
	}

	/* Free memory */
	TARGET void STRUCTURE::Uninit()
	{
		//Normal
		TryDelete(Z_cumu);
		TryDelete(rout);
		TryDelete(cluster);
		TryDelete(clusterb);
		TryDelete(Base);
		TryDelete(Lambda);
		TryDelete(buf);
		TryDelete(bufb);
		TryDelete(buf2);

		//Adm
		if (Z) for (int i = 0; i < N; ++i)
			TryDelete(Z[i]);
		TryDelete(Z);
		TryDelete(ZZ);
		TryDelete(Q);
		TryDelete(Mi);
		TryDelete(Ni);

		//Fmodel
		TryDelete(F);
		TryDelete(f);

		//Loc
		TryDelete(Alpha);
		TryDelete(SumAlpha);
		TryDelete(AlphaLocal);
		TryDelete(Di);
		TryDelete(Eta);
		TryDelete(Gamma);
		TryDelete(kdis);
	}

	/* Perform MCMC */
	TARGET void STRUCTURE::MCMC()
	{
		InitFix();
		InitAdmix();
		InitLocPriori();
		InitFmodel();

		nadmburnin = (!locpriori && !admix) ? nadmburnin : 0;
		nburnin += nadmburnin;
		int ntrep = nburnin + nreps;

		for (m = 0; m < ntrep; ++m)
		{
			binaryq = false;
			UpdateP();//Allele frequency
			UpdateQ();//Individual priori gene proportion
			UpdateLocPriori();//LocPriori r
			UpdateZ();//Individual gene origin
			UpdateAlpha();//LocPriori+Admix
			UpdateLambda();//
			UpdatePA();//
			UpdateF();//
			Record();
		}
		Arrange();
	}
#endif

#ifndef _EXHAPLO
	/* Do nothing */
	TARGET HAPLO_DUMMY_HAPLOTYPE::HAPLO_DUMMY_HAPLOTYPE()
	{

	}

	/* Extract the ith haplotype from an individual */
	TARGET void HAPLO_DUMMY_HAPLOTYPE::ExtractHaplotype(byte vi, IND *ti, int64 st, int64 ed, int nvar, ushort aid, MEMORY &haplo_memory)
	{
		alleleid = aid;
		haplo_memory.Alloc(alleles, nvar);
		int sc = 0;

		for (int64 l = st; l <= ed; ++l)
		{
			int ta = ti->GetGenotype(GetLocId(l), GetLoc(l).GetGtab()).GetAlleleCopy(vi);//fine

			//encounter an missing allele, set all haplotype to missing
			if (ta == 0xFFFF)
			{
				SetFF(alleles, nvar);
				return;
			}
			alleles[sc++] = ta;
		}
	}

	/* Print information for an extracted locus */
	TARGET void HAPLO_DUMMY_HAPLOTYPE::PrintHaplotype(FILE *f1, int64 st, int64 ed)
	{
		fprintf(f1, "%s%d", g_linebreak_val, alleleid);
		for (int64 l = st, sc = 0; l <= ed; ++l)
			fprintf(f1, "%c%s", g_delimiter_val, GetLoc(l).GetAlleleName(alleles[sc++]));
	}
#endif

#ifndef _VCF

	/* Initialize */
	TARGET void Initialize()
	{
		// Initialize variables
		individual_memory = NULL;
		locus_memory = NULL;
		conversion_memory = NULL;
		conversion_memory2 = NULL;
		conversion_string = NULL;
		gd_memory = NULL;

		gdtab = NULL;
		gdlock = NULL;

		genotype_bucket = NULL;
		genotype_size = 0;
		genotype_coffset = 0;
		genotype_offset = NULL;

		alleledepth_bucket = NULL;
		alleledepth_size = 0;
		alleledepth_coffset = 0;
		alleledepth_offset = NULL;

		allele_freq_offset = NULL;
		maxK = 0;
		KT = 0;

		genotype_count_offset = NULL;
		maxG = 0;
		GT = 0;

		load_buf = NULL;
		load_buf_size = 0;
		vcf_header = NULL;

		pop.Clear();
		for (int i = 0; i < reg.size; ++i)
			reg[i].Clear();
		reg.Clear();

		locus = NULL;
		slocus = NULL;
		useslocus = false;

		ainds = NULL;
		apops = NULL;
		aregs = NULL;
		npop = 0;
		lreg = 0;
		SetZero(nreg, N_MAX_REG);
		nregt = 0;
		nregt2 = 0;

		total_pop = NULL;

		nloc = 0;
		nfilter = 0;
		nind = 0;
		nfilterind = 0;
		reassigned = 0;
		rinds = NULL;
		cpop = NULL;
		maxploidy = 0;
		minploidy = 0;

		progress1 = 0;
		progress2 = 0;
		state_lock = NULL;

		haplotype_offset = NULL;
		haplotype_bucket = NULL;
		haplotype_size = 0;
		qslstack.Clear();
		haplotype_locus.Clear();

		convert_file = NULL;
		convert_buf = NULL;
		convert_linesize = 0;

		diversity_buf = NULL;
		SetZero(&diversity_sum, 1);
		diversity_stage = 0;

		SetZero(fst_buf, 6);
		fst_type = 0;

		gdist_buf = NULL;
		gdist_type = 0;

		amova_buf = NULL;

		relatedness_buf = NULL;
		kinship_buf = NULL;
		SetZero(gdindex, N_GD_ESTIMATOR + 1);

		pcoa_matrix = NULL;
		structure_par = NULL;
		structure_totalruns = 0;

		spa_x = NULL;
		spa_n = 0;
		spa_np = 0;

		if (g_input_row * g_input_col == 0)
			Exit("\nError: no input files detected.\n");

		g_filetotlen = 0;
		for (int i = 0; i < g_input_row; ++i)
		{
			for (int j = 0; j < g_input_col; ++j)
			{
				if ((g_filehandle[i][j] = fopen(g_filepath[i][j], "rb")) == NULL)
					Exit("\nError: input file %s does not exist.\n", g_filepath[i][j]);
				uint magic = FGetUshort(g_filehandle[i][j]);
				fclose(g_filehandle[i][j]);

				if (magic == 0x8B1F)//*.gz
				{
					g_format_val = -abs(g_format_val);
					if ((g_filehandle[i][j] = fopen(g_filepath[i][j], "rb")) == NULL)
						Exit("\nError: input file %s does not exist.\n", g_filepath[i][j]);
					fseeko64(g_filehandle[i][j], 0, SEEK_END);
					g_filelen[i][j] = ftello64(g_filehandle[i][j]);
					fclose(g_filehandle[i][j]);

					g_filetotlen += g_filelen[i][j];
					g_filehandle[i][j] = FOpen(g_filepath[i][j], "rb");
				}
				else
				{
					if ((g_filehandle[i][j] = FOpen(g_filepath[i][j], "rb")) == NULL)
						Exit("\nError: input file %s does not exist.\n", g_filepath[i][j]);
					g_filelen[i][j] = GetFileLen(g_filehandle[i][j]);
					g_filetotlen += g_filelen[i][j];
				}

				if (g_filehandle[i][j] == NULL)
					Exit("\nError: cannot open input file %s.\n", g_filepath[i][j]);
			}
		}

		if (g_format_val > 2)
		{
			g_filebuf = new char[64 * 1024];
			setvbuf(g_filehandle[0][0], g_filebuf, _IOFBF, 64 * 1024);
		}
		else
			g_filebuf = NULL;

		genotype_coffset = 0;
		V_BASE_GENOTYPE    = 0x080000000000;
		V_BASE_ALLELEDEPTH = 0x100000000000;
		ALLELE_IDENTIFIER = abs(g_format_val) <= 2;

		FRES_BUF = new char[LINE_BUFFER];
		FRES_NAME = new char[FILE_NAME_LEN];
		omp_set_num_threads(g_nthread_val);
		Eigen::initParallel();

		InitBinomial();
		InitAlpha();
		InitCryptTable();
		InitCryptTable();

		if (g_input_row * g_input_col != 1 && abs(g_format_val) > 2)
			Exit("\nError: multiple input files are only supported for vcf format.\n");

		individual_memory = new MEMORY[1];
		locus_memory = new MEMORY[g_nthread_val];

		//initial missing genotypes
		SetFF(missing_array, N_MAX_PLOIDY);
		ushort *nullgatab = NULL;
		for (int i = 0; i <= N_MAX_PLOIDY; ++i)
		{
			missing_hash[i] = HashGenotype(missing_array, i);
			new(&missing_genotype[i]) GENOTYPE(nullgatab, missing_array, i);//sorted
		}

		POP tp;
		LIST<POP> tlist(NULL);
		pop.Push(tp);
		new(&pop[0]) POP((char*)"DefPop", NULL, 0, 0, 0, 0, true);

		//assign pop 1
		if (g_indfile_b || g_indtext_b)
		{
			int indtextlen = (int)strlen(g_indtext_val);
			char *indtextend = g_indtext_val + indtextlen;
			ReplaceChar(g_indtext_val, '\r', '\n');
			ReplaceChar(g_indtext_val, ' ', '\n');
			char *t1 = g_indtext_val;
			char *t2 = g_indtext_val - 1;
			int rl = -2;
			while (*t1 && t1 < indtextend)
			{
				t1 = t2 + 1;
				while (*t1 == '\n' || *t1 == '\r' || *t1 == ' ' || *t1 == '\t') t1++;
				if (t1 >= indtextend) break;

				if (reg.size == 0 || LwrLineCmp("#REG", t1) == 0)
				{
					//add reg
					rl++;
					if (reg.size)
					{
						t1 = StrNextIdx(t1, "\n", 1) + 1;
						while (*t1 == '\n' || *t1 == '\r' || *t1 == ' ' || *t1 == '\t') t1++;
					}
					POP tr2;
					reg.Push(tlist);
					reg[rl + 1].Push(tr2);
					char *b = new char[15];
					sprintf(b, "DefRegL%d", rl + 2);
					new(&reg[rl + 1][0]) POP(b, NULL, 0, 0, 0, 0, false);
				}

				if (rl == -1)
				{
					t2 = StrNextIdx(t1, ":", 1);
					if (t2 == NULL) break;
					*t2 = '\0';

					tp.name = t1;

					for (int i = 0; i < pop.size; ++i)
						if (strcmp(pop[i].name, t1) == 0)
							Exit("\nError: Two populations have the same name: %s\n", t1);

					t1 = t2 + 1;
					t2 = StrNextIdx(t1, "\n", 1);
					if (t2) *t2 = '\0';
					else t2 = StrNextIdx0(t1, 1);

					tp.ispop = true;
					tp.id = pop.size;
					tp.names = SplitStr(t1, ',', tp.nind);//deleted
					tp.inds = NULL;
					tp.allelefreq = NULL;
					tp.genocount = NULL;
					tp.loc_stat = NULL;
					tp.rid = 0;
					tp.vpop = NULL;
					tp.npop = 0;
					pop.Push(tp);
				}
				else
				{
					t2 = StrNextIdx(t1, ":", 1);
					if (t2 == NULL) break;
					*t2 = '\0';

					POP tr3;
					tr3.ispop = false;
					tr3.id = reg[rl].size;
					tr3.name = t1;

					t1 = t2 + 1;
					t2 = StrNextIdx(t1, "\n", 1);
					if (t2) *t2 = '\0';
					else t2 = StrNextIdx0(t1, 1);

					tr3.names = SplitStr(t1, ',', tr3.npop);//deleted
					tr3.nind = 0;
					tr3.inds = NULL;
					tr3.allelefreq = NULL;
					tr3.genocount = NULL;
					tr3.loc_stat = NULL;
					tr3.vpop = NULL;
					tr3.rid = 0;

					//check subreg/pop name
					for (int rl2 = 0; rl2 <= rl; ++rl2)
						for (int i = 0; i < reg[rl2].size; ++i)
							if (strcmp(reg[rl2][i].name, tr3.name) == 0)
								Exit("\nError: repeat population/region name: %s.\n", tr3.name);

					for (int i = 0; i < pop.size; ++i)
						if (strcmp(pop[i].name, tr3.name) == 0)
							Exit("\nError: repeat population/region name: %s.\n", tr3.name);

					//assign pop
					LIST<POP> &popv = rl == 0 ? pop : reg[rl - 1];
					int npop2 = 0, namelen = 0;

					for (int j2 = 0; j2 < tr3.npop; ++j2)
					{
						int v1 = 0, v2 = 0;
						ParseTwoNumber(tr3.names[j2], v1, v2);

						// #1-#2 format
						if (v1 != -1 && v2 != -1)
						{
							if (v1 <= 0 || v1 >= popv.size)
								Exit("\nError: the first number of indtext %s for region %s exceeds range.", tr3.names[j2], tr3.name);
							if (v2 <= 0 || v2 >= popv.size)
								Exit("\nError: the second number of indtext %s for region %s exceeds range.", tr3.names[j2], tr3.name);

							int vmin = Min(v1, v2), vmax = Max(v1, v2);
							for (int i = vmin; i <= vmax; ++i)
							{
								if (popv[i].rid == tr3.id)
									Exit("\nError: Population/subregion appear twice in Region %s.\n", popv[i].name, tr3.name);
								if (popv[i].rid != 0)
									Exit("\nError: Two regions %s and %s have the same population/subregion: %s.\n", reg[rl][popv[i].rid].name, tr3.name, popv[i].name);
								namelen += (int)strlen(popv[i].name) + 1;
								popv[i].rid = (ushort)tr3.id;
								npop2++;
							}
						}
						// #1 or id format
						else 
						{
							bool found = false;
							for (int i = 0; i < popv.size; ++i)
								if (v1 == i || strcmp(popv[i].name, tr3.names[j2]) == 0)
								{
									if (popv[i].rid == tr3.id)
										Exit("\nError: Population/subregion appear twice in Region %s.\n", popv[i].name, tr3.name);
									if (popv[i].rid != 0)
										Exit("\nError: Two regions %s and %s have the same population/subregion: %s.\n", reg[rl][popv[i].rid].name, tr3.name, popv[i].name);
									namelen += (int)strlen(popv[i].name) + 1;
									popv[i].rid = (ushort)tr3.id;
									npop2++;
									found = true;
									break;
								}

							if (!found)
								Exit("\nError: Can not found population/subregion %s in region %s.\n", tr3.names[j2], tr3.name);
						}
					}

					delete[] tr3.names;
					tr3.names = NULL;
					tr3.npop = npop2;
					reg[rl].Push(tr3);
				}
			}
		}

		NBUF = CALC_THREAD_BUFFER * g_nthread_val;
		state_lock = new atomic<int64>[NBUF];
	}

	/* UnInitialize */
	TARGET void UnInitialize()
	{
		//Exit("\nCalculation complete!\n");

		if (lreg == -1 && npop == 1)
			delete[] total_pop->name;

		delete[] FRES_BUF;
		delete[] FRES_NAME;

		if (useslocus)
		{
			delete[] slocus;
			slocus = NULL;
		}
		else
		{
			delete[] locus;
			locus = NULL;
		}

		if (ainds)
		{
			delete[] ainds;
			ainds = rinds = NULL;
		}

		for (int rl = 0; rl < reg.size; ++rl)
			for (int i = 0; i < reg[rl].size; ++i)
			{
				if (i == 0) delete[] reg[rl][i].name;
				reg[rl][i].Uninitialize();
			}

		if (aregs)
		{
			for (int rl = 0; rl < reg.size; ++rl)
				if (aregs[rl] && (uint)aregs[rl] % 0xFFFFFFFF != 0xCDCDCDCD)
					delete[] aregs[rl];
			delete[] aregs;
		}

		for (int i = 0; i < pop.size; ++i)
			pop[i].Uninitialize();

		if (apops) delete[] apops;

		qslstack.~LIST();
		pop.~LIST();
		for (int rl = 0; rl < reg.size; ++rl)
			reg[rl].~LIST();
		reg.~LIST();
		haplotype_locus.~LIST();

		if (VIRTUAL_MEMORY)
		{
			VUnAllocGenotype();
			VUnAllocAlleleDepth();
		}
		else
		{
			delete[] genotype_bucket;
			if (ploidyinfer) delete[] alleledepth_bucket;
		}
		VIRTUAL_MEMORY = false;

		genotype_size = 0;
		if (genotype_offset) delete[] genotype_offset;

		if (ploidyinfer)
		{
			alleledepth_size = 0;
			delete[] alleledepth_offset;
		}

		if (allele_freq_offset) delete[] allele_freq_offset;
		if (genotype_count_offset)  delete[] genotype_count_offset;

		delete[] cryptTable;
		delete[] individual_memory;
		delete[] locus_memory;
		delete[] state_lock;
	}

	/* Close input files */
	TARGET void CloseInput()
	{
		for (int i = 0; i < g_input_row; ++i)
		{
			for (int j = 0; j < g_input_col; ++j)
			{
				FClose(g_filehandle[i][j]);
				if (abs(g_format_val) == 2)
					delete g_bcfheader[i][j];
			}
			delete[] g_filehandle[i];
			delete[] g_filelen[i];
			if (abs(g_format_val) == 2) delete[] g_bcfheader[i];
		}
		if (abs(g_format_val) == 2) delete[] g_bcfheader;
		if (g_filenamebuf) delete[] g_filenamebuf;
		if (g_filehandle) delete[] g_filehandle;
		if (g_filebuf) delete[] g_filebuf;
		if (g_filelen) delete[] g_filelen;
	}

	/* Assign individual indid to population popid */
	TARGET void AssignIndSub(int indid, int popid, int &namelen, int &nind2)
	{
		if (ainds[indid]->popid != 0)
			Exit("\nError: populations %s and %s have the same individual %s.", pop[ainds[indid]->popid].name, pop[popid].name, ainds[indid]->name);

		ainds[indid]->popid = (ushort)popid;
		namelen += (int)strlen(ainds[indid]->name) + 1;
		nind2++;
	}

	/* Assign individuals to their populations */
	TARGET void AssignInds()
	{
		//Assign all individuals to defpop initially and prepare name to indid hash table
		TABLE<HASH, int> name2id(false, NULL);
		for (int i = 0; i < nind; ++i)
		{
			ainds[i]->popid = 0;
			HASH ha = HashString(ainds[i]->name);
			if (name2id.ContainsKey(ha))
				Exit("\nError: individiuals #%d and #%d have the same name %s.", name2id[ha], i, ainds[i]->name);
			name2id[ha] = i;
		}

		//Parse indtext
		for (int pp = pop.size - 1; pp >= 0; --pp)
		{
			int nind2 = 0;
			int namelen = 0;
			for (int j2 = 0; j2 < pop[pp].nind; ++j2)
			{
				int v1 = 0, v2 = 0;
				ParseTwoNumber(pop[pp].names[j2], v1, v2);

				//id format
				if (v1 == -1 && v2 == -1)
				{
					HASH ha = HashString(pop[pp].names[j2]);
					if (!name2id.ContainsKey(ha))
						Exit("\nError: individiual %s in population %s is not in the genotype data.", pop[pp].names[j2], pop[pp].name);

					AssignIndSub(name2id[ha], pp, namelen, nind2);
				}

				//#1 or id format
				else if (v2 == -1)
				{
					HASH ha = HashString(pop[pp].names[j2]);
					int indid = -1;
					if (name2id.ContainsKey(ha))
						indid = name2id[ha];
					else
						indid = v1 - 1;

					if (indid >= 0 && indid < nind)
						AssignIndSub(indid, pp, namelen, nind2);
					else
						Exit("\nError: can not found individual %s for population %s.", pop[pp].names[j2], pop[pp].name);
				}

				//#1-#2 format
				else
				{
					if (v1 <= 0 || v1 > nind)
						Exit("\nError: the first number of indtext %s for population %s exceeds range.", pop[pp].names[j2], pop[pp].name);
					if (v2 <= 0 || v2 > nind)
						Exit("\nError: the second number of indtext %s for population %s exceeds range.", pop[pp].names[j2], pop[pp].name);

					int vmin = Min(v1, v2), vmax = Max(v1, v2);
					for (int i = vmin - 1; i <= vmax - 1; ++i)
						AssignIndSub(i, pp, namelen, nind2);
				}
			}

			if (pop[pp].nind == 0) for (int i = 0; i < nind; ++i)
				if (ainds[i]->popid == pp)
				{
					ainds[i]->popid = (ushort)pp;
					namelen += (int)strlen(ainds[i]->name) + 1;
					nind2++;
				}

			//Allocate individuals names
			delete[] pop[pp].names;
			pop[pp].names = NULL;
			pop[pp].nind = 0;
		}

		//pop[i].inds will be allocated and copy individuals after filtered in function AssignVInds
	}

	TARGET void LoadFile()
	{
		//Vcf/bcf format
		if (abs(g_format_val) <= 2)
		{
			int64 buflen = 1024 * 64;
			vcf_header = new char[buflen]; vcf_header[0] = '\0';
			int64 title_buflen = buflen, title_len = 0;
			char ***header = new char**[g_input_row];
			char *readbuf = new char[buflen];
			if (abs(g_format_val) == 2)
				g_bcfheader = new BCFHEADER**[g_input_row];

			for (int i = 0; i < g_input_row; ++i)
			{
				header[i] = new char*[g_input_col];
				if (abs(g_format_val) == 2)
					g_bcfheader[i] = new BCFHEADER*[g_input_col];

				for (int j = 0; j < g_input_col; ++j)
				{
					if (abs(g_format_val) == 2)
					{
						FRead(readbuf, 1, 9, g_filehandle[i][j]);
						if (*(uint*)readbuf != 0x02464342) Exit("\nError: Unsupported BCF format.\n");
						int headerlen = *(int*)(readbuf + 5);
						VLA_NEW(headerbuf, char, headerlen);
						FRead(headerbuf, 1, headerlen, g_filehandle[i][j]);
						g_bcfheader[i][j] = new BCFHEADER(headerbuf);
						VLA_DELETE(headerbuf);
						FSeek(g_filehandle[i][j], 9, SEEK_SET);
					}

					bool findhead = false;
					int64 readlen = 0;
					while (!FEof(g_filehandle[i][j]))
					{
						readlen = FGets(readbuf, buflen, g_filehandle[i][j]);
						if (LwrLineCmp("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	", readbuf) == 0)
						{
							if (readbuf[readlen - 1] == '\n') readbuf[--readlen] = '\0';
							if (readbuf[readlen - 1] == '\r') readbuf[--readlen] = '\0';
							header[i][j] = new char[readlen + 1];
							strcpy(header[i][j], readbuf);
							findhead = true;
							if (abs(g_format_val) == 2)
								g_bcfheader[i][j]->nsample = CountChar(header[i][j], '\t') - 8;
							break;
						}
					}

					if (i == 0)
					{
						if (title_len + readlen + 1 >= title_buflen)
						{
							title_buflen <<= 1;
							char *tbuf2 = new char[title_buflen];
							strcpy(tbuf2, vcf_header);
							delete[] vcf_header;
						}
						char *tbuf = StrNextIdx(header[i][j], '\t', 9);
						int64 tlen = readlen - (tbuf - header[i][j]);
						memcpy(vcf_header + title_len, tbuf, tlen); //concatenate individual names
						title_len += tlen;
						vcf_header[title_len] = '\0';
					}

					if ((abs(g_format_val) == 2) && g_bcfheader[i][j]->format_gqid != g_bcfheader[i][0]->format_gqid)
						Exit("\nError: BCF format error, FORMAT GQ tag indices are different among BCFs. \n", g_filename[i][j], g_filename[i][0]);

					if ((abs(g_format_val) == 2) && g_bcfheader[i][j]->format_dpid != g_bcfheader[i][0]->format_dpid)
						Exit("\nError: BCF format error, FORMAT DP tag indices are different among BCFs. \n", g_filename[i][j], g_filename[i][0]);

					if ((abs(g_format_val) == 2) && g_bcfheader[i][j]->format_adid != g_bcfheader[i][0]->format_adid)
						Exit("\nError: BCF format error, FORMAT AD tag indices are different among BCFs. \n", g_filename[i][j], g_filename[i][0]);

					if (i && strcmp(header[i][j], header[0][j]))
						Exit("\nError: VCF/BCF format error, individual identifiers in %s and %s are different. \n", g_filename[i][j], g_filename[i][0]);

					if (!findhead)
						Exit("\nError: VCF/BCF format error, cannot find genotype table header in %s.\n", g_filename[i][j]);

					if (abs(g_format_val) == 2)
						FSeek(g_filehandle[i][j], 1, SEEK_CUR);
				}
			}

			for (int i = 0; i < g_input_row; ++i)
			{
				for (int j = 0; j < g_input_col; ++j)
					delete[] header[i][j];
				delete[] header[i];
			}
			delete[] header;
			delete[] readbuf;

			//Create locus
			RunThreads((abs(g_format_val) == 1) ? &GetVCFLines : &GetBCFLines, NULL, NULL, g_filetotlen, g_filetotlen, "\nLoading VCF/BCF:\n", 1, true, (int)(0.15 * g_progress_val));
			locus = new LOCUS[nloc];
			SetZero(locus, nloc);

			//Create individuals
			vcf_header[title_len + 1] = '\0';
			nind = CountChar(vcf_header, '\t');
			char *title = vcf_header + 1;
			ainds = new IND*[nind];

			genotype_offset = new OFFSET[nloc];
			genotype_size = (((int64)((ceil(nind / 4.0) + 3) * nloc * 1.10) >> 16) << 16) + 65536;
			genotype_bucket = (byte*)V_BASE_GENOTYPE;
			VAllocGenotype(genotype_size);

			if (ploidyinfer)
			{
				alleledepth_offset = new OFFSET[nloc];
				alleledepth_size = (((uint64)((nind + 3) * nloc * 1.10) >> 16) << 16) + 65536;
				alleledepth_bucket = (byte*)V_BASE_ALLELEDEPTH;
				VAllocAlleleDepth(alleledepth_size);
			}

			int indc = 0;
			do
			{
				individual_memory->Alloc(ainds[indc], 1);
				new(ainds[indc]) IND(title, indc);
				indc++;
			} while (*title);
			delete[] vcf_header; vcf_header = (char*)"\0";

			AssignInds();

			progress1 = 0;
			load_buf = new char*[NBUF];
			load_buf_size = new int [NBUF];
			for (int i = 0; i < NBUF; ++i)
			{
				load_buf_size[i] = 64 * 1024;
				load_buf[i] = new char[load_buf_size[i]];
			}

			PROGRESS_VALUE = 0; PROGRESS_TOTAL = nloc; PROGRESS_NOUTPUTED2 = PROGRESS_NOUTPUTED; PROGRESS_NOUTPUTED = 0;

			if (abs(g_format_val) == 1)
				RunThreads(&LoadVCF, &LoadVCFGuard, NULL, nloc, nloc, NULL, g_nthread_val, false, (int)(0.85 * g_progress_val));
			else
				RunThreads(&LoadBCF, &LoadBCFGuard, NULL, nloc, nloc, NULL, g_nthread_val, false, (int)(0.85 * g_progress_val));

			//Use allele depth as allele frequency, each ind represents a population
			if (ad)
			{
				ad++;
				genotype_coffset = 0;
				progress1 = progress2 = 0;
				allele_freq_offset = new LOCN[nloc];
				genotype_count_offset = new LOCN[nloc];
				SetFF(allele_freq_offset, nloc);
				SetFF(genotype_count_offset, nloc);
				KT = GT = maxK = maxG = 0;

				for (int64 l = 0; l < nloc; ++l)
				{
					allele_freq_offset[l] = KT;
					genotype_count_offset[l] = GT;
					KT += locus[l].k;
					GT += locus[l].ngeno;
					maxK = Max((int)locus[l].k, maxK);
					maxG = Max((int)locus[l].ngeno, maxG);
				}

				// use allele frequency only 
				for (int i = 0; i < pop.size; ++i)
				{
					pop[i].AllocFreq();
					delete[] pop[i].loc_stat;
					delete[] pop[i].genocount;
					pop[i].genocount = NULL;
					pop[i].loc_stat = NULL;
				}

				if (abs(g_format_val) == 1)
					RunThreads(&LoadVCF, &LoadVCFGuard, NULL, nloc, nloc, "\nConverting allelic depth into allele frequency:\n", g_nthread_val, true, g_progress_val);
				else
					RunThreads(&LoadBCF, &LoadBCFGuard, NULL, nloc, nloc, "\nConverting allelic depth into allele frequency:\n", g_nthread_val, true, g_progress_val);

				delete[] allele_freq_offset;
				delete[] genotype_count_offset;
				allele_freq_offset = genotype_count_offset = NULL;
			}

			for (int i = 0; i < NBUF; ++i)
				delete[] load_buf[i];
			delete[] load_buf;
			delete[] load_buf_size;
		}
		//Non vcf format: genepop|spagedi|cervus|arlequin|structure|polygene|polyrelatedness
		else
		{
			if (ad) Exit("\nError: non-vcf and non-bcf format are incompatible with allelic depth (-ad) option.\n");
			if (haplotype) Exit("\nError: Cannot extract haplotype for non-vcf and non-bcf format.\n");
			switch (abs(g_format_val))
			{
			case 3:
				RunThreads(&LoadGenepop, NULL, NULL, g_filetotlen << 1, g_filetotlen << 1, "\nLoading genepop file:\n", 1, true);
				break;
			case 4:
				RunThreads(&LoadSpagedi, NULL, NULL, g_filetotlen << 1, g_filetotlen << 1, "\nLoading spagedi file:\n", 1, true);
				break;
			case 5:
				RunThreads(&LoadCervus, NULL, NULL, g_filetotlen << 1, g_filetotlen << 1, "\nLoading cervus file:\n", 1, true);
				break;
			case 6:
				RunThreads(&LoadArlequin, NULL, NULL, g_filetotlen << 1, g_filetotlen << 1, "\nLoading structure file:\n", 1, true);
				break;
			case 7:
				RunThreads(&LoadStructure, NULL, NULL, g_filetotlen << 1, g_filetotlen << 1, "\nLoading structure file:\n", 1, true);
				break;
			case 8:
				RunThreads(&LoadPolyGene, NULL, NULL, g_filetotlen << 1, g_filetotlen << 1, "\nLoading polygene file:\n", 1, true);
				break;
			case 9:
				RunThreads(&LoadPolyRelatedness, NULL, NULL, g_filetotlen << 1, g_filetotlen << 1, "\nLoading polyrelatedness file:\n", 1, true);
				break;
			}
			AssignInds();
		}

		//Finish read
		CloseInput();
		CheckGenotypeId();
		
		//Compress locus into slocus
		MEMORY *nlocus_memory = new MEMORY[g_nthread_val];
		slocus = new SLOCUS[nloc];
		if (haplotype) locus_pos = new uint64[nloc];

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			threadid = omp_get_thread_num();
			if (haplotype) locus_pos[l] = locus[l].pos;
			new(&slocus[l]) SLOCUS(nlocus_memory[threadid], locus[l]);
		}

		delete[] locus;
		delete[] locus_memory;
		locus = NULL;
		locus_memory = nlocus_memory;
		useslocus = true;

		CheckGenotypeId();
	}

	/* Check anisoploid within contigs in haplotype extraction */
	THREAD(CheckAnisoploid)
	{
		AssignPloidyThreadIn();
		nfilterind = 0;

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int i = 0; i < nind; ++i)
		{
			if (ainds[i]->vmin == ainds[i]->vmax)
			{
				nfilterind++;
				PROGRESS_VALUE++;
				continue;
			}

			byte vst = 0;
			for (int64 st = 0, ed = 0; ed < nloc; ++ed)
			{
				//if ed is in a different contig, reset st
				if (strcmp(GetLoc(st).GetChrom(), GetLoc(ed).GetChrom()))
				{
					st = ed;
					vst = 0;
					continue;
				}

				GENOTYPE &gt = ainds[i]->GetGenotype(GetLocId(ed), GetLoc(ed).GetGtab());//fine
				if (gt.Nalleles() == 0) continue;

				int ved = gt.Ploidy();
				if (vst == 0) vst = ved;
				if (vst != ved)
				{
					ainds[i] = NULL;
					break;
				}
			}

			if (ainds[i]) nfilterind++;
			PROGRESS_VALUE++;
		}
	}

	/* Perform haplotype extraction */
	TARGET void CalcHaplotype()
	{
		if (abs(g_format_val) > 2 || !haplotype)
		{
			if (locus_pos)
			{
				delete[] locus_pos;
				locus_pos = NULL;
			}
			haplotype = false;
			return;
		}
		if (ad) Exit("\nError: haplotype extraction (-haplotype) is incompatible with allelic depth (-ad) option.\n");
		if (ploidyinfer) Exit("\nError: haplotype extraction (-haplotype) is incompatible with ploidy inference (-ploidyinfer) option.\n");

		//test xxx debug
		if (false)
		{
			for (int64 l = 0; l < nloc; ++l)
				slocus[l].GetChrom()[1] = '0' + (l % 10);
		}

		//Backup original locus id for sorting locus
		if (useslocus)
		{
			locus_id = new LOCN[nloc];
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
			for (int64 l = 0; l < nloc; ++l)
				locus_id[l] = l;
		}

		//Sort locus according to contig and pos, so as to iterator locus
		QSLPAR par = { 0, nloc - 1 };
		qslstack.Push(par);
		RunThreads(&QSWorker, NULL, NULL, nloc, nloc, "\nSorting loci according to contig/chromosome and position:\n", g_nthread_val, true);

		//Marker anisoploids in the same contig
		RunThreads(&CheckAnisoploid, NULL, NULL, nind + nloc, nind + nloc, "\nChecking ansioploids in the same contig:\n", 1, true);

		if (nfilterind == 0)
			Exit("\nError: all %d individuals are excluded from analysis due to the anisoploid in the same contig, please check data.\n", nind);

		//Remove anisoploids in the same contig
		if (nfilterind != nind)
			RunThreads(&RemoveIndividual, NULL, NULL, nloc, nloc, "\nFilter ansioploids in the same contig to perform haplotype extraction:\n", 1, true);

		//CreateHaplotypeLocus
		haplotype_contig = -1;
		RunThreads(&CreateHaplotypeLocus, NULL, NULL, nloc, nloc, "\nCreate haplotype from VCF file:\n", g_nthread_val, true);

		//Sort extracted locus
		int64 nl = haplotype_locus.size;
		QSLPAR par2 = { 0, nl - 1 };
		qslstack.Push(par2);
		RunThreads(&QSHapWorker, NULL, NULL, nl, nl, "\nSorting extracted loci according to contig/chromosome and position:\n", g_nthread_val, true);

		if (nl == 0) 
			Exit("\nError: Cannot extract haplotype, because all variants are not genotyped in all individuals or the conditions are too strict.\n");
		
		//Calculate new genotype index table offset
		haplotype_offset = new OFFSET[nl];
		genotype_coffset = 0;
		for (int64 l = 0; l < nl; ++l)
		{
			int64 size = CeilLog2(haplotype_locus[l].gsize);
			haplotype_offset[l].offset = genotype_coffset;
			haplotype_offset[l].size = size;
			genotype_coffset += (int64)ceil(size * nind / 8.0) + 3;
		}
		haplotype_size = genotype_coffset;
		haplotype_bucket = new byte[haplotype_size];

		//Create new locus
		MEMORY *omemory = locus_memory;
		locus_memory = new MEMORY[g_nthread_val];
		haplotype_nslocus = useslocus ? new SLOCUS[nl] : NULL;
		haplotype_nlocus = useslocus ? NULL : new LOCUS[nl];

		//Write haplotypes
		OpenResFile("-haplotype", "Extracted haplotype information");
		OpenTempFiles(g_nthread_val, ".haplotype");
		PROGRESS_VALUE = 0; PROGRESS_TOTAL = nl; PROGRESS_NOUTPUTED2 = PROGRESS_NOUTPUTED; PROGRESS_NOUTPUTED = 0;
		RunThreads(&WriteHaplotypeLocus, NULL, NULL, nl, nl, "\nWriting extracted loci:\n", g_nthread_val, true);

		ALLELE_IDENTIFIER = false;

		//Finish, release memory
		nloc = nl;
		if (useslocus)
		{
			delete[] slocus;
			slocus = haplotype_nslocus;

			delete[] locus_id;
			locus_id = NULL;

			delete[] locus_pos;
			locus_pos = NULL;
		}
		else
		{
			delete[] locus;
			locus = haplotype_nlocus;
		}
		delete[] omemory;

		haplotype_locus.~LIST();

		if (VIRTUAL_MEMORY)
			VUnAllocGenotype();
		else
			delete[] genotype_bucket;

		delete[] genotype_offset; genotype_offset = haplotype_offset;
		genotype_size = haplotype_size;
		genotype_bucket = haplotype_bucket;

		JoinTempFiles(g_nthread_val);
		CloseResFile();			

		//Using unphased genotype and no allele name saved in the locus
		haplotype = false;

		//Calculate total number of alleles and genotypes
		KT = GT = maxK = maxG = 0;
		for (int64 l = 0; l < nloc; ++l)
		{
			KT += GetLoc(l).k;
			GT += GetLoc(l).ngeno;
			maxK = Max((int)GetLoc(l).k, maxK);
			maxG = Max((int)GetLoc(l).ngeno, maxG);
		}
	}

	/* Applying individual and diversity filters */
	TARGET void ApplyFilter()
	{
		//1. Applying individual filter
		if (individual_filter)
		{
			RunThreads(&MarkerIndividual, NULL, NULL, nind + nloc, nind, "\nApplying individual filter:\n", g_nthread_val, true, (int)(0.5 * g_progress_val));
			RunThreads(&RemoveIndividual, NULL, NULL, nind + nloc, nloc, NULL, 1, false, (int)(0.5 * g_progress_val));
		}
		CheckGenotypeId();

		//Assign individuals into their vpop
		AssignVInds();
		CheckGenotypeId();

		//2. Apply diversity filter
		if (diversity_filter)
		{
			if (f_pop_b)
			{
				if (strcmp("total", f_pop_val) == 0)
					cpop = total_pop;
				else
				{
					bool find = false;
					for (int i = 0; i < npop; ++i)
						if (strcmp(pop[i].name, f_pop_val) == 0)
						{
							find = true;
							cpop = &pop[i];
						}

					if (!find) Exit("\nError: Cannot find target population %d, check parameter -f_pop.\n", f_pop_val);
				}
			}
			else if (f_region_b)
			{
				bool find = false;
				for (int rl = 0; rl < reg.size - 1; ++rl)
					for (int i = 0; i < reg[rl].size; ++i)
						if (strcmp(reg[rl][i].name, f_region_val) == 0)
						{
							find = true;
							cpop = &reg[rl][i];
						}

				if (!find) Exit("\nError: Cannot find target region %s, check parameter -f_region.\n", f_region_val);
			}
			else
				cpop = total_pop;

			//Temp use allele frequency and genotype count for diversity calculation
			allele_freq_offset = new LOCN[nloc];
			genotype_count_offset = new LOCN[nloc];
			SetFF(allele_freq_offset, nloc);
			SetFF(genotype_count_offset, nloc);

			//Calculate offset
			KT = GT = maxK = maxG = 0;
			for (int64 l = 0; l < nloc; ++l)
			{
				allele_freq_offset[l] = KT;
				genotype_count_offset[l] = GT;
				KT += GetLoc(l).k;
				GT += GetLoc(l).ngeno;
				maxK = Max((int)GetLoc(l).k, maxK);
				maxG = Max((int)GetLoc(l).ngeno, maxG);
			}

			//Calculate freq for cpop and diversity estimation
			cpop->AllocFreq();
			RunThreads(&MarkerDiversity, NULL, NULL, nloc, nloc, "\nMarking locus diversity filter:\n", 1, true);
			cpop->UnAllocFreq();

			//Release temp allele frequency and genotype count
			delete[] allele_freq_offset;
			delete[] genotype_count_offset;
			allele_freq_offset = genotype_count_offset = NULL;
		}
		CheckGenotypeId();

		//Mark monomorphic locus after previous filtering
		TABLE<int, int> *allele_table = new TABLE<int, int>[g_nthread_val];
		for (int i = 0; i < g_nthread_val; ++i)
			new(&allele_table[i]) TABLE<int, int>(true, NULL);

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			if (GetLoc(l).k < 2)
				GetLoc(l).flag_pass = 0;
			else 
			{
				threadid = omp_get_thread_num();
				TABLE<int, int> &atab = allele_table[threadid];
				atab.Clear();

				GENOTYPE *gtab = GetLoc(l).GetGtab();
				int ngeno = GetLoc(l).ngeno;
				bool pass_flag = false;

				for (int gi = 0; gi < ngeno; ++gi)
				{
					GENOTYPE &gt = gtab[gi];
					int nalleles = gt.Nalleles();
					if (nalleles == 0) continue;
					if (nalleles >= 2) { pass_flag = true; break; }
					atab[gt.GetAlleleArray()[0]] = 1;
					if (atab.size > 1) { pass_flag = true; break; }
				}

				GetLoc(l).flag_pass = pass_flag;
			}
		}
		delete[] allele_table;
		CheckGenotypeId();

		//3. Apply locus filter
		nfilter = 0;
		KT = GT = maxK = maxG = 0;
		for (int64 l = 0; l < nloc; ++l)
		{
			if (GetLoc(l).flag_pass)
			{
				nfilter++;
				KT += GetLoc(l).k;
				GT += GetLoc(l).ngeno;
				maxK = Max((int)GetLoc(l).k, maxK);
				maxG = Max((int)GetLoc(l).ngeno, maxG);
			}
		}

		CheckGenotypeId();
		if (nfilter == 0)
			Exit("\nError: all %d loci are excluded from analysis, try less strict filters.\n", nloc);

		if (nfilter != nloc)
			RunThreads(&RemoveLocus, NULL, NULL, nloc, nloc, "\nApplying locus diversity filter:\n", 1, true);

		CheckGenotypeId();
	}

	/* Recursive set vid, vpop, ind for each region */
	TARGET void SetVReg(int rl, int i)
	{
		if (rl >= 0)
		{
			if (reg[rl][i].npop == 0 || reg[rl][i].nind == 0) return;
			reg[rl][i].id = nreg[rl];
			reg[rl][i].vpop = rl >= 1 ? aregs[rl - 1] + nreg[rl - 1] : apops + npop;
			reg[rl][i].inds = rinds + nind;
			reg[rl][i].ind0id = nind;
			aregs[rl][nreg[rl]++] = &reg[rl][i];
		}

		if (rl > 0) for (int j = 0; j < reg[rl - 1].size; ++j)
		{
			if (reg[rl - 1][j].rid != i) continue;
			reg[rl - 1][j].rid = reg[rl][i].id;
			SetVReg(rl - 1, j);
		}
		else for (int j = 0, np = 0; j < pop.size; ++j)
		{
			if (pop[j].rid != i || pop[j].nind == 0) continue;

			apops[npop] = &pop[j];
			pop[j].id = npop++;
			pop[j].rid = reg[rl][i].id;
			pop[j].vpop = NULL;
			pop[j].inds = rinds + nind;
			pop[j].ind0id = nind; 
			nind += pop[j].nind;
		}
	}

	/* Sort individuals by population index to rinds array */
	TARGET void AssignVInds()
	{
		/* Original order is in ind and should not be changed to correctly obtain genotype index */

		//Assign pop/reg
		for (int i = 0; i < pop.size; ++i)
			pop[i].nind = 0;

		for (int rl = 0; rl < reg.size; ++rl)
			for (int i = 0; i < reg[rl].size; ++i)
				reg[rl][i].npop = reg[rl][i].nind = 0;

		//First scan, count inds
		for (int j = 0; j < nind; ++j)
		{
			//default pop or input pop
			int i = (ainds[j]->popid < 1 || ainds[j]->popid >= pop.size) ? 0 : ainds[j]->popid;
			pop[i].nind++;
		}

		//Allocate vreg
		nregt = nregt2 = 0;
		lreg = reg.size;
		aregs = new POP**[reg.size];
		SetZero(nreg, lreg);

		//Add up nind and npop to regions at lay 0
		for (int i = 0; i < pop.size; ++i)
		{
			if (pop[i].nind == 0) continue;
			if (reg.size > 0 && reg[0].size > 0)
			{
				reg[0][pop[i].rid].npop++;
				reg[0][pop[i].rid].nind += pop[i].nind;
			}
		}

		//Add up nind and npop to regions at lay 1+
		for (int rl = 0; rl < reg.size; ++rl)
		{
			for (int i = 0; i < reg[rl].size; ++i)
			{
				if (reg[rl][i].npop > 0)
				{
					nreg[rl]++;  nregt++;
					if (rl < reg.size - 1)
					{
						reg[rl + 1][reg[rl][i].rid].npop++;
						reg[rl + 1][reg[rl][i].rid].nind += reg[rl][i].nind;
					}
				}
			}
			nregt2 += nreg[rl] * nreg[rl];
		}

		//Remove top lay while there is only one region in this lay
		while (lreg && nreg[lreg - 1] == 1)
		{
			nregt--;
			nregt2--;
			lreg--;
			if (lreg == 0) break;
		}

		//Allocate vpop
		npop = 0;
		for (int i = 0; i < pop.size; ++i)
			if (pop[i].nind > 0)
				npop++;
		apops = new POP*[npop];

		//Contruct total pop
		if (lreg == 0 && npop == 1)
		{
			lreg = -1;
			for (int i = 0; i < pop.size; ++i)
				if (pop[i].nind > 0)
					total_pop = &pop[i];
			total_pop->name = new char[6];
			nregt2 = nregt = 0;
		}
		else
			total_pop = &reg[lreg][0];

		sprintf(total_pop->name, "Total");
		total_pop->rid = -1;

		for (int rl = 0; rl <= lreg; ++rl)
		{
			aregs[rl] = new POP*[nreg[rl]];
			nreg[rl] = 0;
		}

		//Allocate rinds
		rinds = new IND*[nind]; 

		//Recursive arrange vpop, vreg and vind
		if (lreg == -1)
		{
			apops[0] = total_pop;
			total_pop->id = 0;
			total_pop->rid = -1;
			total_pop->vpop = NULL;
			total_pop->inds = rinds;
			total_pop->ind0id = 0;
		}
		else
		{
			nind = 0; npop = 0;
			SetVReg(lreg, 0);
		}


		//Set nind = 0
		for (int i = 0; i < pop.size; ++i)
			pop[i].nind = 0;

		//Move ind to pop ind array and assign vpopid
		bool swap = false;
		for (int j = 0; j < nind; ++j)
		{
			//Default pop or input pop
			IND& i = *ainds[j];
			POP& p = pop[(i.popid < 1 || i.popid >= pop.size) ? 0 : i.popid];
			if (i.indid != p.ind0id + p.nind)
				swap = true;
			i.indid = p.ind0id + p.nind;
			p.inds[p.nind++] = ainds[j];
			i.popid = p.id;
		}

		if (!swap) return;

		//swap genotype
		ushort *gtbuf = new ushort[g_nthread_val * nind];

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			threadid = omp_get_thread_num();
			ushort *gtb = gtbuf + threadid * nind;

			GENO_ITERATOR rg(0u, l, true);//OK
			GENO_ITERATOR wg(0u, l, false, genotype_bucket, genotype_offset);//OK

			for (int i = 0; i < nind; ++i)
				gtb[ainds[i]->indid] = rg.Read();

			for (int i = 0; i < nind; ++i)
				wg.Write(gtb[i]);

			wg.FinishWrite();
		}

		delete[] gtbuf;

		if (ploidyinfer)
		{
			//swap allelic depth
			uint *dpbuf = new uint[g_nthread_val * nind * maxK];

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
			for (int64 l = 0; l < nloc; ++l)
			{
				threadid = omp_get_thread_num();
				uint *dpb = dpbuf + threadid * nind * maxK;

				GENO_ITERATOR rd(0u, l,  true, alleledepth_bucket, alleledepth_offset);//OK
				GENO_ITERATOR wd(0u, l, false, alleledepth_bucket, alleledepth_offset);//OK
				int k2 = GetLoc(l).k;

				for (int i = 0; i < nind; ++i)
				{
					uint *dpb2 = dpb + ainds[i]->indid * k2;

					for (int k = 0; k < k2; ++k)
						dpb2[k] = rd.Read();
				}

				for (int i = 0, nend = k2 * nind; i < nend; ++i)
					wd.Write(*dpb++);
			}

			delete[] dpbuf;
		}

		delete[] ainds; ainds = rinds; 
	}

	/* Calculate individual min and max ploidy, and sum ploidy levels */
	TARGET void AssignPloidy()
	{
		RunThreads(&AssignPloidyThread, NULL, NULL, nloc, nloc,
			"\nAssigning individual ploidy:\n", 1, true);

		//Calculate the total number of haplotypes in each population
		for (int i = 0; i < npop; ++i)
			apops[i]->nhaplotypes = 0;

		for (int i = 0; i < nind; ++i)
			apops[ainds[i]->popid]->nhaplotypes += ainds[i]->vmax;

		for (int rl = 0; rl <= lreg; ++rl)
			for (int i = 0; i < nreg[rl]; ++i)
			{
				POP *tp = aregs[rl][i];
				tp->nhaplotypes = 0;
				for (int j = 0; j < tp->npop; ++j)
					tp->nhaplotypes += tp->vpop[j]->nhaplotypes;
			}

		reassigned = true;
	}

	/* Calculate allele frequencies for each population and region for further use */
	TARGET void CalcFreq()
	{
		//Calculate allele frequency and genotype count
		RunThreads(&CalcAlleleFreq, NULL, NULL, nloc * (int64)nind * (lreg + 2), nloc * (int64)nind * (lreg + 2),
			"\nCalculating allele frequencies:\n", 1, true);
	}

	/* Convert into other genotype formats */
	TARGET void ConvertFile()
	{
		if (!convert) return;
		if (ad) Exit("\nError: file convertion (-convert) is incompatible with allelic depth (-ad) option.\n");

		convert_buf = new char*[NBUF];
		SetZero(convert_buf, NBUF);

		bool isfirst = true;
		int ntot = 0;
		for (int i = 1; i <= N_CONVERTER; ++i)
			if (convert_format_val[i])
				ntot += nind;

		//Alloc memory to buffer genotype string
		conversion_memory = new MEMORY[g_nthread_val];
		conversion_string = (LIST<char*>*)malloc(nloc * sizeof(LIST<char*>));

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			threadid = omp_get_thread_num();
			new(&conversion_string[l]) LIST<char*>(&conversion_memory[threadid]);
		}

		for (int i = 1; i <= N_CONVERTER; ++i)
		{
			if (convert_format_val[i] == 0) continue;
			conversion_memory2 = new MEMORY[g_nthread_val];

			switch (i)
			{
			case 1: ConvertGenepop(ntot, isfirst); break;
			case 2: ConvertSpagedi(ntot, isfirst); break;
			case 3: ConvertCervus(ntot, isfirst); break;
			case 4: ConvertArlequin(ntot, isfirst); break;
			case 5: ConvertStructure(ntot, isfirst); break;
			case 6: ConvertPolygene(ntot, isfirst); break;
			case 7: ConvertPolyrelatedness(ntot, isfirst); break;
			}

			delete[] conversion_memory2;
		}

		delete[] conversion_memory;
		free(conversion_string);
		delete[] convert_buf;
	}

	/* Calculate genetic diveristy indices */
	TARGET void CalcDiversity()
	{
		if (!diversity) return;
		OpenResFile("-diversity", "Individual statistics");
		OpenTempFiles(2, ".diversity");

		diversity_buf = new DIVERSITY[NBUF];

		bool isfirst = true;
		int64 ntot = 0;
		if (diversity_level_val[1] || diversity_level_val[4]) ntot += nind * nloc;
		if ((diversity_level_val[2] || diversity_level_val[5]) && lreg > 0) ntot += lreg * nind * nloc;
		if (diversity_level_val[3] || diversity_level_val[6]) ntot += nind * nloc;

		//Total population
		if (diversity_level_val[3] || diversity_level_val[6])
		{
			diversity_stage = 3;
			cpop = total_pop;

			if (diversity_level_val[3]) DIVERSITY::WriteHeader(TEMP_FILES[0]);
			if (diversity_level_val[6]) DIVERSITY::WriteHeader(TEMP_FILES[1]);

			RunThreads(&DiversityThread, &DiversityGuard, NULL, ntot, cpop->nind * nloc,
				"\nCalculating genetic diversity:\n", g_nthread_val, isfirst);
			isfirst = false;

			if (diversity_level_val[3]) diversity_sum.Write(TEMP_FILES[0], "Total");
		}

		//Regions
		if (diversity_level_val[2] || diversity_level_val[5])
		{
			diversity_stage = 2;
			if (diversity_level_val[2]) DIVERSITY::WriteHeader(TEMP_FILES[0]);
			if (diversity_level_val[5]) DIVERSITY::WriteHeader(TEMP_FILES[1]);

			for (int rl = lreg - 1; rl >= 0; --rl)
				for (int i = 0; i < nreg[rl]; ++i)
				{
					cpop = aregs[rl][i];
					RunThreads(&DiversityThread, &DiversityGuard, NULL, ntot, cpop->nind * nloc,
						"\nCalculating genetic diversity:\n", g_nthread_val, isfirst);
					isfirst = false;

					if (diversity_level_val[2]) diversity_sum.Write(TEMP_FILES[0], cpop->name);
				}
		}

		//Populations
		if (diversity_level_val[1] || diversity_level_val[4])
		{
			diversity_stage = 1;
			if (diversity_level_val[1]) DIVERSITY::WriteHeader(TEMP_FILES[0]);
			if (diversity_level_val[4]) DIVERSITY::WriteHeader(TEMP_FILES[1]);

			for (int i = 0; i < npop; ++i)
			{
				cpop = apops[i];
				RunThreads(&DiversityThread, &DiversityGuard, NULL, ntot, cpop->nind * nloc,
					"\nCalculating genetic diversity:\n", g_nthread_val, isfirst);
				isfirst = false;

				if (diversity_level_val[1]) diversity_sum.Write(TEMP_FILES[0], cpop->name);
			}
		}

		JoinTempFiles(2);
		CloseResFile();
		CheckGenotypeId();

		delete[] diversity_buf;
	}

	/* Calculate individual statistics */
	TARGET void CalcIndstat()
	{
		if (!indstat) return;
		if (ad) Exit("\nError: individual statistics (-indstat) is incompatible with allelic depth (-ad) option.\n");
		OpenResFile("-indstat", "Individual statistics");
		OpenTempFiles(g_nthread_val, ".indstat");

		IND::IndividualStatisticsHeader(FRES);
		RunThreads(&IndividualStatisticsThread, NULL, NULL, nind, nind,
			"\nCalculating individual statistics:\n", g_nthread_val, true);

		JoinTempFiles(g_nthread_val);
		CloseResFile();
	}

	/* Calculate genetic differentiation */
	TARGET void CalcDiff()
	{
		if (!fst) return;
		OpenResFile("-fst", "Genetic differentiation");
		OpenTempFiles(1, ".fst");

		bool isfirst = true; int64 ntot = 0;
		if (fst_level_val[1] && lreg > 0) ntot += lreg;
		if (fst_level_val[2]) ntot += 1;
		if (fst_level_val[3] && lreg > 0) ntot += nregt;
		if (fst_level_val[4] && lreg > 0) ntot += (nregt2 - nregt) / 2;
		if (fst_level_val[5]) ntot += npop * (npop - 1) / 2;

		for (fst_type = 1; fst_type <= 5; ++fst_type)
		{
			if (fst_level_val[fst_type] == 0) continue;
			if ((fst_type == 1 || fst_type == 4 || fst_type == 3) && lreg == 0) continue;

			int n1 = 0, n2 = 0;
			switch (fst_type)
			{
			case 1: n1 = lreg >= 0 ? lreg : 0;  n2 = lreg >= 0 ? lreg : 0; break;
			case 2: n1 = 1;     n2 = 1;    break;
			case 3: n1 = nregt; n2 = nregt; break;
			case 4: n1 = nregt2; n2 = (nregt2 - nregt) / 2;  break;
			case 5: n1 = npop * npop; n2 = npop * (npop - 1) / 2; break;
			}

			if (n1 == 0 || n2 == 0) continue;
			fst_buf[fst_type] = new FST[n1];
			SetZero(fst_buf[fst_type], n1);
			RunThreads(&GeneticDifferentiationThread, NULL, NULL, ntot, n2,
				"\nCalculating genetic differentiation:\n", g_nthread_val, isfirst);
			isfirst = false;

			if ((fst_type == 4 || fst_type == 5) && fst_fmt_val[1])
				FST::MatrixPrint(TEMP_FILES[0], fst_buf[fst_type], npop, fst_type);
		}

		FST::ColumnPrint(FRES);
		JoinTempFiles(1);
		CloseResFile();

		for (fst_type = 1; fst_type <= 5; ++fst_type)
		{
			if (fst_level_val[fst_type] == 0) continue;
			if ((fst_type == 1 || fst_type == 4) && lreg == 0) continue;

			int n1 = 0, n2 = 0;
			switch (fst_type)
			{
			case 1: n1 = lreg;  n2 = lreg; break;
			case 2: n1 = 1;     n2 = 1;    break;
			case 3: n1 = nregt; n2 = nregt; break;
			case 4: n1 = nregt2; n2 = (nregt2 - nregt) / 2;  break;
			case 5: n1 = npop * npop; n2 = npop * (npop - 1) / 2; break;
			}

			if (fst_type <= 3)
			{
				for (int ni = 0; ni < n1; ++ni)
					fst_buf[fst_type][ni].Uninitialize();
			}
			else if (fst_type == 4)
			{
				auto fst_buf2 = fst_buf[fst_type];
				for (int rl = lreg - 1; rl >= 0; --rl)
				{
					for (int i = 0; i < nreg[rl]; ++i)
						for (int j = i + 1; j < nreg[rl]; ++j)
							fst_buf2[i * nreg[rl] + j].Uninitialize();
					fst_buf2 += nreg[rl] * nreg[rl];
				}
			}
			else if (fst_type == 5)
			{
				auto fst_buf2 = fst_buf[fst_type];
				for (int i = 0; i < npop; ++i)
					for (int j = i + 1; j < npop; ++j)
						fst_buf2[i * npop + j].Uninitialize();
			}
			delete[] fst_buf[fst_type];
			fst_buf[fst_type] = NULL;
		}
	}

	/* Calculate genetic distance */
	TARGET void CalcDist()
	{
		if (!gdist) return;
		GDIST_METHOD = 1;
		OpenResFile("-gdist", "Genetic distance");

		bool isfirst = true;
		gdist_buf = new GDIST[NBUF];

		int64 ntot = 0;
		int64 n[] = { 0, nind * nind * 2, npop * npop * 2 * 100 };

		for (int k = 1; k <= 2; ++k)
			if (gdist_level_val[k])
				ntot += n[k];

		if (gdist_level_val[3])
			for (int rl = 0; rl < lreg; ++rl)
				ntot += nreg[rl] * nreg[rl] * 2 * 100;

		for (int k = 1; k <= 3; ++k)
		{
			if (gdist_level_val[k] == 0)
				continue;

			//Calculate genetic distance table between any two genotypes
			if (k == 1)
			{
				gd_memory = new MEMORY[g_nthread_val];
				gdtab = new TABLE<HASH, INDGD>[nloc];
				gdlock = new shared_mutex[nloc];

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
				for (int64 l = 0; l < nloc; ++l)
				{
					threadid = omp_get_thread_num();
					int ng = GetLoc(l).ngeno;
					if (ng < N_MAX_GDTAB)
						new(&gdtab[l]) TABLE<HASH, INDGD>(false, &gd_memory[threadid], (int)(Binomial(ng, 2) + ng + 0.5));
				}
			}

			for (int rl = 0; rl < (k <= 2 ? 1 : lreg); ++rl)
			{
				int64 nthis = k <= 2 ? n[k] : nreg[rl] * nreg[rl] * 2 * 100;
				gdist_estimator_val[0] = true;
				OpenTempFiles(N_GD_ESTIMATOR + 1, ".gdist", gdist_estimator_val);
				gdist_type = k + rl;

				SetZero(gdist_buf, NBUF);

				RunThreads(&GeneticDistanceThread, &GeneticDistanceGuard2, &GeneticDistanceGuard1, ntot, nthis,
					"\nCalculating genetic distance:\n", g_nthread_val, isfirst);
				
				isfirst = false;
				JoinTempFiles(N_GD_ESTIMATOR + 1, gdist_estimator_val);
				gdist_estimator_val[0] = false;
			}

			if (k == 1)
			{
				delete[] gdtab;
				delete[] gdlock;
				delete[] gd_memory;
			}
		}
		delete[] gdist_buf;
		GDIST_METHOD = 0;
		CloseResFile();
	}

	/* Calculate analysis of molecular variance */
	TARGET void CalcAMOVA()
	{
		if (!amova) return;
		if (ad) Exit("\nError: AMOVA (-amova) is incompatible with allelic depth (-ad) option.\n");

		OpenResFile("-amova", "Analysis of molecular variances");

		bool isfirst = true;
		int64 nt = 0;
		for (int i1 = 1; i1 <= 3; ++i1) if (amova_method_val[i1])
			for (int i2 = 1; i2 <= (i1 == 3 ? 1 : 2); ++i2) if (amova_mutation_val[i2])
				for (int i3 = 1; i3 <= 2; ++i3) if (amova_ind_val[i3])
				{
					int64 nperm = amova_nperm_val * (amova_test_val == 1);

					switch (i1)
					{
					case 1: nt += (int64)(nloc * 10 + Binomial(total_pop->nhaplotypes, 2) + nperm * BINOMIAL[(i3 == 1) + 2 + lreg][2] * 100 + 0.5); break;
					case 2: 
					{
						if (amova_pseudo_val > 0 && amova_pseudo_val < 10)
							Exit("\nError: In anisoploid AMOVA method, the number of pseudo-permutations should be greater than 10. \n");

						bool PseudoPerm = amova_pseudo_val > 0;
						int64 M = Min((int)Max(pow(amova_nperm_val * 10000, 1.0 / nloc), (double)amova_pseudo_val), amova_nperm_val) * (amova_test_val == 1);
						int64 nperm = amova_nperm_val * (amova_test_val == 1);
						if (!PseudoPerm) M = nperm;
						int npair = BINOMIAL[(i3 == 1) + 2 + lreg][2] + 0.5;
						nt += nloc * M * npair * 50 + nperm * npair;
						break;
					}
					case 3: nt += (int64)(nperm * BINOMIAL[2 + lreg][2] * 500 + 0.5); break;
					}
				}

		PROGRESS_VALUE2 = PROGRESS_VALUE3 = 0;

		for (int i1 = 1; i1 <= 3; ++i1) if (amova_method_val[i1])
			for (int i2 = 1; i2 <= (i1 == 3 ? 1 : 2); ++i2) if (amova_mutation_val[i2])
				for (int i3 = 1; i3 <= 2; ++i3) if (amova_ind_val[i3])
				{
					int nthread = 1;
					int64 nperm = amova_nperm_val * (amova_test_val == 1);

					amova_cmethod_val = i1;
					amova_cmutation_val = i2;
					amova_cind_val = i3;
					amova_buf = new AMOVA();

					int64 nced = 0;
					switch (i1)
					{
					case 1: nced += (int64)(nloc * 10 + Binomial(total_pop->nhaplotypes, 2) + nperm * BINOMIAL[(i3 == 1) + 2 + lreg][2] * 100 + 0.5); break;
					case 2: 
					{
						bool PseudoPerm = amova_pseudo_val > 0;
						int64 M = Min((int)Max(pow(amova_nperm_val * 10000, 1.0 / nloc), (double)amova_pseudo_val), amova_nperm_val) * (amova_test_val == 1);
						int64 nperm = amova_nperm_val * (amova_test_val == 1);
						if (!PseudoPerm) M = nperm;
						int npair = BINOMIAL[(i3 == 1) + 2 + lreg][2] + 0.5;
						nced += nloc * M * npair * 50 + nperm * npair;
						break;
					}
					case 3: nced += (int64)(nperm * BINOMIAL[2 + lreg][2] * 500 + 0.5); break;
					}

					RunThreads(&AMOVAThread, NULL, NULL, nt, nced,
						"\nPerforming analysis of molecular variance:\n", nthread, isfirst);

					isfirst = false;
					int nlay = amova_buf[0].nlay;
					for (int i = 1; i < nthread; ++i)
					{
						Add(amova_buf[0].G, amova_buf[i].G, nlay * nlay);
						Add(amova_buf[0].E, amova_buf[i].E, nlay * nlay);
						Add(amova_buf[0].EF2, amova_buf[i].EF2, nlay * nlay);
						if (amova_buf[i].nSSW) for (int j = 0; j < nlay; ++j)
							Add(amova_buf[0].SSW[j], amova_buf[i].SSW[j], amova_buf[i].nSSW[j]);
					}
					amova_buf[0].PrintAMOVA(FRES);

					delete amova_buf;
				}

		CloseResFile();
	}

	/* Calculate population assignment */
	TARGET void CalcAssignment()
	{
		if (!popas) return;
		if (ad) Exit("\nError: population assignment (-popas) is incompatible with allelic depth (-ad) option.\n");
		OpenResFile("-popas", "Population assignment");
		OpenTempFiles(g_nthread_val, ".popas");

		IND::AssignmentHeader(FRES);
		RunThreads(&PopulationAssignmentThread, NULL, NULL, nind, nind,
			"\nPerforming population assignment:\n", g_nthread_val, true);

		JoinTempFiles(g_nthread_val);
		CloseResFile();
	}

	/* Calculate relatedness coefficient */
	TARGET void CalcRelatedness()
	{
		if (!relatedness) return;
		if (ad) Exit("\nError: relatedness estimation (-relatedness) is incompatible with allelic depth (-ad) option.\n");
		if (relatedness_estimator_val[11]) RELATEDNESS::Huang2015_Initialize();
		OpenResFile("-relatedness", "Relatedness coefficient");

		setNbThreads(1);

		bool isfirst = true;
		relatedness_buf = new RELATEDNESS[NBUF];
		int64 ntot = 0;
		if (relatedness_range_val[1]) for (int i = 0; i < npop; ++i)
			ntot += apops[i]->nind * apops[i]->nind;
		if (relatedness_range_val[2])
			for (int rl = 0; rl < lreg; ++rl)
				for (int i = 0; i < nreg[rl]; ++i)
					ntot += aregs[rl][i]->nind * aregs[rl][i]->nind;
		if (relatedness_range_val[3])
			ntot += total_pop->nind * total_pop->nind;
		ntot <<= 1;

		for (int k = 1; k <= 3; ++k)
		{
			if (relatedness_range_val[k] == 0) continue;
			for (int rl = 0; rl < (k == 2 ? lreg : 1); ++rl)
			{
				int n = k == 1 ? npop : (k == 2 ? nreg[rl] : 1);
				for (int i = 0; i < n; ++i)
				{
					OpenTempFiles(N_RELATEDNESS_ESTIMATOR + 1, ".relatedness");
					cpop = k == 1 ? apops[i] : (k == 2 ? aregs[rl][i] : total_pop);
					SetZero(relatedness_buf, NBUF);

					RunThreads(&RelatednessThread, &RelatednessGuard1, &RelatednessGuard2, ntot, cpop->nind * cpop->nind * 2,
						"\nEstimating relatedness coefficient between individuals:\n", g_nthread_val, isfirst);
					isfirst = false;
					JoinTempFiles(N_RELATEDNESS_ESTIMATOR + 1);
				}
			}
		}

		delete[] relatedness_buf;
		CloseResFile();
		if (relatedness_estimator_val[11]) RELATEDNESS::Huang2015_Uninitialize();
	}

	/* Calculate kinship coefficient */
	TARGET void CalcKinship()
	{
		if (!kinship) return;
		if (ad) Exit("\nError: kinship estimation (-kinship) is incompatible with allelic depth (-ad) option.\n");
		OpenResFile("-kinship", "Kinship coefficient");

		bool isfirst = true;
		kinship_buf = new KINSHIP[NBUF];
		int64 ntot = 0;
		if (kinship_range_val[1]) for (int i = 0; i < npop; ++i)
			ntot += apops[i]->nind * apops[i]->nind;
		if (kinship_range_val[2])
			for (int rl = 0; rl < lreg; ++rl)
				for (int i = 0; i < nreg[rl]; ++i)
					ntot += aregs[rl][i]->nind * aregs[rl][i]->nind;
		if (kinship_range_val[3])
			ntot += total_pop->nind * total_pop->nind;
		ntot <<= 1;

		for (int k = 1; k <= 3; ++k)
		{
			if (kinship_range_val[k] == 0) continue;
			for (int rl = 0; rl < (k == 2 ? lreg : 1); ++rl)
			{
				int n = k == 1 ? npop : (k == 2 ? nreg[rl] : 1);
				for (int i = 0; i < n; ++i)
				{
					OpenTempFiles(N_KINSHIP_ESTIMATOR + 1, ".kinship");
					cpop = k == 1 ? apops[i] : (k == 2 ? aregs[rl][i] : total_pop);
					SetZero(kinship_buf, NBUF);

					RunThreads(&KinshipThread, &KinshipGuard1, &KinshipGuard2, ntot, cpop->nind * cpop->nind * 2,
						"\nEstimating kinship coefficient between individuals:\n", g_nthread_val, isfirst);
					isfirst = false;
					JoinTempFiles(N_KINSHIP_ESTIMATOR + 1);
				}
			}
		}

		delete[] kinship_buf;
		CloseResFile();
	}

	/* Calculate principal coordinate analysis */
	TARGET void CalcPCOA()
	{
		if (!pcoa) return;
		GDIST_METHOD = 2;
		OpenResFile("-pcoa", "Principal coordinate analysis");

		bool isfirst = true;
		int ntot = 0;
		if (pcoa_level_val[1]) ntot += nind * (nind - 1) / 2;
		if (pcoa_level_val[2]) ntot += npop * (npop - 1) / 2;
		if (pcoa_level_val[3])
			for (int rl = 0; rl < lreg; ++rl)
				ntot += nreg[rl] * (nreg[rl] - 1) / 2;

		setNbThreads(g_nthread_val);

		int nestimator = 0;
		for (int k = 1; k <= N_GD_ESTIMATOR; ++k)
			if (pcoa_estimator_val[k])
				gdindex[k] = nestimator++;

		for (int m = 1; m <= 3; ++m)
		{
			if (pcoa_level_val[m] == 0) continue;
			for (int rl = 0; rl < (m < 3 ? 1 : lreg); ++rl)
			{
				PCOA tpcoa;
				gdindex[0] = m + rl;
				int n = m == 1 ? nind : (m == 2 ? npop : nreg[rl]);

				pcoa_matrix = new double[n * n * nestimator];
				SetZero(pcoa_matrix, n * n * nestimator);

				//Calculate genetic distance table between any two genotypes
				if (m == 1)
				{
					gd_memory = new MEMORY[g_nthread_val];
					gdtab = new TABLE<HASH, INDGD>[nloc];
					gdlock = new shared_mutex[nloc];

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
					for (int64 l = 0; l < nloc; ++l)
					{
						threadid = omp_get_thread_num();
						int ng = GetLoc(l).ngeno;
						if (ng < N_MAX_GDTAB)
							new(&gdtab[l]) TABLE<HASH, INDGD>(false, &gd_memory[threadid], (int)(Binomial(ng, 2) + ng + 0.5));
					}
				}

				RunThreads(&PCoAClusteringThread, NULL, NULL, ntot, n * (n - 1) / 2,
					"\nPerforming principal coordinate analysis:\n", g_nthread_val, isfirst);
				isfirst = false;

				for (int k = 1; k <= N_GD_ESTIMATOR; ++k)
					if (pcoa_estimator_val[k])
						tpcoa.PrintPCoA(FRES, pcoa_matrix + gdindex[k] * n * n, n, k, m + rl);

				if (m == 1)
				{
					delete[] gdtab;
					delete[] gdlock;
					delete[] gd_memory;
				}
				delete[] pcoa_matrix;
			}
		}

		CloseResFile();
		GDIST_METHOD = 0;
	}

	/* Calculate hierarchical clustering */
	TARGET void CalcClustering()
	{
		if (!cluster) return;
		GDIST_METHOD = 3;
		OpenResFile("-cluster", "#Hierarchical clustering");

		bool isfirst = true;
		int ntot = 0;
		if (cluster_level_val[1]) ntot += nind * (nind - 1) / 2;
		if (cluster_level_val[2]) ntot += npop * (npop - 1) / 2;
		if (cluster_level_val[3])
			for (int rl = 0; rl < lreg; ++rl)
				ntot += nreg[rl] * (nreg[rl] - 1) / 2;

		int nestimator = 0;
		for (int k = 1; k <= N_GD_ESTIMATOR; ++k)
			if (cluster_estimator_val[k])
				gdindex[k] = nestimator++;
		const char *level[] = { "", "individual", "population", "region level" };

		for (int m = 1; m <= 3; ++m)
		{
			if (cluster_level_val[m] == 0) continue;
			for (int rl = 0; rl < (m == 3 ? lreg : 1); ++rl)
			{
				gdindex[0] = m + rl;
				int n = m == 1 ? nind : (m == 2 ? npop : nreg[rl]);

				//Calculate genetic distance table between any two genotypes
				gd_memory = new MEMORY[g_nthread_val];
				if (m == 1)
				{
					gdtab = new TABLE<HASH, INDGD>[nloc];
					gdlock = new shared_mutex[nloc];

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
					for (int64 l = 0; l < nloc; ++l)
					{
						threadid = omp_get_thread_num();
						int ng = GetLoc(l).ngeno;
						if (ng < N_MAX_GDTAB)
							new(&gdtab[l]) TABLE<HASH, INDGD>(false, &gd_memory[threadid], (int)(Binomial(ng, 2) + ng + 0.5));
					}
				}

				pcoa_matrix = new double[n * n * nestimator];
				SetZero(pcoa_matrix, n * n * nestimator);

				RunThreads(&PCoAClusteringThread, NULL, NULL, ntot, n * (n - 1) / 2,
					"\nPerforming hierarchical clustering:\n", g_nthread_val, isfirst);
				isfirst = false;

				for (int k = 1; k <= N_GD_ESTIMATOR; ++k)
				{
					if (cluster_estimator_val[k] == 0) continue;
					for (int64 l = 1; l <= N_CLUSTER_METHOD; ++l)
					{
						if (cluster_method_val[l] == 0) continue;
						switch (m)
						{
						case 1:
						{
							fprintf(FRES, "#Level: %s, genetic distance: %s, method: %s%s", level[m], GD_ESTIMATOR[k], CLUSTER_METHOD[l], g_linebreak_val);
							HCLUSTERING tcluster(pcoa_matrix + gdindex[k] * n * n, ainds, n, l, gd_memory);
							tcluster.PrintClustering(FRES);
							break;
						}
						case 2:
						{
							fprintf(FRES, "#Level: %s, genetic distance: %s, method: %s%s", level[m], GD_ESTIMATOR[k], CLUSTER_METHOD[l], g_linebreak_val);
							HCLUSTERING tcluster(pcoa_matrix + gdindex[k] * n * n, apops, n, l, gd_memory);
							tcluster.PrintClustering(FRES);
							break;
						}
						case 3:
						{
							for (int rl2 = 0; rl2 < lreg; ++rl2)
							{
								fprintf(FRES, "#Level: %s %d, genetic distance: %s, method: %s%s", level[m], rl2 + 1, GD_ESTIMATOR[k], CLUSTER_METHOD[l], g_linebreak_val);
								HCLUSTERING tcluster(pcoa_matrix + gdindex[k] * n * n, aregs[rl2], n, l, gd_memory);
								tcluster.PrintClustering(FRES);
							}
							break;
						}
						}
					}
				}
				if (m == 1)
				{
					delete[] gdtab;
					delete[] gdlock;
				}
				delete[] gd_memory;
				delete[] pcoa_matrix;
			}
		}

		GDIST_METHOD = 0;
		CloseResFile();
	}

	/* Calculate spatical pattern */
	TARGET void CalcSPA()
	{
		//locus parameter
		if (!spa) return;
		if (spa_level_val == 1 && ad)
			Exit("\nError: SPA analysis at individual level (-spa_level=ind) is incompatible with allelic depth (-ad) option.\n");
		if (spa_level_val == 1 && nind <= spa_dim_val)
			Exit("\nError: the number of individuals should be greater than the dimension of coordinate %d.\n", spa_dim_val);
		if (spa_level_val == 2 && npop <= spa_dim_val)
			Exit("\nError: the number of populations should be greater than the dimension of coordinate %d.\n", spa_dim_val);

		OpenResFile("-spa", "Analysis of spatial structure");
		OpenTempFiles(g_nthread_val, ".spa");

		spa_n = spa_level_val == 1 ? nind : npop;
		spa_np = spa_dim_val + 1;
		spa_x = new double[spa_n * spa_np];
		SetVal(spa_x, -123456789.0, spa_n * spa_np);

		//Load coordinates
		char *t1 = spa_coord_val, *t2 = 0;
		while (*t1)
		{
			while (*t1 == '\n' || *t1 == '\r' || *t1 == ' ' || *t1 == '\t') t1++;
			t2 = StrNextIdx(t1, ":", 1);

			if (spa_level_val == 1)
			{
				bool flag = false;
				for (int i = 0; i < nind; ++i)
				{
					if (memcmp(ainds[i]->name, t1, t2 - t1) == 0 && (int)strlen(ainds[i]->name) == (int)(t2 - t1))
					{
						if (spa_x[ainds[i]->indid * spa_np] != -123456789.0)
							Exit("\nError: the coordinate of individual %s appear twice.\n", t1);
						flag = true;
						t1 = t2 + 1;
						for (int j = 0; j < spa_dim_val; ++j)
						{
							while (*t1 != '-' && *t1 != '.' && !(*t1 >= '0' && *t1 <= '9')) t1++;
							spa_x[ainds[i]->indid * spa_np + j] = ReadDouble(t1);
						}
						break;
					}
				}

				if (!flag)
				{
					*t2 = 0;
					Exit("\nError: cannot find individual %s.\n", t1);
				}
			}
			else
			{
				bool flag = false;
				for (int p = 0; p < npop; ++p)
				{
					if (memcmp(apops[p]->name, t1, t2 - t1) == 0 && (int)strlen(apops[p]->name) == (int)(t2 - t1))
					{
						if (spa_x[apops[p]->id * spa_np] != -123456789.0)
							Exit("\nError: the coordinate of population %s appear twice.\n", t1);
						flag = true;
						t1 = t2 + 1;
						for (int j = 0; j < spa_dim_val; ++j)
						{
							while (*t1 != '-' && *t1 != '.' && !(*t1 >= '0' && *t1 <= '9')) t1++;
							spa_x[apops[p]->id * spa_np + j] = ReadDouble(t1);
						}
						break;
					}
				}

				if (!flag)
				{
					*t2 = 0;
					Exit("\nError: cannot find population %s.\n", t1);
				}
			}
		}

		for (int i = 0; i < spa_n; ++i)
		{
			if (spa_x[i * spa_np] == -123456789.0)
				Exit("\nError: the coordinate of %s %s is abscent.\n", spa_level_val == 1 ? "individual" : "population", ainds[i]->name);
			spa_x[i * spa_np + spa_dim_val] = 1.0;
		}

		RunThreads(&SPAThread, NULL, NULL, nloc, nloc, "\nCalculating SPA:\n", g_nthread_val, true);
		JoinTempFiles(g_nthread_val);
		CloseResFile();
		delete[] spa_x;
	}

	/* Calculate bayesian clustering */
	TARGET void CalcBayesian()
	{
		if (!structure) return;
		if (ad) Exit("\nError: Bayesian clustering (-structure) is incompatible with allelic depth (-ad) option.\n");

		//alloc tasks
		structure_totalruns = (structure_krange_max + 1 - structure_krange_min) * structure_nruns_val;
		structure_par = new SRUNINFO[structure_totalruns];
		int lc = 0;

		for (int i = structure_krange_min; i <= structure_krange_max; ++i)
			for (int j = 1; j <= structure_nruns_val; ++j)
			{
				structure_par[lc].flag.clear();
				structure_par[lc].k = i;
				structure_par[lc].id = lc + 1;
				structure_par[lc].rep = j;
				structure_par[lc].MeanlnL = 0;
				structure_par[lc].VarlnL = 0;
				structure_par[lc].lnPD = 0;
				lc++;
			}

		int nrun = ((structure_locpriori_val == 2 && structure_admix_val == 2) ? structure_nadmburnin_val : 0) + structure_nreps_val + structure_nburnin_val;
		RunThreads(&StructureThread, NULL, NULL, structure_totalruns * nrun, structure_totalruns * nrun,
			"\nPerforming bayesian clustering:\n", g_nthread_val, true);
		STRUCTURE::PrintSummary(structure_par, structure_totalruns);
		delete[] structure_par;
	}

	/* Calculate individual ploidy inference */
	TARGET void CalcPloidyInference()
	{
		if (!ploidyinfer) return;
		if (ad) Exit("\nError: ploidy inference (-ploidyinfer) is incompatible with allelic depth (-ad) option.\n");
		if (haplotype) Exit("\nError: ploidy inference (-ploidyinfer) is incompatible with haplotype extraction (-haplotype) option.\n");
		OpenResFile("-ploidyinfer", "Ploidy inference");
		OpenTempFiles(g_nthread_val, ".ploidyinfer");

		IND::PloidyInferenceHeader(FRES);
		RunThreads(&PloidyInferenceThread, NULL, NULL, nind, nind,
			"\nCalculating individual ploiy level:\n", g_nthread_val, true);

		JoinTempFiles(g_nthread_val);
		CloseResFile();
	}

	/* Calculate various analyses */
	TARGET void Calculate()
	{
		Initialize();

		LoadFile();

		// 1. Convert
		ApplyFilter();//allow ad

		// 3. Haplotype
		CalcHaplotype();//forbid ad

		AssignPloidy();//allow ad

		CalcFreq();//allow ad

		//estimate pes model
		if ((diversity && diversity_model_val[4] == 1) || (indstat && indstat_model_val[4] == 1) || (popas && popas_model_val[4] == 1))
			RunThreads(&GetLocusPESModel, NULL, NULL, nloc, nloc, "\nEstimating double-reduction rate for PES model:\n", 1, true);

		// 4. Convert
		ConvertFile(); //forbid ad, circle buf

		// 5. Genetic diversity
		CalcDiversity(); //allow ad, circle buf

		// 6. Individual statistics
		CalcIndstat(); //forbid ad, circle buf

		// 7. Genetic differentiation
		CalcDiff(); //allow ad

		// 8. Genetic distance
		CalcDist(); //allow ad, circle buf

		// 9. AMOVA
		CalcAMOVA(); //forbid ad

		// 10. Population assignment
		CalcAssignment(); //forbid ad

		// 11. Relatedness coefficient
		CalcRelatedness(); //forbid ad, circle buf

		// 12. Kinship coefficient
		CalcKinship(); //forbid ad, circle buf

		// 13. Principal coordinate analysis
		CalcPCOA(); //allow ad

		// 14. Hierarchical clustering
		CalcClustering(); //allow ad

		//TEST SPA
		CalcSPA();//allow ad

		//15. Bayesian clustering
		CalcBayesian(); //forbid ad

		//16. Ploidy Inference
		CalcPloidyInference(); //forbid ad

		UnInitialize();
	}

	/* Create genotype table for non-vcf/bcf input files */
	TARGET void CreateGenoIndexTable(GENO_ITERATOR *iter)
	{
		genotype_offset = new OFFSET[nloc];
		genotype_coffset = 0;
		int64 size = 0;
		for (int64 l = 0; l < nloc; ++l)
		{
			size = CeilLog2((int)locus[l].ngeno);//in bits 
			genotype_offset[l].offset = genotype_coffset;
			genotype_offset[l].size = size;
			genotype_coffset += ceil(size * nind / 8.0) + 3;//in bytes
		}
		genotype_size = genotype_coffset;
		genotype_bucket = new byte[genotype_size];

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
			new(&iter[l]) GENO_ITERATOR(0, l, false);//OK
	}

	/* 0. Loading thread functions */

	/* load from Genepop input format */
	THREAD(LoadGenepop)
	{
		//load genepop
		int64 buflen = 64 * 1024;
		char *buf = new char[buflen];
		FGets(buf, buflen, g_filehandle[0][0]); //ignore firstline
		bool findhead = false;

		//count locus
		int64 loc_pos = FTell(g_filehandle[0][0]); //backup position
		while (!FEof(g_filehandle[0][0]))
		{
			FGets(buf, buflen, g_filehandle[0][0]);
			if (LwrLineCmp("pop", buf) == 0)
			{
				findhead = true;
				break;
			}
			nloc++;
		}

		if (!findhead)
			Exit("\nError: Cannot parse input file.\n");

		if (nloc == 0)
			Exit("\nError: there are no loci within this file.\n");

		int64 ind_pos = FTell(g_filehandle[0][0]); //backup position
		while (!FEof(g_filehandle[0][0]))
		{
			int64 rlen = FGets(buf, buflen, g_filehandle[0][0]);
			if (LwrLineCmp("pop", buf) && rlen > 4ll)
				nind++;
		}

		if (nind == 0)
			Exit("\nError: there are no individuals within this file.\n");

		ainds = new IND*[nind];
		/*   1   */
		ushort **gatab = new ushort*[nloc];
		GENOTYPE **gtab = new GENOTYPE*[nloc];
		GENO_ITERATOR *iter = new GENO_ITERATOR[nloc];//OK

		locus = new LOCUS[nloc];
		nvcf_memory = new MEMORY[g_nthread_val];
		nvcf_gfid = new TABLE<HASH, uint>[nloc];

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			threadid = omp_get_thread_num();
			new(&nvcf_gfid[l]) TABLE<HASH, GENOTYPE*>(true, &nvcf_memory[threadid]);
		}

		for (int64 j = 1; j >= 0; --j)
		{
			FSeek(g_filehandle[0][0], ind_pos, SEEK_SET);
			for (int i = 0; i < nind; )
			{
				FGets(buf, buflen, g_filehandle[0][0]);
				if (LwrLineCmp("pop", buf) == 0) continue;
				if (j != 0) individual_memory->Alloc(ainds[i], 1);
				/*   2   */
				new(ainds[i]) IND(buf, j, i, gtab, gatab, iter);
				i++;
				PROGRESS_VALUE = FTell(g_filehandle[0][0]) + (j ? 0 : g_filelen[0][0]);
			}

			if (j != 0)
			{
				/*   3   */
				FSeek(g_filehandle[0][0], loc_pos, SEEK_SET);
				for (int64 l = 0; l < nloc; ++l)
				{
					FGets(buf, buflen, g_filehandle[0][0]);
					if (LwrLineCmp("pop", buf) == 0) continue;
					new(&locus[l]) LOCUS(locus_memory[threadid], buf, l, nvcf_gfid[l].size, gtab[l], gatab[l]);
				}
				CreateGenoIndexTable(iter);
			}
		}

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
			iter[l].FinishWrite();

		PROGRESS_VALUE = PROGRESS_CEND;
		IndexAlleleLength();

		/*   4   */
		delete[] iter;
		delete[] gatab;
		delete[] gtab;
		delete[] nvcf_gfid;
		delete[] nvcf_memory;
		delete[] buf;
	}

	/* load from Spagedi input format */
	THREAD(LoadSpagedi)
	{
		//load spagedi
		int64 buflen = 64 * 1024;
		char *buf = new char[buflen];

		int64 tpos = 0;

		do
		{
			tpos = FTell(g_filehandle[0][0]);
			FGets(buf, buflen, g_filehandle[0][0]); //ignore firstline
		} while (buf[0] == '/' && buf[1] == '/');

		char *bufp = buf;
		for (int i = 0; i < 5; ++i)
			genotype_digit = ReadInteger(bufp);
		FGets(buf, buflen, g_filehandle[0][0]); //-3
		int extracol = ReadIntegerKeep(buf) + 3;
		genotype_extracol = extracol;

		//count locus
		int64 loc_pos = FTell(g_filehandle[0][0]); //backup position
		FGets(buf, buflen, g_filehandle[0][0]);
		bufp = buf + strlen(buf);
		while (*bufp == '\n' || *bufp == '\r' || *bufp == '\t' || *bufp == ' ')
			*bufp-- = '\0';
		nloc = CountChar(buf, '\t') + 1 - 2 - extracol;

		if (nloc == 0)
			Exit("\nError: there are no loci within this file.\n");

		//count individual
		int64 ind_pos = FTell(g_filehandle[0][0]); //backup position
		while (!FEof(g_filehandle[0][0]))
		{
			FGets(buf, buflen, g_filehandle[0][0]);
			if (LwrLineCmp("end", buf) == 0) break;
			nind++;
		}

		if (nind == 0)
			Exit("\nError: there are no individuals within this file.\n");

		FSeek(g_filehandle[0][0], ind_pos, SEEK_SET);

		ainds = new IND*[nind];
		/*   1   */
		ushort **gatab = new ushort*[nloc];
		GENOTYPE **gtab = new GENOTYPE*[nloc];
		GENO_ITERATOR *iter = new GENO_ITERATOR[nloc];//OK

		locus = new LOCUS[nloc];
		nvcf_memory = new MEMORY[g_nthread_val];
		nvcf_gfid = new TABLE<HASH, uint>[nloc];
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			threadid = omp_get_thread_num();
			new(&nvcf_gfid[l]) TABLE<HASH, GENOTYPE*>(true, &nvcf_memory[threadid]);
		}

		for (int64 j = 1; j >= 0; --j)
		{
			FSeek(g_filehandle[0][0], ind_pos, SEEK_SET);
			for (int i = 0; i < nind; )
			{
				FGets(buf, buflen, g_filehandle[0][0]);
				if (LwrLineCmp("end", buf) == 0) break;
				if (j != 0) individual_memory->Alloc(ainds[i], 1);
				/*   2   */
				new(ainds[i]) IND(buf, j, i, gtab, gatab, iter);
				i++;
				PROGRESS_VALUE = FTell(g_filehandle[0][0]) + (j ? 0 : g_filetotlen);
			}

			if (j != 0)
			{
				/*   3   */
				FSeek(g_filehandle[0][0], loc_pos, SEEK_SET);
				FGets(buf, buflen, g_filehandle[0][0]);
				char *locus_name = StrNextIdx(buf, '\t', 2 + extracol) + 1;
				for (int64 l = 0; l < nloc; ++l)
				{
					new(&locus[l]) LOCUS(locus_memory[threadid], locus_name, l, nvcf_gfid[l].size, gtab[l], gatab[l]);
					locus_name = StrNextIdx(locus_name, '\t', 1) + 1;
				}
				CreateGenoIndexTable(iter);
			}
		}

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
			iter[l].FinishWrite();

		PROGRESS_VALUE = PROGRESS_CEND;
		IndexAlleleLength();

		/*   4   */
		delete[] iter;
		delete[] gatab;
		delete[] gtab;
		delete[] nvcf_gfid;
		delete[] nvcf_memory;
		delete[] buf;
	}

	/* load from Cervus input format */
	THREAD(LoadCervus)
	{
		//load cervus
		int64 buflen = 64 * 1024;
		char *buf = new char[buflen];
		char *bufp = buf;

		//count locus
		int64 loc_pos = FTell(g_filehandle[0][0]); //backup position
		FGets(buf, buflen, g_filehandle[0][0]);

		bufp = buf + strlen(buf);
		while (*bufp == '\n' || *bufp == '\r' || *bufp == ',' || *bufp == ' ')
			*bufp-- = '\0';
		nloc = CountChar(buf, ',') + 1;
		int extracol = (nloc - ((nloc - 1) / 2) * 2) - 1;
		nloc = (nloc - 1) / 2;
		genotype_extracol = extracol;

		if (nloc == 0)
			Exit("\nError: there are no loci within this file.\n");

		//count individual
		int64 ind_pos = FTell(g_filehandle[0][0]); //backup position
		while (!FEof(g_filehandle[0][0]))
			if (FGets(buf, buflen, g_filehandle[0][0]) > 3) nind++;

		if (nind == 0)
			Exit("\nError: there are no individuals within this file.\n");

		FSeek(g_filehandle[0][0], ind_pos, SEEK_SET);

		ainds = new IND*[nind];
		/*   1   */
		ushort **gatab = new ushort*[nloc];
		GENOTYPE **gtab = new GENOTYPE*[nloc];
		GENO_ITERATOR *iter = new GENO_ITERATOR[nloc];//OK

		locus = new LOCUS[nloc];
		nvcf_memory = new MEMORY[g_nthread_val];
		nvcf_gfid = new TABLE<HASH, uint>[nloc];
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			threadid = omp_get_thread_num();
			new(&nvcf_gfid[l]) TABLE<HASH, GENOTYPE*>(true, &nvcf_memory[threadid]);
		}

		for (int64 j = 1; j >= 0; --j)
		{
			FSeek(g_filehandle[0][0], ind_pos, SEEK_SET);
			for (int i = 0; i < nind; )
			{
				if (FGets(buf, buflen, g_filehandle[0][0]) <= 3)  continue;
				if (j != 0) individual_memory->Alloc(ainds[i], 1);
				/*   2   */
				new(ainds[i]) IND(buf, j, i, gtab, gatab, iter);
				i++;
				PROGRESS_VALUE = FTell(g_filehandle[0][0]) + (j ? 0 : g_filetotlen);
			}
			if (j != 0)
			{
				/*   3   */
				FSeek(g_filehandle[0][0], loc_pos, SEEK_SET);
				FGets(buf, buflen, g_filehandle[0][0]);
				char *locus_name = StrNextIdx(buf, ',', extracol + 1) + 1;

				for (int64 l = 0; l < nloc; ++l)
				{
					char *locus_name2 = StrNextIdx(locus_name, ',', 1);
					*(locus_name2 - 1) = '\0';
					new(&locus[l]) LOCUS(locus_memory[threadid], locus_name, l, nvcf_gfid[l].size, gtab[l], gatab[l]);
					*(locus_name2 - 1) = 'A';
					locus_name = StrNextIdx(locus_name, ',', 2) + 1;
				}
				CreateGenoIndexTable(iter);
			}
		}

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
			iter[l].FinishWrite();

		PROGRESS_VALUE = PROGRESS_CEND;
		IndexAlleleLength();

		/*   4   */
		delete[] iter;
		delete[] gatab;
		delete[] gtab;
		delete[] nvcf_gfid;
		delete[] nvcf_memory;
		delete[] buf;
	}

	/* load from Arlequin input format */
	THREAD(LoadArlequin)
	{
		//load structure
		int64 buflen = 64 * 1024;
		char *buf = new char[buflen];
		char *bufp, *bufp2;

		//count locus
		int64 loc_pos = 0; //backup position
		bool header = false;
		for (;;)
		{
			FGets(buf, buflen, g_filehandle[0][0]);
			bufp = buf;
			while (*bufp == '\t' || *bufp == ' ') bufp++;
			if (LwrLineCmp("SampleData", bufp) == 0)
			{
				header = true;
				loc_pos = FTell(g_filehandle[0][0]);
				break;
			}
		}

		FGets(buf, buflen, g_filehandle[0][0]);
		ReplaceChar(buf, '\t', ' ');
		bufp = buf;
		while (*bufp == ' ') bufp++; //name
		bufp = StrNextIdx(bufp, ' ', 1); while (*bufp == ' ') bufp++; //1
		bufp = StrNextIdx(bufp, ' ', 1); while (*bufp == ' ') bufp++; //genotype
		bufp2 = bufp + strlen(bufp) - 1;
		while (*bufp2 == ' ' || *bufp2 == '\n' || *bufp2 == '\r') *bufp2-- = '\0';
		nloc = CountChar(bufp, ' ') + 1;
		if (nloc == 0) Exit("\nError: there are no loci within this file.\n");

		int64 ind_pos = loc_pos;
		FSeek(g_filehandle[0][0], ind_pos, SEEK_SET);
		while (!FEof(g_filehandle[0][0]))
		{
			//read name
			if (FGets(buf, buflen, g_filehandle[0][0]) == 0) break;
			bufp = buf;
			ReplaceChar(buf, '\t', ' ');
			while (*bufp == ' ') bufp++;
			if (*bufp == '}') for (;;)
			{
				if (FGets(buf, buflen, g_filehandle[0][0]) == 0) break;
				ReplaceChar(buf, '\t', ' ');
				bufp = buf;
				while (*bufp == ' ') bufp++;
				if (LwrLineCmp("SampleData", bufp) == 0) break;
			}
			else
			{
				if (FGets(bufp, buflen, g_filehandle[0][0], buf) == 0) break;
				nind++;
			}
		}

		if (nind == 0)
			Exit("\nError: there are no individuals within this file.\n");

		ainds = new IND*[nind];
		/*   1   */
		ushort **gatab = new ushort*[nloc];
		GENOTYPE **gtab = new GENOTYPE*[nloc];
		GENO_ITERATOR *iter = new GENO_ITERATOR[nloc];//OK

		locus = new LOCUS[nloc];
		nvcf_memory = new MEMORY[g_nthread_val];
		nvcf_gfid = new TABLE<HASH, uint>[nloc];
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			threadid = omp_get_thread_num();
			new(&nvcf_gfid[l]) TABLE<HASH, GENOTYPE*>(true, &nvcf_memory[threadid]);
		}

		for (int64 j = 1; j >= 0; --j)
		{
			FSeek(g_filehandle[0][0], ind_pos, SEEK_SET);
			for (int i = 0; i < nind; )
			{
				//read name
				if (FGets(buf, buflen, g_filehandle[0][0]) == 0) break;
				bufp = buf;
				ReplaceChar(buf, '\t', ' ');
				while (*bufp == ' ') bufp++;
				if (*bufp == '}') for (;;)
				{
					if (FGets(buf, buflen, g_filehandle[0][0]) == 0) break;
					ReplaceChar(buf, '\t', ' ');
					bufp = buf;
					while (*bufp == ' ') bufp++;
					if (LwrLineCmp("SampleData", bufp) == 0) break;
				}
				else
				{
					bufp = buf + strlen(buf);
					if (FGets(bufp, buflen, g_filehandle[0][0], buf) == 0) break;
					if (j != 0) individual_memory->Alloc(ainds[i], 1);
					/*   2   */
					new(ainds[i]) IND(buf, j, i, gtab, gatab, iter);
					i++;
					PROGRESS_VALUE = FTell(g_filehandle[0][0]) + (j ? 0 : g_filetotlen);
				}
			}
			if (j != 0)
			{
				/*   3   */
				FSeek(g_filehandle[0][0], ind_pos, SEEK_SET);
				for (int64 l = 0; l < nloc; ++l)
					new(&locus[l]) LOCUS(locus_memory[threadid], (char*)NULL, l, nvcf_gfid[l].size, gtab[l], gatab[l]);
				CreateGenoIndexTable(iter);
			}
		}

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
			iter[l].FinishWrite();

		PROGRESS_VALUE = PROGRESS_CEND;
		IndexAlleleLength();

		/*   4   */
		delete[] iter;
		delete[] gatab;
		delete[] gtab;
		delete[] nvcf_gfid;
		delete[] nvcf_memory;
		delete[] buf;
	}

	/* load from Structure input format */
	THREAD(LoadStructure)
	{
		//load structure
		int64 buflen = 64 * 1024;
		char *buf = new char[buflen];
		char *bufp = buf;

		//count locus
		int64 loc_pos = 0; //backup position
		FGets(buf, buflen, g_filehandle[0][0]);
		int buflen2 = (int)strlen(buf);
		while (buf[buflen2 - 1] == ' ' || buf[buflen2 - 1] == '\r' || buf[buflen2 - 1] == '\n')
			buf[buflen2-- - 1] = '\0';
		ReplaceChar(buf, '\t', ' ');
		nloc = CountChar(buf, ' ') + 1;

		FGets(buf, buflen, g_filehandle[0][0]);
		ReplaceChar(buf, '\t', ' ');
		int64 nloc2 = CountChar(buf, ' ') + 1;

		bool headerrow = false;
		if (nloc != nloc2)
			headerrow = true;
		else
			nloc = nloc2 - 1 - g_extracol_val;

		if (nloc == 0) Exit("\nError: there are no loci within this file.\n");

		FSeek(g_filehandle[0][0], 0, SEEK_SET);

		//jump 1
		FGets(buf, buflen, g_filehandle[0][0]);
		ReplaceChar(buf, '\t', ' ');
		bufp = buf;

		//count individual
		int64 ind_pos = FTell(g_filehandle[0][0]); //backup position

		while (!FEof(g_filehandle[0][0]))
		{
			//read name
			int64 tpos1 = FTell(g_filehandle[0][0]);
			if (FGetsSmall(buf, 1024, g_filehandle[0][0]) == 0) break;
			ReplaceChar(buf, '\t', ' ');
			bufp = buf;
			while (*bufp != ' ') bufp++; //skip name
			*bufp++ = '\0'; *bufp = '\0';
			FSeek(g_filehandle[0][0], tpos1, SEEK_SET);

			//read lines
			while (!FEof(g_filehandle[0][0]))
			{
				//read one more section
				int64 tpos2 = FTell(g_filehandle[0][0]);
				FGetsSmall(bufp, 1024, g_filehandle[0][0]);
				ReplaceChar(bufp, '\t', ' ');
				FSeek(g_filehandle[0][0], tpos2, SEEK_SET);
				if (LineCmpABterm(buf, bufp, ' ')) break; //not the same individual, break

				//the same individual, add length
				FGets(bufp, buflen, g_filehandle[0][0], buf);
				ReplaceChar(bufp, '\t', ' ');
			}

			//end line
			nind++;
		}

		if (nind == 0)
			Exit("\nError: there are no individuals within this file.\n");

		ainds = new IND*[nind];
		/*   1   */
		ushort **gatab = new ushort*[nloc];
		GENOTYPE **gtab = new GENOTYPE*[nloc];
		GENO_ITERATOR *iter = new GENO_ITERATOR[nloc];//OK

		locus = new LOCUS[nloc];
		nvcf_memory = new MEMORY[g_nthread_val];
		nvcf_gfid = new TABLE<HASH, uint>[nloc];
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			threadid = omp_get_thread_num();
			new(&nvcf_gfid[l]) TABLE<HASH, GENOTYPE*>(true, &nvcf_memory[threadid]);
		}

		for (int64 j = 1; j >= 0; --j)
		{
			FSeek(g_filehandle[0][0], ind_pos, SEEK_SET);
			for (int i = 0; i < nind; )
			{
				//read name
				int64 tpos1 = FTell(g_filehandle[0][0]);
				if (FGetsSmall(buf, 1024, g_filehandle[0][0]) == 0) break;
				ReplaceChar(buf, '\t', ' ');

				bufp = buf;
				while (*bufp != ' ') bufp++; //skip name
				*bufp++ = '\0'; *bufp = '\0';
				FSeek(g_filehandle[0][0], tpos1, SEEK_SET);

				//read lines
				char *bufp2 = bufp;
				while (!FEof(g_filehandle[0][0]))
				{
					//read one more section
					int64 tpos2 = FTell(g_filehandle[0][0]);
					FGetsSmall(bufp2, 1024, g_filehandle[0][0]);
					ReplaceChar(bufp2, '\t', ' ');
					FSeek(g_filehandle[0][0], tpos2, SEEK_SET);
					if (LineCmpABterm(buf, bufp2, ' '))
					{
						*bufp2 = '\0';
						break; //not the same individual, break
					}

					//the same individual, add length
					while (!FEof(g_filehandle[0][0]))
					{
						FGets(bufp2, buflen, g_filehandle[0][0], buf);
						ReplaceChar(bufp2, '\t', ' ');
						while (*bufp2) bufp2++;
						if (*(bufp2 - 1) == '\n')
							break;
					}
				}

				//end line
				if (j != 0) individual_memory->Alloc(ainds[i], 1);
				/*   2   */
				new(ainds[i]) IND(bufp, j, i, gtab, gatab, iter);
				i++;
				PROGRESS_VALUE = FTell(g_filehandle[0][0]) + (j ? 0 : g_filetotlen);
			}
			if (j != 0)
			{
				/*   3   */
				FSeek(g_filehandle[0][0], loc_pos, SEEK_SET);
				FGets(buf, buflen, g_filehandle[0][0]);
				ReplaceChar(buf, '\t', ' ');
				char *locus_name = buf;

				if (headerrow) for (int64 l = 0; l < nloc; ++l)
				{
					char *locus_name2 = StrNextIdx(locus_name, ' ', 1);
					if (locus_name2) *locus_name2 = '\0';
					new(&locus[l]) LOCUS(locus_memory[threadid], locus_name, l, nvcf_gfid[l].size, gtab[l], gatab[l]);
					if (locus_name2) *locus_name2 = ' ';
					locus_name = locus_name2 + 1;
				}
				else for (int64 l = 0; l < nloc; ++l)
					new(&locus[l]) LOCUS(locus_memory[threadid], (char*)NULL, l, nvcf_gfid[l].size, gtab[l], gatab[l]);
				CreateGenoIndexTable(iter);
			}
		}

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
			iter[l].FinishWrite();

		PROGRESS_VALUE = PROGRESS_CEND;
		IndexAlleleLength();

		/*   4   */
		delete[] iter;
		delete[] gatab;
		delete[] gtab;
		delete[] nvcf_gfid;
		delete[] nvcf_memory;
		delete[] buf;
	}

	/* load from PolyGene input format */
	THREAD(LoadPolyGene)
	{
		int64 buflen = 64 * 1024;
		char *buf = new char[64 * 1024];

		int64 tlen = FGets(buf, buflen, g_filehandle[0][0]); //ignore firstline
		while (buf[tlen] == '\t') buf[tlen--] = '\0';

		nloc = CountChar(buf, '\t') - 2;  //ind pop ploidy
		if (nloc == '\0')
			Exit("\nError: there are no loci within this file.\n");

		//Read locus
		char *locus_name = StrNextIdx(buf, '\t', 3) + 1;
		char *locus_name_b = new char[strlen(locus_name) + 1];
		strcpy(locus_name_b, locus_name);
		locus_name = locus_name_b;

		//count individual
		int64 ind_pos = FTell(g_filehandle[0][0]); //backup position
		while (!FEof(g_filehandle[0][0]))
			if (FGets(buf, buflen, g_filehandle[0][0]) > 5)
				nind++;

		if (nind == 0)
			Exit("\nError: there are no individuals within this file.\n");

		FSeek(g_filehandle[0][0], ind_pos, SEEK_SET);
		ainds = new IND*[nind];
		/*   1   */
		ushort **gatab = new ushort*[nloc];
		GENOTYPE **gtab = new GENOTYPE*[nloc];
		GENO_ITERATOR *iter = new GENO_ITERATOR[nloc];//OK

		locus = new LOCUS[nloc];
		nvcf_memory = new MEMORY[g_nthread_val];
		nvcf_gfid = new TABLE<HASH, uint>[nloc];

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			threadid = omp_get_thread_num();
			new(&nvcf_gfid[l]) TABLE<HASH, GENOTYPE*>(true, &nvcf_memory[threadid]);
		}

		for (int64 j = 1; j >= 0; --j)
		{
			ind_pos = FTell(g_filehandle[0][0]);
			for (int i = 0; i < nind; )
			{
				FGets(buf, buflen, g_filehandle[0][0]);
				if (j != 0) individual_memory->Alloc(ainds[i], 1);
				/*   2   */
				new(ainds[i]) IND(buf, j, i, gtab, gatab, iter);
				i++;
				PROGRESS_VALUE = FTell(g_filehandle[0][0]) + (j ? 0 : g_filetotlen);
			}
			FSeek(g_filehandle[0][0], ind_pos, SEEK_SET);
			if (j != 0)
			{
				/*   3   */
				for (int64 l = 0; l < nloc; ++l)
				{
					new(&locus[l]) LOCUS(locus_memory[threadid], locus_name, l, nvcf_gfid[l].size, gtab[l], gatab[l]);
					locus_name = StrNextIdx(locus_name, '\t', 1) + 1;
				}
				CreateGenoIndexTable(iter);
			}
		}

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
			iter[l].FinishWrite();

		PROGRESS_VALUE = PROGRESS_CEND;
		IndexAlleleLength();

		/*   4   */
		delete[] locus_name_b;
		delete[] iter;
		delete[] gatab;
		delete[] gtab;
		delete[] nvcf_gfid;
		delete[] nvcf_memory;
		delete[] buf;
	}

	/* load from PolyRelatedness input format */
	THREAD(LoadPolyRelatedness)
	{
		//load polyrelatedness
		int64 buflen = 64 * 1024;
		char *buf = new char[64 * 1024];

		int64 tpos = 0;
		do
		{
			tpos = FTell(g_filehandle[0][0]);
			FGets(buf, buflen, g_filehandle[0][0]); //ignore firstline
		} while (buf[0] == '/' && buf[1] == '/');

		char *bufp = buf;
		genotype_digit = ReadInteger(bufp);
		genotype_missing = ReadInteger(bufp);
		genotype_missing = ReadInteger(bufp);
		genotype_ambiguous = ReadInteger(bufp);

		//count locus
		bool findhead = false;
		while (!FEof(g_filehandle[0][0]))
		{
			FGets(buf, buflen, g_filehandle[0][0]);
			if (LwrLineCmp("//genotype", buf) == 0)
			{
				findhead = true;
				break;
			}
		}
		int64 loc_pos = FTell(g_filehandle[0][0]); //backup position
		FGets(buf, buflen, g_filehandle[0][0]);

		bufp = buf + strlen(buf);
		while (*bufp == '\n' || *bufp == '\r' || *bufp == '\t' || *bufp == ' ')
			*bufp-- = '\0';
		nloc = CountChar(buf, '\t') + 1 - 2;

		if (nloc == '\0')
			Exit("\nError: there are no loci within this file.\n");

		// alloc memory
		FSeek(g_filehandle[0][0], loc_pos, SEEK_SET);

		//jump 2
		FGets(buf, buflen, g_filehandle[0][0]);
		char *locus_name = StrNextIdx(buf, '\t', 2) + 1, *locus_name_b = NULL;
		if (locus_name != (char*)1)
		{
			locus_name_b = new char[strlen(locus_name) + 1];
			strcpy(locus_name_b, locus_name);
			locus_name = locus_name_b;
		}

		//count individual
		int64 ind_pos = FTell(g_filehandle[0][0]); //backup position
		while (!FEof(g_filehandle[0][0]))
		{
			FGets(buf, buflen, g_filehandle[0][0]);
			if (LwrLineCmp("//end of file", buf) == 0) break;
			nind++;
		}

		if (nind == 0)
			Exit("\nError: there are no individuals within this file.\n");

		ainds = new IND*[nind];
		/*   1   */
		ushort **gatab = new ushort*[nloc];
		GENOTYPE **gtab = new GENOTYPE*[nloc];
		GENO_ITERATOR *iter = new GENO_ITERATOR[nloc];//OK

		locus = new LOCUS[nloc];
		nvcf_memory = new MEMORY[g_nthread_val];
		nvcf_gfid = new TABLE<HASH, uint>[nloc];

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			threadid = omp_get_thread_num();
			new(&nvcf_gfid[l]) TABLE<HASH, GENOTYPE*>(true, &nvcf_memory[threadid]);
		}

		for (int64 j = 1; j >= 0; --j)
		{
			FSeek(g_filehandle[0][0], ind_pos, SEEK_SET);
			for (int i = 0; i < nind; )
			{
				FGets(buf, buflen, g_filehandle[0][0]);
				if (LwrLineCmp("//end of file", buf) == 0) break;
				if (j != 0) individual_memory->Alloc(ainds[i], 1);
				/*   2   */
				new(ainds[i]) IND(buf, j, i, gtab, gatab, iter);
				i++;
				PROGRESS_VALUE = FTell(g_filehandle[0][0]) + (j ? 0 : g_filetotlen);
			}

			if (j != 0)
			{
				/*   3   */
				FSeek(g_filehandle[0][0], loc_pos, SEEK_SET);
				for (int64 l = 0; l < nloc; ++l)
				{
					new(&locus[l]) LOCUS(locus_memory[threadid], locus_name == (char*)1 ? NULL : locus_name, l, nvcf_gfid[l].size, gtab[l], gatab[l]);
					if (locus_name != (char*)1) locus_name = StrNextIdx(locus_name, '\t', 1) + 1;
				}
				CreateGenoIndexTable(iter);
			}
		}

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
			iter[l].FinishWrite();

		PROGRESS_VALUE = PROGRESS_CEND;
		IndexAlleleLength();

		/*   4   */
		if (locus_name_b) delete[] locus_name_b;
		delete[] iter;
		delete[] gatab;
		delete[] gtab;
		delete[] nvcf_gfid;
		delete[] nvcf_memory;
		delete[] buf;
	}

	/* Indexing integer alleles to zero based for non-vcf input, with allele identifier being the size */
	TARGET void IndexAlleleLength()
	{
		//paralleled

		TABLE<ushort, ushort> *aftab = new TABLE<ushort, ushort>[g_nthread_val];
		for (int i = 0; i < g_nthread_val; ++i)
			new(&aftab[i]) TABLE<ushort,ushort>(true, NULL);

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			threadid = omp_get_thread_num();
			TABLE<ushort, ushort> &atab = aftab[threadid]; 
			atab.Clear();

			GENOTYPE *gtab = locus[l].GetGtab();
			int ngenotype = locus[l].ngeno;

			for (int gi = 0; gi < ngenotype; ++gi)
			{
				GENOTYPE &gt = gtab[gi];
				if (gt.Nalleles() == 0) continue;
				int v = gt.Ploidy();
				ushort *als = gt.GetAlleleArray();
				for (int vi = 0; vi < v; ++vi)
					atab.PushIndex(als[vi]);
			}

			atab.GetLength(locus[l]._alen, locus_memory[threadid]);

			locus[l].k = (ushort)atab.size;
			locus[l].flag_alen = 1;
			if (locus[l].k < 2) locus[l].flag_pass = false;

			for (int gi = 0; gi < ngenotype; ++gi)
			{
				GENOTYPE &gt = gtab[gi];
				if (gt.Nalleles() == 0) continue;
				int v = gt.Ploidy();
				ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
				ushort *als = gt.GetAlleleArray();

				//mapping alleles
				for (int vi = 0; vi < v; ++vi)
					alleles[vi] = atab[als[vi]];

				//non-vcf genotypes are unphased, can be sort
				Sort(alleles, v);

				ushort *nullgatab = NULL;

				//sorted, reconstruct alleles, hash is invalid
				new(&gt) GENOTYPE(nullgatab, alleles, v);
			}

		}

		delete[] aftab;
	}

	/* Check input files rows and columns and count the number of variants */
	THREAD(GetBCFLines)
	{
		char *bufj = new char[64 * 1024];
		char *loc_name = new char[64 * 1024];
		VLA_NEW(filep, int64, g_input_col);
		VLA_NEW(jloc, int64, g_input_col);
		VLA_NEW(offsetj, int64, g_input_col);
		int64 totlen = 0;

		for (int i = 0; i < g_input_row; ++i)
		{
			for (int j = 0; j < g_input_col; ++j)
			{
				filep[j] = FTell(g_filehandle[i][j]);
				offsetj[j] = FOffset(g_filehandle[i][j]);
				PROGRESS_VALUE += offsetj[j];
				jloc[j] = 0;
			}

			if (g_input_col == 1) for (int64 nextpos = FTell(g_filehandle[i][0]);;)
			{
				nextpos += 8 + FGetUint(g_filehandle[i][0]) + FGetUint(g_filehandle[i][0]);
				if (FEof(g_filehandle[i][0])) break;
				FSeek(g_filehandle[i][0], nextpos, SEEK_SET);
				PROGRESS_VALUE = offsetj[0] + FOffset(g_filehandle[i][0]);
				jloc[0]++;
			}
			else for (int64 nextpos = FTell(g_filehandle[i][0]);;) for (int j = 0; j < g_input_col; ++j)
			{
				//read variant info
				nextpos += 8 + FGetUint(g_filehandle[i][0]) + FGetUint(g_filehandle[i][0]);
				if (FEof(g_filehandle[i][j]))
				{
					if (j == g_input_col - 1) goto end;
					continue;
				}
				int64 noffset = FOffset(g_filehandle[i][j]);
				PROGRESS_VALUE += noffset - offsetj[j];
				offsetj[j] = noffset;

				char *bufjj = bufj;

				//CHROM
				int chrom = FGetUint(g_filehandle[i][j]);
				AppendString(bufjj, g_bcfheader[i][j]->contig_name[chrom]); *bufjj++ = '\t';

				//POS
				int pos = FGetUint(g_filehandle[i][j]);
				AppendInt(bufjj, pos); *bufjj++ = '\t';

				FSeek(g_filehandle[i][j], 10, SEEK_CUR);
				ushort nallele = FGetUshort(g_filehandle[i][j]);
				FSeek(g_filehandle[i][j], 4, SEEK_CUR);

				//ID
				int typelen = FGetByte(g_filehandle[i][j]) >> 4;
				if (typelen == 0xF) typelen = FGetTypedInt(g_filehandle[i][j]);
				FGet(g_filehandle[i][j], bufjj, typelen); bufjj += typelen;  *bufjj++ = '\t';

				//REF ALT
				for (int a = 0; a < nallele; ++a)
				{
					typelen = FGetByte(g_filehandle[i][j]) >> 4;
					typelen = (typelen == 0xF) ? FGetTypedInt(g_filehandle[i][j]) : typelen;
					FGet(g_filehandle[i][j], bufjj, typelen); bufjj += typelen;  *bufjj++ = '\t';
				}
				*bufjj++ = '\0';
				jloc[j]++;

				if (j != 0)
				{
					if (LineCmp(loc_name, bufj))
						Exit("\nError: variant information in %s and %s mismatch, at lines %d:\n%s\nvs\n%s\n",
							g_filename[i][0], g_filename[i][j], jloc[j], loc_name, bufj);
				}
				else
					memcpy(loc_name, bufj, bufjj - bufj + 1);

				FSeek(g_filehandle[i][j], nextpos, SEEK_SET);
			}

		end:
			for (int j = 0; j < g_input_col; ++j)
			{
				totlen += g_filelen[i][j];
				if (!FEof(g_filehandle[i][j]))
					Exit("\nError: number of variants in %s mismatch that in %s.\n", g_filename[i][j], g_filename[i][0]);
				if (jloc[j] != jloc[0])
					Exit("\nError: number of variants in %s mismatch that in %s.\n", g_filename[i][j], g_filename[i][0]);
			}
			nloc += jloc[0];
			PROGRESS_VALUE = totlen;

			for (int j = 0; j < g_input_col; ++j)
				FSeek(g_filehandle[i][j], filep[j], SEEK_SET);
		}

		delete[] bufj;
		delete[] loc_name;
		VLA_DELETE(filep);
		VLA_DELETE(jloc);
		VLA_DELETE(offsetj);
	}

	/* Process lines from memory */
	THREAD(LoadBCF)
	{
		TABLE<HASH, uint> gfid(false, NULL);

		for (int64 ii = 0; ii < nloc; ii++)
		{
			while (ii >= progress1)
				SLEEP(SLEEP_TIME_TINY);
			int64 avail_val = ii * 4 + 2, calc_val = avail_val + 1;

			if (state_lock[ii % NBUF].compare_exchange_strong(avail_val, calc_val)) {

			gfid.Clear();

			char *str = load_buf[ii % NBUF];
			int64 byteflag = *(int64*)(str + sizeof(int64));
			str += sizeof(int64) * 2;
			char *dst = StrNextIdx(str, '\t', 9) + 1;//GT

			byte type = *dst++;
			int ploidy = type >> 4;
			if (ploidy == 0xF) ploidy = ReadTypedInt(dst);
			if (ploidy >= 10) Exit("\nError: do not support a ploidy level greater than 10.\n");
			type &= 0xF;
			int asize = 0;
			switch (type)
			{
			case 1: asize = 1; break;
			case 2: asize = 2; break;
			case 3: asize = 4; Exit("\nError: do not support 4 byte allele identifier in BCF format, the maximum number of alleles allowed is 65535.\n"); break; break;
			default: Exit("\nError: GT format is should encode as integers.\n"); break;
			}

			bool usedploidy[N_MAX_PLOIDY + 1] = { 0 };
			for (int i = 0; i < nind; ++i)
			{
				char *genostr = dst;
				
				for (int ai = 0; ai < ploidy; ++ai)
					switch (type)
					{
					case 1: 
						if (*(byte*)genostr == (byte)0x81 || *(byte*)genostr == (byte)0x80 || *(byte*)genostr >> 1 == 0) 
							break;
					case 2: 
						if (*(ushort*)genostr == (ushort)0x8001 || *(ushort*)genostr == (ushort)0x8000 || *(ushort*)genostr >> 1 == 0) 
							break;
					case 3: 
						if (*(uint*)genostr == (uint)0x80000001 || *(uint*)genostr == (uint)0x80000000 || *(uint*)genostr >> 1 == 0) 
							break;
					}

				//add missing at ploidy v
				if (!usedploidy[ploidy])
				{
					int tid = gfid.size;
					gfid[missing_hash[ploidy]] = tid;
					usedploidy[ploidy] = true;
				}

				//add gfid if not missing
				if (genostr[0] != '.')
				{
					ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
					bool phase = true;
					int indploidy = 0;
					ReadBCFGenoString(genostr, alleles, phase, indploidy, asize, ploidy, ii, NULL);

					if (haplotype && !phase)
						locus[ii].flag_pass = false;  //Warning

					if (!haplotype) Sort(alleles, indploidy); //phaseed, do not sort alleles

					HASH ha = HashGenotype(alleles, indploidy);
					if (!gfid.ContainsKey(ha))
					{
						//add to hash and alleles array size
						int tid = gfid.size;
						gfid[ha] = tid;
						locus[ii].gasize += indploidy + GetNalleles(alleles, indploidy);
					}
				}

				dst = genostr + ploidy * asize;
			}

			uint64 gtsize = (uint64)CeilLog2(gfid.size);
			Lock(glock2);

			genotype_offset[ii].offset = genotype_coffset;
			genotype_offset[ii].size = gtsize;

			genotype_coffset += (uint64)ceil(gtsize * nind / 8.0) + 3;
			if (genotype_coffset > genotype_size)
			{
				uint64 nsize = ((uint64)(genotype_coffset * nloc / (ii + 1) * 1.1) >> 16) << 16;
				VAllocGenotype(nsize);
				genotype_size += nsize;
			}
			UnLock(glock2);

			////////////////////////////////////////////////////////////////////////////////////

			//alloc in one piece
			GENOTYPE *gtab;  ushort *gatab;
			if (ad != 2)
				new(&locus[ii]) LOCUS(locus_memory[threadid], str, (uint64)byteflag, gfid.size, gtab, gatab);//LoadBCF
			else
				str = StrNextIdx(str, '\t', 9) + 1;

			char *gt = str + 1, *buf1 = str + ploidy * asize * nind + 1, *gq = NULL, *dp = NULL, *adval = NULL;
			int gqlen = 0, dplen = 0, adlen = 0;
			if (locus[ii].gqid != 0xFFFF)
			{
				gq = buf1;
				switch (*gq++ & 0xF)
				{
				case 1: gqlen = 1; break;
				case 2: gqlen = 2; break;
				case 3: gqlen = 4; break;
				}
				buf1 += gqlen * nind + 1;
			}
			if (locus[ii].dpid != 0xFFFF)
			{
				dp = buf1;
				switch (*dp++ & 0xF)
				{
				case 1: dplen = 1; break;
				case 2: dplen = 2; break;
				case 3: dplen = 4; break;
				}
				buf1 += dplen * nind + 1;
			}
			if (locus[ii].adid != 0xFFFF)
			{
				adval = buf1;
				switch (*adval++ & 0xF)
				{
				case 1: adlen = 1; break;
				case 2: adlen = 2; break;
				case 3: adlen = 4; break;
				}
				buf1 += adlen * nind + 1;
			}

			uint *depth = ploidyinfer ? new uint[locus[ii].k * nind] : NULL, *depthbak = depth;

			GENO_ITERATOR iter(0u, ii, false);//OK

			for (int j = 0; j < nind; ++j)
				ainds[j]->AddBCFGenotype(ii, gt, gq, dp, adval, ploidy, asize, gqlen, dplen, adlen, depth, gfid, gtab, gatab, iter);

			iter.FinishWrite();

			if (ploidyinfer)
			{
				//add allele depth
				uint64 adsize = (uint64)CeilLog2((int)locus[ii].maxdepth + 1);
				Lock(glock2);

				alleledepth_offset[ii].offset = alleledepth_coffset;
				alleledepth_offset[ii].size = adsize;

				alleledepth_coffset += ceil(adsize * locus[ii].k * nind / 8.0) + 3;
				if (alleledepth_coffset > alleledepth_size)
				{
					uint64 nsize = ((uint64)(alleledepth_coffset * nloc / (ii + 1) * 1.1) >> 16) << 16;
					VAllocAlleleDepth(nsize);
					alleledepth_size += nsize;
				}
				UnLock(glock2);
				IND::SetAlleleDepth(ii, depthbak, nind * locus[ii].k, 0u);
				delete[] depthbak;
			}

			PROGRESS_VALUE++;

			state_lock[ii % NBUF] = (ii + NBUF) * 4; }
		}
	}

	/* Read lines from BCF file */
	THREAD(LoadBCFGuard)
	{
		int64 &ii = progress1 = 0;
		VLA_NEW(next_line, int64, g_input_col);
		VLA_NEW(backup_offset, int64, g_input_row * g_input_col);
		int buf_len = LINE_BUFFER;
		int bufgt_len = 64 * 1024;
		int bufdp_len = 64 * 1024;
		int bufgq_len = 64 * 1024;
		int bufad_len = 64 * 1024;
		char *buf = new char[buf_len];
		char *bufgt = new char[bufgt_len];
		char *bufgq = new char[bufgq_len];
		char *bufdp = new char[bufdp_len];
		char *bufad = new char[bufad_len];

		for (int i = 0; i < g_input_row; ++i)
		{
			for (int j = 0; j < g_input_col; ++j)
				next_line[j] = backup_offset[i * g_input_col + j] = FTell(g_filehandle[i][j]);

			//read a whole text line across bcfs in the same row
			for (;;)
			{
				while (state_lock[ii % NBUF] != ii * 4)
					SLEEP(SLEEP_TIME_TINY);
				state_lock[ii % NBUF]++;

				double qual = 1e300;
				bool pass = true;
				char *bufj = buf, *bufjgt = bufgt, *bufjgq = bufgq, *bufjdp = bufdp, *bufjad = bufad;
				uint64 byteflag = 0;
				bool breakflag = false;
				for (int j = 0; j < g_input_col; ++j)
				{
					int infolen = FGetUint(g_filehandle[i][j]);
					int64 fmtpos = next_line[j] + 8 + infolen;
					next_line[j] = fmtpos + FGetUint(g_filehandle[i][j]);

					if (FEof(g_filehandle[i][j]))
					{
						breakflag = true;
						break;
					}

					//extend buf
					if (bufj - buf + infolen + 4096 > buf_len)
					{
						int nlen = bufgt_len;
						while (bufj - buf + infolen + 4096 > nlen) nlen <<= 1;
						char *tbuf = new char[nlen];
						SetVal(tbuf, buf, bufj - buf);
						bufj = tbuf + (bufj - buf);
						delete[] buf; buf = tbuf;
						buf_len = nlen;
					}

					//CHROM
					int chromid = FGetUint(g_filehandle[i][j]);
					if (j == 0) { AppendString(bufj, g_bcfheader[i][j]->contig_name[chromid]); *bufj++ = '\t'; }

					//POS
					int pos = FGetUint(g_filehandle[i][j]);
					if (j == 0) { AppendInt(bufj, pos); *bufj++ = '\t'; }

					//RLEN
					int rlen = FGetUint(g_filehandle[i][j]);

					//QUAL
					qual = Min((double)FGetFloat(g_filehandle[i][j]), qual);

					//#INFO
					int ninfo = (int)FGetUshort(g_filehandle[i][j]);

					//#ALLELE
					int nallele = (int)FGetUshort(g_filehandle[i][j]);

					//#SAMPLE
					int nsamp = (int)FGetUshort(g_filehandle[i][j]) | ((int)FGetByte(g_filehandle[i][j]) << 16);

					//#FORMAT
					int nformat = (int)FGetByte(g_filehandle[i][j]);

					//ID
					int typelen = (int)FGetByte(g_filehandle[i][j]) >> 4;

					if (typelen == 0xF)
						typelen = FGetTypedInt(g_filehandle[i][j]);
					if (j == 0)
					{
						FGet(g_filehandle[i][j], bufj, typelen);
						bufj += typelen;
						*bufj++ = '\t';
					}
					else
						FSeek(g_filehandle[i][j], typelen, SEEK_CUR);

					//REF ALT
					for (int a = 0; a < nallele; ++a)
					{
						typelen = FGetByte(g_filehandle[i][j]) >> 4;
						if (typelen == 0xF) typelen = FGetTypedInt(g_filehandle[i][j]);
						if (j == 0)
						{
							FGet(g_filehandle[i][j], bufj, typelen);
							bufj += typelen;
							*bufj++ = a == 0 ? '\t' : ',';
						}
						else
							FSeek(g_filehandle[i][j], typelen, SEEK_CUR);
					}
					if (j == 0) *(bufj - 1) = '\t';


					if (f_qual_b)
					{
						//qual filter
						if (f_qual_b && IsNormal(qual) && qual < f_qual_min && qual > f_qual_max)
							byteflag |= 0x4;
					}

					//FILTER
					int filtertype = (int)FGetByte(g_filehandle[i][j]);
					int numfilter = filtertype >> 4; filtertype &= 0xF;

					if (f_original_b && f_original_val == 1)
					{
						for (int k = 0; k < numfilter; ++k)
						{
							//original filter
							int filterid = 0;
							switch (filtertype)
							{
							case 1: filterid = (int)FGetByte(g_filehandle[i][j]); break;
							case 2: filterid = (int)FGetUshort(g_filehandle[i][j]); break;
							case 3: filterid = (int)FGetUint(g_filehandle[i][j]); break;
							}

							if (filterid != g_bcfheader[i][j]->filter_passidx)
							{
								byteflag |= 0x2;
								pass = false;
							}
						}
					}

					//SKIP INFO
					FSeek(g_filehandle[i][j], fmtpos, SEEK_SET);
					int gtid = g_bcfheader[i][j]->format_gtid, gqid = g_bcfheader[i][j]->format_gqid, dpid = g_bcfheader[i][j]->format_dpid, adid = g_bcfheader[i][j]->format_adid;

					/////////////////////////////////////////////////////////////

					//FORMAT
					for (int a = 0; a < nformat; ++a)
					{
						int fkey = (int)FGetTypedInt(g_filehandle[i][j]);
						int ftype = (int)FGetByte(g_filehandle[i][j]);
						int flen = (ftype >> 4 == 0xF) ? FGetTypedInt(g_filehandle[i][j]) : (ftype >> 4);

						switch (ftype & 0xF)
						{
						case 1: break;
						case 2: flen <<= 1; break;
						case 3: flen <<= 2; break;
						case 5: flen <<= 2; break;
						case 7: break;
						}
						int datalen = g_bcfheader[i][j]->nsample * flen;

						if (fkey == gtid)
						{
							if (bufjgt - bufgt + datalen + 4096 > bufgt_len)
							{
								int nlen = bufgt_len;
								while (bufjgt - bufgt + datalen + 4096 > nlen) nlen <<= 1;
								char *tbuf = new char[nlen];
								SetVal(tbuf, bufgt, bufjgt - bufgt);
								bufjgt = tbuf + (bufjgt - bufgt);
								delete[] bufgt; bufgt = tbuf;
								bufgt_len = nlen;
							}
							if (j == 0)
							{
								*(byte*)bufjgt++ = (byte)ftype;
								if ((ftype >> 4) == 0xF) AppendTypedInt(bufjgt, flen);
							}
							FGet(g_filehandle[i][j], bufjgt, datalen);
							bufjgt += datalen;
						}
						else if (fkey == gqid)
						{
							if (bufjgq - bufgq + datalen + 4096 > bufgq_len)
							{
								int nlen = bufgq_len;
								while (bufjgq - bufgq + datalen + 4096 > nlen) nlen <<= 1;
								char *tbuf = new char[nlen];
								SetVal(tbuf, bufgq, bufjgq - bufgq);
								bufjgq = tbuf + (bufjgq - bufgq);
								delete[] bufgq; bufgq = tbuf;
								bufgq_len = nlen;
							}
							if (j == 0)
							{
								*(byte*)bufjgq++ = (byte)ftype;
								if ((ftype >> 4) == 0xF) AppendTypedInt(bufjgq, flen);
							}
							FGet(g_filehandle[i][j], bufjgq, datalen);
							bufjgq += datalen;
						}
						else if (fkey == dpid)
						{
							if (bufjdp - bufdp + datalen + 4096 > bufdp_len)
							{
								int nlen = bufdp_len;
								while (bufjdp - bufdp + datalen + 4096 > nlen) nlen <<= 1;
								char *tbuf = new char[nlen];
								SetVal(tbuf, bufdp, bufjdp - bufdp);
								bufjdp = tbuf + (bufjdp - bufdp);
								delete[] bufdp; bufdp = tbuf;
								bufdp_len = nlen;
							}
							if (j == 0)
							{
								*(byte*)bufjdp++ = (byte)ftype;
								if ((ftype >> 4) == 0xF) AppendTypedInt(bufjdp, flen);
							}
							FGet(g_filehandle[i][j], bufjdp, datalen);
							bufjdp += datalen;
						}
						else if (fkey == adid)
						{
							if (bufjad - bufad + datalen + 4096 > bufad_len)
							{
								int nlen = bufad_len;
								while (bufjad - bufad + datalen + 4096 > nlen) nlen <<= 1;
								char *tbuf = new char[nlen];
								SetVal(tbuf, bufad, bufjad - bufad);
								bufjad = tbuf + (bufjad - bufad);
								delete[] bufad; bufad = tbuf;
								bufad_len = nlen;
							}
							if (j == 0)
							{
								*(byte*)bufjad++ = (byte)ftype;
								if ((ftype >> 4) == 0xF) AppendTypedInt(bufjad, flen);
							}
							FGet(g_filehandle[i][j], bufjad, datalen);
							bufjad += datalen;
						}
						else
							FSeek(g_filehandle[i][j], datalen, SEEK_CUR);
					}

					FSeek(g_filehandle[i][j], next_line[j], SEEK_SET);
				}
				if (breakflag) break;
				sprintf(bufj, "%f", qual); while (*bufj) bufj++;
				if (pass) AppendString(bufj, "\tPASS\t\tGT");
				else AppendString(bufj, "\t.\t\tGT");
				if (bufjgq - bufgq) AppendString(bufj, ":GQ");
				if (bufjdp - bufdp) AppendString(bufj, ":DP");
				if (bufjad - bufad) AppendString(bufj, ":AD");

				*bufj++ = '\t';
				int64 writelen = sizeof(int64) * 2 + 1 + (bufj - buf) + (bufjgt - bufgt) + (bufjgq - bufgq) + (bufjdp - bufdp) + (bufjad - bufad);

				if (bufjgq > bufgq) byteflag |= 0x08; //hasgq
				if (bufjdp > bufdp) byteflag |= 0x10; //hasdp
				if (bufjad > bufad) byteflag |= 0x20; //hasad

				if (load_buf_size[ii % NBUF] < writelen)
				{
					delete[] load_buf[ii % NBUF];
					load_buf_size[ii % NBUF] = writelen + writelen / 2;
					load_buf[ii % NBUF] = new char[load_buf_size[ii % NBUF]];
				}

				char *writebuf = load_buf[ii % NBUF];
				char *writep = writebuf + sizeof(int64) * 2;
				writebuf[writelen - 1] = '\0';
				*(int64*)(writebuf                ) = writelen;
				*(int64*)(writebuf + sizeof(int64)) = byteflag;

				SetVal(writep, buf, bufj - buf);     writep += bufj - buf;
				SetVal(writep, bufgt, bufjgt - bufgt); writep += bufjgt - bufgt;
				SetVal(writep, bufgq, bufjgq - bufgq); writep += bufjgq - bufgq;
				SetVal(writep, bufdp, bufjdp - bufdp); writep += bufjdp - bufdp;
				SetVal(writep, bufad, bufjad - bufad); writep += bufjad - bufad;

				state_lock[ii % NBUF]++;

				if (ii++ == nloc) break;
			}
		}
		delete[] bufgt;
		delete[] bufgq;
		delete[] bufdp;
		delete[] bufad;
		delete[] buf;
		for (int i = 0; i < g_input_row; ++i)
			for (int j = 0; j < g_input_col; ++j)
				FSeek(g_filehandle[i][j], backup_offset[i * g_input_col + j], SEEK_SET);
		VLA_DELETE(next_line);
		VLA_DELETE(backup_offset);
	}

	/* Check input files rows and columns and count the number of variants */
	THREAD(GetVCFLines)
	{
		char **data = new char*[g_input_col];
		int64 *data_offset = new int64[g_input_col];
		int64 *data_used = new int64[g_input_col];
		for (int j = 0; j < g_input_col; ++j)
			data[j] = new char[64 * 1024 + 1];
		char *locus_name = new char[64 * 1024];
		int64 *uncompress_offset = new int64[g_input_col];
		int64 *num_loc = new int64[g_input_col];
		int64 *compress_offset = new int64[g_input_col];

		for (int i = 0; i < g_input_row; ++i)
		{
			int64 progress = PROGRESS_VALUE;//backup progress
			for (int j = 0; j < g_input_col; ++j)
			{
				uncompress_offset[j] = FTell(g_filehandle[i][j]);
				compress_offset[j] = FOffset(g_filehandle[i][j]);
				PROGRESS_VALUE += compress_offset[j];
				if (g_input_col > 1)
				{
					data_used[j] = FRead(data[j], 1, 64 * 1024, g_filehandle[i][j]);
					data[j][data_used[j]] = '\0';
				}
				data_offset[j] = 0;
				num_loc[j] = 0;
			}

			if (g_input_col == 1)
			{
				int64 tlen = 0;
				while (!FEof(g_filehandle[i][0]))
				{
					tlen = FRead(data[0], 1, 64 * 1024, g_filehandle[i][0]);
					if (g_format_val == -1)
						PROGRESS_VALUE = progress + FOffset(g_filehandle[i][0]);
					else
						PROGRESS_VALUE += tlen;
					data[0][tlen] = '\0';
					num_loc[0] += CountChar(data[0], '\n', tlen);
				}
				if (data[0][tlen - 1] != '\n')
					num_loc[0]--;
			}
			else for (;;) for (int j = 0; j < g_input_col; ++j)
			{
				int64 res = data_used[j] - data_offset[j];
				char *buf = data[j] + data_offset[j];
				char *bufnext = StrNextIdx(buf, '\n', 1, res) + 1;
				if (bufnext == (char*)1)
				{
					memcpy(data[j], buf, res);
					res = data_used[j] = res + FRead(data[j] + res, 1, 64 * 1024 - res, g_filehandle[i][j]);
					data[j][data_used[j]] = '\0';
					buf = data[j];
					data_offset[j] = 0;
					bufnext = StrNextIdx(buf, '\n', 1, res) + 1;
				}
				data_offset[j] = bufnext - data[j];

				if (g_format_val == -1)
				{
					int64 ncompress_offset = FOffset(g_filehandle[i][j]);
					PROGRESS_VALUE += ncompress_offset - compress_offset[j];
					compress_offset[j] = ncompress_offset;
				}
				else
					PROGRESS_VALUE += bufnext - buf;

				if (buf[0] != '\0' && buf[0] != '\r' && buf[0] != '\n')
					num_loc[j]++;
				else
					goto nextstep;

				if (j != 0)
				{
					if (LineCmp(locus_name, buf))
					{
						buf[strlen(locus_name)] = '\0';
						Exit("\nError: variant information in %s and %s mismatch, at lines %d:\n%s\nvs\n%s\n",
							g_filename[i][0], g_filename[i][j], num_loc[j], locus_name, buf);
					}
				}
				else// if (g_input_col > 1)
				{
					int64 infolen = StrNextIdx(buf, '\t', 5, res) - buf;
					memmove(locus_name, buf, infolen);
					locus_name[infolen] = '\0';
				}
			}

		nextstep: 
			for (int j = 0; j < g_input_col; ++j)
			{
				progress += g_filelen[i][j];
				if (!FEof(g_filehandle[i][j]))
					Exit("\nError: number of variants in %s mismatch that in %s.\n", g_filename[i][j], g_filename[i][0]);
				if (num_loc[j] != num_loc[0])
					Exit("\nError: number of variants in %s mismatch that in %s.\n", g_filename[i][j], g_filename[i][0]);
			}

			nloc += num_loc[0];
			PROGRESS_VALUE = progress;

			for (int j = 0; j < g_input_col; ++j)
				FSeek(g_filehandle[i][j], uncompress_offset[j], SEEK_SET);
		}

		for (int j = 0; j < g_input_col; ++j)
			delete[] data[j];
		delete[] data;
		delete[] data_offset;
		delete[] data_used;
		delete[] uncompress_offset;
		delete[] locus_name;
		delete[] num_loc;
		delete[] compress_offset;
	}

	/* Process lines from memory */
	THREAD(LoadVCF)
	{
		TABLE<HASH, uint> gfid(false, NULL);
		bool fmt_filter[1024] = { 0 };

		for (int64 ii = 0; ii < nloc; ii++)
		{
			while (ii >= progress1)
				SLEEP(SLEEP_TIME_TINY); 
			int64 avail_val = ii * 4 + 2, calc_val = avail_val + 1; 
			if (state_lock[ii % NBUF].compare_exchange_strong(avail_val, calc_val)) {
			gfid.Clear();
			int fmt_count = 0;

			char *str = load_buf[ii % NBUF];
			int64 res =      *(int64*)(str);
			int64 byteflag = *(int64*)(str + sizeof(int64));

			str += sizeof(int64) * 2;
			char *src1 = str;
			char *src2 = StrNextIdx(src1, '\t', 7, res - (src1 - str)) + 1;		//INFO
			char *dst = src2;
			src1 = StrNextIdx(src2, '\t', 1, res - (src2 - str)) + 1;			//FORMAT
			src2 = StrNextIdx(src1, '\t', 1, res - (src1 - str));				//GENOTYPES

			*dst++ = '\t';
			fmt_count = 0;
			bool writed = false;
			int gttag = -1;

			while (src1 > (char*)1 && src1 < src2)
			{
				ushort tagid = *(ushort*)src1;
				char tagnext = src1[2];

				if ((tagnext == ':' || tagnext == '\t') && 
					(tagid == 0x5447 || //GT
					 tagid == 0x5044 || //DP
					 tagid == 0x5147 || //GQ
					(tagid == 0x4441 && (ad || ploidyinfer)))) //AD
				{
					if (writed) *dst++ = ':';
					*(ushort*)dst = *(ushort*)src1;
					dst += 2; src1 += 2;
					if (tagid == 0x5447) gttag = fmt_count;
					fmt_filter[fmt_count++] = writed = true;
				}
				else
					fmt_filter[fmt_count++] = false;

				while (*src1 != ':' && *src1 != '\t') 
					src1++; 
				src1++;
			}

			*dst++ = '\t';
			src1 = src2 + 1;
			bool usedploidy[N_MAX_PLOIDY + 1] = { 0 };
			int count = 0;
			while (*src1)
			{
				count++;
				writed = false;
				for (int i = 0; i < fmt_count; ++i)
				{
					if (fmt_filter[i])
					{
						char *genostr = src1;
						if (writed) *dst++ = ':';

						while (*src1 != ':' && *src1 != '\0' && *src1 != '\t' && *src1 != '\n') 
							*dst++ = *src1++;

						if (gttag == i)
						{
							// genotype
							int64 len = src1 - genostr;
							int v = CountChars(genostr, "/|", len) + 1;
							if (v <= 0 || v > N_MAX_PLOIDY)
								Exit("\nError: ploidy level greater than %d in individual %s at %d-th locus.\n", N_MAX_PLOIDY, ainds[count]->name, ii + 1);
							
							//add missing at ploidy v
							if (!usedploidy[v])
							{
								int tid = gfid.size;
								gfid[missing_hash[v]] = tid;
								usedploidy[v] = true;
							}

							//add gfid if not missing
							if (genostr[0] != '.')
							{
								ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
								ReadVCFGenoString(alleles, genostr, v, ii, NULL);
								if (!haplotype) Sort(alleles, v); //phaseed, do not sort alleles
								
								HASH ha = HashGenotype(alleles, v);
								if (!gfid.ContainsKey(ha))
								{
									//add to hash and alleles array size
									int tid = gfid.size;
									gfid[ha] = tid;
									locus[ii].gasize += v + GetNalleles(alleles, v);
								}
							}
						}
						if (*src1 == '\t' || *src1 == '\n')
						{
							src1++;
							break;
						}
						else if (*src1 == '\0') break;
						src1++;
						writed = true;
					}
					else
					{
						//for (register char c = *src1; c != ':' && c != '\0' && c != '\t' && c != '\n'; c = *++src1);
						while (*src1 != ':' && *src1 != '\0' && *src1 != '\t' && *src1 != '\n') src1++;
						if (*src1 == '\t' || *src1 == '\n')
						{
							src1++;
							break;
						}
						else if (*src1 == '\0') break;
						src1++;
					}
				}
				if (count == nind) break;
				*dst++ = '\t';
			}

			*dst++ = '\n';
			*dst++ = '\0';

			uint64 gtsize = (uint64)CeilLog2(gfid.size);
			Lock(glock2);

			genotype_offset[ii].offset = genotype_coffset;
			genotype_offset[ii].size = gtsize;

			genotype_coffset += (uint64)ceil(gtsize * nind / 8.0) + 3;
			if (genotype_coffset > genotype_size)
			{
				uint64 nsize = ((uint64)(genotype_coffset * nloc / (ii + 1) * 1.1) >> 16) << 16;
				VAllocGenotype(nsize);
				genotype_size += nsize;
			}
			UnLock(glock2);

			////////////////////////////////////////////////////////////////////////////////////

			//alloc in one piece
			GENOTYPE *gtab;  ushort *gatab;
			if (ad != 2)
				new(&locus[ii]) LOCUS(locus_memory[threadid], str, (uint64)byteflag, gfid.size, gtab, gatab);//LoadVCF
			else
				str = StrNextIdx(str, '\t', 9) + 1;

			uint *depth = ploidyinfer ? new uint[locus[ii].k * nind] : NULL, *depthbak = depth;

			GENO_ITERATOR iter(0u, ii, false);//OK

			for (int j = 0; j < nind; ++j)
				ainds[j]->AddVCFGenotype(str, ii, depth, gfid, gtab, gatab, iter);

			iter.FinishWrite();

			if (ploidyinfer)
			{
				//add allele depth
				uint64 adsize = (uint64)CeilLog2((int64)locus[ii].maxdepth + 1);
				Lock(glock2);
				alleledepth_offset[ii].offset = alleledepth_coffset;
				alleledepth_offset[ii].size = adsize;

				alleledepth_coffset += (int64)ceil(adsize * locus[ii].k * nind / 8.0) + 3;
				if (alleledepth_coffset > alleledepth_size)
				{
					uint64 nsize = ((uint64)(alleledepth_coffset * nloc / (ii + 1) * 1.1) >> 16) << 16;
					VAllocAlleleDepth(nsize);
					alleledepth_size += nsize;
				}
				UnLock(glock2);
				IND::SetAlleleDepth(ii, depthbak, nind * locus[ii].k, 0u);
				delete[] depthbak;
			}

			PROGRESS_VALUE++;
			
			state_lock[ii % NBUF] = (ii + NBUF) * 4; }
		}
	}

	/* Read lines from VCF file */
	THREAD(LoadVCFGuard)
	{
		int64 &ii = progress1 = 0;
		VLA_NEW(backup_offset, int64, g_input_row * g_input_col);
		VLA_NEW(read_len, int64, g_input_col);
		VLA_NEW(line_array, char*, g_input_col);
		VLA_NEW(data_buffer, char*, g_input_col);
		VLA_NEW(data_offset, int64, g_input_col);
		VLA_NEW(data_size, int64, g_input_col);

		for (int j = 0; j < g_input_col; ++j)
			data_buffer[j] = new char[64 * 1024 + 1];

		for (int i = 0; i < g_input_row; ++i)
		{
			for (int j = 0; j < g_input_col; ++j)
			{
				backup_offset[i * g_input_col + j] = FTell(g_filehandle[i][j]);
				data_size[j] = FRead(data_buffer[j], 1, 64 * 1024, g_filehandle[i][j]);
				data_buffer[j][data_size[j]] = '\0';
				data_offset[j] = 0;  //fix data_offset
			}

			for (;;)
			{
				while (state_lock[ii % NBUF] != ii * 4) 
					SLEEP(SLEEP_TIME_TINY);
				state_lock[ii % NBUF]++;

				char *first_col_line = NULL; 
				int64 byte_flag = 0; 
				int64 first_col_remain_size = 0; 
				int64 write_len = 1 + 2 * sizeof(int64);

				for (int j = 0; j < g_input_col; ++j)
				{
					char *this_line = data_buffer[j] + data_offset[j];			   
					int64 remain_size = data_size[j] - data_offset[j];        
					char *next_line = StrNextIdx(this_line, '\n', 1, remain_size) + 1; 

					//do not have next line, move remaining data_buffer to begin and read a new batch of lines
					if (next_line == (char*)1)
					{
						memcpy(data_buffer[j], this_line, remain_size);
						remain_size = data_size[j] = remain_size + FRead(data_buffer[j] + remain_size, 1, 64 * 1024 - remain_size, g_filehandle[i][j]);
						data_buffer[j][data_size[j]] = '\0';
						this_line = data_buffer[j];
						data_offset[j] = 0;
						next_line = StrNextIdx(this_line, '\n', 1, remain_size) + 1;
					}

					//first col if many vcfs are used
					if (j == 0) 
					{ 
						first_col_line = this_line; 
						first_col_remain_size = remain_size; 
					}

					data_offset[j] = next_line - data_buffer[j];
					read_len[j] = next_line - this_line;

					//complete files in row
					if (read_len[j] == 0 || this_line[0] == NULL)
						break;

					//set line break to \0
					if (this_line[read_len[j] - 1] == '\n') this_line[--read_len[j]] = '\0';
					if (this_line[read_len[j] - 1] == '\r') this_line[--read_len[j]] = '\0';

					//save to a array, and process many vcfs in the same row
					line_array[j] = this_line;

					//test qual filter
					if (f_qual_b)
					{
						char *qualstr = StrNextIdx(this_line, '\t', 5, remain_size) + 1;
						double qual = -1;
						if (*qualstr != '.') qual = ReadDouble(qualstr);
						if (f_qual_b && qual != -1 && qual < f_qual_min && qual > f_qual_max)
							byte_flag |= 0x4;	//fail to pass qual filter
						//avg_qual += qual;
					}

					//test original filter
					if (f_original_b && f_original_val == 1)
					{
						if (LwrLineCmp("pass", StrNextIdx(this_line, '\t', 6, remain_size) + 1))
							byte_flag |= 0x2;	//fail to pass original filter
					}

					//check formats of the vcfs in the same row
					if (j != 0)
					{
						//skip text before genotype
						line_array[j] = StrNextIdx(this_line, '\t', 9, remain_size);
						read_len[j] -= line_array[j] - this_line;
						if (LineCmpAterm(StrNextIdx(first_col_line, '\t', 8, first_col_remain_size) + 1, StrNextIdx(this_line, '\t', 8, remain_size) + 1, '\t'))
						{
							char *fmta = StrNextIdx(first_col_line, '\t', 8, first_col_remain_size) + 1;
							char *fmtb = StrNextIdx(this_line, '\t', 8, remain_size) + 1;
							*StrNextIdx(fmtb, '\t', 1) = '\0';
							*StrNextIdx(fmta, '\t', 1) = '\0';
							*StrNextIdx(first_col_line, '\t', 3) = '\0';
							Exit("\nError: FORMAT field in %s and %s are different, at line %d, variant %s: \n%s\nvs\n%s\n", g_filename[i][0], g_filename[i][j], ii + 1, first_col_line, fmta, fmtb);
						}
					}
					write_len += read_len[j];
				}

				//read complete , break
				if (read_len[0] == 0 || first_col_line[0] == NULL)
					break;

				//avg_qual /= g_input_col;
				//float avg_qual_32 = (float)avg_qual;
				//*(int64*)(locus + ii) |= (((int64)*(int*)&avg_qual_32) << 32);

				if (load_buf_size[ii % NBUF] < write_len)
				{
					delete[] load_buf[ii % NBUF];
					load_buf_size[ii % NBUF] = write_len + write_len / 2;
					load_buf[ii % NBUF] = new char[load_buf_size[ii % NBUF]];
				}

				char *writebuf = load_buf[ii % NBUF];
				*(int64*)(writebuf                ) = write_len;
				*(int64*)(writebuf + sizeof(int64)) = byte_flag; //0 for pass
				writebuf[write_len - 1] = '\0';

				write_len = sizeof(int64) * 2;
				for (int j = 0; j < g_input_col; ++j)
				{
					SetVal(writebuf + write_len, line_array[j], read_len[j]);
					write_len += read_len[j];
				}
				writebuf[write_len++] = '\0';

				state_lock[ii % NBUF]++;

				if (ii++ == nloc) break;
			}
		}

		for (int j = 0; j < g_input_col; ++j)
			delete[] data_buffer[j];

		VLA_DELETE(read_len);
		VLA_DELETE(line_array);
		VLA_DELETE(data_buffer);
		VLA_DELETE(data_offset);
		VLA_DELETE(data_size);

		for (int i = 0; i < g_input_row; ++i)
			for (int j = 0; j < g_input_col; ++j)
				FSeek(g_filehandle[i][j], backup_offset[i * g_input_col + j], SEEK_SET);
		VLA_DELETE(backup_offset);
	}

	/* Calculate allele frequencies for each population and region */
	THREAD(CalcAlleleFreq)
	{
		if (allele_freq_offset)      delete[] allele_freq_offset;
		if (genotype_count_offset)   delete[] genotype_count_offset;

		//Allocate new offset
		allele_freq_offset = new LOCN[nloc];
		genotype_count_offset = new LOCN[nloc];
		SetFF(allele_freq_offset, nloc);
		SetFF(genotype_count_offset, nloc);

		//Calculate new offset and total number of alleles and genotypes
		KT = GT = maxK = maxG = 0;
		for (int64 l = 0; l < nloc; ++l)
		{
			allele_freq_offset[l] = KT;
			genotype_count_offset[l] = GT;
			KT += GetLoc(l).k;
			GT += GetLoc(l).ngeno;
			maxK = Max((int)GetLoc(l).k, maxK);
			maxG = Max((int)GetLoc(l).ngeno, maxG);
		}

		//Allocate memory for allele frequency and genotype count for each population and region
		for (int i = 0; i < pop.size; ++i)
			if (pop[i].nind && ad == 0)
				pop[i].AllocFreq();

		for (int rl = 0; rl < reg.size; ++rl)
			for (int i = 0; i < reg[rl].size; ++i)
				if (reg[rl][i].npop)
					reg[rl][i].AllocFreq();

		if (ad)
		{
			for (int rl = 0; rl < lreg; ++rl)
				for (int i = 0; i < nreg[rl]; ++i)
				{
					POP *cp = aregs[rl][i];
					POP **vp = cp->vpop;
					for (int p = 0; p < cp->npop; ++p)
						Add(cp->allelefreq, vp[p]->allelefreq, KT);
				}

			if (lreg > -1)
			{
				POP *cp = total_pop;
				POP **vp = cp->vpop;
				for (int p = 0; p < cp->npop; ++p)
					Add(cp->allelefreq, vp[p]->allelefreq, KT);
			}
		}

		for (int i = 0; i < npop; ++i)
			apops[i]->CalcFreqGcount();

		for (int rl = 0; rl < lreg; ++rl)
			for (int i = 0; i < nreg[rl]; ++i)
				aregs[rl][i]->CalcFreqGcount();

		if (lreg > -1) total_pop->CalcFreqGcount();
	}

	/* 1. Filter thread functions */

	/* Marker individual filtered or not */
	THREAD(MarkerIndividual)
	{
		atomic<int64> *ntype = new atomic<int64>[nind];
		int *ploidytab = new int[maxG * g_nthread_val];
		int *nallelestab = new int[maxG * g_nthread_val];

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int i = 0; i < nind; ++i)
		{
			ainds[i]->vt = 0;
			ainds[i]->vmin = 100;
			ainds[i]->vmax = 0;
		}

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			threadid = omp_get_thread_num();

			int *ploidy = ploidytab + maxG * threadid;
			int *nalleles = nallelestab + maxG * threadid;

			GENOTYPE *gtab = GetLoc(l).GetGtab();
			int ngeno = GetLoc(l).ngeno;

			for (int gi = 0; gi < ngeno; ++gi)
			{
				GENOTYPE &gt = gtab[gi];
				nalleles[gi] = gt.Nalleles();
				ploidy[gi] = gt.Ploidy();
			}

			GENO_ITERATOR iter(0u, l, true);//OK
			for (int j = 0; j < nind; ++j)
			{
				int gid = iter.Read();
				if (nalleles[gid])
				{
					byte v = ploidy[gid];
					AtomicAdd8(ainds[j]->vt, (int64)v);
					if (v < ainds[j]->vmin) AtomicMin1(ainds[j]->vmin, v);
					if (v > ainds[j]->vmax) AtomicMax1(ainds[j]->vmax, v);
					ntype[j]++;
				}
			}

			PROGRESS_VALUE++;
		}

		delete[] ploidytab;
		delete[] nallelestab;

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int i = 0; i < nind; ++i)
		{
			if (ainds[i] && f_ntype_b && (f_ntype_min > ntype[i] || ntype[i] > f_ntype_max))
				ainds[i] = NULL;

			if (ainds[i] && f_nploidy_b && (f_nploidy_min > minploidy || minploidy > f_nploidy_max))
				ainds[i] = NULL;

			if (ainds[i] && f_nploidy_b && (f_nploidy_min > maxploidy || maxploidy > f_nploidy_max))
				ainds[i] = NULL;

			PROGRESS_VALUE++;
		}

		delete[] ntype;
	}

	/* Marker locus filtered or not */
	THREAD(MarkerDiversity)
	{
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			DIVERSITY d;
			if (GetLoc(l).flag_pass)
			{
				d.CalcDiversity(l);

				if ((d.k == 0) ||
					(f_bmaf_b && d.k == 2 && (d.bmaf < f_bmaf_min || d.bmaf > f_bmaf_max)) ||
					(f_k_b && (d.k < f_k_min || d.k > f_k_max)) ||
					(f_n_b && (d.n < f_n_min || d.n > f_n_max)) ||
					(f_ptype_b && (d.ptype < f_ptype_min || d.ptype > f_ptype_max)) ||
					(f_pval_b && d.pval != -1 && (d.pval < f_pval_min || d.pval > f_pval_max)) ||
					(f_he_b && (d.he < f_he_min || d.he > f_he_max)) ||
					(f_ho_b && (d.ho < f_ho_min || d.ho > f_ho_max)) ||
					(f_pic_b && (d.pic < f_pic_min || d.pic > f_pic_max)) ||
					(f_ae_b && (d.ae < f_ae_min || d.ae > f_ae_max)))
					GetLoc(l).flag_pass = false;
			}
			else
				GetLoc(l).flag_pass = false;

			PROGRESS_VALUE++;
		}
	}

	/* Remove individual fail to pass filter */
	THREAD(RemoveIndividual)
	{
		MEMORY *oindividual_memory = individual_memory;
		individual_memory = new MEMORY[1];

		int nf = 0;
		for (int i = 0; i < nind; ++i)
			if (ainds[i]) nf++;

		if (nf == nind)
		{
			PROGRESS_VALUE = PROGRESS_CEND;
			delete[] individual_memory;
			individual_memory = oindividual_memory;
			return;
		}

		if (nf == 0)
			Exit("\nError: all %d individuals are excluded from analysis due to the anisoploid in the same contig, please check data.\n", nind);

		IND **newind = new IND*[nf];
		for (int i = 0, nid = 0; i < nind; ++i)
			if (ainds[i])
			{
				newind[nid] = new(individual_memory->Alloc(sizeof(IND))) IND(*ainds[i]);
				newind[nid]->indid = nid++;
			}

		//move genotype
		OFFSET *ngenotype_offset = new OFFSET[nloc];
		genotype_coffset = 0;

		OFFSET *nalleledepth_offset = ploidyinfer ? new OFFSET[nloc] : NULL;
		alleledepth_coffset = ploidyinfer ? 0 : alleledepth_coffset;

		for (int64 l = 0; l < nloc; ++l)
		{
			uint64 gtsize = genotype_offset[l].size;
			ngenotype_offset[l].offset = genotype_coffset;
			ngenotype_offset[l].size = gtsize;
			genotype_coffset += (int64)ceil(gtsize * nf / 8.0) + 3;//in bytes
			if (ploidyinfer)
			{
				int k2 = GetLoc(l).k;
				uint64 adsize = alleledepth_offset[l].size;
				nalleledepth_offset[l].offset = alleledepth_coffset;
				nalleledepth_offset[l].size = adsize;
				alleledepth_coffset += (int64)ceil(adsize * k2 * nf / 8.0) + 3;//in bytes
			}
		}

		uint64 ngenotype_size = genotype_coffset;
		byte *ngenotype_bucket = new byte[ngenotype_size];
		uint64 nalleledepth_size = alleledepth_coffset;
		byte *nalleledepth_bucket = ploidyinfer ? new byte[nalleledepth_size] : NULL;

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			threadid = omp_get_thread_num();
			GENO_ITERATOR rg(0u, l, true);//OK
			GENO_ITERATOR wg(0u, l, false, ngenotype_bucket, ngenotype_offset);//OK

			for (int i = 0; i < nind; ++i)
			{
				uint gid = rg.Read();
				if (ainds[i]) wg.Write(gid);
			}

			if (ploidyinfer)
			{
				GENO_ITERATOR rd(0u, l, true, alleledepth_bucket, alleledepth_offset);//OK
				GENO_ITERATOR wd(0u, l, false, nalleledepth_bucket, nalleledepth_offset);//OK
				int k2 = GetLoc(l).k;

				for (int i = 0; i < nind; ++i)
				{
					for (int k = 0; k < k2; ++k)
					{
						uint depth = rd.Read();
						if (ainds[i]) wd.Write(depth);
					}
				}
			}

			PROGRESS_VALUE++;
		}

		delete[] ainds; ainds = newind; nind = nf;
		delete[] genotype_offset; genotype_offset = ngenotype_offset;
		delete[] alleledepth_offset; alleledepth_offset = nalleledepth_offset;

		if (VIRTUAL_MEMORY)
		{
			VUnAllocGenotype();
			VUnAllocAlleleDepth();
		}
		else
		{
			delete[] genotype_bucket;
			delete[] alleledepth_bucket;
		}

		delete[] oindividual_memory;
		genotype_bucket = ngenotype_bucket;
		genotype_size = ngenotype_size;
		alleledepth_bucket = nalleledepth_bucket;
		alleledepth_size = nalleledepth_size;
	}

	/* Remove locus fail to pass filter */
	THREAD(RemoveLocus)
	{
		MEMORY *olocus_memory = locus_memory;
		locus_memory = new MEMORY[g_nthread_val];

		LOCUS *nlocus = useslocus ? NULL : new LOCUS[nfilter];
		SLOCUS *nslocus = useslocus ? new SLOCUS[nfilter] : NULL;

		//genotype id offset
		OFFSET *ngenotype_offset = new OFFSET[nfilter];
		genotype_coffset = 0;

		//allele depth offset
		OFFSET *nalleledepth_offset = NULL;
		if (ploidyinfer) nalleledepth_offset = new OFFSET[nfilter];
		alleledepth_coffset = 0;

		//genotype count, allele frequency
		LOCN *ngenotype_count_offset = NULL, *nallele_freq_offset = NULL;
		int64 *newid = new int64[nloc];

		//calculate new offsets
		for (int64 l = 0, nl = 0; l < nloc; ++l)
		{
			if (!GetLoc(l).flag_pass) continue;

			newid[l] = nl;
			uint64 gtsize = genotype_offset[l].size, gtlinesize = (uint64)ceil(gtsize * nind / 8.0) + 3;
			ngenotype_offset[nl].offset = genotype_coffset;
			ngenotype_offset[nl].size = gtsize;
			genotype_coffset += gtlinesize;

			if (ploidyinfer)
			{
				uint64 adsize = ploidyinfer ? alleledepth_offset[l].size : 0;
				uint64 adlinesize = (uint64)ceil(adsize * GetLoc(l).k * nind / 8.0) + 3;
				nalleledepth_offset[nl].offset = alleledepth_coffset;
				nalleledepth_offset[nl].size = adsize;
				alleledepth_coffset += adlinesize;
			}
			nl++;
		}

		//move genotype id and allele depth table
		byte *ngenotype_bucket = new byte[genotype_coffset];
		byte *nalleledepth_bucket = ploidyinfer ? new byte[alleledepth_coffset] : NULL;

		uint64 *nlocus_pos = NULL;
		if (useslocus && haplotype)
			nlocus_pos = new uint64[nfilter];

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			threadid = omp_get_thread_num();

			if (GetLoc(l).flag_pass)
			{
				int64 nl = newid[l];

				uint64 gtsize = genotype_offset[l].size, gtlinesize = (uint64)ceil(gtsize * nind / 8.0) + 3;
				memcpy(ngenotype_bucket + (ngenotype_offset[nl].offset), genotype_bucket + (genotype_offset[l].offset), gtlinesize);

				if (ploidyinfer)
				{
					uint64 adsize = ploidyinfer ? alleledepth_offset[l].size : 0;
					uint64 adlinesize = (uint64)ceil(adsize * GetLoc(l).k * nind / 8.0) + 3;
					memcpy(nalleledepth_bucket + (nalleledepth_offset[nl].offset), alleledepth_bucket + (alleledepth_offset[l].offset), adlinesize);
				}

				//deep copy locus
				if (useslocus)
				{
					new(&nslocus[nl]) SLOCUS(locus_memory[threadid], slocus[l]);//remove locus
					if (haplotype) nlocus_pos[nl] = locus_pos[l];
				}
				else
					new(&nlocus[nl]) LOCUS(locus_memory[threadid], nl, locus[l]);//remove locus
			}
			PROGRESS_VALUE++;
		}

		delete[] genotype_offset; 
		genotype_offset = ngenotype_offset;

		if (ploidyinfer)
		{
			delete[] alleledepth_offset; 
			alleledepth_offset = nalleledepth_offset;
		}

		if (VIRTUAL_MEMORY)
		{
			VUnAllocGenotype();
			VUnAllocAlleleDepth();
		}
		else
		{
			delete[] genotype_bucket;
			if (ploidyinfer) 
				delete[] alleledepth_bucket;
		}

		genotype_size = genotype_coffset;
		genotype_bucket = ngenotype_bucket;

		alleledepth_size = alleledepth_coffset;
		alleledepth_bucket = nalleledepth_bucket;

		delete[] newid;
		delete[] olocus_memory;

		if (useslocus)
		{
			delete[] slocus;
			slocus = nslocus;

			if (haplotype)
			{
				delete[] locus_pos;
				locus_pos = nlocus_pos;
			}
		}
		else 
		{ 
			delete[] locus; 
			locus = nlocus; 
		}

		nloc = nfilter;
	}

	/* 3. Haplotype extraction */

	/* Quick sort locus by contig and position */
	TARGET void QSLocus(int64 left, int64 right)
	{
		int64 i = left, j = right;

		int64 mid = (left + right) >> 1;
		LOCUS &pivot = locus[mid];
		SLOCUS &spivot = slocus[mid];

		while (left < j || i < right)
		{
			if (useslocus)
			{
				while (strcmp(slocus[i].GetChrom(), spivot.GetChrom()) < 0 || (strcmp(slocus[i].GetChrom(), spivot.GetChrom()) == 0 && GetLocPos(i) < GetLocPos(mid))) i++;
				while (strcmp(slocus[j].GetChrom(), spivot.GetChrom()) > 0 || (strcmp(slocus[j].GetChrom(), spivot.GetChrom()) == 0 && GetLocPos(j) > GetLocPos(mid))) j--;
				if (i <= j)
				{
					Swap(locus_pos[i], locus_pos[j]);
					Swap(locus_id[i], locus_id[j]);
					Swap(slocus[i], slocus[j]);
					i++; j--;
				}
			}
			else 
			{
				while (strcmp(locus[i].GetChrom(), pivot.GetChrom()) < 0 || (strcmp(locus[i].GetChrom(), pivot.GetChrom()) == 0 && locus[i].pos < pivot.pos)) i++;
				while (strcmp(locus[j].GetChrom(), pivot.GetChrom()) > 0 || (strcmp(locus[j].GetChrom(), pivot.GetChrom()) == 0 && locus[j].pos > pivot.pos)) j--;
				if (i <= j)
				{
					Swap(locus[i], locus[j]);
					i++; j--;
				}
			}

			if (i > j)
			{
				if (i == j + 2) PROGRESS_VALUE++;
				if (left < j)
				{
					QSLPAR par = { left, j };
					Lock(glock2);
					qslstack.Push(par);
					UnLock(glock2);
				}
				else if (left == j) PROGRESS_VALUE++;
				if (i < right)
				{
					QSLPAR par = { i, right };
					Lock(glock2);
					qslstack.Push(par);
					UnLock(glock2);
				}
				else if (i == right) PROGRESS_VALUE++;
				return;
			}
		}
	}

	/* Quick sort locus in a contig */
	THREAD(QSWorker)
	{
		QSLPAR par;
		while (PROGRESS_VALUE != PROGRESS_CEND)
		{
			Lock(glock2);
			if (qslstack.size)
			{
				par = qslstack.Pop();
				QSLocus(par.left, par.right);
			}
			UnLock(glock2);
			SLEEP(SLEEP_TIME_TINY);
		}
	}

	/* Quick sort extracted locus by contig and position */
	TARGET void QSHapLocus(int64 left, int64 right)
	{
		int64 i = left, j = right;

		int64 mid = (left + right) >> 1;
		HAPLO_DUMMY_LOCUS &pivot = haplotype_locus[mid];

		while (left < j || i < right)
		{
			while (haplotype_locus[i].chrid < pivot.chrid || (haplotype_locus[i].chrid == pivot.chrid && haplotype_locus[i].stpos < pivot.stpos)) i++;
			while (haplotype_locus[j].chrid > pivot.chrid || (haplotype_locus[j].chrid == pivot.chrid && haplotype_locus[j].stpos > pivot.stpos)) j--;

			if (i <= j)
			{
				Swap(haplotype_locus[i], haplotype_locus[j]);
				i++; j--;
			}

			if (i > j)
			{
				if (i == j + 2) PROGRESS_VALUE++;
				if (left < j)
				{
					QSLPAR par = { left, j };
					Lock(glock2);
					qslstack.Push(par);
					UnLock(glock2);
				}
				else if (left == j) PROGRESS_VALUE++;
				if (i < right)
				{
					QSLPAR par = { i, right };
					Lock(glock2);
					qslstack.Push(par);
					UnLock(glock2);
				}
				else if (i == right) PROGRESS_VALUE++;
				return;
			}
		}
	}

	/* Quick sort extracted locus in a contig */
	THREAD(QSHapWorker)
	{
		QSLPAR par;
		while (PROGRESS_VALUE != PROGRESS_CEND)
		{
			Lock(glock2);
			if (qslstack.size)
			{
				par = qslstack.Pop();
				QSHapLocus(par.left, par.right);
			}
			UnLock(glock2);
			SLEEP(SLEEP_TIME_TINY);
		}
	}

	/* Get number of alleles and genotypes at a dummy locus */
	TARGET double GetDummyK(int64 st, int64 ed, TABLE<HASH, ushort> &hfidx, TABLE<HASH, ushort> &gfidx)
	{
		hfidx.Clear(); gfidx.Clear();
		HASH hash[N_MAX_PLOIDY];
		ushort alleles[N_MAX_PLOIDY];
		bool usedploidy[N_MAX_PLOIDY] = { 0 };
		int ntyped = 0;

		for (int i = 0; i < nind; ++i)
		{
			int v = 0;
			HashHaplotype(ainds[i], st, ed, hash, v);

			if (!usedploidy[v])
			{
				usedploidy[v] = true;
				gfidx.PushIndex(missing_hash[v]);
			}

			if (hash[0] != (HASH)-1)
			{
				for (int vi = 0; vi < v; ++vi)
					alleles[vi] = (ushort)hfidx.PushIndex(hash[vi]);
				Sort(alleles, v); //unphase
				gfidx.PushIndex(HashGenotype(alleles, v));
				ntyped++;
			}
		}

		return ntyped / (double)nind;
	}

	/* Create locus for haplotype extraction */
	THREAD(CreateHaplotypeLocus)
	{
		int nthread = g_nthread_val;
		TABLE<HASH, ushort> hfidx(false, NULL), gfidx(false, NULL);
		LIST<HAPLO_DUMMY_LOCUS> dlocus(NULL);    //dummy locus

		//For each contig, find is begin and end variant
		for (int64 contig_st = -1, chrid = 0, contig_ed = -1; contig_st < nloc; )
		{
			int64 test_val = contig_st, next_val = contig_ed + 1;
			while (!haplotype_contig.compare_exchange_strong(test_val, next_val))
			{
				//jump a contig
				contig_st = contig_ed = contig_ed + 1;
				if (contig_st >= nloc) return;
				while (contig_ed < nloc && strcmp(GetLoc(contig_st).GetChrom(), GetLoc(contig_ed).GetChrom()) == 0) contig_ed++; contig_ed--;
				chrid++;
				test_val = contig_st; next_val = contig_ed + 1;
			}

			//1. Set i = 1 and j = 1 in the beginning.
			contig_st = contig_ed = contig_ed + 1;
			while (contig_ed < nloc && strcmp(GetLoc(contig_st).GetChrom(), GetLoc(contig_ed).GetChrom()) == 0) contig_ed++; contig_ed--;
			chrid++;
			int64 st = contig_st, ed = contig_st;

			//2. If the value of j exceeds the number of variants in this chromosome or contig, then terminate this algorithm. 
		step2:
			if (ed > contig_ed) continue;
			PROGRESS_VALUE = Max((uint64)PROGRESS_VALUE, (uint64)ed);

			//3. Calculate the values of several parameters (the length of each haplotype, the number of variants, the number of alleles and the number of genotypes). 
			int64 clen = GetLocPos(ed) - GetLocPos(st) + 1;
			int64 nvariants = ed - st + 1;
			double gtrate = GetDummyK(st, ed, hfidx, gfidx);

			//4. If any parameter exceeds the upper bound then increase i and go to step 2. 
			if (clen > haplotype_length_max ||
				nvariants > haplotype_variants_max ||
				gtrate < haplotype_ptype_min ||
				hfidx.size > haplotype_alleles_max ||
				gfidx.size > haplotype_genotypes_max)
			{
				st++;
				goto step2;
			}

			//5. If any parameter exceeds the lower bound then increase j and go to step 2. 
			if (clen < haplotype_length_min ||
				nvariants < haplotype_variants_min ||
				gtrate > haplotype_ptype_max ||
				hfidx.size < haplotype_alleles_min ||
				gfidx.size < haplotype_genotypes_min)
			{
				ed++;
				goto step2;
			}

			//6. Combine all variants between the variants iand j, and extract the haplotypes.
			HAPLO_DUMMY_LOCUS tentry { chrid, (int64)GetLocPos(st), st, ed, (int)nvariants, hfidx.size, gfidx.size };
			Lock(glock2);
			haplotype_locus.Push(tentry);
			UnLock(glock2);

			//7. Set i and j to the next applicable variant according to -haplotype_interval, and do to Step 2.
			st = ed = ed + haplotype_interval_val + 1;
			goto step2;
		}
		PROGRESS_VALUE = nloc;
	}

	/* Output locus for haplotype extraction */
	THREAD(WriteHaplotypeLocus)
	{
		int nthread = g_nthread_val;
		int64 newloc = haplotype_locus.size; 

		double nsec = newloc / (double)nthread + 1e-8;
		int64 st1 = (int64)(threadid * nsec), ed1 = (int64)((threadid + 1) * nsec);
		FILE *fout = TEMP_FILES[threadid];
		MEMORY haplo_memory;

		TABLE<HASH, HAPLO_DUMMY_HAPLOTYPE*> hftab(true, NULL);
		TABLE<HASH, TEMP_GENOTYPE> temptab(true, NULL);

		for (int64 l = st1; l < ed1; ++l)
		{
			hftab.Clear();    temptab.Clear();

			HAPLO_DUMMY_LOCUS &dlocus = haplotype_locus[l];
			HAPLO_DUMMY_LOCUS *dnext = l + 1 < newloc ? &haplotype_locus[l + 1] : NULL;

			int64 st = dlocus.st, ed = dlocus.ed, nvar = dlocus.nvar;
			int gasize = 0;
			GENO_ITERATOR wt(0u, l, false, haplotype_bucket, haplotype_offset);//OK

			for (int i = 0; i < nind; ++i)
			{
				HASH haplohash[N_MAX_PLOIDY];
				ushort alleles[N_MAX_PLOIDY] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
				int ploidy = 0;
				
				//get haplotype hash
				HashHaplotype(ainds[i], st, ed, haplohash, ploidy);

				//missing data
				if (haplohash[0] == (HASH)-1)
					SetFF(alleles, ploidy);

				//map hash into allele index
				else for (int ai = 0; ai < ploidy; ++ai)
				{
					HASH ha = haplohash[ai];
					HAPLO_DUMMY_HAPLOTYPE *ht = NULL;

					if ((ht = hftab.Try(ha)) == NULL)
					{
						haplo_memory.Alloc(ht, 1);
						hftab[ha] = ht;
						ht->ExtractHaplotype(ai, ainds[i], st, ed, nvar, (ushort)hftab.size - 1, haplo_memory);
					}
					alleles[ai] = ht->alleleid;
				}

				//unphase
				Sort(alleles, ploidy);

				//calculate genotype hash
				HASH hash = HashGenotype(alleles, ploidy);

				//add missing genotype
				if (!temptab.ContainsKey(missing_hash[ploidy]))
				{
					TEMP_GENOTYPE &tgt = temptab[missing_hash[ploidy]];
					SetVal(tgt.alleles, missing_array, N_MAX_PLOIDY);
					tgt.hash = missing_hash[ploidy];
					tgt.ploidy = ploidy;
					tgt.gid = temptab.size - 1;
					gasize += 0;
				}

				//add genotype
				if (!temptab.ContainsKey(hash))
				{
					TEMP_GENOTYPE &tgt = temptab[hash];
					SetVal(tgt.alleles, alleles, N_MAX_PLOIDY);
					tgt.hash = hash;
					tgt.ploidy = ploidy;
					tgt.gid = temptab.size - 1;
					gasize += ploidy + GetNalleles(alleles, ploidy);
				}

				//write genotype id
				wt.Write(temptab[hash].gid);
			}
			wt.FinishWrite();

			//create locus
			int ngeno = temptab.size;
			LOCUS *loc; SLOCUS *sloc;

			if (ngeno != haplotype_locus[l].gsize)
				Exit("\nError: haplotype extraction find more genotypes than expect.");

#define Loc (useslocus ? sloc : loc)
			if (useslocus)
				sloc = new(&haplotype_nslocus[l]) SLOCUS(locus_memory[threadid], slocus[st], l, ngeno, gasize, temptab);
			else
				loc = new(&haplotype_nlocus[l]) LOCUS(locus_memory[threadid], locus[st], l, ngeno, gasize, temptab);

			//write results
			int hsize = hftab.size;
			Loc->k = (ushort)hsize;
			fprintf(fout, "%s%sLocus:%s", g_linebreak_val, g_linebreak_val, Loc->GetName());
			fprintf(fout, "%s#CHROM:%s", g_linebreak_val, Loc->GetChrom());
			fprintf(fout, "%s#Variants:%lld", g_linebreak_val, nvar);
			fprintf(fout, "%s#Haplotypes:%d", g_linebreak_val, hftab.size);
			fprintf(fout, "%sRange:%lld-%lld", g_linebreak_val,  GetLocPos(dlocus.st), GetLocPos(dlocus.ed));
			fprintf(fout, "%sLength:%lld", g_linebreak_val, GetLocPos(dlocus.ed) - GetLocPos(dlocus.st) + 1);
			
			if (dnext && strcmp(GetLoc(dnext->st).GetChrom(), Loc->GetChrom()) == 0)
				fprintf(fout, "%sDistance to next extracted locus:%lld", g_linebreak_val, GetLocPos(dnext->st) - GetLocPos(dlocus.ed) - 1);
			else
				fprintf(fout, "%sNext extracted locus is in the different contig or chromosome", g_linebreak_val);
			fprintf(fout, "%sHapId", g_linebreak_val);

			for (int64 l2 = st; l2 <= ed; ++l2)
				fprintf(fout, "%c%s", g_delimiter_val, GetLoc(l2).GetName());

			for (int hi = 0; hi < hsize; ++hi)
				hftab(hi)->PrintHaplotype(fout, st, ed);
			haplo_memory.ClearMemory();

			PROGRESS_VALUE++;
		}
	}

	/* Find the optimal PES model for each locus */
	THREAD(GetLocusPESModel)
	{
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			GENOTYPE *gtab = GetLoc(l).GetGtab();
			int ngeno = GetLoc(l).ngeno;

			int model = 0;
			double maxli = -1e300;
			for (int m = 4; m <= N_DRE_MODELT; ++m)
			{
				double li = 0, li2 = 1;

				OpenLog(li, li2);
				for (int p = 0; p < npop; ++p)
				{
					ushort *gcount = cpop->GetGenoCount(l);
					double *freq = apops[p]->GetFreq(l);

					for (int gi = 0; gi < ngeno; ++gi)
					{
						GENOTYPE &gt = gtab[gi];
						if (gt.Nalleles() == 0 || gcount[gi] == 0) continue;

						double gfz = gt.GFZ(m, freq);
						gfz = gfz > 0 ? gfz : 1;

						if (gcount[gi] < 10)
							for (int ci = 0; ci < gcount[gi]; ++ci)
								ChargeLog(li, li2, gfz);
						else
							li += MyLog(gcount[gi]) * gcount[gi];
					}
				}
				CloseLog(li, li2);

				if (li > maxli)
				{
					maxli = li;
					model = m;
				}
			}

			GetLoc(l).pes_model = (byte)model;
			VLA_DELETE(gtab);

			PROGRESS_VALUE++;
		}
	}

	/* Calculate individual min and max ploidy, and sum ploidy levels */
	THREAD(AssignPloidyThread)
	{
		maxploidy = 0;  minploidy = 100;

		VLA_NEW(ploidytab, int, maxG * g_nthread_val);
		VLA_NEW(nallelestab, int, maxG * g_nthread_val);

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int i = 0; i < nind; ++i)
		{
			ainds[i]->vt = 0;
			ainds[i]->vmin = 100;
			ainds[i]->vmax = 0;
		}

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			threadid = omp_get_thread_num();

			int *nalleles = nallelestab + maxG * threadid;
			int *ploidy = ploidytab + maxG * threadid;

			GENOTYPE *gtab = GetLoc(l).GetGtab();
			int ngeno = GetLoc(l).ngeno;

			for (int gi = 0; gi < ngeno; ++gi)
			{
				GENOTYPE &gt = gtab[gi];
				nalleles[gi] = gt.Nalleles();
				ploidy[gi] = gt.Ploidy();
			}

			GENO_ITERATOR iter(0u, l, true);//OK
			for (int j = 0; j < nind; ++j)
			{
				int gid = iter.Read();
				byte v = ploidy[gid];
				if (nalleles[gid]) AtomicAdd8(ainds[j]->vt, (int64)v);
				if (v < ainds[j]->vmin) AtomicMin1(ainds[j]->vmin, v);
				if (v > ainds[j]->vmax) AtomicMax1(ainds[j]->vmax, v);
			}

			PROGRESS_VALUE++;
		}

#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int i = 0; i < nind; ++i)
		{
			if (ainds[i]->vmin < minploidy) AtomicMin1(minploidy, ainds[i]->vmin);
			if (ainds[i]->vmax > maxploidy) AtomicMax1(maxploidy, ainds[i]->vmax);
		}

		VLA_DELETE(ploidytab);
		VLA_DELETE(nallelestab);
	}

	/* 4. Conversion thread functions */

	/* Write convert genepop genotypes in a guard thread */
	THREAD(ConvertGenepopGuard)
	{
		for (int64 &ii = progress1 = 0; (int)ii < nind; ++ii, ++PROGRESS_VALUE)
		{
			GUARD_BEGIN

			if (ii == 0 || rinds[ii]->popid != rinds[ii - 1]->popid)
				fprintf(convert_file, "pop\r\n");
			fwrite(convert_buf[ii % NBUF], (int)strlen(convert_buf[ii % NBUF]), 1, convert_file);

			GUARD_END
		}
	}

	/* Write convert arlequin genotypes in a guard thread */
	THREAD(ConvertArlequinGuard)
	{
		for (int64 &ii = progress1 = 0; (int)ii < nind; ++ii, ++PROGRESS_VALUE)
		{
			GUARD_BEGIN

			if (ii == 0 || rinds[ii]->popid != rinds[ii - 1]->popid)
			{
				ushort popid = rinds[ii]->popid;
				if (ii) fprintf(convert_file, "\t\t}\r\n\r\n");
				fprintf(convert_file, "\t\tSampleName=\"%s\"\r\n\t\tSampleSize=%d\r\n\t\tSampleData={\r\n",
					apops[popid]->name, apops[popid]->nind);
			}
			fwrite(convert_buf[ii % NBUF], (int)strlen(convert_buf[ii % NBUF]), 1, convert_file);

			GUARD_END
		}
	}

	/* Write convert genotypes in a guard thread */
	THREAD(ConvertGuard)
	{
		for (int64 &ii = progress1 = 0; (int)ii < nind; ++ii, ++PROGRESS_VALUE)
		{
			GUARD_BEGIN

			fwrite(convert_buf[ii % NBUF], (int)strlen(convert_buf[ii % NBUF]), 1, convert_file);

			GUARD_END
		}
	}

	/* Convert genotype string */
	TARGET void PrepareGenotypeString(int format)
	{
#pragma omp parallel  for num_threads(g_nthread_val)  schedule(dynamic, 1)
		for (int64 l = 0; l < nloc; ++l)
		{
			threadid = omp_get_thread_num();
			GENOTYPE *gtab = GetLoc(l).GetGtab();
			int ngeno = GetLoc(l).ngeno;

			conversion_string[l].Clear();
			for (int gi = 0; gi < ngeno; ++gi)
			{
				char *str = NULL;
				
				switch (format)
				{
				case 1: str = gtab[gi].GetGenepopStr(); break;
				case 2: str = gtab[gi].GetSpagediStr(); break;
				case 3: str = gtab[gi].GetCervusStr(); break;
				case 4: str = gtab[gi].GetArlequinStr(); break;
				case 5: str = gtab[gi].GetStructureStr(); break;
				case 6: str = gtab[gi].GetPolygeneStr(); break;
				case 7: str = gtab[gi].GetPolyrelatednessStr(); break;
				}
				conversion_string[l].Push(str);
			}
		}
	}

	/* Convert into genepop format */
	TARGET void ConvertGenepop(int ntot, bool &isfirst)
	{
		struct tm *t1;
		time_t tt1;
		time(&tt1);
		t1 = localtime(&tt1);
		char filename[FILE_NAME_LEN];
		convert_file = FOpen(filename, "wb", "%s%s", g_output_val, ".convert.genepop.txt");
		fprintf(convert_file, "Genepop data file created by vcfpop %s on %04d-%02d-%02d %02d:%02d:%02d\r\n",
			VERSION, t1->tm_year + 1900, t1->tm_mon + 1, t1->tm_mday, t1->tm_hour, t1->tm_min, t1->tm_sec);

		for (int64 l = 0; l < nloc; ++l)
			fprintf(convert_file, "%s\r\n", GetLoc(l).GetName());

		convert_linesize = IND_NAME_LEN + 7 * nloc;
		for (int j = 0; j < NBUF; ++j)
			conversion_memory2->Alloc(convert_buf[j], convert_linesize);

		PrepareGenotypeString(1);

		RunThreads(&ConvertGenepopInd, &ConvertGenepopGuard, NULL, ntot, nind,
			"\nConverting population genetics software format:\n", g_nthread_val, isfirst);
		
		isfirst = false;
		fclose(convert_file);
	}

	/* Convert individual genotypes into genepop format in multiple threads */
	THREAD(ConvertGenepopInd)
	{
		for (int64 ii = 0; (int)ii < nind; ++ii)
		{
			THREAD_BEGIN

			char *str = convert_buf[ii % NBUF];
			sprintf(str, "%s,", rinds[ii]->name); while (*str) str++;
			for (int64 l = 0; l < nloc; ++l)
			{
				uint gid = rinds[ii]->GetGenotypeId(l);//fine
				GENOTYPE &gt = GetLoc(l).GetGtab()[gid];
				if (gt.Ploidy() != 2)
					Exit("\nError: Cannot convert genepop format, because the genotype of individual %s at locus %s is not diploid.\n", rinds[ii]->name, GetLoc(l).GetName());
				char *gstr = conversion_string[l][gid];
				do { *str++ = *gstr++; } while (*gstr);
			}
			sprintf(str, "\r\n");

			THREAD_END(2)
		}
	}

	/* Convert into spagedi format */
	TARGET void ConvertSpagedi(int ntot, bool &isfirst)
	{
		struct tm *t1;
		time_t tt1;
		time(&tt1);
		t1 = localtime(&tt1);
		char filename[FILE_NAME_LEN];
		convert_file = FOpen(filename, "wb", "%s%s", g_output_val, ".convert.spagedi.txt");
		fprintf(convert_file, "//spagedi data file created by vcfpop %s on %04d-%02d-%02d %02d:%02d:%02d\r\n",
			VERSION, t1->tm_year + 1900, t1->tm_mon + 1, t1->tm_mday, t1->tm_hour, t1->tm_min, t1->tm_sec);
		fprintf(convert_file, "//#individuals	#categories	#coordinates	#loci	#digits/allele	max ploidy\r\n");
		fprintf(convert_file, "%d	%d	0	%lu	3	%d\r\n-3\r\nind	pop", nind, npop, nloc, maxploidy);

		for (int64 l = 0; l < nloc; ++l)
			fprintf(convert_file, "\t%s", GetLoc(l).GetName());
		fprintf(convert_file, "\r\n");

		convert_linesize = IND_NAME_LEN + (maxploidy * 3 + 1) * nloc;
		for (int j = 0; j < NBUF; ++j)
			conversion_memory2->Alloc(convert_buf[j], convert_linesize);

		PrepareGenotypeString(2);

		RunThreads(&ConvertSpagediInd, &ConvertGuard, NULL, ntot, nind,
			"\nConverting population genetics software format:\n", g_nthread_val, isfirst);
		
		isfirst = false;
		fprintf(convert_file, "END");
		fclose(convert_file);
	}

	/* Convert individual genotypes into spagedi format in multiple threads */
	THREAD(ConvertSpagediInd)
	{
		for (int64 ii = 0; (int)ii < nind; ++ii)
		{
			THREAD_BEGIN

			char *str = convert_buf[ii % NBUF];
			sprintf(str, "%s\t%s", rinds[ii]->name, apops[rinds[ii]->popid]->name); while (*str) str++;
			for (int64 l = 0; l < nloc; ++l)
			{
				char *gstr = conversion_string[l][rinds[ii]->GetGenotypeId(l)];//fine
				do { *str++ = *gstr++; } while (*gstr);
			}
			sprintf(str, "\r\n");

			THREAD_END(2)
		}
	}

	/* Convert into cervus format */
	TARGET void ConvertCervus(int ntot, bool &isfirst)
	{
		char filename[FILE_NAME_LEN];
		convert_file = FOpen(filename, "wb", "%s%s", g_output_val, ".convert.cervus.csv");
		fprintf(convert_file, "ind,pop");

		for (int64 l = 0; l < nloc; ++l)
			fprintf(convert_file, ",%sA,%sB", GetLoc(l).GetName(), GetLoc(l).GetName());
		fprintf(convert_file, "\r\n");

		convert_linesize = IND_NAME_LEN + 8 * nloc;
		for (int j = 0; j < NBUF; ++j)
			conversion_memory2->Alloc(convert_buf[j], convert_linesize);

		PrepareGenotypeString(3);

		RunThreads(&ConvertCervusInd, &ConvertGuard, NULL, ntot, nind,
			"\nConverting population genetics software format:\n", g_nthread_val, isfirst);
		isfirst = false;
		fclose(convert_file);
	}

	/* Convert individual genotypes into cervus format in multiple threads */
	THREAD(ConvertCervusInd)
	{
		for (int64 ii = 0; (int)ii < nind; ++ii)
		{
			THREAD_BEGIN

			char *str = convert_buf[ii % NBUF];
			sprintf(str, "%s,%s", rinds[ii]->name, apops[rinds[ii]->popid]->name); while (*str) str++;
			for (int64 l = 0; l < nloc; ++l)
			{
				uint gid = rinds[ii]->GetGenotypeId(l);//fine
				GENOTYPE &gt = GetLoc(l).GetGtab()[gid];
				if (gt.Ploidy() != 2)
					Exit("\nError: Cannot convert cervus format, because the genotype of individual %s at locus %s is not diploid.\n", rinds[ii]->name, GetLoc(l).GetName());
				char *gstr = conversion_string[l][gid];
				do { *str++ = *gstr++; } while (*gstr);
			}
			sprintf(str, "\r\n");

			THREAD_END(2)
		}
	}

	/* Convert into arlequin format */
	TARGET void ConvertArlequin(int ntot, bool &isfirst)
	{
		struct tm *t1;
		time_t tt1;
		time(&tt1);
		t1 = localtime(&tt1);
		char filename[FILE_NAME_LEN];
		convert_file = FOpen(filename, "wb", "%s%s", g_output_val, ".convert.arlequin.arp");
		fprintf(convert_file, "#arlequin data file created by vcfpop %s on %04d-%02d-%02d %02d:%02d:%02d\r\n",
			VERSION, t1->tm_year + 1900, t1->tm_mon + 1, t1->tm_mday, t1->tm_hour, t1->tm_min, t1->tm_sec);
		fprintf(convert_file, "[Profile]\r\n\tTitle=\"vcfpop\"\r\n\tNbSamples=%d\r\n\tGenotypicData=1\r\n\tGameticPhase=0\r\n\tDataType=MICROSAT\r\n\tLocusSeparator=WHITESPACE\r\n\tMissingData=\"?\"\r\n\r\n# Locus Name", npop);

		for (int64 l = 0; l < nloc; ++l)
			fprintf(convert_file, "\r\n# %s", GetLoc(l).GetName());
		fprintf(convert_file, "\r\n\r\n[Data]\r\n\t[[Samples]]\r\n");

		convert_linesize = IND_NAME_LEN + 9 * nloc;
		for (int j = 0; j < NBUF; ++j)
			conversion_memory2->Alloc(convert_buf[j], convert_linesize);

		PrepareGenotypeString(4);

		RunThreads(&ConvertArlequinInd, &ConvertArlequinGuard, NULL, ntot, nind,
			"\nConverting population genetics software format:\n", g_nthread_val, isfirst);
		isfirst = false;

		fprintf(convert_file, "\t\t}\r\n\r\n");
		fprintf(convert_file, "[[Structure]]\r\n\r\n\t\tStructureName=\"Region\"\r\n\t\tNbGroups=%d\r\n\r\n", nreg[0]);
		for (int i = 0; i < nreg[0]; ++i)
		{
			POP *r = lreg >= 0 ? aregs[0][i] : NULL;
			POP **pops = r->vpop;
			fprintf(convert_file, "\t\tGroup={\r\n");
			for (int j = 0; j < r->npop; ++j)
				fprintf(convert_file, "\t\t\t\"%s\"\r\n", pops[j]->name);
			fprintf(convert_file, "\t\t}\r\n\r\n");
		}
		fclose(convert_file);
	}

	/* Convert individual genotypes into arlequin format in multiple threads */
	THREAD(ConvertArlequinInd)
	{
		for (int64 ii = 0; (int)ii < nind; ++ii)
		{
			THREAD_BEGIN

			char *str = convert_buf[ii % NBUF];
			sprintf(str, "\t\t\t%s 1", rinds[ii]->name); while (*str) str++;
			for (int64 l = 0; l < nloc; ++l)
			{
				uint gid = rinds[ii]->GetGenotypeId(l);//fine
				GENOTYPE &gt = GetLoc(l).GetGtab()[gid];
				if (gt.Ploidy() != 2)
					Exit("\nError: Cannot convert arlequin format, because the genotype of individual %s at locus %s is not diploid.\n", rinds[ii]->name, GetLoc(l).GetName());
				char *gstr = conversion_string[l][gid];//yes
				do { *str++ = *gstr++; } while (*gstr);
			}
			sprintf(str, "\r\n\t\t\t"); while (*str) str++;

			for (int64 l = 0; l < nloc; ++l)
			{
				char *gstr = conversion_string[l][rinds[ii]->GetGenotypeId(l)] + 5;//fine
				do { *str++ = *gstr++; } while (*gstr);
			}
			sprintf(str, "\r\n");

			THREAD_END(2)
		}
	}

	/* Convert into structure format */
	TARGET void ConvertStructure(int ntot, bool &isfirst)
	{
		struct tm *t1;
		time_t tt1;
		time(&tt1);
		t1 = localtime(&tt1);
		char filename[FILE_NAME_LEN];
		convert_file = FOpen(filename, "wb", "%s%s", g_output_val, ".convert.structure.txt");

		for (int64 l = 0; l < nloc; ++l)
			fprintf(convert_file, "%s ", GetLoc(l).GetName());
		FSeek(convert_file, -1, SEEK_CUR);
		fprintf(convert_file, "\r\n");

		convert_linesize = IND_NAME_LEN + (maxploidy * 4) * nloc;
		for (int j = 0; j < NBUF; ++j)
			conversion_memory2->Alloc(convert_buf[j], convert_linesize);

		PrepareGenotypeString(5);

		RunThreads(&ConvertStructureInd, &ConvertGuard, NULL, ntot, nind,
			"\nConverting population genetics software format:\n", g_nthread_val, isfirst);
		isfirst = false;
		fclose(convert_file);
	}

	/* Convert individual genotypes into structure format in multiple threads */
	THREAD(ConvertStructureInd)
	{
		for (int64 ii = 0; (int)ii < nind; ++ii)
		{
			THREAD_BEGIN
			char *str = convert_buf[ii % NBUF];
			int ploidy = maxploidy;

			for (int k = 0; k < ploidy; ++k)
			{
				sprintf(str, "%s %s", rinds[ii]->name, apops[rinds[ii]->popid]->name); while (*str) str++;
				for (int64 l = 0; l < nloc; ++l)
				{
					uint gid = rinds[ii]->GetGenotypeId(l);//fine
					GENOTYPE &gt = GetLoc(l).GetGtab()[gid];
					if (gt.Ploidy() != maxploidy)
						Exit("\nError: Cannot convert structure format, because the ploidy level of individual %s at locus %s is different.\n", rinds[ii]->name, GetLoc(l).GetName());
					sprintf(str, " %s", conversion_string[l][gid] + k * 5); while (*str) str++;
				}
				sprintf(str, "\r\n"); while (*str) str++;
			}
			THREAD_END(2)
		}
	}

	/* Convert into polygene format */
	TARGET void ConvertPolygene(int ntot, bool &isfirst)
	{
		char filename[FILE_NAME_LEN];
		convert_file = FOpen(filename, "wb", "%s%s", g_output_val, ".convert.polygene.txt");
		fprintf(convert_file, "ID\tPop\tPloidy");

		for (int64 l = 0; l < nloc; ++l)
			fprintf(convert_file, "\t%s", GetLoc(l).GetName());
		fprintf(convert_file, "\r\n");

		convert_linesize = IND_NAME_LEN + maxploidy * 4 * nloc;
		for (int j = 0; j < NBUF; ++j)
			conversion_memory2->Alloc(convert_buf[j], convert_linesize);

		PrepareGenotypeString(6);

		RunThreads(&ConvertPolygeneInd, &ConvertGuard, NULL, ntot, nind,
			"\nConverting population genetics software format:\n", g_nthread_val, isfirst);
		isfirst = false;
		fclose(convert_file);
	}

	/* Convert individual genotypes into polygene format in multiple threads */
	THREAD(ConvertPolygeneInd)
	{
		for (int64 ii = 0; (int)ii < nind; ++ii)
		{
			THREAD_BEGIN

			char *str = convert_buf[ii % NBUF];
			int ploidy = rinds[ii]->GetGenotype(0).Ploidy();//fine
			sprintf(str, "%s\t%s\t%d", rinds[ii]->name, apops[rinds[ii]->popid]->name, ploidy); while (*str) str++;

			for (int64 l = 0; l < nloc; ++l)
			{
				uint gid = rinds[ii]->GetGenotypeId(l);//fine
				GENOTYPE &gt = GetLoc(l).GetGtab()[gid];
				if (gt.Ploidy() != ploidy)
					Exit("\nError: Cannot convert polygene format, because the ploidy of individual %s at locus %s is different.\n", rinds[ii]->name, GetLoc(l).GetName());
				char *gstr = conversion_string[l][gid];//yes
				do { *str++ = *gstr++; } while (*gstr);
			}
			sprintf(str, "\r\n");

			THREAD_END(2)
		}
	}

	/* Convert into polyrelatedness format */
	TARGET void ConvertPolyrelatedness(int ntot, bool &isfirst)
	{
		struct tm *t1;
		time_t tt1;
		time(&tt1);
		t1 = localtime(&tt1);
		char filename[FILE_NAME_LEN];
		convert_file = FOpen(filename, "wb", "%s%s", g_output_val, ".convert.polyrelatedness.txt");
		fprintf(convert_file, "//configuration\r\n//#alleledigits(1~4)\t#outputdigits(0~10)\t#missingallele\t#ambiguousallele\t#nthreads(1~64)\r\n3\t8\t000\t999\t2\r\n//genotype\r\nInd\tPop");
		for (int64 l = 0; l < nloc; ++l)
			fprintf(convert_file, "\t%s", GetLoc(l).GetName());
		fprintf(convert_file, "\r\n");

		convert_linesize = IND_NAME_LEN + (maxploidy * 3 + 1) * nloc;
		for (int j = 0; j < NBUF; ++j)
			conversion_memory2->Alloc(convert_buf[j], convert_linesize);

		PrepareGenotypeString(7);

		RunThreads(&ConvertPolyrelatednessInd, &ConvertGuard, NULL, ntot, nind,
			"\nConverting population genetics software format:\n", g_nthread_val, isfirst);
		isfirst = false;

		fprintf(convert_file, "//end of file");
		fclose(convert_file);
	}

	/* Convert individual genotypes into polyrelatedness format in multiple threads */
	THREAD(ConvertPolyrelatednessInd)
	{
		for (int64 ii = 0; (int)ii < nind; ++ii)
		{
			THREAD_BEGIN

			char *str = convert_buf[ii % NBUF];
			sprintf(str, "%s\t%s", rinds[ii]->name, apops[rinds[ii]->popid]->name); while (*str) str++;
			for (int64 l = 0; l < nloc; ++l)
			{
				char *gstr = conversion_string[l][rinds[ii]->GetGenotypeId(l)];//fine
				do { *str++ = *gstr++; } while (*gstr);
			}
			sprintf(str, "\r\n");

			THREAD_END(2)
		}
	}

	/* 5. Diversity thread functions */

	/* Add and write genetic diversity */
	THREAD(DiversityGuard)
	{
		int ni = cpop->nind;
		diversity_sum.Init();

		bool addsum = diversity_level_val[diversity_stage];
		bool writelocus = diversity_level_val[diversity_stage + 3];
		char *popname = cpop == total_pop ? (char*)"Total" : cpop->name;

		for (int64 &ii = progress1 = 0; ii < nloc; ii++)
		{
			while (state_lock[ii % NBUF] >> 1 != ii * 2 + 1) 
				SLEEP(SLEEP_TIME_TINY);


			PROGRESS_VALUE += ni;

			if (addsum)
				diversity_sum.Add(diversity_buf[ii % NBUF]);

			if (writelocus)
				diversity_buf[ii % NBUF].WriteLocus(TEMP_FILES[1], popname);

			GUARD_END
		}
	}

	/* Calculate genetic diversity using multiple threads */
	THREAD(DiversityThread)
	{
		for (int64 ii = 0; ii < nloc; ii++)
		{
			THREAD_BEGIN

			diversity_buf[ii % NBUF].CalcDiversity(ii);

			THREAD_END(2)
		}
	}

	/* 6. Individual Statistics thread function */

	/* Write individual statistics 
	THREAD(IndividualStatisticsGuard)
	{
		for (int64 &ii = progress1 = 0; ii < nind; ii++)
		{
			GUARD_BEGIN

			fprintf(FRES, "%s", indstat_buf[ii % NBUF]);
			PROGRESS_VALUE++;

			GUARD_END
		}
	}*/

	/* Calculate individual statistics using multiple threads */
	THREAD(IndividualStatisticsThread)
	{
		uint nthread = g_nthread_val;
		double nsec = nind / (double)nthread + 1e-8;
		uint st = (uint64)(threadid * nsec), ed = (uint64)((threadid + 1) * nsec);

		for (uint i = st; i < ed; ++i)
		{ 
			ainds[i]->PrintIndividualStatistics(TEMP_FILES[threadid]);
			PROGRESS_VALUE++;
		}
	}

	/* 7. Genetic Differentiation thread function */

	/* Calculate genetic differentiation using multiple threads */
	THREAD(GeneticDifferentiationThread)
	{
		//load ind
		int nthread = g_nthread_val;
		int64 progress = 0;

		FST *fst_buf2 = fst_buf[fst_type];
		if (fst_type == 1)
		{
			//among regions
			for (int rl = 0; rl < lreg; ++rl)
				if (progress++ % nthread == threadid)
					fst_buf2[lreg - 1 - rl].CalcFst(aregs[rl], nreg[rl]);
		}
		if (fst_type == 2)
		{
			//among pops
			if (progress++ % nthread == threadid)
				fst_buf2[0].CalcFst(apops, npop);
		}
		if (fst_type == 3)
		{
			//among pops/regs in region
			for (int rl = lreg - 1, p = 0; rl >= 0; --rl)
				for (int i = 0; i < nreg[rl]; ++i, ++p)
					if (progress++ % nthread == threadid)
						fst_buf2[p].CalcFst(aregs[rl][i]->vpop, aregs[rl][i]->npop);
		}
		if (fst_type == 4)
		{
			//between regions
			for (int rl = lreg - 1; rl >= 0; --rl)
			{
				int n = nreg[rl];
				for (int i = 0; i < n; ++i)
					for (int j = i + 1; j < n; ++j)
						if (progress++ % nthread == threadid)
						{
							fst_buf2[i * n + j].CalcFst(aregs[rl][i], aregs[rl][j]);
							fst_buf2[j * n + i] = fst_buf2[i * n + j];
						}
				fst_buf2 += n * n;
			}
		}
		if (fst_type == 5)
		{
			//between pops
			int n = npop;
			for (int i = 0; i < n; ++i)
				for (int j = i + 1; j < n; ++j)
					if (progress++ % nthread == threadid)
					{
						fst_buf2[i * n + j].CalcFst(apops[i], apops[j]);
						fst_buf2[j * n + i] = fst_buf2[i * n + j];
					}
		}
	}

	/* 8. Genetic Distance thread function */

	/* Write column format genetic distance results in a guard thread */
	THREAD(GeneticDistanceGuard1)
	{
		int64 &ii = progress1 = 0;
		int n = gdist_type == 1 ? nind : (gdist_type == 2 ? npop : nreg[gdist_type - 3]);
		int64 nadd2 = gdist_type == 1 ? 1 : 100;
		if (gdist_fmt_val[2])
			GDIST::ColumnPrintHeader();

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j, ++ii)
			{
				GUARD_BEGIN

				if (j >= i && gdist_fmt_val[2])
					gdist_buf[ii % NBUF].ColumnPrintLine(i, j);

				PROGRESS_VALUE += nadd2;

				GUARD_END
			}
		}
	}

	/* Write matrix format genetic distance results in a guard thread */
	THREAD(GeneticDistanceGuard2)
	{
		int64 &ii = progress2 = 0;
		uint nk = gdist_type == 1 ? N_GD_ESTIMATOR - 2 * N_FST_ESTIMATOR : N_GD_ESTIMATOR;
		uint n = gdist_type == 1 ? nind : (gdist_type == 2 ? npop : nreg[gdist_type - 3]);
		int64 nadd2 = gdist_type == 1 ? 1 : 100;

		if (gdist_fmt_val[1])
			for (uint k = 1; k <= nk; ++k)
				GDIST::MatrixPrintMatrixHeader(k, n);

		for (uint i = 0; i < n; ++i)
		{
			if (gdist_fmt_val[1])
				for (uint k = 1; k <= nk; ++k)
					GDIST::MatrixPrintRowHeader(k, i);

			for (uint j = 0; j < n; ++j, ++ii)
			{
				GUARD_BEGIN

				GDIST &gd = gdist_buf[ii % NBUF];
				if (gdist_fmt_val[1])
					for (uint k = 1; k <= nk; ++k)
						gd.MatrixPrintCell(k);

				PROGRESS_VALUE += nadd2;

				GUARD_END
			}
		}
	}

	/* Calculate genetic distance using multiple threads */
	THREAD(GeneticDistanceThread)
	{
		int64 ii = 0;
		if (gdist_type == 1)
		{
			VLA_NEW(p1, double, maxK);
			VLA_NEW(p2, double, maxK);
			for (int i = 0; i < nind; ++i)
			{
				for (int j = 0; j < nind; ++j, ++ii)
				{
					THREAD_BEGIN2

					gdist_buf[ii % NBUF].CalcGD(ainds[i], ainds[j], p1, p2);

					THREAD_END(1)
				}
			}
			VLA_DELETE(p1);
			VLA_DELETE(p2);
		}
		else if (gdist_type >= 2)
		{
			int rl = gdist_type - 3;
			int n = gdist_type == 2 ? npop : nreg[rl];
			POP **tpop = gdist_type == 2 ? apops : aregs[rl];
			VLA_NEW(p1, double, maxK);
			for (int i = 0; i < n; ++i)
			{
				for (int j = 0; j < n; ++j, ++ii)
				{
					THREAD_BEGIN2

					gdist_buf[ii % NBUF].CalcGD(tpop[i], tpop[j], p1);

					THREAD_END(1)
				}
			}
			VLA_DELETE(p1);
		}
	}

	/* 9. AMOVA thread function */

	/* Calculate AMOVA using multiple threads */
	THREAD(AMOVAThread)
	{
		switch (amova_cmethod_val)
		{
		case 1: amova_buf[threadid].CalcAMOVA_homo(); break;
		case 2: amova_buf[threadid].CalcAMOVA_aniso(); break;
		case 3: amova_buf[threadid].CalcAMOVA_ml(); break;
		}
	}

	/* 10. Population assignment thread function */

	/* Calculate population assignment using multiple threads */
	THREAD(PopulationAssignmentThread)
	{
		//load ind
		double nsec = nind / (double)g_nthread_val + 1e-8;
		int st = (int)(threadid * nsec), ed = (int)((threadid + 1) * nsec);
		for (int i = st; i < ed; ++i)
		{
			ainds[i]->PrintAssignment(TEMP_FILES[threadid]);

			PROGRESS_VALUE++;
		}
	}

	/* 11. Relatedness thread function */

	/* Write column format relatedness coefficient results in a guard thread */
	THREAD(RelatednessGuard1)
	{
		int64 &ii = progress1 = 0;
		int n = cpop->nind;
		if (relatedness_fmt_val[2])
			RELATEDNESS::ColumnPrintHeader();

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j, ++ii)
			{
				GUARD_BEGIN

				RELATEDNESS& re = relatedness_buf[ii % NBUF];
				if (j >= i && relatedness_fmt_val[2])
					re.ColumnPrintLine(i, j);

				PROGRESS_VALUE++;

				GUARD_END
			}
		}
	}

	/* Write matrix format relatedness coefficient results in a guard thread */
	THREAD(RelatednessGuard2)
	{
		int64 &ii = progress2 = 0;
		int n = cpop->nind;

		if (relatedness_fmt_val[1])
			for (int k = 1; k <= N_RELATEDNESS_ESTIMATOR; ++k)
				RELATEDNESS::MatrixPrintMatrixHeader(k, n);

		for (int i = 0; i < n; ++i)
		{
			if (relatedness_fmt_val[1])
				for (int k = 1; k <= N_RELATEDNESS_ESTIMATOR; ++k)
					RELATEDNESS::MatrixPrintRowHeader(k, i);

			for (int j = 0; j < n; ++j, ++ii)
			{
				GUARD_BEGIN

				RELATEDNESS& re = relatedness_buf[ii % NBUF];

				if (relatedness_fmt_val[1])
					for (int k = 1; k <= N_RELATEDNESS_ESTIMATOR; ++k)
						re.MatrixPrintCell(k);

				PROGRESS_VALUE++;

				GUARD_END
			}
		}
	}

	/* Calculate relatedness coefficient using multiple threads */
	THREAD(RelatednessThread)
	{
		//load ind
		int64 ii = 0;
		int ni = cpop->nind;

		if (relatedness_estimator_val[8] || relatedness_estimator_val[9])  Anderson2007_Coef = new double[nloc * 3];
		if (relatedness_estimator_val[11]) Huang2015_Coef = new double[nloc * 9];

		for (int i = 0; i < ni; ++i)
		{
			for (int j = 0; j < ni; ++j, ++ii)
			{
				THREAD_BEGIN2

				relatedness_buf[ii % NBUF].CalcRelatedness(ainds[i], ainds[j]);

				THREAD_END(1)
			}
		}

		if (relatedness_estimator_val[8] || relatedness_estimator_val[9])  delete[] Anderson2007_Coef;
		if (relatedness_estimator_val[11]) delete[] Huang2015_Coef;
	}

	/* 12. Kinship thread function */

	/* Write column format kinship coefficient results in a guard thread */
	THREAD(KinshipGuard1)
	{
		int64 &ii = progress1 = 0;
		int n = cpop->nind;
		if (kinship_fmt_val[2])
			KINSHIP::ColumnPrintHeader();

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j, ++ii)
			{
				GUARD_BEGIN

				KINSHIP& ks = kinship_buf[progress1 % NBUF];
				if (j >= i && kinship_fmt_val[2])
					ks.ColumnPrintLine(i, j);

				PROGRESS_VALUE++;

				GUARD_END
			}
		}
	}

	/* Write matrix format kinship coefficient results in a guard thread */
	THREAD(KinshipGuard2)
	{
		int64 &ii = progress2 = 0;
		int n = cpop->nind;

		if (kinship_fmt_val[1])
			for (int k = 1; k <= N_KINSHIP_ESTIMATOR; ++k)
				KINSHIP::MatrixPrintMatrixHeader(k, n);

		for (int i = 0; i < n; ++i)
		{
			if (kinship_fmt_val[1])
				for (int k = 1; k <= N_KINSHIP_ESTIMATOR; ++k)
					KINSHIP::MatrixPrintRowHeader(k, i);

			for (int j = 0; j < n; ++j, ++progress2)
			{
				GUARD_BEGIN

				KINSHIP& ks = kinship_buf[progress2 % NBUF];
				if (kinship_fmt_val[1])
					for (int k = 1; k <= N_KINSHIP_ESTIMATOR; ++k)
						ks.MatrixPrintCell(k);

				PROGRESS_VALUE++;

				GUARD_END
			}
		}
	}

	/* Calculate kinship coefficient using multiple threads */
	THREAD(KinshipThread)
	{
		int64 ii = 0;
		int ni = cpop->nind;

		for (int i = 0; i < ni; ++i)
		{
			for (int j = 0; j < ni; ++j, ++ii)
			{
				THREAD_BEGIN2

				kinship_buf[ii % NBUF].CalcKinship(ainds[i], ainds[j]);

				THREAD_END(1)
			}
		}
	}

	/* 13. Principal Coordinate Analysis thread function */
	/* 14. Hierarchical clustering */

	/* Calculate genetic distance for PCoA and Hierarchical clustering using multiple threads */
	THREAD(PCoAClusteringThread)
	{
		//load ind
		int nthread = g_nthread_val;
		int type = gdindex[0];
		int n = type == 1 ? nind : (type == 2 ? npop : nreg[type - 3]);
		int64 progress = 0;
		byte *estimator = GDIST_METHOD == 2 ? pcoa_estimator_val : cluster_estimator_val;
		GDIST tbuf;

		VLA_NEW(p1, double, maxK);
		VLA_NEW(p2, double, maxK);
		for (int i = 0; i < n; ++i)
		{
			for (int j = i + 1; j < n; ++j)
			{
				if (progress++ % nthread != threadid) continue;
				switch (type)
				{
				case 1:  tbuf.CalcGD(ainds[i], ainds[j], p1, p2); break;
				case 2:  tbuf.CalcGD(apops[i], apops[j], p1); break;
				case 3:
				default: tbuf.CalcGD(aregs[type - 3][i], aregs[type - 3][j], p1);  break;
				}
				for (int k = 1; k <= N_GD_ESTIMATOR; ++k)
					if (estimator[k])
						*(pcoa_matrix + n * n * gdindex[k] + j * n + i) =
						*(pcoa_matrix + n * n * gdindex[k] + i * n + j) =
						*(&tbuf.Nei1972 + k - 1);

				PROGRESS_VALUE++;
			}
		}
		VLA_DELETE(p1);
		VLA_DELETE(p2);
	}

	/* Test. spatial pattern, do not use */
	TARGET double SPA_Fij(double *a, int i)
	{
		double re = 1.0 / (exp(-SumProd(a, spa_x + i * spa_np, spa_np)) + 1.0);
		if (re < spa_truncate_val) re = spa_truncate_val;
		if (re > spa_truncate_val2) re = spa_truncate_val2;
		return re;
	}

	TARGET double SPA_Likelihood(CPOINT &xx, void **Param)
	{
		double *a = xx.image;
		uint64 l = *(uint64*)Param[0];
		return SPA_Likelihood(l, a);
	}

	TARGET double SPA_Likelihood(uint64 l, double *a)
	{
		// L = sum_i (ad_1 lnfij + ad_2 ln(1-fij))
		double re = 0;
		double f1 = 0, f2 = 0;
		if (spa_level_val == 1)
		{
			GENOTYPE *gtab = GetLoc(l).GetGtab();
			GENO_ITERATOR rt(0u, l, true);

			for (int i = 0; i < nind; ++i)
			{
				GENOTYPE &gt = gtab[rt.Read()];
				if (!spa_valid[i]) continue;
				f1 = SPA_Fij(a, i);
				f2 = log(1.0 - f1);
				f1 = log(f1);

				ushort *als = gt.GetAlleleArray();
				for (int j = 0, v = gt.Ploidy(); j < v; ++j)
					re += als[j] == 0 ? f1 : f2;
			}
		}
		else for (int p = 0; p < npop; ++p)
		{
			if (!spa_valid[p]) continue;
			f1 = SPA_Fij(a, p);
			f2 = log(1.0 - f1);
			f1 = log(f1);

			int nh = apops[p]->loc_stat[l].nhaplo;
			int n1 = (int)(apops[p]->GetFreq(l, 0) * nh + 0.5);
			re += f1 * n1 + f2 * (nh - n1);
		}
		return re;
	}

	TARGET double SPA_Hessian(uint64 l, double *G, double *H, double *a, double *a2)
	{
		double eps = 1e-6;
		double j00 = SPA_Likelihood(l, a);
		for (int j = 0; j < spa_np; ++j)
		{
			SetVal(a2, a, spa_np);
			a2[j] += eps;
			double j01 = SPA_Likelihood(l, a2);
			double d11 = (j01 - j00) / eps;
			G[j] = d11;
			for (int k = j; k < spa_np; ++k)
			{
				SetVal(a2, a, spa_np);
				a2[k] += eps;
				double j10 = SPA_Likelihood(l, a2);

				a2[j] += eps;
				double j11 = SPA_Likelihood(l, a2);
				double d12 = (j11 - j10) / eps;
				H[j * spa_np + k] = H[k * spa_np + j] = (d12 - d11) / eps;
			}
		}
		return j00;
	}

	TARGET void SPA_Locus(uint64 l, double *x, double *f, double *a, double *a2, double *am, double *g, double *h)
	{
		//At least dim + 1 = 3 individuals should be genotyped

		// L = sum_i (ad_1 lnfij + ad_2 ln(1-fij))
		/*
			x1 a1 + y1 a2 + 1 a3 = -ln(1 / p1 - 1)
			x2 a1 + y2 a2 + 1 a3 = -ln(1 / p2 - 1)
			x3 a1 + y3 a2 + 1 a3 = -ln(1 / p3 - 1)
			x A = f
			A = inv(x) * b
		*/

		//Initial Value
		spa_tn = 0;
		if (spa_level_val == 1)
		{
			GENOTYPE *gtab = GetLoc(l).GetGtab();
			GENO_ITERATOR rt(0u, l, true);

			for (int i = 0; i < nind; ++i)
			{
				GENOTYPE &gt = gtab[rt.Read()];
				if (gt.Nalleles() == 0)
				{
					spa_valid[i] = false;
					continue;
				}
				spa_valid[i] = true;
				int ac = 0, v = gt.Ploidy();
				ushort *als = gt.GetAlleleArray();
				for (int j = 0; j < v; ++j)
					if (als[j] == 0)
						ac++;
				f[spa_tn++] = ac / (double)v;
			}
		}
		else for (int p = 0; p < npop; ++p)
		{
			if (apops[p]->loc_stat[l].nhaplo == 0)
			{
				spa_valid[p] = false;
				continue;
			}
			spa_valid[p] = true;
			f[spa_tn++] = apops[p]->GetFreq(l, 0);
		}
		if (spa_tn < spa_np) return;

		//rearrange coordinate
		for (int i = 0, c = 0; i < spa_n; ++i)
		{
			if (!spa_valid[i]) continue;
			for (int j = 0; j < spa_np; ++j)
				x[j * spa_tn + c] = spa_x[i * spa_np + j];
			c++;
		}

		//truncate
		for (int i = 0; i < spa_tn; ++i)
		{
			if (f[i] < spa_truncate_val) f[i] = spa_truncate_val;
			if (f[i] > spa_truncate_val2) f[i] = spa_truncate_val2;
			f[i] = -log(1.0 / f[i] - 1.0);
		}

		
		Map<MatrixXd> B(f, spa_tn, 1);
		//may be wrong, matrix saved in column format
		Map<MatrixXd> X(x, spa_tn, spa_np);

		MatrixXd A = X.householderQr().solve(B);
		for (int i = 0; i < spa_np; ++i)
			a[i] = A(i, 0);

		//Downhill simplex
		int dim = spa_np;
		void *Param[] = { (void*)&l };
		CPOINT xx0 = CPOINT::DownHillSimplex(dim, 0, false, 0.0001, 15, SPA_Likelihood, Param);
		SetVal(am, xx0.image, spa_np);
		am[spa_np + 2] = xx0.li;

		//Newton's method
		Map<MatrixXd> H(h, spa_np, spa_np);//ok, symmetric
		Map<MatrixXd> G(g, spa_np, 1);//ok, column vector
		SetVal(a2, a, spa_np);
		a2[0] += 1.0;
		double eps = 1e-6, likelihood = 0;
		for (int m = 0; m < 100; ++m)
		{
			likelihood = SPA_Hessian(l, g, h, a, a2);
			if (likelihood > am[spa_np + 2])
			{
				SetVal(am, a, spa_np);
				am[spa_np + 2] = likelihood;
			}

			bool flag = true;
			for (int j = 0; j < spa_np; ++j)
				if (abs(a[j] - a2[j]) > eps)
					flag = false;
			SetVal(a2, a, spa_np);
			MatrixXd D = H.lu().solve(G);
			for (int j = 0; j < spa_np; ++j)
				a[j] -= D(j, 0);
			if (flag) break;
		}

		double ex2 = 0, ex = 0;
		for (int i = 0; i < spa_n; ++i)
		{
			if (!spa_valid[i]) continue;
			double p = SPA_Fij(am, i);
			ex2 += p * p;
			ex += p;
		}
		ex2 /= spa_tn;
		ex /= spa_tn;
		am[spa_np] = sqrt((ex2 - ex * ex) * spa_tn);
		am[spa_np + 1] = sqrt((ex2 - ex * ex) * spa_tn / (spa_tn - 1));

		FILE *fout = TEMP_FILES[threadid];
		fprintf(fout, "%s", GetLoc(l).GetName());
		if (spa_odepth_val == 1)
		{
			if (spa_level_val == 1)
			{
				GENOTYPE *gtab = GetLoc(l).GetGtab();
				GENO_ITERATOR rt(0u, l, true);

				for (int i = 0; i < nind; ++i)
				{
					GENOTYPE &gt = gtab[rt.Read()];
					if (!spa_valid[i])
					{
						fprintf(fout, "%c", g_delimiter_val);
						continue;
					}
					ushort *als = gt.GetAlleleArray();
					int v = gt.Ploidy(), n1 = 0;
					for (int j = 0; j < v; ++j)
						if (als[j] == 0) n1++;
					fprintf(fout, "%c", g_delimiter_val);
					fprintf(fout, "%d/%d", n1, v - n1);
				}

			}
			else for (int p = 0; p < npop; ++p)
			{
				if (!spa_valid[p])
				{
					fprintf(fout, "%c", g_delimiter_val);
					continue;
				}
				int nh = apops[p]->loc_stat[l].nhaplo;
				int n1 = (int)(apops[p]->GetFreq(l, 0) * nh + 0.5);
				fprintf(fout, "%c", g_delimiter_val);
				fprintf(fout, "%d/%d", n1, nh - n1);
			}
		}
		if (spa_ofreq_val == 1) for (int i = 0; i < spa_n; ++i)
		{
			if (!spa_valid[i])
			{
				fprintf(fout, "%c", g_delimiter_val);
				continue;
			}
			double p = SPA_Fij(am, i);
			fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, p);
			fprintf(fout, "/");
			WriteReal(fout, 1.0 - p);
		}
		for (int i = 0; i < spa_np + 3; ++i)
		{
			fprintf(fout, "%c", g_delimiter_val);
			WriteReal(fout, am[i]);
		}
		fprintf(fout, "%s", g_linebreak_val);
	}

	/* Calculate spatial pattern using multiple threads */
	THREAD(SPAThread)
	{
		FILE *fout = TEMP_FILES[threadid];
		if (threadid == 0)
		{
			fprintf(fout, "%s%sLocus", g_linebreak_val, g_linebreak_val);
			if (spa_odepth_val == 1)
			{
				if (spa_level_val == 1) for (int i = 0; i < nind; ++i)
				{
					fprintf(fout, "%c", g_delimiter_val);
					fprintf(fout, "%s_depth", ainds[i]->name);
				}
				else for (int p = 0; p < npop; ++p)
				{
					fprintf(fout, "%c", g_delimiter_val);
					fprintf(fout, "%s_depth", apops[p]->name);
				}
			}
			if (spa_ofreq_val == 1)
			{
				if (spa_level_val == 1) for (int i = 0; i < nind; ++i)
				{
					fprintf(fout, "%c", g_delimiter_val);
					fprintf(fout, "%s_freq", ainds[i]->name);
				}
				else for (int p = 0; p < npop; ++p)
				{
					fprintf(fout, "%c", g_delimiter_val);
					fprintf(fout, "%s_freq", apops[p]->name);
				}
			}
			for (int i = 0; i < spa_dim_val; ++i)
				fprintf(fout, "%ca%d", g_delimiter_val, i + 1);
			fprintf(fout, "%cb%cSPA%cSD%cLikelihood%s", g_delimiter_val, g_delimiter_val, g_delimiter_val, g_delimiter_val, g_linebreak_val);
		}

		int nthread = g_nthread_val;
		double nsec = nloc / (double)nthread + 1e-8;
		uint64 st1 = (uint64)(threadid * nsec), ed1 = (uint64)((threadid + 1) * nsec);

		VLA_NEW(x, double, spa_n * spa_np);
		VLA_NEW(f, double, spa_n);
		VLA_NEW(g, double, spa_np);
		VLA_NEW(h, double, spa_np * spa_np);
		VLA_NEW(a, double, spa_dim_val + 3);
		VLA_NEW(a2, double, spa_np);
		VLA_NEW(am, double, spa_dim_val + 4);
		spa_valid = new bool[spa_n];

		for (uint64 l = st1; l < ed1; ++l)
		{
			SPA_Locus(l, x, f, a, a2, am, g, h);
			
			PROGRESS_VALUE++;
		}

		VLA_DELETE(x);
		VLA_DELETE(f);
		VLA_DELETE(g);
		VLA_DELETE(h);
		VLA_DELETE(a);
		VLA_DELETE(a2);
		VLA_DELETE(am);
		delete[] spa_valid;
	}

	/* 15. Bayesian clustering */

	/* Calculate Bayesian clustering using multiple threads */
	THREAD(StructureThread)
	{
		for (int i = 0; i < structure_totalruns; ++i)
		{
			if (structure_par[i].flag.test_and_set()) continue;

			STRUCTURE Structure;
			Structure.ReadPar(structure_par + i);
			Structure.MCMC();
			Structure.PrintStructure();
		}
	}

	/* 16. Ploidy Inference thread function */

	/* Calculate ploidy inference using multiple threads */
	THREAD(PloidyInferenceThread)
	{
		int nthread = g_nthread_val;
		double nsec = nind / (double)nthread + 1e-8;
		int st = (int)(threadid * nsec), ed = (int)((threadid + 1) * nsec);
		for (int i = st; i < ed; ++i)
		{
			ainds[i]->PrintPloidyInference(TEMP_FILES[threadid]);
			PROGRESS_VALUE++;
		}
	}
#endif

#undef extern
#pragma pack(pop)