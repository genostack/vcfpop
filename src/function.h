/* Functions */

#pragma once
#include "vcfpop.h"

#pragma pack(push, 1)

/* Temporatory genotype, will be compressed into GENOTYPE */
struct TEMP_GENOTYPE
{
	HASH hash;								//Genotype hash
	int ploidy;								//Ploidy level
	int gid;								//Genotype id
	ushort alleles[N_MAX_PLOIDY];			//Indexed alleles
};

/* Memory offset for compat bit-wise storage */
struct OFFSET
{
	uint64 offset : 48;						//Offset, some base pointer add this offset can obtain data at a locus
	uint64 size : 16;						//Size in bits of each piece of data, use piece wise storage
};

/* Header class of BCF file */

class BCFHEADER
{
public:
	int format_gtid;						//Index of GT field
	int format_gqid;						//Index of GQ field
	int format_dpid;						//Index of DP field
	int format_adid;						//Index of AD field
	int filter_passidx;						//Pass filter index
	char **contig_name;						//Config names
	int64 contig_size;						//Contig size
	int64 nsample;							//Number of sample

	//##FILTER=<ID=PASS,Description="All filters passed",IDX=0>
	//##contig=<ID=1,assembly=b37,length=249250621,IDX=0>
	//##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype",IDX=1>
	//##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants",IDX=2>

	/* Uninitialize */
	TARGET ~BCFHEADER();

	/* Initialize */
	TARGET BCFHEADER(char *header);
};

/* Individual genetic distance */
struct INDGD
{
public:
	double ABtype;							//Number of loci both genotyped
	double Jx1, Jx2, Jxy;					//Nei1972; Nei1974;
	double Cavalli1967;						//Cavalli1967
	double t1, t2;							//Reynolds1983;
	double Nei1983;							//Nei1983
	double Euclidean;						//Euclidean
	double Goldstein1995;					//Goldstein1995
	double Roger1972;						//Roger1972
};

/* Quick sort locus */
struct QSLPAR
{
	int64 left;								//Left bound of quick sort range
	int64 right;							//Right bound of quick sort range
};

/* Infomation of a run in Bayesian Clustering */
struct SRUNINFO
{
	atomic_flag flag;						//Is task been allocated to a thread
	int k;									//Number of ancestral clusters
	int id;									//Run index
	int rep;								//Replicate index
	double MeanlnL;							//Mean lnLSetAlleleDepth
	double VarlnL;							//Var lnL
	double lnPD;							//Ln Prob of Data
};

struct GENOTYPE
{
public:
			//4 bytes
	uint offset : 24;						//this pointer + offset * 2 is address of alleles array
			//ushort *alleles;						//Alleles copy in ascending order [ploidy]
											//Unique alleles [nalleles] order by dosage descending
	uint patternid : 8;						//Pattern index

	/* Do nothing */
	TARGET GENOTYPE();

	/* Create genotype from alleles and hash */
	TARGET GENOTYPE(ushort *&gatab, ushort *alleles, int ploidy);

	/* Copy from a reference genotype */
	TARGET GENOTYPE(ushort *&gatab, GENOTYPE &ref);

	/* Get allele copy at ith haplotype */
	TARGET ushort GetAlleleCopy(int i);

	/* Get allele array */
	TARGET ushort *GetAlleleArray();

	/* Set allele array */
	TARGET void SetAlleleArray(ushort *alleles);

	/* Get pattern code */
	TARGET uint64 GetPattern();

	/* Number of alleles */
	TARGET int Nalleles();

	/* Ploidy level */
	TARGET int Ploidy();

	/* Crc32 hash */
	TARGET HASH Hash();

	/* Heterozygosity in this genotype */
	TARGET double HIndex();

	/* SS within genotype under IAM model */
	TARGET double SS_IAM();

	/* SS within genotype under IAM model */
	TARGET double SS_SMM(ushort *alen2, int k2);

	/* The multinomial coefficient for HWE/RCS genotype frequency */
	TARGET double HWECoef();

	/* Genotypic frequency for zygotes under inbreeding */
	TARGET double GFZ(int *allele_count, int sum, double f);

	/* Genotypic frequency for zygotes under specific double-reduction model */
	TARGET double GFZ(int DR_MODE, double *f);

	/* Number of copies of target allele */
	TARGET int GetAlleleCount(int a);

	/* Frequency of target allele in this genotype */
	TARGET double GetFreq(int a);

	/* Frequencies of all alleles in this genotype */
	TARGET void GetFreq(double *p, int k2);

	/* Obtain Genepop genotype string */
	TARGET char *GetGenepopStr();

	/* Obtain Arlequin genotype string */
	TARGET char *GetArlequinStr();

	/* Obtain Arlequin genotype string */
	TARGET char *GetStructureStr();

	/* Obtain Spagedi genotype string */
	TARGET char *GetSpagediStr();

	/* Obtain PolyRelatedness genotype string */
	TARGET char *GetPolyrelatednessStr();

	/* Obtain Polygene genotype string */
	TARGET char *GetPolygeneStr();

	/* Obtain Cervus genotype string */
	TARGET char *GetCervusStr();
};

class SLOCUS
{
	//small locus, 12 bytes

public:
	uint64 bits1 : 48;						//Genotype table[ngeno]
											//alen[k + k * k] {for non-vcf SMM distance}
											//chrom \0 name \0 {(allele identifiers \0)[k] for vcf/bcf)}
											//genotype alleles[gasize]
	uint64 k : 16;							//Number of alleles

	uint   ngeno : 22;						//Number of genotypes
	uint   flag_pass : 1;					//0 pass filter
	uint   flag_alen : 1;					//Has alen table?
	uint   pes_model : 8;					//1 for RCS, 2 for PRCS, 3 for CES, 4+ for PES

	/* Initialize */
	TARGET SLOCUS();

	/* Convert from LOCUS */
	TARGET SLOCUS(MEMORY &memory, LOCUS& ref);

	/* Create unphase locus */
	TARGET SLOCUS(MEMORY &memory, SLOCUS& ref, TABLE<HASH, uint> &gitab, ushort *gtmap);

	/* Create locus for haplotype extraction and Chi-square test */
	TARGET SLOCUS(MEMORY &memory, SLOCUS& ref, int64 _id, int _ngeno, int _gasize, TABLE<HASH, TEMP_GENOTYPE> &temptab);

	/* Deep copy from SLOCUS */
	TARGET SLOCUS(MEMORY &memory, SLOCUS& ref);

	/* Get Genotype array */
	TARGET GENOTYPE *GetGtab();

	/* Get end of chrom \0 name \0 {(allele identifiers \0)[k] */
	TARGET char *GetEnd();

	/* Get chrom string */
	TARGET char *GetChrom();

	/* Get locus identifier */
	TARGET char *GetName();

	/* Get allele name for vcf/bcf */
	TARGET char *GetAlleleName(int a);

	/* Get alen array */
	TARGET ushort *GetAlenArray();

	/* Get genotype allele array */
	TARGET ushort *GetGenoAlleleArray();

	/* Get Genotype alleles array size Sum(ploidy+nalleles)*/
	TARGET uint GetGenoAlleleSize();

	/* Get SMM distance */
	TARGET ushort GetSMMDist(int a, int b);

	//TARGET void GetGfTab(TABLE<HASH, GENOTYPE*> &gftab);
};

class LOCUS : public SLOCUS
{
	//+ 32 + 35 bytes
public:
	TABLE<HASH, GENOTYPE*> gftab;			//Genotype table, do not need release
	char *_chrom;							//Chrom \0 name \0 {(allele identifiers \0)[k] for vcf/bcf)}
	ushort *_alen;							//alen[k + k * k] {for non-vcf SMM distance}
	ushort *_genoallele;
	
	uint64 pos;								//Position in chromosome	
	LOCN id;								//Index
	uint maxdepth;							//Max allele depth
	uint gasize;							//Size of genotype alleles sum(ploidy + nalleles)

	ushort gtid;							//Index of GT tag
	ushort gqid;							//Index of GQ tag
	ushort dpid;							//Index of DP tag
	ushort adid;							//Index of AD tag
	ushort format_size;						//Format offset

	byte flag_original : 1;					//1 pass original filter
	byte flag_qual : 1;						//2 qual
	byte flag_indel : 1;					//3 isindel
	byte flag_hasgq : 1;					//4 hasgq
	byte flag_hasdp : 1;					//5 hasdp
	byte flag_hasad : 1;					//7 hasad

	/* Initialize */
	TARGET LOCUS();

	/* Deep Copy Locus, SLOCUS */
	TARGET LOCUS(MEMORY &memory, int64 _id, LOCUS& ref);

	/* Create unphase locus */
	TARGET LOCUS(MEMORY &memory, LOCUS& ref, TABLE<HASH, uint> &gitab, ushort *gtmap);

	/* Create locus for haplotype extraction and Chi-square test, SLOCUS */
	TARGET LOCUS(MEMORY &memory, LOCUS& ref, int64 _id, int _ngeno, int _gasize, TABLE<HASH, TEMP_GENOTYPE> &alstab);

	/* For dummy locus for collapse alleles during testing genotype distributions */
	TARGET LOCUS(MEMORY &memory, SLOCUS& ref);

	/* For non-vcf input, set locus name and id, SLOCUS */
	TARGET LOCUS(MEMORY &memory, char *line, int64 _id, int _ngenotype, GENOTYPE *&gtab, ushort *&gatab);

	/* For vcf input, set locus name and id, SLOCUS */
	TARGET LOCUS(MEMORY &memory, char *&line, uint64 _mask, int _ngenotype, GENOTYPE *&gtab, ushort *&gatab);

	/* Get end of chrom \0 name \0 {(allele identifiers \0)[k] */
	TARGET char *GetEnd();

	/* Get chrom string */
	TARGET char *GetChrom();

	/* Get locus identifier */
	TARGET char *GetName();

	/* Get allele name for vcf/bcf */
	TARGET char *GetAlleleName(int a);

	/* Get alen array */
	TARGET ushort *GetAlenArray();

	/* Get genotype allele array */
	TARGET ushort *GetGenoAlleleArray();

	/* Get Genotype alleles array size Sum(ploidy+nalleles)*/
	TARGET uint GetGenoAlleleSize();

	/* Get index of a target format */
	TARGET ushort GetFormatId(char *base, char *fmt_name, ushort *fmt_offset);
};

class GENO_ITERATOR
{
public:
	uint *pos;								//Current read pointer
	uint64 data;							//Readed bits
	uint64 mask64;							//Mask to extract genotype id bits
	uint size;								//Number of bits a genotype id used
	uint mask;								//Mask to extract genotype id bits
	uint nbits;								//Number of bits remaining in data
	uint nskipbits;							//Nskipbits for first write back 

	/* Do nothing */
	TARGET GENO_ITERATOR();

	/* Initialize reader/writer */
	TARGET GENO_ITERATOR(int indid, int64 l, bool isread, byte *bucket = NULL, OFFSET *offset = NULL);

	/* Get id of next ind (order by indid) */
	TARGET int Read();

	/* Write id of next ind to buffer */
	TARGET void Write(uint gtid);

	/* Write all remaining bits to buffer */
	TARGET void FinishWrite();
};

struct SPF
{
	//Used in ploidy infer
	int count;
	double val1[N_MAX_PLOIDY + 1];	//
	double val2[N_MAX_PLOIDY + 1];	//
};

class IND
{
public:
	char *name;								//Individual identifier
	byte vmax;								//Maximum ploidy level
	int indid;								//Individual id
	int popid;								//Population id
	byte vmin;								//Minimum ploidy level
	int64 vt;								//Sum of allele copies across loci, missing genotype do not account

	/* Initialize */
	TARGET IND();

	/* Create individual for non-vcf input */
	TARGET IND(char *t, bool iscount, int id, GENOTYPE **gtab, ushort **gatab, GENO_ITERATOR *iter);

	/* Create individual for vcf input */
	TARGET IND(char *&title, int id);

	/* Create individual from a reference */
	TARGET IND(IND& ref);

	/* Unnitialize */
	TARGET ~IND();

	/* Set allele sequencing depth, for ad, test */
	TARGET static void SetAlleleDepth(int64 l, uint *depth, int K, int indid);

	/* Set allele sequencing depth, for ad, test */
	TARGET void SetAlleleDepth(int64 l, uint *depth, int K);

	/* Set allele sequencing depth, for ad, test */
	TARGET void SetAlleleDepth(int64 l, uint *depth, int K, OFFSET *_offset, byte *bucket);

	/* Set allele sequencing depth, for ad, test */
	TARGET void GetAlleleDepth(int64 l, uint *depth);

	/* Set individual genotype with default bucket */
	TARGET void SetGenotype(int64 l, uint gid);

	/* Set individual genotype with local bucket */
	TARGET void SetGenotype(int64 l, uint gid, OFFSET *_offset, byte *bucket);

	/* Get index for a pair of genotype */
	TARGET static void GetDyadGenotypeIdx(int &id1, int &id2, int64 l);

	/* Get individual genotype index from default table */
	TARGET int GetGenotypeId(int64 l);

	/* Get individual genotype from default table */
	TARGET GENOTYPE &GetGenotype(int64 l);

	/* Get individual genotype from local table */
	TARGET GENOTYPE &GetGenotype(int64 l, TABLE<HASH, GENOTYPE*> &gftab);

	/* Get individual genotype from local table */
	TARGET GENOTYPE &GetGenotype(int64 l, GENOTYPE *gtab);

	/* Create individual from genepop */
	TARGET void genepop(char *t, bool iscount, GENOTYPE **gtab, ushort **gatab, GENO_ITERATOR *iter);

	/* Create individual from cervus */
	TARGET void cervus(char *t, bool iscount, GENOTYPE **gtab, ushort **gatab, GENO_ITERATOR *iter);

	/* Create individual from spagedi */
	TARGET void spagedi(char *t, bool iscount, GENOTYPE **gtab, ushort **gatab, GENO_ITERATOR *iter);

	/* Create individual from arlequin */
	TARGET void arlequin(char *t, bool iscount, GENOTYPE **gtab, ushort **gatab, GENO_ITERATOR *iter);

	/* Create individual from structure */
	TARGET void structure(char *t, bool iscount, GENOTYPE **gtab, ushort **gatab, GENO_ITERATOR *iter);

	/* Create individual from polyrelatedness */
	TARGET void polyrelatedness(char *t, bool iscount, GENOTYPE **gtab, ushort **gatab, GENO_ITERATOR *iter);

	/* Create individual from polygene */
	TARGET void polygene(char *t, bool iscount, GENOTYPE **gtab, ushort **gatab, GENO_ITERATOR *iter);

	/* Calculate the likelihood of genotype data */
	TARGET double Likelihood(POP *grp, int model, int64 loc, double e);

	/* Calculate the individual kinship coefficient */
	TARGET void Theta(POP *grp, double &f_ritland, double &f_loiselle, double &f_weir, double &t_ritland, double &t_loiselle, double &t_weir, int64 loc = -1);

	/* Write header row for individual statistics */
	TARGET static void IndividualStatisticsHeader(FILE *fout);

	/* Write result row for individual statistics */
	TARGET void PrintIndividualStatistics(FILE *fout);

	/* Write header row for ploidy inference, test */
	TARGET static void PloidyInferenceHeader(FILE *fout);

	/* Write result row for ploidy inference, test */
	TARGET void PrintPloidyInference(FILE *fout);

	/* Calculate the likelihood for ploidy inference */
	TARGET void PloidyInferlnL3(map<int64, SPF> &depth, int v, double f0, double f1, double f2, double &l0, double &l1, double &l2);

	/* Calculate the likelihood for ploidy inference */
	TARGET double PloidyInferlnL(map<int64, SPF> &depth, int v, double f);

	/* Infer the ploidy level from allelic depth distribution */
	TARGET void PloidyInference(int v, double &lnL, double &f, map<int64, SPF> &depth);

	/* average genotypic frequency at a diallelic locus given m copies of A */
	TARGET double AvgG(int v, int m, double f);

	/* Write header row for population assignment */
	TARGET static void AssignmentHeader(FILE *fout);

	/* Write result row for population assignment */
	TARGET void PrintAssignment(FILE *fout);

	/* Read and set genotype from bcf input */
	TARGET void AddBCFGenotype(int64 l, char *&gtstr, char *&gqstr, char *&dpstr, char *&adstr, int vlen, int asize, int gqlen, int dplen, int adlen, uint *&depth, TABLE<HASH, uint> &gfid, GENOTYPE *&gtab, ushort *&gatab, GENO_ITERATOR &iter);

	/* Read and set genotype from vcf input */
	TARGET void AddVCFGenotype(char *&line, int64 l, uint *&depth, TABLE<HASH, uint> &gfid, GENOTYPE *&gtab, ushort *&gatab, GENO_ITERATOR &iter);

	/* Get tag value */
	TARGET char *GetTagValue(char *re, int tagid);
};

/* Locus Diversity */
struct DIVERSITY
{
public:
	int64 l;								//Locus id

	/* Add */
	double bmaf;							//Minor allele freq of biallelic locus
	double ptype;							//Genotype rate
	double pval;							//P-val of genotype distribution test
	double he;								//Expected heterozygosity
	double ho;								//Observed heterozygosity
	double pic;								//Polymorphic information contents
	double ae;								//Effective number of alleles
	double I;								//Shannon¡¯s Information Index

	double fis;								//Inbreeding coefficient
	double g;								//G-statistic in HWE test
	double df;								//Degrees of freedom in HWE test

	int k;									//Number of alleles, int
	int n;									//Number of individuals, int
	int nhaplo;								//Number of allele copies, int
	int v2i;								//Number of within-allele pairs to weight Ho, int

	byte minploidy;							//Min ploidy
	byte maxploidy;							//Max ploidy
	bool varploidy;							//Does ploidy level varying among individuals
	bool unusued;							//Unusued byte for alignment

	/* Multiply */
	double NE1P;							//Exclusion rate without known parent
	double NE2P;							//Exclusion rate with known parent
	double NEPP;							//Exclusion rate for parent pair
	double NEID;							//Exclusion nonrelatives in identity test
	double NESI;							//Exclusion full-sibs in identity test

	/* Initialize */
	TARGET DIVERSITY();

	/* Calculate diveristy indices */
	TARGET void CalcDiversity(int64 _l);

	/* Write header to the result file */
	TARGET static void WriteHeader(FILE *f);

	/* Write diversity of a locus to the result file */
	TARGET void WriteLocus(FILE *f, const char *name);
};

/* Sum of locus diversity */
struct DIVSUM
{
public:
	/* Add */
	double k;								//Number of alleles
	double n;								//Number of individuals
	double nhaplo;							//Number of allele copies
	double bmaf;							//Minor allele freq of biallelic locus
	double ptype;							//Genotype rate
	double pval;							//P-val of genotype distribution test
	double he;								//Expected heterozygosity
	double ho;								//Observed heterozygosity
	double pic;								//Polymorphic information contents
	double ae;								//Effective number of alleles
	double I;								//Shannon¡¯s Information Index
	double fis;								//Inbreeding coefficient
	double minploidy;						//Min ploidy
	double maxploidy;						//Max ploidy

	/* Count */
	int kc;									//Number of alleles
	int nc;									//Number of individuals
	int nhaploc;							//Number of allele copies
	int bmafc;								//Minor allele freq of biallelic locus
	int ptypec;								//Genotype rate
	int pvalc;								//P-val of genotype distribution test
	int hec;								//Expected heterozygosity
	int hoc;								//Observed heterozygosity
	int picc;								//Polymorphic information contents
	int aec;								//Effective number of alleles
	int Ic;									//Shannon¡¯s Information Index
	int fisc;								//Inbreeding coefficient
	int minploidyc;							//Min ploidy
	int maxploidyc;							//Max ploidy

	/* Multiply */
	double NE1P;							//Exclusion rate without known parent
	double NE2P;							//Exclusion rate with known parent
	double NEPP;							//Exclusion rate for parent pair
	double NEID;							//Exclusion nonrelatives in identity test
	double NESI;							//Exclusion full-sibs in identity test

	LOCK lock;

	/* Initialize sum */
	TARGET void Init();

	/* Add loc to the sum */
	TARGET void Add(DIVERSITY& loc);

	/* Write sum */
	TARGET void Write(FILE *f, const char *name);
};

/* Statistics of locus */
struct LOCSTAT
{
public:
	ushort k;								//Number of alleles in target pop
	int nhaplo;								//Number of allele copies
	double a2;								//Sum(pi^2)
	double a3;								//Sum(pi^3)
	double a4;								//Sum(pi^4)
	double lm1;								//Mean difference between an allele and missing allele
	double lm2;								//Mean difference between two missing allele
};

struct POP
{
public:
	int id;									//Population id, will be rearrange using tree strucutre of total population
	char *name;								//Population name
	char **names;							//Individual names, used during parsing indtext
	int nind;								//Number of individuals
	int ind0id;								//Indid of the first individual, used to fast  genotype id iteration in this pop
	IND **inds;								//Individuals

	LOCSTAT *loc_stat;						//Locus diversity and statistics
	double *allelefreq;						//+ allele_freq_offset[l] is the allele frequency array at locus l
	ushort *genocount;						//+ genotype_count_offset[l] is the genotype count array at locus l, consider missing genotype

	int rid;								//Index of the region it belongs

	POP **vpop;								//Subpopulations
	int npop;								//Number of subpopulations
	int nhaplotypes;						//Number of haplotypes, used in amova homoploid method

	bool ispop;								//Is a population or a region

	/* Initialize */
	TARGET POP();

	/* Create a pop */
	TARGET POP(char *_name, char **_names, int _nind, int _regid, int _npop, int _id, bool _ispop);

	/* Get genotype count array */
	TARGET ushort *GetGenoCount(int64 l);

	/* Get genotype count of a genotype */
	TARGET ushort GetGenoCount(int64 l, int id);

	/* Get allele frequencies array */
	TARGET double *GetFreq(int64 l);

	/* Get allele frequency of an allele */
	TARGET double GetFreq(int64 l, int a);

	/* Uninitialize a pop */
	TARGET void Uninitialize();

	/* Uncllocate memory for locstat, allele frequency, genotype count */
	TARGET void UnAllocFreq();

	/* Allocate memory for locstat, allele frequency, genotype count */
	TARGET void AllocFreq();

	/* Move after filter monomorphic locus in bayesian clustering */
	TARGET void MoveFreq(LOCN *nafoffset, LOCN *ngcoffset);

	/* Clear memory for locstat, allele frequency, genotype count */
	TARGET void ClearFreqGcount();

	/* Calculate loc stat, allele frequency, genotype count */
	TARGET void CalcFreqGcount();
};

struct FST
{
public:
	/* Total value for all loci */
	double Nei1973T;						//Gst
	double Weir1984T;						//ANOVA
	double Hudson1992T;						//Mean allele difference
	double Slatkin1995T;					//Allele size difference
	double Hedrick2005T;					//G'st
	double Jost2008T;						//D
	double Huang2021homoT;					//AMOVA homoploid
	double Huang2021anisoT;					//AMOVA anisoploid

	/* Differentiation test by genotype distribution */
	double Genotype_GT;						//G-static
	int Genotype_DFT;						//degrees-of-freedom
	double Genotype_PT;						//P-val

	/* Differentiation test by allele distribution */
	double Allele_GT;						//G-static
	int Allele_DFT;							//Degrees-of-freedom
	double Allele_PT;						//P-val

	/* Value for each locus */
	double *Nei1973;						//Gst
	double *Weir1984;						//ANOVA
	double *Hudson1992;						//Mean allele difference
	double *Slatkin1995;					//Allele size difference
	double *Hedrick2005;					//G'st
	double *Jost2008;						//D
	double *Huang2021_homo;					//ANOVA homoploid
	double *Huang2021_aniso;				//ANOVA anisoploid

	/* Differentiation test by genotype distribution */
	double *Genotype_G;						//G-static
	int *Genotype_DF;						//degrees-of-freedom
	double *Genotype_P;						//P-val

	/* Differentiation test by allele distribution */
	double *Allele_G;						//G-static
	int *Allele_DF;							//Degrees-of-freedom
	double *Allele_P;						//P-val

	/* Fst estimator Warpper */
	TARGET static double FstEstimator(POP **grps, int n, int e, double *each, double *buf);

	/* Estimate Fst and test differentiation for two pops */
	TARGET void CalcFst(POP *a, POP *b);

	/* Estimate Fst and test differentiation for multiple pops */
	TARGET void CalcFst(POP **grps, int n);

	/* Uninitialize */
	TARGET void Uninitialize();

	/* Nei 1973 Fst estimator based on heterozgysotiy */
	TARGET static double Fst_Nei1973(POP **grps, int n, double *each, double *buf);

	/* Nei 1973 Fst estimator based on anova */
	TARGET static double Fst_Weir1984(POP **grps, int n, double *each);

	/* Nei 1973 Fst estimator based on mean allele difference */
	TARGET static double Fst_Hudson1992(POP **grps, int n, double *each, double *buf);

	/* Slatkin 1995 Fst estimator based on allele size */
	TARGET static double Fst_Slatkin1995(POP **grps, int n, double *each, double *buf);

	/* Hedrick 2005 G'st */
	TARGET static double Fst_Hedrick2005(POP **grps, int n, double *each, double *buf);

	/* Jost 2008 D */
	TARGET static double Fst_Jost2008(POP **grps, int n, double *each, double *buf);

	/* Huang 2021 Fst estimator based on multi-level amova */
	TARGET static double Fst_Huang2021_homo(POP **grps, int n, int layer, bool isiam, double *each);

	/* Huang 2021 Fst estimator based on multi-level amova */
	TARGET static double Fst_Huang2021_aniso(POP **grps, int n, int layer, bool isiam, bool sumss, double *each, double *buf);

	/* Write results file in column format */
	TARGET static void ColumnPrint(FILE *fout);

	/* Write results file in matrix format */
	TARGET static void MatrixPrint(FILE *fout, FST *Fst, int n, int type);
};

struct GDIST
{
public:
	double Nei1972;							//Genetic distance estimates
	double Cavalli1967;
	double Reynolds1983;
	double Nei1983;
	double Euclidean;
	double Goldstein1995;
	double Nei1974;
	double Roger1972;

	/* Converted from Fst by Slatkin's transform d = Fst/(1-Fst) */
	double Slatkin_Nei1973;					//Gst
	double Slatkin_Weir1984;				//ANOVA
	double Slatkin_Hudson1992;				//Mean allele difference
	double Slatkin_Slatkin1995;				//Allele size difference
	double Slatkin_Hedrick2005;				//G'st
	double Slatkin_Jost2008;				//D
	double Slatkin_Huang2021_homo;			//AMOVA homoploid
	double Slatkin_Huang2021_aniso;			//AMOVA anisoploid

	/* Converted from Fst by Reynolds's transform d = -ln(1 - Fst)*/
	double Reynolds_Nei1973;				//Gst
	double Reynolds_Weir1984;				//ANOVA
	double Reynolds_Hudson1992;				//Mean allele difference
	double Reynolds_Slatkin1995;			//Allele size difference
	double Reynolds_Hedrick2005;			//G'st
	double Reynolds_Jost2008;				//D
	double Reynolds_Huang2021_homo;			//AMOVA homoploid
	double Reynolds_Huang2021_aniso;		//AMOVA anisoploid

	/* Write column format header row for genetic distance estimation */
	TARGET static void ColumnPrintHeader();

	/* Write column format result row for genetic distance estimation */
	TARGET void ColumnPrintLine(int i, int j);

	/* Write matrix format header for genetic distance estimation */
	TARGET static void MatrixPrintMatrixHeader(int k, int n);

	/* Write matrix format row header for genetic distance estimation */
	TARGET static void MatrixPrintRowHeader(int k, int i);

	/* Write matrix format grid for genetic distance estimation */
	TARGET void MatrixPrintCell(int k);

	/* Use population/region allele frequency as the missing data */
	TARGET void GetMissingFreq(IND* a, GENOTYPE &gt, int64 l, double *pbuf);

	/* Calculate genetic distance between two individuals */
	TARGET void CalcGD(IND* a, IND* b, double *p1, double *p2);

	/* Calculate genetic distance between two populations/regions */
	TARGET void CalcGD(POP *a, POP *b, double *buf);
};

/* Used to permute vessels in AMOVA */
class VESSEL_ITERATOR
{
public:
	int relative_id[N_MAX_REG + 3];			//Index in its nesting vessles (from 0~nsubunits)
	int universal_id[N_MAX_REG + 3];		//Index among all vessles at this lay (from 0~nhaplo)
	VESSEL* trace[N_MAX_REG + 3];			//Current and its higher level vessels
	int lay;								//Current lay

	/* Go to start */
	TARGET void Rewind(int nlay);

	/* Copy from a reference*/
	TARGET void Copy(VESSEL_ITERATOR &ref, int nlay);

	/* Initialize */
	TARGET VESSEL_ITERATOR();

	/* Initialize */
	TARGET VESSEL_ITERATOR(int _lay, VESSEL &root, int nlay);

	/* Uninitialize */
	TARGET ~VESSEL_ITERATOR();

	/* Go to next vessel */
	TARGET void Next(int nlay);

	/* Get haplotype index to calculate genetic distance */
	TARGET int GetHapId();

	/* Get allele to fetch calculate distance */
	TARGET ushort GetAllele();

	/* Get subpopulation in print SS */
	TARGET POP *GetSubpop(int nlay, int tlay);

	/* Get individual in print SS */
	TARGET IND* GetInd(int nlay, int tlay);
};

/* Vessel of genes in AMOVA */
class VESSEL
{
public:
	VESSEL **subunits;						//Vessels nested within this vessel
	int *nhaplos;							//Number of haplotypes at locus l
	int *allelecount;						//KT elements, + allele_freq_offset[l] is the count of alleles at locus l in this vessel
	int nsubunits;							//Number of subunits
	int nhaplo;								//Number of haplotypes, homoploid model
	int hid;								//Haplotype id, -1 for non-allele vessel
	short lay;								//Level
	ushort allele;							//Allele for aniso model

	/* Uninitialize */
	TARGET ~VESSEL();

	/* Initialize */
	TARGET VESSEL();

	/* Deep copy a vessel */
	TARGET VESSEL(VESSEL &r);

	/* Create vessel from population */
	TARGET VESSEL(POP *s, int _lay, int &_hid, int64 loc, int method);

	/* Create vessel from individual */
	TARGET VESSEL(IND* s, int _lay, int &_hid, int64 loc, int method);

	/* Create vessel from haplotype */
	TARGET VESSEL(int _lay, int &_hid, ushort _allele);

	/* Get the allele count array at locus l*/
	TARGET int *GetAlleleCount(int l);

	/* Save all vellels in level fa into an array */
	TARGET void GetVessels(VESSEL **vs, int &nvessels, int fa);

	/* Replace with shuffled vessels */
	TARGET int Replace(VESSEL **vs, int &nvessels, int fa, int method);

	/* Shuffle fa level vessels among fb level vessels */
	TARGET void Shuffle(RNG &rng, int fa, int fb, int method, VESSEL **buf);

	/* Calculate matrix C for maximum-likelihood method */
	TARGET void GetCML(double *C, int64 l, int *tid, double *tw, int Nh, int nlay, double **W);

	/* Initialize matrix C, the coefficient matrix of S = CV */
	TARGET void InitC(double *C, int *tid, int Nh, int nlay);

	/* Calculate matrix C */
	TARGET void GetC(double *tw, int *tid, double *C, int nlay, double **W, int64 l);

	/* Count number of vessels in each hierarchy */
	TARGET void CountVessels(int *count);

	/* Initialize W, W[lay][tid[lay]] = 1 / nhaplo */
	TARGET void InitW(MEMORY &mem, double **&W, int nlay);

	/* Calculate SS for homoploid method */
	TARGET void GetSSHomo(double *SS, double *gd, int Nh, double **W, int nlay, VESSEL_ITERATOR &ve1, VESSEL_ITERATOR &ve2);

	/* Calculate SS for anisoploid method */
	TARGET void GetSSAniso(ushort *hap_bucket, double *SS, bool isiam, int nh, int64 l, double **W, int nlay, VESSEL_ITERATOR &ve1, VESSEL_ITERATOR &ve2);

	/* Calculate SS for anisoploid method */
	TARGET void GetSSAniso(double *SS, bool isiam, int nh, int64 l, double **W, int nlay, VESSEL_ITERATOR &ve1, VESSEL_ITERATOR &ve2);

	/* Calculate variance component matrix V */
	TARGET static void GetV(double *C, double *SS, double *&V, int nlay);

	/* Calculate F-statistics */
	TARGET static void GetF(double *V, double *F, double *vs, int nlay);

	/* Calculate F-statistics */
	TARGET static void GetF(double *Fi, double *F, int nlay);
};

struct AMOVA
{
	//One time
	double *V;								//Variance components at each level, nlay elements
	double *F;								//F-statistics, nlay*nlay elements
	int *G;									//Number of permuations with F'>F, nlay*nlay elements
	int *E;									//Number of permuations with F'=F, nlay*nlay elements

	double *EF2;							//Mean permuated squared F', nlay*nlay elements
	double *EF;								//Mean permuated F', nlay*nlay elements

	double *SS;								//SS within each level, nlay elements
	double *DF;								//Degrees-of-freedom for each level, nlay elements
	double **SSW;							//SS within each vessel
	int *nSSW;								//Size of each SSW, nlay elements

	int nlay;								//Number of levels
	int method;								//Method: 1 homoploid, 2 anisoploid or 4 likelihood
	int npermed;							//Number of permuated iterations
	int Lind;								//Consider individual level?

	/* Initialize */
	TARGET AMOVA();

	/* Extract dummy haplotype for homoploid method */
	TARGET void GetHaplotype(POP *tpop, ushort *bucket);

	/* Perform AMOVA using homoploid method */
	TARGET void CalcAMOVA_homo();

	/* Perform AMOVA using anisoploid method */
	TARGET void CalcAMOVA_aniso();

	/* Calculate likelihood for permuated data */
	TARGET static double Likelihood(CPOINT &xx, void **Param);

	/* Perform AMOVA using maximum-likelihood method */
	TARGET void CalcAMOVA_ml();

	/* Destructor */
	TARGET ~AMOVA();

	/* Write results */
	TARGET void PrintAMOVA(FILE *fout);
};

struct Huang2015ENTRY
{
public:
	int ibs;								//IBS model index
	uint64 pattern;							//Genotype pair pattern
};

struct RELATEDNESS
{
public:
	double Lynch1999;						//Estimates for various relatedness estimators
	double Wang2002;
	double Thomas2010;
	double Li1993;
	double Queller1989;
	double Huang2016A;
	double Huang2016B;
	double Milligan2003;
	double Anderson2007;
	double Huang2014;
	double Huang2015;
	double Ritland1996m;
	double Loiselle1995m;
	double Ritland1996;
	double Loiselle1995;
	double Weir1996;
	int ABtype;								//Number of loci genotyped in both individuals
	int Atype;								//Number of loci genotyped in individual A
	int Btype;								//Number of loci genotyped in individual B

	/* Write header row for relatedness estimation */
	TARGET static void ColumnPrintHeader();

	/* Write result row for relatedness estimation */
	TARGET void ColumnPrintLine(int i, int j);

	/* Write matrix format header for relatedness estimation */
	TARGET static void MatrixPrintMatrixHeader(int k, int n);

	/* Write matrix format row header for relatedness estimation */
	TARGET static void MatrixPrintRowHeader(int k, int i);

	/* Write matrix format grid for relatedness estimation */
	TARGET void MatrixPrintCell(int k);

	/* Calculate relatedness coefficient */
	TARGET void CalcRelatedness(IND* x, IND* y);

	/* Relatedness estimator Warpper */
	TARGET static double RelatednessEstimator(int k, IND* x, IND* y);

	/* Lynch 1999 relatedness estimator */
	TARGET static double R_Lynch1999(IND* x, IND* y);

	/* Wang 2002 relatedness estimator */
	TARGET static double R_Wang2002(IND* x, IND* y);

	/* Thomas 2010 relatedness estimator */
	TARGET static double R_Thomas2010(IND* x, IND* y);

	/* Li 1993 relatedness estimator */
	TARGET static double R_Li1993(IND* x, IND* y);

	/* Queller 1989 relatedness estimator */
	TARGET static double R_Queller1989(IND* x, IND* y);

	/* Huang 2016 relatedness estimator A */
	TARGET static double R_Huang2016A(IND* x, IND* y);

	/* Huang 2016 relatedness estimator B */
	TARGET static double R_Huang2016B(IND* x, IND* y);

	/* Initialize Anderson 2007 relatedness estimator */
	TARGET static void R_AndersonInitialize(IND* x, IND* y);

	/* Milligan 2003 relatedness estimator */
	TARGET static double R_Milligan2003(IND* x, IND* y);

	/* Anderson 2007 relatedness estimator */
	TARGET static double R_Anderson2007(IND* x, IND* y, bool confine);

	/* Calculate Anderson 2007 likelihood */
	TARGET static double L_Anderson(CPOINT &x, void **unusued);

	/* Ritland 1996 kinship estimator, convert into relatedness */
	TARGET static double R_Ritland1996(IND* x, IND* y, bool iscorrect, bool mulv);

	/* Loiselle 1995 kinship estimator, convert into relatedness */
	TARGET static double R_Loiselle1995(IND* x, IND* y, bool iscorrect, bool mulv);

	/* Weir 1996 kinship estimator, convert into relatedness */
	TARGET static double R_Weir1996(IND* x, IND* y, bool mulv);

	/* Huang 2014 relatedness estimator */
	TARGET static double R_Huang2014(IND* x, IND* y);

	/* Huang 2014 relatedness estimator : similarity index */
	TARGET static double S_Index(int *c, int *d, int ploidyx, int ploidyy);

	/* Huang 2014 relatedness estimator : get genotype pattern for reference individual */
	TARGET static int GetRefMode(int *a, int ploidy);

	/* Huang 2014 relatedness estimator : calculate relatedness */
	TARGET static double HuangMoment(GENOTYPE &gx, GENOTYPE &gy, int64 l, double &weight);

	/* Initialize Huang 2015 relatedness estimator */
	TARGET static void Huang2015_Initialize();

	/* Uninitialize Huang 2015 relatedness estimator */
	TARGET static void Huang2015_Uninitialize();

	/* Huang 2015 relatedness estimator */
	TARGET static double R_Huang2015(IND* x, IND* y);

	/* Calculate Huang 2015 likelihood */
	TARGET static double L_Huang2015(CPOINT &xx, void **unusued);

	/* Huang 2015 likelihood estimator: Match genotype-pair pattern and assign alleles */
	TARGET static void Huang2015_MatchAllele(int64 pattern, int *gx, int *gy, int *alleles, int p);
};

struct KINSHIP
{
public:
	double Ritland1996;						//Estimates for various kinship estimators
	double Loiselle1995;
	double Weir1996;
	int ABtype;								//Number of loci genotyped in both individuals
	int Atype;								//Number of loci genotyped in individual A
	int Btype;								//Number of loci genotyped in individual B

	TARGET static void ColumnPrintHeader();

	TARGET void ColumnPrintLine(int i, int j);

	TARGET static void MatrixPrintMatrixHeader(int k, int n);

	TARGET static void MatrixPrintRowHeader(int k, int i);

	TARGET void MatrixPrintCell(int k);

	TARGET void CalcKinship(IND* a, IND* b);
};

struct PCOA
{
	double *D;								//Euclidean distance matrix
	int N;									//Number of objects
	double *U;								//Eigen-vector matrix
	double *V;								//Eigen-value matrix
	double Vt;								//Total variance
	int p;									//Number of expected coordinates
	int estimator;							//PCoA or PCA
	int type;								//Type of objects: 1 individual, 2 pop, 3 reg
	int maxp;								//Number of extracted coordinates

	/* Do nothing */
	TARGET PCOA();

	/* Destructor */
	TARGET ~PCOA();

	/* Perform PCoA */
	TARGET int CalcPCoA(int _maxp);

	/* Print PCOA */
	TARGET void PrintPCoA(FILE *fout, double *d, int n, int _est, int _type);
};

struct HCLUSTER
{
public:
	bool isend;								//Is a leaf node
	char *endname;							//Object name
	double x;								//Node coordinate x
	double y;								//Node coordinate y
	uint *id;								//Objects ids
	uint  idlen;							//Objects size
	HCLUSTER *left;							//Left node
	HCLUSTER *right;						//Right node
};

class HCLUSTERING
{
public:

	int method;								//Clustering method
	double *dori; 							//Original distance matrix
	double *dcur;							//New distance matrix
	double *dnew;							//Current distance matrix
	int nori;								//Dimension of dori
	int ncur;								//Dimension of dcur
	LIST<HCLUSTER*> node;					//Nodes
	MEMORY *memory;							//Memory class

	TARGET HCLUSTERING(double *d, IND **obj, int n, int m, MEMORY *_memory);

	TARGET HCLUSTERING(double *d, POP **obj, int n, int m, MEMORY *_memory);

	TARGET ~HCLUSTERING();

	TARGET void CalcClustering();

	TARGET double FindMinIdx(int &a, int &b);

	TARGET void ReduceMatrix(int _a, int _b);

	TARGET void PrintClustering(FILE *fout, HCLUSTER *c = NULL, double cy = 0);
};

struct SCLUSTER
{
public:
	double *bucket;							//Allele frequency bucket, KT elements

	/* Get allele frequency array */
	TARGET double *GetFreq(int64 l);

	/* Get allele frequency */
	TARGET double GetFreq(int64 l, int allele);

	/* Set allele frequency pointer */
	TARGET void SetFreq(double *_bucket);
};

class STRUCTURE
{
public:
	/* Parameter */
	int S;									//Number of sampling locations
	int K;									//Number of clusters
	int N;									//Number of individuals
	int64 L;								//Number of loci

	/* Model */
	bool admix;								//Admix model
	bool locpriori;							//Locpriori model
	bool fmodel;							//F model

	/* MCMC */
	int nburnin;							//Number of burnin iterations
	int nreps;								//Number of iteration after burnin
	int nthinning;							//Sampling interval
	int nruns;								//Number of independent runs for each K
	int nadmburnin;							//Number of admburnin iterations

	/* Misc */
	double lambda;							//Dirichlet parameter to update the allele frequencies 
	double stdlambda;						//Standard deviation of new lambda
	double maxlambda;						//Maximum of new lambda
	bool inferlambda;						//Updated lambda in each iteration
	bool difflambda;						//Use separate lambda for each cluster
	int diversity;							//Output diversity

	/* Admix */
	double alpha;							//The initial alpha, the priori Dirichlet parameter of admixture proportions Q
	bool inferalpha;						//Update alpha in ADMIX model
	bool diffalpha;							//Use separate alpha for each cluster
	bool uniformalpha;						//Priori distribution for alpha, uniform or gamma alpha
	double stdalpha;						//Standard deviation of uniform priori distribution of alpha
	double maxalpha;						//Maximum of uniform priori distribution of alpha
	double alphapriora;						//One gamma priori distribution parameter
	double alphapriorb;						//The other gamma priori distribution parameter
	int metrofreq;							//Frequency of Metropolis-Hastings update of admixture proportions Q

	/* Locpriori model */
	double r;								//Initial value of r, where r evaluates the informativeness of data for the sampling location
	double maxr;							//Maximum of new r, where r evaluates the informativeness of data for the sampling location
	double epsr;							//Max step value of new r
	double epseta;							//Max step value of new eta for the non-ADMIXTURE model
	double epsgamma;						//Max step value of new gamma for the non-ADMIXTURE model

	/* F model */
	double pmeanf;							//Priori mean F
	double pstdf;							//Priori standard deviation of F
	double stdf;							//Standard deviation of new F
	bool fsame;								//Use the same F in all clusters

	
	double *buf;							//Buffer with length K
	double *bufb;							//Buffer with length K
	double **buf2;							//Buffer with length K

			//Fixed
	SCLUSTER* cluster;						//Clusters, K elements
	SCLUSTER* clusterb;						//Additioanl clusters to write temp allele freq, K elements
	double *Base;							//(2 * K + 3) * KT elements, cluster[K*KT], clusterb[K*KT], PA[KT], PA1[KT], PA2[KT] 

	double *Lambda;							//K, Priori dirichlet parameter for each cluster
	double *Z_cumu;							//N*K, Z_cumu added by Mi in each iteration
	ushort **Z;								//[N] * vt, current culster of allel copy
	ushort **ZZ;							//Loop variable for Z
	double *Mi;								//N*K, Mi[i,k] is the number of allele copies of individuali assigned to cluster k
	double *Q;								//N*K, Q[i,k] is the priori proportion of indidivual i's gene from cluster k
	int *Ni;								//K*KT, Ni[k, k2] is the number of allele copies in each cluster
	double *Alpha;							//K, the priori Dirichlet parameter of admixture proportions Q for each cluster

			//F model, epsilon, ancestral population
	double *f;								//K,   (1-F)/F
	double *F;								//K,    Fst
	SCLUSTER PA, PA1, PA2;

			//LocPriori model
	double *Eta;							//K, eta for each cluster
	double *Gamma;							//S*K
	double *Di;								//S*K, number of individual from s into k

	double *AlphaLocal;						//S*K
	double *SumAlpha;						//S

	RNG rng;								//Random number generator
	int m;									//Current itation id
	bool binaryq;							//Is q only takes from 0 and 1
	SRUNINFO *par2;							//Structure parameter
	int nr;									//Number of recorded records

	double *rout;							//Output buffer
	int rlen;								//Rout size

	int *kdis;								//kdis[k] = number of loci has k alleles


	/* Set all bits to 0 */
	TARGET STRUCTURE();

	/* Write results for a run */
	TARGET void PrintStructure();

	/* Write results summary for all runs */
	TARGET static void PrintSummary(SRUNINFO *sp, int len);

	/* Copy parameters */
	TARGET void ReadPar(SRUNINFO *_par2);

	/* Initialize MCMC */
	TARGET void InitFix();

	/* Initialize MCMC for admix model */
	TARGET void InitAdmix();

	/* Initialize MCMC for locpriori model */
	TARGET void InitLocPriori();

	/* Initialize MCMC for F model */
	TARGET void InitFmodel();

	/* Update allele frequency for all clusters */
	TARGET void UpdateP();

	/* Update a priori ancetral proportion for non-admix model */
	TARGET void UpdateQNoAdmix();

	/* Update a priori ancetral proportion for admix model */
	TARGET void UpdateQAdmix();

	/* Update a priori ancetral proportion by Metropolis-Hastings for admix model*/
	TARGET void UpdateQMetro();

	/* Update a priori ancetral proportion */
	TARGET void UpdateQ();

	/* Update locpriori parameters */
	TARGET void UpdateLocPriori();

	/* Update ancestral proportion for each allele or for each individual */
	TARGET void UpdateZ();

	/* Update Dirichlet parameter alpha (to draw admixture proportion Q) in the admix model */
	TARGET void UpdateAlpha();

	/* Update Dirichlet parameter lambda (to draw allele frequency) */
	TARGET void UpdateLambda();

	/* Update allele frequency of ancestral population for the F model */
	TARGET void UpdatePA();

	/* Update population-specific Fst for the F model */
	TARGET void UpdateF();

	/* Finalize records */
	TARGET void Arrange();

	/* Record updated MCMC parameters */
	TARGET void Record();

	/* Free memory */
	TARGET void Uninit();

	/* Perform MCMC */
	TARGET void MCMC();
};

struct HAPLO_DUMMY_HAPLOTYPE
{
public:
	ushort alleleid;						//Dummy allele index
	ushort *alleles;						//True allele array in [st-ed] variant, alloc in memory class

	/* Do nothing */
	TARGET HAPLO_DUMMY_HAPLOTYPE();

	/* Extract the ith haplotype from an individual */
	TARGET void ExtractHaplotype(byte vi, IND* ti, int64 st, int64 ed, int nvar, ushort aid, MEMORY &haplo_memory);

	/* Print information for an extracted locus */
	TARGET void PrintHaplotype(FILE *f1, int64 st, int64 ed);
};

struct HAPLO_DUMMY_LOCUS
{
	int64 chrid;							//Chrom id
	int64 stpos;							//Pos of start locus
	int64 st;								//Start locus id
	int64 ed;								//End locus id
	int nvar;								//Number of variants
	int hsize;								//Number of haplotypes
	int gsize;								//Number of genotypes
};


#ifndef _VCF

/* Initialize */
TARGET void Initialize();

/* UnInitialize */
TARGET void UnInitialize();

/* Close input files */
TARGET void CloseInput();

/* Assign individual indid to population popid */
TARGET void AssignIndSub(int indid, int popid, int &namelen, int &nind2);

/* Assign individuals to their populations */
TARGET void AssignInds();

/* Load data from input files */
TARGET void LoadFile();

/* Check anisoploid in haplotype extraction */
THREADH(CheckAnisoploid);

/* Perform haplotype extraction */
TARGET void CalcHaplotype();

/* Applying individual and diversity filters */
TARGET void ApplyFilter();

/* Recursive set id, pops, inds for each region */
TARGET void SetVReg(int rl, int i);

/* Sort individuals by population index to rinds array */
TARGET void AssignVInds();

/* Calculate individual min and max ploidy, and sum ploidy levels */
TARGET void AssignPloidy();

/* Calculate allele frequencies for each population and region for further use */
TARGET void CalcFreq();

/* Convert into other genotype formats */
TARGET void ConvertFile();

/* Calculate genetic diveristy indices */
TARGET void CalcDiversity();

/* Calculate individual statistics */
TARGET void CalcIndstat();

/* Calculate genetic differentiation */
TARGET void CalcDiff();

/* Calculate genetic distance */
TARGET void CalcDist();

/* Calculate analysis of molecular variance */
TARGET void CalcAMOVA();

/* Calculate population assignment */
TARGET void CalcAssignment();

/* Calculate relatedness coefficient */
TARGET void CalcRelatedness();

/* Calculate kinship coefficient */
TARGET void CalcKinship();

/* Calculate principal coordinate analysis */
TARGET void CalcPCOA();

/* Calculate hierarchical clustering */
TARGET void CalcClustering();

/* Calculate spatical pattern */
TARGET void CalcSPA();

/* Calculate bayesian clustering */
TARGET void CalcBayesian();

/* Calculate individual ploidy inference */
TARGET void CalcPloidyInference();

/* Calculate various analyses */
TARGET void Calculate();

/* Create genotype table for non-vcf/bcf input files */
TARGET void CreateGenoIndexTable(GENO_ITERATOR *iter = NULL);

/* 0. Loading thread functions */

/* load from Genepop input format */
THREADH(LoadGenepop);

/* load from Spagedi input format */
THREADH(LoadSpagedi);

/* load from Cervus input format */
THREADH(LoadCervus);

/* load from Arlequin input format */
THREADH(LoadArlequin);

/* load from Structure input format */
THREADH(LoadStructure);

/* load from PolyRelatedness input format */
THREADH(LoadPolyRelatedness);

/* load from PolyGene input format */
THREADH(LoadPolyGene);

/* Indexing alleles for non-vcf input, with allele identifier being the size */
TARGET void IndexAlleleLength();

/* Check input files rows and columns and count the number of variants */
THREADH(GetBCFLines);

/* Process lines from memory */
THREADH(LoadBCF);

/* Read lines from BCF file */
THREADH(LoadBCFGuard);

/* Check input files rows and columns and count the number of variants */
THREADH(GetVCFLines);

/* Process lines from memory */
THREADH(LoadVCF);

/* Read lines from VCF file */
THREADH(LoadVCFGuard);

/* Calculate allele frequencies for each population and region */
THREADH(CalcAlleleFreq);

/* 1. Filter thread functions */

/* Marker individual filtered or not */
THREADH(MarkerIndividual);

/* Marker locus filtered or not */
THREADH(MarkerDiversity);

/* Remove individual fail to pass filter */
THREADH(RemoveIndividual);

/* Remove locus fail to pass filter */
THREADH(RemoveLocus);

/* 3. Haplotype extraction */

/* Quick sort locus by contig and position */
TARGET void QSLocus(int64 left, int64 right);

/* Quick sort locus in a contig */
THREADH(QSWorker);

/* Quick sort extracted locus by contig and position */
TARGET void QSHapLocus(int64 left, int64 right);

/* Quick sort extracted locus in a contig */
THREADH(QSHapWorker);

/* Get number of alleles and genotypes at a dummy locus */
TARGET double GetDummyK(uint64 st, uint64 ed, TABLE<HASH, ushort> &hfidx, TABLE<HASH, ushort> &gfidx);

/* Create locus for haplotype extraction */
THREADH(CreateHaplotypeLocus);

/* Output locus for haplotype extraction */
THREADH(WriteHaplotypeLocus);

/* Find the optimal PES model for each locus */
THREADH(GetLocusPESModel);

/* Calculate individual min and max ploidy, and sum ploidy levels */
THREADH(AssignPloidyThread);

/* 4. Conversion thread functions */

/* Write convert genepop genotypes in a guard thread */
THREADH(ConvertGenepopGuard);

/* Write convert arlequin genotypes in a guard thread */
THREADH(ConvertArlequinGuard);

/* Write convert genotypes in a guard thread */
THREADH(ConvertGuard);

/* Convert genotype string */
TARGET void PrepareGenotypeString(int format);

/* Convert into genepop format */
TARGET void ConvertGenepop(int ntot, bool &isfirst);

/* Convert individual genotypes into genepop format in multiple threads */
THREADH(ConvertGenepopInd);

/* Convert into spagedi format */
TARGET void ConvertSpagedi(int ntot, bool &isfirst);

/* Convert individual genotypes into spagedi format in multiple threads */
THREADH(ConvertSpagediInd);

/* Convert into cervus format */
TARGET void ConvertCervus(int ntot, bool &isfirst);

/* Convert individual genotypes into cervus format in multiple threads */
THREADH(ConvertCervusInd);

/* Convert into arlequin format */
TARGET void ConvertArlequin(int ntot, bool &isfirst);

/* Convert individual genotypes into arlequin format in multiple threads */
THREADH(ConvertArlequinInd);

/* Convert into structure format */
TARGET void ConvertStructure(int ntot, bool &isfirst);

/* Convert individual genotypes into structure format in multiple threads */
THREADH(ConvertStructureInd);

/* Convert into polygene format */
TARGET void ConvertPolygene(int ntot, bool &isfirst);

/* Convert individual genotypes into polygene format in multiple threads */
THREADH(ConvertPolygeneInd);

/* Convert into polyrelatedness format */
TARGET void ConvertPolyrelatedness(int ntot, bool &isfirst);

/* Convert individual genotypes into polyrelatedness format in multiple threads */
THREADH(ConvertPolyrelatednessInd);

/* 5. Diversity thread functions */

/* Add and write genetic diversity */
THREADH(DiversityGuard);

/* Calculate genetic diversity using multiple threads */
THREADH(DiversityThread);

/* 6. Individual Statistics thread function */

/* Write individual statistics */
THREADH(IndividualStatisticsGuard);

/* Calculate individual statistics using multiple threads */
THREADH(IndividualStatisticsThread);

/* 7. Genetic Differentiation thread function */

/* Calculate genetic differentiation using multiple threads */
THREADH(GeneticDifferentiationThread);

/* 8. Genetic Distance thread function */

/* Write column format genetic distance results in a guard thread */
THREADH(GeneticDistanceGuard1);

/* Write matrix format genetic distance results in a guard thread */
THREADH(GeneticDistanceGuard2);

/* Calculate genetic distance using multiple threads */
THREADH(GeneticDistanceThread);

/* 9. AMOVA thread function */

/* Calculate AMOVA using multiple threads */
THREADH(AMOVAThread);

/* 10. Population assignment thread function */

/* Calculate population assignment using multiple threads */
THREADH(PopulationAssignmentThread);

/* 11. Relatedness thread function */

/* Write column format relatedness coefficient results in a guard thread */
THREADH(RelatednessGuard1);

/* Write matrix format relatedness coefficient results in a guard thread */
THREADH(RelatednessGuard2);

/* Calculate relatedness coefficient using multiple threads */
THREADH(RelatednessThread);

/* 12. Kinship thread function */

/* Write column format kinship coefficient results in a guard thread */
THREADH(KinshipGuard1);

/* Write matrix format kinship coefficient results in a guard thread */
THREADH(KinshipGuard2);

/* Calculate kinship coefficient using multiple threads */
THREADH(KinshipThread);

/* 13. Principal Coordinate Analysis thread function */
/* 14. Hierarchical clustering */

/* Calculate genetic distance for PCoA and Hierarchical clustering using multiple threads */
THREADH(PCoAClusteringThread);


/* Test. SPA */
TARGET double SPA_Fij(double *a, int i);

TARGET double SPA_Likelihood(CPOINT &xx, void **Param);

TARGET double SPA_Likelihood(uint64 l, double *a);

TARGET double SPA_Hessian(uint64 l, double *G, double *H, double *a, double *a2);

TARGET void SPA_Locus(uint64 l, double *x, double *f, double *a, double *a2, double *am, double *g, double *h);

/* Calculate spatial pattern using multiple threads */
THREADH(SPAThread);

/* 15. Bayesian clustering */

/* Calculate Bayesian clustering using multiple threads */
THREADH(StructureThread);

/* 16. Ploidy Inference thread function */

/* Calculate ploidy inference using multiple threads */
THREADH(PloidyInferenceThread);
#endif


/* Thread-specific variables */

extern _thread MEMORY *amova_memory;				//Memory class for amova vessels
extern _thread double *Anderson2007_Coef;			//Anderson 2007 relatedness estimator coefficients
extern _thread double *Huang2015_Coef;				//Huang 2015 relatedness estimator coefficients
extern _thread int threadid;						//Thread index
extern _thread bool *spa_valid;
extern _thread int spa_tn;

/* Functions */
extern TABLE<int, Huang2015ENTRY> *Huang2015_maps;	//Huang2015_maps[ploidylevel][hash] is a entry saves the ibs modex index and genotype pair pattern

/* Hash */
extern uint *cryptTable;							//Crypt table for calculate hash for genotypes

/* Functions */
extern int NBUF;									//CALC_THREAD_BUFFER * g_nthread_val;
extern MEMORY *individual_memory;					//Individual memory class
extern MEMORY *locus_memory;						//Locus memory class
extern MEMORY *nvcf_memory;							//Locus memory for first round counting ngeno
extern TABLE<HASH, uint> *nvcf_gfid;				//Hash table for counting ngeno 
extern MEMORY *conversion_memory;					//Memory class for conversion_string
extern MEMORY *conversion_memory2;					//Memory class for genotype_string and convert_buf
extern LIST<char*> *conversion_string;				//Genotype string for each genotype for file conversion
extern MEMORY *gd_memory;							//Individual genetic distance memory class

extern TABLE<HASH, INDGD> *gdtab;					//Hash table saves the genetic distance between genotypes
extern shared_mutex* gdlock;						//gdtab lock2
extern double *amova_gd;							//Genetic distance matrix used in AMOVA

extern byte *genotype_bucket;						//Genotype index data bucket
extern uint64 genotype_size;						//Size of bucket in bytes
extern uint64 genotype_coffset;						//Size of used bucket memory
extern OFFSET *genotype_offset;						//Genotype index data at each locus

extern ushort missing_array[N_MAX_PLOIDY];			//Allele array of the missing genotypes
extern GENOTYPE missing_genotype[N_MAX_PLOIDY + 1];	//Missing genotype at different ploidy level
extern HASH missing_hash[N_MAX_PLOIDY + 1];			//Hash of missing genotype

/* Allelic depth, test*/
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
extern IND **rinds;									//Rearranged individuals (by enumerate tree structure of total population)

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

extern double *spa_x;								//Spatial pattern, test, don't use, n * (dim + 1)
extern int spa_n;									//Spatial pattern, test, don't use, n
extern int spa_np;									//Spatial pattern, test, don't use, n

#pragma pack(pop)