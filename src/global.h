/* Constants and Global variables */

#pragma once
#include "vcfpop.h"

#ifdef _WIN64
	#define _thread __declspec(thread)
#else 
	#define _thread __thread
#endif

//#define genoid bits2[3]
#define VERSION						("1.01")
#define DATE						("2022-04-12")

#define N_FUNC						(17)						//Number of main options
#define N_MAX_OPTION				(35)						//Maximum options
#define N_MAX_REG					(10)						//Maximum region levels
#define N_MAX_GDTAB					(20)						//Maximum number of genotypes a locus have to buffer individual genetic distance
#define N_CONVERTER					(7)							//Number of file converters
#define N_GD_ESTIMATOR				(24)						//Number of genetic distance estimators
#define N_INDGD						(11)						//Number of genetic distances in individual gd structure
#define N_CLUSTER_METHOD			(7)							//Number of hierarchical clustering estimators
#define N_FST_ESTIMATOR				(8)							//Number of Fst estimators
#define N_RELATEDNESS_ESTIMATOR 	(16)						//Number of relatedness estimators
#define N_KINSHIP_ESTIMATOR			(3)							//Number of kinship estimators
#define N_DRE_MODEL					(4)							//Number of double-reductio models
#define N_DRE_MODELT				(104)						//Number of double-reductio models with additional 1% gridents for PES models
#define N_MAX_PLOIDY				(10)						//Maximum ploidy level supported
#define N_MAXGAMMALN				(1024)						//Number of gammaln bufferred 
#define N_PATTERN_END				(139)						//Number of patterns, for pid>=139, pid-139 is ploidy level

#define PI							(3.14159265358979323846)
#define LIKELIHOOD_TERM				(1e-8)						//Threshold of likelihood difference to terminate iteration
#define MAX_ITER_DOWNHILL			(1200)						//Max iterations for down-hill simplex
#define MAX_ITER_SVD				(60)						//Max iterations for SVD decompsition
#define EPSILON_SVD					(1e-30)						//Minimum positive value for SVD decompsition
#define DOUBLE_UNDERFLOW			(1e-100)					//Lower bound of real number
#define DOUBLE_OVERFLOW				(1e100)						//Upper bound of real number
#define MIN_FREQ					(1e-10)						//Minimum allele freq to avoid being zero in unify

#define FILE_NAME_LEN				(8192)						//Buffer size for filename
#define IND_NAME_LEN				(256)						//Buffer size for individual identifier
#define CALC_THREAD_BUFFER			(1024)						//Buffer size for kinship, relatedness, distance threads
#define LINE_BUFFER					(1048576)					//Buffer size for a line
#define PATH_LEN					(4096)						//Buffer size for file name

#define SLEEP_TIME_TINY				(5)							//A small sleep time to wait other threads

/* Global Variables */

extern LOCK glock1, glock2;
extern double NA;
extern bool ALLELE_IDENTIFIER;									//No allele name string for non-vcf/bcf input and vcf/bcf that performed haplotype extraction
extern int SIMD_TYPE;											//Single-instruction-multiple-data instruction used
extern int GDIST_METHOD;										//Current GD method, 1 genetic distance, 2 pcoa or 3 hierarchical clustering
extern char EXEDIR[PATH_LEN];									//Executable directory
extern char CURDIR[PATH_LEN];									//Current directory
extern double ALPHA[N_DRE_MODELT + 1][N_MAX_PLOIDY + 1][3];		//Double reduction rates
extern double BINOMIAL[N_MAX_PLOIDY + 1][N_MAX_PLOIDY + 1];		//Binomial coefficients
extern double GAMMALN[N_MAXGAMMALN + 1];						//Logarithm of gamma functions
extern double LOGFRAC[N_MAX_PLOIDY + 1][N_MAX_PLOIDY + 1];		//Logarithm of fractions

extern atomic<int64> PROGRESS_VALUE;							//Progress value 
extern atomic<int64> PROGRESS_VALUE2;							//Progress value 2 
extern atomic<int64> PROGRESS_VALUE3;							//Progress value 3
extern int64 PROGRESS_TOTAL;									//Total tasks in this function
extern int64 PROGRESS_CSTART;									//Start value of this batch
extern int64 PROGRESS_CEND;										//End value of this batch
extern int64 PROGRESS_NOUTPUTED;								//Number of outputed characters
extern int64 PROGRESS_NOUTPUTED2;								//Number of previously outputed characters

extern const char *GD_ESTIMATOR[];									//Genetic distance estimator names
extern const char *CLUSTER_METHOD[]; 									//Hierarchical method names
extern const char *FST_ESTIMATOR[]; 									//Fst estimator names
extern const char *RELATEDNESS_ESTIMATOR[]; 							//Relatedness estimator names
extern const char *KINSHIP_ESTIMATOR[]; 								//Kinship estimator names
extern const char *DRE_MODEL[];  										//Double-reduction model names

/* Virtual Memory */

extern bool VIRTUAL_MEMORY; 									//Allocate continous memory at a fixed address to avoid frequent realloc and move
extern bool BIG_FILE;											//use > 10Gib genotype id memory

extern int64 V_BASE_GENOTYPE;									//Allocate at 0x080000000000 to save genotype index
extern int64 *V_ALLOC_LEN_GENOTYPE;								//Size of each block
extern int64 V_ALLOC_SIZE_GENOTYPE;								//Number of blocks
extern int64 V_ALLOC_CSIZE_GENOTYPE;							//Current block index

extern int64 V_BASE_ALLELEDEPTH;								//Allocate at 0x100000000000 to save genotype index
extern int64 *V_ALLOC_LEN_ALLELEDEPTH;							//Size of each block
extern int64 V_ALLOC_SIZE_ALLELEDEPTH;							//Number of blocks
extern int64 V_ALLOC_CSIZE_ALLELEDEPTH;							//Current block index


/* Temple and results file */

extern char *FRES_NAME;											//Result file name
extern char *FRES_BUF;											//Result buffer
extern tm *FRES_TIME;											//Result time
extern FILE *FRES;												//Result file pointers
extern FILE **TEMP_FILES;										//Temporatory file pointers
extern char **TEMP_BUFS;										//Temporatory buffer
extern char **TEMP_NAMES;										//Temporatory file names


/* Parameters */

/* Parameter file */
extern bool p_b;								extern char *p_val;									extern char *p_val_spar;

/* Global settings */
extern bool g_decimal_b;						extern int g_decimal_val;							extern char g_decimal_str[10];
extern bool g_scientific_b;						extern int g_scientific_val;
extern bool g_nthread_b;						extern int g_nthread_val;
extern bool g_simd_b;							extern int g_simd_val;
extern bool g_benchmark_b;						extern int g_benchmark_val;
extern bool g_seed_b;							extern int g_seed_val;
extern bool g_tmpdir_b;							extern char *g_tmpdir_val;
extern bool g_progress_b;						extern int g_progress_val;

extern bool g_input_b;							extern char *g_input_val;
extern int g_input_row, g_input_col;
extern char ***g_filepath;
extern char ***g_filename;
extern char *g_filenamebuf;
extern int64 g_filetotlen;
extern int64 **g_filelen;
extern FILE ***g_filehandle;
extern BCFHEADER ***g_bcfheader;
extern char *g_filebuf;

extern bool g_format_b;							extern int g_format_val;
extern bool g_extracol_b;						extern int g_extracol_val;
extern bool g_output_b;							extern char *g_output_val;
extern bool g_indfile_b;
extern bool g_indtext_b;						extern char *g_indtext_val;
extern bool g_delimiter_b;						extern char g_delimiter_val;
extern bool g_linebreak_b;						extern char *g_linebreak_val;

/* Filters */
extern bool f_filter;
extern bool f_qual_b;							extern double f_qual_min, f_qual_max;
extern bool f_type_b;							extern int f_type_val;
extern bool f_original_b;						extern int f_original_val;
extern bool f_pop_b;							extern char *f_pop_val;
extern bool f_region_b;							extern char *f_region_val;
extern bool f_bmaf_b;							extern double f_bmaf_min, f_bmaf_max;
extern bool f_k_b;								extern int f_k_min, f_k_max;
extern bool f_n_b;								extern int f_n_min, f_n_max;
extern bool f_ptype_b;							extern double f_ptype_min, f_ptype_max;
extern bool f_pval_b;							extern double f_pval_min, f_pval_max;
extern bool f_model_b;							extern int f_model_val;
extern bool f_he_b;								extern double f_he_min, f_he_max;
extern bool f_ho_b;								extern double f_ho_min, f_ho_max;
extern bool f_pic_b;							extern double f_pic_min, f_pic_max;
extern bool f_ae_b;								extern double f_ae_min, f_ae_max;
extern bool f_dp_b;								extern uint f_dp_min, f_dp_max;
extern bool f_gq_b;								extern int f_gq_min, f_gq_max;
extern bool f_ploidy_b;							extern int f_ploidy_min, f_ploidy_max;
extern bool f_ntype_b;							extern int f_ntype_min, f_ntype_max;
extern bool f_nploidy_b;						extern int f_nploidy_min, f_nploidy_max;

/* Haplotype extraction */
extern bool haplotype;
extern bool haplotype_ptype_b;					extern double haplotype_ptype_min, haplotype_ptype_max;
extern bool haplotype_length_b;					extern int64 haplotype_length_min, haplotype_length_max;
extern bool haplotype_variants_b;				extern int haplotype_variants_min, haplotype_variants_max;
extern bool haplotype_interval_b;				extern int64 haplotype_interval_val;
extern bool haplotype_alleles_b;				extern int haplotype_alleles_min, haplotype_alleles_max;
extern bool haplotype_genotypes_b;				extern int haplotype_genotypes_min, haplotype_genotypes_max;

/* File conversion */
extern bool convert;
extern bool convert_format_b;					extern byte convert_format_val[N_MAX_OPTION];

/* Individual statistics */
extern bool indstat;
extern bool indstat_type_b;						extern byte indstat_type_val[N_MAX_OPTION];
extern bool indstat_model_b;					extern byte indstat_model_val[N_MAX_OPTION];
extern bool indstat_estimator_b;				extern byte indstat_estimator_val[N_MAX_OPTION];
extern bool indstat_ref_b;						extern byte indstat_ref_val[N_MAX_OPTION];
extern bool indstat_locus_b;					extern byte indstat_locus_val[N_MAX_OPTION];

/* Genetic Diversity */
extern bool diversity;
extern bool diversity_level_b;					extern byte diversity_level_val[N_MAX_OPTION];
extern bool diversity_model_b;					extern byte diversity_model_val[N_MAX_OPTION];

/* Genetic differentiation */
extern bool fst;
extern bool fst_level_b;						extern byte fst_level_val[N_MAX_OPTION];
extern bool fst_estimator_b;					extern byte fst_estimator_val[N_MAX_OPTION];
extern bool fst_fmt_b;							extern byte fst_fmt_val[N_MAX_OPTION];
extern bool fst_locus_b;						extern byte fst_locus_val[N_MAX_OPTION];
extern bool fst_test_b;							extern byte fst_test_val[N_MAX_OPTION];

/* Genetic distance */
extern bool gdist;
extern bool gdist_level_b;						extern byte gdist_level_val[N_MAX_OPTION];
extern bool gdist_weightmissing_b;				extern int gdist_weightmissing_val;
extern bool gdist_estimator_b;					extern byte gdist_estimator_val[N_MAX_OPTION];
extern bool gdist_fmt_b;						extern byte gdist_fmt_val[N_MAX_OPTION];

/* Analysis of molecular variances */
extern bool amova;
extern bool amova_method_b;						extern byte amova_method_val[N_MAX_OPTION];					extern int amova_cmethod_val;
extern bool amova_mutation_b;					extern byte amova_mutation_val[N_MAX_OPTION];				extern int amova_cmutation_val;
extern bool amova_ind_b;						extern byte amova_ind_val[N_MAX_OPTION];					extern int amova_cind_val;
extern bool amova_test_b;						extern int amova_test_val;
extern bool amova_nperm_b;						extern int amova_nperm_val;
extern bool amova_pseudo_b;						extern int amova_pseudo_val;
extern bool amova_printss_b;					extern int amova_printss_val;

/* Population assignment */
extern bool popas;
extern bool popas_model_b;						extern byte popas_model_val[N_MAX_OPTION];
extern bool popas_level_b;						extern byte popas_level_val[N_MAX_OPTION];
extern bool popas_error_b;						extern double popas_error_val;

/* Relatedness coefficient estimation */
extern bool relatedness;
extern bool relatedness_range_b;				extern byte relatedness_range_val[N_MAX_OPTION];
extern bool relatedness_fmt_b;					extern byte relatedness_fmt_val[N_MAX_OPTION];
extern bool relatedness_estimator_b;			extern byte relatedness_estimator_val[N_MAX_OPTION];

/* Kinship coefficient estimation */
extern bool kinship;
extern bool kinship_range_b;					extern byte kinship_range_val[N_MAX_OPTION];
extern bool kinship_fmt_b;						extern byte kinship_fmt_val[N_MAX_OPTION];
extern bool kinship_estimator_b;				extern byte kinship_estimator_val[N_MAX_OPTION];

/* Principal coordinate analysis */
extern bool pcoa;
extern bool pcoa_level_b;						extern byte pcoa_level_val[N_MAX_OPTION];
extern bool pcoa_dim_b;							extern int pcoa_dim_val;
extern bool pcoa_estimator_b;					extern byte pcoa_estimator_val[N_MAX_OPTION];

/* Hierarchical clustering */
extern bool cluster;
extern bool cluster_level_b;					extern byte cluster_level_val[N_MAX_OPTION];
extern bool cluster_method_b;					extern byte cluster_method_val[N_MAX_OPTION];
extern bool cluster_estimator_b;				extern byte cluster_estimator_val[N_MAX_OPTION];

/* Bayesian clustering */
extern bool structure;

/* Bayesian clustering: Model */
extern bool structure_admix_b;					extern int structure_admix_val;
extern bool structure_locpriori_b;				extern int structure_locpriori_val;
extern bool structure_f_b;						extern int structure_f_val;

/* Bayesian clustering: MCMC */
extern bool structure_krange_b;					extern int structure_krange_min, structure_krange_max;
extern bool structure_nburnin_b;				extern int structure_nburnin_val;
extern bool structure_nreps_b;					extern int structure_nreps_val;
extern bool structure_nthinning_b;				extern int structure_nthinning_val;
extern bool structure_nruns_b;					extern int structure_nruns_val;
extern bool structure_nadmburnin_b;				extern int structure_nadmburnin_val;

/* Bayesian clustering: Misc */
extern bool structure_lambda_b;					extern double structure_lambda_val;
extern bool structure_stdlambda_b;				extern double structure_stdlambda_val;
extern bool structure_maxlambda_b;				extern double structure_maxlambda_val;
extern bool structure_inferlambda_b;			extern int structure_inferlambda_val;
extern bool structure_difflambda_b;				extern int structure_difflambda_val;
extern bool structure_diversity_b;				extern int structure_diversity_val;

/* Bayesian clustering: Admix */
extern bool structure_alpha_b;					extern double structure_alpha_val;
extern bool structure_stdalpha_b;				extern double structure_stdalpha_val;
extern bool structure_maxalpha_b;				extern double structure_maxalpha_val;
extern bool structure_inferalpha_b;				extern int structure_inferalpha_val;
extern bool structure_diffalpha_b;				extern int structure_diffalpha_val;
extern bool structure_uniformalpha_b;			extern int structure_uniformalpha_val;
extern bool structure_alphapriora_b;			extern double structure_alphapriora_val;
extern bool structure_alphapriorb_b;			extern double structure_alphapriorb_val;
extern bool structure_metrofreq_b;				extern int structure_metrofreq_val;

/* Bayesian clustering: locpriori */
extern bool structure_r_b;						extern double structure_r_val;
extern bool structure_maxr_b;					extern double structure_maxr_val;
extern bool structure_stdr_b;					extern double structure_stdr_val;
extern bool structure_epseta_b;					extern double structure_epseta_val;
extern bool structure_epsgamma_b;				extern double structure_epsgamma_val;

/* Bayesian clustering: fmodel */
extern bool structure_pmeanf_b;					extern double structure_pmeanf_val;
extern bool structure_pstdf_b;					extern double structure_pstdf_val;
extern bool structure_stdf_b;					extern double structure_stdf_val;
extern bool structure_singlef_b;				extern int structure_singlef_val;

/* Ploidy inference from allelic depth distribution, test, do not use */
extern bool ploidyinfer;
extern bool ploidyinfer_type_b;					extern byte ploidyinfer_type_val[N_MAX_OPTION];
extern bool ploidyinfer_histogram_b;			extern int ploidyinfer_histogram_val;
extern bool ploidyinfer_nbins_b;				extern int ploidyinfer_nbins_val;

/* Allelic depth analysis, test, do not use */
extern int ad;

/* Analysis of spatial structure, test, do not use */
extern bool spa;
extern bool spa_dim_b;							extern int spa_dim_val;
extern bool spa_level_b;						extern int spa_level_val;
extern bool spa_odepth_b;						extern int spa_odepth_val;
extern bool spa_ofreq_b;						extern int spa_ofreq_val;
extern bool spa_coord_b;						extern char *spa_coord_val;
extern bool spa_truncate_b;						extern double spa_truncate_val, spa_truncate_val2;//n

/* Filter types */
extern bool diversity_filter;
extern bool genotype_filter;
extern bool individual_filter;
extern bool info_filter;

/* MISC */

extern int argc;
extern char **argv;

extern int genotype_digit;
extern int genotype_extracol;
extern int genotype_missing;
extern int genotype_ambiguous;

/* Functions */
extern TABLE<int, Huang2015ENTRY> *Huang2015_maps;		//Huang2015_maps[ploidylevel][hash] is a entry saves the ibs modex index and genotype pair pattern

/* Hash */
extern uint *cryptTable;								//crypt table for calculate hash for genotypes
