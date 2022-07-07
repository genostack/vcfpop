/* Parameter Functions */

#include "vcfpop.h"

/* Initialize parameters */
TARGET void SetDefaultParameters()
{
	SIMD_TYPE = 3;
	*(uint64*)&NA = 0xFFF8000000000000ull;
	srand((unsigned int)time(NULL));

	InitLock(glock1);
	InitLock(glock2);

	GetCurDir(CURDIR);

	diversity_filter = 0;
	genotype_filter = 0;
	individual_filter = 0;
	info_filter = 0;

	ad = 0;

	spa = false;
	spa_dim_b = false;				spa_dim_val = 2;
	spa_level_b = false;			spa_level_val = 1;
	spa_odepth_b = false;			spa_odepth_val = 2;
	spa_ofreq_b = false;			spa_ofreq_val = 2;
	spa_coord_b = false;			spa_coord_val = 0;
	spa_truncate_b = false;			spa_truncate_val = 0.001;		spa_truncate_val2 = 1.0 - spa_truncate_val;

	f_filter = false;
	f_qual_b = false;				f_qual_min = 0;					f_qual_max = 0;
	f_type_b = false;				f_type_val = 0;
	f_original_b = false;			f_original_val = 0;
	f_pop_b = false;				f_pop_val = (char*)"total";
	f_region_b = false;				f_region_val = 0;
	f_bmaf_b = false;				f_bmaf_min = 0;					f_bmaf_max = 0;
	f_k_b = false;					f_k_min = 0;					f_k_max = 0;
	f_n_b = false;					f_n_min = 0;					f_n_max = 0;
	f_ptype_b = false;				f_ptype_min = 0; 				f_ptype_max = 0;
	f_pval_b = false; 				f_pval_min = 0; 				f_pval_max = 0;
	f_model_b = false; 				f_model_val = 0;
	f_he_b = false; 				f_he_min = 0; 					f_he_max = 0;
	f_ho_b = false; 				f_ho_min = 0; 					f_ho_max = 0;
	f_pic_b = false; 				f_pic_min = 0; 					f_pic_max = 0;
	f_ae_b = false; 				f_ae_min = 0; 					f_ae_max = 0;
	f_dp_b = false; 				f_dp_min = 0; 					f_dp_max = 0;
	f_gq_b = false; 				f_gq_min = 0; 					f_gq_max = 0;
	f_ploidy_b = false; 			f_ploidy_min = 0; 				f_ploidy_max = 0;
	f_ntype_b = false;				f_ntype_min = 0;				f_ntype_max = 0;
	f_nploidy_b = false;			f_nploidy_min = 0;				f_nploidy_max = 0;

	g_decimal_b = false;			g_decimal_val = 5;				sprintf(g_decimal_str, "%%0.%df", 5);
	g_scientific_b = false;			g_scientific_val = 2;
	g_nthread_b = false;			g_nthread_val = 2;
	g_simd_b = false;				g_simd_val = 2;
	g_seed_b = false;				g_seed_val = rand();
	g_tmpdir_b = false;				g_tmpdir_val = CURDIR;
	g_progress_b = false;			g_progress_val = 80;

	g_input_b = false;				g_input_val = (char*)"";
	g_input_row = g_input_col = 0;
	g_filepath = NULL;
	g_filename = NULL;
	g_filenamebuf = NULL;
	g_format_b = false;				g_format_val = 1;
	g_extracol_b = false;			g_extracol_val = 0;
	g_output_b = false;				g_output_val = (char*)"vcfpop.out";
	g_indfile_b = false;
	g_indtext_b = false;			g_indtext_val = 0;
	g_delimiter_b = false;			g_delimiter_val = '\t';
	g_linebreak_b = false;			g_linebreak_val = (char*)"\n";

	haplotype = false;
	haplotype_ptype_b = false;	haplotype_ptype_min = 0.8;	haplotype_ptype_max = 1;
	haplotype_length_b = false;		haplotype_length_min = 1ull;	haplotype_length_max = 1000000ull;
	haplotype_variants_b = false;	haplotype_variants_min = 5;		haplotype_variants_max = 20;
	haplotype_interval_b = false;	haplotype_interval_val = 0ll;
	haplotype_alleles_b = false;	haplotype_alleles_min = 2;		haplotype_alleles_max = 65535;
	haplotype_genotypes_b = false;	haplotype_genotypes_min = 2;	haplotype_genotypes_max = 65535;

	convert = false;
	convert_format_b = false;		SetZero(convert_format_val, N_MAX_OPTION);		convert_format_val[1] = 1;

	indstat = false;
	indstat_type_b = false;			SetZero(indstat_type_val, N_MAX_OPTION);		indstat_type_val[1] = indstat_type_val[2] = 1;
	indstat_model_b = false;		SetZero(indstat_model_val, N_MAX_OPTION);		indstat_model_val[1] = 1;
	indstat_estimator_b = false;	SetZero(indstat_estimator_val, N_MAX_OPTION);	indstat_estimator_val[1] = 1;
	indstat_ref_b = false;			SetZero(indstat_ref_val, N_MAX_OPTION);			indstat_ref_val[3] = 1;
	indstat_locus_b = false;		SetZero(indstat_locus_val, N_MAX_OPTION);		indstat_locus_val[1] = 1;

	diversity = false;
	diversity_level_b = false;		SetZero(diversity_level_val, N_MAX_OPTION);		diversity_level_val[1] = diversity_level_val[4] = 1;
	diversity_model_b = false;		SetZero(diversity_model_val, N_MAX_OPTION);		diversity_model_val[1] = 1;

	fst = false;
	fst_level_b = false;			SetZero(fst_level_val, N_MAX_OPTION);			fst_level_val[1] = 1;
	fst_estimator_b = false;		SetZero(fst_estimator_val, N_MAX_OPTION);		fst_estimator_val[1] = 1;
	fst_fmt_b = false;				SetZero(fst_fmt_val, N_MAX_OPTION);				fst_estimator_val[1] = 1;
	fst_locus_b = false;			SetZero(fst_locus_val, N_MAX_OPTION);			fst_locus_val[1] = 1;
	fst_test_b = false;				fst_test_val[2] = 1;

	gdist = false;
	gdist_level_b = false;			SetZero(gdist_level_val, N_MAX_OPTION);			gdist_level_val[2] = 1;
	gdist_weightmissing_b = false;	gdist_weightmissing_val = 1;
	gdist_estimator_b = false;		SetZero(gdist_estimator_val, N_MAX_OPTION);		gdist_estimator_val[1] = 1;
	gdist_fmt_b = false;			SetZero(gdist_fmt_val, N_MAX_OPTION);			gdist_fmt_val[1] = 1;

	amova = false;
	amova_method_b = false;			amova_method_val[1] = 1;
	amova_mutation_b = false;		amova_mutation_val[1] = 1;
	amova_ind_b = false;			amova_ind_val[1] = 1;
	amova_test_b = false;			amova_test_val = 1;
	amova_nperm_b = false;			amova_nperm_val = 9999;
	amova_pseudo_b = false;			amova_pseudo_val = 50;
	amova_printss_b = false;		amova_printss_val = false;

	popas = false;
	popas_model_b = false;			SetZero(popas_model_val, N_MAX_OPTION);			popas_model_val[1] = 1;
	popas_level_b = false;			SetZero(popas_level_val, N_MAX_OPTION);			popas_level_val[1] = 1;
	popas_error_b = false;			popas_error_val = 0.01;

	relatedness = false;
	relatedness_range_b = false;	SetZero(relatedness_range_val, N_MAX_OPTION);	relatedness_range_val[3] = 1;
	relatedness_fmt_b = false;		SetZero(relatedness_fmt_val, N_MAX_OPTION);		relatedness_fmt_val[1] = 1;
	relatedness_estimator_b = false; SetZero(relatedness_estimator_val, N_MAX_OPTION); relatedness_estimator_val[1] = 1;

	kinship = false;
	kinship_range_b = false;		SetZero(kinship_range_val, N_MAX_OPTION);		kinship_range_val[3] = 1;
	kinship_fmt_b = false;			SetZero(kinship_fmt_val, N_MAX_OPTION);			kinship_fmt_val[1] = 1;
	kinship_estimator_b = false;	SetZero(kinship_estimator_val, N_MAX_OPTION);	kinship_estimator_val[9] = 1;

	pcoa = false;
	pcoa_level_b = false;			SetZero(pcoa_level_val, N_MAX_OPTION);			pcoa_level_val[1] = 1;
	pcoa_dim_b = false;				pcoa_dim_val = 3;
	pcoa_estimator_b = false;		SetZero(pcoa_estimator_val, N_MAX_OPTION);		pcoa_estimator_val[1] = 1;

	cluster = false;
	cluster_level_b = false;		SetZero(cluster_level_val, N_MAX_OPTION);		cluster_level_val[1] = 1;
	cluster_method_b = false;		SetZero(cluster_method_val, N_MAX_OPTION);		cluster_method_val[3] = 1;
	cluster_estimator_b = false;	SetZero(cluster_estimator_val, N_MAX_OPTION);	cluster_estimator_val[1] = 1;

	structure = false;
	structure_admix_b = false;		structure_admix_val = false;
	structure_locpriori_b = false;	structure_locpriori_val = false;
	structure_f_b = false;			structure_f_val = false;
	structure_krange_b = false;		structure_krange_min = 1;						structure_krange_max = 5;

	structure_nburnin_b = false;	structure_nburnin_val = 10000;
	structure_nreps_b = false;		structure_nreps_val = 100000;
	structure_nthinning_b = false;	structure_nthinning_val = 1;
	structure_nruns_b = false;		structure_nruns_val = 1;
	structure_nadmburnin_b = false; structure_nadmburnin_val = 500;

	structure_lambda_b = false;		structure_lambda_val = 1;
	structure_stdlambda_b = false;	structure_stdlambda_val = 0.3;
	structure_maxlambda_b = false;	structure_maxlambda_val = 10;
	structure_inferlambda_b = false;structure_inferlambda_val = 1;
	structure_difflambda_b = false; structure_difflambda_val = 1;
	structure_diversity_b = false;	structure_diversity_val = 2;

	structure_alpha_b = false;		structure_alpha_val = 1;
	structure_stdalpha_b = false;	structure_stdalpha_val = 0.025;
	structure_maxalpha_b = false;	structure_maxalpha_val = 10;
	structure_inferalpha_b = false;	structure_inferalpha_val = 2;
	structure_diffalpha_b = false;	structure_diffalpha_val = 0;
	structure_uniformalpha_b = false; structure_uniformalpha_val = 1;
	structure_alphapriora_b = false; structure_alphapriora_val = 0.025;
	structure_alphapriorb_b = false; structure_alphapriorb_val = 0.001;

	structure_metrofreq_b = false;	structure_metrofreq_val = 10;
	structure_r_b = false;			structure_r_val = 1;
	structure_maxr_b = false;		structure_maxr_val = 20;
	structure_stdr_b = false;		structure_stdr_val = 0.1;
	structure_epseta_b = false;		structure_epseta_val = 0.025;
	structure_epsgamma_b = false;	structure_epsgamma_val = 0.025;

	structure_pmeanf_b = false;		structure_pmeanf_val = 0.01;
	structure_pstdf_b = false;		structure_pstdf_val = 0.05;
	structure_stdf_b = false;		structure_stdf_val = 0.05;
	structure_singlef_b = false;	structure_singlef_val = 2;

	ploidyinfer = false;
	ploidyinfer_type_b = false;		SetZero(ploidyinfer_type_val, N_MAX_OPTION);	ploidyinfer_type_val[2] = ploidyinfer_type_val[4] = 1;
	ploidyinfer_histogram_b = false; ploidyinfer_histogram_val = false;
	ploidyinfer_nbins_b = false;	ploidyinfer_nbins_val = 20;

	diversity_filter = false;
	genotype_filter = false;
	individual_filter = false;
	info_filter = false;

	argc = NULL;
	argv = 0;
}

/* Set parameters */
TARGET void SetParameters(int _argc, char **_argv, bool isinparfile)
{
	if (!isinparfile) SetDefaultParameters();
	argc = _argc;
	argv = _argv;
	for (int i = 0; i < _argc; ++i)
	{
		if (!isinparfile && !LwrLineCmp("-p=", _argv[i]))
		{
			p_b = true;
			p_val = TrimQuote(_argv[i] + 3);
			break;
		}
		else if (!LwrLineCmp("-ad", _argv[i]))
		{
			ad = 1;
		}
		else if (!LwrLineCmp("-spa", _argv[i]))
		{
			if (!LwrStrCmp("-spa", _argv[i]))
				GetParBool(_argv[i], spa);
			else if (!LwrLineCmp("-spa_dim=", _argv[i]))
				GetParInteger(_argv[i], spa_dim_b, spa_dim_val, 1, 15);
			else if (!LwrLineCmp("-spa_level=", _argv[i]))
				GetParString(_argv[i], "ind|pop", spa_level_b, spa_level_val);
			else if (!LwrLineCmp("-spa_odepth=", _argv[i]))
				GetParString(_argv[i], "yes|no", spa_odepth_b, spa_odepth_val);
			else if (!LwrLineCmp("-spa_ofreq=", _argv[i]))
				GetParString(_argv[i], "yes|no", spa_ofreq_b, spa_ofreq_val);
			else if (!LwrLineCmp("-spa_truncate=", _argv[i]))
			{
				GetParDouble(_argv[i], spa_truncate_b, spa_truncate_val, 0.0000000001, 0.01);
				spa_truncate_val2 = 1.0 - spa_truncate_val;
			}
			else if (!LwrLineCmp("-spa_coord=", _argv[i]))
			{
				if (spa_coord_b) Exit("\nError: parameter %s has been assigned twice.\n", _argv[i]);
				spa_coord_b = true;
				spa_coord_val = ReplaceStr(_argv[i] + 11, "#n", "\n");
			}
		}
		else if (!LwrLineCmp("-f_", _argv[i]) || !LwrStrCmp("-f", _argv[i]))
		{
			if (!LwrStrCmp("-f", _argv[i]))
				GetParBool(_argv[i], f_filter);
			else if (!LwrLineCmp("-f_qual=", _argv[i]))
				GetRangeParDouble(_argv[i], f_qual_b, f_qual_min, f_qual_max, 0, 255);
			else if (!LwrLineCmp("-f_type=", _argv[i]))
				GetParString(_argv[i], "snp|indel|both", f_type_b, f_type_val);
			else if (!LwrLineCmp("-f_original=", _argv[i]))
				GetParString(_argv[i], "yes|no", f_original_b, f_original_val);
			else if (!LwrLineCmp("-f_pop=", _argv[i]))
			{
				if (f_pop_b) Exit("\nError: parameter %s has been assigned twice.\n", _argv[i]);
				f_pop_b = true;
				f_pop_val = TrimQuote(_argv[i] + 7);
				if (f_region_b) Exit("\nError: options -f_pop and -f_region are exclusive, delete sum.\n");
			}
			else if (!LwrLineCmp("-f_region=", _argv[i]))
			{
				if (f_region_b) Exit("\nError: parameter %s has been assigned twice.\n", _argv[i]);
				f_region_b = true;
				f_region_val = TrimQuote(_argv[i] + 10);
				if (f_pop_b) Exit("\nError: options -f_pop and -f_region are exclusive, delete sum.\n");
			}
			else if (!LwrLineCmp("-f_bmaf=", _argv[i]))
				GetRangeParDouble(_argv[i], f_bmaf_b, f_bmaf_min, f_bmaf_max, 0, 0.5);
			else if (!LwrLineCmp("-f_k=", _argv[i]))
				GetRangeParInteger(_argv[i], f_k_b, f_k_min, f_k_max, 1, 65535);
			else if (!LwrLineCmp("-f_n=", _argv[i]))
				GetRangeParInteger(_argv[i], f_n_b, f_n_min, f_n_max, 0, 65535);
			else if (!LwrLineCmp("-f_ptype=", _argv[i]))
				GetRangeParDouble(_argv[i], f_ptype_b, f_ptype_min, f_ptype_max, 0, 1);
			else if (!LwrLineCmp("-f_pval=", _argv[i]))
				GetRangeParDouble(_argv[i], f_pval_b, f_pval_min, f_pval_max, 0, 1);
			else if (!LwrLineCmp("-f_model=", _argv[i]))
				GetParString(_argv[i], "rcs|prcs|ces|pes", f_model_b, f_model_val);
			else if (!LwrLineCmp("-f_he=", _argv[i]))
				GetRangeParDouble(_argv[i], f_he_b, f_he_min, f_he_max, 0, 1);
			else if (!LwrLineCmp("-f_ho=", _argv[i]))
				GetRangeParDouble(_argv[i], f_ho_b, f_ho_min, f_ho_max, 0, 1);
			else if (!LwrLineCmp("-f_pic=", _argv[i]))
				GetRangeParDouble(_argv[i], f_pic_b, f_pic_min, f_pic_max, 0, 1);
			else if (!LwrLineCmp("-f_ae=", _argv[i]))
				GetRangeParDouble(_argv[i], f_ae_b, f_ae_min, f_ae_max, 0, 65535);
			else if (!LwrLineCmp("-f_dp=", _argv[i]))
				GetRangeParInteger(_argv[i], f_dp_b, f_dp_min, f_dp_max, 0, 0x7FFFFFFF);
			else if (!LwrLineCmp("-f_gq=", _argv[i]))
				GetRangeParInteger(_argv[i], f_gq_b, f_gq_min, f_gq_max, 0, 100);
			else if (!LwrLineCmp("-f_ploidy=", _argv[i]))
				GetRangeParInteger(_argv[i], f_ploidy_b, f_ploidy_min, f_ploidy_max, 1, N_MAX_PLOIDY);
			else if (!LwrLineCmp("-f_ntype=", _argv[i]))
				GetRangeParInteger(_argv[i], f_ntype_b, f_ntype_min, f_ntype_max, 0, 99999999);
			else if (!LwrLineCmp("-f_nploidy=", _argv[i]))
				GetRangeParInteger(_argv[i], f_nploidy_b, f_nploidy_min, f_nploidy_max, 1, N_MAX_PLOIDY);
			else
				Exit("\nError: Unrecognized parameter: %s\n", _argv[i]);
		}
		else if (!LwrLineCmp("-g_", _argv[i]))
		{
			if (!LwrLineCmp("-g_decimal=", _argv[i]))
			{
				GetParInteger(_argv[i], g_decimal_b, g_decimal_val, 1, 15);
				sprintf(g_decimal_str, g_scientific_val == 1 ? "%%0.%de" : "%%0.%df", g_decimal_val);
			}
			else if (!LwrLineCmp("-g_scientific=", _argv[i]))
			{
				GetParString(_argv[i], "yes|no", g_scientific_b, g_scientific_val);
				sprintf(g_decimal_str, g_scientific_val == 1 ? "%%0.%de" : "%%0.%df", g_decimal_val);
			}
			else if (!LwrLineCmp("-g_nthread=", _argv[i]))
				GetParInteger(_argv[i], g_nthread_b, g_nthread_val, 1, 4096);
			else if (!LwrLineCmp("-g_simd=", _argv[i]))
			{
				GetParString(_argv[i], "mmx|sse|avx|avx512", g_simd_b, g_simd_val);
				SIMD_TYPE = g_simd_val;
			}
			else if (!LwrLineCmp("-g_benchmark=", _argv[i]))
				GetParString(_argv[i], "yes|no", g_benchmark_b, g_benchmark_val);
			else if (!LwrLineCmp("-g_seed=", _argv[i]))
			{
				GetParInteger(_argv[i], g_seed_b, g_seed_val, 0, 0x7FFFFFFF);
				if (g_seed_val == 0) g_seed_val = (int)time(NULL);
			}
			else if (!LwrLineCmp("-g_tmpdir=", _argv[i]))
			{
				if (g_tmpdir_b) Exit("\nError: parameter %s has been assigned twice.\n", _argv[i]);
				g_tmpdir_b = true;
				uint64 tlen = strlen(_argv[i] + 10);
				g_tmpdir_val = new char[tlen + 4];
				strcpy(g_tmpdir_val, TrimQuote(_argv[i] + 10));
#ifdef _WIN64
				if (g_tmpdir_val[tlen - 1] != '\\')
				{
					g_tmpdir_val[tlen] = '\\';
					g_tmpdir_val[tlen + 1] = '\0';
				}
				char *t;
				t = ReplaceStr(g_tmpdir_val, "/", "\\"); strcpy(g_tmpdir_val, t); delete[] t;
				t = ReplaceStr(g_tmpdir_val, "\\\\", "\\"); strcpy(g_tmpdir_val, t); delete[] t;
#else
				if (g_tmpdir_val[tlen - 1] != '/')
				{
					g_tmpdir_val[tlen] = '/';
					g_tmpdir_val[tlen + 1] = '\0';
				}
				char *t;
				t = ReplaceStr(g_tmpdir_val, "\\", "/"); strcpy(g_tmpdir_val, t); delete[] t;
				t = ReplaceStr(g_tmpdir_val, "//", "/"); strcpy(g_tmpdir_val, t); delete[] t;
#endif
			}
			else if (!LwrLineCmp("-g_progress=", _argv[i]))
				GetParInteger(_argv[i], g_progress_b, g_progress_val, 10, 1000000);
			else if (!LwrLineCmp("-g_input=", _argv[i]))
			{
				if (g_input_b) Exit("\nError: parameter %s has been assigned twice.\n", _argv[i]);
				g_input_b = true;
				g_input_val = TrimQuote(_argv[i] + 9);
				g_filenamebuf = new char[strlen(g_input_val) + 1];
				strcpy(g_filenamebuf, g_input_val);

				char **ts = SplitStr(g_filenamebuf, '&', g_input_row);
				g_input_col = CountChar(ts[0], '|') + 1;

				g_filepath = new char**[g_input_row];
				g_filename = new char**[g_input_row];
				g_filehandle = new FILE**[g_input_row];
				g_filelen = new int64*[g_input_row];
				for (int ii = 0; ii < g_input_row; ++ii)
				{
					int tcol = 0;
					g_filepath[ii] = SplitStr(ts[ii], '|', tcol);
					if (tcol != g_input_col) Exit("\nError: multiple input vcf files at each row should have the same number of columns.\n");
					g_filename[ii] = new char*[g_input_col];
					g_filehandle[ii] = new FILE*[g_input_col];
					g_filelen[ii] = new int64[g_input_col];
					for (int j = 0; j < g_input_col; ++j)
					{
						g_filename[ii][j] = g_filepath[ii][j] + strlen(g_filepath[ii][j]);
						while (g_filename[ii][j] > g_filepath[ii][j])
						{
							if (*g_filename[ii][j] == '/' || *g_filename[ii][j] == '\\')
							{
								g_filename[ii][j]++;
								break;
							}
							g_filename[ii][j]--;
						}
					}
				}
				delete[] ts;
			}
			else if (!LwrLineCmp("-g_format=", _argv[i]))
				GetParString(_argv[i], "vcf|bcf|genepop|spagedi|cervus|arlequin|structure|polygene|polyrelatedness", g_format_b, g_format_val);
			else if (!LwrLineCmp("-g_extracol=", _argv[i]))
				GetParInteger(_argv[i], g_extracol_b, g_extracol_val, 0, 4096);
			else if (!LwrLineCmp("-g_output=", _argv[i]))
			{
				if (g_output_b) Exit("\nError: parameter %s has been assigned twice.\n", _argv[i]);
				g_output_b = true;
				g_output_val = TrimQuote(_argv[i] + 10);
			}
			else if (!LwrLineCmp("-g_indfile=", _argv[i]))
			{
				if (g_indfile_b) Exit("\nError: parameter %s has been assigned twice.\n", _argv[i]);
				g_indfile_b = true;
				g_indtext_val = ReadAllText(TrimQuote(_argv[i] + 11));
				if (g_indtext_b) Exit("\nError: options -g_indfile and -g_indtext are exclusive, delete sum.\n");
			}
			else if (!LwrLineCmp("-g_indtext=", _argv[i]))
			{
				if (g_indtext_b) Exit("\nError: parameter %s has been assigned twice.\n", _argv[i]);
				g_indtext_b = true;
				g_indtext_val = ReplaceStr(_argv[i] + 11, "#n", "\n");
				if (g_indfile_b) Exit("\nError: options -g_indfile and -g_indtext are exclusive, delete sum.\n");
			}
			else if (!LwrLineCmp("-g_delimiter=", _argv[i]))
			{
				int tval = 0;
				GetParString(_argv[i], "comma|tab", g_delimiter_b, tval);
				if (tval == 1) g_delimiter_val = ',';
				else g_delimiter_val = '\t';
			}
			else if (!LwrLineCmp("-g_linebreak=", _argv[i]))
			{
				int tval = 0;
				GetParString(_argv[i], "unix|win", g_linebreak_b, tval);
				if (tval == 1) g_linebreak_val = (char*)"\n";
				else g_linebreak_val = (char*)"\r\n";
			}
			else
				Exit("\nError: Unrecognized parameter: %s\n", _argv[i]);
		}
		else if (!LwrLineCmp("-haplotype", _argv[i]))
		{
			if (!LwrStrCmp("-haplotype", _argv[i]))
				GetParBool(_argv[i], haplotype);
			else if (!LwrLineCmp("-haplotype_ptype=", _argv[i]))
				GetRangeParDouble(_argv[i], haplotype_ptype_b, haplotype_ptype_min, haplotype_ptype_max, 0.1, 1);
			else if (!LwrLineCmp("-haplotype_length=", _argv[i]))
				GetRangeParLong(_argv[i], haplotype_length_b, haplotype_length_min, haplotype_length_max, 1ll, 100000000000ll);
			else if (!LwrLineCmp("-haplotype_variants=", _argv[i]))
				GetRangeParInteger(_argv[i], haplotype_variants_b, haplotype_variants_min, haplotype_variants_max, 1, 1000000000);
			else if (!LwrLineCmp("-haplotype_interval_val=", _argv[i]))
				GetParLong(_argv[i], haplotype_interval_b, haplotype_interval_val, 0ll, 100000000000ll);
			else if (!LwrLineCmp("-haplotype_alleles=", _argv[i]))
				GetRangeParInteger(_argv[i], haplotype_alleles_b, haplotype_alleles_min, haplotype_alleles_max, 2, 65535);
			else if (!LwrLineCmp("-haplotype_genotypes=", _argv[i]))
				GetRangeParInteger(_argv[i], haplotype_genotypes_b, haplotype_genotypes_min, haplotype_genotypes_max, 2, 65535);
		}
		else if (!LwrLineCmp("-convert", _argv[i]))
		{
			if (!LwrStrCmp("-convert", _argv[i]))
				GetParBool(_argv[i], convert);
			else if (!LwrLineCmp("-convert_format=", _argv[i]))
				GetParStringMultiSel(_argv[i], "genepop|spagedi|cervus|arlequin|structure|polygene|polyrelatedness", convert_format_b, convert_format_val);
			else
				Exit("\nError: Unrecognized parameter: %s\n", _argv[i]);
		}
		else if (!LwrLineCmp("-indstat", _argv[i]))
		{
			if (!LwrStrCmp("-indstat", _argv[i]))
				GetParBool(_argv[i], indstat);
			else if (!LwrLineCmp("-indstat_type=", _argv[i]))
				GetParStringMultiSel(_argv[i], "hidx|lnl|f|theta", indstat_type_b, indstat_type_val);
			else if (!LwrLineCmp("-indstat_model=", _argv[i]))
				GetParStringMultiSel(_argv[i], "rcs|prcs|ces|pes", indstat_model_b, indstat_model_val);
			else if (!LwrLineCmp("-indstat_estimator=", _argv[i]))
				GetParStringMultiSel(_argv[i], "Ritland1996|Loiselle1995|Weir1996", indstat_estimator_b, indstat_estimator_val);
			else if (!LwrLineCmp("-indstat_ref=", _argv[i]))
				GetParStringMultiSel(_argv[i], "pop|reg|total", indstat_ref_b, indstat_ref_val);
			else if (!LwrLineCmp("-indstat_locus=", _argv[i]))
				GetParStringMultiSel(_argv[i], "all|each", indstat_locus_b, indstat_locus_val);
			else
				Exit("\nError: Unrecognized parameter: %s\n", _argv[i]);
		}
		else if (!LwrLineCmp("-diversity", _argv[i]))
		{
			if (!LwrStrCmp("-diversity", _argv[i]))
				GetParBool(_argv[i], diversity);
			else if (!LwrLineCmp("-diversity_level=", _argv[i]))
				GetParStringMultiSel(_argv[i], "pop|reg|tot|popXloc|regXloc|totXloc", diversity_level_b, diversity_level_val);
			else if (!LwrLineCmp("-diversity_model=", _argv[i]))
				GetParStringMultiSel(_argv[i], "rcs|prcs|ces|pes", diversity_model_b, diversity_model_val);
			else
				Exit("\nError: Unrecognized parameter: %s\n", _argv[i]);
		}
		else if (!LwrLineCmp("-fst", _argv[i]))
		{
			if (!LwrStrCmp("-fst", _argv[i]))
				GetParBool(_argv[i], fst);
			else if (!LwrLineCmp("-fst_level=", _argv[i]))
				GetParStringMultiSel(_argv[i], "regXtot|popXtot|popXreg|reg|pop", fst_level_b, fst_level_val);
			else if (!LwrLineCmp("-fst_estimator=", _argv[i]))
				GetParStringMultiSel(_argv[i], "Nei1973|Weir1984|Hudson1992|Slatkin1995|Hedrick2005|Jost2008|Huang2021_homo|Huang2021_aniso", fst_estimator_b, fst_estimator_val);
			else if (!LwrLineCmp("-fst_fmt=", _argv[i]))
				GetParStringMultiSel(_argv[i], "matrix|table", fst_fmt_b, fst_fmt_val);
			else if (!LwrLineCmp("-fst_locus=", _argv[i]))
				GetParStringMultiSel(_argv[i], "all|each", fst_locus_b, fst_locus_val);
			else if (!LwrLineCmp("-fst_test=", _argv[i]))
				GetParStringMultiSel(_argv[i], "genotype|allele", fst_test_b, fst_test_val);
			else
				Exit("\nError: Unrecognized parameter: %s\n", _argv[i]);
		}
		else if (!LwrLineCmp("-gdist", _argv[i]))
		{
			if (!LwrStrCmp("-gdist", _argv[i]))
				GetParBool(_argv[i], gdist);
			else if (!LwrLineCmp("-gdist_level=", _argv[i]))
				GetParStringMultiSel(_argv[i], "ind|pop|reg", gdist_level_b, gdist_level_val);
			else if (!LwrLineCmp("-gdist_weightmissing=", _argv[i]))
				GetParString(_argv[i], "yes|no", gdist_weightmissing_b, gdist_weightmissing_val);
			else if (!LwrLineCmp("-gdist_estimator=", _argv[i]))
				GetParStringMultiSel(_argv[i], "Nei1972|Cavalli-Sforza1967|Reynolds1983|Nei1983|Euclidean|Goldstein1995|Nei1974|Roger1972|Slatkin_Nei1973|Slatkin_Weir1984|Slatkin_Hudson1992|Slatkin_Slatkin1995|Slatkin_Hedrick2005|Slatkin_Jost2008|Slatkin_Huang2021_homo|Slatkin_Huang2021_aniso|Reynolds_Nei1973|Reynolds_Weir1984|Reynolds_Hudson1992|Reynolds_Slatkin1995|Reynolds_Hedrick2005|Reynolds_Jost2008|Reynolds_Huang2021_homo|Reynolds_Huang2021_aniso", gdist_estimator_b, gdist_estimator_val);
			else if (!LwrLineCmp("-gdist_fmt=", _argv[i]))
				GetParStringMultiSel(_argv[i], "matrix|table", gdist_fmt_b, gdist_fmt_val);
			else
				Exit("\nError: Unrecognized parameter: %s\n", _argv[i]);
		}
		else if (!LwrLineCmp("-amova", _argv[i]))
		{
			if (!LwrStrCmp("-amova", _argv[i]))
				GetParBool(_argv[i], amova);
			else if (!LwrLineCmp("-amova_method=", _argv[i]))
				GetParStringMultiSel(_argv[i], "homoploid|anisoploid|likelihood", amova_method_b, amova_method_val);
			else if (!LwrLineCmp("-amova_mutation=", _argv[i]))
				GetParStringMultiSel(_argv[i], "iam|smm", amova_mutation_b, amova_mutation_val);
			else if (!LwrLineCmp("-amova_ind=", _argv[i]))
				GetParStringMultiSel(_argv[i], "yes|no", amova_ind_b, amova_ind_val);
			else if (!LwrLineCmp("-amova_test=", _argv[i]))
				GetParString(_argv[i], "yes|no", amova_test_b, amova_test_val);
			else if (!LwrLineCmp("-amova_nperm=", _argv[i]))
				GetParInteger(_argv[i], amova_nperm_b, amova_nperm_val, 1, 99999999);
			else if (!LwrLineCmp("-amova_pseudo=", _argv[i]))
				GetParInteger(_argv[i], amova_pseudo_b, amova_pseudo_val, 0, 9999);
			else if (!LwrLineCmp("-amova_printss=", _argv[i]))
				GetParString(_argv[i], "yes|no", amova_printss_b, amova_printss_val);
			else
				Exit("\nError: Unrecognized parameter: %s\n", _argv[i]);
		}
		else if (!LwrLineCmp("-popas", _argv[i]))
		{
			if (!LwrStrCmp("-popas", _argv[i]))
				GetParBool(_argv[i], popas);
			else if (!LwrLineCmp("-popas_model=", _argv[i]))
				GetParStringMultiSel(_argv[i], "rcs|prcs|ces|pes", popas_model_b, popas_model_val);
			else if (!LwrLineCmp("-popas_level=", _argv[i]))
				GetParStringMultiSel(_argv[i], "pop|reg", popas_level_b, popas_level_val);
			else if (!LwrLineCmp("-popas_error=", _argv[i]))
				GetParDouble(_argv[i], popas_error_b, popas_error_val, 0, 0.2);
			else
				Exit("\nError: Unrecognized parameter: %s\n", _argv[i]);
		}
		else if (!LwrLineCmp("-relatedness", _argv[i]))
		{
			if (!LwrStrCmp("-relatedness", _argv[i]))
				GetParBool(_argv[i], relatedness);
			else if (!LwrLineCmp("-relatedness_range=", _argv[i]))
				GetParStringMultiSel(_argv[i], "pop|reg|total", relatedness_range_b, relatedness_range_val);
			else if (!LwrLineCmp("-relatedness_fmt=", _argv[i]))
				GetParStringMultiSel(_argv[i], "matrix|table", relatedness_fmt_b, relatedness_fmt_val);
			else if (!LwrLineCmp("-relatedness_estimator=", _argv[i]))
				GetParStringMultiSel(_argv[i], "Lynch1999|Wang2002|Thomas2010|Li1993|Queller1989|Huang2016A|Huang2016B|Milligan2003|Anderson2007|Huang2014|Huang2015|Ritland1996_modified|Loiselle1995_modified|Ritland1996|Loiselle1995|Weir1996", relatedness_estimator_b, relatedness_estimator_val);
			else
				Exit("\nError: Unrecognized parameter: %s\n", _argv[i]);
		}
		else if (!LwrLineCmp("-kinship", _argv[i]))
		{
			if (!LwrStrCmp("-kinship", _argv[i]))
				GetParBool(_argv[i], kinship);
			else if (!LwrLineCmp("-kinship_range=", _argv[i]))
				GetParStringMultiSel(_argv[i], "pop|reg|total", kinship_range_b, kinship_range_val);
			else if (!LwrLineCmp("-kinship_fmt=", _argv[i]))
				GetParStringMultiSel(_argv[i], "matrix|table", kinship_fmt_b, kinship_fmt_val);
			else if (!LwrLineCmp("-kinship_estimator=", _argv[i]))
				GetParStringMultiSel(_argv[i], "Ritland1996|Loiselle1995|Weir1996", kinship_estimator_b, kinship_estimator_val);
			else
				Exit("\nError: Unrecognized parameter: %s\n", _argv[i]);
		}
		else if (!LwrLineCmp("-pcoa", _argv[i]))
		{
			if (!LwrStrCmp("-pcoa", _argv[i]))
				GetParBool(_argv[i], pcoa);
			else if (!LwrLineCmp("-pcoa_level=", _argv[i]))
				GetParStringMultiSel(_argv[i], "ind|pop|reg", pcoa_level_b, pcoa_level_val);
			else if (!LwrLineCmp("-pcoa_dim=", _argv[i]))
				GetParInteger(_argv[i], pcoa_dim_b, pcoa_dim_val, 1, 4096);
			else if (!LwrLineCmp("-pcoa_estimator=", _argv[i]))
				GetParStringMultiSel(_argv[i], "Nei1972|Cavalli-Sforza1967|Reynolds1983|Nei1983|Euclidean|Goldstein1995|Nei1974|Roger1972|Slatkin_Nei1973|Slatkin_Weir1984|Slatkin_Hudson1992|Slatkin_Slatkin1995|Slatkin_Hedrick2005|Slatkin_Jost2008|Slatkin_Huang2021_homo|Slatkin_Huang2021_aniso|Reynolds_Nei1973|Reynolds_Weir1984|Reynolds_Hudson1992|Reynolds_Slatkin1995|Reynolds_Hedrick2005|Reynolds_Jost2008|Reynolds_Huang2021_homo|Reynolds_Huang2021_aniso", pcoa_estimator_b, pcoa_estimator_val);
			else
				Exit("\nError: Unrecognized parameter: %s\n", _argv[i]);
		}
		else if (!LwrLineCmp("-cluster", _argv[i]))
		{
			if (!LwrStrCmp("-cluster", _argv[i]))
				GetParBool(_argv[i], cluster);
			else if (!LwrLineCmp("-cluster_level=", _argv[i]))
				GetParStringMultiSel(_argv[i], "ind|pop|reg", cluster_level_b, cluster_level_val);
			else if (!LwrLineCmp("-cluster_method=", _argv[i]))
				GetParStringMultiSel(_argv[i], "NEAREST|FURTHEST|UPGMA|WPGMA|UPGMC|WPGMC|WARD", cluster_method_b, cluster_method_val);
			else if (!LwrLineCmp("-cluster_estimator=", _argv[i]))
				GetParStringMultiSel(_argv[i], "Nei1972|Cavalli-Sforza1967|Reynolds1983|Nei1983|Euclidean|Goldstein1995|Nei1974|Roger1972|Slatkin_Nei1973|Slatkin_Weir1984|Slatkin_Hudson1992|Slatkin_Slatkin1995|Slatkin_Hedrick2005|Slatkin_Jost2008|Slatkin_Huang2021_homo|Slatkin_Huang2021_aniso|Reynolds_Nei1973|Reynolds_Weir1984|Reynolds_Hudson1992|Reynolds_Slatkin1995|Reynolds_Hedrick2005|Reynolds_Jost2008|Reynolds_Huang2021_homo|Reynolds_Huang2021_aniso", cluster_estimator_b, cluster_estimator_val);
			else
				Exit("\nError: Unrecognized parameter: %s\n", _argv[i]);
		}
		else if (!LwrLineCmp("-structure", _argv[i]))
		{
			if (!LwrStrCmp("-structure", _argv[i]))
				GetParBool(_argv[i], structure);
			else if (!LwrLineCmp("-structure_admix=", _argv[i]))
				GetParString(_argv[i], "yes|no", structure_admix_b, structure_admix_val);
			else if (!LwrLineCmp("-structure_locpriori=", _argv[i]))
				GetParString(_argv[i], "yes|no", structure_locpriori_b, structure_locpriori_val);
			else if (!LwrLineCmp("-structure_f=", _argv[i]))
				GetParString(_argv[i], "yes|no", structure_f_b, structure_f_val);

			else if (!LwrLineCmp("-structure_krange=", _argv[i]))
				GetRangeParInteger(_argv[i], structure_krange_b, structure_krange_min, structure_krange_max, 1, 1000);
			else if (!LwrLineCmp("-structure_nburnin=", _argv[i]))
				GetParInteger(_argv[i], structure_nburnin_b, structure_nburnin_val, 1000, 10000000);
			else if (!LwrLineCmp("-structure_nreps=", _argv[i]))
				GetParInteger(_argv[i], structure_nreps_b, structure_nreps_val, 1000, 10000000);
			else if (!LwrLineCmp("-structure_nthinning=", _argv[i]))
				GetParInteger(_argv[i], structure_nthinning_b, structure_nthinning_val, 1, 1000);
			else if (!LwrLineCmp("-structure_nruns=", _argv[i]))
				GetParInteger(_argv[i], structure_nruns_b, structure_nruns_val, 1, 1000);

			else if (!LwrLineCmp("-structure_nadmburnin=", _argv[i]))
				GetParInteger(_argv[i], structure_nadmburnin_b, structure_nadmburnin_val, 100, 10000000);
			else if (!LwrLineCmp("-structure_lambda=", _argv[i]))
				GetParDouble(_argv[i], structure_lambda_b, structure_lambda_val, 0, 10000);
			else if (!LwrLineCmp("-structure_stdlambda=", _argv[i]))
				GetParDouble(_argv[i], structure_stdlambda_b, structure_stdlambda_val, 0, 10000);
			else if (!LwrLineCmp("-structure_maxlambda=", _argv[i]))
				GetParDouble(_argv[i], structure_maxlambda_b, structure_maxlambda_val, 0, 10000);
			else if (!LwrLineCmp("-structure_inferlambda=", _argv[i]))
				GetParString(_argv[i], "yes|no", structure_inferlambda_b, structure_inferlambda_val);
			else if (!LwrLineCmp("-structure_difflambda=", _argv[i]))
				GetParString(_argv[i], "yes|no", structure_difflambda_b, structure_difflambda_val);
			else if (!LwrLineCmp("-structure_diversity=", _argv[i]))
				GetParString(_argv[i], "yes|no", structure_diversity_b, structure_diversity_val);

			else if (!LwrLineCmp("-structure_alpha=", _argv[i]))
				GetParDouble(_argv[i], structure_alpha_b, structure_alpha_val, 0, 10000);
			else if (!LwrLineCmp("-structure_inferalpha=", _argv[i]))
				GetParString(_argv[i], "yes|no", structure_inferalpha_b, structure_inferalpha_val);
			else if (!LwrLineCmp("-structure_diffalpha=", _argv[i]))
				GetParString(_argv[i], "yes|no", structure_diffalpha_b, structure_diffalpha_val);
			else if (!LwrLineCmp("-structure_uniformalpha=", _argv[i]))
				GetParString(_argv[i], "yes|no", structure_uniformalpha_b, structure_uniformalpha_val);
			else if (!LwrLineCmp("-structure_stdalpha=", _argv[i]))
				GetParDouble(_argv[i], structure_stdalpha_b, structure_stdalpha_val, 0, 10000);
			else if (!LwrLineCmp("-structure_maxalpha=", _argv[i]))
				GetParDouble(_argv[i], structure_maxalpha_b, structure_maxalpha_val, 0, 10000);
			else if (!LwrLineCmp("-structure_alphapriora=", _argv[i]))
				GetParDouble(_argv[i], structure_alphapriora_b, structure_alphapriora_val, 0, 10000);
			else if (!LwrLineCmp("-structure_alphapriorb=", _argv[i]))
				GetParDouble(_argv[i], structure_alphapriorb_b, structure_alphapriorb_val, 0, 10000);
			else if (!LwrLineCmp("-structure_metrofreq=", _argv[i]))
				GetParInteger(_argv[i], structure_metrofreq_b, structure_metrofreq_val, 0, 1000000);

			else if (!LwrLineCmp("-structure_r=", _argv[i]))
				GetParDouble(_argv[i], structure_r_b, structure_r_val, 0, 10000);
			else if (!LwrLineCmp("-structure_maxr=", _argv[i]))
				GetParDouble(_argv[i], structure_maxr_b, structure_maxr_val, 0, 10000);
			else if (!LwrLineCmp("-structure_stdr=", _argv[i]))
				GetParDouble(_argv[i], structure_stdr_b, structure_stdr_val, 0, 10000);
			else if (!LwrLineCmp("-structure_epseta=", _argv[i]))
				GetParDouble(_argv[i], structure_epseta_b, structure_epseta_val, 0, 10000);
			else if (!LwrLineCmp("-structure_epsgamma=", _argv[i]))
				GetParDouble(_argv[i], structure_epsgamma_b, structure_epsgamma_val, 0, 10000);

			else if (!LwrLineCmp("-structure_pmeanf=", _argv[i]))
				GetParDouble(_argv[i], structure_pmeanf_b, structure_pmeanf_val, 0, 10000);
			else if (!LwrLineCmp("-structure_pstdf=", _argv[i]))
				GetParDouble(_argv[i], structure_pstdf_b, structure_pstdf_val, 0, 10000);
			else if (!LwrLineCmp("-structure_stdf=", _argv[i]))
				GetParDouble(_argv[i], structure_stdf_b, structure_stdf_val, 0, 10000);
			else if (!LwrLineCmp("-structure_singlef=", _argv[i]))
				GetParString(_argv[i], "yes|no", structure_singlef_b, structure_singlef_val);
			else
				Exit("\nError: Unrecognized parameter: %s\n", _argv[i]);
		}
		else if (!LwrLineCmp("-ploidyinfer", _argv[i]))
		{
			if (!LwrStrCmp("-ploidyinfer", _argv[i]))
				GetParBool(_argv[i], ploidyinfer);
			else if (!LwrLineCmp("-ploidyinfer_type=", _argv[i]))
				GetParStringMultiSel(_argv[i], "1|2|3|4|5|6|7|8|9|10", ploidyinfer_type_b, ploidyinfer_type_val);
			else if (!LwrLineCmp("-ploidyinfer_histogram=", _argv[i]))
				GetParString(_argv[i], "yes|no", ploidyinfer_histogram_b, ploidyinfer_histogram_val);
			else if (!LwrLineCmp("-ploidyinfer_nbins=", _argv[i]))
				GetParInteger(_argv[i], ploidyinfer_nbins_b, ploidyinfer_nbins_val, 10, 100);
			else
				Exit("\nError: Unrecognized parameter: %s\n", _argv[i]);
		}
	}

	if (f_filter && (f_bmaf_b || f_k_b || f_n_b || f_ptype_b || f_pval_b || f_he_b || f_ho_b || f_pic_b || f_ae_b))
		diversity_filter = true;
	if (f_filter && (f_dp_b || f_gq_b || f_ploidy_b))
		genotype_filter = true;
	if (f_filter && (f_ntype_b || f_nploidy_b))
		individual_filter = true;
	if (f_filter && (f_qual_b || f_type_b || f_original_b))
		info_filter = true;

	if (!isinparfile && p_b)
	{
		p_val_spar = ReadAllText(p_val);
		char *tpar = new char[strlen(p_val_spar) + strlen(_argv[0]) + 4];
		sprintf(tpar, "\"%s\"\n%s", _argv[0], p_val_spar);
		delete[] p_val_spar;
		p_val_spar = tpar;
		p_val_spar = ReplaceStr(p_val_spar, "\r\n", " "); delete[] tpar; tpar = p_val_spar;
		p_val_spar = ReplaceStr(p_val_spar, "\n", " "); delete[] tpar; tpar = p_val_spar;
		p_val_spar = ReplaceStr(p_val_spar, "\t", " "); delete[] tpar; tpar = p_val_spar;
		_argv = SplitStr(p_val_spar, ' ', _argc);
		for (int i = 0; i < _argc; ++i)
			printf("%s\r\n", _argv[i]);
		SetParameters(_argc, _argv, true);
	}
}

/* Release parameter memory */
TARGET void ReleaseParameters()
{
	if (g_tmpdir_b) delete[] g_tmpdir_val;
	if (g_indtext_val) delete[] g_indtext_val;
	for (int i = 0; i < g_input_row; ++i)
	{
		delete[] g_filepath[i];
		delete[] g_filename[i];
	}
	if (g_filepath) delete[] g_filepath;
	if (g_filename) delete[] g_filename;

	if (p_b)
	{
		delete[] p_val_spar;
		delete[] argv;
	}
}

/* Write parameters to result file */
TARGET void WriteParameters(FILE *f1, char *type, char *prefix)
{
	fprintf(f1, "%s  Parameter: %s%s", prefix, argv[0], g_linebreak_val);
	for (int i = 1; i < argc; ++i)
	{
		if (!LwrLineCmp("-f_", argv[i]) || !LwrLineCmp("-g_", argv[i]) || !LwrLineCmp(type, argv[i]))
			fprintf(f1, "%s    %s%s", prefix, argv[i], g_linebreak_val);
	}
}