import tempfile
import rpy2.robjects
import pandas as pd
def sleuth_wrapper(exp_design_matrix, full_model_string, reduced_model_string, dry_run=False, outfile=None, organism='hsapiens_gene_ensembl'):
    fname_design = tempfile.NamedTemporaryFile().name
    fname_de_table = tempfile.NamedTemporaryFile().name
    exp_design_matrix.to_csv(fname_design)

    lib_loading = """
    library(Biobase)
    library(sleuth)
    """
    #                          host = 'ensembl.org'
    sleuth_load_design_str = f's2c <- read.table("{fname_design}", header = TRUE, stringsAsFactors=FALSE, sep=",")'
    biomart_str = f"""
    mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                             dataset = "{organism}")
    t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                         "external_gene_name", "transcript_biotype"), mart = mart)
    t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                         ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
    """
    prep_str = 'so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, target_mapping = t2g)'
    sleuth_fit1_str = f'so <- sleuth_fit(so, {full_model_string}, "full")'
    sleuth_fit2_str = f'so <- sleuth_fit(so, {reduced_model_string}, "reduced")'
    sleuth_lrt_str = 'so <- sleuth_lrt(so, "reduced", "full")'
    sleuth_table_str = 'sleuth_table <- sleuth_results(so, "reduced:full", "lrt", show_all = FALSE)'
    sleuth_write_str = f'write.csv(sleuth_table, "{fname_de_table}")'

    if dry_run:
        return "\n".join([lib_loading,
                          sleuth_load_design_str,
                          biomart_str,
                          prep_str,
                          sleuth_fit1_str,
                         sleuth_fit2_str,
                         sleuth_lrt_str,
                         sleuth_table_str,
                         sleuth_write_str])

    rpy2.robjects.r(lib_loading)
    rpy2.robjects.r(sleuth_load_design_str)

    print('biomart')
    rpy2.robjects.r(biomart_str)

    print('sleuth prep')
    rpy2.robjects.r(prep_str)

    print('sleuth fit')
    rpy2.robjects.r(sleuth_fit1_str)
    rpy2.robjects.r(sleuth_fit2_str)
    rpy2.robjects.r(sleuth_lrt_str)
    rpy2.robjects.r(sleuth_table_str)
    rpy2.robjects.r(sleuth_write_str)

    print('done')

    if outfile:
        rpy2.robjects.r(f'sleuth_save(so, "{outfile}")')


    return pd.read_csv(fname_de_table)


# for live shiny app:
# sleuth_live(so, options=options(browser = "kfmclient newTab"))
# or options(browser='/usr/bin/brave')
