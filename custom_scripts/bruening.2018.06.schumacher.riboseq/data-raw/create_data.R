library(magrittr)

raw_counts <- system.file(
    "extdata/salmon_merged_gene_counts.csv.gz",
    package = "bruening.2018.06.schumacher.riboseq"
  ) %>%
  readr::read_csv()

gene_counts_norm <-
  raw_counts %>%
  # Filter mitochondrial genes
  dplyr::filter(!(gene_id %in% c("ENSMUSG00000064337", "ENSMUSG00000064339"))) %>%
  # Normalize the values
  dplyr::mutate(
    glp1r_1 = (glp1r_pd_1 + 1) / (glp1r_input_1 + 1),
    glp1r_2 = (glp1r_pd_2 + 1) / (glp1r_input_2 + 1),
    glp1r_3 = (glp1r_pd_3 + 1) / (glp1r_input_3 + 1),
    lepr_2 = (lepr_pd_2 + 1) / (lepr_input_2 + 1),
    lepr_3 = (lepr_pd_3 + 1) / (lepr_input_3 + 1),
    lepr_4 = (lepr_pd_4 + 1) / (lepr_input_4 + 1),
    lepr_5 = (lepr_pd_5 + 1) / (lepr_input_5 + 1)
  ) %>%
  dplyr::select(gene_id, glp1r_1:lepr_5) %>%
  # Replace nan values with 0
  dplyr::mutate_all(~ifelse(is.nan(.), 0, .)) %>%
  # Multiply with constant factor to have high counts
  dplyr::mutate_if(is.numeric, ~ . * 1000) %>%
  # Round values to comply with DESeq2 requirements
  dplyr::mutate_if(is.numeric, round)

# Set WD here, because when we do it in the function call, it uses the wd from rmarkdown
my_wd <- getwd()
nfRNAseqDESeq2::run_differential_expression(
  path_config_json = file.path(my_wd, "data-raw/groups.json"),
  count_data = gene_counts_norm,
  out_path = my_wd,
  biomart_attributes = c("external_gene_name", "gene_biotype")
)

gene_norm_de <- readr::read_csv(file.path(getwd(), "deseq_diff/deseq2_diff.csv")) %>%
  dplyr::rename(gene_id = row) %>%
  # Attach normalized counts
  tidylog::left_join(gene_counts_norm) %>%
  # Attach original counts
  tidylog::left_join(
    raw_counts %>%
      # Alphabetically order the columns
      dplyr::select(order(colnames(.)))
  )

usethis::use_data(gene_norm_de, overwrite = TRUE)

nfRNAseqDESeq2::goterm_analysis_of_all_comparisons(
  deseq2_diff_path = file.path(my_wd, "deseq_diff", "deseq2_diff.csv"),
  out_path = my_wd,
  simplify_ontologies = TRUE,
  do_gse = FALSE
)
