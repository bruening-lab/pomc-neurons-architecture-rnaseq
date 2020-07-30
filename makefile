SHELL=/bin/bash
current_date := $(shell date +'%Y-%m-%d_%H-%M')
rpackage=bruening.2018.06.schumacher.riboseq
SINGULARITY=singularity exec /beegfs/scratch/bruening_scratch/pklemm/singularity/singularity-images/mytidyverse-base-3.6.3-2.simg

# Copy salmon data to package
custom_scripts/${rpackage}/inst/extdata/salmon_merged_gene_counts.csv.gz:
	mkdir -p custom_scripts/${rpackage}/inst/extdata/
	cp nfcore-rnaseq-pipeline/results/salmon/salmon_merged_gene_counts.csv custom_scripts/${rpackage}/inst/extdata/
	gzip custom_scripts/${rpackage}/inst/extdata/salmon_merged_gene_counts.csv
	# Rebuild the project to add the file
	cd custom_scripts/${rpackage} && $(SINGULARITY) Rscript -e 'roxygen2::roxygenise(); devtools::install()'

.PHONY: docs/nfcore_rnaseq_qc
docs/nfcore_rnaseq_qc:
	mkdir -p docs/nfcore_rnaseq_qc
	cp nfcore-rnaseq-pipeline/results/MultiQC/multiqc_report.html docs/nfcore_rnaseq_qc/${current_date}_multiqc.html

.PHONY: normalized_differential_expression
normalized_differential_expression:
	-rm -r release/deseq_normalized
	-mkdir -p release/deseq_normalized
	cd custom_scripts/${rpackage} && $(SINGULARITY) Rscript data-raw/create_data.R
	mv custom_scripts/${rpackage}/deseq_diff release/deseq_normalized
	mv custom_scripts/${rpackage}/differential_expression.nb.html release/deseq_normalized
	mv custom_scripts/${rpackage}/deseq_figures release/deseq_normalized
	zip -r release/deseq_normalized/deseq_figures.zip release/deseq_normalized/deseq_figures
	-mkdir -p docs/deseq_normalized
	cp release/deseq_normalized/differential_expression.nb.html docs/deseq_normalized/${current_date}_deseq_normalized.html
	# GO-terms
	-mkdir -p release/deseq_normalized/goterms
	mv custom_scripts/${rpackage}/glp1r_vs_lepr release/deseq_normalized/goterms
	zip -r release/deseq_normalized/goterms/glp1r_vs_lepr/goterm_plots.zip release/deseq_normalized/goterms/glp1r_vs_lepr/plots
	-mkdir -p docs/deseq_normalized/goterms
	cp release/deseq_normalized/goterms/glp1r_vs_lepr/goterm_report.nb.html docs/deseq_normalized/goterms/${current_date}_glp1r_vs_lepr.html
	# Rebuild the package
	cd custom_scripts/${rpackage} && $(SINGULARITY) Rscript -e 'roxygen2::roxygenise(); devtools::install()'
	# Export raw and normalized counts
	cd custom_scripts/${rpackage} && $(SINGULARITY) Rscript -e 'library(magrittr); ${rpackage}::gene_norm_de %>% dplyr::select(gene_id, external_gene_name, glp1r_1:lepr_pd_5) %>% WriteXLS::WriteXLS(ExcelFileName = "${shell readlink -f .}/release/deseq_normalized/counts.xlsx", AdjWidth = TRUE, AutoFilter = TRUE, BoldHeaderRow = TRUE, FreezeRow = 1, SheetNames = "counts")'
