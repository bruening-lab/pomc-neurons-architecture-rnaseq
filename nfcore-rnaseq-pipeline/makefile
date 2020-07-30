SHELL=/bin/bash
SINGULARITY=singularity exec $(SINGULARITY_IMAGES)/mytidyverse-base-3.6.3-2.simg

all:
	clean \
	results

nextflow:
	curl -s https://get.nextflow.io | bash

results: nextflow mpisf-rnaseq.config raw_data
	# Run nextflow
	./nextflow run nf-core/rnaseq \
		-r 1.4.2 \
		--genome GRCm38 \
		--unStranded \
		--pseudo_aligner salmon \
		-c ${shell readlink -f mpisf-rnaseq.config} \
		--reads '${shell readlink -f raw_data}/*{_R1,_R2}.fastq.gz'

.PHONY: clean
clean:
	-rm -r raw_data
	-rm mpisf-rnaseq.config
	-rm nextflow
	-rm -r .nextflow
	-rm .nextflow.log*
	-rm -r results

.PHONY: clean_work_files
clean_work_files:
	rm -r work
