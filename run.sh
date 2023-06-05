 nextflow run main.nf \
 	--dataprep \
 	--transcriptome gencode.v39.transcripts.fa \
	--transcriptome_gtf gencode.v39.annotation.sorted.gtf \
 	--input_csv input_template.csv \
 	-with-singularity ubuntu:18.04 \
 	-resume
