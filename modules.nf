process Basecall { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/basecall/${base}", mode: 'symlink', overwrite: true
container  "genomicpariscentre/guppy-gpu:latest"
cpus 8
beforeScript 'chmod o+rw .'
label (params.GPU == "ON" ? 'with_gpus': 'with_cpus')
input: 
    tuple val(base), val(condition), val(replicated), file(fast5_dir), val(flowcell), val(kit), file(summary)

output: 
    tuple val(base), file("${base}.basecalled_output/")
    tuple val(base), file("${base}.pass.fastq.gz")
    tuple val(base), file("${fast5_dir}"), file ("${summary}")

script:
"""
#!/bin/bash

ls -lah 

/usr/bin/guppy_basecaller --input_path ${fast5_dir} \
    --save_path ${base}.basecalled_output \
    --device "auto" \
    --flowcell ${flowcell} \
    --kit ${kit} \
    --compress_fastq \
    -r 

cat ${base}.basecalled_output/pass/*.fastq.gz > ${base}.pass.fastq.gz
"""
}

process Align_transcriptome { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/transcriptome_alignment/${base}", mode: 'symlink', overwrite: true
container  "vpeddu/nanopore_metagenomics:latest"
cpus 8
beforeScript 'chmod o+rw .'
input: 
    tuple val(base), file(fastq)
    file transcriptome

output: 
    tuple val(base), file("${base}.sorted.bam"), file("${base}.sorted.bam.bai")

script:
"""
#!/bin/bash

ls -lah 

minimap2 -ax map-ont \
    -uf \
    -t ${task.cpus} \
    --secondary=no \
    ${transcriptome} \
    ${fastq} | samtools view -Sb - | samtools sort -o ${base}.sorted.bam 

samtools index ${base}.sorted.bam
"""
}

process Nanopolish { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/transcriptome_alignment/${base}", mode: 'symlink', overwrite: true
container  "ttubb/nanopolish:latest"
cpus 8
beforeScript 'chmod o+rw .'
input: 
    tuple val(base), file(fastq), file(fast5_dir), file(summary), file(bam), file(bamindex)
    file transcriptome

output: 

script:
"""
#!/bin/bash

echo running ${base}

ls -lah 

nanopolish index -d ${fast5_dir} ${fastq}


nanopolish eventalign --reads ${fastq} \
--bam ${bam} \
--genome ${transcriptome} \
--signal-index \
--scale-events \
--summary ${summary} \
--threads 32 > ${base}.eventalign.txt


"""
}