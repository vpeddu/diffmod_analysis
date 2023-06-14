process Basecall { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/basecall/${base}", mode: 'symlink', overwrite: true
container  "genomicpariscentre/guppy-gpu:latest"
cpus 8
beforeScript 'chmod o+rw .'
label (params.GPU == "ON" ? 'with_gpus': 'with_cpus')
input: 
    tuple val(base), val(condition), val(replicated), path(fast5_dir), val(flowcell), val(kit), file(summary)

output: 
    tuple val(base), file("${base}.basecalled_output/")
    tuple val(base), file("${base}.pass.fastq.gz")
    tuple val(base), file("${fast5_dir}"), file ("${summary}")

script:
"""
#!/bin/bash

ls -lah 

/usr/bin/guppy_basecaller --input_path ${fast5_dir}/ \
    --save_path ${base}.basecalled_output \
    --device "cuda:0" \
    --flowcell ${flowcell} \
    --kit ${kit} \
    --compress_fastq \
    --recursive 

echo "finished basecalling ${base}"

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
    file transcriptome_fasta

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
    ${transcriptome_fasta} \
    ${fastq} | samtools view -Sb - | samtools sort -o ${base}.sorted.bam 

samtools index ${base}.sorted.bam
"""
}

process Nanopolish_index { 
//conda "${baseDir}/env/env.yml"
//publishDir "${params.output}/nanop/${base}", mode: 'symlink', overwrite: true
container  "quay.io/biocontainers/nanopolish:0.14.0--hb24e783_1"
cpus 8
beforeScript 'chmod o+rw .'
input: 
    tuple val(base), file(fastq), file(fast5_dir), file(summary), file(bam), file(bamindex)
    file transcriptome_fasta

output: 
    tuple val(base), file("${base}.pass.fastq.gz.index"), file("${base}.pass.fastq.gz.index.fai"), file("${base}.pass.fastq.gz.index.gzi"), file("${base}.pass.fastq.gz.index.readdb")

script:
"""
#!/bin/bash

echo running ${base}

ls -lah 

echo "starting index on ${base}"

nanopolish index -d ${fast5_dir} ${fastq}

echo "finished index for ${base}"

"""
}

process Nanopolish_eventalign { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/nanopolish_eventalign/${base}", mode: 'symlink', overwrite: true
container  "quay.io/biocontainers/nanopolish:0.14.0--h773013f_3"
cpus 16
beforeScript 'chmod o+rw .'
input: 
    tuple val(base), file(nanopolish_fqindex),file(nanopolish_fai),file(nanopolish_gzi),file(nanopolish_readbb), file(fastq), file(fast5_dir), file(summary), file(bam), file(bamindex)
    file transcriptome_fasta

output: 
    tuple val(base),file("${base}.eventalign.txt")
script:
"""
#!/bin/bash

echo running ${base}

ls -lah 

echo "starting eventalign for ${base} "

# not sure why, but this isn't defined in the ephemeral container
export HDF5_PLUGIN_PATH=/usr/local/hdf5/lib/plugin/

gunzip ${transcriptome_fasta}

nanopolish eventalign --reads ${fastq} \
--bam ${bam} \
--genome ${transcriptome_fasta} \
--signal-index \
--scale-events \
--summary ${summary} \
--threads ${task.cpus} > ${base}.eventalign.txt

"""
}

process Xpore_dataprep { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/xpore_dataprep/${base}", mode: 'symlink', overwrite: true
container  "yuukiiwa/xpore:2.1"
cpus 8
beforeScript 'chmod o+rw .'
input: 
    tuple val(base), file(eventalign)
    file transcriptome_fasta
    file transcriptome_gtf

output: 
    file "${base}.dataprep"
script:
"""
#!/bin/bash

echo running ${base}

ls -lah 

echo "starting dataprep on ${base}"

/usr/local/bin/xpore dataprep \
--eventalign ${eventalign} \
--gtf_or_gff ${transcriptome_gtf} \
--transcript_fasta ${transcriptome_fasta} \
--out_dir ${base}.dataprep \
--n_processes ${task.cpus}
"""
}


process Make_yaml { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/make_yaml/${base}", mode: 'symlink', overwrite: true
container  "yuukiiwa/xpore:2.1"
cpus 8
beforeScript 'chmod o+rw .'
input: 
    file input_csv

output: 
    file "xpore_config.yaml"

script:
"""
#!/usr/bin/env python
import os
import csv
import sys
import yaml

os.listdir()

def csv_to_yaml(csv_path, output_path, output_dir):
    data = {}
    with open(csv_path, 'r') as csv_file:
        reader = csv.reader(csv_file)
        next(reader)  # Skip header row
        for row in reader:
            sample_name = row[0]
            condition = row[1]
            replicate = row[2]
            data_path = str('./' + sample_name + '.dataprep')
            
            if condition not in data:
                data[condition] = {}
            
            data[condition][replicate] = data_path
    
    output = {'data': data, 'out': output_dir}
    
    with open(output_path, 'w') as yaml_file:
        yaml.dump(output, yaml_file, default_flow_style=False)

# Example usage
csv_file_path = "${input_csv}"
output_directory = './xpore_output'
yaml_output_path = 'xpore_config.yaml'
csv_to_yaml(csv_file_path, yaml_output_path, output_directory)
"""
}

process Xpore { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/xpore/${base}", mode: 'symlink', overwrite: true
container  "yuukiiwa/xpore:2.1"
cpus 8
beforeScript 'chmod o+rw .'
input: 
    file yaml
    file datapreps

output: 

script:
"""
#!/bin/bash

ls -lah 

echo "starting Xpore run"

/usr/local/bin/xpore diffmod \
    --config ${yaml} \
    --n_processes ${task.cpus}

/usr/local/bin/xpore postprocessing \
    --diffmod_dir xpore_output
"""
}

process M6anet_dataprep { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/m6anet_dataprep/${base}", mode: 'symlink', overwrite: true
container  "yuukiiwa/m6anet:1.0"
cpus 8
beforeScript 'chmod o+rw .'
input: 
    tuple val(base), file(eventalign)


output: 
    tuple val(base), file("${base}.dataprep")
script:
"""

m6anet-dataprep --eventalign ${eventalign} \\
                --out_dir ${base}.dataprep \\
                --n_processes 4

"""
}

process M6anet { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/m6anet/${base}", mode: 'symlink', overwrite: true
container  "yuukiiwa/m6anet:1.0"
cpus 8
beforeScript 'chmod o+rw .'
input: 
    tuple val(base),file(datapreps)

output: 
    file "${base}.m6anet.output"
script:
"""
m6anet-run_inference \\
    --input_dir ${datapreps} \\
    --out_dir ${base}.m6anet.output  \\
    --n_processes ${task.cpus} \\
    --num_iterations 1000
"""
}
