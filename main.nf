#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { Basecall } from './modules.nf'
include { Align_transcriptome } from './modules.nf'
include { Nanopolish_index } from './modules.nf'
include { Nanopolish_eventalign } from './modules.nf'
include { Xpore_dataprep } from './modules.nf'
include { Make_yaml } from './modules.nf'
include { Xpore } from './modules.nf'
include { M6anet_dataprep } from './modules.nf'
include { M6anet } from './modules.nf'

params.generate_db = false

    workflow{
        if ( params.dataprep ){
        Generate_ch = Channel
            .fromPath(params.input_csv)
            // can't get parser to work with headers
            // workaround for now 
            .splitCsv(header: false, skip:1)
            .map { row -> [row[0], row[1], row[2],row[3], row[4], row[5], file(row[6])] }
        //Generate_ch.view()
        Basecall( 
            Generate_ch
        ) 
        Align_transcriptome(
            Basecall.out[1],
            file(params.transcriptome)
        )
        //Basecall.out[1].join(Basecall.out[2]).join(Align_transcriptome.out).groupTuple().view()
        Nanopolish_index(
            Basecall.out[1].join(Basecall.out[2]).join(Align_transcriptome.out).groupTuple(),
            file(params.transcriptome)
        )
        Nanopolish_eventalign(
            Nanopolish_index.out.join(Basecall.out[1]).join(Basecall.out[2]).join(Align_transcriptome.out).groupTuple(),
            file(params.transcriptome)
        )
        Xpore_dataprep( 
            Nanopolish_eventalign.out,
            file(params.transcriptome),
            file(params.transcriptome_gtf)
        )
        Make_yaml( 
            file(params.input_csv)
        )
        Xpore(
            Make_yaml.out,
            Xpore_dataprep.out.collect()
        )
        M6anet_dataprep( 
            Nanopolish_eventalign.out,
        )
        M6anet( 
            M6anet_dataprep.out
        )
        }
    }