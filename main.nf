#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { Basecall } from './modules.nf'
include { Align_transcriptome } from './modules.nf'
include { Nanopolish } from './modules.nf'
params.generate_db = false

    workflow{
        if ( params.dataprep ){
        Generate_ch = Channel
            .fromPath(params.input_csv)
            // can't get parser to work with headers
            // workaround for now 
            .splitCsv(header: false, skip:1)
            .map { row -> [row[0], row[1], row[2],file(row[3]), row[4], row[5]] }
        //Generate_ch.view()
        Basecall( 
            Generate_ch
        )
        Align_transcriptome(
            Basecall.out[1],
            file(params.transcriptome)
        )
                    Basecall.out[1].join(Basecall.out[2]).join(Align_transcriptome.out).groupTuple().view()
        Nanopolish(
            Basecall.out[1].join(Basecall.out[2]).join(Align_transcriptome.out).groupTuple(),
            file(params.transcriptome)


        )

        }
    }