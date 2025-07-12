#!/usr/bin/env nextflow

// This version of rnaseq-nf shows how it would be written
// in Nextflow 25.04, including the following new features:
//
// - topic channels (final)
// - workflow outputs (third preview)
//

nextflow.preview.output = true

params.reads = null
params.transcriptome = null
params.outdir = "results"
params.multiqc = "$projectDir/multiqc"

workflow {
    main:
    log.info """\
        R N A S E Q - N F   P I P E L I N E
        ===================================
        transcriptome: ${params.transcriptome}
        reads        : ${params.reads}
        outdir       : ${params.outdir}
    """.stripIndent()

    inputs_ch = channel.fromPath(params.reads)
        .splitCsv()
        .map { id, fastq_1, fastq_2 ->
            tuple(id, file(fastq_1, checkIfExists: true), file(fastq_2, checkIfExists: true))
        }

    (samples_ch, index) = RNASEQ( params.transcriptome, inputs_ch )

    multiqc_files_ch = samples_ch
        .flatMap { id, fastqc, quant -> [fastqc, quant] }
        .collect()

    multiqc_report = MULTIQC( multiqc_files_ch, params.multiqc ).report

    versions_ch = channel.topic('versions')
        .unique()
        .map { process, name, version ->
            [
                (process.tokenize(':').last()): [
                    (name): version
                ]
            ]
        }

    workflow.onComplete = {
        log.info(
            workflow.success
                ? "\nDone! Open the following report in your browser --> ${workflow.outputDir}/multiqc_report.html\n"
                : "Oops .. something went wrong"
        )
    }

    publish:
    samples = samples_ch.map { id, fastqc, quant -> [id: id, fastqc: fastqc, quant: quant] }
    index = index
    multiqc_report = multiqc_report
    versions = versions_ch
}

output {
    samples {
        path { sample ->
            sample.fastqc >> "fastqc/${sample.id}"
            sample.quant >> "quant/${sample.id}"
        }
        index {
            path 'samples.json'
        }
    }

    index {
        path '.'
    }

    multiqc_report {
        path '.'
    }

    versions {
        path '.'
        index {
            path 'versions.yml'
        }
    }
}


workflow RNASEQ {
    take:
    transcriptome
    reads_ch

    main:
    index = INDEX(transcriptome).result
    fastqc_ch = FASTQC(reads_ch).result
    quant_ch = QUANT(index, reads_ch).result
    samples_ch = fastqc_ch.join(quant_ch)

    emit:
    samples = samples_ch
    index = index
}


process FASTQC {
    tag "$id"
    conda 'bioconda::fastqc=0.12.1'

    input:
    tuple val(id), path(fastq_1), path(fastq_2)

    output:
    tuple val(id), path("fastqc_${id}_logs"), emit: result
    tuple val("${task.process}"), val('fastqc'), val('0.12.1'), topic: versions

    script:
    """
    fastqc.sh "$id" "$fastq_1 $fastq_2"
    """
}


process INDEX {
    tag "$transcriptome.simpleName"
    conda 'bioconda::salmon=1.10.3'

    input:
    path transcriptome

    output:
    path 'index', emit: result
    tuple val("${task.process}"), val('salmon'), val('1.10.3'), topic: versions

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}


process MULTIQC {
    conda 'bioconda::multiqc=1.27.1'

    input:
    path '*'
    path config

    output:
    path 'multiqc_report.html', emit: report
    tuple val("${task.process}"), val('multiqc'), val('1.27.1'), topic: versions

    script:
    """
    cp $config/* .
    echo "custom_logo: \$PWD/nextflow_logo.png" >> multiqc_config.yaml
    multiqc -n multiqc_report.html .
    """
}


process QUANT {
    tag "$id"
    conda 'bioconda::salmon=1.10.3'

    input:
    path index
    tuple val(id), path(fastq_1), path(fastq_2)

    output:
    tuple val(id), path("quant_${id}"), emit: result
    tuple val("${task.process}"), val('salmon'), val('1.10.3'), topic: versions

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${fastq_1} -2 ${fastq_2} -o quant_$id
    """
}
