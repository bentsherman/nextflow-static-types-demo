#!/usr/bin/env nextflow 

// This version of rnaseq-nf shows how it would be written
// in Nextflow 25.10, including the following new features:
//
// - workflow params
// - workflow outputs (final)
// - worklfow `onComplete:` and `onError:` section
// - type annotations
// - typed process inputs/outputs
//
// NOTES:
//
// - type checking is enabled in the runtime / language server,
//   not as a feature flag
//
// - typed process inputs/outputs can be (mostly) automatically generated
//   from the existing syntax (e.g. as LSP command)
//
// - some workflow logic must be adjusted to satisfy typed processes,
//   which only accept one input channel (e,g, RNASEQ:QUANT)
//
// - in this case, Nextflow should be able to type-check the entire pipeline,
//   but it might not be able to do this in all cases, depending on which
//   operators are used
//

params {
    // The input read-pair files
    reads: Path

    // The input transcriptome file
    transcriptome: Path

    // Directory containing multiqc configuration
    multiqc: Path = "${projectDir}/multiqc"
}

workflow {
    main:
    log.info """\
        R N A S E Q - N F   P I P E L I N E
        ===================================
        transcriptome: ${params.transcriptome}
        reads        : ${params.reads}
        outdir       : ${workflow.outputDir}
    """.stripIndent()

    inputs_ch = channel.of(params.reads)
        .splitCsv()
        .map { id, fastq_1, fastq_2 ->
            tuple(id, file(fastq_1, checkIfExists: true), file(fastq_2, checkIfExists: true))
        }

    (samples_ch, index) = RNASEQ( params.transcriptome, inputs_ch )

    multiqc_files_ch = samples_ch
        .flatMap { id, fastqc, quant -> [fastqc, quant] }
        .collect()

    multiqc_report = MULTIQC( multiqc_files_ch.combine(channel.value(params.multiqc)) )

    versions_ch = channel.topic('versions')
        .unique()
        .map { process, name, version ->
            [
                (process.tokenize(':').last()): [
                    (name): version
                ]
            ]
        }

    publish:
    samples = samples_ch.map { id, fastqc, quant -> [id: id, fastqc: fastqc, quant: quant] }
    index = index
    multiqc_report = multiqc_report
    versions = versions_ch

    onComplete:
    log.info(
        workflow.success
            ? "\nDone! Open the following report in your browser --> ${workflow.outputDir}/multiqc_report.html\n"
            : "Oops .. something went wrong"
    )
}

output {
    samples: Channel<Map> {
        path { sample ->
            sample.fastqc >> "fastqc/${sample.id}"
            sample.quant >> "quant/${sample.id}"
        }
        index {
            path 'samples.json'
        }
    }

    index: Path {
        path '.'
    }

    multiqc_report: Path {
        path '.'
    }

    versions: Map<String,Map> {
        path '.'
        index {
            path 'versions.yml'
        }
    }
}


workflow RNASEQ {
    take:
    transcriptome   : Path
    reads_ch        : Channel<String,Path,Path>

    main:
    index = INDEX( channel.value(transcriptome) )   // Channel<Path>
    fastqc_ch = FASTQC(reads_ch)                    // Channel<String,Path>
    quant_ch = QUANT(index.combine(reads_ch))       // Channel<String,Path>
    samples_ch = fastqc_ch.join(quant_ch)           // Channel<String,Path,Path>

    emit:
    samples : Channel<String,Path,Path> = samples_ch
    index   : Path = index
}


process FASTQC {
    tag "$id"
    conda 'bioconda::fastqc=0.12.1'

    input:
    id      : String
    fastq_1 : Path
    fastq_2 : Path

    output:
    tuple(id, file("fastqc_${id}_logs"))

    topic:
    tuple("${task.process}", 'fastqc', '0.12.1') >> 'versions'

    script:
    """
    fastqc.sh "$id" "$fastq_1 $fastq_2"
    """
}


process INDEX {
    tag "$transcriptome.simpleName"
    conda 'bioconda::salmon=1.10.3'
    
    input:
    transcriptome   : Path

    output:
    result: Path = file('index')

    topic:
    tuple("${task.process}", 'salmon', '1.10.3') >> 'versions'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}


process MULTIQC {
    conda 'bioconda::multiqc=1.27.1'

    input:
    $in1    : Bag<Path>
    config  : Path

    stage:
    stageAs $in1, '*'

    output:
    report: Path = file('multiqc_report.html')

    topic:
    tuple("${task.process}", 'multiqc', '1.27.1') >> 'versions'

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
    index   : Path
    id      : String
    fastq_1 : Path
    fastq_2 : Path

    output:
    result = tuple(id, file("quant_${id}"))

    topic:
    tuple("${task.process}", 'salmon', '1.10.3') >> 'versions'

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${fastq_1} -2 ${fastq_2} -o quant_$id
    """
}
