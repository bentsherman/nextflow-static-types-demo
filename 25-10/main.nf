#!/usr/bin/env nextflow

// This version of rnaseq-nf shows how it would be written
// in Nextflow 25.10, including the following new features:
//
// - workflow params
// - workflow outputs (final)
// - worklfow `onComplete:` and `onError:` section
// - type annotations
// - typed processes
//
// NOTES:
//
// - type checking is enabled in the runtime / language server,
//   not as a feature flag
//
// - typed process inputs/outputs can be automatically migrated
//   from the existing syntax (e.g. as an LSP command)
//
// - in this case, Nextflow should be able to type-check the entire pipeline,
//   but it might not be able to do this in all cases, depending on which
//   operators are used
//

nextflow.preview.types = true

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

    read_pairs_ch = channel.of(params.reads)
        .flatMap { csv ->
            csv.splitCsv(header: true, sep: ',')
        }
        .map { row -> row as List<String> }
        .map { row ->
            tuple(row[0], file(row[1], checkIfExists: true), file(row[2], checkIfExists: true))
        }

    (samples_ch, index) = RNASEQ( read_pairs_ch, params.transcriptome )

    multiqc_files_ch = samples_ch
        .flatMap { id, fastqc, quant -> [fastqc, quant] }
        .collect()

    multiqc_report = MULTIQC( multiqc_files_ch, params.multiqc )

    versions_ch = channel.topic('versions')
        .map { entry -> entry as Tuple<String,String,String> }
        .reduce([:]) { acc, entry ->
            def (process, name, version) = entry
            acc[process.tokenize(':').last()] = [
                (name): version
            ]
            return acc
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
    reads_ch        : Channel<Tuple<String,Path,Path>>
    transcriptome   : Path

    main:
    index = INDEX(transcriptome)            // Value<Path>
    fastqc_ch = FASTQC(reads_ch)            // Channel<Tuple<String,Path>>
    quant_ch = QUANT(reads_ch, index)       // Channel<Tuple<String,Path>>
    samples_ch = fastqc_ch.join(quant_ch)   // Channel<Tuple<String,Path,Path>>

    emit:
    samples : Channel<Tuple<String,Path,Path>> = samples_ch
    index   : Value<Path> = index
}


process FASTQC {
    tag "$id"
    conda 'bioconda::fastqc=0.12.1'

    input:
    (id, fastq_1, fastq_2): Tuple<String, Path, Path>

    output:
    tuple(id, file("fastqc_${id}_logs"))

    topic:
    tuple(task.process, 'fastqc', '0.12.1') >> 'versions'

    script:
    """
    fastqc.sh "${id}" "${fastq_1} ${fastq_2}"
    """
}


process INDEX {
    tag transcriptome.simpleName
    conda 'bioconda::salmon=1.10.3'

    input:
    transcriptome: Path

    output:
    file('index')

    topic:
    tuple(task.process, 'salmon', '1.10.3') >> 'versions'

    script:
    """
    salmon index --threads ${task.cpus} -t ${transcriptome} -i index
    """
}


process MULTIQC {
    conda 'bioconda::multiqc=1.27.1'

    input:
    $in1    : Bag<Path>
    config  : Path

    stage:
    stageAs '*', $in1

    output:
    file('multiqc_report.html')

    topic:
    tuple(task.process, 'multiqc', '1.27.1') >> 'versions'

    script:
    """
    cp ${config}/* .
    echo "custom_logo: \$PWD/nextflow_logo.png" >> multiqc_config.yaml
    multiqc -n multiqc_report.html .
    """
}


process QUANT {
    tag "$id"
    conda 'bioconda::salmon=1.10.3'

    input:
    (id, fastq_1, fastq_2): Tuple<String, Path, Path>
    index: Path

    output:
    tuple(id, file("quant_${id}"))

    topic:
    tuple(task.process, 'salmon', '1.10.3') >> 'versions'

    script:
    """
    salmon quant --threads ${task.cpus} --libType=U -i ${index} -1 ${fastq_1} -2 ${fastq_2} -o quant_${id}
    """
}
