#!/usr/bin/env nextflow

// This version of rnaseq-nf shows how it would be written
// in Nextflow 26.04, including the following new features:
//
// - record types
// - dataflow v2
// - typed processes (v2)
// - strict type checking
//
// NOTES:
//
// - records are just maps, and record types are applied to maps
//   in specific contexts to ensure a minimum schema. tuples are
//   replaced entirely by maps
//
// - the `reads` param is now a list of records instead of a file.
//   the input file should be a JSON or YAML file and is
//   automatically parsed into a list of records
//
// - dataflow v2 provides a new operator library, which is a streamlined
//   version of the v1 operators. workflow logic may need to be changed
//   significantly, depending on how the workflow is written. in this case,
//   the workflow logic is mostly the same
//
// - in dataflow v2, processes are called like an operator closure, typically
//   using the `map` operator. this makes the syntax more consistent and
//   removes the need for syntactic variants such as `|`, `&`, and `.out`
//
// - in dataflow v2, processes can only accept a single input and emit a single
//   output. multiple inputs/outputs are treated as a single record. when calling
//   a process, additional inputs can be supplied as named args
//
// - when type checking is set to 'strict' in the language server / compiler,
//   an error will be reported for any value that can't be type-checked
//

nextflow.enable.dataflow = 'v2'

params {
    // The input read-pair files
    reads: List<FastqPair>

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
        reads        : ${params.reads*.id}
        outdir       : ${workflow.outputDir}
    """.stripIndent()

    inputs_ch = channel.fromList(params.reads)

    (samples_ch, index) = RNASEQ( params.transcriptome, inputs_ch )

    multiqc_files_ch = samples_ch
        .flatMap { sample -> [sample.fastqc, sample.quant] }
        .collect()

    multiqc_report = MULTIQC( multiqc_files_ch, config: params.multiqc )

    versions_ch = channel.topic('versions')
        .unique()
        .map { t: ToolVersion ->
            [
                (t.process.tokenize(':').last()): [
                    (t.name): t.version
                ]
            ]
        }

    publish:
    samples = samples_ch
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
    samples: Channel<Sample> {
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

record ToolVersion {
    process : String
    name    : String
    version : String
}


workflow RNASEQ {
    take:
    transcriptome   : Path
    reads_ch        : Channel<FastqPair>

    main:
    index = INDEX(transcriptome)                    // Future<Path>
    fastqc_ch = reads_ch.map(FASTQC)                // Channel<id: String, fastqc: Path>
    quant_ch = reads_ch.map(QUANT, index: index)    // Channel<id: String, quant: Path>
    samples_ch = fastqc_ch.join(quant_ch, 'id')     // Channel<id: String, fastqc: Path, quant: Path>

    emit:
    samples : Channel<Sample> = samples_ch
    index   : Value<Path> = index
}

record FastqPair {
  id      : String
  fastq_1 : Path
  fastq_2 : Path
}

record Sample {
  id      : String
  fastqc  : Path
  quant   : Path
}


process FASTQC {
    tag "$id"
    conda 'bioconda::fastqc=0.12.1'

    input:
    id      : String
    fastq_1 : Path
    fastq_2 : Path

    output:
    id      : String = id
    fastqc  : Path = file("fastqc_${id}_logs")

    topic:
    [process: task.process, name: 'fastqc', version: '0.12.1'] >> 'versions'

    script:
    """
    fastqc.sh "$id" "$fastq_1 $fastq_2"
    """
}


process INDEX {
    tag "$transcriptome.simpleName"
    conda 'bioconda::salmon=1.10.3'

    input:
    transcriptome: Path

    output:
    file('index')

    topic:
    [process: task.process, name: 'salmon', version: '1.10.3'] >> 'versions'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}


process MULTIQC {
    conda 'bioconda::multiqc=1.27.1'

    input:
    inputs  : Bag<Path>
    config  : Path

    output:
    file('multiqc_report.html')

    topic:
    [process: task.process, name: 'multiqc', version: '1.27.1'] >> 'versions'

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
    id      : String = id
    quant   : QuantResult = quantResult(file("quant_${id}"))

    topic:
    [process: task.process, name: 'salmon', version: '1.10.3'] >> 'versions'

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${fastq_1} -2 ${fastq_2} -o quant_$id
    """
}


record QuantResult {
    aux_info    : Path
    libParams   : Path
    logs        : Path
    cmd_info    : Path
    counts      : Path
    quant       : Path
}


def quantResult(root: Path) -> QuantResult {
    record(
        aux_info:   root / 'aux_info',
        libParams:  root / 'libParams',
        logs:       root / 'logs',
        cmd_info:   root / 'cmd_info.json',
        counts:     root / 'lib_format_counts.json',
        quant:      root / 'quant.sf',
    )
}
