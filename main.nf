#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_main } from "./nevermore/workflows/nevermore"
include { nevermore_align } from "./nevermore/workflows/align"
include { gffquant_flow } from "./nevermore/workflows/gffquant"
include { fastq_input } from "./nevermore/workflows/input"
include { collate_stats } from "./nevermore/modules/stats"

include { bowtie2_build; bowtie2_align } from "./nevermore/modules/align/bowtie2"


process scims {
	conda "/scratch/schudoma/envs/scims"
	tag "${sample.id}"
    cpus 1
    time {4.h * task.attempt}

	input:
	tuple val(sample), path(index_stats)
	tuple val(homogametic_chr), val(heterogametic_chr), path(scaffolds)

	output:
	tuple val(sample), path("${sample.id}/scims/*_results.txt"), emit: results
	// work/1e/0e6c5d680b3c8d1fc4a4a5796d7052/HP51/scims/HP51_results.txt

	script:
	"""
	set -e -o pipefail

	mkdir -p ${sample.id}/scims

	scims call \
    --idxstats_file ${index_stats} \
    --scaffolds ${scaffolds} \
    --homogametic_id ${homogametic_chr} \
    --heterogametic_id ${heterogametic_chr} \
    --output_dir ${sample.id}/scims
	"""
}

process scims_collate {
	tag "Collating..."
	executor "local"
	publishDir "${params.output_dir}", mode: "copy"

	input:
	path(scims_results)

	output:
	path("scims_collated.txt")

	script:
	"""
	awk -v OFS='\\t' 'NR==1 || NFR>1 { print \$0; }' ${scims_results} > scims_collated.txt
	"""
}


process samtools_idxstats {
	container "registry.git.embl.org/schudoma/bowtie2-docker:latest"
	tag "${sample.id}"
    cpus 1
    time {4.h * task.attempt}

	input:
	tuple val(sample), path(bam), path(index)
	
	output:
	tuple val(sample), path("${sample.id}/idxstats/${sample.id}.idxstats.txt"), emit: stats

	script:
	"""
	set -e -o pipefail

	mkdir -p ${sample.id}/idxstats/

	samtools idxstats ${bam} > ${sample.id}/idxstats/${sample.id}.idxstats.txt
	"""
}


workflow {
	def input_dir = (params.input_dir) ? params.input_dir : params.remote_input_dir
	def do_alignment = params.run_gffquant || !params.skip_alignment
	def do_stream = params.gq_stream
	def do_preprocessing = (!params.skip_preprocessing || params.run_preprocessing)

	params.ignore_dirs = ""

	print "PARAMS-MAIN: ${params}"

	fastq_input(
		Channel.fromPath(input_dir + "/**"),
		Channel.of(null)
	)

	fastq_ch = fastq_input.out.fastqs
	
	bowtie2_build(
		Channel.fromPath(params.reference_fasta).map { file -> [ [id: "reference"], file ] }
	)
	nevermore_main(fastq_ch)

	nevermore_main.out.fastqs.dump(pretty: true, tag: "fastqs_ch")

	align_ch = nevermore_main.out.fastqs
		.filter { sample, files -> sample.is_paired }
		.combine(bowtie2_build.out.index.map { _sample, index -> [ index ] })
	
	align_ch.dump(pretty: true, tag: "align_ch")

	counts_ch = nevermore_main.out.readcounts

	bowtie2_align(align_ch)

	samtools_idxstats(
		bowtie2_align.out.bam.join(bowtie2_align.out.bai, by: 0)
	)

	scims_ch = Channel.fromPath(params.scaffolds)
		.map { file -> [ params.homogametic_chr, params.heterogametic_ch, file ] }
	// tuple val(sample), path(index_stats)
	// tuple val(homogametic_chr), val(heterogametic_chr), path(scaffolds)

	scims(
		samtools_idxstats.out.stats,
		[ params.homogametic_chr, params.heterogametic_ch, params.scaffolds ]
	)

	scims_collate(scims.out.results.map { sample, file -> file }.collect())

	if (!do_stream && do_alignment) {
		nevermore_align(nevermore_main.out.fastqs)
		align_ch = nevermore_align.out.alignments
		counts_ch = counts_ch.mix(
			nevermore_align.out.aln_counts
				.map { sample, file -> return file }
				.collect()
		)
	}

	if (do_preprocessing && params.run_qa) {
		collate_stats(counts_ch.collect())		
	}

	if (params.run_gffquant) {

		if (params.gq_stream) {
			gq_input_ch = nevermore_main.out.fastqs
				.map { sample, fastqs ->
				sample_id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
				return tuple(sample_id, [fastqs].flatten())
			}
			.groupTuple()
			.map { sample_id, fastqs -> return tuple(sample_id, [fastqs].flatten()) }
			gq_input_ch.dump(pretty: true, tag: "gq_input_ch")

		} else {

			gq_input_ch = align_ch

		}

		gffquant_flow(gq_input_ch)		

	}

}
