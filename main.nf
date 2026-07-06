#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_main } from "./nevermore/workflows/nevermore"
include { nevermore_align } from "./nevermore/workflows/align"
include { gffquant_flow } from "./nevermore/workflows/gffquant"
include { fastq_input } from "./nevermore/workflows/input"
include { collate_stats } from "./nevermore/modules/stats"

include { bowtie2_build } from "./nevermore/modules/align/bowtie2"


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
		Channel.fromPath(params.reference_fasta).map { file -> [ "reference", file ] }
	)
	nevermore_main(fastq_ch)

	nevermore_main.out.fastqs.dump(pretty: true, tag: "fastqs_ch")

	align_ch = nevermore_main.out.fastqs
		.filter { sample, files -> sample.is_paired }
	
	align_ch.dump(pretty: true, tag: "align_ch")

	counts_ch = nevermore_main.out.readcounts






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
