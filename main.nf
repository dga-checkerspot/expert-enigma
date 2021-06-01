#!/usr/bin/env nextflow

params.reads='s3://algaetranscriptomics/CHK*_R{1,2}_001.fastq.gz'
pairInt='s3://transcriptomepipeline/PairInterleaves.sh'
genome='s3://hic.genome/PGA_scaffolds.fa'


Channel
	.fromFilePairs(params.reads)
	.ifEmpty {error "Cannot find any reads matching: ${params.reads}"}
	.set { read_pairs_ch }

process minimapS31 {

	input:
	tuple val(pair_id), path(reads) from read_pairs_ch
	path genom from genome
	
	output:
	set pair_id, "aln.sorted.bam" into align_ch
	

    """
	minimap2 -ax sr  $genom "${pair_id}_R1_001.fastq" "${pair_id}_R2_001.fastq" > aln.sam
	samtools view -bS aln.sam > aln.bam
	samtools sort aln.bam > aln.sorted.bam
	rm -f aln.sam
    """

}

process merge {

	memory '16G'

	input:
	tuple val(pair_id), path(fileList) from align_ch.collect()
	
	output:
	file 'merged.s.bam' into merged_bam
	
	"""
	echo `ls * ` > dir.txt
	samtools merge -b dir.txt merged.bam
	samtools sort merged.bam > merged.s.bam
	
	"""
}











