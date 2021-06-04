#!/usr/bin/env nextflow

params.reads='s3://algaetranscriptomics/CHK*_R{1,2}_001.fastq.gz'
pairInt='s3://transcriptomepipeline/PairInterleaves.sh'
genome='s3://hic.genome/PGA_scaffolds.fa'
genome2='s3://hic.genome/PGA_scaffolds.fa'
genome3='s3://hic.genome/PGA_scaffolds.fa'
protein='s3://hic.genome/GCF_000733215.1_ASM73321v1_protein.faa'


Channel
	.fromFilePairs(params.reads)
	.ifEmpty {error "Cannot find any reads matching: ${params.reads}"}
	.set { read_pairs_ch }


//First do genomethreader

process gth { 
	
	memory '8G'
	
	input:
	path genom from genome3
	path prot from protein
	
	output: 
	file 'bonafide.gtf' into protHints
	file 'bonafide.gb' into protTrain
	
	"""
	startAlign.pl --genome=$genom --prot=$prot --prg=gth
	gth2gtf.pl align_gth/gth.concat.aln bonafide.gtf
	gff2gbSmallDNA.pl bonafide.gtf genome.fa 10000 bonafide.gb
	"""

}



process Augustus {

	input:
	path genom from genome
	
	output:
	file 'aug.gtf' into abinitio
	
	"""
	augustus --species=chlorella $genom > aug.gtf
	
	"""
	

}


process STARALIGN {

	memory '16G'

	input:
	tuple val(pair_id), path(reads) from read_pairs_ch
	path genom from genome2
	path genes from abinitio
	
	output:
	set pair_id, "STARAligned.sortedByCoord.out.bam" into align_ch
	

    """
    gunzip $reads
    STAR --runMode genomeGenerate --genomeDir /opt --genomeFastaFiles $genom --sjdbGTFfile $genes --sjdbOverhang 99 
    STAR --genomeDir /opt --outFileNamePrefix /opt/STAR --outSAMtype BAM SortedByCoordinate --readFilesIn "${pair_id}_R1_001.fastq" "${pair_id}_R2_001.fastq" 

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







