#!/usr/bin/env nextflow

params.reads='s3://algaetranscriptomics/CHK*_R{1,2}_001.fastq.gz'
perlAnnot='s3://hic.genome/perlAnnotation.sh'
params.geno='s3://hic.genome/ragoo_primary.fasta'

params.prot='s3://hic.genome/*protein.faa'
params.cdna='s3://hic.genome/AWSBatch_transcriptome.fasta'

prot_datasets = Channel.fromPath(params.prot)
prot_datasets.into{Proteins; protein1}

cdna_datasets= Channel.fromPath(params.cdna)
cdna_datasets.into{cdnafile; cdnafile1}

geno=Channel.fromPath(params.geno)

Channel
	.fromFilePairs(params.reads)
	.ifEmpty {error "Cannot find any reads matching: ${params.reads}"}
	.set { read_pairs_ch }


//Do repeatmasker

process repeatMask {
	memory '4G'
	
	input:
	path genome from geno
	
	output:
	file "${genome.baseName}.fasta.masked" into maskedGenome
	
	"""
	RepeatMasker --species arabidopsis -xsmall $genome 
	"""
	
}

maskedGenome.into{genome; genome1 ; genome2 ; genome3 ; genome4 ; genome5 ; genome6}



//Do genomethreader

process gth { 
	
	memory '2G'
	
	input:
	path genom from genome3
	path prot from Proteins
	
	output: 
	file "${prot.baseName}.bonafide.gtf" into protHints
	
	"""
	startAlign.pl --genome=$genom --prot=$prot --prg=gth
	gth2gtf.pl align_gth/gth.concat.aln "${prot.baseName}.bonafide.gtf"
	"""

}

process gtHints {
	memory '2G'
	
	input:
	path gtf from protHints
	
	output:
	file "${gtf.baseName}.hints" into gtHintsFiles
	
	"""
	align2hints.pl --in=$gtf --out=${gtf.baseName}.hints --prg=gth
	"""
}

process gthMerge {
	memory '2G'
	
	input:
	path hints from gtHintsFiles.collect()
	
	output:
	file "gtHints.merge.hints" into mergeGtHints
	
	
	"""
	cat *.hints > gtHints.merge.hints
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

abinitio.into{abinitio1; abinitio2}


process STARALIGN {

	memory '4G'

	input:
	tuple val(pair_id), path(reads) from read_pairs_ch
	path genom from genome2
	path genes from abinitio1
	
	output:
	set pair_id, "STARAligned.sortedByCoord.out.bam" into align_ch
	

    """
    STAR --runMode genomeGenerate --genomeDir /opt --genomeFastaFiles $genom --sjdbGTFfile $genes --sjdbOverhang 99 --genomeSAindexNbases 10
    STAR --genomeDir /opt --outFileNamePrefix STAR --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --readFilesIn "${pair_id}_R1_001.fastq.gz" "${pair_id}_R2_001.fastq.gz" --limitBAMsortRAM 44000000000 

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


process wiggle { 

	memory '4G'
	
	input:
	path merged from merged_bam
	
	output:
	file 'rnaseq.gff' into rnaHints
	
	"""
	bam2wig $merged > rnaseq.wig
	cat rnaseq.wig | wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep --UCSC=unstranded.track --radius=4.5 --pri=4 --strand="." > rnaseq.gff
	"""
}

process cdna {

	memory '4G'
	
	input:
	path transcriptome from cdnafile
	path genome from genome5
	
	output:
	file 'cdna.hints' into cdnaHints
	
	
	"""
	blat -noHead -minIdentity=88 $genome $transcriptome blat_cdna.psl

	blat2hints.pl --in=blat_cdna.psl --out=cdna.hints --minintronlen=35
	"""
}


process AllHints {

	memory '4G'

	input:
	path cdna from cdnaHints
	path rna from rnaHints
	path prot from mergeGtHints
	
	output:
	file 'mergedAllHints.hints' into hintsFile
	
	"""
	cat $cdna $rna $prot > mergedAllHints.hints
	"""

}

process bonafide {

	memory '8G'
	
	input:
	path bam from merged_bam
	path genom from genome4
	path genes from abinitio2
	
	output:
	file 'bonafide.gtf' into RNA_gtf
	file 'bonafide.gb' into RNA_gb
	
	"""
	
	bam2hints --intronsonly --in=$bam --out=introns.gff

	filterIntronsFindStrand.pl $genom introns.gff --score > introns.f.gff

	filterGenemark.pl --genemark=$genes --introns=introns.f.gff

	ln -s aug.f.good.gtf bonafide.gtf

	gff2gbSmallDNA.pl $genes $genom 300 bonafide.gb


	"""
	
}


process runAnnotation {

	cpus 4
	memory '8G'
	
	input:
	path genome from genome6
	path bonafide from RNA_gb
	path hints from hintsFile
	path cdna from cdnafile1
	path script from perlAnnot
	
	output:
	file 'bug_optimized_hints.gff' into gff
	file 'bug.tar.gz' into config_dir
	
	"""
	AUGUSTUS_CONFIG_PATH=/root/miniconda3/config
	export AUGUSTUS_CONFIG_PATH
	new_species.pl --species=bug

	randomSplit.pl $bonafide 3500
	mv bonafide.gb.test test.gb
	mv bonafide.gb.train train.gb
	etraining --species=bug train.gb &> etrain.out
	augustus --species=bug test.gb > test.out
	
	optimize_augustus.pl --species=bug --rounds=2 --kfold=4 --cpus=4 --jg=1 train.gb > optimize.out
	etraining --species=bug train.gb &> etrain.out
	augustus --species=bug test.gb > test.opt.out
	
	blat -noHead -minIdentity=88 $genome $cdna blat_cdna.psl
	blat2hints.pl --in=blat_cdna.psl --out=cdna.hints --minintronlen=35
	augustus --species=chlamydomonas --extrinsicCfgFile=/root/miniconda3/config/extrinsic/extrinsic.E.cfg --hintsfile=cdna.hints --softmasking=off $genome > chlamy_CDNA_hints.gff
	
	chmod 777 $script
	./$script
	
	gff2gbSmallDNA.pl --good=supported.lst chlamy_CDNA_hints.gff $genome 800 chlamy_bonafide.gb
	optimize_augustus.pl --species=bug --rounds=3 chlamy_bonafide.gb --UTR=on --metapars=/root/miniconda3/config/species/bug/bug_metapars.utr.cfg --trainOnlyUtr=1

	augustus --species=bug --extrinsicCfgFile=/root/miniconda3/config/extrinsic/extrinsic.M.RM.E.W.P.cfg --hintsfile=$hints --softmasking=on --UTR=on --print_utr=on --alternatives-from-sampling=true --alternatives-from-evidence=true $genome > bug_optimized_hints.gff
	cp -r /root/miniconda3/config/species/bug ./bug/
	tar -zcvf bug.tar.gz bug
	"""

}





