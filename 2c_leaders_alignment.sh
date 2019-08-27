#! /usr/bin/env bash

#Bowtie used to align tenmers to hg19 with a 10bp window around 	
#ctss. Up to 20 potential locations are allowed for each 10mer.

#modules
#module add apps/gcc/BEDTools/2.25.0
#module add apps/gcc/samtools/1.6
#module add apps/gcc/bowtie2/2.2.6 



#create index 
#bowtie2-build tss_regions.fasta human_tss





for FILE in `find -type f -name '*.5.fasta'` ; do
	a=${FILE}
	sample=$(basename $FILE)
	b=${sample%.fasta}
	b+=.bam

	c=${sample%.fasta}
	c+=.sam
	d=${sample%.fasta}
	d+=.mapped.bam

	#bedtools bamtofastq -i $a -fq $b

	if [ -e $a ]; then
		echo $FILE
		echo fastq exists
		bowtie2 -k 20 -x human_tss_10 -f $FILE -S $c	
	fi

done


	