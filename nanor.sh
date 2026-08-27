#!/usr/bin/bash
#PBS -S /usr/bin/bash
#PBS -q rgq
#PBS -N nanor
#PBS -l select=1:ncpus=8
#PBS -V
#PBS -k oe
#PBS -m abe
#PBS -M 

cd $PBS_O_WORKDIR/ 

## last updated 20240614
## Naohiro Makise, Masahito Kawazu
## 2 inputs
## qsub \
## -v "CASE=$CASE,SAMPLE=$SAMPLE" \
## nanord1.pbs

CURRENT_DIR=$(pwd)

CASE="$CASE"
SAMPLE="$SAMPLE"
DIR="$CURRENT_DIR/$CASE/$SAMPLE"


## parse d1
DAY=d1

OUTDIR="$DIR/$DAY"
PREF="$SAMPLE-$DAY-hacm4"
BAM="$PREF.sort.bam"

echo "CASE="$CASE
echo "SAMPLE="$SAMPLE
echo "DIR="$DIR
echo "OUTDIR="$OUTDIR 
echo "PREF="$PREF
echo "BAM="$BAM

date

cd $DIR

## if d1.bam file does not exist, exit 1
if [[ ! -e $DIR/$PREF.bam ]]; then
    echo $PREF".bam does not exist"
    echo "Stop analysis"
    exit 1
else
    echo $PREF".bam exists!"
fi

###########################

## if mamba cannot be activated, exit 1
source $CURRENT_DIR/.bashrc
micromamba activate clairs-to || echo "Failed to activate micromamba & stop anlysis"
micromamba activate clairs-to || exit 1
echo "micromamba can be activated!" && date
micromamba deactivate

###############################

##################################

## if conda cannot be activated, exit 1
source $CURRENT_DIR/miniconda3/etc/profile.d/conda.sh
conda activate || echo "Filed to activate conda & stop anlysis"
conda activate || exit 1
echo "conda can be activated!"
conda deactivate

################################

## samtools
#### sort bam
	echo "Start samtools" && date
	echo "Start sort $DAY" && date
	/opt/bin/samtools sort -o $PREF.sort.bam $PREF.bam
	echo "Finished sort $DAY" && date

#### index sorted bam
	echo "Start index $DAY" && date
	/opt/bin/samtools index $PREF.sort.bam
	echo "Finished index $DAY" && date

#### stats sorted bam
	echo "Start stats $DAY" && date
	/opt/bin/samtools stats $PREF.sort.bam > $PREF.stats.txt
	echo "Finished stats $DAY" && date

#### extract u1000 bam from sorted bam
	echo "Start extracting u1000 $DAY" && date
	/opt/bin/samtools view -h $BAM | \
	awk 'length($10) <= 1000 || $1 ~ /^@/' | \
	/opt/bin/samtools view -b -o $PREF-u1000.bam
	echo "Finished extracting u1000 $DAY!" && date

#### index u1000 bam
	echo "Start index u1000 $DAY" && date
	/opt/bin/samtools index $PREF-u1000.bam
	echo "Finished index u1000 $DAY" && date

#### stats u1000 bam
	echo "Start stats u1000 $DAY" && date
	/opt/bin/samtools stats $PREF-u1000.bam > $PREF-u1000.stats.txt
	echo "Finished stats u1000 $DAY" && date

echo "Finished samtools $DAY" && date


## parse
#### copy number
###### mosdepth
		mkdir -p $OUTDIR/mosdepth
		cd $OUTDIR/mosdepth
		conda activate mosdepth

		echo "Start mosdepth $DAY" && date
		mosdepth -n \
		--by $CURRENT_DIR/data/20230818_target_only_bed.bed \
		$PREF-gene \
		$DIR/$BAM

		mosdepth -n \
		--by $CURRENT_DIR/data/20230812-target-bed.bed \
		$PREF-target \
		$DIR/$BAM

		mosdepth -n \
		--by $CURRENT_DIR/data/20230812-off-target-bed.bed \
		$PREF-off-target \
		$DIR/$BAM

		conda deactivate
		echo "Finished mosdepth $DAY" && date
		cd $DIR

#### SV
###### SVIM
		echo "Start SVIM $DAY" && date
		conda activate svim

		svim alignment . $BAM $CURRENT_DIR/genomes/hg38_seq/hg38.fa \
		--verbose \
		--sequence_alleles \
		--insertion_sequences \
		--tandem_duplications_as_insertions \
		--read_names

		cd signatures
		cp trans.bed $PREF-trans.bed
		awk -F'\t' '$5 >= 2' trans.bed > $PREF-trans2.bed
		awk -F'\t' '$5 >= 3' trans.bed > $PREF-trans3.bed
		awk -F'\t' '$5 >= 4' trans.bed > $PREF-trans4.bed
		awk -F'\t' '$5 >= 5' trans.bed > $PREF-trans5.bed

		cd ..
		mv variants.vcf signatures/
		mv SVIM_* signatures/
		mv candidates/ signatures/candidates/
		mv signatures/ $OUTDIR/SVIM_signatures/
		conda deactivate
		echo "Finished svim $DAY" && date

#### SNV
## clairs-to
micromamba activate clairs-to

/rgdata/home/makise/ClairS-TO/run_clairs_to \
	-T $DIR/$BAM \
	-R $CURRENT_DIR/genomes/hg38_seq/hg38.fa \
	-o $OUTDIR/clairs-to-target \
	-t 8 \
	-p ont_r10_dorado_hac_4khz \
	-b $CURRENT_DIR/data/20230812-target-bed.bed

gzip -dc $OUTDIR/clairs-to-target/snv.vcf.gz > $OUTDIR/clairs-to-target/$PREF.CSTO.snv.vcf
gzip -dc $OUTDIR/clairs-to-target/indel.vcf.gz > $OUTDIR/clairs-to-target/$PREF.CSTO.indel.vcf

echo "Finished clairs-to target $DAY" && date
micromamba deactivate

###################################

###### snpsift
		echo "Start snpsift CSTO $DAY" && date
		conda activate snpsift
		cd $OUTDIR/clairs-to-target
		
		java -jar $CURRENT_DIR/.conda/pkgs/snpsift-5.2-hdfd78af_0/share/snpsift-5.2-0/SnpSift.jar annotate \
		$CURRENT_DIR/data/clinvar.vcf.gz \
		$PREF.CSTO.snv.vcf > $PREF.CSTO.snv.cl.vcf
		
		java -jar $CURRENT_DIR/.conda/pkgs/snpsift-5.2-hdfd78af_0/share/snpsift-5.2-0/SnpSift.jar annotate \
		$CURRENT_DIR/data/clinvar.vcf.gz \
		$PREF.CSTO.indel.vcf > $PREF.CSTO.indel.cl.vcf
		
		cat $PREF.CSTO.snv.cl.vcf | \
		java -jar $CURRENT_DIR/.conda/pkgs/snpsift-5.2-hdfd78af_0/share/snpsift-5.2-0/SnpSift.jar filter \
		"(CLNSIG =~ 'Pathogenic') | (CLNSIG =~ 'pathogenic') | (CLNSIG =~ 'Uncertain') | (CLNSIG =~ 'uncertain')" \
		> $PREF.CSTO.snv.cl.PU.vcf
		
		cat $PREF.CSTO.indel.cl.vcf | \
		java -jar $CURRENT_DIR/.conda/pkgs/snpsift-5.2-hdfd78af_0/share/snpsift-5.2-0/SnpSift.jar filter \
		"(CLNSIG =~ 'Pathogenic') | (CLNSIG =~ 'pathogenic') | (CLNSIG =~ 'Uncertain') | (CLNSIG =~ 'uncertain')" \
		> $PREF.CSTO.indel.cl.PU.vcf
		
		java -jar $CURRENT_DIR/.conda/pkgs/snpsift-5.2-hdfd78af_0/share/snpsift-5.2-0/SnpSift.jar annotate \
        $CURRENT_DIR/data/CosmicCodingMuts.vcf.gz \
        $PREF.CSTO.snv.cl.PU.vcf > $PREF.CSTO.snv.cl.PU.c.vcf
        
        java -jar $CURRENT_DIR/.conda/pkgs/snpsift-5.2-hdfd78af_0/share/snpsift-5.2-0/SnpSift.jar annotate \
        $CURRENT_DIR/data/CosmicCodingMuts.vcf.gz \
        $PREF.CSTO.indel.cl.PU.vcf > $PREF.CSTO.indel.cl.PU.c.vcf
        
		java -jar $CURRENT_DIR/.conda/pkgs/snpsift-5.2-hdfd78af_0/share/snpsift-5.2-0/SnpSift.jar annotate \
        $CURRENT_DIR/data/00-All.vcf.gz \
        $PREF.CSTO.snv.cl.PU.c.vcf > $PREF.CSTO.snv.cl.PU.cd.vcf
        
        java -jar $CURRENT_DIR/.conda/pkgs/snpsift-5.2-hdfd78af_0/share/snpsift-5.2-0/SnpSift.jar annotate \
        $CURRENT_DIR/data/00-All.vcf.gz \
        $PREF.CSTO.indel.cl.PU.c.vcf > $PREF.CSTO.indel.cl.PU.cd.vcf
		
		conda deactivate
		echo "Finished snpsift CSTO $DAY" && date

###### snpeff
		echo "Start snpEff CSTO $DAY" && date
		
		/opt/bin/snpEff GRCh38.p14 $PREF.CSTO.snv.cl.PU.cd.vcf > $PREF.CSTO.snv.cl.PU.cde.vcf
		/opt/bin/snpEff GRCh38.p14 $PREF.CSTO.indel.cl.PU.cd.vcf > $PREF.CSTO.indel.cl.PU.cde.vcf
		
		echo "FInished snpEff CSTO $DAY" && date
		
#################################
		
#### methyl
###### modkit
		cd $DIR
		echo "Start modkit $DAY" && date
		conda activate modkit

		modkit pileup \
		$BAM \
		$PREF-modkit.bed \
		--ref $CURRENT_DIR/genomes/hg38_seq/hg38.fa \
		--cpg \
		--only-tabs \
		--include-bed $CURRENT_DIR/data/20230827_target_cpg_only_bed.bed

###### modkit mgmt
		echo "Start modkit MGMT $DAY" && date
		
		modkit pileup \
		$BAM \
		$PREF-modkit-mgmt.bed \
		--ref $CURRENT_DIR/genomes/hg38_seq/hg38.fa \
		--cpg \
		--only-tabs \
		--include-bed $CURRENT_DIR/data/20240413_mgmt_cpg_only_bed.bed

		conda deactivate
		echo "Finished modkit $DAY" && date


## R
#### QDNAseq
	echo "Start QDNAseq on R $DAY" && date
	
	R --vanilla --slave --args $CASE $SAMPLE $DAY < $CURRENT_DIR/QDNAseq.R

	echo "Finished QDNAseq on R $DAY" && date

#### circlize
	echo "Start circlize on R $DAY" && date
	
	R --vanilla --slave --args $CASE $SAMPLE $DAY < $CURRENT_DIR/circlize.R

	echo "Finished circlize on R $DAY" && date

#### methyl classify sarc
	echo "Start mSARC on R $DAY" && date
	
	R --vanilla --slave --args $CASE $SAMPLE $DAY < $CURRENT_DIR/msarc.R

	echo "Finished mSARC on R $DAY" && date

echo #Finished all $DAY"

