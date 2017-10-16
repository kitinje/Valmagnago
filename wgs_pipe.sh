# Load modules
module load autoload bwa
module load autoload gatk
module load autoload r

# Environmental variables
Fastq_R1=/pico/scratch/userexternal/mfratell/CORNELIA/analisi/fastq/MADRE_S2_R1_001.fastq
Fastq_R2=/pico/scratch/userexternal/mfratell/CORNELIA/analisi/fastq/MADRE_S2_R2_001.fastq
hg38_fasta=/pico/scratch/userexternal/mfratell/GENOME/iGENOMES/Homo_sapiens/hg38/Sequence/BWAIndex/genome.fa

bwa mem -M -R '@RG\tID:foo\tSM:bar' -t 20 $hg38_fasta $Fastq_R1 $Fastq_R2 > madre_bwa.sam

samtools view -bS -@ 20 madre_bwa.sam | samtools sort -@ 20 -o madre_bwa_sorted.bam
samtools index sofia_sorted.bam

samtools view -@ 20 -b sofia_sorted.bam "chr5:36544862-37398309" > nipbl.bam
samtools sort -@ 20 nipbl.bam -o nipbl_sorted.bam
samtools index nipbl_sorted.bam

GATK_PATH=/pico/scratch/userexternal/mfratell/SW/GATK38/

indels=/pico/scratch/userexternal/mfratell/GENOME/gatk_hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf
dbSNP=/pico/scratch/userexternal/mfratell/GENOME/gatk_hg38/dbsnp150.vcf
G1000_SNP=/pico/scratch/userexternal/mfratell/GENOME/gatk_hg38/1000G_phase1.snps.high_confidence.hg38.vcf
omni_SNP=/pico/scratch/userexternal/mfratell/GENOME/gatk_hg38/1000G_omni2.5.hg38.vcf
hapmap_SNP=/pico/scratch/userexternal/mfratell/GENOME/gatk_hg38/hapmap_3.3.hg38.vcf
hg38_fasta=/pico/scratch/userexternal/mfratell/GENOME/iGENOMES/Homo_sapiens/hg38/Sequence/BWAIndex/genome.fa

cd alla vostra cartella con i dati di padre/madre o sofia. 

#####################################
######    # PRE-PROCESSING #    #####
#####################################

# Mark Duplicates
java -XX:ParallelGCThreads=20 -jar $PICARD_PATH/picard.jar MarkDuplicates INPUT=padre_bwa_sorted.bam OUTPUT=dedup_padre_bwa_sorted.bam METRICS_FILE=metrics.txt

# Re-index bam
java -XX:ParallelGCThreads=20 -jar $PICARD_PATH/picard.jar BuildBamIndex INPUT=dedup_padre_bwa_sorted.bam

# if missing genome.dict
# java -jar /pico/scratch/userexternal/mfratell/SW/picard.jar CreateSequenceDictionary R=genome.fa O=genome.dict

# Replace read groups
java -Xmx500g -XX:ParallelGCThreads=20 -jar /pico/scratch/userexternal/mfratell/SW/picard.jar AddOrReplaceReadGroups I=dedup_padre_bwa_sorted.bam O=dedup_padre_bwa_sorted_RG.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=1
samtools index dedup_padre_bwa_sorted_RG.bam

# Recalibrate Bases
java -Xmx500g -XX:ParallelGCThreads=20 -jar $GATK_PATH/GenomeAnalysisTK.jar -T BaseRecalibrator -R $hg38_fasta -I dedup_padre_bwa_sorted_RG.bam -knownSites $dbSNP -knownSites $indels -o recal_data.table 

# Recalibrate Bases 2nd-pass
 java -Xmx500g -XX:ParallelGCThreads=20  -jar $GATK_PATH/GenomeAnalysisTK.jar -T BaseRecalibrator -R $hg38_fasta -I dedup_padre_bwa_sorted_RG.bam -knownSites $dbSNP -knownSites $indels -BQSR recal_data.table -o post_recal_data.table 
 
 # Generate before/after plots
 # install R-packages ggplots, gplots, reshape, gsalib
 java -Xmx500g -XX:ParallelGCThreads=20  -jar $GATK_PATH/GenomeAnalysisTK.jar -T AnalyzeCovariates -R $hg38_fasta -before recal_data.table -after post_recal_data.table -plots recalibration_plots.pdf
    
 # Apply recalibration to data
 # removed options suchs as ParallelGCThreads and -Xmx to avoid SIGSEGV error
 java -jar $GATK_PATH/GenomeAnalysisTK.jar -T PrintReads -R $hg38_fasta -I dedup_padre_bwa_sorted_RG.bam -BQSR recal_data.table -o padre_recalled.bam 
  
  
#####################################
######  # VARIANT DISCOVERY #  ######
#####################################

java -Xmx500g -XX:ParallelGCThreads=20  -jar $GATK_PATH/GenomeAnalysisTK.jar -T HaplotypeCaller -R $hg38_fasta -I padre_recalled.bam  --genotyping_mode DISCOVERY  -stand_call_conf 30 -o raw_variants.vcf 

#Generating the bamout for a single site or interval
#java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R human_b37_20.fasta -I recalibrated.bam -o hc_variants.vcf -L 20:10255630-10255840 -bamout bamout.bam

#Generating the bamout for multiples intervals/whole genome
#java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R human_b37_20.fasta -I recalibrated.bam -o hc_variants.vcf -bamout bamout.bam

# JOINT GENOTYPE (MULTIPLE SAMPLES)

# Recalibrate variant quality scores - VSQR
java -Xmx500g -XX:ParallelGCThreads=20 -jar $GATK_PATH/GenomeAnalysisTK.jar -T VariantRecalibrator -R $hg38_fasta -input raw_variants.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap_SNP -resource:omni,known=false,training=true,truth=true,prior=12.0 $omni_SNP -resource:1000G,known=false,training=true,truth=false,prior=10.0 $G1000_SNP -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbSNP -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile recalibrate_SNP.recal -tranchesFile recalibrate_SNP.tranches -rscriptFile recalibrate_SNP_plots.R 
