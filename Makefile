export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files

LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'

SAMPLES=14_0152 14_0174 14_0291 StA_BS StA_ES StA_SS
#SAMPLES=test

all: fastqc bwa picard gatk varscan snpeff filtered-variants

include ~/generic/scripts/rna-seq/fastqc.mk

clean:
	rm -rf bwa picard gatk varscan snpeff filtered-variants tmp
	
#----------
# DOWNLOAD SAMPLES
#----------
download:
#	wget --user=gomftp --password='#gomf%74' --limit-rate=700k ftp://ftp.uniklinik-duesseldorf.de/medgen/14_0152_R1.fastq.gz
#	wget --user=gomftp --password='#gomf%74' --limit-rate=700k ftp://ftp.uniklinik-duesseldorf.de/medgen/14_0152_R2.fastq.gz
#	wget --user=gomftp --password='#gomf%74' --limit-rate=700k ftp://ftp.uniklinik-duesseldorf.de/medgen/14_0174_R1.fastq.gz
#	wget --user=gomftp --password='#gomf%74' --limit-rate=700k ftp://ftp.uniklinik-duesseldorf.de/medgen/14_0174_R2.fastq.gz
#	wget --user=gomftp --password='#gomf%74' --limit-rate=700k ftp://ftp.uniklinik-duesseldorf.de/medgen/14_0291_R1.fastq.gz
#	wget --user=gomftp --password='#gomf%74' --limit-rate=700k ftp://ftp.uniklinik-duesseldorf.de/medgen/14_0291_R2.fastq.gz

	wget --user=gomftp --password='#gomf%74' --limit-rate=700k ftp://ftp.uniklinik-duesseldorf.de/medgen/StA_BS_R1.fastq.gz
	wget --user=gomftp --password='#gomf%74' --limit-rate=700k ftp://ftp.uniklinik-duesseldorf.de/medgen/StA_BS_R2.fastq.gz
	wget --user=gomftp --password='#gomf%74' --limit-rate=700k ftp://ftp.uniklinik-duesseldorf.de/medgen/StA_ES_R1.fastq.gz
	wget --user=gomftp --password='#gomf%74' --limit-rate=700k ftp://ftp.uniklinik-duesseldorf.de/medgen/StA_ES_R2.fastq.gz
	wget --user=gomftp --password='#gomf%74' --limit-rate=700k ftp://ftp.uniklinik-duesseldorf.de/medgen/StA_SS_R1.fastq.gz
	wget --user=gomftp --password='#gomf%74' --limit-rate=700k ftp://ftp.uniklinik-duesseldorf.de/medgen/StA_SS_R2.fastq.gz

#-----------	
# ALIGNMENT, SORTING, MARK DUPLICATES, INDEXING
#-----------

.PHONY: bwa
bwa: $(foreach S, $(SAMPLES), bwa/$S.bwa.sorted.dupmarked.bam.bai)
	
bwa/%.bwa.bam: ~/generic/data/broad/human_g1k_v37.fasta ~/oskar/data/fastq/%_R1.fastq.gz ~/oskar/data/fastq/%_R2.fastq.gz
	mkdir -p bwa
	flock -x .lock ~/tools/bwa-0.7.10/bwa mem \
		-t 20 \
		-R '@RG\tID:RG1\tPL:Illumina\tSM:$*' \
		$(word 1,$^) $(word 2,$^) $(word 3,$^) | \
			~/tools/samtools-0.1.19/samtools view -bhS - \
		2>&1 1>$@.part | $(LOG)
	mv $@.part $@

bwa/%.bwa.sorted.bam: bwa/%.bwa.bam
	~/tools/samtools-0.1.19/samtools sort -@ 5 -o $< $* 2>&1 1>$@.part | $(LOG)
	mv $@.part $@
	rm $<

bwa/%.bwa.sorted.dupmarked.bam: bwa/%.bwa.sorted.bam
	mkdir -p picard
	java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar ~/tools/picard-tools-1.114/MarkDuplicates.jar \
		INPUT=$< \
		OUTPUT=$@.part \
		METRICS_FILE=picard/$*.mark_duplicates_metrics \
		VALIDATION_STRINGENCY=LENIENT \
		2>&1 1>$@.part | $(LOG)
	mv $@.part $@
	rm $<

bwa/%.bwa.sorted.dupmarked.bam.bai: bwa/%.bwa.sorted.dupmarked.bam
	rm -f $@
	~/tools/samtools-0.1.19/samtools index $^ $@.part 2>&1 | $(LOG)
	mv $@.part $@

#-----------
# GATK INDEL REALIGNMENT
#-----------

.PHONY: gatk
gatk: $(foreach S, $(SAMPLES), gatk/$S.bwa.sorted.dupmarked.realigned.bam)

gatk/%.intervalList.intervals: bwa/%.bwa.sorted.dupmarked.bam bwa/%.bwa.sorted.dupmarked.bam.bai
	mkdir -p gatk
	java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -Xmx2g -jar ~/tools/GenomeAnalysisTK-3.3.0/GenomeAnalysisTK.jar \
		-T RealignerTargetCreator \
		-R ~/generic/data/broad/human_g1k_v37.fasta \
		-I $< \
		-o $@.part \
		2>&1 | $(LOG)
	mv $@.part $@
   
gatk/%.bwa.sorted.dupmarked.realigned.bam: bwa/%.bwa.sorted.dupmarked.bam gatk/%.intervalList.intervals
	java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -Xmx4g -jar ~/tools/GenomeAnalysisTK-3.3.0/GenomeAnalysisTK.jar \
		-T IndelRealigner \
		-R ~/generic/data/broad/human_g1k_v37.fasta \
		-I $(word 1,$^) \
		-targetIntervals $(word 2,$^) \
		-o $@.part \
		2>&1 | $(LOG)
	mv $@.part $@
	rm -f gatk/%.bwa.sorted.dupmarked.realigned.bam.part.bai
	
gatk/%.bwa.sorted.dupmarked.realigned.bam.bai: gatk/%.bwa.sorted.dupmarked.realigned.bam
	rm -f $@
	~/tools/samtools-0.1.19/samtools index $^ $@.part 2>&1 | $(LOG)
	mv $@.part $@

#-----------	
# PICARD METRICS
#-----------	

.PHONY: picard
picard: $(foreach S, $(SAMPLES), picard/$S.multiplemetrics picard/$S.hs_metrics)

picard/%.multiplemetrics: bwa/%.bwa.sorted.dupmarked.bam bwa/%.bwa.sorted.dupmarked.bam.bai
	mkdir -p picard
	java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar ~/tools/picard-tools-1.114/CollectMultipleMetrics.jar \
		INPUT=$< \
		OUTPUT=picard/$* \
		VALIDATION_STRINGENCY=LENIENT \
		PROGRAM=CollectAlignmentSummaryMetrics \
		PROGRAM=CollectInsertSizeMetrics \
		PROGRAM=QualityScoreDistribution \
		PROGRAM=MeanQualityByCycle \
		2>&1 | $(LOG)
	touch $@

picard/%.hs_metrics: bwa/%.bwa.sorted.dupmarked.bam bwa/%.bwa.sorted.dupmarked.bam.bai
	mkdir -p picard
	~/tools/samtools-0.1.19/samtools view -H $< 2>&1 1> picard/$*.truseq-for-picard.bed | $(LOG)
	gawk 'BEGIN { OFS="\t"} {print $$1,$$2,$$3,$$6,$$4 }' ~/generic/data/illumina/truseq_exome_targeted_regions.hg19.bed >> picard/$*.truseq-for-picard.bed
	java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar ~/tools/picard-tools-1.114/CalculateHsMetrics.jar \
		BAIT_INTERVALS=picard/$*.truseq-for-picard.bed \
		TARGET_INTERVALS=picard/$*.truseq-for-picard.bed \
		INPUT=$< \
		OUTPUT=$@.part \
		REFERENCE_SEQUENCE=~/generic/data/broad/human_g1k_v37.fasta \
		PER_TARGET_COVERAGE=picard/$*.hs_metrics.per_target_coverage.part \
		VALIDATION_STRINGENCY=LENIENT \
		2>&1 | $(LOG)
	mv picard/$*.hs_metrics.per_target_coverage.part picard/$*.hs_metrics.per_target_coverage
	rm picard/$*.truseq-for-picard.bed
	mv $@.part $@

#-----------	
# VARSCAN
#-----------	

.PHONY: varscan
varscan: $(foreach S, $(SAMPLES), varscan/$S.varscan.vcf)

varscan/%.varscan.vcf: gatk/%.bwa.sorted.dupmarked.realigned.bam gatk/%.bwa.sorted.dupmarked.realigned.bam.bai
	mkdir -p varscan
	# -q 10 : remove non-uniquely mapping reads and reads with low mapping quality
	# -F 1024 : remove PCR and optical duplicates
	# -f 2 : only reads mapped in proper pairs
	java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -jar ~/tools/varscan-2.3.6/VarScan.v2.3.6.jar mpileup2cns \
		<(~/tools/samtools-0.1.19/samtools view -u -q 10 -F 1024 -f 2 $(word 1,$^) | ~/tools/samtools-0.1.19/samtools mpileup -f ~/generic/data/broad/human_g1k_v37.fasta -) \
		--min-coverage 8 \
		--min-reads2 4 \
		--min-avg-qual 20 \
		--min-var-freq 0.2 \
		--min-freq-for-hom 0.75 \
		--p-value 1 \
		--strand-filter 1 \
		--output-vcf 1 \
		--variants 1 \
		2>&1 1>$@.part | $(LOG)
	mv $@.part $@

#-----------	
# SNPEFF
#-----------

.PHONY: snpeff
snpeff: $(foreach S, $(SAMPLES), snpeff/$S.varscan.dbsnp.snpeff.dbNSFP.vcf.bgz.tbi)

snpeff/%.varscan.dbsnp.vcf: varscan/%.varscan.vcf ~/generic/data/ncbi/common_no_known_medical_impact_20140826.vcf
	mkdir -p snpeff
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.6; java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -jar SnpSift.jar annotate \
		-v ~/generic/data/ncbi/common_no_known_medical_impact_20140826.vcf \
		<(cat $(PWD)/$< | perl -ne 's/\trs\d+\t/\t.\t/; print $$_;' -) \
		2>&1 1>$(PWD)/$@.part) | $(LOG)
	test -s $@.part
	mv $@.part $@

snpeff/%.varscan.dbsnp.snpeff.vcf: snpeff/%.varscan.dbsnp.vcf
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.6; java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -Xmx4g -jar snpEff.jar -v -lof GRCh37.75 -stats $(PWD)/snpeff/$*.snpeff.summary.html $(PWD)/$< 2>&1 1>$(PWD)/$@.part) | $(LOG)
	mv $@.part $@

snpeff/%.varscan.dbsnp.snpeff.dbNSFP.vcf: snpeff/%.varscan.dbsnp.snpeff.vcf ~/generic/data/dbNSFP-2.6/dbNSFP2.6_variant.tsv.gz
	PWD=$(pwd)
	(cd ~/tools/snpEff-3.6; java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -jar SnpSift.jar dbnsfp \
		-v ~/generic/data/dbNSFP-2.6/dbNSFP2.6_variant.tsv.gz \
		-collapse \
		-f SIFT_pred,SIFT_score,Polyphen2_HVAR_pred,Polyphen2_HVAR_score,SiPhy_29way_logOdds,LRT_pred,LRT_score,MutationTaster_pred,MutationTaster_score,MutationAssessor_pred,MutationAssessor_score,FATHMM_pred,FATHMM_score,RadialSVM_pred,RadialSVM_score,GERP++_RS,1000Gp1_AF,1000Gp1_AFR_AF,1000Gp1_EUR_AF,1000Gp1_AMR_AF,1000Gp1_ASN_AF,ESP6500_AA_AF,ESP6500_EA_AF,Uniprot_acc,Interpro_domain, \
		$(PWD)/$< 2>&1 1>$(PWD)/$@.part) | $(LOG)
	mv $@.part $@

snpeff/%.vcf.bgz.tbi: snpeff/%.vcf
	bgzip -c $< > snpeff/$*.vcf.bgz.part
	mv snpeff/$*.vcf.bgz.part snpeff/$*.vcf.bgz
	~/tools/tabix-0.2.6/tabix -p vcf snpeff/$*.vcf.bgz

#-----------	
# FILTERED VARIANTS LIST
#-----------	

.PHONY: filtered-variants
filtered-variants: $(foreach S, $(SAMPLES), filtered-variants/$S.tsv)

filtered-variants/%.tsv: snpeff/%.varscan.dbsnp.snpeff.dbNSFP.vcf ~/oskar/scripts/filter-variants.pl
	mkdir -p filtered-variants
	perl ~/oskar/scripts/filter-variants.pl \
		$< \
		--patient $* \
		--rmsk-file ~/generic/data/hg19/hg19.rmsk.txt.gz \
		--simpleRepeat-file ~/generic/data/hg19/hg19.simpleRepeat.txt.gz \
		--segdup-file ~/generic/data/hg19/hg19.genomicSuperDups.txt.gz \
		--blacklist-file ~/generic/data/hg19/hg19.wgEncodeDacMapabilityConsensusExcludable.txt.gz \
		--g1k-accessible ~/generic/data/hg19/paired.end.mapping.1000G..pilot.bed.gz \
		--ucscRetro ~/generic/data/hg19/hg19.ucscRetroAli5.txt.gz \
		--remission-variants-file ~/hdall/results/remission-variants.tsv.gz \
		--cosmic-mutation-file ~/generic/data/cosmic/v67/CosmicMutantExport_v67_241013.tsv \
		--evs-file ~/generic/data/evs/ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt.gz \
		--exac-file ~/generic/data/ExAC/ExAC.r0.2.sites.vep.vcf.gz \
		--clinvar-file ~/generic/data/clinvar/clinvar_20140929.vcf.gz \
		2>&1 1> $@.part | $(LOG)
	mv $@.part $@
	
filtered-variants/StA.tsv: snpeff/StA_BS.varscan.dbsnp.snpeff.dbNSFP.vcf snpeff/StA_ES.varscan.dbsnp.snpeff.dbNSFP.vcf.bgz.tbi snpeff/StA_SS.varscan.dbsnp.snpeff.dbNSFP.vcf.bgz.tbi ~/oskar/scripts/filter-variants.pl
	mkdir -p filtered-variants
	perl ~/oskar/scripts/filter-variants.pl \
		$< \
		--patient StA_BS \
		--mother snpeff/StA_ES.varscan.dbsnp.snpeff.dbNSFP.vcf.bgz \
		--father snpeff/StA_SS.varscan.dbsnp.snpeff.dbNSFP.vcf.bgz \
		--rmsk-file ~/generic/data/hg19/hg19.rmsk.txt.gz \
		--simpleRepeat-file ~/generic/data/hg19/hg19.simpleRepeat.txt.gz \
		--segdup-file ~/generic/data/hg19/hg19.genomicSuperDups.txt.gz \
		--blacklist-file ~/generic/data/hg19/hg19.wgEncodeDacMapabilityConsensusExcludable.txt.gz \
		--g1k-accessible ~/generic/data/hg19/paired.end.mapping.1000G..pilot.bed.gz \
		--ucscRetro ~/generic/data/hg19/hg19.ucscRetroAli5.txt.gz \
		--remission-variants-file ~/hdall/results/remission-variants.tsv.gz \
		--cosmic-mutation-file ~/generic/data/cosmic/v67/CosmicMutantExport_v67_241013.tsv \
		--evs-file ~/generic/data/evs/ESP6500SI-V2-SSA137.updatedRsIds.chrAll.snps_indels.txt.gz \
		--exac-file ~/generic/data/ExAC/ExAC.r0.2.sites.vep.vcf.gz \
		--clinvar-file ~/generic/data/clinvar/clinvar_20140929.vcf.gz \
		2>&1 1> $@.part | $(LOG)
	mv $@.part $@
	