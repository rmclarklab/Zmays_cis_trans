# Quantificaiton of allele-specific expression in maize (Zea mays)

## Experimental Design
The B96 and B75 maize inbred lines are highly resistant to multiple herbivores, including the devastating European corn borer (<i> Ostrinia nubilalis </i>) and the generalist two-spotted spider mites (<i>Tetranychus urticae</i>) [(Bui, et al., 2021)](https://doi.org/10.3389/fpls.2021.693088). To further understand the genetic basis underlying their resistant phenotype, leaf tissue from maize plants under uninfested (control) or infested (24 h, <i> T. urticae </i> treatment) conditions was collected from B73, B96, and B75, and the respective F1 plants from crosses of B96 and B75 to B73, where RNA-seq was extracted from. See the publication for this study [ref]().
- 2 independent experiement for B96 vs. B73 and B75 vs. B73 (B73 is relative susceptible line)
- 2 conditions: the steady (control, **C**) and infested (mite herbivory treatmetn, **T**) conditions
- 4 biological replicates 

## Programs version used in this study
**bwa** Version: 0.7.17-r1188 <br>
**GATK** Version: 4.1.4.1 <br>
**samtools** Version: 1.9 (using htslib 1.9) <br>
**python** Version: 3.7.3 (including python pacakages: pysam, pandas, etc. compatible with the python version) <br>

## SNPs calling to enable allele-specific expression analysis
#### We followed best practice recommendations for the [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-) pipeline for variants calling purposes. <br>
1. Generate index for reference genome <br>
- create reference index for bwa mapping <br>
  <code> bwa index -p <prefix> reference.fasta </code>
- build reference directionry using picard <br>
  <code> picard CreateSequenceDictionary --REFERENCE reference.fasta </code>
2. DNA-seq mapping using bwa <br>
  <code> bwa mem -t 20 reference.fasta read_1.fq read_2.fq | samtools view -Su - | samtools sort -@ 20 - -o output.bam </code>
3. mark duplicates for bam file using picard <br>
  <code> picard MarkDuplicates -I input.bam -O output_mark.bam -M output.metrics -ASO coordinate --CREATE_INDEX true --REMOVE_DUPLICATES && samtools index -@ 20 output_mark.bam </code>
4. align indexs on the left side <br>
  <code> gatk LeftAlignIndels -R reference.fasta -I output_mark.bam -O output_mark_leftalign.bam --create-output-bam-index true </code>
5. call SNPs and indels via local de-novo assembly of haplotypes <br>
  <code> gatk HaplotypeCaller -R reference.fasta -I output_mark_leftalign.bam -ERC GVCF -L chromosome_1 -O chromosome_1.g.vcf.gz --sequence-dictionary reference.dict </code>
** Note: haplotypecaller needs to be called for each individual chromosome. <br>
6. Combine all separate g.vcf.gz for individual chromosome (/scaffold) into one g.vcf.gz <br>
  <code> gatk combineGVCFs -R reference.fasta --variant chromosome_1.g.vcf.gz --variant chromosome_2.g.vcf.gz ...... -O all.g.vcf.gz </code>
7. Call genotype file in vcf format <br>
  <code> gatk GenotypeGVCFs -R reference.fasta -V all.g.vcf.gz -O all.vcf.gz </code> 
8. To split SNPs and INDELs in the vcf file <br>
  <code> gatk SelectVariants -R reference.fasta -V all.vcf.gz -select-type-to-include SNP -O all.SNP.vcf.gz </code>
  <code> gatk SelectVariants -R reference.fasta -V all.vcf.gz -select-type-to-include INDEL -O all.INDEL.vcf.gz </code>
9. We performed hard-filtering to ratain high-confident SNPs. <br> 
  Custom python script vcf_pass.py
 for this purposes (criteria: MQ > 40 and homozygous genotype): <br>
  <code> python vcf_pass.py -vcf all.SNP.vcf.gz -R reference.fasta -O SNP_fitered </code>

## RNA-seq mapping against the its (corrected) reference genome and expression estimation
#### We aligned RNA-reads from B73, B96, B75, and the respective F1s, to the one or more of the B73 or SNP-corrected genome sequences. [B73v4](https://www.maizegdb.org/genome/assembly/Zm-B73-REFERENCE-GRAMENE-4.0) was used as a master reference sequence for the SNP-corrected genome reconstruction. <br>
1. Build reference index for star mapping <br>
  <code> STAR --runMode genomeGenerate --genomeDir STAR_index --runThreadN 20 --genomeFastaFiles reference.fasta </code>
2. RNA-seq mapping against the reference genome using STAR <br>
  <code> STAR --genomeDir STAR_index --runThreadN 20 --readFilesIn read_1.fastq.gz read_2.fastq.gz --twopassMode Basic --sjdbOverhang 149 --outFileNamePrefix STAR_map --readFilesCommand zcat --alignIntronMax 60000 --outSAMtype SortedByCoordinate && samtools index STAR_map.bam -@ 40 </code>
3. htseq-count was used to count the reads on gene basis <br>
  <code> htseq-count -f bam -r pos -s reverse -t exon --nonunique none STAR_map.bam reference.gtf > sample_count.txt </code>
4. Scripts developed for the reconstruction of SNP-corrected genomes see restore_genome.py. <br>
  usage: <br>
  <code> python restore_genome.py -R reference.fasta -vcf_snp <snp.vcf> -O ref_SNP_corrected.fa </code>
5. Transcript per million (TPM) was calculated using [RSEM](https://deweylab.github.io/RSEM/README.html) <br>
    - prepare reference for star-mapping based TPM calculation <br>
    <code> rsem-prepare-reference -gtf gtf --star reference.fasta reference_name </code>
    - calculate gene and isoform level expression <br>
    <code> rsem-calculate-expression --paired-end --star-gzipped-read-file --star --strandedness reverse -p 45 read_1.fastq.gz read_2.fastq.gz reference_name output </code>

## Developed pipeline for allele-specific expression (ASE) analyses (step by step)
Preparing files:
- GTF annotation of reference genome (sorted and indexed);
- VCF of SNPs from DNA-seq alignment (after hard-filtering, sorted and indexed);
- BAM of RNA-seq alignment for parental lines and F1 hybrid (sorted and indexed);
1. Since the VCF from DNA-seq alignment includes SNPs not on gene coding region, which is not informative for ASE analysis, we rewrite a new VCF file and only include SNPs on gene coding region based on genome annotation file. <br>
<code> python vcf_coding.py -vcf VCF -gtf GTF -O output </code>
Note: To sort and index VCF file, you need install [bcftools](https://samtools.github.io/bcftools/) <br>
2. Before ASE, we further validated SNPs from RNA-seq alignment of parental lines. <br>
<code> python snp_allele_count.py -vcf VCF -bam parental_BAM -gtf GTF -O parental_count </code> <br>
Output table includes information of allelic read count at SNP sites based on RNA BAM file. 
3. Two parental lines provided RNA-seq for SNP validation, replicates from single genotype can be combined (no matter what condition) for given genotype SNP validation. <br> 
<code> python SNP_filter.py -dir DIR -O output </code> <br>
NOTE: Given DIR should include allelic read count on SNP-basis of RNA-seq alignment from single genotype. If the count for the genotype is reference, indicate "-is_ref" in the command line, otherwise don't use that label. 
From this step, we can sort out SNPs for which RNA-seq inconsistent with DNA-seq. And we will drop these SNPs from further ASE analysis. 
4. ASE on gene-basis for given BAM file. <br>
<code> python ASE.py -v VCF -gtf GTF -bam BAM -O output </code> <br>
NOTE: VCF file should only include SNPs pass filtering of DNA and RNA-based seq alignment. 

## ASE for cis and trans classification <br>
The methodology of cis and trans classification was adapted from previous [publications](https://www.nature.com/articles/s41467-019-13386-w). <br>
    
