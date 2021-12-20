# Zmays_cis_trans

 # Concerted cis and trans effects underpin heightened defense gene expression in multi-herbivore resistant maize lines
#### Data storage location: 

## Experimental Design

#### The B96 and B75 maize inbred lines are highly resistant to multiple herbivores, including the devastating European corn borer (<i> Ostrinia nubilalis </i>) and the generalist two-spotted spider mites (<i> Tetranychus urticae </i>) [ref]("https://www.frontiersin.org/articles/10.3389/fpls.2021.693088/full"). To further understand the genetic basis underlying their resistant phenotype, leaf tissue from maize plants under uninfested (control) or infested (24 h, <i> T. urticae </i> treatment) conditions was collected from B73, B96, and B75, and the respective F1 plants from crosses of B96 and B75 to B73, where RNA-seq was extracted from.
- 2 independent experiement for B96 vs. B73 and B75 vs. B73 (B73 is relative susceptible line)
- 2 conditions: the steady (control, **C**) and infested (treatment, **T**) conditions
- 4 biological replicates 

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
  Custom python script [vcf_pass.py](https://github.com/richardmichaelclark/Maize_cistrans/blob/main/vcf_pass.py)
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
4. Scripts developed for the reconstruction of SNP-corrected genomes see [restore_genome.py](https://github.com/richardmichaelclark/Maize_cistrans/blob/main/restore_genome.py). <br>
  usage: <br>
  <code> python restore_genome.py -R reference.fasta -vcf_snp <snp.vcf> -O ref_SNP_corrected.fa </code>
5. Transcript per million (TPM) was calculated using [RSEM](https://deweylab.github.io/RSEM/README.html) <br>
    - prepare reference for star-mapping based TPM calculation <br>
    <code> rsem-prepare-reference -gtf gtf --star reference.fasta reference_name </code>
    - calculate gene and isoform level expression <br>
    <code> rsem-calculate-expression --paired-end --star-gzipped-read-file --star --strandedness reverse -p 45 read_1.fastq.gz read_2.fastq.gz reference_name output </code>
    
## Differential gene expression analysis by performing pairwise comparisons
#### We employing [DESeq2](https://github.com/mikelove/DESeq2) for the differential gene expression (DGE) analysis by comparing between genotypes by treatment, or between control and treatment by a given genotype. 
#### Inputs prepared for DESeq2:
  - htseq-count for all samples in one folder;
  - sample information table matched to htseq-count files;
#### Please refer to the script [diff_expr_htseq.R](https://github.com/richardmichaelclark/Maize_cistrans/blob/de501d2ef5689686828c9c78a6b6c41a437c963c/diff_expr_htseq.R) for DGE analysis.
#### Please refer to the script [PCA_DESeq2.R](https://github.com/richardmichaelclark/Maize_cistrans/blob/main/PCA_DESeq2.R) for principal component analysis using DESeq2. 

## Developed pipeline for allele-specific expression (ASE) analyses 
#### By taking SNP variants after DNA-based filtering, we performed gene level ASE analyses. 
#### Inputs for ASE analysis:
  - allele-specific read count at allele site (for RNA-based SNP filtering)
  - RNA-mapping file for parental lines and its F1 hybrid (B96 x B73, B75 x B73);
  - Filtered VCF file with SNP variant on gene coding region;
  - Genes with its variants on coding region;
  - referene genome GTF file;
  - RNA-seq mapped BAM file;
#### Please refer to the script [ASE.py](https://github.com/richardmichaelclark/Maize_cistrans/blob/main/ASE.py) for measurement of gene-level allelic expression.
##### usage: 
  <code> python ASE.py -dir count_dir -v vcf -t coding_table -gtf gtf -bam bam -O ASE_result </code>

## Cis-/Trans- control underlying gene expression variation 
#### We used the ASE data to estimate the relative contribution of cis and trans regulatory effects in variation of gene expression. <br>
#### The methodology is adapted from [here](https://github.com/Wendellab/CisTransRegulation#analysis-of-cis-trans-regulation-underlying-cotton-fiber-domestication). Briefly, (1) the combination of cis and trans effects was measured by taking the log2 ratios of expression difference between parental lines (**A**=log2(P2/P1)), (2) cis effects were measured by log2 ratios of the corresponding allelic expression in F1 hybrids (**B**=log2(F1P2/F1P1)), and (3) trans effects were derived by subtracting the F1 allelic log2 ratios from the parental log2 ratios (**A**-**B**). The log2 ratios for comparison were calculated in DESeq2.
#### For the cis/trans classification detail:

#### Inputs: 
- DESeq2 output for differential expression of parental lines (allelic data utilized);
- DESeq2 output for allelic differential expression within F1 
#### Script [cistrans.R](https://github.com/richardmichaelclark/Maize_cistrans/blob/main/cistrans.R) usage: 
<code> Rscript cistrans.R -par par_DESeq.txt -F1 F1_DESeq.txt -O gene_classification </code>

## Classification of inheritance mode in gene expression 
#### The categories for gene expression inheritance were adapted from previous studies (e.g., in [Saccharomyces](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5604594/), [cotton](https://www.nature.com/articles/s41467-019-13386-w)). Our methodology was based on [here](https://github.com/Wendellab/CisTransRegulation#analysis-of-cis-trans-regulation-underlying-cotton-fiber-domestication). [publication](https://github.com/Wendellab/CisTransRegulation#analysis-of-cis-trans-regulation-underlying-cotton-fiber-domestication). On its basis, we further explored inheritance mode into partial dominance. So the final classification includes two files: one with partial dominance; one without partial dominance.
#### Inputs prepared for this analysis:
- read counts include parents and its F1 offspring
- sample table matched to the read count data
#### A universal script was developed in the present study, people who are not good at coding can easily perform their analysis by taking our [code](https://github.com/richardmichaelclark/Maize_cistrans/blob/main/inheritance_mode.R). Usage:
<code> Rscript inheritance_mode.R -count cnt_table.txt -sample sample_table.txt -par Parent1 Parent2 -F1 F1_label -O inheritance </code>
    
## Scripts for results vasulisation:
#### We generated all heatmap figure using [ComplexHeatmap](https://jokergoo.github.io/ComplexHeatmap-reference/book/) in R. Other types of plots were generated using [ggplot2](http://r-statistics.co/Complete-Ggplot2-Tutorial-Part1-With-R-Code.html) in R. 
* [lfc_heatmap.R](https://github.com/richardmichaelclark/Maize_cistrans/blob/main/lfc_heatmap.R): generate heatmap with log2foldchange values labeled based on its p-value.
* [venn.R](https://github.com/richardmichaelclark/Maize_cistrans/blob/main/venn.R): generate Venn Diagram for different gene set
* [regprop_lfc.R](https://github.com/richardmichaelclark/Maize_cistrans/blob/main/regprop_lfc.R): generate bar plot for the composition of each genetic control along with log2 foldchange between parental lines
* [GO_bar.R](https://github.com/richardmichaelclark/Maize_cistrans/blob/main/GO_bar.R): bar plot for top enriched GO terms with its description
* [mode_correlation.Rmd](https://github.com/richardmichaelclark/Maize_cistrans/blob/main/mode_correlation.Rmd): the association between genetic control and gene expression inheritance mode (Chi-squared test of independence)
* [sankey.Rmd](https://github.com/richardmichaelclark/Maize_cistrans/blob/main/sankey.Rmd): generate sankey plot for tracking the gene regulatory mode under control and treatment conditions
* [TPM_gene_horizontal.R](https://github.com/richardmichaelclark/Maize_cistrans/blob/main/TPM_gene_horizontal.R): generate bar plot with TPM value for individual genes in each genotype and condition

## Other scripts:
* [GO_enrich.R](https://github.com/richardmichaelclark/Maize_cistrans/blob/main/GO_enrich.R): perform Gene Ontology (GO) enrichment analysis by employing [clusterProfiler](https://guangchuangyu.github.io/software/clusterProfiler/)
