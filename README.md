# Estimation of paternity share from PoolSeq data

This repository is part of Henry Arenas-Castro's PhD thesis dissertation at the University of Queensland and contains step-by-step instructions to estimate the fertilisation outcome of competitive crosses from PoolSeq data. Maddie E. James assisted with data analysis and Daniel Ortiz-Barrientos and Jan Engelstädter supervised this project.

## Fertilisation null model

In absence of fertilisation biases, it is expected that when two sires simultaneously mate or pollinate a dam, each sire fertilises half of the available ovules. Consequently, the estimation of paternity share in the offspring of competitive crosses may be used to test for deviations in the fertilisation outcome, proven that other technical and biological factors are controlled (see below).

![Alt text](Figures/Figure_FertilisationBias.png?raw=true "Title")

Traditional approaches have estimated the paternity share of the offspring by genotyping each individual. However, this is infeasible in large-scale studies. Here, we develop a method to estimate the paternity share in the offspring of competitive crosses using PoolSeq data. Instead of individually genotyping the offspring, we pooled equal amount of tissue from each individual, homogenise them, extracted DNA, prepared RAD-seq libraries, and sequenced a single sample per cross.

To infer the paternity share of the offspring, we (*i*) identified paternity markers based on the genetic variation of the populations of the parentals and (*ii*) estimated the allelic frequency of these markers in the offspring. 

In the case of competitive crosses between sires from different populations, the genetic markers that are most informative about the paternity of the offpsring are those that both are monomorphic in each population and have different alleles fixed between the two population. These markers should also be monomorphic in the dam population. Depending on whether the dam population shares the same allele with one of the sire populations or not, the expected allelic frequencies in the offspring vary.

![Alt text](Figures/Figure_DiagramPaternityMarkers.png?raw=true "Title")

## Paternity markers

All the parental individuals used in the competitive crosses (N=139) were individually genotyped following James et al. (2021) RAD-seq protocol. Then, we jointly called the SNPs for each set of individuals used either as dams or sires per population to identify the sites that were monomorphic. Finally, we identified those sites that had differentially fixed alleles between the sire populations.

### Reference genome indexing

We used ```bwa``` to index the *Senecio lautus* reference genome version *SPD_CN1K_CtgRN*.

```
bwa index Senecio.contigs.fasta
```

### Reads alignment to the reference genome

For each individual, we aligned their reads to the reference genome with the ```BWA-MEM``` algorithm. The ```*_1.fq.gz``` and ```*_2.fq.gz``` files correspond to the forward and reverse read files for each individual. We further used ```samtools``` to filter for mapping quality ```-q 20```.

```
bwa mem -M -t 12 reference_genome.fasta ind1_1.fq.gz ind1_2.fq.gz | samtools view -q 20 -b | samtools sort -@ 12 -T ind1 > ind1_sorted.bam
```

### Alignment quality stats

We also used ```samtoools``` to produce a summary of alignment stats per individual, including mean depth and the number and percentage of reads mapped to the reference genome.

```
samtools depth ind1_sorted.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print sum/NR,"\t", sqrt(sumsq/NR - (sum/NR)**2)}' > ind1_bamdepth.txt
```

```
samtools flagstat -O tsv ind1_sorted.bam > ind1_bamstats.txt
```

![Alt text](Figures/Figure_SummaryStats.png?raw=true "Title")

### Alignment cleaning

To clean the bam files and add read groups to each individual, we used ```PicardTools```, where ```RGID``` is the read group ID (i.e. the name of the individual), ```RGLB``` is the read group library (i.e. the sequencing lane number), ```RGPU``` is the read group platform unit (i.e. run barcode of the lane), and ```RGSM``` is the read group sample name (i.e. the name of the individual).

```
java -jar picard.jar CleanSam INPUT=ind1_sorted.bam OUTPUT=ind1_sorted_cleaned.bam
```

```
java -jar picard.jar AddOrReplaceReadGroups INPUT=ind1_sorted_cleaned.bam OUTPUT=ind1_sorted_cleaned_rg.bam SORT_ORDER=coordinate RGID=ind1 RGLB=lib1 RGPL=ILLUMINA RGPU=lane1_parentals RGSM=ind1 CREATE_INDEX=True
```

### SNP calling

We used ```FreeBayes``` to jointly call SNPs for each set of dams and sires per population. Joint calling assumes all samples are genetically similar and uses information across all samples to call a SNP. SNPs from these cohorts are then combined into one VCF file. The code below shows an example where ```ind1```, ```ind2```, and ```ind3``` conform the set of dams from ```pop1```.

```
freebayes -f Senecio.contigs.fasta ind1_sorted_cleaned_rg.bam ind2_sorted_cleaned_rg.bam ind3_sorted_cleaned_rg.bam --use-best-n-alleles 4 --report-monomorphic --genotype-qualities > pop1_dams.vcf
```

### SNP filtering

For each cross, or pool, we used ```bcftools view``` to subsample the individuals that effectively served as sires and dams. The ```-s``` flag specifies the name of the individuals to sample from the population's VCF (i.e. only ind1 and ind3 from the step above).

```
bcftools view -s ind1,ind3 pop1_dams.vcf > pop1_dams_ind.vcf
```

Before starting to filter the SNPs, we used ```bgzip``` and ```tabix``` to compress and index each VCF.

```
bgzip -c pop1_dams_ind.vcf > pop1_dams_ind.vcf.gz
tabix -p vcf pop1_dams_ind.vcf.gz
```

To identify and retain the sites that vary among the three parental populations, we used ```bcftools merge``` and ```vcftools``` to merge the VCFs and keep the sites genotyped in at least 50% of individuals with a minimum quality score of 30.

```
bcftools merge pop1_dams_ind.vcf.gz pop2_sires_ind.vcf.gz pop3_sires_ind.vcf.gz -o Pool1.vcf
vcftools --gzvcf Pool1.vcf --max-missing 0.5 --mac 1 --minQ 30 --recode --recode-INFO-all --stdout | gzip -c > Pool1_Q30mac1.vcf.gz
```

We applied a minimum depth per sample of 3 reads (genotypes with less than 3 reads were recoded as missing data).

```
vcftools --gzvcf Pool1_Q30mac1.vcf.gz --minDP 3 --recode --recode-INFO-all --stdout > Pool1_Q30mac1dp3.vcf
```

Sites with high coverage could be multi-copy regions of the genome, so we removed these potential paralogues. We first calculated the mean depth per site after pooling the VCFs of all the individuals used in this study and plot their distribution in ```R```. Then, we visually examined the distribution of mean read depth across all sites and selected a threshold value of 250.

```
vcftools --vcf Pool1_Q30mac1dp3.vcf --site-mean-depth --out mean_depth
```

![Alt text](Figures/Figure_MeanReadDepth.png?raw=true "Title")

```
vcftools --vcf  Pool1_Q30mac1dp3.vcf --recode --recode-INFO-all --out Pool1_Q30mac1dp3MaxDP250 --max-meanDP 250
```

We then filtered for a minimum mean depth of 10.

```
vcftools --vcf Pool1_Q30mac1dp3MaxDP250.recode.vcf --min-meanDP 10 --recode --recode-INFO-all --out Pool1_Q30mac1dp3MaxDP250MinDP10
```

We additionally filtered for an overall missing data, removing sites if they had greater than 20% missing data.

```
vcftools --vcf Pool1_Q30mac1dp3MaxDP250MinDP10.recode.vcf --recode --recode-INFO-all --max-missing 0.8 --out Pool1_Q30mac1dp3MaxDP250MinDP10md80
```







## References

James ME, Arenas-Castro H, Groh JS, Engelstädter J, Allen SL, Ortiz-Barrientos D. 2021. Highly replicated evolution of parapatric ecotypes. Molecular Biology and Evolution, In press.
