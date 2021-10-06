# Estimation of paternity share from PoolSeq data

This repository is part of Henry Arenas-Castro's PhD thesis dissertation at the University of Queensland and contains step-by-step instructions to estimate the fertilisation outcome of competitive crosses from PoolSeq data. Maddie E. James assisted with data analysis and Daniel Ortiz-Barrientos and Jan Engelstädter supervised this project.

## Fertilisation null model

In absence of fertilisation biases, it is expected that when two sires simultaneously mate or pollinate a dam, each sire fertilises half of the available ovules. Consequently, the estimation of paternity share in the offspring may be used to test for deviations in the fertilisation outcome, proven that other technical and biological factors are controlled (see below).

![Alt text](Figures/Figure_FertilisationBias.png?raw=true "Title")

Traditional approaches have estimated the paternity share of the offspring by genotyping each individual. However, this is infeasible in large-scale studies. Here, we develop a method to estimate the paternity share in the offspring of competitive crosses using PoolSeq data. Instead of individually genotyping the offspring, we pooled equal amount of tissue from each individual, homogenise them, extracted DNA, prepared RAD-seq libraries, and sequenced a single sample per cross.

To infer the paternity share of the offspring, we (*i*) identified paternity markers based on the genetic variation of the populations of the parentals and (*ii*) estimated the allelic frequency of these markers in the offspring. 

In the case of competitive crosses between sires from different populations, the genetic markers that are most informative about the paternity of the offpsring are those that both are monomorphic within each sire population and have different alleles fixed in each population. These markers should also be monomorphic in the dam population. Depending on whether the dam population shares the same allele with one of the sire populations or not, the expected allelic frequencies in the offspring vary.

![Alt text](Figures/Figure_DiagramPaternityMarkers.png?raw=true "Title")

## Paternity markers

All the parental individuals used in the competitive crosses (N=139) were individually genotyped following James et al. (2021) RAD-seq protocol. Then, we jointly called the SNPs for each set of individuals used either as dams or sires per population to identify the sites that were monomorphic. Finally, we identify those sites that had alleles differentially fixed between sire populations.

### Reference genome indexing


```
bwa index ${TMPDIR}/Senecio.contigs.fasta
```

Senecio reference genome version: SPD_CN1K_CtgRN


### Reads alignment to the reference genome



```
bwa mem -M -t 12 reference_genome.fasta ind1_1.fq.gz ind1_2.fq.gz | samtools view -q 20 -b | samtools sort -@ 12 -T ind1 > ind1_sorted.bam
```

### Alignment quality stats


```
samtools depth ind1_sorted.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print sum/NR,"\t", sqrt(sumsq/NR - (sum/NR)**2)}' > ind1_bamdepth.txt
```

```
samtools flagstat -O tsv ind1_sorted.bam > ind1_bamstats.txt
```

![Alt text](Figures/Figure_SummaryStats.png?raw=true "Title")


### Alignment cleaning

PicardTools was used to clean bam files (to set soft-clipping for reads beyond end of the reference, and MAPQ to 0 for unmapped reads). PCR duplicates were not marked for removal.

```
java -jar picard.jar CleanSam INPUT=ind1.bam OUTPUT=ind1.cleaned.bam 2>ind1.cleaned.bam.log
```

PicardTools was used to add read groups to each individual. RGID is the read group ID, (i.e. the name of the individual). RGLB is the read group library, (i.e. the sequencing lane number). RGPU is the read group platform unit (i.e. run barcode of the lane). RGSM is the read group sample name (i.e. the name of the individual).

```
java -jar picard.jar AddOrReplaceReadGroups INPUT=ind1.cleaned.bam OUTPUT=ind1.sort.rg.bam SORT_ORDER=coordinate RGID=ind1 RGLB=lib1 RGPL=ILLUMINA RGPU= lane1_parentals RGSM=ind1 CREATE_INDEX=True 2>ind1.sort.rg.log
```

### SNP calling


```
./freebayes -f reference.fasta ind1.sort.rg.bam ind2.sort.rg.bam --use-best-n-alleles 4 --report-monomorphic --genotype-qualities > pop1.vcf
```



## References

James ME, Arenas-Castro H, Groh JS, Engelstädter J, Allen SL, Ortiz-Barrientos D. 2021. Highly replicated evolution of parapatric ecotypes. Molecular Biology and Evolution, In press.
