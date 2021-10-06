# Estimation of paternity share from PoolSeq data

This repository contains step-by-step instructions to estimate the fertilisation outcome of competitive crosses, where two sires from different populations simultaneously mate or pollinate a dam. The paternity share of the offspring is inferred from the allelic frequency of a set of markers differentially fixed in the sire populations and monomorphic in the dam population for either allele. 



## Parental individuals - Calling SNPs



### Indexing the reference genome


```
bwa index ${TMPDIR}/Senecio.contigs.fasta
```

Senecio reference genome version: SPD_CN1K_CtgRN


### Alignmening the parental reads to the reference genome



```
bwa mem -M -t 12 reference_genome.fasta ind1_1.fq.gz ind1_2.fq.gz | samtools view -q 20 -b | samtools sort -@ 12 -T ind1 > ind1_sorted.bam
```

### Stats


```
samtools depth ind1_sorted.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print sum/NR,"\t", sqrt(sumsq/NR - (sum/NR)**2)}' > ind1_bamdepth.txt
```

```
samtools flagstat -O tsv ind1_sorted.bam > ind1_bamstats.txt
```


### 

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





