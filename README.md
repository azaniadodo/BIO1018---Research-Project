# BIO1018: Research Project
Examining The Impact Of Knocking Out The Alzheimer’s Disease Susceptibility Gene TREM2 On Gene Expression In Human iPS-Derived Microglial Cells

Detailed transcript of the code and methodology employed to obtain my results, using the supplied data. All procedures were applied to both the primary and secondary (comparator dataset)

## At the Commandline

## 1. Pre-Alignment Quality Control
To ensure suitability for further RNA-Seq processing, the raw FASTQ sequence data for each sample were quality-assessed using FastQC (v 0.12.1), generating results plots for several quality control metrics

###### Tool Used: FastQC - https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt

#### Step 1: Move all raw fastq files into a new directory called  raw_fastq
    mkdir -p /mnt/data/GCB2025/username/raw_fastq/
#### Step 2: Run FastQC on multiple samples in the  directory (using 4 threads) and send outputs to a new designated directory
    fastqc -t 4 *.fastq.gz --outdir=/mnt/data/GCB2025/username/raw_fastq/FastQC
#### Step 3: Collate FastQC results of each sample using MultiQC and send output to a new designated directory
###### Tool Used: MultiQC - https://multiqc.info/
    multiqc /mnt/data/GCB2025/username/raw_fastq/FastQC


## 2. Alignment to Reference Genome
To facilitate all further downstream processes the raw RNA-sequencing FASTQ data for each sample were aligned to the Homo sapiens (Human) genome assembly GRCh38 (hg38) reference genome using the splice-aware aligner STAR (v2.7.11b). The reference genome was obtained from the UCSC Genome Browser. This process generated multiple output files, notably the BAM output file. 

###### Tool Used: STAR - https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

#### Step 1: Generate Genome Indices
      STAR --runThreadN 16 --runMode genomeGenerate --genomeDir /mnt/data/Reference_genomes/Human/UCSC/STAR --genomeFastaFiles \
      /mnt/data/Reference_genomes/Human/UCSC/hg38.fa --sjdbGTFfile /mnt/data/Reference_genomes/Human/UCSC/hg38.ncbiRefSeq.gtf --sjdbOverhang 100
#### Step 2: Aligne all fastq files within a directory to the reference genome
    for i in *_1.fastq.gz; do 
    STAR --genomeDir /mnt/data/Reference_genomes/Human/UCSC/STAR/ --runThreadN 12 --readFilesIn $i ${i%_1.fastq.gz}_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix ${i%_1.fastq.gz} --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard
    done
#### Step 3:  Move all BAM and related files into a new directory and create a MultiQC summary plot of the results 
      mkdir /mnt/data/GCB2025/username/STAR_files
      mv *[.bam,.out,.tab] /mnt/data/GCB2025/username/STAR_files

      multiqc /mnt/data/GCB2025/username/STAR_files

      
## 3. Post-Alignment Quality Control
###### Tool Used: Samtools - http://www.htslib.org/doc/samtools.html 
###### Tool Used: RSeQC - http://rseqc.sourceforge.net/

### 3.1 Generate Mapping Statistics
For each sample, the number of alignments to the reference genome was counted from the coordinate-sorted BAM output files using Samtools View -c (v1.21+htslib-1.21). The strandedness of the RNA-Seq data was then inferred using the BAM file from one sample using infer_experiment.py (v5.0.4). The mapping statistics for each BAM file were summarised, including the number of duplicated reads, read pairing and total reads mapped using bam_stat.py (v5.0.4)

### 3.2 Assess Mapping Quality
The mapping quality of each alignment file was evaluated through numerous RNA-Seq QC alignment processes. First, the STAR alignments MultiQC summary plot (Appendix X) was analysed for each dataset to assess the proportion of read mapping and overall mapping success. Next, the distribution of read mapping across genomic features (e.g. exons, introns, unannotated regions)  was examined. Based on the alignment counts obtained, a subsample of approximately 200,000 aligned reads was extracted from each coordinate-sorted BAM file using Samtools View -s. The subsampled BAM files were then indexed using Samtools index, producing indexed BAM files (BAI). The BAI files were aligned to the BED12 format of the reference genome to determine the consistency of  5’ to 3’ gene body coverage using geneBody_coverage.py (v5.0.4). 

 #### Step 1: Count number of alignments for each sample
     samtools view -c SampleNameAligned.sortedByCoord.out.bam
 #### Step 2: Index BAM files (using 12 threads)
     for i in *Aligned.sortedByCoord.out.bam; do
     samtools index $i -@ 12
     done
 #### Step 3: Convert GTF to genePred format and then to BED12 format
     /mnt/data/Reference_genomes/Human/UCSC/hg38.ncbiRefSeq.gtf /mnt/data/Reference_genomes/Human/UCSC/hg38.ncbiRefSeq.genePred

     /mnt/data/Reference_genomes/Human/UCSC/hg38.ncbiRefSeq.genePred /mnt/data/Reference_genomes/Human/UCSC/hg38.ncbiRefSeq.bed12
 #### Step 4: Infer Strandness of RNA-seq data (only performed on one BAM file)
     infer_experiment.py -r
     /mnt/data/Reference_genomes/Human/UCSC/hg38.ncbiRefSeq.bed12 -i SampleNameAligned.sortedByCoord.out.bam
#### Step 5:  Create RSeq Directory for output files
    mkdir /mnt/data/GCB2025/username/rseqc
#### Step 6: Summarise mapping statistics of a BAM or SAM file
###### Tool Used: bam_stat.py
    for i in *Aligned.sortedByCoord.out.bam; do
    bam_stat.py -i $i > ~/rseqc/${i%Aligned.sortedByCoord.out.bam}.bamstats.txt
    done    
#### Step 7: Calculate how mapped reads were distributed over genome features (like CDS exon, 5’UTR exon, 3’UTR exon, Intron, Intergenic regions)
###### Tool Used: read_distribution.py
     for i in *Aligned.sortedByCoord.out.bam; do
     read_distribution.py -i $i -r /mnt/data/Reference_genomes/Human/UCSC/hg38.ncbiRefSeq.bed12 >  ~/rseqc/${i%Aligned.sortedByCoord.out.bam}.read_dist.txt
     done 
##### Step 8: Subsample a proportion of the aligned reads, to end up with ~200,000
      for i in *Aligned.sortedByCoord.out.bam; do
      samtools view -s 0.08 -o ${i%.bam}_subset.bam $i
      done
##### Step 9: index subsampled BAM files
    for i in *Aligned.sortedByCoord.out_subset.bam; do
    samtools index $i
    done
#### Step 10:  Count number of alignments in the subsampled BAM files
    for i in *Aligned.sortedByCoord.out_subset.bam; do
    samtools view -c $i
    done
#### Step 11: Calculate the RNA-seq reads coverage over gene body
###### Tool Used: geneBody_coverage.py
    for i in *Aligned.sortedByCoord.out_subset.bam; do
    geneBody_coverage.py -i $i -r /mnt/data/Reference_genomes/Human/UCSC/hg38.ncbiRefSeq.bed12 -o ~/rseqc/${i%Aligned.sortedByCoord.out_subset.bam}
    done
#### Step 12: Collate RSeQC results of each sample for each QC analysis (Steps 6, 7 and 11) using MultiQC 
    multiqc /mnt/data/GCB2025/username/rseqc

    
## 4. Gene Expression Quantification
Each coordinated-sorted BAM file was resorted by name using Samtools Sort -n. Alignments of the name-sorted BAM files to the GTF file format of the reference genome were performed, further providing information on the coordinates of each gene/exon. The mapped reads overlapping genomic features in the reference were assigned to said genomic feature and expression levels (i.e. number of reads) counted using featureCounts (v2.0.6) — an essential preliminary procedure for downstream differential gene expression analysis. This step was altered for a reverse-stranded library to distinguish between genes with different coding strands (antisense and sense).

###### Tool Used: featurecounts - http://subread.sourceforge.net/featureCounts.html
##### Step 1: Create a name-sorted BAM file directory
    mkdir /mnt/data/GCB2025/username/name_sorted_bams
#### Step 2: Create name-sorted versions of BAM files, as is required for featureCounts
    for i in *Aligned.sortedByCoord.out.bam; do
    samtools sort -n -@ 12 -o ~/name_sorted_bams/${i%Aligned.sortedByCoord.out.bam}_namesorted.bam $i
    done
 #### Step 3: Perform featureCounts for a reverse-stranded library
     featureCounts -T 12 -s 2 -p --countReadPairs -C -a /mnt/data/Reference_genomes/Human/UCSC/hg38.ncbiRefSeq.gtf -o featurecounts.txt *namesorted.bam 2> featurecounts.screen-output.log

## In R
## Download Packages
```
        if (!require("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager") }

        if (!require("apeglm")) {
          BiocManager::install("apeglm") }
        library(apeglm)
        
        
        if (!require("DESeq2")) {
          BiocManager::install("DESeq2") }
        library(DESeq2)
        packageVersion("fgsea")
        
        
        if (!require("tidyverse")) {install.packages("tidyverse") }
        library(tidyverse)
        
        
        if (!require("RColorBrewer")) {
          install.packages("RColorBrewer") }
        library(RColorBrewer)
        
        
        if (!require("pheatmap")) {
          install.packages("pheatmap") }
        library(pheatmap)
        
        
        if (!require("ggrepel")) {
          install.packages("ggrepel") }
        library(ggrepel)
        
        
        if (!require("msigdbr")) {
          install.packages("msigdbr") }
        library(msigdbr)
        
        if (!require("tibble")) {
          install.packages("tibble") }
        library(tibble)
        
        if (!require("fgsea")) {
          BiocManager::install("fgsea") }
        library(fgsea)
```
## 5. Differential Gene Expression 
#### For Primary Dataset
### 5.1 Data Normalisation and Variance Analysis
The raw gene/exon counts obtained from the featureCounts output were normalised for gene expression level variability due to covariates or batch effects. The package DESeq2 v1.46.0 was utilised in R (v4.4.2) via R Studio (v2024.12.1+563) to normalise the counts based on differences in sequencing depth and sample composition for each dataset.
To assess the overall between-sample and between-replicate similarity and deduce the major sources of variation within the RNA-Seq data, a principal component analysis (PCA) and a hierarchical clustering analysis (HCA) of the regularised-logarithm transformed (rlog) normalised count data were performed for each dataset, using the plotPCA() and pheatmap() functions, respectively. 

```
counts <- read.table("featurecounts.txt", header = TRUE)
View(counts)
row.names(counts) <- counts$Geneid
head(counts)

counts <- counts[ , -c(1:6) ]
head(counts)

orig_names <- names(counts)
orig_names
names(counts) <- gsub("(PersonA_KO_rep|PersonB_KO_rep|PersonA_WT_rep|PersonB_WT_rep)([1-9]).*" , "\\1\\2" , orig_names)
names(counts)
head(counts)

metadata <- data.frame(genotype=gsub("(PersonA_|PersonB_)(KO|WT)(_rep[1-9])", "\\2", names(counts)), row.names = names(counts))
metadata$individual <-gsub("(Person)(A|B)(_KO|_WT)(_rep[1-9])", "\\2", names(counts)) #adding a column for individual ID
metadata

class(counts)
class(metadata)

# Check that sample names match in both files
all(colnames(counts) %in% rownames(metadata))

all(colnames(counts) == rownames(metadata))

dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ individual + genotype) #accounting for both factors
colData(dds) %>% head
colData(dds)
assay(dds, "counts" ) %>% head
dds$genotype <- factor(dds$genotype, levels = c("WT","KO"))
str(dds$genotype) # added just to check the WT was successfully made to be the reference

# Check number of rows (i.e number of genes) and remove rows/genes where there are zero counts in every sample
dim(dds)
dim(dds[rowSums(counts(dds)) > 0, ])
dds <- dds[rowSums(counts(dds)) > 0, ]

# Perform a regularised-logarithm transformation of the count data to make it more homoskedastic
rld <- rlog(dds, blind=TRUE)

# Compare counts in original DESeqDataSet object (dds) with transformed DESeqDataSet object (rld)
assay(dds) %>% head
assay(rld) %>% head

# Create a PCA plot
p <-plotPCA(rld, intgroup=c("genotype","individual")) + geom_text_repel(aes(label=name), size=2) 
p$layers[[1]]$aes_params$size <- 1
p + guides(color=guide_legend(override.aes = list(size =3))) 

# Create a hierarchical clustering plot
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$genotype, rld$individual,sep = "-")
colnames(sampleDistMatrix) <- paste(rld$genotype, rld$individual,sep = "-") #added column names to help better interpret the plot 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
```

### 5.2 DGE Analaysis
DESeq2 was run on the raw genomic feature counts to identify any differential expression of the genes assigned in the featureCounts file.
To achieve this, a count matrix and tabular metadata dataframe were created. A design formula was specified to measure differential expression based on the effect of interest — genotype; WT/WT-2, KO/KO-2 or R47H — whilst accounting for covariates of (A) individual sequenced (for the primary dataset) and, (B) the week sequencing was performed (for the comparator dataset). 

```
# Run the differential expression pipeline on the raw counts
dds <- DESeq(dds)

# Obtain results for specific comparisons
KOvsWT <-results(dds, contrast=c("genotype","KO","WT"), independentFiltering = TRUE , alpha = 0.05)

summary(KOvsWT)
head(KOvsWT)
table(KOvsWT$padj < 0.05)

KOvsWT[order(KOvsWT$padj),] %>% head

# Write the results to a file
write.table(as.data.frame(KOvsWT[order(KOvsWT$padj),] ), file="KO_vs_WT_dge.txt",sep="\t", quote = FALSE)
```

### 5.3 Visualisation
The normalised count levels of the GOI — TREM2 — were plotted using the ggplot2() function to assess the differences according to genotype (WT, KO and R47H) for each dataset.  The DGE  levels for each gene were collectively assessed and visualised through a volcano plot of log2fc vs P-adj,  labelled with the Top 30 most significantly differentially expressed genes, using the ggplot() function. The magnitude of DGE change was further visualised by generating a heatmap of the Top 30 significantly differentially expressed genes using the pheatmap() function. These processes were performed on the KO vs WT/ KO-2 vs WT-2 DGE level comparison results, using their respective Top 30 subset of significantly differentially expressed genes; these were derived from the DESeq2 results of the primary and comparator datasets and computed into tibble dataframes. Likewise, the same procedure was carried out for the R47H vs WT-2 DGE level comparison tibble dataframe, using its Top 30 significantly differentially expressed genes.  All data visualisation steps were performed in R/Rstudio, using a significance threshold value of Padj < 0.05
```
# Plot expression for a single gene by genotype
plotCounts(dds, gene="TREM2", intgroup="genotype") # are counts higher in WTs or KOs

# Plot expression for a single gene by individual
plotCounts(dds, gene="TREM2", intgroup="individual") #checking if the gene expression patterns are different between person A and B 

# Plot expression for a single gene by individual and genotype
plotCounts(dds, gene="CLU", intgroup=c("individual","genotype") ) #consistency of WT being higher than KO in both individuals

# Identify the names of the top 30 significantly differentially expressed genes by genotype
KOvsWT_sorted <- KOvsWT[order(KOvsWT$padj),]
KOvsWT_top30 <- head(KOvsWT_sorted, n=30)
top30_genes_KOvsWT <- rownames(KOvsWT_top30)
top30_genes_KOvsWT

# Extract rlog-transformed values into a matrix and create a heatmap
rlog.dge <- rld[top30_genes_KOvsWT,] %>% assay
pheatmap(rlog.dge, scale="row", main = "Differential Gene Expression (row-based z-score) for KO vs WT")

# Create a tibble dataframe of results
KOvsWT_tb <- KOvsWT %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble()
KOvsWT_tb[order(KOvsWT_tb$gene),]

# Obtain logical vector where TRUE values denote padj values < 0.05; mutate() adds new variables and preserves existing ones
KOvsWT_tb <- KOvsWT_tb %>% mutate(threshold = padj < 0.05)
KOvsWT_tb <- KOvsWT_tb %>% arrange(padj) %>% mutate(genelabels = "")
KOvsWT_tb$genelabels[1:30] <- KOvsWT_tb $gene[1:30] #labelling top 30 most significant genes
KOvsWT_tb$severity[KOvsWT_tb$log2FoldChange > 0 & KOvsWT_tb$padj < 0.05] <- "Upregulated" 
KOvsWT_tb$severity[KOvsWT_tb$log2FoldChange < 0 & KOvsWT_tb$padj < 0.05] <- "Downregulated"
KOvsWT_tb$severity[KOvsWT_tb$padj > 0.05] <- "Not Significant"
KOvsWT_tb$dataset <- "Original"
KOvsWT_significant <- KOvsWT_tb %>% filter(padj < 0.05) #create subset dataframe of genes with a significant p-adj value 

## Create a Volcano plot
ggplot(KOvsWT_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = severity)) +
  geom_text_repel(aes(label = genelabels,colour = severity)) + 
  scale_color_manual(values = c("blue","grey","red"), limits=c("Downregulated", "Not Significant", "Upregulated")) +
  guides(color=guide_legend(override.aes = list(size =3)))+
  ggtitle("KO vs WT") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))
```

####  For Comparator Dataset

**Replicate procedures in section 5, making minor variations to align with he supplied data.**
- Account for factors of  week and genotype (rather than individual and genotype) during data normalisation and variance analysis steps
- Obtain DGE results for and visualise comparisons of KO-2 vs WT-2, **and** R47H vs WT-2 after running DESeq
- **Optional**: Obtain DGE analysis results, Top 30 Most Significantly Differentially Expressed genes for the KO-2 vs R47H comparison and visualise as a DGE Heatmap for further analysis and understanding of gene variant effects

        rlog.dge2C <- rld2[top30_genes_KO2vsR47H ,c(1:2,4:5,7:8)] %>% assay
          pheatmap(rlog.dge2C, scale="row",  main = "Differential Gene Expression KO-2 vs R47H (row-based z-score)")
  
Comparator Code:
```
counts2 <- read.table("featurecounts_comparator.txt", header = TRUE)
row.names(counts2) <- counts2$Geneid
counts2 <- counts2[ , -c(1:6) ]

orig_names2 <- names(counts2)
names(counts2) <- gsub("(pMac)(_week[7-9])(_KO|_WT|_R47H)(_.*)" , "\\1\\2\\3" , orig_names2)

metadata2 <- data.frame(genotype=gsub("(pMac_)(week[7-9]_)(WT|KO|R47H)", "\\3", names(counts2)), row.names = names(counts2))
metadata2$week <-gsub("(pMac_week)([7-9])(_KO|_WT|_R47H)", "\\2", names(counts2)) #adding a column for week number

dds2 <- DESeqDataSetFromMatrix(countData = counts2, colData = metadata2, design = ~ week + genotype) #accounting for both factors
colData(dds2) %>% head
colData(dds2)
assay(dds2, "counts" ) %>% head
dds2$genotype <- factor(dds2$genotype, levels = c("WT","KO","R47H"))
str(dds2$genotype) # added just to check the WT was successfully made to be the reference

# Check number of rows (i.e number of genes) and remove rows/genes where there are zero counts in every sample
dim(dds2)
dim(dds2[rowSums(counts(dds2)) > 0, ])
dds2 <- dds2[rowSums(counts(dds2)) > 0, ]

# Perform a regularised-logarithm transformation of the count data to make it more homoskedastic
rld2 <- rlog(dds2, blind=TRUE)

# Compare counts in original DESeqDataSet object (dds) with transformed DESeqDataSet object (rld)
assay(dds2) %>% head
assay(rld2) %>% head

# Create a PCA plot
p <-plotPCA(rld2, intgroup=c("genotype","week")) + geom_text_repel(aes(label=name), size=2) 
p$layers[[1]]$aes_params$size <- 1
p + guides(color=guide_legend(override.aes = list(size =3))) #not sure if this looks much better and the legend looks strange

# Create a hierarchical clustering plot
sampleDists <- dist(t(assay(rld2)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld2$genotype, rld2$week,sep = "-")
colnames(sampleDistMatrix) <- paste(rld2$genotype, rld2$week,sep = "-") #added column names to help me interpret the plot better
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# Run the differential expression pipeline on the raw counts
dds2 <- DESeq(dds2)

# Obtain results for specific comparisons
KO2vsWT2 <-results(dds2, contrast=c("genotype","KO","WT"), independentFiltering = TRUE , alpha = 0.05)
R47HvsWT2<-results(dds2, contrast=c("genotype","R47H","WT"), independentFiltering = TRUE , alpha = 0.05)
KO2vsR47H <-results(dds2, contrast=c("genotype","R47H","KO"), independentFiltering = TRUE , alpha = 0.05)

# Write the results to a file
write.table(as.data.frame(KO2vsWT2[order(KO2vsWT2$padj),] ), file="KO_vs_WT_dgeC.txt",sep="\t", quote = FALSE)
write.table(as.data.frame(R47HvsWT2[order(R47HvsWT2$padj),] ), file="R47H_vs_WT_dge.txt",sep="\t", quote = FALSE)
write.table(as.data.frame(KO2vsR47H[order(KO2vsR47H$padj),] ), file="R47H_vs_KO2_dge.txt",sep="\t", quote = FALSE)


# Plot expression for a single gene by genotype
plotCounts(dds2, gene="TREM2", intgroup="genotype") # are counts higher in WTs or KOs or R47H

# Plot expression for a single gene by individual
plotCounts(dds2, gene="TREM2", intgroup="week") #checking if the gene expression patterns are different per week 


# Identify the names of the top 30 significantly differentially expressed genes by genotype
KO2vsWT2_sorted <- KO2vsWT2[order(KO2vsWT2$padj),]
R47HvsWT2_sorted <- R47HvsWT2[order(R47HvsWT2$padj),]
KO2vsR47H_sorted <- KO2vsR47H[order(KO2vsR47H$padj),]
KO2vsWT2_top30 <- head(KO2vsWT2_sorted, n=30)
R47HvsWT2_top30 <- head(R47HvsWT2_sorted, n=30)
KO2vsR47H_top30 <-head(KO2vsR47H_sorted, n=30)
top30_genes_KO2vsWT2 <- rownames(KO2vsWT2_top30)
top30_genes_KO2vsWT2
top30_genes_R47HvsWT2 <- rownames(R47HvsWT2_top30)
top30_genes_R47HvsWT2
top30_genes_KO2vsR47H <-rownames(KO2vsR47H_top30)
top30_genes_KO2vsR47H

# Extract rlog-transformed values into a matrix and generate DGE heatmaps that only show specific comparison data
rlog.dge2 <- rld2[top30_genes_KO2vsWT2,c(1,3,4,6,7,9)] %>% assay
pheatmap(rlog.dge2, scale="row", main = "Differential Gene Expression KO-2 vs WT-2 (row-based z-score)")
rlog.dge2B <- rld2[top30_genes_R47HvsWT2,c(2:3,5:6,8:9)] %>% assay
pheatmap(rlog.dge2B, scale="row",  main = "Differential Gene Expression R47H vs WT-2 (row-based z-score)")

# Create a tibble dataframe of results
KO2vsWT2_tb <- KO2vsWT2 %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble()
R47HvsWT2_tb <- R47HvsWT2 %>% data.frame() %>% rownames_to_column(var ="gene") %>% as_tibble()

#create subset dataframe of annotated genes with a significant p-adj value 
KO2vsWT2_significant <- KO2vsWT2_tb %>% filter(padj <0.05)
KO2vsWT2_annotated <- KO2vsWT2_tb %>% filter(gene != "" | gene != NA) 
R47HvsWT2_significant <- R47HvsWT2_tb %>% filter(padj <0.05)
R47HvsWT2_annotated <- R47HvsWT2_tb %>% filter(gene != "" | gene != NA)

# Obtain logical vector where TRUE values denote padj values < 0.05; mutate() adds new variables and preserves existing ones
KO2vsWT2_tb <- KO2vsWT2_tb %>% mutate(threshold = padj < 0.05)
KO2vsWT2_tb <- KO2vsWT2_tb %>% arrange(padj) %>% mutate(genelabels = "")
KO2vsWT2_tb$genelabels[1:30] <- KO2vsWT2_tb $gene[1:30] #labelling top 30 most significant genes
KO2vsWT2_tb$severity[KO2vsWT2_tb$log2FoldChange > 0 & KO2vsWT2_tb$padj < 0.05] <- "Upregulated" 
KO2vsWT2_tb$severity[KO2vsWT2_tb$log2FoldChange < 0 & KO2vsWT2_tb$padj < 0.05] <- "Downregulated"
KO2vsWT2_tb$severity[KO2vsWT2_tb$padj > 0.05] <- "Not Significant"
KO2vsWT2_tb$dataset <- "Comparator"

R47HvsWT2_tb <- R47HvsWT2_tb %>% mutate(threshold = padj < 0.05)
R47HvsWT2_tb <- R47HvsWT2_tb %>% arrange(padj) %>% mutate(genelabels = "")
R47HvsWT2_tb$genelabels[1:30] <- R47HvsWT2_tb $gene[1:30] #labelling top 30 most significant genesR
R47HvsWT2_tb$severity[R47HvsWT2_tb$log2FoldChange > 0 & R47HvsWT2_tb$padj < 0.05] <- "Upregulated" 
R47HvsWT2_tb$severity[R47HvsWT2_tb$log2FoldChange < 0 & R47HvsWT2_tb$padj < 0.05] <- "Downregulated"
R47HvsWT2_tb$severity[R47HvsWT2_tb$padj > 0.05] <- "Not Significant"
R47HvsWT2_tb$dataset <- "R47H"

## Create a Volcano plot
ggplot(KO2vsWT2_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = severity)) +
  geom_text_repel(aes(label = genelabels,colour = severity)) + #when I colour the labels by severity it changes the legend mark to the letter 'a' 
  scale_color_manual(values = c("blue","grey","red"), limits=c("Downregulated", "Not Significant", "Upregulated")) +
  guides(color=guide_legend(override.aes = list(size =3)))+
  ggtitle("KO-2 vs WT-2 ") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 

ggplot(R47HvsWT2_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = severity)) +
  geom_text_repel(aes(label = genelabels,colour = severity)) + #when I colour the labels by severity it changes the legend mark to the letter 'a' 
  scale_color_manual(values = c("blue","grey","red"), limits=c("Downregulated", "Not Significant", "Upregulated")) +
  guides(color=guide_legend(override.aes = list(size =3)))+
  ggtitle("R47H vs WT-2") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 

```

## 6. Gene Data Filtering
All significantly differentially expressed genes were extracted from each comparison tibble dataframe, i.e. KO vs WT, KO-2 vs WT-2 and R47H vs WT-2, using a significance threshold of P-adj < 0.05. These subset dataframes were combined to create two merged dataframes of significantly differentially expressed genes derived from (1) KO vs WT and KO-2 vs WT-2 comparisons, and (2) KO vs WT and R47H vs WT-2 comparisons. The merged dataframes were filtered again to select only genes that were significantly differentially expressed and showing the same direction of differential expression (+/-log2fc ), i.e. genes significantly upregulated or downregulated in both of the comparison dataframes they were derived from. The two final filtered datasets were searched for known Alzheimer’s susceptibility genes  and DAM-related genes.

```
# Join the two tibbles  and filter for only significant annotated genes (Comparison 1)
KOvsWT_joined_tb <- full_join(x= KOvsWT_significant, y= KO2vsWT2_annotated, join_by(gene)) %>% na.omit(KOvsWT_joined_tb)
KOvsWT_log2FC_tb <- KOvsWT_joined_tb %>% filter(log2FoldChange.x >0 & log2FoldChange.y>0|log2FoldChange.x<0 &log2FoldChange.y<0) #significsntly differentially expressed genes in both comparisons
KOvsWT_log2FC_Padj_tb<- KOvsWT_joined_tb %>% filter(log2FoldChange.x >0 & log2FoldChange.y>0|log2FoldChange.x<0 &log2FoldChange.y<0) %>% filter(padj.y<=0.1) # significantly differentially expressed genes showing the same direction of differential expression (+/-log2fc ) in both comparisons

# Visualise filtered data
ggplot(KOvsWT_joined_tb, aes(x=log2FoldChange.x, y=log2FoldChange.y)) +
  geom_point()+
  xlab("log2FC Primary Dataset (KO vs WT)") + 
  ylab("log2FC Comparator Dataset (KO-2 vs WT-2)") +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 

ggplot(KOvsWT_log2FC_Padj_tb, aes(x=log2FoldChange.x, y=log2FoldChange.y)) +
  geom_point()+
  xlab("log2FC Primary Dataset (KO vs WT)") + 
  ylab("log2FC Comparator Dataset (KO-2 vs WT-2)") +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 

AD_GWASgenes <- read.csv("AD_Loci.csv") #Import Alzheimer's Disease Susceptibility loci
AD_GWAS_justgenes <- AD_GWASgenes %>% select(Name)
colnames(AD_GWAS_justgenes) <- "gene"

#Compare subset
KOvsWT_justgenes <- KOvsWT_log2FC_tb %>% select(gene)
semi_join(x= KOvsWT_justgenes, y=AD_GWAS_justgenes) %>% print(n=26) 
#26 of the AD genes are either upregulated/downregulated in the Original dataset and the comparator dataset

KOvsWT_justgenes_final <- KOvsWT_log2FC_Padj_tb %>% select(gene)
semi_join(x= KOvsWT_justgenes_final, y=AD_GWAS_justgenes) 
#14 of the AD genes are either upregulated/downregulated in the Original dataset and the comparator dataset and  significantly in both datasets

DAM_genes <- read.csv("DAM_genes.csv") #Import Disease-associated microglia loci
semi_join(x= KOvsWT_justgenes, y=DAM_genes) 
#9 of the DAM genes are either upregulated/downregulated in the Original dataset and the comparator dataset

semi_join(x= KOvsWT_justgenes_final, y=DAM_genes) 
#5 of the DAM genes are either upregulated/downregulated in the Original dataset and the comparator dataset and  significantly in  both datasets
```

### **Replicate procedures in section 6 for Comparison 2 (KO vs WT and R47H vs WT-2)**
Comparator Code:
```
KOvsWTvsR47H_joined_tb <- full_join(x= KOvsWT_significant, y= R47HvsWT2_annotated, join_by(gene)) %>% na.omit(KOvsWTvsR47H_joined_tb)
KOvsWTvsR47H_log2FC_tb <- KOvsWTvsR47H_joined_tb %>% filter(log2FoldChange.x >0 & log2FoldChange.y>0|log2FoldChange.x<0 &log2FoldChange.y<0)
KOvsWTvsR47H_log2FC_Padj_tb<- KOvsWTvsR47H_joined_tb %>% filter(log2FoldChange.x >0 & log2FoldChange.y>0|log2FoldChange.x<0 &log2FoldChange.y<0) %>% filter(padj.y<=0.1)
view(KOvsWTvsR47H_log2FC_Padj_tb)

#allsignificant_joined_tb <- full_join(x= KOvsWT_log2FC_Padj_tb, y= KOvsWTvsR47H_log2FC_Padj_tb, join_by(gene)) %>% na.omit(allsignificant_joined_tb)


ggplot(KOvsWTvsR47H_joined_tb, aes(x=log2FoldChange.x, y=log2FoldChange.y)) +
  geom_point()+
  xlab("log2FC Primary Dataset (KO vs WT)") + 
  ylab("log2FC Comparator Dataset (R47H vs WT-2)") +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 

ggplot(KOvsWTvsR47H_log2FC_Padj_tb, aes(x=log2FoldChange.x, y=log2FoldChange.y)) +
  geom_point()+
  xlab("log2FC Primary Dataset (KO vs WT)") + 
  ylab("log2FC Comparator Dataset (R47H vs WT-2)") +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 

all(AD_GWASgenes$Name %in% KOvsWTvsR47H_log2FC_Padj_tb)
all(KOvsWTvsR47H_log2FC_Padj_tb %in% AD_GWASgenes$Name)

KOvsWTvsR47H_justgenes <- KOvsWTvsR47H_log2FC_tb %>% select(gene)
semi_join(x= KOvsWTvsR47H_justgenes, y=AD_GWAS_justgenes) %>% print(n=22)
#22 of the AD genes are either upregulated/downregulated in the Original dataset and the variant dataset

KOvsWTvsR47H_justgenes_final <- KOvsWTvsR47H_log2FC_Padj_tb %>% select(gene)
semi_join(x= KOvsWTvsR47H_justgenes_final, y=AD_GWAS_justgenes) 
#6 of the AD genes are either upregulated/downregulated in the Original dataset and the variant dataset and  significantly in the both datasets

semi_join(x= KOvsWTvsR47H_justgenes, y=DAM_genes) 
#9 of the DAM genes are either upregulated/downregulated in the Original dataset and the variant dataset

semi_join(x= KOvsWTvsR47H_justgenes_final, y=DAM_genes) 
#2 of the DAM genes are either upregulated/downregulated in the Original dataset and the variant dataset and  significantly in  both datasets

```

## 7. Functional Enrichment Analysis
To distinguish which pathways or gene networks the significantly differentially expressed genes were involved in, and their role in the disease phenotype, a functional enrichment analysis, also known as a gene-set enrichment analysis or overrepresentation analysis, was performed on the two final filtered datasets using gProfiler’s g:GOSt tool. 
Each input list of significantly differentially expressed genes was mapped to functional data sources — Gene Ontology: Biological Processes (GO: BP), KEGG, Reactome and WikiPathways — and statistically significantly enriched biological pathways and processes were identified 


#### Step 1: writing final subset of significantly differentially expressed genes for Comparison 1 to a text file
##### In Primary Dataset
```
KOvsWT_2FC_Padj_justgenes <- KOvsWT_log2FC_Padj_tb$gene
write(KOvsWT_2FC_Padj_justgenes, file = "KOvsWT_2FC_Padj_justgenes.txt")
```
##### In Comparator Dataset
```
KOvsWTvsR47H_2FC_Padj_justgenes <- KOvsWTvsR47H_log2FC_Padj_tb$gene
write(KOvsWTvsR47H_2FC_Padj_justgenes, file = "KOvsWTvsR47H_2FC_Padj_justgenes.txt"
```
#### Step 2: Input list of significantly differentially expressed genes into gProfiler’s g:GOSt tool
