# BIO1018---Research-Project
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
