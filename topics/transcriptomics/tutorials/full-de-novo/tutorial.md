---
layout: tutorial_hands_on
enable: false

title: De novo transcriptome assembly, annotation, and differential expression analysis
zenodo_link: 'https://zenodo.org/record/3541678'
questions:
- Which biological questions are addressed by the tutorial?
- Which bioinformatics techniques are important to know for this type of data?
objectives:
- The learning objectives are the goals of the tutorial
- They will be informed by your audience and will communicate to them and to yourself
  what you should focus on during the course
- They are single sentences describing what a learner should be able to do once they
  have completed the tutorial
- You can use Bloom's Taxonomy to write effective learning objectives
time_estimation: 3H
key_points:
- The take-home messages
- They will appear at the end of the tutorial
contributors:
- abretaud
- lecorguille
- r1corre
- xiliu
- lleroi
- alexcorm
- paulineauffret

---


# Introduction
{:.no_toc}

<!-- This is a comment. -->

As a result of the development of novel sequencing technologies, the years between 2008 and 2012 saw a large drop in the cost of sequencing. Per megabase and genome, the cost dropped to 1/100,000th and 1/10,000th of the price, respectively. Prior to this, only transcriptomes of organisms that were of broad interest and utility to scientific research were sequenced; however, these developed in 2010s high-throughput sequencing (also called next-generation sequencing) technologies are both cost- and labor- effective, and the range of organisms studied via these methods is expanding.

Examining non-model organisms can provide novel insights into the mechanisms underlying the "diversity of fascinating morphological innovations" that have enabled the abundance of life on planet Earth. In animals and plants, the "innovations" that cannot be examined in common model organisms include mimicry, mutualism, parasitism, and asexual reproduction. De novo transcriptome assembly is often the preferred method to studying non-model organisms, since it is cheaper and easier than building a genome, and reference-based methods are not possible without an existing genome. The transcriptomes of these organisms can thus reveal novel proteins and their isoforms that are implicated in such unique biological phenomena.

[(source)](https://en.wikipedia.org/wiki/De_novo_transcriptome_assembly)

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Read cleaning (20 minutes)

Known sequencing biases:
- Unknown nucleotides (Ns)
- Bad quality nucleotides
- Hexamers biases (Illumina. Now corrected ?)

Why do we need to correct those?
- To remove a lot of sequencing errors (detrimental to the vast majority of assemblers)
- Because most de-bruijn graph based assemblers can’t handle unknown nucleotides

## Get data

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Create a new history for this tutorial
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the 12 `fq.gz` into a `List of Pairs` collection named `fastq_raw`
>    - Option 1: from a shared data library (ask your instructor)
>    - Option 2: from Zenodo using the URLs given below
>
>      [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3541678.svg)](https://doi.org/10.5281/zenodo.3541678)
>
>    ```
>    https://zenodo.org/record/3541678/files/A1_left.fq.gz
>    https://zenodo.org/record/3541678/files/A1_right.fq.gz
>    https://zenodo.org/record/3541678/files/A2_left.fq.gz
>    https://zenodo.org/record/3541678/files/A2_right.fq.gz
>    https://zenodo.org/record/3541678/files/A3_left.fq.gz
>    https://zenodo.org/record/3541678/files/A3_right.fq.gz
>    https://zenodo.org/record/3541678/files/B1_left.fq.gz
>    https://zenodo.org/record/3541678/files/B1_right.fq.gz
>    https://zenodo.org/record/3541678/files/B2_left.fq.gz
>    https://zenodo.org/record/3541678/files/B2_right.fq.gz
>    https://zenodo.org/record/3541678/files/B3_left.fq.gz
>    https://zenodo.org/record/3541678/files/B3_right.fq.gz
>    ```
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md collection=true collection_type="List of Pairs" collection_name="fastq_raw" %}
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
> 3. Rename the datasets
> 4. Check that the datatype
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="datatypes" %}
>
> 5. Add to each database a tag corresponding to ...
>
>    {% snippet faqs/galaxy/datasets_add_tag.md %}
>
{: .hands_on}

## Quality control

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **FastQC** {% icon tool %} with the following parameters:
>   - *"Short read data from your current history"*: `fastq_raw` (collection)
>
{: .hands_on}

<!-- ## Quality control with **MultiQC** - step 2/2

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **MultiQC** {% icon tool %} with the following parameters:
>    - In *"Results"*:
>        - {% icon param-repeat %} *"Insert Results"*
>            - *"Which tool was used generate logs?"*: `FastQC`
>                - In *"FastQC output"*:
>                    - *"Type of FastQC output?"*: `Raw data`
>                    - *"FastQC output"*: `data XX, data XX, and others (flattened)`
>
>    > ### {% icon comment %} Comment
>    >
>    > We agree that it's not comfortable. The wrapper of MultiQC must be improved
>    {: .comment}
>
{: .hands_on} -->

## Sequencing error correction with **Rcorrector**  

To do so, we can use [Rcorrector (Song et al., 2015)](https://github.com/mourisl/Rcorrector), which is a kmer-based error correction method for RNA-seq data, based on the path search algorithm. 

Rcorrector distinguishes among solid and non-solid k-mers as the basis for its correction algorithm. A solid k-mer is one that passes a given count threshold and therefore can be trusted to be correct. Rcorrector uses a flexible threshold for solid k-mers, which is calculated for each k-mer within each read sequence. At run time, Rcorrector scans the read sequence and, at each position, decides whether the next k-mer and each of its alternatives are solid and therefore represent valid continuations of the path. The path with the smallest number of differences from the read sequence, representing the likely transcript of origin, is then used to correct k-mers in the original read.

![Rcorrector_algo](../../images/full-de-novo/Rcorrector_path-search-algo.png)   

Path extension in Rcorrector. Four possible path continuations at the AGTC k-mer (k=4) in the De Bruijn graph for the r= AAGTCATAA read sequence. Numbers in the vertices represent k-mer counts. The first (top) path corresponds to the original read’s representation in the De Bruijn graph. The extension is pruned after the first step, AGTC →GTCA, as the count M(GTCA)=4 falls below the local cutoff (determined based on the maximum k-mer count (494) of the four possible successors of AGTC). The second path (yellow) has higher k-mer counts but it introduces four corrections, changing the read into AAGTCCGTC. The third path (blue) introduces only two corrections, to change the sequence into AAGTCGTTA, and is therefore chosen to correct the read. The fourth (bottom) path is pruned as the k-mer count for GTCT does not pass the threshold. Paths 2 and 3 are likely to indicate paralogs and/or splice variants of this gene.   
   

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Rcorrector** {% icon tool %} with the following parameters: 
>    - *"Is this library paired- or single-end?"*: `Paired-end (as collection)`
>        - *"FastQ file R1 (left)"*: `R1`
>        - *"FastQ file R2 (right)"*: `R2`
>        - *"Filter uncorrectable reads"*: `Yes`
>    - *"Additional options"*: `No`   
> 2. **Rename the resulting datasets**
>    - `RNA-seq Rcorrector on XX ->  XX_corrected`  
> 
{: .hands_on}

## rRNA removal with **Bowtie2**  

One of the main source of contamination of RNA-seq samples is **ribosomal RNA**. Indeed, ~90% of total RNA correspond to rRNA. Before sequencing, ribodepletion and polyA selection are common mmethods  to clean the samples, but it does not filter out all rRNA. After sequencing, removal rRNA reads from raw reads and detect rRNA transcripts are usefull processes.

To do so, we use [Bowtie 2](https://github.com/BenLangmead/bowtie2) whiwh is an ultrafast and memory-efficient tool for aligning sequencing reads to reference sequences. 
Here, the reference will be Silva database.

> ### {% icon hands_on %} Hands-on: Task description
> 1. **Bowtie2** {% icon tool %} with the following parameters:
>    - *"Is this single or paired library"*: `Paired-end`
>        - *"FASTA/Q file #1"*: `R1_corrected`
>        - *"FASTA/Q file #2"*: `R2_corrected`
>        - *"Write unaligned reads (in fastq format) to separate file(s)"*: `Yes`
>        - *"Write aligned reads (in fastq format) to separate file(s)"*: `No`
>        - *"Do you want to set paired-end options?"*: `No`
>    - *"Will you select a reference genome from your history or use a built-in index?"*: `Use a built-in genome index`
>        - *"Select reference genome"*: `Silva Ribosomal Database`
>    - *"Set read groups information?"*: `Do not set`
>    - *"Select analysis mode"*: `1:Default setting only`
>        - *"Do you want to use presets?"*: `Very sensitive end-to-end (--very-sensitive)`
>    - *"Do you want to tweak SAM/BAM Options?"*: `No`
>    - *"Save the bowtie2 mapping statistics to the history"*: `Yes`
> 2. **Rename the resulting datasets**
>    - `Bowtie2 on XX: alignments -> `     
>
{: .hands_on}

## Read cleaning with **Trimmomatic** 



> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Trimmomatic** {% icon tool %} with the following parameters:
>    - *"Single-end or paired-end reads?"*: `Paired-end (as collection)`
>    - *"Select FASTQ dataset collection with R1/R2 pair"*: `Raw Data`
>    - *"Perform initial ILLUMINACLIP step?"*: `Yes`
>    - *"Adapter sequences to use"*: `TruSeq3 (additional seqs) (paired-ended, for MiSeq and HiSeq)`
>    - In *"Trimmomatic Operation"*:
>        - {% icon param-repeat %} *"Insert Trimmomatic Operation"*
>            - *"Select Trimmomatic operation to perform"*: `Cut bases off end of a read, if below a threshold quality (TRAILING)`
>        - {% icon param-repeat %} *"Insert Trimmomatic Operation"*
t_toy_dataset](../../images/full-de-novo/ExN50_plot_toy_dataset.png)
t](../../images/full-de-novo/ExN50_plot.png)
>            - *"Select Trimmomatic operation to perform"*: `Cut bases off start of a read, if below a threshold quality (LEADING)`
>        - {% icon param-repeat %} *"Insert Trimmomatic Operation"*
>            - *"Select Trimmomatic operation to perform"*: `Sliding window trimming (SLIDINGWINDOW)`
>        - {% icon param-repeat %} *"Insert Trimmomatic Operation"*
>            - *"Select Trimmomatic operation to perform"*: `Drop reads with average quality lower than a specific level (AVGQUAL)`
>                - *"Minimum length of reads to be kept"*: `25`
>        - {% icon param-repeat %} *"Insert Trimmomatic Operation"*
>            - *"Select Trimmomatic operation to perform"*: `Drop reads below a specified length (MINLEN)`
>                - *"Minimum length of reads to be kept"*: `50`
>    - *"Output trimmomatic log messages?"*: `Yes`
> 2. **Rename** the Dataset Collection
>    - `Trimmomatic on collection XX: paired` -> `fastqc_cleaned`
>
>    > ### {% icon comment %} Comment
>    >
>    > You can check the Trimmomatic log files to get the number of read before and after the cleaning
>    > ```
>    > Input Read Pairs: 10000
>    > Both Surviving: 8804 (88.04%)
>    > Forward Only Surviving: 491 (4.91%)
>    > Reverse Only Surviving: 456 (4.56%) Dropped: 249 (2.49%)
>    > ```
>    {: .comment}
>
>    {% snippet faqs/galaxy/collections_rename.md %}
>
{: .hands_on}

## Quality control after cleaning

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **FastQC** {% icon tool %} with the following parameters:
>   - *"Short read data from your current history"*: `fastqc_cleaned` (collection)
>
{: .hands_on}

# Assembly (120 minutes - computing)

We will use *Trinity*, a de novo transcriptome assembler for short sequencing reads. 
*Trinity* is the most widely used de novo transcriptome assembler and in continuous development since several years.
All information about Trinity assembler are here [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)

## Assembly with **Trinity**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Trinity** {% icon tool %} with the following parameters:
>    - *"Are you pooling sequence datasets?"*: `Yes`
>        - *"Paired or Single-end data?"*: `Paired-end collection`
>            - *"Strand specific data"*: `No`
>    - *"Run in silico normalization of reads"*: `No`
>    - In *"Additional Options"*:
>        - *"Use the genome guided mode?"*: `No`
> 2. **Rename** the Trinity output
>    - `Trinity on data 52, data 51, and others: Assembled Transcripts` -> `transcriptome_raw.fasta`
>
>    {% snippet faqs/galaxy/datasets_rename.md %}
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. TODO
> 2. TODO
>
> > ### {% icon solution %} Solution
> >
> > 1. TODO
> > 2. TODO
> >
> {: .solution}
>
{: .question}


> ### {% icon comment %} Try it on!
> *rnaSPAdes* is a most recent assembler and can outperform *Trinity* results most of the time but not always. You can do
> the de novo assembly with **rnaSPAdes** {% icon tool %} and compare the results!
> 
> 
> > ### {% icon hands_on %} Hands-on: Task description
> >
> > 1. **rnaSPAdes** {% icon tool %} with the following parameters:
> >   - *"Single-end or paired-end short-reads"*: `Paired-end: individual datasets`
> >   - *"Select optional output file(s)"*: `Transcripts`
> >
> {: .hands_on}
{: .comment}

> ### {% icon tip %} Benchmarking
> [De novo transcriptome assembly: A comprehensive cross-species comparison of short-read RNA-Seq assemblers](https://academic.oup.com/gigascience/article/8/5/giz039/5488105)
> 
{: .comment}

# Assembly assessment / cleaning

## Checking of the assembly statistics

TODO: ***Trinity Statistics*** displays the summary statistics for a fasta file.

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Trinity Statistics** {% icon tool %} with the following parameters:
>    - *"Trinity assembly"*: `transcriptome_raw.fasta`
>
>    > ### {% icon comment %} Comment
>    > This step, even with this toy dataset, will take around 2 hours
>    {: .comment}
>
{: .hands_on}

> ### {% icon question %} Questions
>
> 1. TODO
>
> > ### {% icon solution %} Solution
> >
> > 1. TODO
> >
> {: .solution}
>
{: .question}


## Remapping on the raw transcriptome

TODO: presentation and aim of this part

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Align reads and estimate abundance** {% icon tool %} with the following parameters:
>    - *"Transcripts"*: `transcriptome_raw.fasta`
>    - *"Paired or Single-end data?"*: `Paired`
>        - *"Left/Forward strand reads"* -> `Multiple datasets`
>            - Click on the *Folder* button at the right
>                - *Type to Search*: `left`
>                - Select the 6 `Trimmomatic on ..._left.fq.gz`
>        - *"Right/Reverse strand reads"* -> `Multiple datasets`
>            - Click on the *Folder* button at the right
>                - *Type to Search*: `right`
>                - Select the 6 `Trimmomatic on ..._left.fq.gz`
>        - *"Strand specific data"*: `Yes`
>    - *"Abundance estimation method"*: `Salmon`
>    - In *"Additional Options"*:
>        - *"Trinity assembly?"*: `Yes`
> 2. **Rename** the 6 `* isoforms counts` :(
>    - Check in the information panel (**i** icon) the lineage of your file (ex: `A1_left.fq.gz` ... )
>    - Rename the datasets: `A1_raw`, `A2_raw`, `A3_raw`, `B1_raw`, `B2_raw`, `B3_raw`.
>
>    > ### {% icon comment %} Comment
>    >
>    > If you check at the Standard Error messages of your outputs. You can get the `Mapping rate`
>    > 1. Click on one dataset
>    > 2. Click on the little **i** icon
>    > 3. Click on *Tool Standard Error:	stderr*
>    > ```
>    > [2019-11-14 15:44:21.500] [jointLog] [info] Mapping rate = 44.4358%
>    > ```
>    {: .comment}
>
>    > ### {% icon comment %} Comment
>    >
>    > At this stage, you can now delete some useless datasets
>    > - `Trimmomatic on collection XX: unpaired`
>    > - `Align reads and estimate abundance on *: genes counts`
>    > Note that the dataset are just hidden. You can delete them permanently and make some room in the history options (the little wheel icon)
>    {: .comment}
>
>
{: .hands_on}

## Merge the mapping tables and compute normalizations

TODO: presentation and aim of this part

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Build expression matrix** {% icon tool %} with the following parameters:
>    - *"Abundance estimates"*: `A1_raw`, `A2_raw`, `A3_raw`, `B1_raw`, `B2_raw`, `B3_raw`
>    - *"Abundance estimation method"*: `Salmon`
>
{: .hands_on}

> ### {% icon question %} Questions
>
> What are the three tables?
>
> > ### {% icon solution %} Solution
> >
> > 1. `estimated RNA-Seq fragment isoform counts (raw counts)``
> > 2. `matrix of isoform TPM expression values (not cross-sample normalized)`
> > 3. `matrix of TMM-normalized expression values`
> >
> {: .solution}
>
{: .question}

## Compute contig Ex90N50 statistic and Ex90 transcript count

TODO: presentation and aim of this part

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Compute contig Ex90N50 statistic and Ex90 transcript count** {% icon tool %} with the following parameters:
>    - *"Expression matrix"*: `Build expression matrix: matrix of TMM-normalized expression values`
>    - *"Transcripts"*: `transcriptome_raw.fasta`
> 2. Click on the visulization icon on the dataset `Compute contig Ex90N50 statistic and Ex90 transcript count: ExN50 statistics`
>    1. **Scatterplot - Creates a 2D-scatterplot from tabular datapoints**
>    2. *"X Column"*: select the Columns `1`
>    3. *"Y Column"*: select the Columns `2`
>
{: .hands_on}

### Definition

Ex90N50 values are computed as usual N50 but limited to the top most highly expressed transcripts that represent 90% of the total normalized expression data. 

### What we get
![ExN50_plot_toy_dataset](../../images/full-de-novo/ExN50_plot_toy_dataset.png)

### What we should get with a real dataset
![ExN50_plot](../../images/full-de-novo/ExN50_plot.png)
[(source)](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats)

## Transcriptome annotation completeness

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Busco** {% icon tool %} with the following parameters:
>    - *"Sequence to analyse"*: `transcriptome_raw.fasta`
>    - *"Mode"*: `transcriptome`
>    - *"Lineage"*: `eukaryota_odb9`
>
{: .hands_on}

![RNASeq samples quality check Graphs](../../images/full-de-novo/rnaseq_samples_quality_check.png)

## Filter low expression transcripts

TODO: presentation and aim of this part

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Filter low expression transcripts** {% icon tool %} with the following parameters:
>    - *"Trinity assembly"*: `transcriptome_raw.fasta`
>    - *"Expression matrix"*: `Build expression matrix: matrix of isoform TPM expression values (not cross-sample normalized)`
>    - *"Minimum expression level required across any sample"*: `1.0`
>    - *"Isoform filtering method"*: `Keep all isoforms above a minimum percent of dominant expression`
>        - *"Minimum percent of dominant isoform expression"*: `1`
>
>    > ### {% icon comment %} Comment
>    >
>    > If you check at the Standard Error messages of your outputs. You can get the `Retained` rate
>    > 1. Click on one dataset
>    > 2. Click on the little **i** icon
>    > 3. Click on *Tool Standard Error:	stderr*
>    > ```
>    > 	Retained 2096 / 2102 = 99.71% of total transcripts.
>    > ```
>    {: .comment}
>
> 2. **Rename** the output
>    - `Filter low expression transcripts on data 42 and data 14: filtered low expression transcripts` -> `transcriptome_filtered.fasta`
>
{: .hands_on}

## Checking of the assembly statistics after cleaning

TODO: presentation and aim of this part

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Trinity Statistics** {% icon tool %} with the following parameters:
>    - *"Trinity assembly"*: `transcriptome_filtered.fasta`
>
{: .hands_on}


# Annotation
## Generate gene to transcript map

TODO: presentation and aim of this part

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Generate gene to transcript map** {% icon tool %} with the following parameters:
>    - *"Trinity assembly"*: `transcriptome_filtered.fasta`
>
{: .hands_on}

## Peptide prediction

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **TransDecoder** {% icon tool %} with the following parameters:
>    - *"Transcripts"*: `transcriptome_filtered.fasta`
>    - In *"Training Options"*:
>        - *"Select the training method"*: `Train with the top longest ORFs`
>
{: .hands_on}

## Similarity search

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Diamond** {% icon tool %} with the following parameters:
>    - *"What do you want to align?"*: `Align amino acid query sequences (blastp)`
>    - *"Input query file in FASTA or FASTQ format"*: `TransDecoder on data XXX: pep`
>    - *"Select a reference database"*: `Uniprot Swissprot`
>    - *"Format of output file"*: `BLAST Tabular`
>    - In *"Method to restrict the number of hits?"*: `Maximum number of target sequences`
>        - *"The maximum number of target sequence per query to report alignments for"*: `1`
> 3. **Rename** the Diamond output
>    - `Diamond on data XXX` -> `Diamond (blastp)`
> 2. **Diamond** {% icon tool %} with the following parameters:
>    - *"What do you want to align?"*: `Align DNA query sequences (blastx)`
>    - *"Input query file in FASTA or FASTQ format"*: `transcriptome_filtered.fasta`
>    - *"Select a reference database"*: `Uniprot Swissprot`
>    - *"Format of output file"*: `BLAST Tabular`
>    - In *"Method to restrict the number of hits?"*: `Maximum number of target sequences`
>        - *"The maximum number of target sequence per query to report alignments for"*: `1`
> 4. **Rename** the Diamond output
>    - `Diamond on data XXX` -> `Diamond (blastx)`
>
>    > ### {% icon comment %} Comment
>    >
>    > Note that you can both use **Diamond** {% icon tool %} or the **NCBI BLAST+ blastp** {% icon tool %} and **NCBI BLAST+ blast** {% icon tool %}
>    {: .comment}
>
{: .hands_on}

## Find signal peptides

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **SignalP 3.0** {% icon tool %} with the following parameters:
>    - *"Fasta file of protein sequences"*: `TransDecoder on data XXX: pep`
>
{: .hands_on}

## Find transmembrane domains

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **TMHMM 2.0** {% icon tool %} with the following parameters:
>    - *"FASTA file of protein sequences"*: `TransDecoder on data XXX: pep`
>
{: .hands_on}

## Search again profile database

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **hmmscan** {% icon tool %} with the following parameters:
>    - *"Sequence file"*: `TransDecoder on data XXX: pep`
>
{: .hands_on}

## Transcriptome annotation using **Trinotate**

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Trinotate** {% icon tool %} with the following parameters:
>    - *"Transcripts"*: `transcriptome_filtered.fasta`
>    - *"Peptides"*: `TransDecoder on data XXX: pep`
>    - *"Genes to transcripts map"*: `Generate gene to transcript map on data XXX: Genes to transcripts map`
>    - *"BLASTP: Peptides vs Uniprot.SwissProt"*: `Diamond (blastp)`
>    - *"BLASTX: Transcripts vs Uniprot.SwissProt"*: `Diamond (blastx)`
>    - *"HMMER hmmscan: Peptides vs PFAM"*: `Table of per-domain hits from HMM matches of TransDecoder on data XXX: pep against the profile database`
>    - *"TMHMM on Peptides"*: `TMHMM results`
>    - *"SignalP on Peptides"*: `SignalP euk results`
>    - *"Let Galaxy downloading the Trinotate Pre-generated Resource SQLite database"*: `Yes`
>
{: .hands_on}

# Differential Expression (DE) Analysis

TODO: presentation and aim of this part

## Remapping on the filtered transcriptome using

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Align reads and estimate abundance** {% icon tool %} with the following parameters:
>    - *"Transcripts"*: `transcriptome_filtered.fasta`
>    - *"Paired or Single-end data?"*: `Paired`
>        - *"Left/Forward strand reads"* -> `Multiple datasets`
>            - Click on the *Folder* button at the right
>                - *Type to Search*: `left`
>                - Select the 6 `Trimmomatic on ..._left.fq.gz`
>        - *"Right/Reverse strand reads"* -> `Multiple datasets`
>            - Click on the *Folder* button at the right
>                - *Type to Search*: `right`
>                - Select the 6 `Trimmomatic on ..._left.fq.gz`
>        - *"Strand specific data"*: `Yes`
>    - *"Abundance estimation method"*: `Salmon`
>    - In *"Additional Options"*:
>        - *"Trinity assembly?"*: `Yes`
> 2. **Rename** the 6 `* isoforms counts` :(
>    - Check in the information panel (**i** icon) the lineage of your file (ex: `A1_left.fq.gz` ... )
>    - Rename the datasets: `A1`, `A2`, `A3`, `B1`, `B2`, `B3`.
>
>    > ### {% icon comment %} Comment
>    >
>    > If you check at the Standard Error messages of your outputs. You can get the `Mapping rate`
>    > 1. Click on one dataset
>    > 2. Click on the little **i** icon
>    > 3. Click on *Tool Standard Error:	stderr*
>    > ```
>    > [2019-11-14 15:44:21.500] [jointLog] [info] Mapping rate = 44.4358%
>    > ```
>    {: .comment}
>
>    > ### {% icon comment %} Comment
>    >
>    > At this stage, you can now delete some useless datasets
>    > - `Align reads and estimate abundance on *: genes counts`
>    > Note that the dataset are just hidden. You can delete them permanently and make some room in the history options (the little wheel icon)
>    {: .comment}
>
{: .hands_on}

## Merge the mapping tables and compute a TMM normalization

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Build expression matrix** {% icon tool %} with the following parameters:
>    - *"Abundance estimates"*: `A1`, `A2`, `A3`, `B1`, `B2`, `B3`
>    - *"Abundance estimation method"*: `Salmon`
> 2. **Describe samples and replicates**  {% icon tool %} with the following parameters:
>    - *"Samples"*
>        - *"1: Samples"*:
>            - *"Full sample name"*: `A1`
>            - *"Condition"*: `A`
>        - *"2: Samples"*:
>            - *"Full sample name"*: `A2`
>            - *"Condition"*: `A`
>        - ...:
>        - *"6: Samples"*:
>            - *"Full sample name"*: `B3`
>            - *"Condition"*: `B`
>
{: .hands_on}

## RNASeq samples quality check

> ### {% icon hands_on %} Hands-on: Task description
> 1. **RNASeq samples quality check** {% icon tool %} with the following parameters:
>    - *"Expression matrix"*: `Build expression matrix: estimated RNA-Seq fragment isoform counts (raw counts)`
>    - *"Samples description"*: `Describe samples`
>
{: .hands_on}

## Differential expression analysis

> ### {% icon hands_on %} Hands-on: Task description
> 1. **Differential expression analysis** {% icon tool %} with the following parameters:
>    - *"Expression matrix"*: `Build expression matrix: estimated RNA-Seq fragment isoform counts (raw counts)`
>    - *"Sample description"*: `Describe samples` (the last one)
>    - *"Differential analysis method"*: `DESeq2`
>
{: .hands_on}

## Extract and cluster differentially expressed transcripts

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Extract and cluster differentially expressed transcripts** {% icon tool %} with the following parameters:
>    - In *"Additional Options"*:
>        - *"Expression matrix"*: `Build expression matrix: estimated RNA-Seq fragment isoform counts (raw counts)`
>        - *"Sample description"*: `Describe samples`
>        - *"Differential expression results"*: `Differential expression results on data XXX and data XXX`
>        - *"p-value cutoff for FDR"*: `1`
>        - *"Run GO enrichment analysis"*: `No`
>
>    > ### {% icon comment %} Comment
>    >
>    > *"p-value cutoff for FDR"*: `1`
>    > Don't do this at home! It's because we have a Toy Dataset. The cutoff should be around `0.001`
>    {: .comment}
>
{: .hands_on}

## Partition genes into expression clusters

> ### {% icon hands_on %} Hands-on: Task description
>
> 1. **Partition genes into expression clusters** {% icon tool %} with the following parameters:
>    - *"RData file"*: `Extract and cluster differentially expressed transcripts: RData file`
>    - *"Method for partitioning genes into clusters"*: `Cut tree based on x percent of max(height) of tree`
>
{: .hands_on}

# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.
