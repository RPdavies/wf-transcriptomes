# wf-transcriptomes

This repository contains a [nextflow](https://www.nextflow.io/) workflow
for assembly and annotation of transcripts from Oxford Nanopore cDNA or direct RNA reads.





## Introduction

This workflow identifies RNA isoforms using either cDNA or direct RNA (dRNA) 
Oxford Nanopore reads.

### Preprocesing
cDNA reads are initially preprocessed by [pychopper](https://github.com/epi2me-labs/pychopper) 
for the identification of full-length reads, as well as trimming and orientation correction (This step is omitted for 
 direct RNA reads).


### Transcript assembly

#### Reference-aided transcript assembly approach
* Full length reads are mapped to a supplied reference genome using [minimap2](https://github.com/lh3/minimap2)
* Transcripts are assembled by [stringtie](http://ccb.jhu.edu/software/stringtie) 
in long read mode (with or without a guide reference annotation) to generate the GFF annotation.
* The annotation generated by the pipeline is compared to the reference annotation. 
using [gffcompare](http://ccb.jhu.edu/software/stringtie/gffcompare.shtml)

#### de novo-based transcript assembly (experimental!)
* Sequence clusters are generated using [isONclust2](https://github.com/nanoporetech/isONclust2)
  * If a reference genome is supplied, cluster quality metrics are determined by comparing    
  with clusters generated from a minimap2 alignment.
* A consensus sequence for each cluster is generated using [spoa](https://github.com/rvaser/spoa)
* Three rounds of polishing using racon and minimap2 to give a final polished CDS for each gene.
* Full-length reads are then mapped to these polished CDS.
* Transcripts are assembled by stringtie as for the reference-based approach.
* __Note__: This approach is currently not supported with direct RNA reads.

### Fusion gene detection
Fusion gene detection is performed using [JAFFA](https://github.com/Oshlack/JAFFA), with the JAFFAL extension for use 
with ONT long reads. 

### Differential expression analysis

Differential gene expression (DGE) and differential transcript usage (DTU) analyses aim to identify genes and/or transcripts that show statistically altered expression patterns in a studied biological system. The results of the differential analyses are presented in a quantitative format and therefore the degree of change (up or down regulation) between experimental conditions can be calculated for each gene identified.

These differential analyses work by taking a “snapshot” of mRNA abundance and calculating the relative levels of transcripts and isoforms. In this context, expression corresponds to the number of messenger RNAs (mRNA) measured from each gene isoform within the organism / tissue / culture being investigated. In order to determine expression levels across the whole genome, sequence data specifically targeting the mRNA molecules can be generated.

Oxford Nanopore Technologies provides a number of sequencing solutions to allow users to generate the required snapshot of gene expression. This can be achieved by both sequencing the mRNA directly, or via a complementary DNA (cDNA) proxy. In contrast to short read sequencing technologies, entire mRNA transcripts can be captured as single reads. The example data provided with this tutorial is from a study based on the PCR-cDNA kit. This is a robust choice for performing differential transcript usage studies. This kit is suitable for preparation of sequence libraries from low mRNA input quantities. The cDNA population is enriched through PCR with low bias; an important prerequisite for the subsequent statistical analysis.

[Workflow-transcriptomes](https://github.com/epi2me-labs/wf-transcriptomes) includes a subworkflow for DGE and DTU. The first step involves using either a reference alignment or _de novo_ assembly approach to create a set of mRNA sequences per sample. These are merged into a non-redundant transcriptome using [stringtie merge](http://ccb.jhu.edu/software/stringtie). The reads are then aligned to the transcriptome using minimap2 in a splice-aware manner. [Salmon](https://github.com/COMBINE-lab/salmon) is used for transcript quantification, giving per transcript counts and then the following R packages are used for analysis.

### Pre-filtering of quantitative data using DRIMSeq
DRIMSeq (Nowicka and Robinson (2016)) is used to filter the transcript count data from the salmon analysis. The filter step will be used to select for genes and transcripts that satisfy rules for the number of samples in which a gene or transcript must be observed and minimum threshold levels for the number of observed reads. The parameters used for filtering are defined in the config.yaml file. The default parameters defined for this analysis include
* min_samps_gene_expr = 3 - a transcript must be mapped to a gene in at least this minimum number of samples for the gene be included in the analysis
*	min_samps_feature_expr = 1 - a transcript must be mapped to an isoform in at least this this minimum number of samples for the gene isoform to be included in the analysis
*	min_gene_expr = 10 - the minimum number of total mapped sequence reads for a gene to be considered expressed
*	min_feature_expr = 3 - the minimum number of total mapped sequence reads for a gene isoform to be considered

### edgeR based differential expression analysis
+A statistical analysis is first performed using edgeR (Robinson, McCarthy, and Smyth (2010), McCarthy et al. (2012)) to identify the subset of differentially expressed genes. The filtered list of gene counts is used as input. A normalisation factor is calculated for each sequence library (using the default TMM method - please see McCarthy et al. (2012) for further details). The defined experimental design is used to calculate estimates of dispersion for each of the gene features. Statistical tests are calculated using the contrasts defined in the experimental design. The differentially expressed genes are corrected for false discovery (fdr) using the method of Benjamini & Hochberg (Benjamini and Hochberg (1995))

### Differential transcript usage using DEXSeq
Differential transcript usage analysis is performed using the R DEXSeq package (Reyes et al. (2013)). Similar to the edgeR package, DEXSeq estimates the variance between the biological replicates and applies generalised linear models for the statistical testing. The key difference is that the DEXSeq method looks for differences at the exon count level. DEXSeq uses the filtered transcript count data prepared earlier in this analysis. 

### StageR stage-wise analysis of DGE and DTU
The final component of this isoform analysis is a stage-wise statistical test using the R software package `stageR` (Van den Berge and Clement (2018)). stageR uses (1) the raw p-values for DTU from the DEXSeq analysis in the previous section and (2) a false-discovery corrected set of p-values from testing whether individual genes contain at least one exon showing DTU. A hierarchical two-stage statistical testing evaluates the set of genes for DTU.

## Running the workflow
For the differential expression analysis section you should have at least 3 repeats for each sample. 
Your fastq data will need to be organised in to 6 directories that represent 3 repeats for each condition. You may also need to provide a condition sheet. 


## Analysis 
Differential gene expression is sensitive to the input data quantity and quality.  There should be equivalence between samples in the number of sequence reads, mapped reads and quality scores. The sequence and alignment summary plots in the report can be used to assess these metrics. There is also a table that shows the transcript per million(TPM) calculated from the salmon counts. TPM normalizes the data for gene length and then sequencing depth, and makes it easier to compare across samples compared to counts.

### Workflow inputs
- Directory containing cDNA/direct RNA reads. Or a directory containing subdirectories each with reads from different samples
  (in fastq/fastq.gz format)
- Reference genome in fasta format (required for reference-based assembly).
- Optional reference annotation in GFF2/3 format (extensions allowed are .gtf(.gz), .gff(.gz), .gff3(.gz)) (required for differential expression analysis `--de_analysis`). Only annotation files from [Encode](https://www.encodeproject.org), [Ensembl](https://www.ensembl.org/index.html) and [NCBI](https://www.ncbi.nlm.nih.gov/) are supported.
- For fusion detection, JAFFAL reference files (see Quickstart) 




## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop),
[Singularity](https://sylabs.io/singularity/) to provide isolation of
the required software. Each method is automated out-of-the-box provided
either docker or singularity is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit out website](https://labs.epi2me.io/wfindex).


### Workflow options

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-transcriptomes --help
```

to see the options for the workflow.

**Download demonstration data**

A small test dataset is provided for the purposes of testing the workflow software. It consists of reads, reference,
and annotations from human chromosome 20 only.
It can be downloaded using:
```shell
wget -O test_data.tar.gz https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-isoforms/wf-isoforms_test_data.tar.gz 
tar -xzvf  test_data.tar.gz
```

**Example execution of a workflow for reference-based transcript assembly and fusion detection**
```
OUTPUT=~/output;
nextflow run epi2me-labs/wf-transcriptomes \
  --fastq ERR6053095_chr20.fastq \
  --ref_genome chr20/hg38_chr20.fa \
  --ref_annotation chr20/gencode.v22.annotation.chr20.gtf \
  --jaffal_refBase chr20/ \
  --jaffal_genome hg38_chr20 \
  --jaffal_annotation "genCode22" \
  --out_dir outdir -w workspace_dir
```

**Example workflow for denovo transcript assembly**
```
OUTPUT=~/output
nextflow run epi2me-labs/wf-transcriptomes \
  --fastq test_data/fastq \
  --transcriptome_source denovo \
  --out_dir ${OUTPUT} \
  -w ${OUTPUT}/workspace \
  --sample sample_id
```
A full list of options can be seen in nextflow_schema.json. 
Parameters can be specified either in a config like `parameter = value` or on the command line like `--parameter value`.
Below are some commonly used parameters in the format used in config files.

Select how the transcriptome used for analysis should be prepared:

- To create a reference transcriptome using an existing reference genome `--transcriptome_source reference-guided` (default)
- Use a a supplied transcriptome `--transcriptome_source precomputed"`
- Gnerate transcriptome via the denovo pipeline `--transcriptome_source denovo"` 


To run the workflow with direct RNA reads `--direct_rna true` (this just skips the pychopper step).

Pychopper and minimap2 can take options via `--minimap2_opts` and `--pychopper_opts`, for example:

- When using the SIRV synthetic test data  
  - `--minimap2_opts '-uf --splice-flank=no'`
- pychopper needs to know which cDNA synthesis kit used, which can be specified with
  - SQK-PCS109: `--pychopper_opts '-k PCS109'` (default)
  - SQK-PCS110: `--pychopper_opts '-k PCS110'`
  - SQK-PCS111: `--pychopper_opts '-k PCS111'`
- pychopper can use one of two available backends for identifying primers in the raw reads
  - nhmmscan `--pychopper opts '-m phmm'` 
  - edlib `--pychopper opts '-m edlib'`

__Note__: edlib is set by default in the config as it's quite a lot faster. However, it may be less sensitive than nhmmscan. 

### Fusion detection

JAFFAL from the [JAFFA](https://github.com/Oshlack/JAFFA)
package is used to identify potential fusion transcripts.  

In order to use JAFFAL, reference files must first be downloaded.
To use pre-processed hg38 genome and GENCODE v22 annotation files (as used in the JAFFAL paper)
do:
```shell
mkdir jaffal_data_dir
cd jaffal_data_dir/
sh path/to/wf-transcriptomes/subworkflows/JAFFAL/download_jaffal_references.sh
````
Then the path to the directory containing the downloaded reference data must be specified with 
`--jaffal_refBase`.


**Using alternative genome and annotation files**

These should be prepared as described
[here](https://github.com/Oshlack/JAFFA/wiki/FAQandTroubleshooting#how-can-i-generate-the-reference-files-for-a-non-supported-genome).

The resulting JAFFAL reference files will look something like `hg38_genCode22.fa`. The following options enable JAFFAL to find these
files:

`--jaffal_genome reference_genome_name` optional (default: `hg38`)  
`--jaffal_annotation jaffal_annotation_prefix` optional (default: `genCode22`) 


__Note__: JAFFAL is not currently working on Mac M1 (osx-arm64 architecture).

### Differential Expression

Differential Expression requires at least 2 replicates of each sample to compare (but we recommend three). You can see an example sample_sheet.csv below.

**Example workflow for differential expression transcript assembly**

#### Sample sheet condition column
The sample sheet should be a comma separated values file (.csv) and include at least three columns named `barcode`, `alias` and `condition`.
- Each `barcode` should refer to a directory of the same name in the input FASTQ directory (in the example below `barcode01` to `barcode06` reflect the `test_data` directory).
- The `alias` column allows you to rename each barcode to an alias that will be used in the report and other output files.
- The condition column will need to contain one of two keys to indicate the two samples being compared. for example: treated/untreated, sample/control etc.

In the default `sample_sheet.csv` available in the test_data directory we have used the following.

eg. sample_sheet.csv
```
barcode,alias,condition
barcode01,sample01,untreated
barcode02,sample02,untreated
barcode03,sample03,untreated
barcode04,sample04,treated
barcode05,sample05,treated
barcode06,sample06,treated
```

You will also need to provide a reference genome and a reference annotation file.
Here is an example cmd to run the workflow. First you will need to download the data with wget. 
eg.
```
wget -O differential_expression.tar.gz https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-isoforms/differential_expression.tar.gz && tar -xzvf differential_expression.tar.gz
nextflow run epi2me-labs/wf-transcriptomes \
  --fastq  differential_expression/differential_expression_fastq \
  --transcriptome-source precomputed \
  --de_analysis \
  --ref_genome differential_expression/hg38_chr20.fa \
  --ref_annotation differential_expression/gencode.v22.annotation.chr20.gff \
  --direct_rna --minimap2_index_opts '-k 15' \
  --ref_transcriptome differential_expression/ref_transcriptome.fasta \
  --sample_sheet test_data/sample_sheet.csv
```
You can also run the differential expression section of the workflow on its own by providing a reference transcriptome and setting the `--transcriptome_source` to precomputed.
eg.
```
nextflow run epi2me-labs/wf-transcriptomes \
  --fastq  differential_expression/differential_expression_fastq \
  --de_analysis \
  --transcriptome_source precomputed \
  --ref_genome differential_expression/hg38_chr20.fa \
  --ref_annotation differential_expression/gencode.v22.annotation.chr20.gtf \
  --direct_rna --minimap2_index_opts '-k 15' \
  --ref_transcriptome differential_expression/ref_transcriptome.fasta \
  --sample_sheet test_data/sample_sheet.csv
```

## Workflow outputs
* an HTML report document detailing the primary findings of the workflow.
* for each sample:
  * [gffcomapre](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml) output directories
  * read_aln_stats.tsv - alignment summary statistics
  * transcriptome.fas - the assembled transcriptome
  * merged_transcritptome.fas - annotated, assembled transcriptome
  * [jaffal](https://github.com/Oshlack/JAFFA) ooutput directories
  
### Fusion detection outputs
in `${out_dir}/jaffal_output_${sample_id}` you will find:
* jaffa_results.csv - the csv results summary file 
* jaffa_results.fasta - fusion transcritpt sequences

### Differential Expression outputs
* `de_analysis/results_dge.tsv` and `de_analysis/results_dge.pdf`- results of `edgeR` differential gene expression analysis.
* `de_analysis/results_dtu_gene.tsv`, `de_analysis/results_dtu_transcript.tsv` and `de_analysis/results_dtu.pdf` - results of differential transcript usage by `DEXSeq`.
* `de_analysis/results_dtu_stageR.tsv` - results of the `stageR` analysis of the `DEXSeq` output.
* `de_analysis/dtu_plots.pdf` - DTU results plot based on the `stageR` results and filtered counts.
* `de_analysis/all_gene_counts.tsv` - Gene counts generated by the `salmon` tool before filtering.
* `de_analysis/de_transcript_counts.tsv` - Transcript counts generated by the `salmon` tool before filtering. 
* `de_analysis/de_tpm_transcript_counts.tsv` - To facilitate comparisons across samples this file shows transcript per million (TPM) of the raw counts. 
* `de_analysis/all_counts_filtered.tsv` - Transcript counts filtered with input criteria. Used for DE analysis.
* `final_non_redundant_transcriptome.fasta` - Transcripts that were used for differential expression analysis including novel transcripts with the identifiers used for DE analysis. 

### References

* Benjamini, Yoav, and Yosef Hochberg. 1995. “Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing.” Journal of the Royal Statistical Society. Series B (Methodological) 57 (1): 289–300. http://www.jstor.org/stable/2346101.
* McCarthy, Davis J., Chen, Yunshun, Smyth, and Gordon K. 2012. “Differential Expression Analysis of Multifactor Rna-Seq Experiments with Respect to Biological Variation.” Nucleic Acids Research 40 (10): 4288–97.
* Nowicka, Malgorzata, and Mark D. Robinson. 2016. “DRIMSeq: A Dirichlet-Multinomial Framework for Multivariate Count Outcomes in Genomics [Version 2; Referees: 2 Approved].” F1000Research 5 (1356). https://doi.org/10.12688/f1000research.8900.2.
* Patro, Robert, Geet Duggal, Michael I Love, Rafael A Irizarry, and Carl Kingsford. 2017. “Salmon Provides Fast and Bias-Aware Quantification of Transcript Expression.” Nature Methods 14 (March). https://doi.org/10.1038/nmeth.4197.
* Robinson, Mark D, Davis J McCarthy, and Gordon K Smyth. 2010. “EdgeR: A Bioconductor Package for Differential Expression Analysis of Digital Gene Expression Data.” Bioinformatics 26 (1): 139–40.
* Love, Michael I., et al. Swimming Downstream: Statistical Analysis of Differential Transcript Usage Following Salmon Quantification. 7:952, F1000Research, 14 Sept. 2018. f1000research.com, https://f1000research.com/articles/7-952



## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [Singularity](https://sylabs.io/singularity/)
* [racon](https://github.com/isovic/racon)
* [spoa](https://github.com/rvaser/spoa)
* [inONclust](https://github.com/ksahlin/isONclust)
* [isONclust2](https://github.com/nanoporetech/isONclust2)