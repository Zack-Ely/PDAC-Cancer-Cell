This repository contains scripts and guides for analyses related to Freed-Pastor, Lambert, & Ely et al (2021).

We have included two directories for scRNA-Seq, one related to analysis of our original mouse data
and one related to our analysis of published human PDAC data from Peng et al (2019). Both directories
contain scripts that can be used to help reproduce our analysis and R workspace for convenient
exploration of our processed data. We have included scripts describing how our PAGODA analysis was
generated. A PAGODA workspace is forthcoming, but note the human and mouse workspaces have these
modules pre-loaded into the relevant Seurat objects. We used Seurat V3.2.2 in our analyses, and
note this version may be required to reproduce the exact same UMAP morphology, clustering annotations,
etc. These can be easily retrieved from our workspaces if one is using a later version of Seurat.

We have also included a directory with example scripts that can be used as references for our neoepitope
prediction workflow. This workflow was not developed to be a bioinformatics package. Some intermediate steps such
as file reformatting, etc were performed on the command line. However, the manuscript methods and these example scripts
should be more than sufficient for anyone with bioinformatics experience to produce their own version of the workflow. 
However, should anyone have trouble with any steps, please feel free to write to us. We will be happy to help. We are 
considering developing a workflow wrapper script and pasting it here as well.

In our neoepitope analysis, two separate cohorts from TCGA and the DFCI PancSeq study were obtained and analyzed. 
The scripts in this repository mostly refer to the TCGA analysis, but the same parameters were applied to 
the PancSeq cohort, with a couple exceptions (see manuscript for details). 
The lists of each cohort's sample IDs used for manuscript figures are also included in this repository.

Also note that some of these commands assume the employed package path has already been implemented in a .bashrc file.

Please refer to the manuscript methods for additional detail and specific package versions.

Paths in all scripts will need to be modified if downloaded. Also, all tools with the versions specified in the manuscript must
be installed.

Also note that some TCGA tumor samples are labeled with the "01B" designation and some normal samples with the (e.g.,) "11A" designation. 
Make sure this is accounted for if analyzing TCGA samples.


General steps for neoepitope prediction workflow:

1: run seq2hla. Same parameters for both TCGA and DFCI.

2: run optitype for all tumor and matched normal samples. Note the example script is for TCGA patients. Refer to manuscript for PancSeq parameter notes.

3: Collate and compare OptiType calls between matched normal and tumor settings. The consistency for allele calls for TCGA
and PancSeq cohorts is reported in our manuscript methods.

4: Collate and compare calls between OptiType and seq2HLA. We collated these into one matrix with the following
format:
	sample_ID,seq2HLA_A_allele_1,A_1_confidence_score,seq2HLA_A_allele_2,A_2_confidence_score,seq2HLA_B_allele_1,B_1_confidence_score,seq2HLA_B_allele_2,B_2_confidence_score,seq2HLA_C_allele_1,C_1_confidence_score,seq2HLA_C_allele_2,C_2_confidence_score,OptiType_A_allele_1,OptiType_A_allele_2,OptiType_B_allele_1,OptiType_B_allele_2,OptiType_C_allele_1,OptiType_C_allele_2
	Example line:
	sample_1,A*68:02,0.0002468358,A*32:01,0.01223943,B*14:02,0.0,B*40:02,0.00430766,C*02:02,0.05207772,C*08:02',0.2209488,A*32:01,A*68:02,B*40:02,B*14:02,C*08:02,C*02:02
	We've included an example python script that was used to evaluate consistency and resolve inconsistencies between seq2HLA and OptiType calls: analyze_hla_typing.py
	Note that readers will need to modify the paths within script if they choose to use it.
	The final HLA matrix used as input into pvacTools is formatted as follows:
	sample_ID,A*31:01,A*26:01,B*35:02,B*38:01,C*04:01,C*04:01.1,E*01:01,E*01:01
		Note that the HLA-E calls only come from seq2HLA and are merged in separately.
 
5: parse the pan-patient MAF file from TCGA into patient-specific maf files. Then convert them into vcf files with the appropriate column header line.
We included a couple of example scripts in this repository.

6: run Manta. An example script is provided. 

7: run Scalpel. Note the input involves a bed file obtained as described in the manuscript. We have included it in this repository
	Note from the manuscript that 13 TCGA samples fail to run through scalpel due to excessive read buildup at certain loci.
	To circumvent this issue, we downsampled the tumor bam files of these samples prior to scalpel, starting with a 50% downsampling
	followed by decrements of 10 until success. We include an example command in this repository.

8: run Strelka. Note the input requires output from Manta. Also note that strelka reports variant allele depth, etc in a 
unique format. We ran 2 subsequent steps after strelka to subset out the PASS variants and then reformat the depth statistics
so that they would be compatible with downstream tools. We arbitrarily use these stats from Strelka as opposed to Scalpel for final variant info.

9: merge variant calls. This starts with taking the intersection of PASSed variants from both Scalpel and Strelka (i.e., only those indels that
were called PASS by both Scalpel and Strelka). This could be done in various ways, but an example bash script is included. After this, the
intersected indel variants are merged with the SNVs. See example script for reference.

10: annotate variant consequence and revise SNV list. Refer to manuscript and documentation for Ensemble VEP and pVACTools for setting
up Ensembl VEP. We have included 2 scripts for reference. The second script revises the list of SNVs but subsetting mostly PASSed variants
but it also retains certain non-PASS designations for some genes that frequently had these designations but also have suspected or confirmed
roles in pancreatic cancer. For example, several KRAS G12 mutations are retained this way. The resulting revised and annotated VCF files
are the ones used as input for neoantigen prediction.

11: Run neoepitope predictions with pVACTools. Filtering parameters are included in the example script, 1_run_neo_pred.sh, but note
that pVACTools outputs all raw epitope scores. We performed the filtering on these results after aggregating the raw results
across all samples and deriving a csv file with all eptiopes predicted with a median affinity <1000 nM. We include an example
custom python script to filter this matrix based on the parameters mentioned in our manuscript. The resulting filtered matrix 
was used to enumerate predicted neoepitopes, etc. 
