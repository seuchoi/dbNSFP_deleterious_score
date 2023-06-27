# dbNSFP_deleterious_score
Defining deleterious score using dbNSFP database. This approach was introducted in ```Jurgens, Choi, Morrill et al (2022). Analysis of rare genetic variation underlying cardiometabolic diseases and traits among 200,000 individuals in the UK Biobank. Nat Genet.``` https://www.nature.com/articles/s41588-021-01011-w.

## Software versions
```
R version 4.2.1
R.utils_2.12.2
R.oo_1.25.0
R.methodsS3_1.8.2
tidyr_1.2.1
data.table_1.14.6
```

## Input file
The input file is VEP annotated file using LOFTEE and dbnsfp. Currently, this function support dbsnfp 4.2, 4.3, and 4.4. VEP annotated file also requires gnomAD continental frequency, BIOTYPE (we will select the protein coding transcripts)

### Example to annotate variants using VEP with LOFTEE and dbNSFP
```
vep -i input.vcf \
-o output.vcf.gz \
--format vcf --compress_output gzip --assembly GRCh38 --species homo_sapiens --offline \
--cache --dir_cache /home/jupyter/.vep/ \
--no_stats --minimal --everything --allele_number --show_ref_allele \
--canonical --tab --buffer_size 5000 \
--force_overwrite --symbol \
--dir_plugins /home/jupyter/.vep/Plugins \
--fasta /home/jupyter/.vep/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
--plugin LoF,loftee_path:/home/jupyter/.vep/Plugins,human_ancestor_fa:/home/jupyter/loftee_data/human_ancestor.fa.gz,gerp_bigwig:/home/jupyter/loftee_data/gerp_conservation_scores.homo_sapiens.GRCh38.bw,conservation_file:/home/jupyter/loftee_data/loftee.sql \
--plugin dbNSFP,/home/jupyter/dbNSFP4.2a.gz,SIFT_pred,SIFT4G_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,LRT_pred,MutationTaster_pred,MutationAssessor_pred,FATHMM_pred,PROVEAN_pred,VEST4_rankscore,MetaSVM_pred,MetaLR_pred,M-CAP_pred,REVEL_rankscore,MutPred_rankscore,MVP_rankscore,MPC_rankscore,PrimateAI_pred,DEOGEN2_pred,BayesDel_addAF_pred,BayesDel_noAF_pred,ClinPred_pred,LIST-S2_pred,Aloft_pred,Aloft_Confidence,CADD_phred,DANN_rankscore,fathmm-MKL_coding_pred,fathmm-XF_coding_pred,Eigen-phred_coding,Eigen-PC-phred_coding
```

## Usage
Clone from github and load R packages
```
# clone
git clone https://github.com/seuchoi/dbNSFP_deleterious_score.git
# run R
R
# load packages and repository
library(data.table)
library(tidyr)
library(R.utils)
sourceDirectory("dbNSFP_deleterious_score")
```
Use appropriate inputs. The command line will generate the .RData file for selected types. This .RData file will be a input file for rare variant analysis using GENESIS. mintools = 7 is default. However, this deleterious score will be useful for missense variants. Recommend to set mintools = 0 for hclof and hclof_noflag variants.

Available variant types are:

  - hclof: high-confidence loss-of-function variants
  - hclof_noflag: high-confidence loss-of-function variants without any flag
  - missense: missense variants
  - synonymous: synonymous variants

```
# infile = output from VEP annotation (gziped file)
# types = multiple options are available
# dbnsfp = select your dbnsfp version (only one)
# mintools = minimum number of tools to create deleterious score. If a variant were annotated with less than minimum number of tools, the variant will not be included in the grouping file  
vep_grouping(infile,types=c("hclof","hclof_noflag","missense","synonymous"),dbnsfp=c(4.2,4.3,4.4),mintools=7)
```
