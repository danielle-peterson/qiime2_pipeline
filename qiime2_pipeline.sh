#!/bin/bash
# exit when any command fails
# NOTE: Get conda env file
# Test inputs:
## C0047-7E-1A_S66_L001_R1_001.fastq.gz
## C0047-7E-1A_S66_L001_R2_001.fastq.gz
## C1009-1F-1A_S237_L001_R1_001.fastq.gz
## C1009-1F-1A_S237_L001_R2_001.fastq.gz
# Test Mapping file (1 row per sample w/metadata)
## test_data_mapping.csv

set -e

# TODO: Don't hard code
. /Users/danielle/bin/qiime2_config

echo "$project_name"
mkdir qiime_outputs
cd qiime_outputs

mkdir reads_qza

# ultimately change this to echo "$data"
####### how to fix this ###########
# NOTE: change hard-coded paths to inputs
qiime tools import --type SampleData[PairedEndSequencesWithQuality] \
                   --input-path /Users/danielle/Documents/thesis/16S/qiime2/test_data \
                   --output-path reads_qza/reads.qza \
                   --input-format CasavaOneEightSingleLanePerSampleDirFmt

qiime cutadapt trim-paired --i-demultiplexed-sequences reads_qza/reads.qza \
                    --p-cores 4 \
                    --p-front-f ^GTGYCAGCMGCCGCGGTAA \
                    --p-front-r ^CCGYCAATTYMTTTRAGTTT \
                    --o-trimmed-sequences reads_qza/reads_trimmed.qza

# "Denoising" is the process of correcting or removing noisy reads
# exports denoising_stats.qza, table.qza, representative_sequences.qza
# gives information on how many reads were retained at each step of the DADA2 pipeline
qiime dada2 denoise-paired --i-demultiplexed-seqs reads_qza/reads_trimmed.qza \
                           --p-trunc-len-f 270 \
                           --p-trunc-len-r 210 \
                           --p-n-threads 4 \
                           --output-dir dada2_output

# convert stats QZA to tsv
qiime tools export --input-path dada2_output/denoising_stats.qza \
                   --output-path dada2_output/denoising


qiime feature-table summarize --i-table dada2_output/table.qza --o-visualization dada2_output/table_summary.qzv

qiime feature-table filter-features --i-table dada2_output/table.qza \
                                    --p-min-frequency 20 \
                                    --p-min-samples 1 \
                                    --o-filtered-table dada2_output/table_filt.qza

qiime feature-table filter-seqs --i-data dada2_output/representative_sequences.qza \
                                --i-table dada2_output/table_filt.qza \
                                --o-filtered-data dada2_output/rep_seqs_filt.qza

qiime feature-table summarize --i-table dada2_output/table_filt.qza --o-visualization dada2_output/table_filt_summary.qzv

qiime tools export --input-path dada2_output/table_filt_summary.qzv \
                   --output-path dada2_output/filt_summary

# finding p-max-depth to make rarefaction tree
## Log
head dada2_output/filt_summary/sample-frequency-detail.csv

# TODO: Fix hard-coded pmax - max number of reads to use, take lowest from files?
pmax=34917

########### Build quick phylogeny with FastTree
# make multiple-sequence alignment using MAFFT

mkdir tree_out

qiime alignment mafft --i-sequences dada2_output/rep_seqs_filt.qza \
                      --p-n-threads 4 \
                      --o-alignment tree_out/rep_seqs_filt_aligned.qza

qiime alignment mask --i-alignment tree_out/rep_seqs_filt_aligned.qza \
                      --o-masked-alignment tree_out/rep_seqs_filt_aligned_masked.qza

############# Running FastTree
# Finally FastTree can be run on this masked multiple-sequence alignment:

qiime phylogeny fasttree --i-alignment tree_out/rep_seqs_filt_aligned_masked.qza \
                         --p-n-threads 4 \
                         --o-tree tree_out/rep_seqs_filt_aligned_masked_tree

# Add root to tree
# FastTree returns an unrooted tree.
# One basic way to add a root to a tree is to add it add it at the midpoint
# of the largest tip-to-tip distance in the tree, which is done with this command:

qiime phylogeny midpoint-root --i-tree tree_out/rep_seqs_filt_aligned_masked_tree.qza \
                              --o-rooted-tree tree_out/rep_seqs_filt_aligned_masked_tree_rooted.qz

# Generate rarefaction curves
# A key quality control step is to plot rarefaction curves for all of your samples to determine if you performed sufficient sequencing.
# The below command will generate these plots (X is a placeholder for the maximum depth in your dataset).

# TODO: set p-max-depth to variable

qiime diversity alpha-rarefaction --i-table dada2_output/table_filt.qza \
                                  --p-max-depth 34917 \
                                  --p-steps 20 \
                                  --i-phylogeny tree_out/rep_seqs_filt_aligned_masked_tree_rooted.qz.qza \
                                  --o-visualization rarefaction_curves_eachsample.qzv

# extract p-max-depth by extracting table_filt_summary.qzv

# # TODO: Could decorate with metadata using `-m-metadata-file`
# # rarefaction curve with metadata
# qiime diversity alpha-rarefaction \
#                ....
#                --m-metadata-file $mapping_file \
#                ....

#### Assign taxonomy
# You can assign taxonomy to your ASVs using a Naive-Bayes approach implemented
# in the scikit learn Python library and the SILVA database. Note that we have
# trained classifiers for a few different amplicon regions already (which are
# available in the /home/shared/taxa_classifiers folder), but you will need to
# generate your own if your region of interest isn't there. The classifier
# filename below is for the V6/V8 B969F-BA1406R region. The regions that we have
# trained classifiers for are:

# Using Pretrained classifier
# TODO: Take out hard-coded path
# add back classifier variable
qiime feature-classifier classify-sklearn --i-reads dada2_output/rep_seqs_filt.qza \
                                          --i-classifier /Users/danielle/Documents/thesis/16S/qiime2/classifiers/silva-132-99-nb-classifier_pretrained.qza \
                                          --p-n-jobs -2 \
                                          --output-dir taxa

qiime tools export --input-path taxa/classification.qza --output-path taxa

# A more useful output is the interactive stacked bar-charts of the taxonomic
# abundances across samples, which can be output with this command:
# TODO: Take out hard-coded path
qiime taxa barplot --i-table dada2_output/table_filt.qza \
                   --i-taxonomy taxa/classification.qza \
                   --o-visualization taxa/taxa_barplot.qzv\
                   --m-metadata-file /Users/danielle/Documents/thesis/16S/qiime2/mappingfiles/test_map.tsv

# Calculating diversity metrics and generating ordination plots Common alpha and
# beta-diversity metrics can be calculated with a single command in QIIME2. In
# addition, ordination plots (such as PCoA plots for weighted UniFrac distances)
# will be generated automatically as well. This command will also rarefy all
# samples to the sample sequencing depth before calculating these metrics (X is
# a placeholder for the lowest reasonable sample depth; samples with depth below
# this cut-off will be excluded).

# TODO: Take out hard-coded path
qiime diversity core-metrics-phylogenetic --i-table dada2_output/table_filt.qza \
                                          --i-phylogeny tree_out/rep_seqs_filt_aligned_masked_tree_rooted.qz.qza \
                                          --p-sampling-depth 34917 \
                                          --m-metadata-file /Users/danielle/Documents/thesis/16S/qiime2/mappingfiles/test_map.tsv \
                                          --p-n-jobs 2 \
                                          --output-dir diversity

# Identifying differentially abundant features with ANCOM ANCOM is one method to
# test for difference in the relative abundance of features between sample
# groupings. It is a compositional approach that makes no assumptions about
# feature distributions. However, it requires that all features have non-zero
# abundances so a pseudocount of 1 first needs to be added:

qiime composition add-pseudocount --i-table dada2_output/table_filt.qza \
                                  --o-composition-table dada2_output/table_filt_pseudocount.qza

# qiime composition ancom --i-table dada2_output/table_filt_pseudocount.qza \
#                         --m-metadata-file $mapping_file \
#                         --m-metadata-column 'Method' \
#                         --output-dir ancom_output

# qiime tools export dada2_output/table_filt.qza --output-dir dada2_output_exported
# qiime tools export dada2_output/rep_seqs_filt.qza --output-dir dada2_output_exported
