#!/bin/bash
# Test inputs found in `/lovelace/echo/sequencing/16S/rawfastq/`:
## test_data/C0047-7E-1A_S66_L001_R1_001.fastq.gz
## test_data/C0047-7E-1A_S66_L001_R2_001.fastq.gz
## test_data/C1009-1F-1A_S237_L001_R1_001.fastq.gz
## test_data/C1009-1F-1A_S237_L001_R2_001.fastq.gz
# Test Mapping file (1 row per sample w/metadata)
## test_data_mapping.csv

# exit when any command fails
set -e

# TODO: Don't hard code
. qiime2_config

echo $Q2_PROJECT
mkdir -p $Q2_OUT/reads_qza

# ultimately change this to echo "$data"
####### how to fix this ###########
# NOTE: change hard-coded paths to inputs

qiime tools import \
    --type SampleData[PairedEndSequencesWithQuality] \
    --input-path   $Q2_IN \
    --output-path  $Q2_OUT/reads_qza/reads.qza \
    --input-format CasavaOneEightSingleLanePerSampleDirFmt

qiime cutadapt trim-paired \
    --i-demultiplexed-sequences $Q2_OUT/reads_qza/reads.qza \
    --p-cores 4 \
    --p-front-f ^GTGYCAGCMGCCGCGGTAA \
    --p-front-r ^CCGYCAATTYMTTTRAGTTT \
    --o-trimmed-sequences $Q2_OUT/reads_qza/reads_trimmed.qza

# "Denoising" is the process of correcting or removing noisy reads
# exports denoising_stats.qza, table.qza, representative_sequences.qza
# gives information on how many reads were retained at each step of the DADA2 pipeline
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs $Q2_OUT/reads_qza/reads_trimmed.qza \
    --p-trunc-len-f 270 \
    --p-trunc-len-r 210 \
    --p-n-threads 4 \
    --output-dir $Q2_OUT/dada2_output

# convert stats QZA to tsv
qiime tools export \
    --input-path  $Q2_OUT/dada2_output/denoising_stats.qza \
    --output-path $Q2_OUT/dada2_output/denoising


qiime feature-table summarize \
    --i-table $Q2_OUT/dada2_output/table.qza \
    --o-visualization $Q2_OUT/dada2_output/table_summary.qzv

qiime feature-table filter-features \
    --i-table $Q2_OUT/dada2_output/table.qza \
    --p-min-frequency 20 \
    --p-min-samples 1 \
    --o-filtered-table $Q2_OUT/dada2_output/table_filt.qza

qiime feature-table filter-seqs \
    --i-data $Q2_OUT/dada2_output/representative_sequences.qza \
    --i-table $Q2_OUT/dada2_output/table_filt.qza \
    --o-filtered-data $Q2_OUT/dada2_output/rep_seqs_filt.qza

qiime feature-table summarize \
    --i-table $Q2_OUT/dada2_output/table_filt.qza \
    --o-visualization $Q2_OUT/dada2_output/table_filt_summary.qzv

qiime tools export \
    --input-path $Q2_OUT/dada2_output/table_filt_summary.qzv \
    --output-path $Q2_OUT/dada2_output/filt_summary

# finding p-max-depth to make rarefaction tree
## Log
head $Q2_OUT/dada2_output/filt_summary/sample-frequency-detail.csv

########## Build quick phylogeny with FastTree
make multiple-sequence alignment using MAFFT

mkdir -p $Q2_OUT/tree_out

qiime alignment mafft \
    --i-sequences $Q2_OUT/dada2_output/rep_seqs_filt.qza \
    --p-n-threads 4 \
    --o-alignment $Q2_OUT/tree_out/rep_seqs_filt_aligned.qza

qiime alignment mask \
    --i-alignment $Q2_OUT/tree_out/rep_seqs_filt_aligned.qza \
    --o-masked-alignment $Q2_OUT/tree_out/rep_seqs_filt_aligned_masked.qza

############# Running FastTree
# Finally FastTree can be run on this masked multiple-sequence alignment:

qiime phylogeny fasttree \
    --i-alignment $Q2_OUT/tree_out/rep_seqs_filt_aligned_masked.qza \
    --p-n-threads 4 \
    --o-tree $Q2_OUT/tree_out/rep_seqs_filt_aligned_masked_tree

# Add root to tree
# FastTree returns an unrooted tree.
# One basic way to add a root to a tree is to add it add it at the midpoint
# of the largest tip-to-tip distance in the tree, which is done with this command:

qiime phylogeny midpoint-root \
    --i-tree $Q2_OUT/tree_out/rep_seqs_filt_aligned_masked_tree.qza \
    --o-rooted-tree $Q2_OUT/tree_out/rep_seqs_filt_aligned_masked_tree_rooted.qz

# Generate rarefaction curves
# A key quality control step is to plot rarefaction curves for all of your samples to determine if you performed sufficient sequencing.
# The below command will generate these plots (X is a placeholder for the maximum depth in your dataset).

# TODO: Fix hard-coded pmax - max number of reads to use, take lowest from files?
PMAX=$(cut -d, -f2 qiime2_output/dada2_output/filt_summary/sample-frequency-detail.csv | sort -n | head -1 | sed -E 's/^([0-9]+)\..+/\1/')
echo $PMAX

qiime diversity alpha-rarefaction \
    --i-table $Q2_OUT/dada2_output/table_filt.qza \
    --p-max-depth $PMAX \
    --p-steps 20 \
    --i-phylogeny $Q2_OUT/tree_out/rep_seqs_filt_aligned_masked_tree_rooted.qz.qza \
    --o-visualization $Q2_OUT/rarefaction_curves_eachsample.qzv

# extract p-max-depth by extracting table_filt_summary.qzv

# # TODO: Could decorate with metadata using `-m-metadata-file`
# # rarefaction curve with metadata
# qiime diversity alpha-rarefaction \
#                ....
#    --m-metadata-file $mapping_file \
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
qiime feature-classifier classify-sklearn \
    --i-reads $Q2_OUT/dada2_output/rep_seqs_filt.qza \
    --i-classifier $Q2_CLASSIFIER \
    --p-n-jobs -2 \
    --output-dir taxa

qiime tools export \
    --input-path $Q2_OUT/taxa/classification.qza --output-path taxa

# A more useful output is the interactive stacked bar-charts of the taxonomic
# abundances across samples, which can be output with this command:
# TODO: Take out hard-coded path
qiime taxa barplot \
    --i-table $Q2_OUT/dada2_output/table_filt.qza \
    --i-taxonomy $Q2_OUT/taxa/classification.qza \
    --o-visualization $Q2_OUT/taxa/taxa_barplot.qzv\
    --m-metadata-file $Q2_MAPPING

# Calculating diversity metrics and generating ordination plots Common alpha and
# beta-diversity metrics can be calculated with a single command in QIIME2. In
# addition, ordination plots (such as PCoA plots for weighted UniFrac distances)
# will be generated automatically as well. This command will also rarefy all
# samples to the sample sequencing depth before calculating these metrics (X is
# a placeholder for the lowest reasonable sample depth; samples with depth below
# this cut-off will be excluded).

# TODO: Take out hard-coded path
qiime diversity core-metrics-phylogenetic \
    --i-table $Q2_OUT/dada2_output/table_filt.qza \
    --i-phylogeny $Q2_OUT/tree_out/rep_seqs_filt_aligned_masked_tree_rooted.qz.qza \
    --p-sampling-depth $PMAX \
    --m-metadata-file $Q2_MAPPING \
    --p-n-jobs 2 \
    --output-dir $Q2_OUT/diversity

# Identifying differentially abundant features with ANCOM ANCOM is one method to
# test for difference in the relative abundance of features between sample
# groupings. It is a compositional approach that makes no assumptions about
# feature distributions. However, it requires that all features have non-zero
# abundances so a pseudocount of 1 first needs to be added:

qiime composition add-pseudocount \
    --i-table $Q2_OUT/dada2_output/table_filt.qza \
    --o-composition-table $Q2_OUT/dada2_output/table_filt_pseudocount.qza

# qiime composition ancom --i-table $Q2_OUT/dada2_output/table_filt_pseudocount.qza \
#    --m-metadata-file $mapping_file \
#    --m-metadata-column 'Method' \
#    --output-dir ancom_output

# qiime tools export $Q2_OUT/dada2_output/table_filt.qza --output-dir $Q2_OUT/dada2_output_exported
# qiime tools export $Q2_OUT/dada2_output/rep_seqs_filt.qza --output-dir $Q2_OUT/dada2_output_exported
