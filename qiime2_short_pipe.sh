# QIIME2 pipeline by Peter Leary 

# Prior to running this pipeline, it is assumed that the user has generated FastQC/MultiQC reports and performed appropriate filtering/trimming of adaptors.

# Options used in the following steps relate to MiSeq v2 500 bp paired end 16S rRNA gene sequencing data.



# First, make a directory to save the QZA file into
mkdir q2_output
mkdir q2_output/demux

# Run the QIIME import command 
qiime tools import \
  --input-path Seqs \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path q2_output/demux/demux.qza \
  --type 'SampleData[PairedEndSequencesWithQuality]'

# Generate interactive quality reports
qiime demux summarize \
  --i-data q2_output/demux/demux.qza \
  --o-visualization q2_output/demux/demux.qzv

# Export the interactive quality report 
qiime tools export \
  --input-path q2_output/demux/demux.qzv \
  --output-path html/demux_plot

# Run the main DADA2 pipeline
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs q2_output/demux/demux.qza \
  --output-dir q2_output/dada2 \
  --p-trim-left-f 13 \
  --p-trim-left-r 13 \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 200 \
  --p-n-threads 10
  
# Create a QZV summarising the ASV counts across the samples 
qiime feature-table summarize \
  --i-table q2_output/dada2/table.qza \
  --o-visualization q2_output/dada2/table.qzv \
  --m-sample-metadata-file meta.txt

# Create an interactive ASV FASTA file 
qiime feature-table tabulate-seqs \
  --i-data q2_output/dada2/representative_sequences.qza \
  --o-visualization q2_output/dada2/representative_sequences.qzv
  
# Export the interactive ASV table
qiime tools export --input-path q2_output/dada2/table.qzv \
  --output-path html/asv_table

# Export the interactive ASV representative sequences file
qiime tools export --input-path q2_output/dada2/representative_sequences.qzv \
  --output-path html/rep_seqs

# Classify ASV taxonomy 
qiime feature-classifier classify-sklearn \
  --i-classifier ../EMP_v4/Silva-v138-V4-classifier-2020_2.qza \
  --i-reads q2_output/dada2/representative_sequences.qza \
  --output-dir q2_output/taxonomy \
  --p-n-jobs -11

# Create table of taxonomic classification confidence 
qiime metadata tabulate \
  --m-input-file q2_output/taxonomy/classification.qza \
  --o-visualization q2_output/taxonomy/classification.qzv 

# Export the interactive taxonomic classification file 
qiime tools export \
  --input-path q2_output/taxonomy/classification.qzv \
  --output-path html/classification 

# Generate interactive taxonomic bar charts 
qiime taxa barplot \
  --i-table q2_output/dada2/table.qza \
  --i-taxonomy q2_output/taxonomy/classification.qza \
  --m-metadata-file meta.txt \
  --o-visualization q2_output/taxonomy/taxa-bar-plots.qzv

# Export interactive taxonomic bar charts 
qiime tools export --input-path q2_output/taxonomy/taxa-bar-plots.qzv \
  --output-path html/taxa_bar_plots

# Multiple sequence alignment of the representative ASV sequences 
qiime alignment mafft \
  --i-sequences q2_output/dada2/representative_sequences.qza \
  --output-dir q2_output/aligned \
  --p-n-threads 10

# Mask highly-variable regions in alignments 
qiime alignment mask \
  --i-alignment q2_output/aligned/alignment.qza \
  --o-masked-alignment q2_output/aligned/masked-alignment.qza

# Create a phylogenetic tree of aligned ASVs 
qiime phylogeny fasttree \
  --i-alignment q2_output/aligned/masked-alignment.qza \
  --o-tree q2_output/aligned/unrooted-tree.qza \
  --p-n-threads 10

# Root the tree between the two ASVs with the greatest distance 
qiime phylogeny midpoint-root \
  --i-tree q2_output/aligned/unrooted-tree.qza \
  --o-rooted-tree q2_output/aligned/rooted-tree.qza

# Filter out any ASVs with a total count of < 100 and that are in < 3 samples 
qiime feature-table filter-features \
    --i-table q2_output/dada2/table.qza \
    --p-min-frequency 100 \
    --p-min-samples 3 \
    --o-filtered-table q2_output/dada2/filtered-table-100.qza

# Filter out the ASV rep. seqs. the same as above 
qiime feature-table filter-seqs \
    --i-data q2_output/dada2/representative_sequences.qza \
    --i-table q2_output/dada2/filtered-table-100.qza \
    --o-filtered-data q2_output/dada2/filtered-rep-seqs-100.qza

# Run the filtered ASV table through PICRUSt2 
qiime picrust2 full-pipeline \
   --i-table q2_output/dada2/filtered-table-100.qza \
   --i-seq q2_output/dada2/filtered-rep-seqs-100.qza \
   --output-dir q2_output/picrust2 \
   --p-threads 30 \
   --p-hsp-method pic \
   --p-max-nsti 2

# Create a directory of 'useful' data 
mkdir useful_data 
qiime tools export --input-path q2_output/dada2/table.qza --output-path useful_data/asv_table
qiime tools export --input-path q2_output/dada2/representative-sequences.qza --output-path useful_data/rep-seqs
qiime tools export --input-path q2_output/aligned/rooted-tree.qza --output-path useful_data/rooted-tree
qiime tools export --input-path q2_output/taxonomy/classification.qza --output-path useful_data/taxonomy
qiime tools export --input-path q2_output/picrust2/ko_metagenome.qza --output-path useful_data/picrust2

biom add-metadata -i useful_data/asv_table/feature-table.biom -o useful_data/asv_table/feature-table-tax.biom \
    --observation-metadata-fp useful_data/taxonomy/taxonomy.tsv \
    --observation-header OTUID,taxonomy --sc-separated taxonomy

biom convert -i useful_data/asv_table/feature-table-tax.biom -o useful_data/asv_table/otu_table_w_tax.txt \
    --to-tsv --table-type="OTU table" --header-key taxonomy
    
biom convert -i useful_data/picrust2/feature-table.biom -o useful_data/picrust2/ko_metagenome.txt \
    --to-tsv --table-type="OTU table"
