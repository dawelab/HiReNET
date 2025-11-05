# HiReNET – Higher-order Repeat Network Exploration Tool

## Introduction:
HiReNET is a graph-based pipeline for detecting and analyzing higher-order repeats (HORs) in genomic sequences. The pipeline has been applied to both Arabidopsis thaliana and maize, providing insights into centromeric and knob repeat structure.
## Installation:
<pre> 
# Clone the Repository
git clone https://github.com/dawelab/HiReNET.git
cd HiReNET

# Build and install
make install

# Add to PATH
echo 'export PATH="$HOME/HiReNET/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

# Verify
HiReNET --help
</pre>
## Dependencies:
HiReNET relies on common bioinformatics tools. Please make sure these are installed and available in your PATH:
- HMMER (3.4)
- BLAST+ (2.16.0)
- BLAT (3.7)
- SeqKit (2.9.0)
- BEDTools (2.31.1)
- bioawk (1.0)
- MAFFT (7.526)
- Clustal-Omega (1.2.4)
- EMBOSS (6.6.0)

R with packages:
- dplyr
- tidyr
- stringr
- purrr
- rlang
- data.table
- ggplot2
- cowplot
- igraph
- caret
- readr

## Quick Start:
<pre>
# Run for mutiple chromosomes or sequences
HiReNET runALL \
  --chr "chr1,chr2,chr3,chr4,chr5" \
  --prefix AthCEN178 \
  --consensus AthCEN178_consensus.fasta \
  --seq Ey15.fasta \
  --variant AthCEN178_consensus_variant.fasta
	
# Run for a single chromosome or sequence
HiReNET runALL \
  --chr "chr1" \
  --prefix AthCEN178 \
  --consensus AthCEN178_consensus.fasta \
  --seq Ey15.fasta \
  --variant AthCEN178_consensus_variant.fasta
	
</pre>
## Usage
see [wiki](https://github.com/dawelab/HiReNET/wiki)
### Part 1: Prepare profile hidden Markov models (phmm’s) for repeat monomer identification. 
Before starting this pipeline, whole-genome sequencing short reads should be analyzed with RepeatExplorer TAREAN to generate consensus repeats and repeat variants for each repeat type. RepeatExplorer is available through the public Galaxy server at https://repeatexplorer-elixir.cerit-sc.cz/. Using the outputs from RepeatExplorer TAREAN, you can then select the repeat of interest and generate profile hidden Markov models (pHMMs). 

<img width="568" height="325" alt="image" src="https://github.com/user-attachments/assets/7b58d4d6-c40f-4d3e-bfb8-5cbb6f6f23f2" />

<pre>
Usage:
# Step 1 — Build profile HMM 
HiReNET getphmm -i AthCEN178_consensus_variant.fasta -o phmm -p AthCEN178
</pre>

### Part 2: Identify repeat monomers and arrange monomers in the customized bins.
Repeat monomers can be identified using profile hidden Markov models (pHMMs) and then organized into customized bins. A bin size of 10 kb is typically a good starting point.

<img width="569" height="340" alt="image" src="https://github.com/user-attachments/assets/a6d9bc54-c2e8-417e-8299-d9dc159dc528" />

<pre>
Usage:
	
# Step 2 — Find repeat arrays
HiReNET arrayfind -g Ey15.fasta -c AthCEN178_consensus.fasta -o AthCEN178_arrayout -p AthCEN178
# Step 3 — Detect monomers
HiReNET monomerfind \
  --arrays-dir AthCEN178_arrayout \
  --chrom-dir  AthCEN178_arrayout/split_seq \
  --outdir     AthCEN178_monomerout \
  --prefix     AthCEN178 \
  --hmm        AthCEN178_phmm/AthCEN178.hmm \
  --chr "chr1,chr2,chr3,chr4,chr5"
# Step 4 — Arrange monomers into bins
HiReNET arrangemonomer \
  --arrays-dir AthCEN178_arrayout \
  --genomic-bed-dir AthCEN178_monomerout \
  --monomer-dir AthCEN178_monomerout/AthCEN178_monomers \
  --outdir AthCEN178_arrangemonomer_10kb \
  --prefix AthCEN178 \
  --bin 10000 \
  --chr "chr1,chr2,chr3,chr4,chr5"
</pre>

### Part 3: Classify repeat bins into three classes (Order, HOR, Disorder) using the pre-trained LDA model.
Within each bin, monomers are compared in an all-to-all manner using BLAT. The resulting output is processed to calculate Jaccard index score, which are used to construct a network. Network structure and monomer information are then combined into a feature table that serves as input for a pre-trained LDA model. Each bin is classified with the LDA model, after which adjacent bins sharing the same class and threshold are merged, and their monomers are rearranged into the merged bins.

<img width="568" height="676" alt="image" src="https://github.com/user-attachments/assets/9d65b6a6-a3d8-495f-a5ba-84ede92f0912" />

<pre>
Usage: 
	
# Step 5 — Compare monomers (self + consensus)
HiReNET comparemonomer \
  --bins-dir AthCEN178_arrangemonomer_10kb/AthCEN178_bin_monomers \
  --outdir AthCEN178_comparemonomers \
  --consensus AthCEN178_consensus.fasta

# Step 6 — Predict HOR classes
# R packages should be installed.
HiReNET classprediction \
  --blatsub AthCEN178_comparemonomers/blat_output_sub \
  --outdir AthCEN178_classpred_out \
  --prefix AthCEN178 \
  --bin 10000 \
  --plot

HiReNET classprediction \
  --blatsub AthCEN178_comparemonomers/blat_output_sub \
  --outdir AthCEN178_classpred_out2 \
  --prefix AthCEN178 \
  --bin 10000 
</pre>


### Part 4: Annotate local HOR patterns through kmer-based analysis for each bin.
For each merged bin, monomers are compared in an all-to-all manner again and a network is constructed again using the optimal similarity threshold. Monomers are then grouped based on network communities, and higher-order repeat (HOR) patterns are identified within each merged bin.

<img width="569" height="690" alt="image" src="https://github.com/user-attachments/assets/7b1f152e-a22e-4c6b-8333-e3f5662e12c9" />

<pre>
Usage:
	
# Step 7 — Rearrange monomers by class
HiReNET rearrangemonomers \
  --bins AthCEN178_classpred_out/AthCEN178_fin_bins_combined.txt \
  --class HOR \
  --prefix AthCEN178 \
  --monomer-dir AthCEN178_monomerout/AthCEN178_monomers \
  --outdir AthCEN178_rearrange_monomers_mergebin \
  --chr "chr1,chr2,chr3,chr4,chr5"
# Step 8 — Compare rearranged monomers
HiReNET comparemonomer \
  --bins-dir AthCEN178_rearrange_monomers_mergebin/re_arrange_monomers \
  --outdir AthCEN178_compare_rearrangemonomers
# Step 9 — Build HOR network for meregd HOR bins
HiReNET networkHOR \
  --blatsub AthCEN178_compare_rearrangemonomers/blat_output_sub \
  --bins AthCEN178_classpred_out/AthCEN178_fin_bins_combined.txt \
  --coor AthCEN178_rearrange_monomers_mergebin/AthCEN178_monomer_bed_inbin.txt \
  --outdir AthCEN178_network_HOR_mergebin

</pre>


### Part 5: Find shared HOR patterns on the chromosome level or the genome level.
In each merged HOR bin, monomers with the same label are extracted to generate consensus HOR monomers. Consensus HOR monomers that share the same threshold are combined, and these are compared in an all-to-all manner across thresholds ranging from 0.90 to 0.99. Monomers are then relabeled, and shared HOR patterns are identified for each threshold.

<img width="567" height="708" alt="image" src="https://github.com/user-attachments/assets/58454a62-43cf-4620-8518-9a69a93510c9" />

<pre>
Usage:
	
# Step 10 — Arrange HOR monomers for consensus 
HiReNET arrangeHORmonomer \
  --groupdir AthCEN178_network_HOR_mergebin \
  --monomer-dir AthCEN178_monomerout/AthCEN178_monomers \
  --outdir AthCEN178_network_mergebin_consensus

# Step 11 — Build consensus HORs per chromosome
HiReNET consensusHORmonomer \
  --outdir AthCEN178_network_mergebin_consensus \
  --threads 10 \
  --chroms "chr1"

HiReNET compareConsensus \
  --chr "chr1" \
  --consensdir AthCEN178_network_mergebin_consensus/all_recluster_consensus_monomer \
  --outdir AthCEN178_compare_consensusHOR_chr1

HiReNET sharedHOR \
  --chr "chr1" \
  --datadir AthCEN178_compare_consensusHOR_chr1/blat_sub \
  --outdir AthCEN178_shared_out_chr1 \
  --letter AthCEN178_network_HOR_mergebin/mergebin_string_output
  --plotv V2
Step 11 — Use loop to build consensus HORs per chromosome 
for chr in chr2 chr3 chr4 chr5; do
  HiReNET consensusHORmonomer \
    --outdir AthCEN178_network_mergebin_consensus \
    --threads 10 \
    --chroms "$chr"

  HiReNET compareConsensus \
    --chr "$chr" \
    --consensdir AthCEN178_network_mergebin_consensus/all_recluster_consensus_monomer \
    --outdir AthCEN178_compare_consensusHOR_${chr}

  HiReNET sharedHOR \
    --chr "$chr" \
    --datadir AthCEN178_compare_consensusHOR_${chr}/blat_sub \
    --outdir AthCEN178_shared_out_${chr} \
    --letter AthCEN178_network_HOR_mergebin/mergebin_string_output
    --plotv V2
done

# Step 12 (optional) — Multi-chromosome shared HOR 
HiReNET sharedHOR_more \
  --chr "chr1,chr2,chr3,chr4,chr5" \
  --datadir AthCEN178_compare_consensusHOR_chr1/blat_sub \
  --outdir AthCEN178_shared_out_all \
  --letter AthCEN178_network_HOR_mergebin/mergebin_string_outputs
</pre>








