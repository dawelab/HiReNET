# HiReNET – Higher-order Repeat Network Exploration Tool

## Introduction:
HiReNET is a graph-based pipeline for detecting and analyzing higher-order repeats (HORs) in genomic sequences. The pipeline has been applied to both Arabidopsis thaliana and maize, providing insights into centromeric and knob repeat structure.
## Installation:
<pre> 
# Clone the Repository
git clone https://github.com/<yourname>/HiReNET.git
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

## Usage:
### Step 1: Prepare profile hidden Markov models (phmm’s) for repeat monomer identification. 
Before starting this pipeline, whole-genome sequencing short reads should be analyzed with RepeatExplorer TAREAN to generate consensus repeats and repeat variants for each repeat type. RepeatExplorer is available through the public Galaxy server at https://repeatexplorer-elixir.cerit-sc.cz/. Using the outputs from RepeatExplorer TAREAN, you can then select the repeat of interest and generate profile hidden Markov models (pHMMs). 

<img width="568" height="325" alt="image" src="https://github.com/user-attachments/assets/7b58d4d6-c40f-4d3e-bfb8-5cbb6f6f23f2" />

<pre>
Usage:
	
HiReNET getphmm -i test_con_variant.fasta -o phmm -p test
</pre>

### Step 2: Identify repeat monomers and arrange monomers in the customized bins.
Repeat monomers can be identified using profile hidden Markov models (pHMMs) and then organized into customized bins. A bin size of 10 kb is typically a good starting point.

<img width="569" height="340" alt="image" src="https://github.com/user-attachments/assets/a6d9bc54-c2e8-417e-8299-d9dc159dc528" />

<pre>
Usage:
	
HiReNET arrayfind -g genome.fasta -c consensus.fasta -o arrayout -p test
HiReNET monomerfind \
	--arrays-dir arrayout \
	--outdir monomerout \
	--prefix test --hmm phmm/test.hmm \
	--chr "chr1,chr2"
HiReNET arrangemonomer \
	--monomer-dir monomerout/test_monomers \
	--outdir arrangemonomer_10kb \
	--prefix test \
	--bin 10000 \
	--chr "chr1,chr2"
</pre>

### Step 3: Classify repeat bins into three classes (Order, HOR, Disorder) using the pre-trained LDA model.
Within each bin, monomers are compared in an all-to-all manner using BLAT. The resulting output is processed to calculate Jaccard index score, which are used to construct a network. Network structure and monomer information are then combined into a feature table that serves as input for a pre-trained LDA model. Each bin is classified with the LDA model, after which adjacent bins sharing the same class and threshold are merged, and their monomers are rearranged into the merged bins.

<img width="568" height="676" alt="image" src="https://github.com/user-attachments/assets/9d65b6a6-a3d8-495f-a5ba-84ede92f0912" />

<pre>
Usage: 
	
HiReNET comparemonomer \
	--bins-dir arrangemonomer_10kb/test_bin_monomers \
	--outdir comparemonomers_wcon \
	--consensus consensus.fasta 
	
HiReNET comparemonomer \
	--bins-dir arrangemonomer_10kb/test_bin_monomers \
	--outdir comparemonomers2_ncon

HiReNET classprediction \
	--blatsub comparemonomers_wcon/blat_output_sub \
	--outdir classpred_out \
	--prefix test \
	--bin 10000
	
HiReNET classprediction \
	--blatsub comparemonomers_wcon/blat_output_sub \
	--outdir classpred_out2 \
	--prefix test \
	--bin 10000 \
	--plot
</pre>


### Step 4: Annotate local HOR patterns through kmer-based analysis for each bin.
For each merged bin, monomers are compared in an all-to-all manner again and a network is constructed again using the optimal similarity threshold. Monomers are then grouped based on network communities, and higher-order repeat (HOR) patterns are identified within each merged bin.

<img width="569" height="690" alt="image" src="https://github.com/user-attachments/assets/7b1f152e-a22e-4c6b-8333-e3f5662e12c9" />

<pre>
Usage:
	
HiReNET rearrangemonomers \
  --bins classpred_out/test_fin_bins_combined.txt \
  --class HOR \
  --prefix test \
  --monomer-dir monomerout/test_monomers \
  --outdir rearrange_monomers_mergebin \
  --chr "chr1,chr2"

HiReNET comparemonomer \
	--bins-dir rearrange_monomers_mergebin/re_arrange_monomers \
	--outdir compare_rearrangemonomers
	
HiReNET networkHOR \
	--blatsub compare_rearrangemonomers/blat_output_sub \
	--bins classpred_out/AthCEN178_fin_bins_combined.txt \
	--coor rearrange_monomers_mergebin/test_monomer_bed_inbin.txt 
	--outdir network_HOR_mergebin

</pre>


### Step 5: Find shared HOR patterns on the chromosome level or the genome level.
In each merged HOR bin, monomers with the same label are extracted to generate consensus HOR monomers. Consensus HOR monomers that share the same threshold are combined, and these are compared in an all-to-all manner across thresholds ranging from 0.90 to 0.99. Monomers are then relabeled, and shared HOR patterns are identified for each threshold.

<img width="567" height="708" alt="image" src="https://github.com/user-attachments/assets/58454a62-43cf-4620-8518-9a69a93510c9" />

<pre>
Usage:
	
HiReNET arrangeHORmonomer \
	--groupdir network_HOR_mergebin \
	--monomer-dir monomerout/test_monomers \
	--outdir network_mergebin_consensus
	
HiReNET consensusHORmonomer \
	--outdir network_mergebin_consensus \
	--threads 4 \
	--chroms "chr1,chr2"
	
HiReNET compareConsensus \
	--chr "chr1" --consensdir network_mergebin_consensus/all_recluster_consensus_monomer \
	--outdir compare_consensusHOR_chr1
	
HiReNET sharedHOR\
	--chr "chr1" --datadir compare_consensusHOR_chr1/blat_sub \
	--outdir shared_out_chr1 \
	--letter network_HOR_mergebin/mergebin_string_outputs  
</pre>








