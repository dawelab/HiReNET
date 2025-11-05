
                                            <img width="400" height="600" alt="HiReNET_logo" src="https://github.com/user-attachments/assets/dd31a097-a47f-40e4-a76c-0deecac53ebb" />

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
## Docs
Please see the [Wiki](https://github.com/dawelab/HiReNET/wiki) for detailed documentation.
## Citation
