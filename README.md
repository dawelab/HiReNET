<p align="center">
<img width="1989" height="1824" alt="logo_v2" src="https://github.com/user-attachments/assets/e736dbe9-d1b7-454f-871b-13f51ec4131c" />

</p>

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
<table>
<tr><th align="center">Bioinformatics Tools</th><th align="center">R Packages</th></tr>
<tr><td align="center">

- HMMER (3.4)  
- BLAST+ (2.16.0)  
- BLAT (3.7)  
- SeqKit (2.9.0)  
- BEDTools (2.31.1)  
- bioawk (1.0)  
- MAFFT (7.526)  
- Clustal-Omega (1.2.4)  
- EMBOSS (6.6.0)  

</td><td align="center">

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

</td></tr>
</table>

## Input Data:
Input sequences must be provided in FASTA format, with each sequence written on a single line. Chromosome names should use chr1, chr2, chr3, rather than chr01, chr02, chr03. If the provided sequences are not derived directly from a genome assembly, name them sequentially as seq1, seq2, seq3.
<pre>
	>chr1
	CCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCT
	>chr2
	AGTAATTCTAGAGCTAATACGTGCAACAAACCCCGACTTATGGAAGGGACGCATTTATTA
	...
	>chr10
	CCCTAAACCCTAAACCCTAAACCTAAACCCTAACTAAACCCTAAACCCTAAACCCTAAAC
	>chr11
	GGATCCGTGGTTTCGCGTATCGGCATGGTCGGGAGCTTTTATCTCGGTCTTGTCGTGCGC
	...
	or 
	>seq1
	CCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCT
	>seq2
	AGTAATTCTAGAGCTAATACGTGCAACAAACCCCGACTTATGGAAGGGACGCATTTATTA
	...
	>seq10
	CCCTAAACCCTAAACCCTAAACCTAAACCCTAACTAAACCCTAAACCCTAAACCCTAAAC
	>seq11
	GGATCCGTGGTTTCGCGTATCGGCATGGTCGGGAGCTTTTATCTCGGTCTTGTCGTGCGC
	...
</pre>
## Quick Start:
<pre>
# Run for a single chromosome or sequence. 
gunzip -k data/test.fasta.gz     # The test.fasta.gz only contains chr1 sequences.
HiReNET runALL \
  --chr chr1 \
  --prefix test \
  --consensus data/AthCEN178_consensus.fasta \
  --seq data/test.fasta \
  --variant data/AthCEN178_consensus_variant.fasta
	
# Run for mutiple chromosomes or sequences. Ey15.fasta is Arabidopsis genome sequences from https://doi.org/10.1038/s41586-023-06062-z.
HiReNET runALL \
  --chr chr1,chr2,chr3,chr4,chr5 \
  --prefix AthCEN178 \
  --consensus data/AthCEN178_consensus.fasta \
  --seq Ey15.fasta \
  --variant data/AthCEN178_consensus_variant.fasta	
</pre>
## Docs
Please see the [Wiki](https://github.com/dawelab/HiReNET/wiki/HiReNET-wiki) for detailed documentation.
## Citation
