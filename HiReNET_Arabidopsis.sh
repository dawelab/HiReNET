#!/bin/bash

#======== run local HOR analysis for each chromosome  ======== 

./getphmm.sh -i AthCEN178_consensus_variant.fasta -o AthCEN178_phmm -p AthCEN178

./arrayfind.sh -g Ey15.fasta -c AthCEN178_consensus.fasta -o AthCEN178_arrayout -p AthCEN178

./monomerfind.sh \
  --arrays-dir AthCEN178_arrayout \
  --chrom-dir  AthCEN178_arrayout/split_seq \
  --outdir     AthCEN178_monomerout \
  --prefix     AthCEN178 \
  --hmm        AthCEN178_phmm/AthCEN178.hmm \
  --chr "chr1,chr2,chr3,chr4,chr5"

./arrangemonomer.sh \
  --arrays-dir AthgCEN178_arrayout \
  --genomic-bed-dir AthCEN178_monomerout \
  --monomer-dir AthCEN178_monomerout/AthCEN178_monomers \
  --outdir AthCEN178_arrangemonomer_10kb \
  --prefix AthCEN178 \
  --bin 10000 \
  --chr "chr1,chr2,chr3,chr4,chr5"

./comparemonomer.sh --bins-dir  AthCEN178_arrangemonomer_10kb/AthCEN178_bin_monomers --outdir AthCEN178_comparemonomers --consensus AthCEN178_consensus.fasta

./classprediction.sh --blatsub AthCEN178_comparemonomers/blat_output_sub --outdir AthCEN178_classpred_out --prefix AthCEN178_ --bin 10000 --plot


./rearrangemonomers.sh \
  --bins AthCEN178_classpred_out/AthCEN178_fin_bins_combined.txt \
  --class HOR \
  --prefix AthCEN178 \
  --monomer-dir AthCEN178_monomerout/AthCEN178_monomers \
  --outdir AthCEN178_rearrange_monomers_mergebin \
  --chr "chr1,chr2,chr3,chr4,chr5"

./comparemonomer.sh --bins-dir AthCEN178_rearrange_monomers_mergebin/re_arrange_monomers --outdir AthCEN178_compare_rearrangemonomers


./networkHOR.sh --blatsub AthCEN178_compare_rearrangemonomers/blat_output_sub --bins AthCEN178_classpred_out/AthCEN178_fin_bins_combined.txt --coor AthCEN178_rearrange_monomers_mergebin/AthCEN178_monomer_bed_inbin.txt --outdir AthCEN178_network_HOR_mergebin

./arrangeHORmonomer.sh --groupdir AthCEN178_network_HOR_mergebin --monomer-dir AthCEN178_monomerout/AthCEN178_monomers --outdir AthCEN178_network_mergebin_consensus 
./consensusHORmonomer.sh --outdir AthCEN178_network_mergebin_consensus --threads 10 --chroms "chr1" 
./compareConsensus.sh --chr "chr1" --consensdir AthCEN178_network_mergebin_consensus/all_recluster_consensus_monomer --outdir AthCEN178_compare_consensusHOR_chr1


./sharedHOR_cus.sh --chr "chr1" --datadir AthCEN178_compare_consensusHOR_chr1/blat_sub --outdir AthCEN178_shared_out_chr1 --letter AthCEN178_network_HOR_mergebin/mergebin_string_outputs 

#================================================================
# Run for all chromosomes Arabidopsis
./arrayfind.sh -g Ey15.fasta -c AthCEN178_consensus.fasta -o AthCEN178_arrayout -p AthCEN178

./monomerfind.sh \
  --arrays-dir AthCEN178_arrayout \
  --chrom-dir  AthCEN178_arrayout/split_seq \
  --outdir     AthCEN178_monomerout \
  --prefix     AthCEN178 \
  --hmm        AthCEN178_phmm/*.hmm \
  --chr "chr1,chr2,chr3,chr4,chr5"

./arrangemonomer.sh \
  --arrays-dir AthCEN178_arrayout \
  --genomic-bed-dir AthCEN178_monomerout \
  --monomer-dir AthCEN178_monomerout/AthCEN178_monomers \
  --outdir AthCEN178_arrangemonomer_10kb \
  --prefix AthCEN178 \
  --bin 10000 \
  --chr "chr1,chr2,chr3,chr4,chr5"

./comparemonomer.sh --bins-dir  AthCEN178_arrangemonomer_10kb/AthCEN178_bin_monomers --outdir AthCEN178_comparemonomers --consensus AthCEN178_consensus.fasta

./classprediction.sh --blatsub AthCEN178_comparemonomers/blat_output_sub --outdir AthCEN178_classpred_out --prefix AthCEN178 --bin 10000

./rearrangemonomers.sh \
  --bins AthCEN178_classpred_out/AthCEN178_fin_bins_combined.txt \
  --class HOR \
  --prefix AthCEN178 \
  --monomer-dir AthCEN178_monomerout/AthCEN178_monomers \
  --outdir AthCEN178_rearrange_monomers_mergebin \
  --chr "chr1,chr2,chr3,chr4,chr5"

./comparemonomer.sh --bins-dir AthCEN178_rearrange_monomers_mergebin/re_arrange_monomers --outdir AthCEN178_compare_rearrangemonomers

./networkHOR.sh --blatsub AthCEN178_compare_rearrangemonomers/blat_output_sub --bins AthCEN178_classpred_out/AthCEN178_fin_bins_combined.txt --coor AthCEN178_rearrange_monomers_mergebin/AthCEN178_monomer_bed_inbin.txt --outdir AthCEN178_network_HOR_mergebin

./arrangeHORmonomer.sh --groupdir AthCEN178_network_HOR_mergebin --monomer-dir AthCEN178_monomerout/AthCEN178_monomers --outdir AthCEN178_network_mergebin_consensus 

# single chromosome
#chr1
./consensusHORmonomer.sh --outdir AthCEN178_network_mergebin_consensus --threads 10 --chroms "chr1"
./compareConsensus.sh --chr "chr1" --consensdir AthCEN178_network_mergebin_consensus/all_recluster_consensus_monomer --outdir AthCEN178_compare_consensusHOR_chr1
./sharedHOR_cus.sh --chr "chr1" --datadir AthCEN178_compare_consensusHOR_chr1/blat_sub --outdir AthCEN178_shared_out_chr1 --letter AthCEN178_network_HOR_mergebin/mergebin_string_outputs  

#chr2
./consensusHORmonomer.sh --outdir AthCEN178_network_mergebin_consensus --threads 10 --chroms "chr2"
./compareConsensus.sh --chr "chr2" --consensdir AthCEN178_network_mergebin_consensus/all_recluster_consensus_monomer --outdir AthCEN178_compare_consensusHOR_chr2
./sharedHOR_cus.sh --chr "chr2" --datadir AthCEN178_compare_consensusHOR_chr2/blat_sub --outdir AthCEN178_shared_out_chr2 --letter AthCEN178_network_HOR_mergebin/mergebin_string_outputs  


#chr3
./consensusHORmonomer.sh --outdir AthCEN178_network_mergebin_consensus --threads 10 --chroms "chr3"
./compareConsensus.sh --chr "chr3" --consensdir AthCEN178_network_mergebin_consensus/all_recluster_consensus_monomer --outdir AthCEN178_compare_consensusHOR_chr3
./sharedHOR_cus.sh --chr "chr3" --datadir AthCEN178_compare_consensusHOR_chr3/blat_sub --outdir AthCEN178_shared_out_chr3 --letter AthCEN178_network_HOR_mergebin/mergebin_string_outputs  

#chr4
./consensusHORmonomer.sh --outdir AthCEN178_network_mergebin_consensus --threads 10 --chroms "chr4"
./compareConsensus.sh --chr "chr4" --consensdir AthCEN178_network_mergebin_consensus/all_recluster_consensus_monomer --outdir AthCEN178_compare_consensusHOR_chr4
./sharedHOR_cus.sh --chr "chr4" --datadir AthCEN178_compare_consensusHOR_chr4/blat_sub --outdir AthCEN178_shared_out_chr4 --letter AthCEN178_network_HOR_mergebin/mergebin_string_outputs  


#chr5
./consensusHORmonomer.sh --outdir AthCEN178_network_mergebin_consensus --threads 10 --chroms "chr5"
./compareConsensus.sh --chr "chr5" --consensdir AthCEN178_network_mergebin_consensus/all_recluster_consensus_monomer --outdir AthCEN178_compare_consensusHOR_chr5
./sharedHOR_cus.sh --chr "chr5" --datadir AthCEN178_compare_consensusHOR_chr5/blat_sub --outdir AthCEN178_shared_out_chr5 --letter AthCEN178_network_HOR_mergebin/mergebin_string_outputs  

