`run.R` showcases how to use the functions to generate desired plots. It takes several sources to summarize expression and mutation profiles of a given set of interested genes. The script also read in UBC1/2 differentially expressed TFs from Ian & Leo to generate the summary. 

Current inputs are: 
- Cavalli 2017 array-based RNA;
- bulk RNA-seq from Liams Hendrikse; 
- PCAWG SNV/INDEL/SVs processed & liftOver from hg19 to hg38 by Xinghan Sun.

Current output plots include:
- RNA-expression by MB groups in array/bulk RNA-seq 
- Mutation (SNV/INDEL/SV) count in the order of 1kb upstream TSS, exon, gene (exon excluded). Each tumour will be only counted once; i.e. if one mutation spans multiple classes (e.g. exon and intron), the plot will only count as 1 exon mutation. NOTE: for transposition SV, I only considered breakpoints.
- Survival using median expression to classify patients.

# Data preprocessing

- bulkRNA Cavalli:
The expression matrix and metadata were downloaded from Cavalli et al 2017 (Affimetrix Human Gene 1.1 ST Array profiling of 763 primary medullobalstoma samples (G3, G4, SHH, WNT) used for identification of Medullobastoma subtypes). Since the matrix was already log2 normalized, it was used to plot as is. 

- bulkRNA Hendrikse:
RNA-seq based bulk expression. 486 G3 & 4 tumours were left after removing 3 samples (2 because of unmatched metadata ID, 1 as it's the only SHH tumour). Count matrix was processed using DESeq2 to generate DESeqDataSet object with normalized counts, which was then used to generate log2(1+normalized counts) for the violin plots. 

- PCAWG SNV/INDEL:
PCAWG .consensus.20160830.somatic.snv_mnv.vcf.gz for MB samples were downloaded and converted to bed format files using a homemade python script (https://github.com/RealPaiLab/MB_noncodingMut/blob/main/Script/snakemake_workflow/workflow/scripts/PCAWG_mutation/merge_mutations.py). The script identifies mutation type (insertion, deletion, snv, long fragment substitution: e.g. AGACAGA -> TTTT) and generate appropriate starts and ends. As PCAWG was in hg19, the converted bed was liftover to hg38 build using cmd LiftOver tool (https://github.com/RealPaiLab/MB_noncodingMut/blob/main/Script/snakemake_workflow/workflow/rules/mutation/generatePcawgMutation.smk). 

- PCAWG SV:
PCAWG .pcawg_consensus_1.6.161116.somatic.sv.bedpe.gz for MB samples were downloaded and converted to bed format files using a homemade python script (https://github.com/RealPaiLab/MB_noncodingMut/blob/main/Script/snakemake_workflow/workflow/scripts/PCAWG_mutation/merge_svs.py). For SV type "DEL", "DUP", "t2tINV", "h2hINV", the start and end of the breakpoints were converted to bed format; for SV type "TRA" (translocation), the two breakpoints were converted as if they were two SNVs. As PCAWG was in hg19, the converted bed was liftover to hg38 build using cmd LiftOver tool (https://github.com/RealPaiLab/MB_noncodingMut/blob/main/Script/snakemake_workflow/workflow/rules/mutation/generatePcawgMutation.smk). 

- Mutations affecting a gene:
Each gene was divided into three regions - 1kb upstream TSS (strand aware), Exon, and Gene (exons excluded) based on gencode.v44.basic.annotation.gtf. In order to correctly show number of tumours mutated in the given gene, mutations overlapping multiple regions were ordered as 1kb Upstream TSS > Exon > Gene (exons excluded); i.e. if mutations are found in all three regions, the tumour will only be counted once and recorded as mutated in 1kb Upstream TSS. 
