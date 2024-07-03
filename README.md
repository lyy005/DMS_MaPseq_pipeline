![image](https://github.com/lyy005/DMS_MaPseq_pipeline/assets/5472908/c8329096-942c-4dd3-8ba7-ccf155c7a937)# DMS_MaPseq_pipeline
We used RNA framework (https://github.com/dincarnato/RNAFramework) in our study and followed the pipeline from Marinus and Incarnato. 2021. RNA Framework for Assaying the Structure of RNAs by High-Throughput Sequencing. RNA Bioinformatics. https://link.springer.com/protocol/10.1007/978-1-0716-1307-8_5

## 1. Our workflow
The RNA Framework pipeline can be used for DMS-seq OR DMS-MaPseq data. But the parameters are different. So when running the pipeline, we set the parameters for DMS-MaPseq as suggested. If your data is not DMS-MaPseq, please refer to Marinus and Incarnato. 2021 OR Zubradt et al. 2016. DMS-MaPseq for genome-wide or targeted RNA structure probing in vivo. Nature. 

<img width="649" alt="image" src="https://github.com/lyy005/DMS_MaPseq_pipeline/assets/5472908/df2c42e0-fc73-480f-b615-236567dae590">

## 2. rf-index, build reference indexes
"Bowtie v1 can only perform ungapped read alignment, thus it is only suitable for the analysis of RT stop-based experiments (i.e., DMS-seq [7], Structure-seq [10], SHAPE-seq [11], CIRS-seq [12]). It is rather advisable to use Bowtie v2 for the analysis of mutational profiling (MaP) experiments (i.e., SHAPE-MaP [2], DMS-MaPseq [4, 5]), as a substantial part of the mutational infor- mation of these experiments is recorded within sequencing reads in the form of insertions and deletions." (Marinus and Incarnato. 2021. RNA Bioinformatics)

So we made bowtie2 indexes for the reference genome. RNA Framework (rf) can automatically download and index reference genomes from UCSC database. However, as UCSC database does not have Arabidopsis thaliana TAIR10 genome, we built the indexes manually. We first downloaded the transcript data of Arabidopsis thaliana (Arabidopsis_thaliana.TAIR10.cdna.all.fa from EMBL). 
https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/arabidopsis_thaliana/cdna/

  # ucsc does not have TAIR10, so have to build it manually
  # use transcripts instead of genome sequences
  # ref: https://rnaframework-docs.readthedocs.io/en/latest/rf-index/

  # sort the reference sequences: 
  awk 'BEGIN{RS=">"} NR>1 {gsub("\n", "\t"); print ">"$0}' transcripts.fasta | LC_ALL=C sort -t ' ' -k 2,2 | awk '{sub("\t", "\n"); gsub("\t", ""); print $0}' > reference_sorted.fasta

  # build indexes on the sorted reference
  bowtie2-build reference_sorted.fasta reference_sorted


## Citation:
Li P, Li J, Ma Y, Ma L, Wang X, and Li Y, 2024, RNA Structural and Transcriptional Changes in Genes Associated with Redox Homeostasis are Responsible for DNA Damage. 

