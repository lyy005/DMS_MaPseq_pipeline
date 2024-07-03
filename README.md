![image](https://github.com/lyy005/DMS_MaPseq_pipeline/assets/5472908/fdb1cce8-9b00-4ddb-be4e-6904ab5666a6)![image](https://github.com/lyy005/DMS_MaPseq_pipeline/assets/5472908/c8329096-942c-4dd3-8ba7-ccf155c7a937)# DMS_MaPseq_pipeline
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

    # build indexes on the sorted reference using Bowtie2 version 2.5.1
    bowtie2-build reference_sorted.fasta reference_sorted

## 3. rf-map, mapping of DMS-MaPseq reads
    # The “-ca3” parameter defines the adapter sequence to clip at the 30 end of reads. With paired-end experiments, a 50 adapter sequence can also be provided via the “-ca5” parameter (this sequence will be automatically reverse-complemented).

    # Therefore, also quality trimming of 50 end can be performed through the “-cq5” parameter.

    # we will use rf-map and Bowtie v2 (enabled by the “-b2” flag):

    # Since we are using Bowtie v2 for read mapping, it is important to pick the right reference index folder (note the “_bt2” suffix).

    # Even though RNA Framework comes with a lot of built-in options, specific mapping parameters can be provided to Bowtie through the “-mp” parameter. In this example, we are directly invoking Bowtie v2 with the “--very-sensitive-local” preset, that causes Bowtie to extensively look for the top-scoring local alignment.

    rf-map --overwrite -p 8 -wt 5 -b2 -cq5 20 -ca3 TCGTATGCCGTCTTCTGCTTG -mp '--very-sensitive-local' -bi reference_sorted -o rf_map_mapseq --keep-logs *.fq.gz

## 4. rf-count, calculates per-base mutation counts (DMS-MaPseq)
    #  -m / --count-mutations flag to make rf-count perform mutation counting

    rf-count --fast --overwrite -p 8 -wt 5 -r --count-mutations --mutation-map -f reference_sorted.fasta -o rf_count_seq ./rf_map_mapseq/*.bam

## 5. rf-norm, Reactivity Normalization of DMS-MaPseq
Raw counts computed by rf-count need to be normalized in order to use them for data-driven RNA folding. This is performed by therf-norm tool in two steps: calculation of raw scores, followed by normalization of base reactivities to values ranging from 0 to 1 (or greater, depending on the normalization method).

When an untreated control is available, this can be used to calculate background mutation frequencies, that will be then sub-tracted from mutation rates in the DMS treated sample. As the Zubradt et al., 2017 scoring scheme does not provide the possibility to account for an untreated control, the Siegfried et al., 2014 scoring scheme [2] (originally introduced for the analysis of SHAPE-MaP data) can be used:

With the Siegfried et al., 2014 scoring method (“-sm 3”), box-plot normalization is recommended (“-nm 3”). The RC file for the untreated control sample is provided through the parameter “-u”.

Optionally, when available, a denatured control sample can also be provided to account for maximum per-base reactivities, through the parameter “-d”. We do not need to worry about this as we do not have a de-natured control sample. 

More about denatured control samples: 
Incarnato et al. 2018. RNA Framework: an all-in-one toolkit for the analysis of RNA structures and post-transcriptional modifications. Nucleic Acids Research. 

Here is the code for rf-norm with a control sample: 

Compared to default settings, we changed -n from 1000 to 5, -ec from 10 to 0 to recover more transcripts

-n: Transcript positions with read coverage behind this threshold will be reported as NaN in the reactivity profile (>0, Default: 10)

-ec (--median-coverage): Discards any transcript with median coverage below this threshold (>=0, Default: 0)

-rb: reactive bases, A and C

    for f in ./rf_count_seq/DMS*.rc 
    do

    rf-norm -rb AC -sm 3 -nm 3 -ec 0 -n 5 -i ./rf_count_seq/index.rci -t $f -u ./rf_count_seq/CK.rc --processors 10 --overwrite -o ./rf_norm/

    done

## 6. rf-correlate, check for correlation between samples
For experiments containing multiple biological replicates, tran- script level (and experiment level) pairwise Pearson correlations can be assessed with the rf-correlate tool, by
    # check correlation
    rf-correlate -m 0.1 ./DMS1_vs_CK_norm/ DMS2_vs_CK_norm/ --overwrite -o DMS1_vs_DMS2

## Citation:
Li P, Li J, Ma Y, Ma L, Wang X, and Li Y, 2024, RNA Structural and Transcriptional Changes in Genes Associated with Redox Homeostasis are Responsible for DNA Damage. 

