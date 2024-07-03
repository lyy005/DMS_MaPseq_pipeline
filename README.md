# DMS_MaPseq_pipeline
We used RNA framework version 2.8.6 (https://github.com/dincarnato/RNAFramework) in our study and followed the pipeline from Marinus and Incarnato. 2021. RNA Framework for Assaying the Structure of RNAs by High-Throughput Sequencing. RNA Bioinformatics. https://link.springer.com/protocol/10.1007/978-1-0716-1307-8_5

## 1. Our workflow
The RNA Framework pipeline can be used for DMS-seq OR DMS-MaPseq data. But the parameters are different. So when running the pipeline, we set the parameters for DMS-MaPseq as suggested. If your data is not DMS-MaPseq, please refer to Marinus and Incarnato. 2021 OR Zubradt et al. 2016. DMS-MaPseq for genome-wide or targeted RNA structure probing in vivo. Nature. 

<img width="649" alt="image" src="https://github.com/lyy005/DMS_MaPseq_pipeline/assets/5472908/df2c42e0-fc73-480f-b615-236567dae590">

## 2. Our samples: 
MMS is for methyl methanesulfonate

DMS is for dimethyl sulfate

We have 8 samples in total:

3 DMS treated samples: DMS1, DMS2, DMS3

1 control sample for DMS treated sample: CKH2O

3 DMS and MMS treated samples: MD1, MD2, MD3

1 control sample for DMS+MMS treated sample: MMSCK


## 3. rf-index, build reference indexes
"Bowtie v1 can only perform ungapped read alignment, thus it is only suitable for the analysis of RT stop-based experiments (i.e., DMS-seq [7], Structure-seq [10], SHAPE-seq [11], CIRS-seq [12]). It is rather advisable to use Bowtie v2 for the analysis of mutational profiling (MaP) experiments (i.e., SHAPE-MaP [2], DMS-MaPseq [4, 5]), as a substantial part of the mutational information of these experiments is recorded within sequencing reads in the form of insertions and deletions." (Marinus and Incarnato. 2021. RNA Bioinformatics)

So we made bowtie2 indexes for the reference genome. RNA Framework (rf) can automatically download and index reference genomes from UCSC database. However, as UCSC database does not have Arabidopsis thaliana TAIR10 genome, we built the indexes manually. We first downloaded the transcript data of Arabidopsis thaliana (Arabidopsis_thaliana.TAIR10.cdna.all.fa from EMBL). 
https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/arabidopsis_thaliana/cdna/

    # ucsc does not have TAIR10, so have to build it manually
    # use transcripts instead of genome sequences
    # ref: https://rnaframework-docs.readthedocs.io/en/latest/rf-index/

    # sort the reference sequences: 
    awk 'BEGIN{RS=">"} NR>1 {gsub("\n", "\t"); print ">"$0}' transcripts.fasta | LC_ALL=C sort -t ' ' -k 2,2 | awk '{sub("\t", "\n"); gsub("\t", ""); print $0}' > reference_sorted.fasta

    # build indexes on the sorted reference using Bowtie2 version 2.5.1
    bowtie2-build reference_sorted.fasta reference_sorted

## 4. rf-map, mapping of DMS-MaPseq reads
    # The “-ca3” parameter defines the adapter sequence to clip at the 30 end of reads. With paired-end experiments, a 50 adapter sequence can also be provided via the “-ca5” parameter (this sequence will be automatically reverse-complemented).

    # Therefore, quality trimming of 50 end can be performed through the “-cq5” parameter.

    # we will use rf-map and Bowtie v2 (enabled by the “-b2” flag):

    # Since we are using Bowtie v2 for read mapping, it is important to pick the right reference index folder (note the “_bt2” suffix).

    # Even though RNA Framework comes with a lot of built-in options, specific mapping parameters can be provided to Bowtie through the “-mp” parameter. In this example, we are directly invoking Bowtie v2 with the “--very-sensitive-local” preset, that causes Bowtie to extensively look for the top-scoring local alignment.

    rf-map --overwrite -p 8 -wt 5 -b2 -cq5 20 -ca3 TCGTATGCCGTCTTCTGCTTG -mp '--very-sensitive-local' -bi reference_sorted -o rf_map_mapseq --keep-logs *.fq.gz
    for i in 
    do
     # map paired-end reads, add read groups
     bowtie2 -x reference_sorted -1 $i\_R1.fastq.gz -2 $i\_R2.fastq.gz -S ./step1_bowtie2_unique/$i\.sam --threads 10 --very-sensitive-local --rg-id $i\ --rg SM:$i
     
     # keep uniquely mapped reads. Ref: https://www.biostars.org/p/101533/
     grep -E "@|NM:" ./step1_bowtie2_unique/$i\.sam | grep -v "XS:" > ./step1_bowtie2_unique/$i\.unique.sam
     samtools view -h -b -S ./step1_bowtie2_unique/$i\.unique.sam -o ./step1_bowtie2_unique/$i\.unique.bam

     # sort and make BAM files
     samtools sort ./step1_bowtie2_unique/${Array[$SLURM_ARRAY_TASK_ID]}.unique.bam -o ./step1_bowtie2_unique/${Array[$SLURM_ARRAY_TASK_ID]}.unique.sorted.bam --threads 10
     samtools index ./step1_bowtie2_unique/${Array[$SLURM_ARRAY_TASK_ID]}.unique.sorted.bam
     
    done

## 5. rf-count, calculates per-base mutation counts (DMS-MaPseq)
    #  -m / --count-mutations flag to make rf-count perform mutation counting

    rf-count --fast --overwrite -p 8 -wt 5 -r --count-mutations -f reference_sorted.fasta -o step2_rf_count_noMM ./step1_bowtie2/*.sorted.bam

## 6. rf-norm, Reactivity Normalization of DMS-MaPseq
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

    # for DMS samples
    for f in ./step2_rf_count_noMM/DMS*.rc 
    do
     rf-norm -rb AC -sm 3 -nm 3 -ec 0 -n 0 -i ./step2_rf_count_noMM/index.rci -t $f -u ./step2_rf_count_noMM/CKH2O.sorted.rc --processors 10 --overwrite -o ./step3_rf_norm_DMS_noMM/
    done

    # for MD samples
    for f in ./step2_rf_count_noMM/MD*.rc
    do
     rf-norm -rb AC -sm 3 -nm 3 -ec 0 -n 0 -i ./step2_rf_count_noMM/index.rci -t $f -u ./step2_rf_count_noMM/MMSCK.sorted.rc --processors 10 --overwrite -o ./step3_rf_norm_MD_noMM/
    done

## 7. rf-correlate, check for correlation between samples
For experiments containing multiple biological replicates, transcript level (and experiment level) pairwise Pearson correlations can be assessed with the rf-correlate tool, by
    
    # check correlation
    # DMS samples
    rf-correlate -m 0.1 ./DMS1_vs_CK_norm/ DMS2_vs_CK_norm/ --overwrite -o DMS1_vs_DMS2
    rf-correlate -m 0.1 ./DMS1_vs_CK_norm/ DMS3_vs_CK_norm/ --overwrite -o DMS1_vs_DMS3
    rf-correlate -m 0.1 ./DMS2_vs_CK_norm/ DMS3_vs_CK_norm/ --overwrite -o DMS2_vs_DMS3

    # MD samples
    rf-correlate -m 0.1 MD1.sorted_vs_MMSCK.sorted_norm/ MD2.sorted_vs_MMSCK.sorted_norm/ --overwrite -o MD1_vs_MD2
    rf-correlate -m 0.1 MD1.sorted_vs_MMSCK.sorted_norm/ MD3.sorted_vs_MMSCK.sorted_norm/ --overwrite -o MD1_vs_MD3
    rf-correlate -m 0.1 MD2.sorted_vs_MMSCK.sorted_norm/ MD3.sorted_vs_MMSCK.sorted_norm/ --overwrite -o MD2_vs_MD3

## 8. rf-combine, combined biological replicates
    
    rf-combine -m 0.1 -o DMS_MaPseq_merge DMS1.sorted_vs_CKH2O.sorted_norm/ DMS2.sorted_vs_CKH2O.sorted_norm/ DMS3.sorted_vs_CKH2O.sorted_norm/
    rf-combine -m 0.1 -o MD_MaPseq_merge MD1.sorted_vs_MMSCK.sorted_norm/ MD2.sorted_vs_MMSCK.sorted_norm/ MD3.sorted_vs_MMSCK.sorted_norm/
    
## 9. rf-fold, fold RNAs

    rf-fold -sl 2.4 -in -0.2 -md 600 -nlp -dp -sh -g DMS_MaPseq_merge/
    rf-fold -sl 2.4 -in -0.2 -md 600 -nlp -dp -sh -g MD_MaPseq_merge/
    
## Citation:
Li P, Li J, Ma Y, Ma L, Wang X, and Li Y, 2024, RNA Structural and Transcriptional Changes in Genes Associated with Redox Homeostasis are Responsible for DNA Damage. 

