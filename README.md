# MemPanG23
MemPanG23 Assembly &amp; Evaluation

* [Part 1 Assembly](https://github.com/arangrhie/MemPanG23/blob/main/Part1_Assembly.md)
* Part 2 Evaluation

## Tools we will use for this course
1. Understand my genome - [Meryl v1.4 (recommend to install with gcc 10+)](https://github.com/marbl/meryl)
2. Create an assembly - [Verkko v1.3.1 (recommend to install through conda)](https://github.com/marbl/verkko)
3. Check the assembly graph on [BandageNG (local)](https://github.com/asl/BandageNG)
4. Label chromosomes, find telomeres - [MashMap v2.0](https://github.com/marbl/mashmap), [seqtk v1.4](https://github.com/lh3/seqtk)
5. Align long-reads back to the assembly - [Winnowmap2 v2.03](https://github.com/marbl/winnowmap), Samtools v1.17
6. Collect anomalies
   [Merqury (clone latest from github; has dependencies on R v4.3.0 (argparse, ggplot2, and scales), Bedtools v2.30.0, Java Runtime Environment)](https://github.com/marbl/merqury)
   [T2T-Polish (clone latest from github; no additional installation required for this course)](https://github.com/arangrhie/T2T-Polish)
   [seqrequester (clone latest on github and install with gcc 10+)](https://github.com/marbl/seqrequester)
7. Browse through the data with [IGV v2.16.1 (local)](https://software.broadinstitute.org/software/igv/download). I prefer the [IGV and IGV Command Line Version](https://data.broadinstitute.org/igv/projects/downloads/2.16/IGV_2.16.1.zip) to easily manipulate memory limit.
8. Browse through some real issues on HG002 diploid assembly with IGV (local, if time permits)

## Where is the data?
All the input and results are available under [a_thal](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/scratch/a_thal/).  
Feel free to check the results here and compare to yours.

Note: Raw input sequence data (`*.fq.gz` files) is available under /opt/assembly_data/

## Some environment variables to set up for convenience
```
export MERQURY=/opt/merqury
export T2T_Polish=/opt/T2T-Polish
```

We will work in `day3_assembly_evaluation`.
```
mkdir -p ~/day3_assembly_evaluation
cd day3_assembly_evaluation
```
