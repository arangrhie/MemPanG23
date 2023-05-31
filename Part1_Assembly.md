# Part1 Assembly

## Where is the data?
All the input and results are available under [a_thal](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/scratch/a_thal/).  
Feel free to check the results here and compare to yours.

Note: part of the input sequencing data is in your system at /opt/assembly_data/

## Some environment variables to set up for convenience
```
export MERQURY=/opt/merqury
export T2T-Polish=/opt/T2T-Polish
```

We will work in `day3_assembly_evaluation`.
```
mkdir -p ~/day3_assembly_evaluation
cd day3_assembly_evaluation
```

## Understand your genome
We will start working with A. thaliana Col-0 strain. There are 4 other strains available for this training course on [a_thal](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/scratch/a_thal/):
* Col-0
* Cvi-0
* Ler-0
* Tanz
* Col-0_x_Cvi-0: a synthetic diploid Col-0 and Cvi-0 cross hybrid F1, using 25x hifi from each Col-0 and Cvi-0 dataset

All datasets were downloaded from [Wlodzimierz et al., Nature 2023](https://www.nature.com/articles/s41586-023-06062-z) and downsampled to 50x hifi. Illumina data for Col-0 has been obtained from [SRR12136403](https://www.ebi.ac.uk/ena/browser/view/SRR12136403).

### Inpsect your genome with GenomeScope2
Before assembling your genome, it's always better to familiarize with your genome.  
The easy way is to count k-mers (short chunks of sequences of length k), and estimate genome features from this. 
We are using `meryl` for counting k-mers, which we can re-use in `merqury`. 


```
ln -s /opt/assembly_data/Col-0.hifi.q20.50x.fq.gz
hifi=Col-0.hifi.q20.50x

# This may take ~10 minutes, feel free to skip and use the $hifi.meryl
meryl count k=21 output $hifi.meryl threads=24 memory=48g $hifi.fq.gz
```
Now let's' get the histogram
```
meryl histogram $hifi.meryl > $hifi.hist
```

Download the histogram locally, and drop it in [GenomeScope2](http://qb.cshl.edu/genomescope/genomescope2.0/).
This will provide you an estimate of the genome size, repetitivenes, level of heterozygosity, which will be the key aspects that will influence your genome assembly.  
Depending on the level of heterozygosity, we will see more or less 'bubbles' in the final assembly graph.

Here are the output of each a. thaliana strain:
* [Col-0](http://genomescope.org/genomescope2.0/analysis.php?code=j75jEqlj2w1l7MxohBY7)
* [Cvi-0](http://genomescope.org/genomescope2.0/analysis.php?code=PVe8io4UzbY6GU8uuH9Z)
* [Ler-0](http://genomescope.org/genomescope2.0/analysis.php?code=QFnDfLPV1WAJubssM4N2)
* [Tanz](http://genomescope.org/genomescope2.0/analysis.php?code=qXhHVqDmTTkNmde8HtKn)

_Compare each results. What do you see?_

### Can I use the Illumina data from the same species, but different sample?
Sometimes, there are limits on samples, and as an alternate source, sequencing set has to be obtained from a different sample but same species.  
As a sanity check, it's a good practice to do this with a tool like [Mash](https://github.com/VGP/vgp-assembly/tree/master/pipeline/mash) in case you have multiple sequencing runs.  
For this simple course, we will use the k-mers from HiFi and Illumina to see how many of them overlap. The Illumina data was obtained from a different study, using different samples.

```
hifi=Col-0.hifi.q20.50x
illumina=Col-0.illumina.k21
meryl intersect $illumina.meryl $hifi.meryl output Col-0.illumina_and_hifi.50x.meryl
$MERQURY/plot/to_hist_for_plotting.sh Col-0.illumina.k21.meryl Illumina Col-0.illumina_and_hifi.50x.meryl Illumina_and_HiFi > Illumina_and_HiFi.hist
Rscript /opt/merqury/plot/plot_spectra_cn.R -f Illumina_and_HiFi.hist -o Illumina_and_HiFi -t line
```

The `meryl intersect` is intersecting the k-mers found of `$hifi.meryl` with `$illumina.meryl`, keeping the counts of `$illumina.meryl`(the first meryl db provided).  
In the output `Illumina_and_HiFi.png`, the original k-mer counts are shown in black, the overlapping k-mers in red.  

_Would you agree the two datasets are coming from the same individual?  
_At home: Try this on k-mers obtained from hifi reads of Cvi-0, and compare against the Illumina k-mers of Col-0_

## Assemble with Verkko

Verkko understands what has been already done, and tries to run based on what's missing.  
For the sake of time, we will only run the last consensus step, assuming all directories from 0- to 6-layout are present.

```
verkko -d asm \
  --hifi Col-0.hifi.q20.50x.fq.gz \
  --nano Col-0.ont_R10.gt_20kb.fastq.gz \
  --threads 24
```
This will create a directory (or assumes a directory) `asm`, and places all the output and intermediate files in there.

The complete output of Col-0 is available at https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/scratch/a_thal/Col-0/asm/

## Assembly graph
Download the `*.noseq.gfa` file, and open on bandage.  
Note the gfa loaded is a `unitig graph`, not the `consensus graph`.

_Download the `*.noseq.gfa` from other assemblies and compare._

Explore the output files.

The graph contains hifi read coverage by default for each unitig node. Add in the ONT coverage by downloading and loading `assembly.ont-coverage.csv` on bandage.

## Assemble with Verkko in trio mode
This time, let's try to understand how a trio-mode assembly works.

First, we need to generate hapmers. Again, for the sake of time, we will skip this step.  
For those who are interested, this step counts k-mers, using k=30, in homopolymer compressed space.
```
meryl count -c k=30 output Col-0_hifi.50x.meryl Col-0.hifi.q20.50x.fq.gz
meryl count -c k=30 output Cvi-0_hifi.50x.meryl Cvi-0.hifi.q20.50x.fq.gz
meryl count -c k=30 output Child_hifi.50x.meryl Col-0.hifi.q20.25x.fq.gz Cvi-0.hifi.q20.25x.fq.gz
```

Now, we will run `hapmers.sh`, to obtain the `inherited` maternal (Col-0) and paternal (Cvi-0) strain.
_Can you interpret the spectra-cn plot?

Now, we are (pretending) we are running verkko in trio mode.
```
verkko -d trio \
  --hifi Col-0.hifi.q20.50x.fq.gz \
  --nano Col-0.ont_R10.gt_20kb.fastq.gz \
  --Col-0_hifi_50x.k30.hapmer.meryl Cvi-0_hifi_50x.k30.hapmer.meryl trio (note to myself: add this later)
  --threads 24
```
The trio mode output of verkko can be found [here](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/scratch/a_thal/Col_x_Cvi/trio).

Download the Col-0_x_Cvi-0 trio version `noseq.gfa`, and `unitig-popped-unitig-normal-connected-tip.colors.csv` from [6-rukki](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/scratch/a_thal/Col_x_Cvi/trio/6-rukki/).

Load the `*.colors.csv` and see the haplotype identified nodes.
_Compare this before and after applying `rukki`.


