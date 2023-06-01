# Part1 Assembly

Learning objectives
* Perform a light qc of the data before beginning a genome assembly project
* Run a haplotype-aware assembly using trio data
* Interprete a genome assembly graph
* Manually fix a path and re-build the consensus

## Where is the data?
All the input and results are available under [a_thal](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/scratch/a_thal/).  
Feel free to check the results here and compare to yours.

Note: Raw input sequence data (`*.fq.gz` files) is available under /opt/assembly_data/

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

All datasets were downloaded from [Wlodzimierz et al., Nature 2023](https://www.nature.com/articles/s41586-023-06062-z) and downsampled to 50x hifi. ONT data for Cvi-0 was about 20x. Ler-0 and Tanz ONT data was downsampled to ~50x. All the ONT read length N50 are about ~30 kb . Illumina data for Col-0 has been obtained from [SRR12136403](https://www.ebi.ac.uk/ena/browser/view/SRR12136403).

Let's begin with Col-0.

```
cd Col-0

# Symlink the hifi and ont data
ln -s /opt/assembly_data/Col-0.hifi_50x.fq.gz
ln -s /opt/assembly_data/Col-0.ont_R10.gt_20kb.fastq.gz
```
_What's the read length N50 and coverage for the ONT data?_

There are multiple ways to obtain this. Here, we will use `[seqrequester](https://github.com/marbl/seqrequester)` for convenience.  
This can be run in the background. Check back later when it's done.
```
seqrequester summarize /opt/assembly_data/Col-0.ont_R10.gt_20kb.fastq.gz > Col-0.ont_R10.gt_20kb.summary &
```


### Inpsect your genome with GenomeScope2
Before assembling your genome, it's always better to familiarize with your genome.  
The easy way is to count k-mers (short chunks of sequences of length k), and estimate genome features from this. 
We are using `[meryl](https://github.com/marbl/meryl)` for counting k-mers, which we can re-use later for evaluating the assembly with `[merqury](https://github.com/marbl/merqury)`. 

```
cd meryl

ln -s /opt/assembly_data/Col-0.hifi.q20.50x.fq.gz
hifi=Col-0.hifi.q20.50x

# This may take ~10 minutes, feel free to skip and use the $hifi.meryl
meryl count k=21 output $hifi.meryl threads=24 memory=48g $hifi.fq.gz
```

Now let's get the histogram.
```
meryl histogram Col-0.hifi_50x.k21.meryl > Col-0.hifi_50x.k21.hist
```

Download the histogram locally, and drop it in [GenomeScope2](http://qb.cshl.edu/genomescope/genomescope2.0/).
This will provide an estimate of the genome size, repetitivenes, level of heterozygosity, which are the key aspects that will influence your genome assembly.  
For example, depending on the level of heterozygosity, we will see more or less 'bubbles' in the final assembly graph.
The genome size estimate is useful if the genome size is unknown, although it might be slightly over estimating the genome size due to excessive copies of chloroplasts in plant genomes or mitochondrial genomes, or some other unknown contaminants.

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
hifi=Col-0.hifi_50x.k21
illumina=Col-0.illumina.k21

meryl intersect $illumina.meryl $hifi.meryl output Col-0.illumina_and_hifi.50x.meryl
$MERQURY/plot/to_hist_for_plotting.sh Col-0.illumina.k21.meryl Illumina Col-0.illumina_and_hifi.50x.meryl Illumina_and_HiFi > Illumina_and_HiFi.hist
Rscript $MERQURY/plot/plot_spectra_cn.R -f Illumina_and_HiFi.hist -o Illumina_and_HiFi -t line
```

The `meryl intersect` is intersecting the k-mers found of `$hifi.meryl` with `$illumina.meryl`, keeping the counts of `$illumina.meryl`(the first meryl db provided).  
In the output `Illumina_and_HiFi.png`, the original k-mer counts from Illumina are shown in black, while the overlapping k-mers found also in the HiFi reads are in red.  

See `Rscript $MERQURY/plot/plot_spectra_cn.R --help` for more available options. The `-m` and `-n` features are useful sometimes if the histogram does not follow a usual "drop then peak and fall" model. [Here is an actual example.](https://github.com/marbl/merqury/issues/27)

_Would you agree the two datasets are coming from the same individual?_  
_At home: Try this on k-mers obtained from hifi reads of Cvi-0, and compare against the Illumina k-mers of Col-0_

## Assemble with Verkko
Verkko understands what has been already done, and tries to run based on what's missing.  
For the sake of time, we will only run the last 7-consensus step, assuming all directories from 0- to 6-layout are present.

```
cd ~/day3_assembly_evaluation/Col-0/

verkko -d asm \
  --hifi Col-0.hifi.q20.50x.fq.gz \
  --nano Col-0.ont_R10.gt_20kb.fastq.gz \
  --threads 24
```
This will create a directory (or assumes a directory) `asm`, and place all the output and intermediate files from verkko.
This step takes ~30 minutes, feel free to take a break or proceed by downloading the outputs described below from https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/scratch/a_thal/Col-0/asm/.

## Assembly graph
Download the `*.noseq.gfa` file, and open on your local bandageNG (or bandage, hereby just called bandage).  
Note the gfa loaded is a `unitig graph`, not the final `consensus graph`.

_Download the `*.noseq.gfa` from all other assemblies and compare._

Explore the output files.

The graph contains hifi read coverage by default for each unitig node. Add in the ONT coverage by downloading and loading `assembly.ont-coverage.csv` on bandage.

## Assemble with Verkko in trio mode
This time, let's try to understand how a trio-mode assembly works.
We will work on a synthetic Col-0 x Cvi-0 child, generated from 25x hifi and ONT data.
Let's move to `Col-0_x_Cvi-0`.
```
cd ~/day3_assembly_evaluation/Col-0_x_Cvi-0
```

First, we need to generate hapmers. Again, for the sake of time, we will skip this step.  
For those who are interested, this step counts k-mers, using k=30, in homopolymer compressed space.
```
cd meryl

# Below is shown to give an idea how the *.meryl folders were generated.
meryl count compress k=30 output Col-0.hifi_50x.hc.k30.meryl /path/to/Col-0.hifi.q20.50x.fq.gz
meryl count compress k=30 output Cvi-0.hifi_50x.hc.k30.meryl /path/to/Cvi-0.hifi.q20.50x.fq.gz
meryl count compress k=30 output Child.hifi_50x.hc.k30.meryl /path/to/Col-0.hifi.q20.25x.fq.gz /path/to/Cvi-0.hifi.q20.25x.fq.gz

# Because we aren't running the above lines, let's link the pre-run meryl dbs.
ln -s ~/day3_assembly_evaluation/Col-0/meryl/Col-0.hifi_50x.hc.k30.meryl
ln -s ~/day3_assembly_evaluation/Col-0/meryl/Cvi-0.hifi_50x.hc.k30.meryl
```

Now, we will run `hapmers.sh`, to obtain the `inherited` maternal (Col-0) and paternal (Cvi-0) strain.
```
mkdir hapmer_hc && cd hapmer_hc
ln -s ../Col-0.hifi_50x.hc.k30.meryl Col-0.meryl
ln -s ../Cvi-0.hifi_50x.hc.k30.meryl Cvi-0.meryl
ln -s ../Child.hifi_50x.hc.k30.meryl Child.meryl
$MERQURY/trio/hapmers.sh Col-0.meryl Cvi-0.meryl Child.meryl
```

_Can you interpret the `inherited_hapmers.ln.png` plot?_

Now, we are (pretending) we are running verkko in trio mode.
```
cd ~/day3_assembly_evaluation/Col-0_x_Cvi-0/

ln -s /opt/assembly_data/Col-0.hifi.q20.25x.fq.gz
ln -s /opt/assembly_data/Cvi-0.hifi.q20.25x.fq.gz
ln -s /opt/assembly_data/Col-0.ont_R10.gt_20kb.fastq.gz
ln -s /opt/assembly_data/Cvi-0.ont.gt_20kb.fastq.gz
ln -s ../meryl/hapmer_hc/Col-0.hapmer.meryl
ln -s ../meryl/hapmer_hc/Cvi-0.hapmer.meryl

# The entire assembly will take ~6 hours, let's skip this for now.
verkko -d trio \
  --hifi Col-0.hifi.q20.25x.fq.gz       Cvi-0.hifi.q20.25x.fq.gz \
  --nano Col-0.ont_R10.gt_20kb.fastq.gz Cvi-0.ont.gt_20kb.fastq.gz \
  --hap-kmers \
      Col-0.hapmer.meryl \
      Cvi-0.hapmer.meryl \
      trio \
  --threads 24
```

The trio mode output of verkko can be also found [here](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/scratch/a_thal/Col_x_Cvi/trio).

Download the `Col-0_x_Cvi-0` trio version `*.noseq.gfa`, and most importantly, the `unitig-popped-unitig-normal-connected-tip.colors.csv` from [6-rukki](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/scratch/a_thal/Col_x_Cvi/trio/6-rukki/). _Note: the colors.csv file will be available at the top level of verkko in later releases._

Load the `*.colors.csv` on bandage and see the haplotype identified nodes.
_Compare this before and after applying `rukki`._

### Advanced: Fix a path and re-do the consensus
Let's pull `utig4-928`. Place `utig4-928` in `Find Nodes` on bandage.

Pull out the path containig `utig4-928` used in rukki.

```
TRIO_PATH_GAF=trio/6-rukki/unitig-popped-unitig-normal-connected-tip.paths.gaf

head -n1 $TRIO_PATH_GAF > new_path.gaf
grep 928 trio/6-rukki/unitig-popped-unitig-normal-connected-tip.paths.gaf >> new_path.gaf
```

This is how `new_path.gaf` looks like:
>name	path	assignment
>Col-0_hifi_50x.k30.hapmer_from_utig4-928	<utig4-927[N25536N]>utig4-928	COL-0_HIFI_50X.K30.HAPMER

We need to edit the current path to replace the `[N25536N]` to the newly walked loop.
On bandage, you can select nodes in order and export them as a "path".
The new path will be `utig4-928-, utig4-925-, utig4-926+, utig4-925-, utig4-927+`. 

Unfortunately, verkko is using its own `gaf` format.
Replace the `-` or `+` to `<` and `>`, and place them in front of the node.  
Manually replace the path to `<utig4-928<utig4-925>utig4-926<utig4-925>utig4-927`.

Now, let's run verkko consensus using the `new_path.gaf`. This will generate only 1 contig as we are providing only 1 path.  
Note that we are re-using output from `trio`. The new consensus will be generated in `new_path`.
```
verkko -d new_path \
  --assembly trio \
  --paths new_path.gaf \
  --hifi Col-0.hifi.q20.25x.fq.gz Cvi-0.hifi.q20.25x.fq.gz \
  --nano Col-0.ont_R10.gt_20kb.fastq.gz Cvi-0.ont.gt_20kb.fastq.gz \
  --hap-kmers Col-0_hifi_50x.k30.hapmer.meryl Cvi-0_hifi_50x.k30.hapmer.meryl trio \
  --threads 24
```
## Closing remarks
Congratulations! You reached the end of this activity.

A lot of bug fixes have been made since the verkko v1.3.1 release. Note that the path fix we manually did will be fixed in the later release of verkko. Specifically, the example gap was caused by an internal coverage threshold, which requred >4 ONT reads to be mapped with MBG.

Still, there might be some edge cases where the graph can answer why the assembly did not went well.
We are happy to get feedbacks. What you think is a 'weird case' could be something _very_ useful for us to find and fix any unseen bugs in verkko. Feel free to report them as [issues](https://github.com/marbl/verkko/issues).

## Thanks to
* Sergey Koren
* Brian Walenz
* Brandon Pickett
* Steven Solar
