Thien’s casual notes on dealing with various transcriptomic analyses in
R (for the BMEX paper)
================
miti
2024-03-04

- [Introduction](#introduction)
- [Step Zero: Downloading and processing RNA-seq data from NCBI Sequence
  Read
  Archive](#step-zero-downloading-and-processing-rna-seq-data-from-ncbi-sequence-read-archive)
  - [Prerequisites](#prerequisites)
  - [Prefetching and processing SRA
    files](#prefetching-and-processing-sra-files)
  - [Mapping and quantifiying reads from FASTQ
    data](#mapping-and-quantifiying-reads-from-fastq-data)
- [Cell-type deconvolution from RNA-seq
  data](#cell-type-deconvolution-from-rna-seq-data)
  - [Introduction](#introduction-1)
  - [Initial setup](#initial-setup)
  - [Running the deconvolution
    scripts](#running-the-deconvolution-scripts)
  - [Preliminary deconvolution
    results](#preliminary-deconvolution-results)
- [Consensus co-expression network
  analysis](#consensus-co-expression-network-analysis)

## Introduction

Hi. This is my notes on performing various kinds of analyses with
RNA-seq data in R. I do not guarantee good or understandable writing in
this, because this is basically my raw thoughts for documentation
purposes, and possibly for future Brain Health Lab younglings.

Future Thien might make it into a proper tutorial.

## Step Zero: Downloading and processing RNA-seq data from NCBI Sequence Read Archive

I wish it was as simple as clicking a download button, but science is
never simple.

### Prerequisites

- Your preferred command-line environment. I use Windows Subsystem for
  Linux (WSL) because Linux and bash is commonly used in scientific
  computing, so if you don’t have experience might as well learn it now.

- Basic experience with using the terminal.

- **A Ton** of free disk space (~500GBs or more), preferably on fast
  storage such as NVME SSDs because these operations are I/O heavy. HDDs
  might struggle.

- Lots and lots of RAM on your PC (minimum 32GB). If you don’t have this
  much you can ask for permission to use the Brain Health Lab server.

- ## Go to NCBI SRA and search for the dataset of your choice. For reference, here are the BioProject accessions that we are going to use:

### Prefetching and processing SRA files

Install the SRA-tools package from NCBI at the following link:

Go to the BioProject page and fetch the SRA accession list by going to
the Run selector and select download Accession list. A text file called
`SRR_Acc_List.txt` will be downloaded

*By default, `prefetch` will download files to the Home directory. You
might want to change this to whatever directory you want, or use the
current working directory. Do this by running `vdb-config -i` and change
it through the interactive menu.*

Run the command to prefetch all files from the downloaded SRA accession
list file: `prefetch --option-file ./SRR_Acc_List.txt`

This might take over 24 hours (for me it took 18 hours with the Toden
dataset). Just leave it running in the background while you do other
productive or unproductive things.

You will obtain a bunch of .sra files. We will now proceed to converting
them into the format that we need: FASTQ. To do this, we use the
`fasterq-dump` command. In our case, we have 338 SRA files so it will be
nice to set up a script to process them all at once. This step will take
up a lot of space (~2-3 times the original downloaded SRA data). I don’t
have enough disk space to process them all.

I found a package called `parallel-fastq-dump` that seems to be very
helpful.

### Mapping and quantifiying reads from FASTQ data

**To be written**, basically feed it into our original pipeline.

*Sidenote: The popular standard package for alignment seems to be
`salmon` or `kallisto` nowadays. I did not know this when I did the
original cf-RNA paper. But for reproducibility sake I will use our
original pipeline which used `rsubread`.*

## Cell-type deconvolution from RNA-seq data

### Introduction

Some people figured out how to extract cell-type data from RNA-seq data,
and applied it to plasma cell-free transcriptome. This is done by
referring to a whole body cell atlas called Tabula Sapiens. Disclaimer:
I don’t really know anything about the math behind deconvolution and
SVR.

We are going to do the same thing to our data. It may or may not produce
anything worthwhile, but let’s try anyway.

See reference

### Initial setup

- Install anaconda:
- Clone the github repository: sevahn/deconvolution
- Create the conda environment:

The authors recommend CPM-normalized values

### Running the deconvolution scripts

Navigate to the deconvolution tutorial folder

Edit the `sh_1.py` script file

Now, you can try to execute the resulting scripts, but Linux won’t let
you do that yet. You need to set the files to be executable. This
command will set all .sh files in the folder to executable
`chmod +x ./*.sh`

We can finally run the scripts. You can execute them all at once by
doing a for loop `for i in *.sh ; do ./$i ; done`

Edit `merge2.py` similarly to point to the correct data files. Finally,
run `merge2.py` to get the final fractions located in the fractions csv
file.

### Preliminary deconvolution results

Now that we have obtained the deconvolved fractions, let’s move on to
visualizing the results. We will import the resulting csv file back into
R. We will tidy up the data to make analysis easier using the tidyverse
packages.

To visualize the data, we will use possibly the greatest plotting
package ever to be invented: ggplot2. We want to make some boxplots to
see the full data distribution (I took a casual look at the data table
and saw a lot of variance between samples, so I don’t want to just
report the average). Let’s pick the top cell types by average fraction
for easier visualization.

We can see that the data is extremely varied. A lot of zero datapoints.

Cool! Next step is probably to perform proper statistical analyses to
determine if there is a statistically significant difference in
composition between control vs. AD plasma. Also, Tabula Sapiens does not
cover all cell types, most notably excluding brain cells. Human Protein
Atlas (HPA) RNA consensus dataset. Or we can narrow it down using
single-cell RNA sequencing datasets.

## Consensus co-expression network analysis

To be done. I have to process the hundreds of gigabytes of sequencing
data from Toden first.