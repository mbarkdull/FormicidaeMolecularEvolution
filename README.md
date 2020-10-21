Comparative Genomic Analysis Workflow
================
Megan Barkdull

## 1\. Introduction

This repository hosts the workflow for a comparative genomics analysis.
The general overview
is:

<img src="README_files/figure-gfm/unnamed-chunk-1-1.png" style="display: block; margin: auto;" />

### Using this workflow:

To use this workflow, simply clone this repository onto your own
machine. You will run all Bash scripts in the main directory associated
with this repository (i.e., in `/FormicidaeMolecularEvolution/`)

## Downloading data:

For a project of this kind, you will need transcript data from your
species of interest. The Bash script `DataDownload` can download these
data for you. This script requires an input file, which should be a
tab-delimited text file with:

  - Download URLs in the first column
  - The file names that you want your data downloaded to in the second
    column, formatted like:
    `FOURLETTERTAXONABBREVIATION_transcripts.fasta`.

See `./scripts/inputurls.txt` for an example.

Once you have an input URLs file, simply run the script by executing the
following command in a Bash shell:

`./scripts/DataDownload ./scripts/inputurls.txt`

This will create a new directory, `./RawData`, containing the downloaded
transcript files.

## Cleaning the raw data:

The raw data is likely to have several features that will make future
steps difficult or annoying. To solve this problem, you’ll want to run
the script `DataCleaning`. This will remove any special characters from
gene names (for example, \_ or (), which Orthofinder will change,
causing errors) and add your four-letter taxon abbreviation to the
beginning of each gene name. To run this script, simply use the command:

`./scripts/DataCleaning ./scripts/inputurls.txt`

This will create a new directory, `./CleanedData`, that contains the
cleaned transcript files.

For the step where you must convert the outputs of OrthoFinder to be
inputs for RERConverge, please see the [Comparative Genomics
repository](https://github.com/mbarkdull/ComparativeGenomics)
