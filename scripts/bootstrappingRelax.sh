#!/bin/bash

# Select a set of 100 random orthogroups that will be used for all 50 iterations:
ls 10_1_RelaxResults/workerReproductionQueens/*.json > allOrthogroups.txt
shuf -n 100 allOrthogroups.txt > subsetOrthogroupFiles.txt
cut -d'_' -f3 subsetOrthogroupFiles.txt > subsetOrthogroups.txt
cut -d'/' -f3 subsetOrthogroups.txt > subsetOrthogroupNumbers.txt
rm subsetOrthogroupFiles.txt
rm subsetOrthogroups.txt

echo Finished getting random orthogroups at `date` >> bootstrappingLog.txt

# Copy trees for randomly selected orthogroups to a directory:
rm -R ./bootstrappingTreeFiles
mkdir ./bootstrappingTreeFiles
while read -r line;
do
  export rawTreeFile=/workdir/mb2337/FormicidaeMolecularEvolution/allGeneTrees/$line"_tree.txt"
  cp $rawTreeFile ./bootstrappingTreeFiles
done < subsetOrthogroupNumbers.txt

# Trim out the previous foreground species:
Rscript ./scripts/trimmingForegroundSpeciesRelaxBootstraps.R $3

echo Selecting $2 foreground species >> bootstrappingLog.txt

# Run this fifty times (randomly labelling foreground, running RELAX):
i=0
while [ $i -ne 50 ]
do
  i=$(($i+1))
  echo "$i"

  # Note the start time:
  echo Starting $i round of bootstrapping at `date` >> bootstrappingLog.txt

  # Get a random list of X foreground species:
  shuf -n $2 ./scripts/inputurls_full.txt > subsetSpeciesFull.txt
  cut -d',' -f4 subsetSpeciesFull.txt > subsetSpecies.txt
  rm subsetSpeciesFull.txt
  echo Finished getting random species at `date` >> bootstrappingLog.txt

  # Label trees with your random species:
  echo "Labelling trees with foreground species"
  export LD_LIBRARY_PATH=
  rm -R ./9_1_LabelledPhylogenies/bootstrapping/
  Rscript ./scripts/LabellingPhylogeniesHYPHY.R /workdir/mb2337/FormicidaeMolecularEvolution/bootstrappingTreeFiles/trimmed subsetSpecies.txt bootstrapping

  echo Finished labelling trees with random species at `date` >> bootstrappingLog.txt

  # Run RELAX on those orthogroups/trees
  # Export required paths:
  export LD_LIBRARY_PATH=/usr/local/gcc-7.3.0/lib64:/usr/local/gcc-7.3.0/lib

  # Make a directory for BUSTED[S] outputs:
  mkdir ./10_1_RelaxResults
  mkdir ./10_1_RelaxResults/bootstrapping

  # Add HyPhy to your path:
  module load gcc/7.3.0
  export PATH=/programs/hyphy-20210923/bin:$PATH
  cd ./10_1_RelaxResults
  git clone https://github.com/veg/hyphy-analyses.git
  cd ../

  # List our sequence files for the boostrap run:
  #echo "List our sequence files for the boostrap run:"
  rm sequencefileList.txt
  while read -r line;
  do
    export removedStopFile=/workdir/mb2337/FormicidaeMolecularEvolution/8_2_RemovedStops/sequences/cleaned_$line"_cds.fasta"
    echo $removedStopFile >> sequencefileList.txt
  done < subsetOrthogroupNumbers.txt

  # Split the list of allowable input files into a specified number of chunks:
  echo "Split the list of allowable input files into a specified number of chunks:"
  export chunkNumber="l/"$1
  split --number=$chunkNumber --additional-suffix=.txt -d sequencefileList.txt bootstrapping_relaxFileList

  # Create a file listing those chunks:
  echo "Create a file listing those chunks:"
  ls bootstrapping_relaxFileList* > bootstrapping_chunkList.txt

  # Create a holder for the chunks:
  export batchSize=$1
  export currentBatch=0
  export batchFileNames=()

  # while reading each line in our list of chunked files,
  #echo "Running RELAX"
  while read -r line;
  do
    export batchFile=$line
    batchFileNames+=($batchFile)

    if [ ${#batchFileNames[@]} -eq $batchSize ]; then
      for batchFile in ${batchFileNames[@]} ; do
        echo $batchFile
        sleep 10 &
        ./scripts/SingleRELAXRun $batchFile  9_1_LabelledPhylogenies/bootstrapping/ bootstrapping &
      done
      wait
      batchFileNames=()
    fi
  done < bootstrapping_chunkList.txt

  mv 10_1_RelaxResults/bootstrapping/ "10_1_RelaxResults/"$i"bootstrapping/"
  mv subsetSpecies.txt $i"subsetSpecies.txt"
done
