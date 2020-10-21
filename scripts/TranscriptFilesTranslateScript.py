#!/usr/bin/env python3
# Importing the required modules.
import sys
import os
from Bio import Seq

# Reads in a fasta file and puts the contents into a dictionary
# such that the gene name is the key and the gene sequence is the key value, and returns the dictionary.
def readFasta(infname):
  outputdick = {}
  handle = open(infname, "r")
  key = ""
  value = ""
  for line in handle:
    if line[0] == ">":
      if len(value) > 0:
        outputdick.update({key: value})
      key = line[1:len(line)].strip()
      value = ""
    else:
      value += line.strip()
  handle.close()
  return outputdick

def run(fastaFileNamesFilename):
  fastaFileNamesHandle = open(fastaFileNamesFilename, "r")
  for fastafileName in fastaFileNamesHandle:
    # print(fastafileName + "\n")
    fastafileName = fastafileName.strip()
    # Seq.translate("xxx", cds = True)
    # os.system("translate %s %s" % fastafileName, true)
    dic = readFasta(fastafileName)
    outputDic = {}
    for key in dic.keys():
      seq = dic[key]
      trans = Seq.translate(seq)
      outputDic[key] = trans


    #outputFileName = fastafileName.replace(".fasta", "translated.fasta")
    outputFileName = fastafileName.replace("cleaned", "translated")
    outputHandle = open(outputFileName, "w")
    for key in outputDic.keys():
      outputHandle.write(">" + key + "\n")
      outputHandle.write(outputDic[key] + "\n")

    outputHandle.close()
    print("")
    # fastaFileHandle=open(fastafileName, "r")
    # geneNames = ""
    # seq = ""
    # trans = ""
    # for line in fastaFileHandle:
    #   print(line)
    #   if line[0] == ">":
    #     if len(seq) > 0:
    #       trans = Seq.translate(seq, cds=True)
    #       print(geneNames)
    #       print(seq)
    #     geneNames = line[1:len(line)].strip()
    #     trans = ""
    #     seq = ""
    #   else:
    #     seq = seq+line.strip()
    #
    # fastaFileHandle.close()

  fastaFileNamesHandle.close()
  # print("done with run")



if __name__ == "__main__":
  if len(sys.argv) < 2:
      sys.stdout.write("usage: %s fastaFileNamesFilename\n" % sys.argv[0])
      sys.exit(1)
  # parse 2 command-line arguments #
  fname = sys.argv[1]
  run(fname)
