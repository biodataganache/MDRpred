#! /usr/bin/env python

# 
# MDRPred version 1.0
#   Prediction of multi-drug resistance transporters from protein sequence
#   Requires PyPro (https://code.google.com/p/protpy/wiki/propy) and
#            biopython (http://biopython.org/)
#
#   Usage: MDRPred.py [pattern filename] [input fasta filename] [output filename]
#   pattern filename: A tab-delimited file that includes a PyPro CTD designation string
#                     (column 1; NA if none) and a regular expression (column 2)
#   fasta filename:   Input sequences to be searched
#   output filename:  Output file. Format will be tab-delimited file with (column order):
#                         1.   sequence identifier 
#                         2.   match count
#                         3-38. Match score from each input pattern
#
#
#   For the MDRpred patterns as described in the paper (http://f1000research.com/articles/4-60/v1)
#        sequences with 36 matches are the highest confidence (positive predictive value ~0.3,
#        that is, 1 out of 3 positive predictions was a true positive in our hands)
#
#   Limitations:
#       The method as currently implemented does not predict if a protein is a transporter in general.
#       That is, a positive prediction by MDRpred should be examined by other means to determine if
#       it's likely to be a transporter-type protein. Spurious results might be output when applied
#       to non-transporter proteins. This will be fixed in future releases.

import os, sys, random, re
import argparse
from Bio import SeqIO

sys.path.append(os.path.join('.','..','..'))
sys.path.append("../../")

import propy.PyPro

def RegexSearch(r,M):
  try:
    Matches = [m for m in re.finditer(r,str(M)) if m!='']
    result = []
    for Match in Matches:
      result.append((Match.start(), Match.end(), 1))
    return(result)
  except Exception, e:
    return []

def AveragePhysioProperty(Property,r,M):
    Matches = [m for m in re.finditer(r,str(M)) if m != '']
    if len(Matches) == 0:
      return []
    else:
      result = []
      for Match in Matches:
        result.append((Match.start(), Match.end(), propy.PyPro.GetProDes(M[Match.start():Match.end()]).GetCTD()[Property]))
      return result

def parse_fasta(filename):
  id_list = []
  sequence_list = []
  handle = open(filename, "r")
  for record in SeqIO.parse(handle, "fasta"):
    id_list.append(record.id)
    sequence_list.append(str(record.seq))
  handle.close()

  return id_list, sequence_list
    
def main(infile, patternfile, outfile):
  models = []

  # parse the patternfile in to models
  # FIXME: no error checking is done for file format
  with open(patternfile, 'r') as fin:
    for line in fin:
      line = line.strip()
      fields = line.split("\t")

      if line != '':
        models.append(fields)

  Labels, Dataset = parse_fasta(infile)

  results_table = []

  if outfile:
    outfilehandle = open(outfile, "w")
  
  i = 0
  for this in range(len(Dataset)):
    seq = Dataset[this]
    thisthing = Labels[this]
    
    s = 0

    results_matrix = ""

    result = [0]*len(seq)

    hits = [0]*len(models)

    l = 0
    for model in models:
      result_n = [0]*len(seq)
      if (model[0] == 'NA'):
        r = RegexSearch(model[1], seq)
      else:
        r = AveragePhysioProperty(model[0], model[1], seq)

      if r:
        hits[l] = 1

      l += 1

      for (start,end,score) in r:
        for j in range(start, end):
          thiss = 0
          if score: thiss = 1

          result_n[j] += thiss
          result[j] += thiss

    #results_matrix = "%s%d\t%d\n" % (results_matrix, start, end)
    for k in result_n:
      st = "+"
      if k==0: st = "."
      results_matrix = "%s%s" % (results_matrix, st)

    results_matrix = "%s\n" % results_matrix

    #print "Sequence %d\tMax: %.1f\tLabel: %s" % (i,max(result),Labels[i])

    # normalize this to between 0 and 1
    if max(result) > 0:
      for q in range(0, len(result)):
        result[q] = result[q]/float(len(models))

    out = "%d\t%s\t%d" % (i, thisthing, sum(hits))
    for m in hits:
      out = "%s\t%d" % (out, m)

    #print out
    #print results_matrix
    
    out = "%s\t" % out
    
    #sys.stdout.write("%d\t" % i)
    for this in result:
      x = int(this*10)
      
      if x == 10: x = "+"
      out = "%s%s" % (out, x)

    out = "%s\n" % out
    if outfile:
      outfilehandle.write(out)
    else:
      print out
    
    i+=1
    results_table.append(result)

  if outfile:
    outfilehandle.close()
    

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='MDRpred v1.0 prediction of multidrug resistance transporters')
  parser.add_argument("-i", "--infile", help="FASTA format file to query")
  parser.add_argument("-p", "--patternfile", help="Pattern file for search", default="PILGram_PATTERNS_MDRpred.txt")
  parser.add_argument("-o", "--outfile", help="Output filename")

  arguments = parser.parse_args()
  main(arguments.infile, arguments.patternfile, arguments.outfile)
