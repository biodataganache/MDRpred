# MDRpred
Code to predict multi-drug resistance transporters.

# 
# MDRPred version 1.0
#   Prediction of multi-drug resistance transporters from protein sequence
#   Requires PyPro (https://code.google.com/p/protpy/wiki/propy) and
#            biopython (http://biopython.org/)
#
#   Usage: MDRPred.py -p [pattern filename] -i [input fasta filename] -o [output filename] -h
#   pattern filename: A tab-delimited file that includes a PyPro CTD designation string
#                     (column 1; NA if none) and a regular expression (column 2)
#   fasta filename:   Input sequences to be searched
#   output filename:  Output file. Format will be tab-delimited file with (column order):
#                         1.   sequence identifier 
#                         2.   match count
#                         3-38. Match score from each input pattern
#                         39.  Match score for each sequence location. 0-9+ indicates the range of
#                              0-100% of the models matching that location.
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

