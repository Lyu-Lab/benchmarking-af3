#!/usr/bin/env python
"""Extract ordered scores from OUTDOCK.

Michael Mysinger 200907 Created combine.py
Ryan Coleman 201006 edited heavily for newdb, now called one_extract.py
"""

import os
import sys
from optparse import OptionParser
import mmmutils
import extract_all

def get_scores_reader(initer, savelimit=None):
  '''reads an OUTDOCK file, newdb formatted, returns lists:
  returns full lists from outdock, 1 per successfully docked pose
  the second line is what happens if the ligand is not scored'''
  for line in initer:
    spl = line.strip().split()
    #if len(spl) == 19 or len(spl) == 20: #new OUTDOCK one-line file format.
    if len(spl) == extract_all.scoreCol: #new OUTDOCK one-line file format.
      if spl[0].isdigit() and len(spl[1]) >= 3 and spl[3].isdigit(): 
        if savelimit is None or \
            float(spl[extract_all.scoreCol - 1]) < savelimit: #passes score filter
          #uses scoreCol -1 since scoreCol is meant for use after the directory
          #column is prepended to each line by extract_all. so -1 here.
          yield spl

def get_scores_f(inf, outf, savelimit=None):
  mmmutils.write_splits(outf, get_scores_reader(inf, savelimit=savelimit), \
      raw_file=True)

def get_scores(infile=os.path.join('.', 'OUTDOCK'), outfile=None, \
    savelimit=None):
  infile = mmmutils.flex_open(infile)
  if outfile is None:
    outfile = sys.stdout
  else:
    outfile = open(outfile, 'w')
  get_scores_f(infile, outfile, savelimit=savelimit)

def main(argv):
  description = "Extract ordered scores from OUTDOCK."
  usage = "%prog [options]"
  version = "%prog *version 201006* created by Ryan Coleman"
  parser = OptionParser(usage=usage, description=description,
                          version=version)
  parser.set_defaults(infile=os.path.join('.', 'OUTDOCK'), outfile=None)
  parser.add_option("-i", "--infile",
           help="input file (default: %default)")  
  parser.add_option("-o", "--outfile",
           help="output file (default: stdout)")
  options, args = parser.parse_args(args=argv[1:])
  if len(args):
    parser.error("program takes no positional arguments.\n" +
                     "  Use --help for more information.")
  get_scores(infile=options.infile, outfile=options.outfile)
  return 0

if __name__ == '__main__':
  sys.exit(main(sys.argv))

