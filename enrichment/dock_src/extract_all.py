#!/bin/env python
"""run extract in each dir in dirlist in input dir. cat them together. 
make separate files with sorted scores and sorted then uniqued scores.
Ryan Coleman 2011-2012

Modifed by Trent Balius, 2014-2015 to support GIST.
"""

import os
import sys
import operator
import string
from optparse import OptionParser
import collections
import math
import bsddb
import marshal
import subprocess

# Common modules
import mmmutils
import checkdir
import one_extract
# mol2extend and pdb are only needed for pose-preference helpers not used by enrich.py

#constants for this module
outFileName = 'extract_all.txt'
outFileNameBdb = 'extract_all.txt.bdb'
sortFileName = 'extract_all.sort.txt'
uniqFileName = 'extract_all.sort.uniq.txt'
outdockName = 'OUTDOCK'
scoreCol = 21 #constants about extract_all.txt file, where total score is
recScoreCol = 18 #receptor energy totaled and saved here
otherRecScoreCols = [12, 13, 14, 15, 16, 18, 19, 20] #these plus recScoreCol equals scoreCol
zincCol = 2 #where the ZINC id is
receptorCol = 3 #where the receptor code column is at
molnumCol = 1 #where the molecule number is
dirCol = 0 #the directory the docking run was in
rankCol = 10 #the column the rank is in
heavyAtomCol = 7 #the column the heavy atom count is in
identicalCheckCols = [8, 9, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21] #columns that are identical if the pose is identical

kT = 0.593

def readIds(idfilename, receptor=None, part=None):
  '''reads a list of ids, puts it in a list, returns it.
  idfilename can be an extract_all.* file, in which case you can use the
  receptor and part options to limit the ids saved & returned.
  idList returned is unique codes in order they are read in.'''
  idList = []
  idSet = set() #want to keep track of added and order, so just do both
  try:
    print(("reading list of identifiers from file:", idfilename))
    idFile = open(idfilename, 'r')
    for line in idFile:
      tokens = string.split(line)
      if len(tokens) > 10: #outdock formatted file
        if (receptor is None) or (tokens[receptorCol] in receptor):
          if (part is None) or (part in string.split(tokens[receptorCol], '.')):
            id = string.strip(str(tokens[2]))
            if id not in idSet:
              idList.append(id)
              idSet.add(id)
      else: #just codes probably
        if len(tokens) > 0: #allow empty lines
          id = string.strip(tokens[0])
          if id not in idSet:
            idList.append(str(id))
            idSet.add(str(id))
  except StopIteration:
    pass
  return idList

def removeDuplicates(scoreList):
  '''list of score lines as input. output is a list of score lines, with 
  duplicate poses/scores removed'''
  newList = []
  for aScore in scoreList:
    okayToAdd = True
    for otherScore in newList:
      matches = True
      for checkCol in identicalCheckCols:
        if aScore[checkCol] != otherScore[checkCol]: 
          matches = False
          break
      if matches:
        okayToAdd = False
        break
    if okayToAdd:
      newList.append(aScore)
  return newList

def energiesToPreference(inList, outList):
  '''reads the energy column, does e^(-E/kT) for each, sums it for in and 
  divides by the sum for in+out to return a preference'''
  inScores = []
  outScores = []
  for aScore in inList:
    inScores.append(aScore[scoreCol])
  for aScore in outList:
    outScores.append(aScore[scoreCol])
  sumIn = 0.0
  for inScore in inScores:
    sumIn += math.exp(-inScore / kT)
  sumOut = sumIn #start with total for in
  for outScore in outScores:
    sumOut += math.exp(-outScore / kT)
  return sumIn / sumOut

def getReceptors(scoreList):
  '''for each line, grab out receptor col, put in set'''
  receptors = set()
  for aScore in scoreList:
    receptors.update(string.split(aScore[receptorCol]))
  receptorList = list(receptors)
  receptorList.sort() #sorted list
  return receptorList

def getParts(scoreList):
  '''for each line, grab out the receptor col, look at the parts, put in set'''
  parts = set()
  for aScore in scoreList:
    try:
      moreParts = string.split(aScore[receptorCol], '.')
      parts.update(moreParts)
    except AttributeError:
      parts.update([aScore[receptorCol]]) #means there are no parts, just invariant
  partList = list(parts) #turn into list for sorting
  partList.sort() 
  return partList

def informationContent(inputList):
  '''does Sum x log x for all x in inputList'''
  total = 0.0
  for data in inputList:
    total += data * math.log(data)
  return total

def makeAllMol2(scoreList, scorePose):
  '''using the ligands in scorelist order, use scorepose dict to make a 
  combined mol2 file'''
  mol2textIn = []
  for position in range(len(scoreList)):
    try:
      thisPoseText = scorePose[tuple(scoreList[position])]
      mol2textIn += thisPoseText
    except KeyError:
      pass #not all scores may have been extracted if they aren't needed. fine.
  import mol2extend
  allMols = mol2extend.Mol2(mol2text=mol2textIn)
  return allMols

def partPreferencesBlackBox(scoreList, scorePose, ndist=10., dof=1., \
    clusterLimit=300, protrmsd=True):
  '''another scheme for getting part preferences. based on Ytreberg & Zuckerman.
   PNAS 2008, 102 no 23. thanks to Gabe Rocklin for pointing it out.
  basic idea is count things less if they are closer to other sampled points.
  use RMSD of ligands as distance.
  scoreList input is a list of all the scores for one ligand.
  output is a dict of part -> %age in that part by energy.
  ndist is how many samples to 'combine' to reweight and dof is the degrees
  of freedom.'''
  parts = getParts(scoreList)
  startRmsdTable = collections.defaultdict(dict)
  if protrmsd: #use lig + prot rmsd, not just lig
    import pdb
    recToPdb = {}
    recRecRmsd = collections.defaultdict(dict)
    recs = getReceptors(scoreList)
    for rec in recs:
      pdbDataIn = pdb.pdbData(rec + '.pdb')
      recToPdb[rec] = pdbDataIn
    for recCount, rec in enumerate(recs):
      recRecRmsd[rec][rec] = 0.0 #self rmsd always 0.0
      for otherRec in recs[recCount+1:]:
        rmsd = recToPdb[rec].calcRMSD(recToPdb[otherRec])
        recRecRmsd[rec][otherRec] = rmsd
        recRecRmsd[otherRec][rec] = rmsd
    for poseCount, poseList in enumerate(scoreList):
      for poseOtherCount in range(poseCount + 1, len(scoreList)): #second part
        thisRmsd = recRecRmsd[poseList[receptorCol]][ \
            scoreList[poseOtherCount][receptorCol]] 
        startRmsdTable[poseCount][poseOtherCount] = thisRmsd
        startRmsdTable[poseOtherCount][poseCount] = thisRmsd
  #make one mol2 molecule with each pose in it
  allMols = makeAllMol2(scoreList, scorePose)
  #we have to force the advanced RMSD to take care of different protonation
  #states, but this also nicely takes care of symmetry in rings and terminal 
  #groups, etc.
  if clusterLimit is not None and \
      len(allMols.atomXyz) > clusterLimit: #want to divisively cluster first
    if protrmsd:
      rmsdTable = allMols.getRMSDtable(forceRedo=True, advanced=True, \
          clusterLimit=clusterLimit, startRmsdTable=startRmsdTable)
    else:
      rmsdTable = allMols.getRMSDtable(forceRedo=True, advanced=True, \
          clusterLimit=clusterLimit)
  else: #normal RMSD, can do all pairs in reasonable time
    if protrmsd: #add prot rmsd to lig
      rmsdTable = allMols.getRMSDtable(forceRedo=True, advanced=True, \
          startRmsdTable=startRmsdTable)
    else:
      rmsdTable = allMols.getRMSDtable(forceRedo=True, advanced=True)
  #possible correct ndist to be higher or lower. can't be lower than #0s in any
  # row, can't be higher than the length of any row.
  minN = 0
  maxN = sys.maxsize
  for position in range(len(scoreList)):
    distances = list(rmsdTable[position].values())
    minN = max(minN, distances.count(0))
    maxN = min(maxN, len(distances))
  if ndist < minN:
    ndist = float(minN)
  if ndist > maxN:
    ndist = float(maxN)
  #print minN, maxN, ndist
  uncorrectedWeight = []
  correctedWeight = []
  for position in range(len(scoreList)):
    boltzEnergy = math.exp(-scoreList[position][scoreCol] / kT)
    uncorrectedWeight.append(boltzEnergy)
    distances = list(rmsdTable[position].values())
    distances.sort()
    rhyp = distances[int(ndist)]
    pobs = ndist/(rhyp**dof)
    #print boltzEnergy, rhyp, pobs       #debugging!
    correctedWeight.append(boltzEnergy / pobs)
  partPrefs = {}
  for aPart in parts:
    inList = []
    for position, aScore in enumerate(scoreList):
      try:
        if aPart in string.split(aScore[receptorCol], '.'):
          inList.append(correctedWeight[position])
      except AttributeError: #no flexible parts
        pass
    partPrefs[aPart] = sum(inList) / sum(correctedWeight)
  return partPrefs, ndist, dof

def partPreferences(scoreList):
  '''foreach flexible receptor part, calculate the partition function for each.
  use this part / all other parts and try to determine identical poses so no 
  double counting happens. would be better inside DOCK with all poses.
  scoreList input is a list of all the scores for one ligand.
  output is a dict of part -> %age in that part by energy.'''
  #first gather all the reccodes, break them up, get the parts.
  parts = getParts(scoreList)
  #now for each part, split scoreList in twain, for ones containing it and not
  #do the math, add the %age to the dictionary
  partPrefs = {} #where data ends up
  for aPart in parts:
    inList = []
    outList = []
    for aScore in scoreList:
      try:
        if aPart in string.split(aScore[receptorCol], '.'):
          inList.append(aScore)
        else:
          outList.append(aScore)
      except AttributeError: #no flexible parts
        if aPart == aScore[receptorCol]:
          inList.append(aScore)
        else:
          outList.append(aScore)
    #prune duplicates, if they have the same energy and look to be the same pose
    inNoDuplicates = removeDuplicates(inList)
    outNoDuplicates = removeDuplicates(outList)
    pref = energiesToPreference(inNoDuplicates, outNoDuplicates)
    partPrefs[aPart] = pref
  return partPrefs

def rankScores(scoreToCheck, indir='.', whichFileName=None, forceit=False, \
               cacheRanked=None):
  '''caches the scores from the ranking file, returns the rank of the code.
  works based on energy since a lot of poses will not be in the uniq file since
  every pose is saved (and may be desired). so this is a global rank estimated
  by score.'''
  if (cacheRanked is None) or forceit: #read from file if not loaded
    rankScores = get_scores_all(indir, whichFileName, forceit)
    cacheRanked = []
    for oneScore in rankScores:
      cacheRanked.append(float(oneScore[scoreCol]))
      #cacheRanked.append(float(oneScore[-1])) #this is a total hack 
  returnRank = 0
  for count, score in enumerate(cacheRanked):
    if scoreToCheck >= score:
      returnRank = count
    else:
      break #we can quit, we found the appropriate rank
  return returnRank + 1, cacheRanked

def get_scores_all(indir='.', whichFileName=None, forceit=False, 
        receptors=None, part=None):
  """Read in combined scores or generate them as needed."""
  if whichFileName is None: #default is uniqFileName
    scorename = os.path.join(indir, uniqFileName)
  else:
    scorename = os.path.join(indir, whichFileName)
  if forceit or not os.path.exists(scorename):
    extract_all(indir=indir)
  try:
    scores = readExtract(scorename, recList=receptors, part=part)
  except IOError:
    scores = None
  return scores

def str2intOrFloat(input):
  '''helper method, processes string into ints or floats if possible'''
  try:
    out = int(input)
  except ValueError:
    try:
      out = float(input)
    except ValueError:
      out = input
  return out

def getId2Scores(indir='.', idlist=None, recList=None, part=None):
  '''read extracted scores, make a map from zincid -> score tokens.
  if idlist is given, only keep those ligands, discard the rest.
  if part is given, only return receptors with that part.
    first #idLimit ids (actually returns one extra due to laziness).
  caches scores found in extract_all.txt.bdb using marshal & bsddb'''
  okay = "OKAY_EXTRACT" #use as key to detect if we should trust the bsddb or not
  id2scoresCache = bsddb.hashopen(os.path.join(indir, outFileNameBdb)) #opens new
  #if no okay key present
  #print "idlist", idlist
  if (not (okay in id2scoresCache) or \
      (not marshal.loads(id2scoresCache[okay]))): #or if key is not true
    #clear the database by closing, deleting the file, and reopening
    id2scoresCache.close()
    os.unlink(os.path.join(indir, outFileNameBdb)) #delete the db cache
    id2scoresCache = bsddb.hashopen(os.path.join(indir, outFileNameBdb)) 
    #print "deleting and reopening", os.path.join(indir, outFileNameBdb)
    #fresh empty one opened
  id2scoresCache[okay] = marshal.dumps(False) #set to false until everything read
  scores = []
  if idlist is None: #have to get every single one, this will likely crash
    scores = getExtractedScores(indir)
    idset = set()
    for aScore in scores:
      idset.add(aScore[zincCol])
    idlist = list(idset) #now set for alter
  else:
    newScores = None
    for oneId in idlist:
      if str(oneId) not in id2scoresCache: #not already read, so read it
        #print "missing key", oneId
        if newScores is None:
          newScores = getExtractedScores(indir, idlist) #read all at once
        for aScore in newScores:
          if aScore[zincCol] == oneId:
            scores.append(aScore)
        id2scoresCache[str(oneId)] = marshal.dumps([]) # clear this key
  #have to write all scores to id2scoresCache now
  #next line caches the last value read, so repeated reading from disk 
  # happens less and is less costly
  lastId, lastValue = None, None
  for aScore in scores:
    valuesList = []
    for scoreToken in aScore:
      valuesList.append(str2intOrFloat(scoreToken))
    thisId = aScore[zincCol]
    aValue = [] #init to empty list if nothing there yet
    if lastId == thisId:
      aValue = lastValue
    elif str(thisId) in id2scoresCache:
      aValue = marshal.loads(id2scoresCache[str(thisId)])
    aValue.append(valuesList)
    id2scoresCache[str(thisId)] = marshal.dumps(aValue)
    lastId = thisId #cache last item only
    lastValue = aValue #cache last item only, keeps from repeated disk reads
  id2scoresCache.sync() #sync now, just in case
  id2scores = collections.defaultdict(list) #return dictionary
  for oneId in idlist:
    allScores = marshal.loads(id2scoresCache[str(oneId)])
    for scoreTokens in allScores: #scoreTokens is one score line
      if (recList is None) or (scoreTokens[receptorCol] in recList):
        if (part is None) or \
            (part in string.split(scoreTokens[receptorCol], '.')):
          id2scores[str(scoreTokens[zincCol])].append(scoreTokens)
  id2scoresCache[okay] = marshal.dumps(True) #set to true, meaning all desired read
  id2scoresCache.sync()
  id2scoresCache.close()
  return id2scores

def getExtractedScores(indir='.', idlist=None, recList=None, idLimit=None):
  """read extracted scores or if they aren't there, call extract_all.
  if idlist present, only read in scores with matching ids.
  if recList present, only read in scores with matching receptor code.
  if idLimit is none, return all, otherwise only return the first # of codes."""
  filename = os.path.join(indir, outFileName)
  if (not os.path.exists(filename)):
    extract_all(indir=indir)
  try:
    scores = readExtract(filename, idlist, recList, idLimit=idLimit)
  except IOError:
    scores = None
  return scores

def readExtract(inFileName, idlist=None, recList=None, part=None, idLimit=None):
  """read an extract_all.sort.uniq.txt format file, return as list of tokens.
  if idlist present, only read in scores with matching ids.
  if recList is present, only read in scores with matching receptor codes.
  if part is present, only read in scores with matching receptor part.
  if idLimit is none, return all, otherwise only return the first # of codes."""
  print(("reading extract info from", inFileName))
  outList = []
  infile = open(inFileName, 'r')
  if idLimit is not None:
    keyCount = set()
  try:
    for line in infile:
      tokens = line.split()
      if (idlist is None) or (tokens[zincCol] in idlist):
        if (recList is None) or (tokens[receptorCol] in recList):
          if (part is None) or (part in string.split(tokens[receptorCol], '.')):
            outList.append(tokens)
            if idLimit is not None:
              keyCount.add(str(tokens[zincCol]))
              if len(keyCount) > idLimit:
                break #get out of the for loop altogether
  except StopIteration:
    pass
  infile.close()
  return outList

def read_scores_rescore_write(inFileName, outFileName, newScoresDict, \
     dockmultiplier=1.0):
  '''reads scores from inFileName, writes to outFileName after rescoring
  using scores from newScoresDict to rescore each flexible part.'''
  infile = open(inFileName, 'r')
  outfile = open(outFileName, 'w')
  try:
    for line in infile:
      tokens = line.split()
      parts = string.split(tokens[receptorCol], '.')
      newRecEnergy = 0.0
      for aPart in parts:
        newRecEnergy += newScoresDict[aPart]
      #newRecEnergy has the new score
      newTotalScore = newRecEnergy #start with this
      for otherScoreCol in otherRecScoreCols: #add all others
        newTotalScore += dockmultiplier * float(tokens[otherScoreCol])
        tokens[otherScoreCol] = str(dockmultiplier * \
            float(tokens[otherScoreCol]))
      #replace in tokens
      tokens[scoreCol] = str(newTotalScore)
      tokens[recScoreCol] = str(newRecEnergy)
      #write to outfile by joining strings with spaces
      outfile.write(string.join(tokens, " ") + '\n')
  except StopIteration:
    pass
  infile.close()
  outfile.close()

def read_scores_write(inFileName, outfileName, \
        idlist=None, recList=None, part=None):
  """read an extract_all.sort.uniq.txt format file, write to outfile
  if idlist present, only read in scores with matching ids.
  if recList is present, only read in scores with matching receptor codes.
  if part is present, only read in scores with matching receptor part"""
  print(("reading extract info from", inFileName))
  infile = open(inFileName, 'r')
  outfile = open(outfileName, 'w')
  try:
    for line in infile:
      tokens = line.split()
      if (idlist is None) or (tokens[zincCol] in idlist):
        if (recList is None) or (tokens[receptorCol] in recList):
          if (part is None) or (part in string.split(tokens[receptorCol], '.')):
            outfile.write(string.join(tokens, " ") + '\n')
  except StopIteration:
    pass
  infile.close()
  outfile.close()

def extract_all(indir='.', savelimit=None, doneflag=False):
  """check if done, if all done, run extract, combine extract.
  new method is run the ones you can run that are finished, but don't combine
  until they are all done."""
  allDone = True  
  for subdir in mmmutils.read_dirlist(indir): 
    if (checkdir.docheckdir(subdir) or doneflag):
      infile = os.path.join(subdir, outdockName)
      outfile = os.path.join(subdir, outFileName) 
      if not os.path.exists(outfile): #only do this if the file isn't already 
        one_extract.get_scores(infile, outfile, savelimit=savelimit) #created by this call
    else:
      allDone = False
  #print "doneflag: ", doneflag
  #print "allDone: ", allDone
  if not (allDone or doneflag):
    print("Error! The above jobs are not done, use --done to override!")
    return False
  allfilename = os.path.join(indir, outFileName)
  allfile = open(allfilename, 'w')
  for subdir in mmmutils.read_dirlist(indir):
    infileName = os.path.join(subdir, outFileName)
    infile = open(infileName, 'r')
    try:
      for data in infile:
        writeOut = True
        if savelimit is not None:
          thisScore = float(string.split(data)[scoreCol - 1]) #scorecol only works
              #after this file is written and an additional column is added at 
              #the beginning, so use -1 here
          if thisScore > savelimit:
            writeOut = False #don't save if over the savelimit
        if writeOut:
          allfile.write(os.path.basename(subdir) + '/ ')  #always append / here
          allfile.write(data)
    except StopIteration:
      pass
    infile.close()
  allfile.close()
  try:
    os.remove(os.path.join(indir, sortFileName))
  except OSError:
    pass #no big deal
  #add 1 to scoreCol since sort is 1-indexed not 0
  sortify(os.path.join(indir, outFileName), os.path.join(indir, sortFileName))
  try:
    os.remove(os.path.join(indir, uniqFileName))
  except OSError:
    pass #no big deal
  #add 1 to zincCol again since awk is 1-indexed
  uniqueify_mod(os.path.join(indir, sortFileName), \
      os.path.join(indir, uniqFileName))
  print(("number of ligands extracted:", \
      subprocess.getoutput("wc -l " + os.path.join(indir, uniqFileName))))
  return True

def write_scores(outputFileName, scores):
  '''writes scores, list of lists, in space delimited format, out to a file'''
  outfile = open(outputFileName, 'w')
  for score in scores:
    outfile.write(string.join(score, " ") + '\n')
  outfile.close()

def sortify(inputFile, outputFile):
  '''uses unix sort to sort by total score, reading in takes too much memory'''
  os.popen("sort -k " + str(scoreCol + 1) + " -n " + inputFile + 
           " > " + outputFile)

def uniqueify(inputFile, outputFile):
  '''uses awk to keep only the first uniq zinc code for each entry'''
  os.popen("awk '!($" + str(zincCol + 1) + "in azinc) {azinc[$" + \
       str(zincCol + 1) + "] ; print}' " + inputFile + " > " + outputFile)

#def uniqueify_zinc_with_prot(inputFile, outputFile):
def uniqueify_mod(inputFile, outputFile):
  ''' keeps only the first uniq zinc code regradless of prot for each entry'''
  '''this function assumes that the list is aready sorted '''
  #os.popen("awk '!($" + str(zincCol + 1) + "in azinc) {azinc[$" + \
  #     str(zincCol + 1) + "] ; print}' " + inputFile + " > " + outputFile)
  infileh = open(inputFile,'r')
  outfileh = open(outputFile,'w')
  zinc_dic = {}
  for line in infileh:
      splitline = line.split() 
      compoundName = splitline[zincCol] # compound name may have this formate zinc_code.prot
      zinccode = compoundName.split('.')[0] # just get zinc. 
      if not (zinccode in zinc_dic): 
         outfileh.write(line)
         zinc_dic[zinccode] = 1
      else:
         zinc_dic[zinccode] = zinc_dic[zinccode] + 1

  # print how meny times it appears in the file.
  #for key in zinc_dic.keys():
  #    print key, zinc_dic[key]

  infileh.close()
  outfileh.close()    

def main(argv):
  description = "Run extract on all directories in dirlist in input directory"
  usage = "%prog [options]"
  version = "%prog *version 201004* created by Ryan Coleman"
  parser = OptionParser(usage=usage, description=description,
                          version=version)
  parser.set_defaults(indir='.')
  parser.add_option("-i", "--indir",
           help="check results inside this directory (default: %default)")  
  parser.add_option("-d", "--done",\
           action="store_true", dest="done", \
           default=False,\
           help="this flag will override the done check.  thing will be processed even if docking is not completed (default: %default)")
  parser.add_option('-s', '--savelimit', type='float', action='store', \
      dest='savelimit', default=None, \
      help="only extract poses with scores better than this. makes prospective extracts faster, don't use retrospectively on small sets (default: %default)")
  options, args = parser.parse_args(args=argv[1:])
  if len(args):
    parser.error("program takes no positional arguments.\n" +
                     "  Use --help for more information.")
  ranOkay = extract_all(indir=options.indir, savelimit=options.savelimit, doneflag = options.done)
  return not ranOkay

if __name__ == '__main__':
  sys.exit(main(sys.argv))
