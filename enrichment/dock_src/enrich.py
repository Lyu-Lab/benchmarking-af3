#!/usr/bin/env python
"""Calculate enrichment and ROC curves.

Output:
<curvetype>.txt - Enrichment/ROC curve against the entire database
<curvetype>_own.txt - Enrichment/ROC curve against its own decoys
         <% searched> <% found>

Michael Mysinger 200706 Created
Michael Mysinger 200707 Split into enrich and combine
Michael Mysinger 200801 Update combine changes and add outdir
Michael Mysinger 200802 Account for all docked molecules and add skip_decoys
Michael Mysinger 200902 Create cleaner interface for reuse by enrich_outliers
Michael Mysinger 201005 Improve LogAUC calculation using analytical area integral
Ryan    Coleman  201006 Use extract_all instead of combine
Ryan    Coleman  201108 add ability to just get the logAUC value out
Trent   Balius   201305 LOGAUC_MAX and LOGAUC_MIN were added as bound for LOGAUC calculation
                        These variable make it easer to modified bounds.
Trent   Balius   201503 Enable the enrichment script to read in outdock with "zinc code . prot" 
                        as name and only use the zinc code.   
"""

import os
import sys
import math
from optparse import OptionParser
import extract_all
import mmmutils

# Input paths and filenames
DUD = '/raid7/people/rgc/dude_db2/defs'
#DUD = '/raid7/people/rgc/dud4/defs'
LIGAND_DIR = "ligands"
DECOY_DIR = "decoys"
LIGAND_EXT = ".name"
DECOY_EXT = ".name"

# constants
#RANDOM_LOGAUC = (1-0.001)/math.log(10)/math.log10(1.0/0.001)
## if you modify also change in plots.py        
LOGAUC_MAX = 1.0   ## this should not change
LOGAUC_MIN = 0.001 ## this you may want to change if you database is large and you have strong early enrichment. 
RANDOM_LOGAUC = (LOGAUC_MAX-LOGAUC_MIN)/math.log(10)/math.log10(LOGAUC_MAX/LOGAUC_MIN)
ROC = "roc.txt"
ROC_OWN = "roc_own.txt"
ENRICH = "enrich.txt"
ENRICH_OWN = "enrich_own.txt"

# Types
CURVE_TYPES = [(ROC, ROC_OWN),
               (ENRICH, ENRICH_OWN)]
ROC_CURVE_TYPES = [(ROC, ROC_OWN)]

def regen_dirs(indirs, forceit=False, skip_decoys=False, ligfile=None, \
               decfile=None, noEnrich=True):
    """Check if curves exist, regenerate as needed. (function used externally)"""
    if noEnrich:
        check_curve_types = [CURVE_TYPES[0]]
    else:
        check_curve_types = CURVE_TYPES
    if forceit:
        regen = set(indirs)
    else:
        all_missing = set()
        decoy_missing = set()
        for indir in indirs:
                for curvetype in check_curve_types:
                    alltype = curvetype[0]
                    filename = os.path.join(indir, alltype)
                    if not os.path.exists(filename):
                        all_missing.add(indir)
                    if not skip_decoys:
                        owntype = curvetype[1]
                        filename = os.path.join(indir, owntype)
                        if not os.path.exists(filename):
                            decoy_missing.add(indir)
        if skip_decoys:
            regen = all_missing
        elif not all_missing and len(decoy_missing) == len(indirs):
            print("Automatically skipping decoy curve generation.")
            print(( "  Rerun enrich.py without --skip-decoys on all inputs " + 
                    "to enable." ))
            skip_decoys = True
            regen = all_missing
        else:
            regen = all_missing | decoy_missing
    for indir in regen:
        print(("Enriching %s" % indir))
        gen_curves(indir=indir, forceit=forceit, skip_decoys=skip_decoys, \
            ligfile=ligfile, decfile=decfile)
    return skip_decoys

def read_dot_name(filename):
    """Read ids from .name file. (ZINC12345678 -> C12345678)"""
    f = open(filename, 'r')
    idstr = f.readline().split()[0]
    ids = [idstr] + [x.split()[0] for x in f]
    f.close()
    return ids

def write_points(filename, points, auc):
    """Write curve points to filename."""
    f = open(filename, 'w')
    f.write("#AUC\t%.2f\tLogAUC\t%.2f\n" % (auc[0]*100, auc[1]*100))
    [f.write("%.4f\t%.4f\n" % (x[0], x[1])) for x in points]
    f.close()

def get_names(indir='.', dud=DUD, ligfile=None, decfile=None):
    """Read ligand and decoy file names."""
    headtail = os.path.split(os.path.split(os.path.abspath(indir))[0])[1]
    enzyme = headtail.replace('_semi', '').replace('_auto', '')
    lfile = os.path.join(dud, LIGAND_DIR, enzyme + LIGAND_EXT)
    dfile = os.path.join(dud, DECOY_DIR, enzyme + DECOY_EXT)
    if ligfile is not None:
        lfile = ligfile
    if decfile is not None:
        dfile = decfile
    return lfile, dfile

def select_own(ligands, decoys, scores):
    """Select ligand ids and decoy ids from full ranked ids."""
    #scores format is full OUTDOCK line
    selected = set(ligands)
    selected.update(decoys)
    results = []
    for scoreline in scores:
        #id = scoreline[extract_all.zincCol] #refer to correct column always
        id = scoreline[extract_all.zincCol].split('.')[0] #refer to correct column always
                                                          # maybe in this form: zinccode.prot
        #print id 
        if id in selected:
            results.append(scoreline)
            #print scoreline
    return results

def enrich(dbsize, ligands, scores, nbins=10000):
    """Calculate enrichment curve given ranked ids and ligand ids."""
    ndata = len(scores)
    binsize = ndata/nbins + 1
    ligs = frozenset(ligands)
    nligs = len(ligs)
    count = 0
    results = []
    for i in range(ndata):
        if i % binsize == 0:
            results.append([i, count])
        # scores[i][extract_all.zincCol] maybe in this form: zinccode.prot
        # split on  "." take only the zinccode
        if scores[i][extract_all.zincCol].split('.')[0] in ligs: #i is the sorted rank, 0 is the id
            count += 1
    results.append([ndata, count])
    results.append([dbsize, nligs])
    results = [[x[0]*100.0/dbsize, x[1]*100.0/nligs] for x in results]
    return results

def roc(dbsize, ligands, scores, nbins=10000):
    """Calculate ROC curve given ranked ids and ligand ids."""
    ndata = len(scores)
    binsize = ndata/nbins + 1
    ligs = frozenset(ligands)
    nligs = len(ligs)
    count = 0
    results = []
    for i in range(ndata):
        if i % binsize == 0:
            results.append([i-count, count])
        # scores[i][extract_all.zincCol] maybe in this form: zinccode.prot
        # split on  "." take only the zinccode
        if scores[i][extract_all.zincCol].split('.')[0] in ligs: #i is the sorted rank, zincCol is the id
            count += 1
    results.append([ndata - count, count])
    ndecs = dbsize - nligs
    results.append([ndecs, nligs])
    results = [[x[0]*100.0/ndecs, x[1]*100.0/nligs] for x in results]
    return results

def interpolate_curve(points):
    """Interpolate curve where needed to get better graph visualization."""
    #    1) add interpolated 0.1 point (done)
    #    2) when both x and y change, and x changes by more than some 
    #           log space value, interpolate intermediate points 
    #           (XXX - still need to implement this!)
    #
    # interpolate in linear space to get y-value at x = 0.1% 
    i = 0
    while i < len(points) and points[i][0] < 0.1:
        i += 1
    slope = (points[i][1] - points[i-1][1])/(points[i][0] - points[i-1][0])
    intercept = points[i][1] - slope * points[i][0]
    point_one =  [0.100001, (slope * 0.100001 + intercept)]
    npoints = [x for x in points]
    npoints.insert(i, point_one)
    return npoints

def auc(points):
    """Calulate the area under the curve using trapezoid rule."""
    return sum((p[0]-lp[0])/100*(lp[1]+(p[1]-lp[1])/2.0)/100
                for p, lp in zip(points[1:], points[:-1]))

def logauc(points):
    """Compute semilog x AUC minus the perfectly random semilog AUC."""
    # assumes we have previously interpolated to get y-value at x = 0.1% 
    # generate new points array clamped between 0.1% and 100%
    #npoints = [[x[0]/100, x[1]/100] for x in points if 0.1 <= x[0] <= 100.0]
    npoints = [[x[0]/100, x[1]/100] for x in points if (LOGAUC_MIN*100) <= x[0] <= (LOGAUC_MAX*100)]
    area = 0.0
    for p, lp in zip(npoints[1:], npoints[:-1]):
        if p[0] - lp[0] < 0.000001:
            continue
        # segment area computed as integral of log transformed equation
        intercept = p[1] - (p[1]-lp[1])/(p[0]-lp[0]) * p[0]
        area += (p[1]-lp[1])/math.log(10) + intercept*(math.log10(p[0])-math.log10(lp[0]))
    # This normalization makes logauc independent of the log base
    #return area/math.log10(1.0/0.001) - RANDOM_LOGAUC
    return area/math.log10(LOGAUC_MAX/LOGAUC_MIN) - RANDOM_LOGAUC

def doaucs(points):
    """Append AUC and LogAUC calculations to aucfile."""
    nauc = auc(points)
    lauc = logauc(points)
    return (nauc, lauc)

def read_ids(idtype, idfile):
    """Read ids from idfile of type idtype.""" 
    print(("Using %s file: %s" % (idtype, idfile)))
    ids = read_dot_name(idfile)
    print(("%d %s read in." % (len(ids), idtype)))
    return ids

def just_roc_own_logAUC(outdir, dbsize, scores, ligands, ownsize=None, \
                        selected=None):
    """gets just the roc_own adjusted logAUC value. no writing files."""
    curvefunc = roc
    own = curvefunc(ownsize, ligands, selected)
    interp = interpolate_curve(own)
    auc, logAUC = doaucs(interp)
    return logAUC

def write_curves(outdir, dbsize, scores, ligands, ownsize=None, selected=None, 
                 skip_decoys=False, noEnrich=True):
    which_curves = CURVE_TYPES
    if noEnrich:
        which_curves = ROC_CURVE_TYPES
    for curvetype in which_curves:
        fname, ofname = curvetype
        if fname == ROC:
            curvefunc = roc
        else:
            curvefunc = enrich
        total = curvefunc(dbsize, ligands, scores)
        interp = interpolate_curve(total)
        write_points(os.path.join(outdir, fname), interp, doaucs(interp))
        if not skip_decoys:
            own = curvefunc(ownsize, ligands, selected)
            interp = interpolate_curve(own)
            write_points(os.path.join(outdir, ofname), interp, doaucs(interp))

def just_gen_curves(indir='.', outdir=None, dud=DUD, forceit=False,
                    ligfile=None, decfile=None, skip_decoys=False, 
                    receptors=[], part=None):
    """Actually just calculate enrichment and ROC curves, no writing."""
    if outdir is None:
      outdir = indir
    if receptors == []:
      receptors = None #later functions expect this to be None if all receptors
    # Combine files or read scores from disk 
    if receptors is None and part is None: #normal case, any receptor is fine
      scores = extract_all.get_scores_all(indir=indir, forceit=forceit)
    else:
      tempFile = os.path.join(outdir, extract_all.sortFileName)
      extract_all.read_scores_write( \
          os.path.join(indir, extract_all.sortFileName), tempFile, \
          recList=receptors, part=part) 
      scores = None #dump memory since it can be huge
      uniqRecFile = os.path.join(outdir, extract_all.uniqFileName)
      extract_all.uniqueify(tempFile, uniqRecFile)
      scores = extract_all.get_scores_all(indir=outdir, 
          whichFileName=extract_all.uniqFileName, forceit=forceit, 
          receptors=receptors, part=part)
    if scores is None:
        return False
    print(("total unique scores read:", len(scores)))
    lfile, dfile = get_names(indir=indir, dud=dud, 
                             ligfile=ligfile, decfile=decfile)
    dbsize = len(scores)
    if 0 == dbsize: #no ligands??
        print(("no ligands found in directory:", indir))
        return False
    ligands = read_ids("ligands", lfile)
    print(("%d of %d (%.1f%%) molecules in overall plot" % (len(scores), 
               dbsize, 100.0*len(scores)/dbsize)))
    if skip_decoys:
        selected = None
        ownsize = None
    else:
        decoys = read_ids("decoys", dfile)
        selected = select_own(ligands, decoys, scores)
        # XXX - should ownsize include bumped ligands and decoys?
        ownsize = len(set(ligands) | set(decoys))
        print(("%d of %d (%.1f%%) selected for own plot." % (len(selected), 
                   ownsize, 100.0*len(selected)/ownsize)))
    return dbsize, scores, ligands, ownsize, selected, outdir

def gen_curves(indir='.', outdir=None, dud=DUD, forceit=False, noEnrich=True, \
        ligfile=None, decfile=None, skip_decoys=False, receptors=[], \
        part=None):
    """Calculate enrichment and ROC curves. using ligfile and decfile (or whole
    database for background. dud controls auto-target parameters. noEnrich
    skips the obsolete enrich curve. receptors & part restrict data to 
    only some flexible receptors"""
    result = just_gen_curves(indir=indir, outdir=outdir, dud=dud, \
        forceit=forceit, ligfile=ligfile, decfile=decfile, \
        skip_decoys=skip_decoys, receptors=receptors, part=part)
    if result:
      dbsize, scores, ligands, ownsize, selected, outdir = result
      write_curves(outdir, dbsize, scores, ligands, ownsize, selected, 
                   skip_decoys=skip_decoys, noEnrich=noEnrich)
      return True
    else:
      return False

def get_roc_own_logAUC(indir='.', outdir=None, dud=DUD, forceit=False, 
               ligfile=None, decfile=None, skip_decoys=False):
    """Calculate just roc curve of own decoys, return logAUC"""
    result = just_gen_curves( \
        indir=indir, outdir=outdir, dud=dud, forceit=forceit, ligfile=ligfile, \
        decfile=decfile, skip_decoys=skip_decoys) 
    if result:
      dbsize, scores, ligands, ownsize, selected, outdir = result
      logAUC = just_roc_own_logAUC(outdir, dbsize, scores, ligands, ownsize, 
                                 selected)
      return logAUC
    else:
      return None

def enrich_parser(description, version, usage="%prog [options]"):  
    parser = OptionParser(usage=usage, description=description,
                          version=version)
    parser.set_defaults(indir='.', outdir=None, dud=DUD, force=False,
        ligfile=None, decfile=None, skip_decoys=False)
    parser.add_option("-i", "--indir",
           help="input directory for extract files (default: %default)")  
    parser.add_option("-o", "--outdir",
           help="output directory (default: --indir)")
    parser.add_option("-l", "--ligand-file", dest='ligfile',
           help="ligand id file (default: pull from DUD)")
    parser.add_option("-d", "--decoy-file", dest='decfile',
           help="decoy id file (default: pull from DUD)")
    parser.add_option("-s", "--skip-decoys", action="store_true",
           help="skip reading decoys and making own curves")
    parser.add_option("-b", "--base", dest='dud', 
           help="base DUD directory (default: %default)")
    parser.add_option("-f", "--force-it", action="store_true",
           dest="force", help="force extract_all to run again")
    parser.add_option("-r", "--receptor", type="string", action="append", \
           dest="receptor", default=[], \
           help="which receptor(s) to get poses for (default: All)")
    parser.add_option("-p", "--part", type="string", action="store", \
           dest="part", default=None, \
           help="which receptor part to get poses for (default: All)")
    return parser

def main(argv):
    description = "Calculate enrichment and ROC curves."
    version = "%prog *version 201203* created by Michael Mysinger, edits Ryan Coleman"
    parser = enrich_parser(description, version)
    options, args = parser.parse_args(args=argv[1:])
    if len(args) != 0:
        parser.error("program takes no positional arguments.\n" + 
                     "  Use --help for more information.")
    if options.outdir is not None and (not os.path.exists(options.outdir)): #create the outdir if necessary
        os.mkdir(options.outdir)
    passed = gen_curves(indir=options.indir, outdir=options.outdir, 
        dud=options.dud, forceit=options.force, ligfile=options.ligfile, 
        decfile=options.decfile, skip_decoys=options.skip_decoys, 
        receptors=options.receptor, part=options.part)
    return not passed

if __name__ == '__main__':
    sys.exit(main(sys.argv))
