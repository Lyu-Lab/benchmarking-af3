#!/bin/env python
"""Plot enrichment and ROC curves into images.

Michael Mysinger 200710 Created
Michael Mysinger 200803 Add skip_decoys and normal options
Michael Mysinger 200902 Enable custom text on plots
"""

import matplotlib  # must import first
matplotlib.use('Agg')  # allows you to not have an x-server running
#these lines must be first, if pylab is imported first it ruins this

import os
import sys
from optparse import OptionParser
import pylab
import mmmutils
import enrich

RANDOM_FILE = "random.txt" 

def get_types(filename):
    """Get type information based on the curve filename."""
    if filename == enrich.ROC or filename == enrich.ROC_OWN:
        plottype = "ROC"
        searchtype = "Decoys Found"
    else:
        plottype = "Enrichment"
        searchtype = "Ranked Hits Searched"
    if filename == enrich.ROC_OWN or filename == enrich.ENRICH_OWN:
        dbtype = "Own Decoys"
    else:
        dbtype = "Whole Database"
    return plottype, searchtype, dbtype

def read_points(filename, normal=False):
    """Read header and curve points from filename."""
    splits = list(mmmutils.read_splits(filename))
    if normal:
        auc = splits[0][1]
    else:
        auc = splits[0][3]
    points = [[float(j) for j in i] for i in splits[1:]]
    x, y = list(zip(*points))
    return x, y, auc

def plot_curves(filename, indirs, labels=[], colors=[], target=None, outdir=".",
                guess=False, skip_decoys=False, normal=False,
                title=None, xlabel=None, ylabel=None, hide_auc=False,
                legend_font_size=14, legend_location='best'):
    """Plot input enrichment curves to generate an image."""
    plottype, searchtype, dbtype = get_types(filename)
    if skip_decoys:
        dbtype = ""
    else:
        dbtype = " - " + dbtype
    if not labels:
        labels = []
        for indir in indirs:
            dname = os.path.dirname(os.path.abspath(indir))
            labels.append(os.path.basename(dname))
    while len(colors) < len(indirs): #set colors to be None if none present
        colors.append(None)
    for indir, label, color in zip(indirs, labels, colors):
        fname = os.path.join(indir, filename)
        try:
            x, y, auc = read_points(fname, normal=normal)
        except IOError:
            print(("No file present, skipping in plot:", fname))
            continue #try to make a plot, even an incomplete plot
        if hide_auc:
            ltext = label.replace('_',' ')
        else:
            ltext = "%s: %s" % (label.replace('_',' '), auc)
        if normal:
            if color is not None:
                pylab.plot(x, y, label=ltext, color=color, linewidth=2)
            else:
                pylab.plot(x, y, label=ltext, linewidth=2)
        else:
            if color is not None:
                pylab.semilogx(x, y, label=ltext, color=color, linewidth=2)
            else:
                pylab.semilogx(x, y, label=ltext, linewidth=2)
    lfont = pylab.matplotlib.font_manager.FontProperties(size=legend_font_size)
    pylab.legend(loc=legend_location.replace('_',' '), prop=lfont)
    rfn = os.path.join(os.path.dirname(__file__), RANDOM_FILE)
    x, y, auc = read_points(rfn)
    if normal:
        pylab.plot(x, y, 'k--', label='Random')
        pylab.axis([0, 100, 0, 100])
    else:
        pylab.semilogx(x, y, 'k--', label='Random')
        pylab.axis([0.1, 100, 0, 100])
    prefix = ""
    if title:
        pylab.title(title)
    elif guess:
        head = os.path.dirname(os.path.abspath(indirs[0]))
        tail = os.path.basename(head)
        pylab.title("%s %s Plot%s" % (tail.upper(), plottype, dbtype))
        prefix = tail.lower() + "_"
    elif target:
        pylab.title("%s %s Plot%s" % (target, plottype, dbtype))
        prefix = target.lower() + "_"
    else:
        pylab.title("%s Plot%s" % (plottype, dbtype))
    if xlabel:
        pylab.xlabel(xlabel)
    else:
        pylab.xlabel("%% %s" % searchtype)
    if ylabel:
        pylab.ylabel(ylabel)
    else:
        pylab.ylabel("% Ligands Found")
    plotname = prefix + os.path.splitext(filename)[0] + '.png'
    figure = pylab.gcf()
    figure.set_size_inches(8., 8.) #square figures
    figure.savefig(os.path.join(outdir, plotname))
    pylab.clf()

def gen_plots(indirs, labels=[], colors=[], target=None, outdir=".",
        guess=False, forceit=False, skip_decoys=False, normal=False,
        enrich_flag=False, hide_auc=False, 
        title=None, xlabel=None, ylabel=None, 
        legend_font_size=14, legend_location='best', only_own=False, \
        ligfile=None, decfile=None):
    """Generate all ROC and Enrichment plots."""
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    enrich.regen_dirs(indirs, forceit=forceit, \
            skip_decoys=skip_decoys, ligfile=ligfile, decfile=decfile)
    if enrich_flag:
        curvetypes = enrich.CURVE_TYPES
    else:
        curvetypes = [enrich.CURVE_TYPES[0]]
    for curvetype in curvetypes:
        if skip_decoys:
            runtypes = [curvetype[0]] #roc only, no decoys
        elif only_own:
            runtypes = [curvetype[1]] #roc_own only
        else:
            runtypes = curvetype #roc & roc_own, database background and decoys
        for filename in runtypes:
            plot_curves(filename, indirs, labels=labels, colors=colors,
                        target=target, outdir=outdir, guess=guess,
                        skip_decoys=skip_decoys, normal=normal, 
                        title=title, xlabel=xlabel, ylabel=ylabel,
                        hide_auc=hide_auc, legend_font_size=legend_font_size,
                        legend_location=legend_location)

def main(argv):
    """Parse arguments."""
    description = "Plot enrichment and ROC curves into images."
    usage = "%prog [options]"
    version = "%prog: version 200902 - created by Michael Mysinger"
    parser = OptionParser(usage=usage, description=description,
                          version=version)
    parser.set_defaults(indir=[], label=[], color=[], target=None, outdir=".", 
                        guess=False, force=False, skip_decoys=False, 
                        normal=False, enrich_flag=False, hide_auc=False, 
                        title=None, xlabel=None, ylabel=None, 
                        legend_font_size=14, legend_location='best', \
                        ligfile=None, decfile=None)
    parser.add_option("-i", "--indir", action="append", 
           help="input directory, with multiple -i options adding additional curves.") 
    parser.add_option("-l", "--label", action="append", 
           help="custom curve label (default: parent of each indir)") 
    parser.add_option("-c", "--color", action="append", dest="color", \
           help="custom color (default: random color scheme)")
    parser.add_option("-o", "--outdir",
           help="output directory (default: %default)")
    parser.add_option("-s", "--skip-decoys", action="store_true",
           help="skip reading decoys and making own decoys plots")
    parser.add_option("-e", "--enrich", action="store_true",
        help="output legacy enrichment plot (note ROC is always superior)")
    parser.add_option("-t", "--target",
           help="use this target name in filename and title text")
    parser.add_option("--custom-title", dest="title",
           help="custom title text")
    parser.add_option("-x", "--x-label", dest="xlabel", 
           help="custom x-axis label")
    parser.add_option("-y", "--y-label", dest="ylabel", 
           help="custom y-axis label")
    parser.add_option("-f", "--force-it", action="store_true",
           dest="force", help="force enrich and combine to run again")
    parser.add_option("-n", "--normal", action="store_true",
           help="use a normal axis and report the normal AUC")
    parser.add_option("-a", "--auc", action="store_true",
           help="hide AUC/LogAUC in legend")
    parser.add_option("--legend-font-size", type="int",
           help="font size for legend text (default: %default)")
    parser.add_option("--legend-location",
           help="location for legend i.e. upper_left, lower_right (default: %default)")
    parser.add_option("-g", "--guess", action="store_true",
            help="guess target from first indir path")
    parser.add_option("--ligand-file", dest='ligfile',
           help="ligand id file (default: pull from DUD)")
    parser.add_option("-d", "--decoy-file", dest='decfile',
           help="decoy id file (default: pull from DUD)")
    options, args = parser.parse_args(args=argv[1:])
    if len(args):
        parser.error("program takes no positional arguments.\n" +
                     "  Use --help for more information.")
    elif not options.indir:
        parser.error("must specify at least one input directory using -i.\n" +
                     "  Use --help for more information.")
    elif options.label and len(options.label) != len(options.indir):
        parser.error("if custom labels are specified, \n" +
                     "there must be one label for each input file.\n" +
                     "  Use --help for more information.")
    gen_plots(indirs=options.indir, labels=options.label, colors=options.color,
        target=options.target, outdir=options.outdir, guess=options.guess,
        forceit=options.force, skip_decoys=options.skip_decoys, 
        normal=options.normal, enrich_flag=options.enrich, 
        hide_auc=options.auc, title=options.title, 
        xlabel=options.xlabel, ylabel=options.ylabel,
        legend_font_size=options.legend_font_size, 
        legend_location=options.legend_location, 
        only_own=True, ligfile=options.ligfile, decfile=options.decfile)
    return 0

if __name__ == "__main__":
    pylab.matplotlib.use('Agg')  # allows you to not have an x-server running
    sys.exit(main(sys.argv))
