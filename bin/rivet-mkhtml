#! /usr/bin/env python

"""\
%(prog)s [options] <yodafile1> [<yodafile2> <yodafile3>...] [PLOT:Key1=Val1:...]

Make web pages from histogram files written out by Rivet.  You can specify
multiple Monte Carlo YODA files to be compared in the same syntax as for
rivet-cmphistos, i.e. including plotting options.

Reference data, analysis metadata, and plot style information should be found
automatically (if not, set the RIVET_ANALYSIS_PATH or similar variables
appropriately).

Any existing output directory will be overwritten.

ENVIRONMENT:
 * RIVET_ANALYSIS_PATH: list of paths to be searched for analysis plugin libraries
 * RIVET_DATA_PATH: list of paths to be searched for data files
"""

from __future__ import print_function

import rivet, sys, os
rivet.util.check_python_version()
rivet.util.set_process_name(os.path.basename(__file__))
COMMAND = " ".join([os.path.basename(sys.argv[0])] + sys.argv[1:])

import glob, shutil
from subprocess import Popen,PIPE


import argparse
parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument("YODAFILES", nargs="+", help="data files to compare")
parser.add_argument("-o", "--outputdir", dest="OUTPUTDIR",
                    default="./rivet-plots", help="directory for Web page output")
parser.add_argument("-c", "--config", dest="CONFIGFILES", action="append", default=["~/.make-plots"],
                    help="plot config file(s) to be used with rivet-cmphistos")
parser.add_argument("-n", "--num-threads", metavar="NUMTHREADS", dest="NUMTHREADS", type=int,
                    default=None, help="request make-plots to use a specific number of threads")
parser.add_argument("--ignore-missing", dest="IGNORE_MISSING", action="store_true",
                    default=False, help="ignore missing YODA files")
parser.add_argument("-i", "--ignore-unvalidated", dest="IGNORE_UNVALIDATED", action="store_true",
                    default=False, help="ignore unvalidated analyses")
# parser.add_argument("--ref", "--refid", dest="REF_ID",
#                   default=None, help="ID of reference data set (file path for non-REF data)")
parser.add_argument("--no-rivet-refs", dest="RIVETREFS", action="store_false",
                        default=True, help="don't use Rivet reference data files")
parser.add_argument("--dry-run", help="don't actually do any plotting or HTML building", dest="DRY_RUN",
                    action="store_true", default=False)
parser.add_argument("--no-cleanup", dest="NO_CLEANUP", action="store_true", default=False,
                    help="keep plotting temporary directory")
parser.add_argument("--no-subproc", dest="NO_SUBPROC", action="store_true", default=False,
                    help="don't use subprocesses to render the plots in parallel -- useful for debugging")
parser.add_argument("--pwd", dest="PATH_PWD", action="store_true", default=False,
                    help="append the current directory (pwd) to the analysis/data search paths (cf. $RIVET_ANALYSIS_PATH)")

stygroup = parser.add_argument_group("Style options")
stygroup.add_argument("-t", "--title", dest="TITLE",
                      default="Plots from Rivet analyses", help="title to be displayed on the main web page")
stygroup.add_argument("--reftitle", dest="REFTITLE",
                      default="Data", help="legend entry for reference data")
stygroup.add_argument("--no-plottitle", dest="NOPLOTTITLE", action="store_true",
                      default=False, help="don't show the plot title on the plot "
                      "(useful when the plot description should only be given in a caption)")
stygroup.add_argument("-s", "--single", dest="SINGLE", action="store_true",
                      default=False, help="display plots on single webpage.")
stygroup.add_argument("--no-ratio", dest="SHOW_RATIO", action="store_false",
                      default=True, help="don't draw a ratio plot under each main plot.")
stygroup.add_argument("--errs", "--mcerrs", "--mc-errs", dest="MC_ERRS", action="store_true",
                      default=False, help="plot error bars.")
stygroup.add_argument("--offline", dest="OFFLINE", action="store_true",
                      default=False, help="generate HTML that does not use external URLs.")
stygroup.add_argument("--pdf", dest="VECTORFORMAT", action="store_const", const="PDF",
                      default="PDF", help="use PDF as the vector plot format.")
stygroup.add_argument("--ps", dest="VECTORFORMAT", action="store_const", const="PS",
                      default="PDF", help="use PostScript as the vector plot format. DEPRECATED")
stygroup.add_argument("--booklet", dest="BOOKLET", action="store_true",
                      default=False, help="create booklet (currently only available for PDF with pdftk or pdfmerge).")
stygroup.add_argument("--font", dest="OUTPUT_FONT", choices="palatino,cm,times,helvetica,minion".split(","),
                      default="palatino", help="choose the font to be used in the plots")
stygroup.add_argument("--palatino", dest="OUTPUT_FONT", action="store_const", const="palatino", default="palatino",
                      help="use Palatino as font (default). DEPRECATED: Use --font")
stygroup.add_argument("--cm", dest="OUTPUT_FONT", action="store_const", const="cm", default="palatino",
                      help="use Computer Modern as font. DEPRECATED: Use --font")
stygroup.add_argument("--times", dest="OUTPUT_FONT", action="store_const", const="times", default="palatino",
                      help="use Times as font. DEPRECATED: Use --font")
stygroup.add_argument("--helvetica", dest="OUTPUT_FONT", action="store_const", const="helvetica", default="palatino",
                      help="use Helvetica as font. DEPRECATED: Use --font")
stygroup.add_argument("--minion", dest="OUTPUT_FONT", action="store_const", const="minion", default="palatino",
                      help="use Adobe Minion Pro as font. DEPRECATED: Use --font")
stygroup.add_argument("--remove-options", help="remove options label from legend", dest="REMOVE_OPTIONS",
                    action="store_true", default=False)

selgroup = parser.add_argument_group("Selective plotting")
selgroup.add_argument("-m", "--match", action="append", dest="PATHPATTERNS", default=[],
                      help="only write out histograms whose $path/$name string matches any of these regexes")
selgroup.add_argument("-M", "--unmatch", action="append", dest="PATHUNPATTERNS", default=[],
                      help="exclude histograms whose $path/$name string matches any of these regexes")
selgroup.add_argument("--ana-match", action="append", dest="ANAPATTERNS", default=[],
                      help="only write out histograms from analyses whose name matches any of these regexes")
selgroup.add_argument("--ana-unmatch", action="append", dest="ANAUNPATTERNS", default=[],
                      help="exclude histograms from analyses whose name matches any of these regexes")
selgroup.add_argument("--no-weights", help="prevent multiweights from being plotted", dest="NO_WEIGHTS",
                    action="store_true", default=False)

vrbgroup = parser.add_argument_group("Verbosity control")
vrbgroup.add_argument("-v", "--verbose", help="add extra debug messages", dest="VERBOSE",
                      action="store_true", default=False)

args = parser.parse_args()
yodafiles = args.YODAFILES

## Add pwd to search paths
if args.PATH_PWD:
    rivet.addAnalysisLibPath(os.path.abspath("."))
    rivet.addAnalysisDataPath(os.path.abspath("."))


## Check that there are some arguments!
if not yodafiles:
    print("Error: You need to specify some YODA files to be plotted!")
    sys.exit(1)


## Make output directory
if not args.DRY_RUN:
    if os.path.exists(args.OUTPUTDIR) and not os.path.realpath(args.OUTPUTDIR)==os.getcwd():
        import shutil
        shutil.rmtree(args.OUTPUTDIR)
    try:
        os.makedirs(args.OUTPUTDIR)
    except:
        print("Error: failed to make new directory '%s'" % args.OUTPUTDIR)
        sys.exit(1)

## Get set of analyses involved in the runs
plotarg = None
analyses = set()
blocked_analyses = set()
import yoda
for yodafile in yodafiles:
    if yodafile.startswith("PLOT:"):
        plotarg = yodafile
        continue
    yodafilepath = os.path.abspath(yodafile.split(":")[0])
    if not os.access(yodafilepath, os.R_OK):
        print("Error: cannot read from %s" % yodafilepath)
        if args.IGNORE_MISSING:
            continue
        else:
            sys.exit(2)

    try:
        ## Note: we use -m/-M flags here as well as when calling rivet-cmphistos, to potentially speed this initial loading
        analysisobjects = yoda.read(yodafilepath, patterns=args.PATHPATTERNS, unpatterns=args.PATHUNPATTERNS)
    except IOError as e:
        print("File reading error: ", e.strerror)
        sys.exit(1)

    for path, ao in analysisobjects.items():
        ## Make a path object and ensure the path is in standard form
        try:
            aop = rivet.AOPath(path)
        except Exception as e:
            #print(e)
            print("Found analysis object with non-standard path structure:", path, "... skipping")
            continue

        ## We don't plot data objects with path components hidden by an underscore prefix
        if aop.istmp() or aop.israw():
            continue

        ## Identify analysis/histo name parts
        analysis = "ANALYSIS"
        if aop.basepathparts(keepref=False):
            analysis = rivet.stripOptions(aop.basepathparts(keepref=False)[0]) #< TODO: for compatibility with rivet-cmphistos... generalise?
            #analysis = "_".join(aop.dirnameparts(keepref=False)[:-1]) #< TODO: would this be nicer? Currently incompatible with rivet-cmphistos

        ## Optionally veto on analysis name pattern matching
        if analysis in analyses.union(blocked_analyses):
            continue
        import re
        matched = True
        if args.ANAPATTERNS:
            matched = False
            for patt in args.ANAPATTERNS:
                if re.match(patt, analysis) is not None:
                    matched = True
                    break
        if matched and args.ANAUNPATTERNS:
            for patt in args.ANAUNPATTERNS:
                if re.match(patt, analysis) is not None:
                    matched = False
                    break
        if matched:
            analyses.add(analysis)
        else:
            blocked_analyses.add(analysis)


## Sort analyses: group ascending by analysis name (could specialise grouping by collider), then
## descending by year, and finally descending by bibliographic archive ID code (INSPIRE first).
def anasort(name):
    rtn = (1, name)
    if name.startswith("MC"):
        rtn = (99999999, name)
    else:
        stdparts = name.split("_")
        try:
            year = int(stdparts[1])
            rtn = (0, stdparts[0], -year, 0)
            idcode = (0 if stdparts[2][0] == "I" else 1e10) - int(stdparts[2][1:])
            rtn = (0, stdparts[0], -year, idcode)
            if len(stdparts) > 3:
                rtn += stdparts[3:]
        except:
            pass
    return rtn

analyses = sorted(analyses, key=anasort)

## Uncomment to test analysis ordering on index page
# print(analyses)
# sys.exit(0)


## Run rivet-cmphistos to get plain .dat files from .yoda
## We do this here since it also makes the necessary directories
ch_cmd = ["rivet-cmphistos"]
if args.MC_ERRS:
    ch_cmd.append("--mc-errs")
if not args.SHOW_RATIO:
    ch_cmd.append("--no-ratio")
if not args.RIVETREFS:
    ch_cmd.append("--no-rivet-refs")
if args.NOPLOTTITLE:
    ch_cmd.append("--no-plottitle")
if args.REMOVE_OPTIONS:
    ch_cmd.append("--remove-options")
if args.NO_WEIGHTS:
    ch_cmd.append("--no-weights")
# if args.REF_ID is not None:
#     ch_cmd.append("--refid=%s" % os.path.abspath(args.REF_ID))
if args.REFTITLE:
    ch_cmd.append("--reftitle=%s" % args.REFTITLE )
if args.PATHPATTERNS:
    for patt in args.PATHPATTERNS:
        ch_cmd += ["-m", patt] #"'"+patt+"'"]
if args.PATHUNPATTERNS:
    for patt in args.PATHUNPATTERNS:
        ch_cmd += ["-M", patt] #"'"+patt+"'"]
ch_cmd.append("--hier-out")
# TODO: Need to be able to override this: provide a --plotinfodir cmd line option?
ch_cmd.append("--plotinfodir=%s" % os.path.abspath("../"))
for af in yodafiles:
    yodafilepath = os.path.abspath(af.split(":")[0])
    if af.startswith("PLOT:"):
        yodafilepath = "PLOT"
    elif not os.access(yodafilepath, os.R_OK):
        continue
    newarg = yodafilepath
    if ":" in af:
        newarg += ":" + af.split(":", 1)[1]
    # print(newarg)
    ch_cmd.append(newarg)

## Pass rivet-mkhtml -c args to rivet-cmphistos
for configfile in args.CONFIGFILES:
    configfile = os.path.abspath(os.path.expanduser(configfile))
    if os.access(configfile, os.R_OK):
        ch_cmd += ["-c", configfile]

if args.VERBOSE:
    ch_cmd.append("--verbose")
    print("Calling rivet-cmphistos with the following command:")
    print(" ".join(ch_cmd))

## Run rivet-cmphistos in a subdir, after fixing any relative paths in Rivet env vars
if not args.DRY_RUN:
    for var in ("RIVET_ANALYSIS_PATH", "RIVET_DATA_PATH", "RIVET_REF_PATH", "RIVET_INFO_PATH", "RIVET_PLOT_PATH"):
        if var in os.environ:
            abspaths = [os.path.abspath(p) for p in os.environ[var].split(":")]
            os.environ[var] = ":".join(abspaths)
    subproc = Popen(ch_cmd, cwd=args.OUTPUTDIR, stdout=PIPE, stderr=PIPE)
    out, err = subproc.communicate()
    retcode = subproc.returncode
    if args.VERBOSE or retcode != 0:
        print('Output from rivet-cmphistos:\n', out.decode("utf-8"))
    if err :
        print('Errors from rivet-cmphistos:\n', err.decode("utf-8"))
    if retcode != 0:
        print('Crash in rivet-cmphistos code = ', retcode, ' exiting')
        exit(retcode)



## Write web page containing all (matched) plots
## Make web pages first so that we can load it locally in
## a browser to view the output before all plots are made
if not args.DRY_RUN:

    style = '''<style>
      html { font-family: sans-serif; }
      img { border: 0; }
      a { text-decoration: none; font-weight: bold; }
    </style>
    '''

    ## Include MathJax configuration
    script = ''
    if not args.OFFLINE:
        # WAS: src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"
        script = '''\
        <script type="text/x-mathjax-config">
        MathJax.Hub.Config({
          tex2jax: {inlineMath: [["$","$"]]}
        });
        </script>
        <script type="text/javascript" async
          src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
        </script>
        '''
    ## Add local JavaScript
    script += "\n\n" + '''\
        <script type="text/javascript">
            function filterPlots(ana, patt) {
                var regex = new RegExp(patt);
                var sec = document.getElementById('ana_'+ana);
                var plots = sec.getElementsByClassName(\'plot\');
                var i; for (i = 0; i < plots.length; ++i) {
                    // var viz = (regex.test(plots[i].id)) ? 'visible' : 'hidden';
                    // plots[i].style.visibility = viz;
                    var dis = (regex.test(plots[i].id)) ? 'block' : 'none';
                    plots[i].style.display = dis;
                }
                return false;
           }
        </script>
        '''
    ## A helper function for metadata LaTeX -> HTML conversion
    from rivet.util import htmlify

    ## Timestamp and command HTML fragments to be used on each page:
    import datetime
    timestamp = '<p>Generated at %s</p>' % datetime.datetime.now().strftime("%A, %d. %B %Y %I:%M%p")
    command = '<p>Created with command: <pre>%s</pre></p>' % COMMAND

    index = open(os.path.join(args.OUTPUTDIR, "index.html"), "w")
    index.write('<html>\n<head>\n<title>%s</title>\n%s</head>\n<body>' % (args.TITLE, style + script))
    if args.BOOKLET and args.VECTORFORMAT == "PDF":
        index.write('<h2><a href="booklet.pdf">%s</a></h2>\n\n' % args.TITLE)
    else:
        index.write('<h2>%s</h2>\n\n' % args.TITLE)

    if args.SINGLE and len(analyses) > 1:
        ## Write table of contents
        index.write('<ul>\n')
        for analysis in analyses:
            summary = analysis
            ana = rivet.AnalysisLoader.getAnalysis(analysis)
            if ana:
                summary = "%s (%s)" % (ana.summary(), analysis)
                if args.IGNORE_UNVALIDATED and ana.status() != "VALIDATED":
                    continue
            index.write('<li><a href="#%s">%s</a>\n' % (analysis, htmlify(summary)) )
        index.write('</ul>\n')

    official_routines = rivet.stdAnalysisNames()
    for analysis in analyses:
        references = []
        summary = htmlify(analysis)
        description, inspireid, spiresid = None, None, None

        if analysis.find("_I") > 0:
            inspireid = analysis[analysis.rfind('_I')+2:len(analysis)]
        elif analysis.find("_S") > 0:
            spiresid = analysis[analysis.rfind('_S')+2:len(analysis)]

        ana = rivet.AnalysisLoader.getAnalysis(analysis)
        if ana:
            if ana.summary():
                summary = htmlify("%s (%s)" % (ana.summary(), analysis))
            references = ana.references()
            description = htmlify(ana.description())
            spiresid = ana.spiresId()
            if args.IGNORE_UNVALIDATED and ana.status().upper() != "VALIDATED":
                continue

        try:
            if args.SINGLE:
                index.write('<section id="ana_{}">\n'.format(analysis))
                index.write('\n<h3 style="clear:left; padding-top:2em;"><a name="%s">%s</a></h3>\n' % (analysis, summary))
            else:
                index.write('\n<h3><a href="%s/index.html" style="text-decoration:none;">%s</a></h3>\n' % (analysis, summary))
        except UnicodeEncodeError as ue:
            print("Unicode error in analysis description for " + analysis + ": " + str(ue))

        reflist = []
        if inspireid:
            reflist.append('<a href="https://inspirehep.net/literature/%s">Inspire</a>' % inspireid)
            reflist.append('<a href="http://hepdata.net/record/ins%s">HepData</a>' % inspireid)
        elif spiresid:
        # elif ana.spiresId():
            reflist.append('<a href="https://inspirehep.net/literature?q=%s">Inspire</a>' % spiresid)
            reflist.append('<a href="http://hepdata.cedar.ac.uk/view/irn%s">HepData</a>' % spiresid)
        if analysis in official_routines:
            reflist.append('<a href="https://rivet.hepforge.org/analyses/%s.html">Analysis reference</a>' % analysis)
        for i, ref in enumerate(references):
            if "arXiv" in ref:
                try:
                    arxivID = ref.split(' ')[0].split(':')[-1]
                    if arxivID:
                        references[i] = '<a href="https://arxiv.org/abs/%s">arXiv:%s</a>' % (arxivID, arxivID)
                except:
                    pass
        reflist += references
        index.write('<p>%s</p>\n' % " &#124; ".join(reflist))

        if description:
            try:
                index.write('<p style="font-size:smaller;">%s</p>\n' % description)
            except UnicodeEncodeError as ue:
                print("Unicode error in analysis description for " + analysis + ": " + str(ue))

        anapath = os.path.join(args.OUTPUTDIR, analysis)
        if not args.SINGLE:
            if not os.path.exists(anapath):
                os.makedirs(anapath)
            anaindex = open(os.path.join(anapath, "index.html"), 'w')
            anaindex.write('<html>\n<head>\n<title>%s &ndash; %s</title>\n%s</head>\n<body>\n' %
                           (htmlify(args.TITLE), analysis, style + script))
            anaindex.write('<section id="ana_{}">\n'.format(analysis))
            anaindex.write('<h3>%s</h3>\n' % htmlify(analysis))
            anaindex.write('<p><a href="../index.html">Back to index</a></p>\n')
            if description:
                try:
                    anaindex.write('<p>\n  %s\n</p>\n' % description)
                except UnicodeEncodeError as ue:
                    print("Unicode error in analysis description for " + analysis + ": " + str(ue))
        else:
            anaindex = index

        ## JS filtering
        anaindex.write('<form onsubmit="var patt = document.getElementById(\'patt_{ana}\').value; filterPlots(\'{ana}\', patt); return false;">\n'.format(ana=analysis))
        anaindex.write('<span>Filter plots:&nbsp;</span>')
        anaindex.write('<input id="patt_{ana}" type="search">\n'.format(ana=analysis))
        anaindex.write('<input type="submit" value="Filter">\n')
        # anaindex.write('<span style="background-color:00cc00; border:2px solid #009933; padding:2px; border-radius:2px;"\n onclick="var patt = document.getElementById(\'patt_{ana}\').value; filterPlots(\'{ana}\', patt);">Filter</span>\n'.format(ana=analysis))
        anaindex.write('</form>\n\n')

        datfiles = glob.glob("%s/*.dat" % anapath)
        #print(datfiles)
        # anaindex.write('<div style="float:none; overflow:auto; width:100%">\n')
        for datfile in sorted(datfiles):
            obsname = os.path.basename(datfile).replace(".dat", "")
            pngfile = obsname+".png"
            vecfile = obsname+"."+args.VECTORFORMAT.lower()
            srcfile = obsname+".dat"
            if args.SINGLE:
                pngfile = os.path.join(analysis, pngfile)
                vecfile = os.path.join(analysis, vecfile)
                srcfile = os.path.join(analysis, srcfile)

            anaindex.write('  <div class="plot" id="plot_{o}" style="float:left; font-size:smaller; font-weight:bold;">\n'.format(o=obsname))
            anaindex.write('    <a href="#%s-%s">&#9875;</a><a href="%s">&#8984</a> %s:<br/>\n' %
                           (analysis, obsname, srcfile, os.path.splitext(vecfile)[0]) )
            anaindex.write('    <a name="%s-%s"><a href="%s">\n' % (analysis, obsname, vecfile) )
            anaindex.write('      <img src="%s">\n' % pngfile )
            anaindex.write('    </a></a>\n')
            anaindex.write('  </div>\n')
        # anaindex.write('</div>\n')

        if args.SINGLE:
            index.write('</section>\n\n')
        else:
            anaindex.write('</section>\n\n')
            anaindex.write('<footer style="clear:both; margin-top:3em; padding-top:3em">\n%s\n</footer>\n\n' % timestamp)
            anaindex.write('</body>\n</html>')
            anaindex.close()

    index.write('<footer style="clear:both; margin-top:3em; padding-top:3em">\n%s\n%s\n</footer>\n\n' % (timestamp, command))
    index.write('</body>\n</html>')
    index.close()

# http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

## Run make-plots on all generated .dat files
# sys.exit(0)
mp_cmd = ["make-plots"]
for configfile in args.CONFIGFILES:
    configfile = os.path.abspath(os.path.expanduser(configfile))
    if os.access(configfile, os.R_OK):
        mp_cmd.append("--config=%s" % configfile)
if args.NUMTHREADS:
    mp_cmd.append("--num-threads=%d" % args.NUMTHREADS)
if args.NO_CLEANUP:
    mp_cmd.append("--no-cleanup")
if args.NO_SUBPROC:
    mp_cmd.append("--no-subproc")
if args.VECTORFORMAT == "PDF":
    mp_cmd.append("--pdfpng")
elif args.VECTORFORMAT == "PS":
    mp_cmd.append("--pspng")
if args.OUTPUT_FONT:
    mp_cmd.append("--font=%s" % args.OUTPUT_FONT)
# if args.OUTPUT_FONT.upper() == "PALATINO":
#     mp_cmd.append("--palatino")
# if args.OUTPUT_FONT.upper() == "CM":
#     mp_cmd.append("--cm")
# elif args.OUTPUT_FONT.upper() == "TIMES":
#     mp_cmd.append("--times")
# elif args.OUTPUT_FONT.upper() == "HELVETICA":
#     mp_cmd.append("--helvetica")
# elif args.OUTPUT_FONT.upper() == "MINION":
#     mp_cmd.append("--minion")
datfiles = []
for analysis in analyses:
    anapath = os.path.join(args.OUTPUTDIR, analysis)
    #print(anapath)
    anadatfiles = glob.glob("%s/*.dat" % anapath)
    datfiles += sorted(anadatfiles)
if datfiles:
    mp_cmd += datfiles
    if args.VERBOSE:
        mp_cmd.append("--verbose")
        print("Calling make-plots with the following options:")
        print(" ".join(mp_cmd))
    if not args.DRY_RUN:
        Popen(mp_cmd).wait()
        if args.BOOKLET and args.VECTORFORMAT=="PDF":
            if which("pdftk") is not None:
                bookletcmd = ["pdftk"]
                for analysis in analyses:
                    anapath = os.path.join(args.OUTPUTDIR, analysis)
                    bookletcmd += sorted(glob.glob("%s/*.pdf" % anapath))
                bookletcmd += ["cat", "output", "%s/booklet.pdf" % args.OUTPUTDIR]
                print(bookletcmd)
                Popen(bookletcmd).wait()
            elif which("pdfmerge") is not None:
                bookletcmd = ["pdfmerge"]
                for analysis in analyses:
                    anapath = os.path.join(args.OUTPUTDIR, analysis)
                    bookletcmd += sorted(glob.glob("%s/*.pdf" % anapath))
                bookletcmd += ["%s/booklet.pdf" % args.OUTPUTDIR]
                print(bookletcmd)
                Popen(bookletcmd).wait()
            else:
                print("Neither pdftk nor pdfmerge available --- not booklet output possible")
