#! /usr/bin/env python

## Get output filename
OUTNAME = "selectedanalyses"
import sys
if len(sys.argv) < 2:
    pass
    #print "Using output name '%s'" % OUTNAME
else:
    OUTNAME = sys.argv[1]


## Get input paths to allow rivet module to be imported from the src dir
import os, re, glob
pybuild = os.path.abspath(os.path.join(os.getcwd(), "..", "pyext", "build"))
dirs = []
for d in os.listdir(pybuild):
    if re.match(r"lib\..*-.*-%d\.%d" % (sys.version_info[0], sys.version_info[1]), d):
        dirs.append(os.path.join(pybuild, d))
sys.path = dirs + sys.path
try:
    os.environ["LD_LIBRARY_PATH"] = os.environ["LD_LIBRARY_PATH"] + ":" + \
        os.path.abspath(os.path.join(os.getcwd(), "..", "src", ".libs"))
except:
    pass
try:
    os.environ["DYLD_LIBRARY_PATH"] = os.environ["DYLD_LIBRARY_PATH"] + ":" + \
        os.path.abspath(os.path.join(os.getcwd(), "..", "src", ".libs"))
except:
    pass
anadirs = glob.glob(os.path.join(os.getcwd(), "..", "src", "Analyses", ".libs"))
#print anadirs
os.environ["RIVET_ANALYSIS_PATH"] = ":".join(anadirs)


## Change dlopen status to GLOBAL for Rivet lib
try:
    import ctypes
    sys.setdlopenflags(sys.getdlopenflags() | ctypes.RTLD_GLOBAL)
except:
    import dl
    sys.setdlopenflags(sys.getdlopenflags() | dl.RTLD_GLOBAL)
import rivet


def texify(s):
    t = s \
        .replace(r"&", r"\&") \
        .replace(r"\\&", r"\&") \
        .replace(r"#", r"\#") \
        .replace(r"->", r"\ensuremath{\to}") \
        .replace(r"pT", r"\pT") \
        .replace(r"sqrt(s)", r"\ensuremath{\sqrt{s}}")
        # .replace(r"_", r"\_") \
        # .replace(r"^", r"") \
    return t

## Build analysis pages
analyses = ["ALEPH_1996_S3196992",
            "DELPHI_1996_S3430090",
            "OPAL_2004_S6132243",
            "SLD_2004_S5693039",
            "CDF_2001_S4751469",
            "D0_2008_S7719523",
            "ALICE_2011_S8945144",
            "ATLAS_2012_I1094568",
            "CMS_2011_S8957746",
            "LHCB_2011_I919315",
            "LHCF_2012_I1115479",
            "TOTEM_2012_I1115294",
            "UA1_1990_S2044935",
            "UA5_1982_S875503",
            "H1_2000_S4129130",
            "STAR_2006_S6500200",
            "ARGUS_1993_S2653028",
            "BABAR_2007_S7266081",
            "BELLE_2006_S6265367",
            "CLEO_2004_S5809304",
            "JADE_1998_S3612880",
            "PDG_HADRON_MULTIPLICITIES",
            "TASSO_1990_S2148048" ]
pages = []
## Use list(...) ctor for 2.3 compatibility
bib = {}
for aname in sorted(list(analyses)):
    #print "Handling analysis '%s'" % aname
    page = ""
    safe_aname = aname.replace(r"_", r"\_")
    ana = rivet.AnalysisLoader.getAnalysis(aname)
    subtitle = "\\subsection{%s:\\\\ %s}\n" % (safe_aname, ana.summary())
    if ana.bibKey() and ana.bibTeX():
        bib[ana.bibKey()] = "%% (%s)\n" % aname + ana.bibTeX()
        citetex = r"\cite{%s}" % ana.bibKey()
        subtitle = "\\subsection[%s]{%s:\\\\ %s}\n" % (safe_aname, safe_aname + "\," + citetex, ana.summary())
    page += subtitle + "\n"

    for para in ana.description().split("\n\n"):
        page += "\n\\noindent " + para + "\n\n"

    if ana.requiredBeams():
        def pid_to_str(pid):
            if pid == 11:
                return "$e^-$"
            elif pid == -11:
                return "$e^+$"
            elif pid == 2212:
                return "$p$"
            elif pid == -2212:
                return "$\\bar{p}$"
            elif pid == 10000:
                return "$*$"
            else:
                return str(pid)
        beamstrs = []
        for bp in ana.requiredBeams():
            beamstrs.append(pid_to_str(bp[0]) + "\\," + pid_to_str(bp[1]))
        page += "\\noindent\\textsc{Beams:} %s \\newline\n" % ", ".join(beamstrs)
    if ana.requiredEnergies():
        page += "\\textsc{Energies:} %s GeV \\newline\n" % \
            ", ".join(["(%0.1f, %0.1f)" % (epair[0], epair[1]) for epair in ana.requiredEnergies()])
    if ana.experiment():
        page += "\\textsc{Experiment:} %s" % ana.experiment()
        if ana.collider():
            page += " (%s)" % ana.collider()
        page += "\\newline\n"
    if ana.inspireId():
        spiresbase = "http://inspire-hep.net/record"
        page += "\\textsc{Inspire ID:} \\href{%s+%s}{%s}\\newline\n" % \
            (spiresbase, ana.inspireId(), ana.inspireId())
    elif ana.spiresId():
        spiresbase = "http://inspire-hep.net/search?p=find+key"
        page += "\\textsc{Spires ID:} \\href{%s+%s}{%s}\\newline\n" % \
            (spiresbase, ana.spiresId(), ana.spiresId())
    page += "\\textsc{Status:} %s\\newline\n" % ana.status()

    if ana.authors():
        page += "\\textsc{Authors:}\n \\penalty 100\n"
        page += "\\begin{itemize}\n"
        for a in ana.authors():
            s = a
            import re
            if re.search(".* <.*@.*>", a):
                name = " ".join(a.split()[:-1])
                email = a.split()[-1].replace("<", "").replace(">", "")
                #s = "\\href{mailto:%s}{%s}" % (email, name)
                s = "%s $\\langle\,$\\href{mailto:%s}{%s}$\,\\rangle$" % (name, email, email)
            page += "  \\item %s\n" % s
        page += "\\end{itemize}\n"
    else:
        page += "\\textsc{No authors listed}\\\\ \n"


    if False and ana.references():
        page += "\\textsc{References:}\n \\penalty 100\n"
        page += "\\begin{itemize}\n"
        for r in ana.references():
            if r.startswith("arXiv:"):
                code = r.split()[0].replace("arXiv:", "")
                url = "http://arxiv.org/abs/" + code
                page += "  \\item %s \\href{%s}{%s}\n" % ("arXiv:", url, code)
            elif r.startswith("doi:"):
                code = r.replace("doi:", "")
                url = "http://dx.doi.org/" + code
                page += "  \\item %s \\href{%s}{%s}\n" % ("DOI:", url, code)
            else:
                page += "  \\item %s\n" % r
        page += "\\end{itemize}\n"


    if ana.runInfo():
        page += "\\textsc{Run details:}\n \\penalty 100\n"
        infos = ana.runInfo().split(" * ")
        #print ana.runInfo(), "->", infos
        page += "\\begin{itemize}\n"
        for i in infos:
            if i:
                page += "\n  \\item %s" % i
        page += "\\end{itemize}\n"
    else:
        page += "\\textsc{No run details listed}\\\\ \n"

    page += "\n"

    page = texify(page)

    pages.append(page)


## Write out LaTeX
prefix = """\
\\makeatletter
\\renewcommand{\\d}[1]{\\ensuremath{\\mathrm{#1}}}
\\let\\old@eta\\eta
\\renewcommand{\\eta}{\\ensuremath{\\old@eta}\\xspace}
\\let\\old@phi\\phi
\\renewcommand{\\phi}{\\ensuremath{\\old@phi}\\xspace}
\\providecommand{\\pT}{\\ensuremath{p_\\perp}\\xspace}
\\providecommand{\\pTmin}{\\ensuremath{p_\\perp^\\text{min}}\\xspace}
\\makeatother

Each Rivet release is accompanied by a standard library of analyses
implementing currently a total of 250 experimental measurements or Monte-Carlo validation
studies. The full listing of these is beyond the scope of this publication, but
it is available both online at \url{http://rivet.hepforge.org/analyses} and as
a part of the manual coming with each release of Rivet in the \kbd{doc/}
sub-directory. Here, we only want to show-case a selection of analyses spanning the full
spectrum of experiments from LEP over HERA to Tevatron and the LHC and
demonstrating the versatility of the Rivet framework.

For each of the 250 analyses, in addition to a brief summary one can find
information about
the collider at which the measurement was made, references to the original
publications, status and authors of the Rivet implementation as well as run
details necessary for comparing a Monte-Carlo prediction with the data.

{\scriptsize
"""

body = ""
for page in pages:
    body = body + page + "\n"

outstr = prefix + body + "}\n"

## Write out to TeX and BibTeX files
f = open("%s.tex" % OUTNAME, "w")
f.write("%auto-ignore\n")
f.write(outstr)
f.close()
f = open("%s.bib" % OUTNAME, "w")
#
bibentries = "\n\n".join(["%% %s\n%s" % (k,b) for k,b in bib.iteritems()])
f.write(bibentries + "\n")
f.close()
