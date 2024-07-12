# cython: embedsignature=True
# cython: c_string_type=str, c_string_encoding=utf8
# distutils: language = c++

import sys

cimport rivet as c
from cython.operator cimport dereference as deref
# Need to be careful with memory management -- perhaps use the base object that
# we used in YODA?

cdef extern from "<utility>" namespace "std" nogil:
    cdef c.unique_ptr[c.Analysis] move(c.unique_ptr[c.Analysis])

# ## Write a string to a file
# ## The file argument can either be a file object, filename, or special "-" reference to stdout
# def _str_to_file(s, file_or_filename):
#     if hasattr(file_or_filename, 'write'):
#         file_or_filename.write(s)
#     elif file_or_filename == "-":
#         sys.stdout.write(s)
#     else:
#         with open(file_or_filename, "w") as f:
#             f.write(s)


cdef void _make_iss(c.istringstream &iss, bytes bs):
    iss.str(bs)


cdef class AnalysisHandler:
    """The main controller class for coordinating analysis setup, execution, and outputs"""

    cdef c.AnalysisHandler *_ptr

    def __cinit__(self):
        self._ptr = new c.AnalysisHandler()

    def __del__(self):
        del self._ptr

    def setIgnoreBeams(self, ignore=True):
        "Choose whether or not to ignore analyses' declared beam requirements [default=False]"
        self._ptr.setIgnoreBeams(ignore)

    def skipMultiWeights(self, ignore=True):
        "Choose whether to ignore weight streams in addition to the nominal [default=False]"
        self._ptr.skipMultiWeights(ignore)

    def selectMultiWeights(self, patterns=""):
        "Choose a subset of variation-weight stream names to consider, by regex pattern"
        self._ptr.selectMultiWeights(patterns)

    def deselectMultiWeights(self, patterns=""):
        "Choose a subset of variation-weight stream names to NOT consider, by regex pattern"
        self._ptr.deselectMultiWeights(patterns)

    def setNominalWeightName(self, name=""):
        "Declare which weight-stream name to treat as the nominal [default=Nominal|Default||0]"
        self._ptr.setNominalWeightName(name)

    def setWeightCap(self, double maxWeight):
        "Set a maximum absolute weight value to use in events with anomalously high weights"
        self._ptr.setWeightCap(maxWeight)

    def setNLOSmearing(self, double smear):
        "Enable the histogram-fill smearing mechanism to tame NLO counter-event non cancellations"
        self._ptr.setNLOSmearing(smear)

    def addAnalysis(self, name):
        "Register an analysis for execution to the AH collection, by name"
        self._ptr.addAnalysis(name)
        return self

    def analysisNames(self):
        "Get the list of registered analyses on this AH"
        anames = self._ptr.analysisNames()
        return [ a for a in anames ]

    def stdAnalysisNames(self):
        "Get the list of registered analyses on this AH, by official names only (no aliases)"
        anames = self._ptr.stdAnalysisNames()
        return [ a for a in anames ]

    # def analysis(self, aname):
    #     cdef c.Analysis* ptr = self._ptr.analysis(aname)
    #     cdef Analysis pyobj = Analysis.__new__(Analysis)
    #     if not ptr:
    #         return None
    #     pyobj._ptr = ptr
    #     return pyobj

    def readData(self, name_or_stream, fmt="yoda", preload=True):
        "Preload histogram data from the provided file name or handle"
        cdef c.istringstream iss
        if type(name_or_stream) is str:
            self._ptr.readData_FILE(name_or_stream, preload)
        else:
            _make_iss(iss, name_or_stream)
            self._ptr.readData_ISTR(iss, fmt, preload)

    def writeData(self, file_or_filename, fmt="yoda"):
        "Write histogram data to the provided file name or handle"
        cdef c.ostringstream oss
        if type(file_or_filename) is str:
            self._ptr.writeData_FILE(file_or_filename)
        else:
            self._ptr.writeData_OSTR(oss, fmt)
            file_or_filename.write(oss.str())

    def nominalCrossSection(self):
        "Get the current nominal cross-section value from the ongoing event run"
        return self._ptr.nominalCrossSection()

    def finalize(self):
        "Perform the finalising operations on all registered analyses"
        self._ptr.finalize()

    def dump(self, name, period):
        "Declare to dump the current status of this AH's histograms to file every <period> events"
        self._ptr.dump(name, period)

    def mergeYodas(self, filelist, delopts, addopts, matches, unmatches, equiv):
        "Access to the API call for merging multiple YODA files correctly, including finalization. Mainly for rivet-merge"
        filelist  = [ f for f in filelist ]
        delopts   = [ d for d in delopts  ]
        addopts   = [ d for d in addopts ]
        matches   = [ d for d in matches ]
        unmatches = [ d for d in unmatches ]
        self._ptr.mergeYodas(filelist, delopts, addopts, matches, unmatches, equiv)

    def merge(self, AnalysisHandler other):
        "Combine analysis data in-memory with another AH object"
        self._ptr.merge(other._ptr[0])


cdef class Run:
    """Object for coordinating a run of events fed to the analyses"""
    cdef c.Run *_ptr

    def __cinit__(self, AnalysisHandler h):
        self._ptr = new c.Run(h._ptr[0])

    def __del__(self):
        del self._ptr

    def setCrossSection(self, double x):
        """\
        Manually set the generated cross-section corresponding to the incoming events

        Overrides any cross-sections provided within the event stream itself.
        """
        self._ptr.setCrossSection(x)
        return self

    def setListAnalyses(self, choice):
        "Tell the analysis handler to list the enabled analyses"
        self._ptr.setListAnalyses(choice)
        return self

    def init(self, name, weight=1.0):
        "Call the init() step on the AnalysisHandler and its analyses"
        return self._ptr.init(name, weight)

    def openFile(self, name, weight=1.0):
        "Open a new event file, with an optional file-level multiplicative weight"
        return self._ptr.openFile(name, weight)

    def readEvent(self):
        "Read the next event from the input event stream"
        return self._ptr.readEvent()

    # def skipEvent(self):
    #     return self._ptr.skipEvent()

    def numEvents(self):
        "Return the number of events processed in this run"
        return self._ptr.numEvents()

    def processEvent(self):
        "Call the analyze() step on the AnalysisHandler and its analyses, for the current event"
        return self._ptr.processEvent()

    def finalize(self):
        "Call the finalize() step on the AnalysisHandler and its analyses"
        return self._ptr.finalize()


cdef class Analysis:
    """\
    The logic of a particular Rivet event analysis routine

    Analysis objects cannot be directly instantiated in the Python interface:
    add them to an analysis run via the named AnalysisHandler.addAnalysis()
    method. This class exists as a method for accessing analysis metadata.
    """

    cdef c.unique_ptr[c.Analysis] _ptr

    def __init__(self):
        raise RuntimeError('This class cannot be instantiated')

    def name(self):
        "Get the analysis-routine name"
        return deref(self._ptr).name()

    def requiredBeams(self):
        "Get the beam configuration required by this analysis"
        return deref(self._ptr).requiredBeams()

    def requiredEnergies(self):
        "Get the beam-energy configuration required by this analysis"
        return deref(self._ptr).requiredEnergies()

    def keywords(self):
        "Get the list of physics keywords for this analysis"
        kws = deref(self._ptr).keywords()
        return [ k for k in kws ]

    def validation(self):
        "Get the validation status of this analysis"
        vld = deref(self._ptr).validation()
        return [ k for k in vld ]

    def reentrant(self):
        "Get whether the analysis is re-entrant, i.e. finalize() can be re-run in post-processing"
        return deref(self._ptr).reentrant()

    def authors(self):
        "Get the list of analysis-routine authors"
        auths = deref(self._ptr).authors()
        return [ a for a in auths ]

    def bibKey(self):
        "Get the BibTeX bibliography key for the corresponding experiment paper"
        return deref(self._ptr).bibKey()

    def bibTeX(self):
        "Get the BibTeX bibliography entry for the corresponding experiment paper"
        return deref(self._ptr).bibTeX()

    def references(self):
        "Get the list of bibliography references for this routine"
        refs = deref(self._ptr).references()
        return [ r for r  in refs ]

    def collider(self):
        "Get the name of the collider on which the corresponding experimental analysis was performed"
        return deref(self._ptr).collider()

    def summary(self):
        "Get a short, one-line description of the analysis"
        return deref(self._ptr).summary()

    def description(self):
        "Get a long description of the analysis methods and context. Often the experimental-paper abstract"
        return deref(self._ptr).description()

    def experiment(self):
        "Get the name of the experiment by which the original analysis was performed"
        return deref(self._ptr).experiment()

    def inspireId(self):
        "Get the Inspire-HEP ID code of the original paper"
        return deref(self._ptr).inspireId()

    def spiresId(self):
        "Get the SPIRES ID code of the original paper [deprecated]"
        return deref(self._ptr).spiresId()

    def runInfo(self):
        "Get information about the MC run conditions required to use this analysis"
        return deref(self._ptr).runInfo()

    def status(self):
        "Get the indicated usability status of this analysis routine"
        return deref(self._ptr).status()

    def warning(self):
        "Get any warning strings indicated for this analysis routine"
        return deref(self._ptr).warning()

    def year(self):
        "Get the year in which the experimental paper was published"
        return deref(self._ptr).year()

    def luminosity(self):
        "Get the corresponding integrated luminosity of the experimental analysis, in picobarns"
        return deref(self._ptr).luminosity()

    def luminosityfb(self):
        "Get the corresponding integrated luminosity of the experimental analysis, in femtobarns"
        return deref(self._ptr).luminosityfb()

    def refMatch(self):
        "A regex for positively filtering matching datasets from the corresponding HepData record"
        return deref(self._ptr).refMatch()

    def refUnmatch(self):
        "A regex for negatively filtering out non-matching datasets from the corresponding HepData record"
        return deref(self._ptr).refUnmatch()

    def writerDoublePrecision(self):
        "Get whether the histogram writer needs to write in double precision for a run containing this analysis"
        return deref(self._ptr).writerDoublePrecision()

    def refFile(self):
        "Get the name of the corresponding reference-data file"
        return deref(self._ptr).refFile()

    def refData(self, asdict=True, patterns=None, unpatterns=None):
        """\
        Get this analysis' reference data, cf. yoda.read()

        NB. There's also a C++ version of this, but this wrapping is nicer for Python.
        """
        import yoda
        return yoda.read(self.refFile(), asdict, patterns, unpatterns)


#cdef object
LEVELS = dict(TRACE = 0, DEBUG = 10, INFO = 20,
              WARN = 30, WARNING = 30, ERROR = 40,
              CRITICAL = 50, ALWAYS = 50)


cdef class AnalysisLoader:
    """\
    Mechanism for finding and loading analyses from Rivet*.so plugin files
    """

    @staticmethod
    def analysisNames():
        "Get the list of available analysis names, not including aliases"
        names = c.AnalysisLoader_analysisNames()
        return [ n for n in names ]

    @staticmethod
    def allAnalysisNames():
        "Get the list of available analysis names, including aliases"
        names = c.AnalysisLoader_allAnalysisNames()
        return [ n for n in names ]

    @staticmethod
    def stdAnalysisNames():
        "Get the list of built-in analysis names"
        names = c.AnalysisLoader_stdAnalysisNames()
        return [ n for n in names ]

    @staticmethod
    def analysisNameAliases():
        "Get the alias -> canonical dict of analysis-name aliases"
        anames = c.AnalysisLoader_analysisNameAliases()
        return { a.first : a.second for a in anames }

    @staticmethod
    def getAnalysis(name):
        "Get a Python wrapper for a named analysis (metadata access only)"
        # try:
        #   name = name
        # except AttributeError:
        #   pass
        cdef c.unique_ptr[c.Analysis] ptr = c.AnalysisLoader_getAnalysis(name)
        cdef Analysis pyobj = Analysis.__new__(Analysis)
        if not ptr:
            return None
        pyobj._ptr = move(ptr)
        # Create python object
        return pyobj

    @staticmethod
    def analysisPlugins():
        "Get the list of active analysis plugin paths"
        ppaths = c.AnalysisLoader_analysisPlugins()
        return [p for p in ppaths]

    # @staticmethod
    # def searchAnalysisPlugins():
    #     "Search for and return the list of analysis plugin paths"
    #     ppaths = c.AnalysisLoader_searchAnalysisPlugins()
    #     return [p for p in ppaths]

    @staticmethod
    def setAnalysisPlugins(ppaths):
        "Set the list of analysis-plugin paths"
        c.AnalysisLoader_setAnalysisPlugins(ppaths)

    @staticmethod
    def loadFromAnalysisPlugins():
        "Load analysis builders from the active analysis plugin libs"
        c.AnalysisLoader_loadFromAnalysisPlugins()


## Convenience versions in main rivet namespace
def analysisNames():
    "Get the list of available analysis names, not including aliases"
    return AnalysisLoader.analysisNames()

def allAnalysisNames():
    "Get the list of available analysis names, including aliases"
    return AnalysisLoader.allAnalysisNames()

def stdAnalysisNames():
    "Get the list of built-in analysis names"
    return AnalysisLoader.stdAnalysisNames()

def analysisNameAliases():
    "Get the list of analysis-name aliases"
    return AnalysisLoader.analysisNameAliases()

def getAnalysis(name):
    "Get a Python wrapper for a named analysis (metadata access only)"
    return AnalysisLoader.getAnalysis(name)

def analysisPlugins():
    "Get the list of active analysis plugin paths"
    return AnalysisLoader.analysisPlugins()

def setAnalysisPlugins(ppaths):
    "Set the list of analysis-plugin paths"
    return AnalysisLoader.setAnalysisPlugins(ppaths)


## Path functions
def getAnalysisLibPaths():
    "Get the list of paths to search for analysis plugin libraries"
    ps = c.getAnalysisLibPaths()
    return [ p for p in ps ]

def setAnalysisLibPaths(xs):
    "Set the list of paths to search for analysis plugin libraries"
    bs = [ x for x in xs ]
    c.setAnalysisLibPaths(bs)

def addAnalysisLibPath(path):
    "Add to the list of paths to search for analysis plugin libraries"
    c.addAnalysisLibPath(path)


def setAnalysisDataPaths(xs):
    "Get the list of paths to search for analysis data files"
    bs = [ x for x in xs ]
    c.setAnalysisDataPaths(bs)

def addAnalysisDataPath(path):
    "Add to the list of paths to search for analysis data files"
    c.addAnalysisDataPath(path)

def getAnalysisDataPaths():
    "Add multiple paths to the list of paths to search for analysis data files"
    ps = c.getAnalysisDataPaths()
    return [ p for p in ps ]

def findAnalysisDataFile(q):
    "Find the first match to a named analysis data file in the search paths"
    f = c.findAnalysisDataFile(q)
    return f

def getAnalysisRefPaths():
    "Get the list of paths to search for analysis reference-data files"
    ps = c.getAnalysisRefPaths()
    return [ p for p in ps ]

def findAnalysisRefFile(q):
    "Find the first match to a named analysis reference-data file in the search paths"
    f = c.findAnalysisRefFile(q)
    return f


def getAnalysisInfoPaths():
    "Get the list of paths to search for analysis info files"
    ps = c.getAnalysisInfoPaths()
    return [ p for p in ps ]

def findAnalysisInfoFile(q):
    "Find the first match to a named analysis info file in the search paths"
    f = c.findAnalysisInfoFile(q)
    return f

def getAnalysisPlotPaths():
    "Get the list of paths to search for analysis plot-style files"
    ps = c.getAnalysisPlotPaths()
    return [ p for p in ps ]

def findAnalysisPlotFile(q):
    "Find the first match to a named analysis plot-style file in the search paths"
    f = c.findAnalysisPlotFile(q)
    return f


def version():
    "Get the Rivet library version"
    return c.version()

def setLogLevel(name, level):
    "Set the log level for a named logger hierarchy"
    c.setLogLevel(name, level)
