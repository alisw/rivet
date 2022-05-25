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
#     s = s.decode('utf-8')
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
    cdef c.AnalysisHandler *_ptr

    def __cinit__(self):
        self._ptr = new c.AnalysisHandler()

    def __del__(self):
        del self._ptr

    def setIgnoreBeams(self, ignore=True):
        self._ptr.setIgnoreBeams(ignore)

    def skipMultiWeights(self, ignore=True):
        self._ptr.skipMultiWeights(ignore)

    def selectMultiWeights(self, patterns=""):
        self._ptr.selectMultiWeights(patterns.encode('utf-8'))

    def deselectMultiWeights(self, patterns=""):
        self._ptr.deselectMultiWeights(patterns.encode('utf-8'))

    def setNominalWeightName(self, name=""):
        self._ptr.setNominalWeightName(name.encode('utf-8'))

    def setWeightCap(self, double maxWeight):
        self._ptr.setWeightCap(maxWeight)

    def setNLOSmearing(self, double smear):
        self._ptr.setNLOSmearing(smear)

    def addAnalysis(self, name):
        self._ptr.addAnalysis(name.encode('utf-8'))
        return self

    def analysisNames(self):
        anames = self._ptr.analysisNames()
        return [ a.decode('utf-8') for a in anames ]

    def stdAnalysisNames(self):
        anames = self._ptr.stdAnalysisNames()
        return [ a.decode('utf-8') for a in anames ]

    # def analysis(self, aname):
    #     cdef c.Analysis* ptr = self._ptr.analysis(aname)
    #     cdef Analysis pyobj = Analysis.__new__(Analysis)
    #     if not ptr:
    #         return None
    #     pyobj._ptr = ptr
    #     return pyobj

    def readData(self, name_or_stream, fmt="yoda", preload=True):
        cdef c.istringstream iss
        if type(name_or_stream) is str:
            self._ptr.readData_FILE(name_or_stream.encode('utf-8'), preload)
        else:
            _make_iss(iss, name_or_stream)
            self._ptr.readData_ISTR(iss, fmt.encode('utf-8'), preload)

    def writeData(self, file_or_filename, fmt="yoda"):
        cdef c.ostringstream oss
        if type(file_or_filename) is str:
            self._ptr.writeData_FILE(file_or_filename.encode('utf-8'))
        else:
            self._ptr.writeData_OSTR(oss, fmt.encode('utf-8'))
            file_or_filename.write(oss.str().decode('utf-8'))


    def nominalCrossSection(self):
        return self._ptr.nominalCrossSection()

    def finalize(self):
        self._ptr.finalize()

    def dump(self, name, period):
        self._ptr.dump(name.encode('utf-8'), period)

    def mergeYodas(self, filelist, delopts, addopts, matches, unmatches, equiv):
        filelist  = [ f.encode('utf-8') for f in filelist ]
        delopts   = [ d.encode('utf-8') for d in delopts  ]
        addopts   = [ d.encode('utf-8') for d in addopts ]
        matches   = [ d.encode('utf-8') for d in matches ]
        unmatches = [ d.encode('utf-8') for d in unmatches ]
        self._ptr.mergeYodas(filelist, delopts, addopts, matches, unmatches, equiv)

    def merge(self, AnalysisHandler other):
        self._ptr.merge(other._ptr[0])


cdef class Run:
    cdef c.Run *_ptr

    def __cinit__(self, AnalysisHandler h):
        self._ptr = new c.Run(h._ptr[0])

    def __del__(self):
        del self._ptr

    def setCrossSection(self, double x):
        self._ptr.setCrossSection(x)
        return self

    def setListAnalyses(self, choice):
        self._ptr.setListAnalyses(choice)
        return self

    def init(self, name, weight=1.0):
        return self._ptr.init(name.encode('utf-8'), weight)

    def openFile(self, name, weight=1.0):
        return self._ptr.openFile(name.encode('utf-8'), weight)

    def readEvent(self):
        return self._ptr.readEvent()

    # def skipEvent(self):
    #     return self._ptr.skipEvent()

    def numEvents(self):
        return self._ptr.numEvents()

    def processEvent(self):
        return self._ptr.processEvent()

    def finalize(self):
        return self._ptr.finalize()


cdef class Analysis:
    cdef c.unique_ptr[c.Analysis] _ptr

    def __init__(self):
        raise RuntimeError('This class cannot be instantiated')

    def requiredBeams(self):
        return deref(self._ptr).requiredBeams()

    def requiredEnergies(self):
        return deref(self._ptr).requiredEnergies()

    def keywords(self):
        kws = deref(self._ptr).keywords()
        return [ k.decode('utf-8') for k in kws ]

    def validation(self):
        vld = deref(self._ptr).validation()
        return [ k.decode('utf-8') for k in vld ]

    def reentrant(self):
        return deref(self._ptr).reentrant()

    def authors(self):
        auths = deref(self._ptr).authors()
        return [ a.decode('utf-8') for a in auths ]

    def bibKey(self):
        return deref(self._ptr).bibKey().decode('utf-8')

    def name(self):
        return deref(self._ptr).name().decode('utf-8')

    def bibTeX(self):
        return deref(self._ptr).bibTeX().decode('utf-8')

    def references(self):
        refs = deref(self._ptr).references()
        return [ r.decode('utf-8') for r  in refs ]

    def collider(self):
        return deref(self._ptr).collider().decode('utf-8')

    def description(self):
        return deref(self._ptr).description().decode('utf-8')

    def experiment(self):
        return deref(self._ptr).experiment().decode('utf-8')

    def inspireId(self):
        return deref(self._ptr).inspireId().decode('utf-8')

    def spiresId(self):
        return deref(self._ptr).spiresId().decode('utf-8')

    def runInfo(self):
        return deref(self._ptr).runInfo().decode('utf-8')

    def status(self):
        return deref(self._ptr).status().decode('utf-8')

    def warning(self):
        return deref(self._ptr).warning().decode('utf-8')

    def summary(self):
        return deref(self._ptr).summary().decode('utf-8')

    def year(self):
        return deref(self._ptr).year().decode('utf-8')

    def luminosity(self):
        return deref(self._ptr).luminosity()

    def luminosityfb(self):
        return deref(self._ptr).luminosityfb()

    def refMatch(self):
        return deref(self._ptr).refMatch().decode('utf-8')

    def refUnmatch(self):
        return deref(self._ptr).refUnmatch().decode('utf-8')

    def writerDoublePrecision(self):
        return deref(self._ptr).writerDoublePrecision().decode('utf-8')

    def refFile(self):
        #return findAnalysisRefFile(self.name() + ".yoda")
        return deref(self._ptr).refFile().decode('utf-8')

    def refData(self, asdict=True, patterns=None, unpatterns=None):
        """Get this analysis' reference data, cf. yoda.read()
        NB. There's also a C++ version of this, but this wrapping is nicer for Python"""
        import yoda
        return yoda.read(self.refFile(), asdict, patterns, unpatterns)


#cdef object
LEVELS = dict(TRACE = 0, DEBUG = 10, INFO = 20,
              WARN = 30, WARNING = 30, ERROR = 40,
              CRITICAL = 50, ALWAYS = 50)


cdef class AnalysisLoader:

    @staticmethod
    def analysisNames():
        names = c.AnalysisLoader_analysisNames()
        return [ n.decode('utf-8') for n in names ]

    @staticmethod
    def allAnalysisNames():
        names = c.AnalysisLoader_allAnalysisNames()
        return [ n.decode('utf-8') for n in names ]

    @staticmethod
    def stdAnalysisNames():
        names = c.AnalysisLoader_stdAnalysisNames()
        return [ n.decode('utf-8') for n in names ]

    @staticmethod
    def analysisNameAliases():
        anames = c.AnalysisLoader_analysisNameAliases()
        return { a.first.decode('utf-8') : a.second.decode('utf-8') for a in anames }

    @staticmethod
    def getAnalysis(name):
        try:
          name = name.encode('utf-8')
        except AttributeError:
          pass
        cdef c.unique_ptr[c.Analysis] ptr = c.AnalysisLoader_getAnalysis(name)
        cdef Analysis pyobj = Analysis.__new__(Analysis)
        if not ptr:
            return None
        pyobj._ptr = move(ptr)
        # Create python object
        return pyobj


## Convenience versions in main rivet namespace
def analysisNames():
    return AnalysisLoader.analysisNames()

def allAnalysisNames():
    return AnalysisLoader.allAnalysisNames()

def stdAnalysisNames():
    return AnalysisLoader.stdAnalysisNames()

def analysisNameAliases():
    return AnalysisLoader.analysisNameAliases()

def getAnalysis(name):
    return AnalysisLoader.getAnalysis(name.encode('utf-8'))


## Path functions
def getAnalysisLibPaths():
    ps = c.getAnalysisLibPaths()
    return [ p.decode('utf-8') for p in ps ]

def setAnalysisLibPaths(xs):
    bs = [ x.encode('utf-8') for x in xs ]
    c.setAnalysisLibPaths(bs)

def addAnalysisLibPath(path):
    c.addAnalysisLibPath(path.encode('utf-8'))


def setAnalysisDataPaths(xs):
    bs = [ x.encode('utf-8') for x in xs ]
    c.setAnalysisDataPaths(bs)

def addAnalysisDataPath(path):
    c.addAnalysisDataPath(path.encode('utf-8'))

def getAnalysisDataPaths():
    ps = c.getAnalysisDataPaths()
    return [ p.decode('utf-8') for p in ps ]

def findAnalysisDataFile(q):
    f = c.findAnalysisDataFile(q.encode('utf-8'))
    return f.decode('utf-8')

def getAnalysisRefPaths():
    ps = c.getAnalysisRefPaths()
    return [ p.decode('utf-8') for p in ps ]

def findAnalysisRefFile(q):
    f = c.findAnalysisRefFile(q.encode('utf-8'))
    return f.decode('utf-8')


def getAnalysisInfoPaths():
    ps = c.getAnalysisInfoPaths()
    return [ p.decode('utf-8') for p in ps ]

def findAnalysisInfoFile(q):
    f = c.findAnalysisInfoFile(q.encode('utf-8'))
    return f.decode('utf-8')

def getAnalysisPlotPaths():
    ps = c.getAnalysisPlotPaths()
    return [ p.decode('utf-8') for p in ps ]

def findAnalysisPlotFile(q):
    f = c.findAnalysisPlotFile(q.encode('utf-8'))
    return f.decode('utf-8')

def version():
    return c.version().decode('utf-8')

def setLogLevel(name, level):
    c.setLogLevel(name.encode('utf-8'), level)
