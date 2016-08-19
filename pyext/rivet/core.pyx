# distutils: language = c++

cimport rivet as c
from cython.operator cimport dereference as deref
# Need to be careful with memory management -- perhaps use the base object that
# we used in YODA?

cdef extern from "<utility>" namespace "std" nogil:
    cdef c.unique_ptr[c.Analysis] move(c.unique_ptr[c.Analysis])

cdef class AnalysisHandler:
    cdef c.AnalysisHandler *_ptr

    def __cinit__(self):
        self._ptr = new c.AnalysisHandler()

    def __del__(self):
        del self._ptr

    def setIgnoreBeams(self, ignore=True):
        self._ptr.setIgnoreBeams(ignore)

    def addAnalysis(self, name):
        self._ptr.addAnalysis(name)
        return self

    def analysisNames(self):
        anames = self._ptr.analysisNames()
        return [a for a in anames]

    # def analysis(self, aname):
    #     cdef c.Analysis* ptr = self._ptr.analysis(aname)
    #     cdef Analysis pyobj = Analysis.__new__(Analysis)
    #     if not ptr:
    #         return None
    #     pyobj._ptr = ptr
    #     return pyobj

    def writeData(self, name):
        self._ptr.writeData(name)

    def crossSection(self):
        return self._ptr.crossSection()

    def finalize(self):
        self._ptr.finalize()


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
        return self._ptr.init(name, weight)

    def openFile(self, name, weight=1.0):
        return self._ptr.openFile(name, weight)

    def readEvent(self):
        return self._ptr.readEvent()

    def skipEvent(self):
        return self._ptr.skipEvent()

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

    def authors(self):
        return deref(self._ptr).authors()

    def bibKey(self):
        return deref(self._ptr).bibKey()

    def name(self):
        return deref(self._ptr).name()

    def bibTeX(self):
        return deref(self._ptr).bibTeX()

    def references(self):
        return deref(self._ptr).references()

    def collider(self):
        return deref(self._ptr).collider()

    def description(self):
        return deref(self._ptr).description()

    def experiment(self):
        return deref(self._ptr).experiment()

    def inspireId(self):
        return deref(self._ptr).inspireId()

    def spiresId(self):
        return deref(self._ptr).spiresId()

    def runInfo(self):
        return deref(self._ptr).runInfo()

    def status(self):
        return deref(self._ptr).status()

    def summary(self):
        return deref(self._ptr).summary()

    def year(self):
        return deref(self._ptr).year()


#cdef object
LEVELS = dict(TRACE = 0, DEBUG = 10, INFO = 20,
              WARN = 30, WARNING = 30, ERROR = 40,
              CRITICAL = 50, ALWAYS = 50)


cdef class AnalysisLoader:
    @staticmethod
    def analysisNames():
        return c.AnalysisLoader_analysisNames()

    @staticmethod
    def getAnalysis(name):
        cdef c.unique_ptr[c.Analysis] ptr = c.AnalysisLoader_getAnalysis(name)
        cdef Analysis pyobj = Analysis.__new__(Analysis)
        if not ptr:
            return None
        pyobj._ptr = move(ptr)
        # Create python object
        return pyobj


def getAnalysisLibPaths():
    return c.getAnalysisLibPaths()

def setAnalysisLibPaths(xs):
    c.setAnalysisLibPaths(xs)

def addAnalysisLibPath(path):
    c.addAnalysisLibPath(path)


def setAnalysisDataPaths(xs):
    c.setAnalysisDataPaths(xs)

def addAnalysisDataPath(path):
    c.addAnalysisDataPath(path)

def getAnalysisDataPaths():
    return c.getAnalysisDataPaths()

def findAnalysisDataFile(q):
    return c.findAnalysisDataFile(q)


def getAnalysisRefPaths():
    return c.getAnalysisRefPaths()

def findAnalysisRefFile(q):
    return c.findAnalysisRefFile(q)


def getAnalysisInfoPaths():
    return c.getAnalysisInfoPaths()

def findAnalysisInfoFile(q):
    return c.findAnalysisInfoFile(q)


def getAnalysisPlotPaths():
    return c.getAnalysisPlotPaths()

def findAnalysisPlotFile(q):
    return c.findAnalysisPlotFile(q)


def version():
    return c.version()

def setLogLevel(name, level):
    c.setLogLevel(name, level)
