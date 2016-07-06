cimport rivet as c
# Need to be careful with memory management -- perhaps use the base object that
# we used in YODA?

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

    def processEvent(self):
        return self._ptr.processEvent()

    def finalize(self):
        return self._ptr.finalize()


cdef class Analysis:
    cdef c.Analysis *_ptr

    def __init__(self):
        raise RuntimeError('This class cannot be instantiated')
    def __del__(self):
        del self._ptr

    def requiredBeams(self):
        return self._ptr.requiredBeams()

    def requiredEnergies(self):
        return self._ptr.requiredEnergies()

    def authors(self):
        return self._ptr.authors()

    def bibKey(self):
        return self._ptr.bibKey()

    def name(self):
        return self._ptr.name()

    def bibTeX(self):
        return self._ptr.bibTeX()

    def references(self):
        return self._ptr.references()

    def collider(self):
        return self._ptr.collider()

    def description(self):
        return self._ptr.description()

    def experiment(self):
        return self._ptr.experiment()

    def inspireId(self):
        return self._ptr.inspireId()

    def spiresId(self):
        return self._ptr.spiresId()

    def runInfo(self):
        return self._ptr.runInfo()

    def status(self):
        return self._ptr.status()

    def summary(self):
        return self._ptr.summary()

    def year(self):
        return self._ptr.year()


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
        cdef c.Analysis* ptr = c.AnalysisLoader_getAnalysis(name)
        cdef Analysis pyobj = Analysis.__new__(Analysis)
        if not ptr:
            return None
        pyobj._ptr = ptr
        # Create python object
        return pyobj


def addAnalysisLibPath(path):
    c.addAnalysisLibPath(path)

def findAnalysisRefFile(q):
    return c.findAnalysisRefFile(q)

def getAnalysisPlotPaths():
    return c.getAnalysisPlotPaths()

def getAnalysisRefPaths():
    return c.getAnalysisRefPaths()

def getAnalysisLibPaths():
    return c.getAnalysisLibPaths()

def setAnalysisLibPaths(xs):
    c.setAnalysisLibPaths(xs)

def version():
    return c.version()

def setLogLevel(name, level):
    c.setLogLevel(name, level)
