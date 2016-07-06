from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libcpp cimport bool
from libcpp.string cimport string

ctypedef int PdgId
ctypedef pair[PdgId,PdgId] PdgIdPair

cdef extern from "Rivet/AnalysisHandler.hh" namespace "Rivet":
    cdef cppclass AnalysisHandler:
        void setIgnoreBeams(bool)
        AnalysisHandler& addAnalysis(string)
        void writeData(string&)
        double crossSection()
        void finalize()

cdef extern from "Rivet/Run.hh" namespace "Rivet":
    cdef cppclass Run:
        Run(AnalysisHandler)
        Run& setCrossSection(double) # For chaining?
        Run& setListAnalyses(bool)
        bool init(string, double) # $2=1.0
        bool openFile(string, double) # $2=1.0
        bool readEvent()
        bool processEvent()
        bool finalize()

cdef extern from "Rivet/Analysis.hh" namespace "Rivet":
    cdef cppclass Analysis:
        vector[PdgIdPair]& requiredBeams()
        vector[pair[double, double]] requiredEnergies()
        vector[string] authors()
        vector[string] references()
        string name()
        string bibTeX()
        string bibKey()
        string collider()
        string description()
        string experiment()
        string inspireId()
        string spiresId()
        string runInfo()
        string status()
        string summary()
        string year()

# Might need to translate the following errors, although I believe 'what' is now
# preserved. But often, we need the exception class name.
#Error
#RangeError
#LogicError
#PidError
#InfoError
#WeightError
#UserError

cdef extern from "Rivet/AnalysisLoader.hh":
    vector[string] AnalysisLoader_analysisNames "Rivet::AnalysisLoader::analysisNames" ()
    Analysis* AnalysisLoader_getAnalysis "Rivet::AnalysisLoader::getAnalysis" (string)

cdef extern from "Rivet/Tools/RivetPaths.hh" namespace "Rivet":
    void addAnalysisLibPath(string)
    string findAnalysisRefFile(string)
    vector[string] getAnalysisPlotPaths()
    vector[string] getAnalysisRefPaths()
    vector[string] getAnalysisLibPaths()
    void setAnalysisLibPaths(vector[string])

cdef extern from "Rivet/Rivet.hh" namespace "Rivet":
    string version()

cdef extern from "Rivet/Tools/Logging.hh":
    void setLogLevel "Rivet::Log::setLevel" (string, int)
