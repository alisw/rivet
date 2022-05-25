from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.memory cimport unique_ptr

ctypedef int PdgId
ctypedef pair[PdgId,PdgId] PdgIdPair

cdef extern from "<sstream>" namespace "std":
    cdef cppclass ostringstream:
        ostringstream()
        string& str()
    cdef cppclass istringstream:
        istringstream()
        string& str(string&)

cdef extern from "Rivet/AnalysisHandler.hh" namespace "Rivet":
    cdef cppclass AnalysisHandler:
        void setIgnoreBeams(bool)
        void skipMultiWeights(bool)
        void selectMultiWeights(string)
        void deselectMultiWeights(string)
        void setNominalWeightName(string)
        void setWeightCap(double)
        void setNLOSmearing(double)
        AnalysisHandler& addAnalysis(string)
        vector[string] analysisNames() const
        vector[string] stdAnalysisNames() const
        # Analysis* analysis(string)
        void writeData_FILE "writeData" (string&) except +
        void writeData_OSTR "writeData" (ostringstream&, string&) except +
        void readData_FILE "readData" (string&, bool) except +
        void readData_ISTR "readData" (istringstream&, string&, bool) except +
        double nominalCrossSection()
        void finalize()
        void dump(string, int)
        void mergeYodas(vector[string]&, vector[string]&, vector[string]&, vector[string]&, vector[string]&, bool)
        void merge(AnalysisHandler&)

cdef extern from "Rivet/Run.hh" namespace "Rivet":
    cdef cppclass Run:
        Run(AnalysisHandler)
        Run& setCrossSection(double) # For chaining?
        Run& setListAnalyses(bool)
        bool init(string, double) except + # $2=1.0
        bool openFile(string, double) except + # $2=1.0
        bool readEvent() except +
        #bool skipEvent() except +
        bool processEvent() except +
        bool finalize() except +
        size_t numEvents()

cdef extern from "Rivet/Analysis.hh" namespace "Rivet":
    cdef cppclass Analysis:
        vector[PdgIdPair]& requiredBeams()
        vector[pair[double, double]] requiredEnergies()
        vector[string] authors()
        vector[string] references()
        vector[string] keywords()
        vector[string] validation()
        bool reentrant()
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
        string warning()
        string summary()
        string year()
        double luminosity()
        double luminosityfb()
        string refFile()
        string refMatch()
        string refUnmatch()
        string writerDoublePrecision()


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
    vector[string] AnalysisLoader_allAnalysisNames "Rivet::AnalysisLoader::allAnalysisNames" ()
    map[string,string] AnalysisLoader_analysisNameAliases "Rivet::AnalysisLoader::analysisNameAliases" ()
    vector[string] AnalysisLoader_stdAnalysisNames "Rivet::AnalysisLoader::stdAnalysisNames" ()
    unique_ptr[Analysis] AnalysisLoader_getAnalysis "Rivet::AnalysisLoader::getAnalysis" (string)

cdef extern from "Rivet/Tools/RivetPaths.hh" namespace "Rivet":
    vector[string] getAnalysisLibPaths()
    void setAnalysisLibPaths(vector[string])
    void addAnalysisLibPath(string)

    vector[string] getAnalysisDataPaths()
    void setAnalysisDataPaths(vector[string])
    void addAnalysisDataPath(string)
    string findAnalysisDataFile(string)

    vector[string] getAnalysisRefPaths()
    string findAnalysisRefFile(string)

    vector[string] getAnalysisInfoPaths()
    string findAnalysisInfoFile(string)

    vector[string] getAnalysisPlotPaths()
    string findAnalysisPlotFile(string)

cdef extern from "Rivet/Rivet.hh" namespace "Rivet":
    string version()

cdef extern from "Rivet/Tools/Logging.hh":
    void setLogLevel "Rivet::Log::setLevel" (string, int)
