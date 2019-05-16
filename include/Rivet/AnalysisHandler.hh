// -*- C++ -*-
#ifndef RIVET_RivetHandler_HH
#define RIVET_RivetHandler_HH

#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Particle.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Tools/RivetYODA.hh"

namespace Rivet {


  // Forward declaration and smart pointer for Analysis
  class Analysis;
  typedef std::shared_ptr<Analysis> AnaHandle;


  // // Needed to make smart pointers compare equivalent in the STL set
  // struct CmpAnaHandle {
  //   bool operator() (const AnaHandle& a, const AnaHandle& b) const {
  //     return a.get() < b.get();
  //   }
  // };


  /// A class which handles a number of analysis objects to be applied to
  /// generated events. An {@link Analysis}' AnalysisHandler is also responsible
  /// for handling the final writing-out of histograms.
  class AnalysisHandler {
  public:

    /// @name Constructors and destructors. */
    //@{

    /// Preferred constructor, with optional run name.
    AnalysisHandler(const string& runname="");

    /// @brief Destructor
    /// The destructor is not virtual, as this class should not be inherited from.
    ~AnalysisHandler();

    //@}


  private:

    /// Get a logger object.
    Log& getLog() const;


  public:

    /// @name Run properties
    //@{

    /// Get the name of this run.
    string runName() const { return _runname; }


    /// Get the number of events seen. Should only really be used by external
    /// steering code or analyses in the finalize phase.
    size_t numEvents() const { return _eventcounter.numEntries(); }

    /// @brief Access the sum of the event weights seen
    ///
    /// This is the weighted equivalent of the number of events. It should only
    /// be used by external steering code or analyses in the finalize phase.
    double sumW() const { return _eventcounter.sumW(); }
    /// Access to the sum of squared-weights
    double sumW2() const { return _eventcounter.sumW2(); }

    /// @brief Compatibility alias for sumOfWeights
    ///
    /// @deprecated Prefer sumW
    double sumOfWeights() const { return sumW(); }

    /// @brief Set the sum of weights
    ///
    /// This is useful if Rivet is steered externally and
    /// the analyses are run for a sub-contribution of the events
    /// (but of course have to be normalised to the total sum of weights)
    ///
    /// @todo What about the sumW2 term? That needs to be set coherently. Need a
    /// new version, with all three N,sumW,sumW2 numbers (or a counter)
    /// supplied.
    ///
    /// @deprecated Weight sums are no longer tracked this way...
    void setSumOfWeights(const double& sum) {
      //_sumOfWeights = sum;
      throw Error("Can't set sum of weights independently, since it's now tracked by a Counter. "
                  "Please contact the Rivet authors if you need this.");
    }


    /// Is cross-section information required by at least one child analysis?
    /// @deprecated Should no-longer be an issue: does any generator not write the cross-section?
    bool needCrossSection() const;
    /// Whether the handler knows about a cross-section.
    /// @deprecated Should no-longer be an issue: does any generator not write the cross-section?
    bool hasCrossSection() const;

    /// Get the cross-section known to the handler.
    double crossSection() const { return _xs; }

    /// Set the cross-section for the process being generated.
    AnalysisHandler& setCrossSection(double xs, double xserr=0);


    /// Set the beam particles for this run
    AnalysisHandler& setRunBeams(const ParticlePair& beams) {
      _beams = beams;
      MSG_DEBUG("Setting run beams = " << beams << " @ " << sqrtS()/GeV << " GeV");
      return *this;
    }

    /// Get the beam particles for this run, usually determined from the first event.
    const ParticlePair& beams() const { return _beams; }

    /// Get beam IDs for this run, usually determined from the first event.
    /// @deprecated Use standalone beamIds(ah.beams()), to clean AH interface
    PdgIdPair beamIds() const;

    /// Get energy for this run, usually determined from the first event.
    /// @deprecated Use standalone sqrtS(ah.beams()), to clean AH interface
    double sqrtS() const;

    /// Setter for _ignoreBeams
    void setIgnoreBeams(bool ignore=true);

    //@}


    /// @name Handle analyses
    //@{

    /// Get a list of the currently registered analyses' names.
    std::vector<std::string> analysisNames() const;

    /// Get the collection of currently registered analyses.
    const std::map<std::string, AnaHandle>& analysesMap() const {
      return _analyses;
    }

    /// Get the collection of currently registered analyses.
    std::vector<AnaHandle> analyses() const {
      std::vector<AnaHandle> rtn;
      rtn.reserve(_analyses.size());
      for (const auto& apair : _analyses) rtn.push_back(apair.second);
      return rtn;
    }

    /// Get a registered analysis by name.
    AnaHandle analysis(const std::string& analysisname) {
      try {
        return _analyses[analysisname];
      } catch (...) {
        throw LookupError("No analysis named '" + analysisname + "' registered in AnalysisHandler");
      }
    }

    /// Add an analysis to the run list by object
    AnalysisHandler& addAnalysis(Analysis* analysis);

    /// @brief Add an analysis to the run list using its name.
    ///
    /// The actual Analysis to be used will be obtained via
    /// AnalysisLoader::getAnalysis(string).  If no matching analysis is found,
    /// no analysis is added (i.e. the null pointer is checked and discarded.
    AnalysisHandler& addAnalysis(const std::string& analysisname);

    /// @brief Add an analysis with a map of analysis options.
    AnalysisHandler& addAnalysis(const std::string& analysisname, std::map<string, string> pars);

    /// @brief Add analyses to the run list using their names.
    ///
    /// The actual {@link Analysis}' to be used will be obtained via
    /// AnalysisHandler::addAnalysis(string), which in turn uses
    /// AnalysisLoader::getAnalysis(string). If no matching analysis is found
    /// for a given name, no analysis is added, but also no error is thrown.
    AnalysisHandler& addAnalyses(const std::vector<std::string>& analysisnames);


    /// Remove an analysis from the run list using its name.
    AnalysisHandler& removeAnalysis(const std::string& analysisname);

    /// Remove analyses from the run list using their names.
    AnalysisHandler& removeAnalyses(const std::vector<std::string>& analysisnames);

    //@}


    /// @name Main init/execute/finalise
    //@{

    /// Initialize a run, with the run beams taken from the example event.
    void init(const GenEvent& event);

    /// @brief Analyze the given \a event by reference.
    ///
    /// This function will call the AnalysisBase::analyze() function of all
    /// included analysis objects.
    void analyze(const GenEvent& event);

    /// @brief Analyze the given \a event by pointer.
    ///
    /// This function will call the AnalysisBase::analyze() function of all
    /// included analysis objects, after checking the event pointer validity.
    void analyze(const GenEvent* event);

    /// Finalize a run. This function calls the AnalysisBase::finalize()
    /// functions of all included analysis objects.
    void finalize();

    //@}


    /// @name Histogram / data object access
    //@{

    /// Add a vector of analysis objects to the current state.
    void addData(const std::vector<AnalysisObjectPtr>& aos);

    /// Read analysis plots into the histo collection (via addData) from the named file.
    void readData(const std::string& filename);

    /// Get all analyses' plots as a vector of analysis objects.
    std::vector<AnalysisObjectPtr> getData(bool includeorphans = false,
                                           bool includetmps = false,
                                           bool usefinalized = true) const;

    /// Write all analyses' plots (via getData) to the named file.
    void writeData(const std::string& filename) const;

    /// Tell Rivet to dump intermediate result to a file named @a
    /// dumpfile every @a period'th event. If @period is not positive,
    /// no dumping will be done.
    void dump(string dumpfile, int period) {
      _dumpPeriod = period;
      _dumpFile = dumpfile;
    }

    /// Take the vector of yoda files and merge them together using
    /// the cross section and weight information provided in each
    /// file. Each file in @a aofiles is assumed to have been produced
    /// by Rivet. By default the files are assumed to contain
    /// different processes (or the same processs but mutually
    /// exclusive cuts), but if @a equiv if ture, the files are
    /// assumed to contain output of completely equivalent (but
    /// statistically independent) Rivet runs. The corresponding
    /// analyses will be loaded and their analysis objects will be
    /// filled with the merged result. finalize() will be run on each
    /// relevant anslysis. The resulting YODA file can then be rwitten
    /// out by writeData(). If delopts is non-empty, it is assumed to
    /// contain names different options to be merged into the same
    /// analysis objects.
    void mergeYodas(const vector<string> & aofiles,
                    const vector<string> & delopts = vector<string>(),
                    bool equiv = false);

    /// Helper function to strip specific options from data object paths.
    void stripOptions(AnalysisObjectPtr ao,
                      const vector<string> & delopts) const;

    //@}

  private:

    /// The collection of Analysis objects to be used.
    std::map<std::string, AnaHandle> _analyses;

    /// A vector of pre-loaded object which do not have a valid
    /// Analysis plugged in.
    vector<AnalysisObjectPtr> _orphanedPreloads;

    /// A vector containing copies of analysis objects after
    /// finalize() has been run.
    vector<AnalysisObjectPtr> _finalizedAOs;

    /// @name Run properties
    //@{

    /// Run name
    std::string _runname;

    /// Event counter
    Counter _eventcounter;

    /// Cross-section known to AH
    double _xs, _xserr;

    /// Beams used by this run.
    ParticlePair _beams;

    /// Flag to check if init has been called
    bool _initialised;

    /// Flag whether input event beams should be ignored in compatibility check
    bool _ignoreBeams;

    /// Determines how often Rivet runs finalize() and writes the
    /// result to a YODA file.
    int _dumpPeriod;

    /// The name of a YODA file to which Rivet periodically dumps
    /// results.
    string _dumpFile;

    /// Flag to indicate periodic dumping is in progress
    int _dumping;

    //@}


  private:

    /// The assignment operator is private and must never be called.
    /// In fact, it should not even be implemented.
    AnalysisHandler& operator=(const AnalysisHandler&);

    /// The copy constructor is private and must never be called.  In
    /// fact, it should not even be implemented.
    AnalysisHandler(const AnalysisHandler&);

  };


}

#endif
