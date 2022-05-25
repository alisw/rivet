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


  /// @brief The key class for coordination of Analysis objects and the event loop
  ///
  /// A class which handles a number of analysis objects to be applied to
  /// generated events. An {@link Analysis}' AnalysisHandler is also responsible
  /// for handling the final writing-out of histograms.
  class AnalysisHandler {
  public:

    /// Preferred constructor, with optional run name.
    AnalysisHandler(const string& runname="");

    /// The copy constructor is deleted, so it can never be called.
    AnalysisHandler(const AnalysisHandler&) = delete;

    /// The assignment operator is deleted, so it can never be called.
    AnalysisHandler& operator=(const AnalysisHandler&) = delete;

    /// The destructor is not virtual, as this class should not be inherited from.
    ~AnalysisHandler();


    /// @name Run properties and weights
    /// @{

    /// Get the name of this run.
    string runName() const;

    /// Get the number of events seen. Should only really be used by external
    /// steering code or analyses in the finalize phase.
    ///
    /// N.B. This only reports the count for the last collapsed event group 
    /// and hence ignores any additional sub-events seen so far.
    size_t numEvents() const { 
      const double N = _eventCounter.get()->_getPersistent(defaultWeightIndex())->numEntries();
      return  size_t(N + 0.5 - (N<0)); // round to nearest integer
    }

    /// @brief Access the sum of the event weights seen
    ///
    /// This is the weighted equivalent of the number of events. It should only
    /// be used by external steering code or analyses in the finalize phase.
    double sumW() const { return _eventCounter->sumW(); }
    /// Access to the sum of squared-weights
    double sumW2() const { return _eventCounter->sumW2(); }

    /// Names of event weight categories
    const vector<string>& weightNames() const { return _weightNames; }

    /// Indices of the weights in the original weight matrix
    //const vector<size_t> weightIndices() const { return _weightIndices; }

    /// Are any of the weights non-numeric?
    size_t numWeights() const { return _weightNames.size(); }

    /// Are any of the weights non-numeric?
    bool haveNamedWeights() const;

    /// Set the weight names from a GenEvent
    void setWeightNames(const GenEvent& ge);

    /// Get the index of the nominal weight-stream
    size_t defaultWeightIndex() const { return _rivetDefaultWeightIdx; }

    /// Set the weight cap
    void setWeightCap(const double maxWeight) { _weightCap = maxWeight; }

    /// Set the relative width of the NLO smearing window.
    void setNLOSmearing(double frac) { _NLOSmearing = frac; }

    /// Setter for _skipWeights
    void skipMultiWeights(bool ignore=false);

    /// Setter for _matchWeightNames
    void selectMultiWeights(std::string patterns="");

    /// Setter for _unmatchWeightNames
    void deselectMultiWeights(std::string patterns="");

    /// Setter for _nominalWeightName
    void setNominalWeightName(std::string name="");

    /// @}


    /// @name Cross-sections
    /// @{

    /// Get the cross-section known to the handler
    Scatter1DPtr crossSection() const { return _xs; }

    /// Set all cross-sections for the process being generated specifically (preferred)
    void setCrossSection(const vector<pair<double,double>>& xsecs, bool isUserSupplied = false);

    /// Set all cross-sections for the process being generated, based on nominal weight
    void setCrossSection(const pair<double, double>& xsec, bool isUserSupplied=false);

    /// Set the cross-section for the process being generated (alternative signature)
    void setCrossSection(double xsec, double xsecerr, bool isUserSupplied=false) {
      setCrossSection({xsec, xsecerr}, isUserSupplied);
    }

    /// Get the nominal cross-section
    double nominalCrossSection() const {
      _xs.get()->setActiveWeightIdx(_rivetDefaultWeightIdx);
      const YODA::Scatter1D::Points& ps = _xs->points();
      if (ps.size() != 1) {
        string errMsg = "value missing when requesting nominal cross-section";
        throw Error(errMsg);
      }
      double xs = ps[0].x();
      _xs.get()->unsetActiveWeight();
      return xs;
    }

    /// @}


    /// @name Beams
    /// @{

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

    /// Option to disable AH beam-consistency checks
    void checkBeams(bool check=true) { setIgnoreBeams(!check); }
    /// Opposite of checkBeams()
    /// @deprecated Use checkBeams()
    void setIgnoreBeams(bool ignore=true);

    /// @}


    /// @name Analysis handling
    /// @{

    /// Get a list of the currently registered analyses' names.
    std::vector<std::string> analysisNames() const;

    /// Get a list of the official analysis names for this release.
    std::vector<std::string> stdAnalysisNames() const;

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
      if ( _analyses.find(analysisname) == _analyses.end() )
        throw LookupError("No analysis named '" + analysisname + "' registered in AnalysisHandler");
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

    /// @}


    /// @name Main init/execute/finalise
    /// @{

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

    /// @}


    /// @name Histogram / data object access
    /// @{

    /// @brief Read analysis plots into the histo collection from the given stream
    ///
    /// Use the @a fmt flag to specify the YODA output format (yoda, yoda.gz, yoda.h5, ...)
    void readData(std::istream& istr, const string& fmt, bool preload = true);

    /// Read analysis plots into the histo collection (via addData) from the named file.
    void readData(const std::string& filename, bool preload = true);

    /// Get all YODA analysis objects (across all weights, optionally including RAW)
    vector<YODA::AnalysisObjectPtr> getYodaAOs(bool includeraw=false) const;

    /// Get a pointer to a preloaded yoda object with the given path,
    /// or null if path is not found.
    const YODA::AnalysisObjectPtr getPreload(string path) const {
      auto it = _preloads.find(path);
      if ( it == _preloads.end() ) return nullptr;
      return it->second;
    }

    /// @brief Write all analyses' plots (via getData) to the given stream.
    ///
    /// Use the @a fmt flag to specify the YODA output format (yoda, yoda.gz, yoda.h5, ...)
    void writeData(std::ostream& ostr, const string& fmt) const;

    /// Write all analyses' plots (via getData) to the named file.
    void writeData(const string& filename) const;

    /// @brief Configure the AnalysisObject dump rate and destination.
    ///
    /// Tell Rivet to dump intermediate result to a file named @a
    /// dumpfile every @a period'th event. If @a period is not positive,
    /// no dumping will be done.
    void setAODump(const string& dumpfile, int period) {
      dump(dumpfile, period);
    }
    /// @brief Configure the AnalysisObject dump rate and destination.
    void setNoAODump() {
      setAODump("DUMMY", -1);
    }
    /// Alias for setAODump()
    /// @deprecated Prefer setAODump()
    void dump(const string& dumpfile, int period) {
      _dumpPeriod = period;
      _dumpFile = dumpfile;
    }

    /// @brief Merge the vector of YODA files, using the cross-section and weight information provided in each.
    ///
    /// Each file in @a aofiles is assumed to have been produced by Rivet. By
    /// default the files are assumed to contain different processes (or the
    /// same processs but mutually exclusive cuts), but if @a equiv if true, the
    /// files are assumed to contain output of completely equivalent (but
    /// statistically independent) Rivet runs. The corresponding analyses will
    /// be loaded and their analysis objects will be filled with the merged
    /// result. finalize() will be run on each relevant analysis. The resulting
    /// YODA file can then be written out by writeData().
    ///
    /// If @a delopts is non-empty, it is assumed to contain names of different
    /// options to be merged into the same analysis objects.
    ///

    void mergeYodas(const vector<string>& aofiles,
                    const vector<string>& delopts=vector<string>(),
                    const vector<string>& addopts=vector<string>(),
                    const vector<string>& matches=vector<string>(),
                    const vector<string>& unmatches=vector<string>(),
                    bool equiv=false);

    /// A method to merge another AnalysisHandler into the current one
    void merge(AnalysisHandler &other);

    /// @}


    /// @name Processing stage
    /// @{

    /// Indicate which Rivet stage we're in.
    /// At the moment, only INIT is used to enable booking.
    enum class Stage { OTHER, INIT, FINALIZE };

    /// Return the current processing stage.
    Stage stage() const { return _stage; }

    /// @}


  private:

    /// @name Internal helper functions
    /// @{

    /// Get a logger object.
    Log& getLog() const;

    /// Get all multi-weight Rivet analysis object wrappers.
    vector<MultiweightAOPtr> getRivetAOs() const;

    /// Helper function to strip specific options from data object paths.
    void stripOptions(YODA::AnalysisObjectPtr ao, const vector<string>& delopts) const;

    /// After all subevents in an event group have been processed, push
    /// all histo fills to the relevant histograms.
    void pushToPersistent();

    /// @brief Merge the AO map @a newaos into @a allaos
    void mergeAOS(map<string, YODA::AnalysisObjectPtr> &allaos,
                  map<string, YODA::AnalysisObject*> &newaos, 
                  map<string, pair<double, double>> &allxsecs,
                  const vector<string>& delopts=vector<string>(),
                  const vector<string>& optAnas=vector<string>(),
                  const vector<string>& optKeys=vector<string>(),
                  const vector<string>& optVals=vector<string>(),
                  bool equiv=false,
                  const bool overwrite_xsec = false,
                  const double user_xsec = 1.0);

     
    /// @brief A method to prepare a re-entrant run for a given set of analysis objects
    ///
    /// The @a unscale parameter multiplies fillable objects with sumW/xsec to counteract
    /// the cross-section scaling in finalize() when merging different processes (non-equiv)
    void loadAOs(const map<string, YODA::AnalysisObjectPtr>& allAOs, const bool unscale = false);

    /// @}


  private:

    /// Current handler processing stage.
    Stage _stage = Stage::OTHER;

    /// The collection of Analysis objects to be used.
    std::map<std::string, AnaHandle> _analyses;

    /// A vector of pre-loaded object which do not have a valid Analysis plugged in.
    ///
    /// @todo Rename to _preloadedAOs for consistency
    map<string,YODA::AnalysisObjectPtr> _preloads;

    /// A vector containing copies of analysis objects after finalize() has been run.
    vector<YODA::AnalysisObjectPtr> _finalizedAOs;


    /// @name Run properties
    /// @{

    /// Weight names
    std::vector<std::string> _weightNames;
    std::vector<std::valarray<double> > _subEventWeights;
    //size_t _numWeightTypes; // always == WeightVector.size()

    /// Weight indices
    std::vector<size_t> _weightIndices;

    /// Run name
    std::string _runname;

    /// Event counter
    CounterPtr _eventCounter;

    /// Cross-section known to AH
    Scatter1DPtr _xs;

    /// Nominal cross-section
    std::pair<double,double> _userxs;

    /// Beams used by this run.
    ParticlePair _beams;

    /// Flag to check if init has been called
    bool _initialised;

    /// Flag whether input event beams should be ignored in compatibility check
    bool _ignoreBeams;

    /// Flag to check if multiweights should be included
    bool _skipWeights;

    /// String of weight names (or regex) to select multiweights
    std::string _matchWeightNames;

    /// String of weight names (or regex) to veto multiweights
    std::string _unmatchWeightNames;

    /// String giving the nominal weight name
    std::string _nominalWeightName;

    /// weight cap value
    double _weightCap;

    /// The relative width of the NLO smearing window.
    ///
    /// @todo Improve & standardise name
    double _NLOSmearing;

    /// Current event number
    int _eventNumber;

    /// The index in the (original) weight vector for the nominal weight stream
    size_t _defaultWeightIdx;

    /// The index in the (possibly pruned) weight vector for the nominal weight stream
    size_t _rivetDefaultWeightIdx;

    /// How often Rivet runs finalize() and writes the result to a YODA file.
    int _dumpPeriod;

    /// The name of a YODA file to which Rivet periodically dumps results.
    string _dumpFile;

    /// Flag to indicate periodic dumping is in progress
    bool _dumping;

    /// @}

  };


}

#endif
