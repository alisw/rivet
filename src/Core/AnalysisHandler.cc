// -*- C++ -*-
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/ParticleName.hh"
#include "Rivet/Tools/BeamConstraint.hh"
#include "Rivet/Tools/RivetPaths.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/Beam.hh"
#include "YODA/IO.h"
#include <iostream>
#include <regex>

using std::cout;
using std::cerr;

namespace Rivet {


  AnalysisHandler::AnalysisHandler(const string& runname)
    : _runname(runname), _userxs{NAN, NAN},
      _initialised(false), _ignoreBeams(false),
      _skipWeights(false), _matchWeightNames(""),
      _unmatchWeightNames(""), _weightCap(0.),
      _NLOSmearing(0.), _defaultWeightIdx(0), 
      _rivetDefaultWeightIdx(0), _dumpPeriod(0), _dumping(false)
  {  }


  AnalysisHandler::~AnalysisHandler() {
      static bool printed = false;
    // Print out MCnet boilerplate
    if (!printed && getLog().getLevel() <= 20) {
      cout << endl
           << "The MCnet usage guidelines apply to Rivet: see http://www.montecarlonet.org/GUIDELINES" << endl
           << "Please acknowledge Rivet in results made using it, and cite https://arxiv.org/abs/1912.05451" << endl;
           // << "https://arxiv.org/abs/1003.0694" << endl;
      printed = true;
    }
  }


  Log& AnalysisHandler::getLog() const {
    return Log::getLog("Rivet.AnalysisHandler");
  }


  /// http://stackoverflow.com/questions/4654636/how-to-determine-if-a-string-is-a-number-with-c
  namespace {
    bool is_number(const std::string& s) {
      std::string::const_iterator it = s.begin();
      while (it != s.end() && std::isdigit(*it)) ++it;
      return !s.empty() && it == s.end();
    }
  }


  /// Check if any of the weightnames is not a number
  bool AnalysisHandler::haveNamedWeights() const {
    bool dec=false;
    for (unsigned int i=0;i<_weightNames.size();++i) {
      string s = _weightNames[i];
      if (!is_number(s)) {
        dec=true;
        break;
      }
    }
    return dec;
  }


  void AnalysisHandler::init(const GenEvent& ge) {
    if (_initialised)
      throw UserError("AnalysisHandler::init has already been called: cannot re-initialize!");

    /// @todo Should the Rivet analysis objects know about weight names?

    setRunBeams(Rivet::beams(ge));
    MSG_DEBUG("Initialising the analysis handler");
    _eventNumber = ge.event_number();

    setWeightNames(ge);
    if (_skipWeights)
        MSG_INFO("Only using nominal weight. Variation weights will be ignored.");
    else if (haveNamedWeights())
        MSG_INFO("Using named weights");
    else
        MSG_INFO("NOT using named weights. Using first weight as nominal weight");

    _eventCounter = CounterPtr(weightNames(), Counter("_EVTCOUNT"));

    // Set the cross section based on what is reported by this event.
    if ( ge.cross_section() ) {
      setCrossSection(HepMCUtils::crossSection(ge));
    }

    // Check that analyses are beam-compatible, and remove those that aren't
    const size_t num_anas_requested = analysisNames().size();
    vector<string> anamestodelete;
    for (const AnaHandle a : analyses()) {
      if (!_ignoreBeams && !a->isCompatible(beams())) {
        //MSG_DEBUG(a->name() << " requires beams " << a->requiredBeams() << " @ " << a->requiredEnergies() << " GeV");
        anamestodelete.push_back(a->name());
      }
    }
    for (const string& aname : anamestodelete) {
      MSG_WARNING("Analysis '" << aname << "' is incompatible with the provided beams: removing");
      removeAnalysis(aname);
    }
    if (num_anas_requested > 0 && analysisNames().empty()) {
      cerr << "All analyses were incompatible with the first event's beams\n"
           << "Exiting, since this probably wasn't intentional!" << endl;
      exit(1);
    }

    // Warn if any analysis' status is not unblemished
    for (const AnaHandle a : analyses()) {
      if ( a->info().preliminary() ) {
        MSG_WARNING("Analysis '" << a->name() << "' is preliminary: be careful, it may change and/or be renamed!");
      } else if ( a->info().obsolete() ) {
        MSG_WARNING("Analysis '" << a->name() << "' is obsolete: please update!");
      } else if (( a->info().unvalidated() ) ) {
        MSG_WARNING("Analysis '" << a->name() << "' is unvalidated: be careful, it may be broken!");
      }
    }

    // Initialize the remaining analyses
    _stage = Stage::INIT;
    for (AnaHandle a : analyses()) {
      MSG_DEBUG("Initialising analysis: " << a->name());
      try {
        // Allow projection registration in the init phase onwards
        a->_allowProjReg = true;
        a->init();
        //MSG_DEBUG("Checking consistency of analysis: " << a->name());
        //a->checkConsistency();
      } catch (const Error& err) {
        cerr << "Error in " << a->name() << "::init method: " << err.what() << endl;
        exit(1);
      }
      MSG_DEBUG("Done initialising analysis: " << a->name());
    }
    _stage = Stage::OTHER;
    _initialised = true;
    MSG_DEBUG("Analysis handler initialised");
  }


  void AnalysisHandler::setWeightNames(const GenEvent& ge) {
    _weightNames = HepMCUtils::weightNames(ge);
    if (_weightNames.empty()) {
      _weightNames.push_back("");
      _rivetDefaultWeightIdx = _defaultWeightIdx = 0;
      _weightIndices = { 0 };
      return;
    } 

    // Find default weights, starting with the preferred empty "" name
    size_t nDefaults = 0;
    _weightIndices.clear();
    for (size_t i = 0, N = _weightNames.size(); i < N; ++i) {
      _weightIndices.push_back(i);
      if (_weightNames[i] == "") {
        if (nDefaults == 0)  _rivetDefaultWeightIdx = _defaultWeightIdx = i;
        ++nDefaults;
      }
    }
    // If no weights with the preferred name, look for acceptable alternatives
    if (nDefaults == 0) {
      for (size_t i = 0, N = _weightNames.size(); i < N; ++i) {
        const string W = toUpper(_weightNames[i]);
        if (W == "WEIGHT" || W == "0" || W == "DEFAULT" || W == "NOMINAL") {
          if (nDefaults == 0) {
            _weightNames[i] = "";
            _rivetDefaultWeightIdx = _defaultWeightIdx = i;
          }
          ++nDefaults;
        }
      }
    }
    // Warn user that no nominal weight could be identified
    if (nDefaults == 0) {
      MSG_WARNING("Could not identify nominal weight. Will continue assuming variations-only run.");
    }
    // Warn if multiple weight names were acceptable alternatives
    if (nDefaults > 1) {
      MSG_WARNING("Found more than " << nDefaults << " default weight candidates. Will use: " << _weightNames[_defaultWeightIdx]);
    }
    if (_skipWeights)  {
      // If running in single-weight mode, remove all bar the nominal weight 
      _weightIndices = { _defaultWeightIdx };
      _weightNames = { _weightNames[_defaultWeightIdx] };
      _rivetDefaultWeightIdx = 0;
    }
    else {
      // check if weight name matches a supplied string/regex
      // and then (de-)select accordingly
      // deseleection takes precedence over selection
      if (_matchWeightNames != "") {
        MSG_DEBUG("Select weight names that match pattern \"" << _matchWeightNames << "\"");
        // compile regex from each string in the comma-separated list
        vector<std::regex> patterns;
        for (const string& pattern : split(_matchWeightNames, ",")) { 
          patterns.push_back( std::regex(pattern) );
        }
        // check which weights match supplied weight-name pattern
        vector<string> selected_subset; _weightIndices.clear();
        for (size_t i = 0, N = _weightNames.size(); i < N; ++i) {
          if (i == _defaultWeightIdx) {
            // default weight cannot be "unselected"
            _rivetDefaultWeightIdx = _weightIndices.size();
            _weightIndices.push_back(i);
            selected_subset.push_back(_weightNames[i]);
            MSG_DEBUG("Selected nominal weight: " << _weightNames[i]);
            continue;
          }
          for (const std::regex& re : patterns) {
            if ( std::regex_match(_weightNames[i], re) ) {
              _weightIndices.push_back(i);
              selected_subset.push_back(_weightNames[i]);
              MSG_DEBUG("Selected variation weight: " << _weightNames[i]);
              break;
            }
          }
        }
        _weightNames = selected_subset;
      }
     if (_unmatchWeightNames != "") {
        MSG_DEBUG("Deselect weight names that match pattern \"" << _unmatchWeightNames << "\"");
        // compile regex from each string in the comma-separated list
        vector<std::regex> patterns;
        for (const string& pattern : split(_unmatchWeightNames, ",")) {
          patterns.push_back( std::regex(pattern) ); 
        }
        // check which weights match supplied weight-name pattern
        vector<string> selected_subset; _weightIndices.clear();
        for (size_t i = 0, N = _weightNames.size(); i < N; ++i) {
          if (i == _defaultWeightIdx) {
            // default weight cannot be vetoed
            _rivetDefaultWeightIdx = _weightIndices.size();
            _weightIndices.push_back(i);
            selected_subset.push_back(_weightNames[i]);
            MSG_DEBUG("Selected nominal weight: " << _weightNames[i]);
            continue;
          }
          bool skip = false;
          for (const std::regex& re : patterns) {
            if ( std::regex_match(_weightNames[i], re) ) { skip = true; break; }
          }
          if (skip) continue;
          _weightIndices.push_back(i);
          selected_subset.push_back(_weightNames[i]);
          MSG_DEBUG("Selected variation weight: " << _weightNames[i]);
        }
        _weightNames = selected_subset;
      }
    }
    // done (de-)selecting weights, add useful debug messages:
    MSG_DEBUG("Default weight name: \"" <<  _weightNames[_rivetDefaultWeightIdx] << "\"");
    MSG_DEBUG("Default weight index (Rivet): " << _rivetDefaultWeightIdx);
    MSG_DEBUG("Default weight index (overall): " << _defaultWeightIdx);
  }


  void AnalysisHandler::analyze(const GenEvent& ge) {
    // Call init with event as template if not already initialised
    if (!_initialised) init(ge);
    assert(_initialised);

    // Ensure that beam details match those from the first event (if we're checking beams)
    if ( !_ignoreBeams ) {
      const PdgIdPair beams = Rivet::beamIds(ge);
      const double sqrts = Rivet::sqrtS(ge);
      if (!compatible(beams, _beams) || !fuzzyEquals(sqrts, sqrtS())) {
        cerr << "Event beams mismatch: "
             << PID::toBeamsString(beams) << " @ " << sqrts/GeV << " GeV" << " vs. first beams "
             << this->beams() << " @ " << this->sqrtS()/GeV << " GeV" << endl;
        exit(1);
      }
    }

    // Create the Rivet event wrapper
    /// @todo Filter/normalize the event here
    bool strip = ( getEnvParam("RIVET_STRIP_HEPMC", string("NOOOO") ) != "NOOOO" );
    Event event(ge, strip);

    // set the cross section based on what is reported by this event.
    // if no cross section
    if ( ge.cross_section() ) setCrossSection(HepMCUtils::crossSection(ge));

    // Won't happen for first event because _eventNumber is set in init()
    if (_eventNumber != ge.event_number()) {

      pushToPersistent();

      _eventNumber = ge.event_number();

    }


    MSG_TRACE("starting new sub event");
    _eventCounter.get()->newSubEvent();

    for (const AnaHandle& a : analyses()) {
        for (auto ao : a->analysisObjects()) {
            ao.get()->newSubEvent();
        }
    }

    _subEventWeights.push_back(pruneWeights(event.weights()));
    if (_weightCap != 0.) {
      MSG_DEBUG("Implementing weight cap using a maximum |weight| = " << _weightCap << " for latest subevent.");
      size_t lastSub = _subEventWeights.size() - 1;
      for (size_t i = 0; i < _subEventWeights[lastSub].size(); ++i) {
        if (abs(_subEventWeights[lastSub][i]) > _weightCap) {
          _subEventWeights[lastSub][i] = sign(_subEventWeights[lastSub][i]) * _weightCap;
        }
      }
    }
    MSG_DEBUG("Analyzing subevent #" << _subEventWeights.size() - 1 << ".");

    _eventCounter->fill();
    // Run the analyses
    for (AnaHandle a : analyses()) {
      MSG_TRACE("About to run analysis " << a->name());
      try {
        a->analyze(event);
      } catch (const Error& err) {
        cerr << "Error in " << a->name() << "::analyze method: " << err.what() << endl;
        exit(1);
      }
      MSG_TRACE("Finished running analysis " << a->name());
    }

    if ( _dumpPeriod > 0 && numEvents() > 0 && numEvents()%_dumpPeriod == 0 ) {
      MSG_DEBUG("Dumping intermediate results to " << _dumpFile << ".");
      _dumping = numEvents()/_dumpPeriod;
      finalize();
      writeData(_dumpFile);
      _dumping = 0;
    }

  }


  void AnalysisHandler::analyze(const GenEvent* ge) {
    if (ge == nullptr) {
      MSG_ERROR("AnalysisHandler received null pointer to GenEvent");
      //throw Error("AnalysisHandler received null pointer to GenEvent");
    }
    analyze(*ge);
  }


  void AnalysisHandler::pushToPersistent() {
    if ( _subEventWeights.empty() ) return;
    MSG_TRACE("AnalysisHandler::analyze(): Pushing _eventCounter to persistent.");
    _eventCounter.get()->pushToPersistent(_subEventWeights);
    for (const AnaHandle& a : analyses()) {
      for (auto ao : a->analysisObjects()) {
        MSG_TRACE("AnalysisHandler::analyze(): Pushing " << a->name()
                  << "'s " << ao->name() << " to persistent.");
        ao.get()->pushToPersistent(_subEventWeights, _NLOSmearing);
      }
      MSG_TRACE("AnalysisHandler::analyze(): finished pushing "
                << a->name() << "'s objects to persistent.");
    }
    _subEventWeights.clear();
  }


  void AnalysisHandler::finalize() {
    if (!_initialised) return;
    MSG_DEBUG("Finalising analyses");

    _stage = Stage::FINALIZE;

    // First push all analyses' objects to persistent and final
    MSG_TRACE("AnalysisHandler::finalize(): Pushing analysis objects to persistent.");
    pushToPersistent();

    // Copy all histos to finalize versions.
    _eventCounter.get()->pushToFinal();
    _xs.get()->pushToFinal();
    for (const AnaHandle& a : analyses())
      for (auto ao : a->analysisObjects())
        ao.get()->pushToFinal();

    for (AnaHandle a : analyses()) {
      if ( _dumping && !a->info().reentrant() )  {
        if ( _dumping == 1 )
          MSG_DEBUG("Skipping finalize in periodic dump of " << a->name() << " as it is not declared re-entrant.");
        continue;
      }
      for (size_t iW = 0; iW < numWeights(); iW++) {
        _eventCounter.get()->setActiveFinalWeightIdx(iW);
        _xs.get()->setActiveFinalWeightIdx(iW);
        for (auto ao : a->analysisObjects())
          ao.get()->setActiveFinalWeightIdx(iW);
        try {
          MSG_TRACE("running " << a->name() << "::finalize() for weight " << iW << ".");
          a->finalize();
        } catch (const Error& err) {
          cerr << "Error in " << a->name() << "::finalize method: " << err.what() << '\n';
          exit(1);
        }
      }
    }

    // Print out number of events processed
    if (!_dumping) {
      const int nevts = numEvents();
      MSG_DEBUG("Processed " << nevts << " event" << (nevts != 1 ? "s" : ""));
    }

    _stage = Stage::OTHER;

  }


  AnalysisHandler& AnalysisHandler::addAnalysis(const string& analysisname, std::map<string, string> pars) {
     // Make an option handle.
    std::string parHandle = "";
    for (map<string, string>::iterator par = pars.begin(); par != pars.end(); ++par) {
      parHandle +=":";
      parHandle += par->first + "=" + par->second;
    }
    return addAnalysis(analysisname + parHandle);
  }


  AnalysisHandler& AnalysisHandler::addAnalysis(const string& analysisname) {
    // Check for a duplicate analysis
    /// @todo Might we want to be able to run an analysis twice, with different params?
    ///       Requires avoiding histo tree clashes, i.e. storing the histos on the analysis objects.
    string ananame = analysisname;
    vector<string> anaopt = split(analysisname, ":");
    if ( anaopt.size() > 1 ) ananame = anaopt[0];
    AnaHandle analysis( AnalysisLoader::getAnalysis(ananame) );
    if (analysis.get() != 0) { // < Check for null analysis.
      MSG_DEBUG("Adding analysis '" << analysisname << "'");
      map<string,string> opts;
      for ( int i = 1, N = anaopt.size(); i < N; ++i ) {
        vector<string> opt = split(anaopt[i], "=");
        if ( opt.size() != 2 ) {
          MSG_WARNING("Error in option specification. Skipping analysis " << analysisname);
          return *this;
        }
        if ( !analysis->info().validOption(opt[0], opt[1]) )
          MSG_WARNING("Setting the option '" << opt[0] << "' to '"
                      << opt[1] << "' for " << analysisname
                      << " has not been declared in the info file "
                      << " and may be ignored in the analysis.");
        opts[opt[0]] = opt[1];
      }
      for ( auto opt: opts) {
        analysis->_options[opt.first] = opt.second;
        analysis->_optstring += ":" + opt.first + "=" + opt.second;
      }
      for (const AnaHandle& a : analyses()) {
        if (a->name() == analysis->name() ) {
          MSG_WARNING("Analysis '" << analysisname << "' already registered: skipping duplicate");
          return *this;
        }
      }
      analysis->_analysishandler = this;
      _analyses[analysisname] = analysis;
    } else {
      MSG_WARNING("Analysis '" << analysisname << "' not found.");
    }
    // MSG_WARNING(_analyses.size());
    // for (const AnaHandle& a : _analyses) MSG_WARNING(a->name());
    return *this;
  }


  AnalysisHandler& AnalysisHandler::removeAnalysis(const string& analysisname) {
    MSG_DEBUG("Removing analysis '" << analysisname << "'");
    if (_analyses.find(analysisname) != _analyses.end()) _analyses.erase(analysisname);
    // }
    return *this;
  }


  void AnalysisHandler::stripOptions(YODA::AnalysisObjectPtr ao,
                                     const vector<string> & delopts) const {
    string path = ao->path();
    string ananame = split(path, "/")[0];
    vector<string> anaopts = split(ananame, ":");
    for ( int i = 1, N = anaopts.size(); i < N; ++i )
      for ( auto opt : delopts )
        if ( opt == "*" || anaopts[i].find(opt + "=") == 0 )
          path.replace(path.find(":" + anaopts[i]), (":" + anaopts[i]).length(), "");
    ao->setPath(path);
  }


  void AnalysisHandler::mergeYodas(const vector<string> & aofiles,
                                   const vector<string> & delopts, 
                                   const vector<string> & addopts,
                                   bool equiv) {

    // Convenience typedef
    typedef multimap<string, YODA::AnalysisObjectPtr> AOMap;

    // Store all found weights here
    set<string> foundWeightNames;

    // Store all found analyses
    set<string> foundAnalyses;

    // Store all analysis objects here
    vector<AOMap> allaos;

    // Parse option adding.
    vector<string> optAnas;
    vector<string> optKeys;
    vector<string> optVals;
    for (string addopt : addopts) {
      size_t pos1 = addopt.find(":");
      size_t pos2 = addopt.find("=");
      if (pos1 == string::npos || pos2 == string::npos || pos2 < pos1) {
        MSG_WARNING("Malformed analysis option: "+addopt+". Format as ANA:OPT=VAL");
        continue;
      }
      optAnas.push_back(addopt.substr(0, pos1));
      optKeys.push_back(addopt.substr(pos1 +1, pos2 - pos1 - 1));
      optVals.push_back(addopt.substr(pos2 +1 , addopt.size() - pos2 - 1));
    }

    // Go through all files and collect information
    for ( auto file : aofiles ) {
      allaos.push_back(AOMap());
      AOMap & aomap = allaos.back();
      vector<YODA::AnalysisObject*> aos_raw;
      try {
        YODA::read(file, aos_raw);
      }
      catch (...) { //< YODA::ReadError&
        throw UserError("Unexpected error in reading file: " + file);
      }
      for (YODA::AnalysisObject* aor : aos_raw) {
        YODA::AnalysisObjectPtr ao(aor);
        AOPath path(ao->path());
        if ( !path )
          throw UserError("Invalid path name in file: " + file);
        if ( !path.isRaw() ) continue;

        foundWeightNames.insert(path.weight());
        // Now check if any options should be removed
        for ( string delopt : delopts )
          if ( path.hasOption(delopt) ) path.removeOption(delopt);
        // ...or added
        for (size_t i = 0; i < optAnas.size(); ++i) {
          if (path.path().find(optAnas[i]) != string::npos ) {
            path.setOption(optKeys[i], optVals[i]);
            path.fixOptionString();   
          }
        }
        path.setPath();
        if ( path.analysisWithOptions() != "" )
          foundAnalyses.insert(path.analysisWithOptions());
        aomap.insert(make_pair(path.path(), ao));
      }
    }

    // Now make analysis handler aware of the weight names present
    _weightNames.clear();
    _rivetDefaultWeightIdx = _defaultWeightIdx = 0;
    for ( string name : foundWeightNames ) _weightNames.push_back(name);

    // Then we create and initialize all analyses
    for ( string ananame : foundAnalyses ) addAnalysis(ananame);
    _stage = Stage::INIT;
    for (AnaHandle a : analyses() ) {
      MSG_TRACE("Initialising analysis: " << a->name());
      if ( !a->info().reentrant() )
        MSG_WARNING("Analysis " << a->name() << " has not been validated to have "
                    << "a reentrant finalize method. The merged result is unpredictable.");
      try {
        // Allow projection registration in the init phase onwards
        a->_allowProjReg = true;
        a->init();
      } catch (const Error& err) {
        cerr << "Error in " << a->name() << "::init method: " << err.what() << endl;
        exit(1);
      }
      MSG_TRACE("Done initialising analysis: " << a->name());
    }
    _stage = Stage::OTHER;
    _initialised = true;

    // Now get all booked analysis objects
    vector<MultiweightAOPtr> raos;
    for (AnaHandle a : analyses()) {
      for (const auto & ao : a->analysisObjects()) {
        raos.push_back(ao);
      }
    }

    // Collect global weights and cross sections and fix scaling for all files
    _eventCounter = CounterPtr(weightNames(), Counter("_EVTCOUNT"));
    _xs = Scatter1DPtr(weightNames(), Scatter1D("_XSEC"));
    for (size_t iW = 0; iW < numWeights(); iW++) {
      _eventCounter.get()->setActiveWeightIdx(iW);
      _xs.get()->setActiveWeightIdx(iW);
      YODA::Counter & sumw = *_eventCounter;
      YODA::Scatter1D & xsec = *_xs;
      vector<YODA::Scatter1DPtr> xsecs;
      vector<YODA::CounterPtr> sows;
      for ( auto & aomap : allaos ) {
        auto xit = aomap.find(xsec.path());
        if ( xit != aomap.end() )
          xsecs.push_back(dynamic_pointer_cast<YODA::Scatter1D>(xit->second));
        else
          xsecs.push_back(YODA::Scatter1DPtr());
        xit = aomap.find(sumw.path());
        if ( xit != aomap.end() )
          sows.push_back(dynamic_pointer_cast<YODA::Counter>(xit->second));
        else
          sows.push_back(YODA::CounterPtr());
      }
      double xs = 0.0, xserr = 0.0;
      for ( int i = 0, N = sows.size(); i < N; ++i ) {
        if ( !sows[i] || !xsecs[i] ) continue;
        double xseci = xsecs[i]->point(0).x();
        double xsecerri = sqr(xsecs[i]->point(0).xErrAvg());
        sumw += *sows[i];
        double effnent = sows[i]->effNumEntries();
        xs += (equiv? effnent: 1.0)*xseci;
        xserr += (equiv? sqr(effnent): 1.0)*xsecerri;
      }
      vector<double> scales(sows.size(), 1.0);
      if ( equiv ) {
        xs /= sumw.effNumEntries();
        xserr = sqrt(xserr)/sumw.effNumEntries();
      } else {
        xserr = sqrt(xserr);
        for ( int i = 0, N = sows.size(); i < N; ++i )
          scales[i] = (sumw.sumW()/sows[i]->sumW())*
           (xsecs[i]->point(0).x()/xs);
      }
      xsec.reset();
      xsec.addPoint(Point1D(xs, xserr));

      // Go through alla analyses and add stuff to their analysis objects;
      for (AnaHandle a : analyses()) {
        for (const auto & ao : a->analysisObjects()) {
          ao.get()->setActiveWeightIdx(iW);
          YODA::AnalysisObjectPtr yao = ao.get()->activeYODAPtr();
          for ( int i = 0, N = sows.size(); i < N; ++i ) {
            if ( !sows[i] || !xsecs[i] ) continue;
            auto range = allaos[i].equal_range(yao->path());
            for ( auto aoit = range.first; aoit != range.second; ++aoit ) {
              if ( !addaos(yao, aoit->second, scales[i]) )
                MSG_WARNING("Cannot merge objects with path " << yao->path()
                            <<" of type " << yao->annotation("Type") );
            }
          }
	  a->rawHookIn(yao);
          ao.get()->unsetActiveWeight();
        }
      }
      _eventCounter.get()->unsetActiveWeight();
      _xs.get()->unsetActiveWeight();
    }

    // Finally we just have to finalize all analyses, leaving to the
    // controlling program to write it out some yoda-file.
    finalize();

  }


  void AnalysisHandler::readData(const string& filename) {
    try {
      /// @todo Use new YODA SFINAE to fill the smart ptr vector directly
      vector<YODA::AnalysisObject*> aos_raw;
      YODA::read(filename, aos_raw);
      for (YODA::AnalysisObject* aor : aos_raw)
        _preloads[aor->path()] = YODA::AnalysisObjectPtr(aor);
    } catch (...) { //< YODA::ReadError&
      throw UserError("Unexpected error in reading file: " + filename);
    }
  }


  vector<MultiweightAOPtr> AnalysisHandler::getRivetAOs() const {
      vector<MultiweightAOPtr> rtn;

      for (AnaHandle a : analyses()) {
          for (const auto & ao : a->analysisObjects()) {
              rtn.push_back(ao);
          }
      }
      rtn.push_back(_eventCounter);
      rtn.push_back(_xs);
      return rtn;
  }


  vector<YODA::AnalysisObjectPtr> AnalysisHandler::getYodaAOs(bool includeraw) const {
    vector<YODA::AnalysisObjectPtr> output;

    // First get all multiweight AOs
    vector<MultiweightAOPtr> raos = getRivetAOs();
    output.reserve(raos.size() * numWeights() * (includeraw ? 2 : 1));

    // Identify an index ordering so that default weight is written out first
    vector<size_t> order = { _rivetDefaultWeightIdx };
    for ( size_t  i = 0; i < numWeights(); ++i ) {
      if ( i != _rivetDefaultWeightIdx )  order.push_back(i);
    }

    // Then we go through all finalized AOs one weight at a time
    for (size_t iW : order ) {
      for ( auto rao : raos ) {
        rao.get()->setActiveFinalWeightIdx(iW);
        if ( rao->path().find("/TMP/") != string::npos ) continue;
        output.push_back(rao.get()->activeYODAPtr());
      }
    }

    // Analyses can make changes neccesary for merging to RAW objects
    // before writing.
    for (size_t iW : order)
      for (auto a : analyses()) a->rawHookOut(raos, iW);
    
    // Finally the RAW objects.
    if (includeraw) {
      for (size_t iW : order ) {
        for ( auto rao : raos ) {
          rao.get()->setActiveWeightIdx(iW);
          output.push_back(rao.get()->activeYODAPtr());
        }
      }
    }

    return output;
  }


  void AnalysisHandler::writeData(const string& filename) const {

    const vector<YODA::AnalysisObjectPtr> output = getYodaAOs(true);
    try {
      YODA::write(filename, output.begin(), output.end());
    } catch (...) { //< YODA::WriteError&
      throw UserError("Unexpected error in writing file: " + filename);
    }

  }


  string AnalysisHandler::runName() const {
    return _runname;
  }


  size_t AnalysisHandler::numEvents() const {
    return _eventCounter->numEntries();
  }


  std::vector<std::string> AnalysisHandler::analysisNames() const {
    std::vector<std::string> rtn;
    for (AnaHandle a : analyses()) {
      rtn.push_back(a->name());
    }
    return rtn;
  }


  std::vector<std::string> AnalysisHandler::stdAnalysisNames() const {
    // std::vector<std::string> rtn;
    // const string anadatpath = findAnalysisDataFile("analyses.dat");
    // if (fileexists(anadatpath)) {
    //   std::ifstream anadat(anadatpath);
    //   string ananame;
    //   while (anadat >> ananame) rtn += ananame;
    // }
    // return rtn;
    return AnalysisLoader::stdAnalysisNames();
  }


  AnalysisHandler& AnalysisHandler::addAnalyses(const std::vector<std::string>& analysisnames) {
    for (const string& aname : analysisnames) {
      //MSG_DEBUG("Adding analysis '" << aname << "'");
      addAnalysis(aname);
    }
    return *this;
  }


  AnalysisHandler& AnalysisHandler::removeAnalyses(const std::vector<std::string>& analysisnames) {
    for (const string& aname : analysisnames) removeAnalysis(aname);
    return *this;
  }


  void AnalysisHandler::setCrossSection(pair<double,double> xsec, bool isUserSupplied) {
    // Update the user xsec
    if (isUserSupplied) _userxs = xsec;

    // If not setting the user xsec, and a user xsec is already set, exit early
    if (!isUserSupplied && notNaN(_userxs.first)) return;

    // Otherwise, update the xs scatter: xs_var = xs_nom * (sumW_var/sumW_nom)
    _xs = Scatter1DPtr(weightNames(), Scatter1D("_XSEC"));
    _eventCounter.get()->setActiveWeightIdx(_rivetDefaultWeightIdx);
    const double nomwgt = sumW();
    for (size_t iW = 0; iW < numWeights(); ++iW) {
      _eventCounter.get()->setActiveWeightIdx(iW);
      const double s = sumW() / nomwgt;
      _xs.get()->setActiveWeightIdx(iW);
      _xs->addPoint(xsec.first*s, xsec.second*s);
    }
    _eventCounter.get()->unsetActiveWeight();
    _xs.get()->unsetActiveWeight();
  }


  AnalysisHandler& AnalysisHandler::addAnalysis(Analysis* analysis) {
    analysis->_analysishandler = this;
    // _analyses.insert(AnaHandle(analysis));
    _analyses[analysis->name()] = AnaHandle(analysis);
    return *this;
  }


  PdgIdPair AnalysisHandler::beamIds() const {
    return Rivet::beamIds(beams());
  }


  double AnalysisHandler::sqrtS() const {
    return Rivet::sqrtS(beams());
  }


  void AnalysisHandler::setIgnoreBeams(bool ignore) {
    _ignoreBeams=ignore;
  }


  void AnalysisHandler::skipMultiWeights(bool ignore) {
    _skipWeights = ignore;
  }

  void AnalysisHandler::selectMultiWeights(std::string patterns) {
    _matchWeightNames = patterns;
  }

  void AnalysisHandler::deselectMultiWeights(std::string patterns) {
    _unmatchWeightNames = patterns;
  }


  std::valarray<double> AnalysisHandler::pruneWeights(const std::valarray<double>& weights) {
    if (_weightIndices.size() == weights.size())  return weights;
    std::valarray<double> acceptedWeights(_weightIndices.size());
    for (size_t i = 0; i < _weightIndices.size(); ++i) {
      acceptedWeights[i] = weights[_weightIndices[i]];
    }
    return acceptedWeights;
  }

}
