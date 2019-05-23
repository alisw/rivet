// -*- C++ -*-
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/ParticleName.hh"
#include "Rivet/Tools/BeamConstraint.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/Beam.hh"
#include "YODA/IO.h"

namespace Rivet {


  AnalysisHandler::AnalysisHandler(const string& runname)
    : _runname(runname),
      _eventcounter("/_EVTCOUNT"),
      _xs(NAN), _xserr(NAN),
      _initialised(false), _ignoreBeams(false), _dumpPeriod(0), _dumping(0)
  {  }


  AnalysisHandler::~AnalysisHandler()
  {  }


  Log& AnalysisHandler::getLog() const {
    return Log::getLog("Rivet.Analysis.Handler");
  }


  void AnalysisHandler::init(const GenEvent& ge) {
    if (_initialised)
      throw UserError("AnalysisHandler::init has already been called: cannot re-initialize!");

    setRunBeams(Rivet::beams(ge));
    MSG_DEBUG("Initialising the analysis handler");
    _eventcounter.reset();

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
      if (toUpper(a->status()) == "PRELIMINARY") {
        MSG_WARNING("Analysis '" << a->name() << "' is preliminary: be careful, it may change and/or be renamed!");
      } else if (toUpper(a->status()) == "OBSOLETE") {
        MSG_WARNING("Analysis '" << a->name() << "' is obsolete: please update!");
      } else if (toUpper(a->status()).find("UNVALIDATED") != string::npos) {
        MSG_WARNING("Analysis '" << a->name() << "' is unvalidated: be careful, it may be broken!");
      }
    }

    // Initialize the remaining analyses
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
    _initialised = true;
    MSG_DEBUG("Analysis handler initialised");
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
    Event event(ge);

    // Weights
    /// @todo Drop this / just report first weight when we support multiweight events
    _eventcounter.fill(event.weight());
    MSG_DEBUG("Event #" << _eventcounter.numEntries() << " weight = " << event.weight());

    // Cross-section
    #ifdef HEPMC_HAS_CROSS_SECTION
    if (ge.cross_section()) {
      _xs = ge.cross_section()->cross_section();
      _xserr = ge.cross_section()->cross_section_error();
    }
    #endif

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

    if ( _dumpPeriod > 0 && numEvents()%_dumpPeriod == 0 ) {
      MSG_INFO("Dumping intermediate results to " << _dumpFile << ".");
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


  void AnalysisHandler::finalize() {
    if (!_initialised) return;

    // First we make copies of all analysis objects.
    map<string,AnalysisObjectPtr> backupAOs;
    for (auto ao : getData(false, true, false) )
      backupAOs[ao->path()] = AnalysisObjectPtr(ao->newclone());

    // Now we run the (re-entrant) finalize() functions for all analyses.
    MSG_INFO("Finalising analyses");
    for (AnaHandle a : analyses()) {
      a->setCrossSection(_xs);
      try {
        if ( !_dumping || a->info().reentrant() )  a->finalize();
        else if ( _dumping == 1 )
          MSG_INFO("Skipping finalize in periodic dump of " << a->name()
                   << " as it is not declared reentrant.");
      } catch (const Error& err) {
        cerr << "Error in " << a->name() << "::finalize method: " << err.what() << endl;
        exit(1);
      }
    }

    // Now we copy all analysis objects to the list of finalized
    // ones, and restore the value to their original ones.
    _finalizedAOs.clear();
    for ( auto ao : getData(false, false, false) )
      _finalizedAOs.push_back(AnalysisObjectPtr(ao->newclone()));
    for ( auto ao : getData(false, true, false) ) {
      // TODO: This should be possible to do in a nicer way, with a flag etc.
      if (ao->path().find("/FINAL") != std::string::npos) continue;
      auto aoit = backupAOs.find(ao->path());
      if ( aoit == backupAOs.end() ) {
        AnaHandle ana = analysis(split(ao->path(), "/")[0]);
        if ( ana ) ana->removeAnalysisObject(ao->path());
      } else
        copyao(aoit->second, ao);
    }

    // Print out number of events processed
    const int nevts = _eventcounter.numEntries();
    MSG_INFO("Processed " << nevts << " event" << (nevts != 1 ? "s" : ""));

    // // Delete analyses
    // MSG_DEBUG("Deleting analyses");
    // _analyses.clear();

    // Print out MCnet boilerplate
    cout << endl;
    cout << "The MCnet usage guidelines apply to Rivet: see http://www.montecarlonet.org/GUIDELINES" << endl;
    cout << "Please acknowledge plots made with Rivet analyses, and cite arXiv:1003.0694 (http://arxiv.org/abs/1003.0694)" << endl;
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
        if ( !analysis->info().validOption(opt[0], opt[1]) ) {
          MSG_WARNING("Cannot set option '" << opt[0] << "' to '" << opt[1]
                      << "'. Skipping analysis " << analysisname);
          return *this;
        }
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
    // std::shared_ptr<Analysis> toremove;
    // for (const AnaHandle a : _analyses) {
    //   if (a->name() == analysisname) {
    //     toremove = a;
    //     break;
    //   }
    // }
    // if (toremove.get() != 0) {
    MSG_DEBUG("Removing analysis '" << analysisname << "'");
    if (_analyses.find(analysisname) != _analyses.end()) _analyses.erase(analysisname);
    // }
    return *this;
  }


  /////////////////////////////


  void AnalysisHandler::addData(const std::vector<AnalysisObjectPtr>& aos) {
    for (const AnalysisObjectPtr ao : aos) {
      string path = ao->path();
      if ( path.substr(0, 5) != "/RAW/" ) {
        _orphanedPreloads.push_back(ao);
        continue;
      }

      path = path.substr(4);
      ao->setPath(path);
      if (path.size() > 1) { // path > "/"
        try {
          const string ananame =  split(path, "/")[0];
          AnaHandle a = analysis(ananame);
          a->addAnalysisObject(ao); /// @todo Need to statistically merge...
        } catch (const Error& e) {
          MSG_TRACE("Adding analysis object " << path <<
                    " to the list of orphans.");
          _orphanedPreloads.push_back(ao);
        }
      }
    }
  }

  void AnalysisHandler::stripOptions(AnalysisObjectPtr ao,
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




  void AnalysisHandler::
  mergeYodas(const vector<string> & aofiles, const vector<string> & delopts, bool equiv) {
    vector< vector<AnalysisObjectPtr> > aosv;
    vector<double> xsecs;
    vector<double> xsecerrs;
    vector<CounterPtr> sows;
    set<string> ananames;
     _eventcounter.reset();

    // First scan all files and extract analysis objects and add the
    // corresponding anayses..
    for ( auto file : aofiles ) {
      Scatter1DPtr xsec;
      CounterPtr sow;

      // For each file make sure that cross section and sum-of-weights
      // objects are present and stor all RAW ones in a vector;
      vector<AnalysisObjectPtr> aos;
      try {
        /// @todo Use new YODA SFINAE to fill the smart ptr vector directly
        vector<YODA::AnalysisObject*> aos_raw;
        YODA::read(file, aos_raw);
        for (AnalysisObject* aor : aos_raw) {
          AnalysisObjectPtr ao = AnalysisObjectPtr(aor);
          if ( ao->path().substr(0, 5) != "/RAW/" ) continue;
          ao->setPath(ao->path().substr(4));
          if ( ao->path() == "/_XSEC" )
            xsec = dynamic_pointer_cast<Scatter1D>(ao);
          else if ( ao->path() == "/_EVTCOUNT" )
            sow = dynamic_pointer_cast<Counter>(ao);
          else {
            stripOptions(ao, delopts);
            string ananame = split(ao->path(), "/")[0];
            if ( ananames.insert(ananame).second ) addAnalysis(ananame);
            aos.push_back(ao);
          }
        }
        if ( !xsec || !sow ) {
          MSG_ERROR( "Error in AnalysisHandler::mergeYodas: The file " << file
                     << " did not contain weights and cross section info.");
          exit(1);
        }
        xsecs.push_back(xsec->point(0).x());
	xsecerrs.push_back(sqr(xsec->point(0).xErrAvg()));
        _eventcounter += *sow;
        sows.push_back(sow);
        aosv.push_back(aos);
      } catch (...) { //< YODA::ReadError&
        throw UserError("Unexpected error in reading file: " + file);
      }
    }

    // Now calculate the scale to be applied for all bins in a file
    // and get the common cross section and sum of weights.
    _xs = _xserr = 0.0;
    for ( int i = 0, N = sows.size(); i < N; ++i ) {
      double effnent = sows[i]->effNumEntries();
      _xs += (equiv? effnent: 1.0)*xsecs[i];
      _xserr += (equiv? sqr(effnent): 1.0)*xsecerrs[i];
    }

    vector<double> scales(sows.size(), 1.0);
    if ( equiv ) {
      _xs /= _eventcounter.effNumEntries();
      _xserr = sqrt(_xserr)/_eventcounter.effNumEntries();
    } else {
      _xserr = sqrt(_xserr);
      for ( int i = 0, N = sows.size(); i < N; ++i )
        scales[i] = (_eventcounter.sumW()/sows[i]->sumW())*(xsecs[i]/_xs);
    }

    // Initialize the analyses allowing them to book analysis objects.
    for (AnaHandle a : analyses()) {
      MSG_DEBUG("Initialising analysis: " << a->name());
      if ( !a->info().reentrant() )
        MSG_WARNING("Analysis " << a->name() << " has not been validated to have "
                    << "a reentrant finalize method. The result is unpredictable.");
      try {
        // Allow projection registration in the init phase onwards
        a->_allowProjReg = true;
        cerr << "sqrtS " << sqrtS() << endl;
        a->init();
        //MSG_DEBUG("Checking consistency of analysis: " << a->name());
        //a->checkConsistency();
      } catch (const Error& err) {
        cerr << "Error in " << a->name() << "::init method: " << err.what() << endl;
        exit(1);
      }
      MSG_DEBUG("Done initialising analysis: " << a->name());
    }
    _initialised = true;
    // Get a list of all anaysis objects to handle.
    map<string,AnalysisObjectPtr> current;
    for ( auto ao : getData(false, true, false) ) current[ao->path()] = ao;
    // Go through all objects to be merged and add them to current
    // after appropriate scaling.
    for ( int i = 0, N = aosv.size(); i < N; ++i)
      for ( auto ao : aosv[i] ) {
        if ( ao->path() == "/_XSEC" || ao->path() == "_EVTCOUNT" ) continue;
	auto aoit = current.find(ao->path());
        if ( aoit == current.end() ) {
          MSG_WARNING("" << ao->path() << " was not properly booked.");
          continue;
        }
        if ( !addaos(aoit->second, ao, scales[i]) )
          MSG_WARNING("Cannot merge objects with path " << ao->path()
                      <<" of type " << ao->annotation("Type") );
      }
    // Now we can simply finalize() the analysis, leaving the
    // controlling program to write it out some yoda-file.
    finalize();

  }


  void AnalysisHandler::readData(const string& filename) {
    vector<AnalysisObjectPtr> aos;
    try {
      /// @todo Use new YODA SFINAE to fill the smart ptr vector directly
      vector<YODA::AnalysisObject*> aos_raw;
      YODA::read(filename, aos_raw);
      for (AnalysisObject* aor : aos_raw) aos.push_back(AnalysisObjectPtr(aor));
    } catch (...) { //< YODA::ReadError&
      throw UserError("Unexpected error in reading file: " + filename);
    }
    if (!aos.empty()) addData(aos);
  }


  vector<AnalysisObjectPtr> AnalysisHandler::
  getData(bool includeorphans, bool includetmps, bool usefinalized) const {
    vector<AnalysisObjectPtr> rtn;
    // Event counter
    rtn.push_back( make_shared<Counter>(_eventcounter) );
    // Cross-section + err as scatter
    YODA::Scatter1D::Points pts; pts.insert(YODA::Point1D(_xs, _xserr));
    rtn.push_back( make_shared<Scatter1D>(pts, "/_XSEC") );
    // Analysis histograms
    vector<AnalysisObjectPtr> aos;
    if (usefinalized)
      aos = _finalizedAOs;
    else {
      for (const AnaHandle a : analyses()) {
        // MSG_WARNING(a->name() << " " << aos.size());
        for (const AnalysisObjectPtr ao : a->analysisObjects()) {
          aos.push_back(ao);
        }
      }
    }
    for (const AnalysisObjectPtr ao : aos) {
      // Exclude paths from final write-out if they contain a "TMP" layer (i.e. matching "/TMP/")
      /// @todo This needs to be much more nuanced for re-entrant histogramming
      if ( !includetmps && ao->path().find("/TMP/" ) != string::npos) continue;
      rtn.push_back(ao);
    }
    // Sort histograms alphanumerically by path before write-out
    sort(rtn.begin(), rtn.end(), [](AnalysisObjectPtr a, AnalysisObjectPtr b) {return a->path() < b->path();});
    if ( includeorphans )
      rtn.insert(rtn.end(), _orphanedPreloads.begin(), _orphanedPreloads.end());
    return rtn;
  }


  void AnalysisHandler::writeData(const string& filename) const {
    vector<AnalysisObjectPtr> out = _finalizedAOs;
    set<string> finalana;
    for ( auto ao : out) finalana.insert(ao->path());
    out.reserve(2*out.size());
    vector<AnalysisObjectPtr> aos = getData(false, true, false);

    if ( _dumping ) {
      for ( auto ao : aos ) {
        if ( finalana.find(ao->path()) == finalana.end() )
          out.push_back(AnalysisObjectPtr(ao->newclone()));
      }
    }

    for ( auto ao : aos ) {
      ao = AnalysisObjectPtr(ao->newclone());
      ao->setPath("/RAW" + ao->path());
      out.push_back(ao);
    }

    try {
      YODA::write(filename, out.begin(), out.end());
    } catch (...) { //< YODA::WriteError&
      throw UserError("Unexpected error in writing file: " + filename);
    }
  }


  std::vector<std::string> AnalysisHandler::analysisNames() const {
    std::vector<std::string> rtn;
    for (AnaHandle a : analyses()) {
      rtn.push_back(a->name());
    }
    return rtn;
  }


  AnalysisHandler& AnalysisHandler::addAnalyses(const std::vector<std::string>& analysisnames) {
    for (const string& aname : analysisnames) {
      //MSG_DEBUG("Adding analysis '" << aname << "'");
      addAnalysis(aname);
    }
    return *this;
  }


  AnalysisHandler& AnalysisHandler::removeAnalyses(const std::vector<std::string>& analysisnames) {
    for (const string& aname : analysisnames) {
      removeAnalysis(aname);
    }
    return *this;
  }


  bool AnalysisHandler::needCrossSection() const {
    bool rtn = false;
    for (const AnaHandle a : analyses()) {
      if (!rtn) rtn = a->needsCrossSection();
      if (rtn) break;
    }
    return rtn;
  }


  AnalysisHandler& AnalysisHandler::setCrossSection(double xs, double xserr) {
    _xs = xs;
    _xserr = xserr;
    return *this;
  }


  bool AnalysisHandler::hasCrossSection() const {
    return (!std::isnan(crossSection()));
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


}
