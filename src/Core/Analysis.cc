// -*- C++ -*-
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/AnalysisInfo.hh"
#include "Rivet/Tools/BeamConstraint.hh"
#include "Rivet/Projections/ImpactParameterProjection.hh"
#include "Rivet/Projections/GeneratedPercentileProjection.hh"
#include "Rivet/Projections/UserCentEstimate.hh"
#include "Rivet/Projections/CentralityProjection.hh"

namespace Rivet {


  Analysis::Analysis(const string& name)
    : _crossSection(-1.0),
      _gotCrossSection(false),
      _analysishandler(NULL)
  {
    ProjectionApplier::_allowProjReg = false;
    _defaultname = name;

    unique_ptr<AnalysisInfo> ai = AnalysisInfo::make(name);
    assert(ai);
    _info = move(ai);
    assert(_info);
  }

  double Analysis::sqrtS() const {
    return handler().sqrtS();
  }

  const ParticlePair& Analysis::beams() const {
    return handler().beams();
  }

  const PdgIdPair Analysis::beamIds() const {
    return handler().beamIds();
  }


  const string Analysis::histoDir() const {
    /// @todo Cache in a member variable
    string _histoDir;
    if (_histoDir.empty()) {
      _histoDir = "/" + name();
      if (handler().runName().length() > 0) {
        _histoDir = "/" + handler().runName() + _histoDir;
      }
      replace_all(_histoDir, "//", "/"); //< iterates until none
    }
    return _histoDir;
  }


  const string Analysis::histoPath(const string& hname) const {
    const string path = histoDir() + "/" + hname;
    return path;
  }


  const string Analysis::histoPath(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) const {
    return histoDir() + "/" + mkAxisCode(datasetId, xAxisId, yAxisId);
  }


  const string Analysis::mkAxisCode(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) const {
    stringstream axisCode;
    axisCode << "d";
    if (datasetId < 10) axisCode << 0;
    axisCode << datasetId;
    axisCode << "-x";
    if (xAxisId < 10) axisCode << 0;
    axisCode << xAxisId;
    axisCode << "-y";
    if (yAxisId < 10) axisCode << 0;
    axisCode << yAxisId;
    return axisCode.str();
  }


  Log& Analysis::getLog() const {
    string logname = "Rivet.Analysis." + name();
    return Log::getLog(logname);
  }


  ///////////////////////////////////////////


  size_t Analysis::numEvents() const {
    return handler().numEvents();
  }

  double Analysis::sumW() const {
    return handler().sumW();
  }

  double Analysis::sumW2() const {
    return handler().sumW2();
  }


  ///////////////////////////////////////////


  bool Analysis::isCompatible(const ParticlePair& beams) const {
    return isCompatible(beams.first.pid(),  beams.second.pid(),
                        beams.first.energy(), beams.second.energy());
  }


  bool Analysis::isCompatible(PdgId beam1, PdgId beam2, double e1, double e2) const {
    PdgIdPair beams(beam1, beam2);
    pair<double,double> energies(e1, e2);
    return isCompatible(beams, energies);
  }


  // bool Analysis::beamIDsCompatible(const PdgIdPair& beams) const {
  //   bool beamIdsOk = false;
  //   for (const PdgIdPair& bp : requiredBeams()) {
  //     if (compatible(beams, bp)) {
  //       beamIdsOk =  true;
  //       break;
  //     }
  //   }
  //   return beamIdsOk;
  // }


  // /// Check that the energies are compatible (within 1% or 1 GeV, whichever is larger, for a bit of UI forgiveness)
  // bool Analysis::beamEnergiesCompatible(const pair<double,double>& energies) const {
  //   /// @todo Use some sort of standard ordering to improve comparisons, esp. when the two beams are different particles
  //   bool beamEnergiesOk = requiredEnergies().size() > 0 ? false : true;
  //   typedef pair<double,double> DoublePair;
  //   for (const DoublePair& ep : requiredEnergies()) {
  //     if ((fuzzyEquals(ep.first, energies.first, 0.01) && fuzzyEquals(ep.second, energies.second, 0.01)) ||
  //         (fuzzyEquals(ep.first, energies.second, 0.01) && fuzzyEquals(ep.second, energies.first, 0.01)) ||
  //         (abs(ep.first - energies.first) < 1*GeV && abs(ep.second - energies.second) < 1*GeV) ||
  //         (abs(ep.first - energies.second) < 1*GeV && abs(ep.second - energies.first) < 1*GeV)) {
  //       beamEnergiesOk =  true;
  //       break;
  //     }
  //   }
  //   return beamEnergiesOk;
  // }


  // bool Analysis::beamsCompatible(const PdgIdPair& beams, const pair<double,double>& energies) const {
  bool Analysis::isCompatible(const PdgIdPair& beams, const pair<double,double>& energies) const {
    // First check the beam IDs
    bool beamIdsOk = false;
    for (const PdgIdPair& bp : requiredBeams()) {
      if (compatible(beams, bp)) {
        beamIdsOk =  true;
        break;
      }
    }
    if (!beamIdsOk) return false;
    // Next check that the energies are compatible (within 1% or 1 GeV, whichever is larger, for a bit of UI forgiveness)
    
    /// @todo Use some sort of standard ordering to improve comparisons, esp. when the two beams are different particles
    bool beamEnergiesOk = requiredEnergies().size() > 0 ? false : true;
    typedef pair<double,double> DoublePair;
    for (const DoublePair& ep : requiredEnergies()) {
      if ((fuzzyEquals(ep.first, energies.first, 0.01) && fuzzyEquals(ep.second, energies.second, 0.01)) ||
          (fuzzyEquals(ep.first, energies.second, 0.01) && fuzzyEquals(ep.second, energies.first, 0.01)) ||
          (abs(ep.first - energies.first) < 1*GeV && abs(ep.second - energies.second) < 1*GeV) ||
          (abs(ep.first - energies.second) < 1*GeV && abs(ep.second - energies.first) < 1*GeV)) {
        beamEnergiesOk =  true;
        break;
      }
    }
    return beamEnergiesOk;
  }


  ///////////////////////////////////////////


  Analysis& Analysis::setCrossSection(double xs) {
    _crossSection = xs;
    _gotCrossSection = true;
    return *this;
  }

  double Analysis::crossSection() const {
    if (!_gotCrossSection || std::isnan(_crossSection)) {
      string errMsg = "You did not set the cross section for the analysis " + name();
      throw Error(errMsg);
    }
    return _crossSection;
  }

  double Analysis::crossSectionPerEvent() const {
    const double sumW = sumOfWeights();
    assert(sumW != 0.0);
    return _crossSection / sumW;
  }



  ////////////////////////////////////////////////////////////
  // Histogramming


  void Analysis::_cacheRefData() const {
    if (_refdata.empty()) {
      MSG_TRACE("Getting refdata cache for paper " << name());
      _refdata = getRefData(getRefDataName());
    }
  }

  vector<AnalysisObjectPtr> Analysis::getAllData(bool includeorphans) const{
    return handler().getData(includeorphans, false, false);
  }

  CounterPtr Analysis::bookCounter(const string& cname,
                                   const string& title) {
    return addOrGetCompatAO(make_shared<Counter>(histoPath(cname), title));
  }


  CounterPtr Analysis::bookCounter(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId,
                                   const string& title) {
    return bookCounter(mkAxisCode(datasetId, xAxisId, yAxisId), title);
  }


  Histo1DPtr Analysis::bookHisto1D(const string& hname,
                                   size_t nbins, double lower, double upper,
                                   const string& title,
                                   const string& xtitle,
                                   const string& ytitle) {
    Histo1DPtr hist = make_shared<Histo1D>(nbins, lower, upper, histoPath(hname), title);
    hist->setTitle(title);
    hist->setAnnotation("XLabel", xtitle);
    hist->setAnnotation("YLabel", ytitle);
    return addOrGetCompatAO(hist);
  }


  Histo1DPtr Analysis::bookHisto1D(const string& hname,
                                   const vector<double>& binedges,
                                   const string& title,
                                   const string& xtitle,
                                   const string& ytitle) {
    Histo1DPtr hist = make_shared<Histo1D>(binedges, histoPath(hname), title);
    hist->setTitle(title);
    hist->setAnnotation("XLabel", xtitle);
    hist->setAnnotation("YLabel", ytitle);
    return addOrGetCompatAO(hist);
  }


  Histo1DPtr Analysis::bookHisto1D(const string& hname,
                                   const initializer_list<double>& binedges,
                                   const string& title,
                                   const string& xtitle,
                                   const string& ytitle) {
    return bookHisto1D(hname, vector<double>{binedges}, title, xtitle, ytitle);
  }


  Histo1DPtr Analysis::bookHisto1D(const string& hname,
                                   const Scatter2D& refscatter,
                                   const string& title,
                                   const string& xtitle,
                                   const string& ytitle) {
    Histo1DPtr hist = make_shared<Histo1D>(refscatter, histoPath(hname));
    if (hist->hasAnnotation("IsRef")) hist->rmAnnotation("IsRef");
    if (hist->hasAnnotation("ErrorBreakdown")) hist->rmAnnotation("ErrorBreakdown");
    if (hist->hasAnnotation("Variations")) hist->rmAnnotation("Variations");
    hist->setTitle(title);
    hist->setAnnotation("XLabel", xtitle);
    hist->setAnnotation("YLabel", ytitle);
    return addOrGetCompatAO(hist);
  }


  Histo1DPtr Analysis::bookHisto1D(const string& hname,
                                   const string& title,
                                   const string& xtitle,
                                   const string& ytitle) {
    const Scatter2D& refdata = refData(hname);
    return bookHisto1D(hname, refdata, title, xtitle, ytitle);
  }


  Histo1DPtr Analysis::bookHisto1D(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId,
                                   const string& title,
                                   const string& xtitle,
                                   const string& ytitle) {
    const string axisCode = mkAxisCode(datasetId, xAxisId, yAxisId);
    return bookHisto1D(axisCode, title, xtitle, ytitle);
  }


  /// @todo Add booking methods which take a path, titles and *a reference Scatter from which to book*


  /////////////////


  Histo2DPtr Analysis::bookHisto2D(const string& hname,
                                   size_t nxbins, double xlower, double xupper,
                                   size_t nybins, double ylower, double yupper,
                                   const string& title,
                                   const string& xtitle,
                                   const string& ytitle,
                                   const string& ztitle)
  {
    const string path = histoPath(hname);
    Histo2DPtr hist = make_shared<Histo2D>(nxbins, xlower, xupper, nybins, ylower, yupper, path, title);
    hist->setAnnotation("XLabel", xtitle);
    hist->setAnnotation("YLabel", ytitle);
    hist->setAnnotation("ZLabel", ztitle);
    return addOrGetCompatAO(hist);
  }


  Histo2DPtr Analysis::bookHisto2D(const string& hname,
                                   const vector<double>& xbinedges,
                                   const vector<double>& ybinedges,
                                   const string& title,
                                   const string& xtitle,
                                   const string& ytitle,
                                   const string& ztitle)
  {
    const string path = histoPath(hname);
    Histo2DPtr hist = make_shared<Histo2D>(xbinedges, ybinedges, path, title);
    hist->setAnnotation("XLabel", xtitle);
    hist->setAnnotation("YLabel", ytitle);
    hist->setAnnotation("ZLabel", ztitle);
    return addOrGetCompatAO(hist);
  }


  Histo2DPtr Analysis::bookHisto2D(const string& hname,
                                   const initializer_list<double>& xbinedges,
                                   const initializer_list<double>& ybinedges,
                                   const string& title,
                                   const string& xtitle,
                                   const string& ytitle,
                                   const string& ztitle)
  {
    return bookHisto2D(hname, vector<double>{xbinedges}, vector<double>{ybinedges},
                       title, xtitle, ytitle, ztitle);
  }


  Histo2DPtr Analysis::bookHisto2D(const string& hname,
                                   const Scatter3D& refscatter,
                                   const string& title,
                                   const string& xtitle,
                                   const string& ytitle,
                                   const string& ztitle) {
    const string path = histoPath(hname);
    Histo2DPtr hist( new Histo2D(refscatter, path) );
    if (hist->hasAnnotation("IsRef")) hist->rmAnnotation("IsRef");
    if (hist->hasAnnotation("ErrorBreakdown")) hist->rmAnnotation("ErrorBreakdown");
    if (hist->hasAnnotation("Variations")) hist->rmAnnotation("Variations");
    hist->setTitle(title);
    hist->setAnnotation("XLabel", xtitle);
    hist->setAnnotation("YLabel", ytitle);
    hist->setAnnotation("ZLabel", ztitle);
    return addOrGetCompatAO(hist);
  }


  Histo2DPtr Analysis::bookHisto2D(const string& hname,
                                   const string& title,
                                   const string& xtitle,
                                   const string& ytitle,
                                   const string& ztitle) {
    const Scatter3D& refdata = refData<Scatter3D>(hname);
    return bookHisto2D(hname, refdata, title, xtitle, ytitle, ztitle);
  }


  Histo2DPtr Analysis::bookHisto2D(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId,
                                   const string& title,
                                   const string& xtitle,
                                   const string& ytitle,
                                   const string& ztitle) {
    const string axisCode = mkAxisCode(datasetId, xAxisId, yAxisId);
    return bookHisto2D(axisCode, title, xtitle, ytitle, ztitle);
  }


  /////////////////


  Profile1DPtr Analysis::bookProfile1D(const string& hname,
                                       size_t nbins, double lower, double upper,
                                       const string& title,
                                       const string& xtitle,
                                       const string& ytitle) {
    const string path = histoPath(hname);
    Profile1DPtr prof = make_shared<Profile1D>(nbins, lower, upper, path, title);
    prof->setAnnotation("XLabel", xtitle);
    prof->setAnnotation("YLabel", ytitle);
    return addOrGetCompatAO(prof);
  }


  Profile1DPtr Analysis::bookProfile1D(const string& hname,
                                       const vector<double>& binedges,
                                       const string& title,
                                       const string& xtitle,
                                       const string& ytitle) {
    const string path = histoPath(hname);
    Profile1DPtr prof = make_shared<Profile1D>(binedges, path, title);
    prof->setAnnotation("XLabel", xtitle);
    prof->setAnnotation("YLabel", ytitle);
    return addOrGetCompatAO(prof);
  }


  Profile1DPtr Analysis::bookProfile1D(const string& hname,
                                       const initializer_list<double>& binedges,
                                       const string& title,
                                       const string& xtitle,
                                       const string& ytitle)
  {
    return bookProfile1D(hname, vector<double>{binedges}, title, xtitle, ytitle);
  }


  Profile1DPtr Analysis::bookProfile1D(const string& hname,
                                       const Scatter2D& refscatter,
                                       const string& title,
                                       const string& xtitle,
                                       const string& ytitle) {
    const string path = histoPath(hname);
    Profile1DPtr prof = make_shared<Profile1D>(refscatter, path);
    if (prof->hasAnnotation("IsRef")) prof->rmAnnotation("IsRef");
    if (prof->hasAnnotation("ErrorBreakdown")) prof->rmAnnotation("ErrorBreakdown");
    if (prof->hasAnnotation("Variations")) prof->rmAnnotation("Variations");
    prof->setTitle(title);
    prof->setAnnotation("XLabel", xtitle);
    prof->setAnnotation("YLabel", ytitle);
    return addOrGetCompatAO(prof);
  }


  Profile1DPtr Analysis::bookProfile1D(const string& hname,
                                       const string& title,
                                       const string& xtitle,
                                       const string& ytitle) {
    const Scatter2D& refdata = refData(hname);
    return bookProfile1D(hname, refdata, title, xtitle, ytitle);
  }


  Profile1DPtr Analysis::bookProfile1D(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId,
                                       const string& title,
                                       const string& xtitle,
                                       const string& ytitle) {
    const string axisCode = mkAxisCode(datasetId, xAxisId, yAxisId);
    return bookProfile1D(axisCode, title, xtitle, ytitle);
  }


  ///////////////////



  Profile2DPtr Analysis::bookProfile2D(const string& hname,
                                   size_t nxbins, double xlower, double xupper,
                                   size_t nybins, double ylower, double yupper,
                                   const string& title,
                                   const string& xtitle,
                                   const string& ytitle,
                                   const string& ztitle)
  {
    const string path = histoPath(hname);
    Profile2DPtr prof = make_shared<Profile2D>(nxbins, xlower, xupper, nybins, ylower, yupper, path, title);
    prof->setAnnotation("XLabel", xtitle);
    prof->setAnnotation("YLabel", ytitle);
    prof->setAnnotation("ZLabel", ztitle);
    return addOrGetCompatAO(prof);
  }


  Profile2DPtr Analysis::bookProfile2D(const string& hname,
                                   const vector<double>& xbinedges,
                                   const vector<double>& ybinedges,
                                   const string& title,
                                   const string& xtitle,
                                   const string& ytitle,
                                   const string& ztitle)
  {
    const string path = histoPath(hname);
    Profile2DPtr prof = make_shared<Profile2D>(xbinedges, ybinedges, path, title);
    prof->setAnnotation("XLabel", xtitle);
    prof->setAnnotation("YLabel", ytitle);
    prof->setAnnotation("ZLabel", ztitle);
    return addOrGetCompatAO(prof);
  }


  Profile2DPtr Analysis::bookProfile2D(const string& hname,
                                       const initializer_list<double>& xbinedges,
                                       const initializer_list<double>& ybinedges,
                                       const string& title,
                                       const string& xtitle,
                                       const string& ytitle,
                                       const string& ztitle)
  {
    return bookProfile2D(hname, vector<double>{xbinedges}, vector<double>{ybinedges},
                         title, xtitle, ytitle, ztitle);
  }


  Profile2DPtr Analysis::bookProfile2D(const string& hname,
                                       const Scatter3D& refscatter,
                                       const string& title,
                                       const string& xtitle,
                                       const string& ytitle,
                                       const string& ztitle) {
    const string path = histoPath(hname);
    Profile2DPtr prof( new Profile2D(refscatter, path) );
    if (prof->hasAnnotation("IsRef")) prof->rmAnnotation("IsRef");
    if (prof->hasAnnotation("ErrorBreakdown")) prof->rmAnnotation("ErrorBreakdown");
    if (prof->hasAnnotation("Variations")) prof->rmAnnotation("Variations");
    prof->setTitle(title);
    prof->setAnnotation("XLabel", xtitle);
    prof->setAnnotation("YLabel", ytitle);
    prof->setAnnotation("ZLabel", ztitle);
    return addOrGetCompatAO(prof);
  }


  Profile2DPtr Analysis::bookProfile2D(const string& hname,
                                       const string& title,
                                       const string& xtitle,
                                       const string& ytitle,
                                       const string& ztitle) {
    const Scatter3D& refdata = refData<Scatter3D>(hname);
    return bookProfile2D(hname, refdata, title, xtitle, ytitle, ztitle);
  }


  Profile2DPtr Analysis::bookProfile2D(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId,
                                       const string& title,
                                       const string& xtitle,
                                       const string& ytitle,
                                       const string& ztitle) {
    const string axisCode = mkAxisCode(datasetId, xAxisId, yAxisId);
    return bookProfile2D(axisCode, title, xtitle, ytitle, ztitle);
  }


  /////////////////


  Scatter2DPtr Analysis::bookScatter2D(unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId,
                                       bool copy_pts,
                                       const string& title,
                                       const string& xtitle,
                                       const string& ytitle) {
    const string axisCode = mkAxisCode(datasetId, xAxisId, yAxisId);
    return bookScatter2D(axisCode, copy_pts, title, xtitle, ytitle);
  }


  Scatter2DPtr Analysis::bookScatter2D(const string& hname,
                                       bool copy_pts,
                                       const string& title,
                                       const string& xtitle,
                                       const string& ytitle) {
    Scatter2DPtr s;
    const string path = histoPath(hname);
    if (copy_pts) {
      const Scatter2D& refdata = refData(hname);
      s = make_shared<Scatter2D>(refdata, path);
      for (Point2D& p : s->points()) p.setY(0, 0);
    } else {
      s = make_shared<Scatter2D>(path);
    }
    if (s->hasAnnotation("IsRef")) s->rmAnnotation("IsRef");
    if (s->hasAnnotation("ErrorBreakdown")) s->rmAnnotation("ErrorBreakdown");
    if (s->hasAnnotation("Variations")) s->rmAnnotation("Variations");
    s->setTitle(title);
    s->setAnnotation("XLabel", xtitle);
    s->setAnnotation("YLabel", ytitle);
    return addOrGetCompatAO(s);
  }


  Scatter2DPtr Analysis::bookScatter2D(const string& hname,
                                       size_t npts, double lower, double upper,
                                       const string& title,
                                       const string& xtitle,
                                       const string& ytitle) {
    const string path = histoPath(hname);
    Scatter2DPtr s = make_shared<Scatter2D>(path);
    const double binwidth = (upper-lower)/npts;
    for (size_t pt = 0; pt < npts; ++pt) {
      const double bincentre = lower + (pt + 0.5) * binwidth;
      s->addPoint(bincentre, 0, binwidth/2.0, 0);
    }
    s->setTitle(title);
    s->setAnnotation("XLabel", xtitle);
    s->setAnnotation("YLabel", ytitle);
    return addOrGetCompatAO(s);
  }

  Scatter2DPtr Analysis::bookScatter2D(const string& hname,
                                       const vector<double>& binedges,
                                       const string& title,
                                       const string& xtitle,
                                       const string& ytitle) {
    const string path = histoPath(hname);
    Scatter2DPtr s = make_shared<Scatter2D>(path);
    for (size_t pt = 0; pt < binedges.size()-1; ++pt) {
      const double bincentre = (binedges[pt] + binedges[pt+1]) / 2.0;
      const double binwidth = binedges[pt+1] - binedges[pt];
      s->addPoint(bincentre, 0, binwidth/2.0, 0);
    }
    s->setTitle(title);
    s->setAnnotation("XLabel", xtitle);
    s->setAnnotation("YLabel", ytitle);
    return addOrGetCompatAO(s);
  }

  Scatter2DPtr Analysis::bookScatter2D(const Scatter2DPtr scPtr, 
    const std::string& path, const std::string& title, 
    const std::string& xtitle, const std::string& ytitle ) {
    
    Scatter2DPtr s = make_shared<Scatter2D>(*scPtr);
    s->setPath(path);
    s->setTitle(title);
    s->setAnnotation("XLabel",xtitle);
    s->setAnnotation("YLabel",ytitle);
    return addOrGetCompatAO(s);
  
  }

  /////////////////////


  void Analysis::divide(CounterPtr c1, CounterPtr c2, Scatter1DPtr s) const {
    const string path = s->path();
    *s = *c1 / *c2;
    s->setPath(path);
  }

  void Analysis::divide(const Counter& c1, const Counter& c2, Scatter1DPtr s) const {
    const string path = s->path();
    *s = c1 / c2;
    s->setPath(path);
  }


  void Analysis::divide(Histo1DPtr h1, Histo1DPtr h2, Scatter2DPtr s) const {
    const string path = s->path();
    *s = *h1 / *h2;
    s->setPath(path);
  }

  void Analysis::divide(const Histo1D& h1, const Histo1D& h2, Scatter2DPtr s) const {
    const string path = s->path();
    *s = h1 / h2;
    s->setPath(path);
  }


  void Analysis::divide(Profile1DPtr p1, Profile1DPtr p2, Scatter2DPtr s) const {
    const string path = s->path();
    *s = *p1 / *p2;
    s->setPath(path);
  }

  void Analysis::divide(const Profile1D& p1, const Profile1D& p2, Scatter2DPtr s) const {
    const string path = s->path();
    *s = p1 / p2;
    s->setPath(path);
  }


  void Analysis::divide(Histo2DPtr h1, Histo2DPtr h2, Scatter3DPtr s) const {
    const string path = s->path();
    *s = *h1 / *h2;
    s->setPath(path);
  }

  void Analysis::divide(const Histo2D& h1, const Histo2D& h2, Scatter3DPtr s) const {
    const string path = s->path();
    *s = h1 / h2;
    s->setPath(path);
  }


  void Analysis::divide(Profile2DPtr p1, Profile2DPtr p2, Scatter3DPtr s) const {
    const string path = s->path();
    *s = *p1 / *p2;
    s->setPath(path);
  }

  void Analysis::divide(const Profile2D& p1, const Profile2D& p2, Scatter3DPtr s) const {
    const string path = s->path();
    *s = p1 / p2;
    s->setPath(path);
  }


  /// @todo Counter and Histo2D efficiencies and asymms


  void Analysis::efficiency(Histo1DPtr h1, Histo1DPtr h2, Scatter2DPtr s) const {
    const string path = s->path();
    *s = YODA::efficiency(*h1, *h2);
    s->setPath(path);
  }

  void Analysis::efficiency(const Histo1D& h1, const Histo1D& h2, Scatter2DPtr s) const {
    const string path = s->path();
    *s = YODA::efficiency(h1, h2);
    s->setPath(path);
  }


  void Analysis::asymm(Histo1DPtr h1, Histo1DPtr h2, Scatter2DPtr s) const {
    const string path = s->path();
    *s = YODA::asymm(*h1, *h2);
    s->setPath(path);
  }

  void Analysis::asymm(const Histo1D& h1, const Histo1D& h2, Scatter2DPtr s) const {
    const string path = s->path();
    *s = YODA::asymm(h1, h2);
    s->setPath(path);
  }


  void Analysis::scale(CounterPtr cnt, double factor) {
    if (!cnt) {
      MSG_WARNING("Failed to scale counter=NULL in analysis " << name() << " (scale=" << factor << ")");
      return;
    }
    if (std::isnan(factor) || std::isinf(factor)) {
      MSG_WARNING("Failed to scale counter=" << cnt->path() << " in analysis: " << name() << " (invalid scale factor = " << factor << ")");
      factor = 0;
    }
    MSG_TRACE("Scaling counter " << cnt->path() << " by factor " << factor);
    try {
      cnt->scaleW(factor);
    } catch (YODA::Exception& we) {
      MSG_WARNING("Could not scale counter " << cnt->path());
      return;
    }
  }


  void Analysis::normalize(Histo1DPtr histo, double norm, bool includeoverflows) {
    if (!histo) {
      MSG_WARNING("Failed to normalize histo=NULL in analysis " << name() << " (norm=" << norm << ")");
      return;
    }
    MSG_TRACE("Normalizing histo " << histo->path() << " to " << norm);
    try {
      const double hint = histo->integral(includeoverflows);
      if (hint == 0)  MSG_WARNING("Skipping histo with null area " << histo->path());
      else            histo->normalize(norm, includeoverflows);
    } catch (YODA::Exception& we) {
      MSG_WARNING("Could not normalize histo " << histo->path());
      return;
    }
  }


  void Analysis::scale(Histo1DPtr histo, double factor) {
    if (!histo) {
      MSG_WARNING("Failed to scale histo=NULL in analysis " << name() << " (scale=" << factor << ")");
      return;
    }
    if (std::isnan(factor) || std::isinf(factor)) {
      MSG_WARNING("Failed to scale histo=" << histo->path() << " in analysis: " << name() << " (invalid scale factor = " << factor << ")");
      factor = 0;
    }
    MSG_TRACE("Scaling histo " << histo->path() << " by factor " << factor);
    try {
      histo->scaleW(factor);
    } catch (YODA::Exception& we) {
      MSG_WARNING("Could not scale histo " << histo->path());
      return;
    }
  }


  void Analysis::normalize(Histo2DPtr histo, double norm, bool includeoverflows) {
    if (!histo) {
      MSG_ERROR("Failed to normalize histo=NULL in analysis " << name() << " (norm=" << norm << ")");
      return;
    }
    MSG_TRACE("Normalizing histo " << histo->path() << " to " << norm);
    try {
      const double hint = histo->integral(includeoverflows);
      if (hint == 0)  MSG_WARNING("Skipping histo with null area " << histo->path());
      else            histo->normalize(norm, includeoverflows);
    } catch (YODA::Exception& we) {
      MSG_WARNING("Could not normalize histo " << histo->path());
      return;
    }
  }


  void Analysis::scale(Histo2DPtr histo, double factor) {
    if (!histo) {
      MSG_ERROR("Failed to scale histo=NULL in analysis " << name() << " (scale=" << factor << ")");
      return;
    }
    if (std::isnan(factor) || std::isinf(factor)) {
      MSG_ERROR("Failed to scale histo=" << histo->path() << " in analysis: " << name() << " (invalid scale factor = " << factor << ")");
      factor = 0;
    }
    MSG_TRACE("Scaling histo " << histo->path() << " by factor " << factor);
    try {
      histo->scaleW(factor);
    } catch (YODA::Exception& we) {
      MSG_WARNING("Could not scale histo " << histo->path());
      return;
    }
  }


  void Analysis::integrate(Histo1DPtr h, Scatter2DPtr s) const {
    // preserve the path info
    const string path = s->path();
    *s = toIntegralHisto(*h);
    s->setPath(path);
  }

  void Analysis::integrate(const Histo1D& h, Scatter2DPtr s) const {
    // preserve the path info
    const string path = s->path();
    *s = toIntegralHisto(h);
    s->setPath(path);
  }


  /// @todo 2D versions of integrate... defined how, exactly?!?


  //////////////////////////////////


  void Analysis::addAnalysisObject(AnalysisObjectPtr ao) {
    _analysisobjects.push_back(ao);
  }


  void Analysis::removeAnalysisObject(const string& path) {
    for (vector<AnalysisObjectPtr>::iterator it = _analysisobjects.begin();  it != _analysisobjects.end(); ++it) {
      if ((*it)->path() == path) {
        _analysisobjects.erase(it);
        break;
      }
    }
  }


  void Analysis::removeAnalysisObject(AnalysisObjectPtr ao) {
    for (vector<AnalysisObjectPtr>::iterator it = _analysisobjects.begin();  it != _analysisobjects.end(); ++it) {
      if (*it == ao) {
        _analysisobjects.erase(it);
        break;
      }
    }
 }

const CentralityProjection &
Analysis::declareCentrality(const SingleValueProjection &proj,
                            string calAnaName, string calHistName,
                            const string projName, bool increasing) {

  CentralityProjection cproj;
  
  // Select the centrality variable from option. Use REF as default.
  // Other selections are "GEN", "IMP" and "USR" (USR only in HEPMC 3).
  string sel = getOption<string>("cent","REF");
  set<string> done;

  if ( sel == "REF" ) {
    Scatter2DPtr refscat;
    auto refmap = getRefData(calAnaName);
    if ( refmap.find(calHistName) != refmap.end() )
      refscat =
        dynamic_pointer_cast<Scatter2D>(refmap.find(calHistName)->second);

    if ( !refscat ) {
      MSG_WARNING("No reference calibration histogram for " <<
                  "CentralityProjection " << projName << " found " <<
                  "(requested histogram " << calHistName << " in " <<
                  calAnaName << ")");
    }
    else {
      MSG_INFO("Found calibration histogram " << sel << " " << refscat->path());
      cproj.add(PercentileProjection(proj, refscat, increasing), sel);
    }
  }
  else if ( sel == "GEN" ) {
    Histo1DPtr genhist;
    string histpath = "/" + calAnaName + "/" + calHistName;
    for ( AnalysisObjectPtr ao : handler().getData(true, false, false) ) {
      if ( ao->path() == histpath )
        genhist = dynamic_pointer_cast<Histo1D>(ao);
    }
    if ( !genhist || genhist->numEntries() <= 1 ) {
      MSG_WARNING("No generated calibration histogram for " <<
               "CentralityProjection " << projName << " found " <<
               "(requested histogram " << calHistName << " in " <<
               calAnaName << ")");
    }
    else {
      MSG_INFO("Found calibration histogram " << sel << " " << genhist->path());
      cproj.add(PercentileProjection(proj, genhist, increasing), sel);
    }
  }
  else if ( sel == "IMP" ) {
    Histo1DPtr imphist =
      getAnalysisObject<Histo1D>(calAnaName, calHistName + "_IMP");
    if ( !imphist || imphist->numEntries() <= 1 ) {
      MSG_WARNING("No impact parameter calibration histogram for " <<
               "CentralityProjection " << projName << " found " <<
               "(requested histogram " << calHistName << "_IMP in " <<
               calAnaName << ")");
    }
    else {
      MSG_INFO("Found calibration histogram " << sel << " " << imphist->path());
      cproj.add(PercentileProjection(ImpactParameterProjection(),
                                     imphist, true), sel);
    }
  }
  else if ( sel == "USR" ) {
#if HEPMC_VERSION_CODE >= 3000000
    Histo1DPtr usrhist =
      getAnalysisObject<Histo1D>(calAnaName, calHistName + "_USR");
    if ( !usrhist || usrhist->numEntries() <= 1 ) {
      MSG_WARNING("No user-defined calibration histogram for " <<
               "CentralityProjection " << projName << " found " <<
               "(requested histogram " << calHistName << "_USR in " <<
               calAnaName << ")");
      continue;
    }
    else {
      MSG_INFO("Found calibration histogram " << sel << " " << usrhist->path());
      cproj.add((UserCentEstimate(), usrhist, true), sel);
     }
#else
      MSG_WARNING("UserCentEstimate is only available with HepMC3.");
#endif
    }
  else if ( sel == "RAW" ) {
#if HEPMC_VERSION_CODE >= 3000000
    cproj.add(GeneratedCentrality(), sel);
#else
    MSG_WARNING("GeneratedCentrality is only available with HepMC3.");
#endif
  }
    else
      MSG_WARNING("'" << sel << "' is not a valid PercentileProjection tag.");

  if ( cproj.empty() )
    MSG_WARNING("CentralityProjection " << projName
                << " did not contain any valid PercentileProjections.");

  return declare(cproj, projName);
  
}

}
