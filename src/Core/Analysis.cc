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
#include <regex>

namespace Rivet {


  Analysis::Analysis(const string& name)
    : _analysishandler(nullptr)
  {
    ProjectionApplier::_allowProjReg = false;
    _defaultname = name;

    unique_ptr<AnalysisInfo> ai = AnalysisInfo::make(name);
    assert(ai);
    _info = move(ai);
    assert(_info);
  }


  double Analysis::sqrtS() const {
    double sqrts = handler().sqrtS();
    if (sqrts <= 0) {
      MSG_DEBUG("Suspicious beam energy. You're probably running rivet-merge. Fetching beam energy from option.");
      sqrts = getOption<double>("ENERGY", 0);
    }
    return sqrts;
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
    std::stringstream axisCode;
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


  bool Analysis::isCompatibleWithSqrtS(const float energy, float tolerance) const {
    return fuzzyEquals(sqrtS()/GeV, energy, tolerance);
  }

  ///////////////////////////////////////////

  double Analysis::crossSection() const {
    const YODA::Scatter1D::Points& ps = handler().crossSection()->points();
    if (ps.size() != 1) {
      string errMsg = "cross section missing for analysis " + name();
      throw Error(errMsg);
    }
    return ps[0].x();
  }

  double Analysis::crossSectionError() const {
    const YODA::Scatter1D::Points& ps = handler().crossSection()->points();
    if (ps.size() != 1) {
      string errMsg = "cross section missing for analysis " + name();
      throw Error(errMsg);
    }
    return ps[0].xErrAvg();
  }

  double Analysis::crossSectionPerEvent() const {
    return crossSection()/sumW();
  }

  double Analysis::crossSectionErrorPerEvent() const {
    return crossSectionError()/sumW();
  }

  ////////////////////////////////////////////////////////////
  // Histogramming


  void Analysis::_cacheRefData() const {
    if (_refdata.empty()) {
      MSG_TRACE("Getting refdata cache for paper " << name());
      _refdata = getRefData(getRefDataName());
    }
  }

  // vector<YODA::AnalysisObjectPtr> Analysis::getAllData(bool includeorphans) const{
  //   return handler().getData(includeorphans, false, false);
  // }


  CounterPtr& Analysis::book(CounterPtr& ctr, const string& cname) {
    // const string path = histoPath(cname);
    // ctr = CounterPtr(handler().weightNames(), Counter(path, title));
    // ctr = addAnalysisObject(ctr);
    // return ctr;
    return ctr = registerAO( Counter(histoPath(cname)) );
  }


  CounterPtr& Analysis::book(CounterPtr& ctr, unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) {
    const string axisCode = mkAxisCode(datasetId, xAxisId, yAxisId);
    return book(ctr, axisCode);
  }




  Histo1DPtr& Analysis::book(Histo1DPtr& histo, const string& hname, size_t nbins, double lower, double upper) {
    const string path = histoPath(hname);

    Histo1D hist = Histo1D(nbins, lower, upper, path);
    _setWriterPrecision(path, hist);

    // histo = Histo1DPtr(handler().weightNames(), hist);
    // histo = addAnalysisObject(histo);
    // return histo;
    return histo = registerAO(hist);
  }

  Histo1DPtr& Analysis::book(Histo1DPtr& histo, const string& hname, const initializer_list<double>& binedges) {
    return book(histo, hname, vector<double>{binedges});
  }

  Histo1DPtr& Analysis::book(Histo1DPtr& histo, const string& hname, const vector<double>& binedges) {
    const string path = histoPath(hname);

    Histo1D hist = Histo1D(binedges, path);
    _setWriterPrecision(path, hist);

    // histo = Histo1DPtr(handler().weightNames(), hist);
    // histo = addAnalysisObject(histo);
    // return histo;
    return histo = registerAO(hist);
  }

  Histo1DPtr& Analysis::book(Histo1DPtr& histo, const string& hname) {
    const Scatter2D& refdata = refData(hname);
    return book(histo, hname, refdata);
  }


  Histo1DPtr& Analysis::book(Histo1DPtr& histo, unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) {
    const string axisCode = mkAxisCode(datasetId, xAxisId, yAxisId);
    return book(histo, axisCode);
  }

  Histo1DPtr& Analysis::book(Histo1DPtr& histo, const string& hname, const Scatter2D& refscatter) {
    const string path = histoPath(hname);

    Histo1D hist = Histo1D(refscatter, path);
    for (const string& a : hist.annotations()) {
      if (a != "Path")  hist.rmAnnotation(a);
    }
    _setWriterPrecision(path, hist);

    // histo = Histo1DPtr(handler().weightNames(), hist);
    // histo = addAnalysisObject(histo);
    // return histo;
    return histo = registerAO(hist);
  }


  /////////////////


  Histo2DPtr& Analysis::book(Histo2DPtr& h2d,const string& hname,
                             size_t nxbins, double xlower, double xupper,
                             size_t nybins, double ylower, double yupper) {
    const string path = histoPath(hname);

    Histo2D hist(nxbins, xlower, xupper, nybins, ylower, yupper, path);
    _setWriterPrecision(path, hist);

    // h2d = Histo2DPtr(handler().weightNames(), hist);
    // h2d = addAnalysisObject(h2d);
    // return h2d;
    return h2d = registerAO(hist);
  }

  Histo2DPtr& Analysis::book(Histo2DPtr& h2d,const string& hname,
                             const initializer_list<double>& xbinedges,
                             const initializer_list<double>& ybinedges) {
    return book(h2d, hname, vector<double>{xbinedges}, vector<double>{ybinedges});
  }

  Histo2DPtr& Analysis::book(Histo2DPtr& h2d,const string& hname,
                             const vector<double>& xbinedges,
                             const vector<double>& ybinedges) {
    const string path = histoPath(hname);

    Histo2D hist(xbinedges, ybinedges, path);
    _setWriterPrecision(path, hist);

    // h2d = Histo2DPtr(handler().weightNames(), hist);
    // h2d = addAnalysisObject(h2d);
    // return h2d;
    return h2d = registerAO(hist);
  }


  Histo2DPtr& Analysis::book(Histo2DPtr& histo, const string& hname, const Scatter3D& refscatter) {
    const string path = histoPath(hname);

    Histo2D hist = Histo2D(refscatter, path);
    for (const string& a : hist.annotations()) {
      if (a != "Path")  hist.rmAnnotation(a);
    }
    _setWriterPrecision(path, hist);

    // histo = Histo2DPtr(handler().weightNames(), hist);
    // histo = addAnalysisObject(histo);
    // return histo;
    return histo = registerAO(hist);
  }


  Histo2DPtr& Analysis::book(Histo2DPtr& histo, const string& hname) {
    const Scatter3D& refdata = refData<Scatter3D>(hname);
    return book(histo, hname, refdata);
  }


  Histo2DPtr& Analysis::book(Histo2DPtr& histo, unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) {
    const string axisCode = mkAxisCode(datasetId, xAxisId, yAxisId);
    return book(histo, axisCode);
  }


  /////////////////


  Profile1DPtr& Analysis::book(Profile1DPtr& p1d,const string& hname, size_t nbins, double lower, double upper) {
    const string path = histoPath(hname);

    Profile1D prof(nbins, lower, upper, path);
    _setWriterPrecision(path, prof);

    // p1d = Profile1DPtr(handler().weightNames(), prof);
    // p1d = addAnalysisObject(p1d);
    // return p1d;
    return p1d = registerAO(prof);
  }


  Profile1DPtr& Analysis::book(Profile1DPtr& p1d,const string& hname, const initializer_list<double>& binedges) {
    return book(p1d, hname, vector<double>{binedges});
  }

  Profile1DPtr& Analysis::book(Profile1DPtr& p1d, const string& hname, const vector<double>& binedges) {
    const string path = histoPath(hname);

    Profile1D prof(binedges, path);
    _setWriterPrecision(path, prof);

    // p1d = Profile1DPtr(handler().weightNames(), prof);
    // p1d = addAnalysisObject(p1d);
    // return p1d;
    return p1d = registerAO(prof);
  }

  Profile1DPtr& Analysis::book(Profile1DPtr& p1d, const string& hname, const Scatter2D& refscatter) {
    const string path = histoPath(hname);

    Profile1D prof(refscatter, path);
    for (const string& a : prof.annotations()) {
      if (a != "Path")  prof.rmAnnotation(a);
    }
    _setWriterPrecision(path, prof);

    // p1d = Profile1DPtr(handler().weightNames(), prof);
    // p1d = addAnalysisObject(p1d);
    // return p1d;
    return p1d = registerAO(prof);
  }


  Profile1DPtr& Analysis::book(Profile1DPtr& p1d,const string& hname) {
    const Scatter2D& refdata = refData(hname);
    return  book(p1d, hname, refdata);
  }


  Profile1DPtr& Analysis::book(Profile1DPtr& p1d,unsigned int datasetId, unsigned int xAxisId, unsigned int yAxisId) {
    const string axisCode = mkAxisCode(datasetId, xAxisId, yAxisId);
    return book(p1d, axisCode);
  }


  ///////////////////


  Profile2DPtr& Analysis::book(Profile2DPtr& p2d, const string& hname,
                               size_t nxbins, double xlower, double xupper,
                               size_t nybins, double ylower, double yupper) {
    const string path = histoPath(hname);

    Profile2D prof(nxbins, xlower, xupper, nybins, ylower, yupper, path);
    _setWriterPrecision(path, prof);

    // p2d = Profile2DPtr(handler().weightNames(), prof);
    // p2d = addAnalysisObject(p2d);
    // return p2d;
    return p2d = registerAO(prof);
  }


  Profile2DPtr& Analysis::book(Profile2DPtr& p2d, const string& hname,
                               const initializer_list<double>& xbinedges,
                               const initializer_list<double>& ybinedges) {
    return book(p2d, hname, vector<double>{xbinedges}, vector<double>{ybinedges});
  }


  Profile2DPtr& Analysis::book(Profile2DPtr& p2d, const string& hname,
                               const vector<double>& xbinedges,
                               const vector<double>& ybinedges) {
    const string path = histoPath(hname);

    Profile2D prof(xbinedges, ybinedges, path);
    _setWriterPrecision(path, prof);

    // p2d = Profile2DPtr(handler().weightNames(), prof);
    // p2d = addAnalysisObject(p2d);
    // return p2d;
    return p2d = registerAO(prof);
  }


  /// @todo REINSTATE

  // Profile2DPtr Analysis::book(Profile2DPtr& prof,const string& hname,
  //                                      const Scatter3D& refscatter) {
  //   const string path = histoPath(hname);

  //   /// @todo Add no-metadata argument to YODA copy constructors
  //   Profile2D prof(refscatter, path);
  //   if (prof.hasAnnotation("IsRef")) prof.rmAnnotation("IsRef");

  //   p2d = Profile2DPtr(handler().weightNames(), prof);
  //   p2d = addAnalysisObject(p2d);
  //   return p2d;
  // }


  // Profile2DPtr Analysis::book(Profile2DPtr& prof, const string& hname) {
  //   const Scatter3D& refdata = refData<Scatter3D>(hname);
  //   return book(prof, hname, refdata);
  // }


  /// @todo Should be able to book Scatter1Ds


  ///////////////


  Scatter2DPtr& Analysis::book(Scatter2DPtr& s2d, unsigned int datasetId,
                               unsigned int xAxisId, unsigned int yAxisId, bool copy_pts) {
    const string axisCode = mkAxisCode(datasetId, xAxisId, yAxisId);
    return book(s2d, axisCode, copy_pts);
  }


  Scatter2DPtr& Analysis::book(Scatter2DPtr& s2d, const string& hname, bool copy_pts) {
    const string path = histoPath(hname);

    Scatter2D scat;
    if (copy_pts) {
      const Scatter2D& refdata = refData(hname);
      scat = Scatter2D(refdata, path);
      for (Point2D& p : scat.points()) p.setY(0, 0);
      for (const string& a : scat.annotations()) {
        if (a != "Path")  scat.rmAnnotation(a);
      }
    } else {
      scat = Scatter2D(path);
    }
    _setWriterPrecision(path, scat);

    // s2d = Scatter2DPtr(handler().weightNames(), scat);
    // s2d = addAnalysisObject(s2d);
    // return s2d;
    return s2d = registerAO(scat);
  }


  Scatter2DPtr& Analysis::book(Scatter2DPtr& s2d, const string& hname, size_t npts, double lower, double upper) {
    const string path = histoPath(hname);

    Scatter2D scat(path);
    const double binwidth = (upper-lower)/npts;
    for (size_t pt = 0; pt < npts; ++pt) {
      const double bincentre = lower + (pt + 0.5) * binwidth;
      scat.addPoint(bincentre, 0, binwidth/2.0, 0);
    }
    _setWriterPrecision(path, scat);

    // s2d = Scatter2DPtr(handler().weightNames(), scat);
    // s2d = addAnalysisObject(s2d);
    // return s2d;
    return s2d = registerAO(scat);
  }

  Scatter2DPtr& Analysis::book(Scatter2DPtr& s2d, const string& hname, const vector<double>& binedges) {
    const string path = histoPath(hname);

    Scatter2D scat(path);
    for (size_t pt = 0; pt < binedges.size()-1; ++pt) {
      const double bincentre = (binedges[pt] + binedges[pt+1]) / 2.0;
      const double binwidth = binedges[pt+1] - binedges[pt];
      scat.addPoint(bincentre, 0, binwidth/2.0, 0);
    }
    _setWriterPrecision(path, scat);

    // s2d = Scatter2DPtr(handler().weightNames(), scat);
    // s2d = addAnalysisObject(s2d);
    // return s2d;
    return s2d = registerAO(scat);
  }

  Scatter2DPtr& Analysis::book(Scatter2DPtr& s2d, const string& hname, const Scatter2D& refscatter) {
    const string path = histoPath(hname);

    Scatter2D scat(refscatter, path);
    for (const string& a : scat.annotations()) {
      if (a != "Path")  scat.rmAnnotation(a);
    }
    _setWriterPrecision(path, scat);

    return s2d = registerAO(scat);
  }


  ///////////////


  Scatter3DPtr& Analysis::book(Scatter3DPtr& s3d, unsigned int datasetId, unsigned int xAxisId,
                               unsigned int yAxisId, unsigned int zAxisId, bool copy_pts) {
    const string axisCode = mkAxisCode(datasetId, xAxisId, yAxisId);
    return book(s3d, axisCode, copy_pts);
  }


  Scatter3DPtr& Analysis::book(Scatter3DPtr& s3d, const string& hname, bool copy_pts) {
    const string path = histoPath(hname);

    Scatter3D scat;
    if (copy_pts) {
      const Scatter3D& refdata = refData<Scatter3D>(hname);
      scat = Scatter3D(refdata, path);
      for (Point3D& p : scat.points()) p.setZ(0, 0);
      for (const string& a : scat.annotations()) {
        if (a != "Path")  scat.rmAnnotation(a);
      }
    } else {
      scat = Scatter3D(path);
    }
    _setWriterPrecision(path, scat);

    // s3d = Scatter3DPtr(handler().weightNames(), scat);
    // s3d = addAnalysisObject(s3d);
    // return s3d;
    return s3d = registerAO(scat);
  }


  Scatter3DPtr& Analysis::book(Scatter3DPtr& s3d, const string& hname,
                               size_t xnpts, double xlower, double xupper,
                               size_t ynpts, double ylower, double yupper) {
    const string path = histoPath(hname);

    Scatter3D scat(path);
    const double xbinwidth = (xupper-xlower)/xnpts;
    const double ybinwidth = (yupper-ylower)/ynpts;
    for (size_t xpt = 0; xpt < xnpts; ++xpt) {
      const double xbincentre = xlower + (xpt + 0.5) * xbinwidth;
      for (size_t ypt = 0; ypt < ynpts; ++ypt) {
        const double ybincentre = ylower + (ypt + 0.5) * ybinwidth;
        scat.addPoint(xbincentre, ybincentre, 0, 0.5*xbinwidth, 0.5*ybinwidth, 0);
      }
    }
    _setWriterPrecision(path, scat);

    // s3d = Scatter3DPtr(handler().weightNames(), scat);
    // s3d = addAnalysisObject(s3d);
    // return s3d;
    return s3d = registerAO(scat);
  }

  Scatter3DPtr& Analysis::book(Scatter3DPtr& s3d, const string& hname,
                               const vector<double>& xbinedges,
                               const vector<double>& ybinedges) {
    const string path = histoPath(hname);

    Scatter3D scat(path);
    for (size_t xpt = 0; xpt < xbinedges.size()-1; ++xpt) {
      const double xbincentre = (xbinedges[xpt] + xbinedges[xpt+1]) / 2.0;
      const double xbinwidth = xbinedges[xpt+1] - xbinedges[xpt];
      for (size_t ypt = 0; ypt < ybinedges.size()-1; ++ypt) {
        const double ybincentre = (ybinedges[ypt] + ybinedges[ypt+1]) / 2.0;
        const double ybinwidth = ybinedges[ypt+1] - ybinedges[ypt];
        scat.addPoint(xbincentre, ybincentre, 0, 0.5*xbinwidth, 0.5*ybinwidth, 0);
      }
    }
    _setWriterPrecision(path, scat);

    // s3d = Scatter3DPtr(handler().weightNames(), scat);
    // s3d = addAnalysisObject(s3d);
    // return s3d;
    return s3d = registerAO(scat);
  }

  Scatter3DPtr& Analysis::book(Scatter3DPtr& s3d, const string& hname, const Scatter3D& refscatter) {
    const string path = histoPath(hname);

    Scatter3D scat(refscatter, path);
    for (const string& a : scat.annotations()) {
      if (a != "Path")  scat.rmAnnotation(a);
    }
    _setWriterPrecision(path, scat);

    return s3d = registerAO(scat);
  }


  ///////////////



  /////////////////////

  void Analysis::barchart(Histo1DPtr h, Scatter2DPtr s, bool usefocus) const {
    const string path = s->path();
    *s = mkScatter(*h, usefocus, false); //< do NOT divide by bin width cf. a differential dsigma/dX histogram
    s->setPath(path);
  }


  void Analysis::barchart(Histo2DPtr h, Scatter3DPtr s, bool usefocus) const {
    const string path = s->path();
    *s = mkScatter(*h, usefocus, false); //< do NOT divide by bin area cf. a differential d^2sigma/dXdY histogram
    s->setPath(path);
  }


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


  void Analysis::scale(CounterPtr cnt, Analysis::CounterAdapter factor) {
    if (!cnt) {
      MSG_WARNING("Failed to scale counter=NULL in analysis " << name() << " (scale=" << double(factor) << ")");
      return;
    }
    if (std::isnan(double(factor)) || std::isinf(double(factor))) {
      MSG_WARNING("Failed to scale counter=" << cnt->path() << " in analysis: " << name() << " (invalid scale factor = " << double(factor) << ")");
      factor = 0;
    }
    MSG_TRACE("Scaling counter " << cnt->path() << " by factor " << double(factor));
    try {
      cnt->scaleW(factor);
    } catch (YODA::Exception& we) {
      MSG_WARNING("Could not scale counter " << cnt->path());
      return;
    }
  }


  void Analysis::normalize(Histo1DPtr histo, Analysis::CounterAdapter norm, bool includeoverflows) {
    if (!histo) {
      MSG_WARNING("Failed to normalize histo=NULL in analysis " << name() << " (norm=" << double(norm) << ")");
      return;
    }
    MSG_TRACE("Normalizing histo " << histo->path() << " to " << double(norm));
    try {
      const double hint = histo->integral(includeoverflows);
      if (hint == 0)  MSG_DEBUG("Skipping histo with null area " << histo->path());
      else            histo->normalize(norm, includeoverflows);
    } catch (YODA::Exception& we) {
      MSG_WARNING("Could not normalize histo " << histo->path());
      return;
    }
  }


  void Analysis::scale(Histo1DPtr histo, Analysis::CounterAdapter factor) {
    if (!histo) {
      MSG_WARNING("Failed to scale histo=NULL in analysis " << name() << " (scale=" << double(factor) << ")");
      return;
    }
    if (std::isnan(double(factor)) || std::isinf(double(factor))) {
      MSG_WARNING("Failed to scale histo=" << histo->path() << " in analysis: " << name() << " (invalid scale factor = " << double(factor) << ")");
      factor = 0;
    }
    MSG_TRACE("Scaling histo " << histo->path() << " by factor " << double(factor));
    try {
      histo->scaleW(factor);
    } catch (YODA::Exception& we) {
      MSG_WARNING("Could not scale histo " << histo->path());
      return;
    }
  }


  void Analysis::normalize(Histo2DPtr histo, Analysis::CounterAdapter norm, bool includeoverflows) {
    if (!histo) {
      MSG_ERROR("Failed to normalize histo=NULL in analysis " << name() << " (norm=" << double(norm) << ")");
      return;
    }
    MSG_TRACE("Normalizing histo " << histo->path() << " to " << double(norm));
    try {
      const double hint = histo->integral(includeoverflows);
      if (hint == 0)  MSG_DEBUG("Skipping histo with null area " << histo->path());
      else            histo->normalize(norm, includeoverflows);
    } catch (YODA::Exception& we) {
      MSG_WARNING("Could not normalize histo " << histo->path());
      return;
    }
  }


  void Analysis::scale(Histo2DPtr histo, Analysis::CounterAdapter factor) {
    if (!histo) {
      MSG_ERROR("Failed to scale histo=NULL in analysis " << name() << " (scale=" << double(factor) << ")");
      return;
    }
    if (std::isnan(double(factor)) || std::isinf(double(factor))) {
      MSG_ERROR("Failed to scale histo=" << histo->path() << " in analysis: " << name() << " (invalid scale factor = " << double(factor) << ")");
      factor = 0;
    }
    MSG_TRACE("Scaling histo " << histo->path() << " by factor " << double(factor));
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

}
/// @todo 2D versions of integrate... defined how, exactly?!?


//////////////////////////////////

// namespace {
//   void errormsg(std::string name) {
//     // #ifdef HAVE_BACKTRACE
//     //      void * buffer[4];
//     //      backtrace(buffer, 4);
//     //      backtrace_symbols_fd(buffer, 4 , 1);
//     // #endif
//     std::cerr << name << ": Can't book objects outside of init().\n";
//     assert(false);
//   }
// }


namespace Rivet {


  // void Analysis::addAnalysisObject(const MultiweightAOPtr& ao) {
  //   if (handler().stage() == AnalysisHandler::Stage::INIT) {
  //     _analysisobjects.push_back(ao);
  //   }
  //   else {
  //     errormsg(name());
  //   }
  // }


  void Analysis::removeAnalysisObject(const string& path) {
    for (auto it = _analysisobjects.begin();
         it != _analysisobjects.end(); ++it) {
      if ((*it)->path() == path) {
        _analysisobjects.erase(it);
        break;
      }
    }
  }

  void Analysis::removeAnalysisObject(const MultiweightAOPtr& ao) {
    for (auto it = _analysisobjects.begin();  it != _analysisobjects.end(); ++it) {
      if ((*it) == ao) {
        _analysisobjects.erase(it);
        break;
      }
    }
  }

  const CentralityProjection&
  Analysis::declareCentrality(const SingleValueProjection &proj,
                              string calAnaName, string calHistName,
                              const string projName, bool increasing) {

    CentralityProjection cproj;

    // Select the centrality variable from option. Use REF as default.
    // Other selections are "GEN", "IMP", "RAW" and "USR" (USR only in HEPMC 3).
    string sel = getOption<string>("cent","REF");
    set<string> done;

    if ( sel == "REF" ) {
      YODA::Scatter2DPtr refscat;
      map<string, YODA::AnalysisObjectPtr> refmap;
      try {
        refmap = getRefData(calAnaName);
      } catch (...) {
        MSG_ERROR("No reference calibration file for CentralityProjection "
                   << calAnaName << " found.\nDid you mean to generate one yourself?\n"
                   << "Once generated, you can preload the calibration file using the "
                   << "-p flag\nand steer the routine logic using the 'cent' option "
                   << "with the appropriate value (e.g. =GEN).");
        exit(1);

      }
      if ( refmap.find(calHistName) != refmap.end() )
        refscat = dynamic_pointer_cast<Scatter2D>(refmap.find(calHistName)->second);

      if ( !refscat ) {
        MSG_WARNING("No reference calibration histogram for " <<
                    "CentralityProjection " << projName << " found " <<
                    "(requested histogram " << calHistName << " in " <<
                    calAnaName << ")");
      }
      else {
        MSG_INFO("Found calibration histogram " << sel << " " << refscat->path());
        cproj.add(PercentileProjection(proj, *refscat, increasing), sel);
      }
    }
    else if ( sel == "GEN" ) {
      YODA::Histo1DPtr genhists =
        getPreload<Histo1D>("/" + calAnaName + "/" + calHistName);
      if ( !genhists || genhists->numEntries() <= 1 ) {
        MSG_WARNING("No generated calibration histogram for " <<
                    "CentralityProjection " << projName << " found " <<
                    "(requested histogram " << calHistName << " in " <<
                    calAnaName << ").\nDid you mean to generate one " <<
		    "yourself and forgot? Generate with " << calAnaName <<
		    " and preload the calibration file using the -p flag.");
      }
      else {
        MSG_INFO("Found calibration histogram " << sel << " " << genhists->path());
        cproj.add(PercentileProjection(proj, *genhists, increasing), sel);
      }
    }
    else if ( sel == "IMP" ) {
      YODA::Histo1DPtr imphists =
        getPreload<Histo1D>("/" + calAnaName + "/" + calHistName + "_IMP");
      if ( !imphists || imphists->numEntries() <= 1 ) {
        MSG_WARNING("No impact parameter calibration histogram for " <<
                    "CentralityProjection " << projName << " found " <<
                    "(requested histogram " << calHistName << "_IMP in " <<
                    calAnaName << ").\nDid you mean to generate one " <<
		    "yourself and forgot? Generate with " << calAnaName <<
		    " and preload the calibration file using the -p flag.");
      }
      else {
        MSG_INFO("Found calibration histogram " << sel << " " << imphists->path());
        cproj.add(PercentileProjection(ImpactParameterProjection(),
                                       *imphists, true), sel);
      }
    }
    else if ( sel == "USR" ) {
      #if HEPMC_VERSION_CODE >= 3000000
      YODA::Histo1DPtr usrhists =
        getPreload<Histo1D>("/" + calAnaName + "/" + calHistName + "_USR");
      if ( !usrhists || usrhists->numEntries() <= 1 ) {
        MSG_WARNING("No user-defined calibration histogram for " <<
                    "CentralityProjection " << projName << " found " <<
                    "(requested histogram " << calHistName << "_USR in " << calAnaName << ")");
        continue;
      } else {
        MSG_INFO("Found calibration histogram " << sel << " " << usrhists->path());
        cproj.add((UserCentEstimate(), usrhists*, true), sel);
      }
      #else
      MSG_ERROR("UserCentEstimate is only available with HepMC3.");
      #endif
    }
    else if ( sel == "RAW" ) {
      #if HEPMC_VERSION_CODE >= 3000000 || defined(RIVET_ENABLE_HEPMC_20610)
      cproj.add(GeneratedPercentileProjection(), sel);
      #else
      MSG_ERROR("GeneratedCentrality is only available with HepMC3 and HepMC 2.06.10.");
      #endif
    }
    else MSG_ERROR("'" << sel << "' is not a valid PercentileProjection tag.");

    if ( cproj.empty() ) {
      MSG_WARNING("CentralityProjection " << projName
                  << " did not contain any valid PercentileProjections.");
    }

    return declare(cproj, projName);
  }


  vector<string> Analysis::_weightNames() const {
    return handler().weightNames();
  }

  YODA::AnalysisObjectPtr Analysis::_getPreload(string path) const {
    return handler().getPreload(path);
  }

  size_t Analysis::defaultWeightIndex() const {
    return handler().defaultWeightIndex();
  }

  MultiweightAOPtr Analysis::_getOtherAnalysisObject(const std::string & ananame, const std::string& name) {
    std::string path = "/" + ananame + "/" + name;
    const auto& ana = handler().analysis(ananame);
    return ana->getAnalysisObject(name); //< @todo includeorphans check??
  }

  void Analysis::_checkBookInit() const {
    if (handler().stage() != AnalysisHandler::Stage::INIT) {
      MSG_ERROR("Can't book objects outside of init()");
      throw UserError(name() + ": Can't book objects outside of init().");
    }
  }

  bool Analysis::_inInit() const {
    return handler().stage() == AnalysisHandler::Stage::INIT;
  }

  bool Analysis::_inFinalize() const {
    return handler().stage() == AnalysisHandler::Stage::FINALIZE;
  }

  template <typename YODAT>
  void Analysis::_setWriterPrecision(const string& path, YODAT& yao) {
    const string re = _info->writerDoublePrecision();
    if (re != "") {
      std::smatch match;
      const bool needsDP = std::regex_search(path, match, std::regex(re));
      if (needsDP)  yao.template setAnnotation("WriterDoublePrecision", "1");
    }
  }

}
