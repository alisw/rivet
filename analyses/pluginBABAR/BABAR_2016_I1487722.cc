// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B0 -> D*- pi+ pi+ pi-
  class BABAR_2016_I1487722 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2016_I1487722);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      DecayedParticles B0(ufs);
      B0.addStable(PID::PI0);
      B0.addStable( 413);
      B0.addStable(-413);
      B0.addStable( 423);
      B0.addStable(-423);
      B0.addStable(PID::PI0);
      declare(B0, "B0");
      // histograms
      book(_h,1,1,1);
      // efficiency
      const Scatter2D& ref = refData(2,1,1);
      _edges.push_back(ref.points()[0].xMin());
      for(auto p : ref.points() ) {
	_edges.push_back(p.xMax());
	_eff.push_back(p.y());
      }
    }

    double eff(double mass) {
      if(mass<_edges[0] || mass>_edges.back()) return 0.;
      for(unsigned int ix=0;ix<_eff.size();++ix) {
	if(mass<_edges[ix+1]) return _eff[ix];
      }
      return 0.;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { -413,1}, { 211,2}, {-211,1} };
      static const map<PdgId,unsigned int> & mode1CC = { {  413,1}, {-211,2}, { 211,1} };
      DecayedParticles B0 = apply<DecayedParticles>(event, "B0");
      for(unsigned int ix=0;ix<B0.decaying().size();++ix) {
	int sign = B0.decaying()[ix].pid()/B0.decaying()[ix].abspid();
	if ( (sign== 1 && B0.modeMatches(ix,4,mode1) ) ||
	     (sign==-1 && B0.modeMatches(ix,4,mode1CC) ) ) {
	  FourMomentum ptotal;
	  for(const Particle & p : B0.decayProducts()[ix].at( sign*211) ) {
	    ptotal+=p.momentum();
	  }
	  for(const Particle & p : B0.decayProducts()[ix].at(-sign*211) ) {
	    ptotal+=p.momentum();
	  }
	  double m = ptotal.mass();
	  _h->fill(m, eff(m));
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    vector<double> _eff;
    vector<double> _edges;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2016_I1487722);

}
