// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief K0 and K*+ spectra
  class TASSO_1990_I284251 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TASSO_1990_I284251);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      declare(Sphericity(cfs), "Sphericity");


      // Book histograms
      _ih=-1; _iy=-1;
      if(isCompatibleWithSqrtS(14.8)) {
	_ih=1;
      }
      else if (isCompatibleWithSqrtS(21.5)) {
	_ih=2;
      }
      else if (isCompatibleWithSqrtS(34.5)) {
	_ih=0;
	_iy=3;
      }
      else if (isCompatibleWithSqrtS(35.0)) {
	_ih=0;
	_iy=2;
      }
      else if (isCompatibleWithSqrtS(42.6)) {
	_ih=0;
	_iy=1;
      }
      else
	MSG_ERROR("Beam energy " << sqrtS() << " not supported!");

      if(_ih==0) {
	book(_h_K0_x, 1,1,_iy);
	if(_iy!=3) {
	  book(_p_K0_S_1, 5,1,2*_iy-1);
	  book(_p_K0_S_2,"TMP/p_K0_S_2",refData(_ih+5,1,2*_iy));
	}
	book(_h_Kstar_x, 8,1,_iy);
	if(_iy==2) {
	  book(_p_Kstar_S_1, 10,1,1);
	  book(_p_Kstar_S_2,"TMP/p_Kstar_S_2",refData(10,1,2));
	}
      }
      else {
	book(_h_K0_x  , _ih+1,1,1);
	book(_p_K0_S_1, _ih+5,1,1);
	book(_p_K0_S_2,"TMP/p_K0_S_2",refData(_ih+5,1,2));
      }
      book(_n_K0   ,"/TMP/nK0"   );
      book(_n_Kstar,"/TMP/nKstar");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      const size_t numParticles = cfs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
				   beams.second.p3().mod() ) / 2.0;
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");

      unsigned int nK0(0),nKstar(0);
      UnstableParticles ufs = apply<UnstableParticles>(event,"UFS");
      for(const Particle & p : ufs.particles(Cuts::abspid==323 or Cuts::pid==130 or Cuts::pid==310)) {
	double xE = p.E()/meanBeamMom;
	double modp = p.p3().mod();
	double beta = modp/p.E();
	if(abs(p.pid())==323) {
	  if(_h_Kstar_x!=Histo1DPtr()) _h_Kstar_x->fill(xE,1./beta);
	  ++nKstar;
	}
	else {
	  if(_h_K0_x!=Histo1DPtr()) _h_K0_x->fill(xE,1./beta);
	  ++nK0;
	}
      }
      _n_K0->fill(nK0);
      _n_Kstar->fill(nKstar);
      double sphere = sphericity.sphericity();
      if(_p_K0_S_1!=Profile1DPtr()) {
	_p_K0_S_1->fill(sphere,nK0);
	_p_K0_S_2->fill(sphere,cfs.particles().size());
      }
      if(_p_Kstar_S_1!=Profile1DPtr()) {
	_p_Kstar_S_1->fill(sphere,nKstar);		    
	_p_Kstar_S_2->fill(sphere,cfs.particles().size());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_K0_x, sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
      if(_h_Kstar_x!=Histo1DPtr())
	scale(_h_Kstar_x, sqr(sqrtS())*crossSection()/nanobarn/sumOfWeights());
      if(_p_K0_S_1!=Profile1DPtr()) {
	Scatter2DPtr temp;
	if(_ih==0)
	  book(temp,5,1,2*_iy);
	else
	  book(temp,_ih+5,1,2);
	divide(_p_K0_S_1,_p_K0_S_2,temp);
      }
      if(_p_Kstar_S_1!=Profile1DPtr()) {
	Scatter2DPtr temp;
	book(temp,10,1,2);
	divide(_p_Kstar_S_1,_p_Kstar_S_2,temp);
      }
      // K0 mult
      scale(_n_K0   ,1./sumOfWeights());
      Scatter2D temphisto(refData(4, 1, 1));
      Scatter2DPtr mult;
      book(mult, 4, 1, 1);
      for (size_t b = 0; b < temphisto.numPoints(); b++) {
      	const double x  = temphisto.point(b).x();
      	pair<double,double> ex = temphisto.point(b).xErrs();
      	pair<double,double> ex2 = ex;
     	if(ex2.first ==0.) ex2. first=0.0001;
     	if(ex2.second==0.) ex2.second=0.0001;
      	if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second))
       	  mult   ->addPoint(x, _n_K0->val(), ex, make_pair(_n_K0->err(),_n_K0->err()));
	else
	  mult   ->addPoint(x, 0., ex, make_pair(0.,.0));
      }
      // K*= mult
      scale(_n_Kstar,1./sumOfWeights());
      Scatter2D temphisto2(refData(9, 1, 1));
      Scatter2DPtr mult2;
      book(mult2,9, 1, 1);
      for (size_t b = 0; b < temphisto2.numPoints(); b++) {
      	const double x  = temphisto2.point(b).x();
      	pair<double,double> ex = temphisto2.point(b).xErrs();
      	pair<double,double> ex2 = ex;
     	if(ex2.first ==0.) ex2. first=0.0001;
     	if(ex2.second==0.) ex2.second=0.0001;
      	if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second))
       	  mult2   ->addPoint(x, _n_Kstar->val(), ex, make_pair(_n_Kstar->err(),_n_Kstar->err()));
	else
	  mult2   ->addPoint(x, 0., ex, make_pair(0.,.0));
      }
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_K0_x, _h_Kstar_x;
    Profile1DPtr _p_K0_S_1, _p_K0_S_2, _p_Kstar_S_1, _p_Kstar_S_2;
    CounterPtr _n_K0,_n_Kstar;
    int _ih,_iy;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(TASSO_1990_I284251);


}
