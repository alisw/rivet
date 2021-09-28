// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief J/Psi chi_c charged particle multiplicities
  class BESIII_2020_I1800404 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BESIII_2020_I1800404);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::pid==443 or Cuts::pid==445 or Cuts::pid==10441 or Cuts::pid==20443),"UFS");
      // histograms
      book(_h_chi0,1,1,1);
      book(_h_chi1,1,1,2);
      book(_h_chi2,1,1,3);
      book(_h_jpsi[0],2,1,1);
      book(_h_jpsi[1],2,1,2);
    }

    void findChildren(const Particle & p, int & nCharged) {
      for( const Particle &child : p.children()) {
	if(child.children().empty()) {
	  if(PID::isCharged(child.pid())) ++nCharged;
	}
	else
	  findChildren(child,nCharged);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over particles
      for( const Particle & p : apply<UnstableParticles>(event, "UFS").particles()) {
	// skip radiative modes
	if(p.children().size()==2) {
	  if(p.children()[0].pid()==PID::GAMMA) {
	    if(p.children()[1].pid()==PID::GAMMA ||
	       p.children()[1].pid()==PID::JPSI) continue;
	  }
	  if(p.children()[1].pid()==PID::GAMMA) {
	    if(p.children()[0].pid()==PID::GAMMA ||
	       p.children()[0].pid()==PID::JPSI) continue;
	  }
	}
	// get the charged particle multiplicity
	int nCharged(0);
	findChildren(p,nCharged);
	if(p.pid()==PID::JPSI) {
	  _h_jpsi[0]->fill(nCharged);
	  _h_jpsi[1]->fill(nCharged);
	}
	else if(p.pid()==10441)
	  _h_chi0->fill(nCharged);
	else if(p.pid()==20443)
	  _h_chi1->fill(nCharged);
	else if(p.pid()==445)
	  _h_chi2->fill(nCharged);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // percentage and bin width
      normalize(_h_jpsi[0],200.);
      normalize(_h_jpsi[1],200.);
      normalize(_h_chi0   ,200.);
      normalize(_h_chi1   ,200.);
      normalize(_h_chi2   ,200.);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_jpsi[2],_h_chi0,_h_chi1,_h_chi2;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(BESIII_2020_I1800404);

}
