// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B->D semileptonic decays
  class BELLE_2020_I1796822 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2020_I1796822);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // Book histograms
      book(_h_mpipi,1,1,1);
      book(_h_q2   ,2,1,1);
      book(_nB,"/TMP/nB");
    }

    void findChildren(const Particle & p, unsigned int &ncount,
		      Particles & pi, Particles & ell, Particles & nu) {
      _nB->fill();
      for (const Particle &child : p.children()) {
	if(child.children().empty()) {
	  if(child.abspid()==211) {
	    ++ncount;
	    pi.push_back(child);
	  }
	  else if(child.abspid()==11 || child.abspid()==13) {
	    ++ncount;
	    ell.push_back(child);
	  }
	  else if(child.abspid()==12 || child.abspid()==14) {
	    ++ncount;
	    nu.push_back(child);
	  }
	  else if(child.pid()!=22)
	    ++ncount;
	}
	// veto gamma gamma decaying mesons and K0
	else if(child.pid()==111 || child.pid()==221 || child.pid()==331 ||
		child.pid()==130 || child.pid()==310) {
	  ++ncount;
	}
	else
	  findChildren(child,ncount,pi,ell,nu);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==PID::BPLUS)) {
	Particles pi,ell,nu;
	unsigned int ncount = 0;
	findChildren(p,ncount,pi,ell,nu);
	// check the decay
	// 4 outgoing
	if(ncount!=4) continue;
	// including pi+ pi-
	if(pi.size()!=2 || pi[0].pid() != -pi[1].pid()) continue;
	// and ell- nubar or ell+ nu
	if(ell.size()!=1 || nu.size()!=1) continue;
	int inu = ell[0].abspid()+1;
	if(ell[0].pid()>0) inu *=-1;
	if(nu[0].pid()!=inu) continue;
	// fill histos
	// m pipi
	FourMomentum ppi = pi[0].momentum()+pi[1].momentum();
	_h_mpipi->fill(ppi.mass());
	// q2
	FourMomentum q = p.momentum() - ppi;
	_h_q2->fill(q.mass2());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_mpipi,1e5/ *_nB);
      scale(_h_q2   ,1e5/ *_nB);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_mpipi,_h_q2;
    CounterPtr _nB;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2020_I1796822);

}
