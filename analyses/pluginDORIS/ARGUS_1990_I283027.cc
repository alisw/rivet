// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief anti-deuteron spectrum in upslion decays
  class ARGUS_1990_I283027 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1990_I283027);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_h_p,1,1,1);
      book(_w_ups,"TMP/w_ups");
    }


    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& deut) {
      // deuteron id code
      static const long id = -1000010020;
      for(const Particle & p: mother.children()) {
	if(p.pid() == id) {
	  deut.push_back(p);
	}
	else if(!p.children().empty())
	  findDecayProducts(p, deut);
      }
    }

    
    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle & ups : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==553) ) {
	_w_ups->fill();
	Particles deut;
	findDecayProducts(ups, deut);
	if(deut.empty()) continue;
	LorentzTransform boost;
	if (ups.p3().mod() > 1*MeV)
	  boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
	for(const Particle& p : deut) {
	  double mom = boost.transform(p.momentum()).p3().mod();
	  _h_p->fill(mom);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // direct upsilon decays, i.e. to gg gamma and ggg so need to renormalize
      // factor and values from arxiv:0612019
      double Bmumu=0.0249, rhad = 3.56;
      double fact = 1.-(3.+rhad)*Bmumu;
      scale(_h_p, 1e5/fact / *_w_ups);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_p;
    CounterPtr _w_ups;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1990_I283027);

}
