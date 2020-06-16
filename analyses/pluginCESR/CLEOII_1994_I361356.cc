// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Sigma_c 0,++ spectra
  class CLEOII_1994_I361356 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CLEOII_1994_I361356);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(), "UFS");
      // book histos
      book(_b_Sigma_pp,1,1,1);
      book(_b_Sigma_0 ,1,1,2);
      book(_h_Sigma_pp,2,1,1);
      book(_h_Sigma_0 ,2,1,2);
      book(_c_ups,"TMP/c_ups");
    }

    void findDecayProducts(Particle parent, Particles & Sigma) {
      for(const Particle & p : parent.children()) {
        int id = abs(p.pid());
	if(id==4112 || id==4222) {
	  Sigma.push_back(p);
	}
	else if(!p.children().empty()) {
	  findDecayProducts(p,Sigma);
	}
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Find the upsilons
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==300553)) {
	_c_ups->fill();
        // Find the decay products we want
	Particles Sigma;
        findDecayProducts(p, Sigma);
	if(Sigma.empty()) continue;
        LorentzTransform boost;
        if (p.p3().mod() > 1*MeV)
          boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());	
	for(const Particle & sig : Sigma) {
	  double mom = boost.transform(sig.momentum()).vector3().mod();
	  if(sig.abspid()==4222) {
	    _h_Sigma_pp->fill(mom);
	    _b_Sigma_pp->fill(0.5);
	  }
	  else {
	    _h_Sigma_0->fill(mom);
	    _b_Sigma_0->fill(0.5);
	  }
	}
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_Sigma_0 , 0.5/ *_c_ups);
      scale(_h_Sigma_pp, 0.5/ *_c_ups);
      scale(_b_Sigma_0 , 0.5/ *_c_ups);
      scale(_b_Sigma_pp, 0.5/ *_c_ups);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_Sigma_0, _h_Sigma_pp;
    Histo1DPtr _b_Sigma_0, _b_Sigma_pp;
    CounterPtr _c_ups;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(CLEOII_1994_I361356);

}
