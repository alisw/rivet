// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief p, Lambda, Lambda_c and Xi spectra at Upsilon(4s)
  class CLEO_1992_I315181 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CLEO_1992_I315181);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_c_ups,"/TMP/c_ups");
      book(_h_p       ,1,1,1);
      book(_h_Lambda_c,2,1,1);
      book(_h_Lambda  ,3,1,1);
      book(_h_Xi      ,4,1,1);
    }

    void findDecayProducts(Particle parent, Particles & protons, Particles & lambda_c, Particles & lambda, Particles & xi) {
      for(const Particle & p : parent.children()) {
	if (p.abspid() == PID::PROTON) {
	  protons.push_back(p);
	  continue;
	}
	else if(p.abspid()== PID::LAMBDA) {
	  lambda.push_back(p);
	  continue;
	}
	else if(p.abspid()== PID::XIMINUS) {
	  xi.push_back(p);
	}
	else if(p.abspid()== 4122) {
	  lambda_c.push_back(p);
	}
	if (!p.children().empty())
          findDecayProducts(p,protons,lambda_c,lambda,xi);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Find the upsilons
      for (const Particle& ups : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==300553)) {
	_c_ups->fill();
	Particles protons,lambda_c,lambda,xi;
	findDecayProducts(ups,protons,lambda_c,lambda,xi);
        LorentzTransform boost;
        if (ups.p3().mod() > 1*MeV)
          boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
	for(const Particle p : protons) {
          double pcm = boost.transform(p.momentum()).p3().mod();
	  _h_p->fill(pcm);
	}
	for(const Particle p : lambda_c) {
          double pcm = boost.transform(p.momentum()).p3().mod();
	  _h_Lambda_c->fill(pcm);
	}
	for(const Particle p : lambda) {
          double pcm = boost.transform(p.momentum()).p3().mod();
	  _h_Lambda->fill(pcm);
	}
	for(const Particle p : xi) {
          double pcm = boost.transform(p.momentum()).p3().mod();
	  _h_Xi->fill(pcm);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // br for lambda_c mode used (p K- pi+)
      double br=0.0623;
      scale(_h_p       ,0.5   / *_c_ups);
      scale(_h_Lambda_c,0.5*br/ *_c_ups);
      scale(_h_Lambda  ,0.5   / *_c_ups);
      scale(_h_Xi      ,0.5   / *_c_ups);
    }

    ///@}


    /// @name Histograms
    ///@{
    CounterPtr _c_ups;
    Histo1DPtr _h_p,_h_Lambda_c,_h_Lambda,_h_Xi;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(CLEO_1992_I315181);

}
