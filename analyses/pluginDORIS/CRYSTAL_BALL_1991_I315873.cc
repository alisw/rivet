// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief direct photon spectrum in upslion decay
  class CRYSTAL_BALL_1991_I315873 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CRYSTAL_BALL_1991_I315873);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(), "UFS");
      book(_h_gamma[0],2,1,1);
      book(_h_gamma[1],2,1,2);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Find the Upsilons among the unstables
      for (const Particle& ups : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==553)) {
	unsigned int nhadron(0);
	Particles photons;
	for(const Particle & child : ups.children()) {
	  if(PID::isHadron(child.pid()))
	    ++nhadron;
	  else if(child.pid()==PID::PHOTON)
	    photons.push_back(child);
	}
	if(nhadron!=0 && !photons.empty()) {
	  LorentzTransform boost;
	  if (ups.p3().mod() > 1*MeV)
	    boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
	  for(const Particle & gamma:photons) {
	    double z = 2.*boost.transform(gamma.momentum()).E()/ups.mass();
	    _h_gamma[0]->fill(z);
	    _h_gamma[1]->fill(z);
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_gamma[0],1.,false);
      normalize(_h_gamma[1],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_gamma[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CRYSTAL_BALL_1991_I315873);

}
