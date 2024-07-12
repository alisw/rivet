// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/Random.hh"

namespace Rivet {


  /// @brief Direct photon spectrum in upslion decay
  class CLEO_1986_I220652 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_1986_I220652);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::pid==553), "UFS");
      // histos
      book(_h,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Find the Upsilons among the unstables
      for (const Particle& ups : apply<UnstableParticles>(event, "UFS").particles()) {
      	unsigned int nhadron(0);
      	Particles photons;
      	for(const Particle & child : ups.children()) {
      	  if(PID::isHadron(child.pid()))
      	    ++nhadron;
      	  else if(child.pid()==PID::PHOTON)
      	    photons.push_back(child);
      	}
	if(nhadron==0 || photons.empty()) continue;
	LorentzTransform boost;
	if (ups.p3().mod() > 1*MeV)
	  boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
	for(const Particle & gamma:photons) {
	  double E = boost.transform(gamma.momentum()).E();
	  E = randnorm(E,0.21*sqrt(E));
	  _h->fill(2.*E/ups.mass());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h,1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_1986_I220652);

}
