// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief CP asymetry in B -> X gamma
  class BABAR_2015_I1337783 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2015_I1337783);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::pid==300553), "UFS");
      // profile hist for asymmetry
      book(_p,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Loop over upslion(4s)
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
      	// boost to rest frame
      	LorentzTransform cms_boost;
      	if (p.p3().mod() > 1*MeV)
      	  cms_boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	// loop over children
      	for(const Particle & p2 : p.children(Cuts::abspid==521 || Cuts::abspid==511)) {
	  Particle bottom=p2;
	  while(bottom.children()[0].abspid()==p2.abspid())
	    bottom = bottom.children()[0];
	  bool charm = false;
	  FourMomentum pgamma(0.,0.,0.,0.);
	  unsigned int ngamma = 0;
	  for (const Particle & child : bottom.children()) {
	    if (child.pid() == PID::PHOTON) {
	      ngamma += 1;
	      pgamma += child.momentum();
	    }
	    else if(PID::isCharmHadron(child.pid()))
	      charm = true;
	  }
	  if (ngamma != 1 || charm ) continue;
	  double Egamma = cms_boost.transform(pgamma).E();
	  double wgt = bottom.pid()<0 ? 100. : -100.;
	  for(const auto & bin : _p->bins())
	    if(bin.xMin()<Egamma) _p->fill(bin.xMid(),wgt);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
    }

    /// @}


    /// @name Histograms
    /// @{
    Profile1DPtr _p;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2015_I1337783);

}
