// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief J/psi -> 3 gamma
  class BESIII_2013_I1126137 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2013_I1126137);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::pid==443), "UFS");
      // histos
      book(_h,1,1,1);
      book(_c,"TMP/nJPsi");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Find the J/psi among the unstables
      for (const Particle& psi : apply<UnstableParticles>(event, "UFS").particles()) {
	_c->fill();
	unsigned int nhadron(0);
	Particles photons;
	for(const Particle & child : psi.children()) {
	  if(PID::isHadron(child.pid()))
	    ++nhadron;
	  else if(child.pid()==PID::PHOTON)
	    photons.push_back(child);
	}
	if(nhadron==0 && photons.size()==3) {
	  LorentzTransform boost;
	  if (psi.p3().mod() > 1*MeV)
	    boost = LorentzTransform::mkFrameTransformFromBeta(psi.momentum().betaVec());
	  for(const Particle & gamma:photons)
	    _h->fill(2.*boost.transform(gamma.momentum()).E()/psi.mass());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h, 1e6/3. / *_c);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    CounterPtr _c;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2013_I1126137);

}
