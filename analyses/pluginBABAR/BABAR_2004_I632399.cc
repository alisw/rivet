// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief phi spectrum at 4S
  class BABAR_2004_I632399 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2004_I632399);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_b_phi, 1, 1, 1);
      book(_h_phi, 2, 1, 1);
      book(_c_4S, "/TMP/N4S");
    }

    void findDecayProducts(const Particle & p, Particles & phi) {
      for(const Particle & child : p.children()) {
	if(child.pid()==333)
	  phi.push_back(child);
	else if(!child.children().empty())
	  findDecayProducts(child,phi);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Find the Upsilons among the unstables
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particles upsilons = ufs.particles(Cuts::pid==300553);
      for (const Particle& ups : upsilons) {
        LorentzTransform cms_boost;
        if (ups.p3().mod() > 1*MeV)
          cms_boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
	Particles phis;
        findDecayProducts(ups, phis);
	_c_4S->fill();

	for(const Particle & phi : phis) {
          FourMomentum p2 = cms_boost.transform(phi.momentum());
	  _h_phi->fill(p2.p3().mod());
	  _b_phi->fill(0.5);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      if(_c_4S->effNumEntries()!=0) {
	scale(_h_phi   , 0.5/ *_c_4S);
	scale(_b_phi   , 50./ *_c_4S );
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_phi,_b_phi;
    CounterPtr _c_4S;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2004_I632399);

}
