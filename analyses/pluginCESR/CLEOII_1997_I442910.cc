// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Xi_c0 Xi_c+ spectrum in upsilon(4s) decays
  class CLEOII_1997_I442910 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CLEOII_1997_I442910);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(), "UFS");
      // book histos
      book(_h_Xi_c0   ,1,1,1);
      book(_h_Xi_cPlus,2,1,1);
    }

    void findDecayProducts(Particle parent, Particles & Xic) {
      for(const Particle & p : parent.children()) {
        int id = abs(p.pid());
	if(id==4132 || id==4232) {
	  Xic.push_back(p);
	}
	else if(!p.children().empty()) {
	  findDecayProducts(p,Xic);
	}
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Find the upsilons
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==300553)) {
        // Find the decay products we want
	Particles Xic;
        findDecayProducts(p, Xic);
	if(Xic.empty()) continue;
        LorentzTransform boost;
        if (p.p3().mod() > 1*MeV)
          boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	
	for(const Particle & xi : Xic) {
	  double Emax = sqrt(0.25*sqr(p.mass())-sqr(xi.mass()));
	  double xp = boost.transform(xi.momentum()).vector3().mod()/Emax;
	  if(xi.abspid()==4132)
	    _h_Xi_c0->fill(xp);
	  else
	    _h_Xi_cPlus->fill(xp);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // normalize to the data as norm in hepdata is weird
      normalize(_h_Xi_c0   ,0.144);
      normalize(_h_Xi_cPlus,0.453);

    }

    //@}

    /// @name Histograms
    //@{
    Histo1DPtr _h_Xi_c0, _h_Xi_cPlus;
    //@}

  };


  DECLARE_RIVET_PLUGIN(CLEOII_1997_I442910);

}
