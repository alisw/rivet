// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Sphericity.hh"

namespace Rivet {


  /// @brief Event shapes at Upsilon(1S)
  class ARGUS_1986_I227324 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1986_I227324);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(), "UFS");
      declare(ChargedFinalState(), "CFS");
      const FinalState fs;
      declare(Thrust(fs)    ,"Thrust");
      declare(Sphericity(fs),"Sphericity");
      // histograms
      if(isCompatibleWithSqrtS(9.98,1e-2)) {
        book(_h_T_cont ,2, 1, 2);
        book(_h_S_cont ,1, 1, 2);
      }
      book(_h_T_Ups ,2, 1, 1);
      book(_h_S_Ups ,1, 1, 1);
    }

    /// Recursively walk the decay tree to find the stable decay products of @a p
    void findDecayProducts(Particle mother, Particles& charged, Particles & neutral) {
      for(const Particle & p: mother.children()) {
	if(!p.children().empty())
	  findDecayProducts(p, charged, neutral);
	else {
	  if(isCharged(p))
	    charged.push_back(p);
	  else
	    neutral.push_back(p);
	}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Find the Upsilons among the unstables
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particles upsilons = ufs.particles(Cuts::pid==553 or Cuts::pid==100553);
      if (upsilons.empty() && _h_T_cont) {
	Particles charged = apply<ChargedFinalState>(event, "CFS").particles(); 
	// at least 6 charged particles
	if(charged.size()<6) vetoEvent;
	// cut on high momentum particles
	unsigned int nHigh(0);
	for(const Particle & p : charged) {
	  if(p.momentum().p3().mod()>2.5) ++nHigh;
	}
	if(nHigh>1) vetoEvent;
        MSG_DEBUG("No Upsilons found => continuum event");
	Thrust thrust = apply<Thrust>(event, "Thrust");
	_h_T_cont->fill(thrust.thrust());
	Sphericity sphericity = apply<Sphericity>(event, "Sphericity");
	_h_S_cont->fill(sphericity.sphericity());
      }
      else {
        for (const Particle& ups : upsilons) {
          LorentzTransform boost;
          if (ups.p3().mod() > 1*MeV)
            boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
          // Find the decay products we want
          Particles charged,neutral;
	  // 6 charged particles
          findDecayProducts(ups, charged, neutral);
	  if(charged.size()<6) continue;
	  // at most 1 |p|>2.5
	  vector<FourMomentum> mom;
	  mom.reserve(neutral.size()+charged.size());
	  unsigned int nHigh(0);
	  for(const Particle & p : charged) {
	    mom.push_back(boost.transform(p.momentum()));
	    if(mom.back().p3().mod()>2.5) ++nHigh;
	  }
	  if(nHigh>1) continue;
	  for(const Particle & p : neutral) {
	    mom.push_back(boost.transform(p.momentum()));
	  }
	  Thrust thrust;
	  thrust.calc(mom);
	  _h_T_Ups->fill(thrust.thrust());
	  Sphericity sphericity;
	  sphericity.calc(mom);
	  _h_S_Ups->fill(sphericity.sphericity());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      if(_h_T_cont) {
	normalize(_h_T_cont);
	normalize(_h_S_cont);
      }
      if(_h_T_Ups->numEntries()!=0.) {
	normalize(_h_T_Ups);
	normalize(_h_S_Ups);
      }
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_T_Ups,_h_T_cont;
    Histo1DPtr _h_S_Ups,_h_S_cont;
    //@}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1986_I227324);

}
