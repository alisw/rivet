// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Thrust.hh"

namespace Rivet {


  /// @brief Lambda0 production at 29 GeV
  class MARKII_1985_I209198 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MARKII_1985_I209198);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      const FinalState fs;
      declare(cfs, "FS");
      declare(Thrust(fs), "Thrust");
      //Histograms
      book(_h_spect,2,1,1);
      book(_h_pT   ,3,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // 5 charged particles
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      if(cfs.particles().size()<5) vetoEvent;
      // thrust
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      // lambdas
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==3122)) {
	const double xp = 2.*p.E()/sqrtS();
      	Vector3 mom3 = p.p3();
	const double beta = mom3.mod() / p.E();
        const double pTin  = dot(mom3, thrust.thrustMajorAxis());
        const double pTout = dot(mom3, thrust.thrustMinorAxis());
	double pT2 = sqr(pTin)+sqr(pTout);
	_h_spect->fill(xp,1./beta);
	_h_pT   ->fill(pT2);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_spect, sqr(sqrtS())*crossSection()/nanobarn/sumOfWeights());
      normalize(_h_pT);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_spect,_h_pT;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(MARKII_1985_I209198);

}
