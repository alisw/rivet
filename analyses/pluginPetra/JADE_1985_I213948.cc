// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief  gamma, pi0 and eta spectra at 14, 22.5 and 34.4 GeV
  class JADE_1985_I213948 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(JADE_1985_I213948);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // find the beam energy
      int ioff=-1;
      if(fuzzyEquals(sqrtS()/GeV,34.5,1e-3)) {
      	ioff=0;
      }
      else if(fuzzyEquals(sqrtS()/GeV,22.5,1e-3)) {
      	ioff=1;
      }
      else if(fuzzyEquals(sqrtS()/GeV,14.0,1e-3)) {
      	ioff=2;
      }
      else
      	MSG_ERROR("Beam energy " << sqrtS() << " not supported!");
      // book histos
      book(_h_gamma,ioff+1,1,1);
      book(_h_pi0  ,ioff+4,1,1);
      if(ioff==0)
	book(_h_eta,7,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      // gamma
      const FinalState& fs = apply<FinalState>(event, "FS");
      for (const Particle& p : fs.particles(Cuts::pid==22)) {
      	double xE = p.E()/meanBeamMom;
      	_h_gamma->fill(xE);
      }
      // pi0, eta
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::pid==111 or Cuts::pid==221)) {
      	double xE = p.E()/meanBeamMom;
      	if(p.pid()==111)
      	  _h_pi0->fill(xE);
      	else if(_h_eta)
      	  _h_eta->fill(xE);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_gamma, crossSection()*sqr(sqrtS())/microbarn/sumOfWeights());
      scale(_h_pi0  , crossSection()*sqr(sqrtS())/microbarn/sumOfWeights());
      if(_h_eta) scale(_h_eta  , crossSection()*sqr(sqrtS())/microbarn/sumOfWeights());
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_gamma, _h_pi0, _h_eta;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(JADE_1985_I213948);

}
