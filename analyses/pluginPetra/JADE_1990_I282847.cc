// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief gamma, pi0 and eta spectra at 35 and 44 GeV
  class JADE_1990_I282847 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(JADE_1990_I282847);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      int ioff=-1;
      if(fuzzyEquals(sqrtS()/GeV,35.,1e-3)) {
      	ioff=0;
      }
      else if(fuzzyEquals(sqrtS()/GeV,44.,1e-3)) {
      	ioff=1;
      }
      else
      	MSG_ERROR("Beam energy " << sqrtS() << " not supported!");
      // Book histograms
      book(_h_gamma, 1+ioff, 1, 1);
      book(_h_pi0  , 3+ioff, 1, 1);
      if(ioff==0) book(_h_eta, 5+ioff, 1, 1);
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
      	double modp = p.p3().mod();
      	double xE = p.E()/meanBeamMom;
      	double beta = modp/p.E();
      	if(p.pid()==111)
      	  _h_pi0->fill(xE,1./beta);
      	else if(_h_eta)
      	  _h_eta->fill(xE,1./beta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_gamma, crossSection()*sqr(sqrtS())/microbarn/sumOfWeights());
      scale(_h_pi0  , crossSection()*sqr(sqrtS())/microbarn/sumOfWeights());
      if(_h_eta) scale(_h_eta  , crossSection()*sqr(sqrtS())/microbarn/sumOfWeights());
    }

    //@}

    /// @name Histograms
    //@{
    Histo1DPtr _h_gamma, _h_pi0, _h_eta;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(JADE_1990_I282847);


}
