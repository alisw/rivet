// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class TOPAZ_1995_I381900 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TOPAZ_1995_I381900);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      const ChargedFinalState cfs;
      declare(cfs, "CFS");

      // Book histograms
      book(_h_charged,1, 1, 1);
      book(_h_pi     ,2, 1, 1);
      book(_h_Kp     ,2, 1, 2);
      book(_h_proton ,2, 1, 3);
      book(_h_K0     ,3, 1, 1);
      book(_wSum,"TMP/wSum");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      int nch = cfs.particles().size();
      if(nch<5) vetoEvent;
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      _wSum->fill();

      // neutral kaons
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for(const Particle & p : ufs.particles(Cuts::pid==130 || Cuts::pid==310)) {
	double xi = -log(p.p3().mod()/meanBeamMom);
	_h_K0->fill(xi);
      }
      // charged particles
      for(const Particle & p : cfs.particles()) {
	double xi = -log(p.p3().mod()/meanBeamMom);
	_h_charged->fill(xi);
	int id = abs(p.pid());
	if(id==211)
	  _h_pi->fill(xi);
	else if(id==321)
	  _h_Kp->fill(xi);
	else if(id==2212)
	  _h_proton->fill(xi);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_charged,1./ *_wSum);
      scale(_h_pi     ,1./ *_wSum);
      scale(_h_Kp     ,1./ *_wSum);
      scale(_h_proton ,1./ *_wSum);
      scale(_h_K0     ,1./ *_wSum);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_charged,_h_pi,_h_Kp,_h_proton,_h_K0;
    CounterPtr _wSum;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(TOPAZ_1995_I381900);


}
