// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Event.hh"

namespace Rivet {


  class EHS_1988_I265504 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(EHS_1988_I265504);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      declare(ChargedFinalState(), "CFS");
      declare(Beam(),"Beam");

      switch ( beamIds().first ) {
      case PID::PIPLUS:
        book(_h_cpos_xF ,1, 1, 1);
        book(_h_cpos_eta ,3, 1, 1);
        book(_h_cpos_pT2 ,5, 1, 1);
        book(_h_cneg_xF ,2, 1, 1);
        book(_h_cneg_eta ,4, 1, 1);
        book(_h_cneg_pT2 ,6, 1, 1);
        break;

      case PID::KPLUS:
        book(_h_cpos_xF ,1, 1, 2);
        book(_h_cpos_eta ,3, 1, 2);
        book(_h_cpos_pT2 ,5, 1, 2);
        book(_h_cneg_xF ,2, 1, 2);
        book(_h_cneg_eta ,4, 1, 2);
        book(_h_cneg_pT2 ,6, 1, 2);
        break;

      case PID::PROTON:
        book(_h_cpos_xF ,1, 1, 3);
        book(_h_cpos_eta ,3, 1, 3);
        book(_h_cpos_pT2 ,5, 1, 3);
        book(_h_cneg_xF ,2, 1, 3);
        book(_h_cneg_eta ,4, 1, 3);
        book(_h_cneg_pT2 ,6, 1, 3);
        break;
      }

      // Calculate boost from lab to CM frame
      _beamboost = cmsTransform( beams() );
      MSG_DEBUG("Boost vector: " << _beamboost );

      // Transform beam into CMS frame
      Particle _beam_cm = beams().first;
      _beam_cm.transformBy(_beamboost);
      // Beam momentum in CM frame defines Feynman-x
      _pz_max = _beam_cm.pz();

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const FinalState& fs = apply<FinalState>(event, "CFS");
      for (const Particle& p: fs.particles()) {
        // Only interested in pi- or positively charged
        if (p.charge() < 0 && p.pid() != PID::PIMINUS) continue;
        // Slow proton cut: reject lab momenta < 1.2GeV
        if (p.pid() == PID::PROTON && p.p() < 1.2*GeV) continue;
        // Transform to cm frame
        const FourMomentum pcm = _beamboost.transform(p);
        const double xF = pcm.pz()/_pz_max;

        if (p.charge() > 0) {
          _h_cpos_xF->fill( xF );
          _h_cpos_pT2->fill( p.pT2() );
          _h_cpos_eta->fill( pcm.eta() );
        } else if (p.pid() == PID::PIMINUS) {
          _h_cneg_xF->fill( xF );
          _h_cneg_pT2->fill( p.pT2() );
          _h_cneg_eta->fill( pcm.eta() );
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSection()/millibarn/sumOfWeights();
      scale(_h_cpos_xF, sf); scale(_h_cpos_pT2, sf);
      scale(_h_cpos_eta, sf); scale(_h_cneg_xF, sf);
      scale(_h_cneg_eta, sf); scale(_h_cneg_pT2, sf);
    }

    //@}


    /// @name Histograms
    //@{
    LorentzTransform _beamboost;
    double _pz_max;
    Histo1DPtr _h_cpos_xF, _h_cpos_eta, _h_cpos_pT2;
    Histo1DPtr _h_cneg_xF, _h_cneg_eta, _h_cneg_pT2;
    //@}

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(EHS_1988_I265504);

}
