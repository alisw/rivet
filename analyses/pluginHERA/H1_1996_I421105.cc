// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/UnstableParticles.hh"
namespace Rivet {


  /// @brief Inclusive D0 and D*+- production in deep inelastic e p scattering at HERA (H1)
  class H1_1996_I421105 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(H1_1996_I421105);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
    
      
      declare(DISKinematics(), "Kinematics");
      declare(UnstableParticles(), "UFS");
 
      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 4.9);
      declare(fs, "FS");

      // Book histograms
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)



      book(_h["p_tD*_norm"], 4, 1, 1);
      book(_h["p_tD*"], 4, 1, 2);
      book(_h["p_tD0_norm"], 5, 1, 1);
      book(_h["p_tD0"], 5, 1, 2);
      book(_h["xD_D*_norm"], 6, 1, 1);
      book(_h["xD_D*"], 6, 1, 2);
      book(_h["xD_D0_norm"], 7, 1, 1);
      book(_h["xD_D0"], 7, 1, 2);

    }


/// Perform the per-event analysis
    void analyze(const Event& event) {

      /// @todo Do the event by event analysis here
      const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");
      const LorentzTransform hcmboost = dk.boostHCM();

      // Get the DIS kinematics
      double y = dk.y();
      double W2 = dk.W2()/GeV2;
      double Q2 = dk.Q2()/GeV;

     bool cut ;
     cut = Q2 > 10 && Q2 < 100 && y > 0.01 && y < 0.7 ;
     
     if (! cut ) vetoEvent ;
     bool etacut ;
     for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
        etacut = abs(p.momentum().eta()) < 1.5 ;
        const FourMomentum hcmMom = hcmboost.transform(p.momentum());
        double p_D ;
        double x_D = 0 ;
        p_D = std::sqrt( hcmMom.px()*hcmMom.px() + hcmMom.py()*hcmMom.py() + hcmMom.pz()*hcmMom.pz() );
        x_D = 2.*p_D/sqrt(W2);
        if (p.abspid() == 421) {
           _h["p_tD0"]->fill(p.momentum().pT()/GeV);
           _h["p_tD0_norm"]->fill(p.momentum().pT()/GeV);
           if (etacut ) _h["xD_D0"]->fill(x_D);
           if (etacut ) _h["xD_D0_norm"]->fill(x_D);
        }
          
        if (p.abspid() == 413) {
          _h["p_tD*"]->fill(p.momentum().pT()/GeV);
          _h["p_tD*_norm"]->fill(p.momentum().pT()/GeV);
          // x_D is defined for the D0
          if (etacut ) _h["xD_D*"]->fill(x_D);
          if (etacut ) _h["xD_D*_norm"]->fill(x_D);
        } 
    }

  }
    /// Normalise histograms etc., after the run
    void finalize() {
        
    
     scale(_h["p_tD*"], crossSection()/nanobarn/sumW()); 
     scale(_h["p_tD0"], crossSection()/nanobarn/sumW()); 
     normalize(_h["p_tD*_norm"]);
     normalize(_h["p_tD0_norm"]);
     
     scale(_h["xD_D*"], crossSection()/nanobarn/sumW()); 
     scale(_h["xD_D0"], crossSection()/nanobarn/sumW()); 
     normalize(_h["xD_D*_norm"]);
     normalize(_h["xD_D0_norm"]);


    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    ///@}

  };


  RIVET_DECLARE_PLUGIN(H1_1996_I421105);

}
