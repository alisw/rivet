// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Thrust.hh"

namespace Rivet {


  /// @brief ALEPH measurement of quark-to-photon fragmentation function
  class ALEPH_1996_S3196992 : public Analysis {
  public:

    /// Constructor
    ALEPH_1996_S3196992()
      : Analysis("ALEPH_1996_S3196992")
    {    }


    /// @name Analysis methods
    //@{

    void init() {
      // Set up projections
      FinalState fs;
      declare(FastJets(fs, FastJets::DURHAM, 0.7), "DurhamJets");
      IdentifiedFinalState ifs; //(Cuts::pT > 0);
      ifs.acceptId(PID::PHOTON);
      declare(ifs, "Photons");
      declare(Thrust(fs), "Thrust");
      declare(ChargedFinalState(), "CFS");

      // Book histograms
      book(_h_z_2jet_001 ,1, 1, 1);
      book(_h_z_2jet_006 ,2, 1, 1);
      book(_h_z_2jet_01  ,3, 1, 1);
      book(_h_z_2jet_033 ,4, 1, 1);
      book(_h_z_3jet_001 ,5, 1, 1);
      book(_h_z_3jet_006 ,6, 1, 1);
      book(_h_z_3jet_01  ,7, 1, 1);
      book(_h_z_4jet_001 ,8, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      if (apply<FinalState>(event, "CFS").particles().size() < 2) vetoEvent;

      const Particles allphotons = apply<IdentifiedFinalState>(event, "Photons").particles();
      Particles photons;
      for (const Particle& photon : allphotons)
        if (fabs(cos(photon.theta())) < 0.95 && photon.E() > 5.0*GeV)
          photons.push_back(photon);
      if (photons.size() < 1) vetoEvent;

      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      if (fabs(cos(thrust.thrustAxis().theta()))>0.95) vetoEvent;

      const FastJets& durjet = apply<FastJets>(event, "DurhamJets");

      for (const Particle& photon : photons) {

        PseudoJets jets_001 = durjet.clusterSeq()->exclusive_jets_ycut(0.01);
        for (const fastjet::PseudoJet& jet : jets_001) {
          if (particleInJet(photon, jet)) {
            double zgamma = photon.E()/jet.E();
            if (jets_001.size() == 2) _h_z_2jet_001->fill(zgamma);
            else if (jets_001.size() == 3) _h_z_3jet_001->fill(zgamma);
            else if (jets_001.size() > 3) _h_z_4jet_001->fill(zgamma);
            break;
          }
        }

        PseudoJets jets_006 = durjet.clusterSeq()->exclusive_jets_ycut(0.06);
        for (const fastjet::PseudoJet& jet : jets_006) {
          if (particleInJet(photon, jet)) {
            double zgamma = photon.E()/jet.E();
            if (jets_006.size() == 2) _h_z_2jet_006->fill(zgamma);
            else if (jets_006.size() == 3) _h_z_3jet_006->fill(zgamma);
            break;
          }
        }

        PseudoJets jets_01 = durjet.clusterSeq()->exclusive_jets_ycut(0.1);
        for (const fastjet::PseudoJet& jet : jets_01) {
          if (particleInJet(photon, jet)) {
            double zgamma = photon.E()/jet.E();
            if (jets_01.size() == 2) _h_z_2jet_01->fill(zgamma);
            else if (jets_01.size() == 3) _h_z_3jet_01->fill(zgamma);
            break;
          }
        }

        PseudoJets jets_033 = durjet.clusterSeq()->exclusive_jets_ycut(0.33);
        for (const fastjet::PseudoJet& jet : jets_033) {
          if (particleInJet(photon, jet)) {
            double zgamma = photon.E()/jet.E();
            if (jets_033.size() == 2) _h_z_2jet_033->fill(zgamma);
            break;
          }
        }

      }
    }


    bool particleInJet(const Particle& p, const fastjet::PseudoJet& jet) {
      for (const fastjet::PseudoJet& jetpart : jet.constituents()) {
        if (fuzzyEquals(jetpart.E(), p.E()) &&
            fuzzyEquals(jetpart.px(), p.px()) &&
            fuzzyEquals(jetpart.py(), p.py()) &&
            fuzzyEquals(jetpart.pz(), p.pz())) {
          return true;
        }
      }
      return false;
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_z_2jet_001, 1000.0/sumOfWeights());
      scale(_h_z_2jet_006, 1000.0/sumOfWeights());
      scale(_h_z_2jet_01, 1000.0/sumOfWeights());
      scale(_h_z_2jet_033, 1000.0/sumOfWeights());
      scale(_h_z_3jet_001, 1000.0/sumOfWeights());
      scale(_h_z_3jet_006, 1000.0/sumOfWeights());
      scale(_h_z_3jet_01, 1000.0/sumOfWeights());
      scale(_h_z_4jet_001, 1000.0/sumOfWeights());
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_z_2jet_001, _h_z_2jet_006, _h_z_2jet_01, _h_z_2jet_033;
    Histo1DPtr _h_z_3jet_001, _h_z_3jet_006, _h_z_3jet_01;
    Histo1DPtr _h_z_4jet_001;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALEPH_1996_S3196992);


}
