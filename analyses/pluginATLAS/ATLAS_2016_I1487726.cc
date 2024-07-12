// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


class ATLAS_2016_I1487726 : public Analysis {

    public:

    /// Constructor
    /// @brief Collinear W emissions at 8 TeV
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2016_I1487726);

    public:

        /// @name Analysis methods
        //@{

        /// Book histograms and initialise projections before the run
        void init() {

            _mode = 0;
            if ( getOption("LMODE") == "EL" )  _mode = 1;

            // These really should include non-prompt leptons/photons
            FinalState mufs(Cuts::abspid == PID::MUON);
            FinalState elfs(Cuts::abspid == PID::ELECTRON);
            FinalState phs(Cuts::abspid == PID::PHOTON);

            Cut lep_fid = (Cuts::abseta < 2.4 && Cuts::pT >= 25*GeV);
            DressedLeptons dlep(phs, _mode? elfs : mufs, 0.1, lep_fid, true);
            declare(dlep, "DressedLeptons");

            FastJets fj(FinalState(), FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
            declare(fj, "AntiKt4Jets");

            book(h_mu_jet_dr,          2, 1, 1);
            book(h_mu_jet_dr_pt500600, 4, 1, 1);
            book(h_mu_jet_dr_pt650,    5, 1, 1);
        }


        /// Perform the per-event analysis
        void analyze(const Event& event) {

          const vector<DressedLepton> leptons = apply<DressedLeptons>(event, "DressedLeptons").dressedLeptons();
          const Jets jets = apply<FastJets>(event, "AntiKt4Jets").jetsByPt(Cuts::pT >= 100*GeV && Cuts::abseta <= 2.1);

          if (leptons.size() != 1)       vetoEvent;
          if (jets.size() < 1)           vetoEvent;
          if (jets[0].pt() < 500.0*GeV)  vetoEvent;

          // find closest jet to the lepton.
          Jet jet;
          double drmin = 999;
          for (const Jet &j : jets) {
            double dr = deltaR(leptons[0], j);
            if (dr < drmin) {
                drmin = dr;
                jet = j;
            }
          }

          h_mu_jet_dr->fill(drmin);
          if (jets[0].pT() > 650*GeV)  h_mu_jet_dr_pt650->fill(drmin);
          else if (jets[0].pT() > 500*GeV && jets[0].pT() < 600*GeV) {
            h_mu_jet_dr_pt500600->fill(drmin);
          }
        }

        /// Normalise histograms etc., after the run
        void finalize() {
          const double sf = crossSection() / femtobarn / sumOfWeights();
          scale(h_mu_jet_dr, sf);
          scale(h_mu_jet_dr_pt500600, sf);
          scale(h_mu_jet_dr_pt650, sf);
        }

        //@}

    protected:

        size_t _mode;

    private:

        /// @name Histograms
        //@{
        Histo1DPtr h_mu_jet_dr;
        Histo1DPtr h_mu_jet_dr_pt500600;
        Histo1DPtr h_mu_jet_dr_pt650;
        //@}
  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2016_I1487726);
}

