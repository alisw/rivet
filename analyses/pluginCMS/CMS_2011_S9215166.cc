// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


    class CMS_2011_S9215166 : public Analysis {
        public:

            /// Constructor
            CMS_2011_S9215166() : Analysis("CMS_2011_S9215166") {  }


            void init() {
                const FinalState fs((Cuts::etaIn(-6.0, 6.0)));
                declare(fs, "FS");
                declare(FastJets(fs, FastJets::ANTIKT, 0.5), "Jets");

                VetoedFinalState fsv(fs);
                fsv.vetoNeutrinos();
                fsv.addVetoPairDetail(PID::MUON, 0.0*GeV, 99999.9*GeV);
                declare(fsv, "fsv");

                // For the MB ND selection
                const ChargedFinalState fschrgd((Cuts::etaIn(-6.0,6.0)));
                declare(fschrgd, "fschrgd");
                VetoedFinalState fschrgdv(fschrgd);
                fschrgdv.vetoNeutrinos();
                declare(fschrgdv, "fschrgdv");

                if (fuzzyEquals(sqrtS()/GeV, 900, 1E-3)) {
                    book(_hist_mb      ,1, 1, 1); // energy flow in MB, 0.9 TeV
                    book(_hist_dijet ,2, 1, 1); // energy flow in dijet events, 0.9 TeV
                } else if (fuzzyEquals(sqrtS()/GeV, 7000, 1E-3)) {
                    book(_hist_mb      ,3, 1, 1); // energy flow in MB, 7 TeV
                    book(_hist_dijet ,4, 1, 1); // energy flow in dijet events, 7 TeV
                }

                book(_weightMB, "/tmp/weightMB");
                book(_weightDiJet, "/tmp/weightDijet");
            }


            void analyze(const Event& event) {
                // Skip if the event is empty
                const FinalState& fsv = apply<FinalState>(event, "fsv");
                if (fsv.empty()) vetoEvent;

                // Veto diffractive topologies according to defined hadron level
                double count_chrg_forward = 0;
                double count_chrg_backward = 0;
                const FinalState& fschrgdv = apply<FinalState>(event, "fschrgdv");
                for (const Particle& p : fschrgdv.particles()) {
                    if (3.9 < p.eta() && p.eta() < 4.4) count_chrg_forward++;
                    if (-4.4 < p.eta() && p.eta() < -3.9) count_chrg_backward++;
                }
                if (count_chrg_forward == 0 || count_chrg_backward == 0) vetoEvent;
                /// @todo "Diffractive" veto should really also veto dijet events?


                // MINIMUM BIAS EVENTS
                _weightMB->fill();
                for (const Particle& p: fsv.particles()) {
                    _hist_mb->fill(p.abseta(), p.E()/GeV);
                }


                // DIJET EVENTS
                double PTCUT = -1.0;
                if (fuzzyEquals(sqrtS()/GeV, 900, 1E-3)) PTCUT = 8.0*GeV;
                else if (fuzzyEquals(sqrtS()/GeV, 7000, 1E-3)) PTCUT = 20.0*GeV;
                const FastJets& jetpro = apply<FastJets>(event, "Jets");
                const Jets jets = jetpro.jetsByPt(PTCUT);
                if (jets.size() >= 2) {
                    // eta cut for the central jets
                    if (fabs(jets[0].eta()) < 2.5 && fabs(jets[1].eta()) < 2.5) {
                        // Back to back condition of the jets
                        const double diffphi = deltaPhi(jets[1].phi(), jets[0].phi());
                        if (diffphi-PI < 1.0) {
                            _weightDiJet->fill();
                            for (const Particle& p: fsv.particles()) {
                                _hist_dijet->fill(p.abseta(), p.E()/GeV);
                            }
                        }
                    }
                }

            }


            void finalize() {
                scale(_hist_mb   , 0.5/_weightMB->sumW());
                scale(_hist_dijet, 0.5/_weightDiJet->sumW());
            }


        private:

            Histo1DPtr _hist_mb, _hist_dijet;
            CounterPtr _weightMB, _weightDiJet;

    };


    // Hook for the plugin system
    DECLARE_RIVET_PLUGIN(CMS_2011_S9215166);

}
