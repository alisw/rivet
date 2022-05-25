// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief D0 differential Z/\f$ \gamma^* \f$ + jet + \f$ X \f$ cross sections
  ///
  /// @author Gavin Hesketh, Andy Buckley, Frank Siegert
  class D0_2008_S7863608 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(D0_2008_S7863608);


    /// @name Analysis methods
    /// @{

    /// Book histograms
    void init() {
      /// @todo These clustering arguments look odd: are they ok?
      Cut cut = Cuts::abseta < 1.7 && Cuts::pT > 15*GeV;
      ZFinder zfinder(FinalState(), cut, PID::MUON, 65*GeV, 115*GeV, 0.2, ZFinder::ClusterPhotons::NONE, ZFinder::AddPhotons::YES);
      declare(zfinder, "ZFinder");

      FastJets conefinder(zfinder.remainingFinalState(), FastJets::D0ILCONE, 0.5);
      declare(conefinder, "ConeFinder");

      book(_sum_of_weights_inclusive, "sum_of_weights_inclusive");

      book(_h_jet_pT_cross_section ,1, 1, 1);
      book(_h_jet_pT_normalised ,1, 1, 2);
      book(_h_jet_y_cross_section ,2, 1, 1);
      book(_h_jet_y_normalised ,2, 1, 2);
      book(_h_Z_pT_cross_section ,3, 1, 1);
      book(_h_Z_pT_normalised ,3, 1, 2);
      book(_h_Z_y_cross_section ,4, 1, 1);
      book(_h_Z_y_normalised ,4, 1, 2);
      book(_h_total_cross_section ,5, 1, 1);
    }


    // Do the analysis
    void analyze(const Event& e) {
      const ZFinder& zfinder = apply<ZFinder>(e, "ZFinder");
      if (zfinder.bosons().size()==1) {
        _sum_of_weights_inclusive->fill();
        const JetAlg& jetpro = apply<JetAlg>(e, "ConeFinder");
        const Jets& jets = jetpro.jetsByPt(20*GeV);
        Jets jets_cut;
        for (const Jet& j : jets) {
          if (j.abseta() < 2.8) {
            jets_cut.push_back(j);
          }
        }

        // Return if there are no jets:
        if(jets_cut.size()<1) {
          MSG_DEBUG("Skipping event " << numEvents() << " because no jets pass cuts ");
          vetoEvent;
        }

        const FourMomentum Zmom = zfinder.bosons()[0].momentum();

        // In jet pT
        _h_jet_pT_cross_section->fill( jets_cut[0].pT());
        _h_jet_pT_normalised->fill( jets_cut[0].pT());
        _h_jet_y_cross_section->fill( fabs(jets_cut[0].rapidity()));
        _h_jet_y_normalised->fill( fabs(jets_cut[0].rapidity()));

        // In Z pT
        _h_Z_pT_cross_section->fill(Zmom.pT());
        _h_Z_pT_normalised->fill(Zmom.pT());
        _h_Z_y_cross_section->fill(Zmom.absrap());
        _h_Z_y_normalised->fill(Zmom.absrap());

        _h_total_cross_section->fill(1960);
      }
    }


    /// Finalize
    void finalize() {
      const double invlumi = crossSection()/sumOfWeights();
      scale(_h_total_cross_section, invlumi);
      scale(_h_jet_pT_cross_section, invlumi);
      scale(_h_jet_y_cross_section, invlumi);
      scale(_h_Z_pT_cross_section, invlumi);
      scale(_h_Z_y_cross_section, invlumi);

      double factor=1/dbl(*_sum_of_weights_inclusive);
      if (_sum_of_weights_inclusive->val() == 0) factor = 0;
      scale(_h_jet_pT_normalised, factor);
      scale(_h_jet_y_normalised, factor);
      scale(_h_Z_pT_normalised, factor);
      scale(_h_Z_y_normalised, factor);
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _h_jet_pT_cross_section;
    Histo1DPtr _h_jet_y_cross_section;
    Histo1DPtr _h_Z_pT_cross_section;
    Histo1DPtr _h_Z_y_cross_section;
    Histo1DPtr _h_total_cross_section;
    Histo1DPtr _h_jet_pT_normalised;
    Histo1DPtr _h_jet_y_normalised;
    Histo1DPtr _h_Z_pT_normalised;
    Histo1DPtr _h_Z_y_normalised;
    /// @}

    CounterPtr _sum_of_weights_inclusive;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(D0_2008_S7863608, D0_2008_I792812);

}
