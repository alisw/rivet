// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// D0 Z+jets angular distributions
  class D0_2009_S8349509 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(D0_2009_S8349509);


    /// @name Analysis methods
    /// @{

    /// Book histograms
    void init() {
      Cut cut = Cuts::abseta < 1.7 && Cuts::pT > 15*GeV;
      ZFinder zfinder(FinalState(), cut, PID::MUON, 65*GeV, 115*GeV, 0.2, ZFinder::ClusterPhotons::NONE, ZFinder::AddPhotons::YES);
      declare(zfinder, "ZFinder");

      FastJets conefinder(zfinder.remainingFinalState(), FastJets::D0ILCONE, 0.5);
      declare(conefinder, "ConeFinder");

      book(_h_dphi_jet_Z25 ,1, 1, 1);
      book(_h_dphi_jet_Z45 ,2, 1, 1);

      book(_h_dy_jet_Z25 ,3, 1, 1);
      book(_h_dy_jet_Z45 ,4, 1, 1);

      book(_h_yboost_jet_Z25 ,5, 1, 1);
      book(_h_yboost_jet_Z45 ,6, 1, 1);

      book(_h_dphi_jet_Z25_xs ,1, 1, 2);
      book(_h_dphi_jet_Z45_xs ,2, 1, 2);

      book(_h_dy_jet_Z25_xs ,3, 1, 2);
      book(_h_dy_jet_Z45_xs ,4, 1, 2);

      book(_h_yboost_jet_Z25_xs ,5, 1, 2);
      book(_h_yboost_jet_Z45_xs ,6, 1, 2);

      book(_inclusive_Z_sumofweights, "_inclusive_Z_sumofweights");
    }


    void analyze(const Event& event) {
      const ZFinder& zfinder = apply<ZFinder>(event, "ZFinder");
      if (zfinder.bosons().size() == 1) {
        // count inclusive sum of weights for histogram normalisation
        _inclusive_Z_sumofweights->fill();

        const FourMomentum& zmom = zfinder.bosons()[0].momentum();
        if (zmom.pT() < 25*GeV) vetoEvent;

        Jets jets;
        for (const Jet& j : apply<JetAlg>(event, "ConeFinder").jetsByPt(20*GeV)) {
          if (j.abseta() < 2.8) {
            jets.push_back(j);
            break;
          }
        }

        // Return if there are no jets:
        if (jets.size() < 1) {
          MSG_DEBUG("Skipping event " << numEvents() << " because no jets pass cuts ");
          vetoEvent;
        }

        const FourMomentum& jetmom = jets[0].momentum();
        const double yZ = zmom.rapidity();
        const double yjet = jetmom.rapidity();
        const double dphi = deltaPhi(zmom, jetmom);
        const double dy = deltaRap(zmom, jetmom);
        const double yboost = fabs(yZ+yjet)/2;

        if (zmom.pT() > 25*GeV) {
          _h_dphi_jet_Z25->fill(dphi);
          _h_dy_jet_Z25->fill(dy);
          _h_yboost_jet_Z25->fill(yboost);
          _h_dphi_jet_Z25_xs->fill(dphi);
          _h_dy_jet_Z25_xs->fill(dy);
          _h_yboost_jet_Z25_xs->fill(yboost);
        }
        if (zmom.pT() > 45*GeV) {
          _h_dphi_jet_Z45->fill(dphi);
          _h_dy_jet_Z45->fill(dy);
          _h_yboost_jet_Z45->fill(yboost);
          _h_dphi_jet_Z45_xs->fill(dphi);
          _h_dy_jet_Z45_xs->fill(dy);
          _h_yboost_jet_Z45_xs->fill(yboost);
        }
      }

    }


    void finalize() {
      if (_inclusive_Z_sumofweights->val() == 0) return;
      scale(_h_dphi_jet_Z25, 1/ *_inclusive_Z_sumofweights);
      scale(_h_dphi_jet_Z45, 1/ *_inclusive_Z_sumofweights);
      scale(_h_dy_jet_Z25, 1/ *_inclusive_Z_sumofweights);
      scale(_h_dy_jet_Z45, 1/ *_inclusive_Z_sumofweights);
      scale(_h_yboost_jet_Z25, 1/ *_inclusive_Z_sumofweights);
      scale(_h_yboost_jet_Z45, 1/ *_inclusive_Z_sumofweights);

      scale(_h_dphi_jet_Z25_xs, crossSectionPerEvent());
      scale(_h_dphi_jet_Z45_xs, crossSectionPerEvent());
      scale(_h_dy_jet_Z25_xs, crossSectionPerEvent());
      scale(_h_dy_jet_Z45_xs, crossSectionPerEvent());
      scale(_h_yboost_jet_Z25_xs, crossSectionPerEvent());
      scale(_h_yboost_jet_Z45_xs, crossSectionPerEvent());
    }

    /// @}


  private:

    /// @name Histograms (normalised)
    /// @{
    Histo1DPtr _h_dphi_jet_Z25;
    Histo1DPtr _h_dphi_jet_Z45;

    Histo1DPtr _h_dy_jet_Z25;
    Histo1DPtr _h_dy_jet_Z45;

    Histo1DPtr _h_yboost_jet_Z25;
    Histo1DPtr _h_yboost_jet_Z45;
    /// @}

    /// @name Histograms (absolute cross sections)
    /// @{
    Histo1DPtr _h_dphi_jet_Z25_xs;
    Histo1DPtr _h_dphi_jet_Z45_xs;

    Histo1DPtr _h_dy_jet_Z25_xs;
    Histo1DPtr _h_dy_jet_Z45_xs;

    Histo1DPtr _h_yboost_jet_Z25_xs;
    Histo1DPtr _h_yboost_jet_Z45_xs;
    /// @}

    CounterPtr _inclusive_Z_sumofweights;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(D0_2009_S8349509, D0_2009_I826756);

}
