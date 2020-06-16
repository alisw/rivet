// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Thrust.hh"

namespace Rivet {

  /// Rivet analysis class for CMS_2011_S8957746 dataset
  class CMS_2011_S8957746 : public Analysis {
  public:

    /// Constructor
    CMS_2011_S8957746()
      : Analysis("CMS_2011_S8957746") {  }


    /// Initialization, called once before running
    void init() {
      // Projections
      const FastJets jets(FinalState((Cuts::etaIn(-5.0, 5.0))), FastJets::ANTIKT, 0.5);
      declare(jets, "Jets");

      // Book histograms
      book(_hist_T_90  ,1, 1, 1);
      book(_hist_m_90  ,2, 1, 1);
      book(_hist_T_125 ,3, 1, 1);
      book(_hist_m_125 ,4, 1, 1);
      book(_hist_T_200 ,5, 1, 1);
      book(_hist_m_200 ,6, 1, 1);
    }


    void analyze(const Event& event) {
      const double weight = 1.0;
      const Jets& jets = apply<FastJets>(event, "Jets").jetsByPt(30.0*GeV);
      if (jets.size() < 2 ||
          fabs(jets[0].eta()) >= 1.3 ||
          fabs(jets[1].eta()) >= 1.3 ||
          jets[0].pT() < 90*GeV) {
        vetoEvent;
      }
      std::vector<Vector3> momenta;
      for (const Jet& j : jets) {
        if (j.abseta() < 1.3) {
          Vector3 mom = j.p3();
          mom.setZ(0.0);
          momenta.push_back(mom);
        }
      }
      if (momenta.size() == 2) {
        // We need to use a ghost so that Thrust.calc() doesn't return 1.
        momenta.push_back(Vector3(1e-10*MeV, 0., 0.));
      }
      Thrust thrust;
      thrust.calc(momenta);

      // The lowest bin also includes the underflow:
      const double T = max(log(1-thrust.thrust()), -12.0);
      const double M = max(log(thrust.thrustMajor()), -6.0);
      if (jets[0].pT()/GeV > 200) {
        _hist_T_200->fill(T, weight);
        _hist_m_200->fill(M, weight);
      } else if (jets[0].pT()/GeV > 125) {
        _hist_T_125->fill(T, weight);
        _hist_m_125->fill(M, weight);
      } else if (jets[0].pT()/GeV > 90) {
        _hist_T_90->fill(T, weight);
        _hist_m_90->fill(M, weight);
      }
    }


    void finalize() {
      normalize(_hist_T_90);
      normalize(_hist_m_90);
      normalize(_hist_T_125);
      normalize(_hist_m_125);
      normalize(_hist_T_200);
      normalize(_hist_m_200);
    }


  private:

    Histo1DPtr _hist_T_90;
    Histo1DPtr _hist_m_90;
    Histo1DPtr _hist_T_125;
    Histo1DPtr _hist_m_125;
    Histo1DPtr _hist_T_200;
    Histo1DPtr _hist_m_200;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2011_S8957746);

}
