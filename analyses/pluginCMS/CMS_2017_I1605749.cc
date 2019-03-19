// -*- C++ -*-
// Rivet framework
#include "Rivet/Analysis.hh"

// Projections
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  using namespace Cuts;

  class CMS_2017_I1605749 : public Analysis {
  public:

    // Constructor
    CMS_2017_I1605749()
      : Analysis("CMS_2017_I1605749")
    {    }

    // Book histograms and initialise projections before the run
    void init() {
      // Projections
      const FinalState fs(-5.0, 5.0, 0.0*GeV);
      addProjection(FastJets(fs, FastJets::ANTIKT, 0.5), "Jets");

      // Jet Charge Histos
      for (int i = 1; i <= 18; i++) {
        _h_Charge[i - 1] = bookHisto1D(i, 1, 1);
      }
    }

    // Perform the per-event analysis
    void analyze(const Event& event) {
      const Jets& jets = applyProjection<FastJets>(event, "Jets").jetsByPt(10.0*GeV);

      if (jets.size() < 2) vetoEvent;

      double leadingpt = jets[0].pt()/GeV;
      double subleadingpt = jets[1].pt()/GeV;

      if (jets.size() < 2 ||
          jets[0].abseta() >= 1.5 ||
          jets[1].abseta() >= 1.5 ||
          leadingpt < 400.0 || subleadingpt < 100.0) {
        vetoEvent;
      }

      const double weight = event.weight();

      vector<Particle> constituents1 = jets[0].constituents();
      std::vector<double> numerator(9, 0), denominator(9, 0);

      double t_jetcharge1, t_jetcharge1k6, t_jetcharge1k3;
      double t_jetchargeL1, t_jetchargeL1k6, t_jetchargeL1k3;
      double t_jetchargeT1, t_jetchargeT1k6, t_jetchargeT1k3;

      denominator[0] = leadingpt;
      denominator[1] = std::pow(leadingpt, 0.6);
      denominator[2] = std::pow(leadingpt, 0.3);

      if (constituents1.size() > 0) {
        for (unsigned j = 0; j < constituents1.size(); j++) {
          if (std::abs(constituents1[j].pdgId()) > 9 &&
              std::abs(constituents1[j].pdgId())!= 21) {
            if (constituents1[j].pt() > 1*GeV) {
              double charge = constituents1[j].charge();
              double mom = constituents1[j].pt();
              double dotproduct = constituents1[j].p3().dot(jets[0].p3()) / jets[0].p();
              double crossproduct = constituents1[j].p3().cross(jets[0].p3()).mod() / jets[0].p();

              numerator[0] += (mom * charge);
              numerator[1] += ((std::pow(mom, 0.6)) * charge);
              numerator[2] += ((std::pow(mom, 0.3)) * charge);

              numerator[3] += (dotproduct * charge);
              numerator[4] += ((std::pow(dotproduct, 0.6)) * charge);
              numerator[5] += ((std::pow(dotproduct, 0.3)) * charge);

              denominator[3] += dotproduct;
              denominator[4] += (std::pow(dotproduct, 0.6));
              denominator[5] += (std::pow(dotproduct, 0.3));

              numerator[6] += (crossproduct * charge);
              numerator[7] += ((std::pow(crossproduct, 0.6)) * charge);
              numerator[8] += ((std::pow(crossproduct, 0.3)) * charge);

              denominator[6] += crossproduct;
              denominator[7] += (std::pow(crossproduct, 0.6));
              denominator[8] += (std::pow(crossproduct, 0.3));
            }
          }
        }
      }

      t_jetcharge1    = (denominator[0] > 0) ? numerator[0] / denominator[0] : 0;
      t_jetcharge1k6  = (denominator[1] > 0) ? numerator[1] / denominator[1] : 0;
      t_jetcharge1k3  = (denominator[2] > 0) ? numerator[2] / denominator[2] : 0;
      t_jetchargeL1   = (denominator[3] > 0) ? numerator[3] / denominator[3] : 0;
      t_jetchargeL1k6 = (denominator[4] > 0) ? numerator[4] / denominator[4] : 0;
      t_jetchargeL1k3 = (denominator[5] > 0) ? numerator[5] / denominator[5] : 0;
      t_jetchargeT1   = (denominator[6] > 0) ? numerator[6] / denominator[6] : 0;
      t_jetchargeT1k6 = (denominator[7] > 0) ? numerator[7] / denominator[7] : 0;
      t_jetchargeT1k3 = (denominator[8] > 0) ? numerator[8] / denominator[8] : 0;

      _h_Charge[0]->fill(t_jetcharge1, weight);
      _h_Charge[1]->fill(t_jetcharge1k6, weight);
      _h_Charge[2]->fill(t_jetcharge1k3, weight);
      _h_Charge[3]->fill(t_jetchargeL1, weight);
      _h_Charge[4]->fill(t_jetchargeL1k6, weight);
      _h_Charge[5]->fill(t_jetchargeL1k3, weight);
      _h_Charge[6]->fill(t_jetchargeT1, weight);
      _h_Charge[7]->fill(t_jetchargeT1k6, weight);
      _h_Charge[8]->fill(t_jetchargeT1k3, weight);

      if (leadingpt > 400 && leadingpt < 700) {
        _h_Charge[9]->fill(t_jetcharge1k6, weight);
        _h_Charge[12]->fill(t_jetchargeL1k6, weight);
        _h_Charge[15]->fill(t_jetchargeT1k6, weight);
      } else if (leadingpt > 700 && leadingpt < 1000) {
        _h_Charge[10]->fill(t_jetcharge1k6, weight);
        _h_Charge[13]->fill(t_jetchargeL1k6, weight);
        _h_Charge[16]->fill(t_jetchargeT1k6, weight);
      } else if (leadingpt > 1000 && leadingpt < 1800) {
        _h_Charge[11]->fill(t_jetcharge1k6, weight);
        _h_Charge[14]->fill(t_jetchargeL1k6, weight);
        _h_Charge[17]->fill(t_jetchargeT1k6, weight);
      }
    }

    // Normalise histograms etc., after the run
    void finalize() {
      for (int j = 0; j < 18; j++) {
        normalize(_h_Charge[j]);
        for (size_t i = 0; i <  _h_Charge[j]-> numBins(); i++) {
          _h_Charge[j]->bin(i).scaleW(1.0 / _h_Charge[j]->bin(i).width());
        }
      }
    }

  private:
    Histo1DPtr _h_Charge[18];
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2017_I1605749);
}
