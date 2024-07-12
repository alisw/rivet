// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/DISFinalState.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  class MC_DIS_CALIB : public SingleValueProjection {
  public:
    MC_DIS_CALIB() {
      declare(DISFinalState(DISFinalState::BoostFrame::HCM),"DISFS");
    }
    virtual ~MC_DIS_CALIB() {}

    DEFAULT_RIVET_PROJ_CLONE(MC_DIS_CALIB);
  
  protected:
    void project(const Event& e) {
      clear();
      int nch = 0;
      for (auto p : apply<DISFinalState>(e, "DISFS").particles())
        if (p.isHadron() && p.isCharged() && abs(p.eta()) > 3.) ++nch;
      set(nch);
    }
    virtual CmpState compare(const Projection& p) const {
      return mkNamedPCmp(p, "DISFS");
    }
  };

  class MC_DIS_PERC : public Analysis {
  public:
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_DIS_PERC);

    void init() {
      declare(MC_DIS_CALIB(),"PERC");
      book(_fm, "FM",25,0,25);
    }
    void analyze(const Event& event) {
      _fm->fill(apply<MC_DIS_CALIB>(event, "PERC")());
    }
    void finalize() {
      _fm->normalize();
    }
    Histo1DPtr _fm;
  
  };

  DECLARE_RIVET_PLUGIN(MC_DIS_PERC);

  class MC_DIS_Mod : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_DIS_Mod);
    
    /// Initialize histograms, projections etc.
    void init() {
      DISLepton lepton(options());
      declare(lepton, "Lepton");
      declare(DISKinematics(lepton), "Kinematics");
      /*string boostSystem = getOption<string>("boost","hcm");
      if (boostSystem == "breit")
        declare(DISFinalState(DISFinalState::BoostFrame::BREIT),"DISFS");
      else
        declare(DISFinalState(DISFinalState::BoostFrame::HCM),"DISFS");
      */

      const DISFinalState& disfs = declare(DISFinalState(DISFinalState::BoostFrame::HCM),"DISFS");
      declare(FastJets(disfs, fastjet::JetAlgorithm::kt_algorithm, fastjet::RecombinationScheme::pt_scheme, 1.0,
        JetAlg::Muons::ALL, JetAlg::Invisibles::NONE, nullptr), "JETS");
      declareCentrality(MC_DIS_CALIB(), "MC_DIS_PERC", "FM", "PERC");
      percentileBins = {20, 40, 60, 100};

      // Book histograms
      book(_hist_Q2, "Q2",logspace(100,0.1, 1000.0));
      book(_hist_y, "y",100,0.,1.);
      book(_hist_x, "xBj",logspace(100,0.00001, 1.0));
      book(_hist_W2, "W2", logspace(100, 0.1, 10000.));
      book(_hist_multeta1, "multeta1", 30, -6, 6);      
      book(_hist_multeta2, "multeta2", 30, -6, 6);
      book(sow1, "_sow1");      
      book(sow2, "_sow2");
      book(_hist_njets, "njets", 5, 0, 5);
      book(_hist_jetET, "jetET", 20, 4, 60);      

      for (size_t i = 0; i < percentileBins.size(); ++i) {
	 double pb = percentileBins[i];
	 book(histEtaPerc[pb], "histEtaCent-"+toString(pb), 30, -6, 6);
	 book(sowPerc[pb], "_sow"+toString(pb));
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");
      const DISFinalState& disfs = apply<DISFinalState>(event, "DISFS");
      if ( dk.failed() ) return;
      double x  = dk.x();
      double y  = dk.y();
      double Q2 = dk.Q2();
      double W2 = dk.W2();
      _hist_Q2->fill(Q2);
      _hist_y->fill(y);
      _hist_x->fill(x);
      _hist_W2->fill(W2);

      const CentralityProjection& cent = apply<CentralityProjection>(event, "PERC");
      double c = cent();
      auto hItr = histEtaPerc.upper_bound(c);
      if (hItr == histEtaPerc.end()) return;
      auto sItr = sowPerc.upper_bound(c);
      if (sItr == sowPerc.end()) return;
      sItr->second->fill();
      if (x < 1e-3) sow1->fill();
      else sow2->fill();
      for (auto p : disfs.particles()) {
	if (p.isHadron() && p.isCharged()) {
	  hItr->second->fill(p.eta());
          if (x < 1e-3) _hist_multeta1->fill(p.eta());
	  else _hist_multeta2->fill(p.eta());
	}
      }
      Jets jets = apply<FastJets>(event, "JETS").jets(Cuts::Et > 4. * GeV, cmpMomByEt);
      _hist_njets->fill(jets.size());
      if (jets.size() == 0) return;
      _hist_jetET->fill(jets[0].Et());
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_hist_Q2); // normalize to unity
      normalize(_hist_y); // normalize to unity
      scale(_hist_multeta1, 1./sow1->sumW());
      scale(_hist_multeta2, 1./sow2->sumW());

      for (size_t i = 0; i < percentileBins.size(); ++i) {
	 double pb = percentileBins[i];
	 histEtaPerc[pb]->scaleW(1./sowPerc[pb]->sumW());
      }
      normalize(_hist_njets);
      scale(_hist_jetET, crossSection()/picobarn/sumOfWeights());
    }

    /// The histograms.
    Histo1DPtr _hist_Q2, _hist_y, _hist_x, _hist_W2;
    Histo1DPtr _hist_multeta1, _hist_multeta2;
    CounterPtr sow1, sow2;
    map<double, Histo1DPtr> histEtaPerc;
    map<double, CounterPtr> sowPerc;
    vector<double> percentileBins;
    Histo1DPtr _hist_njets, _hist_jetET;
  };

  DECLARE_RIVET_PLUGIN(MC_DIS_Mod);
}
