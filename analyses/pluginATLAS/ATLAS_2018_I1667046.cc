// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/Cutflow.hh"

namespace Rivet {


  /// Search for R-parity-violating SUSY in multi-jet final states at 13 TeV
  class ATLAS_2018_I1667046 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2018_I1667046);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Projections
      const FinalState fs (Cuts::abseta < 4.9);
      declare("SmallRJ", FastJets(fs, FastJets::ANTIKT, 0.4));
      declare("LargeRJ", FastJets(fs, FastJets::ANTIKT, 1.0));

      // Book histograms
      book(_h_sigmaM, "sigmaM",   50, 200, 2000);
      book(_h_modeta, "ModEta12", 42,   0,  4.2);

      // Cutflows
      _flows.addCutflow("CutFlow1",
                        {"NJet >= 4 ", "Delta12 < 1.4", "PJet1 > 400 GeV", "M SumJ > 1.0 ",
                            "NbJet > 0", "M SumJ > 1.0  & NbJet > 0"});
      _flows.addCutflow("CutFlow2",
                        {"NJet >= 4 ", "Delta12 < 1.4", "NJet >= 5 ", "M SumJ > 0.8 ",
                            "NbJet > 0", "M SumJ > 0.8  & NbJet > 0"});
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      _flows.fillinit();

      // Trim large-R jets and apply cuts
      Jets LRJ_old = apply<FastJets>(event, "LargeRJ").jetsByPt(Cuts::abseta < 4.9);
      Jets LRJJ;
      fastjet::Filter trimmer(fastjet::JetDefinition(fastjet::kt_algorithm,0.2), fastjet::SelectorPtFractionMin(0.05));
      for (Jet& j: LRJ_old) LRJJ.push_back(trimmer(j));
      Jets LRJ = ifilter_select(LRJJ, Cuts::abseta < 2.0 && Cuts::pT > 200*GeV);
      if (LRJ.size() < 4) vetoEvent;
      LRJ = sortByPt(LRJ); // sorting for constructing sigmaM

      // Small R jets need to pass some cuts + need to be BTagged.
      const Jets SRJ = apply<FastJets>(event, "SmallRJ").jetsByPt(Cuts::abseta < 2.5 && Cuts::pT > 50*GeV);
      Jets BT_SRJ = filter_select(SRJ, hasBTag(Cuts::pT > 5*GeV));
      const int tagg = (BT_SRJ.size() == 0) ? 1 : 0; //if there are no B-TAGGED jets are present, tagg = 1.

      // Now to find B-MATCHED Large R Jets: extra step, not useful for SR regions!
      const Jets BM_LRJ = selectIfAnyDeltaRLess(LRJ, BT_SRJ, 1.0);


      // Now to build observables SigmaM and deltaEta
      // Add mass of leading four large R jets,
      double sigmaM = 0.0;
      for (const Jet& j : head(LRJ, 4)) sigmaM += j.mass();
      _h_sigmaM->fill(sigmaM);

      // Build deta between two leading large-R jets
      const double delta_eta = fabs(deltaEta(LRJ[0], LRJ[1]));
      _h_modeta->fill(delta_eta);

      // CutFlow1
      if (LRJ.size() >= 4) {
        _flows["CutFlow1"].fill(1);
        if (delta_eta < 1.4) {
          _flows["CutFlow1"].fill(2);
          if (LRJ[0].pT() > 400*GeV) {
            _flows["CutFlow1"].fill(3);
            if (sigmaM > 1000*GeV)
              _flows["CutFlow1"].fill(4);
            if (tagg == 0) {
              _flows["CutFlow1"].fill(5);
              if (sigmaM > 1000*GeV){
                _flows["CutFlow1"].fill(6);
              }
            } //end of btagg loop
          } //end of pT loop
        } //end of delta loop
      } //end of Njet>4 loop

      // CutFlow2
      if (LRJ.size() >= 4) {
        _flows["CutFlow2"].fill(1);
        if (delta_eta < 1.4) {
          _flows["CutFlow2"].fill(2);
          if (LRJ.size() >= 5) {
            _flows["CutFlow2"].fill(3);
            if (sigmaM > 800*GeV)
              _flows["CutFlow2"].fill(4);
            if (tagg == 0) {
              _flows["CutFlow2"].fill(5);
              if (sigmaM > 800*GeV){_flows["CutFlow2"].fill(6);
              }
            } //end of btagg loop
          } //Njet>5 loop
        } //end of delta loop
      } //Njet>4 loop

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double expected = 36.1*crossSection()/femtobarn;
      normalize(_h_sigmaM, expected/sumOfWeights());
      normalize(_h_modeta, expected/sumOfWeights());
      // _flows.scale(99.7/numEvents());
      MSG_INFO(_flows);
    }

    //@}


  private:

    /// @name Histograms
    Histo1DPtr _h_sigmaM, _h_modeta;

    // Cutflows
    Cutflows _flows;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2018_I1667046);


}
