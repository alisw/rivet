// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// Analysis looking at various distributions of final state particles
  class MC_FSPARTICLES : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_FSPARTICLES);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Projections
      FinalState fs(Cuts::abseta < 5 && Cuts::pT > 500*MeV);
      declare(fs, "FS");
      declare(ChargedFinalState(fs), "CFS");

      // Histograms
      /// @todo Choose E/pT ranged based on input energies... can't do anything about kin. cuts, though
      book(_histMult   , "Mult", 100, -0.5, 199.5);
      book(_histMultCh , "MultCh", 100, -0.5, 199.5);

      book(_histPt   , "Pt", 300, 0, 30);
      book(_histPtCh , "PtCh", 300, 0, 30);

      book(_histE   , "E", 100, 0, 200);
      book(_histECh , "ECh", 100, 0, 200);

      book(_histEtaSumEt , "EtaSumEt", 25, 0, 5);

      book(_histEta    , "Eta", 50, -5, 5);
      book(_histEtaCh  , "EtaCh", 50, -5, 5);
      book(_tmphistEtaPlus, "TMP/EtaPlus", 25, 0, 5);
      book(_tmphistEtaMinus, "TMP/EtaMinus", 25, 0, 5);
      book(_tmphistEtaChPlus, "TMP/EtaChPlus", 25, 0, 5);
      book(_tmphistEtaChMinus, "TMP/EtaChMinus", 25, 0, 5);

      book(_histRapidity    , "Rapidity", 50, -5, 5);
      book(_histRapidityCh  , "RapidityCh", 50, -5, 5);
      book(_tmphistRapPlus, "TMP/RapPlus", 25, 0, 5);
      book(_tmphistRapMinus, "TMP/RapMinus", 25, 0, 5);
      book(_tmphistRapChPlus, "TMP/RapChPlus", 25, 0, 5);
      book(_tmphistRapChMinus, "TMP/RapChMinus", 25, 0, 5);

      book(_histPhi    , "Phi", 50, 0, TWOPI);
      book(_histPhiCh  , "PhiCh", 50, 0, TWOPI);

      book(_histEtaPMRatio , "EtaPMRatio");
      book(_histEtaChPMRatio , "EtaChPMRatio");
      book(_histRapidityPMRatio , "RapidityPMRatio");
      book(_histRapidityChPMRatio , "RapidityChPMRatio");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Charged + neutral final state
      const FinalState& fs = apply<FinalState>(event, "FS");
      MSG_DEBUG("Total multiplicity = " << fs.size());
      _histMult->fill(fs.size());
      for (const Particle& p : fs.particles()) {
        _histEta->fill(p.eta());
        _histEtaSumEt->fill(p.abseta(), p.Et());
        (p.eta() > 0 ? _tmphistEtaPlus : _tmphistEtaMinus)->fill(p.abseta());
        //
        _histRapidity->fill(p.rap());
        (p.rap() > 0 ? _tmphistRapPlus : _tmphistRapMinus)->fill(p.absrap());
        //
        _histPt->fill(p.pT()/GeV);
        _histE->fill(p.E()/GeV);
        _histPhi->fill(p.phi());
      }

      // Same for the charged FS particles only
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      MSG_DEBUG("Total charged multiplicity = " << cfs.size());
      _histMultCh->fill(cfs.size());
      for (const Particle& p : cfs.particles()) {
        _histEtaCh->fill(p.eta());
        (p.eta() > 0 ? _tmphistEtaChPlus : _tmphistEtaChMinus)->fill(p.abseta());
        //
        _histRapidityCh->fill(p.rap());
        (p.rap() > 0 ? _tmphistRapChPlus : _tmphistRapChMinus)->fill(p.absrap());
        //
        _histPtCh->fill(p.pT()/GeV);
        _histECh->fill(p.E()/GeV);
        _histPhiCh->fill(p.phi());
      }

    }


    /// Finalize
    void finalize() {
      normalize(_histMult); normalize(_histEta); normalize(_histRapidity); 
      normalize(_histPt); normalize(_histE); normalize(_histPhi);
      normalize(_histMultCh); normalize(_histEtaCh); normalize(_histRapidityCh); 
      normalize(_histPtCh); normalize(_histECh); normalize(_histPhiCh);
      divide(_tmphistEtaPlus, _tmphistEtaMinus, _histEtaPMRatio);
      divide(_tmphistEtaChPlus, _tmphistEtaChMinus, _histEtaChPMRatio);
      divide(_tmphistRapPlus, _tmphistRapMinus, _histRapidityPMRatio);
      divide(_tmphistRapChPlus, _tmphistRapChMinus, _histRapidityChPMRatio);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _histMult, _histEta, _histRapidity, _histPt, _histE, _histPhi;
    Histo1DPtr _histMultCh,  _histEtaCh, _histRapidityCh, _histPtCh, _histECh, _histPhiCh;
    Profile1DPtr _histEtaSumEt;
    Scatter2DPtr _histEtaPMRatio, _histEtaChPMRatio, _histRapidityPMRatio, _histRapidityChPMRatio;
    //@}

    /// @name Temporary histos used to calculate +/- rapidity ratio plots
    //@{
    Histo1DPtr _tmphistEtaPlus, _tmphistEtaMinus, _tmphistEtaChPlus, _tmphistEtaChMinus;
    Histo1DPtr _tmphistRapPlus, _tmphistRapMinus, _tmphistRapChPlus, _tmphistRapChMinus;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_FSPARTICLES);

}
