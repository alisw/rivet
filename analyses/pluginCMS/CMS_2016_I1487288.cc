#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/WFinder.hh"

namespace Rivet {


  /// @brief WZ production cross section in pp collisions at 7 and 8 TeV
  class CMS_2016_I1487288 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2016_I1487288);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs(Cuts::abseta < 4.9);

      FastJets fj(fs, FastJets::ANTIKT, 0.5, JetAlg::Muons::ALL, JetAlg::Invisibles::DECAY);
      declare(fj, "Jets");

      ZFinder zeeFinder(fs, Cuts::abseta < 2.5 && Cuts::pT > 20*GeV, PID::ELECTRON, 71*GeV, 111*GeV);
      declare(zeeFinder, "Zee");

      ZFinder zmumuFinder(fs, Cuts::abseta < 2.4 && Cuts::pT > 20*GeV, PID::MUON, 71*GeV, 111*GeV);
      declare(zmumuFinder, "Zmumu");

      WFinder weFinder(fs, Cuts::abseta < 2.5 && Cuts::pT > 20*GeV, PID::ELECTRON, 60*GeV, 100*GeV, 30*GeV);
      declare(weFinder, "We");

      WFinder wmuFinder(fs, Cuts::abseta < 2.4 && Cuts::pT > 20*GeV, PID::MUON, 60*GeV, 100*GeV, 30*GeV);
      declare(wmuFinder, "Wmu");

      book(_h_ZpT, "d03-x01-y01");
      book(_h_Njet, "d04-x01-y01", {-0.5, 0.5, 1.5, 2.5, 3.5}); ///< @todo Ref data has null bin widths
      book(_h_JpT, "d05-x01-y01");

      MSG_WARNING("\033[91;1mLIMITED VALIDITY - check info file for details!\033[m");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Find Z -> l+ l-
      const ZFinder& zeeFS = apply<ZFinder>(event, "Zee");
      const ZFinder& zmumuFS = apply<ZFinder>(event, "Zmumu");
      const Particles zlls = zeeFS.bosons() + zmumuFS.bosons();
      if (zlls.empty()) vetoEvent;

      // Next find the W
      const WFinder& weFS = apply<WFinder>(event, "We");
      const WFinder& wmuFS = apply<WFinder>(event, "Wmu");
      const Particles wls = weFS.bosons() + wmuFS.bosons();
      if (wls.empty()) vetoEvent;


      // If more than one Z candidate, use the one with Mll nearest to MZ
      const Particles zlls_mz = sortBy(zlls, [](const Particle& a, const Particle& b){
          return fabs(a.mass() - 91.2*GeV) < fabs(b.mass() - 91.2*GeV); });
      const Particle& Z = zlls_mz.front();
      // const bool isZee = any(Z.constituents(), hasAbsPID(PID::ELECTRON));

      // If more than one Z candidate, use the one with Mll nearest to MZ
      const Particles wls_mw = sortBy(wls, [](const Particle& a, const Particle& b){
          return fabs(a.mass() - 80.4*GeV) < fabs(b.mass() - 80.4*GeV); });
      const Particle& W = wls_mw.front();
      // const bool isWe = any(W.constituents(), hasAbsPID(PID::ELECTRON));

      // Isolate W and Z charged leptons from each other
      for (const Particle& lw : W.constituents()) {
        if (lw.charge3() == 0) continue;
        for (const Particle& lz : Z.constituents()) {
          if (deltaR(lw, lz) < 0.1) vetoEvent;
        }
      }

      // Fill Z pT histogram
      _h_ZpT->fill(Z.pT()/GeV);


      // Isolate jets from W and Z charged leptons
      const Particles wzleps = filter_select(W.constituents()+Z.constituents(), isChLepton);
      const Jets& jets = apply<FastJets>("Jets", event).jetsByPt(Cuts::pT > 30*GeV and Cuts::abseta < 2.5);
      const Jets isojets = discardIfAnyDeltaRLess(jets, wzleps, 0.5);

      // Fill jet histograms
      _h_Njet->fill(isojets.size());
      if (!isojets.empty()) _h_JpT->fill(isojets[0].pT()/GeV);

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // Total cross-section is corrected for BR(W->lnu), BR(Z->ll), leptonic-tau fraction f_tau = 6.5-7.6%,
      // and unpublished detector/acceptance signal efficiencies epsilon_sig. Fix to published value: valid for shape comparison only
      const double xsec8tev = 24.09; // picobarn;
      normalize(_h_ZpT,  xsec8tev);
      normalize(_h_Njet, xsec8tev);
      normalize(_h_JpT,  xsec8tev);
    }


  private:

    /// Histogram
    Histo1DPtr _h_ZpT, _h_Njet, _h_JpT;

  };



  DECLARE_RIVET_PLUGIN(CMS_2016_I1487288);

}
