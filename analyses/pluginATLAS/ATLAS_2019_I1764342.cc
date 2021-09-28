// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"

namespace Rivet {

  /// @brief Z(ll)y cross-section at 13 TeV
  class ATLAS_2019_I1764342 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2019_I1764342);

    /// @name Analysis methods
    //@{
    /// Book histograms and initialise projections before the run
    void init() {
      // Prompt photons
      const PromptFinalState photon_fs(Cuts::abspid == PID::PHOTON && Cuts::pT > 30*GeV && Cuts::abseta < 2.37);
      declare(photon_fs, "Photons");

      // Prompt leptons
      const PromptFinalState bareelectron_fs = Cuts::abspid == PID::ELECTRON;
      const PromptFinalState baremuon_fs = Cuts::abspid == PID::MUON;

      // Dressed leptons
      const FinalState allphoton_fs(Cuts::abspid == PID::PHOTON); // photons used for lepton dressing
      const Cut leptoncut = Cuts::pT > 25*GeV && Cuts::abseta < 2.47;
      const DressedLeptons dressedelectron_fs(allphoton_fs, bareelectron_fs, 0.1, leptoncut, true); // use *all* photons for lepton dressing
      const DressedLeptons dressedmuon_fs(allphoton_fs, baremuon_fs, 0.1, leptoncut, true); // use *all* photons for lepton dressing

      declare(dressedelectron_fs, "Electrons");
      declare(dressedmuon_fs, "Muons");
      
      // FS excluding the leading photon
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(photon_fs);
      vfs.addVetoOnThisFinalState(dressedmuon_fs);
      vfs.addVetoOnThisFinalState(InvisibleFinalState());
      declare(vfs, "isolatedFS");

      // Histograms
      book(_hist_EgammaT,     2, 1, 1); // dSigma / dE^gamma_T 
      book(_hist_etagamma,    3, 1, 1);
      book(_hist_mZgamma,     4, 1, 1); // dSigma / dm^{Zgamma}
      book(_hist_EZgammaT,    5, 1, 1);
      book(_hist_dPhiZgamma,  6, 1, 1);
      book(_hist_ETbyMZgamma, 7, 1, 1);
    }


/// Perform the per-event analysis
   void analyze(const Event& event) {
     // Get objects
     vector<DressedLepton> electrons = apply<DressedLeptons>(event, "Electrons").dressedLeptons();
     vector<DressedLepton> muons = apply<DressedLeptons>(event, "Muons").dressedLeptons();
     const Particles& photons = apply<PromptFinalState>(event, "Photons").particlesByPt();

     if (photons.empty())  vetoEvent;
     if (electrons.size() < 2 && muons.size() < 2)  vetoEvent;
     vector<DressedLepton> lep;
     // Sort the dressed leptons by pt
     if (electrons.size() >= 2) {
       lep.push_back(electrons[0]);
       lep.push_back(electrons[1]);
     } else {
       lep.push_back(muons[0]);
       lep.push_back(muons[1]);
     }
     if(lep[0].Et() < 30)  vetoEvent;
     double mll = (lep[0].momentum() + lep[1].momentum()).mass();
     if (mll < 40*GeV) vetoEvent;

     vector<Particle> selectedPh;
     Particles fs = apply<VetoedFinalState>(event, "isolatedFS").particles();
     for (const Particle& ph : photons){
       // check photon isolation
       double coneEnergy(0.0);
       for (const Particle& p : fs) {
         if ( deltaR(ph, p) < 0.2 )  coneEnergy += p.Et();
       }
       if (coneEnergy / ph.Et() > 0.07 )  continue;
       if (deltaR(ph, lep[0]) < 0.4) continue;
       if (deltaR(ph, lep[1]) < 0.4) continue;
       selectedPh.push_back(ph);
     }

     if(selectedPh.size()<1) vetoEvent;
     double mlly = (lep[0].momentum() + lep[1].momentum() + selectedPh[0].momentum()).mass();
     if(mll + mlly <= 182*GeV) vetoEvent;

     double ptlly = (lep[0].momentum() + lep[1].momentum() + selectedPh[0].momentum()).pT();
     double dphilly = deltaPhi((lep[0].momentum() + lep[1].momentum()).phi(), selectedPh[0].momentum().phi());

     // Fill plots
     _hist_EgammaT->fill(selectedPh[0].pT()/GeV);
     _hist_etagamma->fill(selectedPh[0].abseta());
     _hist_mZgamma->fill(mlly/GeV);
     _hist_EZgammaT->fill(ptlly/GeV);
     _hist_dPhiZgamma->fill(dphilly/pi);
     _hist_ETbyMZgamma->fill(ptlly/mlly);
   } // end of analysis

   /// Normalise histograms etc., after the run
   void finalize() {
      const double sf = crossSection()/femtobarn/sumOfWeights();
      scale(_hist_EgammaT, sf);
      scale(_hist_etagamma, sf);
      scale(_hist_mZgamma, sf);
      scale(_hist_EZgammaT, sf);
      scale(_hist_dPhiZgamma, sf/pi);
      scale(_hist_ETbyMZgamma, sf);
   }


  protected:

    // Data members like post-cuts event weight counters go here
    size_t _mode;

  private:

    /// Histograms
    Histo1DPtr _hist_EgammaT;
    Histo1DPtr _hist_etagamma;
    Histo1DPtr _hist_mZgamma;
    Histo1DPtr _hist_EZgammaT;
    Histo1DPtr _hist_dPhiZgamma;
    Histo1DPtr _hist_ETbyMZgamma;
  }; 

   // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2019_I1764342);

}
