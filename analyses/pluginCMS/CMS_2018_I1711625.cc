// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"

namespace Rivet {


  /// Drell-Yan differential cross section measurement @ 13 TeV
  class CMS_2018_I1711625 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2018_I1711625);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // prompt final state electrons
      const PromptFinalState el_pfs = PromptFinalState(Cuts::abspid == PID::ELECTRON);
      declare(el_pfs, "PromptFinalStateElectrons");

      // prompt final state muons
      const PromptFinalState mu_pfs = PromptFinalState(Cuts::abspid == PID::MUON);
      declare(mu_pfs, "PromptFinalStateMuons");

      // dressed leptons
      const FinalState photon_fs = FinalState(Cuts::abspid == PID::PHOTON);

      const DressedLeptons mu_dressed(photon_fs, mu_pfs, 0.1, Cuts::open());
      declare(mu_dressed, "DressedMuons");

      book(_h_massMuMu, 3, 1, 1); /// muon channel result in full-phase space @ dressed level
      book(_h_massMuMuFiducial, 5, 1, 1); /// muon channel result in fiducial region @ post-FSR level
      book(_h_massEEFiducial, 6, 1, 1); /// electron channel result in fiducial region @ post-FSR level
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const DressedLeptons muons_dressed = apply<DressedLeptons>(event, "DressedMuons");
      bool filled_mu = FillHistogram_DressedLepton(muons_dressed, 13);
      if ( filled_mu ) {

        const PromptFinalState muons_PFS = apply<PromptFinalState>(event, "PromptFinalStateMuons");
        FillHistogram_PFSLepton(muons_PFS, 13);
      }
      else { // electron channel

        const PromptFinalState electrons_PFS = apply<PromptFinalState>(event, "PromptFinalStateElectrons");
        FillHistogram_PFSLepton(electrons_PFS, 11);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_massMuMu, crossSection()/picobarn/sumOfWeights()); /// norm to cross section
      scale(_h_massMuMuFiducial, crossSection()/picobarn/sumOfWeights()); /// norm to cross section
      scale(_h_massEEFiducial, crossSection()/picobarn/sumOfWeights()); /// norm to cross section

    }

    //@}


  private:


    /// @name Histograms
    //@{
    Histo1DPtr _h_massMuMu;
    Histo1DPtr _h_massMuMuFiducial;
    Histo1DPtr _h_massEEFiducial;
    //@}


    // select two opposite sign leptons with highest pT & fill the histogram for full-phase space diff. x-section
    bool FillHistogram_DressedLepton(DressedLeptons leptons_dressed, int leptonID) {
      bool filled = false;

      vector< DressedLepton > vec_dressedLepByPt = leptons_dressed.dressedLeptons();

      int nLepton_dressed = (int)vec_dressedLepByPt.size();
      if ( nLepton_dressed >= 2 ) {
        int index_lepton1 = -1;
        int index_lepton2 = -1;
        FindDressedLeptonPair_HighestPt(vec_dressedLepByPt, index_lepton1, index_lepton2);

        if ( index_lepton1 != -1 && index_lepton2 != -1 ) {
          DressedLepton lepton1_dressed = vec_dressedLepByPt[index_lepton1];
          DressedLepton lepton2_dressed = vec_dressedLepByPt[index_lepton2];

          const FourMomentum pVec_diLep = lepton1_dressed.mom() + lepton2_dressed.mom();
          double mass = pVec_diLep.mass();

          // // fill histograms
          if ( leptonID == 13 ) _h_massMuMu->fill(mass/GeV);

          filled = true;
        }
      } // end of if( nLepton_dressed >= 2 )

      return filled;
    }


    void FindDressedLeptonPair_HighestPt(vector<DressedLepton>& vec_dressedLepByPt, int& index_lepton1, int& index_lepton2) {
      // 1st lepton: lepton with highest pT
      int nLepton_dressed = int(vec_dressedLepByPt.size());
      for (int i=0; i<nLepton_dressed; ++i) {
        auto& lepton = vec_dressedLepByPt[i]; // decreasing order of pT
        if ( lepton.isLepton() ) {
          index_lepton1 = i;
          break;
        }
      }

      // if no lepton is found in the leptons_dressed
      if ( index_lepton1 < 0 ) {
        index_lepton1 = -1;
        index_lepton2 = -1;
        return;
      }

      // 2nd lepton: lepton with highest-pT among the leptons with the opposite sign with 1st lepton
      int pdgID_lepton1 = vec_dressedLepByPt[index_lepton1].pid();
      for (int i=index_lepton1+1; i<nLepton_dressed; ++i) { // starting after lepton1
        auto& lepton = vec_dressedLepByPt[i];
        if ( lepton.isLepton() && lepton.pid() == (-1)*pdgID_lepton1 ) {
          index_lepton2 = i;
          break;
        }
      }
    }


    // select two opposite sign leptons with highest pT & fill the histogram for the fiducial diff. x-section
    void FillHistogram_PFSLepton(PromptFinalState leptons_PFS, int leptonID) {
      vector< Particle > vec_PFSLepByPt = leptons_PFS.particlesByPt();

      int nLepton_PFS = int(vec_PFSLepByPt.size());
      if ( nLepton_PFS >= 2 ) {
        int index_lepton1 = -1;
        int index_lepton2 = -1;
        FindPFSLeptonPair_HighestPtWithinAcc(vec_PFSLepByPt, leptonID, index_lepton1, index_lepton2);
        if ( index_lepton1 != -1 && index_lepton2 != -1 ) {
          Particle lepton1_PFS = leptons_PFS.particlesByPt()[index_lepton1];
          Particle lepton2_PFS = leptons_PFS.particlesByPt()[index_lepton2];

          const FourMomentum pVec_diLep = lepton1_PFS.mom() + lepton2_PFS.mom();
          double mass = pVec_diLep.mass();

          if      ( leptonID == 11 ) _h_massEEFiducial->fill(mass/GeV);
          else if ( leptonID == 13 ) _h_massMuMuFiducial->fill(mass/GeV);
        }
      } // end of if ( nLepton_PFS >= 2 )
    }


    void FindPFSLeptonPair_HighestPtWithinAcc(vector<Particle>& vec_PFSLepByPt, int pdgID, int& index_lepton1, int& index_lepton2) {
      double pTCut_lead = 0;
      if ( pdgID == 11 ) pTCut_lead = 30.0;
      if ( pdgID == 13 ) pTCut_lead = 22.0;
      double pTCut_sub = 10.0; // same for both channel

      double etaCut_lead = 0;
      if ( pdgID == 11 ) etaCut_lead = 2.5;
      if ( pdgID == 13 ) etaCut_lead = 2.4;
      double etaCut_sub = etaCut_lead;

      int nLepton = int(vec_PFSLepByPt.size());
      for (int i=0; i<nLepton; ++i) {
        auto& lepton = vec_PFSLepByPt[i];
        if ( lepton.isLepton() && lepton.pT() > pTCut_lead && lepton.abseta() < etaCut_lead ) {
          if ( pdgID == 11 ) { // electron channel: check ECAL gap
            if ( !(lepton.abseta() > 1.4442 && lepton.abseta() < 1.566) ) {
              index_lepton1 = i;
              break;
            }
          }
          else { // muon channel
            index_lepton1 = i;
            break;
          }
        }
      } // end of lepton iteration

      // if no lepton is found in the leptons_PFS
      if ( index_lepton1 < 0 ) {
        index_lepton1 = -1;
        index_lepton2 = -1;
        return;
      }

      // 2nd lepton: lepton with highest-pT among the leptons with the opposite sign with 1st lepton
      int pdgID_lepton1 = vec_PFSLepByPt[index_lepton1].pid();
      for (int i=index_lepton1+1; i<nLepton; ++i) { // starting after lepton1
        auto& lepton = vec_PFSLepByPt[i];
        if ( lepton.isLepton() && lepton.pid() == (-1)*pdgID_lepton1 && 
             lepton.pT() > pTCut_sub && lepton.abseta() < etaCut_sub ) {
          if ( pdgID == 11 ) { // electron channel: check ECAL gap
            if ( !(lepton.abseta() > 1.4442 && lepton.abseta() < 1.566) ) {
              index_lepton2 = i;
              break;
            }
          }
          else { // muon channel
            index_lepton2 = i;
            break;
          }
        }
      } // end of lepton iteration
    }
  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(CMS_2018_I1711625);


}
