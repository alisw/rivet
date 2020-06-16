// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"

namespace Rivet {

  /// @brief H(125)->ZZ->4l at 8 TeV
  class ATLAS_2014_I1310835 : public Analysis {
  public:

    /// Default constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2014_I1310835);

    void init() {
      const FinalState fs(Cuts::abseta < 5.0);

      PromptFinalState photons(Cuts::abspid == PID::PHOTON);

      PromptFinalState bare_el(Cuts::abspid == PID::ELECTRON);

      PromptFinalState bare_mu(Cuts::abspid == PID::MUON);

      // Selection: lepton selection
      Cut etaranges_el = Cuts::abseta < 2.47 && Cuts::pT > 7*GeV; 
      DressedLeptons electron_sel4l(photons, bare_el, 0.1, etaranges_el, false);
      declare(electron_sel4l, "electrons");
 
      Cut etaranges_mu = Cuts::abseta < 2.7 && Cuts::pT > 6*GeV;
      DressedLeptons muon_sel4l(photons, bare_mu, 0.1, etaranges_mu, false);
      declare(muon_sel4l, "muons");

      FastJets jetpro(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetpro, "jet");

      // Book histos
      book(_h_pt          , 1, 1, 1);
      book(_h_rapidity    , 2, 1, 1);
      book(_h_m34         , 3, 1, 1);
      book(_h_costheta    , 4, 1, 1);
      book(_h_njets       , 5, 1, 1);
      book(_h_leadingjetpt, 6, 1, 1);

    }



    /// Do the analysis
    void analyze(const Event& e) {
      
      ////////////////////////////////////////////////////////////////////
      // preselection of leptons for ZZ-> llll final state
      ////////////////////////////////////////////////////////////////////

      const vector<DressedLepton>& mu_sel4l = applyProjection<DressedLeptons>(e, "muons").dressedLeptons();
      const vector<DressedLepton>& el_sel4l = applyProjection<DressedLeptons>(e, "electrons").dressedLeptons();

      vector<DressedLepton> leptonsFS_sel4l;
      leptonsFS_sel4l.insert( leptonsFS_sel4l.end(), mu_sel4l.begin(), mu_sel4l.end() );
      leptonsFS_sel4l.insert( leptonsFS_sel4l.end(), el_sel4l.begin(), el_sel4l.end() );

      /////////////////////////////////////////////////////////////////////////////
      /// H->ZZ->4l pairing
      /////////////////////////////////////////////////////////////////////////////
 
      size_t el_p = 0;
      size_t el_n = 0;
      size_t mu_p = 0; 
      size_t mu_n = 0;
      
      for (const Particle& l : leptonsFS_sel4l) {
        if (l.abspid() == PID::ELECTRON) {
          if (l.pid() < 0)  ++el_n;
          if (l.pid() > 0)  ++el_p;
        }
        else if (l.abspid() == PID::MUON) {
          if (l.pid() < 0)  ++mu_n;
          if (l.pid() > 0)  ++mu_p;
        }
      }
            
      bool pass_sfos = ( (el_p >=2 && el_n >=2) || (mu_p >=2 && mu_n >=2) || (el_p >=1 && el_n >=1 && mu_p >=1 && mu_n >=1) );
      
      if (!pass_sfos)  vetoEvent;

      Zstate Z1, Z2, Zcand;
      size_t n_parts = leptonsFS_sel4l.size();
      size_t l1_index = 0;
      size_t l2_index = 0;

      // determine Z1 first
      double min_mass_diff = -1;
      for (size_t i = 0; i < n_parts; ++i) {
        for (size_t j = 0; j < n_parts; ++j) {
          if (i >= j)  continue;

          if (leptonsFS_sel4l[i].pid() != -1*leptonsFS_sel4l[j].pid())  continue; //only pair SFOS leptons

          Zcand = Zstate( ParticlePair(leptonsFS_sel4l[i], leptonsFS_sel4l[j]) );
          double mass_diff = fabs( Zcand.mom().mass() - 91.1876 );
         
          if (min_mass_diff == -1 || mass_diff < min_mass_diff) {
            min_mass_diff = mass_diff;
            Z1 = Zcand;
            l1_index = i;
            l2_index = j;
          }
        }
      }

      //determine Z2 second
      min_mass_diff = -1;
      for (size_t i = 0; i < n_parts; ++i) {
        if (i == l1_index || i == l2_index)  continue;
        for (size_t j = 0; j < n_parts; ++j) {
          if (j == l1_index || j == l2_index || i >= j)  continue;

          if (leptonsFS_sel4l[i].pid() != -1*leptonsFS_sel4l[j].pid())  continue; // only pair SFOS leptons

          Zcand = Zstate( ParticlePair(leptonsFS_sel4l[i], leptonsFS_sel4l[j]) );
          double mass_diff = fabs( Zcand.mom().mass() - 91.1876 );

          if (min_mass_diff == -1 || mass_diff < min_mass_diff) {
            min_mass_diff = mass_diff;
            Z2 = Zcand;
          }
        }
      }

      Particles leptons_sel4l;
      leptons_sel4l.push_back(Z1.first);
      leptons_sel4l.push_back(Z1.second);
      leptons_sel4l.push_back(Z2.first);
      leptons_sel4l.push_back(Z2.second);

      ////////////////////////////////////////////////////////////////////////////
      // Kinematic Requirements
      ///////////////////////////////////////////////////////////////////////////
      
      //leading lepton pT requirement
      std::vector<double> lepton_pt;
      for (const Particle& i : leptons_sel4l) lepton_pt.push_back(i.pT() / GeV);
      std::sort(lepton_pt.begin(), lepton_pt.end(), [](const double pT1, const double pT2) { return pT1 > pT2; });
      
      if (!(lepton_pt[0] > 20*GeV && lepton_pt[1] > 15*GeV && lepton_pt[2] > 10*GeV))  vetoEvent;
      
      //invariant mass requirements
      if (!(inRange(Z1.mom().mass(), 50*GeV, 106*GeV) && inRange(Z2.mom().mass(), 12*GeV, 115*GeV)))  vetoEvent;
      
      //lepton separation requirements
      for (unsigned int i = 0; i < 4; ++i) {
        for (unsigned int j = 0; j < 4; ++j) {
          if (i >= j) continue;
          double dR = deltaR(leptons_sel4l[i], leptons_sel4l[j]);
          bool sameflavor = leptons_sel4l[i].abspid() == leptons_sel4l[j].abspid();

          if ( sameflavor && dR < 0.1)  vetoEvent;
          if (!sameflavor && dR < 0.2)  vetoEvent;
        }
      }

      // J/Psi veto requirement
      for (unsigned int i = 0; i < 4; ++i) {
        for (unsigned int j = 0; j < 4; ++j) {
          if (i >= j) continue;
          if ( leptons_sel4l[i].pid() != -1*leptons_sel4l[j].pid() )  continue;
          if ((leptons_sel4l[i].momentum() + leptons_sel4l[j].momentum()).mass() <= 5*GeV)  vetoEvent;
        }
      }
 
      // 4-lepton invariant mass requirement
      double m4l = (Z1.mom() + Z2.mom()).mass();
      if (!(inRange(m4l, 118*GeV, 129*GeV)))  vetoEvent;
  
  
      ////////////////////////////////////////////////////////////////////////////
      // Higgs observables
      ///////////////////////////////////////////////////////////////////////////
      FourMomentum Higgs = Z1.mom() + Z2.mom();

      double H4l_pt       = Higgs.pt()/GeV; 
      double H4l_rapidity = Higgs.absrap(); 
      LorentzTransform HRF_boost;
      //HRF_boost.mkFrameTransformFromBeta(Higgs.betaVec());
      HRF_boost.setBetaVec(- Higgs.betaVec());
      FourMomentum Z1_in_HRF = HRF_boost.transform( Z1.mom() );
      double H4l_costheta = fabs(cos( Z1_in_HRF.theta())); 
      double H4l_m34      = Z2.mom().mass()/GeV;
      
      ////////////////////////////////////////////////////////////////////////////
      // Jet observables
      ///////////////////////////////////////////////////////////////////////////

      Jets jets;
      for (const Jet& jet : applyProjection<FastJets>(e, "jet").jetsByPt(Cuts::pT > 30*GeV && Cuts::absrap < 4.4)) {
        bool overlaps = false;
        for (const Particle& lep : leptonsFS_sel4l) {
          if (lep.abspid() != PID::ELECTRON)  continue;
          const double dR = deltaR(lep, jet);
          if (dR < 0.2) { overlaps = true; break; }
        }
        if (!overlaps) jets += jet;
      }
      size_t n_jets = jets.size();
      if (n_jets > 3)  n_jets = 3;

      std::vector<double> jet_pt;
      for (const Jet& i : jets) jet_pt.push_back(i.pT()/GeV);

      double leading_jet_pt = n_jets? jet_pt[0] : 0.;

      ////////////////////////////////////////////////////////////////////////////
      // End of H->ZZ->llll selection: now fill histograms
      ////////////////////////////////////////////////////////////////////////////


      _h_pt->fill(H4l_pt);
      _h_rapidity->fill(H4l_rapidity);
      _h_costheta->fill(H4l_costheta);
      _h_m34->fill(H4l_m34);
      _h_njets->fill(n_jets + 1);
      _h_leadingjetpt->fill(leading_jet_pt);


    }


    /// Generic Z candidate
    struct Zstate : public ParticlePair {
      Zstate() { }
      Zstate(ParticlePair _particlepair) : ParticlePair(_particlepair) { }
      FourMomentum mom() const { return first.momentum() + second.momentum(); }
      operator FourMomentum() const { return mom(); }
    };

    /// Finalize
    void finalize() {

      const double norm = crossSection()/sumOfWeights()/femtobarn;
      std::cout << "xsec: " << crossSection() << '\n';
      std::cout << "sumw: " << sumOfWeights() << '\n';
      std::cout << "femb: " << femtobarn << '\n';
      std::cout << "norm: " << norm << '\n';

      scale(_h_pt, norm);
      scale(_h_rapidity, norm);
      scale(_h_costheta, norm);
      scale(_h_m34, norm);
      scale(_h_njets, norm);
      scale(_h_leadingjetpt, norm);
    }


  private:

    Histo1DPtr _h_pt, _h_rapidity, _h_costheta;
    Histo1DPtr _h_m34, _h_njets, _h_leadingjetpt;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1310835);

}
