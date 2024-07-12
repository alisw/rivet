// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/DISFinalState.hh"
#include <cmath>

namespace Rivet {


  /// @brief DIS dijets in the breit frame
  class ZEUS_2010_I875006 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ZEUS_2010_I875006);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(DISKinematics(), "Kinematics");

      // All final state particles boosted to Breit frame then clustered 
      //using FastJet KT algorithm with jet radius parameter 1      
      const DISFinalState DISfs(DISFinalState::BoostFrame::BREIT);
      FastJets DISjetfs(DISfs, FastJets::KT, 1.0);
      declare(DISjetfs, "DISjets");

      // Book histograms
      // specify custom binning

      //Non-Grouped histgrams
      book(_h_Q2, 1,1,1);
      book(_h_XBj, 2,1,1);
      book(_h_Et, 3,1,1);
      book(_h_Mjj, 4,1,1);
      book(_h_Eta, 5,1,1);
      book(_h_Zeta, 6,1,1);

      //Zeta values seperated into Q2 ranges
      book(_h_ZetaQ2[0], 7,1,1);
      book(_h_ZetaQ2[1], 8,1,1);
      book(_h_ZetaQ2[2], 9,1,1);
      book(_h_ZetaQ2[3], 10,1,1);
      book(_h_ZetaQ2[4], 11,1,1);
      book(_h_ZetaQ2[5], 12,1,1);

      //Transverse jet energy seperated into Q2 ranges
      book(_h_EtQ2[0], 13,1,1);
      book(_h_EtQ2[1], 14,1,1);
      book(_h_EtQ2[2], 15,1,1);
      book(_h_EtQ2[3], 16,1,1);
      book(_h_EtQ2[4], 17,1,1);
      book(_h_EtQ2[5], 18,1,1);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      
      //First Lorentz invariant quantities in Lab frame
      DISKinematics dis = apply<DISKinematics>(event, "Kinematics");
      double Q2 = dis.Q2();
      double xbj = dis.x();
      double y = dis.y();

      //Perform required cut on Q2 and y
      if (!inRange(Q2, 125*GeV2, 20000*GeV2)) vetoEvent;
      if (!inRange(y, 0.2, 0.6)) vetoEvent;
	
      //Get Lorentz transforms for Breit Boost and Lab Boost	
      const LorentzTransform breitboost = dis.boostBreit();
      const LorentzTransform labboost = breitboost.inverse();	

      //Get jets clustered in Breit frame 
      Jets jets = apply<FastJets>(event, "DISjets").jetsByPt();
	
      //Boost jets to lab frame
      for(std::vector<int>::size_type i=0; i<jets.size(); i++){
        jets[i].transformBy(labboost);
      }
    
      //Cut on Pseudorapidity in lab frame 
      const int orientation = dis.orientation();
      vector<Jet> cutJets;
      for(std::vector<int>::size_type i=0; i<jets.size(); i++){
  	double etaJet = jets[i].eta()*orientation;
        if(etaJet < 2.5 && etaJet > -1){
		cutJets.push_back(jets[i]);
	}
      }
          
      //veto event if only single jet
      if(cutJets.size()<2){
          vetoEvent;
      }
      
      //Boost jets to Breit frame 
      for(std::vector<int>::size_type i=0; i<cutJets.size(); i++){
      	cutJets[i].transformBy(breitboost);
      }
      
      //Sort jets by et in descending order
      std::sort(cutJets.begin(), cutJets.end(),[](const Jet& j1, const Jet& j2){
	 return j1.Et()>j2.Et();
      }); 

      //Ensure two hardest jets have Et>8GeV in Breit frame
      const Jet& jet1 = cutJets[0];
      const Jet& jet2 = cutJets[1];

      if(jet1.Et()<8*GeV || jet2.Et()<8*GeV){
      	  vetoEvent;
      }
      
      //Extract required quantities in Breit frame
      //Dijet mean transverse energy 
      const double dijetEt = (jet1.Et() + jet2.Et())/2;      

      //Invariant dijet mass of hardest transverse jets > 20GeV
      const double Mjj =  FourMomentum(jet1.momentum() + jet2.momentum()).mass();
      if(Mjj<20*GeV){
          vetoEvent;
      } 

      const double eta1 = orientation*jet1.eta();
      const double eta2 = orientation*jet2.eta();
      const double etastar = abs(eta1 - eta2)/2;

      const double logZeta = log10(xbj*(1 + pow(Mjj,2)/Q2));

      //Fill histograms
      _h_Q2->fill(Q2);
      _h_XBj->fill(xbj);
      _h_Et->fill(dijetEt);
      _h_Mjj->fill(Mjj);
      _h_Eta->fill(etastar);
      _h_Zeta->fill(logZeta);

      //Fill histograms for different Q2 ranges
      if(Q2>125*GeV2 && Q2<=250*GeV2){
          _h_ZetaQ2[0]->fill(logZeta);
          _h_EtQ2[0]->fill(dijetEt);
      }else if (Q2>250*GeV2 && Q2<=500*GeV2){
          _h_ZetaQ2[1]->fill(logZeta);
          _h_EtQ2[1]->fill(dijetEt);
      }else if (Q2>500*GeV2 && Q2<=1000*GeV2){
          _h_ZetaQ2[2]->fill(logZeta);
          _h_EtQ2[2]->fill(dijetEt);
      }else if (Q2>1000*GeV2 && Q2<=2000*GeV2){
          _h_ZetaQ2[3]->fill(logZeta);
          _h_EtQ2[3]->fill(dijetEt);
      }else if (Q2>2000*GeV2 && Q2<=5000*GeV2){
          _h_ZetaQ2[4]->fill(logZeta);
          _h_EtQ2[4]->fill(dijetEt);
      }else if (Q2>5000*GeV2 && Q2<=20000*GeV2){
          _h_ZetaQ2[5]->fill(logZeta);
          _h_EtQ2[5]->fill(dijetEt);
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      //Calculate scaling factor from cross section
      const double sf = crossSection()/picobarn/sumW(); //Scale factor with cuts

      scale(_h_Q2, sf);
      scale(_h_XBj, sf);
      scale(_h_Et, sf);
      scale(_h_Mjj, sf);
      scale(_h_Eta, sf);
      scale(_h_Zeta, sf);

      for(int i = 0; i<6;i++){
          scale(_h_ZetaQ2[i], sf);
          scale(_h_EtQ2[i],sf);
      } 
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_Q2;	
    Histo1DPtr _h_XBj;
    Histo1DPtr _h_Et;
    Histo1DPtr _h_Mjj;
    Histo1DPtr _h_Eta;	
    Histo1DPtr _h_Zeta;
    Histo1DPtr _h_ZetaQ2[6];
    Histo1DPtr _h_EtQ2[6];	
    ///@}
  };


  RIVET_DECLARE_PLUGIN(ZEUS_2010_I875006);

}
