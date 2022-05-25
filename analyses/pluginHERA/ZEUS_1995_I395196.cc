// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/DISLepton.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/DISRapidityGap.hh"
namespace Rivet {


  /// @brief Neutral strange particle production in deep inelastic scattering at HERA (ZEUS)
  class ZEUS_1995_I395196 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ZEUS_1995_I395196);


    void init() {
     
      declare(DISLepton(), "Lepton");
      declare(DISKinematics(), "Kinematics");
      declare(DISRapidityGap(), "Rapidity Gap"); 
    
      const Cut cut = Cuts::abseta < 1.3;

      const FinalState fs(cut);
      declare(fs, "FS");

      const ChargedFinalState cfs(cut);
      declare(cfs, "CFS");

       // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      book(_h["pT_kaon"], 1, 1, 1);
      book(_h["eta_kaon"], 2, 1, 1);
      book(_h["pT_lambda"], 3, 1, 1);
      book(_h["eta_lambda"], 4, 1, 1);
      book(_h_multK0_0,"TMP/mult_0", refData(5,1,1));
      book(_h_multK0_1,"TMP/mult_1", refData(5,1,1));
      book(_h_multK0_2,"TMP/mult_2", refData(6,1,1));
      book(_h_multK0_3,"TMP/mult_3", refData(6,1,1));
      book(_h_scatratio, 6, 1, 1);
      book(_h["K0_NRG_data_pT"], 7, 1, 1);
      book(_h["K0_LRG_data_pT"], 8, 1, 1);
      book(_h["K0_NRG_data_eta"], 9, 1, 1);
      book(_h["K0_LRG_data_eta"], 10, 1, 1);
      book(_h_scat,5,1,1);
      book(_c["dis"],"TMP/Nevt_after_cuts");

      
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const FinalState& fs = apply<FinalState>(event, "FS");
      const DISKinematics& dk = applyProjection<DISKinematics>(event, "Kinematics");
      
      const DISRapidityGap& g = applyProjection<DISRapidityGap>(event, "Rapidity Gap");
  
      const ChargedFinalState& charged = apply<ChargedFinalState>(event, "CFS");
      const size_t numPartcharged = charged.particles().size();
      //const size_t numPart = fs.particles().size();      
      //_c["charged"] -> fill(numPartcharged);  
      //_c["all"] -> fill(numPart);
      const size_t numParticles = fs.particles().size();
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      
      double rgap = g.gap();   
      // Get the DIS kinematics
      double xbj  = dk.x();
      double ybj = dk.y();
      double Q2 = dk.Q2()/GeV;
      double W = sqrt(dk.W2()/GeV);
      bool cut = Q2 >10 && Q2<640 && xbj>0.0003 && xbj<0.01 && ybj>0.04 && ybj<1.0;
      if (!cut) vetoEvent;
      _h_multK0_1 -> fill(Q2);
      _h_multK0_2 -> fill(Q2,numPartcharged);     
    
      _c["dis"] -> fill();
      int kaon=0;
      int lambda=0;

      for(const Particle& p : fs.particles()){
          const double eta= p.eta();
          const double pT = p.pT()/GeV;
	  const int pid = abs(p.pid());
          //const double ybj= (p.E()-p.pz())/(2*27.5);
          if (pid == 310 || pid == 130) {  //K0S
             //cout << " pid " << pid << " eta " << eta << endl;
             if (pT>0.5  && pT<4.0){
                kaon++ ;
                //fill histograms related to the kaons in here. 
                _h["pT_kaon"] -> fill(pT,0.5/pT);
                _h["eta_kaon"] -> fill(eta);
                _h_multK0_0 -> fill(Q2);
                _h_multK0_3 -> fill(Q2);
                if(rgap<1.5 && W>140.0 ) {
                   _h["K0_LRG_data_pT"] -> fill(pT,0.5/pT);
                   _h["K0_LRG_data_eta"] -> fill(eta);
                   //cout<< abs(eta) <<endl;
                }
                else if(rgap>1.5 && W>140.0) {
                   _h["K0_NRG_data_pT"] -> fill(pT,0.5/pT);
                   _h["K0_NRG_data_eta"] -> fill(eta);
                } 
             }
          }
          else if (pid==3122){ // Lambda 
             if (pT>0.5  && pT<3.5){
                lambda++ ;
                //fill histograms related to the lambdas  in here. 
                _h["pT_lambda"] -> fill(pT,0.5/pT);
                _h["eta_lambda"] -> fill(eta);
             }
          }
     }
     
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      divide(_h_multK0_0, _h_multK0_1, _h_scat);   
      divide(_h_multK0_3, _h_multK0_2, _h_scatratio);
      

      //cout<< "#of kaons per events"<< kaon/numEvents() <<endl;
      //cout<< "Num mean charged p multiplicity"<< *_c["charged"]<< endl;      
      scale(_h["pT_kaon"],1./ *_c["dis"]);
      scale(_h["eta_kaon"],1./ *_c["dis"]);
      scale(_h["pT_lambda"],1./ *_c["dis"]);
      scale(_h["eta_lambda"],1./ *_c["dis"]); 
      scale(_h["K0_LRG_data_pT"],1./ *_c["dis"]);
      scale(_h["K0_NRG_data_pT"],1./ *_c["dis"]);
      scale(_h["K0_LRG_data_eta"],1./ *_c["dis"]);
      scale(_h["K0_NRG_data_eta"],1./ *_c["dis"]);

    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    ///@}
   private:
     Scatter2DPtr _h_scat, _h_scatratio ;
     Histo1DPtr _h_multK0_0 , _h_multK0_1 ,_h_multK0_2 ,_h_multK0_3 ;
  };


  RIVET_DECLARE_PLUGIN(ZEUS_1995_I395196);

}


