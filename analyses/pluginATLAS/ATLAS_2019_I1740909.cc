#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  /// @brief jet fragmentation at 13 TeV
  class ATLAS_2019_I1740909: public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2019_I1740909);
    
    /// Book cuts and projections
    void init() {

      const FinalState bare_MU(Cuts::abspid == PID::MUON);

      VetoedFinalState jetinput;
      jetinput.addVetoOnThisFinalState(bare_MU);
      jetinput.addVetoOnThisFinalState(InvisibleFinalState());

      FastJets jetpro(jetinput, FastJets::ANTIKT, 0.4);
      declare(jetpro, "Jets");

      book(_p["nch_jetpt_F"] , 1, 1, 1);
      book(_p["nch_jetpt_C"] , 2, 1, 1);
      book(_p["nch_jetpt_B"] , 9, 1, 1);
    
      for (size_t i_bin = 0; i_bin < 14; ++i_bin) {
        book(_h["nch_B"+to_str(i_bin)] , 13 + i_bin, 1, 1); 
        book(_hr["r_B"+to_str(i_bin)] , 27 + i_bin, 1, 1); 
        book(_h["zeta_B"+to_str(i_bin)] , 41 + i_bin, 1, 1); 
        book(_h["pTrel_B"+to_str(i_bin)] , 55 + i_bin, 1, 1); 
        book(_h["nch_F"+to_str(i_bin)] , 69 + i_bin, 1, 1); 
        book(_hr["r_F"+to_str(i_bin)] , 83 + i_bin, 1, 1); 
        book(_h["zeta_F"+to_str(i_bin)] , 97 + i_bin, 1, 1); 
        book(_h["pTrel_F"+to_str(i_bin)] , 111 + i_bin, 1, 1);  
        book(_h["nch_C"+to_str(i_bin)] , 125 + i_bin, 1, 1); 
        book(_hr["r_C"+to_str(i_bin)] , 139 + i_bin, 1, 1); 
        book(_h["zeta_C"+to_str(i_bin)] , 153 + i_bin, 1, 1); 
        book(_h["pTrel_C"+to_str(i_bin)] , 167 + i_bin, 1, 1);          
        }
    }

    void analyze(const Event& event) {
    
      //Init
      double fnch=0;
      double cnch=0;
      double fzval=0;
      double czval=0;
      double frval=0;
      double crval=0;
      double ftval=0;
      double ctval=0;

      // Event selection 
      Jets m_goodJets = apply<JetAlg>(event, "Jets").jetsByPt(Cuts::pT > 100*GeV && Cuts::abseta < 2.1);
      if (m_goodJets.size() < 2) vetoEvent;
      if (fabs(1.0 - m_goodJets[0].pT()/m_goodJets[1].pT()) > 0.5)  vetoEvent;
      // Decide forward or central   
      bool check = m_goodJets[0].abseta() < m_goodJets[1].abseta();
      int pos_f = int(check);
      int pos_c = int(!check);
      
      // Calculate obs, separately for central and fwd, also get bin
      double fpt = m_goodJets[pos_f].pT();
      size_t f_bin = GetJetBin(fpt);
      double cpt = m_goodJets[pos_c].pT();
      size_t c_bin = GetJetBin(cpt);

      for (const Particle& p : m_goodJets[pos_f].particles()) {
        if (p.pT() < 0.5*GeV)  continue;
        if (p.charge() != 0) {
          ++fnch;
	                
          fzval = p.pT() / m_goodJets[pos_f].pt();        
          ftval = p.pT()*sin(p.phi()-m_goodJets[pos_f].phi());                
          frval = deltaR(m_goodJets[pos_f], p);
               
          _hr["r_F"+to_str(f_bin)]->fill(frval);
          _hr["r_B"+to_str(f_bin)]->fill(frval);                               
          _h["zeta_F"+to_str(f_bin)]->fill(fzval);
          _h["zeta_B"+to_str(f_bin)]->fill(fzval);     
          _h["pTrel_F"+to_str(f_bin)]->fill(ftval);
          _h["pTrel_B"+to_str(f_bin)]->fill(ftval);
        }    
      }
	               
	               
	    for (const Particle& p : m_goodJets[pos_c].particles()) {            
        if (p.pT() < 0.5*GeV)  continue;
        if (p.charge() != 0) {
          ++cnch;
          czval = p.pT() / m_goodJets[pos_c].pt();
          ctval = p.pT()*sin(p.phi()-m_goodJets[pos_c].phi());
          crval = deltaR(m_goodJets[pos_c], p);
                       
          _hr["r_C"+to_str(c_bin)]->fill(crval);
          _hr["r_B"+to_str(c_bin)]->fill(crval);
          _h["zeta_C"+to_str(c_bin)]->fill(czval);
          _h["zeta_B"+to_str(c_bin)]->fill(czval);
          _h["pTrel_C"+to_str(c_bin)]->fill(ctval);
          _h["pTrel_B"+to_str(c_bin)]->fill(ctval);
        }
      }

      if (fnch > 63)  fnch = 63;
      if (cnch > 63)  cnch = 63;
       
       //Fill nchg histo
       
      _p["nch_jetpt_F"]->fill(fpt,fnch);
      _p["nch_jetpt_C"]->fill(cpt,cnch);
      _p["nch_jetpt_B"]->fill(fpt,fnch);
      _p["nch_jetpt_B"]->fill(cpt,cnch);
            
      _h["nch_F"+to_str(f_bin)]->fill(fnch);
      _h["nch_C"+to_str(c_bin)]->fill(cnch);
      _h["nch_B"+to_str(f_bin)]->fill(fnch);
      _h["nch_B"+to_str(c_bin)]->fill(cnch);
    }
 
    
    void finalize() {
    
      // For r only
      /// @todo Replace with barchart()
      for (auto& hist : _hr) {
        for(size_t i=0; i < hist.second->numBins(); ++i) {
          double x = hist.second->bin(i).xMid();
          double bW = hist.second->bin(i).xWidth();
          hist.second->bin(i).scaleW(bW/(2.0*M_PI*x)); 
        }  
      }

      // The rest
      /// @todo Replace with barchart()
      for (auto& hist : _h) {                
        for (size_t i=0; i < hist.second->numBins(); ++i) {
          double bW = hist.second->bin(i).xWidth();
          hist.second->bin(i).scaleW(bW); 
        }
      }
                          
      for (size_t i_bin = 0; i_bin < 14; ++i_bin) {
      
        double sfB =  _h["nch_B"+to_str(i_bin)]->sumW();                         
        if (sfB) {
          scale(_h["zeta_B"+to_str(i_bin)],2.0/sfB);
          scale(_h["pTrel_B"+to_str(i_bin)],2.0/sfB);
          scale(_hr["r_B"+to_str(i_bin)],2.0/sfB);                      
          scale(_h["nch_B"+to_str(i_bin)], 2.0/sfB);
        }
        
        double sfF =  _h["nch_F"+to_str(i_bin)]->sumW();                         
        if (sfF) {
          scale(_h["zeta_F"+to_str(i_bin)],1.0/sfF);
          scale(_h["pTrel_F"+to_str(i_bin)],1.0/sfF);
          scale(_hr["r_F"+to_str(i_bin)],1.0/sfF);                      
          scale(_h["nch_F"+to_str(i_bin)], 1.0/sfF);
        }
        
        double sfC =  _h["nch_C"+to_str(i_bin)]->sumW();                         
        if (sfC) {
          scale(_h["zeta_C"+to_str(i_bin)],1.0/sfC);
          scale(_h["pTrel_C"+to_str(i_bin)],1.0/sfC);
          scale(_hr["r_C"+to_str(i_bin)],1.0/sfC);                      
          scale(_h["nch_C"+to_str(i_bin)], 1.0/sfC);                 
        }
      }
    }

  private:
    
    size_t GetJetBin(const double jetpt){
      size_t i_bin = 0;
      if (inRange(jetpt,100,200)) i_bin=0;
      if (inRange(jetpt,200,300)) i_bin=1;
      if (inRange(jetpt,300,400)) i_bin=2;
      if (inRange(jetpt,400,500)) i_bin=3;
      if (inRange(jetpt,500,600)) i_bin=4;
      if (inRange(jetpt,600,700)) i_bin=5;
      if (inRange(jetpt,700,800)) i_bin=6;
      if (inRange(jetpt,800,900)) i_bin=7;
      if (inRange(jetpt,900,1000)) i_bin=8;
      if (inRange(jetpt,1000,1200)) i_bin=9;
      if (inRange(jetpt,1200,1400)) i_bin=10;
      if (inRange(jetpt,1400,1600)) i_bin=11;
      if (inRange(jetpt,1600,2000)) i_bin=12;
      if (inRange(jetpt,2000,2500)) i_bin=13;   
      if(jetpt < 100) i_bin=0;
      if(jetpt > 2500) i_bin=13; 
      return i_bin;
    } 
  
    map<string, Histo1DPtr> _h;
    map<string, Histo1DPtr> _hr;
    map<string, Profile1DPtr> _p; 
  };
  
  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2019_I1740909);
  
}
