// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/LundGenerator.hh"

namespace Rivet {

  /// @brief Lund jet plane with charged particles
  class ATLAS_2020_I1790256: public Analysis {
  public:

    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2020_I1790256);

    /// @name Analysis methods
    //@{

    void init() {
    
      //Projections
      FinalState fs(Cuts::abseta < 4.5); 
      FastJets jet4(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jet4, "Jets");
      
      ChargedFinalState tracks(Cuts::pT > 0.5*GeV && Cuts::abseta < 2.5);
      declare(tracks, "tracks");    
  
      book(_h_lundplane, 1,1,1);  
     
      _h_vs.resize(13);
      for (size_t i = 0; i < _h_vs.size(); ++i) {
        book(_h_vs[i] , i+3 , 1, 1); 
      }
            
      _h_hs.resize(19);
      for (size_t i = 0; i < _h_hs.size(); ++i) {
        book(_h_hs[i], i+16, 1, 1); 
      }        
        
        
      book(_njets, "_njets");  
                      
    }

    void analyze(const Event& event) {
 
      const Jets jets = apply<JetAlg>(event, "Jets").jetsByPt(Cuts::pT > 300*GeV && Cuts::abseta < 2.1);  
         
      if (jets.size() < 2)  vetoEvent;
      if (jets[0].pT() < 675*GeV)  vetoEvent;
      
      if ( (jets[0].pT()/jets[1].pT()) > 1.5 ) vetoEvent;

       _njets->fill(2);

       const Particles& tracks = apply<ChargedFinalState>(event, "tracks").particlesByPt();
 
       Particles intracks1;
       Particles intracks2;

       const Jet& j1 = jets[0];
       const Jet& j2 = jets[1];


      for (const Particle& p : tracks) {
        const double dr = deltaR(j1, p, PSEUDORAPIDITY);
        if (dr > 0.4) continue;
	      if (abs(p.pid()) == 13) continue;
        intracks1.push_back(p);
      }

      for (const Particle& p : tracks) {
        const double dr = deltaR(j2, p, PSEUDORAPIDITY);
        if (dr > 0.4) continue;
        if (abs(p.pid()) == 13) continue;
        intracks2.push_back(p);
      }

      JetDefinition tjet1_def(fastjet::cambridge_algorithm, 10);
      ClusterSequence tjet1_cs(intracks1, tjet1_def);     
      vector<PseudoJet> tjets1 = fastjet::sorted_by_pt(tjet1_cs.inclusive_jets(0.0));  
     
      JetDefinition tjet2_def(fastjet::cambridge_algorithm, 10);
      ClusterSequence tjet2_cs(intracks2, tjet2_def);     
      vector<PseudoJet> tjets2 = fastjet::sorted_by_pt(tjet2_cs.inclusive_jets(0.0));    
    
      if (tjets1.size() < 1 || tjets2.size() < 1) vetoEvent;
    
      fjcontrib::LundGenerator lund;
      vector<fjcontrib::LundDeclustering> declusts1 = lund(tjets1[0]);
      for (size_t idecl = 0; idecl < declusts1.size(); ++idecl) {
        pair<double,double> coords = declusts1[idecl].lund_coordinates();
        double X = -0.9163 + coords.first;
        double Y = - log(declusts1[idecl].z());
            
        if (X > 0 && X < 4.33 && Y > log(1/0.5)  && Y < 8.6*log(1/0.5) ){
                 
          _h_lundplane->fill(X, Y);

          double hdiv = (double)4.33/(double)13;    
          size_t i = floor(X/hdiv);
          _h_vs[i]->fill(Y);
                
          double vdiv = (8.6*log(1/0.5) - log(1/0.5))/(double)19; 
          size_t j = floor((Y - log(1/0.5))/vdiv);
          _h_hs[j]->fill(X);
                 
        }
      }
  
      vector<fjcontrib::LundDeclustering> declusts2 = lund(tjets2[0]);
      for (size_t idecl = 0; idecl < declusts2.size(); ++idecl) {
        pair<double,double> coords = declusts2[idecl].lund_coordinates();

        double X = -0.9163 + coords.first;
        double Y = - log(declusts2[idecl].z());

        if (X > 0 && X < 4.33 && Y > log(1/0.5)  && Y < 8.6*log(1/0.5) ) {

          _h_lundplane->fill(X, Y);

          double hdiv = (double)4.33/(double)13;    
          size_t i = floor(X/hdiv);
          _h_vs[i]->fill(Y);
                
          double vdiv = (8.6*log(1/0.5) - log(1/0.5))/(double)19; 
          size_t j = floor((Y - log(1/0.5))/vdiv);
          _h_hs[j]->fill(X);
        }
      }
    }


    void finalize() {
    
      double area = _njets->sumW();
      scale(_h_lundplane, 1/area);
      scale(_h_vs, 1/(area*0.333));
      scale(_h_hs, 1/(area*0.277)); 

    }

  private:

   
    Histo2DPtr _h_lundplane;
    vector<Histo1DPtr> _h_vs, _h_hs;
    CounterPtr _njets;
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2020_I1790256);
}



    

