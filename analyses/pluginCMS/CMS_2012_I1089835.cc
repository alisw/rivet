// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {
  
  /// @brief Inclusive b-jet production in pp collisions at 7 TeV
  class CMS_2012_I1089835 : public Analysis {
  public:
    
    /// @name Constructors etc.
    //@{
    
    /// Constructor
    CMS_2012_I1089835()
      : Analysis("CMS_2012_I1089835")
    {        
      
    }
    
    //@}
    
    
  public:
    
    /// @name Analysis methods
    //@{
    
    /// Book histograms and initialise projections before the run
    void init() {
      
      const FinalState cnfs;
      declare(cnfs, "FS");
      declare(FastJets(cnfs, FastJets::ANTIKT, 0.5), "Jets");
      
      book(_h_dsigdpty05, 4, 1, 1);
      book(_h_dsigdpty10, 5, 1, 1);
      book(_h_dsigdpty15, 6, 1, 1);
      book(_h_dsigdpty20, 7, 1, 1);
      book(_h_dsigdpty22, 8, 1, 1);
      
    }
    
    
    /// Perform the per-event analysis
    void analyze(const Event& event) {
      
      const FastJets& fastjets = apply<FastJets>(event, "Jets"); 
      const Jets jets = fastjets.jetsByPt(10.);
      
      for (const Jet& j : jets) {
        
        const double ptB= j.pT();
        const double yB= j.rapidity();
        
        if (j.bTagged()) {
          
          if( fabs(yB) < 0.5) { _h_dsigdpty05->fill( ptB, 1.0 );}
          else if( fabs(yB) >= 0.5 && fabs(yB) < 1.0) { _h_dsigdpty10->fill( ptB, 1.0 );}
          else if( fabs(yB) >= 1.0 && fabs(yB) < 1.5) { _h_dsigdpty15->fill( ptB, 1.0 );}
          else if( fabs(yB) >= 1.5 && fabs(yB) < 2.0) { _h_dsigdpty20->fill( ptB, 1.0 );}
          else if( fabs(yB) >= 2.0 && fabs(yB) < 2.2) { _h_dsigdpty22->fill( ptB, 1.0 );}
        }
      }
    }
    
    /// Normalise histograms etc., after the run
    void finalize() {
      
      double invlumi = crossSection()/picobarn/sumOfWeights();
      
      scale(_h_dsigdpty05, invlumi); 
      scale(_h_dsigdpty10, invlumi); 
      scale(_h_dsigdpty15, invlumi); 
      scale(_h_dsigdpty20, invlumi); 
      scale(_h_dsigdpty22, invlumi/0.4); 
      
    }
    
    //@}
    
    
  private:
    
    /// @name Histograms
    //@{
    Histo1DPtr _h_dsigdpty05;
    Histo1DPtr _h_dsigdpty10;
    Histo1DPtr _h_dsigdpty15;
    Histo1DPtr _h_dsigdpty20;
    Histo1DPtr _h_dsigdpty22;
    //@}
    
    
  };
  
  
  
  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2012_I1089835);

}
