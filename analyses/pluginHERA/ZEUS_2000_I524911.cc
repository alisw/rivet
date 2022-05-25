// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/DISLepton.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Measurement of azimuthal asymmetries in deep inelastic scattering (ZEUS)
  class ZEUS_2000_I524911 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ZEUS_2000_I524911);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(DISLepton(), "Lepton");
      declare(DISKinematics(), "Kinematics");

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance

      const ChargedFinalState cfs;
      declare(cfs, "CFS");

      book(_h["A1"], 1, 1, 1);
      book(_h["A2"], 1, 1, 2);
      book(_h["A3"], 1, 1, 3);
      book(_h["A4"], 1, 1, 4);
                
      book(_p["cosphi"],2, 1, 1) ;
      book(_p["cos2phi"],2, 1, 2) ;
   
// counter pointer to store the no. of events           
      book(_Nevt_after_cuts, "TMP/Nevt_after_cuts");

      

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
    
    //const FinalState& fsall = apply<FinalState>(event, "FS");
    const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");

    const DISKinematics& dk = applyProjection<DISKinematics>(event, "Kinematics");
    const DISLepton& dl = applyProjection<DISLepton>(event,"Lepton");

    double x = dk.x();
    double y = dk.y();
    const double Q2 = dk.Q2();
    
     double PT[] = { 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0};
    // Extract the particles other than the lepton

      if(x<0.01||x>0.1) vetoEvent;
      if(y<0.2||y>0.8) vetoEvent;
      if(Q2 > 7220 || Q2 < 180) vetoEvent;

      _Nevt_after_cuts -> fill();
        
      Particles particles;
      particles.reserve(cfs.particles().size());
      
      ConstGenParticlePtr dislepGP = dl.out().genParticle();
      for (const Particle& p : cfs.particles()) {
          ConstGenParticlePtr loopGP = p.genParticle();

          if (loopGP == dislepGP) continue;
          particles.push_back(p);
      }
    
    const LorentzTransform hcmboost = dk.boostHCM();
    for (size_t ip1 = 0; ip1 < particles.size(); ++ip1) {
        const Particle& p = particles[ip1];
       
       // calculate zh        
        double zh = 2.*x/Q2* (dk.beamHadron().E()*p.momentum().E() - dk.beamHadron().pz()*p.momentum().pz()) ;
       // cout << " zh " << zh << endl;
       // Boost to hcm
       
       if (zh < 0.2 ) continue ;

         const FourMomentum hcmMom = hcmboost.transform(p.momentum());      
                  
         const double phi =mapAngleMPiToPi(hcmMom.phi())/degree ;
         
         
         
// Filling histograms with values of cos(phi) and cos(2phi) wrt the corresponding momentum cuts

         for (size_t i = 0; i < 8; ++i) {
            if(hcmMom.pT() > PT[i] ) { 
             _p["cosphi"]->fill(i+1,cos(hcmMom.phi()));
             _p["cos2phi"]->fill(i+1,cos(2.*hcmMom.phi()));
             }
         }
         
         
                    
         if(hcmMom.pT() > PT[1] ) { _h["A1"] -> fill(phi); }          
                     
         if(hcmMom.pT() > PT[3] ) { _h["A2"] -> fill(phi); }
         
         if(hcmMom.pT() > PT[5] ) { _h["A3"] -> fill(phi); }
                   
         if(hcmMom.pT() > PT[7] ) { _h["A4"] -> fill(phi); }
        
     }
        
  
  }



    /// Normalise histograms etc., after the run
    void finalize() {

      // correct binwidth in degree to correct for binning from degree to rad by: binwidth/(2PI/10.)
      double norm = dbl(*_Nevt_after_cuts) ;
   //   cout << " Nev " << norm << " bin_width= " <<_h["A2"]->bin(0).xWidth() << endl;
      double degTOrad_width = _h["A1"]->bin(0).xWidth()*10./2./M_PI ;
      if (norm > 1 ) {
         scale(_h["A1"], degTOrad_width/norm); 
         scale(_h["A2"], degTOrad_width/norm); 
         scale(_h["A3"], degTOrad_width/norm);
         scale(_h["A4"], degTOrad_width/norm); 
      }     
      
    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    CounterPtr _Nevt_after_cuts;
    ///@}
    


  };


  DECLARE_RIVET_PLUGIN(ZEUS_2000_I524911);

}
