// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Measurement of D* meson cross-sections at HERA
  class H1_1999_I481112 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(H1_1999_I481112);


    /// @name Analysis methods
    ///@{
    

    /// Book histograms and initialise projections before the run
    void init() {
     
      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 1.5);
      declare(fs, "fs");
      // The final-state particles declared above are clustered using FastJet with
    
      //Initialize quantities needed for cuts
      declare(DISKinematics(), "Kinematics");
      declare(UnstableParticles(), "DStars");
      

      Histo1DPtr dummy;
     
      book(_h["211"],2,1,1);
      book(_h["311"],3,1,1);
      book(_h["411"],4,1,1);
      book(_h["511"],5,1,1);
      book(_h["611"],6,1,1);
      book(_h["rap194"], 7,1,1);
      book(_h["pt194"], 8,1,1);
      book(_h["rap88"], 9,1,1);
      book(_h["pt88"], 10,1,1);
      _hpt.add(2.5, 3.5, book(dummy, 11,1,1));
      _hpt.add(3.5, 5.0, book(dummy, 11,1,2));
      _hpt.add(5.0, 10.0, book(dummy, 11,1,3));
      book(_h["1211"],12,1,1);
      book(_h["1212"],12,1,2);
      book(_h["1311"],13,1,1);
    }
    
    
    /// Perform the per-event analysis
    void analyze(const Event& event) {
     
      
      const DISKinematics& dk = applyProjection<DISKinematics>(event, "Kinematics");
      
      
      bool isDIS = false;
      bool ETAG44 = false;
      bool ETAG33 = false;

      const double y = dk.y();
      const double Q2 = dk.Q2();

      if(Q2 > 2 && Q2 <100 &&  y > 0.05 && y < 0.7) isDIS = true;
      if(Q2 < 0.009 && y > 0.02 && y <0.32  ) ETAG44 = true;
      if(Q2 < 0.01 &&  y > 0.29 && y <0.62 ) ETAG33 = true;
      if(isDIS == false && ETAG44 == false && ETAG33 == false) vetoEvent;
      
      
      //Creating array of D*
      Particles unstables ;
      if (isDIS ) { unstables = apply<ParticleFinder>(event, "DStars").particles(Cuts::pT > 1.5*GeV && Cuts::absetaIn(0,1.5));}
      else { unstables = apply<ParticleFinder>(event, "DStars").particles(Cuts::pT > 2.0*GeV && Cuts::absrapIn(0,1.5));}
      const Particles dstars = filter_select(unstables, [](const Particle& p){return p.abspid() == PID::DSPLUS;});
      
      if(dstars.empty() ) vetoEvent;
      MSG_DEBUG("D*" << dstars.size());
      
      const Particle& dstar = dstars.front();
      // boosting the system
      const LorentzTransform hcmboost = dk.boostHCM();
      const FourMomentum hcmMom = hcmboost.transform(dstar.momentum());

      //discriminate between dis and photoprod, and between ETA33 and ETA 44

      //kinematics quantities
      const double E = dstar.E();
      const double p_z = dstar.pz();
      //std::cout<<"y: "<<y<<endl;
      if(y<0.02) vetoEvent;
      const double m2 = 2.25; // charm mass^2
      const double E_e = dk.beamLepton().E();
      const double z = (E - p_z)/(2*y*E_e);
      
     
      const double M2 = (1.44*hcmMom.pT2() + m2)/(z*(1-z));
      const double x_g = (M2 + Q2)/(y*dk.s());
      // std::cout<<"s:  "<<dk.s()<<endl;
      // std::cout<<"M2: "<<M2<<endl;
      //  std::cout<<"z: "<<z<<endl;
      
      const double y_capp = dstar.rapidity();
      const double W = sqrt(dk.W2());
      //std::cout<<"x_g: "<<x_g<<endl;
      
      //perform the cuts
      if(isDIS == true){
         _h["211"]->fill(dstar.pT());
         _h["411"]->fill(dstar.eta());
         _h["511"]->fill(Q2);
         _h["611"]->fill(log10(x_g));
		//boosting to the hcm frame
         _h["311"]->fill(hcmMom.pT());
      }
      if(ETAG33 == true && abs(y_capp) < 1.5 && dstar.pT() > 2.5*GeV  ) {
		
         _h["rap194"]->fill(y_capp);
         _h["pt194"] ->fill(dstar.pT());
         _hpt.fill(dstar.pT(), y_capp);	
	
	if( W > 173 &&  W<273){
	   _h["1211"]->fill(log10(x_g));
	  
	}
	if( W > 130 &&  W<230){
	   _h["1212"]->fill(log10(x_g));
	  
	}
	
	
      }
      if(ETAG44 == true && abs(y_capp) < 1.5 && dstar.pT()> 2) { 

	   _h["pt88"]->fill(dstar.pT());
	   _h["rap88"]->fill(y_capp);
	   _h["1311"]->fill(log10(x_g));
	
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      // conversion factors from ep to gamma p xsections (as given in publication)
      const double F_etag33 = 0.0128;
      const double F_etag44 = 0.0838;

      double norm = crossSection()/nanobarn/sumW() ;
      scale( _h["211"], norm);
      
      scale(_h["311"], norm);  
      scale(_h["411"], norm);
      scale(_h["511"], norm);
      scale(_h["611"], norm);
      double norm_mub = crossSection()/microbarn/sumW();

      scale(_h["rap194"],norm_mub/F_etag33 );
      scale(_h["pt194"], norm_mub/F_etag33 );
      scale(_h["rap88"], norm_mub/F_etag44 );
      scale(_h["pt88"],  norm_mub/F_etag44 );
         
      _hpt.scale( norm/F_etag33, this);
      scale(_h["1211"], norm_mub/F_etag33);
      scale(_h["1212"], norm_mub/F_etag33);
      scale(_h["1311"], norm_mub/F_etag44);
    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    BinnedHistogram  _hpt; 
    ///@}


  };


  RIVET_DECLARE_PLUGIN(H1_1999_I481112);

}
