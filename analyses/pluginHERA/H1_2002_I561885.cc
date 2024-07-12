// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Measurement of D*+- meson production in deep inelastic scattering at HERA (H1)
  class H1_2002_I561885 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(H1_2002_I561885);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
    
      
      declare(DISKinematics(), "Kinematics");
      declare(UnstableParticles(), "Dstars");
      //Cuts::abspid == PID::DSTARPLUS

      // Initialise and register projections

      
      Histo1DPtr dummy; //Introducing


      // Book histograms
        book(_h["W_GeV"], 2, 1, 1);
        book(_h["p_tD*"], 3, 1, 1);
        book(_h["logx"], 4, 1, 1);
        book(_h["etaD*"], 5, 1, 1);
        book(_h["Q2"], 6, 1, 1);
        book(_h["Z_D*"], 7, 1, 1);
        
        _h_Q2eta.add(-1.5,-0.5, book(dummy,8,1,1)); 
        _h_Q2eta.add(-0.5,0.5, book(dummy,8,1,2)); 
        _h_Q2eta.add(0.5,1.5, book(dummy,8,1,3)); 


        _h_Q2pt.add(1.5,4., book(dummy,9,1,1)); 
        _h_Q2pt.add(4.,10., book(dummy,9,1,2)); 


        _h_eta1.add(0,0.25, book(dummy,10,1,1));
        _h_eta1.add(0.25,0.50, book(dummy,10,1,2));
        _h_eta1.add(0.5,1, book(dummy,10,1,3));
         
         
        _h_eta2.add(1.5,2.5, book(dummy,11,1,1));
        _h_eta2.add(2.5,4.0, book(dummy,11,1,2));
        _h_eta2.add(4.0,10.0,book(dummy,11,1,3));
         
         
        _h_zD1.add(1.5,2.5, book(dummy,12,1,1));
        _h_zD1.add(2.5,4.0, book(dummy,12,1,2));
        _h_zD1.add(4.0,10.0, book(dummy,12,1,3));
                 
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
          
      const DISKinematics& kin = applyProjection<DISKinematics>(event, "Kinematics");
      
       // Q2 and inelasticity cuts
      if (!inRange(kin.Q2(), 1.0*GeV2, 100*GeV2)) vetoEvent;
      if (!inRange(kin.y(), 0.05, 0.7)) vetoEvent;
      
        // D* reconstruction
      // const Particles unstables = apply<ParticleFinder>(event, "Dstars").particles(Cuts::pT > 1.5*GeV && Cuts::abseta < 1.5);
      const Particles unstables = apply<ParticleFinder>(event, "Dstars").particles(Cuts::pT > 1.5*GeV );
      const Particles dstars = filter_select(unstables, [](const Particle& p){ return p.abspid() == PID::DSTARPLUS; });
      if (dstars.empty()) vetoEvent;
      // MSG_DEBUG("#D* = " << dstars.size());
      //const Particle& dstar = dstars.front();
      for(const Particle& dstar : dstars){
      const double zD = (dstar.E() - dstar.pz()) / (2*kin.beamLepton().E()*kin.y());
  
  
  // Single-differential histograms
    //     cout << " ' D* found "<< dstar.pid() << " " <<  kin.Q2() << " " << kin.x() << " " << dstar.pT() << " " << dstar.eta() << " " << dstar.rapidity() << endl;
        _h["p_tD*"]->fill(dstar.pT()/GeV);
        _h["etaD*"]->fill(dstar.eta());
        _h["Z_D*"]->fill(zD/GeV);
        _h["Q2"]->fill(kin.Q2()/GeV2);
        _h["W_GeV"]->fill(sqrt(kin.W2())/GeV);
        _h["logx"]->fill(log10(kin.x()));
  // Double-differential (y,Q2) histograms

        _h_Q2eta.fill(dstar.eta(),kin.Q2());
        _h_Q2pt.fill(dstar.pT(),kin.Q2());

        _h_eta1.fill(zD, dstar.eta());
        _h_eta2.fill(dstar.pT(), dstar.eta());
        _h_zD1.fill(dstar.pT(), zD/GeV);

        }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
    

     const double sf = crossSection()/nanobarn/sumOfWeights();
      scale(_h["p_tD*"], sf);
      scale(_h["etaD*"], sf);
      scale(_h["Z_D*"], sf);
      scale(_h["Q2"], sf);
      scale(_h["W_GeV"], sf);
      scale(_h["logx"], sf);
      
      _h_Q2eta.scale(sf, this );
      _h_Q2pt.scale(sf, this );
      _h_eta1.scale(sf, this );
      _h_eta2.scale(sf, this );
      _h_zD1.scale(sf, this );
        

    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    BinnedHistogram _h_eta1,_h_eta2,_h_zD1, _h_Q2eta, _h_Q2pt;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(H1_2002_I561885);

}
