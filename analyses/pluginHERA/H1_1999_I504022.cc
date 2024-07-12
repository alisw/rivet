// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

#include "Rivet/Projections/DISKinematics.hh"

namespace Rivet {


  /// @brief Forward pi0 meson production at HERA (H1)
  class H1_1999_I504022 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(H1_1999_I504022);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // The basic final-state projection:
 
      declare(FinalState(Cuts::abseta < 7 ), "FS");
      declare(DISKinematics(), "Kinematics");
      declare(UnstableParticles(), "UFS");


      // Book histograms
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      
      book(_h["p_T>2.5&x"], 1, 1, 1);
      book(_h["Q:2.0-4.5-eta"], 2, 1, 1);
      book(_h["Q:2.0-4.5-p_T"], 3, 1, 1);
      book(_h["Q:4.5-15.0-x"], 4, 1, 1);
      book(_h["Q:4.5-15.0-eta"], 5, 1, 1);
      book(_h["Q:4.5-15.0-p_T"], 6, 1, 1);
      book(_h["Q:15.0-70.0-x"], 7, 1, 1);
      book(_h["Q:15.0-70.0-eta"], 8, 1, 1);
      book(_h["Q:15.0-70.0-p_T"], 9, 1, 1);
      book(_h["p_T>2.5&Q"], 10, 1, 1);
      book(_h["p_T>3.5&x"], 11, 1, 1);
      book(_h["p_T>3.5&Q"], 12, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      /// @todo Do the event by event analysis here

        const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");

        // Get the DIS kinematics
        double xbj  = dk.x();
        double ybj = dk.y();
        double Q2 = dk.Q2()/GeV2;
        
        // Q2 and inelasticity cuts
        //cout << " after xbj " << xbj << endl;
        if (!inRange(ybj, 0.1, 0.6)) vetoEvent;
        //cout << " after ybj " << ybj << endl;
        if (!inRange(Q2, 2.0*GeV2, 70.0*GeV2)) vetoEvent;
        //cout << " after Q2 " << Q2 << endl;

      const FinalState& fs = apply<FinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();
      //cout << " Num all     final state particles " << numParticles << endl;
      //proton beam energy: 820 GeV
      double e_proton = dk.beamHadron().E()/GeV;
      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }

      // Extracting the pi0 
      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");
      //Get the hadronic CMS kinematics
      const LorentzTransform hcmboost = dk.boostHCM();

      for (const Particle& p : ufs.particles(Cuts::pid==PID::PI0)) {
        //cout << " Pid = " << p.pid() << endl;
        //Get the LAB kinematics
        double theta = p.theta();
        double eta = p.momentum().pseudorapidity();
        // double pT = p.momentum().pT()/GeV; 
        //Boost hcm
        const FourMomentum hcmMom = hcmboost.transform(p.momentum());
        double pThcm =hcmMom.pT();

        double e_pi0 = p.E()/GeV;
        double x_pi0_proton = e_pi0/e_proton;

        // epi0/e_proton, theta and pThcm cuts
        if(x_pi0_proton < 0.01) continue;
        //cout << " theta " << theta << " in deg " << theta/degree << endl;
        if (!inRange(theta/degree, 5, 25)) continue;
        if(pThcm < 2.5*GeV) continue;
        
        //Three cuts for Q2:
        if (Q2 > 2.0*GeV2 && Q2 < 4.5*GeV2){
          _h["Q:2.0-4.5-eta"]->fill(eta);
          _h["Q:2.0-4.5-p_T"]->fill(pThcm);
        } 
        if (Q2 > 4.5*GeV2 && Q2 < 15.0*GeV2){
          _h["Q:4.5-15.0-x"]->fill(xbj);
          _h["Q:4.5-15.0-eta"]->fill(eta);
          _h["Q:4.5-15.0-p_T"]->fill(pThcm);
        }  
        if (Q2 > 15.0*GeV2 && Q2 < 70.0*GeV2){
          _h["Q:15.0-70.0-x"]->fill(xbj);
          _h["Q:15.0-70.0-eta"]->fill(eta); 
          _h["Q:15.0-70.0-p_T"]->fill(pThcm);
        }

        //Two cuts for p_T:
        if (pThcm > 2.5*GeV){
          _h["p_T>2.5&x"]->fill(xbj);
          _h["p_T>2.5&Q"]->fill(Q2);
        }
        if (pThcm > 3.5*GeV){
          _h["p_T>3.5&x"]->fill(xbj);
          _h["p_T>3.5&Q"]->fill(Q2);
        }
        
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double sn = crossSection()/nanobarn/sumW();
      const double sp = crossSection()/picobarn/sumW();
      scale(_h["p_T>2.5&x"], sn);
      scale(_h["Q:2.0-4.5-eta"], sp);
      scale(_h["Q:2.0-4.5-p_T"], sp);
      scale(_h["Q:4.5-15.0-x"], sn);
      scale(_h["Q:4.5-15.0-eta"], sp);
      scale(_h["Q:4.5-15.0-p_T"], sp);
      scale(_h["Q:15.0-70.0-x"], sn);
      scale(_h["Q:15.0-70.0-eta"], sp);
      scale(_h["Q:15.0-70.0-p_T"], sp);
      scale(_h["p_T>2.5&Q"], sp);
      scale(_h["p_T>3.5&x"], sn);
      scale(_h["p_T>3.5&Q"], sp);
    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(H1_1999_I504022);

}
