// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/DISLepton.hh"

namespace Rivet {


  /// @brief A Study of the Fragmentation of Quarks in ep Collisions at HERA (H1)
  class H1_1995_I394793 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(H1_1995_I394793);


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
      const FinalState fs(Cuts::abseta < 4.9);
      declare(fs, "FS");
      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      
      // book a counter
      book(_Nevt_after_cuts, "TMP/Nevt_after_cuts");
      book(_Nevt_afterfwd_cuts, "TMP/Nevt_afterfwd_cuts");
      book(_Nevt_afterh_cuts, "TMP/Nevt_afterh_cuts");
      book(_Nevt_afterhfwd_cuts, "TMP/Nevt_afterhfwd_cuts");


      // Book histograms
      // specify custom binning
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      book(_h["costh_lowQ"], 1, 1, 1);
      book(_h["costh_highQ"], 1, 1, 2);     
      book(_h["costh_lowQ_noEfwd"], 2, 1, 1);
      book(_h["costh_highQ_noEfwd"], 2, 1, 2);
      book(_h["xp_posCharge_lowQ"], 3, 1, 1);
      book(_h["xp_negCharge_lowQ"], 3, 1, 2);
      book(_h["xp_posCharge_highQ"], 3, 1, 3);
      book(_h["xp_negCharge_highQ"], 3, 1, 4);
      book(_h["ksi_lowQ"], 4, 1, 1);
      book(_h["ksi_highQ"], 4, 1, 2);
      book(_s["Mult_vrs_Q2"], 5, 1, 1);
      book(_s["Mult_vrs_Q2_noEfwd"], 6, 1, 1);
      
      book(_h["Mult_vrs_Q2_nchrg"],"TMP/Mult_vrs_Q2_nchrg", refData(5,1,1));
      book(_h["Mult_vrs_Q2_noEfwd_nchrg"],"TMP/Mult_vrs_Q2_noEfwd_nchrg", refData(6,1,1));
      book(_h["Mult_vrs_Q2_count"],"TMP/Mult_vrs_Q2_count", refData(5,1,1));
      book(_h["Mult_vrs_Q2_noEfwd_count"],"TMP/Mult_vrs_Q2_noEfwd_count", refData(6,1,1));

    }
 

    /// Perform the per-event analysis
    void analyze(const Event& event) {

     const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
    
   
    //DIS kinematics
      const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");
      if ( dk.failed() ) vetoEvent;
      double y   = dk.y();
      double w2  = dk.W2();
      double Q2  = dk.Q2();
      
      
      
      bool cut ;
      
      cut = Q2 > 12 && y < 0.6 && w2 > 3000 ;
      
      
      if ( !cut ) vetoEvent ; 
      
      const DISLepton& dl = applyProjection<DISLepton>(event,"Lepton");
      if ( dl.failed() ) vetoEvent;
      /*
      cout << "  scattered lepton angle " << 180.- dl.out().momentum().angle(dl.in().momentum())/degree << endl;
      cout << " in  lepton " << dl.in().momentum() << endl;
      cout << " out lepton " << dl.out().momentum() << endl;
      */
      
      const FinalState& fs = apply<FinalState>(event, "FS");
      Particles particles; particles.reserve(fs.size());
      ConstGenParticlePtr dislepGP = dl.out().genParticle();
      for (const Particle& p : cfs.particles()) {
        ConstGenParticlePtr loopGP = p.genParticle();
        if (loopGP == dislepGP) continue;
        particles.push_back(p);
      }


      double efwd = 0.;

      for (const Particle& p : particles) {
        const double th = 180. - p.momentum().angle(dl.in().momentum())/degree;

        if (inRange(th, 4.4, 15.0)) { 
           efwd += p.E();
           //cout << " angle " << th << " pid " << p.pid() << " Efwd = " << efwd << endl;
           } 
      }

      bool evcut[2];
      evcut[0] =   efwd > 0.5;
      
       // fill the counter
      _Nevt_after_cuts -> fill();
      if (Q2 > 100 )  _Nevt_afterh_cuts -> fill();
       
      if ( evcut[0] && (Q2 < 80 ) ) _Nevt_afterfwd_cuts -> fill();
      if ( evcut[0] && (Q2 > 100 ) ) _Nevt_afterhfwd_cuts -> fill();

      double n_charg = 0; 

      // Boost to Breit
      const LorentzTransform breitboost = dk.boostBreit();
       
      for (size_t ip1 = 0; ip1 < particles.size(); ++ip1) {
        const Particle& p = particles[ip1];
        const FourMomentum BreMom = breitboost.transform(p.momentum());
       // cout << BreMom.pz() << endl;
        double x = cos(BreMom.theta());
        
        if (Q2 < 80 ) _h["costh_lowQ"] ->fill(x);
        if (Q2 > 100 ) _h["costh_highQ"]  ->fill(x);
        
        if (Q2 < 80  && evcut[0]) _h["costh_lowQ_noEfwd"] ->fill(x);
        if (Q2 > 100 && evcut[0]) _h["costh_highQ_noEfwd"]  ->fill(x);


        if ( BreMom.pz() > 0. ) continue;
        double pcal= sqrt(BreMom.px2() + BreMom.py2()+ BreMom.pz2()) ;
        double xp = 2*pcal/(sqrt(Q2));
        double xi = log(1/xp);
        
        double charge = p.charge() ;
       // cout << " charge " << charge << endl;
       
        if (charge > 0 ) {
          if (Q2 < 80 ) _h["xp_posCharge_lowQ"] -> fill(xp);
          if (Q2 > 100 ) _h["xp_posCharge_highQ"] -> fill(xp);
        } else {
         if (Q2 < 80 ) _h["xp_negCharge_lowQ"] -> fill(xp);
         if (Q2 > 100 ) _h["xp_negCharge_highQ"] -> fill(xp);
        }
    
        if (Q2 < 80 ) _h["ksi_lowQ"] -> fill(xi);
        if (Q2 > 100 ) _h["ksi_highQ"] -> fill (xi);
 
        n_charg = n_charg + 1; 

     }
     _h["Mult_vrs_Q2_nchrg"] -> fill(Q2,n_charg) ;
     _h["Mult_vrs_Q2_count"] -> fill(Q2) ;
     if ( evcut[0]) _h["Mult_vrs_Q2_noEfwd_nchrg"] -> fill(Q2,n_charg) ;
     if ( evcut[0]) _h["Mult_vrs_Q2_noEfwd_count"] -> fill(Q2) ;
}

    /// Normalise histograms etc., after the run
    void finalize() {
     
   
      normalize(_h["xp_posCharge_lowQ"]);
      if(dbl(*_Nevt_afterh_cuts)>0) scale(_h["xp_posCharge_highQ"], 1.0/ *_Nevt_afterh_cuts);

      normalize(_h["xp_negCharge_lowQ"]);
      if(dbl(*_Nevt_afterh_cuts)>0) scale(_h["xp_negCharge_highQ"], 1.0/ *_Nevt_afterh_cuts);

      scale(_h["costh_lowQ"], 1.0/ *_Nevt_after_cuts);
      if(dbl(*_Nevt_afterh_cuts)>0) scale(_h["costh_highQ"], 1.0/ *_Nevt_afterh_cuts);
      //cout << " after fwd cuts " << dbl(*_Nevt_afterfwd_cuts) << endl;
      
      if(dbl(*_Nevt_afterfwd_cuts)>0) scale(_h["costh_lowQ_noEfwd"], 1.0/ *_Nevt_afterfwd_cuts);
      if(dbl(*_Nevt_afterhfwd_cuts)>0) scale(_h["costh_highQ_noEfwd"], 1.0/ *_Nevt_afterhfwd_cuts);
      
      scale(_h["ksi_lowQ"], 1.0/ *_Nevt_after_cuts);
      if(dbl(*_Nevt_afterh_cuts)>0) scale(_h["ksi_highQ"], 1.0/ *_Nevt_afterh_cuts);
      //cout << " Nevt " << dbl(*_Nevt_after_cuts) << endl;

      divide(_h["Mult_vrs_Q2_nchrg"], _h["Mult_vrs_Q2_count"], _s["Mult_vrs_Q2"]);
      divide(_h["Mult_vrs_Q2_noEfwd_nchrg"], _h["Mult_vrs_Q2_noEfwd_count"], _s["Mult_vrs_Q2_noEfwd"]);
    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    map<string, Scatter2DPtr> _s;
    CounterPtr _Nevt_after_cuts;
    CounterPtr _Nevt_afterfwd_cuts;
    CounterPtr _Nevt_afterh_cuts;
    CounterPtr _Nevt_afterhfwd_cuts;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(H1_1995_I394793);

}
