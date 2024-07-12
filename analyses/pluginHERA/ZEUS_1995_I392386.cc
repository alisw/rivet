// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"


namespace Rivet {

const vector<double> QEdges {10., 20., 40., 80., 160., 320.};
const vector<double> xEdges {0.6e-3,1.2e-3,2.4e-3,1.0e-2,5.0e-2};

  /// @brief Charged particle multiplicity and momentum spectra in Breit frame at ZEUS
  class ZEUS_1995_I392386 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ZEUS_1995_I392386);


    /// @name Analysis methods
    ///@{
      
    /// Book histograms and initialise projections before the run
    void init() {
        
      // Initialise and register projections
        declare(DISKinematics(), "Kinematics");
        
      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
        const ChargedFinalState fs;
        declare(fs, "FS");


      // Book histograms
      // specify custom binning
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      
        for(int iQ = 0; iQ < 11; ++iQ) {
          book(_Nevt_after_cuts_Q[iQ], "TMP/Nevt_after_cuts_Q"+ to_string(iQ));     
        }

        
        book(_h["mult1"],"TMP/mult 1", refData(1, 1, 1)); // Multiplicity
        book(_h["mult1_Q"],"TMP/mult 1_1", refData(1, 1, 1));
        book(_h["mult2"],"TMP/mult 2", refData(2, 1, 1));
        book(_h["mult2_Q"],"TMP/mult 2_1", refData(2, 1, 1));
        book(_h["mult3"], "TMP/mult 3", refData(3, 1, 1));
        book(_h["mult3_Q"], "TMP/mult 3_1",refData(3, 1, 1));
        book(_h["mult4"], "TMP/mult 4", refData(4, 1, 1));
        book(_h["mult4_Q"],"TMP/mult 4_1", refData(4, 1, 1));
        
        book(_h["mom1"], "TMP/mult 5", refData(5, 1, 1)); // Momentum spectra
        book(_h["mom1_Q"],"TMP/mult 5_1", refData(5, 1, 1));
        book(_h["mom2"], "TMP/mult 6", refData(6, 1, 1));
        book(_h["mom2_Q"],"TMP/mult 6_1", refData(6, 1, 1));
        book(_h["mom3"], "TMP/mult 7", refData(7, 1, 1));
        book(_h["mom3_Q"],"TMP/mult 7_1", refData(7, 1, 1));
        book(_h["mom4"], "TMP/mult 8", refData(8, 1, 1));
        book(_h["mom4_Q"], "TMP/mult 8_1", refData(8, 1, 1));
        
        book(_h_mult1, 1,1,1);
        book(_h_mult2, 2,1,1);
        book(_h_mult3, 3,1,1);
        book(_h_mult4, 4,1,1);
        book(_h_mom1, 5,1,1);
        book(_h_mom2, 6,1,1);
        book(_h_mom3, 7,1,1);
        book(_h_mom4, 8,1,1);
        
        book(_h["nch1"], 9, 1, 1);  // Multiplicity
        book(_h["nch2"], 10, 1, 1);
        book(_h["nch3"], 10, 1, 2);
        book(_h["nch4"], 10, 1, 3);
        book(_h["nch5"], 11, 1, 1);
        book(_h["nch6"], 11, 1, 2);
        book(_h["nch7"], 11, 1, 3);
        book(_h["nch8"], 11, 1, 4);
        book(_h["nch9"], 12, 1, 1);
        book(_h["nch10"], 12, 1, 2);
        
        book(_h["loginvmom1"], 13, 1, 1); // Momentum spectra
        book(_h["loginvmom2"], 14, 1, 1);
        book(_h["loginvmom3"], 14, 1, 2);
        book(_h["loginvmom4"], 14, 1, 3);
        book(_h["loginvmom5"], 15, 1, 1);
        book(_h["loginvmom6"], 15, 1, 2);
        book(_h["loginvmom7"], 15, 1, 3);
        book(_h["loginvmom8"], 15, 1, 4);
        book(_h["loginvmom9"], 16, 1, 1);
        book(_h["loginvmom10"], 16, 1, 2);
        
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
        const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "FS");
        const DISKinematics& dk = applyProjection<DISKinematics>(event, "Kinematics");

        double xbj = dk.x(); // momentum fraction
        double Q2 = dk.Q2()/GeV; // momentum transfer
        const LorentzTransform Breitboost = dk.boostBreit();
        
        // Multiplicity counters
        int n911(0), n1011(0), n1012(0), n1013(0), n1111(0), n1112(0), n1113(0), n1114(0), n1211(0), n1212(0);


        if(0.6e-3<xbj && xbj<1.2e-3) {
            if(10<Q2 && Q2<20) {
                    _Nevt_after_cuts_Q[1] -> fill();
            }
            if(10<Q2 && Q2<20) {
                  _Nevt_after_cuts_Q[2] -> fill();
            }
            if(20<Q2 && Q2<40) {
                  _Nevt_after_cuts_Q[3] -> fill();
            }
            if(40<Q2 && Q2<80) {
                  _Nevt_after_cuts_Q[4] -> fill();
            }
        }
        if(2.4e-3<xbj && xbj<1.0e-2) {
            if(20<Q2 && Q2<40) {
                  _Nevt_after_cuts_Q[5] -> fill();
            }
            if(40<Q2 && Q2<80){
                  _Nevt_after_cuts_Q[6] -> fill();
            }
            if(80<Q2 && Q2<160) {
                  _Nevt_after_cuts_Q[7] -> fill();
            }
            if(160<Q2 && Q2<320) {
                 _Nevt_after_cuts_Q[8] -> fill();
            }
        }
        if(1.0e-2<xbj && xbj<5.0e-2) {
            if(320<Q2 && Q2<640) {
                 _Nevt_after_cuts_Q[9] -> fill();
            }
            if(640<Q2 && Q2<1280) {
                 _Nevt_after_cuts_Q[10] -> fill();
            }
        }
            

        for (const Particle& p : cfs.particles()) {
            //??? calculating ln(1/x_p) ??? --> part to ask
            const FourMomentum BrMom = Breitboost.transform(p.momentum());
            double pp = sqrt(BrMom.px()*BrMom.px() + BrMom.py()*BrMom.py() + BrMom.pz()*BrMom.pz() );
            double xp = 2*pp/(sqrt(Q2));
            const double logInvScaledMom = log(1/xp);
            
            if ( BrMom.pz() > 0. ) continue;
            
            if(0.6e-3<xbj && xbj<1.2e-3) {
                _h["mom1"] ->fill(Q2, logInvScaledMom);
                _h["mom1_Q"] ->fill(Q2);
                if(10<Q2 && Q2<20) {
                    _h["loginvmom1"] ->fill(logInvScaledMom);
                    ++n911;
                }
            }
            
          if(1.2e-3<xbj && xbj<2.4e-3) {
                _h["mom2"] ->fill(Q2, logInvScaledMom);
                _h["mom2_Q"] ->fill(Q2);
              if(10<Q2 && Q2<20) {
                  _h["loginvmom2"] ->fill(logInvScaledMom);
                  ++n1011;
                  }
              if(20<Q2 && Q2<40) {
                  _h["loginvmom3"] ->fill(logInvScaledMom);
                  ++n1012;
                  }
              if(40<Q2 && Q2<80) {
                  _h["loginvmom4"] ->fill(logInvScaledMom);
                  ++n1013;
                  }
            }
            
            if(2.4e-3<xbj && xbj<1.0e-2) {
                _h["mom3"] ->fill(Q2, logInvScaledMom);
                _h["mom3_Q"] ->fill(Q2);
                if(20<Q2 && Q2<40) {
                    _h["loginvmom5"] ->fill(logInvScaledMom);
                    ++n1111;
                    }
                if(40<Q2 && Q2<80){
                    _h["loginvmom6"] ->fill(logInvScaledMom);
                    ++n1112;
                    }
                if(80<Q2 && Q2<160) {
                    _h["loginvmom7"] ->fill(logInvScaledMom);
                    ++n1113;
                    }
                if(160<Q2 && Q2<320) {
                    _h["loginvmom8"] ->fill(logInvScaledMom);
                    ++n1114;
                    }
            }
            
            if(1.0e-2<xbj && xbj<5.0e-2) {
                _h["mom4"] ->fill(Q2,logInvScaledMom);
                _h["mom4_Q"] ->fill(Q2);
                if(320<Q2 && Q2<640) {
                    _h["loginvmom9"] ->fill(logInvScaledMom);
                    ++n1211;
                    }
                    if(640<Q2 && Q2<1280) {
                    _h["loginvmom10"] ->fill(logInvScaledMom);
                    ++n1212;
                    }
            }
            
            }
        
        if(0.6e-3<xbj && xbj<1.2e-3) {
            if(10<Q2 && Q2<20) {
                _h["mult1"] ->fill(Q2, n911);
                _h["mult1_Q"] ->fill(Q2);
                _h["nch1"] ->fill(n911);
            }
        }
        
      if(1.2e-3<xbj && xbj<2.4e-3) {
          if(10<Q2 && Q2<80) {
            _h["mult2"] ->fill(Q2, n1011+n1012+n1013);
            _h["mult2_Q"] ->fill(Q2);
              if(10<Q2 && Q2<20) {
                  _h["nch2"] ->fill(n1011); }
              if(20<Q2 && Q2<40) {
                  _h["nch3"] ->fill(n1012);}
              if(40<Q2 && Q2<80) {
                  _h["nch4"] ->fill(n1013);}
          }
        }
        
        if(2.4e-3<xbj && xbj<1.0e-2) {
            _h["mult3"] ->fill(Q2,n1111+n1112+n1113+n1114);
            _h["mult3_Q"] ->fill(Q2);
            if(20<Q2 && Q2<40) {
                _h["nch5"] ->fill(n1111);}
            if(40<Q2 && Q2<80) {
                _h["nch6"] ->fill(n1112);}
            if(80<Q2 && Q2<160) {
                _h["nch7"] ->fill(n1113);}
            if(160<Q2 && Q2<320) {
                _h["nch8"] ->fill(n1114);}

        }
        
        if(1.0e-2<xbj && xbj<5.0e-2) {
            _h["mult4"] ->fill(Q2, n1211+n1212);
            _h["mult4_Q"] ->fill(Q2);
            if(320<Q2 && Q2<640) {
                _h["nch9"] ->fill(n1211);}
            if(640<Q2 && Q2<1280) {
                _h["nch10"] ->fill(n1212); }
        }
        
    }

      
    /// Normalise histograms etc., after the run
    void finalize() {
   //   normalize(_h["XXXX"]); // normalize to unity
        divide(_h["mult1"],_h["mult1_Q"],_h_mult1); // Q versus Multiplicity
        divide(_h["mult2"],_h["mult2_Q"],_h_mult2);
        divide(_h["mult3"],_h["mult3_Q"],_h_mult3);
        divide(_h["mult4"],_h["mult4_Q"],_h_mult4);
        divide(_h["mom1"],_h["mom1_Q"],_h_mom1); // Q versus Momentum
        divide(_h["mom2"],_h["mom2_Q"],_h_mom2);
        divide(_h["mom3"],_h["mom3_Q"],_h_mom3);
        divide(_h["mom4"],_h["mom4_Q"],_h_mom4);
        
        
        normalize(_h["nch1"]); //multiplicity
        normalize(_h["nch2"]);
        normalize(_h["nch3"]);
        normalize(_h["nch4"]);
        normalize(_h["nch5"]);
        normalize(_h["nch6"]);
        normalize(_h["nch7"]);
        normalize(_h["nch8"]);
        normalize(_h["nch9"]);
        normalize(_h["nch10"]);

        if(dbl(*_Nevt_after_cuts_Q[1])>0 ) scale(_h["loginvmom1"],1./ *_Nevt_after_cuts_Q[1]);  //momentum
        if(dbl(*_Nevt_after_cuts_Q[2])>0 ) scale(_h["loginvmom2"],1./ *_Nevt_after_cuts_Q[2]);
        if(dbl(*_Nevt_after_cuts_Q[3])>0 ) scale(_h["loginvmom3"],1./ *_Nevt_after_cuts_Q[3] );
        if(dbl(*_Nevt_after_cuts_Q[4])>0 ) scale(_h["loginvmom4"],1./ *_Nevt_after_cuts_Q[4] );
        if(dbl(*_Nevt_after_cuts_Q[5])>0 ) scale(_h["loginvmom5"],1./ *_Nevt_after_cuts_Q[5] );
        if(dbl(*_Nevt_after_cuts_Q[6])>0 ) scale(_h["loginvmom6"],1./ *_Nevt_after_cuts_Q[6]);
        if(dbl(*_Nevt_after_cuts_Q[7])>0 ) scale(_h["loginvmom7"],1./ *_Nevt_after_cuts_Q[7]);
        if(dbl(*_Nevt_after_cuts_Q[8])>0 ) scale(_h["loginvmom8"],1./ *_Nevt_after_cuts_Q[8]);
        if(dbl(*_Nevt_after_cuts_Q[9])>0 ) scale(_h["loginvmom9"],1./ *_Nevt_after_cuts_Q[9] );
        if(dbl(*_Nevt_after_cuts_Q[10])>0 ) scale(_h["loginvmom10"],1./ *_Nevt_after_cuts_Q[10] );


        
     // scale(_h["ZZZZ"], crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts)

    }

    ///@}
      

    /// @name Histograms
    ///@{
      
      Scatter2DPtr _h_mult1;
      Scatter2DPtr _h_mult2;
      Scatter2DPtr _h_mult3;
      Scatter2DPtr _h_mult4;
      Scatter2DPtr _h_mom1;
      Scatter2DPtr _h_mom2;
      Scatter2DPtr _h_mom3;
      Scatter2DPtr _h_mom4;
      
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    CounterPtr _Nevt_after_cuts_Q[11];
    BinnedHistogram _h_invmom1,_h_invmom2,_h_invmom3;

    ///@}


  };


  RIVET_DECLARE_PLUGIN(ZEUS_1995_I392386);

}
