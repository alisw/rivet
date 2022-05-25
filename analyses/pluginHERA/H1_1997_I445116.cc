// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/DISLepton.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief Evolution of e p fragmentation and multiplicity distributions in the Breit frame (H1)
  class H1_1997_I445116 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(H1_1997_I445116);


    /// @name Analysis methods
    /// @{
    const vector<double> QEdges {3.17, 3.915, 4.72, 6.13, 7.635, 8.85, 10.81, 13.415, 16.365, 21.745, 33.835, 50.745};
    const vector<double> AvgEdges {3.13, 3.875, 4.675, 6.065, 7.55, 8.76, 9.785, 14.675, 16.25, 21.45, 33.06, 49.26};
    const vector<double> xp_range{0.02, 0.05, 0.10, 0.20, 0.3, 0.4, 0.5, 0.7};
    const size_t iPmax = 7;

    /// Book histograms and initialise projections before the run
    void init() {

      declare(DISKinematics(), "Kinematics");
      const DISLepton dl;
      declare(ChargedFinalState(dl.remainingFinalState()), "CFS");

      book(_Nevt_after_cuts_Qlow, "TMP/Nevt_after_cuts_Qlow");
      book(_Nevt_after_cuts_QHigh, "TMP/Nevt_after_cuts_QHigh");

      book(_h["xp"], 1, 1, 1);
      book(_h["xpQgt"], 1, 1, 2);
      book(_h["xi"], 2, 1, 1);
      book(_h["xiQgt"], 2, 1, 2);

      Histo1DPtr dummy;
      for (size_t iQ = 0; iQ < 11; ++iQ) {
        _h_Q2_xp.add(QEdges[iQ], QEdges[iQ+1],book(dummy,"TMP/xpQ"+ to_string(iQ), xp_range));
        book(_Nevt_after_cuts_Q[iQ], "TMP/Nevt_after_cuts_Q"+ to_string(iQ));
      }

      for (size_t iP = 0 ; iP < iPmax; ++iP) {
        book(_s["Qxp"+to_string(iP)], 3+iP, 1, 1);
      }

      book(_s["Avg1"], 10,1,1);
      book(_s["Avg2"], 11,1,1);

      book(_h["QT1"], "TMP/QT1",refData(10,1,1));
      book(_h["QT2"], "TMP/QT2",refData(11,1,1));

      book(_h["AvgT1"], "TMP/AvgT1",refData(10,1,1));
      book(_h["AvgT2"], "TMP/AvgT2",refData(11,1,1));

      for (size_t iE = 0; iE < 8; ++iE) {
        book(_h["E"+to_string(iE)], 12+iE, 1, 1);
        book(_Nevt_after_cuts_E[iE], "TMP/Nevt_after_cuts_E"+to_string(iE));
      }

      for (size_t iN = 0; iN < 6; ++iN) {
        book(_h["N"+to_string(iN)], 20+iN, 1, 1);
        book(_Nevt_after_cuts_N[iN], "TMP/Nevt_after_cuts_N"+to_string(iN));
      }

      book(_h["MeanTest1"], "TMP/Meantest1",20,5.,6.14);
      book(_h["MeanTest2"], "TMP/Meantest2",50,19.,8000.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      Particles particles = cfs.particles();
      const size_t numParticles = particles.size();

      const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");
      double Q2 = dk.Q2()/GeV2;
      double y= dk.y();
      double x= dk.x();
      double W2= dk.W2();
      double Q = sqrt(Q2);

      if (y < 0.05 or y > 0.6 ) vetoEvent;
      if (W2 < 4400) vetoEvent ;

      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }

      for (int iQ = 0; iQ < 11; ++iQ) {
        if (inRange(sqrt(Q2), QEdges[iQ],QEdges[iQ+1])) {
          _Nevt_after_cuts_Q[iQ] -> fill();
        }
      }

      if (Q > 5 && Q < 6.14)  _Nevt_after_cuts_E[0]->fill();
      if (Q >19 && Q < 8000)  _Nevt_after_cuts_E[1]->fill();
      if (Q2 >12 && Q2 <  15) _Nevt_after_cuts_E[2]->fill();
      if (Q2 >15 && Q2 <  25) _Nevt_after_cuts_E[3]->fill();
      if (Q2 >20 && Q2 <  40) _Nevt_after_cuts_E[4]->fill();
      if (Q2 >40 && Q2 <  60) _Nevt_after_cuts_E[5]->fill();
      if (Q2 >60 && Q2 <  80) _Nevt_after_cuts_E[6]->fill();
      if (Q2 >80 && Q2 < 100) _Nevt_after_cuts_E[7]->fill();

      if (Q2 >  12 && Q2 <  30 && x>6e-4 && x<2e-3)  _Nevt_after_cuts_N[0]->fill();
      if (Q2 >  12 && Q2 <  30 && x>2e-3 && x<1e-2)  _Nevt_after_cuts_N[1]->fill();
      if (Q2 >  30 && Q2 <  80 && x>6e-4 && x<2e-3)  _Nevt_after_cuts_N[2]->fill();
      if (Q2 >  30 && Q2 <  80 && x>2e-3 && x<1e-2)  _Nevt_after_cuts_N[3]->fill();
      if (Q2 > 100 && Q2 < 500 && x>2e-3 && x<1e-2)  _Nevt_after_cuts_N[4]->fill();
      if (Q2 > 100 && Q2 < 500 && x>1e-2 && x<2e-1)  _Nevt_after_cuts_N[5]->fill();

      if ( Q2 >12 && Q2 < 100 )  {
        _Nevt_after_cuts_Qlow->fill();
      }

      if ( Q2 >100 && Q2 < 8000 )  {
        _Nevt_after_cuts_QHigh->fill();
      }

      double multi=0;
      double multiperevent1 = 0.0;
      double multiperevent2 = 0.0;
      double multiperevent3 = 0.0;
      double multiperevent4 = 0.0;
      double multiperevent5 = 0.0;
      double multiperevent6 = 0.0;


      const LorentzTransform Breitboost = dk.boostBreit();


      for (size_t ip1 = 0; ip1 < particles.size(); ++ip1) {
        const Particle& p = particles[ip1];

        const FourMomentum BreMom = Breitboost.transform(p.momentum());

        if ( BreMom.pz() > 0. ) continue;
        double pcal= std::sqrt(BreMom.px2() + BreMom.py2()+ BreMom.pz2());
        double xp = 2*pcal/(sqrt(Q2));

        double E = std::sqrt((.27*.27)+(pcal*pcal));
        double dp = 4*M_PI*pcal*pcal;
        double dE = pcal;
        double factor = dE/dp;

        if ( Q2 >12 && Q2 < 100 )  {
          _h["xp"] -> fill(xp);
          _h["xi"] -> fill(log(1/(xp)));
        }

        if ( Q2 >100 && Q2 < 8000 )  {
          _h["xpQgt"] -> fill(xp);
          _h["xiQgt"] -> fill(log(1/(xp)));
        }

        if (Q  >  5 && Q < 6.14)  _h["E0"]->fill(E, factor);
        if (Q  > 19 && Q < 8000)  _h["E1"]->fill(E, factor);
        if (Q2 > 12 && Q2 <  15)  _h["E2"]->fill(E, factor);
        if (Q2 > 15 && Q2 <  25)  _h["E3"]->fill(E, factor);
        if (Q2 > 20 && Q2 <  40)  _h["E4"]->fill(E, factor);
        if (Q2 > 40 && Q2 <  60)  _h["E5"]->fill(E, factor);
        if (Q2 > 60 && Q2 <  80)  _h["E6"]->fill(E, factor);
        if (Q2 > 80 && Q2 < 100)  _h["E7"]->fill(E, factor);

        _h_Q2_xp.fill(sqrt(Q2),xp);

        multi = multi+1;

        if (Q2 >  12 && Q2 <  30 && x>6e-4 && x<2e-3)   multiperevent1 = multiperevent1+1;
        if (Q2 >  12 && Q2 <  30 && x>2e-3 && x<1e-2)   multiperevent2 = multiperevent2+1;
        if (Q2 >  30 && Q2 <  80 && x>6e-4 && x<2e-3)   multiperevent3 = multiperevent3+1;
        if (Q2 >  30 && Q2 <  80 && x>2e-3 && x<1e-2)   multiperevent4 = multiperevent4+1;
        if (Q2 > 100 && Q2 < 500 && x>2e-3 && x<1e-2)   multiperevent5 = multiperevent5+1;
        if (Q2 > 100 && Q2 < 500 && x>1e-2 && x<2e-1)   multiperevent6 = multiperevent6+1;

      }

      if (Q2 >  12 && Q2 <  30 && x>6e-4 && x<2e-3)   _h["N0"]->fill(multiperevent1);
      if (Q2 >  12 && Q2 <  30 && x>2e-3 && x<1e-2)   _h["N1"]->fill(multiperevent2);
      if (Q2 >  30 && Q2 <  80 && x>6e-4 && x<2e-3)   _h["N2"]->fill(multiperevent3);
      if (Q2 >  30 && Q2 <  80 && x>2e-3 && x<1e-2)   _h["N3"]->fill(multiperevent4);
      if (Q2 > 100 && Q2 < 500 && x>2e-3 && x<1e-2)   _h["N4"]->fill(multiperevent5);
      if (Q2 > 100 && Q2 < 500 && x>1e-2 && x<2e-1)   _h["N5"]->fill(multiperevent6);

      if (Q2 > 12) {
       _h["AvgT2"] -> fill(sqrt(Q2), multi);
       _h["QT2"] -> fill(sqrt(Q2));
      }

      _h["AvgT1"] -> fill(sqrt(Q2), multi);
      _h["QT1"] -> fill(sqrt(Q2));

      _h["MeanTest1"] -> fill(sqrt(Q2));
      _h["MeanTest2"] -> fill(sqrt(Q2));

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      MSG_DEBUG("Nevt Qlow " << dbl(*_Nevt_after_cuts_Qlow));
      scale(_h["xp"], 1.0/ *_Nevt_after_cuts_Qlow);
      scale(_h["xi"], 1.0/ *_Nevt_after_cuts_Qlow);
      scale(_h["xpQgt"], 1.0/ *_Nevt_after_cuts_QHigh);
      scale(_h["xiQgt"], 1.0/ *_Nevt_after_cuts_QHigh);

      for(int iE=0 ; iE< 8 ; ++iE){
        scale(_h["E"+to_string(iE)], 1.0/ *_Nevt_after_cuts_E[iE]);
      }

      for(int iN=0 ; iN< 8 ; ++iN){
        scale(_h["N"+to_string(iN)], 1.0/ *_Nevt_after_cuts_N[iN]);
      }

      divide(_h["AvgT1"], _h["QT1"],_s["Avg1"] );
      divide(_h["AvgT2"], _h["QT2"],_s["Avg2"] );



      int iQ = 0;
      double mean;
      for (Histo1DPtr histo : _h_Q2_xp.histos()) {
        double Qerr = (QEdges[iQ+1] - QEdges[iQ])/2. ;
        double Qmid = (QEdges[iQ+1] + QEdges[iQ])/2. ;
        double Nev = dbl(*_Nevt_after_cuts_Q[iQ]) ;
        if (Nev != 0) scale(histo, 1./Nev);

        for (size_t iP = 0; iP < iPmax; ++iP) {
          mean = histo->bin(iP).height() ;
          double mean_err = mean/100;
          _s["Qxp"+to_string(iP)]->addPoint(Qmid, mean, Qerr, mean_err);
        }
        ++iQ;
      }

      const double x1 = _h["MeanTest1"]->xMean(false);
      const double x2 = _h["MeanTest2"]->xMean(false);
      MSG_DEBUG("Mean of low Q = " << x1);
      MSG_DEBUG("Mean of High Q = " << x2);
    }

    /// @}


  private:

    /// @name Histograms
    /// @{

    CounterPtr _Nevt_after_cuts_Qlow,_Nevt_after_cuts_QHigh,_Nevt_after_cuts_E[8],_Nevt_after_cuts_N[6];
    CounterPtr _Nevt_after_cuts_Q[12];
    map<string, CounterPtr> _Nevt;
    map<string, Histo1DPtr> _h;
    map<string, Scatter2DPtr> _s;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    Scatter2DPtr _h_pt_06_ratio;

    BinnedHistogram _h_Q2_xp,_h_Avg,_h_xipeak,_h_xiwidth;

    /// @}

  };


  RIVET_DECLARE_PLUGIN(H1_1997_I445116);

}
