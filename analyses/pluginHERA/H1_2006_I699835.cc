// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/DISFinalState.hh"
#include "Rivet/Projections/Thrust.hh"

namespace Rivet {



  /// @brief Measurement of Event Shape Variables in Deep-Inelastic Scattering at HERA (H1)
  class H1_2006_I699835 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(H1_2006_I699835);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      //cerr << endl << endl << "Initializing --------------------------------" << endl << endl;
      
      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 4.9);
      declare(fs, "FS");

      const DISFinalState DISfs(DISFinalState::BoostFrame::BREIT);

      const FinalState DISfsCut(DISfs, Cuts::eta < 0);

      declare(Thrust(DISfsCut), "ThrustCut");
      
      declare(DISKinematics(), "Kinematics");
      
      //Book histograms for different Q ranges:
      book(_Nevt_after_cuts, "TMP/Nevt_after_cuts");
      
      Histo1DPtr dummy;
      for(int iQ = 0; iQ < 7; ++iQ) {
         // cout << " iQ " << iQ << " " << QEdges[iQ] << " " << QEdges[iQ+1] << endl;
        _h_tauc.add(QEdges[iQ], QEdges[iQ+1], book(dummy,1+iQ,1,1));
        _h_tau.add(QEdges[iQ], QEdges[iQ+1], book(dummy,8+iQ,1,1));
        _h_B.add(QEdges[iQ], QEdges[iQ+1], book(dummy,15+iQ,1,1));
        _h_rho.add(QEdges[iQ], QEdges[iQ+1], book(dummy,22+iQ,1,1));
        book(_Nevt_after_cuts_Q[iQ], "TMP/Nevt_after_cuts_Q"+ to_string(iQ));     
      }
      
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      
      
      const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");
      
      // The kinematic region covered by the analysis is defined by ranges of Q² and y
      if (dk.Q2() < 196 or dk.Q2() > 40000 or dk.y() < 0.1 or dk.y() > 0.7) {
      	// cerr << "Event out of the kinematic region covered by the analysis" << endl;
      	vetoEvent;
      }

      double Q=sqrt(dk.Q2());
      _Nevt_after_cuts -> fill();
       for(int iQ = 0; iQ < 7; ++iQ) {
          if (inRange(Q, QEdges[iQ],QEdges[iQ+1])) {
         //  cout << " Q " << Q << " " << iQ << endl;
          _Nevt_after_cuts_Q[iQ] -> fill(); }
       }

      
      // Boost to hadronic Breit frame
      const LorentzTransform breitboost = dk.boostBreit();
      
      const FinalState& fs = apply<FinalState>(event, "FS");
      
      /*Calculate event shape variables:
      thrust_num is \sum |pz_h|
      thrust_den is \sum |p_h|
      b_num is \sum |pt_h|
      sumE is the sum of energies
      (All sums run through the particles in the current hemisphere)
      */
      double thrust_num, thrust_den, b_num, sumE;
      thrust_num = thrust_den = b_num = sumE = 0;
      Vector3 sumMom;
      
      for (const Particle& p : fs.particles()) {
      	// Boost to Breit frame
        const FourMomentum BreitMom = breitboost.transform(p.momentum());
      	if (BreitMom.eta() < 0) {
      		thrust_num += abs(BreitMom.pz());
      		thrust_den += BreitMom.p();
      		b_num += abs(BreitMom.pt());
      		sumMom.operator+=(BreitMom.p3());
      		sumE += BreitMom.E();
      	}
      }
      
      //The energy in the current hemisphere must exceed a certain value.
      if (sumE <= Q/10.0) {
      	// cerr << "Energy in the current hemisphere too low" << endl;
      	vetoEvent;
      }
      
      /* Comment from A. Galvan:
      Thrust here is with respect to the z axis (tau in the paper), while the
      one from Rivet projection is with respect to the maximum thrust 
      axis (1 - tau_c in the paper)
      */
      
      double thrust_mine = 1.0 - ((double)thrust_num)/((double)thrust_den);
      double b_mine = ((double)b_num)/(2.0*(double)thrust_den);
      double rho_num = thrust_den*thrust_den - sumMom.dot(sumMom);
      double rho_mine = ((double)rho_num)/(4.0*(double)thrust_den*(double)thrust_den);
      
      const Thrust& thrCut = applyProjection<Thrust>(event, "ThrustCut");

      //Fill histograms:
      _h_tauc.fill(Q, 1.0 - thrCut.thrust());
      _h_tau.fill(Q, thrust_mine);
      _h_B.fill(Q, b_mine);
      _h_rho.fill(Q, rho_mine);
      
      /* Comment from A. Galvan:
       As for the C-parameter, my results did not fit the reference data at
       all. The formula given in the paper is a bit ambiguous, because there is
       a sum that runs through all pairs of particles h,h' and I was not sure
       whether each pair should be counted twice (h,h' and h',h) or not, or if
       the pair of a particle with itself (hh) should be considered. 
       That is why I tried all of the possibilities, but none of them worked.
      */
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      //cerr << endl << endl << "Finalizing --------------------------------" << endl << endl;
      
      // cerr << "Cross section: " << crossSection() << endl;
      
      // double Nev = dbl(*_Nevt_after_cuts) ;
      // cout << " N total " << Nev << endl; 


      int iQ=0;
      for (Histo1DPtr histo : _h_tauc.histos()) { 
          double Nev = dbl(*_Nevt_after_cuts_Q[iQ]) ;
          // cout << " Nev " << Nev << " " << iQ << endl;
          if (Nev != 0) scale(histo, 1./Nev);
          vector<YODA::HistoBin1D>& bins = histo -> bins();
          for (auto & b : bins) b.scaleW(b.xWidth());
          ++iQ;
      }
      
      iQ=0;
      for (Histo1DPtr histo : _h_tau.histos()) { 
          double Nev = dbl(*_Nevt_after_cuts_Q[iQ]) ;
          //cout << " Nev " << Nev << " " << iQ << endl;
          if (Nev != 0) scale(histo, 1./Nev);
          vector<YODA::HistoBin1D>& bins = histo -> bins();
          for (auto & b : bins) b.scaleW(b.xWidth());
          ++iQ;
      }
      
      iQ = 0;
      for (Histo1DPtr histo : _h_B.histos()) { 
          double Nev = dbl(*_Nevt_after_cuts_Q[iQ]) ;
          // cout << " Nev " << Nev << " " << iQ << endl;
          if (Nev != 0) scale(histo, 1./Nev);
          vector<YODA::HistoBin1D>& bins = histo -> bins();
          for (auto & b : bins) b.scaleW(b.xWidth());
          ++iQ;
      }
      
      iQ = 0;
      for (Histo1DPtr histo : _h_rho.histos()) { 
          double Nev = dbl(*_Nevt_after_cuts_Q[iQ]) ;
          // cout << " Nev " << Nev << " " << iQ << endl;
          if (Nev != 0) scale(histo, 1./Nev);
          vector<YODA::HistoBin1D>& bins = histo -> bins();
          for (auto & b : bins) b.scaleW(b.xWidth());
          ++iQ;
      }

    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    BinnedHistogram _h_tauc, _h_tau, _h_B, _h_rho;
    CounterPtr _Nevt_after_cuts;
    CounterPtr _Nevt_after_cuts_Q[7];
    ///@}
    const vector<double> QEdges {14., 16., 20., 30., 50., 70., 100., 200.};


  };


  RIVET_DECLARE_PLUGIN(H1_2006_I699835);

}
