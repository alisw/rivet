// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief Charged particle multiplicities in deep inelastic scattering at HERA (H1)
  class H1_1996_I422230 : public Analysis {
  public:

    const vector<double> WEdges {80., 115, 150., 185., 220.};


    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(H1_1996_I422230);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      const DISLepton disl;
      declare(disl, "Lepton");
      declare(DISKinematics(), "Kinematics");

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs;
      declare(fs,"FS");
      const ChargedFinalState cfs(disl.remainingFinalState());
      declare(cfs,"CFS");



      //binned histograms to count for the multiplicity
      for (size_t ix = 0; ix < 4; ++ix) {
        book(_Nevt_after_cuts[ix], "TMP/Nevt_after_cuts"+ to_string(ix));
      }
      Histo1DPtr dummy;
      for (int iW = 0; iW < 4; ++iW) {
        // cout << " iW " << iW << " " << WEdges[iW+1] << " " << WEdges[iW] << endl;
        _h_mult1.add(WEdges[iW], WEdges[iW+1],book(dummy,iW+1,1,1));
        _h_mult2.add(WEdges[iW], WEdges[iW+1],book(dummy,iW+1,1,2));
        _h_mult3.add(WEdges[iW], WEdges[iW+1],book(dummy,iW+1,1,3));
        _h_mult4.add(WEdges[iW], WEdges[iW+1],book(dummy,iW+1,1,4));
        _h_mult_all.add(WEdges[iW], WEdges[iW+1],book(dummy,"TMP/dummy"+ to_string(iW),refData(4,1,4)));
        _h_mult10_all.add(WEdges[iW], WEdges[iW+1],book(dummy,"TMP/dummy1"+ to_string(iW),refData(4,1,4)));
        _h_mult11_all.add(WEdges[iW], WEdges[iW+1],book(dummy,"TMP/dummy2"+ to_string(iW),refData(4,1,4)));
        _h_mult12_all.add(WEdges[iW], WEdges[iW+1],book(dummy,"TMP/dummy3"+ to_string(iW),refData(4,1,4)));
      }
      //histograms for the statistical moments and the mean

      book(_h_mean0,5,1,1);
      book(_h_D2_0,5,1,2);
      book(_h_D3_0,5,1,3);
      book(_h_D4_0,5,1,4);
      book(_h_C2_0,5,1,5);
      book(_h_C3_0,5,1,6);
      book(_h_C4_0,5,1,7);
      book(_h_R2_0,5,1,8);
      book(_h_R3_0,5,1,9);

      book(_h_mean12,6,1,1);
      book(_h_D2_12,6,1,2);
      book(_h_D3_12,6,1,3);
      book(_h_D4_12,6,1,4);
      book(_h_C2_12,6,1,5);
      book(_h_C3_12,6,1,6);
      book(_h_C4_12,6,1,7);
      book(_h_R2_12,6,1,8);
      book(_h_R3_12,6,1,9);
      book(_h_K3_12,6,1,10);

      book(_h_mean13,7,1,1);
      book(_h_D2_13,7,1,2);
      book(_h_D3_13,7,1,3);
      book(_h_D4_13,7,1,4);
      book(_h_C2_13,7,1,5);
      book(_h_C3_13,7,1,6);
      book(_h_C4_13,7,1,7);
      book(_h_R2_13,7,1,8);
      book(_h_R3_13,7,1,9);
      book(_h_K3_13,7,1,10);

      book(_h_mean14,8,1,1);
      book(_h_D2_14,8,1,2);
      book(_h_D3_14,8,1,3);
      book(_h_D4_14,8,1,4);
      book(_h_C2_14,8,1,5);
      book(_h_C3_14,8,1,6);
      book(_h_C4_14,8,1,7);
      book(_h_R2_14,8,1,8);
      book(_h_R3_14,8,1,9);

      book(_h_mean15,9,1,1);
      book(_h_D2_15,9,1,2);
      book(_h_D3_15,9,1,3);
      book(_h_D4_15,9,1,4);
      book(_h_C2_15,9,1,5);
      book(_h_C3_15,9,1,6);
      book(_h_C4_15,9,1,7);
      book(_h_R2_15,9,1,8);
      book(_h_R3_15,9,1,9);

      book(_h_mean23,10,1,1);
      book(_h_D2_23,10,1,2);
      book(_h_D3_23,10,1,3);
      book(_h_D4_23,10,1,4);
      book(_h_C2_23,10,1,5);
      book(_h_C3_23,10,1,6);
      book(_h_C4_23,10,1,7);
      book(_h_R2_23,10,1,8);
      book(_h_R3_23,10,1,9);
      book(_h_K3_23,10,1,10);

      book(_h_mean34,11,1,1);
      book(_h_D2_34,11,1,2);
      book(_h_D3_34,11,1,3);
      book(_h_D4_34,11,1,4);
      book(_h_C2_34,11,1,5);
      book(_h_C3_34,11,1,6);
      book(_h_C4_34,11,1,7);
      book(_h_R2_34,11,1,8);
      book(_h_R3_34,11,1,9);
      book(_h_K3_34,11,1,10);

      book(_h_mean45,12,1,1);
      book(_h_D2_45,12,1,2);
      book(_h_D3_45,12,1,3);
      book(_h_D4_45,12,1,4);
      book(_h_C2_45,12,1,5);
      book(_h_C3_45,12,1,6);
      book(_h_C4_45,12,1,7);
      book(_h_R2_45,12,1,8);
      book(_h_R3_45,12,1,9);
      book(_h_K3_45,12,1,10);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      //apply the final state of all particles and the one for only charged particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event,"CFS");
      const Particles& particles = cfs.particles();

      //because it does not make sense to have a collision with the numparticles is less than two,we use the vetoEvent so that if there is an event like this it does not run the analysis and goes to the next one

      if (numParticles<2) {
        MSG_DEBUG("Failed leptonic cut");
        vetoEvent;
      }

      //apply DIS kinematics in the event
      const DISKinematics& dk = apply<DISKinematics>(event,"Kinematics");
      //const DISLepton& dl = apply<DISLepton>(event,"Lepton");

      //get the DIS Kinematics but not in a loop because they are variables that descrube the type of event not the particles
      double Q2 = dk.Q2();
      double W2 = dk.W2();
      double W = std::sqrt(W2);
      bool cut1 = W<220. && W>80.;
      if (!cut1) vetoEvent;
      bool cut = Q2<1000. && Q2>10. && W>80. && W<220.;
      if (!cut) vetoEvent;

      double Efwd = 0. ;
      for (size_t ip1 = 0; ip1 < particles.size(); ++ip1) {
        const Particle& p = particles[ip1];
        double theta = p.theta()/degree ;
        if ( inRange(theta,4.4,15.) ) {
          Efwd = Efwd + p.E() ;
        }
      }

      bool cut_fwd = Efwd > 0.5 && dk.beamLepton().E() > 12. ;
      if (!cut_fwd) vetoEvent ;

      for (int iW = 0; iW < 4; ++iW) {
        if (inRange(W, WEdges[iW],WEdges[iW+1])) {
          _Nevt_after_cuts[iW] -> fill();
        }
      }


      //boost to the hadronic centre of mass frame
      const LorentzTransform hcmboost = dk.boostHCM();
      double kall = 0.;
      double k1 = 0.;
      double k2 = 0.;
      double k3 = 0.;
      double k4 = 0.;
      double k10 = 0.;
      double k11 = 0.;
      double k12 = 0.;
      for (size_t ip1 = 0; ip1 < particles.size(); ++ip1) {
        const Particle& p = particles[ip1];
        const FourMomentum hcmMom = hcmboost.transform(p.momentum());
        const double etahcm_char = hcmMom.eta();
        if (etahcm_char>0.) {
          kall = kall+1;
        }
        if (etahcm_char>1. && etahcm_char<2.) {
          k1 = k1+1;
        }
        if (etahcm_char>1. && etahcm_char<3.) {
          k2 = k2+1;
        }
        if (etahcm_char>1. && etahcm_char<4.) {
          k3 = k3+1;
        }
        if (etahcm_char>1. && etahcm_char<5.) {
          k4 = k4+1;
        }
        if (etahcm_char>2. && etahcm_char<3.) {
          k10 = k10+1;
        }
        if (etahcm_char>3. && etahcm_char<4.) {
          k11 = k11+1;
        }
        if (etahcm_char>4. && etahcm_char<5.) {
          k12 = k12+1;
        }

      }
      // cout<<k1<<endl;
      _h_mult_all.fill(W,kall);
      _h_mult10_all.fill(W,k10);
      _h_mult11_all.fill(W,k11);
      _h_mult12_all.fill(W,k12);
      _h_mult1.fill(W,k1);
      _h_mult2.fill(W,k2);
      _h_mult3.fill(W,k3);
      _h_mult4.fill(W,k4);

    }

    /// Normalise histograms etc., after the run
    void finalize() {
      // cout << _h_mult1.histo() << endl;
      int iW = 0 ;
      double iq  ;
      double mean, dispersion, cq ,R2,R3,K3 ;
      for (Histo1DPtr histo : _h_mult_all.histos()) {
        double Werr = (WEdges[iW+1] - WEdges[iW])/2. ;
        //cout << " iW " << iW << " Edges " << WEdges[iW+1] << " " << WEdges[iW] << " Werr = " << Werr << endl;
        double Wmid = (WEdges[iW+1] + WEdges[iW])/2. ;
        //cout << " new histo: W = "  << Wmid << " " << histo->name() << endl;
        iq = 2 ;
        _histo_to_moments(histo, iq , mean, dispersion,cq,R2,R3,K3) ;
        //cout << " mean " << mean << " dispersion " << dispersion << " for iq = " << iq << endl;

        // just to have some values, needs to be corrected

        _h_mean0->addPoint(Wmid, mean, Werr, dispersion/2);
        _h_D2_0->addPoint(Wmid, dispersion, Werr, dispersion/mean);
        _h_C2_0->addPoint(Wmid, cq, Werr, cq/mean);
        _h_R2_0->addPoint(Wmid, R2, Werr, R2/mean);
        _h_R3_0->addPoint(Wmid,R3,Werr,R3/mean);
        iq = 3 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq , mean, dispersion,cq,R2,R3,K3) ;
        _h_D3_0->addPoint(Wmid, dispersion, Werr, dispersion/mean);
        _h_C3_0->addPoint(Wmid, cq, Werr, cq/mean);
        iq = 4 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq , mean, dispersion,cq,R2,R3,K3) ;
        _h_D4_0->addPoint(Wmid, dispersion, Werr, dispersion/mean);
        _h_C4_0->addPoint(Wmid, cq, Werr, cq/mean);

        ++iW;
      }
      iW = 0 ;
      mean = 0.;
      dispersion = 0.;
      cq = 0.;
      R2 = 0;
      R3 = 0.;
      K3 = 0.;
      for (Histo1DPtr histo : _h_mult1.histos()) {
        // use scale to get Prob in %, and use counter to count events after cuts
        //cout << " Nevt " << dbl(*_Nevt_after_cuts[iW]) << endl;
        scale(histo, 100.0/ *_Nevt_after_cuts[iW]);
        //
        double Werr = (WEdges[iW+1] - WEdges[iW])/2. ;
        double Wmid = (WEdges[iW+1] + WEdges[iW])/2. ;
        iq = 2 ;
        _histo_to_moments(histo, iq , mean, dispersion,cq,R2,R3,K3) ;

        _h_mean12->addPoint(Wmid, mean, Werr, dispersion/2);
        _h_D2_12->addPoint(Wmid, dispersion, Werr, dispersion/mean);
        _h_C2_12->addPoint(Wmid, cq, Werr, cq/mean);
        _h_R2_12->addPoint(Wmid, R2, Werr, R2/mean);
        _h_R3_12->addPoint(Wmid,R3,Werr,R3/mean);
        _h_K3_12->addPoint(Wmid,K3,Werr,K3/mean);
        iq = 3 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq , mean, dispersion,cq,R2,R3,K3) ;
        _h_D3_12->addPoint(Wmid, dispersion, Werr, dispersion/mean);
        _h_C3_12->addPoint(Wmid, cq, Werr, cq/mean);
        iq = 4 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq , mean, dispersion,cq,R2,R3,K3) ;
        _h_D4_12->addPoint(Wmid, dispersion, Werr, dispersion/mean);
        _h_C4_12->addPoint(Wmid, cq, Werr, cq/mean);
        ++iW;
      }
      iW = 0 ;
      mean = 0.;
      dispersion = 0.;
      cq = 0.;
      R2 = 0;
      R3 = 0.;
      K3 = 0.;

      for (Histo1DPtr histo : _h_mult2.histos()) {
        // use scale to get Prob in %, and use counter to count events after cuts
        //cout << " Nevt " << dbl(*_Nevt_after_cuts[0]) << endl;
        scale(histo, 100.0/ *_Nevt_after_cuts[iW]);
        //
        double Werr = (WEdges[iW+1] - WEdges[iW])/2. ;
        double Wmid = (WEdges[iW+1] + WEdges[iW])/2. ;
        iq = 2 ;
        _histo_to_moments(histo, iq , mean, dispersion,cq,R2,R3,K3) ;

        _h_mean13->addPoint(Wmid, mean, Werr, dispersion/2);
        _h_D2_13->addPoint(Wmid, dispersion, Werr, dispersion/mean);
        _h_C2_13->addPoint(Wmid, cq, Werr, cq/mean);
        _h_R2_13->addPoint(Wmid, R2, Werr, R2/mean);
        _h_R3_13->addPoint(Wmid,R3,Werr,R3/mean);
        _h_K3_13->addPoint(Wmid,K3,Werr,K3/mean);
        iq = 3 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq , mean, dispersion,cq,R2,R3,K3) ;
        _h_D3_13->addPoint(Wmid, dispersion, Werr, dispersion/mean);
        _h_C3_13->addPoint(Wmid, cq, Werr, cq/mean);
        iq = 4 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq , mean, dispersion,cq,R2,R3,K3) ;
        _h_D4_13->addPoint(Wmid, dispersion, Werr, dispersion/mean);
        _h_C4_13->addPoint(Wmid, cq, Werr, cq/mean);
        ++iW;
      }
      iW = 0 ;
      mean = 0.;
      dispersion = 0.;
      cq = 0.;
      R2 = 0;
      R3 = 0.;
      K3 = 0.;

      for (Histo1DPtr histo : _h_mult3.histos()) {
        // use scale to get Prob in %, and use counter to count events after cuts
        //cout << " Nevt " << dbl(*_Nevt_after_cuts[0]) << endl;
        scale(histo, 100.0/ *_Nevt_after_cuts[iW]);
        //
        double Werr = (WEdges[iW+1] - WEdges[iW])/2. ;
        double Wmid = (WEdges[iW+1] + WEdges[iW])/2. ;
        iq = 2 ;
        _histo_to_moments(histo, iq , mean, dispersion,cq,R2,R3,K3) ;
        _h_mean14->addPoint(Wmid, mean, Werr, dispersion/2);
        _h_D2_14->addPoint(Wmid, dispersion, Werr, dispersion/mean);
        _h_C2_14->addPoint(Wmid, cq, Werr, cq/mean);
        _h_R2_14->addPoint(Wmid, R2, Werr, R2/mean);
        _h_R3_14->addPoint(Wmid,R3,Werr,R3/mean);
        iq = 3 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq , mean, dispersion,cq,R2,R3,K3) ;
        _h_D3_14->addPoint(Wmid, dispersion, Werr, dispersion/mean);
        _h_C3_14->addPoint(Wmid, cq, Werr, cq/mean);
        iq = 4 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq , mean, dispersion,cq,R2,R3,K3) ;
        _h_D4_14->addPoint(Wmid, dispersion, Werr, dispersion/mean);
        _h_C4_14->addPoint(Wmid, cq, Werr, cq/mean);
        ++iW;
      }
      iW = 0 ;
      mean = 0.;
      dispersion = 0.;
      cq = 0.;
      R2 = 0;
      R3 = 0.;
      K3 = 0.;

      for (Histo1DPtr histo : _h_mult4.histos()) {
        // use scale to get Prob in %, and use counter to count events after cuts
        //cout << " Nevt " << dbl(*_Nevt_after_cuts[0]) << endl;
        scale(histo, 100.0/ *_Nevt_after_cuts[iW]);
        //
        double Werr = (WEdges[iW+1] - WEdges[iW])/2. ;
        double Wmid = (WEdges[iW+1] + WEdges[iW])/2. ;
        iq = 2 ;
        _histo_to_moments(histo, iq , mean, dispersion,cq,R2,R3,K3) ;
        _h_mean15->addPoint(Wmid, mean, Werr, dispersion/2);
        _h_D2_15->addPoint(Wmid, dispersion, Werr, dispersion/mean);
        _h_C2_15->addPoint(Wmid, cq, Werr, cq/mean);
        _h_R2_15->addPoint(Wmid, R2, Werr, R2/mean);
        _h_R3_15->addPoint(Wmid,R3,Werr,R3/mean);
        iq = 3 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq , mean, dispersion,cq,R2,R3,K3) ;
        _h_D3_15->addPoint(Wmid, dispersion, Werr, dispersion/mean);
        _h_C3_15->addPoint(Wmid, cq, Werr, cq/mean);
        iq = 4 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq , mean, dispersion,cq,R2,R3,K3) ;
        _h_D4_15->addPoint(Wmid, dispersion, Werr, dispersion/mean);
        _h_C4_15->addPoint(Wmid, cq, Werr, cq/mean);
        ++iW;
      }
      iW = 0 ;
      mean = 0.;
      dispersion = 0.;
      cq = 0.;
      R2 = 0;
      R3 = 0.;
      K3 = 0.;

      for (Histo1DPtr histo : _h_mult10_all.histos()) {
        // use scale to get Prob in %, and use counter to count events after cuts
        //cout << " Nevt " << dbl(*_Nevt_after_cuts[0]) << endl;
        scale(histo, 100.0/ *_Nevt_after_cuts[iW]);
        //
        double Werr = (WEdges[iW+1] - WEdges[iW])/2. ;
        double Wmid = (WEdges[iW+1] + WEdges[iW])/2. ;
        iq = 2 ;
        _histo_to_moments(histo, iq , mean, dispersion,cq,R2,R3,K3) ;
        _h_mean23->addPoint(Wmid, mean, Werr, dispersion/2);
        _h_D2_23->addPoint(Wmid, dispersion, Werr, dispersion/mean);
        _h_C2_23->addPoint(Wmid, cq, Werr, cq/mean);
        _h_R2_23->addPoint(Wmid, R2, Werr,R2/mean);
        _h_R3_23->addPoint(Wmid,R3,Werr,R3/mean);
        _h_K3_23->addPoint(Wmid,K3,Werr,K3/mean);
        iq = 3 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq , mean, dispersion,cq,R2,R3,K3) ;
        _h_D3_23->addPoint(Wmid, dispersion, Werr, dispersion/mean);
        _h_C3_23->addPoint(Wmid, cq, Werr, cq/mean);
        iq = 4 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq , mean, dispersion,cq,R2,R3,K3) ;
        _h_D4_23->addPoint(Wmid, dispersion, Werr, dispersion/mean);
        _h_C4_23->addPoint(Wmid, cq, Werr, cq/mean);
        ++iW;
      }

      iW = 0 ;
      mean = 0.;
      dispersion = 0.;
      cq = 0.;
      R2 = 0;
      R3 = 0.;
      K3 = 0.;

      for (Histo1DPtr histo : _h_mult11_all.histos()) {
        // use scale to get Prob in %, and use counter to count events after cuts
        //cout << " Nevt " << dbl(*_Nevt_after_cuts[0]) << endl;
        scale(histo, 100.0/ *_Nevt_after_cuts[iW]);
        //
        double Werr = (WEdges[iW+1] - WEdges[iW])/2. ;
        double Wmid = (WEdges[iW+1] + WEdges[iW])/2. ;
        iq = 2 ;
        _histo_to_moments(histo, iq , mean, dispersion,cq,R2,R3,K3) ;
        _h_mean34->addPoint(Wmid, mean, Werr, dispersion/2);
        _h_D2_34->addPoint(Wmid, dispersion, Werr, dispersion/mean);
        _h_C2_34->addPoint(Wmid, cq, Werr, cq/mean);
        _h_R2_34->addPoint(Wmid, R2, Werr, R2/mean);
        _h_R3_34->addPoint(Wmid,R3,Werr,R3/mean);
        _h_K3_34->addPoint(Wmid,K3,Werr,K3/mean);
        iq = 3 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq , mean, dispersion,cq,R2,R3,K3) ;
        _h_D3_34->addPoint(Wmid, dispersion, Werr, dispersion/mean);
        _h_C3_34->addPoint(Wmid, cq, Werr, cq/mean);
        iq = 4 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq , mean, dispersion,cq,R2,R3,K3) ;
        _h_D4_34->addPoint(Wmid, dispersion, Werr, dispersion/mean);
        _h_C4_34->addPoint(Wmid, cq, Werr, cq/mean);
        ++iW;
      }
      iW = 0 ;
      mean = 0.;
      dispersion = 0.;
      cq = 0.;
      R2 = 0;
      R3 = 0.;
      K3 = 0.;

      for (Histo1DPtr histo : _h_mult12_all.histos()) {
        // use scale to get Prob in %, and use counter to count events after cuts
        //cout << " Nevt " << dbl(*_Nevt_after_cuts[0]) << endl;
        scale(histo, 100.0/ *_Nevt_after_cuts[iW]);
        //
        double Werr = (WEdges[iW+1] - WEdges[iW])/2. ;
        double Wmid = (WEdges[iW+1] + WEdges[iW])/2. ;
        iq = 2 ;
        _histo_to_moments(histo, iq , mean, dispersion,cq,R2,R3,K3) ;
        _h_mean45->addPoint(Wmid, mean, Werr, dispersion/2);
        _h_D2_45->addPoint(Wmid, dispersion, Werr, dispersion/mean);
        _h_C2_45->addPoint(Wmid, cq, Werr, cq/mean);
        _h_R2_45->addPoint(Wmid, R2, Werr, R2/mean);
        _h_R3_45->addPoint(Wmid,R3,Werr,R3/mean);
        // cout<<"R2= "<<R2<<"R3= "<<R3<<"K3= "<<K3<<endl;
        _h_K3_45->addPoint(Wmid,K3,Werr,K3/mean);
        iq = 3 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq , mean, dispersion,cq,R2,R3,K3) ;
        _h_D3_45->addPoint(Wmid, dispersion, Werr, dispersion/mean);
        _h_C3_45->addPoint(Wmid, cq, Werr, cq/mean);
        iq = 4 ;
        dispersion = 0 ;
        cq = 0;
        _histo_to_moments(histo, iq , mean, dispersion,cq,R2,R3,K3) ;
        _h_D4_45->addPoint(Wmid, dispersion, Werr, dispersion/mean);
        _h_C4_45->addPoint(Wmid, cq, Werr, cq/mean);
        ++iW;
      }

    }


    inline void _histo_to_moments(Histo1DPtr histo_input, double iq, double & mean, double & dispersion, double & cq,double &R2,double &R3, double &K3) {

      //cout << " histo mean = " << histo_input->xMean() << " variance " << histo_input->xVariance() << endl;
      double mysumWX = 0. ;
      double mysumW2X = 0. ;
      double mysumWX2 = 0. ;
      double mysumW2 = 0. ;
      double mysumW = 0. ;

      // cout << histo_input->name() << endl;
      if (histo_input->effNumEntries() == 0 || histo_input->sumW() == 0) {
        MSG_WARNING("Requested mean of a distribution with no net fill weights");
      } else {
        // loop to calcualte mean
        for (size_t b = 0; b < histo_input->numBins(); ++b) { // loop over points
          mysumWX  += histo_input->bin(b).height()      *  histo_input->bin(b).xMid() ;
          mysumW2X += sqr(histo_input->bin(b).height()) *  histo_input->bin(b).xMid() ;
          mysumWX2 += histo_input->bin(b).height()      *  sqr(histo_input->bin(b).xMid()) ;
          mysumW2  += sqr(histo_input->bin(b).height()) ;
          mysumW   += histo_input->bin(b).height() ;
        }
        mean = mysumWX/mysumW ;


        // loop to calculate dispersion (variance)
        double var = 0.;
        for (size_t b = 0; b < histo_input->numBins(); ++b) { // loop over points
          double xval = histo_input->bin(b).xMid() ;
          double weight = histo_input->bin(b).height() ;
          var = var + weight * pow((xval - mean),iq) ;
        }

        var = var/mysumW ;

        dispersion = pow(var,1./iq) ;
        double c = 0.;
        for (size_t b = 0; b < histo_input->numBins(); ++b) { // loop over points
          double xval = histo_input->bin(b).xMid() ;
          double weight = histo_input->bin(b).height() ;
          c = c+pow(xval,iq)*weight;
        }
        cq = c/(mysumW*pow(mean,iq));
        double r2 = 0.;
        double r3 = 0.;
        for (size_t b = 0; b < histo_input->numBins(); ++b) { // loop over points
          double xval = histo_input->bin(b).xMid() ;
          double weight = histo_input->bin(b).height() ;
          r2 = r2+xval*(xval-1)*weight;
          r3 = r3+xval*(xval-1)*(xval-2)*weight;
        }
        R2 = r2/(mysumW*pow(mean,2));
        R3 = r3/(mysumW*pow(mean,3));
        K3 = R3 - 3*R2 + 2;
      }
    }

    ///@}


  private:

    /// @name Histograms
    ///@{

    BinnedHistogram _h_mult1;
    BinnedHistogram _h_mult2;
    BinnedHistogram _h_mult3;
    BinnedHistogram _h_mult4;
    BinnedHistogram _h_mult_all;
    BinnedHistogram _h_mult10_all;
    BinnedHistogram _h_mult11_all;
    BinnedHistogram _h_mult12_all;

    Histo1DPtr _h2_W5;
    Histo1DPtr _h2_W6;
    Histo1DPtr _h2_W7;
    Histo1DPtr _h2_W8;
    Histo1DPtr _h2_W9;
    Histo1DPtr _h2_W10;
    Histo1DPtr _h2_W11;
    Histo1DPtr _h2_W12;
    //Histo1DPtr _h2_W5_R2;
    Scatter2DPtr _h_mean0;
    Scatter2DPtr _h_D2_0;
    Scatter2DPtr _h_D3_0;
    Scatter2DPtr _h_D4_0;
    Scatter2DPtr _h_C2_0;
    Scatter2DPtr _h_C3_0;
    Scatter2DPtr _h_C4_0;
    Scatter2DPtr _h_R2_0;
    Scatter2DPtr _h_R3_0;
    Scatter2DPtr _h_mean12;
    Scatter2DPtr _h_D2_12;
    Scatter2DPtr _h_D3_12;
    Scatter2DPtr _h_D4_12;
    Scatter2DPtr _h_C2_12;
    Scatter2DPtr _h_C3_12;
    Scatter2DPtr _h_C4_12;
    Scatter2DPtr _h_R2_12;
    Scatter2DPtr _h_R3_12;
    Scatter2DPtr _h_K3_12;
    Scatter2DPtr _h_mean13;
    Scatter2DPtr _h_D2_13;
    Scatter2DPtr _h_D3_13;
    Scatter2DPtr _h_D4_13;
    Scatter2DPtr _h_C2_13;
    Scatter2DPtr _h_C3_13;
    Scatter2DPtr _h_C4_13;
    Scatter2DPtr _h_R2_13;
    Scatter2DPtr _h_R3_13;
    Scatter2DPtr _h_K3_13;
    Scatter2DPtr _h_mean14;
    Scatter2DPtr _h_D2_14;
    Scatter2DPtr _h_D3_14;
    Scatter2DPtr _h_D4_14;
    Scatter2DPtr _h_C2_14;
    Scatter2DPtr _h_C3_14;
    Scatter2DPtr _h_C4_14;
    Scatter2DPtr _h_R2_14;
    Scatter2DPtr _h_R3_14;
    Scatter2DPtr _h_mean15;
    Scatter2DPtr _h_D2_15;
    Scatter2DPtr _h_D3_15;
    Scatter2DPtr _h_D4_15;
    Scatter2DPtr _h_C2_15;
    Scatter2DPtr _h_C3_15;
    Scatter2DPtr _h_C4_15;
    Scatter2DPtr _h_R2_15;
    Scatter2DPtr _h_R3_15;
    Scatter2DPtr _h_mean23;
    Scatter2DPtr _h_D2_23;
    Scatter2DPtr _h_D3_23;
    Scatter2DPtr _h_D4_23;
    Scatter2DPtr _h_C2_23;
    Scatter2DPtr _h_C3_23;
    Scatter2DPtr _h_C4_23;
    Scatter2DPtr _h_R2_23;
    Scatter2DPtr _h_R3_23;
    Scatter2DPtr _h_K3_23;
    Scatter2DPtr _h_mean34;
    Scatter2DPtr _h_D2_34;
    Scatter2DPtr _h_D3_34;
    Scatter2DPtr _h_D4_34;
    Scatter2DPtr _h_C2_34;
    Scatter2DPtr _h_C3_34;
    Scatter2DPtr _h_C4_34;
    Scatter2DPtr _h_R2_34;
    Scatter2DPtr _h_R3_34;
    Scatter2DPtr _h_K3_34;
    Scatter2DPtr _h_mean45;
    Scatter2DPtr _h_D2_45;
    Scatter2DPtr _h_D3_45;
    Scatter2DPtr _h_D4_45;
    Scatter2DPtr _h_C2_45;
    Scatter2DPtr _h_C3_45;
    Scatter2DPtr _h_C4_45;
    Scatter2DPtr _h_R2_45;
    Scatter2DPtr _h_R3_45;
    Scatter2DPtr _h_K3_45;

    CounterPtr _Nevt_after_cuts[4];

    ///@}

  };


  RIVET_DECLARE_PLUGIN(H1_1996_I422230);

}
