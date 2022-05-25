// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Math/Constants.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"

namespace Rivet {


  /// @brief H1 energy flow and charged particle spectra
  ///
  /// @author Peter Richardson
  ///
  /// Based on the HZTOOL analysis HZ99091
  class H1_2000_S4129130 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(H1_2000_S4129130);


    /// @name Analysis methods
    //@{

    /// Initialise projections and histograms
    void init() {
      // Projections
      declare(DISLepton(), "Lepton");
      declare(DISKinematics(), "Kinematics");
      declare(FinalState(), "FS");

      // Histograms and weight vectors for low Q^2 a
      _histETLowQa.resize(17);
      for (size_t ix = 0; ix < 17; ++ix) {
        book(_histETLowQa[ix], ix+1, 1, 1);
        book(_weightETLowQa[ix], "TMP/ETLowQa" + to_string(ix));
      }

      // Histograms and weight vectors for high Q^2 a
      _histETHighQa.resize(7);
      for (size_t ix = 0; ix < 7; ++ix) {
        book(_histETHighQa[ix], ix+18, 1, 1);
        book(_weightETHighQa[ix], "TMP/ETHighQa" + to_string(ix));
      }

      // Histograms and weight vectors for low Q^2 b
      _histETLowQb.resize(5);
      for (size_t ix = 0; ix < 5; ++ix) {
        book(_histETLowQb[ix], ix+25, 1, 1);
        book(_weightETLowQb[ix], "TMP/ETLowQb" + to_string(ix));
      }

      // Histograms and weight vectors for high Q^2 b
      _histETHighQb.resize(5);
      for (size_t ix = 0; ix < 3; ++ix) {
        book(_histETHighQb[ix], 30+ix, 1, 1);
        book(_weightETHighQb[ix], "TMP/ETHighQb" + to_string(ix));
      }

      // Histograms for the averages
      book(_histAverETCentral ,33,  1, 1);
      book(_histAverETFrag ,34,  1, 1);
    }


    /// Analyze each event
    void analyze(const Event& event) {

      // DIS kinematics
      const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");
      if ( dk.failed() ) vetoEvent;
      double q2  = dk.Q2();
      double x   = dk.x();
      double y   = dk.y();
      double w2  = dk.W2();

      // Kinematics of the scattered lepton
      const DISLepton& dl = apply<DISLepton>(event,"Lepton");
      if ( dl.failed() ) vetoEvent;
      const FourMomentum leptonMom = dl.out();
      const double enel = leptonMom.E();
      const double thel = 180 - leptonMom.angle(dl.in().mom())/degree;

      // Extract the particles other than the lepton
      const FinalState& fs = apply<FinalState>(event, "FS");
      Particles particles; particles.reserve(fs.size());
      ConstGenParticlePtr dislepGP = dl.out().genParticle(); ///< @todo Is the GenParticle stuff necessary? (Not included in Particle::==?)
      for(const Particle& p: fs.particles()) {
        ConstGenParticlePtr loopGP = p.genParticle();
        if (loopGP == dislepGP) continue;
        particles.push_back(p);
      }

      // Cut on the forward energy
      double efwd = 0.;
      for (const Particle& p : particles) {
        const double th = 180 - p.angle(dl.in())/degree;
        if (inRange(th, 4.4, 15.0)) efwd += p.E();
      }
      // There are four possible selections for events
      bool evcut[4];
      // Low  Q2 selection a
      evcut[0] = enel/GeV > 12. && w2 >= 4400.*GeV2 && efwd/GeV > 0.5 && inRange(thel,157.,176.);
      // Low  Q2 selection b
      evcut[1] = enel/GeV > 12. && inRange(y,0.3,0.5);
      // High Q2 selection a
      evcut[2] = inRange(thel,12.,150.) && inRange(y,0.05,0.6) && w2 >= 4400.*GeV2 && efwd > 0.5;
      // High Q2 selection b
      evcut[3] = inRange(thel,12.,150.) && inRange(y,0.05,0.6) && inRange(w2,27110.*GeV2,45182.*GeV2);

      // Veto if fails all cuts
      /// @todo Can we use all()?
      if (! (evcut[0] || evcut[1] || evcut[2] || evcut[3]) ) vetoEvent;

      // Find the bins
      int bin[4] = {-1,-1,-1,-1};
      // For the low Q2 selection a)
      if (q2 > 2.5*GeV2 && q2 <= 5.*GeV2) {
        if (x > 0.00005 && x <= 0.0001 ) bin[0] = 0;
        if (x > 0.0001  && x <= 0.0002 ) bin[0] = 1;
        if (x > 0.0002  && x <= 0.00035) bin[0] = 2;
        if (x > 0.00035 && x <= 0.0010 ) bin[0] = 3;
      }
      else if (q2 > 5.*GeV2 && q2 <= 10.*GeV2) {
        if (x > 0.0001  && x <= 0.0002 ) bin[0] = 4;
        if (x > 0.0002  && x <= 0.00035) bin[0] = 5;
        if (x > 0.00035 && x <= 0.0007 ) bin[0] = 6;
        if (x > 0.0007  && x <= 0.0020 ) bin[0] = 7;
      }
      else if (q2 > 10.*GeV2 && q2 <= 20.*GeV2) {
        if (x > 0.0002 && x <= 0.0005) bin[0] = 8;
        if (x > 0.0005 && x <= 0.0008) bin[0] = 9;
        if (x > 0.0008 && x <= 0.0015) bin[0] = 10;
        if (x > 0.0015 && x <= 0.040 ) bin[0] = 11;
      }
      else if (q2 > 20.*GeV2 && q2 <= 50.*GeV2) {
        if (x > 0.0005 && x <= 0.0014) bin[0] = 12;
        if (x > 0.0014 && x <= 0.0030) bin[0] = 13;
        if (x > 0.0030 && x <= 0.0100) bin[0] = 14;
      }
      else if (q2 > 50.*GeV2 && q2 <= 100.*GeV2) {
        if (x >0.0008 && x <= 0.0030) bin[0] = 15;
        if (x >0.0030 && x <= 0.0200) bin[0] = 16;
      }
      // check in one of the bins
      evcut[0] &= bin[0] >= 0;
      // For the low Q2 selection b)
      if (q2 > 2.5*GeV2 && q2 <= 5.  *GeV2) bin[1] = 0;
      if (q2 > 5. *GeV2 && q2 <= 10. *GeV2) bin[1] = 1;
      if (q2 > 10.*GeV2 && q2 <= 20. *GeV2) bin[1] = 2;
      if (q2 > 20.*GeV2 && q2 <= 50. *GeV2) bin[1] = 3;
      if (q2 > 50.*GeV2 && q2 <= 100.*GeV2) bin[1] = 4;
      // check in one of the bins
      evcut[1] &= bin[1] >= 0;
      // for the high Q2 selection a)
      if (q2 > 100.*GeV2 && q2 <= 400.*GeV2) {
        if (x > 0.00251 && x <= 0.00631) bin[2] = 0;
        if (x > 0.00631 && x <= 0.0158 ) bin[2] = 1;
        if (x > 0.0158  && x <= 0.0398 ) bin[2] = 2;
      }
      else if (q2 > 400.*GeV2 && q2 <= 1100.*GeV2) {
        if (x > 0.00631 && x <= 0.0158 ) bin[2] = 3;
        if (x > 0.0158  && x <= 0.0398 ) bin[2] = 4;
        if (x > 0.0398  && x <= 1.     ) bin[2] = 5;
      }
      else if (q2 > 1100.*GeV2 && q2 <= 100000.*GeV2) {
        if (x > 0. && x <= 1.) bin[2] = 6;
      }
      // check in one of the bins
      evcut[2] &= bin[2] >= 0;
      // for the high Q2 selection b)
      if      (q2 > 100.*GeV2 && q2 <= 220.*GeV2) bin[3] = 0;
      else if (q2 > 220.*GeV2 && q2 <= 400.*GeV2) bin[3] = 1;
      else if (q2 > 400.              ) bin[3] = 2;
      // check in one of*GeV the bins
      evcut[3] &= bin[3] >= 0;

      // Veto if fails all cuts after bin selection
      /// @todo Can we use all()?
      if (! (evcut[0] || evcut[1] || evcut[2] || evcut[3])) vetoEvent;

      // Increment the count for normalisation
      if (evcut[0]) _weightETLowQa [bin[0]]->fill();
      if (evcut[1]) _weightETLowQb [bin[1]]->fill();
      if (evcut[2]) _weightETHighQa[bin[2]]->fill();
      if (evcut[3]) _weightETHighQb[bin[3]]->fill();

      // Boost to hadronic CoM
      const LorentzTransform hcmboost = dk.boostHCM();

      // Loop over the particles
      double etcent = 0;
      double etfrag = 0;
      for (const Particle& p : particles) {
        // Boost momentum to CMS
        const FourMomentum hcmMom = hcmboost.transform(p.momentum());
        double et = fabs(hcmMom.Et());
        double eta = hcmMom.eta();
        // Averages in central and forward region
        if (fabs(eta) < .5 ) etcent += et;
        if (eta > 2 && eta <= 3.) etfrag += et;
        // Histograms of Et flow
        if (evcut[0]) _histETLowQa [bin[0]]->fill(eta, et);
        if (evcut[1]) _histETLowQb [bin[1]]->fill(eta, et);
        if (evcut[2]) _histETHighQa[bin[2]]->fill(eta, et);
        if (evcut[3]) _histETHighQb[bin[3]]->fill(eta, et);
      }
      // Fill histograms for the average quantities
      if (evcut[1] || evcut[3]) {
        _histAverETCentral->fill(q2, etcent);
        _histAverETFrag   ->fill(q2, etfrag);
      }
    }


    // Finalize
    void finalize() {
      // Normalization of the Et distributions
      /// @todo Simplify by using normalize() instead? Are all these being normalized to area=1?
      for (size_t ix = 0; ix < 17; ++ix) if (_weightETLowQa[ix]->val()  != 0) scale(_histETLowQa[ix],  1/ *_weightETLowQa[ix]);
      for (size_t ix = 0; ix <  7; ++ix) if (_weightETHighQa[ix]->val() != 0) scale(_histETHighQa[ix], 1/ *_weightETHighQa[ix]);
      for (size_t ix = 0; ix <  5; ++ix) if (_weightETLowQb[ix]->val()  != 0) scale(_histETLowQb[ix],  1/ *_weightETLowQb[ix]);
      for (size_t ix = 0; ix <  3; ++ix) if (_weightETHighQb[ix]->val() != 0) scale(_histETHighQb[ix], 1/ *_weightETHighQb[ix]);
    }

    //@}


  private:

    /// @name Histograms
    /// @{
    vector<Histo1DPtr> _histETLowQa;
    vector<Histo1DPtr> _histETHighQa;
    vector<Histo1DPtr> _histETLowQb;
    vector<Histo1DPtr> _histETHighQb;
    Profile1DPtr _histAverETCentral;
    Profile1DPtr _histAverETFrag;
    /// @}

    /// @name Storage of weights for normalisation
    /// @{
    array<CounterPtr,17> _weightETLowQa;
    array<CounterPtr, 7> _weightETHighQa;
    array<CounterPtr, 5> _weightETLowQb;
    array<CounterPtr, 3> _weightETHighQb;
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(H1_2000_S4129130, H1_2000_I503947);

}
