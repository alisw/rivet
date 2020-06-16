// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Sphericity.hh"

namespace Rivet {


  /// Rivet analysis class for ATLAS min bias event shapes
  class ATLAS_2012_I1124167 : public Analysis {
  public:

    /// Constructor
    ATLAS_2012_I1124167()
      : Analysis("ATLAS_2012_I1124167") {  }


    /// Initialization, called once before running
    void init() {
      // Projections
      ChargedFinalState cfs(Cuts::abseta < 2.5 && Cuts::pT > 0.5*GeV);
      declare(cfs, "CFS");

      // Book histograms
      book(_hist_T_05_25 ,1,1,1);
      book(_hist_T_05    ,2,1,1);
      book(_hist_T_25_50 ,1,1,2);
      book(_hist_T_25    ,2,1,2);
      book(_hist_T_50_75 ,1,1,3);
      book(_hist_T_50    ,2,1,3);
      book(_hist_T_75_100,1,1,4);
      book(_hist_T_75    ,2,1,4);
      book(_hist_T_100   ,2,1,5);

      book(_hist_TM_05_25 ,3,1,1);
      book(_hist_TM_05    ,4,1,1);
      book(_hist_TM_25_50 ,3,1,2);
      book(_hist_TM_25    ,4,1,2);
      book(_hist_TM_50_75 ,3,1,3);
      book(_hist_TM_50    ,4,1,3);
      book(_hist_TM_75_100,3,1,4);
      book(_hist_TM_75    ,4,1,4);
      book(_hist_TM_100   ,4,1,5);

      book(_hist_S_05_25 ,5,1,1);
      book(_hist_S_05    ,6,1,1);
      book(_hist_S_25_50 ,5,1,2);
      book(_hist_S_25    ,6,1,2);
      book(_hist_S_50_75 ,5,1,3);
      book(_hist_S_50    ,6,1,3);
      book(_hist_S_75_100,5,1,4);
      book(_hist_S_75    ,6,1,4);
      book(_hist_S_100   ,6,1,5);


      book(_hist_T_N  ,7,1,1);
      book(_hist_TM_N ,7,1,2);
      book(_hist_S_N  ,7,1,3);

      book(_hist_T_S  ,8,1,1);
      book(_hist_TM_S ,8,1,2);
      book(_hist_S_S  ,8,1,3);
    }


    void analyze(const Event& event) {
      const double weight = 1.0;

      // CFS projection and particles
      const Particles& particles500 = apply<ChargedFinalState>(event, "CFS").particlesByPt();

      // Require at least 6 charged particles
      if (particles500.size() < 6) vetoEvent;

      // Preparation for Thrust calculation
      vector<Vector3> momenta;

      // Counters
      double num500 = 0;
      double ptSum500 = 0;

      double pTlead = particles500[0].pT()/GeV;

      // Loop over particles
      for (const Particle& p : particles500) {
        num500 += 1;
        ptSum500 += p.pT()/GeV;

        // Transverse Thrust calculation requires p_z to be set to 0
        Vector3 mom = p.p3();
        mom.setZ(0.0);
        momenta.push_back(mom);
      }

      // If only 2 particles, we need to use a ghost so that Thrust.calc() doesn't return 1.
      if (momenta.size() == 2) {
        momenta.push_back(Vector3(1e-10*MeV, 0., 0.));
      }

      // Actual thrust calculation
      Thrust thrust;
      thrust.calc(momenta);

      const double T  = 1.0 - thrust.thrust();
      const double TM = thrust.thrustMajor();

      Sphericity sphericity;
      sphericity.calc(momenta);

      double S = sphericity.transSphericity();
      if ( std::isnan(S) )  S = -1.0; // put this in the underflow bin

      // Fill histos, most inclusive first

      // pTlead > 0.5
      _hist_T_05->fill(T , weight);
      _hist_TM_05->fill(TM, weight);
      _hist_S_05->fill(S , weight);

      // pTlead 0.5 - 2.5
      if (pTlead <= 2.5) {
        _hist_T_05_25->fill(T , weight);
        _hist_TM_05_25->fill(TM, weight);
        _hist_S_05_25->fill(S , weight);
      }

      // pTlead > 2.5
      if (pTlead > 2.5) {
        _hist_T_25->fill(T , weight);
        _hist_TM_25->fill(TM, weight);
        _hist_S_25->fill(S , weight);
      }

      // pTlead 2.5 - .5
      if (inRange(pTlead, 2.5, 5.0)) {
        _hist_T_25_50->fill(T , weight);
        _hist_TM_25_50->fill(TM, weight);
        _hist_S_25_50->fill(S , weight);
      }

      // pTlead > 5
      if (pTlead > 5) {
        _hist_T_50->fill(T , weight);
        _hist_TM_50->fill(TM, weight);
        _hist_S_50->fill(S , weight);
      }

      // pTlead 5 - 7.5
      if (inRange(pTlead, 5.0, 7.5)) {
        _hist_T_50_75->fill(T , weight);
        _hist_TM_50_75->fill(TM, weight);
        _hist_S_50_75->fill(S , weight);
      }

      // pTlead > 7.5
      if (pTlead > 7.5) {
        _hist_T_75->fill(T , weight);
        _hist_TM_75->fill(TM, weight);
        _hist_S_75->fill(S , weight);
      }

      // pTlead 7.5 - 10
      if (inRange(pTlead, 7.5, 10)) {
        _hist_T_75_100->fill(T , weight);
        _hist_TM_75_100->fill(TM, weight);
        _hist_S_75_100->fill(S , weight);
      }

      // pTlead > 10
      if (pTlead > 10) {
        _hist_T_100->fill(T , weight);
        _hist_TM_100->fill(TM, weight);
        _hist_S_100->fill(S , weight);
      }


      // Profiles Nch vs. ES
      _hist_T_N->fill(num500, T, weight);
      _hist_TM_N->fill(num500, TM, weight);
      _hist_S_N->fill(num500, S, weight);

      // Profiles pTsum vs. ES
      _hist_T_S->fill(ptSum500, T, weight);
      _hist_TM_S->fill(ptSum500, TM, weight);
      _hist_S_S->fill(ptSum500, S, weight);
    }


    void finalize() {
      normalize(_hist_T_05_25);
      normalize(_hist_T_05);
      normalize(_hist_T_25_50);
      normalize(_hist_T_25);
      normalize(_hist_T_50_75);
      normalize(_hist_T_50);
      normalize(_hist_T_75_100);
      normalize(_hist_T_75);
      normalize(_hist_T_100);

      normalize(_hist_TM_05_25);
      normalize(_hist_TM_05);
      normalize(_hist_TM_25_50);
      normalize(_hist_TM_25);
      normalize(_hist_TM_50_75);
      normalize(_hist_TM_50);
      normalize(_hist_TM_75_100);
      normalize(_hist_TM_75);
      normalize(_hist_TM_100);

      normalize(_hist_S_05_25);
      normalize(_hist_S_05);
      normalize(_hist_S_25_50);
      normalize(_hist_S_25);
      normalize(_hist_S_50_75);
      normalize(_hist_S_50);
      normalize(_hist_S_75_100);
      normalize(_hist_S_75);
      normalize(_hist_S_100);
    }


  private:

    Histo1DPtr _hist_T_05_25, _hist_T_05, _hist_T_25_50, _hist_T_25, _hist_T_50_75, _hist_T_50, _hist_T_75_100, _hist_T_75, _hist_T_100;
    Histo1DPtr _hist_TM_05_25, _hist_TM_05, _hist_TM_25_50, _hist_TM_25, _hist_TM_50_75, _hist_TM_50, _hist_TM_75_100, _hist_TM_75, _hist_TM_100;
    Histo1DPtr _hist_S_05_25, _hist_S_05, _hist_S_25_50, _hist_S_25, _hist_S_50_75, _hist_S_50, _hist_S_75_100, _hist_S_75, _hist_S_100;
    Profile1DPtr _hist_T_N, _hist_TM_N, _hist_S_N;
    Profile1DPtr _hist_T_S, _hist_TM_S, _hist_S_S;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1124167);

}
