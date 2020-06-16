// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  namespace {

    inline double calcz(const Jet& j, const Particle& p) {
      const double num = j.p3().dot(p.p3());
      const double den = j.p3().mod2();
      return num/den;
    }

    inline double calcptrel(const Jet& j, const Particle& p) {
      const double num = j.p3().cross(p.p3()).mod();
      const double den = j.p3().mod();
      return num/den;
    }

    inline double calcr(const Jet& j, const Particle& p) {
      return deltaR(j.rapidity(), j.phi(), p.rapidity(), p.phi());
    }

    // For annulus area kludge
    /// @todo Improve somehow... need normalisation *without* bin width factors!
    inline double calcrweight(const Jet& j, const Particle& p) {
      size_t nBins_r = 26;
      double bins_r[] = { 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10,
                          0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30,
                          0.35, 0.40, 0.45, 0.50, 0.55, 0.60 };
      double r = calcr(j,p);
      for (size_t bin = 0 ; bin < nBins_r ; bin++) {
        if (r < bins_r[bin+1]) {
          double up = bins_r[bin+1];
          double down = bins_r[bin];
          return ((up-down)/(M_PI*(up*up-down*down)));
        }
      }
      return 1.0;
    }
  }


  class ATLAS_2011_I919017 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    ATLAS_2011_I919017()
      : Analysis("ATLAS_2011_I919017")
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      ChargedFinalState cfs((Cuts::etaIn(-2.5, 2.5) && Cuts::pT >=  0.3*GeV));
      FastJets trkjets04(cfs, FastJets::ANTIKT, 0.4);
      FastJets trkjets06(cfs, FastJets::ANTIKT, 0.6);
      declare(trkjets04, "Jets04");
      declare(trkjets06, "Jets06");

      // Book histograms
      book(_h_pt04_00_05 ,1, 1, 1);
      book(_h_pt06_00_05 ,2, 1, 1);
      book(_h_N04_00_05_04_06 ,1, 2, 1+5);
      book(_h_N06_00_05_04_06 ,2, 2, 1+5);
      book(_h_N04_00_05_06_10 ,1, 2, 2+5);
      book(_h_N06_00_05_06_10 ,2, 2, 2+5);
      book(_h_N04_00_05_10_15 ,1, 2, 3+5);
      book(_h_N06_00_05_10_15 ,2, 2, 3+5);
      book(_h_N04_00_05_15_24 ,1, 2, 4+5);
      book(_h_N06_00_05_15_24 ,2, 2, 4+5);
      book(_h_N04_00_05_24_40 ,1, 2, 5+5);
      book(_h_N06_00_05_24_40 ,2, 2, 5+5);
      book(_h_z04_00_05_04_06 ,1, 3, 1+5);
      book(_h_z06_00_05_04_06 ,2, 3, 1+5);
      book(_h_z04_00_05_06_10 ,1, 3, 2+5);
      book(_h_z06_00_05_06_10 ,2, 3, 2+5);
      book(_h_z04_00_05_10_15 ,1, 3, 3+5);
      book(_h_z06_00_05_10_15 ,2, 3, 3+5);
      book(_h_z04_00_05_15_24 ,1, 3, 4+5);
      book(_h_z06_00_05_15_24 ,2, 3, 4+5);
      book(_h_z04_00_05_24_40 ,1, 3, 5+5);
      book(_h_z06_00_05_24_40 ,2, 3, 5+5);
      book(_h_ptrel04_00_05_04_06 ,1, 4, 1+5);
      book(_h_ptrel06_00_05_04_06 ,2, 4, 1+5);
      book(_h_ptrel04_00_05_06_10 ,1, 4, 2+5);
      book(_h_ptrel06_00_05_06_10 ,2, 4, 2+5);
      book(_h_ptrel04_00_05_10_15 ,1, 4, 3+5);
      book(_h_ptrel06_00_05_10_15 ,2, 4, 3+5);
      book(_h_ptrel04_00_05_15_24 ,1, 4, 4+5);
      book(_h_ptrel06_00_05_15_24 ,2, 4, 4+5);
      book(_h_ptrel04_00_05_24_40 ,1, 4, 5+5);
      book(_h_ptrel06_00_05_24_40 ,2, 4, 5+5);
      book(_h_rdA04_00_05_04_06 ,1, 5, 1+5);
      book(_h_rdA06_00_05_04_06 ,2, 5, 1+5);
      book(_h_rdA04_00_05_06_10 ,1, 5, 2+5);
      book(_h_rdA06_00_05_06_10 ,2, 5, 2+5);
      book(_h_rdA04_00_05_10_15 ,1, 5, 3+5);
      book(_h_rdA06_00_05_10_15 ,2, 5, 3+5);
      book(_h_rdA04_00_05_15_24 ,1, 5, 4+5);
      book(_h_rdA06_00_05_15_24 ,2, 5, 4+5);
      book(_h_rdA04_00_05_24_40 ,1, 5, 5+5);
      book(_h_rdA06_00_05_24_40 ,2, 5, 5+5);

      book(_h_pt04_05_10 ,1, 1, 2);
      book(_h_pt06_05_10 ,2, 1, 2);
      book(_h_N04_05_10_04_06 ,1, 2, 1+10);
      book(_h_N06_05_10_04_06 ,2, 2, 1+10);
      book(_h_N04_05_10_06_10 ,1, 2, 2+10);
      book(_h_N06_05_10_06_10 ,2, 2, 2+10);
      book(_h_N04_05_10_10_15 ,1, 2, 3+10);
      book(_h_N06_05_10_10_15 ,2, 2, 3+10);
      book(_h_N04_05_10_15_24 ,1, 2, 4+10);
      book(_h_N06_05_10_15_24 ,2, 2, 4+10);
      book(_h_N04_05_10_24_40 ,1, 2, 5+10);
      book(_h_N06_05_10_24_40 ,2, 2, 5+10);
      book(_h_z04_05_10_04_06 ,1, 3, 1+10);
      book(_h_z06_05_10_04_06 ,2, 3, 1+10);
      book(_h_z04_05_10_06_10 ,1, 3, 2+10);
      book(_h_z06_05_10_06_10 ,2, 3, 2+10);
      book(_h_z04_05_10_10_15 ,1, 3, 3+10);
      book(_h_z06_05_10_10_15 ,2, 3, 3+10);
      book(_h_z04_05_10_15_24 ,1, 3, 4+10);
      book(_h_z06_05_10_15_24 ,2, 3, 4+10);
      book(_h_z04_05_10_24_40 ,1, 3, 5+10);
      book(_h_z06_05_10_24_40 ,2, 3, 5+10);
      book(_h_ptrel04_05_10_04_06 ,1, 4, 1+10);
      book(_h_ptrel06_05_10_04_06 ,2, 4, 1+10);
      book(_h_ptrel04_05_10_06_10 ,1, 4, 2+10);
      book(_h_ptrel06_05_10_06_10 ,2, 4, 2+10);
      book(_h_ptrel04_05_10_10_15 ,1, 4, 3+10);
      book(_h_ptrel06_05_10_10_15 ,2, 4, 3+10);
      book(_h_ptrel04_05_10_15_24 ,1, 4, 4+10);
      book(_h_ptrel06_05_10_15_24 ,2, 4, 4+10);
      book(_h_ptrel04_05_10_24_40 ,1, 4, 5+10);
      book(_h_ptrel06_05_10_24_40 ,2, 4, 5+10);
      book(_h_rdA04_05_10_04_06 ,1, 5, 1+10);
      book(_h_rdA06_05_10_04_06 ,2, 5, 1+10);
      book(_h_rdA04_05_10_06_10 ,1, 5, 2+10);
      book(_h_rdA06_05_10_06_10 ,2, 5, 2+10);
      book(_h_rdA04_05_10_10_15 ,1, 5, 3+10);
      book(_h_rdA06_05_10_10_15 ,2, 5, 3+10);
      book(_h_rdA04_05_10_15_24 ,1, 5, 4+10);
      book(_h_rdA06_05_10_15_24 ,2, 5, 4+10);
      book(_h_rdA04_05_10_24_40 ,1, 5, 5+10);
      book(_h_rdA06_05_10_24_40 ,2, 5, 5+10);

      book(_h_pt04_10_15 ,1, 1, 3);
      book(_h_pt06_10_15 ,2, 1, 3);
      book(_h_N04_10_15_04_06 ,1, 2, 1+15);
      book(_h_N06_10_15_04_06 ,2, 2, 1+15);
      book(_h_N04_10_15_06_10 ,1, 2, 2+15);
      book(_h_N06_10_15_06_10 ,2, 2, 2+15);
      book(_h_N04_10_15_10_15 ,1, 2, 3+15);
      book(_h_N06_10_15_10_15 ,2, 2, 3+15);
      book(_h_N04_10_15_15_24 ,1, 2, 4+15);
      book(_h_N06_10_15_15_24 ,2, 2, 4+15);
      book(_h_N04_10_15_24_40 ,1, 2, 5+15);
      book(_h_N06_10_15_24_40 ,2, 2, 5+15);
      book(_h_z04_10_15_04_06 ,1, 3, 1+15);
      book(_h_z06_10_15_04_06 ,2, 3, 1+15);
      book(_h_z04_10_15_06_10 ,1, 3, 2+15);
      book(_h_z06_10_15_06_10 ,2, 3, 2+15);
      book(_h_z04_10_15_10_15 ,1, 3, 3+15);
      book(_h_z06_10_15_10_15 ,2, 3, 3+15);
      book(_h_z04_10_15_15_24 ,1, 3, 4+15);
      book(_h_z06_10_15_15_24 ,2, 3, 4+15);
      book(_h_z04_10_15_24_40 ,1, 3, 5+15);
      book(_h_z06_10_15_24_40 ,2, 3, 5+15);
      book(_h_ptrel04_10_15_04_06 ,1, 4, 1+15);
      book(_h_ptrel06_10_15_04_06 ,2, 4, 1+15);
      book(_h_ptrel04_10_15_06_10 ,1, 4, 2+15);
      book(_h_ptrel06_10_15_06_10 ,2, 4, 2+15);
      book(_h_ptrel04_10_15_10_15 ,1, 4, 3+15);
      book(_h_ptrel06_10_15_10_15 ,2, 4, 3+15);
      book(_h_ptrel04_10_15_15_24 ,1, 4, 4+15);
      book(_h_ptrel06_10_15_15_24 ,2, 4, 4+15);
      book(_h_ptrel04_10_15_24_40 ,1, 4, 5+15);
      book(_h_ptrel06_10_15_24_40 ,2, 4, 5+15);
      book(_h_rdA04_10_15_04_06 ,1, 5, 1+15);
      book(_h_rdA06_10_15_04_06 ,2, 5, 1+15);
      book(_h_rdA04_10_15_06_10 ,1, 5, 2+15);
      book(_h_rdA06_10_15_06_10 ,2, 5, 2+15);
      book(_h_rdA04_10_15_10_15 ,1, 5, 3+15);
      book(_h_rdA06_10_15_10_15 ,2, 5, 3+15);
      book(_h_rdA04_10_15_15_24 ,1, 5, 4+15);
      book(_h_rdA06_10_15_15_24 ,2, 5, 4+15);
      book(_h_rdA04_10_15_24_40 ,1, 5, 5+15);
      book(_h_rdA06_10_15_24_40 ,2, 5, 5+15);

      book(_h_pt04_15_19 ,1, 1, 4);
      book(_h_pt06_15_19 ,2, 1, 4);
      book(_h_N04_15_19_04_06 ,1, 2, 1+20);
      book(_h_N06_15_19_04_06 ,2, 2, 1+20);
      book(_h_N04_15_19_06_10 ,1, 2, 2+20);
      book(_h_N06_15_19_06_10 ,2, 2, 2+20);
      book(_h_N04_15_19_10_15 ,1, 2, 3+20);
      book(_h_N06_15_19_10_15 ,2, 2, 3+20);
      book(_h_N04_15_19_15_24 ,1, 2, 4+20);
      book(_h_N06_15_19_15_24 ,2, 2, 4+20);
      book(_h_N04_15_19_24_40 ,1, 2, 5+20);
      book(_h_N06_15_19_24_40 ,2, 2, 5+20);
      book(_h_z04_15_19_04_06 ,1, 3, 1+20);
      book(_h_z06_15_19_04_06 ,2, 3, 1+20);
      book(_h_z04_15_19_06_10 ,1, 3, 2+20);
      book(_h_z06_15_19_06_10 ,2, 3, 2+20);
      book(_h_z04_15_19_10_15 ,1, 3, 3+20);
      book(_h_z06_15_19_10_15 ,2, 3, 3+20);
      book(_h_z04_15_19_15_24 ,1, 3, 4+20);
      book(_h_z06_15_19_15_24 ,2, 3, 4+20);
      book(_h_z04_15_19_24_40 ,1, 3, 5+20);
      book(_h_z06_15_19_24_40 ,2, 3, 5+20);
      book(_h_ptrel04_15_19_04_06 ,1, 4, 1+20);
      book(_h_ptrel06_15_19_04_06 ,2, 4, 1+20);
      book(_h_ptrel04_15_19_06_10 ,1, 4, 2+20);
      book(_h_ptrel06_15_19_06_10 ,2, 4, 2+20);
      book(_h_ptrel04_15_19_10_15 ,1, 4, 3+20);
      book(_h_ptrel06_15_19_10_15 ,2, 4, 3+20);
      book(_h_ptrel04_15_19_15_24 ,1, 4, 4+20);
      book(_h_ptrel06_15_19_15_24 ,2, 4, 4+20);
      book(_h_ptrel04_15_19_24_40 ,1, 4, 5+20);
      book(_h_ptrel06_15_19_24_40 ,2, 4, 5+20);
      book(_h_rdA04_15_19_04_06 ,1, 5, 1+20);
      book(_h_rdA06_15_19_04_06 ,2, 5, 1+20);
      book(_h_rdA04_15_19_06_10 ,1, 5, 2+20);
      book(_h_rdA06_15_19_06_10 ,2, 5, 2+20);
      book(_h_rdA04_15_19_10_15 ,1, 5, 3+20);
      book(_h_rdA06_15_19_10_15 ,2, 5, 3+20);
      book(_h_rdA04_15_19_15_24 ,1, 5, 4+20);
      book(_h_rdA06_15_19_15_24 ,2, 5, 4+20);
      book(_h_rdA04_15_19_24_40 ,1, 5, 5+20);
      book(_h_rdA06_15_19_24_40 ,2, 5, 5+20);

      book(_h_N04_00_19_04_06 ,1, 2, 1+0);
      book(_h_N06_00_19_04_06 ,2, 2, 1+0);
      book(_h_N04_00_19_06_10 ,1, 2, 2+0);
      book(_h_N06_00_19_06_10 ,2, 2, 2+0);
      book(_h_N04_00_19_10_15 ,1, 2, 3+0);
      book(_h_N06_00_19_10_15 ,2, 2, 3+0);
      book(_h_N04_00_19_15_24 ,1, 2, 4+0);
      book(_h_N06_00_19_15_24 ,2, 2, 4+0);
      book(_h_N04_00_19_24_40 ,1, 2, 5+0);
      book(_h_N06_00_19_24_40 ,2, 2, 5+0);
      book(_h_z04_00_19_04_06 ,1, 3, 1+0);
      book(_h_z06_00_19_04_06 ,2, 3, 1+0);
      book(_h_z04_00_19_06_10 ,1, 3, 2+0);
      book(_h_z06_00_19_06_10 ,2, 3, 2+0);
      book(_h_z04_00_19_10_15 ,1, 3, 3+0);
      book(_h_z06_00_19_10_15 ,2, 3, 3+0);
      book(_h_z04_00_19_15_24 ,1, 3, 4+0);
      book(_h_z06_00_19_15_24 ,2, 3, 4+0);
      book(_h_z04_00_19_24_40 ,1, 3, 5+0);
      book(_h_z06_00_19_24_40 ,2, 3, 5+0);
      book(_h_ptrel04_00_19_04_06 ,1, 4, 1+0);
      book(_h_ptrel06_00_19_04_06 ,2, 4, 1+0);
      book(_h_ptrel04_00_19_06_10 ,1, 4, 2+0);
      book(_h_ptrel06_00_19_06_10 ,2, 4, 2+0);
      book(_h_ptrel04_00_19_10_15 ,1, 4, 3+0);
      book(_h_ptrel06_00_19_10_15 ,2, 4, 3+0);
      book(_h_ptrel04_00_19_15_24 ,1, 4, 4+0);
      book(_h_ptrel06_00_19_15_24 ,2, 4, 4+0);
      book(_h_ptrel04_00_19_24_40 ,1, 4, 5+0);
      book(_h_ptrel06_00_19_24_40 ,2, 4, 5+0);
      book(_h_rdA04_00_19_04_06 ,1, 5, 1+0);
      book(_h_rdA06_00_19_04_06 ,2, 5, 1+0);
      book(_h_rdA04_00_19_06_10 ,1, 5, 2+0);
      book(_h_rdA06_00_19_06_10 ,2, 5, 2+0);
      book(_h_rdA04_00_19_10_15 ,1, 5, 3+0);
      book(_h_rdA06_00_19_10_15 ,2, 5, 3+0);
      book(_h_rdA04_00_19_15_24 ,1, 5, 4+0);
      book(_h_rdA06_00_19_15_24 ,2, 5, 4+0);
      book(_h_rdA04_00_19_24_40 ,1, 5, 5+0);
      book(_h_rdA06_00_19_24_40 ,2, 5, 5+0);

      book(_sumofweights04, "_sumofweights04");
      book(_sumofweights06, "_sumofweights06");
      book(_numjets04_00_05_04_06, "_numjets04_00_05_04_06");
      book(_numjets04_00_05_06_10, "_numjets04_00_05_06_10");
      book(_numjets04_00_05_10_15, "_numjets04_00_05_10_15");
      book(_numjets04_00_05_15_24, "_numjets04_00_05_15_24");
      book(_numjets04_00_05_24_40, "_numjets04_00_05_24_40");
      book(_numjets06_00_05_04_06, "_numjets06_00_05_04_06");
      book(_numjets06_00_05_06_10, "_numjets06_00_05_06_10");
      book(_numjets06_00_05_10_15, "_numjets06_00_05_10_15");
      book(_numjets06_00_05_15_24, "_numjets06_00_05_15_24");
      book(_numjets06_00_05_24_40, "_numjets06_00_05_24_40");
      book(_numjets04_05_10_04_06, "_numjets04_05_10_04_06");
      book(_numjets04_05_10_06_10, "_numjets04_05_10_06_10");
      book(_numjets04_05_10_10_15, "_numjets04_05_10_10_15");
      book(_numjets04_05_10_15_24, "_numjets04_05_10_15_24");
      book(_numjets04_05_10_24_40, "_numjets04_05_10_24_40");
      book(_numjets06_05_10_04_06, "_numjets06_05_10_04_06");
      book(_numjets06_05_10_06_10, "_numjets06_05_10_06_10");
      book(_numjets06_05_10_10_15, "_numjets06_05_10_10_15");
      book(_numjets06_05_10_15_24, "_numjets06_05_10_15_24");
      book(_numjets06_05_10_24_40, "_numjets06_05_10_24_40");
      book(_numjets04_10_15_04_06, "_numjets04_10_15_04_06");
      book(_numjets04_10_15_06_10, "_numjets04_10_15_06_10");
      book(_numjets04_10_15_10_15, "_numjets04_10_15_10_15");
      book(_numjets04_10_15_15_24, "_numjets04_10_15_15_24");
      book(_numjets04_10_15_24_40, "_numjets04_10_15_24_40");
      book(_numjets06_10_15_04_06, "_numjets06_10_15_04_06");
      book(_numjets06_10_15_06_10, "_numjets06_10_15_06_10");
      book(_numjets06_10_15_10_15, "_numjets06_10_15_10_15");
      book(_numjets06_10_15_15_24, "_numjets06_10_15_15_24");
      book(_numjets06_10_15_24_40, "_numjets06_10_15_24_40");
      book(_numjets04_15_19_04_06, "_numjets04_15_19_04_06");
      book(_numjets04_15_19_06_10, "_numjets04_15_19_06_10");
      book(_numjets04_15_19_10_15, "_numjets04_15_19_10_15");
      book(_numjets04_15_19_15_24, "_numjets04_15_19_15_24");
      book(_numjets04_15_19_24_40, "_numjets04_15_19_24_40");
      book(_numjets06_15_19_04_06, "_numjets06_15_19_04_06");
      book(_numjets06_15_19_06_10, "_numjets06_15_19_06_10");
      book(_numjets06_15_19_10_15, "_numjets06_15_19_10_15");
      book(_numjets06_15_19_15_24, "_numjets06_15_19_15_24");
      book(_numjets06_15_19_24_40, "_numjets06_15_19_24_40");
      book(_numjets04_00_19_04_06, "_numjets04_00_19_04_06");
      book(_numjets04_00_19_06_10, "_numjets04_00_19_06_10");
      book(_numjets04_00_19_10_15, "_numjets04_00_19_10_15");
      book(_numjets04_00_19_15_24, "_numjets04_00_19_15_24");
      book(_numjets04_00_19_24_40, "_numjets04_00_19_24_40");
      book(_numjets06_00_19_04_06, "_numjets06_00_19_04_06");
      book(_numjets06_00_19_06_10, "_numjets06_00_19_06_10");
      book(_numjets06_00_19_10_15, "_numjets06_00_19_10_15");
      book(_numjets06_00_19_15_24, "_numjets06_00_19_15_24");
      book(_numjets06_00_19_24_40, "_numjets06_00_19_24_40");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const Jets& jets04 = apply<JetAlg>(event, "Jets04").jets();
      if (!jets04.empty()) {
        _sumofweights04->fill();
        for (const Jet& j : jets04) {
          const double jetpt = j.pT();
          if (j.absrap() < 0.5) {
            _h_pt04_00_05->fill(jetpt/GeV);
            if (inRange(jetpt/GeV, 4., 6.)) {
              _numjets04_00_05_04_06->fill();
              _h_N04_00_05_04_06->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z04_00_05_04_06->fill(calcz(j,p));
                _h_ptrel04_00_05_04_06->fill(calcptrel(j,p));
                _h_rdA04_00_05_04_06->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 6., 10.)) {
              _numjets04_00_05_06_10->fill();
              _h_N04_00_05_06_10->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z04_00_05_06_10->fill(calcz(j,p));
                _h_ptrel04_00_05_06_10->fill(calcptrel(j,p));
                _h_rdA04_00_05_06_10->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 10., 15.)) {
              _numjets04_00_05_10_15->fill();
              _h_N04_00_05_10_15->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z04_00_05_10_15->fill(calcz(j,p));
                _h_ptrel04_00_05_10_15->fill(calcptrel(j,p));
                _h_rdA04_00_05_10_15->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 15., 24.)) {
              _numjets04_00_05_15_24->fill();
              _h_N04_00_05_15_24->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z04_00_05_15_24->fill(calcz(j,p));
                _h_ptrel04_00_05_15_24->fill(calcptrel(j,p));
                _h_rdA04_00_05_15_24->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 24., 40.)) {
              _numjets04_00_05_24_40->fill();
              _h_N04_00_05_24_40->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z04_00_05_24_40->fill(calcz(j,p));
                _h_ptrel04_00_05_24_40->fill(calcptrel(j,p));
                _h_rdA04_00_05_24_40->fill(calcr(j,p),calcrweight(j,p));
              }
            }
          }
          if (j.absrap() > 0.5 && j.absrap() < 1.0) {
            _h_pt04_05_10->fill(jetpt/GeV);
            if (inRange(jetpt/GeV, 4., 6.)) {
              _numjets04_05_10_04_06->fill();
              _h_N04_05_10_04_06->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z04_05_10_04_06->fill(calcz(j,p));
                _h_ptrel04_05_10_04_06->fill(calcptrel(j,p));
                _h_rdA04_05_10_04_06->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 6., 10.)) {
              _numjets04_05_10_06_10->fill();
              _h_N04_05_10_06_10->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z04_05_10_06_10->fill(calcz(j,p));
                _h_ptrel04_05_10_06_10->fill(calcptrel(j,p));
                _h_rdA04_05_10_06_10->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 10., 15.)) {
              _numjets04_05_10_10_15->fill();
              _h_N04_05_10_10_15->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z04_05_10_10_15->fill(calcz(j,p));
                _h_ptrel04_05_10_10_15->fill(calcptrel(j,p));
                _h_rdA04_05_10_10_15->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 15., 24.)) {
              _numjets04_05_10_15_24->fill();
              _h_N04_05_10_15_24->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z04_05_10_15_24->fill(calcz(j,p));
                _h_ptrel04_05_10_15_24->fill(calcptrel(j,p));
                _h_rdA04_05_10_15_24->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 24., 40.)) {
              _numjets04_05_10_24_40->fill();
              _h_N04_05_10_24_40->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z04_05_10_24_40->fill(calcz(j,p));
                _h_ptrel04_05_10_24_40->fill(calcptrel(j,p));
                _h_rdA04_05_10_24_40->fill(calcr(j,p),calcrweight(j,p));
              }
            }
          }
          if (j.absrap() > 1.0 && j.absrap() < 1.5) {
            _h_pt04_10_15->fill(jetpt/GeV);
            if (inRange(jetpt/GeV, 4., 6.)) {
              _numjets04_10_15_04_06->fill();
              _h_N04_10_15_04_06->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z04_10_15_04_06->fill(calcz(j,p));
                _h_ptrel04_10_15_04_06->fill(calcptrel(j,p));
                _h_rdA04_10_15_04_06->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 6., 10.)) {
              _numjets04_10_15_06_10->fill();
              _h_N04_10_15_06_10->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z04_10_15_06_10->fill(calcz(j,p));
                _h_ptrel04_10_15_06_10->fill(calcptrel(j,p));
                _h_rdA04_10_15_06_10->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 10., 15.)) {
              _numjets04_10_15_10_15->fill();
              _h_N04_10_15_10_15->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z04_10_15_10_15->fill(calcz(j,p));
                _h_ptrel04_10_15_10_15->fill(calcptrel(j,p));
                _h_rdA04_10_15_10_15->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 15., 24.)) {
              _numjets04_10_15_15_24->fill();
              _h_N04_10_15_15_24->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z04_10_15_15_24->fill(calcz(j,p));
                _h_ptrel04_10_15_15_24->fill(calcptrel(j,p));
                _h_rdA04_10_15_15_24->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 24., 40.)) {
              _numjets04_10_15_24_40->fill();
              _h_N04_10_15_24_40->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z04_10_15_24_40->fill(calcz(j,p));
                _h_ptrel04_10_15_24_40->fill(calcptrel(j,p));
                _h_rdA04_10_15_24_40->fill(calcr(j,p),calcrweight(j,p));
              }
            }
          }
          if (j.absrap() > 1.5 && j.absrap() < 1.9) {
            _h_pt04_15_19->fill(jetpt/GeV);
            if (inRange(jetpt/GeV, 4., 6.)) {
              _numjets04_15_19_04_06->fill();
              _h_N04_15_19_04_06->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z04_15_19_04_06->fill(calcz(j,p));
                _h_ptrel04_15_19_04_06->fill(calcptrel(j,p));
                _h_rdA04_15_19_04_06->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 6., 10.)) {
              _numjets04_15_19_06_10->fill();
              _h_N04_15_19_06_10->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z04_15_19_06_10->fill(calcz(j,p));
                _h_ptrel04_15_19_06_10->fill(calcptrel(j,p));
                _h_rdA04_15_19_06_10->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 10., 15.)) {
              _numjets04_15_19_10_15->fill();
              _h_N04_15_19_10_15->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z04_15_19_10_15->fill(calcz(j,p));
                _h_ptrel04_15_19_10_15->fill(calcptrel(j,p));
                _h_rdA04_15_19_10_15->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 15., 24.)) {
              _numjets04_15_19_15_24->fill();
              _h_N04_15_19_15_24->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z04_15_19_15_24->fill(calcz(j,p));
                _h_ptrel04_15_19_15_24->fill(calcptrel(j,p));
                _h_rdA04_15_19_15_24->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 24., 40.)) {
              _numjets04_15_19_24_40->fill();
              _h_N04_15_19_24_40->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z04_15_19_24_40->fill(calcz(j,p));
                _h_ptrel04_15_19_24_40->fill(calcptrel(j,p));
                _h_rdA04_15_19_24_40->fill(calcr(j,p),calcrweight(j,p));
              }
            }
          } // 1.5 < rapidity < 1.9
          if (j.absrap() < 1.9) {
            if (inRange(jetpt/GeV, 4., 6.)) {
              _numjets04_00_19_04_06->fill();
              _h_N04_00_19_04_06->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z04_00_19_04_06->fill(calcz(j,p));
                _h_ptrel04_00_19_04_06->fill(calcptrel(j,p));
                _h_rdA04_00_19_04_06->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 6., 10.)) {
              _numjets04_00_19_06_10->fill();
              _h_N04_00_19_06_10->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z04_00_19_06_10->fill(calcz(j,p));
                _h_ptrel04_00_19_06_10->fill(calcptrel(j,p));
                _h_rdA04_00_19_06_10->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 10., 15.)) {
              _numjets04_00_19_10_15->fill();
              _h_N04_00_19_10_15->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z04_00_19_10_15->fill(calcz(j,p));
                _h_ptrel04_00_19_10_15->fill(calcptrel(j,p));
                _h_rdA04_00_19_10_15->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 15., 24.)) {
              _numjets04_00_19_15_24->fill();
              _h_N04_00_19_15_24->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z04_00_19_15_24->fill(calcz(j,p));
                _h_ptrel04_00_19_15_24->fill(calcptrel(j,p));
                _h_rdA04_00_19_15_24->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 24., 40.)) {
              _numjets04_00_19_24_40->fill();
              _h_N04_00_19_24_40->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z04_00_19_24_40->fill(calcz(j,p));
                _h_ptrel04_00_19_24_40->fill(calcptrel(j,p));
                _h_rdA04_00_19_24_40->fill(calcr(j,p),calcrweight(j,p));
              }
            }
          } // 0.0 < rapidity < 1.9
        } // each jet
      } // jets04 not empty

      const Jets& jets06 = apply<JetAlg>(event, "Jets06").jets();
      if (!jets06.empty()) {
        _sumofweights06->fill();
        for (const Jet& j : jets06) {
          const double jetpt = j.pT();
          if (j.absrap() < 0.5) {
            _h_pt06_00_05->fill(jetpt/GeV);
            if (inRange(jetpt/GeV, 4., 6.)) {
              _numjets06_00_05_04_06->fill();
              _h_N06_00_05_04_06->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z06_00_05_04_06->fill(calcz(j,p));
                _h_ptrel06_00_05_04_06->fill(calcptrel(j,p));
                _h_rdA06_00_05_04_06->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 6., 10.)) {
              _numjets06_00_05_06_10->fill();
              _h_N06_00_05_06_10->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z06_00_05_06_10->fill(calcz(j,p));
                _h_ptrel06_00_05_06_10->fill(calcptrel(j,p));
                _h_rdA06_00_05_06_10->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 10., 15.)) {
              _numjets06_00_05_10_15->fill();
              _h_N06_00_05_10_15->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z06_00_05_10_15->fill(calcz(j,p));
                _h_ptrel06_00_05_10_15->fill(calcptrel(j,p));
                _h_rdA06_00_05_10_15->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 15., 24.)) {
              _numjets06_00_05_15_24->fill();
              _h_N06_00_05_15_24->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z06_00_05_15_24->fill(calcz(j,p));
                _h_ptrel06_00_05_15_24->fill(calcptrel(j,p));
                _h_rdA06_00_05_15_24->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 24., 40.)) {
              _numjets06_00_05_24_40->fill();
              _h_N06_00_05_24_40->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z06_00_05_24_40->fill(calcz(j,p));
                _h_ptrel06_00_05_24_40->fill(calcptrel(j,p));
                _h_rdA06_00_05_24_40->fill(calcr(j,p),calcrweight(j,p));
              }
            }
          }
          if (j.absrap() > 0.5 && j.absrap() < 1.0) {
            _h_pt06_05_10->fill(jetpt/GeV);
            if (inRange(jetpt/GeV, 4., 6.)) {
              _numjets06_05_10_04_06->fill();
              _h_N06_05_10_04_06->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z06_05_10_04_06->fill(calcz(j,p));
                _h_ptrel06_05_10_04_06->fill(calcptrel(j,p));
                _h_rdA06_05_10_04_06->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 6., 10.)) {
              _numjets06_05_10_06_10->fill();
              _h_N06_05_10_06_10->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z06_05_10_06_10->fill(calcz(j,p));
                _h_ptrel06_05_10_06_10->fill(calcptrel(j,p));
                _h_rdA06_05_10_06_10->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 10., 15.)) {
              _numjets06_05_10_10_15->fill();
              _h_N06_05_10_10_15->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z06_05_10_10_15->fill(calcz(j,p));
                _h_ptrel06_05_10_10_15->fill(calcptrel(j,p));
                _h_rdA06_05_10_10_15->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 15., 24.)) {
              _numjets06_05_10_15_24->fill();
              _h_N06_05_10_15_24->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z06_05_10_15_24->fill(calcz(j,p));
                _h_ptrel06_05_10_15_24->fill(calcptrel(j,p));
                _h_rdA06_05_10_15_24->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 24., 40.)) {
              _numjets06_05_10_24_40->fill();
              _h_N06_05_10_24_40->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z06_05_10_24_40->fill(calcz(j,p));
                _h_ptrel06_05_10_24_40->fill(calcptrel(j,p));
                _h_rdA06_05_10_24_40->fill(calcr(j,p),calcrweight(j,p));
              }
            }
          }
          if (j.absrap() > 1.0 && j.absrap() < 1.5) {
            _h_pt06_10_15->fill(jetpt/GeV);
            if (inRange(jetpt/GeV, 4., 6.)) {
              _numjets06_10_15_04_06->fill();
              _h_N06_10_15_04_06->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z06_10_15_04_06->fill(calcz(j,p));
                _h_ptrel06_10_15_04_06->fill(calcptrel(j,p));
                _h_rdA06_10_15_04_06->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 6., 10.)) {
              _numjets06_10_15_06_10->fill();
              _h_N06_10_15_06_10->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z06_10_15_06_10->fill(calcz(j,p));
                _h_ptrel06_10_15_06_10->fill(calcptrel(j,p));
                _h_rdA06_10_15_06_10->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 10., 15.)) {
              _numjets06_10_15_10_15->fill();
              _h_N06_10_15_10_15->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z06_10_15_10_15->fill(calcz(j,p));
                _h_ptrel06_10_15_10_15->fill(calcptrel(j,p));
                _h_rdA06_10_15_10_15->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 15., 24.)) {
              _numjets06_10_15_15_24->fill();
              _h_N06_10_15_15_24->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z06_10_15_15_24->fill(calcz(j,p));
                _h_ptrel06_10_15_15_24->fill(calcptrel(j,p));
                _h_rdA06_10_15_15_24->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 24., 40.)) {
              _numjets06_10_15_24_40->fill();
              _h_N06_10_15_24_40->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z06_10_15_24_40->fill(calcz(j,p));
                _h_ptrel06_10_15_24_40->fill(calcptrel(j,p));
                _h_rdA06_10_15_24_40->fill(calcr(j,p),calcrweight(j,p));
              }
            }
          }
          if (j.absrap() > 1.5 && j.absrap() < 1.9) {
            _h_pt06_15_19->fill(jetpt/GeV);
            if (inRange(jetpt/GeV, 4., 6.)) {
              _numjets06_15_19_04_06->fill();
              _h_N06_15_19_04_06->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z06_15_19_04_06->fill(calcz(j,p));
                _h_ptrel06_15_19_04_06->fill(calcptrel(j,p));
                _h_rdA06_15_19_04_06->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 6., 10.)) {
              _numjets06_15_19_06_10->fill();
              _h_N06_15_19_06_10->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z06_15_19_06_10->fill(calcz(j,p));
                _h_ptrel06_15_19_06_10->fill(calcptrel(j,p));
                _h_rdA06_15_19_06_10->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 10., 15.)) {
              _numjets06_15_19_10_15->fill();
              _h_N06_15_19_10_15->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z06_15_19_10_15->fill(calcz(j,p));
                _h_ptrel06_15_19_10_15->fill(calcptrel(j,p));
                _h_rdA06_15_19_10_15->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 15., 24.)) {
              _numjets06_15_19_15_24->fill();
              _h_N06_15_19_15_24->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z06_15_19_15_24->fill(calcz(j,p));
                _h_ptrel06_15_19_15_24->fill(calcptrel(j,p));
                _h_rdA06_15_19_15_24->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 24., 40.)) {
              _numjets06_15_19_24_40->fill();
              _h_N06_15_19_24_40->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z06_15_19_24_40->fill(calcz(j,p));
                _h_ptrel06_15_19_24_40->fill(calcptrel(j,p));
                _h_rdA06_15_19_24_40->fill(calcr(j,p),calcrweight(j,p));
              }
            }
          } // 1.5 < rapidity < 1.9
          if (j.absrap() < 1.9) {
            if (inRange(jetpt/GeV, 4., 6.)) {
              _numjets06_00_19_04_06->fill();
              _h_N06_00_19_04_06->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z06_00_19_04_06->fill(calcz(j,p));
                _h_ptrel06_00_19_04_06->fill(calcptrel(j,p));
                _h_rdA06_00_19_04_06->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 6., 10.)) {
              _numjets06_00_19_06_10->fill();
              _h_N06_00_19_06_10->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z06_00_19_06_10->fill(calcz(j,p));
                _h_ptrel06_00_19_06_10->fill(calcptrel(j,p));
                _h_rdA06_00_19_06_10->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 10., 15.)) {
              _numjets06_00_19_10_15->fill();
              _h_N06_00_19_10_15->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z06_00_19_10_15->fill(calcz(j,p));
                _h_ptrel06_00_19_10_15->fill(calcptrel(j,p));
                _h_rdA06_00_19_10_15->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 15., 24.)) {
              _numjets06_00_19_15_24->fill();
              _h_N06_00_19_15_24->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z06_00_19_15_24->fill(calcz(j,p));
                _h_ptrel06_00_19_15_24->fill(calcptrel(j,p));
                _h_rdA06_00_19_15_24->fill(calcr(j,p),calcrweight(j,p));
              }
            }
            if (inRange(jetpt/GeV, 24., 40.)) {
              _numjets06_00_19_24_40->fill();
              _h_N06_00_19_24_40->fill(j.particles().size());
              for (const Particle& p : j.particles()) {
                _h_z06_00_19_24_40->fill(calcz(j,p));
                _h_ptrel06_00_19_24_40->fill(calcptrel(j,p));
                _h_rdA06_00_19_24_40->fill(calcr(j,p),calcrweight(j,p));
              }
            }
          }
        } // each jet
      } // jets06 not empty

    } // end of event


    /// Normalise histograms etc., after the run
    void finalize() {

      // pT histos: d2sigma_jet / deta dpT
      const double xsec = crossSection()/microbarn;
      safeinvscale(_h_pt04_00_05, _sumofweights04->val()*(2*0.5)/xsec);
      safeinvscale(_h_pt06_00_05, _sumofweights06->val()*(2*0.5)/xsec);
      safeinvscale(_h_pt04_05_10, _sumofweights04->val()*(2*0.5)/xsec);
      safeinvscale(_h_pt06_05_10, _sumofweights06->val()*(2*0.5)/xsec);
      safeinvscale(_h_pt04_10_15, _sumofweights04->val()*(2*0.5)/xsec);
      safeinvscale(_h_pt06_10_15, _sumofweights06->val()*(2*0.5)/xsec);
      safeinvscale(_h_pt04_15_19, _sumofweights04->val()*(2*0.4)/xsec);
      safeinvscale(_h_pt06_15_19, _sumofweights06->val()*(2*0.4)/xsec);

      // N histos: 1/N_jet dN_jet / dN^{ch}_jet
      safeinvscale(_h_N04_00_05_04_06, _numjets04_00_05_04_06->val());
      safeinvscale(_h_N06_00_05_04_06, _numjets06_00_05_04_06->val());
      safeinvscale(_h_N04_00_05_06_10, _numjets04_00_05_06_10->val());
      safeinvscale(_h_N06_00_05_06_10, _numjets06_00_05_06_10->val());
      safeinvscale(_h_N04_00_05_10_15, _numjets04_00_05_10_15->val());
      safeinvscale(_h_N06_00_05_10_15, _numjets06_00_05_10_15->val());
      safeinvscale(_h_N04_00_05_15_24, _numjets04_00_05_15_24->val());
      safeinvscale(_h_N06_00_05_15_24, _numjets06_00_05_15_24->val());
      safeinvscale(_h_N04_00_05_24_40, _numjets04_00_05_24_40->val());
      safeinvscale(_h_N06_00_05_24_40, _numjets06_00_05_24_40->val());
      safeinvscale(_h_N04_05_10_04_06, _numjets04_05_10_04_06->val());
      safeinvscale(_h_N06_05_10_04_06, _numjets06_05_10_04_06->val());
      safeinvscale(_h_N04_05_10_06_10, _numjets04_05_10_06_10->val());
      safeinvscale(_h_N06_05_10_06_10, _numjets06_05_10_06_10->val());
      safeinvscale(_h_N04_05_10_10_15, _numjets04_05_10_10_15->val());
      safeinvscale(_h_N06_05_10_10_15, _numjets06_05_10_10_15->val());
      safeinvscale(_h_N04_05_10_15_24, _numjets04_05_10_15_24->val());
      safeinvscale(_h_N06_05_10_15_24, _numjets06_05_10_15_24->val());
      safeinvscale(_h_N04_05_10_24_40, _numjets04_05_10_24_40->val());
      safeinvscale(_h_N06_05_10_24_40, _numjets06_05_10_24_40->val());
      safeinvscale(_h_N04_10_15_04_06, _numjets04_10_15_04_06->val());
      safeinvscale(_h_N06_10_15_04_06, _numjets06_10_15_04_06->val());
      safeinvscale(_h_N04_10_15_06_10, _numjets04_10_15_06_10->val());
      safeinvscale(_h_N06_10_15_06_10, _numjets06_10_15_06_10->val());
      safeinvscale(_h_N04_10_15_10_15, _numjets04_10_15_10_15->val());
      safeinvscale(_h_N06_10_15_10_15, _numjets06_10_15_10_15->val());
      safeinvscale(_h_N04_10_15_15_24, _numjets04_10_15_15_24->val());
      safeinvscale(_h_N06_10_15_15_24, _numjets06_10_15_15_24->val());
      safeinvscale(_h_N04_10_15_24_40, _numjets04_10_15_24_40->val());
      safeinvscale(_h_N06_10_15_24_40, _numjets06_10_15_24_40->val());
      safeinvscale(_h_N04_15_19_04_06, _numjets04_15_19_04_06->val());
      safeinvscale(_h_N06_15_19_04_06, _numjets06_15_19_04_06->val());
      safeinvscale(_h_N04_15_19_06_10, _numjets04_15_19_06_10->val());
      safeinvscale(_h_N06_15_19_06_10, _numjets06_15_19_06_10->val());
      safeinvscale(_h_N04_15_19_10_15, _numjets04_15_19_10_15->val());
      safeinvscale(_h_N06_15_19_10_15, _numjets06_15_19_10_15->val());
      safeinvscale(_h_N04_15_19_15_24, _numjets04_15_19_15_24->val());
      safeinvscale(_h_N06_15_19_15_24, _numjets06_15_19_15_24->val());
      safeinvscale(_h_N04_15_19_24_40, _numjets04_15_19_24_40->val());
      safeinvscale(_h_N06_15_19_24_40, _numjets06_15_19_24_40->val());
      safeinvscale(_h_N04_00_19_04_06, _numjets04_00_19_04_06->val());
      safeinvscale(_h_N06_00_19_04_06, _numjets06_00_19_04_06->val());
      safeinvscale(_h_N04_00_19_06_10, _numjets04_00_19_06_10->val());
      safeinvscale(_h_N06_00_19_06_10, _numjets06_00_19_06_10->val());
      safeinvscale(_h_N04_00_19_10_15, _numjets04_00_19_10_15->val());
      safeinvscale(_h_N06_00_19_10_15, _numjets06_00_19_10_15->val());
      safeinvscale(_h_N04_00_19_15_24, _numjets04_00_19_15_24->val());
      safeinvscale(_h_N06_00_19_15_24, _numjets06_00_19_15_24->val());
      safeinvscale(_h_N04_00_19_24_40, _numjets04_00_19_24_40->val());
      safeinvscale(_h_N06_00_19_24_40, _numjets06_00_19_24_40->val());
      
      // z histos: 1/N_jet dN_track / dz_track->val()
      safeinvscale(_h_z04_00_05_04_06, _numjets04_00_05_04_06->val());
      safeinvscale(_h_z06_00_05_04_06, _numjets06_00_05_04_06->val());
      safeinvscale(_h_z04_00_05_06_10, _numjets04_00_05_06_10->val());
      safeinvscale(_h_z06_00_05_06_10, _numjets06_00_05_06_10->val());
      safeinvscale(_h_z04_00_05_10_15, _numjets04_00_05_10_15->val());
      safeinvscale(_h_z06_00_05_10_15, _numjets06_00_05_10_15->val());
      safeinvscale(_h_z04_00_05_15_24, _numjets04_00_05_15_24->val());
      safeinvscale(_h_z06_00_05_15_24, _numjets06_00_05_15_24->val());
      safeinvscale(_h_z04_00_05_24_40, _numjets04_00_05_24_40->val());
      safeinvscale(_h_z06_00_05_24_40, _numjets06_00_05_24_40->val());
      safeinvscale(_h_z04_05_10_04_06, _numjets04_05_10_04_06->val());
      safeinvscale(_h_z06_05_10_04_06, _numjets06_05_10_04_06->val());
      safeinvscale(_h_z04_05_10_06_10, _numjets04_05_10_06_10->val());
      safeinvscale(_h_z06_05_10_06_10, _numjets06_05_10_06_10->val());
      safeinvscale(_h_z04_05_10_10_15, _numjets04_05_10_10_15->val());
      safeinvscale(_h_z06_05_10_10_15, _numjets06_05_10_10_15->val());
      safeinvscale(_h_z04_05_10_15_24, _numjets04_05_10_15_24->val());
      safeinvscale(_h_z06_05_10_15_24, _numjets06_05_10_15_24->val());
      safeinvscale(_h_z04_05_10_24_40, _numjets04_05_10_24_40->val());
      safeinvscale(_h_z06_05_10_24_40, _numjets06_05_10_24_40->val());
      safeinvscale(_h_z04_10_15_04_06, _numjets04_10_15_04_06->val());
      safeinvscale(_h_z06_10_15_04_06, _numjets06_10_15_04_06->val());
      safeinvscale(_h_z04_10_15_06_10, _numjets04_10_15_06_10->val());
      safeinvscale(_h_z06_10_15_06_10, _numjets06_10_15_06_10->val());
      safeinvscale(_h_z04_10_15_10_15, _numjets04_10_15_10_15->val());
      safeinvscale(_h_z06_10_15_10_15, _numjets06_10_15_10_15->val());
      safeinvscale(_h_z04_10_15_15_24, _numjets04_10_15_15_24->val());
      safeinvscale(_h_z06_10_15_15_24, _numjets06_10_15_15_24->val());
      safeinvscale(_h_z04_10_15_24_40, _numjets04_10_15_24_40->val());
      safeinvscale(_h_z06_10_15_24_40, _numjets06_10_15_24_40->val());
      safeinvscale(_h_z04_15_19_04_06, _numjets04_15_19_04_06->val());
      safeinvscale(_h_z06_15_19_04_06, _numjets06_15_19_04_06->val());
      safeinvscale(_h_z04_15_19_06_10, _numjets04_15_19_06_10->val());
      safeinvscale(_h_z06_15_19_06_10, _numjets06_15_19_06_10->val());
      safeinvscale(_h_z04_15_19_10_15, _numjets04_15_19_10_15->val());
      safeinvscale(_h_z06_15_19_10_15, _numjets06_15_19_10_15->val());
      safeinvscale(_h_z04_15_19_15_24, _numjets04_15_19_15_24->val());
      safeinvscale(_h_z06_15_19_15_24, _numjets06_15_19_15_24->val());
      safeinvscale(_h_z04_15_19_24_40, _numjets04_15_19_24_40->val());
      safeinvscale(_h_z06_15_19_24_40, _numjets06_15_19_24_40->val());
      safeinvscale(_h_z04_00_19_04_06, _numjets04_00_19_04_06->val());
      safeinvscale(_h_z06_00_19_04_06, _numjets06_00_19_04_06->val());
      safeinvscale(_h_z04_00_19_06_10, _numjets04_00_19_06_10->val());
      safeinvscale(_h_z06_00_19_06_10, _numjets06_00_19_06_10->val());
      safeinvscale(_h_z04_00_19_10_15, _numjets04_00_19_10_15->val());
      safeinvscale(_h_z06_00_19_10_15, _numjets06_00_19_10_15->val());
      safeinvscale(_h_z04_00_19_15_24, _numjets04_00_19_15_24->val());
      safeinvscale(_h_z06_00_19_15_24, _numjets06_00_19_15_24->val());
      safeinvscale(_h_z04_00_19_24_40, _numjets04_00_19_24_40->val());
      safeinvscale(_h_z06_00_19_24_40, _numjets06_00_19_24_40->val());

      // ptrel histos: 1/N_jet dN_track / dptrel_track
      safeinvscale(_h_ptrel04_00_05_04_06, _numjets04_00_05_04_06->val());
      safeinvscale(_h_ptrel06_00_05_04_06, _numjets06_00_05_04_06->val());
      safeinvscale(_h_ptrel04_00_05_06_10, _numjets04_00_05_06_10->val());
      safeinvscale(_h_ptrel06_00_05_06_10, _numjets06_00_05_06_10->val());
      safeinvscale(_h_ptrel04_00_05_10_15, _numjets04_00_05_10_15->val());
      safeinvscale(_h_ptrel06_00_05_10_15, _numjets06_00_05_10_15->val());
      safeinvscale(_h_ptrel04_00_05_15_24, _numjets04_00_05_15_24->val());
      safeinvscale(_h_ptrel06_00_05_15_24, _numjets06_00_05_15_24->val());
      safeinvscale(_h_ptrel04_00_05_24_40, _numjets04_00_05_24_40->val());
      safeinvscale(_h_ptrel06_00_05_24_40, _numjets06_00_05_24_40->val());
      safeinvscale(_h_ptrel04_05_10_04_06, _numjets04_05_10_04_06->val());
      safeinvscale(_h_ptrel06_05_10_04_06, _numjets06_05_10_04_06->val());
      safeinvscale(_h_ptrel04_05_10_06_10, _numjets04_05_10_06_10->val());
      safeinvscale(_h_ptrel06_05_10_06_10, _numjets06_05_10_06_10->val());
      safeinvscale(_h_ptrel04_05_10_10_15, _numjets04_05_10_10_15->val());
      safeinvscale(_h_ptrel06_05_10_10_15, _numjets06_05_10_10_15->val());
      safeinvscale(_h_ptrel04_05_10_15_24, _numjets04_05_10_15_24->val());
      safeinvscale(_h_ptrel06_05_10_15_24, _numjets06_05_10_15_24->val());
      safeinvscale(_h_ptrel04_05_10_24_40, _numjets04_05_10_24_40->val());
      safeinvscale(_h_ptrel06_05_10_24_40, _numjets06_05_10_24_40->val());
      safeinvscale(_h_ptrel04_10_15_04_06, _numjets04_10_15_04_06->val());
      safeinvscale(_h_ptrel06_10_15_04_06, _numjets06_10_15_04_06->val());
      safeinvscale(_h_ptrel04_10_15_06_10, _numjets04_10_15_06_10->val());
      safeinvscale(_h_ptrel06_10_15_06_10, _numjets06_10_15_06_10->val());
      safeinvscale(_h_ptrel04_10_15_10_15, _numjets04_10_15_10_15->val());
      safeinvscale(_h_ptrel06_10_15_10_15, _numjets06_10_15_10_15->val());
      safeinvscale(_h_ptrel04_10_15_15_24, _numjets04_10_15_15_24->val());
      safeinvscale(_h_ptrel06_10_15_15_24, _numjets06_10_15_15_24->val());
      safeinvscale(_h_ptrel04_10_15_24_40, _numjets04_10_15_24_40->val());
      safeinvscale(_h_ptrel06_10_15_24_40, _numjets06_10_15_24_40->val());
      safeinvscale(_h_ptrel04_15_19_04_06, _numjets04_15_19_04_06->val());
      safeinvscale(_h_ptrel06_15_19_04_06, _numjets06_15_19_04_06->val());
      safeinvscale(_h_ptrel04_15_19_06_10, _numjets04_15_19_06_10->val());
      safeinvscale(_h_ptrel06_15_19_06_10, _numjets06_15_19_06_10->val());
      safeinvscale(_h_ptrel04_15_19_10_15, _numjets04_15_19_10_15->val());
      safeinvscale(_h_ptrel06_15_19_10_15, _numjets06_15_19_10_15->val());
      safeinvscale(_h_ptrel04_15_19_15_24, _numjets04_15_19_15_24->val());
      safeinvscale(_h_ptrel06_15_19_15_24, _numjets06_15_19_15_24->val());
      safeinvscale(_h_ptrel04_15_19_24_40, _numjets04_15_19_24_40->val());
      safeinvscale(_h_ptrel06_15_19_24_40, _numjets06_15_19_24_40->val());

      safeinvscale(_h_ptrel04_00_19_04_06, _numjets04_00_19_04_06->val());
      safeinvscale(_h_ptrel06_00_19_04_06, _numjets06_00_19_04_06->val());
      safeinvscale(_h_ptrel04_00_19_06_10, _numjets04_00_19_06_10->val());
      safeinvscale(_h_ptrel06_00_19_06_10, _numjets06_00_19_06_10->val());
      safeinvscale(_h_ptrel04_00_19_10_15, _numjets04_00_19_10_15->val());
      safeinvscale(_h_ptrel06_00_19_10_15, _numjets06_00_19_10_15->val());
      safeinvscale(_h_ptrel04_00_19_15_24, _numjets04_00_19_15_24->val());
      safeinvscale(_h_ptrel06_00_19_15_24, _numjets06_00_19_15_24->val());
      safeinvscale(_h_ptrel04_00_19_24_40, _numjets04_00_19_24_40->val());
      safeinvscale(_h_ptrel06_00_19_24_40, _numjets06_00_19_24_40->val());

      // r histos: 1/N_jet dN_track / dA
      safeinvscale(_h_rdA04_00_05_04_06, _numjets04_00_05_04_06->val());
      safeinvscale(_h_rdA06_00_05_04_06, _numjets06_00_05_04_06->val());
      safeinvscale(_h_rdA04_00_05_06_10, _numjets04_00_05_06_10->val());
      safeinvscale(_h_rdA06_00_05_06_10, _numjets06_00_05_06_10->val());
      safeinvscale(_h_rdA04_00_05_10_15, _numjets04_00_05_10_15->val());
      safeinvscale(_h_rdA06_00_05_10_15, _numjets06_00_05_10_15->val());
      safeinvscale(_h_rdA04_00_05_15_24, _numjets04_00_05_15_24->val());
      safeinvscale(_h_rdA06_00_05_15_24, _numjets06_00_05_15_24->val());
      safeinvscale(_h_rdA04_00_05_24_40, _numjets04_00_05_24_40->val());
      safeinvscale(_h_rdA06_00_05_24_40, _numjets06_00_05_24_40->val());
      safeinvscale(_h_rdA04_05_10_04_06, _numjets04_05_10_04_06->val());
      safeinvscale(_h_rdA06_05_10_04_06, _numjets06_05_10_04_06->val());
      safeinvscale(_h_rdA04_05_10_06_10, _numjets04_05_10_06_10->val());
      safeinvscale(_h_rdA06_05_10_06_10, _numjets06_05_10_06_10->val());
      safeinvscale(_h_rdA04_05_10_10_15, _numjets04_05_10_10_15->val());
      safeinvscale(_h_rdA06_05_10_10_15, _numjets06_05_10_10_15->val());
      safeinvscale(_h_rdA04_05_10_15_24, _numjets04_05_10_15_24->val());
      safeinvscale(_h_rdA06_05_10_15_24, _numjets06_05_10_15_24->val());
      safeinvscale(_h_rdA04_05_10_24_40, _numjets04_05_10_24_40->val());
      safeinvscale(_h_rdA06_05_10_24_40, _numjets06_05_10_24_40->val());
      safeinvscale(_h_rdA04_10_15_04_06, _numjets04_10_15_04_06->val());
      safeinvscale(_h_rdA06_10_15_04_06, _numjets06_10_15_04_06->val());
      safeinvscale(_h_rdA04_10_15_06_10, _numjets04_10_15_06_10->val());
      safeinvscale(_h_rdA06_10_15_06_10, _numjets06_10_15_06_10->val());
      safeinvscale(_h_rdA04_10_15_10_15, _numjets04_10_15_10_15->val());
      safeinvscale(_h_rdA06_10_15_10_15, _numjets06_10_15_10_15->val());
      safeinvscale(_h_rdA04_10_15_15_24, _numjets04_10_15_15_24->val());
      safeinvscale(_h_rdA06_10_15_15_24, _numjets06_10_15_15_24->val());
      safeinvscale(_h_rdA04_10_15_24_40, _numjets04_10_15_24_40->val());
      safeinvscale(_h_rdA06_10_15_24_40, _numjets06_10_15_24_40->val());
      safeinvscale(_h_rdA04_15_19_04_06, _numjets04_15_19_04_06->val());
      safeinvscale(_h_rdA06_15_19_04_06, _numjets06_15_19_04_06->val());
      safeinvscale(_h_rdA04_15_19_06_10, _numjets04_15_19_06_10->val());
      safeinvscale(_h_rdA06_15_19_06_10, _numjets06_15_19_06_10->val());
      safeinvscale(_h_rdA04_15_19_10_15, _numjets04_15_19_10_15->val());
      safeinvscale(_h_rdA06_15_19_10_15, _numjets06_15_19_10_15->val());
      safeinvscale(_h_rdA04_15_19_15_24, _numjets04_15_19_15_24->val());
      safeinvscale(_h_rdA06_15_19_15_24, _numjets06_15_19_15_24->val());
      safeinvscale(_h_rdA04_15_19_24_40, _numjets04_15_19_24_40->val());
      safeinvscale(_h_rdA06_15_19_24_40, _numjets06_15_19_24_40->val());

      safeinvscale(_h_rdA04_00_19_04_06, _numjets04_00_19_04_06->val());
      safeinvscale(_h_rdA06_00_19_04_06, _numjets06_00_19_04_06->val());
      safeinvscale(_h_rdA04_00_19_06_10, _numjets04_00_19_06_10->val());
      safeinvscale(_h_rdA06_00_19_06_10, _numjets06_00_19_06_10->val());
      safeinvscale(_h_rdA04_00_19_10_15, _numjets04_00_19_10_15->val());
      safeinvscale(_h_rdA06_00_19_10_15, _numjets06_00_19_10_15->val());
      safeinvscale(_h_rdA04_00_19_15_24, _numjets04_00_19_15_24->val());
      safeinvscale(_h_rdA06_00_19_15_24, _numjets06_00_19_15_24->val());
      safeinvscale(_h_rdA04_00_19_24_40, _numjets04_00_19_24_40->val());
      safeinvscale(_h_rdA06_00_19_24_40, _numjets06_00_19_24_40->val());
    }

    //@}


  private:

    void safeinvscale(Histo1DPtr h, double denom) {
      if (denom != 0) {
        scale(h, 1.0/denom);
      } else {
        normalize(h, 0);
      }
    }


    /// Event weights
    CounterPtr _sumofweights04, _sumofweights06;


    /// Jet counters
    CounterPtr _numjets04_00_05_04_06, _numjets04_00_05_06_10, _numjets04_00_05_10_15, _numjets04_00_05_15_24, _numjets04_00_05_24_40;
    CounterPtr _numjets06_00_05_04_06, _numjets06_00_05_06_10, _numjets06_00_05_10_15, _numjets06_00_05_15_24, _numjets06_00_05_24_40;
    CounterPtr _numjets04_05_10_04_06, _numjets04_05_10_06_10, _numjets04_05_10_10_15, _numjets04_05_10_15_24, _numjets04_05_10_24_40;
    CounterPtr _numjets06_05_10_04_06, _numjets06_05_10_06_10, _numjets06_05_10_10_15, _numjets06_05_10_15_24, _numjets06_05_10_24_40;
    CounterPtr _numjets04_10_15_04_06, _numjets04_10_15_06_10, _numjets04_10_15_10_15, _numjets04_10_15_15_24, _numjets04_10_15_24_40;
    CounterPtr _numjets06_10_15_04_06, _numjets06_10_15_06_10, _numjets06_10_15_10_15, _numjets06_10_15_15_24, _numjets06_10_15_24_40;
    CounterPtr _numjets04_15_19_04_06, _numjets04_15_19_06_10, _numjets04_15_19_10_15, _numjets04_15_19_15_24, _numjets04_15_19_24_40;
    CounterPtr _numjets06_15_19_04_06, _numjets06_15_19_06_10, _numjets06_15_19_10_15, _numjets06_15_19_15_24, _numjets06_15_19_24_40;
    CounterPtr _numjets04_00_19_04_06, _numjets04_00_19_06_10, _numjets04_00_19_10_15, _numjets04_00_19_15_24, _numjets04_00_19_24_40;
    CounterPtr _numjets06_00_19_04_06, _numjets06_00_19_06_10, _numjets06_00_19_10_15, _numjets06_00_19_15_24, _numjets06_00_19_24_40;


  private:

    /// @name Histograms
    //@{

    Histo1DPtr _h_pt04_00_05, _h_pt06_00_05;
    Histo1DPtr _h_N04_00_05_04_06, _h_N06_00_05_04_06;
    Histo1DPtr _h_N04_00_05_06_10, _h_N06_00_05_06_10;
    Histo1DPtr _h_N04_00_05_10_15, _h_N06_00_05_10_15;
    Histo1DPtr _h_N04_00_05_15_24, _h_N06_00_05_15_24;
    Histo1DPtr _h_N04_00_05_24_40, _h_N06_00_05_24_40;
    Histo1DPtr _h_z04_00_05_04_06, _h_z06_00_05_04_06;
    Histo1DPtr _h_z04_00_05_06_10, _h_z06_00_05_06_10;
    Histo1DPtr _h_z04_00_05_10_15, _h_z06_00_05_10_15;
    Histo1DPtr _h_z04_00_05_15_24, _h_z06_00_05_15_24;
    Histo1DPtr _h_z04_00_05_24_40, _h_z06_00_05_24_40;
    Histo1DPtr _h_ptrel04_00_05_04_06, _h_ptrel06_00_05_04_06;
    Histo1DPtr _h_ptrel04_00_05_06_10, _h_ptrel06_00_05_06_10;
    Histo1DPtr _h_ptrel04_00_05_10_15, _h_ptrel06_00_05_10_15;
    Histo1DPtr _h_ptrel04_00_05_15_24, _h_ptrel06_00_05_15_24;
    Histo1DPtr _h_ptrel04_00_05_24_40, _h_ptrel06_00_05_24_40;
    Histo1DPtr _h_rdA04_00_05_04_06, _h_rdA06_00_05_04_06;
    Histo1DPtr _h_rdA04_00_05_06_10, _h_rdA06_00_05_06_10;
    Histo1DPtr _h_rdA04_00_05_10_15, _h_rdA06_00_05_10_15;
    Histo1DPtr _h_rdA04_00_05_15_24, _h_rdA06_00_05_15_24;
    Histo1DPtr _h_rdA04_00_05_24_40, _h_rdA06_00_05_24_40;

    Histo1DPtr _h_pt04_05_10, _h_pt06_05_10;
    Histo1DPtr _h_N04_05_10_04_06, _h_N06_05_10_04_06;
    Histo1DPtr _h_N04_05_10_06_10, _h_N06_05_10_06_10;
    Histo1DPtr _h_N04_05_10_10_15, _h_N06_05_10_10_15;
    Histo1DPtr _h_N04_05_10_15_24, _h_N06_05_10_15_24;
    Histo1DPtr _h_N04_05_10_24_40, _h_N06_05_10_24_40;
    Histo1DPtr _h_z04_05_10_04_06, _h_z06_05_10_04_06;
    Histo1DPtr _h_z04_05_10_06_10, _h_z06_05_10_06_10;
    Histo1DPtr _h_z04_05_10_10_15, _h_z06_05_10_10_15;
    Histo1DPtr _h_z04_05_10_15_24, _h_z06_05_10_15_24;
    Histo1DPtr _h_z04_05_10_24_40, _h_z06_05_10_24_40;
    Histo1DPtr _h_ptrel04_05_10_04_06, _h_ptrel06_05_10_04_06;
    Histo1DPtr _h_ptrel04_05_10_06_10, _h_ptrel06_05_10_06_10;
    Histo1DPtr _h_ptrel04_05_10_10_15, _h_ptrel06_05_10_10_15;
    Histo1DPtr _h_ptrel04_05_10_15_24, _h_ptrel06_05_10_15_24;
    Histo1DPtr _h_ptrel04_05_10_24_40, _h_ptrel06_05_10_24_40;
    Histo1DPtr _h_rdA04_05_10_04_06, _h_rdA06_05_10_04_06;
    Histo1DPtr _h_rdA04_05_10_06_10, _h_rdA06_05_10_06_10;
    Histo1DPtr _h_rdA04_05_10_10_15, _h_rdA06_05_10_10_15;
    Histo1DPtr _h_rdA04_05_10_15_24, _h_rdA06_05_10_15_24;
    Histo1DPtr _h_rdA04_05_10_24_40, _h_rdA06_05_10_24_40;

    Histo1DPtr _h_pt04_10_15, _h_pt06_10_15;
    Histo1DPtr _h_N04_10_15_04_06, _h_N06_10_15_04_06;
    Histo1DPtr _h_N04_10_15_06_10, _h_N06_10_15_06_10;
    Histo1DPtr _h_N04_10_15_10_15, _h_N06_10_15_10_15;
    Histo1DPtr _h_N04_10_15_15_24, _h_N06_10_15_15_24;
    Histo1DPtr _h_N04_10_15_24_40, _h_N06_10_15_24_40;
    Histo1DPtr _h_z04_10_15_04_06, _h_z06_10_15_04_06;
    Histo1DPtr _h_z04_10_15_06_10, _h_z06_10_15_06_10;
    Histo1DPtr _h_z04_10_15_10_15, _h_z06_10_15_10_15;
    Histo1DPtr _h_z04_10_15_15_24, _h_z06_10_15_15_24;
    Histo1DPtr _h_z04_10_15_24_40, _h_z06_10_15_24_40;
    Histo1DPtr _h_ptrel04_10_15_04_06, _h_ptrel06_10_15_04_06;
    Histo1DPtr _h_ptrel04_10_15_06_10, _h_ptrel06_10_15_06_10;
    Histo1DPtr _h_ptrel04_10_15_10_15, _h_ptrel06_10_15_10_15;
    Histo1DPtr _h_ptrel04_10_15_15_24, _h_ptrel06_10_15_15_24;
    Histo1DPtr _h_ptrel04_10_15_24_40, _h_ptrel06_10_15_24_40;
    Histo1DPtr _h_rdA04_10_15_04_06, _h_rdA06_10_15_04_06;
    Histo1DPtr _h_rdA04_10_15_06_10, _h_rdA06_10_15_06_10;
    Histo1DPtr _h_rdA04_10_15_10_15, _h_rdA06_10_15_10_15;
    Histo1DPtr _h_rdA04_10_15_15_24, _h_rdA06_10_15_15_24;
    Histo1DPtr _h_rdA04_10_15_24_40, _h_rdA06_10_15_24_40;

    Histo1DPtr _h_pt04_15_19, _h_pt06_15_19;
    Histo1DPtr _h_N04_15_19_04_06, _h_N06_15_19_04_06;
    Histo1DPtr _h_N04_15_19_06_10, _h_N06_15_19_06_10;
    Histo1DPtr _h_N04_15_19_10_15, _h_N06_15_19_10_15;
    Histo1DPtr _h_N04_15_19_15_24, _h_N06_15_19_15_24;
    Histo1DPtr _h_N04_15_19_24_40, _h_N06_15_19_24_40;
    Histo1DPtr _h_z04_15_19_04_06, _h_z06_15_19_04_06;
    Histo1DPtr _h_z04_15_19_06_10, _h_z06_15_19_06_10;
    Histo1DPtr _h_z04_15_19_10_15, _h_z06_15_19_10_15;
    Histo1DPtr _h_z04_15_19_15_24, _h_z06_15_19_15_24;
    Histo1DPtr _h_z04_15_19_24_40, _h_z06_15_19_24_40;
    Histo1DPtr _h_ptrel04_15_19_04_06, _h_ptrel06_15_19_04_06;
    Histo1DPtr _h_ptrel04_15_19_06_10, _h_ptrel06_15_19_06_10;
    Histo1DPtr _h_ptrel04_15_19_10_15, _h_ptrel06_15_19_10_15;
    Histo1DPtr _h_ptrel04_15_19_15_24, _h_ptrel06_15_19_15_24;
    Histo1DPtr _h_ptrel04_15_19_24_40, _h_ptrel06_15_19_24_40;
    Histo1DPtr _h_rdA04_15_19_04_06, _h_rdA06_15_19_04_06;
    Histo1DPtr _h_rdA04_15_19_06_10, _h_rdA06_15_19_06_10;
    Histo1DPtr _h_rdA04_15_19_10_15, _h_rdA06_15_19_10_15;
    Histo1DPtr _h_rdA04_15_19_15_24, _h_rdA06_15_19_15_24;
    Histo1DPtr _h_rdA04_15_19_24_40, _h_rdA06_15_19_24_40;

    Histo1DPtr _h_N04_00_19_04_06, _h_N06_00_19_04_06;
    Histo1DPtr _h_N04_00_19_06_10, _h_N06_00_19_06_10;
    Histo1DPtr _h_N04_00_19_10_15, _h_N06_00_19_10_15;
    Histo1DPtr _h_N04_00_19_15_24, _h_N06_00_19_15_24;
    Histo1DPtr _h_N04_00_19_24_40, _h_N06_00_19_24_40;
    Histo1DPtr _h_z04_00_19_04_06, _h_z06_00_19_04_06;
    Histo1DPtr _h_z04_00_19_06_10, _h_z06_00_19_06_10;
    Histo1DPtr _h_z04_00_19_10_15, _h_z06_00_19_10_15;
    Histo1DPtr _h_z04_00_19_15_24, _h_z06_00_19_15_24;
    Histo1DPtr _h_z04_00_19_24_40, _h_z06_00_19_24_40;
    Histo1DPtr _h_ptrel04_00_19_04_06, _h_ptrel06_00_19_04_06;
    Histo1DPtr _h_ptrel04_00_19_06_10, _h_ptrel06_00_19_06_10;
    Histo1DPtr _h_ptrel04_00_19_10_15, _h_ptrel06_00_19_10_15;
    Histo1DPtr _h_ptrel04_00_19_15_24, _h_ptrel06_00_19_15_24;
    Histo1DPtr _h_ptrel04_00_19_24_40, _h_ptrel06_00_19_24_40;
    Histo1DPtr _h_rdA04_00_19_04_06, _h_rdA06_00_19_04_06;
    Histo1DPtr _h_rdA04_00_19_06_10, _h_rdA06_00_19_06_10;
    Histo1DPtr _h_rdA04_00_19_10_15, _h_rdA06_00_19_10_15;
    Histo1DPtr _h_rdA04_00_19_15_24, _h_rdA06_00_19_15_24;
    Histo1DPtr _h_rdA04_00_19_24_40, _h_rdA06_00_19_24_40;

    //@}

  };


  // This global object acts as a hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2011_I919017);

}
