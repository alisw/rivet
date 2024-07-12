// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {

  class ALICE_2012_I1181770 : public Analysis {
  public:

    ALICE_2012_I1181770()
      : Analysis("ALICE_2012_I1181770")
    {    }


    void init() {
      // Projection setup
      declare(ChargedFinalState(), "CFS");

      // Book (energy-specific) histograms
      int isqrts = -1;
      if (isCompatibleWithSqrtS(900)) isqrts = 1;
      else if (isCompatibleWithSqrtS(2760)) isqrts = 2;
      else if (isCompatibleWithSqrtS(7000)) isqrts = 3;
      assert(isqrts > 0);

      book(_h_frac_sd_inel, 1, 1, isqrts);
      book(_h_frac_dd_inel, 2, 1, isqrts);
      book(_h_xsec_sd     , 3, 1, isqrts);
      book(_h_xsec_dd     , 4, 1, isqrts);
      book(_h_xsec_inel   , 5, 1, isqrts);
    }


    void analyze(const Event& event) {
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      if (cfs.size() < 2) vetoEvent; // need at least two particles to calculate gaps

      // Fill INEL plots for each event
      _h_xsec_inel->fill(sqrtS()/GeV);

      // Identify particles with most positive/most negative rapidities
      const Particles particlesByRap = cfs.particles(cmpMomByRap);
      const Particle pslowest = particlesByRap.front();
      const Particle pfastest = particlesByRap.back();

      // Find gap sizes
      const Particles particlesByEta = cfs.particles(cmpMomByEta); // sorted from minus to plus
      const size_t num_particles = particlesByEta.size();
      vector<double> gaps;
      for (size_t ip = 1; ip < num_particles; ++ip) {
        const Particle& p1 = particlesByEta[ip-1];
        const Particle& p2 = particlesByEta[ip];
        const double gap = p2.eta() - p1.eta();
        assert(gap >= 0);
        gaps.push_back(gap);
      }

      // First, last, and largest gaps
      const double gapmax = *max_element(gaps.begin(), gaps.end());
      const double gapbwd = gaps.front();
      const double gapfwd = gaps.back();

      // Mx calculation
      FourMomentum p4lead;
      if (pslowest.pid() == PID::PROTON && pfastest.pid() == PID::PROTON) {
        p4lead = (fabs(pslowest.rapidity()) > fabs(pfastest.rapidity())) ? pslowest.momentum() : pfastest.momentum();
      } else if (pslowest.pid() == PID::PROTON) {
        p4lead = pslowest.momentum();
      } else if (pfastest.pid() == PID::PROTON) {
        p4lead = pfastest.momentum();
      }
      const double Mx = sqrt( (sqrtS()-p4lead.E()-p4lead.p3().mod()) * (sqrtS()-p4lead.E()+p4lead.p3().mod()) );

      // Fill SD (and escape) if Mx is sufficiently low
      if (Mx < 200*GeV) {
        _h_xsec_sd->fill(sqrtS()/GeV);
        return;
      }

      // Also remove SD-like events in NSD events
      if (fuzzyEquals(gapbwd, gapmax) || fuzzyEquals(gapfwd, gapmax)) vetoEvent;

      // Fill DD plots
      if (gapmax > 3) _h_xsec_dd->fill(sqrtS()/GeV);
    }


    void finalize() {

      // get the ratio plots: SD/inel, DD/inel
      divide(_h_xsec_sd , _h_xsec_inel, _h_frac_sd_inel);
      divide(_h_xsec_sd , _h_xsec_inel, _h_frac_dd_inel);

      const double scaling = crossSection()/millibarn/sumOfWeights();
      scale(_h_xsec_sd,   scaling);
      scale(_h_xsec_dd,   scaling);
      scale(_h_xsec_inel, scaling);

    }

  private:

    Scatter2DPtr _h_frac_sd_inel;
    Scatter2DPtr _h_frac_dd_inel;
    Histo1DPtr   _h_xsec_sd;
    Histo1DPtr   _h_xsec_dd;
    Histo1DPtr   _h_xsec_inel;

  };

  // Hook for the plugin system
  RIVET_DECLARE_PLUGIN(ALICE_2012_I1181770);

}
