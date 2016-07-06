// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief SFM charged multiplicities in NSD and inelastic minbias events
  class SFM_1984_S1178091 : public Analysis {
  public:

    /// Constructor
    SFM_1984_S1178091() : Analysis("SFM_1984_S1178091") {
      _sumW = 0;
      _sumWDiff = 0;
    }


    /// @name Analysis methods
    //@{

    void init() {
      // Projections
      /// @todo Corrected to full phase space?
      addProjection(ChargedFinalState(-10,10,250*MeV), "FS");

      // Histograms
      if (fuzzyEquals(sqrtS()/GeV, 30.4, 1E-1)) {
        _hist_multiplicity_inel = bookHisto1D(1, 1, 1);
        _hist_multiplicity_nsd = bookHisto1D(2, 1, 1);
      } else if (fuzzyEquals(sqrtS(), 44.5, 1E-1)) {
        _hist_multiplicity_inel = bookHisto1D(1, 1, 2);
        _hist_multiplicity_nsd = bookHisto1D(2, 1, 2);
      } else if (fuzzyEquals(sqrtS(), 52.2, 1E-1)) {
        _hist_multiplicity_inel = bookHisto1D(1, 1, 3);
        _hist_multiplicity_nsd = bookHisto1D(2, 1, 3);
      } else if (fuzzyEquals(sqrtS(), 62.2, 1E-1)) {
        _hist_multiplicity_inel = bookHisto1D(1, 1, 4);
        _hist_multiplicity_nsd = bookHisto1D(2, 1, 4);
      }

    }


    // Analyse each event
    void analyze(const Event& event) {
      const double weight = event.weight();
      const ChargedFinalState& fs = applyProjection<ChargedFinalState>(event, "FS");
      //const size_t numParticles = fs.particles().size();
      size_t N=0;

      /// @todo Any trigger?
      //
      // Trigger
      //if (numParticles <1 ) vetoEvent;

      // Decide whether event is of diffractive type or not
      // @todo It is not so clear in the paper how this distinction is made.
      // They seem to require either exactly one particle with Feynman x larger
      // than 0.8 to call an event diffractive or that there are no tracks
      // reconstructed in either of the two hemispheres. For the latter
      // they require in addition also the number of charged particles
      // to be smaller than 8.
      int n_left(0), n_right(0), n_large_x(0);
      foreach (const Particle& p, fs.particles()) {
        // Calculate the particles' Feynman x
        if (p.pT() <=3*GeV) {
          N++;
          const double x_feyn = 2.0 * fabs(p.pz())/sqrtS();
          if (x_feyn > 0.8 ) n_large_x += 1;

          // Pseudorapidity
          const double eta = p.eta();
          if (eta > 0.0) n_right += 1;
          else if (eta < 0.0) n_left += 1;
        }
      }
      MSG_DEBUG("N_left: " << n_left << ", "
                << "N_right: " << n_right << ", "
                << "N_large_x: " << n_large_x);

      if ( N < 1 ) vetoEvent;
      // Not sure about the "=="!
      // @todo Not sure about the "== 1", the paper says no charged particle
      // that was reconstructed so the incoming protons must run down the beam
      // pipe. Since we look a the complete final state here no particle being
      // reconstructed should be equal to one particle (proton) in each
      // hemisphere.  The "< 8" is also not certain.
      const bool isDiffractive = (n_large_x == 1) ||
        (( (n_left + n_right) == 1) && (N < 7) );
        //((n_left == 1 || n_right == 1) && N < 8 );
        //((n_left == 1 || n_right == 1) && numParticles < 8 );

      // Increment weight counters
      _sumW += weight;
      _sumWDiff += weight;

      // Fill histos of charged multiplicity distributions
      //_hist_multiplicity_inel->fill(numParticles, weight);
      //if (!isDiffractive) _hist_multiplicity_nsd->fill(numParticles, weight);
      _hist_multiplicity_inel->fill(N, weight);
      if (!isDiffractive) _hist_multiplicity_nsd->fill(N, weight);
    }


    void finalize() {
      scale(_hist_multiplicity_inel, 1.0/_sumWDiff);
      scale(_hist_multiplicity_nsd, 1.0/_sumW );
    }

    //@}


  private:

    /// @name Weight counters
    //@{
    double _sumW, _sumWDiff;
    //@}

    /// @name Histograms
    //@{
    Histo1DPtr _hist_multiplicity_inel;
    Histo1DPtr _hist_multiplicity_nsd;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(SFM_1984_S1178091);

}
