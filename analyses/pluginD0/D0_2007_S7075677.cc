// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// @brief Measurement of D0 Run II Z \f$ p_\perp \f$ diff cross-section shape
  /// @author Andy Buckley
  /// @author Gavin Hesketh
  /// @author Frank Siegert
  class D0_2007_S7075677 : public Analysis {
  public:

    /// Default constructor.
    D0_2007_S7075677()
      : Analysis("D0_2007_S7075677")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      ZFinder zfinder(FinalState(), Cuts::open(), PID::ELECTRON,
                      71*GeV, 111*GeV, 0.2, ZFinder::ClusterPhotons::NODECAY, ZFinder::AddPhotons::YES);
      declare(zfinder, "ZFinder");

      book(_h_yZ ,1, 1, 1);
    }


    /// Do the analysis
    void analyze(const Event & e) {
      const ZFinder& zfinder = apply<ZFinder>(e, "ZFinder");
      if (zfinder.bosons().size() == 1) {
        const Particles& el(zfinder.constituents());
        if (el[0].pT() > 25*GeV || el[1].pT() > 25*GeV) {
          _h_yZ->fill(fabs(zfinder.bosons()[0].rapidity()));
        }
      } else {
        MSG_DEBUG("No unique lepton pair found.");
      }
    }


    // Finalize
    void finalize() {
      // Data seems to have been normalized for the avg of the two sides
      // (+ve & -ve rapidity) rather than the sum, hence the 0.5:
      normalize(_h_yZ, 0.5);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_yZ;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(D0_2007_S7075677);

}
