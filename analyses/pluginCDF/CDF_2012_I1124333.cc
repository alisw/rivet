// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {

  /// @ CDF Run II Z \f$ p_\perp \f$ in Drell-Yan events
  /// @author Simone Amoroso
  class CDF_2012_I1124333 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CDF_2012_I1124333);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      ///  Initialise and register projections here
      ZFinder zfinder(FinalState(), Cuts::open(), PID::ELECTRON, 66*GeV, 116*GeV, 0.0, ZFinder::ClusterPhotons::NONE);
      declare(zfinder, "ZFinder");


      ///  Book histograms here, e.g.:
      //book(      _hist_z_xs ,1, 1, 1);
      book(_hist_zpt, 2, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      /// @todo Do the event by event analysis here
      const ZFinder& zfinder = apply<ZFinder>(event, "ZFinder");
      if (zfinder.bosons().size() != 1) {
        MSG_DEBUG("Num e+ e- pairs found = " << zfinder.bosons().size());
        vetoEvent;
      }
      const FourMomentum& pZ = zfinder.bosons()[0].momentum();
      if (pZ.mass2() < 0) {
        MSG_DEBUG("Negative Z mass**2 = " << pZ.mass2()/GeV2 << "!");
        vetoEvent;
      }

      MSG_DEBUG("Dilepton mass = " << pZ.mass()/GeV << " GeV");
      _hist_zpt->fill(pZ.pT());
      //      _hist_z_xs->fill(1);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_hist_zpt, crossSection()/picobarn/sumOfWeights());
    }

    //@}


  private:
    /// @name Histograms
    //@{
    Histo1DPtr _hist_zpt;
    //    Histo1DPtr _hist_z_xs;
    //@}


  };



  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(CDF_2012_I1124333);

}
