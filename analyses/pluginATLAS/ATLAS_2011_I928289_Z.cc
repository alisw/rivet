// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  class ATLAS_2011_I928289_Z : public Analysis {
  public:

    /// Constructor
    ATLAS_2011_I928289_Z()
      : Analysis("ATLAS_2011_I928289_Z")
    {

    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs;

      Cut cut = (Cuts::pT >= 20.0*GeV);

      ZFinder zfinder_ee_bare(   fs, cut, PID::ELECTRON, 66.0*GeV, 116.0*GeV, 0.0, ZFinder::ClusterPhotons::NODECAY, ZFinder::AddPhotons::NO);
      ZFinder zfinder_ee_dressed(fs, cut, PID::ELECTRON, 66.0*GeV, 116.0*GeV, 0.1, ZFinder::ClusterPhotons::NODECAY, ZFinder::AddPhotons::NO);
      ZFinder zfinder_mm_bare(   fs, cut, PID::MUON    , 66.0*GeV, 116.0*GeV, 0.0, ZFinder::ClusterPhotons::NODECAY, ZFinder::AddPhotons::NO);
      ZFinder zfinder_mm_dressed(fs, cut, PID::MUON    , 66.0*GeV, 116.0*GeV, 0.1, ZFinder::ClusterPhotons::NODECAY, ZFinder::AddPhotons::NO);

      declare(zfinder_ee_bare   , "ZFinder_ee_bare"   );
      declare(zfinder_ee_dressed, "ZFinder_ee_dressed");
      declare(zfinder_mm_bare   , "ZFinder_mm_bare"   );
      declare(zfinder_mm_dressed, "ZFinder_mm_dressed");

      // y(Z) cross-section dependence
      book(_h_Z_y_ee_bare     ,1, 1, 1);
      book(_h_Z_y_ee_dressed  ,1, 1, 2);
      book(_h_Z_y_mm_bare     ,1, 1, 3);
      book(_h_Z_y_mm_dressed  ,1, 1, 4);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ZFinder& zfinder_ee_bare     = apply<ZFinder>(event, "ZFinder_ee_bare"   );
      const ZFinder& zfinder_ee_dressed  = apply<ZFinder>(event, "ZFinder_ee_dressed");
      const ZFinder& zfinder_mm_bare     = apply<ZFinder>(event, "ZFinder_mm_bare"   );
      const ZFinder& zfinder_mm_dressed  = apply<ZFinder>(event, "ZFinder_mm_dressed");

      fillPlots1D(zfinder_ee_bare   , _h_Z_y_ee_bare);
      fillPlots1D(zfinder_ee_dressed, _h_Z_y_ee_dressed);
      fillPlots1D(zfinder_mm_bare   , _h_Z_y_mm_bare);
      fillPlots1D(zfinder_mm_dressed, _h_Z_y_mm_dressed);

    }


    void fillPlots1D(const ZFinder& zfinder, Histo1DPtr hist) {
      if (zfinder.bosons().size() != 1) return;
      const FourMomentum zmom = zfinder.bosons()[0].momentum();
      hist->fill(zmom.absrap());
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // Print summary info
      const double xs_pb(crossSection() / picobarn);
      const double sumw(sumOfWeights());

      // Normalise, scale and otherwise manipulate histograms here
      const double sf(0.5 * xs_pb / sumw); // 0.5 accounts for rapidity bin width
      scale(_h_Z_y_ee_bare   , sf);
      scale(_h_Z_y_ee_dressed, sf);
      scale(_h_Z_y_mm_bare   , sf);
      scale(_h_Z_y_mm_dressed, sf);

    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_Z_y_ee_bare;
    Histo1DPtr _h_Z_y_ee_dressed;
    Histo1DPtr _h_Z_y_mm_bare;
    Histo1DPtr _h_Z_y_mm_dressed;
    //@}

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2011_I928289_Z);

}
