// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include <fstream>

#define _USE_MATH_DEFINES
#include <cmath>

namespace Rivet {

  /// @brief ALICE PbPb at 2.76 TeV multiplicity at mid-rapidity
  class ALICE_2010_I880049 : public Analysis {

  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2010_I880049);

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Declare centrality projection
      declareCentrality(ALICE::V0MMultiplicity(),
        "ALICE_2015_PBPBCentrality", "V0M", "V0M");

      // Trigger projections
      declare(ChargedFinalState((Cuts::eta > 2.8 && Cuts::eta < 5.1) &&
        Cuts::pT > 0.1*GeV), "VZERO1");
      declare(ChargedFinalState((Cuts::eta > -3.7 && Cuts::eta < -1.7) &&
        Cuts::pT > 0.1*GeV), "VZERO2");
      declare(ChargedFinalState(Cuts::abseta < 1. && Cuts::pT > 0.15*GeV),
        "SPD");

      // Charged, primary particles with |eta| < 0.5 and pT > 50 MeV
      declare(ALICE::PrimaryParticles(Cuts::abseta < 0.5 &&
        Cuts::pT > 50*MeV && Cuts::abscharge > 0), "APRIM");

      // Histograms and variables initialization
      _histNchVsCentr = bookProfile1D(1, 1, 1);
      _histNpartVsCentr = bookProfile1D(1, 1, 2);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      // Charged, primary particles with at least pT = 50 MeV
      // in eta range of |eta| < 0.5
      Particles chargedParticles =
        applyProjection<ALICE::PrimaryParticles>(event,"APRIM").particles();

      // Trigger projections
      const ChargedFinalState& vz1 =
        applyProjection<ChargedFinalState>(event,"VZERO1");
      const ChargedFinalState& vz2 =
        applyProjection<ChargedFinalState>(event,"VZERO2");
      const ChargedFinalState& spd =
        applyProjection<ChargedFinalState>(event,"SPD");
      int fwdTrig = (vz1.particles().size() > 0 ? 1 : 0);
      int bwdTrig = (vz2.particles().size() > 0 ? 1 : 0);
      int cTrig = (spd.particles().size() > 0 ? 1 : 0);

      if (fwdTrig + bwdTrig + cTrig < 2) vetoEvent;

      const CentralityProjection& centrProj =
        apply<CentralityProjection>(event, "V0M");
      double centr = centrProj();
      if (centr > 80.)
        vetoEvent;

      // Calculate number of charged particles and fill histogram
      double nch = chargedParticles.size();
      _histNchVsCentr->fill(centr, nch, weight);

      // Attempt to extract Npart form GenEvent.
      // TODO: Unclear how to handle this in HepMC3
      const HepMC::HeavyIon* hi = event.genEvent()->heavy_ion();
      if (hi && hi->is_valid()) {
        _histNpartVsCentr->fill(centr, hi->Npart_proj() + hi->Npart_targ(),
          weight);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

    }

    //@}

  private:

    /// @name Histograms
    //@{
    Profile1DPtr _histNchVsCentr;
    Profile1DPtr _histNpartVsCentr;
    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2010_I880049);

}
