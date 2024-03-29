// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include <fstream>

#define _USE_MATH_DEFINES
#include <math.h>

namespace Rivet {

  /// @brief ALICE PbPb at 2.76 TeV R_AA analysis.
  class ALICE_2012_I1127497 : public Analysis {

  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2012_I1127497);

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Declare centrality projection
      declareCentrality(ALICE::V0MMultiplicity(),
        "ALICE_2015_PBPBCentrality", "V0M", "V0M");

      // Charged, primary particles with |eta| < 0.5 and pT > 150 MeV
      declare(ALICE::PrimaryParticles(Cuts::abseta < 0.5 &&
        Cuts::pT > 150*MeV && Cuts::abscharge > 0), "APRIM");

      // Loop over all histograms
      for (size_t ihist = 0; ihist < NHISTOS; ++ihist) {

        // Initialize PbPb objects
        _histNch[PBPB][ihist] = bookHisto1D(ihist+1, 1, 1);

        std::string nameCounterPbPb = "counter.pbpb." + std::to_string(ihist);
        _counterSOW[PBPB][ihist] = bookCounter(nameCounterPbPb,
          "Sum of weights counter for PbPb");

        std::string nameCounterNcoll = "counter.ncoll." + std::to_string(ihist);
        _counterNcoll[ihist] = bookCounter(nameCounterNcoll,
          "Ncoll counter for PbPb");

        // Initialize pp objects. In principle, only one pp histogram would be
        // needed since centrality does not make any difference here. However,
        // in some cases in this analysis the binning differ from each other,
        // so this is easy-to-implement way to account for that.
        std::string namePP = _histNch[PBPB][ihist]->name() + "-pp";
        // The binning is taken from the reference data
        _histNch[PP][ihist] = bookHisto1D(namePP, refData(ihist+1, 1, 1));

        std::string nameCounterpp = "counter.pp." + std::to_string(ihist);
        _counterSOW[PP][ihist] = bookCounter(nameCounterpp,
          "Sum of weights counter for pp");

      }

      // Centrality regions keeping boundaries for a certain region.
      // Note, that some regions overlap with other regions.
      _centrRegions.clear();
      _centrRegions = {{0., 5.},   {5., 10.},  {10., 20.},
                       {20., 30.}, {30., 40.}, {40., 50.},
                       {50., 60.}, {60., 70.}, {70., 80.},
                       {0., 10.},  {0., 20.},  {20., 40.},
                       {40., 60.}, {40., 80.}, {60., 80.}};

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      // Charged, primary particles with at least pT = 150 MeV
      // in eta range of |eta| < 0.5
      Particles chargedParticles =
        applyProjection<ALICE::PrimaryParticles>(event,"APRIM").particlesByPt();

      // Check type of event. This may not be a perfect way to check for the
      // type of event as there might be some weird conditions hidden inside.
      // For example some HepMC versions check if number of hard collisions
      // is equal to 0 and assign 'false' in that case, which is usually wrong.
      // This might be changed in the future
      const HepMC::HeavyIon* hi = event.genEvent()->heavy_ion();
      if (hi && hi->is_valid()) {

        // Prepare centrality projection and value
        const CentralityProjection& centrProj =
          apply<CentralityProjection>(event, "V0M");
        double centr = centrProj();
        // Veto event for too large centralities since those are not used
        // in the analysis at all
        if ((centr < 0.) || (centr > 80.))
          vetoEvent;

        // Fill PbPb histograms and add weights based on centrality value
        for (size_t ihist = 0; ihist < NHISTOS; ++ihist) {
          if (inRange(centr, _centrRegions[ihist].first, _centrRegions[ihist].second)) {
            _counterSOW[PBPB][ihist]->fill(weight);
            _counterNcoll[ihist]->fill(event.genEvent()->heavy_ion()->Ncoll(), weight);
            foreach (const Particle& p, chargedParticles) {
              float pT = p.pT()/GeV;
              if (pT < 50.) {
                double pTAtBinCenter = _histNch[PBPB][ihist]->binAt(pT).xMid();
                _histNch[PBPB][ihist]->fill(pT, weight/pTAtBinCenter);
              }
            }
          }
        }

      }
      else {

        // Fill all pp histograms and add weights
        for (size_t ihist = 0; ihist < NHISTOS; ++ihist) {
          _counterSOW[PP][ihist]->fill(weight);
          foreach (const Particle& p, chargedParticles) {
            float pT = p.pT()/GeV;
            if (pT < 50.) {
              double pTAtBinCenter = _histNch[PP][ihist]->binAt(pT).xMid();
              _histNch[PP][ihist]->fill(pT, weight/pTAtBinCenter);
            }
          }
        }

      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // Right scaling of the histograms with their individual weights.
      for (size_t itype = 0; itype < EVENT_TYPES; ++itype ) {
        for (size_t ihist = 0; ihist < NHISTOS; ++ihist) {
          if (_counterSOW[itype][ihist]->sumW() > 0.) {
            scale(_histNch[itype][ihist],
              (1. / _counterSOW[itype][ihist]->sumW() / 2. / M_PI));
          }
        }
      }

      // Postprocessing of the histograms
      for (size_t ihist = 0; ihist < NHISTOS; ++ihist) {
        // If there are entires in histograms for both beam types
        if (_histNch[PP][ihist]->numEntries() > 0 && _histNch[PBPB][ihist]->numEntries() > 0) {
          // Initialize and fill R_AA histograms
          _histRAA[ihist] = bookScatter2D(ihist+16, 1, 1);
          divide(_histNch[PBPB][ihist], _histNch[PP][ihist], _histRAA[ihist]);
          // Scale by Ncoll. Unfortunately some generators does not provide
          // Ncoll value (eg. JEWEL), so the following scaling will be done
          // only if there are entries in the counters
          double ncoll = _counterNcoll[ihist]->sumW();
          double sow = _counterSOW[PBPB][ihist]->sumW();
          if (ncoll > 1e-6 && sow > 1e-6)
            _histRAA[ihist]->scaleY(1. / (ncoll / sow));

        }
      }

    }

    //@}

  private:

    static const int NHISTOS = 15;
    static const int EVENT_TYPES = 2;
    static const int PP = 0;
    static const int PBPB = 1;

    /// @name Histograms
    //@{
    Histo1DPtr _histNch[EVENT_TYPES][NHISTOS];
    CounterPtr _counterSOW[EVENT_TYPES][NHISTOS];
    CounterPtr _counterNcoll[NHISTOS];
    Scatter2DPtr _histRAA[NHISTOS];
    //@}

    std::vector<std::pair<double, double>> _centrRegions;

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2012_I1127497);


}
