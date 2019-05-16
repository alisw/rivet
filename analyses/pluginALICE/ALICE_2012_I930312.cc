// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"

namespace Rivet {

  /// @brief ALICE PbPb at 2.76 TeV azimuthal di-hadron correlations
  class ALICE_2012_I930312 : public Analysis {

  public:

    // Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2012_I930312);

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Declare centrality projection
      declareCentrality(ALICE::V0MMultiplicity(),
        "ALICE_2015_PBPBCentrality", "V0M", "V0M");

      // Projection for trigger particles: charged, primary particles
      // with |eta| < 1.0 and 8 < pT < 15 GeV/c
      declare(ALICE::PrimaryParticles(Cuts::abseta < 1.0 && Cuts::abscharge > 0
        && Cuts::ptIn(8.*GeV, 15.*GeV)), "APRIMTrig");

      // pT bins edges
      vector<double> ptBins = { 3., 4., 6., 8., 10. };

      // Projections for associated particles: charged, primary particles
      // with |eta| < 1.0 and different pT bins
      for (int ipt = 0; ipt < PT_BINS; ++ipt) {
        Cut cut = Cuts::abseta < 1.0 && Cuts::abscharge > 0 &&
          Cuts::ptIn(ptBins[ipt]*GeV, ptBins[ipt+1]*GeV);
        declare(ALICE::PrimaryParticles(cut), "APRIMAssoc" + toString(ipt));
      }

      // Create event strings
      vector<string> evString = { "pp", "central", "peripheral" };

      // Initialize trigger counters and yield histograms
      string title = "Per trigger particle yield";
      string xtitle = "$\\Delta\\eta$ (rad)";
      string ytitle =
        "$1 / N_{trig} {\\rm d}N_{assoc} / {\\rm d}\\Delta\\eta$ (rad$^-1$)";
      for (int itype = 0; itype < EVENT_TYPES; ++itype) {
        _counterTrigger[itype] = bookCounter("counter." + toString(itype));
        for (int ipt = 0; ipt < PT_BINS; ++ipt) {
          string name = "yield." + evString[itype] + ".pt" + toString(ipt);
          _histYield[itype][ipt] = bookHisto1D(name, 36,
            -0.5*M_PI, 1.5*M_PI, title, xtitle, ytitle);
        }
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      // Trigger particles
      Particles trigParticles =
        applyProjection<ALICE::PrimaryParticles>(event,"APRIMTrig").particles();

      // Associated particles
      Particles assocParticles[PT_BINS];
      for (int ipt = 0; ipt < PT_BINS; ++ipt) {
        string pname = "APRIMAssoc" + toString(ipt);
        assocParticles[ipt] =
          applyProjection<ALICE::PrimaryParticles>(event,pname).particles();
      }

      // Check type of event. This may not be a perfect way to check for the
      // type of event as there might be some weird conditions hidden inside.
      // For example some HepMC versions check if number of hard collisions
      // is equal to 0 and assign 'false' in that case, which is usually wrong.
      // This might be changed in the future
      int ev_type = 0; // pp
      const HepMC::HeavyIon* hi = event.genEvent()->heavy_ion();
      if (hi && hi->is_valid()) {
        // Prepare centrality projection and value
        const CentralityProjection& centrProj =
          apply<CentralityProjection>(event, "V0M");
        double centr = centrProj();
        // Set the flag for the type of the event
        if (centr > 0.0 && centr < 5.0)
          ev_type = 1; // PbPb, central
        else if (centr > 60.0 && centr < 90.0)
          ev_type = 2; // PbPb, peripherial
        else
          vetoEvent; // PbPb, other, this is not used in the analysis at all
      }

      // Fill trigger histogram for a proper event type
      _counterTrigger[ev_type]->fill(trigParticles.size());

      // Loop over trigger particles
      for (const Particle& trigParticle : trigParticles) {
        // For each pt bin
        for (int ipt = 0; ipt < PT_BINS; ++ipt) {
          // Loop over associated particles
          for (const Particle& assocParticle : assocParticles[ipt]) {
            // If associated and trigger particle are not the same particles.
            if (!isSame(trigParticle, assocParticle)) {
              // Test trigger particle.
              if (trigParticle.pt() > assocParticle.pt()) {
                // Calculate delta phi in range (-0.5*PI, 1.5*PI).
                double dPhi = deltaPhi(trigParticle, assocParticle, true);
                if (dPhi < -0.5 * M_PI) dPhi += 2 * M_PI;
                // Fill yield histogram for calculated delta phi
                _histYield[ev_type][ipt]->fill(dPhi, weight);
              }
            }
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // Check for the reentrant finalize
      bool pp_available = false, PbPb_available = false;
      for (int itype = 0; itype < EVENT_TYPES; ++itype) {
        for (int ipt = 0; ipt < PT_BINS; ++ipt) {
          if (_histYield[itype][ipt]->numEntries() > 0)
            itype == 0 ? pp_available = true : PbPb_available = true;
        }
      }
      // Skip postprocessing if pp or PbPb histograms are not available
      if (!(pp_available && PbPb_available))
        return;

      // Initialize IAA and ICP histograms
      _histIAA[0] = bookScatter2D(1, 1, 1);
      _histIAA[1] = bookScatter2D(2, 1, 1);
      _histIAA[2] = bookScatter2D(5, 1, 1);
      _histIAA[3] = bookScatter2D(3, 1, 1);
      _histIAA[4] = bookScatter2D(4, 1, 1);
      _histIAA[5] = bookScatter2D(6, 1, 1);

      // Initialize background-subtracted yield histograms
      for (int itype = 0; itype < EVENT_TYPES; ++itype) {
        for (int ipt = 0; ipt < PT_BINS; ++ipt) {
          string newname = _histYield[itype][ipt]->name() + ".nobkg";
          string newtitle = _histYield[itype][ipt]->title() +
            ", background subtracted";
          _histYieldNoBkg[itype][ipt] =
            bookHisto1D(newname, 36, -0.5*M_PI, 1.5*M_PI, newtitle,
            _histYield[itype][ipt]->annotation("XLabel"),
            _histYield[itype][ipt]->annotation("YLabel"));
        }
      }

      // Variable for near and away side peak integral calculation
      double integral[EVENT_TYPES][PT_BINS][2] = { { {0.0} } };

      // Variables for background calculation
      double bkg = 0.0;
      double bkgErr[EVENT_TYPES][PT_BINS] = { {0.0} };

      // Variables for integration error calculation
      double norm[EVENT_TYPES] = {0.0};
      double numEntries[EVENT_TYPES][PT_BINS][2] = { { {0.0} } };
      int numBins[EVENT_TYPES][PT_BINS][2] = { { {0} } };

      // For each event type
      for (int itype = 0; itype < EVENT_TYPES; ++itype) {
        // Get counter
        CounterPtr counter = _counterTrigger[itype];
        // For each pT range
        for (int ipt = 0; ipt < PT_BINS; ++ipt) {

          // Get yield histograms
          Histo1DPtr hYield = _histYield[itype][ipt];
          Histo1DPtr hYieldNoBkg = _histYieldNoBkg[itype][ipt];

          // Check if histograms are fine
          if (counter->sumW() == 0 || hYield->numEntries() == 0) {
            MSG_WARNING("There are no entries in one of the histograms");
            continue;
          }

          // Scale yield histogram
          norm[itype] = 1. / counter->sumW();
          scale(hYield, norm[itype]);

          // Calculate background
          double sum = 0.0;
          int nbins = 0;
          for (size_t ibin = 0; ibin < hYield->numBins(); ++ibin) {
            double xmid = hYield->bin(ibin).xMid();
            if (inRange(xmid, -0.5 * M_PI, -0.5 * M_PI + 0.4) ||
                inRange(xmid, 0.5 * M_PI - 0.4, 0.5 * M_PI + 0.4) ||
                inRange(xmid, 1.5 * M_PI - 0.4, 1.5 * M_PI)) {
              sum += hYield->bin(ibin).sumW();
              nbins += 1;
            }
          }
          if (nbins == 0) {
            MSG_WARNING("Failed to estimate background!");
            continue;
          }
          bkg = sum / nbins;

          // Calculate background error
          sum = 0.0;
          nbins = 0;
          for (size_t ibin = 0; ibin < hYield->numBins(); ++ibin) {
            double xmid = hYield->bin(ibin).xMid();
            if (inRange(xmid, 0.5 * M_PI - 0.4, 0.5 * M_PI + 0.4)) {
              sum += (hYield->bin(ibin).sumW() - bkg) *
                     (hYield->bin(ibin).sumW() - bkg);
              nbins++;
            }
          }
          if (nbins < 2) {
            MSG_WARNING("Failed to estimate background error!");
            continue;
          }
          bkgErr[itype][ipt] = sqrt(sum / (nbins - 1));

          // Fill histograms with removed background
          for (size_t ibin = 0; ibin < hYield->numBins(); ++ibin) {
            hYieldNoBkg->fillBin(ibin, hYield->bin(ibin).sumW() - bkg);
          }

          // Integrate near-side yield
          size_t lowerBin = hYield->binIndexAt(-0.7 + 0.02);
          size_t upperBin = hYield->binIndexAt( 0.7 - 0.02) + 1;
          nbins = upperBin - lowerBin;
          numBins[itype][ipt][NEAR] = nbins;
          integral[itype][ipt][NEAR] =
            hYield->integralRange(lowerBin, upperBin) - nbins * bkg;
          numEntries[itype][ipt][NEAR] =
            hYield->integralRange(lowerBin, upperBin) * counter->sumW();

          // Integrate away-side yield
          lowerBin = hYield->binIndexAt(M_PI - 0.7 + 0.02);
          upperBin = hYield->binIndexAt(M_PI + 0.7 - 0.02) + 1;
          nbins = upperBin - lowerBin;
          numBins[itype][ipt][AWAY] = nbins;
          integral[itype][ipt][AWAY] =
            hYield->integralRange(lowerBin, upperBin) - nbins * bkg;
          numEntries[itype][ipt][AWAY] =
            hYield->integralRange(lowerBin, upperBin) * counter->sumW();

        }
      }

      // Variables for IAA/ICP plots
      double yval[2] = { 0.0, 0.0 };
      double yerr[2] = { 0.0, 0.0 };
      double xval[PT_BINS] = { 3.5, 5.0, 7.0, 9.0 };
      double xerr[PT_BINS] = { 0.5, 1.0, 1.0, 1.0 };

      int types1[3] = {1, 2, 1};
      int types2[3] = {0, 0, 2};

      // Fill IAA/ICP plots for near and away side peak
      for (int ihist = 0; ihist < 3; ++ihist) {
        int type1 = types1[ihist];
        int type2 = types2[ihist];
        double norm1 = norm[type1];
        double norm2 = norm[type2];
        for (int ipt = 0; ipt < PT_BINS; ++ipt) {
          double bkgErr1 = bkgErr[type1][ipt];
          double bkgErr2 = bkgErr[type2][ipt];
          for (int ina = 0; ina < 2; ++ina) {
            double integ1 = integral[type1][ipt][ina];
            double integ2 = integral[type2][ipt][ina];
            double numEntries1 = numEntries[type1][ipt][ina];
            double numEntries2 = numEntries[type2][ipt][ina];
            double numBins1 = numBins[type1][ipt][ina];
            double numBins2 = numBins[type2][ipt][ina];
            yval[ina] = integ1 / integ2;
            yerr[ina] = norm1 * norm1 * numEntries1 +
              norm2 * norm2 * numEntries2 * integ1 * integ1 / (integ2 * integ2) +
              numBins1 * numBins1 * bkgErr1 * bkgErr1 +
              numBins2 * numBins2 * bkgErr2 * bkgErr2 * integ1 * integ1 / (integ2 * integ2);
            yerr[ina] = sqrt(yerr[ina])/integ2;
          }
          _histIAA[ihist]->addPoint(xval[ipt], yval[NEAR], xerr[ipt], yerr[NEAR]);
          _histIAA[ihist + 3]->addPoint(xval[ipt], yval[AWAY], xerr[ipt], yerr[AWAY]);
        }
      }

    }

    //@}

  private:

    static const int PT_BINS = 4;
    static const int EVENT_TYPES = 3;
    static const int NEAR = 0;
    static const int AWAY = 1;

    /// @name Histograms
    //@{
    Histo1DPtr _histYield[EVENT_TYPES][PT_BINS];
    Histo1DPtr _histYieldNoBkg[EVENT_TYPES][PT_BINS];
    CounterPtr _counterTrigger[EVENT_TYPES];
    Scatter2DPtr _histIAA[6];
    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2012_I930312);

}
