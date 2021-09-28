// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ParisiTensor.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include <cmath>

namespace Rivet {


  /// @brief OPAL event shapes and moments at 91, 133, 177, and 197 GeV
  /// @author Andy Buckley
  class OPAL_2004_S6132243 : public Analysis {
  public:

    /// Constructor
    OPAL_2004_S6132243()
      : Analysis("OPAL_2004_S6132243"),
        _isqrts(-1)
    {
      //
    }


    /// @name Analysis methods
    //@{

    /// Energies: 91, 133, 177 (161-183), 197 (189-209) => index 0..4
    int getHistIndex(double sqrts) {
      int ih = -1;
      if (inRange(sqrts/GeV, 89.9, 91.5)) {
        ih = 0;
      } else if (fuzzyEquals(sqrts/GeV, 133)) {
        ih = 1;
      } else if (fuzzyEquals(sqrts/GeV, 177)) { // (161-183)
        ih = 2;
      } else if (fuzzyEquals(sqrts/GeV, 197)) { // (189-209)
        ih = 3;
      } else {
        stringstream ss;
        ss << "Invalid energy for OPAL_2004 analysis: "
           << sqrts/GeV << " GeV != 91, 133, 177, or 197 GeV";
        throw Error(ss.str());
      }
      assert(ih >= 0);
      return ih;
    }


    void init() {
      // Projections
      declare(Beam(), "Beams");
      const FinalState fs;
      declare(fs, "FS");
      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      declare(FastJets(fs, FastJets::DURHAM, 0.7), "DurhamJets");
      declare(Sphericity(fs), "Sphericity");
      declare(ParisiTensor(fs), "Parisi");
      const Thrust thrust(fs);
      declare(thrust, "Thrust");
      declare(Hemispheres(thrust), "Hemispheres");

      // Get beam energy index
      _isqrts = getHistIndex(sqrtS());

      // Book histograms
      book(_hist1MinusT[_isqrts]    ,1, 1, _isqrts+1);
      book(_histHemiMassH[_isqrts]  ,2, 1, _isqrts+1);
      book(_histCParam[_isqrts]     ,3, 1, _isqrts+1);
      book(_histHemiBroadT[_isqrts] ,4, 1, _isqrts+1);
      book(_histHemiBroadW[_isqrts] ,5, 1, _isqrts+1);
      book(_histY23Durham[_isqrts]  ,6, 1, _isqrts+1);
      book(_histTMajor[_isqrts]     ,7, 1, _isqrts+1);
      book(_histTMinor[_isqrts]     ,8, 1, _isqrts+1);
      book(_histAplanarity[_isqrts] ,9, 1, _isqrts+1);
      book(_histSphericity[_isqrts] ,10, 1, _isqrts+1);
      book(_histOblateness[_isqrts] ,11, 1, _isqrts+1);
      book(_histHemiMassL[_isqrts]  ,12, 1, _isqrts+1);
      book(_histHemiBroadN[_isqrts] ,13, 1, _isqrts+1);
      book(_histDParam[_isqrts]     ,14, 1, _isqrts+1);
      //
      book(_hist1MinusTMom[_isqrts]    ,15, 1, _isqrts+1);
      book(_histHemiMassHMom[_isqrts]  ,16, 1, _isqrts+1);
      book(_histCParamMom[_isqrts]     ,17, 1, _isqrts+1);
      book(_histHemiBroadTMom[_isqrts] ,18, 1, _isqrts+1);
      book(_histHemiBroadWMom[_isqrts] ,19, 1, _isqrts+1);
      book(_histY23DurhamMom[_isqrts]  ,20, 1, _isqrts+1);
      book(_histTMajorMom[_isqrts]     ,21, 1, _isqrts+1);
      book(_histTMinorMom[_isqrts]     ,22, 1, _isqrts+1);
      book(_histSphericityMom[_isqrts] ,23, 1, _isqrts+1);
      book(_histOblatenessMom[_isqrts] ,24, 1, _isqrts+1);
      book(_histHemiMassLMom[_isqrts]  ,25, 1, _isqrts+1);
      book(_histHemiBroadNMom[_isqrts] ,26, 1, _isqrts+1);

      book(_sumWTrack2, "_sumWTrack2");
      book(_sumWJet3, "_sumWJet3");
    }


    void analyze(const Event& event) {
      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      if (cfs.size() < 2) vetoEvent;

      _sumWTrack2->fill();

      // Thrusts
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      _hist1MinusT[_isqrts]->fill(1-thrust.thrust());
      _histTMajor[_isqrts]->fill(thrust.thrustMajor());
      _histTMinor[_isqrts]->fill(thrust.thrustMinor());
      _histOblateness[_isqrts]->fill(thrust.oblateness());
      for (int n = 1; n <= 5; ++n) {
        _hist1MinusTMom[_isqrts]->fill(n, pow(1-thrust.thrust(), n));
        _histTMajorMom[_isqrts]->fill(n, pow(thrust.thrustMajor(), n));
        _histTMinorMom[_isqrts]->fill(n, pow(thrust.thrustMinor(), n));
        _histOblatenessMom[_isqrts]->fill(n, pow(thrust.oblateness(), n));
      }

      // Jets
      const FastJets& durjet = apply<FastJets>(event, "DurhamJets");
      if (durjet.clusterSeq()) {
        _sumWJet3->fill();
        const double y23 = durjet.clusterSeq()->exclusive_ymerge_max(2);
        if (y23>0.0) {
          _histY23Durham[_isqrts]->fill(y23);
          for (int n = 1; n <= 5; ++n) {
            _histY23DurhamMom[_isqrts]->fill(n, pow(y23, n));
          }
        }
      }

      // Sphericities
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");
      const double sph = sphericity.sphericity();
      const double apl = sphericity.aplanarity();
      _histSphericity[_isqrts]->fill(sph);
      _histAplanarity[_isqrts]->fill(apl);
      for (int n = 1; n <= 5; ++n) {
        _histSphericityMom[_isqrts]->fill(n, pow(sph, n));
      }

      // C & D params
      const ParisiTensor& parisi = apply<ParisiTensor>(event, "Parisi");
      const double cparam = parisi.C();
      const double dparam = parisi.D();
      _histCParam[_isqrts]->fill(cparam);
      _histDParam[_isqrts]->fill(dparam);
      for (int n = 1; n <= 5; ++n) {
        _histCParamMom[_isqrts]->fill(n, pow(cparam, n));
      }

      // Hemispheres
      const Hemispheres& hemi = apply<Hemispheres>(event, "Hemispheres");
      // The paper says that M_H/L are scaled by sqrt(s), but scaling by E_vis is the way that fits the data...
      const double hemi_mh = hemi.scaledMhigh();
      const double hemi_ml = hemi.scaledMlow();
      /// @todo This shouldn't be necessary... what's going on? Memory corruption suspected :(
      // if (std::isnan(hemi_ml)) {
      //   MSG_ERROR("NaN in HemiL! Event = " << numEvents());
      //   MSG_ERROR(hemi.M2low() << ", " << hemi.E2vis());
      // }
      if (!std::isnan(hemi_mh) && !std::isnan(hemi_ml)) {
        const double hemi_bmax = hemi.Bmax();
        const double hemi_bmin = hemi.Bmin();
        const double hemi_bsum = hemi.Bsum();
        _histHemiMassH[_isqrts]->fill(hemi_mh);
        _histHemiMassL[_isqrts]->fill(hemi_ml);
        _histHemiBroadW[_isqrts]->fill(hemi_bmax);
        _histHemiBroadN[_isqrts]->fill(hemi_bmin);
        _histHemiBroadT[_isqrts]->fill(hemi_bsum);
        for (int n = 1; n <= 5; ++n) {
          // if (std::isnan(pow(hemi_ml, n))) MSG_ERROR("NaN in HemiL moment! Event = " << numEvents());
          _histHemiMassHMom[_isqrts]->fill(n, pow(hemi_mh, n));
          _histHemiMassLMom[_isqrts]->fill(n, pow(hemi_ml, n));
          _histHemiBroadWMom[_isqrts]->fill(n, pow(hemi_bmax, n));
          _histHemiBroadNMom[_isqrts]->fill(n, pow(hemi_bmin, n));
          _histHemiBroadTMom[_isqrts]->fill(n, pow(hemi_bsum, n));
        }
      }
    }


    void finalize() {
      scale(_hist1MinusT[_isqrts], 1.0 / *_sumWTrack2);
      scale(_histTMajor[_isqrts], 1.0 / *_sumWTrack2);
      scale(_histTMinor[_isqrts], 1.0 / *_sumWTrack2);
      scale(_histOblateness[_isqrts], 1.0 / *_sumWTrack2);
      scale(_histSphericity[_isqrts], 1.0 / *_sumWTrack2);
      scale(_histAplanarity[_isqrts], 1.0 / *_sumWTrack2);
      scale(_histHemiMassH[_isqrts], 1.0 / *_sumWTrack2);
      scale(_histHemiMassL[_isqrts], 1.0 / *_sumWTrack2);
      scale(_histHemiBroadW[_isqrts], 1.0 / *_sumWTrack2);
      scale(_histHemiBroadN[_isqrts], 1.0 / *_sumWTrack2);
      scale(_histHemiBroadT[_isqrts], 1.0 / *_sumWTrack2);
      scale(_histCParam[_isqrts], 1.0 / *_sumWTrack2);
      scale(_histDParam[_isqrts], 1.0 / *_sumWTrack2);
      scale(_histY23Durham[_isqrts], 1.0 / *_sumWJet3);
      //
      scale(_hist1MinusTMom[_isqrts], 1.0 / *_sumWTrack2);
      scale(_histTMajorMom[_isqrts], 1.0 / *_sumWTrack2);
      scale(_histTMinorMom[_isqrts], 1.0 / *_sumWTrack2);
      scale(_histOblatenessMom[_isqrts], 1.0 / *_sumWTrack2);
      scale(_histSphericityMom[_isqrts], 1.0 / *_sumWTrack2);
      scale(_histHemiMassHMom[_isqrts], 1.0 / *_sumWTrack2);
      scale(_histHemiMassLMom[_isqrts], 1.0 / *_sumWTrack2);
      scale(_histHemiBroadWMom[_isqrts], 1.0 / *_sumWTrack2);
      scale(_histHemiBroadNMom[_isqrts], 1.0 / *_sumWTrack2);
      scale(_histHemiBroadTMom[_isqrts], 1.0 / *_sumWTrack2);
      scale(_histCParamMom[_isqrts], 1.0 / *_sumWTrack2);
      scale(_histY23DurhamMom[_isqrts], 1.0 / *_sumWJet3);
    }

    //@}


  private:

    /// Beam energy index for histograms
    int _isqrts;

    /// @name Counters of event weights passing the cuts
    //@{
    CounterPtr _sumWTrack2, _sumWJet3;
    //@}

    /// @name Event shape histos at 4 energies
    //@{
    Histo1DPtr _hist1MinusT[4];
    Histo1DPtr _histHemiMassH[4];
    Histo1DPtr _histCParam[4];
    Histo1DPtr _histHemiBroadT[4];
    Histo1DPtr _histHemiBroadW[4];
    Histo1DPtr _histY23Durham[4];
    Histo1DPtr _histTMajor[4];
    Histo1DPtr _histTMinor[4];
    Histo1DPtr _histAplanarity[4];
    Histo1DPtr _histSphericity[4];
    Histo1DPtr _histOblateness[4];
    Histo1DPtr _histHemiMassL[4];
    Histo1DPtr _histHemiBroadN[4];
    Histo1DPtr _histDParam[4];
    //@}

    /// @name Event shape moment histos at 4 energies
    //@{
    Histo1DPtr _hist1MinusTMom[4];
    Histo1DPtr _histHemiMassHMom[4];
    Histo1DPtr _histCParamMom[4];
    Histo1DPtr _histHemiBroadTMom[4];
    Histo1DPtr _histHemiBroadWMom[4];
    Histo1DPtr _histY23DurhamMom[4];
    Histo1DPtr _histTMajorMom[4];
    Histo1DPtr _histTMinorMom[4];
    Histo1DPtr _histSphericityMom[4];
    Histo1DPtr _histOblatenessMom[4];
    Histo1DPtr _histHemiMassLMom[4];
    Histo1DPtr _histHemiBroadNMom[4];
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(OPAL_2004_S6132243);

}
