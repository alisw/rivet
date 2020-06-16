// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ParisiTensor.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /**
   * @brief DELPHI event shapes and identified particle spectra
   * @author Andy Buckley
   * @author Hendrik Hoeth
   *
   * This is the paper which was used for the original PROFESSOR MC tuning
   * study. It studies a wide range of e+ e- event shape variables, differential
   * jet rates in the Durham and JADE schemes, and incorporates identified
   * particle spectra, from other LEP analyses.
   *
   * @par Run conditions
   *
   * @arg LEP1 beam energy: \f$ \sqrt{s} = \f$ 91.2 GeV
   * @arg Run with generic QCD events.
   * @arg No \f$ p_\perp^\text{min} \f$ cutoff is required
   */
  class DELPHI_1996_S3430090 : public Analysis {
  public:

    /// Constructor
    DELPHI_1996_S3430090()
      : Analysis("DELPHI_1996_S3430090")
    {    }


    /// @name Analysis methods
    //@{

    void init() {
      declare(Beam(), "Beams");
      // Don't try to introduce a pT or eta cut here. It's all corrected
      // back. (See Section 2 of the paper.)
      const ChargedFinalState cfs;
      declare(cfs, "FS");
      declare(UnstableParticles(), "UFS");
      declare(FastJets(cfs, FastJets::JADE, 0.7), "JadeJets");
      declare(FastJets(cfs, FastJets::DURHAM, 0.7), "DurhamJets");
      declare(Sphericity(cfs), "Sphericity");
      declare(ParisiTensor(cfs), "Parisi");
      const Thrust thrust(cfs);
      declare(thrust, "Thrust");
      declare(Hemispheres(thrust), "Hemispheres");

      book(_histPtTIn, 1, 1, 1);
      book(_histPtTOut,2, 1, 1);
      book(_histPtSIn, 3, 1, 1);
      book(_histPtSOut,4, 1, 1);

      book(_histRapidityT, 5, 1, 1);
      book(_histRapidityS, 6, 1, 1);
      book(_histScaledMom, 7, 1, 1);
      book(_histLogScaledMom, 8, 1, 1);

      book(_histPtTOutVsXp ,9,  1, 1);
      book(_histPtVsXp ,10, 1, 1);

      book(_hist1MinusT, 11, 1, 1);
      book(_histTMajor, 12, 1, 1);
      book(_histTMinor, 13, 1, 1);
      book(_histOblateness, 14, 1, 1);

      book(_histSphericity, 15, 1, 1);
      book(_histAplanarity, 16, 1, 1);
      book(_histPlanarity, 17, 1, 1);

      book(_histCParam, 18, 1, 1);
      book(_histDParam, 19, 1, 1);

      book(_histHemiMassH, 20, 1, 1);
      book(_histHemiMassL, 21, 1, 1);
      book(_histHemiMassD, 22, 1, 1);

      book(_histHemiBroadW, 23, 1, 1);
      book(_histHemiBroadN, 24, 1, 1);
      book(_histHemiBroadT, 25, 1, 1);
      book(_histHemiBroadD, 26, 1, 1);

      // Binned in y_cut
      book(_histDiffRate2Durham, 27, 1, 1);
      book(_histDiffRate2Jade, 28, 1, 1);
      book(_histDiffRate3Durham, 29, 1, 1);
      book(_histDiffRate3Jade, 30, 1, 1);
      book(_histDiffRate4Durham, 31, 1, 1);
      book(_histDiffRate4Jade, 32, 1, 1);

      // Binned in cos(chi)
      book(_histEEC, 33, 1, 1);
      book(_histAEEC, 34, 1, 1);

      book(_histMultiCharged, 35, 1, 1);

      book(_histMultiPiPlus, 36, 1, 1);
      book(_histMultiPi0, 36, 1, 2);
      book(_histMultiKPlus, 36, 1, 3);
      book(_histMultiK0, 36, 1, 4);
      book(_histMultiEta, 36, 1, 5);
      book(_histMultiEtaPrime, 36, 1, 6);
      book(_histMultiDPlus, 36, 1, 7);
      book(_histMultiD0, 36, 1, 8);
      book(_histMultiBPlus0, 36, 1, 9);

      book(_histMultiF0, 37, 1, 1);

      book(_histMultiRho, 38, 1, 1);
      book(_histMultiKStar892Plus, 38, 1, 2);
      book(_histMultiKStar892_0, 38, 1, 3);
      book(_histMultiPhi, 38, 1, 4);
      book(_histMultiDStar2010Plus, 38, 1, 5);

      book(_histMultiF2, 39, 1, 1);
      book(_histMultiK2Star1430_0, 39, 1, 2);

      book(_histMultiP, 40, 1, 1);
      book(_histMultiLambda0, 40, 1, 2);
      book(_histMultiXiMinus, 40, 1, 3);
      book(_histMultiOmegaMinus, 40, 1, 4);
      book(_histMultiDeltaPlusPlus, 40, 1, 5);
      book(_histMultiSigma1385Plus, 40, 1, 6);
      book(_histMultiXi1530_0, 40, 1, 7);
      book(_histMultiLambdaB0, 40, 1, 8);

      book(_weightedTotalPartNum,"_TotalPartNum");
      book(_passedCutWeightSum, "_passedCutWeightSum");
      book(_passedCut3WeightSum, "_passedCut3WeightSum");
      book(_passedCut4WeightSum, "_passedCut4WeightSum");
      book(_passedCut5WeightSum, "_passedCut5WeightSum");

    }



    void analyze(const Event& e) {
      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = apply<FinalState>(e, "FS");
      const size_t numParticles = fs.particles().size();
      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");

      _passedCutWeightSum->fill();
      _weightedTotalPartNum->fill(numParticles);

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(e, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      // Thrusts
      MSG_DEBUG("Calculating thrust");
      const Thrust& thrust = apply<Thrust>(e, "Thrust");
      _hist1MinusT->fill(1 - thrust.thrust());
      _histTMajor->fill(thrust.thrustMajor());
      _histTMinor->fill(thrust.thrustMinor());
      _histOblateness->fill(thrust.oblateness());

      // Jets
      const FastJets& durjet = apply<FastJets>(e, "DurhamJets");
      const FastJets& jadejet = apply<FastJets>(e, "JadeJets");
      if (numParticles >= 3) {
        _passedCut3WeightSum->fill();
        if (durjet.clusterSeq()) _histDiffRate2Durham->fill(durjet.clusterSeq()->exclusive_ymerge_max(2));
        if (jadejet.clusterSeq()) _histDiffRate2Jade->fill(jadejet.clusterSeq()->exclusive_ymerge_max(2));
      }
      if (numParticles >= 4) {
        _passedCut4WeightSum->fill();
        if (durjet.clusterSeq()) _histDiffRate3Durham->fill(durjet.clusterSeq()->exclusive_ymerge_max(3));
        if (jadejet.clusterSeq()) _histDiffRate3Jade->fill(jadejet.clusterSeq()->exclusive_ymerge_max(3));
      }
      if (numParticles >= 5) {
        _passedCut5WeightSum->fill();
        if (durjet.clusterSeq()) _histDiffRate4Durham->fill(durjet.clusterSeq()->exclusive_ymerge_max(4));
        if (jadejet.clusterSeq()) _histDiffRate4Jade->fill(jadejet.clusterSeq()->exclusive_ymerge_max(4));
      }

      // Sphericities
      MSG_DEBUG("Calculating sphericity");
      const Sphericity& sphericity = apply<Sphericity>(e, "Sphericity");
      _histSphericity->fill(sphericity.sphericity());
      _histAplanarity->fill(sphericity.aplanarity());
      _histPlanarity->fill(sphericity.planarity());

      // C & D params
      MSG_DEBUG("Calculating Parisi params");
      const ParisiTensor& parisi = apply<ParisiTensor>(e, "Parisi");
      _histCParam->fill(parisi.C());
      _histDParam->fill(parisi.D());

      // Hemispheres
      MSG_DEBUG("Calculating hemisphere variables");
      const Hemispheres& hemi = apply<Hemispheres>(e, "Hemispheres");
      _histHemiMassH->fill(hemi.scaledM2high());
      _histHemiMassL->fill(hemi.scaledM2low());
      _histHemiMassD->fill(hemi.scaledM2diff());
      _histHemiBroadW->fill(hemi.Bmax());
      _histHemiBroadN->fill(hemi.Bmin());
      _histHemiBroadT->fill(hemi.Bsum());
      _histHemiBroadD->fill(hemi.Bdiff());

      // Iterate over all the charged final state particles.
      double Evis = 0.0;
      double Evis2 = 0.0;
      MSG_DEBUG("About to iterate over charged FS particles");
      for (const Particle& p : fs.particles()) {
        // Get momentum and energy of each particle.
        const Vector3 mom3 = p.p3();
        const double energy = p.E();
        Evis += energy;

        // Scaled momenta.
        const double mom = mom3.mod();
        const double scaledMom = mom/meanBeamMom;
        const double logInvScaledMom = -std::log(scaledMom);
        _histLogScaledMom->fill(logInvScaledMom);
        _histScaledMom->fill(scaledMom);

        // Get momenta components w.r.t. thrust and sphericity.
        const double momT = dot(thrust.thrustAxis(), mom3);
        const double momS = dot(sphericity.sphericityAxis(), mom3);
        const double pTinT = dot(mom3, thrust.thrustMajorAxis());
        const double pToutT = dot(mom3, thrust.thrustMinorAxis());
        const double pTinS = dot(mom3, sphericity.sphericityMajorAxis());
        const double pToutS = dot(mom3, sphericity.sphericityMinorAxis());
        const double pT = sqrt(pow(pTinT, 2) + pow(pToutT, 2));
        _histPtTIn->fill(fabs(pTinT/GeV));
        _histPtTOut->fill(fabs(pToutT/GeV));
        _histPtSIn->fill(fabs(pTinS/GeV));
        _histPtSOut->fill(fabs(pToutS/GeV));
        _histPtVsXp->fill(scaledMom, fabs(pT/GeV));
        _histPtTOutVsXp->fill(scaledMom, fabs(pToutT/GeV));

        // Calculate rapidities w.r.t. thrust and sphericity.
        const double rapidityT = 0.5 * std::log((energy + momT) / (energy - momT));
        const double rapidityS = 0.5 * std::log((energy + momS) / (energy - momS));
        _histRapidityT->fill(fabs(rapidityT));
        _histRapidityS->fill(fabs(rapidityS));
        MSG_TRACE(fabs(rapidityT) << " " << scaledMom/GeV);
      }
      Evis2 = Evis*Evis;

      // (A)EEC
      // Need iterators since second loop starts at current outer loop iterator, i.e. no "for" here!
      for (Particles::const_iterator p_i = fs.particles().begin(); p_i != fs.particles().end(); ++p_i) {
        for (Particles::const_iterator p_j = p_i; p_j != fs.particles().end(); ++p_j) {
          if (p_i == p_j) continue;
          const Vector3 mom3_i = p_i->momentum().p3();
          const Vector3 mom3_j = p_j->momentum().p3();
          const double energy_i = p_i->momentum().E();
          const double energy_j = p_j->momentum().E();
          const double cosij = dot(mom3_i.unit(), mom3_j.unit());
          const double eec = (energy_i*energy_j) / Evis2;
          _histEEC->fill(cosij, eec);
          if (cosij < 0)
            _histAEEC->fill( cosij,  eec);
          else
            _histAEEC->fill(-cosij, -eec);
        }
      }

      _histMultiCharged->fill(_histMultiCharged->bin(0).xMid(), numParticles);


      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableFinalState>(e, "UFS");

      for (const Particle& p : ufs.particles()) {
        int id = p.abspid();
        switch (id) {
        case 211:
          _histMultiPiPlus->fill(_histMultiPiPlus->bin(0).xMid());
          break;
        case 111:
          _histMultiPi0->fill(_histMultiPi0->bin(0).xMid());
          break;
        case 321:
          _histMultiKPlus->fill(_histMultiKPlus->bin(0).xMid());
          break;
        case 130:
        case 310:
          _histMultiK0->fill(_histMultiK0->bin(0).xMid());
          break;
        case 221:
          _histMultiEta->fill(_histMultiEta->bin(0).xMid());
          break;
        case 331:
          _histMultiEtaPrime->fill(_histMultiEtaPrime->bin(0).xMid());
          break;
        case 411:
          _histMultiDPlus->fill(_histMultiDPlus->bin(0).xMid());
          break;
        case 421:
          _histMultiD0->fill(_histMultiD0->bin(0).xMid());
          break;
        case 511:
        case 521:
        case 531:
          _histMultiBPlus0->fill(_histMultiBPlus0->bin(0).xMid());
          break;
        case 9010221:
          _histMultiF0->fill(_histMultiF0->bin(0).xMid());
          break;
        case 113:
          _histMultiRho->fill(_histMultiRho->bin(0).xMid());
          break;
        case 323:
          _histMultiKStar892Plus->fill(_histMultiKStar892Plus->bin(0).xMid());
          break;
        case 313:
          _histMultiKStar892_0->fill(_histMultiKStar892_0->bin(0).xMid());
          break;
        case 333:
          _histMultiPhi->fill(_histMultiPhi->bin(0).xMid());
          break;
        case 413:
          _histMultiDStar2010Plus->fill(_histMultiDStar2010Plus->bin(0).xMid());
          break;
        case 225:
          _histMultiF2->fill(_histMultiF2->bin(0).xMid());
          break;
        case 315:
          _histMultiK2Star1430_0->fill(_histMultiK2Star1430_0->bin(0).xMid());
          break;
        case 2212:
          _histMultiP->fill(_histMultiP->bin(0).xMid());
          break;
        case 3122:
          _histMultiLambda0->fill(_histMultiLambda0->bin(0).xMid());
          break;
        case 3312:
          _histMultiXiMinus->fill(_histMultiXiMinus->bin(0).xMid());
          break;
        case 3334:
          _histMultiOmegaMinus->fill(_histMultiOmegaMinus->bin(0).xMid());
          break;
        case 2224:
          _histMultiDeltaPlusPlus->fill(_histMultiDeltaPlusPlus->bin(0).xMid());
          break;
        case 3114:
          _histMultiSigma1385Plus->fill(_histMultiSigma1385Plus->bin(0).xMid());
          break;
        case 3324:
          _histMultiXi1530_0->fill(_histMultiXi1530_0->bin(0).xMid());
          break;
        case 5122:
          _histMultiLambdaB0->fill(_histMultiLambdaB0->bin(0).xMid());
          break;
        }
      }
    }



    // Finalize
    void finalize() {
      // Normalize inclusive single particle distributions to the average number
      // of charged particles per event.
      const double avgNumParts = dbl(*_weightedTotalPartNum / *_passedCutWeightSum);

      normalize(_histPtTIn, avgNumParts);
      normalize(_histPtTOut, avgNumParts);
      normalize(_histPtSIn, avgNumParts);
      normalize(_histPtSOut, avgNumParts);

      normalize(_histRapidityT, avgNumParts);
      normalize(_histRapidityS, avgNumParts);

      normalize(_histLogScaledMom, avgNumParts);
      normalize(_histScaledMom, avgNumParts);

      scale(_histEEC, 1.0 / *_passedCutWeightSum);
      scale(_histAEEC, 1.0 / *_passedCutWeightSum);
      scale(_histMultiCharged, 1.0 / *_passedCutWeightSum);

      scale(_histMultiPiPlus, 1.0 / *_passedCutWeightSum);
      scale(_histMultiPi0, 1.0 / *_passedCutWeightSum);
      scale(_histMultiKPlus, 1.0 / *_passedCutWeightSum);
      scale(_histMultiK0, 1.0 / *_passedCutWeightSum);
      scale(_histMultiEta, 1.0 / *_passedCutWeightSum);
      scale(_histMultiEtaPrime, 1.0 / *_passedCutWeightSum);
      scale(_histMultiDPlus, 1.0 / *_passedCutWeightSum);
      scale(_histMultiD0, 1.0 / *_passedCutWeightSum);
      scale(_histMultiBPlus0, 1.0 / *_passedCutWeightSum);

      scale(_histMultiF0, 1.0 / *_passedCutWeightSum);

      scale(_histMultiRho, 1.0 / *_passedCutWeightSum);
      scale(_histMultiKStar892Plus, 1.0 / *_passedCutWeightSum);
      scale(_histMultiKStar892_0, 1.0 / *_passedCutWeightSum);
      scale(_histMultiPhi, 1.0 / *_passedCutWeightSum);
      scale(_histMultiDStar2010Plus, 1.0 / *_passedCutWeightSum);

      scale(_histMultiF2, 1.0 / *_passedCutWeightSum);
      scale(_histMultiK2Star1430_0, 1.0 / *_passedCutWeightSum);

      scale(_histMultiP, 1.0 / *_passedCutWeightSum);
      scale(_histMultiLambda0, 1.0 / *_passedCutWeightSum);
      scale(_histMultiXiMinus, 1.0 / *_passedCutWeightSum);
      scale(_histMultiOmegaMinus, 1.0 / *_passedCutWeightSum);
      scale(_histMultiDeltaPlusPlus, 1.0 / *_passedCutWeightSum);
      scale(_histMultiSigma1385Plus, 1.0 / *_passedCutWeightSum);
      scale(_histMultiXi1530_0, 1.0 / *_passedCutWeightSum);
      scale(_histMultiLambdaB0, 1.0 / *_passedCutWeightSum);

      scale(_hist1MinusT, 1.0 / *_passedCutWeightSum);
      scale(_histTMajor, 1.0 / *_passedCutWeightSum);
      scale(_histTMinor, 1.0 / *_passedCutWeightSum);
      scale(_histOblateness, 1.0 / *_passedCutWeightSum);

      scale(_histSphericity, 1.0 / *_passedCutWeightSum);
      scale(_histAplanarity, 1.0 / *_passedCutWeightSum);
      scale(_histPlanarity, 1.0 / *_passedCutWeightSum);

      scale(_histHemiMassD, 1.0 / *_passedCutWeightSum);
      scale(_histHemiMassH, 1.0 / *_passedCutWeightSum);
      scale(_histHemiMassL, 1.0 / *_passedCutWeightSum);

      scale(_histHemiBroadW, 1.0 / *_passedCutWeightSum);
      scale(_histHemiBroadN, 1.0 / *_passedCutWeightSum);
      scale(_histHemiBroadT, 1.0 / *_passedCutWeightSum);
      scale(_histHemiBroadD, 1.0 / *_passedCutWeightSum);

      scale(_histCParam, 1.0 / *_passedCutWeightSum);
      scale(_histDParam, 1.0 / *_passedCutWeightSum);

      scale(_histDiffRate2Durham, 1.0 / *_passedCut3WeightSum);
      scale(_histDiffRate2Jade, 1.0 / *_passedCut3WeightSum);
      scale(_histDiffRate3Durham, 1.0 / *_passedCut4WeightSum);
      scale(_histDiffRate3Jade, 1.0 / *_passedCut4WeightSum);
      scale(_histDiffRate4Durham, 1.0 / *_passedCut5WeightSum);
      scale(_histDiffRate4Jade, 1.0 / *_passedCut5WeightSum);
    }

    //@}


  private:

    /// Store the weighted sums of numbers of charged / charged+neutral
    /// particles - used to calculate average number of particles for the
    /// inclusive single particle distributions' normalisations.
    CounterPtr _weightedTotalPartNum;

    /// @name Sums of weights past various cuts
    //@{
    CounterPtr _passedCutWeightSum;
    CounterPtr _passedCut3WeightSum;
    CounterPtr _passedCut4WeightSum;
    CounterPtr _passedCut5WeightSum;
    //@}

    /// @name Histograms
    //@{
    Histo1DPtr _histPtTIn;
    Histo1DPtr _histPtTOut;
    Histo1DPtr _histPtSIn;
    Histo1DPtr _histPtSOut;

    Histo1DPtr _histRapidityT;
    Histo1DPtr _histRapidityS;

    Histo1DPtr _histScaledMom, _histLogScaledMom;

    Profile1DPtr   _histPtTOutVsXp, _histPtVsXp;

    Histo1DPtr _hist1MinusT;
    Histo1DPtr _histTMajor;
    Histo1DPtr _histTMinor;
    Histo1DPtr _histOblateness;

    Histo1DPtr _histSphericity;
    Histo1DPtr _histAplanarity;
    Histo1DPtr _histPlanarity;

    Histo1DPtr _histCParam;
    Histo1DPtr _histDParam;

    Histo1DPtr _histHemiMassD;
    Histo1DPtr _histHemiMassH;
    Histo1DPtr _histHemiMassL;

    Histo1DPtr _histHemiBroadW;
    Histo1DPtr _histHemiBroadN;
    Histo1DPtr _histHemiBroadT;
    Histo1DPtr _histHemiBroadD;

    Histo1DPtr _histDiffRate2Durham;
    Histo1DPtr _histDiffRate2Jade;
    Histo1DPtr _histDiffRate3Durham;
    Histo1DPtr _histDiffRate3Jade;
    Histo1DPtr _histDiffRate4Durham;
    Histo1DPtr _histDiffRate4Jade;

    Histo1DPtr _histEEC, _histAEEC;

    Histo1DPtr _histMultiCharged;

    Histo1DPtr _histMultiPiPlus;
    Histo1DPtr _histMultiPi0;
    Histo1DPtr _histMultiKPlus;
    Histo1DPtr _histMultiK0;
    Histo1DPtr _histMultiEta;
    Histo1DPtr _histMultiEtaPrime;
    Histo1DPtr _histMultiDPlus;
    Histo1DPtr _histMultiD0;
    Histo1DPtr _histMultiBPlus0;

    Histo1DPtr _histMultiF0;

    Histo1DPtr _histMultiRho;
    Histo1DPtr _histMultiKStar892Plus;
    Histo1DPtr _histMultiKStar892_0;
    Histo1DPtr _histMultiPhi;
    Histo1DPtr _histMultiDStar2010Plus;

    Histo1DPtr _histMultiF2;
    Histo1DPtr _histMultiK2Star1430_0;

    Histo1DPtr _histMultiP;
    Histo1DPtr _histMultiLambda0;
    Histo1DPtr _histMultiXiMinus;
    Histo1DPtr _histMultiOmegaMinus;
    Histo1DPtr _histMultiDeltaPlusPlus;
    Histo1DPtr _histMultiSigma1385Plus;
    Histo1DPtr _histMultiXi1530_0;
    Histo1DPtr _histMultiLambdaB0;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(DELPHI_1996_S3430090);

}
