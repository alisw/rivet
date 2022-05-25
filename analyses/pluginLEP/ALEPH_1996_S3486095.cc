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


  /// @brief ALEPH QCD study with event shapes and identified particles
  ///
  /// @author Holger Schulz
  class ALEPH_1996_S3486095 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(ALEPH_1996_S3486095);


    /// @name Analysis methods
    /// @{

    void init() {
      // Set up projections
      declare(Beam(), "Beams");
      const ChargedFinalState cfs;
      declare(cfs, "FS");
      declare(UnstableParticles(), "UFS");
      declare(FastJets(cfs, FastJets::DURHAM, 0.7), "DurhamJets");
      declare(Sphericity(cfs), "Sphericity");
      declare(ParisiTensor(cfs), "Parisi");
      const Thrust thrust(cfs);
      declare(thrust, "Thrust");
      declare(Hemispheres(thrust), "Hemispheres");

      // Book histograms
      book(_histSphericity   ,1, 1, 1);
      book(_histAplanarity   ,2, 1, 1);

      book(_hist1MinusT      ,3, 1, 1);
      book(_histTMinor       ,4, 1, 1);

      book(_histY3           ,5, 1, 1);
      book(_histHeavyJetMass ,6, 1, 1);
      book(_histCParam       ,7, 1, 1);
      book(_histOblateness   ,8, 1, 1);

      book(_histScaledMom    ,9, 1, 1);
      book(_histRapidityT    ,10, 1, 1);

      book(_histPtSIn        ,11, 1, 1);
      book(_histPtSOut       ,12, 1, 1);

      book(_histLogScaledMom ,17, 1, 1);

      book(_histChMult       ,18, 1, 1);
      book(_histMeanChMult   ,19, 1, 1);

      book(_histMeanChMultRapt05,20, 1, 1);
      book(_histMeanChMultRapt10,21, 1, 1);
      book(_histMeanChMultRapt15,22, 1, 1);
      book(_histMeanChMultRapt20,23, 1, 1);


      // Particle spectra
      book(_histMultiPiPlus        ,25, 1, 1);
      book(_histMultiKPlus         ,26, 1, 1);
      book(_histMultiP             ,27, 1, 1);
      book(_histMultiPhoton        ,28, 1, 1);
      book(_histMultiPi0           ,29, 1, 1);
      book(_histMultiEta           ,30, 1, 1);
      book(_histMultiEtaPrime      ,31, 1, 1);
      book(_histMultiK0            ,32, 1, 1);
      book(_histMultiLambda0       ,33, 1, 1);
      book(_histMultiXiMinus       ,34, 1, 1);
      book(_histMultiSigma1385Plus ,35, 1, 1);
      book(_histMultiXi1530_0      ,36, 1, 1);
      book(_histMultiRho           ,37, 1, 1);
      book(_histMultiOmega782      ,38, 1, 1);
      book(_histMultiKStar892_0    ,39, 1, 1);
      book(_histMultiPhi           ,40, 1, 1);

      book(_histMultiKStar892Plus  ,43, 1, 1);

      // Mean multiplicities
      book(_histMeanMultiPi0           ,44, 1,  2);
      book(_histMeanMultiEta           ,44, 1,  3);
      book(_histMeanMultiEtaPrime      ,44, 1,  4);
      book(_histMeanMultiK0            ,44, 1,  5);
      book(_histMeanMultiRho           ,44, 1,  6);
      book(_histMeanMultiOmega782      ,44, 1,  7);
      book(_histMeanMultiPhi           ,44, 1,  8);
      book(_histMeanMultiKStar892Plus  ,44, 1,  9);
      book(_histMeanMultiKStar892_0    ,44, 1, 10);
      book(_histMeanMultiLambda0       ,44, 1, 11);
      book(_histMeanMultiSigma0        ,44, 1, 12);
      book(_histMeanMultiXiMinus       ,44, 1, 13);
      book(_histMeanMultiSigma1385Plus ,44, 1, 14);
      book(_histMeanMultiXi1530_0      ,44, 1, 15);
      book(_histMeanMultiOmegaOmegaBar ,44, 1, 16);
      book(_weightedTotalPartNum, "/TMP/TotalPartNum");

      book(_weightedTotalPartNum, "/TMP/weightedTotalPartNum");
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
      _histTMinor->fill(thrust.thrustMinor());
      _histOblateness->fill(thrust.oblateness());

      // Jets
      MSG_DEBUG("Calculating differential jet rate plots:");
      const FastJets& durjet = apply<FastJets>(e, "DurhamJets");
      if (durjet.clusterSeq()) {
        double y3 = durjet.clusterSeq()->exclusive_ymerge_max(2);
        if (y3>0.0) _histY3->fill(-1. * std::log(y3));
      }

      // Sphericities
      MSG_DEBUG("Calculating sphericity");
      const Sphericity& sphericity = apply<Sphericity>(e, "Sphericity");
      _histSphericity->fill(sphericity.sphericity());
      _histAplanarity->fill(sphericity.aplanarity());

      // C param
      MSG_DEBUG("Calculating Parisi params");
      const ParisiTensor& parisi = apply<ParisiTensor>(e, "Parisi");
      _histCParam->fill(parisi.C());

      // Hemispheres
      MSG_DEBUG("Calculating hemisphere variables");
      const Hemispheres& hemi = apply<Hemispheres>(e, "Hemispheres");
      _histHeavyJetMass->fill(hemi.scaledM2high());

      // Iterate over all the charged final state particles.
      double Evis = 0.0;
      double rapt05 = 0.;
      double rapt10 = 0.;
      double rapt15 = 0.;
      double rapt20 = 0.;
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
        const double pTinS = dot(mom3, sphericity.sphericityMajorAxis());
        const double pToutS = dot(mom3, sphericity.sphericityMinorAxis());
        _histPtSIn->fill(fabs(pTinS/GeV));
        _histPtSOut->fill(fabs(pToutS/GeV));

        // Calculate rapidities w.r.t. thrust.
        const double rapidityT = 0.5 * std::log((energy + momT) / (energy - momT));
        _histRapidityT->fill(fabs(rapidityT));
        if (std::fabs(rapidityT) <= 0.5)  {
            rapt05 += 1.0;
        }
        if (std::fabs(rapidityT) <= 1.0)  {
            rapt10 += 1.0;
        }
        if (std::fabs(rapidityT) <= 1.5) {
            rapt15 += 1.0;
        }
        if (std::fabs(rapidityT) <= 2.0)  {
            rapt20 += 1.0;
        }

      }

      _histChMult->fill(numParticles);

      _histMeanChMultRapt05->fill(_histMeanChMultRapt05->bin(0).xMid(), rapt05);
      _histMeanChMultRapt10->fill(_histMeanChMultRapt10->bin(0).xMid(), rapt10);
      _histMeanChMultRapt15->fill(_histMeanChMultRapt15->bin(0).xMid(), rapt15);
      _histMeanChMultRapt20->fill(_histMeanChMultRapt20->bin(0).xMid(), rapt20);
      _histMeanChMult->fill(_histMeanChMult->bin(0).xMid(), numParticles);


      //// Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(e, "UFS");
      for (Particles::const_iterator p = ufs.particles().begin(); p != ufs.particles().end(); ++p) {
        const Vector3 mom3 = p->momentum().p3();
        int id = abs(p->pid());
        const double mom = mom3.mod();
        const double energy = p->momentum().E();
        const double scaledMom = mom/meanBeamMom;
        const double scaledEnergy = energy/meanBeamMom;  // meanBeamMom is approximately beam energy
        switch (id) {
           case 22:
              _histMultiPhoton->fill(-1.*std::log(scaledMom));
              break;
           case -321:
           case 321:
              _histMultiKPlus->fill(scaledMom);
              break;
           case 211:
           case -211:
              _histMultiPiPlus->fill(scaledMom);
              break;
           case 2212:
           case -2212:
              _histMultiP->fill(scaledMom);
              break;
           case 111:
              _histMultiPi0->fill(scaledMom);
              _histMeanMultiPi0->fill(_histMeanMultiPi0->bin(0).xMid());
              break;
           case 221:
              if (scaledMom >= 0.1) {
                _histMultiEta->fill(scaledEnergy);
                _histMeanMultiEta->fill(_histMeanMultiEta->bin(0).xMid());
              }
              break;
           case 331:
              if (scaledMom >= 0.1) {
                _histMultiEtaPrime->fill(scaledEnergy);
                _histMeanMultiEtaPrime->fill(_histMeanMultiEtaPrime->bin(0).xMid());
              }
              break;
           case 130: //klong
           case 310: //kshort
              _histMultiK0->fill(scaledMom);
              _histMeanMultiK0->fill(_histMeanMultiK0->bin(0).xMid());
              break;
           case 113:
              _histMultiRho->fill(scaledMom);
              _histMeanMultiRho->fill(_histMeanMultiRho->bin(0).xMid());
              break;
           case 223:
              _histMultiOmega782->fill(scaledMom);
              _histMeanMultiOmega782->fill(_histMeanMultiOmega782->bin(0).xMid());
              break;
           case 333:
              _histMultiPhi->fill(scaledMom);
              _histMeanMultiPhi->fill(_histMeanMultiPhi->bin(0).xMid());
              break;
           case 313:
           case -313:
              _histMultiKStar892_0->fill(scaledMom);
              _histMeanMultiKStar892_0->fill(_histMeanMultiKStar892_0->bin(0).xMid());
              break;
           case 323:
           case -323:
              _histMultiKStar892Plus->fill(scaledEnergy);
              _histMeanMultiKStar892Plus->fill(_histMeanMultiKStar892Plus->bin(0).xMid());
              break;
           case 3122:
           case -3122:
              _histMultiLambda0->fill(scaledMom);
              _histMeanMultiLambda0->fill(_histMeanMultiLambda0->bin(0).xMid());
              break;
           case 3212:
           case -3212:
              _histMeanMultiSigma0->fill(_histMeanMultiSigma0->bin(0).xMid());
              break;
           case 3312:
           case -3312:
              _histMultiXiMinus->fill(scaledEnergy);
              _histMeanMultiXiMinus->fill(_histMeanMultiXiMinus->bin(0).xMid());
              break;
           case 3114:
           case -3114:
           case 3224:
           case -3224:
              _histMultiSigma1385Plus->fill(scaledEnergy);
              _histMeanMultiSigma1385Plus->fill(_histMeanMultiSigma1385Plus->bin(0).xMid());
              break;
           case 3324:
           case -3324:
              _histMultiXi1530_0->fill(scaledEnergy);
              _histMeanMultiXi1530_0->fill(_histMeanMultiXi1530_0->bin(0).xMid());
              break;
           case 3334:
              _histMeanMultiOmegaOmegaBar->fill(_histMeanMultiOmegaOmegaBar->bin(0).xMid());
              break;
        }
      }

    }



    /// Finalize
    void finalize() {
      // Normalize inclusive single particle distributions to the average number
      // of charged particles per event.
      const double avgNumParts = _weightedTotalPartNum->sumW() / sumOfWeights();

      normalize(_histPtSIn, avgNumParts);
      normalize(_histPtSOut, avgNumParts);

      normalize(_histRapidityT, avgNumParts);
      normalize(_histY3);

      normalize(_histLogScaledMom, avgNumParts);
      normalize(_histScaledMom, avgNumParts);

      // particle spectra
      scale(_histMultiPiPlus        ,1./sumOfWeights());
      scale(_histMultiKPlus         ,1./sumOfWeights());
      scale(_histMultiP             ,1./sumOfWeights());
      scale(_histMultiPhoton        ,1./sumOfWeights());
      scale(_histMultiPi0           ,1./sumOfWeights());
      scale(_histMultiEta           ,1./sumOfWeights());
      scale(_histMultiEtaPrime      ,1./sumOfWeights());
      scale(_histMultiK0            ,1./sumOfWeights());
      scale(_histMultiLambda0       ,1./sumOfWeights());
      scale(_histMultiXiMinus       ,1./sumOfWeights());
      scale(_histMultiSigma1385Plus ,1./sumOfWeights());
      scale(_histMultiXi1530_0      ,1./sumOfWeights());
      scale(_histMultiRho           ,1./sumOfWeights());
      scale(_histMultiOmega782      ,1./sumOfWeights());
      scale(_histMultiKStar892_0    ,1./sumOfWeights());
      scale(_histMultiPhi           ,1./sumOfWeights());

      scale(_histMultiKStar892Plus  ,1./sumOfWeights());

      // event shape
      normalize(_hist1MinusT);
      normalize(_histTMinor);
      normalize(_histOblateness);

      normalize(_histSphericity);
      normalize(_histAplanarity);
      normalize(_histHeavyJetMass);
      normalize(_histCParam);


      // mean multiplicities
      scale(_histChMult              , 2.0/sumOfWeights()); // taking into account the binwidth of 2
      scale(_histMeanChMult          , 1.0/sumOfWeights());
      scale(_histMeanChMultRapt05    , 1.0/sumOfWeights());
      scale(_histMeanChMultRapt10    , 1.0/sumOfWeights());
      scale(_histMeanChMultRapt15    , 1.0/sumOfWeights());
      scale(_histMeanChMultRapt20    , 1.0/sumOfWeights());


      scale(_histMeanMultiPi0          , 1.0/sumOfWeights());
      scale(_histMeanMultiEta          , 1.0/sumOfWeights());
      scale(_histMeanMultiEtaPrime     , 1.0/sumOfWeights());
      scale(_histMeanMultiK0           , 1.0/sumOfWeights());
      scale(_histMeanMultiRho          , 1.0/sumOfWeights());
      scale(_histMeanMultiOmega782     , 1.0/sumOfWeights());
      scale(_histMeanMultiPhi          , 1.0/sumOfWeights());
      scale(_histMeanMultiKStar892Plus , 1.0/sumOfWeights());
      scale(_histMeanMultiKStar892_0   , 1.0/sumOfWeights());
      scale(_histMeanMultiLambda0      , 1.0/sumOfWeights());
      scale(_histMeanMultiSigma0       , 1.0/sumOfWeights());
      scale(_histMeanMultiXiMinus      , 1.0/sumOfWeights());
      scale(_histMeanMultiSigma1385Plus, 1.0/sumOfWeights());
      scale(_histMeanMultiXi1530_0     , 1.0/sumOfWeights());
      scale(_histMeanMultiOmegaOmegaBar, 1.0/sumOfWeights());
    }

    /// @}


  private:

    /// Store the weighted sums of numbers of charged / charged+neutral
    /// particles - used to calculate average number of particles for the
    /// inclusive single particle distributions' normalisations.
    CounterPtr _weightedTotalPartNum;

    /// @name Histograms
    /// @{
    Histo1DPtr _histSphericity;
    Histo1DPtr _histAplanarity;

    Histo1DPtr _hist1MinusT;
    Histo1DPtr _histTMinor;

    Histo1DPtr _histY3;
    Histo1DPtr _histHeavyJetMass;
    Histo1DPtr _histCParam;
    Histo1DPtr _histOblateness;

    Histo1DPtr _histScaledMom;
    Histo1DPtr _histRapidityT;

    Histo1DPtr _histPtSIn;
    Histo1DPtr _histPtSOut;

    Histo1DPtr _histJetRate2Durham;
    Histo1DPtr _histJetRate3Durham;
    Histo1DPtr _histJetRate4Durham;
    Histo1DPtr _histJetRate5Durham;

    Histo1DPtr _histLogScaledMom;

    Histo1DPtr _histChMult;

    Histo1DPtr _histMultiPiPlus;
    Histo1DPtr _histMultiKPlus;
    Histo1DPtr _histMultiP;
    Histo1DPtr _histMultiPhoton;
    Histo1DPtr _histMultiPi0;
    Histo1DPtr _histMultiEta;
    Histo1DPtr _histMultiEtaPrime;
    Histo1DPtr _histMultiK0;
    Histo1DPtr _histMultiLambda0;
    Histo1DPtr _histMultiXiMinus;
    Histo1DPtr _histMultiSigma1385Plus;
    Histo1DPtr _histMultiXi1530_0;
    Histo1DPtr _histMultiRho;
    Histo1DPtr _histMultiOmega782;
    Histo1DPtr _histMultiKStar892_0;
    Histo1DPtr _histMultiPhi;
    Histo1DPtr _histMultiKStar892Plus;

    // mean multiplicities
    Histo1DPtr _histMeanChMult;
    Histo1DPtr _histMeanChMultRapt05;
    Histo1DPtr _histMeanChMultRapt10;
    Histo1DPtr _histMeanChMultRapt15;
    Histo1DPtr _histMeanChMultRapt20;

    Histo1DPtr _histMeanMultiPi0;
    Histo1DPtr _histMeanMultiEta;
    Histo1DPtr _histMeanMultiEtaPrime;
    Histo1DPtr _histMeanMultiK0;
    Histo1DPtr _histMeanMultiRho;
    Histo1DPtr _histMeanMultiOmega782;
    Histo1DPtr _histMeanMultiPhi;
    Histo1DPtr _histMeanMultiKStar892Plus;
    Histo1DPtr _histMeanMultiKStar892_0;
    Histo1DPtr _histMeanMultiLambda0;
    Histo1DPtr _histMeanMultiSigma0;
    Histo1DPtr _histMeanMultiXiMinus;
    Histo1DPtr _histMeanMultiSigma1385Plus;
    Histo1DPtr _histMeanMultiXi1530_0;
    Histo1DPtr _histMeanMultiOmegaOmegaBar;
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(ALEPH_1996_S3486095, ALEPH_1996_I428072);

}
