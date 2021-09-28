// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Implementation of PDG hadron multiplicities
  /// @author Hendrik Hoeth
  class PDG_HADRON_MULTIPLICITIES : public Analysis {
  public:

    /// Constructor
    PDG_HADRON_MULTIPLICITIES() : Analysis("PDG_HADRON_MULTIPLICITIES")
    {
    }


    /// @name Analysis methods
    //@{

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

      MSG_DEBUG("sqrt(s) = " << sqrtS()/GeV << " GeV");

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableFinalState>(e, "UFS");

      if (sqrtS()/GeV >= 9.5 && sqrtS()/GeV <= 10.5) {
        for (const Particle& p : ufs.particles()) {
          const PdgId id = p.abspid();
          switch (id) {
             case 211:
                _histMeanMultiPiPlus->fill(_histMeanMultiPiPlus->bin(0).xMid());
                break;
             case 111:
                _histMeanMultiPi0->fill(_histMeanMultiPi0->bin(0).xMid());
                break;
             case 321:
                _histMeanMultiKPlus->fill(_histMeanMultiKPlus->bin(0).xMid());
                break;
             case 130:
             case 310:
                _histMeanMultiK0->fill(_histMeanMultiK0->bin(0).xMid());
                break;
             case 221:
                _histMeanMultiEta->fill(_histMeanMultiEta->bin(0).xMid());
                break;
             case 331:
                _histMeanMultiEtaPrime->fill(_histMeanMultiEtaPrime->bin(0).xMid());
                break;
             case 411:
                _histMeanMultiDPlus->fill(_histMeanMultiDPlus->bin(0).xMid());
                break;
             case 421:
                _histMeanMultiD0->fill(_histMeanMultiD0->bin(0).xMid());
                break;
             case 431:
                _histMeanMultiDPlus_s->fill(_histMeanMultiDPlus_s->bin(0).xMid());
                break;
             case 9010221:
                _histMeanMultiF0_980->fill(_histMeanMultiF0_980->bin(0).xMid());
                break;
             case 113:
                _histMeanMultiRho770_0->fill(_histMeanMultiRho770_0->bin(0).xMid());
                break;
             case 223:
                _histMeanMultiOmega782->fill(_histMeanMultiOmega782->bin(0).xMid());
                break;
             case 323:
                _histMeanMultiKStar892Plus->fill(_histMeanMultiKStar892Plus->bin(0).xMid());
                break;
             case 313:
                _histMeanMultiKStar892_0->fill(_histMeanMultiKStar892_0->bin(0).xMid());
                break;
             case 333:
                _histMeanMultiPhi1020->fill(_histMeanMultiPhi1020->bin(0).xMid());
                break;
             case 413:
                _histMeanMultiDStar2010Plus->fill(_histMeanMultiDStar2010Plus->bin(0).xMid());
                break;
             case 423:
                _histMeanMultiDStar2007_0->fill(_histMeanMultiDStar2007_0->bin(0).xMid());
                break;
             case 433:
                _histMeanMultiDStar_s2112Plus->fill(_histMeanMultiDStar_s2112Plus->bin(0).xMid());
                break;
             case 443:
                _histMeanMultiJPsi1S->fill(_histMeanMultiJPsi1S->bin(0).xMid());
                break;
             case 225:
                _histMeanMultiF2_1270->fill(_histMeanMultiF2_1270->bin(0).xMid());
                break;
             case 2212:
                _histMeanMultiP->fill(_histMeanMultiP->bin(0).xMid());
                break;
             case 3122:
                _histMeanMultiLambda->fill(_histMeanMultiLambda->bin(0).xMid());
                break;
             case 3212:
                _histMeanMultiSigma0->fill(_histMeanMultiSigma0->bin(0).xMid());
                break;
             case 3312:
                _histMeanMultiXiMinus->fill(_histMeanMultiXiMinus->bin(0).xMid());
                break;
             case 2224:
                _histMeanMultiDelta1232PlusPlus->fill(_histMeanMultiDelta1232PlusPlus->bin(0).xMid());
                break;
             case 3114:
                _histMeanMultiSigma1385Minus->fill(_histMeanMultiSigma1385Minus->bin(0).xMid());
                _histMeanMultiSigma1385PlusMinus->fill(_histMeanMultiSigma1385PlusMinus->bin(0).xMid());
                break;
             case 3224:
                _histMeanMultiSigma1385Plus->fill(_histMeanMultiSigma1385Plus->bin(0).xMid());
                _histMeanMultiSigma1385PlusMinus->fill(_histMeanMultiSigma1385PlusMinus->bin(0).xMid());
                break;
             case 3324:
                _histMeanMultiXi1530_0->fill(_histMeanMultiXi1530_0->bin(0).xMid());
                break;
             case 3334:
                _histMeanMultiOmegaMinus->fill(_histMeanMultiOmegaMinus->bin(0).xMid());
                break;
             case 4122:
                _histMeanMultiLambda_c_Plus->fill(_histMeanMultiLambda_c_Plus->bin(0).xMid());
                break;
             case 4222:
             case 4112:
                _histMeanMultiSigma_c_PlusPlus_0->fill(_histMeanMultiSigma_c_PlusPlus_0->bin(0).xMid());
                break;
             case 3124:
                _histMeanMultiLambda1520->fill(_histMeanMultiLambda1520->bin(0).xMid());
                break;
          }
        }
      }

      if (sqrtS()/GeV >= 29 && sqrtS()/GeV <= 35) {
        for (const Particle& p : ufs.particles()) {
          const PdgId id = p.abspid();
          switch (id) {
             case 211:
                _histMeanMultiPiPlus->fill(_histMeanMultiPiPlus->bin(0).xMid());
                break;
             case 111:
                _histMeanMultiPi0->fill(_histMeanMultiPi0->bin(0).xMid());
                break;
             case 321:
                _histMeanMultiKPlus->fill(_histMeanMultiKPlus->bin(0).xMid());
                break;
             case 130:
             case 310:
                _histMeanMultiK0->fill(_histMeanMultiK0->bin(0).xMid());
                break;
             case 221:
                _histMeanMultiEta->fill(_histMeanMultiEta->bin(0).xMid());
                break;
             case 331:
                _histMeanMultiEtaPrime->fill(_histMeanMultiEtaPrime->bin(0).xMid());
                break;
             case 411:
                _histMeanMultiDPlus->fill(_histMeanMultiDPlus->bin(0).xMid());
                break;
             case 421:
                _histMeanMultiD0->fill(_histMeanMultiD0->bin(0).xMid());
                break;
             case 431:
                _histMeanMultiDPlus_s->fill(_histMeanMultiDPlus_s->bin(0).xMid());
                break;
             case 9010221:
                _histMeanMultiF0_980->fill(_histMeanMultiF0_980->bin(0).xMid());
                break;
             case 113:
                _histMeanMultiRho770_0->fill(_histMeanMultiRho770_0->bin(0).xMid());
                break;
             case 323:
                _histMeanMultiKStar892Plus->fill(_histMeanMultiKStar892Plus->bin(0).xMid());
                break;
             case 313:
                _histMeanMultiKStar892_0->fill(_histMeanMultiKStar892_0->bin(0).xMid());
                break;
             case 333:
                _histMeanMultiPhi1020->fill(_histMeanMultiPhi1020->bin(0).xMid());
                break;
             case 413:
                _histMeanMultiDStar2010Plus->fill(_histMeanMultiDStar2010Plus->bin(0).xMid());
                break;
             case 423:
                _histMeanMultiDStar2007_0->fill(_histMeanMultiDStar2007_0->bin(0).xMid());
                break;
             case 225:
                _histMeanMultiF2_1270->fill(_histMeanMultiF2_1270->bin(0).xMid());
                break;
             case 325:
                _histMeanMultiK2Star1430Plus->fill(_histMeanMultiK2Star1430Plus->bin(0).xMid());
                break;
             case 315:
                _histMeanMultiK2Star1430_0->fill(_histMeanMultiK2Star1430_0->bin(0).xMid());
                break;
             case 2212:
                _histMeanMultiP->fill(_histMeanMultiP->bin(0).xMid());
                break;
             case 3122:
                _histMeanMultiLambda->fill(_histMeanMultiLambda->bin(0).xMid());
                break;
             case 3312:
                _histMeanMultiXiMinus->fill(_histMeanMultiXiMinus->bin(0).xMid());
                break;
             case 3114:
                _histMeanMultiSigma1385Minus->fill(_histMeanMultiSigma1385Minus->bin(0).xMid());
                _histMeanMultiSigma1385PlusMinus->fill(_histMeanMultiSigma1385PlusMinus->bin(0).xMid());
                break;
             case 3224:
                _histMeanMultiSigma1385Plus->fill(_histMeanMultiSigma1385Plus->bin(0).xMid());
                _histMeanMultiSigma1385PlusMinus->fill(_histMeanMultiSigma1385PlusMinus->bin(0).xMid());
                break;
             case 3334:
                _histMeanMultiOmegaMinus->fill(_histMeanMultiOmegaMinus->bin(0).xMid());
                break;
             case 4122:
                _histMeanMultiLambda_c_Plus->fill(_histMeanMultiLambda_c_Plus->bin(0).xMid());
                break;
          }
        }
      }

      if (sqrtS()/GeV >= 89.5 && sqrtS()/GeV <= 91.8) {
        for (const Particle& p : ufs.particles()) {
          const PdgId id = p.abspid();
          switch (id) {
             case 211:
                _histMeanMultiPiPlus->fill(_histMeanMultiPiPlus->bin(0).xMid());
                break;
             case 111:
                _histMeanMultiPi0->fill(_histMeanMultiPi0->bin(0).xMid());
                break;
             case 321:
                _histMeanMultiKPlus->fill(_histMeanMultiKPlus->bin(0).xMid());
                break;
             case 130:
             case 310:
                _histMeanMultiK0->fill(_histMeanMultiK0->bin(0).xMid());
                break;
             case 221:
                _histMeanMultiEta->fill(_histMeanMultiEta->bin(0).xMid());
                break;
             case 331:
                _histMeanMultiEtaPrime->fill(_histMeanMultiEtaPrime->bin(0).xMid());
                break;
             case 411:
                _histMeanMultiDPlus->fill(_histMeanMultiDPlus->bin(0).xMid());
                break;
             case 421:
                _histMeanMultiD0->fill(_histMeanMultiD0->bin(0).xMid());
                break;
             case 431:
                _histMeanMultiDPlus_s->fill(_histMeanMultiDPlus_s->bin(0).xMid());
                break;
             case 511:
                _histMeanMultiBPlus_B0_d->fill(_histMeanMultiBPlus_B0_d->bin(0).xMid());
                break;
             case 521:
                _histMeanMultiBPlus_B0_d->fill(_histMeanMultiBPlus_B0_d->bin(0).xMid());
                _histMeanMultiBPlus_u->fill(_histMeanMultiBPlus_u->bin(0).xMid());
                break;
             case 531:
                _histMeanMultiB0_s->fill(_histMeanMultiB0_s->bin(0).xMid());
                break;
             case 9010221:
                _histMeanMultiF0_980->fill(_histMeanMultiF0_980->bin(0).xMid());
                break;
             case 9000211:
                _histMeanMultiA0_980Plus->fill(_histMeanMultiA0_980Plus->bin(0).xMid());
                break;
             case 113:
                _histMeanMultiRho770_0->fill(_histMeanMultiRho770_0->bin(0).xMid());
                break;
             case 213:
                _histMeanMultiRho770Plus->fill(_histMeanMultiRho770Plus->bin(0).xMid());
                break;
             case 223:
                _histMeanMultiOmega782->fill(_histMeanMultiOmega782->bin(0).xMid());
                break;
             case 323:
                _histMeanMultiKStar892Plus->fill(_histMeanMultiKStar892Plus->bin(0).xMid());
                break;
             case 313:
                _histMeanMultiKStar892_0->fill(_histMeanMultiKStar892_0->bin(0).xMid());
                break;
             case 333:
                _histMeanMultiPhi1020->fill(_histMeanMultiPhi1020->bin(0).xMid());
                break;
             case 413:
                _histMeanMultiDStar2010Plus->fill(_histMeanMultiDStar2010Plus->bin(0).xMid());
                break;
             case 433:
                _histMeanMultiDStar_s2112Plus->fill(_histMeanMultiDStar_s2112Plus->bin(0).xMid());
                break;
             case 513:
             case 523:
             case 533:
                _histMeanMultiBStar->fill(_histMeanMultiBStar->bin(0).xMid());
                break;
             case 443:
                _histMeanMultiJPsi1S->fill(_histMeanMultiJPsi1S->bin(0).xMid());
                break;
             case 100443:
                _histMeanMultiPsi2S->fill(_histMeanMultiPsi2S->bin(0).xMid());
                break;
             case 553:
                _histMeanMultiUpsilon1S->fill(_histMeanMultiUpsilon1S->bin(0).xMid());
                break;
             case 20223:
                _histMeanMultiF1_1285->fill(_histMeanMultiF1_1285->bin(0).xMid());
                break;
             case 20333:
                _histMeanMultiF1_1420->fill(_histMeanMultiF1_1420->bin(0).xMid());
                break;
             case 445:
                _histMeanMultiChi_c1_3510->fill(_histMeanMultiChi_c1_3510->bin(0).xMid());
                break;
             case 225:
                _histMeanMultiF2_1270->fill(_histMeanMultiF2_1270->bin(0).xMid());
                break;
             case 335:
                _histMeanMultiF2Prime1525->fill(_histMeanMultiF2Prime1525->bin(0).xMid());
                break;
             case 315:
                _histMeanMultiK2Star1430_0->fill(_histMeanMultiK2Star1430_0->bin(0).xMid());
                break;
             case 515:
             case 525:
             case 535:
                _histMeanMultiBStarStar->fill(_histMeanMultiBStarStar->bin(0).xMid());
                break;
             case 10433:
             case 20433:
                _histMeanMultiDs1Plus->fill(_histMeanMultiDs1Plus->bin(0).xMid());
                break;
             case 435:
                _histMeanMultiDs2Plus->fill(_histMeanMultiDs2Plus->bin(0).xMid());
                break;
             case 2212:
                _histMeanMultiP->fill(_histMeanMultiP->bin(0).xMid());
                break;
             case 3122:
                _histMeanMultiLambda->fill(_histMeanMultiLambda->bin(0).xMid());
                break;
             case 3212:
                _histMeanMultiSigma0->fill(_histMeanMultiSigma0->bin(0).xMid());
                break;
             case 3112:
                _histMeanMultiSigmaMinus->fill(_histMeanMultiSigmaMinus->bin(0).xMid());
                _histMeanMultiSigmaPlusMinus->fill(_histMeanMultiSigmaPlusMinus->bin(0).xMid());
                break;
             case 3222:
                _histMeanMultiSigmaPlus->fill(_histMeanMultiSigmaPlus->bin(0).xMid());
                _histMeanMultiSigmaPlusMinus->fill(_histMeanMultiSigmaPlusMinus->bin(0).xMid());
                break;
             case 3312:
                _histMeanMultiXiMinus->fill(_histMeanMultiXiMinus->bin(0).xMid());
                break;
             case 2224:
                _histMeanMultiDelta1232PlusPlus->fill(_histMeanMultiDelta1232PlusPlus->bin(0).xMid());
                break;
             case 3114:
                _histMeanMultiSigma1385Minus->fill(_histMeanMultiSigma1385Minus->bin(0).xMid());
                _histMeanMultiSigma1385PlusMinus->fill(_histMeanMultiSigma1385PlusMinus->bin(0).xMid());
                break;
             case 3224:
                _histMeanMultiSigma1385Plus->fill(_histMeanMultiSigma1385Plus->bin(0).xMid());
                _histMeanMultiSigma1385PlusMinus->fill(_histMeanMultiSigma1385PlusMinus->bin(0).xMid());
                break;
             case 3324:
                _histMeanMultiXi1530_0->fill(_histMeanMultiXi1530_0->bin(0).xMid());
                break;
             case 3334:
                _histMeanMultiOmegaMinus->fill(_histMeanMultiOmegaMinus->bin(0).xMid());
                break;
             case 4122:
                _histMeanMultiLambda_c_Plus->fill(_histMeanMultiLambda_c_Plus->bin(0).xMid());
                break;
             case 5122:
                _histMeanMultiLambda_b_0->fill(_histMeanMultiLambda_b_0->bin(0).xMid());
                break;
             case 3124:
                _histMeanMultiLambda1520->fill(_histMeanMultiLambda1520->bin(0).xMid());
                break;
          }
        }
      }

      if (sqrtS()/GeV >= 130 && sqrtS()/GeV <= 200) {
        for (const Particle& p : ufs.particles()) {
          const PdgId id = p.abspid();
          switch (id) {
             case 211:
                _histMeanMultiPiPlus->fill(_histMeanMultiPiPlus->bin(0).xMid());
                break;
             case 321:
                _histMeanMultiKPlus->fill(_histMeanMultiKPlus->bin(0).xMid());
                break;
             case 130:
             case 310:
                _histMeanMultiK0->fill(_histMeanMultiK0->bin(0).xMid());
                break;
             case 2212:
                _histMeanMultiP->fill(_histMeanMultiP->bin(0).xMid());
                break;
             case 3122:
                _histMeanMultiLambda->fill(_histMeanMultiLambda->bin(0).xMid());
                break;
          }
        }
      }

    }



    void init() {
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      if (sqrtS()/GeV >= 9.5 && sqrtS()/GeV <= 10.5) {
        book(_histMeanMultiPiPlus             , 1, 1, 1);
        book(_histMeanMultiPi0                , 2, 1, 1);
        book(_histMeanMultiKPlus              , 3, 1, 1);
        book(_histMeanMultiK0                 , 4, 1, 1);
        book(_histMeanMultiEta                , 5, 1, 1);
        book(_histMeanMultiEtaPrime           , 6, 1, 1);
        book(_histMeanMultiDPlus              , 7, 1, 1);
        book(_histMeanMultiD0                 , 8, 1, 1);
        book(_histMeanMultiDPlus_s            , 9, 1, 1);
        book(_histMeanMultiF0_980             ,13, 1, 1);
        book(_histMeanMultiRho770_0           ,15, 1, 1);
        book(_histMeanMultiOmega782           ,17, 1, 1);
        book(_histMeanMultiKStar892Plus       ,18, 1, 1);
        book(_histMeanMultiKStar892_0         ,19, 1, 1);
        book(_histMeanMultiPhi1020            ,20, 1, 1);
        book(_histMeanMultiDStar2010Plus      ,21, 1, 1);
        book(_histMeanMultiDStar2007_0        ,22, 1, 1);
        book(_histMeanMultiDStar_s2112Plus    ,23, 1, 1);
        book(_histMeanMultiJPsi1S             ,25, 1, 1);
        book(_histMeanMultiF2_1270            ,31, 1, 1);
        book(_histMeanMultiP                  ,38, 1, 1);
        book(_histMeanMultiLambda             ,39, 1, 1);
        book(_histMeanMultiSigma0             ,40, 1, 1);
        book(_histMeanMultiXiMinus            ,44, 1, 1);
        book(_histMeanMultiDelta1232PlusPlus  ,45, 1, 1);
        book(_histMeanMultiSigma1385Minus     ,46, 1, 1);
        book(_histMeanMultiSigma1385Plus      ,47, 1, 1);
        book(_histMeanMultiSigma1385PlusMinus ,48, 1, 1);
        book(_histMeanMultiXi1530_0           ,49, 1, 1);
        book(_histMeanMultiOmegaMinus         ,50, 1, 1);
        book(_histMeanMultiLambda_c_Plus      ,51, 1, 1);
        book(_histMeanMultiSigma_c_PlusPlus_0 ,53, 1, 1);
        book(_histMeanMultiLambda1520         ,54, 1, 1);
      }

      if (sqrtS()/GeV >= 29 && sqrtS()/GeV <= 35) {
        book(_histMeanMultiPiPlus             , 1, 1, 2);
        book(_histMeanMultiPi0                , 2, 1, 2);
        book(_histMeanMultiKPlus              , 3, 1, 2);
        book(_histMeanMultiK0                 , 4, 1, 2);
        book(_histMeanMultiEta                , 5, 1, 2);
        book(_histMeanMultiEtaPrime           , 6, 1, 2);
        book(_histMeanMultiDPlus              , 7, 1, 2);
        book(_histMeanMultiD0                 , 8, 1, 2);
        book(_histMeanMultiDPlus_s            , 9, 1, 2);
        book(_histMeanMultiF0_980             ,13, 1, 2);
        book(_histMeanMultiRho770_0           ,15, 1, 2);
        book(_histMeanMultiKStar892Plus       ,18, 1, 2);
        book(_histMeanMultiKStar892_0         ,19, 1, 2);
        book(_histMeanMultiPhi1020            ,20, 1, 2);
        book(_histMeanMultiDStar2010Plus      ,21, 1, 2);
        book(_histMeanMultiDStar2007_0        ,22, 1, 2);
        book(_histMeanMultiF2_1270            ,31, 1, 2);
        book(_histMeanMultiK2Star1430Plus     ,33, 1, 1);
        book(_histMeanMultiK2Star1430_0       ,34, 1, 1);
        book(_histMeanMultiP                  ,38, 1, 2);
        book(_histMeanMultiLambda             ,39, 1, 2);
        book(_histMeanMultiXiMinus            ,44, 1, 2);
        book(_histMeanMultiSigma1385Minus     ,46, 1, 2);
        book(_histMeanMultiSigma1385Plus      ,47, 1, 2);
        book(_histMeanMultiSigma1385PlusMinus ,48, 1, 2);
        book(_histMeanMultiOmegaMinus         ,50, 1, 2);
        book(_histMeanMultiLambda_c_Plus      ,51, 1, 2);
      }

      if (sqrtS()/GeV >= 89.5 && sqrtS()/GeV <= 91.8) {
        book(_histMeanMultiPiPlus             , 1, 1, 3);
        book(_histMeanMultiPi0                , 2, 1, 3);
        book(_histMeanMultiKPlus              , 3, 1, 3);
        book(_histMeanMultiK0                 , 4, 1, 3);
        book(_histMeanMultiEta                , 5, 1, 3);
        book(_histMeanMultiEtaPrime           , 6, 1, 3);
        book(_histMeanMultiDPlus              , 7, 1, 3);
        book(_histMeanMultiD0                 , 8, 1, 3);
        book(_histMeanMultiDPlus_s            , 9, 1, 3);
        book(_histMeanMultiBPlus_B0_d         ,10, 1, 1);
        book(_histMeanMultiBPlus_u            ,11, 1, 1);
        book(_histMeanMultiB0_s               ,12, 1, 1);
        book(_histMeanMultiF0_980             ,13, 1, 3);
        book(_histMeanMultiA0_980Plus         ,14, 1, 1);
        book(_histMeanMultiRho770_0           ,15, 1, 3);
        book(_histMeanMultiRho770Plus         ,16, 1, 1);
        book(_histMeanMultiOmega782           ,17, 1, 2);
        book(_histMeanMultiKStar892Plus       ,18, 1, 3);
        book(_histMeanMultiKStar892_0         ,19, 1, 3);
        book(_histMeanMultiPhi1020            ,20, 1, 3);
        book(_histMeanMultiDStar2010Plus      ,21, 1, 3);
        book(_histMeanMultiDStar_s2112Plus    ,23, 1, 2);
        book(_histMeanMultiBStar              ,24, 1, 1);
        book(_histMeanMultiJPsi1S             ,25, 1, 2);
        book(_histMeanMultiPsi2S              ,26, 1, 1);
        book(_histMeanMultiUpsilon1S          ,27, 1, 1);
        book(_histMeanMultiF1_1285            ,28, 1, 1);
        book(_histMeanMultiF1_1420            ,29, 1, 1);
        book(_histMeanMultiChi_c1_3510        ,30, 1, 1);
        book(_histMeanMultiF2_1270            ,31, 1, 3);
        book(_histMeanMultiF2Prime1525        ,32, 1, 1);
        book(_histMeanMultiK2Star1430_0       ,34, 1, 2);
        book(_histMeanMultiBStarStar          ,35, 1, 1);
        book(_histMeanMultiDs1Plus            ,36, 1, 1);
        book(_histMeanMultiDs2Plus            ,37, 1, 1);
        book(_histMeanMultiP                  ,38, 1, 3);
        book(_histMeanMultiLambda             ,39, 1, 3);
        book(_histMeanMultiSigma0             ,40, 1, 2);
        book(_histMeanMultiSigmaMinus         ,41, 1, 1);
        book(_histMeanMultiSigmaPlus          ,42, 1, 1);
        book(_histMeanMultiSigmaPlusMinus     ,43, 1, 1);
        book(_histMeanMultiXiMinus            ,44, 1, 3);
        book(_histMeanMultiDelta1232PlusPlus  ,45, 1, 2);
        book(_histMeanMultiSigma1385Minus     ,46, 1, 3);
        book(_histMeanMultiSigma1385Plus      ,47, 1, 3);
        book(_histMeanMultiSigma1385PlusMinus ,48, 1, 3);
        book(_histMeanMultiXi1530_0           ,49, 1, 2);
        book(_histMeanMultiOmegaMinus         ,50, 1, 3);
        book(_histMeanMultiLambda_c_Plus      ,51, 1, 3);
        book(_histMeanMultiLambda_b_0         ,52, 1, 1);
        book(_histMeanMultiLambda1520         ,54, 1, 2);
      }

      if (sqrtS()/GeV >= 130 && sqrtS()/GeV <= 200) {
        book(_histMeanMultiPiPlus            , 1, 1, 4);
        book(_histMeanMultiKPlus             , 3, 1, 4);
        book(_histMeanMultiK0                , 4, 1, 4);
        book(_histMeanMultiP                 ,38, 1, 4);
        book(_histMeanMultiLambda            ,39, 1, 4);
      }
    }



    // Finalize
    void finalize() {
      if (sqrtS()/GeV >= 9.5 && sqrtS()/GeV <= 10.5) {
        scale(_histMeanMultiPiPlus            , 1.0/sumOfWeights());
        scale(_histMeanMultiPi0               , 1.0/sumOfWeights());
        scale(_histMeanMultiKPlus             , 1.0/sumOfWeights());
        scale(_histMeanMultiK0                , 1.0/sumOfWeights());
        scale(_histMeanMultiEta               , 1.0/sumOfWeights());
        scale(_histMeanMultiEtaPrime          , 1.0/sumOfWeights());
        scale(_histMeanMultiDPlus             , 1.0/sumOfWeights());
        scale(_histMeanMultiD0                , 1.0/sumOfWeights());
        scale(_histMeanMultiDPlus_s           , 1.0/sumOfWeights());
        scale(_histMeanMultiF0_980            , 1.0/sumOfWeights());
        scale(_histMeanMultiRho770_0          , 1.0/sumOfWeights());
        scale(_histMeanMultiOmega782          , 1.0/sumOfWeights());
        scale(_histMeanMultiKStar892Plus      , 1.0/sumOfWeights());
        scale(_histMeanMultiKStar892_0        , 1.0/sumOfWeights());
        scale(_histMeanMultiPhi1020           , 1.0/sumOfWeights());
        scale(_histMeanMultiDStar2010Plus     , 1.0/sumOfWeights());
        scale(_histMeanMultiDStar2007_0       , 1.0/sumOfWeights());
        scale(_histMeanMultiDStar_s2112Plus   , 1.0/sumOfWeights());
        scale(_histMeanMultiJPsi1S            , 1.0/sumOfWeights());
        scale(_histMeanMultiF2_1270           , 1.0/sumOfWeights());
        scale(_histMeanMultiP                 , 1.0/sumOfWeights());
        scale(_histMeanMultiLambda            , 1.0/sumOfWeights());
        scale(_histMeanMultiSigma0            , 1.0/sumOfWeights());
        scale(_histMeanMultiXiMinus           , 1.0/sumOfWeights());
        scale(_histMeanMultiDelta1232PlusPlus , 1.0/sumOfWeights());
        scale(_histMeanMultiSigma1385Minus    , 1.0/sumOfWeights());
        scale(_histMeanMultiSigma1385Plus     , 1.0/sumOfWeights());
        scale(_histMeanMultiSigma1385PlusMinus, 1.0/sumOfWeights());
        scale(_histMeanMultiXi1530_0          , 1.0/sumOfWeights());
        scale(_histMeanMultiOmegaMinus        , 1.0/sumOfWeights());
        scale(_histMeanMultiLambda_c_Plus     , 1.0/sumOfWeights());
        scale(_histMeanMultiSigma_c_PlusPlus_0, 1.0/sumOfWeights());
        scale(_histMeanMultiLambda1520        , 1.0/sumOfWeights());
      }

      if (sqrtS()/GeV >= 29 && sqrtS()/GeV <= 35) {
        scale(_histMeanMultiPiPlus            , 5.0/sumOfWeights());
        scale(_histMeanMultiPi0               , 5.0/sumOfWeights());
        scale(_histMeanMultiKPlus             , 5.0/sumOfWeights());
        scale(_histMeanMultiK0                , 5.0/sumOfWeights());
        scale(_histMeanMultiEta               , 5.0/sumOfWeights());
        scale(_histMeanMultiEtaPrime          , 5.0/sumOfWeights());
        scale(_histMeanMultiDPlus             , 5.0/sumOfWeights());
        scale(_histMeanMultiD0                , 5.0/sumOfWeights());
        scale(_histMeanMultiDPlus_s           , 5.0/sumOfWeights());
        scale(_histMeanMultiF0_980            , 5.0/sumOfWeights());
        scale(_histMeanMultiRho770_0          , 5.0/sumOfWeights());
        scale(_histMeanMultiKStar892Plus      , 5.0/sumOfWeights());
        scale(_histMeanMultiKStar892_0        , 5.0/sumOfWeights());
        scale(_histMeanMultiPhi1020           , 5.0/sumOfWeights());
        scale(_histMeanMultiDStar2010Plus     , 5.0/sumOfWeights());
        scale(_histMeanMultiDStar2007_0       , 5.0/sumOfWeights());
        scale(_histMeanMultiF2_1270           , 5.0/sumOfWeights());
        scale(_histMeanMultiK2Star1430Plus    , 5.0/sumOfWeights());
        scale(_histMeanMultiK2Star1430_0      , 5.0/sumOfWeights());
        scale(_histMeanMultiP                 , 5.0/sumOfWeights());
        scale(_histMeanMultiLambda            , 5.0/sumOfWeights());
        scale(_histMeanMultiXiMinus           , 5.0/sumOfWeights());
        scale(_histMeanMultiSigma1385Minus    , 5.0/sumOfWeights());
        scale(_histMeanMultiSigma1385Plus     , 5.0/sumOfWeights());
        scale(_histMeanMultiSigma1385PlusMinus, 5.0/sumOfWeights());
        scale(_histMeanMultiOmegaMinus        , 5.0/sumOfWeights());
        scale(_histMeanMultiLambda_c_Plus     , 5.0/sumOfWeights());
      }

      if (sqrtS()/GeV >= 89.5 && sqrtS()/GeV <= 91.8) {
        scale(_histMeanMultiPiPlus            , 1.0/sumOfWeights());
        scale(_histMeanMultiPi0               , 1.0/sumOfWeights());
        scale(_histMeanMultiKPlus             , 1.0/sumOfWeights());
        scale(_histMeanMultiK0                , 1.0/sumOfWeights());
        scale(_histMeanMultiEta               , 1.0/sumOfWeights());
        scale(_histMeanMultiEtaPrime          , 1.0/sumOfWeights());
        scale(_histMeanMultiDPlus             , 1.0/sumOfWeights());
        scale(_histMeanMultiD0                , 1.0/sumOfWeights());
        scale(_histMeanMultiDPlus_s           , 1.0/sumOfWeights());
        scale(_histMeanMultiBPlus_B0_d        , 1.0/sumOfWeights());
        scale(_histMeanMultiBPlus_u           , 1.0/sumOfWeights());
        scale(_histMeanMultiB0_s              , 1.0/sumOfWeights());
        scale(_histMeanMultiF0_980            , 1.0/sumOfWeights());
        scale(_histMeanMultiA0_980Plus        , 1.0/sumOfWeights());
        scale(_histMeanMultiRho770_0          , 1.0/sumOfWeights());
        scale(_histMeanMultiRho770Plus        , 1.0/sumOfWeights());
        scale(_histMeanMultiOmega782          , 1.0/sumOfWeights());
        scale(_histMeanMultiKStar892Plus      , 1.0/sumOfWeights());
        scale(_histMeanMultiKStar892_0        , 1.0/sumOfWeights());
        scale(_histMeanMultiPhi1020           , 1.0/sumOfWeights());
        scale(_histMeanMultiDStar2010Plus     , 1.0/sumOfWeights());
        scale(_histMeanMultiDStar_s2112Plus   , 1.0/sumOfWeights());
        scale(_histMeanMultiBStar             , 1.0/sumOfWeights());
        scale(_histMeanMultiJPsi1S            , 1.0/sumOfWeights());
        scale(_histMeanMultiPsi2S             , 1.0/sumOfWeights());
        scale(_histMeanMultiUpsilon1S         , 1.0/sumOfWeights());
        scale(_histMeanMultiF1_1285           , 1.0/sumOfWeights());
        scale(_histMeanMultiF1_1420           , 1.0/sumOfWeights());
        scale(_histMeanMultiChi_c1_3510       , 1.0/sumOfWeights());
        scale(_histMeanMultiF2_1270           , 1.0/sumOfWeights());
        scale(_histMeanMultiF2Prime1525       , 1.0/sumOfWeights());
        scale(_histMeanMultiK2Star1430_0      , 1.0/sumOfWeights());
        scale(_histMeanMultiBStarStar         , 1.0/sumOfWeights());
        scale(_histMeanMultiDs1Plus           , 1.0/sumOfWeights());
        scale(_histMeanMultiDs2Plus           , 1.0/sumOfWeights());
        scale(_histMeanMultiP                 , 1.0/sumOfWeights());
        scale(_histMeanMultiLambda            , 1.0/sumOfWeights());
        scale(_histMeanMultiSigma0            , 1.0/sumOfWeights());
        scale(_histMeanMultiSigmaMinus        , 1.0/sumOfWeights());
        scale(_histMeanMultiSigmaPlus         , 1.0/sumOfWeights());
        scale(_histMeanMultiSigmaPlusMinus    , 1.0/sumOfWeights());
        scale(_histMeanMultiXiMinus           , 1.0/sumOfWeights());
        scale(_histMeanMultiDelta1232PlusPlus , 1.0/sumOfWeights());
        scale(_histMeanMultiSigma1385Minus    , 1.0/sumOfWeights());
        scale(_histMeanMultiSigma1385Plus     , 1.0/sumOfWeights());
        scale(_histMeanMultiSigma1385PlusMinus, 1.0/sumOfWeights());
        scale(_histMeanMultiXi1530_0          , 1.0/sumOfWeights());
        scale(_histMeanMultiOmegaMinus        , 1.0/sumOfWeights());
        scale(_histMeanMultiLambda_c_Plus     , 1.0/sumOfWeights());
        scale(_histMeanMultiLambda_b_0        , 1.0/sumOfWeights());
        scale(_histMeanMultiLambda1520        , 1.0/sumOfWeights());
      }

      if (sqrtS()/GeV >= 130 && sqrtS()/GeV <= 200) {
        scale(_histMeanMultiPiPlus           , 70.0/sumOfWeights());
        scale(_histMeanMultiKPlus            , 70.0/sumOfWeights());
        scale(_histMeanMultiK0               , 70.0/sumOfWeights());
        scale(_histMeanMultiP                , 70.0/sumOfWeights());
        scale(_histMeanMultiLambda           , 70.0/sumOfWeights());
      }
    }

    //@}


  private:

    Histo1DPtr _histMeanMultiPiPlus;
    Histo1DPtr _histMeanMultiPi0;
    Histo1DPtr _histMeanMultiKPlus;
    Histo1DPtr _histMeanMultiK0;
    Histo1DPtr _histMeanMultiEta;
    Histo1DPtr _histMeanMultiEtaPrime;
    Histo1DPtr _histMeanMultiDPlus;
    Histo1DPtr _histMeanMultiD0;
    Histo1DPtr _histMeanMultiDPlus_s;
    Histo1DPtr _histMeanMultiBPlus_B0_d;
    Histo1DPtr _histMeanMultiBPlus_u;
    Histo1DPtr _histMeanMultiB0_s;
    Histo1DPtr _histMeanMultiF0_980;
    Histo1DPtr _histMeanMultiA0_980Plus;
    Histo1DPtr _histMeanMultiRho770_0;
    Histo1DPtr _histMeanMultiRho770Plus;
    Histo1DPtr _histMeanMultiOmega782;
    Histo1DPtr _histMeanMultiKStar892Plus;
    Histo1DPtr _histMeanMultiKStar892_0;
    Histo1DPtr _histMeanMultiPhi1020;
    Histo1DPtr _histMeanMultiDStar2010Plus;
    Histo1DPtr _histMeanMultiDStar2007_0;
    Histo1DPtr _histMeanMultiDStar_s2112Plus;
    Histo1DPtr _histMeanMultiBStar;
    Histo1DPtr _histMeanMultiJPsi1S;
    Histo1DPtr _histMeanMultiPsi2S;
    Histo1DPtr _histMeanMultiUpsilon1S;
    Histo1DPtr _histMeanMultiF1_1285;
    Histo1DPtr _histMeanMultiF1_1420;
    Histo1DPtr _histMeanMultiChi_c1_3510;
    Histo1DPtr _histMeanMultiF2_1270;
    Histo1DPtr _histMeanMultiF2Prime1525;
    Histo1DPtr _histMeanMultiK2Star1430Plus;
    Histo1DPtr _histMeanMultiK2Star1430_0;
    Histo1DPtr _histMeanMultiBStarStar;
    Histo1DPtr _histMeanMultiDs1Plus;
    Histo1DPtr _histMeanMultiDs2Plus;
    Histo1DPtr _histMeanMultiP;
    Histo1DPtr _histMeanMultiLambda;
    Histo1DPtr _histMeanMultiSigma0;
    Histo1DPtr _histMeanMultiSigmaMinus;
    Histo1DPtr _histMeanMultiSigmaPlus;
    Histo1DPtr _histMeanMultiSigmaPlusMinus;
    Histo1DPtr _histMeanMultiXiMinus;
    Histo1DPtr _histMeanMultiDelta1232PlusPlus;
    Histo1DPtr _histMeanMultiSigma1385Minus;
    Histo1DPtr _histMeanMultiSigma1385Plus;
    Histo1DPtr _histMeanMultiSigma1385PlusMinus;
    Histo1DPtr _histMeanMultiXi1530_0;
    Histo1DPtr _histMeanMultiOmegaMinus;
    Histo1DPtr _histMeanMultiLambda_c_Plus;
    Histo1DPtr _histMeanMultiLambda_b_0;
    Histo1DPtr _histMeanMultiSigma_c_PlusPlus_0;
    Histo1DPtr _histMeanMultiLambda1520;

    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(PDG_HADRON_MULTIPLICITIES);

}
