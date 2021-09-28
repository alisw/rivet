// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Thrust.hh"

#define I_KNOW_THE_INITIAL_QUARKS_PROJECTION_IS_DODGY_BUT_NEED_TO_USE_IT
#include "Rivet/Projections/InitialQuarks.hh"

namespace Rivet {


  /// @brief SLD flavour-dependent fragmentation paper
  /// @author Peter Richardson
  class SLD_1999_S3743934 : public Analysis {
  public:

    /// Constructor
    SLD_1999_S3743934()
      : Analysis("SLD_1999_S3743934"),
         _multPiPlus(4),_multKPlus(4),_multK0(4),
         _multKStar0(4),_multPhi(4),
         _multProton(4),_multLambda(4)
    {    }


    /// @name Analysis methods
    //@{

    void analyze(const Event& e) {
      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = apply<FinalState>(e, "FS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed ncharged cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed ncharged cut");

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(e, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      int flavour = 0;
      const InitialQuarks& iqf = apply<InitialQuarks>(e, "IQF");

      // If we only have two quarks (qqbar), just take the flavour.
      // If we have more than two quarks, look for the highest energetic q-qbar pair.
      /// @todo Can we make this based on hadron flavour instead?
      Particles quarks;
      if (iqf.particles().size() == 2) {
        flavour = iqf.particles().front().abspid();
        quarks = iqf.particles();
      } else {
        map<int, Particle > quarkmap;
        for (const Particle& p : iqf.particles()) {
          if (quarkmap.find(p.pid()) == quarkmap.end()) quarkmap[p.pid()] = p;
          else if (quarkmap[p.pid()].E() < p.E()) quarkmap[p.pid()] = p;
        }
        double maxenergy = 0.;
        for (int i = 1; i <= 5; ++i) {
          double energy(0.);
          if (quarkmap.find( i) != quarkmap.end())
            energy += quarkmap[ i].E();
          if (quarkmap.find(-i) != quarkmap.end())
            energy += quarkmap[-i].E();
          if (energy > maxenergy)
            flavour = i;
        }
        if (quarkmap.find(flavour) != quarkmap.end())
          quarks.push_back(quarkmap[flavour]);
        if (quarkmap.find(-flavour) != quarkmap.end())
          quarks.push_back(quarkmap[-flavour]);
      }
      switch (flavour) {
      case PID::DQUARK:
      case PID::UQUARK:
      case PID::SQUARK:
        _SumOfudsWeights->fill();
        break;
      case PID::CQUARK:
        _SumOfcWeights->fill();
        break;
      case PID::BQUARK:
        _SumOfbWeights->fill();
        break;
      }
      // thrust axis for projections
      Vector3 axis = apply<Thrust>(e, "Thrust").thrustAxis();
      double dot(0.);
      if (!quarks.empty()) {
        dot = quarks[0].p3().dot(axis);
        if (quarks[0].pid() < 0) dot *= -1;
      }

      for (const Particle& p : fs.particles()) {
        const double xp = p.p3().mod()/meanBeamMom;
        // if in quark or antiquark hemisphere
        bool quark = p.p3().dot(axis)*dot > 0.;
        _h_XpChargedN->fill(xp);
        _temp_XpChargedN1->fill(xp);
        _temp_XpChargedN2->fill(xp);
        _temp_XpChargedN3->fill(xp);
        int id = p.abspid();
        // charged pions
        if (id == PID::PIPLUS) {
          _h_XpPiPlusN->fill(xp);
          _multPiPlus[0]->fill();
          switch (flavour) {
          case PID::DQUARK:
          case PID::UQUARK:
          case PID::SQUARK:
            _multPiPlus[1]->fill();
            _h_XpPiPlusLight->fill(xp);
            if( ( quark && p.pid()>0 ) || ( !quark && p.pid()<0 ))
              _h_RPiPlus->fill(xp);
            else
              _h_RPiMinus->fill(xp);
            break;
          case PID::CQUARK:
            _multPiPlus[2]->fill();
            _h_XpPiPlusCharm->fill(xp);
            break;
          case PID::BQUARK:
            _multPiPlus[3]->fill();
            _h_XpPiPlusBottom->fill(xp);
            break;
          }
        }
        else if (id == PID::KPLUS) {
          _h_XpKPlusN->fill(xp);
          _multKPlus[0]->fill();
          switch (flavour) {
          case PID::DQUARK:
          case PID::UQUARK:
          case PID::SQUARK:
            _multKPlus[1]->fill();
            _temp_XpKPlusLight->fill(xp);
            _h_XpKPlusLight->fill(xp);
            if( ( quark && p.pid()>0 ) || ( !quark && p.pid()<0 ))
              _h_RKPlus->fill(xp);
            else
              _h_RKMinus->fill(xp);
            break;
         break;
          case PID::CQUARK:
            _multKPlus[2]->fill();
            _h_XpKPlusCharm->fill(xp);
            _temp_XpKPlusCharm->fill(xp);
            break;
          case PID::BQUARK:
            _multKPlus[3]->fill();
            _h_XpKPlusBottom->fill(xp);
            break;
          }
        }
        else if (id == PID::PROTON) {
          _h_XpProtonN->fill(xp);
          _multProton[0]->fill();
          switch (flavour) {
          case PID::DQUARK:
          case PID::UQUARK:
          case PID::SQUARK:
            _multProton[1]->fill();
            _temp_XpProtonLight->fill(xp);
            _h_XpProtonLight->fill(xp);
            if( ( quark && p.pid()>0 ) || ( !quark && p.pid()<0 ))
              _h_RProton->fill(xp);
            else
              _h_RPBar  ->fill(xp);
            break;
         break;
          case PID::CQUARK:
            _multProton[2]->fill();
            _temp_XpProtonCharm->fill(xp);
            _h_XpProtonCharm->fill(xp);
            break;
          case PID::BQUARK:
            _multProton[3]->fill();
            _h_XpProtonBottom->fill(xp);
            break;
          }
        }
      }

      const UnstableParticles& ufs = apply<UnstableFinalState>(e, "UFS");
      for (const Particle& p : ufs.particles()) {
        const double xp = p.p3().mod()/meanBeamMom;
        // if in quark or antiquark hemisphere
        bool quark = p.p3().dot(axis)*dot>0.;
        int id = p.abspid();
        if (id == PID::LAMBDA) {
          _multLambda[0]->fill();
          _h_XpLambdaN->fill(xp);
          switch (flavour) {
          case PID::DQUARK:
          case PID::UQUARK:
          case PID::SQUARK:
            _multLambda[1]->fill();
            _h_XpLambdaLight->fill(xp);
            if( ( quark && p.pid()>0 ) || ( !quark && p.pid()<0 ))
              _h_RLambda->fill(xp);
            else
              _h_RLBar  ->fill(xp);
            break;
          case PID::CQUARK:
            _multLambda[2]->fill();
            _h_XpLambdaCharm->fill(xp);
            break;
          case PID::BQUARK:
            _multLambda[3]->fill();
            _h_XpLambdaBottom->fill(xp);
            break;
          }
        }
        else if (id == 313) {
          _multKStar0[0]->fill();
          _h_XpKStar0N->fill(xp);
          switch (flavour) {
          case PID::DQUARK:
          case PID::UQUARK:
          case PID::SQUARK:
            _multKStar0[1]->fill();
            _temp_XpKStar0Light->fill(xp);
            _h_XpKStar0Light->fill(xp);
            if ( ( quark && p.pid()>0 ) || ( !quark && p.pid()<0 ))
              _h_RKS0   ->fill(xp);
            else
              _h_RKSBar0->fill(xp);
            break;
            break;
          case PID::CQUARK:
            _multKStar0[2]->fill();
            _temp_XpKStar0Charm->fill(xp);
            _h_XpKStar0Charm->fill(xp);
            break;
          case PID::BQUARK:
            _multKStar0[3]->fill();
            _h_XpKStar0Bottom->fill(xp);
            break;
          }
        }
        else if (id == 333) {
          _multPhi[0]->fill();
          _h_XpPhiN->fill(xp);
          switch (flavour) {
          case PID::DQUARK:
          case PID::UQUARK:
          case PID::SQUARK:
            _multPhi[1]->fill();
            _h_XpPhiLight->fill(xp);
            break;
          case PID::CQUARK:
            _multPhi[2]->fill();
            _h_XpPhiCharm->fill(xp);
            break;
          case PID::BQUARK:
            _multPhi[3]->fill();
            _h_XpPhiBottom->fill(xp);
            break;
          }
        }
        else if (id == PID::K0S || id == PID::K0L) {
          _multK0[0]->fill();
          _h_XpK0N->fill(xp);
          switch (flavour) {
          case PID::DQUARK:
          case PID::UQUARK:
          case PID::SQUARK:
            _multK0[1]->fill();
            _h_XpK0Light->fill(xp);
            break;
          case PID::CQUARK:
            _multK0[2]->fill();
            _h_XpK0Charm->fill(xp);
            break;
          case PID::BQUARK:
            _multK0[3]->fill();
            _h_XpK0Bottom->fill(xp);
            break;
          }
        }
      }
    }


    void init() {
      // Projections
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      declare(InitialQuarks(), "IQF");
      declare(Thrust(FinalState()), "Thrust");

      book(_temp_XpChargedN1 ,"TMP/XpChargedN1", refData( 1, 1, 1));
      book(_temp_XpChargedN2 ,"TMP/XpChargedN2", refData( 2, 1, 1));
      book(_temp_XpChargedN3 ,"TMP/XpChargedN3", refData( 3, 1, 1));

      book(_h_XpPiPlusN      , 1, 1, 2);
      book(_h_XpKPlusN       , 2, 1, 2);
      book(_h_XpProtonN      , 3, 1, 2);
      book(_h_XpChargedN     , 4, 1, 1);
      book(_h_XpK0N          , 5, 1, 1);
      book(_h_XpLambdaN      , 7, 1, 1);
      book(_h_XpKStar0N      , 8, 1, 1);
      book(_h_XpPhiN         , 9, 1, 1);

      book(_h_XpPiPlusLight  ,10, 1, 1);
      book(_h_XpPiPlusCharm  ,10, 1, 2);
      book(_h_XpPiPlusBottom ,10, 1, 3);
      book(_h_XpKPlusLight   ,12, 1, 1);
      book(_h_XpKPlusCharm   ,12, 1, 2);
      book(_h_XpKPlusBottom  ,12, 1, 3);
      book(_h_XpKStar0Light  ,14, 1, 1);
      book(_h_XpKStar0Charm  ,14, 1, 2);
      book(_h_XpKStar0Bottom ,14, 1, 3);
      book(_h_XpProtonLight  ,16, 1, 1);
      book(_h_XpProtonCharm  ,16, 1, 2);
      book(_h_XpProtonBottom ,16, 1, 3);
      book(_h_XpLambdaLight  ,18, 1, 1);
      book(_h_XpLambdaCharm  ,18, 1, 2);
      book(_h_XpLambdaBottom ,18, 1, 3);
      book(_h_XpK0Light      ,20, 1, 1);
      book(_h_XpK0Charm      ,20, 1, 2);
      book(_h_XpK0Bottom     ,20, 1, 3);
      book(_h_XpPhiLight     ,22, 1, 1);
      book(_h_XpPhiCharm     ,22, 1, 2);
      book(_h_XpPhiBottom    ,22, 1, 3);

      book(_temp_XpKPlusCharm   ,"TMP/XpKPlusCharm", refData(13, 1, 1));
      book(_temp_XpKPlusLight   ,"TMP/XpKPlusLight", refData(13, 1, 1));
      book(_temp_XpKStar0Charm  ,"TMP/XpKStar0Charm", refData(15, 1, 1));
      book(_temp_XpKStar0Light  ,"TMP/XpKStar0Light", refData(15, 1, 1));
      book(_temp_XpProtonCharm  ,"TMP/XpProtonCharm", refData(17, 1, 1));
      book(_temp_XpProtonLight  ,"TMP/XpProtonLight", refData(17, 1, 1));

      book(_h_RPiPlus  , 26, 1, 1);
      book(_h_RPiMinus , 26, 1, 2);
      book(_h_RKS0     , 28, 1, 1);
      book(_h_RKSBar0  , 28, 1, 2);
      book(_h_RKPlus   , 30, 1, 1);
      book(_h_RKMinus  , 30, 1, 2);
      book(_h_RProton  , 32, 1, 1);
      book(_h_RPBar    , 32, 1, 2);
      book(_h_RLambda  , 34, 1, 1);
      book(_h_RLBar    , 34, 1, 2);

      book(_s_Xp_PiPl_Ch      , 1, 1, 1);
      book(_s_Xp_KPl_Ch       , 2, 1, 1);
      book(_s_Xp_Pr_Ch        , 3, 1, 1);
      book(_s_Xp_PiPlCh_PiPlLi, 11, 1, 1);
      book(_s_Xp_PiPlBo_PiPlLi, 11, 1, 2);
      book(_s_Xp_KPlCh_KPlLi  , 13, 1, 1);
      book(_s_Xp_KPlBo_KPlLi  , 13, 1, 2);
      book(_s_Xp_KS0Ch_KS0Li  , 15, 1, 1);
      book(_s_Xp_KS0Bo_KS0Li  , 15, 1, 2);
      book(_s_Xp_PrCh_PrLi    , 17, 1, 1);
      book(_s_Xp_PrBo_PrLi    , 17, 1, 2);
      book(_s_Xp_LaCh_LaLi    , 19, 1, 1);
      book(_s_Xp_LaBo_LaLi    , 19, 1, 2);
      book(_s_Xp_K0Ch_K0Li    , 21, 1, 1);
      book(_s_Xp_K0Bo_K0Li    , 21, 1, 2);
      book(_s_Xp_PhiCh_PhiLi  , 23, 1, 1);
      book(_s_Xp_PhiBo_PhiLi  , 23, 1, 2);

      book(_s_PiM_PiP   , 27, 1, 1);
      book(_s_KSBar0_KS0, 29, 1, 1);
      book(_s_KM_KP     , 31, 1, 1);
      book(_s_Pr_PBar   , 33, 1, 1);
      book(_s_Lam_LBar  , 35, 1, 1);

      book(_SumOfudsWeights, "_SumOfudsWeights");
      book(_SumOfcWeights, "_SumOfcWeights");
      book(_SumOfbWeights, "_SumOfbWeights");

      for ( size_t i=0; i<4; ++i) {
      	book(_multPiPlus[i], "_multPiPlus_"+to_str(i));
      	book(_multKPlus[i], "_multKPlus_"+to_str(i));
      	book(_multK0[i], "_multK0_"+to_str(i));
      	book(_multKStar0[i], "_multKStar0_"+to_str(i));
      	book(_multPhi[i], "_multPhi_"+to_str(i));
      	book(_multProton[i], "_multProton_"+to_str(i));
      	book(_multLambda[i], "_multLambda_"+to_str(i));
      }

      book(tmp1, 24, 1, 1, true);
      book(tmp2, 24, 1, 2, true);
      book(tmp3, 24, 1, 3, true);
      book(tmp4, 24, 1, 4, true);
      book(tmp5, 25, 1, 1, true);
      book(tmp6, 25, 1, 2, true);
      book(tmp7, 24, 2, 1, true);
      book(tmp8, 24, 2, 2, true);
      book(tmp9, 24, 2, 3, true);
      book(tmp10, 24, 2, 4, true);
      book(tmp11, 25, 2, 1, true);
      book(tmp12, 25, 2, 2, true);
      book(tmp13, 24, 3, 1, true);
      book(tmp14, 24, 3, 2, true);
      book(tmp15, 24, 3, 3, true);
      book(tmp16, 24, 3, 4, true);
      book(tmp17, 25, 3, 1, true);
      book(tmp18, 25, 3, 2, true);
      book(tmp19, 24, 4, 1, true);
      book(tmp20, 24, 4, 2, true);
      book(tmp21, 24, 4, 3, true);
      book(tmp22, 24, 4, 4, true);
      book(tmp23, 25, 4, 1, true);
      book(tmp24, 25, 4, 2, true);
      book(tmp25, 24, 5, 1, true);
      book(tmp26, 24, 5, 2, true);
      book(tmp27, 24, 5, 3, true);
      book(tmp28, 24, 5, 4, true);
      book(tmp29, 25, 5, 1, true);
      book(tmp30, 25, 5, 2, true);
      book(tmp31, 24, 6, 1, true);
      book(tmp32, 24, 6, 2, true);
      book(tmp33, 24, 6, 3, true);
      book(tmp34, 24, 6, 4, true);
      book(tmp35, 25, 6, 1, true);
      book(tmp36, 25, 6, 2, true);
      book(tmp37, 24, 7, 1, true);
      book(tmp38, 24, 7, 2, true);
      book(tmp39, 24, 7, 3, true);
      book(tmp40, 24, 7, 4, true);
      book(tmp41, 25, 7, 1, true);
      book(tmp42, 25, 7, 2, true);
    }


    /// Finalize
    void finalize() {
      // Get the ratio plots sorted out first
      divide(_h_XpPiPlusN, _temp_XpChargedN1, _s_Xp_PiPl_Ch);
      divide(_h_XpKPlusN, _temp_XpChargedN2, _s_Xp_KPl_Ch);
      divide(_h_XpProtonN, _temp_XpChargedN3, _s_Xp_Pr_Ch);
      divide(_h_XpPiPlusCharm, _h_XpPiPlusLight, _s_Xp_PiPlCh_PiPlLi);
      _s_Xp_PiPlCh_PiPlLi->scaleY(dbl(*_SumOfudsWeights / *_SumOfcWeights));
      divide(_h_XpPiPlusBottom, _h_XpPiPlusLight, _s_Xp_PiPlBo_PiPlLi);
       _s_Xp_PiPlBo_PiPlLi->scaleY(dbl(*_SumOfudsWeights / *_SumOfbWeights));
      divide(_temp_XpKPlusCharm , _temp_XpKPlusLight, _s_Xp_KPlCh_KPlLi);
      _s_Xp_KPlCh_KPlLi->scaleY(dbl(*_SumOfudsWeights / *_SumOfcWeights));
      divide(_h_XpKPlusBottom, _h_XpKPlusLight, _s_Xp_KPlBo_KPlLi);
       _s_Xp_KPlBo_KPlLi->scaleY(dbl(*_SumOfudsWeights / *_SumOfbWeights));
      divide(_temp_XpKStar0Charm, _temp_XpKStar0Light, _s_Xp_KS0Ch_KS0Li);
      _s_Xp_KS0Ch_KS0Li->scaleY(dbl(*_SumOfudsWeights / *_SumOfcWeights));
      divide(_h_XpKStar0Bottom, _h_XpKStar0Light, _s_Xp_KS0Bo_KS0Li);
      _s_Xp_KS0Bo_KS0Li->scaleY(dbl(*_SumOfudsWeights / *_SumOfbWeights));
      divide(_temp_XpProtonCharm, _temp_XpProtonLight, _s_Xp_PrCh_PrLi);
      _s_Xp_PrCh_PrLi->scaleY(dbl(*_SumOfudsWeights / *_SumOfcWeights));
      divide(_h_XpProtonBottom, _h_XpProtonLight, _s_Xp_PrBo_PrLi);
      _s_Xp_PrBo_PrLi->scaleY(dbl(*_SumOfudsWeights / *_SumOfbWeights));
      divide(_h_XpLambdaCharm, _h_XpLambdaLight, _s_Xp_LaCh_LaLi);
      _s_Xp_LaCh_LaLi->scaleY(dbl(*_SumOfudsWeights / *_SumOfcWeights));
      divide(_h_XpLambdaBottom, _h_XpLambdaLight, _s_Xp_LaBo_LaLi);
      _s_Xp_LaBo_LaLi->scaleY(dbl(*_SumOfudsWeights / *_SumOfbWeights));
      divide(_h_XpK0Charm, _h_XpK0Light, _s_Xp_K0Ch_K0Li);
      _s_Xp_K0Ch_K0Li->scaleY(dbl(*_SumOfudsWeights / *_SumOfcWeights));
      divide(_h_XpK0Bottom, _h_XpK0Light, _s_Xp_K0Bo_K0Li);
      _s_Xp_K0Bo_K0Li->scaleY(dbl(*_SumOfudsWeights / *_SumOfbWeights));
      divide(_h_XpPhiCharm, _h_XpPhiLight, _s_Xp_PhiCh_PhiLi);
      _s_Xp_PhiCh_PhiLi->scaleY(dbl(*_SumOfudsWeights / *_SumOfcWeights));
      divide(_h_XpPhiBottom, _h_XpPhiLight, _s_Xp_PhiBo_PhiLi);
      _s_Xp_PhiBo_PhiLi->scaleY(dbl(*_SumOfudsWeights / *_SumOfbWeights));

      // Then the leading particles
      divide(*_h_RPiMinus - *_h_RPiPlus, *_h_RPiMinus + *_h_RPiPlus, _s_PiM_PiP);
      divide(*_h_RKSBar0 - *_h_RKS0, *_h_RKSBar0 + *_h_RKS0, _s_KSBar0_KS0);
      divide(*_h_RKMinus - *_h_RKPlus, *_h_RKMinus + *_h_RKPlus, _s_KM_KP);
      divide(*_h_RProton - *_h_RPBar, *_h_RProton + *_h_RPBar, _s_Pr_PBar);
      divide(*_h_RLambda - *_h_RLBar, *_h_RLambda + *_h_RLBar, _s_Lam_LBar);

      // Then the rest
      scale(_h_XpPiPlusN,      1/sumOfWeights());
      scale(_h_XpKPlusN,       1/sumOfWeights());
      scale(_h_XpProtonN,      1/sumOfWeights());
      scale(_h_XpChargedN,     1/sumOfWeights());
      scale(_h_XpK0N,          1/sumOfWeights());
      scale(_h_XpLambdaN,      1/sumOfWeights());
      scale(_h_XpKStar0N,      1/sumOfWeights());
      scale(_h_XpPhiN,         1/sumOfWeights());
      scale(_h_XpPiPlusLight,  1 / *_SumOfudsWeights);
      scale(_h_XpPiPlusCharm,  1 / *_SumOfcWeights);
      scale(_h_XpPiPlusBottom, 1 / *_SumOfbWeights);
      scale(_h_XpKPlusLight,   1 / *_SumOfudsWeights);
      scale(_h_XpKPlusCharm,   1 / *_SumOfcWeights);
      scale(_h_XpKPlusBottom,  1 / *_SumOfbWeights);
      scale(_h_XpKStar0Light,  1 / *_SumOfudsWeights);
      scale(_h_XpKStar0Charm,  1 / *_SumOfcWeights);
      scale(_h_XpKStar0Bottom, 1 / *_SumOfbWeights);
      scale(_h_XpProtonLight,  1 / *_SumOfudsWeights);
      scale(_h_XpProtonCharm,  1 / *_SumOfcWeights);
      scale(_h_XpProtonBottom, 1 / *_SumOfbWeights);
      scale(_h_XpLambdaLight,  1 / *_SumOfudsWeights);
      scale(_h_XpLambdaCharm,  1 / *_SumOfcWeights);
      scale(_h_XpLambdaBottom, 1 / *_SumOfbWeights);
      scale(_h_XpK0Light,      1 / *_SumOfudsWeights);
      scale(_h_XpK0Charm,      1 / *_SumOfcWeights);
      scale(_h_XpK0Bottom,     1 / *_SumOfbWeights);
      scale(_h_XpPhiLight,     1 / *_SumOfudsWeights);
      scale(_h_XpPhiCharm ,    1 / *_SumOfcWeights);
      scale(_h_XpPhiBottom,    1 / *_SumOfbWeights);
      scale(_h_RPiPlus,        1 / *_SumOfudsWeights);
      scale(_h_RPiMinus,       1 / *_SumOfudsWeights);
      scale(_h_RKS0,           1 / *_SumOfudsWeights);
      scale(_h_RKSBar0,        1 / *_SumOfudsWeights);
      scale(_h_RKPlus,         1 / *_SumOfudsWeights);
      scale(_h_RKMinus,        1 / *_SumOfudsWeights);
      scale(_h_RProton,        1 / *_SumOfudsWeights);
      scale(_h_RPBar,          1 / *_SumOfudsWeights);
      scale(_h_RLambda,        1 / *_SumOfudsWeights);
      scale(_h_RLBar,          1 / *_SumOfudsWeights);

      // Multiplicities
      double avgNumPartsAll, avgNumPartsLight,avgNumPartsCharm, avgNumPartsBottom;
      // pi+/-
      // all
      avgNumPartsAll = dbl(*_multPiPlus[0])/sumOfWeights();
      tmp1->point(0).setY(avgNumPartsAll);
      // light
      avgNumPartsLight = dbl(*_multPiPlus[1] / *_SumOfudsWeights);
      tmp2->point(0).setY(avgNumPartsLight);
      // charm
      avgNumPartsCharm = dbl(*_multPiPlus[2] / *_SumOfcWeights);
      tmp3->point(0).setY(avgNumPartsCharm);
      // bottom
      avgNumPartsBottom = dbl(*_multPiPlus[3] / *_SumOfbWeights);
      tmp4->point(0).setY(avgNumPartsBottom);
      // charm-light
      tmp5->point(0).setY(avgNumPartsCharm - avgNumPartsLight);
      // bottom-light
      tmp6->point(0).setY(avgNumPartsBottom - avgNumPartsLight);
      // K+/-
      // all
      avgNumPartsAll = dbl(*_multKPlus[0])/sumOfWeights();
      tmp7->point(0).setY(avgNumPartsAll);
      // light
      avgNumPartsLight = dbl(*_multKPlus[1] / *_SumOfudsWeights);
      tmp8->point(0).setY(avgNumPartsLight);
      // charm
      avgNumPartsCharm = dbl(*_multKPlus[2] / *_SumOfcWeights);
      tmp9->point(0).setY(avgNumPartsCharm);
      // bottom
      avgNumPartsBottom = dbl(*_multKPlus[3] / *_SumOfbWeights);
      tmp10->point(0).setY(avgNumPartsBottom);
      // charm-light
      tmp11->point(0).setY(avgNumPartsCharm - avgNumPartsLight);
      // bottom-light
      tmp12->point(0).setY(avgNumPartsBottom - avgNumPartsLight);
      // K0
      // all
      avgNumPartsAll = dbl(*_multK0[0])/sumOfWeights();
      tmp13->point(0).setY(avgNumPartsAll);
      // light
      avgNumPartsLight = dbl(*_multK0[1] / *_SumOfudsWeights);
      tmp14->point(0).setY(avgNumPartsLight);
      // charm
      avgNumPartsCharm = dbl(*_multK0[2] / *_SumOfcWeights);
      tmp15->point(0).setY(avgNumPartsCharm);
      // bottom
      avgNumPartsBottom = dbl(*_multK0[3] / *_SumOfbWeights);
      tmp16->point(0).setY(avgNumPartsBottom);
      // charm-light
      tmp17->point(0).setY(avgNumPartsCharm - avgNumPartsLight);
      // bottom-light
      tmp18->point(0).setY(avgNumPartsBottom - avgNumPartsLight);
      // K*0
      // all
      avgNumPartsAll = dbl(*_multKStar0[0])/sumOfWeights();
      tmp19->point(0).setY(avgNumPartsAll);
      // light
      avgNumPartsLight = dbl(*_multKStar0[1] / *_SumOfudsWeights);
      tmp20->point(0).setY(avgNumPartsLight);
      // charm
      avgNumPartsCharm = dbl(*_multKStar0[2] / *_SumOfcWeights);
      tmp21->point(0).setY(avgNumPartsCharm);
      // bottom
      avgNumPartsBottom = dbl(*_multKStar0[3] / *_SumOfbWeights);
      tmp22->point(0).setY(avgNumPartsBottom);
      // charm-light
      tmp23->point(0).setY(avgNumPartsCharm - avgNumPartsLight);
      // bottom-light
      tmp24->point(0).setY(avgNumPartsBottom - avgNumPartsLight);
      // phi
      // all
      avgNumPartsAll = dbl(*_multPhi[0])/sumOfWeights();
      tmp25->point(0).setY(avgNumPartsAll);
      // light
      avgNumPartsLight = dbl(*_multPhi[1] / *_SumOfudsWeights);
      tmp26->point(0).setY(avgNumPartsLight);
      // charm
      avgNumPartsCharm = dbl(*_multPhi[2] / *_SumOfcWeights);
      tmp27->point(0).setY(avgNumPartsCharm);
      // bottom
      avgNumPartsBottom = dbl(*_multPhi[3] / *_SumOfbWeights);
      tmp28->point(0).setY(avgNumPartsBottom);
      // charm-light
      tmp29->point(0).setY(avgNumPartsCharm - avgNumPartsLight);
      // bottom-light
      tmp30->point(0).setY(avgNumPartsBottom - avgNumPartsLight);
      // p
      // all
      avgNumPartsAll = dbl(*_multProton[0])/sumOfWeights();
      tmp31->point(0).setY(avgNumPartsAll);
      // light
      avgNumPartsLight = dbl(*_multProton[1] / *_SumOfudsWeights);
      tmp32->point(0).setY(avgNumPartsLight);
      // charm
      avgNumPartsCharm = dbl(*_multProton[2] / *_SumOfcWeights);
      tmp33->point(0).setY(avgNumPartsCharm);
      // bottom
      avgNumPartsBottom = dbl(*_multProton[3] / *_SumOfbWeights);
      tmp34->point(0).setY(avgNumPartsBottom);
      // charm-light
      tmp35->point(0).setY(avgNumPartsCharm - avgNumPartsLight);
      // bottom-light
      tmp36->point(0).setY(avgNumPartsBottom - avgNumPartsLight);
      // Lambda
      // all
      avgNumPartsAll = dbl(*_multLambda[0])/sumOfWeights();
      tmp37->point(0).setY(avgNumPartsAll);
      // light
      avgNumPartsLight = dbl(*_multLambda[1] / *_SumOfudsWeights);
      tmp38->point(0).setY(avgNumPartsLight);
      // charm
      avgNumPartsCharm = dbl(*_multLambda[2] / *_SumOfcWeights);
      tmp39->point(0).setY(avgNumPartsCharm);
      // bottom
      avgNumPartsBottom = dbl(*_multLambda[3] / *_SumOfbWeights);
      tmp40->point(0).setY(avgNumPartsBottom);
      // charm-light
      tmp41->point(0).setY(avgNumPartsCharm - avgNumPartsLight);
      // bottom-light
      tmp42->point(0).setY(avgNumPartsBottom - avgNumPartsLight);
    }

    //@}


  private:

    /// Store the weighted sums of numbers of charged / charged+neutral
    /// particles. Used to calculate average number of particles for the
    /// inclusive single particle distributions' normalisations.
    CounterPtr _SumOfudsWeights, _SumOfcWeights, _SumOfbWeights;
    vector<CounterPtr> _multPiPlus, _multKPlus, _multK0,
      _multKStar0, _multPhi, _multProton, _multLambda;

    Histo1DPtr _h_XpPiPlusSig, _h_XpPiPlusN;
    Histo1DPtr _h_XpKPlusSig, _h_XpKPlusN;
    Histo1DPtr _h_XpProtonSig, _h_XpProtonN;
    Histo1DPtr _h_XpChargedN;
    Histo1DPtr _h_XpK0N, _h_XpLambdaN;
    Histo1DPtr _h_XpKStar0N, _h_XpPhiN;
    Histo1DPtr _h_XpPiPlusLight, _h_XpPiPlusCharm, _h_XpPiPlusBottom;
    Histo1DPtr _h_XpKPlusLight, _h_XpKPlusCharm, _h_XpKPlusBottom;
    Histo1DPtr _h_XpKStar0Light, _h_XpKStar0Charm, _h_XpKStar0Bottom;
    Histo1DPtr _h_XpProtonLight, _h_XpProtonCharm, _h_XpProtonBottom;
    Histo1DPtr _h_XpLambdaLight, _h_XpLambdaCharm, _h_XpLambdaBottom;
    Histo1DPtr _h_XpK0Light, _h_XpK0Charm, _h_XpK0Bottom;
    Histo1DPtr _h_XpPhiLight, _h_XpPhiCharm, _h_XpPhiBottom;

    Histo1DPtr _temp_XpChargedN1, _temp_XpChargedN2, _temp_XpChargedN3;
    Histo1DPtr _temp_XpKPlusCharm , _temp_XpKPlusLight;
    Histo1DPtr _temp_XpKStar0Charm, _temp_XpKStar0Light;
    Histo1DPtr _temp_XpProtonCharm, _temp_XpProtonLight;

    Histo1DPtr _h_RPiPlus, _h_RPiMinus;
    Histo1DPtr _h_RKS0, _h_RKSBar0;
    Histo1DPtr _h_RKPlus, _h_RKMinus;
    Histo1DPtr _h_RProton, _h_RPBar;
    Histo1DPtr _h_RLambda, _h_RLBar;

    Scatter2DPtr _s_Xp_PiPl_Ch, _s_Xp_KPl_Ch,  _s_Xp_Pr_Ch;
    Scatter2DPtr _s_Xp_PiPlCh_PiPlLi, _s_Xp_PiPlBo_PiPlLi;
    Scatter2DPtr _s_Xp_KPlCh_KPlLi, _s_Xp_KPlBo_KPlLi;
    Scatter2DPtr _s_Xp_KS0Ch_KS0Li, _s_Xp_KS0Bo_KS0Li;
    Scatter2DPtr _s_Xp_PrCh_PrLi, _s_Xp_PrBo_PrLi;
    Scatter2DPtr _s_Xp_LaCh_LaLi, _s_Xp_LaBo_LaLi;
    Scatter2DPtr _s_Xp_K0Ch_K0Li, _s_Xp_K0Bo_K0Li;
    Scatter2DPtr _s_Xp_PhiCh_PhiLi, _s_Xp_PhiBo_PhiLi;

    Scatter2DPtr _s_PiM_PiP, _s_KSBar0_KS0, _s_KM_KP, _s_Pr_PBar, _s_Lam_LBar;

    //@}
      Scatter2DPtr tmp1;
      Scatter2DPtr tmp2;
      Scatter2DPtr tmp3;
      Scatter2DPtr tmp4;
      Scatter2DPtr tmp5;
      Scatter2DPtr tmp6;
      Scatter2DPtr tmp7;
      Scatter2DPtr tmp8;
      Scatter2DPtr tmp9;
      Scatter2DPtr tmp10;
      Scatter2DPtr tmp11;
      Scatter2DPtr tmp12;
      Scatter2DPtr tmp13;
      Scatter2DPtr tmp14;
      Scatter2DPtr tmp15;
      Scatter2DPtr tmp16;
      Scatter2DPtr tmp17;
      Scatter2DPtr tmp18;
      Scatter2DPtr tmp19;
      Scatter2DPtr tmp20;
      Scatter2DPtr tmp21;
      Scatter2DPtr tmp22;
      Scatter2DPtr tmp23;
      Scatter2DPtr tmp24;
      Scatter2DPtr tmp25;
      Scatter2DPtr tmp26;
      Scatter2DPtr tmp27;
      Scatter2DPtr tmp28;
      Scatter2DPtr tmp29;
      Scatter2DPtr tmp30;
      Scatter2DPtr tmp31;
      Scatter2DPtr tmp32;
      Scatter2DPtr tmp33;
      Scatter2DPtr tmp34;
      Scatter2DPtr tmp35;
      Scatter2DPtr tmp36;
      Scatter2DPtr tmp37;
      Scatter2DPtr tmp38;
      Scatter2DPtr tmp39;
      Scatter2DPtr tmp40;
      Scatter2DPtr tmp41;
      Scatter2DPtr tmp42;
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(SLD_1999_S3743934);

}
