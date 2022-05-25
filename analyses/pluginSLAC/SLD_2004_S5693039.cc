// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Thrust.hh"

#define I_KNOW_THE_INITIAL_QUARKS_PROJECTION_IS_DODGY_BUT_NEED_TO_USE_IT
#include "Rivet/Projections/InitialQuarks.hh"

namespace Rivet {


  /// @brief SLD flavour-dependent fragmentation paper
  ///
  /// @author Peter Richardson
  class SLD_2004_S5693039 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SLD_2004_S5693039);


    /// @name Analysis methods
    //@{

    void analyze(const Event& e) {
      // First, veto on leptonic events by requiring at least 2 charged FS particles
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
      Particles quarks;
      if (iqf.particles().size() == 2) {
        flavour = iqf.particles().front().abspid();
        quarks = iqf.particles();
      }
      else {
        map<int, Particle > quarkmap;
        for (const Particle& p : iqf.particles()) {
          if (quarkmap.find(p.pid())==quarkmap.end())
            quarkmap[p.pid()] = p;
          else if (quarkmap[p.pid()].E() < p.E())
            quarkmap[p.pid()] = p;
        }
        double maxenergy = 0.;
        for (int i = 1; i <= 5; ++i) {
          double energy(0.);
          if(quarkmap.find( i)!=quarkmap.end())
            energy += quarkmap[ i].E();
          if(quarkmap.find(-i)!=quarkmap.end())
            energy += quarkmap[-i].E();
          if (energy > maxenergy) flavour = i;
        }
        if(quarkmap.find( flavour)!=quarkmap.end())
          quarks.push_back(quarkmap[ flavour]);
        if(quarkmap.find(-flavour)!=quarkmap.end())
          quarks.push_back(quarkmap[-flavour]);
      }

      // total multiplicities
      switch (flavour) {
      case PID::DQUARK:
      case PID::UQUARK:
      case PID::SQUARK:
        _weightLight  ->fill();
        _weightedTotalChargedPartNumLight  ->fill(numParticles);
        break;
      case PID::CQUARK:
        _weightCharm  ->fill();
        _weightedTotalChargedPartNumCharm  ->fill(numParticles);
        break;
      case PID::BQUARK:
        _weightBottom ->fill();
        _weightedTotalChargedPartNumBottom ->fill(numParticles);
        break;
      }
      // thrust axis for projections
      Vector3 axis = apply<Thrust>(e, "Thrust").thrustAxis();
      double dot(0.);
      if(!quarks.empty()) {
        dot = quarks[0].p3().dot(axis);
        if(quarks[0].pid()<0) dot *= -1.;
      }
      // spectra and individual multiplicities
      for (const Particle& p : fs.particles()) {
        double pcm = p.p3().mod();
        const double xp = pcm/meanBeamMom;

        // if in quark or antiquark hemisphere
        bool quark = p.p3().dot(axis)*dot>0.;

        _h_PCharged ->fill(pcm     );
        // all charged
        switch (flavour) {
        case PID::DQUARK:
        case PID::UQUARK:
        case PID::SQUARK:
          _h_XpChargedL->fill(xp);
          break;
        case PID::CQUARK:
          _h_XpChargedC->fill(xp);
          break;
        case PID::BQUARK:
          _h_XpChargedB->fill(xp);
          break;
        }

        int id = p.abspid();
        // charged pions
        if (id == PID::PIPLUS) {
          _h_XpPiPlus->fill(xp);
          _h_XpPiPlusTotal->fill(xp);
          switch (flavour) {
          case PID::DQUARK:
          case PID::UQUARK:
          case PID::SQUARK:
            _h_XpPiPlusL->fill(xp);
            _h_NPiPlusL->fill(sqrtS());
            if( ( quark && p.pid()>0 ) || ( !quark && p.pid()<0 ))
              _h_RPiPlus->fill(xp);
            else
              _h_RPiMinus->fill(xp);
            break;
          case PID::CQUARK:
            _h_XpPiPlusC->fill(xp);
            _h_NPiPlusC->fill(sqrtS());
            break;
          case PID::BQUARK:
            _h_XpPiPlusB->fill(xp);
            _h_NPiPlusB->fill(sqrtS());
            break;
          }
        }
        else if (id == PID::KPLUS) {
          _h_XpKPlus->fill(xp);
          _h_XpKPlusTotal->fill(xp);
          switch (flavour) {
          case PID::DQUARK:
          case PID::UQUARK:
          case PID::SQUARK:
            _h_XpKPlusL->fill(xp);
            _h_NKPlusL->fill(sqrtS());
            if( ( quark && p.pid()>0 ) || ( !quark && p.pid()<0 ))
              _h_RKPlus->fill(xp);
            else
              _h_RKMinus->fill(xp);
            break;
          case PID::CQUARK:
            _h_XpKPlusC->fill(xp);
            _h_NKPlusC->fill(sqrtS());
            break;
          case PID::BQUARK:
            _h_XpKPlusB->fill(xp);
            _h_NKPlusB->fill(sqrtS());
            break;
          }
        }
        else if (id == PID::PROTON) {
          _h_XpProton->fill(xp);
          _h_XpProtonTotal->fill(xp);
          switch (flavour) {
          case PID::DQUARK:
          case PID::UQUARK:
          case PID::SQUARK:
            _h_XpProtonL->fill(xp);
            _h_NProtonL->fill(sqrtS());
            if( ( quark && p.pid()>0 ) || ( !quark && p.pid()<0 ))
              _h_RProton->fill(xp);
            else
              _h_RPBar  ->fill(xp);
            break;
          case PID::CQUARK:
            _h_XpProtonC->fill(xp);
            _h_NProtonC->fill(sqrtS());
            break;
          case PID::BQUARK:
            _h_XpProtonB->fill(xp);
            _h_NProtonB->fill(sqrtS());
            break;
          }
        }
      }
    }


    void init() {
      // Projections
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(InitialQuarks(), "IQF");
      declare(Thrust(FinalState()), "Thrust");

      // Book histograms
      book(_h_PCharged   , 1, 1, 1);
      book(_h_XpPiPlus   , 2, 1, 2);
      book(_h_XpKPlus    , 3, 1, 2);
      book(_h_XpProton   , 4, 1, 2);
      book(_h_XpPiPlusTotal , 2, 2, 2);
      book(_h_XpKPlusTotal  , 3, 2, 2);
      book(_h_XpProtonTotal , 4, 2, 2);
      book(_h_XpPiPlusL  , 5, 1, 1);
      book(_h_XpPiPlusC  , 5, 1, 2);
      book(_h_XpPiPlusB  , 5, 1, 3);
      book(_h_XpKPlusL   , 6, 1, 1);
      book(_h_XpKPlusC   , 6, 1, 2);
      book(_h_XpKPlusB   , 6, 1, 3);
      book(_h_XpProtonL  , 7, 1, 1);
      book(_h_XpProtonC  , 7, 1, 2);
      book(_h_XpProtonB  , 7, 1, 3);
      book(_h_XpChargedL , 8, 1, 1);
      book(_h_XpChargedC , 8, 1, 2);
      book(_h_XpChargedB , 8, 1, 3);

      book(_h_NPiPlusL  , 5, 2, 1);
      book(_h_NPiPlusC  , 5, 2, 2);
      book(_h_NPiPlusB  , 5, 2, 3);
      book(_h_NKPlusL   , 6, 2, 1);
      book(_h_NKPlusC   , 6, 2, 2);
      book(_h_NKPlusB   , 6, 2, 3);
      book(_h_NProtonL  , 7, 2, 1);
      book(_h_NProtonC  , 7, 2, 2);
      book(_h_NProtonB  , 7, 2, 3);

      book(_h_RPiPlus  , 9, 1, 1);
      book(_h_RPiMinus , 9, 1, 2);
      book(_h_RKPlus   ,10, 1, 1);
      book(_h_RKMinus  ,10, 1, 2);
      book(_h_RProton  ,11, 1, 1);
      book(_h_RPBar    ,11, 1, 2);

      // Ratios: used as target of divide() later
      book(_s_PiM_PiP,  9, 1, 3);
      book(_s_KM_KP  , 10, 1, 3);
      book(_s_Pr_PBar, 11, 1, 3);

      book(_weightedTotalChargedPartNumLight, "_weightedTotalChargedPartNumLight");
      book(_weightedTotalChargedPartNumCharm, "_weightedTotalChargedPartNumCharm");
      book(_weightedTotalChargedPartNumBottom, "_weightedTotalChargedPartNumBottom");
      book(_weightLight, "_weightLight");
      book(_weightCharm, "_weightCharm");
      book(_weightBottom, "_weightBottom");

      book(tmp1, 8, 2, 1, true);
      book(tmp2, 8, 2, 2, true);
      book(tmp3, 8, 2, 3, true);
      book(tmp4, 8, 3, 2, true);
      book(tmp5, 8, 3, 3, true);


    }


    /// Finalize
    void finalize() {

      // Multiplicities
      /// @todo Include errors
      const double avgNumPartsLight = _weightedTotalChargedPartNumLight->val() / _weightLight->val();
      const double avgNumPartsCharm = _weightedTotalChargedPartNumCharm->val() / _weightCharm->val();
      const double avgNumPartsBottom = _weightedTotalChargedPartNumBottom->val() / _weightBottom->val();
      tmp1->point(0).setY(avgNumPartsLight);
      tmp2->point(0).setY(avgNumPartsCharm);
      tmp3->point(0).setY(avgNumPartsBottom);
      tmp4->point(0).setY(avgNumPartsCharm - avgNumPartsLight);
      tmp5->point(0).setY(avgNumPartsBottom - avgNumPartsLight);

      // Do divisions
      divide(*_h_RPiMinus - *_h_RPiPlus, *_h_RPiMinus + *_h_RPiPlus, _s_PiM_PiP);
      divide(*_h_RKMinus - *_h_RKPlus, *_h_RKMinus + *_h_RKPlus, _s_KM_KP);
      divide(*_h_RProton - *_h_RPBar, *_h_RProton + *_h_RPBar, _s_Pr_PBar);

      // Scale histograms
      scale(_h_PCharged,      1./sumOfWeights());
      scale(_h_XpPiPlus,      1./sumOfWeights());
      scale(_h_XpKPlus,       1./sumOfWeights());
      scale(_h_XpProton,      1./sumOfWeights());
      scale(_h_XpPiPlusTotal, 1./sumOfWeights());
      scale(_h_XpKPlusTotal,  1./sumOfWeights());
      scale(_h_XpProtonTotal, 1./sumOfWeights());
      scale(_h_XpPiPlusL,     1. / *_weightLight);
      scale(_h_XpPiPlusC,     1. / *_weightCharm);
      scale(_h_XpPiPlusB,     1. / *_weightBottom);
      scale(_h_XpKPlusL,      1. / *_weightLight);
      scale(_h_XpKPlusC,      1. / *_weightCharm);
      scale(_h_XpKPlusB,      1. / *_weightBottom);
      scale(_h_XpProtonL,     1. / *_weightLight);
      scale(_h_XpProtonC,     1. / *_weightCharm);
      scale(_h_XpProtonB,     1. / *_weightBottom);

      scale(_h_XpChargedL, 1. / *_weightLight);
      scale(_h_XpChargedC, 1. / *_weightCharm);
      scale(_h_XpChargedB, 1. / *_weightBottom);

      scale(_h_NPiPlusL, 1. / *_weightLight);
      scale(_h_NPiPlusC, 1. / *_weightCharm);
      scale(_h_NPiPlusB, 1. / *_weightBottom);
      scale(_h_NKPlusL,  1. / *_weightLight);
      scale(_h_NKPlusC,  1. / *_weightCharm);
      scale(_h_NKPlusB,  1. / *_weightBottom);
      scale(_h_NProtonL, 1. / *_weightLight);
      scale(_h_NProtonC, 1. / *_weightCharm);
      scale(_h_NProtonB, 1. / *_weightBottom);

      // Paper suggests this should be 0.5/weight but it has to be 1.0 to get normalisations right...
      scale(_h_RPiPlus,  1. / *_weightLight);
      scale(_h_RPiMinus, 1. / *_weightLight);
      scale(_h_RKPlus,   1. / *_weightLight);
      scale(_h_RKMinus,  1. / *_weightLight);
      scale(_h_RProton,  1. / *_weightLight);
      scale(_h_RPBar,    1. / *_weightLight);

      // convert ratio to %
      _s_PiM_PiP->scaleY(100.);
      _s_KM_KP  ->scaleY(100.);
      _s_Pr_PBar->scaleY(100.);
    }

    //@}


  private:

    Scatter2DPtr tmp1, tmp2, tmp3, tmp4, tmp5;

    /// Multiplicities
    CounterPtr _weightedTotalChargedPartNumLight, _weightedTotalChargedPartNumCharm, _weightedTotalChargedPartNumBottom;

    /// Weights
    CounterPtr _weightLight, _weightCharm, _weightBottom;

    // Histograms
    /// @{
    Histo1DPtr _h_PCharged;
    Histo1DPtr _h_XpPiPlus, _h_XpKPlus, _h_XpProton;
    Histo1DPtr _h_XpPiPlusTotal, _h_XpKPlusTotal, _h_XpProtonTotal;
    Histo1DPtr _h_XpPiPlusL, _h_XpPiPlusC, _h_XpPiPlusB;
    Histo1DPtr _h_XpKPlusL, _h_XpKPlusC, _h_XpKPlusB;
    Histo1DPtr _h_XpProtonL, _h_XpProtonC, _h_XpProtonB;
    Histo1DPtr _h_XpChargedL, _h_XpChargedC, _h_XpChargedB;
    Histo1DPtr _h_NPiPlusL, _h_NPiPlusC, _h_NPiPlusB;
    Histo1DPtr _h_NKPlusL, _h_NKPlusC, _h_NKPlusB;
    Histo1DPtr _h_NProtonL, _h_NProtonC, _h_NProtonB;
    Histo1DPtr _h_RPiPlus, _h_RPiMinus, _h_RKPlus;
    Histo1DPtr _h_RKMinus, _h_RProton, _h_RPBar;
    Scatter2DPtr _s_PiM_PiP, _s_KM_KP, _s_Pr_PBar;
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(SLD_2004_S5693039, SLD_2004_I630327);

}
