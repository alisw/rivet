// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

#define I_KNOW_THE_INITIAL_QUARKS_PROJECTION_IS_DODGY_BUT_NEED_TO_USE_IT
#include "Rivet/Projections/InitialQuarks.hh"

namespace Rivet {


  /// @brief OPAL flavour-dependent fragmentation paper
  /// @author Hendrik Hoeth
  class OPAL_1998_S3780481 : public Analysis {
  public:

    /// Constructor
    OPAL_1998_S3780481() : Analysis("OPAL_1998_S3780481") {
    }


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

      _weightedTotalPartNum->fill(numParticles);

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(e, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      int flavour = 0;
      const InitialQuarks& iqf = apply<InitialQuarks>(e, "IQF");

      // If we only have two quarks (qqbar), just take the flavour.
      // If we have more than two quarks, look for the highest energetic q-qbar pair.
      /// @todo Yuck... does this *really* have to be quark-based?!?
      if (iqf.particles().size() == 2) {
        flavour = iqf.particles().front().abspid();
      } else {
        map<int, double> quarkmap;
        for (const Particle& p : iqf.particles()) {
          if (quarkmap[p.pid()] < p.E()) {
            quarkmap[p.pid()] = p.E();
          }
        }
        double maxenergy = 0.;
        for (int i = 1; i <= 5; ++i) {
          if (quarkmap[i]+quarkmap[-i] > maxenergy) {
            flavour = i;
          }
        }
      }

      switch (flavour) {
      case 1:
      case 2:
      case 3:
        _SumOfudsWeights->fill();
        break;
      case 4:
        _SumOfcWeights->fill();
        break;
      case 5:
        _SumOfbWeights->fill();
        break;
      }

      for (const Particle& p : fs.particles()) {
        const double xp = p.p3().mod()/meanBeamMom;
        const double logxp = -std::log(xp);
        _histXpall->fill(xp);
        _histLogXpall->fill(logxp);
        _histMultiChargedall->fill(_histMultiChargedall->bin(0).xMid());
        switch (flavour) {
          /// @todo Use PDG code enums
        case PID::DQUARK:
        case PID::UQUARK:
        case PID::SQUARK:
          _histXpuds->fill(xp);
          _histLogXpuds->fill(logxp);
          _histMultiChargeduds->fill(_histMultiChargeduds->bin(0).xMid());
          break;
        case PID::CQUARK:
          _histXpc->fill(xp);
          _histLogXpc->fill(logxp);
          _histMultiChargedc->fill(_histMultiChargedc->bin(0).xMid());
          break;
        case PID::BQUARK:
          _histXpb->fill(xp);
          _histLogXpb->fill(logxp);
          _histMultiChargedb->fill(_histMultiChargedb->bin(0).xMid());
          break;
        }
      }

    }


    void init() {
      // Projections
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(InitialQuarks(), "IQF");

      // Book histos
      book(_histXpuds           ,1, 1, 1);
      book(_histXpc             ,2, 1, 1);
      book(_histXpb             ,3, 1, 1);
      book(_histXpall           ,4, 1, 1);
      book(_histLogXpuds        ,5, 1, 1);
      book(_histLogXpc          ,6, 1, 1);
      book(_histLogXpb          ,7, 1, 1);
      book(_histLogXpall        ,8, 1, 1);
      book(_histMultiChargeduds ,9, 1, 1);
      book(_histMultiChargedc   ,9, 1, 2);
      book(_histMultiChargedb   ,9, 1, 3);
      book(_histMultiChargedall ,9, 1, 4);
      // Counters
      book(_weightedTotalPartNum, "_TotalPartNum");
      book(_SumOfudsWeights, "_udsWeights");
      book(_SumOfcWeights, "_cWeights");
      book(_SumOfbWeights, "_bWeights");
    }


    /// Finalize
    void finalize() {
      const double avgNumParts = dbl(*_weightedTotalPartNum) / sumOfWeights();
      normalize(_histXpuds    , avgNumParts);
      normalize(_histXpc      , avgNumParts);
      normalize(_histXpb      , avgNumParts);
      normalize(_histXpall    , avgNumParts);
      normalize(_histLogXpuds , avgNumParts);
      normalize(_histLogXpc   , avgNumParts);
      normalize(_histLogXpb   , avgNumParts);
      normalize(_histLogXpall , avgNumParts);

      scale(_histMultiChargeduds, 1.0/ *_SumOfudsWeights);
      scale(_histMultiChargedc  , 1.0/ *_SumOfcWeights);
      scale(_histMultiChargedb  , 1.0/ *_SumOfbWeights);
      scale(_histMultiChargedall, 1.0/sumOfWeights());
    }

    //@}


  private:

    /// Store the weighted sums of numbers of charged / charged+neutral
    /// particles - used to calculate average number of particles for the
    /// inclusive single particle distributions' normalisations.
    CounterPtr _weightedTotalPartNum;

    CounterPtr _SumOfudsWeights;
    CounterPtr _SumOfcWeights;
    CounterPtr _SumOfbWeights;

    Histo1DPtr _histXpuds;
    Histo1DPtr _histXpc;
    Histo1DPtr _histXpb;
    Histo1DPtr _histXpall;
    Histo1DPtr _histLogXpuds;
    Histo1DPtr _histLogXpc;
    Histo1DPtr _histLogXpb;
    Histo1DPtr _histLogXpall;
    Histo1DPtr _histMultiChargeduds;
    Histo1DPtr _histMultiChargedc;
    Histo1DPtr _histMultiChargedb;
    Histo1DPtr _histMultiChargedall;

    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(OPAL_1998_S3780481);

}
