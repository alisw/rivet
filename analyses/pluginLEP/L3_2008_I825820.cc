// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include "Rivet/Projections/ParisiTensor.hh"
#include "Rivet/Projections/FastJets.hh"

#define I_KNOW_THE_INITIAL_QUARKS_PROJECTION_IS_DODGY_BUT_NEED_TO_USE_IT
#include "Rivet/Projections/InitialQuarks.hh"

namespace Rivet {


  /// @brief Event shapes at 197 GeV
  class L3_2008_I825820 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(L3_2008_I825820);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Projections to use
      const FinalState FS;
      declare(FS, "FS");
      const ChargedFinalState CFS;
      declare(CFS, "CFS");
      const Thrust thrust(FS);
      declare(thrust, "Thrust");
      declare(ParisiTensor(FS), "Parisi");
      declare(Hemispheres(thrust), "Hemispheres");
      declare(InitialQuarks(), "InitialQuarks");
      declare(FastJets(FS, FastJets::JADE, 0.7), "Jets");

      // histograms
      book(_h_T         ,1,1,1);
      book(_h_T_udsc    ,1,1,2);
      book(_h_T_bottom  ,1,1,3);
      book(_h_rho       ,2,1,1);
      book(_h_rho_udsc  ,2,1,2);
      book(_h_rho_bottom,2,1,3);
      book(_h_B_T       ,3,1,1);
      book(_h_B_T_udsc  ,3,1,2);
      book(_h_B_T_bottom,3,1,3);
      book(_h_B_W       ,4,1,1);
      book(_h_B_W_udsc  ,4,1,2);
      book(_h_B_W_bottom,4,1,3);
      book(_h_C         ,5,1,1);
      book(_h_C_udsc    ,5,1,2);
      book(_h_C_bottom  ,5,1,3);
      book(_h_y23       ,6,1,1);
      book(_h_y23_udsc  ,6,1,2);
      book(_h_y23_bottom,6,1,3);
      book(_sumW_udsc, "_sumW_udsc");
      book(_sumW_b, "_sumW_b");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      int flavour = 0;
      const InitialQuarks& iqf = apply<InitialQuarks>(event, "InitialQuarks");
      Particles quarks;
      if ( iqf.particles().size() == 2 ) {
        flavour = iqf.particles().front().abspid();
        quarks  = iqf.particles();
      } else {
        map<int, Particle> quarkmap;
        for (const Particle& p : iqf.particles()) {
          if (quarkmap.find(p.pid()) == quarkmap.end()) quarkmap[p.pid()] = p;
          else if (quarkmap[p.pid()].E() < p.E()) quarkmap[p.pid()] = p;
        }
        double max_energy = 0.;
        for (int i = 1; i <= 5; ++i) {
          double energy = 0.;
          if (quarkmap.find(i) != quarkmap.end())
            energy += quarkmap[ i].E();
          if (quarkmap.find(-i) != quarkmap.end())
            energy += quarkmap[-i].E();
          if (energy > max_energy)
            flavour = i;
        }
        if (quarkmap.find(flavour) != quarkmap.end())
          quarks.push_back(quarkmap[flavour]);
        if (quarkmap.find(-flavour) != quarkmap.end())
          quarks.push_back(quarkmap[-flavour]);
      }

      // Flavour label
      int iflav = (flavour == PID::DQUARK || flavour == PID::UQUARK ||
                   flavour == PID::SQUARK || flavour == PID::CQUARK) ?
        1 : (flavour == PID::BQUARK) ? 5 : 0;

      // Update weight sums
      if (iflav == 1) {
        _sumW_udsc->fill();
      } else if (iflav == 5) {
        _sumW_b->fill();
      }

      // Thrust
      const Thrust& thrust = applyProjection<Thrust>(event, "Thrust");
      if (iflav == 1) {
        _h_T_udsc->fill(thrust.thrust());
      } else if (iflav == 5) {
        _h_T_bottom->fill(thrust.thrust());
      }
      _h_T->fill(thrust.thrust());


      // The hemisphere variables
      const Hemispheres& hemisphere = applyProjection<Hemispheres>(event, "Hemispheres");
      if (iflav == 1) {
        _h_rho_udsc->fill(hemisphere.scaledM2high());
        _h_B_T_udsc->fill(hemisphere.Bsum());
        _h_B_W_udsc->fill(hemisphere.Bmax());
      } else if (iflav == 5) {
        _h_rho_bottom->fill(hemisphere.scaledM2high());
        _h_B_T_bottom->fill(hemisphere.Bsum());
        _h_B_W_bottom->fill(hemisphere.Bmax());
      }
      _h_rho->fill(hemisphere.scaledM2high());
      _h_B_T->fill(hemisphere.Bsum());
      _h_B_W->fill(hemisphere.Bmax());

      const ParisiTensor& parisi = applyProjection<ParisiTensor>(event, "Parisi");
      if (iflav == 1) {
        _h_C_udsc->fill(parisi.C());
      } else if (iflav == 5) {
        _h_C_bottom->fill(parisi.C());
      }
      _h_C->fill(parisi.C());
      // y_23
      const FastJets& durjet = apply<FastJets>(event, "Jets");
      const double y23 = durjet.clusterSeq()->exclusive_ymerge_max(2);
      if (iflav == 1) {
        _h_y23_udsc->fill(y23);
      } else if (iflav == 5) {
        _h_y23_bottom->fill(y23);
      }
      _h_y23->fill(y23);
    }

    void finalize() {
      scale(_h_T_udsc  , 1./ *_sumW_udsc);
      scale(_h_T_bottom, 1./ *_sumW_b   );
      normalize(_h_T);
      scale(_h_rho_udsc  , 1./ *_sumW_udsc);
      scale(_h_rho_bottom, 1./ *_sumW_b   );
      normalize(_h_rho);
      scale(_h_B_T_udsc  , 1./ *_sumW_udsc);
      scale(_h_B_T_bottom, 1./ *_sumW_b   );
      normalize(_h_B_T);
      scale(_h_B_W_udsc  , 1./ *_sumW_udsc);
      scale(_h_B_W_bottom, 1./ *_sumW_b   );
      normalize(_h_B_W);
      scale(_h_C_udsc  , 1./ *_sumW_udsc);
      scale(_h_C_bottom, 1./ *_sumW_b   );
      normalize(_h_C);
      scale(_h_y23_udsc  , 1./ *_sumW_udsc);
      scale(_h_y23_bottom, 1./ *_sumW_b   );
      normalize(_h_y23);
    }
    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_T  ,_h_T_udsc  , _h_T_bottom;
    Histo1DPtr _h_rho,_h_rho_udsc, _h_rho_bottom;
    Histo1DPtr _h_B_T,_h_B_T_udsc, _h_B_T_bottom;
    Histo1DPtr _h_B_W,_h_B_W_udsc, _h_B_W_bottom;
    Histo1DPtr _h_C  ,_h_C_udsc  , _h_C_bottom;
    Histo1DPtr _h_y23,_h_y23_udsc, _h_y23_bottom;
    CounterPtr _sumW_udsc, _sumW_b;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(L3_2008_I825820);

}
