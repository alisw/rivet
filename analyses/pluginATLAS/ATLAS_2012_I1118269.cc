// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Particle.hh"

namespace Rivet {

  class ATLAS_2012_I1118269 : public Analysis {
  public:

    ATLAS_2012_I1118269() : Analysis("ATLAS_2012_I1118269")
    {  }

    void init() {
      book(_h_sigma_vs_pt  ,1, 1, 1);
      book(_h_sigma_vs_eta ,2, 1, 1);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles bhadrons;
      for(ConstGenParticlePtr p: HepMCUtils::particles(event.genEvent())) {

        if (!( PID::isHadron( p->pdg_id() ) && PID::hasBottom( p->pdg_id() )) ) continue;

        ConstGenVertexPtr dv = p->end_vertex();

        /// @todo In future, convert to use built-in 'last B hadron' function
        bool hasBdaughter = false;
        if ( PID::isHadron( p->pdg_id() ) && PID::hasBottom( p->pdg_id() )) { // b-hadron selection
          if (dv) {
            /// @todo particles_out_const_iterator is deprecated in HepMC3
            for(ConstGenParticlePtr pp: HepMCUtils::particles(dv, Relatives::CHILDREN)){
              if ( PID::isHadron( pp->pdg_id() ) && PID::hasBottom( pp->pdg_id()) ) {
                hasBdaughter = true;
              }
            }
          }
        }
        if (hasBdaughter) continue;

        bhadrons += Particle(*p);
      }

      for (const Particle& particle : bhadrons) {
        double eta = particle.eta();
        double pt = particle.pT();

        if (!(inRange(eta, -2.5, 2.5))) continue;
        if (pt < 9.*GeV) continue;

        _h_sigma_vs_pt->fill(pt);
        _h_sigma_vs_eta->fill(fabs(eta));

      }

    }


    void finalize() {
      scale(_h_sigma_vs_pt,  crossSection()/nanobarn/sumOfWeights());
      scale(_h_sigma_vs_eta, crossSection()/microbarn/sumOfWeights());
    }


  private:

    Histo1DPtr _h_sigma_vs_pt;
    Histo1DPtr _h_sigma_vs_eta;

  };


  // Hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1118269);

}
