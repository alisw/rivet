// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Particle.hh"


namespace Rivet {

  class LHCB_2010_I867355 : public Analysis {
  public:

    LHCB_2010_I867355() : Analysis("LHCB_2010_I867355")
    {  }

    void init() {

      //@ Results are presented for two different fragmentation functions, LEP and Tevatron. Therefore, we have two sets of histograms.
      book(_h_sigma_vs_eta_lep ,1, 1, 1);
      book(_h_sigma_vs_eta_tvt ,1, 1, 2);
      book(_h_sigma_total_lep  ,2, 1, 1);
      book(_h_sigma_total_tvt  ,2, 1, 2);

    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles bhadrons;
      for(ConstGenParticlePtr p: HepMCUtils::particles(event.genEvent())) {
        if (!( PID::isHadron( p->pdg_id() ) && PID::hasBottom( p->pdg_id() )) ) continue;

        ConstGenVertexPtr dv = p->end_vertex();

        bool hasBdaughter = false;
        if ( PID::isHadron( p->pdg_id() ) && PID::hasBottom( p->pdg_id() )) { // selecting b-hadrons
          if (dv) {
            for (ConstGenParticlePtr pp: HepMCUtils::particles(dv, Relatives::CHILDREN)){
              if (PID::isHadron(pp->pdg_id() ) && PID::hasBottom(pp->pdg_id() )) {
                hasBdaughter = true;
              }
            }
          }
        }
        if (hasBdaughter) continue; // continue if the daughter is another b-hadron

        bhadrons += Particle(*p);
      }

      for (const Particle& particle : bhadrons) {

        // take fabs() to use full statistics and then multiply weight by 0.5 because LHCb is single-sided
        double eta = fabs(particle.eta());

        _h_sigma_vs_eta_lep->fill( eta, 0.5 );
        _h_sigma_vs_eta_tvt->fill( eta, 0.5 );

        _h_sigma_total_lep->fill( eta, 0.5 ); // histogram for full kinematic range
        _h_sigma_total_tvt->fill( eta, 0.5 ); // histogram for full kinematic range

      }

    }


    void finalize() {
      double norm = crossSection()/microbarn/sumOfWeights();
      double binwidth = 4.;  // integrated over full rapidity space from 2 to 6.

      // to get the avergae of b and bbar, we scale with 0.5
      scale(_h_sigma_vs_eta_lep, 0.5*norm);
      scale(_h_sigma_vs_eta_tvt, 0.5*norm);
      scale(_h_sigma_total_lep, 0.5*norm*binwidth);
      scale(_h_sigma_total_tvt, 0.5*norm*binwidth);
    }


  private:

    Histo1DPtr _h_sigma_total_lep;
    Histo1DPtr _h_sigma_total_tvt;
    Histo1DPtr _h_sigma_vs_eta_lep;
    Histo1DPtr _h_sigma_vs_eta_tvt;

  };


  // Hook for the plugin system
  DECLARE_RIVET_PLUGIN(LHCB_2010_I867355);

}
