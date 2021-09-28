// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief BaBar pion, kaon and proton production in the continuum
  /// @author Peter Richardson
  class BABAR_2013_I1238276 : public Analysis {
  public:

    BABAR_2013_I1238276()
      : Analysis("BABAR_2013_I1238276")
    { }


    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");

      book(_histPion_no_dec   ,1,1,1);
      book(_histKaon_no_dec   ,1,1,2);
      book(_histProton_no_dec ,1,1,3);
      book(_histPion_dec      ,2,1,1);
      book(_histKaon_dec      ,2,1,2);
      book(_histProton_dec    ,2,1,3);
    }


    void analyze(const Event& e) {
      // Loop through charged FS particles and look for charmed mesons/baryons
      const ChargedFinalState& fs = apply<ChargedFinalState>(e, "FS");

      const Beam beamproj = apply<Beam>(e, "Beams");
      const ParticlePair& beams = beamproj.beams();
      const FourMomentum mom_tot = beams.first.momentum() + beams.second.momentum();
      const LorentzTransform cms_boost = LorentzTransform::mkFrameTransformFromBeta(mom_tot.betaVec());
      MSG_DEBUG("CMS Energy sqrt s = " << beamproj.sqrtS());

      for (const Particle& p : fs.particles()) {
        // check if prompt or not
        ConstGenParticlePtr pmother = p.genParticle();
        ConstGenVertexPtr ivertex = pmother->production_vertex();
        bool prompt = true;
        while (ivertex) {
          vector<ConstGenParticlePtr> inparts = HepMCUtils::particles(ivertex, Relatives::PARENTS);
          int n_inparts = inparts.size();
          if (n_inparts < 1) break;
          pmother = inparts[0]; // first mother particle
          int mother_pid = abs(pmother->pdg_id());
          if (mother_pid==PID::K0S || mother_pid==PID::LAMBDA) {
            prompt = false;
            break;
          }
          else if (mother_pid<6) {
            break;
          }
          ivertex = pmother->production_vertex();
        }

        // momentum in CMS frame
        const double mom = cms_boost.transform(p.momentum()).vector3().mod();
        const int PdgId = p.abspid();
        MSG_DEBUG("pdgID = " << PdgId << " Momentum = " << mom);
        switch (PdgId) {
        case PID::PIPLUS:
          if(prompt) _histPion_no_dec->fill(mom);
          _histPion_dec   ->fill(mom);
          break;
        case PID::KPLUS:
          if(prompt) _histKaon_no_dec->fill(mom);
          _histKaon_dec   ->fill(mom);
          break;
        case PID::PROTON:
          if(prompt) _histProton_no_dec->fill(mom);
          _histProton_dec   ->fill(mom);
        default :
          break;
        }
      }
    }


    void finalize() {
      scale(_histPion_no_dec  ,1./sumOfWeights());
      scale(_histKaon_no_dec  ,1./sumOfWeights());
      scale(_histProton_no_dec,1./sumOfWeights());
      scale(_histPion_dec     ,1./sumOfWeights());
      scale(_histKaon_dec     ,1./sumOfWeights());
      scale(_histProton_dec   ,1./sumOfWeights());
    }


  private:

    //@{
    // Histograms for continuum data (sqrt(s) = 10.52 GeV)
    // no K_S and Lambda decays
    Histo1DPtr _histPion_no_dec;
    Histo1DPtr _histKaon_no_dec;
    Histo1DPtr _histProton_no_dec;
    // including decays
    Histo1DPtr _histPion_dec;
    Histo1DPtr _histKaon_dec;
    Histo1DPtr _histProton_dec;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BABAR_2013_I1238276);

}
