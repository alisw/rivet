// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Babar charmonium spectra
  /// @author Peter Richardson
  class BABAR_2003_I593379 : public Analysis {
  public:

    BABAR_2003_I593379()
      : Analysis("BABAR_2003_I593379")
    { }


    void analyze(const Event& e) {
      // Find the charmonia
      Particles upsilons;
      // First in unstable final state
      const UnstableParticles& ufs = apply<UnstableParticles>(e, "UFS");
      for (const Particle& p : ufs.particles())
        if (p.pid() == 300553) upsilons.push_back(p);
      // Then in whole event if fails
      if (upsilons.empty()) {
        for(ConstGenParticlePtr p: HepMCUtils::particles(e.genEvent())) {
          if (p->pdg_id() != 300553) continue;
          ConstGenVertexPtr pv = p->production_vertex();
          bool passed = true;
          if (pv) {
            for(ConstGenParticlePtr pp: HepMCUtils::particles(pv, Relatives::PARENTS)){
              if ( p->pdg_id() == pp->pdg_id() ) {
                passed = false;
                break;
              }
            }
          }
          if (passed) upsilons.push_back(Particle(p));
        }
      }

      // Find upsilons
      for (const Particle& p : upsilons) {
        _weightSum->fill();
        // Find the charmonium resonances
        /// @todo Use Rivet::Particles
        vector<ConstGenParticlePtr> allJpsi, primaryJpsi, Psiprime, all_chi_c1, all_chi_c2, primary_chi_c1, primary_chi_c2;
        findDecayProducts(p.genParticle(), allJpsi, primaryJpsi, Psiprime,
                          all_chi_c1, all_chi_c2, primary_chi_c1, primary_chi_c2);
        const LorentzTransform cms_boost = LorentzTransform::mkFrameTransformFromBeta(p.mom().betaVec());
        for (size_t i = 0; i < allJpsi.size(); i++) {
          const double pcm = cms_boost.transform(FourMomentum(allJpsi[i]->momentum())).p();
          _hist_all_Jpsi->fill(pcm);
        }
        _mult_JPsi->fill(10.58, double(allJpsi.size()));
        for (size_t i = 0; i < primaryJpsi.size(); i++) {
          const double pcm = cms_boost.transform(FourMomentum(primaryJpsi[i]->momentum())).p();
          _hist_primary_Jpsi->fill(pcm);
        }
        _mult_JPsi_direct->fill(10.58, double(primaryJpsi.size()));
        for (size_t i=0; i<Psiprime.size(); i++) {
          const double pcm = cms_boost.transform(FourMomentum(Psiprime[i]->momentum())).p();
          _hist_Psi_prime->fill(pcm);
        }
        _mult_Psi2S->fill(10.58, double(Psiprime.size()));
        for (size_t i = 0; i < all_chi_c1.size(); i++) {
          const double pcm = cms_boost.transform(FourMomentum(all_chi_c1[i]->momentum())).p();
          _hist_chi_c1->fill(pcm);
        }
        _mult_chi_c1->fill(10.58, double(all_chi_c1.size()));
        _mult_chi_c1_direct->fill(10.58, double(primary_chi_c1.size()));
        for (size_t i = 0; i < all_chi_c2.size(); i++) {
          const double pcm = cms_boost.transform(FourMomentum(all_chi_c2[i]->momentum())).p();
          _hist_chi_c2->fill(pcm);
        }
        _mult_chi_c2->fill(10.58, double(all_chi_c2.size()));
        _mult_chi_c2_direct->fill(10.58, double(primary_chi_c2.size()));
      }
    } // analyze


    void finalize() {
      scale(_hist_all_Jpsi     , 0.5*0.1 / *_weightSum);
      scale(_hist_chi_c1       , 0.5*0.1 / *_weightSum);
      scale(_hist_chi_c2       , 0.5*0.1 / *_weightSum);
      scale(_hist_Psi_prime    , 0.5*0.1 / *_weightSum);
      scale(_hist_primary_Jpsi , 0.5*0.1 / *_weightSum);
      scale(_mult_JPsi         , 0.5*100. / *_weightSum);
      scale(_mult_JPsi_direct  , 0.5*100. / *_weightSum);
      scale(_mult_chi_c1       , 0.5*100. / *_weightSum);
      scale(_mult_chi_c1_direct, 0.5*100. / *_weightSum);
      scale(_mult_chi_c2       , 0.5*100. / *_weightSum);
      scale(_mult_chi_c2_direct, 0.5*100. / *_weightSum);
      scale(_mult_Psi2S        , 0.5*100. / *_weightSum);
    } // finalize


    void init() {
      declare(UnstableParticles(), "UFS");

      book(_mult_JPsi          ,1, 1, 1);
      book(_mult_JPsi_direct   ,1, 1, 2);
      book(_mult_chi_c1        ,1, 1, 3);
      book(_mult_chi_c1_direct ,1, 1, 4);
      book(_mult_chi_c2        ,1, 1, 5);
      book(_mult_chi_c2_direct ,1, 1, 6);
      book(_mult_Psi2S         ,1, 1, 7);
      book(_hist_all_Jpsi      ,6, 1, 1);
      book(_hist_chi_c1        ,7, 1, 1);
      book(_hist_chi_c2        ,7, 1, 2);
      book(_hist_Psi_prime     ,8, 1, 1);
      book(_hist_primary_Jpsi  ,10, 1, 1);

      book(_weightSum, "TMP/weightSum");
    } // init

  private:

    //@{
    // count of weights
    CounterPtr _weightSum;
    /// Histograms
    Histo1DPtr _hist_all_Jpsi;
    Histo1DPtr _hist_chi_c1;
    Histo1DPtr _hist_chi_c2;
    Histo1DPtr _hist_Psi_prime;
    Histo1DPtr _hist_primary_Jpsi;

    Histo1DPtr _mult_JPsi;
    Histo1DPtr _mult_JPsi_direct;
    Histo1DPtr _mult_chi_c1;
    Histo1DPtr _mult_chi_c1_direct;
    Histo1DPtr _mult_chi_c2;
    Histo1DPtr _mult_chi_c2_direct;
    Histo1DPtr _mult_Psi2S;
    //@}

    void findDecayProducts(ConstGenParticlePtr p,
                           vector<ConstGenParticlePtr>& allJpsi,
                           vector<ConstGenParticlePtr>& primaryJpsi,
                           vector<ConstGenParticlePtr>& Psiprime,
                           vector<ConstGenParticlePtr>& all_chi_c1, vector<ConstGenParticlePtr>& all_chi_c2,
                           vector<ConstGenParticlePtr>& primary_chi_c1, vector<ConstGenParticlePtr>& primary_chi_c2) {
      ConstGenVertexPtr dv = p->end_vertex();
      bool isOnium = false;
      /// @todo Use better looping
      for (ConstGenParticlePtr pp: HepMCUtils::particles(dv, Relatives::PARENTS)){
        int id = pp->pdg_id();
        id = id%1000;
        id -= id%10;
        id /= 10;
        if (id==44) isOnium = true;
      }
      /// @todo Use better looping
      for (ConstGenParticlePtr pp: HepMCUtils::particles(dv, Relatives::CHILDREN)){
        int id = pp->pdg_id();
        if (id==100443) {
          Psiprime.push_back(pp);
        }
        else if (id==20443) {
          all_chi_c1.push_back(pp);
          if (!isOnium) primary_chi_c1.push_back(pp);
        }
        else if (id==445) {
          all_chi_c2.push_back(pp);
          if (!isOnium) primary_chi_c2.push_back(pp);
        }
        else if (id==443) {
          allJpsi.push_back(pp);
          if (!isOnium) primaryJpsi.push_back(pp);
        }
        if (pp->end_vertex()) {
          findDecayProducts(pp, allJpsi, primaryJpsi, Psiprime, all_chi_c1, all_chi_c2, primary_chi_c1, primary_chi_c2);
        }
      }
    }

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BABAR_2003_I593379);

}
