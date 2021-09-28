// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// CMS cross-section and angular correlations in Z boson + b-hadrons events at 7 TeV
  class CMS_2013_I1256943 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2013_I1256943);


    /// Add projections and book histograms
    void init() {
      book(_sumW, "sumW");
      book(_sumW50, "sumW50");
      book(_sumWpT, "sumWpT");

      FinalState fs(Cuts::abseta < 2.4 && Cuts::pT > 20*GeV);
      declare(fs, "FS");

      UnstableParticles ufs(Cuts::abseta < 2 && Cuts::pT > 15*GeV);
      declare(ufs, "UFS");

      Cut zetacut = Cuts::abseta < 2.4;

      ZFinder zfindermu(fs, zetacut, PID::MUON, 81.0*GeV, 101.0*GeV, 0.1, ZFinder::ClusterPhotons::NONE, ZFinder::AddPhotons::YES, 91.2*GeV);
      declare(zfindermu, "ZFinderMu");

      ZFinder zfinderel(fs, zetacut, PID::ELECTRON, 81.0*GeV, 101.0*GeV, 0.1, ZFinder::ClusterPhotons::NONE, ZFinder::AddPhotons::YES, 91.2*GeV);
      declare(zfinderel, "ZFinderEl");


      // Histograms in non-boosted region of Z pT
      book(_h_dR_BB ,1, 1, 1);
      book(_h_dphi_BB ,2, 1, 1);
      book(_h_min_dR_ZB ,3, 1, 1);
      book(_h_A_ZBB ,4, 1, 1);

      // Histograms in boosted region of Z pT (pT > 50 GeV)
      book(_h_dR_BB_boost ,5, 1, 1);
      book(_h_dphi_BB_boost ,6, 1, 1);
      book(_h_min_dR_ZB_boost ,7, 1, 1);
      book(_h_A_ZBB_boost ,8, 1, 1);

      book(_h_min_ZpT ,9,1,1);
    }


    /// Do the analysis
    void analyze(const Event& e) {

      const UnstableParticles& ufs = apply<UnstableFinalState>(e, "UFS");
      const ZFinder& zfindermu = apply<ZFinder>(e, "ZFinderMu");
      const ZFinder& zfinderel = apply<ZFinder>(e, "ZFinderEl");

      // Look for a Z --> mu+ mu- event in the final state
      if (zfindermu.empty() && zfinderel.empty()) vetoEvent;

      const Particles& z = !zfindermu.empty() ? zfindermu.bosons() : zfinderel.bosons();
      const bool is_boosted = ( z[0].pT() > 50*GeV );

      // Loop over the unstable particles
      vector<FourMomentum> Bmom;
      for (const Particle& p : ufs.particles()) {
        const PdgId pid = p.pid();

        // Look for particles with a bottom quark
        if (PID::hasBottom(pid)) {

          bool good_B = false;
          ConstGenParticlePtr pgen = p.genParticle();
          ConstGenVertexPtr vgen = pgen -> end_vertex();

          // Loop over the decay products of each unstable particle, looking for a b-hadron pair
          /// @todo Avoid HepMC API
          for (ConstGenParticlePtr it: HepMCUtils::particles(vgen, Relatives::CHILDREN)){
            // If the particle produced has a bottom quark do not count it and go to the next loop cycle.
            if (!( PID::hasBottom( it->pdg_id() ) ) ) {
              good_B = true;
              continue;
            } else {
              good_B = false;
              break;
            }
          }
          if (good_B ) Bmom.push_back( p.momentum() );
        }
        else continue;
      }

      // If there are more than two B's in the final state veto the event
      if (Bmom.size() != 2 ) vetoEvent;

      // Calculate the observables
      double dphiBB = deltaPhi(Bmom[0], Bmom[1]);
      double dRBB = deltaR(Bmom[0], Bmom[1]);

      const FourMomentum& pZ = z[0].momentum();
      const bool closest_B = ( deltaR(pZ, Bmom[0]) < deltaR(pZ, Bmom[1]) );
      const double mindR_ZB = closest_B ? deltaR(pZ, Bmom[0]) : deltaR(pZ, Bmom[1]);
      const double maxdR_ZB = closest_B ? deltaR(pZ, Bmom[1]) : deltaR(pZ, Bmom[0]);
      const double AZBB = ( maxdR_ZB - mindR_ZB ) / ( maxdR_ZB + mindR_ZB );

      // Fill the histograms in the non-boosted region
      _h_dphi_BB->fill(dphiBB);
      _h_dR_BB->fill(dRBB);
      _h_min_dR_ZB->fill(mindR_ZB);
      _h_A_ZBB->fill(AZBB);
      _sumW->fill();
      _sumWpT->fill();

      // Fill the histograms in the boosted region
      if (is_boosted) {
        _sumW50->fill();
        _h_dphi_BB_boost->fill(dphiBB);
        _h_dR_BB_boost->fill(dRBB);
        _h_min_dR_ZB_boost->fill(mindR_ZB);
        _h_A_ZBB_boost->fill(AZBB);
      }

      // Fill Z pT (cumulative) histogram
      _h_min_ZpT->fill(0);
      if (pZ.pT() > 40*GeV ) {
        _sumWpT->fill();
        _h_min_ZpT->fill(40);
      }
      if (pZ.pT() > 80*GeV ) {
        _sumWpT->fill();
        _h_min_ZpT->fill(80);
      }
      if (pZ.pT() > 120*GeV ) {
        _sumWpT->fill();
        _h_min_ZpT->fill(120);
      }

      Bmom.clear();
    }


    /// Finalize
    void finalize() {

      // Normalize excluding overflow bins (d'oh)
      normalize(_h_dR_BB, 0.7*crossSection()*dbl(*_sumW)/sumOfWeights(), false);  // d01-x01-y01
      normalize(_h_dphi_BB, 0.53*crossSection()*dbl(*_sumW)/sumOfWeights(), false);   // d02-x01-y01
      normalize(_h_min_dR_ZB, 0.84*crossSection()*dbl(*_sumW)/sumOfWeights(), false); // d03-x01-y01
      normalize(_h_A_ZBB, 0.2*crossSection()*dbl(*_sumW)/sumOfWeights(), false);  // d04-x01-y01

      normalize(_h_dR_BB_boost, 0.84*crossSection()*dbl(*_sumW50)/sumOfWeights(), false); // d05-x01-y01
      normalize(_h_dphi_BB_boost, 0.63*crossSection()*dbl(*_sumW50)/sumOfWeights(), false);   // d06-x01-y01
      normalize(_h_min_dR_ZB_boost, 1*crossSection()*dbl(*_sumW50)/sumOfWeights(), false);    // d07-x01-y01
      normalize(_h_A_ZBB_boost, 0.25*crossSection()*dbl(*_sumW50)/sumOfWeights(), false); // d08-x01-y01

      normalize(_h_min_ZpT, 40*crossSection()*dbl(*_sumWpT)/sumOfWeights(), false);   // d09-x01-y01
    }


  private:

    /// @name Weight counters
    //@{
    CounterPtr _sumW, _sumW50, _sumWpT;
    //@}

    /// @name Histograms
    //@{
    Histo1DPtr _h_dphi_BB, _h_dR_BB, _h_min_dR_ZB, _h_A_ZBB;
    Histo1DPtr _h_dphi_BB_boost, _h_dR_BB_boost, _h_min_dR_ZB_boost, _h_A_ZBB_boost, _h_min_ZpT;
    //@}

  };


  DECLARE_RIVET_PLUGIN(CMS_2013_I1256943);

}
