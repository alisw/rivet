// -*- C++ -*-
#include "Rivet/Analysis.hh"

#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"


namespace Rivet {


  /// ATLAS Z+jets in pp at 7 TeV
  class ATLAS_2011_I945498 : public Analysis {
  public:

    /// Constructor
    ATLAS_2011_I945498()
      : Analysis("ATLAS_2011_I945498")
    {    }


    /// Book histograms and initialise projections before the run
    void init() {

      // Variable initialisation
      _isZeeSample = false;
      _isZmmSample = false;
      for (size_t chn = 0; chn < 3; ++chn) {
        book(weights_nj0[chn], "weights_nj0_" + to_str(chn));
        book(weights_nj1[chn], "weights_nj1_" + to_str(chn));
        book(weights_nj2[chn], "weights_nj2_" + to_str(chn));
        book(weights_nj3[chn], "weights_nj3_" + to_str(chn));
        book(weights_nj4[chn], "weights_nj4_" + to_str(chn));
      }

      // Set up projections
	  FinalState fs;
      ZFinder zfinder_mu(fs, Cuts::abseta < 2.4 && Cuts::pT > 20*GeV, PID::MUON, 66*GeV, 116*GeV, 0.1, ZFinder::ClusterPhotons::NODECAY);
      declare(zfinder_mu, "ZFinder_mu");

      Cut cuts = (Cuts::abseta < 1.37 || Cuts::absetaIn(1.52, 2.47)) && Cuts::pT > 20*GeV;

      ZFinder zfinder_el(fs, cuts, PID::ELECTRON, 66*GeV, 116*GeV, 0.1, ZFinder::ClusterPhotons::NODECAY);
      declare(zfinder_el, "ZFinder_el");

	  Cut cuts25_20 = Cuts::abseta < 2.5 && Cuts::pT > 20*GeV;
      // For combined cross-sections (combined phase space + dressed level)
      ZFinder zfinder_comb_mu(fs, cuts25_20, PID::MUON, 66.0*GeV, 116.0*GeV, 0.1, ZFinder::ClusterPhotons::NODECAY);
      declare(zfinder_comb_mu, "ZFinder_comb_mu");
      ZFinder zfinder_comb_el(fs, cuts25_20, PID::ELECTRON, 66.0*GeV, 116.0*GeV, 0.1, ZFinder::ClusterPhotons::NODECAY);
      declare(zfinder_comb_el, "ZFinder_comb_el");

      // Define veto FS in order to prevent Z-decay products entering the jet algorithm
      VetoedFinalState remfs;
      remfs.addVetoOnThisFinalState(zfinder_el);
      remfs.addVetoOnThisFinalState(zfinder_mu);
      VetoedFinalState remfs_comb;
      remfs_comb.addVetoOnThisFinalState(zfinder_comb_el);
      remfs_comb.addVetoOnThisFinalState(zfinder_comb_mu);

      FastJets jets(remfs, FastJets::ANTIKT, 0.4);
      jets.useInvisibles();
      declare(jets, "jets");
      FastJets jets_comb(remfs_comb, FastJets::ANTIKT, 0.4);
      jets_comb.useInvisibles();
      declare(jets_comb, "jets_comb");

      // 0=el, 1=mu, 2=comb
      for (size_t chn = 0; chn < 3; ++chn) {
        book(_h_njet_incl[chn]  ,1, 1, chn+1);
        book(_h_njet_ratio[chn] ,2, 1, chn+1);
        book(_h_ptjet[chn]      ,3, 1, chn+1);
        book(_h_ptlead[chn]     ,4, 1, chn+1);
        book(_h_ptseclead[chn]  ,5, 1, chn+1);
        book(_h_yjet[chn]       ,6, 1, chn+1);
        book(_h_ylead[chn]      ,7, 1, chn+1);
        book(_h_yseclead[chn]   ,8, 1, chn+1);
        book(_h_mass[chn]       ,9, 1, chn+1);
        book(_h_deltay[chn]     ,10, 1, chn+1);
        book(_h_deltaphi[chn]   ,11, 1, chn+1);
        book(_h_deltaR[chn]     ,12, 1, chn+1);
      }
    }


    // Jet selection criteria universal for electron and muon channel
    /// @todo Replace with a Cut passed to jetsByPt
    Jets selectJets(const ZFinder* zf, const FastJets* allJets) {
      const FourMomentum l1 = zf->constituents()[0].momentum();
      const FourMomentum l2 = zf->constituents()[1].momentum();
      Jets jets;
      for (const Jet& jet : allJets->jetsByPt(30*GeV)) {
        const FourMomentum jmom = jet.momentum();
        if (jmom.absrap() < 4.4 &&
            deltaR(l1, jmom) > 0.5  && deltaR(l2, jmom) > 0.5) {
          jets.push_back(jet);
        }
      }
      return jets;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      vector<const ZFinder*> zfs;
      zfs.push_back(& (apply<ZFinder>(event, "ZFinder_el")));
      zfs.push_back(& (apply<ZFinder>(event, "ZFinder_mu")));
      zfs.push_back(& (apply<ZFinder>(event, "ZFinder_comb_el")));
      zfs.push_back(& (apply<ZFinder>(event, "ZFinder_comb_mu")));

      vector<const FastJets*> fjs;
      fjs.push_back(& (apply<FastJets>(event, "jets")));
      fjs.push_back(& (apply<FastJets>(event, "jets_comb")));

      // Determine what kind of MC sample this is
      const bool isZee = (zfs[0]->bosons().size() == 1) || (zfs[2]->bosons().size() == 1);
      const bool isZmm = (zfs[1]->bosons().size() == 1) || (zfs[3]->bosons().size() == 1);
      if (isZee) _isZeeSample = true;
      if (isZmm) _isZmmSample = true;

      // Require exactly one electronic or muonic Z-decay in the event
      bool isZeemm = ( (zfs[0]->bosons().size() == 1 && zfs[1]->bosons().size() != 1) ||
                       (zfs[1]->bosons().size() == 1 && zfs[0]->bosons().size() != 1) );
      bool isZcomb = ( (zfs[2]->bosons().size() == 1 && zfs[3]->bosons().size() != 1) ||
                       (zfs[3]->bosons().size() == 1 && zfs[2]->bosons().size() != 1) );
      if (!isZeemm && !isZcomb) vetoEvent;

      vector<int> zfIDs;
      vector<int> fjIDs;
      if (isZeemm) {
        int chn = zfs[0]->bosons().size() == 1 ? 0 : 1;
        zfIDs.push_back(chn);
        fjIDs.push_back(0);
      }
      if (isZcomb) {
        int chn = zfs[2]->bosons().size() == 1 ? 2 : 3;
        zfIDs.push_back(chn);
        fjIDs.push_back(1);
      }

      for (size_t izf = 0; izf < zfIDs.size(); ++izf) {
        int zfID = zfIDs[izf];
        int fjID = fjIDs[izf];

        int chn = zfID;
        if (zfID == 2 || zfID == 3) chn = 2;

        Jets jets = selectJets(zfs[zfID], fjs[fjID]);

        switch (jets.size()) {
        case 0:
          weights_nj0[chn]->fill();
          break;
        case 1:
          weights_nj0[chn]->fill();
          weights_nj1[chn]->fill();
          break;
        case 2:
          weights_nj0[chn]->fill();
          weights_nj1[chn]->fill();
          weights_nj2[chn]->fill();
          break;
        case 3:
          weights_nj0[chn]->fill();
          weights_nj1[chn]->fill();
          weights_nj2[chn]->fill();
          weights_nj3[chn]->fill();
          break;
        default: // >= 4
          weights_nj0[chn]->fill();
          weights_nj1[chn]->fill();
          weights_nj2[chn]->fill();
          weights_nj3[chn]->fill();
          weights_nj4[chn]->fill();
        }

        // Require at least one jet
        if (jets.empty()) continue;

        // Fill jet multiplicities
        for (size_t ijet = 1; ijet <= jets.size(); ++ijet) {
          _h_njet_incl[chn]->fill(ijet);
        }

        // Loop over selected jets, fill inclusive jet distributions
        for (size_t ijet = 0; ijet < jets.size(); ++ijet) {
          _h_ptjet[chn]->fill(jets[ijet].pT()/GeV);
          _h_yjet [chn]->fill(fabs(jets[ijet].rapidity()));
        }

        // Leading jet histos
        const double ptlead   = jets[0].pT()/GeV;
        const double yabslead = fabs(jets[0].rapidity());
        _h_ptlead[chn]->fill(ptlead);
        _h_ylead [chn]->fill(yabslead);

        if (jets.size() >= 2) {
          // Second jet histos
          const double pt2ndlead   = jets[1].pT()/GeV;
          const double yabs2ndlead = fabs(jets[1].rapidity());
          _h_ptseclead[chn] ->fill(pt2ndlead);
          _h_yseclead [chn] ->fill(yabs2ndlead);

          // Dijet histos
          const double deltaphi = fabs(deltaPhi(jets[1], jets[0]));
          const double deltarap = fabs(jets[0].rapidity() - jets[1].rapidity()) ;
          const double deltar   = fabs(deltaR(jets[0], jets[1], RAPIDITY));
          const double mass     = (jets[0].momentum() + jets[1].momentum()).mass();
          _h_mass    [chn] ->fill(mass/GeV);
          _h_deltay  [chn] ->fill(deltarap);
          _h_deltaphi[chn] ->fill(deltaphi);
          _h_deltaR  [chn] ->fill(deltar);
        }
      }
    }


    /// @name Ratio calculator util functions
    //@{

    /// Calculate the ratio, being careful about div-by-zero
    double ratio(double a, double b) {
      return (b != 0) ? a/b : 0;
    }

    /// Calculate the ratio error, being careful about div-by-zero
    double ratio_err(double a, double b) {
      return (b != 0) ? sqrt(a/sqr(b) + sqr(a)/(b*b*b)) : 0;
    }

    //@}


    void finalize() {
      // Fill ratio histograms
      for (size_t chn = 0; chn < 3; ++chn) {
        _h_njet_ratio[chn]->addPoint(1, ratio(weights_nj1[chn]->val(), weights_nj0[chn]->val()), 0.5, ratio_err(weights_nj1[chn]->val(), weights_nj0[chn]->val()));
        _h_njet_ratio[chn]->addPoint(2, ratio(weights_nj2[chn]->val(), weights_nj1[chn]->val()), 0.5, ratio_err(weights_nj2[chn]->val(), weights_nj1[chn]->val()));
        _h_njet_ratio[chn]->addPoint(3, ratio(weights_nj3[chn]->val(), weights_nj2[chn]->val()), 0.5, ratio_err(weights_nj3[chn]->val(), weights_nj2[chn]->val()));
        _h_njet_ratio[chn]->addPoint(4, ratio(weights_nj4[chn]->val(), weights_nj3[chn]->val()), 0.5, ratio_err(weights_nj4[chn]->val(), weights_nj3[chn]->val()));
      }

      // Scale other histos
      for (size_t chn = 0; chn < 3; ++chn) {
        // For ee and mumu channels: normalize to Njet inclusive cross-section
        double xs = crossSectionPerEvent()/picobarn;
        if (chn != 2 && weights_nj0[chn]->val() != 0.)  xs = 1.0 / weights_nj0[chn]->val();
        // For inclusive MC sample(ee/mmu channels together) we want the single-lepton-flavor xsec
        if (_isZeeSample && _isZmmSample) xs *= 0.5;

        // Special case histogram: always not normalized
        scale(_h_njet_incl[chn], (chn < 2) ? crossSectionPerEvent()/picobarn : xs);

        scale(_h_ptjet[chn]    , xs);
        scale(_h_ptlead[chn]   , xs);
        scale(_h_ptseclead[chn], xs);
        scale(_h_yjet[chn]     , xs);
        scale(_h_ylead[chn]    , xs);
        scale(_h_yseclead[chn] , xs);
        scale(_h_deltaphi[chn] , xs);
        scale(_h_deltay[chn]   , xs);
        scale(_h_deltaR[chn]   , xs);
        scale(_h_mass[chn]     , xs);
      }

    }

    //@}


  private:

    bool _isZeeSample;
    bool _isZmmSample;

    CounterPtr weights_nj0[3];
    CounterPtr weights_nj1[3];
    CounterPtr weights_nj2[3];
    CounterPtr weights_nj3[3];
    CounterPtr weights_nj4[3];

    Scatter2DPtr _h_njet_ratio[3];
    Histo1DPtr _h_njet_incl[3];
    Histo1DPtr _h_ptjet[3];
    Histo1DPtr _h_ptlead[3];
    Histo1DPtr _h_ptseclead[3];
    Histo1DPtr _h_yjet[3];
    Histo1DPtr _h_ylead[3];
    Histo1DPtr _h_yseclead[3];
    Histo1DPtr _h_deltaphi[3];
    Histo1DPtr _h_deltay[3];
    Histo1DPtr _h_deltaR[3];
    Histo1DPtr _h_mass[3];

  };


  DECLARE_RIVET_PLUGIN(ATLAS_2011_I945498);


}
