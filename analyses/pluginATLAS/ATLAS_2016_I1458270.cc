// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/SmearedParticles.hh"
#include "Rivet/Projections/SmearedJets.hh"
#include "Rivet/Projections/SmearedMET.hh"
#include "Rivet/Tools/Cutflow.hh"

namespace Rivet {


  /// @brief ATLAS 0-lepton SUSY search with 3.2/fb of 13 TeV pp data
  class ATLAS_2016_I1458270 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2016_I1458270);


    /// @name Analysis methods
    //@{

    // method to turn Hist1D into Scatter... so we can write this out witout dividing by bin width
    // since the HEPData entry corresponding to this does not divide the refData by bin width!
    // Have requested they update their HEPData entry but until they do so, we use this workaround.
    Scatter2DPtr convertToScatterWithoutBinWidthDivision(Histo1DPtr input , Scatter2DPtr output ){
      for (size_t b = 0; b < input->numBins(); ++b) {
            const double x   = input->bin(b).xMid();
            const double ex  = input->bin(b).xWidth()/2.;
            const double val = input->bin(b).sumW();
            const double err = input->bin(b).relErr() * val;
            output->addPoint(x, val, ex, err);
       }
       return output;
    }

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      FinalState calofs(Cuts::abseta < 4.8);
      FastJets fj(calofs, FastJets::ANTIKT, 0.4);
      declare(fj, "TruthJets");
      declare(SmearedJets(fj, JET_SMEAR_ATLAS_RUN2, JET_BTAG_ATLAS_RUN2_MV2C20), "RecoJets");

      MissingMomentum mm(calofs);
      declare(mm, "TruthMET");
      declare(SmearedMET(mm, MET_SMEAR_ATLAS_RUN2), "RecoMET");

      PromptFinalState es(Cuts::abseta < 2.47 && Cuts::abspid == PID::ELECTRON, true, true);
      declare(es, "TruthElectrons");
      declare(SmearedParticles(es, ELECTRON_RECOEFF_ATLAS_RUN2, ELECTRON_SMEAR_ATLAS_RUN2), "RecoElectrons");

      PromptFinalState mus(Cuts::abseta < 2.7 && Cuts::abspid == PID::MUON, true);
      declare(mus, "TruthMuons");
      declare(SmearedParticles(mus, MUON_EFF_ATLAS_RUN2, MUON_SMEAR_ATLAS_RUN2), "RecoMuons");


      // Book histograms/counters
      book(_h_2jl, "2jl");
      book(_h_2jm, "2jm");
      book(_h_2jt, "2jt");
      book(_h_4jt, "4jt");
      book(_h_5j , "5j");
      book(_h_6jm, "6jm");
      book(_h_6jt, "6jt");

      book(_hMeff_2jl, 4,1,1);
      book(_hMeff_2jm, 5,1,1);
      book(_hMeff_2jt, 6,1,1);
      book(_hMeff_4jt, 7,1,1);
      book(_hMeff_5j , 8,1,1);
      book(_hMeff_6jm, 9,1,1);
      book(_hMeff_6jt, 10,1,1);

      book(_h_temp_Meff_2jl, "_temp_Meff_2jl",refData( 4,1,1));
      book(_h_temp_Meff_2jm, "_temp_Meff_2jm",refData( 5,1,1));
      book(_h_temp_Meff_2jt, "_temp_Meff_2jt",refData( 6,1,1));
      book(_h_temp_Meff_4jt, "_temp_Meff_4jt",refData( 7,1,1));
      book(_h_temp_Meff_5j , "_temp_Meff_5j" ,refData( 8,1,1));
      book(_h_temp_Meff_6jm, "_temp_Meff_6jm",refData( 9,1,1));
      book(_h_temp_Meff_6jt, "_temp_Meff_6jt",refData(10,1,1));




      // Book cut-flows
      const vector<string> cuts2j = {"Pre-sel+MET+pT1", "Njet", "Dphi_min(j,MET)", "pT2", "MET/sqrtHT", "m_eff(incl)"};
      _flows.addCutflow("2jl", cuts2j);
      _flows.addCutflow("2jm", cuts2j);
      _flows.addCutflow("2jt", cuts2j);
      const vector<string> cutsXj = {"Pre-sel+MET+pT1", "Njet", "Dphi_min(j,MET)", "pT2", "pT4", "Aplanarity", "MET/m_eff(Nj)", "m_eff(incl)"};
      _flows.addCutflow("4jt", cutsXj);
      _flows.addCutflow("5j",  cutsXj);
      _flows.addCutflow("6jm", cutsXj);
      _flows.addCutflow("6jt", cutsXj);

    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      _flows.fillinit();

      // Same MET cut for all signal regions
      //const Vector3 vmet = -apply<MissingMomentum>(event, "TruthMET").vectorEt();
      const Vector3 vmet = -apply<SmearedMET>(event, "RecoMET").vectorEt();
      const double met = vmet.mod();
      if (met < 200*GeV) vetoEvent;

      // Get baseline electrons, muons, and jets
      Particles elecs = apply<ParticleFinder>(event, "RecoElectrons").particles(Cuts::pT > 10*GeV);
      Particles muons = apply<ParticleFinder>(event, "RecoMuons").particles(Cuts::pT > 10*GeV);
      Jets jets = apply<JetAlg>(event, "RecoJets").jetsByPt(Cuts::pT > 20*GeV && Cuts::abseta < 2.8); ///< @todo Pile-up subtraction

      // Jet/electron/muons overlap removal and selection
      // Remove any |eta| < 2.8 jet within dR = 0.2 of a baseline electron
      for (const Particle& e : elecs)
        ifilter_discard(jets, deltaRLess(e, 0.2, RAPIDITY));
      // Remove any electron or muon with dR < 0.4 of a remaining (Nch > 3) jet
      for (const Jet& j : jets) {
        /// @todo Add track efficiency random filtering
        ifilter_discard(elecs, deltaRLess(j, 0.4, RAPIDITY));
        if (j.particles(Cuts::abscharge > 0 && Cuts::pT > 500*MeV).size() >= 3)
          ifilter_discard(muons, deltaRLess(j, 0.4, RAPIDITY));
      }
      // Discard the softer of any electrons within dR < 0.05
      for (size_t i = 0; i < elecs.size(); ++i) {
        const Particle& e1 = elecs[i];
        /// @todo Would be nice to pass a "tail view" for the filtering, but awkward without range API / iterator guts
        ifilter_discard(elecs, [&](const Particle& e2){ return e2.pT() < e1.pT() && deltaR(e1,e2) < 0.05; });
      }

      // Loose electron selection
      ifilter_select(elecs, ParticleEffFilter(ELECTRON_EFF_ATLAS_RUN2_LOOSE));

      // Veto the event if there are any remaining baseline leptons
      if (!elecs.empty()) vetoEvent;
      if (!muons.empty()) vetoEvent;

      // Signal jets have pT > 50 GeV
      const Jets jets50 = filter_select(jets, Cuts::pT > 50*GeV);
      if (jets50.size() < 2) vetoEvent;
      vector<double> jetpts; transform(jets, jetpts, pT);
      vector<double> jetpts50; transform(jets50, jetpts50, pT);
      const double j1pt = jetpts50[0];
      const double j2pt = jetpts50[1];
      if (j1pt < 200*GeV) vetoEvent;

      // Construct multi-jet observables
      const double ht = sum(jetpts, 0.0);
      const double met_sqrt_ht = met / sqrt(ht);
      const double meff_incl = sum(jetpts50, met);

      // Get dphis between MET and jets
      vector<double> dphimets50; transform(jets50, dphimets50, deltaPhiWRT(vmet));
      const double min_dphi_met_3 = min(head(dphimets50, 3));
      MSG_DEBUG(dphimets50 << ", " << min_dphi_met_3);

      // Jet aplanarity
      Sphericity sph; sph.calc(jets);
      const double aplanarity = sph.aplanarity();


      // Fill SR counters
      // 2-jet SRs
      if (_flows["2jl"].filltail({true, true, min_dphi_met_3 > 0.8, j2pt > 200*GeV,
            met_sqrt_ht > 15*sqrt(GeV), meff_incl > 1200*GeV})) _h_2jl->fill();
      if (_flows["2jm"].filltail({j1pt > 300*GeV, true, min_dphi_met_3 > 0.4, j2pt > 50*GeV,
            met_sqrt_ht > 15*sqrt(GeV), meff_incl > 1600*GeV})) _h_2jm->fill();
      if (_flows["2jt"].filltail({true, true, min_dphi_met_3 > 0.8, j2pt > 200*GeV,
            met_sqrt_ht > 20*sqrt(GeV), meff_incl > 2000*GeV})) _h_2jt->fill();

      // Fill SR Meff Histo1Ds
      // 2-jet SRs
      if ((min_dphi_met_3 > 0.8) && (j2pt > 200*GeV) && (met_sqrt_ht > 15*sqrt(GeV))) _h_temp_Meff_2jl->fill(meff_incl);
      if ((j1pt > 300*GeV)  && (min_dphi_met_3 > 0.4) && (j2pt > 50*GeV) && (met_sqrt_ht > 15*sqrt(GeV))) _h_temp_Meff_2jm->fill(meff_incl);
      if ((min_dphi_met_3 > 0.8) && (j2pt > 200*GeV ) && (met_sqrt_ht > 20*sqrt(GeV))) _h_temp_Meff_2jt->fill(meff_incl);

      // Upper multiplicity SRs
      const double j4pt = jets50.size() > 3 ? jetpts50[3] : -1;
      const double j5pt = jets50.size() > 4 ? jetpts50[4] : -1;
      const double j6pt = jets50.size() > 5 ? jetpts50[5] : -1;
      const double meff_4 = jets50.size() > 3 ? sum(head(jetpts50, 4), met) : -1;
      const double meff_5 = jets50.size() > 4 ? meff_4 + jetpts50[4] : -1;
      const double meff_6 = jets50.size() > 5 ? meff_5 + jetpts50[5] : -1;
      const double met_meff_4 = met / meff_4;
      const double met_meff_5 = met / meff_5;
      const double met_meff_6 = met / meff_6;
      const double min_dphi_met_more = jets50.size() > 3 ? min(tail(dphimets50, -3)) : -1;


      if (_flows["4jt"].filltail({true, jets50.size() >= 4, min_dphi_met_3 > 0.4 && min_dphi_met_more > 0.2,
            jetpts[1] > 100*GeV, j4pt > 100*GeV, aplanarity > 0.04, met_meff_4 > 0.20, meff_incl > 2200*GeV}))
        _h_4jt->fill();
      if (_flows["5j"].filltail({true, jets50.size() >= 5, min_dphi_met_3 > 0.4 && min_dphi_met_more > 0.2,
            jetpts[1] > 100*GeV, j4pt > 100*GeV && j5pt > 50*GeV, aplanarity > 0.04, met_meff_5 > 0.25, meff_incl > 1600*GeV}))
        _h_5j->fill();
      if (_flows["6jm"].filltail({true, jets50.size() >= 6, min_dphi_met_3 > 0.4 && min_dphi_met_more > 0.2,
            jetpts[1] > 100*GeV, j4pt > 100*GeV && j6pt > 50*GeV, aplanarity > 0.04, met_meff_6 > 0.25, meff_incl > 1600*GeV}))
        _h_6jm->fill();
      if (_flows["6jt"].filltail({true, jets50.size() >= 6, min_dphi_met_3 > 0.4 && min_dphi_met_more > 0.2,
            jetpts[1] > 100*GeV, j4pt > 100*GeV && j6pt > 50*GeV, aplanarity > 0.04, met_meff_6 > 0.20, meff_incl > 2000*GeV}))
        _h_6jt->fill();

      // Fill SR Meff Histo1Ds
      // Upper multiplicity SRs
      if (((jets50.size() >= 4) && (min_dphi_met_3 > 0.4) && (min_dphi_met_more > 0.2) &&  (jetpts[1] > 100*GeV) && (j4pt > 100*GeV) && (aplanarity > 0.04) && (met_meff_4 > 0.20))) _h_temp_Meff_4jt->fill(meff_incl);
      if (((jets50.size() >= 5) && (min_dphi_met_3 > 0.4) && (min_dphi_met_more > 0.2 ) &&
            (jetpts[1] > 100*GeV) && (j4pt > 100*GeV) && (j5pt > 50*GeV) && (aplanarity > 0.04) && (met_meff_5 > 0.25))) _h_temp_Meff_5j->fill(meff_incl);
      if (((jets50.size() >= 6) && (min_dphi_met_3 > 0.4) && (min_dphi_met_more > 0.2) &&
            (jetpts[1] > 100*GeV) && (j4pt > 100*GeV) && (j6pt > 50*GeV) && (aplanarity > 0.04) && (met_meff_6 > 0.25))) _h_temp_Meff_6jm->fill(meff_incl);
      if (((jets50.size() >= 6) && (min_dphi_met_3 > 0.4) && (min_dphi_met_more > 0.2) &&
            (jetpts[1] > 100*GeV) && (j4pt > 100*GeV) && (j6pt > 50*GeV) && (aplanarity > 0.04) && (met_meff_6 > 0.20))) _h_temp_Meff_6jt->fill(meff_incl);

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double sf = 3.2*crossSection()/femtobarn/sumOfWeights();
      scale(_h_2jl, sf); scale(_h_2jm, sf); scale(_h_2jt, sf);
      scale(_h_4jt, sf); scale(_h_5j, sf);
      scale(_h_6jm, sf); scale(_h_6jt, sf);

      scale(_h_temp_Meff_2jl, sf); scale(_h_temp_Meff_2jm, sf); scale(_h_temp_Meff_2jt, sf);
      scale(_h_temp_Meff_4jt, sf); scale(_h_temp_Meff_5j, sf);
      scale(_h_temp_Meff_6jm, sf); scale(_h_temp_Meff_6jt, sf);


      // the HEPData entry corresponding to this does not divide their distributions
      // by bin width... so to avoid this we need to convert to Scatter2D which is not divided by bw
      _hMeff_2jl = convertToScatterWithoutBinWidthDivision(_h_temp_Meff_2jl,_hMeff_2jl);
      _hMeff_2jm = convertToScatterWithoutBinWidthDivision(_h_temp_Meff_2jm,_hMeff_2jm);
      _hMeff_2jt = convertToScatterWithoutBinWidthDivision(_h_temp_Meff_2jt,_hMeff_2jt);
      _hMeff_4jt = convertToScatterWithoutBinWidthDivision(_h_temp_Meff_4jt,_hMeff_4jt);
      _hMeff_5j  = convertToScatterWithoutBinWidthDivision(_h_temp_Meff_5j ,_hMeff_5j ) ;
      _hMeff_6jm = convertToScatterWithoutBinWidthDivision(_h_temp_Meff_6jm,_hMeff_6jm);
      _hMeff_6jt = convertToScatterWithoutBinWidthDivision(_h_temp_Meff_6jt,_hMeff_6jt);
      MSG_INFO("CUTFLOWS:\n\n" << _flows);

    }

    //@}


    private:

    /// @name Histograms
    //@{
    CounterPtr _h_2jl, _h_2jm, _h_2jt;
    CounterPtr _h_4jt, _h_5j;
    CounterPtr _h_6jm, _h_6jt;

    Scatter2DPtr  _hMeff_2jl, _hMeff_2jm, _hMeff_2jt;
    Scatter2DPtr  _hMeff_4jt, _hMeff_5j;
    Scatter2DPtr  _hMeff_6jm, _hMeff_6jt;

    Histo1DPtr _h_temp_Meff_2jl, _h_temp_Meff_2jm, _h_temp_Meff_2jt;
    Histo1DPtr _h_temp_Meff_4jt, _h_temp_Meff_5j;
    Histo1DPtr _h_temp_Meff_6jm, _h_temp_Meff_6jt;
    //@}

    /// Cut-flows
    Cutflows _flows;

  };



  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2016_I1458270);


}
