// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/CentralityProjection.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Tools/AliceCommon.hh"

namespace Rivet {


  /// @brief Strangeness enhancement in pp 7 TeV by ALICE.
  class ALICE_2016_I1471838 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2016_I1471838);

    int profileIndex(vector<double> cBins, double c) {
      int index = 100;
      if (c > 0 && c <= cBins[0]) return cBins.size() - 1;
      for (size_t i = 0; i < cBins.size() - 1; ++i) {
        if (c > cBins[i] && c <= cBins[i + 1]) {
	  index = i;
	  break;
	}
      }
      return max(0, int(cBins.size() - index - 2));
    }

    /// Book histograms and initialise projections before the run
    void init() {
      // Centrality projection.
      declareCentrality(ALICE::V0MMultiplicity(),
        "ALICE_2015_PPCentrality","V0M","V0M");
      // Central primary particles
      declare(ChargedFinalState(Cuts::abseta < 1.0),"PP");
      declare(ALICE::PrimaryParticles(Cuts::absrap < 0.5),"PPy");
      centralityBins = {1.,5.,10.,15.,20., 30., 40., 50., 70., 100.};
      centralityBinsOmega = {5.,15.,30.,50.,100.};
      // Book histograms
      for (int i = 0; i < 10; ++i) {
        book(K0SpT[centralityBins[i]], i+1,1,1);
        book(LambdapT[centralityBins[i]], i+11,1,1);
        book(XipT[centralityBins[i]], i+21,1,1);
        book(sow[centralityBins[i]], "sow_" + toString(i));
      }
      for (int i = 0; i < 5; ++i) {
        book(OmegapT[centralityBinsOmega[i]], i+31,1,1);
        book(sowOmega[centralityBinsOmega[i]], "sowO_" + toString(i));
      }
      book(piYield, 40,1,1);
      book(pYield, 41,1,1);
      book(kYield, 42,1,1);
      book(lambdaYield, 43,1,1);
      book(xiYield, 44,1,1);
      book(omegaYield, 45,1,1);
      book(piRebinned, "/piRebinned", refData(45,1,1));

      // Make the ratios
      book(kpi, 36, 1, 1, true);
      book(ppi, 47, 1, 1, true);
      book(lpi, 37, 1, 1, true);
      book(xpi, 38, 1, 1, true);
      book(opi, 39, 1, 1, true);
      book(lk, 46, 1, 1, true);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (apply<ChargedFinalState>(event,"PP").particles().size() < 1) vetoEvent;
      const ALICE::PrimaryParticles& prim = apply<ALICE::PrimaryParticles>(event,"PPy");
      const CentralityProjection& cent = apply<CentralityProjection>(event,"V0M");
      double c  = cent();
      // Find the correct histograms
      auto kptItr = K0SpT.upper_bound(c);
      if (kptItr == K0SpT.end()) return;
      auto lptItr = LambdapT.upper_bound(c);
      if (lptItr == LambdapT.end()) return;
      auto xptItr = XipT.upper_bound(c);
      if (xptItr == XipT.end()) return;
      auto optItr = OmegapT.upper_bound(c);
      if (optItr == OmegapT.end()) return;
      // Fill the sow.
      auto sowItr = sow.upper_bound(c);
      if (sowItr == sow.end()) return;
      auto sowOmegaItr = sowOmega.upper_bound(c);
      if (sowOmegaItr == sowOmega.end()) return;
      sowItr->second->fill();
      sowOmegaItr->second->fill();
      // Fill the pt histograms and count yields.
      int npi = 0, npr = 0, nk = 0;
      int nla = 0, nxi = 0, nom = 0;
      for (auto p : prim.particles()) {
        const double pT = p.pT();
	const int pid = abs(p.pid());
	if (pid == 211) ++npi;
	else if (pid == 2212) ++npr;
	else if (pid == 310) {
	  kptItr->second->fill(pT);
	  ++nk;
	}
	else if (pid == 3122) {
	  lptItr->second->fill(pT);
	  ++nla;
	}
	else if (pid == 3312) {
	  xptItr->second->fill(pT);
	  ++nxi;
	}
	else if (pid == 3334) {
	  optItr->second->fill(pT);
	  ++nom;
	}
      }
      // Fill the profiles of yields.
      int index = profileIndex(centralityBins,c);
      piYield->fillBin(index, double(npi));
      pYield->fillBin(index, double(npr));
      kYield->fillBin(index, double(nk));
      lambdaYield->fillBin(index, double(nla));
      xiYield->fillBin(index, double(nxi));
      index = profileIndex(centralityBinsOmega, c);
      omegaYield->fillBin(index, double(nom));
      piRebinned->fillBin(index,double(npi));
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // Normalize the spectra
      for (int i = 0; i < 10; ++i) {
        K0SpT[centralityBins[i]]->scaleW(1./sow[centralityBins[i]]->sumW());
        XipT[centralityBins[i]]->scaleW(1./sow[centralityBins[i]]->sumW());
        LambdapT[centralityBins[i]]->scaleW(1./sow[centralityBins[i]]->sumW());
      }
      for (int i = 0; i < 5; ++i) {
        OmegapT[centralityBinsOmega[i]]->scaleW(1./sowOmega[centralityBinsOmega[i]]->sumW());
      }

      divide(kYield, piYield, kpi);
      kpi->scaleY(2.);
      divide(pYield, piYield, ppi);
      divide(lambdaYield, piYield, lpi);
      divide(xiYield, piYield, xpi);
      divide(omegaYield, piRebinned, opi);
      divide(lambdaYield, kYield, lk);
      lk->scaleY(0.5);
    }

    //@}


    /// @name Histograms
    //@{
    // Histograms ordered in centrality classes
    vector<double> centralityBins;
    vector<double> centralityBinsOmega;

    // pT spectra
    map<double, Histo1DPtr> K0SpT;
    map<double, Histo1DPtr> LambdapT;
    map<double, Histo1DPtr> XipT;
    map<double, Histo1DPtr> OmegapT;
    map<double, CounterPtr> sow;
    map<double, CounterPtr> sowOmega;

    // Total yields
    Profile1DPtr piYield;
    Profile1DPtr pYield;
    Profile1DPtr kYield;
    Profile1DPtr lambdaYield;
    Profile1DPtr xiYield;
    Profile1DPtr omegaYield;
    Profile1DPtr piRebinned;

    // Ratios
    Scatter2DPtr kpi;
    Scatter2DPtr ppi;
    Scatter2DPtr lpi;
    Scatter2DPtr xpi;
    Scatter2DPtr opi;
    Scatter2DPtr lk;
    //@}
  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ALICE_2016_I1471838);


}
