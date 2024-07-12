// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/CentralityProjection.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Tools/Cuts.hh"

namespace Rivet {


  /// @brief Identified particles in p--Pb @ 5 TeV
  class ALICE_2014_I1244523 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2014_I1244523);


    /// @name Analysis methods
    //@{

    int profileIndex(vector<double> cBins, double c) {
      int index = 100;
      if (c > 0 && c <= cBins[0]) return cBins.size() - 1;
      for (size_t i = 0; i < cBins.size() - 1; ++i) {
        if (c > cBins[i] && c <= cBins[i + 1]) {
	  index = i;
	  break;
	}
      }
      // Catch low fluctuation.
      return max(0, int(cBins.size() - index - 2));
    }

    void scaleHisto(Histo1DPtr h) {
      vector<YODA::HistoBin1D>& bins = h->bins();
      for (vector<YODA::HistoBin1D>::iterator b = bins.begin(); b != bins.end(); ++b) {
        b->scaleW(1./b->width()/b->xMid());
      }
    }

    /// Book histograms and initialise projections before the run
    void init() {
      // The centrality projection.
      declareCentrality(ALICE::V0AMultiplicity(),
           "ALICE_2015_PPBCentrality", "V0A", "V0A");

      // Define the cuts for the analysis:
      // pPb Collision has a centre of mass system shift of +0.465
      // They study -0.5 < yCoM < 0.0 -> -0.035 < y < 0.465
      const Cut& cut = Cuts::rap < 0.035 && Cuts::rap > -0.465;
      //const Cut& cut = Cuts::rap > -0.035 && Cuts::rap < 0.465;
      const ALICE::PrimaryParticles fs(cut);
      declare(fs,"FS");

      // The event trigger.
      declare(ALICE::V0AndTrigger(), "V0-AND");

      // The centrality bins
      centralityBins = {5.,10.,20.,40.,60.,80.,100.};

      for (int i = 0; i < 4; ++i) {
       // First we book the invariant spectra.
        book(_histPipT[centralityBins[i]], 1, 1, 1 + i);
        if (i < 3) book(_histPipT[centralityBins[i + 4]], 2, 1, 1 + i);
        book(_histKpT[centralityBins[i]], 3, 1, 1 + i);
        if (i < 3) book(_histKpT[centralityBins[i + 4]], 4, 1, 1 + i);
        book(_histK0SpT[centralityBins[i]], 5, 1, 1 + i);
        if (i < 3) book(_histK0SpT[centralityBins[i + 4]], 6, 1, 1 + i);
        book(_histProtonpT[centralityBins[i]], 7, 1, 1 + i);
        if (i < 3) book(_histProtonpT[centralityBins[i + 4]], 8, 1, 1 + i);
        book(_histLambdapT[centralityBins[i]], 9, 1, 1 + i);
        if (i < 3) book(_histLambdapT[centralityBins[i + 4]], 10, 1, 1 + i);
        // The associated sow counters.
        book(_sow[centralityBins[i]], "TMP/sow" + toString(i));
        if (i < 3) book(_sow[centralityBins[i + 4]], "TMP/sow" + toString(i + 4));
      	// Then the pi spectra going into the centrality dependent pT ratios.
        book(_tmpPi4KpT[centralityBins[i]], "TMP/NPi4K" + toString(i), refData(11, 1, 1 + i));
        if (i < 3) book(_tmpPi4KpT[centralityBins[i + 4]], "TMP/NPi4K" + toString(i + 4), refData(12, 1, 1 + i));
        book(_tmpPi4PpT[centralityBins[i]], "TMP/NPi4P" + toString(i), refData(13, 1, 1 + i));
        if (i < 3) book(_tmpPi4PpT[centralityBins[i + 4]], "TMP/NPi4P" + toString(i + 4), refData(14, 1, 1 + i));
        book(_tmpK4LpT[centralityBins[i]], "TMP/NK4L" + toString(i), refData(15, 1, 1 + i));
        if (i < 3) book(_tmpK4LpT[centralityBins[i + 4]], "TMP/NK4L" + toString(i + 4), refData(16, 1, 1 + i));
	// Then the rest of the spectra going into the cent. dep't pT ratios.
        book(_tmpKpT[centralityBins[i]], "TMP/NK" + toString(i), refData(11, 1, 1 + i));
        if (i < 3) book(_tmpKpT[centralityBins[i + 4]], "TMP/NK" + toString(i + 4), refData(12, 1, 1 + i));
        book(_tmpProtonpT[centralityBins[i]], "TMP/NP" + toString(i), refData(13, 1, 1 + i));
        if (i < 3) book(_tmpProtonpT[centralityBins[i + 4]], "TMP/NP" + toString(i + 4), refData(14, 1, 1 + i));
        book(_tmpLambdapT[centralityBins[i]], "TMP/NL" + toString(i), refData(15, 1, 1 + i));
        if (i < 3) book(_tmpLambdapT[centralityBins[i + 4]], "TMP/NL" + toString(i + 4), refData(16, 1, 1 + i));
        // Then the centrality dependent pT ratios.
        book(_ratioKPi[centralityBins[i]], 11, 1, 1 + i, true);
        if (i < 3) book(_ratioKPi[centralityBins[i + 4]], 12, 1, 1 + i, true);
        book(_ratioPPi[centralityBins[i]], 13, 1, 1 + i, true);
        if (i < 3) book(_ratioPPi[centralityBins[i + 4]], 14, 1, 1 + i, true);
        book(_ratioLK[centralityBins[i]], 15, 1, 1 + i, true);
        if (i < 3) book(_ratioLK[centralityBins[i + 4]], 16, 1, 1 + i, true);
      }

      // Mean pT vs. multiplicity class.
      book(_histLambdaMeanpT, 17, 1, 1);
      book(_histProtonMeanpT, 18, 1, 1);
      book(_histK0SMeanpT,    19, 1, 1);
      book(_histKMeanpT,      20, 1, 1);
      book(_histPiMeanpT,     21, 1, 1);

      // Yield ratios.
      book(_histKtoPiYield,      22, 1, 1, true);
      book(_histProtontoPiYield, 22, 1, 2, true);
      book(_histLambdatoPiYield, 22, 1, 3, true);

      book(_histKYield,      "TMP/KY", refData(22,1,1));
      book(_histProtonYield, "TMP/PrY",refData(22,1,2));
      book(_histLambdaYield, "TMP/LY", refData(22,1,3));
      book(_histPiYield,     "TMP/PiY",refData(22,1,1));
      book(_histPi4LYield,   "TMP/PiLY",refData(22,1,3)); // HepData entry is wrong -- look in the paper.

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Event trigger.
      if (!apply<ALICE::V0AndTrigger>(event, "V0-AND")() ) vetoEvent;
      // Centrality
      const CentralityProjection& cent = apply<CentralityProjection>(event,"V0A");
      double c = cent();
      // Find the index for the profiles.
      int index = profileIndex(centralityBins, c);
      // Find the correct histograms
      // all the pion histos
      auto pi1Itr = _histPipT.upper_bound(c);
      // Test the first one.
      if (pi1Itr == _histPipT.end()) return;
      auto pi2Itr = _tmpPi4KpT.upper_bound(c);
      auto pi3Itr = _tmpPi4PpT.upper_bound(c);
      // Then the rest
      auto kItr = _histKpT.upper_bound(c);
      auto k0Itr = _histK0SpT.upper_bound(c);
      auto krItr = _tmpKpT.upper_bound(c);
      auto klItr = _tmpK4LpT.upper_bound(c);
      auto pItr = _histProtonpT.upper_bound(c);
      auto prItr = _tmpProtonpT.upper_bound(c);
      auto lItr = _histLambdapT.upper_bound(c);
      auto lrItr = _tmpLambdapT.upper_bound(c);
      // And the sow
      auto sowItr = _sow.upper_bound(c);
      sowItr->second->fill();


      const ALICE::PrimaryParticles& fs =
        apply<ALICE::PrimaryParticles>(event,"FS");
      // Count number of particles for yields.
      int npi = 0, nk = 0, np = 0, nlam = 0;
      for(auto p : fs.particles()) {
	  const double pT = p.pT();
	  const int pid = abs(p.pid());
	  const double nW = 1 / M_PI / pT; // Dividing and multiplying by 2 because dy.
	  if (pid == 211) { // pi+/-
	    ++npi;
	    pi1Itr->second->fill(pT, nW);
	    pi2Itr->second->fill(pT);
	    pi3Itr->second->fill(pT);
	    _histPiMeanpT->fillBin(index, pT);
	  }
	  else if (pid == 321) { // K +/-
	    ++nk;
	    kItr->second->fill(pT, nW);
	    krItr->second->fill(pT);
	    _histKMeanpT->fillBin(index, pT);
	  }
	  else if (pid == 310) { // K0S
	    k0Itr->second->fill(pT, nW);
	    klItr->second->fill(pT);
	    _histK0SMeanpT->fillBin(index, pT);
	  }
	  else if (pid == 2212) { // p + pbar
	    ++np;
	    pItr->second->fill(pT, nW);
	    prItr->second->fill(pT);
	    _histProtonMeanpT->fillBin(index, pT);
	  }
	  else if (pid == 3122) { // Lambda + Lambdabar
	    ++nlam;
	    lItr->second->fill(pT, nW);
	    lrItr->second->fill(pT);
	    _histLambdaMeanpT->fillBin(index, pT);
	  }
        }
      // Fill the yield profiles.
      _histKYield->fillBin(index, double(nk));
      _histPi4LYield->fillBin(index, double(npi));
      _histProtonYield->fillBin(index, double(np));
      _histPiYield->fillBin(index, double(npi));
      _histLambdaYield->fillBin(index, double(nlam));
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      // Loop over centrality classes.
      for (int i = 0; i < 7; i++){

         // Normalize the spectra.
        _histPipT[centralityBins[i]]->scaleW(1./_sow[centralityBins[i]]->sumW());
        _histKpT[centralityBins[i]]->scaleW(1./_sow[centralityBins[i]]->sumW());
        _histK0SpT[centralityBins[i]]->scaleW(1./_sow[centralityBins[i]]->sumW());
        _histProtonpT[centralityBins[i]]->scaleW(1./_sow[centralityBins[i]]->sumW());
        _histLambdapT[centralityBins[i]]->scaleW(1./_sow[centralityBins[i]]->sumW());

	// Make the pT ratios.
        divide(_tmpKpT[centralityBins[i]], _tmpPi4KpT[centralityBins[i]],
	  _ratioKPi[centralityBins[i]]);
        divide(_tmpProtonpT[centralityBins[i]], _tmpPi4PpT[centralityBins[i]],
	  _ratioPPi[centralityBins[i]]);
        divide(_tmpLambdapT[centralityBins[i]], _tmpK4LpT[centralityBins[i]],
	  _ratioLK[centralityBins[i]]);
      }

      divide(_histKYield,      _histPiYield,  _histKtoPiYield);
      divide(_histProtonYield, _histPiYield,  _histProtontoPiYield);
      divide(_histLambdaYield, _histPi4LYield, _histLambdatoPiYield);

    }

    //@}

private:
    vector<double> centralityBins;
    // pT spectra (separated by multiplicity classes)
    map<double, Histo1DPtr> _histPipT;
    map<double, Histo1DPtr> _histKpT;
    map<double, Histo1DPtr> _histK0SpT;
    map<double, Histo1DPtr> _histProtonpT;
    map<double, Histo1DPtr> _histLambdapT;

    // Associated sum of weights.
    map<double, CounterPtr> _sow;

    // pT spectra for ratios.
    map<double, Histo1DPtr> _tmpPi4KpT;
    map<double, Histo1DPtr> _tmpPi4PpT;
    map<double, Histo1DPtr> _tmpK4LpT;
    map<double, Histo1DPtr> _tmpKpT;
    map<double, Histo1DPtr> _tmpProtonpT;
    map<double, Histo1DPtr> _tmpLambdapT;

    // The acual ratios.
    map<double, Scatter2DPtr> _ratioKPi;
    map<double, Scatter2DPtr> _ratioPPi;
    map<double, Scatter2DPtr> _ratioLK;

    // Mean pT vs. Multiplicity
    Profile1DPtr       _histKMeanpT;
    Profile1DPtr       _histK0SMeanpT;
    Profile1DPtr       _histProtonMeanpT;
    Profile1DPtr       _histLambdaMeanpT;
    Profile1DPtr       _histPiMeanpT;

    // Total yields
    Profile1DPtr        _histKYield;
    Profile1DPtr        _histProtonYield;
    Profile1DPtr        _histLambdaYield;
    Profile1DPtr        _histPiYield;
    Profile1DPtr        _histPi4LYield;

    // Yield ratios.
    Scatter2DPtr       _histKtoPiYield;
    Scatter2DPtr       _histProtontoPiYield;
    Scatter2DPtr       _histLambdatoPiYield;

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ALICE_2014_I1244523);


}
