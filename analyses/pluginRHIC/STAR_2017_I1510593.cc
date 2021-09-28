// -*- C++ -*-
#include <complex>
#include <iostream>
#include <string>
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ImpactParameterProjection.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Tools/Percentile.hh"
#include "Rivet/Tools/RHICCommon.hh"
#include "Rivet/Projections/HepMCHeavyIon.hh"

namespace Rivet {


  /// pT distributions, ratios and production yields of hadrons in STAR
  class STAR_2017_I1510593 : public Analysis {
  public:

    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2017_I1510593);

    string coStr(int i, int j, int k) {
      return "/TMP/d" + toString(i) + "x" + toString(j) + "y" + toString(k);
    }


    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declareCentrality(STAR_BES_Centrality(), "STAR_BES_CALIB", "CMULT", "CMULT");

      // The observed particles.
      declare(ChargedFinalState(Cuts::abseta < 0.5 &&
                                Cuts::absrap < 0.1 && Cuts::pT > 0.2), "CFS");

      // Access the HepMC heavy ion info
      declare(HepMCHeavyIon(), "HepMC");


      // Energy bins
      energies = {7.7, 11.5, 19.6, 27.0, 39.0};
      for (size_t i = 0, N = energies.size(); i < N; ++i) {
        if (fuzzyEquals(sqrtS() / 197. / GeV, energies[i])) enebin = i;
      }

      // Centrality bins
      centralities = {5, 10, 20, 30, 40, 50, 60, 70, 80};

      // Energy bins for Fig. 25
      enebinfig = -1;
      if (fuzzyEquals(sqrtS() / 197. / GeV, energies[0])) enebinfig = 0;
      if (fuzzyEquals(sqrtS() / 197. / GeV, energies[4])) enebinfig = 1;

      // Book all histograms for all energies in order to use re-entrant finalize
      /// @todo Raw arrays would be a *lot* easier to read here (and N_cent is fixed)
      _h_dpT_Piplus = vector<vector<Histo1DPtr> >(energies.size(), vector<Histo1DPtr>(centralities.size()));
      _h_dpT_Pi = vector<vector<Histo1DPtr> >(energies.size(), vector<Histo1DPtr>(centralities.size()));
      _h_dpT_Kaonplus = vector<vector<Histo1DPtr> >(energies.size(), vector<Histo1DPtr>(centralities.size()));
      _h_dpT_Kaon = vector<vector<Histo1DPtr> >(energies.size(), vector<Histo1DPtr>(centralities.size()));
      _h_dpT_Proton = vector<vector<Histo1DPtr> >(energies.size(), vector<Histo1DPtr>(centralities.size()));
      _h_dpT_AntiProton = vector<vector<Histo1DPtr> >(energies.size(), vector<Histo1DPtr>(centralities.size()));
      _wght_PiPlus = vector<vector<CounterPtr> >(energies.size(), vector<CounterPtr>(centralities.size()));
      _wght_Pi = vector<vector<CounterPtr> >(energies.size(), vector<CounterPtr>(centralities.size()));
      _wght_KaonPlus = vector<vector<CounterPtr> >(energies.size(), vector<CounterPtr>(centralities.size()));
      _wght_Kaon = vector<vector<CounterPtr> >(energies.size(), vector<CounterPtr>(centralities.size()));
      _wght_Proton = vector<vector<CounterPtr> >(energies.size(), vector<CounterPtr>(centralities.size()));
      _wght_AntiProton = vector<vector<CounterPtr> >(energies.size(), vector<CounterPtr>(centralities.size()));
      _h_npart_PiPlus = vector<Histo1DPtr>(energies.size());
      _h_npart_PiMinus = vector<Histo1DPtr>(energies.size());
      _h_npart_KaPlus = vector<Histo1DPtr>(energies.size());
      _h_npart_KaMinus = vector<Histo1DPtr>(energies.size());
      _h_npart_Proton = vector<Histo1DPtr>(energies.size());
      _h_npart_AntiProton = vector<Histo1DPtr>(energies.size());
      _wght_npart_PiPlus = vector<CounterPtr>(energies.size());
      _wght_npart_PiMinus = vector<CounterPtr>(energies.size());
      _wght_npart_KaonPlus = vector<CounterPtr>(energies.size());
      _wght_npart_KaonMinus = vector<CounterPtr>(energies.size());
      _wght_npart_Proton = vector<CounterPtr>(energies.size());
      _wght_npart_AntiProton = vector<CounterPtr>(energies.size());
      _h_npart_pT_PiPlus = vector<Profile1DPtr>(energies.size());
      _h_npart_pT_PiMinus = vector<Profile1DPtr>(energies.size());
      _h_npart_pT_KaPlus = vector<Profile1DPtr>(energies.size());
      _h_npart_pT_KaMinus = vector<Profile1DPtr>(energies.size());
      _h_npart_pT_Proton = vector<Profile1DPtr>(energies.size());
      _h_npart_pT_AntiProton = vector<Profile1DPtr>(energies.size());

      _h_npart_Piratio = vector<Profile1DPtr>(energies.size());
      _h_npart_Karatio = vector<Profile1DPtr>(energies.size());
      _h_npart_Pratio = vector<Profile1DPtr>(energies.size());
      _h_npart_KaPi = vector<Profile1DPtr>(energies.size());
      _h_npart_AntiPPi = vector<Profile1DPtr>(energies.size());
      _h_npart_KaPiplus = vector<Profile1DPtr>(energies.size());
      _h_npart_PPiplus = vector<Profile1DPtr>(energies.size());

      for (size_t j = 0, N = energies.size(); j < N; ++j) {
        for (size_t i = 0, M = centralities.size(); i < M; ++i) {
          // Book [energy][centrality] histograms.
          book(_h_dpT_Piplus[j][i], 12 + j, 1, 1 + i);
          book(_h_dpT_Pi[j][i], 12 + j, 2, 1 + i);
          book(_h_dpT_Kaonplus[j][i], 12 + j, 3, 1 + i);
          book(_h_dpT_Kaon[j][i], 12 + j, 4, 1 + i);
          book(_h_dpT_Proton[j][i], 12 + j, 5, 1 + i);
          book(_h_dpT_AntiProton[j][i], 12 + j, 6, 1 + i);
          // Book ditto sum of weights.
          book(_wght_PiPlus[j][i], coStr(12 + j, 1, 1 + i));
          book(_wght_Pi[j][i], coStr(12 + j, 2, 1 + i));
          book(_wght_KaonPlus[j][i], coStr(12 + j, 3, 1 + i));
          book(_wght_Kaon[j][i], coStr(12 + j, 4, 1 + i));
          book(_wght_Proton[j][i], coStr(12 + j, 5, 1 + i));
          book(_wght_AntiProton[j][i], coStr(12 + j, 6, 1 + i));
        }
      }

      /// Booking npart histograms
      for (size_t i = 0, N = energies.size(); i < N; ++i) {
        book(_h_npart_PiPlus[i],17, 1, i+1);
        book(_h_npart_PiMinus[i],17, 2, i+1);
        book(_h_npart_KaPlus[i],17, 3, i+1);
        book(_h_npart_KaMinus[i],17, 4, i+1);
        book(_h_npart_Proton[i],17, 5, i+1);
        book(_h_npart_AntiProton[i],17, 6, i+1);
        // ...and the weights.
        book(_wght_npart_PiPlus[i],coStr(17, 1, i+1));
        book(_wght_npart_PiMinus[i],coStr(17, 2, i+1));
        book(_wght_npart_KaonPlus[i],coStr(17, 3, i+1));
        book(_wght_npart_KaonMinus[i],coStr(17, 4, i+1));
        book(_wght_npart_Proton[i],coStr(17, 5, i+1));
        book(_wght_npart_AntiProton[i],coStr(17, 6, i+1));
        // ... and the profiles.
        book(_h_npart_pT_PiPlus[i], 18, 1, i+1);
        book(_h_npart_pT_PiMinus[i], 18, 2, i+1);
        book(_h_npart_pT_KaPlus[i], 18, 3, i+1);
        book(_h_npart_pT_KaMinus[i], 18, 4, i+1);
        book(_h_npart_pT_Proton[i], 18, 5, i+1);
        book(_h_npart_pT_AntiProton[i], 18, 6, i+1);

        book(_h_npart_Piratio[i], 19, 1, i+1);
        book(_h_npart_Karatio[i], 19, 2, i+1);
        book(_h_npart_Pratio[i], 19, 3, i+1);
        book(_h_npart_KaPi[i], 20, 1, i+1);
        book(_h_npart_AntiPPi[i], 20, 2, i+1);
        book(_h_npart_KaPiplus[i], 20, 3, i+1);
        book(_h_npart_PPiplus[i], 20, 4, i+1);
      }

      book(_h_snn_npart_PiPlus, 21, 1, 1);
      book(_h_snn_npart_PiMinus, 21, 1, 2);
      book(_h_snn_npart_KaPlus, 21, 2, 1);
      book(_h_snn_npart_KaMinus, 21, 2, 2);
      book(_h_snn_npart_Proton, 21, 3, 1);
      book(_h_snn_npart_AntiProton, 21, 3, 2);

      book(_h_snn_mt_PiPlus, 22, 1, 1);
      book(_h_snn_mt_PiMinus, 22, 1, 2);
      book(_h_snn_mt_KaPlus, 22, 2, 1);
      book(_h_snn_mt_KaMinus, 22, 2, 2);
      book(_h_snn_mt_Proton, 22, 3, 1);
      book(_h_snn_mt_AntiProton, 22, 3, 2);

      book(_h_snn_Piratio, 23, 1, 1);
      book(_h_snn_Karatio, 23, 2, 1);
      book(_h_snn_Pratio, 23, 3, 1);
      book(_h_snn_KaPiplus, 24, 1, 1);
      book(_h_snn_KaPiminus, 24, 1, 2);

      _h_yields = vector<Profile1DPtr>(2);
      _h_ratios = vector<Profile1DPtr>(2);
      for (size_t i = 0; i < 2; ++i) {
        book(_h_yields[i], 25, i + 1, 1);
        book(_h_ratios[i], 25, i + 3, 1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& cfs = applyProjection<ChargedFinalState>(event, "CFS");
      // Require at least two charged particles for the analysis to
      // make sense. No further triggers are described in the paper.
      const Particles& particles = cfs.particles();
      nprtcl = particles.size();
      if (nprtcl < 2) return;

      /// Determine the centrality
      const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
      const double c = cent();

      /// Determine the impact parameter
      const HepMCHeavyIon & hi = apply<HepMCHeavyIon>(event, "HepMC");
      const double Npart = hi.Npart_targ();

      /// Determine the centrality bin
      cenbin = (c < 5) ? 0 : c / 10 + 1;

      /// Initializing for each event
      for (size_t i = 0; i < 10; ++i) nparts[i] = 0;
      for (size_t i = 0, N = energies.size(); i < N; ++i) {
        nPi[i] = 0;
        nPiPlus[i] = 0;
        nKaon[i] = 0;
        nKaonPlus[i] = 0;
        nProton[i] = 0;
        nAntiProton[i] = 0;
      }

      /// Loop over all charged particles of the CFS
      for (const Particle& p : cfs.particles()) {
        double pT = p.pT()/GeV;
        double mass = p.mass()/GeV;
        double mTm = sqrt(pT * pT + mass * mass) - mass;
        if (p.absrap() < 0.1) {
          const PdgId id = p.pid();
          switch (id) {
          case 211:
            if (c < 80) {
              _h_dpT_Piplus[enebin][cenbin]->fill(pT, 1. / pT);
              _h_npart_PiPlus[enebin]->fill(Npart, 1. / (0.2 * 0.5 * Npart));
              _h_npart_pT_PiPlus[enebin]->fill(Npart, pT, 5);
            }
            if (c < 5) {
              ++nparts[0];
              _h_snn_npart_PiPlus->fill(energies[enebin], 1.0 / (0.2 * 0.5 * Npart));
              _h_snn_mt_PiPlus->fillBin(enebin, mTm);
            }
            ++nPiPlus[enebin];
            break;
          case -211:
            if (c < 80) {
              _h_dpT_Pi[enebin][cenbin]->fill(pT, 1.0 / pT);
              _h_npart_PiMinus[enebin]->fill(Npart, 1.0 / (0.2 * 0.5 * Npart));
              _h_npart_pT_PiMinus[enebin]->fill(Npart, pT, 5);
            }
            if (c < 5) {
              ++nparts[1];
              _h_snn_npart_PiMinus->fillBin(enebin, 1.0 / (0.2 * 0.5 * Npart));
              _h_snn_mt_PiMinus->fillBin(enebin, mTm);
            }
            ++nPi[enebin];
            break;
          case 321:
            if (c < 80) {
              _h_dpT_Kaonplus[enebin][cenbin]->fill(pT, 1.0 / pT);
              _h_npart_KaPlus[enebin]->fill(Npart, 1.0 /(0.2 * 0.5 * Npart));
              _h_npart_pT_KaPlus[enebin]->fill(Npart, pT, 5);
            }
            if (c < 5) {
              ++nparts[2];
              _h_snn_npart_KaPlus->fillBin(enebin, 1.0 / (0.2 * 0.5 * Npart));
              _h_snn_mt_KaPlus->fillBin(enebin, mTm);
            }
            ++nKaonPlus[enebin];
            break;
          case -321:
            if (c < 80) {
              _h_dpT_Kaon[enebin][cenbin]->fill(pT, 1.0 / pT);
              _h_npart_KaMinus[enebin]->fill(Npart, 1.0 / (0.2 * 0.5 * Npart));
              _h_npart_pT_KaMinus[enebin]->fill(Npart, pT, 5);
            }
            if (c < 5) {
              ++nparts[3];
              _h_snn_npart_KaMinus->fillBin(enebin, 1.0 / (0.2 * 0.5 * Npart));
              _h_snn_mt_KaMinus->fillBin(enebin, mTm);
            }
            ++nKaon[enebin];
            break;
          case 2212:
            if (c < 80) {
              _h_dpT_Proton[enebin][cenbin]->fill(pT, 1.0 / pT);
              _h_npart_Proton[enebin]->fill(Npart, 1.0 /(0.2 * 0.5 * Npart));
              _h_npart_pT_Proton[enebin]->fill(Npart, pT, 5);
            }
            if (c < 5) {
              ++nparts[4];
              _h_snn_npart_Proton->fillBin(enebin, 1.0 / (0.2 * 0.5 * Npart));
              _h_snn_mt_Proton->fillBin(enebin, mTm);
            }
            ++nProton[enebin];
            break;
          case -2212:
            if (c < 80) {
              _h_dpT_AntiProton[enebin][cenbin]->fill(pT, 1.0 / pT);
              _h_npart_AntiProton[enebin]->fill(Npart, 1.0 / (0.2 * 0.5 * Npart));
              _h_npart_pT_AntiProton[enebin]->fill(Npart, pT, 5);
            }
            if (c < 5) {
              ++nparts[5];
              _h_snn_npart_AntiProton->fillBin(enebin, 1.0 / (0.2 * 0.5 * Npart));
              _h_snn_mt_AntiProton->fillBin(enebin, mTm);
            }
            ++nAntiProton[enebin];
            break;
          case 3122:
            if (c < 5) ++nparts[6];
            break;
          case -3122:
            if (c < 5) ++nparts[7];
            break;
          case 3312:
            if (c < 5) ++nparts[8];
            break;
          case -3312:
            if (c < 5) ++nparts[9];
            break;
          }
        }
      }

      /// Particle Ratios
      //"if( > 1e-6)" because "> 0" or "!= 0" can cause errors
      if (nPiPlus[enebin] > 1e-6) {
        _h_npart_Piratio[enebin]->fill(Npart, nPi[enebin] / nPiPlus[enebin], 5);
        _h_npart_KaPiplus[enebin]->fill(Npart, nKaonPlus[enebin] / nPiPlus[enebin], 5);
        _h_npart_PPiplus[enebin]->fill(Npart, nProton[enebin] / nPiPlus[enebin], 5);
      }

      if (nPi[enebin] > 1e-6) {
        _h_npart_KaPi[enebin]->fill(Npart, nKaon[enebin] / nPi[enebin], 5);
        _h_npart_AntiPPi[enebin]->fill(Npart, nAntiProton[enebin] / nPi[enebin], 5);
      }

      if (nKaonPlus[enebin] > 1e-6)
        _h_npart_Karatio[enebin]->fill(Npart, nKaon[enebin] / nKaonPlus[enebin], 5);

      if (nProton[enebin] > 1e-6)
        _h_npart_Pratio[enebin]->fill(Npart, nAntiProton[enebin] / nProton[enebin], 5);

      /// Particle Yields
      if (enebinfig == 0 || enebinfig == 1) {
        for (size_t i = 0; i < 10; i++) {
          if (nparts[i] > 1e-6)
            _h_yields[enebinfig]->fill(i + 1, nparts[i], 5);
        }
        if (nparts[0] > 1e-6)
          _h_ratios[enebinfig]->fill(1, nparts[1] / nparts[0], 5);
        if (nparts[2] > 1e-6)
          _h_ratios[enebinfig]->fill(2, nparts[3] / nparts[2], 5);
        if (nparts[4] > 1e-6)
          _h_ratios[enebinfig]->fill(3, nparts[5] / nparts[4], 5);
        if (nparts[6] > 1e-6)
          _h_ratios[enebinfig]->fill(4, nparts[7] / nparts[6], 5);
        if (nparts[8] > 1e-6)
          _h_ratios[enebinfig]->fill(5, nparts[9] / nparts[8], 5);
        if (nparts[1] > 1e-6) {
          _h_ratios[enebinfig]->fill(6, nparts[3] / nparts[1], 5);
          _h_ratios[enebinfig]->fill(7, nparts[5] / nparts[1], 5);
          _h_ratios[enebinfig]->fill(8, nparts[6] / nparts[1], 5);
          _h_ratios[enebinfig]->fill(9, nparts[9] / nparts[1], 5);
        }
      }

      if (nparts[0] > 1e-6) {
        _h_snn_Piratio->fill(energies[enebin], nparts[1] / nparts[0], 5);
        _h_snn_KaPiplus->fill(energies[enebin], nparts[2] / nparts[0], 5);
      }

      if (nparts[1] > 1e-6)
        _h_snn_KaPiminus->fill(energies[enebin],nparts[3] / nparts[1], 5);
      if (nparts[2] > 1e-6)
        _h_snn_Karatio->fill(energies[enebin],nparts[3] / nparts[2], 5);
      if (nparts[4] > 1e-6)
        _h_snn_Pratio->fill(energies[enebin],nparts[5] / nparts[4], 5);

      /// Sum the weight of the event
      if (c < 80) {
        _wght_Pi[enebin][cenbin]->fill();
        _wght_PiPlus[enebin][cenbin]->fill();
        _wght_Kaon[enebin][cenbin]->fill();
        _wght_KaonPlus[enebin][cenbin]->fill();
        _wght_Proton[enebin][cenbin]->fill();
        _wght_AntiProton[enebin][cenbin]->fill();
        _wght_npart_PiPlus[enebin]->fill();
        _wght_npart_PiMinus[enebin]->fill();
        _wght_npart_KaonPlus[enebin]->fill();
        _wght_npart_KaonMinus[enebin]->fill();
        _wght_npart_Proton[enebin]->fill();
        _wght_npart_AntiProton[enebin]->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      /// Normalisation
      for (size_t j = 0; j < 5; ++j) {
        for (size_t i = 0; i < 9; ++i) {
          if (_h_dpT_Pi[j][i]->integral() != 0 && _wght_Pi[j][i]->sumW() != 0)
            scale(_h_dpT_Pi[j][i], 1. / (TWOPI * 0.2 * _wght_Pi[j][i]->sumW()));
          if (_h_dpT_Piplus[j][i]->integral() != 0 && _wght_PiPlus[j][i]->sumW() != 0)
            scale(_h_dpT_Piplus[j][i], 1. / (TWOPI * 0.2 * _wght_PiPlus[j][i]->sumW()));
          if (_h_dpT_Kaon[j][i]->integral() != 0 && _wght_Kaon[j][i]->sumW() != 0)
            scale(_h_dpT_Kaon[j][i], 1. / (TWOPI * 0.2 * _wght_Kaon[j][i]->sumW()));
          if (_h_dpT_Kaonplus[j][i]->integral() != 0 && _wght_KaonPlus[j][i]->sumW() != 0)
            scale(_h_dpT_Kaonplus[j][i], 1. / (TWOPI * 0.2 * _wght_KaonPlus[j][i]->sumW()));
          if (_h_dpT_AntiProton[j][i]->integral() != 0 && _wght_Proton[j][i]->sumW() != 0)
            scale(_h_dpT_AntiProton[j][i], 1. / (TWOPI * 0.2 * _wght_Proton[j][i]->sumW()));
          if (_h_dpT_Proton[j][i]->integral() != 0 && _wght_AntiProton[j][i]->sumW() != 0)
            scale(_h_dpT_Proton[j][i], 1. / (TWOPI * 0.2 * _wght_AntiProton[j][i]->sumW()));
        }
      }

      /// Filling the bins with a value (here out of the defined
      /// screening range of the plot) when it has not been filled by
      /// anything, otherwise it won't want to plot
      for (size_t j = 0, N = energies.size(); j < N; ++j) {
        for (size_t i = 0, M = _h_npart_PiPlus[j]->numBins(); i < M; ++i)
          if (_h_npart_PiPlus[j]->bin(i).numEntries() == 0)
            _h_npart_PiPlus[j]->fillBin(i, -0.1);

        for (size_t i = 0, M = _h_npart_PiMinus[j]->numBins(); i < M; ++i)
          if (_h_npart_PiMinus[j]->bin(i).numEntries() == 0)
            _h_npart_PiMinus[j]->fillBin(i, -0.1);

        for (size_t i = 0, M = _h_npart_KaPlus[j]->numBins(); i < M; ++i)
          if (_h_npart_KaPlus[j]->bin(i).numEntries() == 0)
            _h_npart_KaPlus[j]->fillBin(i, -0.1);

        for (size_t i = 0, M = _h_npart_KaMinus[j]->numBins(); i < M; ++i)
          if (_h_npart_KaMinus[j]->bin(i).numEntries() == 0)
            _h_npart_KaMinus[j]->fillBin(i, -0.1);

        for (size_t i = 0, M = _h_npart_Proton[j]->numBins(); i < M; ++i)
          if (_h_npart_Proton[j]->bin(i).numEntries() == 0)
            _h_npart_Proton[j]->fillBin(i, -0.1);

        for (size_t i = 0, M = _h_npart_AntiProton[j]->numBins(); i < M; ++i)
          if (_h_npart_AntiProton[j]->bin(i).numEntries() == 0)
            _h_npart_AntiProton[j]->fillBin(i, -0.1);

        for (size_t i = 0, M = _h_npart_pT_PiPlus[j]->numBins(); i < M; ++i)
          if (_h_npart_pT_PiPlus[j]->bin(i).numEntries() == 0)
            _h_npart_pT_PiPlus[j]->fillBin(i, -0.1);

        for (size_t i = 0, M = _h_npart_pT_PiMinus[j]->numBins(); i < M; ++i)
          if (_h_npart_pT_PiMinus[j]->bin(i).numEntries() == 0)
            _h_npart_pT_PiMinus[j]->fillBin(i, -0.1);

        for (size_t i = 0, M = _h_npart_pT_KaPlus[j]->numBins(); i < M; ++i)
          if (_h_npart_pT_KaPlus[j]->bin(i).numEntries() == 0)
            _h_npart_pT_KaPlus[j]->fillBin(i, -0.1);

        for (size_t i = 0, M = _h_npart_pT_KaMinus[j]->numBins(); i < M; ++i)
          if (_h_npart_pT_KaMinus[j]->bin(i).numEntries() == 0)
            _h_npart_pT_KaMinus[j]->fillBin(i, -0.1);

        for (size_t i = 0, M = _h_npart_pT_Proton[j]->numBins(); i < M; ++i)
          if (_h_npart_pT_Proton[j]->bin(i).numEntries() == 0)
            _h_npart_pT_Proton[j]->fillBin(i, -0.1);

        for (size_t i = 0, M = _h_npart_pT_AntiProton[j]->numBins(); i < M; ++i)
          if (_h_npart_pT_AntiProton[j]->bin(i).numEntries() == 0)
            _h_npart_pT_AntiProton[j]->fillBin(i, -0.1);
      }

      for (size_t j = 0; j < 5; ++j)
        for (size_t i = 0; i < 9; ++i) {
          if (_h_npart_Piratio[j]->bin(i).numEntries() == 0)
            _h_npart_Piratio[j]->fillBin(i, -0.1);

          if (_h_npart_Karatio[j]->bin(i).numEntries() == 0)
            _h_npart_Karatio[j]->fillBin(i, -0.1);

          if (_h_npart_Pratio[j]->bin(i).numEntries() == 0)
            _h_npart_Pratio[j]->fillBin(i, -0.1);

          if (_h_npart_KaPi[j]->bin(i).numEntries() == 0)
            _h_npart_KaPi[j]->fillBin(i, -0.1);

          if (_h_npart_AntiPPi[j]->bin(i).numEntries() == 0)
            _h_npart_AntiPPi[j]->fillBin(i, -0.1);

          if (_h_npart_KaPiplus[j]->bin(i).numEntries() == 0)
            _h_npart_KaPiplus[j]->fillBin(i, -0.1);

          if (_h_npart_PPiplus[j]->bin(i).numEntries() == 0)
            _h_npart_PPiplus[j]->fillBin(i, -0.1);
        }


      for (size_t j = 0; j < 2; ++j) {
        for (size_t i = 0, N = _h_ratios[j]->numBins(); i < N; ++i)
          if (_h_ratios[j]->bin(i).numEntries() == 0)
            _h_ratios[j]->fillBin(i, -0.1);
        for (size_t i = 0, N = _h_yields[j]->numBins(); i < N; ++i)
          if (_h_yields[j]->bin(i).numEntries() == 0)
            _h_yields[j]->fillBin(i, -0.1);
      }

      for (size_t i = 0, N = energies.size(); i < N; ++i) {
        if (_h_snn_npart_PiPlus->bin(i).numEntries() == 0)
          _h_snn_npart_PiPlus->fillBin(i, -0.1);

        if (_h_snn_npart_PiMinus->bin(i).numEntries() == 0)
          _h_snn_npart_PiMinus->fillBin(i, -0.1);

        if (_h_snn_npart_KaPlus->bin(i).numEntries() == 0)
          _h_snn_npart_KaPlus->fillBin(i, -0.1);

        if (_h_snn_npart_KaMinus->bin(i).numEntries() == 0)
          _h_snn_npart_KaMinus->fillBin(i, -0.1);

        if (_h_snn_npart_Proton->bin(i).numEntries() == 0)
          _h_snn_npart_Proton->fillBin(i, -0.1);

        if (_h_snn_npart_AntiProton->bin(i).numEntries() == 0)
          _h_snn_npart_AntiProton->fillBin(i, -0.1);

        if (_h_snn_mt_PiPlus->bin(i).numEntries() == 0)
          _h_snn_mt_PiPlus->fillBin(i, -0.1);

        if (_h_snn_mt_PiMinus->bin(i).numEntries() == 0)
          _h_snn_mt_PiMinus->fillBin(i, -0.1);

        if (_h_snn_mt_KaPlus->bin(i).numEntries() == 0)
          _h_snn_mt_KaPlus->fillBin(i, -0.1);

        if (_h_snn_mt_KaMinus->bin(i).numEntries() == 0)
          _h_snn_mt_KaMinus->fillBin(i, -0.1);

        if (_h_snn_mt_Proton->bin(i).numEntries() == 0)
          _h_snn_mt_Proton->fillBin(i, -0.1);

        if (_h_snn_mt_AntiProton->bin(i).numEntries() == 0)
          _h_snn_mt_AntiProton->fillBin(i, -0.1);

        if (_h_snn_KaPiplus->bin(i).numEntries() == 0)
          _h_snn_KaPiplus->fillBin(i, -0.1);

        if (_h_snn_KaPiminus->bin(i).numEntries() == 0)
          _h_snn_KaPiminus->fillBin(i, -0.1);

        if (_h_snn_Piratio->bin(i).numEntries() == 0)
          _h_snn_Piratio->fillBin(i, -0.1);

        if (_h_snn_Karatio->bin(i).numEntries() == 0)
          _h_snn_Karatio->fillBin(i, -0.1);

        if (_h_snn_Pratio->bin(i).numEntries() == 0)
          _h_snn_Pratio->fillBin(i, -0.1);
      }
    }


  private:

    /// @name Histograms
    /// @{
    vector<vector<Histo1DPtr>> _h_dpT_Pi;
    vector<vector<Histo1DPtr>> _h_dpT_Piplus;
    vector<vector<Histo1DPtr>> _h_dpT_Kaon;
    vector<vector<Histo1DPtr>> _h_dpT_Kaonplus;
    vector<vector<Histo1DPtr>> _h_dpT_AntiProton;
    vector<vector<Histo1DPtr>> _h_dpT_Proton;

    vector<vector<CounterPtr>> _wght_Pi;
    vector<vector<CounterPtr>> _wght_PiPlus;
    vector<vector<CounterPtr>> _wght_Kaon;
    vector<vector<CounterPtr>> _wght_KaonPlus;
    vector<vector<CounterPtr>> _wght_Proton;
    vector<vector<CounterPtr>> _wght_AntiProton;

    vector<Profile1DPtr> _h_npart_Piratio;
    vector<Profile1DPtr> _h_npart_Karatio;
    vector<Profile1DPtr> _h_npart_Pratio;
    vector<Profile1DPtr> _h_npart_KaPi;
    vector<Profile1DPtr> _h_npart_AntiPPi;
    vector<Profile1DPtr> _h_npart_KaPiplus;
    vector<Profile1DPtr> _h_npart_PPiplus;
    vector<Profile1DPtr> _h_yields;
    vector<Profile1DPtr> _h_ratios;

    vector<Histo1DPtr> _h_npart_PiPlus;
    vector<Histo1DPtr> _h_npart_PiMinus;
    vector<Histo1DPtr> _h_npart_KaPlus;
    vector<Histo1DPtr> _h_npart_KaMinus;
    vector<Histo1DPtr> _h_npart_Proton;
    vector<Histo1DPtr> _h_npart_AntiProton;

    vector<CounterPtr> _wght_npart_PiPlus;
    vector<CounterPtr> _wght_npart_PiMinus;
    vector<CounterPtr> _wght_npart_KaonPlus;
    vector<CounterPtr> _wght_npart_KaonMinus;
    vector<CounterPtr> _wght_npart_Proton;
    vector<CounterPtr> _wght_npart_AntiProton;

    vector<Profile1DPtr> _h_npart_pT_PiPlus;
    vector<Profile1DPtr> _h_npart_pT_PiMinus;
    vector<Profile1DPtr> _h_npart_pT_KaPlus;
    vector<Profile1DPtr> _h_npart_pT_KaMinus;
    vector<Profile1DPtr> _h_npart_pT_Proton;
    vector<Profile1DPtr> _h_npart_pT_AntiProton;

    Histo1DPtr _h_snn_npart_PiPlus;
    Histo1DPtr _h_snn_npart_PiMinus;
    Histo1DPtr _h_snn_npart_KaPlus;
    Histo1DPtr _h_snn_npart_KaMinus;
    Histo1DPtr _h_snn_npart_Proton;
    Histo1DPtr _h_snn_npart_AntiProton;
    Profile1DPtr _h_snn_mt_PiPlus;
    Profile1DPtr _h_snn_mt_PiMinus;
    Profile1DPtr _h_snn_mt_KaPlus;
    Profile1DPtr _h_snn_mt_KaMinus;
    Profile1DPtr _h_snn_mt_Proton;
    Profile1DPtr _h_snn_mt_AntiProton;
    Profile1DPtr _h_snn_KaPiplus;
    Profile1DPtr _h_snn_KaPiminus;
    Profile1DPtr _h_snn_Piratio;
    Profile1DPtr _h_snn_Karatio;
    Profile1DPtr _h_snn_Pratio;
    /// @}

    /// @name Variables
    /// @{
    vector<double> energies;
    vector<double> centralities;
    int cenbin, enebin = 0, enebinfig = 0;
    double nprtcl, nPi[5], nPiPlus[5], nKaon[5], nKaonPlus[5], nProton[5],
      nAntiProton[5];
    // The following vector contains the counters for all particles used in
    // Fig. 25. In the right order : pi+, pi-, K+, K-, p, Antip, Lambda,
    // AntiLambda, Xi, AntiXi
    double nparts[10];
    /// @}

  };


  DECLARE_RIVET_PLUGIN(STAR_2017_I1510593);

}
