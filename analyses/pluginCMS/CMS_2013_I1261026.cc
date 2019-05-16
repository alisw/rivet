// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// Jet and underlying event properties as a function of particle multiplicity
  class CMS_2013_I1261026 : public Analysis {
  public:

    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2013_I1261026);

    void init() {
      FastJets jetpro(ChargedFinalState(-2.4, 2.4, 0.25*GeV), FastJets::ANTIKT, 0.5);
      declare(jetpro, "Jets");

      const ChargedFinalState cfs(-2.4, 2.4, 0.25*GeV);
      declare(cfs, "CFS250");

      // For min bias trigger
      const ChargedFinalState cfsBSCplus(3.23, 4.65, 500*MeV);
      declare(cfsBSCplus, "cfsBSCplus");

      const ChargedFinalState cfsBSCminus(-4.65, -3.23, 500*MeV);
      declare(cfsBSCminus, "cfsBSCminus");

      // Histograms:
      _h_AllTrkMeanPt            = bookProfile1D(1, 1, 1);
      _h_SoftTrkMeanPt           = bookProfile1D(2, 1, 1);
      _h_IntrajetTrkMeanPt       = bookProfile1D(3, 1, 1);
      _h_IntrajetLeaderTrkMeanPt = bookProfile1D(4, 1, 1);
      _h_MeanJetPt               = bookProfile1D(5, 1, 1);
      _h_JetRate5GeV             = bookProfile1D(6, 1, 1);
      _h_JetRate30GeV            = bookProfile1D(7, 1, 1);

      for (int ihist = 0; ihist < 5; ++ihist) {
        _h_JetSpectrum[ihist] = bookHisto1D(ihist+8, 1, 1);
        _h_JetStruct[ihist]   = bookHisto1D(ihist+13, 1, 1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // MinBias trigger
      const ChargedFinalState& cfsBSCplus = apply<ChargedFinalState>(event, "cfsBSCplus");
      if (cfsBSCplus.empty()) vetoEvent;
      const ChargedFinalState& cfsBSCminus = apply<ChargedFinalState>(event, "cfsBSCminus");
      if (cfsBSCminus.empty()) vetoEvent;

      const ChargedFinalState& cfsp = apply<ChargedFinalState>(event, "CFS250");
      if (cfsp.empty()) vetoEvent;

      const FastJets& jetpro = apply<FastJets>(event, "Jets");
      const Jets& jets = jetpro.jetsByPt(5.0*GeV);

      const int mult = cfsp.size();

      int multbin[6] = { 10, 30, 50, 80, 110, 140 };
      for (int ibin = 0; ibin < 5; ++ibin) {
        if (mult > multbin[ibin] && mult <= multbin[ibin + 1]) {
          eventDecomp(event, mult, ibin, weight);
          unsigned int jetCounter5GeV(0), jetCounter30GeV(0);
          for (size_t ijets = 0; ijets < jets.size(); ++ijets) {
            if (jets[ijets].abseta() < 1.9) {
              _h_JetSpectrum[ibin]->fill(jets[ijets].pT()/GeV, weight);
              _h_MeanJetPt->fill(mult, jets[ijets].pT()/GeV, weight);
              if (jets[ijets].pT() > 5*GeV)  ++jetCounter5GeV;
              if (jets[ijets].pT() > 30*GeV)  ++jetCounter30GeV;
            }
          }
          _h_JetRate5GeV->fill( mult,  jetCounter5GeV, weight);
          _h_JetRate30GeV->fill(mult, jetCounter30GeV, weight);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (size_t i = 0; i < 5; ++i) {
        normalize(_h_JetSpectrum[i],  4.0);
        normalize(_h_JetStruct[i]  , 0.08);
      }

    }

    void eventDecomp(const Event& event, int mult, size_t ibin, double weight) {

      struct TrkInJet { double pt; double eta; double phi; double R; };
      TrkInJet jetConstituents[100][100]; //1-st index - the number of the jet, 2-nd index - track in the jet
      TrkInJet jetsEv[100];
      size_t j[100];
      size_t jCount = 0;

      for (size_t i = 0; i < 100; ++i) {
        j[i] = 0;
        jetsEv[i].pt = 0;
        jetsEv[i].eta = 0;
        jetsEv[i].phi = 0;
        for (size_t k = 0; k < 100; ++k) {
          jetConstituents[i][k].pt = 0;
          jetConstituents[i][k].phi = 0;
          jetConstituents[i][k].eta = 0;
          jetConstituents[i][k].R = 0;
        }
      }

      const FastJets& jetpro = apply<FastJets>(event, "Jets");
      const Jets& jets = jetpro.jetsByPt(5.0*GeV);

      // Start event decomp

      for (size_t ijets = 0; ijets < jets.size(); ++ijets) {
        jetsEv[ijets].pt = jets[ijets].pT();
        jetsEv[ijets].eta = jets[ijets].eta();
        jetsEv[ijets].phi = jets[ijets].phi();
        ++jCount;
      }

      const ChargedFinalState& cfsp = apply<ChargedFinalState>(event, "CFS250");
      foreach (const Particle& p, cfsp.particles()) {
        _h_AllTrkMeanPt->fill(mult, p.pT()/GeV, weight);
        int flag = 0;
        for (size_t i = 0; i < jCount; ++i) {
          const double delta_phi = deltaPhi(jetsEv[i].phi, p.phi());
          const double delta_eta = jetsEv[i].eta - p.eta();
          const double R = sqrt(delta_phi * delta_phi + delta_eta * delta_eta);
          if (R <= 0.5) {
            ++flag;
            jetConstituents[i][j[i]].pt = p.pT();
            jetConstituents[i][j[i]].R = R;
            ++j[i];
          }
        }
        if (flag == 0) _h_SoftTrkMeanPt->fill(mult, p.pT()/GeV, weight);
      }

      for (size_t i = 0; i < jCount; ++i) {
        double ptInjetLeader = 0;
        if (!inRange(jetsEv[i].eta, -1.9, 1.9)) continue; // only fully reconstructed jets for internal jet studies
        for (size_t k = 0; k < j[i]; ++k) {
          _h_IntrajetTrkMeanPt->fill(mult, jetConstituents[i][k].pt , weight);
          _h_JetStruct[ibin]->fill(jetConstituents[i][k].R, jetConstituents[i][k].pt/jetsEv[i].pt);
          if (ptInjetLeader < jetConstituents[i][k].pt) ptInjetLeader = jetConstituents[i][k].pt;
        }
        if (ptInjetLeader != 0) _h_IntrajetLeaderTrkMeanPt->fill(mult, ptInjetLeader, weight);
      }

    }


  private:

    Profile1DPtr _h_AllTrkMeanPt, _h_SoftTrkMeanPt;
    Profile1DPtr _h_IntrajetTrkMeanPt, _h_IntrajetLeaderTrkMeanPt;
    Profile1DPtr _h_MeanJetPt;
    Profile1DPtr _h_JetRate5GeV, _h_JetRate30GeV;

    Histo1DPtr _h_JetSpectrum[5];
    Histo1DPtr _h_JetStruct[5];

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2013_I1261026);

}
