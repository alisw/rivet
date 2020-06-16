// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// Rivet analysis class for ATLAS 2010 minimum bias analysis
  class ATLAS_2010_S8918562 : public Analysis {
  public:

    /// Helper for collectively filling Nch, pT, eta, and pT vs. Nch histograms
    void fillPtEtaNch(const ChargedFinalState& cfs, const int nchcut, const string& label) {
      // Get number of particles and skip if event fails cut
      const int nch = cfs.size();
      if (nch < nchcut) return;

      // Fill nch
      _h[label + "_nch"]->fill(nch);
      // Loop over particles, fill pT, eta and ptnch
      for (const Particle& p : cfs.particles()) {
        const double pt = p.pT();
        _h[label + "_pt"]->fill(pt/GeV, 1.0/pt);
        _h[label + "_eta"]->fill(p.eta());
        if (_p[label + "_ptnch"])  _p[label + "_ptnch"]->fill(nch, pt/GeV);
      }
    }


    /// Default constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2010_S8918562);


    /// Initialization, called once before running
    void init() {
      // Projections
      const ChargedFinalState cfs100(Cuts::abseta < 2.5 && Cuts::pT > 100*MeV);
      declare(cfs100, "CFS100");
      const ChargedFinalState cfs500(Cuts::abseta < 2.5 && Cuts::pT > 500*MeV);
      declare(cfs500, "CFS500");
      const ChargedFinalState cfs2500(Cuts::abseta < 2.5 && Cuts::pT > 2500*MeV);
      declare(cfs2500, "CFS2500");

      // Book histograms
      if (fuzzyEquals(sqrtS()/GeV, 900)) {
        book(_h["pt100_nch2_nch"],   18, 1, 1);
        book(_h["pt100_nch2_pt"],    11, 1, 1);
        book(_h["pt100_nch2_eta"],    4, 1, 1);
        book(_p["pt100_nch2_ptnch"], 24, 1, 1);

        book(_h["pt100_nch20_nch"], 34, 1, 1);
        book(_h["pt100_nch20_pt"],  30, 1, 1);
        book(_h["pt100_nch20_eta"], 26, 1, 1);

        book(_h["pt500_nch1_nch"],   15, 1, 1);
        book(_h["pt500_nch1_pt"],     8, 1, 1);
        book(_h["pt500_nch1_eta"],    1, 1, 1);
        book(_p["pt500_nch1_ptnch"], 22, 1, 1);

        book(_h["pt500_nch6_nch"], 20, 1, 1);
        book(_h["pt500_nch6_pt"],  13, 1, 1);
        book(_h["pt500_nch6_eta"],  6, 1, 1);

        book(_h["pt2500_nch1_nch"],   36, 1, 1);
        book(_h["pt2500_nch1_pt"],    32, 1, 1);
        book(_h["pt2500_nch1_eta"],   28, 1, 1);
        book(_p["pt2500_nch1_ptnch"], 38, 1, 1);

      } else if (fuzzyEquals(sqrtS()/GeV, 2360)) {

        book(_h["pt500_nch1_nch"], 16, 1, 1);
        book(_h["pt500_nch1_pt"],   9, 1, 1);
        book(_h["pt500_nch1_eta"],  2, 1, 1);
        _p["pt500_nch1_ptnch"] = nullptr;

      } else if (fuzzyEquals(sqrtS()/GeV, 7000)) {

        book(_h["pt100_nch2_nch"],   19, 1, 1);
        book(_h["pt100_nch2_pt"],    12, 1, 1);
        book(_h["pt100_nch2_eta"],    5, 1, 1);
        book(_p["pt100_nch2_ptnch"], 25, 1, 1);

        book(_h["pt100_nch20_nch"], 35, 1, 1);
        book(_h["pt100_nch20_pt"],  31, 1, 1);
        book(_h["pt100_nch20_eta"], 27, 1, 1);

        book(_h["pt500_nch1_nch"],   17, 1, 1);
        book(_h["pt500_nch1_pt"],    10, 1, 1);
        book(_h["pt500_nch1_eta"],    3, 1, 1);
        book(_p["pt500_nch1_ptnch"], 23, 1, 1);

        book(_h["pt500_nch6_nch"], 21, 1, 1);
        book(_h["pt500_nch6_pt"],  14, 1, 1);
        book(_h["pt500_nch6_eta"],  7, 1, 1);

        book(_h["pt2500_nch1_nch"],   37, 1, 1);
        book(_h["pt2500_nch1_pt"],    33, 1, 1);
        book(_h["pt2500_nch1_eta"],   29, 1, 1);
        book(_p["pt2500_nch1_ptnch"], 39, 1, 1);

      } else {
        throw LogicError("The ATLAS_2010_S8918562 analysis is only valid for sqrt(s) = 900, 2360 and 7000 GeV!");
      }

    }


    void analyze(const Event& event) {
      // 100 GeV final states
      if (!fuzzyEquals(sqrtS()/GeV, 2360)) {
        const ChargedFinalState& cfs100 = apply<ChargedFinalState>(event, "CFS100");
        // nch>=2
        fillPtEtaNch(cfs100, 2, "pt100_nch2");
        // nch>=20
        fillPtEtaNch(cfs100, 20, "pt100_nch20");
      }

      // 500 GeV final states
      const ChargedFinalState& cfs500 = apply<ChargedFinalState>(event, "CFS500");
      // nch>=1
      fillPtEtaNch(cfs500, 1, "pt500_nch1");
      // nch>=6
      if (!fuzzyEquals(sqrtS()/GeV, 2360)) {
        fillPtEtaNch(cfs500, 6, "pt500_nch6");
      }

      // 2500 GeV final states
      if (!fuzzyEquals(sqrtS()/GeV, 2360)) {
        const ChargedFinalState& cfs2500 = apply<ChargedFinalState>(event, "CFS2500");
        // nch>=1
        fillPtEtaNch(cfs2500, 1, "pt2500_nch1");
      }

    }


    void finalize() {

      double sf = safediv(1.0, _h["pt500_nch1_nch"]->integral(true), 1.0);
      scale(_h["pt500_nch1_nch"], sf);
      scale(_h["pt500_nch1_pt"],  sf/TWOPI/5);
      scale(_h["pt500_nch1_eta"], sf);

      if (!fuzzyEquals(sqrtS()/GeV, 2360)) {
        sf = safediv(1.0, _h["pt100_nch2_nch"]->integral(true), 1.0);
        scale(_h["pt100_nch2_nch"], sf);
        scale(_h["pt100_nch2_pt"],  sf/TWOPI/5);
        scale(_h["pt100_nch2_eta"], sf);

        sf = safediv(1.0, _h["pt100_nch20_nch"]->integral(true), 1.0);
        scale(_h["pt100_nch20_nch"], sf);
        scale(_h["pt100_nch20_pt"],  sf/TWOPI/5);
        scale(_h["pt100_nch20_eta"], sf);

        sf = safediv(1.0, _h["pt500_nch6_nch"]->integral(true), 1.0);
        scale(_h["pt500_nch6_nch"], sf);
        scale(_h["pt500_nch6_pt"],  sf/TWOPI/5);
        scale(_h["pt500_nch6_eta"], sf);

        sf = safediv(1.0, _h["pt2500_nch1_nch"]->integral(true), 1.0);
        scale(_h["pt2500_nch1_nch"], sf);
        scale(_h["pt2500_nch1_pt"],  sf/TWOPI/5);
        scale(_h["pt2500_nch1_eta"], sf);
      }
    }


  private:

    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2010_S8918562);

}
