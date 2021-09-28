// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"

namespace Rivet {


  class ATLAS_2015_I1351916 : public Analysis {
  public:

    /// Constructor
    ATLAS_2015_I1351916(const string name="ATLAS_2015_I1351916", size_t mode=0,
                        const string ref_data="ATLAS_2015_I1351916") : Analysis(name) {
      _mode = mode; // pick electron channel by default
      setRefDataName(ref_data);
    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs;

      IdentifiedFinalState bareleptons(fs);
      bareleptons.acceptIdPair(_mode? PID::MUON : PID::ELECTRON);

      const Cut cuts = (_mode == 0) ? (Cuts::pT > 25*GeV && Cuts::abseta < 4.9) : (Cuts::pT > 20*GeV && Cuts::abseta < 2.47);
      DressedLeptons leptons(fs, bareleptons, 0.1, cuts, true);
      declare(leptons, "leptons");


      // Book dummy histograms for heterogeneous merging
      const Scatter2D& ref = refData(_mode? 4 : 2, 1, 2);
      book(_h["NCC_pos"], "_ncc_pos", ref);
      book(_h["NCC_neg"], "_ncc_neg", ref);
      book(_s["CC"], _mode ? 4 : 2, 1, 2, true);

      if (_mode == 0) { // electron-channel only
        const Scatter2D& ref_cf = refData(3, 1, 2);
        book(_h["NCF_pos"], "_ncf_pos", ref_cf);
        book(_h["NCF_neg"], "_ncf_neg", ref_cf);
        book(_s["CF"], 3, 1, 2, true);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& e) {

      // Get and cut on dressed leptons
      const vector<DressedLepton>& leptons = apply<DressedLeptons>(e, "leptons").dressedLeptons();
      if (leptons.size() != 2) vetoEvent; // require exactly two leptons
      if (leptons[0].charge3() * leptons[1].charge3() > 0) vetoEvent; // require opposite charge

      // Identify lepton vs antilepton
      const Particle& lpos = leptons[(leptons[0].charge3() > 0) ? 0 : 1];
      const Particle& lneg = leptons[(leptons[0].charge3() < 0) ? 0 : 1];

      string label = "N";
      if (_mode == 1) {// electron channel
        label += "CC"; // only central-central for muons
      } else { // electron channel
        const double eta1 = lpos.abseta();
        const double eta2 = lneg.abseta();
        if ( (eta1 < 2.47 && inRange(eta2, 2.5, 4.9)) || (eta2 < 2.47 && inRange(eta1, 2.5, 4.9)) )
          label += "CF"; // central-forward
        else if (eta1 < 2.47 && eta2 < 2.47)
          label += "CC"; // central-central
        else vetoEvent; // ain't no forward-forward
      }

      const double cosThetaStar = cosCollinsSoper(lneg, lpos);
      const double mll = (lpos.mom() + lneg.mom()).mass();
      label += cosThetaStar < 0.0?  "_neg" : "_pos";
      _h[label]->fill(mll/GeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSectionPerEvent() / picobarn;
      for (const auto& key_hist : _h) scale(key_hist.second, sf);
      divide(*_h["NCC_pos"] - *_h["NCC_neg"], *_h["NCC_pos"] + *_h["NCC_neg"], _s["CC"]);
      if (!_mode)  divide(*_h["NCF_pos"] - *_h["NCF_neg"], *_h["NCF_pos"] + *_h["NCF_neg"], _s["CF"]);
    }


    // Cosine of the decay angle in the Collins-Soper frame
    double cosCollinsSoper(const FourMomentum& l1, const FourMomentum& l2) {
      const FourMomentum ll = l1 + l2;
      const double nom  = (l1.E() + l1.pz()) * (l2.E() - l2.pz()) - (l1.E() - l1.pz()) * (l2.E() + l2.pz());
      const double denom = ll.mass() * sqrt( sqr(ll.mass()) + sqr(ll.pt()) );
      return sign(ll.pz()) * safediv(nom, denom); // protect against division by zero, you never know...
    }

    //@}


  protected:

    /// Electron or muon mode = 0 or 1, for use by derived _EL, _MU analysis classes
    size_t _mode;


  private:

    /// Histograms
    map<string, Histo1DPtr> _h;
    /// Asymmetries
    map<string, Scatter2DPtr> _s;

  };



  class ATLAS_2015_I1351916_EL : public ATLAS_2015_I1351916 {
  public:
    ATLAS_2015_I1351916_EL() : ATLAS_2015_I1351916("ATLAS_2015_I1351916_EL", 0) { }
  };


  class ATLAS_2015_I1351916_MU : public ATLAS_2015_I1351916 {
  public:
    ATLAS_2015_I1351916_MU() : ATLAS_2015_I1351916("ATLAS_2015_I1351916_MU", 1) { }
  };


  DECLARE_RIVET_PLUGIN(ATLAS_2015_I1351916);
  DECLARE_RIVET_PLUGIN(ATLAS_2015_I1351916_EL);
  DECLARE_RIVET_PLUGIN(ATLAS_2015_I1351916_MU);

}
