// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// @brief Forward Z production cross-section in pp collisions at 13 TeV
  class LHCB_2021_I1990313 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2021_I1990313);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      ZFinder zmumufinder(FinalState(), Cuts::absetaIn(2.0, 4.5) && Cuts::pT > 20 * GeV, PID::MUON, 60*GeV, 120*GeV);
      declare(zmumufinder, "ZmumuFinder");

      // Book histograms
      // specify custom binning
      book(_h_sigma_vs_y,    18, 1, 1);
      book(_h_sigma_vs_pt,   19, 1, 1);
      book(_h_sigma_vs_phi,  20, 1, 1);
      double ylow, yhigh;
      for (int i = 1; i < 6; ++i) {
      	ylow = 1.5 + 0.5*(double)i;
      	yhigh = ylow + 0.5;
      	{Histo1DPtr tmp; _h_sigma_vs_ypt.add (ylow, yhigh, book(tmp, 21, 1, i) );};
      	{Histo1DPtr tmp; _h_sigma_vs_yphi.add(ylow, yhigh, book(tmp, 22, 1, i) );};
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Retrieve dressed leptons, sorted by pT
      const ZFinder& zmumufinder = apply<ZFinder>(event, "ZmumuFinder");
      if(zmumufinder.empty()) vetoEvent;
      if(zmumufinder.bosons().size() > 1)
        MSG_WARNING("Found multiple (" << zmumufinder.bosons().size() << ") Z -> mu+ mu- decays!");

      // Z momenta
      FourMomentum zmumu = zmumufinder.bosons()[0].momentum();
      if(zmumufinder.constituentLeptons().size() < 2) vetoEvent;

      const Particle& muon_p = zmumufinder.constituentLeptons()[0];
      const Particle& muon_m = zmumufinder.constituentLeptons()[1];

      const double diffphi = deltaPhi(muon_p, muon_m);
      const double diffpsd = deltaEta(muon_p, muon_m);
      const double accphi = M_PI - diffphi;
      const double angular = tan(accphi/2) / cosh(diffpsd/2);

      _h_sigma_vs_y->fill(zmumu.rapidity());
      _h_sigma_vs_pt->fill(zmumu.pT()/GeV);
      _h_sigma_vs_phi->fill(angular);
      _h_sigma_vs_ypt.fill(zmumu.rapidity(), zmumu.pT()/GeV);
      _h_sigma_vs_yphi.fill(zmumu.rapidity(), angular);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double xs = crossSection()/picobarn;
      double scale_f = xs/sumOfWeights()/2.;
      _h_sigma_vs_y->scaleW(scale_f);
      _h_sigma_vs_pt->scaleW(scale_f);
      _h_sigma_vs_phi->scaleW(scale_f);
      for (Histo1DPtr h : _h_sigma_vs_ypt.histos())  h->scaleW(scale_f*2.);
      for (Histo1DPtr h : _h_sigma_vs_yphi.histos()) h->scaleW(scale_f*2.);
    }

    ///@}

    Histo1DPtr _h_sigma_vs_y, _h_sigma_vs_pt, _h_sigma_vs_phi;
    BinnedHistogram _h_sigma_vs_ypt, _h_sigma_vs_yphi;


  };


  RIVET_DECLARE_PLUGIN(LHCB_2021_I1990313);

}
