// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// @brief MC validation analysis for higgs [-> tau tau] events
  class MC_HINC : public Analysis {
  public:

    /// Default constructor
    MC_HINC()
      : Analysis("MC_HINC")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      Cut cut = Cuts::abseta < 3.5 && Cuts::pT > 25*GeV;
      /// @todo Urk, abuse! Need explicit HiggsFinder and TauFinder?
      ZFinder hfinder(FinalState(), cut, PID::TAU, 115*GeV, 135*GeV, 0.0, ZFinder::ClusterPhotons::NONE, ZFinder::AddPhotons::NO, 125*GeV);
      declare(hfinder, "Hfinder");
      book(_h_H_mass ,"H_mass", 50, 119.7, 120.3);
      book(_h_H_pT ,"H_pT", logspace(100, 1.0, 0.5*(sqrtS()>0.?sqrtS():14000.)/GeV));
      book(_h_H_pT_peak ,"H_pT_peak", 25, 0.0, 25.0);
      book(_h_H_y ,"H_y", 40, -4, 4);
      book(_h_H_phi ,"H_phi", 25, 0.0, TWOPI);
      book(_h_lepton_pT ,"lepton_pT", logspace(100, 10.0, 0.25*(sqrtS()>0.?sqrtS():14000.)/GeV));
      book(_h_lepton_eta ,"lepton_eta", 40, -4, 4);
    }


    /// Do the analysis
    void analyze(const Event & e) {
      const ZFinder& hfinder = apply<ZFinder>(e, "Hfinder");
      if (hfinder.bosons().size() != 1) vetoEvent;
      const double weight = 1.0;

      FourMomentum hmom(hfinder.bosons()[0].momentum());
      _h_H_mass->fill(hmom.mass()/GeV, weight);
      _h_H_pT->fill(hmom.pT()/GeV, weight);
      _h_H_pT_peak->fill(hmom.pT()/GeV, weight);
      _h_H_y->fill(hmom.rapidity(), weight);
      _h_H_phi->fill(hmom.phi(), weight);
      for (const Particle& l : hfinder.constituents()) {
        _h_lepton_pT->fill(l.pT()/GeV, weight);
        _h_lepton_eta->fill(l.eta(), weight);
      }
    }


    /// Finalize
    void finalize() {
      const double xsec = crossSection()/picobarn;
      normalize(_h_H_mass, xsec);
      normalize(_h_H_pT, xsec);
      normalize(_h_H_pT_peak, xsec);
      normalize(_h_H_y, xsec);
      normalize(_h_H_phi, xsec);
      normalize(_h_lepton_pT, xsec);
      normalize(_h_lepton_eta, xsec);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_H_mass;
    Histo1DPtr _h_H_pT;
    Histo1DPtr _h_H_pT_peak;
    Histo1DPtr _h_H_y;
    Histo1DPtr _h_H_phi;
    Histo1DPtr _h_lepton_pT;
    Histo1DPtr _h_lepton_eta;
    //@}

  };



  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MC_HINC);

}
