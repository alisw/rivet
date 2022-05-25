// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Particle.hh"

namespace Rivet {


  class CMS_2010_PAS_QCD_10_024 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CMS_2010_PAS_QCD_10_024() : Analysis("CMS_2010_PAS_QCD_10_024"),
		       _weight_pt05_eta08(0.), _weight_pt10_eta08(0.),
		       _weight_pt05_eta24(0.), _weight_pt10_eta24(0.) {  }


    void init() {
      declare(ChargedFinalState((Cuts::etaIn(-0.8, 0.8) && Cuts::pT >=  0.5*GeV)), "CFS_08_05");
      declare(ChargedFinalState((Cuts::etaIn(-0.8, 0.8) && Cuts::pT >=  1.0*GeV)), "CFS_08_10");
      declare(ChargedFinalState((Cuts::etaIn(-2.4, 2.4) && Cuts::pT >=  0.5*GeV)), "CFS_24_05");
      declare(ChargedFinalState((Cuts::etaIn(-2.4, 2.4) && Cuts::pT >=  1.0*GeV)), "CFS_24_10");

      size_t offset = 0;
      if (isCompatibleWithSqrtS(7000)) offset = 0;
      if (isCompatibleWithSqrtS( 900)) offset = 4;
      book(_hist_dNch_deta_pt05_eta08 ,1+offset, 1, 1);
      book(_hist_dNch_deta_pt10_eta08 ,2+offset, 1, 1);
      book(_hist_dNch_deta_pt05_eta24 ,3+offset, 1, 1);
      book(_hist_dNch_deta_pt10_eta24 ,4+offset, 1, 1);
    }


    void analyze(const Event& event) {
      const double weight = 1.0;
      const ChargedFinalState& cfs_08_05 = apply<ChargedFinalState>(event, "CFS_08_05");
      const ChargedFinalState& cfs_08_10 = apply<ChargedFinalState>(event, "CFS_08_10");
      const ChargedFinalState& cfs_24_05 = apply<ChargedFinalState>(event, "CFS_24_05");
      const ChargedFinalState& cfs_24_10 = apply<ChargedFinalState>(event, "CFS_24_10");

      // Plot distributions
      if(!cfs_08_05.particles().empty()) _weight_pt05_eta08 += weight;
      if(!cfs_24_05.particles().empty()) _weight_pt05_eta24 += weight;
      for (const Particle& p : cfs_24_05.particles()) {
        _hist_dNch_deta_pt05_eta24->fill(p.eta(), weight);
        if(!cfs_08_05.particles().empty())
	  _hist_dNch_deta_pt05_eta08->fill(p.eta(), weight);
      }
      if(!cfs_08_10.particles().empty()) _weight_pt10_eta08 += weight;
      if(!cfs_24_10.particles().empty()) _weight_pt10_eta24 += weight;
      for (const Particle& p : cfs_24_10.particles()) {
        _hist_dNch_deta_pt10_eta24->fill(p.eta(), weight);
	if(!cfs_08_10.particles().empty())
	  _hist_dNch_deta_pt10_eta08->fill(p.eta(), weight);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_hist_dNch_deta_pt05_eta08,1./_weight_pt05_eta08);
      scale(_hist_dNch_deta_pt10_eta08,1./_weight_pt10_eta08);
      scale(_hist_dNch_deta_pt05_eta24,1./_weight_pt05_eta24);
      scale(_hist_dNch_deta_pt10_eta24,1./_weight_pt10_eta24);
    }


  private:

    Histo1DPtr _hist_dNch_deta_pt05_eta08;
    Histo1DPtr _hist_dNch_deta_pt10_eta08;
    Histo1DPtr _hist_dNch_deta_pt05_eta24;
    Histo1DPtr _hist_dNch_deta_pt10_eta24;
    double _weight_pt05_eta08,_weight_pt10_eta08,_weight_pt05_eta24,_weight_pt10_eta24;
  };


  // Hook for the plugin system
  RIVET_DECLARE_PLUGIN(CMS_2010_PAS_QCD_10_024);

}
