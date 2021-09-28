// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {



  /// @brief MC validation analysis for truth-MET measurement
  /// @todo Add plots for MET based on prompt invisibles
  class MC_MET : public Analysis {
  public:

    MC_MET()
      : Analysis("MC_MET")
    {    }


    void init() {
      FinalState inclfs;
      FinalState calofs(Cuts::abseta < 5);
      declare(MissingMomentum(inclfs), "InclMET");
      declare(MissingMomentum(calofs), "CaloMET");

      book(_h_met_incl ,"met_incl", logspace(30, 1, 150));
      book(_h_met_calo ,"met_calo", logspace(30, 1, 150));
      book(_h_set_incl ,"set_incl", logspace(30, 1, sqrtS()/GeV/2));
      book(_h_set_calo ,"set_calo", logspace(30, 1, sqrtS()/GeV/2));
    }


    void analyze(const Event& event) {
      const double weight = 1.0;

      const MissingMomentum& mmincl = apply<MissingMomentum>(event, "InclMET");
      _h_met_incl->fill(mmincl.met()/GeV, weight);
      _h_set_incl->fill(mmincl.set()/GeV, weight);

      const MissingMomentum& mmcalo = apply<MissingMomentum>(event, "CaloMET");
      _h_met_calo->fill(mmcalo.met()/GeV, weight);
      _h_set_calo->fill(mmcalo.set()/GeV, weight);

    }


    void finalize() {
      normalize(_h_met_incl); normalize(_h_set_incl);
      normalize(_h_met_calo); normalize(_h_set_calo);
    }


  private:

    Histo1DPtr _h_met_incl, _h_set_incl;
    Histo1DPtr _h_met_calo, _h_set_calo;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_MET);


}
