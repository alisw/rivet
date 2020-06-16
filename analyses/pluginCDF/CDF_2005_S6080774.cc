// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include <array>

namespace Rivet {


  /// @brief CDF diff cross-sections for prompt di-photon production
  class CDF_2005_S6080774 : public Analysis {
  public:

    /// Constructor
    CDF_2005_S6080774()
      : Analysis("CDF_2005_S6080774")
    {    }


    /// @name Analysis methods
    //@{

    void init() {
      FinalState fs;
      declare(fs, "FS");

      IdentifiedFinalState ifs(Cuts::abseta < 0.9 && Cuts::pT > 13*GeV);
      ifs.acceptId(PID::PHOTON);
      declare(ifs, "IFS");

      for (size_t yAxisId=0; yAxisId<4; ++yAxisId) {
        book(_h_m_PP[yAxisId],    1, 1, yAxisId + 1);
        book(_h_pT_PP[yAxisId],   2, 1, yAxisId + 1);
        book(_h_dphi_PP[yAxisId], 3, 1, yAxisId + 1);
      }
    }


    void analyze(const Event& event) {
      Particles photons = apply<IdentifiedFinalState>(event, "IFS").particlesByPt();
      if (photons.size() < 2 || photons[0].pT() < 14.0*GeV) {
        vetoEvent;
      }

      // Isolate photons with ET_sum in cone
      Particles isolated_photons;
      Particles fs = apply<FinalState>(event, "FS").particles();
      for (const Particle& photon : photons) {
        FourMomentum mom_in_cone;
        double eta_P = photon.eta();
        double phi_P = photon.phi();
        for (const Particle& p : fs) {
          if (deltaR(eta_P, phi_P, p.eta(), p.phi()) < 0.4) {
            mom_in_cone += p.momentum();
          }
        }
        if (mom_in_cone.Et()-photon.Et() < 1.0*GeV) {
          isolated_photons.push_back(photon);
        }
      }

      if (isolated_photons.size() != 2) {
        vetoEvent;
      }

      FourMomentum mom_PP = isolated_photons[0].momentum() + isolated_photons[1].momentum();
      for (size_t i=0; i<4; ++i) {
        _h_m_PP[i]->fill(mom_PP.mass());
        _h_pT_PP[i]->fill(mom_PP.pT());
        _h_dphi_PP[i]->fill(mapAngle0ToPi(isolated_photons[0].phi()-
                                          isolated_photons[1].phi())/M_PI);
      }
    }


    void finalize() {
      for (size_t i=0; i<4; ++i) {
        scale(_h_m_PP[i], crossSection()/sumOfWeights());
        scale(_h_pT_PP[i], crossSection()/sumOfWeights());
        scale(_h_dphi_PP[i], crossSection()/M_PI/sumOfWeights());
      }
    }

    //@}


  private:

    /// @name Histograms
    //@{
    std::array<Histo1DPtr,4> _h_m_PP;
    std::array<Histo1DPtr,4> _h_pT_PP;
    std::array<Histo1DPtr,4> _h_dphi_PP;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CDF_2005_S6080774);

}
