// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief D0 azimuthal correlation of jets widely separated in rapidity
  class D0_1996_S3324664 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    D0_1996_S3324664() : Analysis("D0_1996_S3324664")
    {    }


    /// @name Analysis methods
    //@{

    void init() {
      const FinalState fs;
      declare(fs, "FS");
      /// @todo Use correct jet algorithm
      declare(FastJets(fs, FastJets::D0ILCONE, 0.7), "ConeJets");

      book(_h_deta ,1, 1, 1);
      {Histo1DPtr tmp; _h_dphi.add(0.0, 2.0, book(tmp, 2, 1, 1));}
      {Histo1DPtr tmp; _h_dphi.add(2.0, 4.0, book(tmp, 2, 1, 2));}
      {Histo1DPtr tmp; _h_dphi.add(4.0, 6.0, book(tmp, 2, 1, 3));}
      book(_h_cosdphi_deta ,3, 1, 1);
    }


    void analyze(const Event& event) {
      const double weight = 1.0;
      Jets jets = apply<FastJets>(event, "ConeJets").jets(Cuts::Et > 20*GeV && Cuts::abseta<3, cmpMomByEt);

      if (jets.size() < 2) vetoEvent;

      FourMomentum minjet = jets[0].momentum();
      FourMomentum maxjet = jets[1].momentum();
      double mineta = minjet.eta();
      double maxeta = maxjet.eta();

      for (const Jet& jet : jets) {
        double eta = jet.eta();
        if (eta < mineta) {
          minjet = jet.momentum();
          mineta = eta;
        } else if (eta > maxeta) {
          maxjet = jet.momentum();
          maxeta = eta;
        }
      }

      if (minjet.Et() < 50*GeV && maxjet.Et() < 50.0*GeV) vetoEvent;

      double deta = maxjet.eta()-minjet.eta();
      double dphi = mapAngle0To2Pi(maxjet.phi()-minjet.phi());

      _h_deta->fill(deta, weight);
      _h_dphi.fill(deta, 1.0-dphi/M_PI, weight);
      _h_cosdphi_deta->fill(deta, cos(M_PI-dphi), weight);
    }


    void finalize() {
      // Normalised to #events
      normalize(_h_deta, 8830.); // fixed norm OK

      // Normalied to 1/(4pi)
      for (Histo1DPtr histo : _h_dphi.histos()) {
        normalize(histo, 1./(4.*M_PI));
      }

    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_deta;
    BinnedHistogram _h_dphi;
    Profile1DPtr _h_cosdphi_deta;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(D0_1996_S3324664);

}
