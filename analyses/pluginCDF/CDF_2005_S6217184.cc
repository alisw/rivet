// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/JetShape.hh"

namespace Rivet {


  /// @brief CDF Run II jet shape analysis
  ///
  /// @author Andy Buckley
  class CDF_2005_S6217184 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CDF_2005_S6217184);


    /// @name Analysis methods
    /// @{

    void init() {
      // Set up projections
      const FinalState fs(Cuts::abseta < 2);
      declare(fs, "FS");
      FastJets fj(fs, FastJets::CDFMIDPOINT, 0.7, JetAlg::Muons::ALL, JetAlg::Invisibles::ALL);
      declare(fj, "Jets");

      // Register a jet shape projection and histogram for each pT bin
      for (size_t i = 0; i < 6; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          const size_t k = i*3 + j;
          stringstream ss; ss << "JetShape" << k;
          const string pname = ss.str();
          _jsnames_pT[k] = pname;
          const JetShape jsp(fj, 0.0, 0.7, 7, PTEDGES[k], PTEDGES[k+1], 0.1, 0.7, RAPIDITY);
          declare(jsp, pname);
          book(_profhistRho_pT[k] ,i+1, 1, j+1);
          book(_profhistPsi_pT[k] ,6+i+1, 1, j+1);
        }
      }

      // Final histo
      book(_profhistPsi_vs_pT, 13, 1, 1, true);
    }


    /// Do the analysis
    void analyze(const Event& evt) {

      // Get jets and require at least one to pass pT and y cuts
      const Jets jets = apply<FastJets>(evt, "Jets")
        .jetsByPt(Cuts::ptIn(PTEDGES.front()*GeV, PTEDGES.back()*GeV) && Cuts::absrap < 0.7);
      MSG_DEBUG("Jet multiplicity before cuts = " << jets.size());
      if (jets.size() == 0) {
        MSG_DEBUG("No jets found in required pT and rapidity range");
        vetoEvent;
      }

      // Calculate and histogram jet shapes
      for (size_t ipt = 0; ipt < 18; ++ipt) {
        const JetShape& jsipt = apply<JetShape>(evt, _jsnames_pT[ipt]);
        for (size_t ijet = 0; ijet < jsipt.numJets(); ++ijet) {
          for (size_t rbin = 0; rbin < jsipt.numBins(); ++rbin) {
            const double r_rho = jsipt.rBinMid(rbin);
            MSG_DEBUG(ipt << " " << rbin << " (" << r_rho << ") " << jsipt.diffJetShape(ijet, rbin));
            /// @note Bin width Jacobian factor of 0.7/0.1 = 7 in the differential shapes plot
            _profhistRho_pT[ipt]->fill(r_rho/0.7, (0.7/0.1)*jsipt.diffJetShape(ijet, rbin));
            const double r_Psi = jsipt.rBinMax(rbin);
            _profhistPsi_pT[ipt]->fill(r_Psi/0.7, jsipt.intJetShape(ijet, rbin));
          }
        }
      }

    }


    // Finalize
    void finalize() {
      // Construct final 1-Psi(0.3/0.7) profile from Psi profiles
      for (size_t i = 0; i < PTEDGES.size()-1; ++i) {
        // Get entry for rad_Psi = 0.2 bin
        /// @todo Not a great handling of empty bins! Skip point, or set NaN values?!
        ProfileBin1D& bi = _profhistPsi_pT[i]->bin(2);
        const double y  = (bi.effNumEntries() > 0) ? bi.mean() : 0;
        const double ey = (bi.effNumEntries() > 1) ? bi.stdErr() : 0;
        _profhistPsi_vs_pT->point(i).setY(y, ey);
      }
    }

    /// @}


  private:

    /// Jet \f$ p_\perp\f$ bins.
    const array<double, 19> PTEDGES =
      {{ 37.0, 45.0, 55.0, 63.0, 73.0, 84.0, 97.0, 112.0, 128.0, 148.0,
         166.0, 186.0, 208.0, 229.0, 250.0, 277.0, 304.0, 340.0, 380.0 }};

    /// JetShape projection name for each \f$p_\perp\f$ bin.
    array<string, 18> _jsnames_pT;


    /// @name Histograms
    /// @{
    Profile1DPtr _profhistRho_pT[18];
    Profile1DPtr _profhistPsi_pT[18];
    Scatter2DPtr _profhistPsi_vs_pT;
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CDF_2005_S6217184, CDF_2005_I682179);

}
