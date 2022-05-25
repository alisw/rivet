#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// Drell-Yan dimuon absolute cross-sections in 800 GeV pp and pd collisions
  class NUSEA_2003_I613362 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(NUSEA_2003_I613362);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Projections
      const FinalState fs;
      declare(fs, "FS");
      Cut cut = Cuts::etaIn(-10.,10.);
      ZFinder zfinder(fs, cut, PID::MUON, 4.0*GeV, 100.0*GeV, 0.1, ZFinder::ClusterPhotons::NONE );
      declare(zfinder, "ZFinder");

      // Booking histograms
      // hydrogen d01-d16
      Histo1DPtr dummy;
      _hist_M_xF.add(-0.05, 0.05, book(dummy,1, 1, 1));
      _hist_M_xF.add( 0.05, 0.10, book(dummy,2, 1, 1));
      _hist_M_xF.add( 0.10, 0.15, book(dummy,3, 1, 1));
      _hist_M_xF.add( 0.15, 0.20, book(dummy,4, 1, 1));
      _hist_M_xF.add( 0.20, 0.25, book(dummy,5, 1, 1));
      _hist_M_xF.add( 0.25, 0.30, book(dummy,6, 1, 1));
      _hist_M_xF.add( 0.30, 0.35, book(dummy,7, 1, 1));
      _hist_M_xF.add( 0.35, 0.40, book(dummy,8, 1, 1));
      _hist_M_xF.add( 0.40, 0.45, book(dummy,9, 1, 1));
      _hist_M_xF.add( 0.45, 0.50, book(dummy,10, 1, 1));
      _hist_M_xF.add( 0.50, 0.55, book(dummy,11, 1, 1));
      _hist_M_xF.add( 0.55, 0.60, book(dummy,12, 1, 1));
      _hist_M_xF.add( 0.60, 0.65, book(dummy,13, 1, 1));
      _hist_M_xF.add( 0.65, 0.70, book(dummy,14, 1, 1));
      _hist_M_xF.add( 0.70, 0.75, book(dummy,15, 1, 1));
      _hist_M_xF.add( 0.75, 0.80, book(dummy,16, 1, 1));

      // deuterium d17-d32

      // hydrogen d40
      _hist_pT_M.add(4.2, 5.2, book(dummy,40, 1, 1));
      _hist_pT_M.add(5.2, 6.2, book(dummy,40, 1, 2));
      _hist_pT_M.add(6.2, 7.2, book(dummy,40, 1, 3));
      _hist_pT_M.add(7.2, 8.7, book(dummy,40, 1, 4));
      _hist_pT_M.add(10.85, 12.85, book(dummy,40, 1, 5));

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double sqrts_tol = 10.;
      if (!isCompatibleWithSqrtS(38.8, sqrts_tol)) {
        MSG_ERROR("Incorrect beam energy used: " << sqrtS()/GeV);
        throw Error("Unexpected sqrtS ! Only 38.8 GeV is supported");
      }

      // Muons
      const ZFinder& zfinder = applyProjection<ZFinder>(event, "ZFinder");
      if (zfinder.particles().size() <= 0) vetoEvent;

      double Zmass = zfinder.bosons()[0].momentum().mass()/GeV;
      double Zpt   = zfinder.bosons()[0].momentum().pT()/GeV;
      double Zpl   = zfinder.bosons()[0].momentum().pz()/GeV;
      double ZE    = zfinder.bosons()[0].momentum().E();
      double xF = 2.*Zpl/sqrtS();

      // Filling dimuon mass in bins of xF
      _hist_M_xF.fill(xF, Zmass/GeV, pow(Zmass,3));

      // Filling pT in bins of Zmass
      if ( xF > -0.05 && xF <= 0.15 ) {
        // Include here all factors which are run-dependent for later scaling
        if (Zpt > 0) _hist_pT_M.fill(Zmass,Zpt, 1./2./Zpt*2.*ZE/sqrtS());
      }

      MSG_DEBUG("Dimuon pT = "<< Zpt<<"   Dimuon E = ");
      MSG_DEBUG("DiMuon mass " << Zmass/GeV);
      MSG_DEBUG("DiMuon pT "<<Zpt);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // xf bin width = 0.2, x-section in picobarn
      double scalefactor=crossSection()/picobarn/(sumOfWeights() * M_PI *0.2 );
      _hist_pT_M.scale(scalefactor, this);

      // x-section is quoted in nanobarn
      _hist_M_xF.scale(crossSection()/nanobarn/sumOfWeights(), this);
    }

    /// @}


  private:

    /// @name Histograms
    ///@{
    BinnedHistogram _hist_pT_M, _hist_M_xF;
    Histo1DPtr  _h_m_DiMuon ;
    Histo1DPtr  _h_pT_DiMuon;
    Histo1DPtr  _h_eta_DiMuon;
    Histo1DPtr  _h_y_DiMuon;
    Histo1DPtr  _h_phi_DiMuon;
    Histo1DPtr _h_yDiff_DiMuon;
    Histo1DPtr _h_dPhi_DiMuon;
    Histo1DPtr _h_xF_DiMuon;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(NUSEA_2003_I613362);

}
