// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
//#include "Rivet/ParticleName.hh"

namespace Rivet {

  


  /// Underlying event activity in the Drell-Yan process at 13 TeV
  class CMS_2017_I1635889 : public Analysis {
  public:

    /// Constructor
    CMS_2017_I1635889()
      : Analysis("CMS_2017_I1635889")
    {   }


    /// Initialization
    void init() {

      /// @note Using a bare muon Z (but with a clustering radius!?)
      Cut cut = Cuts::abseta < 2.4 && Cuts::pT > 10*GeV;
      ZFinder zfinder(FinalState(), cut, PID::MUON, 81*GeV, 101*GeV, 0.2, ZFinder::NOCLUSTER);
      addProjection(zfinder, "ZFinder");

      ChargedFinalState cfs(zfinder.remainingFinalState());
      addProjection(cfs, "cfs");

      _h_Nchg_towards_pTmumu                 = bookProfile1D(1, 1, 1);
      _h_Nchg_transverse_pTmumu              = bookProfile1D(2, 1, 1);
      _h_Nchg_away_pTmumu                    = bookProfile1D(3, 1, 1);
      _h_pTsum_towards_pTmumu                = bookProfile1D(4, 1, 1);
      _h_pTsum_transverse_pTmumu             = bookProfile1D(5, 1, 1);
      _h_pTsum_away_pTmumu                   = bookProfile1D(6, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      const ZFinder& zfinder = applyProjection<ZFinder>(event, "ZFinder");

      if (zfinder.bosons().size() != 1) vetoEvent;
      if (zfinder.constituents()[0].pT()<20 && zfinder.constituents()[1].pT()<20)vetoEvent;
      //std::cout<<"pt[0] = "<<zfinder.constituents()[0].pT()<<"pt[1] = "<<zfinder.constituents()[1].pT()<<std::endl;
      double Zpt = zfinder.bosons()[0].pT()/GeV;
      double Zphi = zfinder.bosons()[0].phi();
      //double Zmass = zfinder.bosons()[0].mass()/GeV;

     Particles particles = applyProjection<ChargedFinalState>(event, "cfs").particlesByPt(Cuts::pT > 0.5*GeV && Cuts::abseta <2.0);

      int nTowards = 0;
      int nTransverse = 0;
      int nAway = 0;
      double ptSumTowards = 0;
      double ptSumTransverse = 0;
      double ptSumAway = 0;

      foreach (const Particle& p, particles) {
        double dphi = fabs(deltaPhi(Zphi, p.phi()));
        double pT = p.pT();

        if ( dphi < M_PI/3 ) {
          nTowards++;
          ptSumTowards += pT;
        } else if ( dphi < 2.*M_PI/3 ) {
          nTransverse++;
          ptSumTransverse += pT;
        } else {
          nAway++;
          ptSumAway += pT;
        }

      } // Loop over particles


      const double area = 8./3.*M_PI;
        _h_Nchg_towards_pTmumu->         fill(Zpt, 1./area * nTowards, weight);
        _h_Nchg_transverse_pTmumu->      fill(Zpt, 1./area * nTransverse, weight);
        _h_Nchg_away_pTmumu->            fill(Zpt, 1./area * nAway, weight);
        _h_pTsum_towards_pTmumu->        fill(Zpt, 1./area * ptSumTowards, weight);
        _h_pTsum_transverse_pTmumu->     fill(Zpt, 1./area * ptSumTransverse, weight);
        _h_pTsum_away_pTmumu->           fill(Zpt, 1./area * ptSumAway, weight);


    }


    /// Normalise histograms etc., after the run
    void finalize() {
    }


  private:


    /// @name Histogram objects
    //@{
    Profile1DPtr _h_Nchg_towards_pTmumu;
    Profile1DPtr _h_Nchg_transverse_pTmumu;
    Profile1DPtr _h_Nchg_away_pTmumu;
    Profile1DPtr _h_pTsum_towards_pTmumu;
    Profile1DPtr _h_pTsum_transverse_pTmumu;
    Profile1DPtr _h_pTsum_away_pTmumu;
    //@}

  };


  // Hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2017_I1635889);

}
