// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief CMS strange particle spectra (Ks, Lambda, Cascade) in pp at 900 and 7000 GeV
  ///
  /// @author Kevin Stenson
  class CMS_2011_S8978280 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2011_S8978280);


    void init() {
      UnstableParticles ufs(Cuts::absrap < 2);
      declare(ufs, "UFS");
      int beamEnergy = -1;
      if (isCompatibleWithSqrtS(900.))  beamEnergy = 1;
      else if (isCompatibleWithSqrtS(7000.))  beamEnergy = 2;
      else {
        MSG_WARNING("Could not decipher beam energy. For rivet-merge set -a CMS_2011_S8978280:energy=OPT, where OPT is 900 or 7000 (GeV is implied).");
      }
      
      // Particle distributions versus rapidity and transverse momentum
      if (beamEnergy == 1){
        book(_h_dNKshort_dy  ,1, 1, 1);
        book(_h_dNKshort_dpT ,2, 1, 1);
        book(_h_dNLambda_dy  ,3, 1, 1);
        book(_h_dNLambda_dpT ,4, 1, 1);
        book(_h_dNXi_dy      ,5, 1, 1);
        book(_h_dNXi_dpT     ,6, 1, 1);
        //
        book(_h_LampT_KpT , 7, 1, 1);
        book(_h_XipT_LampT, 8, 1, 1);
        book(_h_Lamy_Ky   , 9, 1, 1);
        book(_h_Xiy_Lamy  , 10, 1, 1);

      } else if (beamEnergy == 2){
        book(_h_dNKshort_dy  ,1, 1, 2);
        book(_h_dNKshort_dpT ,2, 1, 2);
        book(_h_dNLambda_dy  ,3, 1, 2);
        book(_h_dNLambda_dpT ,4, 1, 2);
        book(_h_dNXi_dy      ,5, 1, 2);
        book(_h_dNXi_dpT     ,6, 1, 2);
        //
        book(_h_LampT_KpT , 7, 1, 2);
        book(_h_XipT_LampT, 8, 1, 2);
        book(_h_Lamy_Ky   , 9, 1, 2);
        book(_h_Xiy_Lamy  , 10, 1, 2);
      } else {
        MSG_WARNING("Could not initialize properly.");
      }
    }


    void analyze(const Event& event) {

      const UnstableParticles& parts = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : parts.particles()) {
        switch (p.abspid()) {
        case PID::K0S:
          _h_dNKshort_dy->fill(p.absrap());
          _h_dNKshort_dpT->fill(p.pT()/GeV);
          break;

        case PID::LAMBDA:
          // Lambda should not have Cascade or Omega ancestors since they should not decay. But just in case...
          if ( !( p.hasAncestor(3322) || p.hasAncestor(-3322) || p.hasAncestor(3312) || p.hasAncestor(-3312) || p.hasAncestor(3334) || p.hasAncestor(-3334) ) ) {
            _h_dNLambda_dy->fill(p.absrap());
            _h_dNLambda_dpT->fill(p.pT()/GeV);
          }
          break;

        case PID::XIMINUS:
          // Cascade should not have Omega ancestors since it should not decay.  But just in case...
          if ( !( p.hasAncestor(3334) || p.hasAncestor(-3334) ) ) {
            _h_dNXi_dy->fill(p.absrap());
            _h_dNXi_dpT->fill(p.pT()/GeV);
          }
          break;
        }

      }
    }


    void finalize() {
      divide(_h_dNLambda_dpT,_h_dNKshort_dpT, _h_LampT_KpT);
      divide(_h_dNXi_dpT,_h_dNLambda_dpT, _h_XipT_LampT);
      divide(_h_dNLambda_dy,_h_dNKshort_dy, _h_Lamy_Ky);
      divide(_h_dNXi_dy,_h_dNLambda_dy, _h_Xiy_Lamy);
      const double normpT = 1.0/sumOfWeights();
      const double normy = 0.5*normpT; // Accounts for using |y| instead of y
      scale(_h_dNKshort_dy, normy);
      scale(_h_dNKshort_dpT, normpT);
      scale(_h_dNLambda_dy, normy);
      scale(_h_dNLambda_dpT, normpT);
      scale(_h_dNXi_dy, normy);
      scale(_h_dNXi_dpT, normpT);
    }


  private:

    /// @name Particle distributions versus rapidity and transverse momentum
    /// @{
    Histo1DPtr _h_dNKshort_dy, _h_dNKshort_dpT, _h_dNLambda_dy, _h_dNLambda_dpT, _h_dNXi_dy, _h_dNXi_dpT;
    Scatter2DPtr _h_LampT_KpT, _h_XipT_LampT, _h_Lamy_Ky, _h_Xiy_Lamy;
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CMS_2011_S8978280, CMS_2011_I890166);

}
