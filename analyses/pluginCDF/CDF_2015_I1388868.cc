// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief CDF leading track underlying event at 300, 900 and 1960 GeV
  /// @author Orestes Tumbarell Aranda (Havana), Hannes Jung (DESY)
  class CDF_2015_I1388868 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CDF_2015_I1388868);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Energy selection
      int isqrts = -1;
      if (isCompatibleWithSqrtS(300)) {
        isqrts = 2;
      } else if (isCompatibleWithSqrtS(900)) {
        isqrts = 1;
      } else if (isCompatibleWithSqrtS(1960)) {
        isqrts = 0;
      } else {
        throw UserError("Unexpected sqrtS ! Only 300, 900, 1960 GeV is supported by CDF_2015_I1388868");
      }
      MSG_DEBUG("CDF Tevatron UE: running with " << sqrtS()/GeV);

      // Book projection
      const ChargedFinalState cfs(Cuts::abseta < 0.8 && Cuts::pT > 0.5*GeV);
      declare(cfs, "Tracks");

      // Book profile histos
      book(_NchgPDFden1,  8*isqrts+4,1,1);
      book(_NchgPMNden1,  8*isqrts+2,1,1);
      book(_NchgPMXden1,  8*isqrts+1,1,1);
      book(_NchgPden1,    8*isqrts+3,1,1);
      book(_PTsumPDFden1, 8*isqrts+8,1,1);
      book(_PTsumPMNden1, 8*isqrts+6,1,1);
      book(_PTsumPMXden1, 8*isqrts+5,1,1);
      book(_PTsumPden1,   8*isqrts+7,1,1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Require at least one track in the event with pT >= 0.5 GeV
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "Tracks");
      if (cfs.empty()) vetoEvent;
      const Particles trks = cfs.particlesByPt();

      // Get lead track
      const Particle p_lead = trks[0];
      const double philead = p_lead.phi();
      const double ptlead  = p_lead.pT();

      // Loop over tracks and compute variables
      double NchgP1 = 0, NchgP2 = 0, PTsumP1 = 0, PTsumP2 = 0;
      for (const Particle& p : trks) {

        // Region definition -- if not in transverse region, ignore
        const double dphi = mapAngle0To2Pi(p.phi() - philead);
        if (!inRange(dphi, PI/3, 2*PI/3) && !inRange(dphi, 4*PI/3, 5*PI/3)) continue;

        // Transverse region 1
        if (inRange(dphi, PI/3, 2*PI/3)) {
          NchgP1 += 1;
          PTsumP1 += p.pT();
        }
        // Transverse region 2
        else if (inRange(dphi, 4*PI/3, 5*PI/3)) {
          NchgP2 += 1;
          PTsumP2 += p.pT();
        }
      }

      // Calculate total variables
      const double NchgPtot = (NchgP1 + NchgP2)/2;
      const double NchgPmax = max(NchgP1,NchgP2);
      const double NchgPmin = min(NchgP1,NchgP2);
      const double PTsumPtot = (PTsumP1 + PTsumP2)/2;
      const double PTsumPmax = max(PTsumP1,PTsumP2);
      const double PTsumPmin = min(PTsumP1,PTsumP2);
      //
      const double PTsumPMXden = PTsumPmax/AREA;
      const double PTsumPMNden = PTsumPmin/AREA;
      const double NchgPMXden = NchgPmax/AREA;
      const double NchgPMNden = NchgPmin/AREA;
      //
      const double NchgPDFden = NchgPMXden - NchgPMNden;
      const double PTsumPDFden = PTsumPMXden - PTsumPMNden;

      // Fill histograms
      _NchgPden1  ->fill(ptlead/GeV, NchgPtot/AREA);
      _NchgPMXden1->fill(ptlead/GeV, NchgPmax/AREA);
      _NchgPMNden1->fill(ptlead/GeV, NchgPmin/AREA);
      _NchgPDFden1->fill(ptlead/GeV, NchgPDFden );
      _PTsumPden1  ->fill(ptlead/GeV, PTsumPtot/AREA);
      _PTsumPMXden1->fill(ptlead/GeV, PTsumPmax/AREA);
      _PTsumPMNden1->fill(ptlead/GeV, PTsumPmin/AREA);
      _PTsumPDFden1->fill(ptlead/GeV, PTsumPDFden );
    }

    //@}


    /// eta-phi area of the transverse region
    constexpr static double AREA = 2*0.8 * M_PI/3;

    /// Histograms
    Profile1DPtr _NchgPden1, _NchgPMXden1,_NchgPMNden1,_NchgPDFden1,_PTsumPden1,_PTsumPMXden1,_PTsumPMNden1,_PTsumPDFden1;

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(CDF_2015_I1388868);

}
