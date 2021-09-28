
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {
 
  /// CMS azimuthal decorrelations in back-to-back dijet events at 13 TeV 
  class CMS_2019_I1719955 : public Analysis {
  public:
    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2019_I1719955);
 

    /// Book projections and histograms 
    void init() {
      const FinalState fs;
      declare(FastJets(fs, FastJets::ANTIKT, 0.4), "ANTIKT");
       Histo1DPtr dummy;

      _h_deltaPhi_2J.add( 200.,  300., book(dummy,1,1,1));
      _h_deltaPhi_2J.add( 300.,  400., book(dummy,2,1,1));
      _h_deltaPhi_2J.add( 400.,  500., book(dummy,3,1,1));
      _h_deltaPhi_2J.add( 500.,  600., book(dummy,4,1,1));
      _h_deltaPhi_2J.add( 600.,  700., book(dummy,5,1,1));
      _h_deltaPhi_2J.add( 700.,  800., book(dummy,6,1,1));
      _h_deltaPhi_2J.add( 800., 1000., book(dummy,7,1,1));
      _h_deltaPhi_2J.add( 1000.,1200., book(dummy,8,1,1));
      _h_deltaPhi_2J.add( 1200.,4000., book(dummy,9,1,1));

      _h_deltaPhi_3J.add( 200.,  300., book(dummy,10,1,1));
      _h_deltaPhi_3J.add( 300.,  400., book(dummy,11,1,1));
      _h_deltaPhi_3J.add( 400.,  500., book(dummy,12,1,1));
      _h_deltaPhi_3J.add( 500.,  600., book(dummy,13,1,1));
      _h_deltaPhi_3J.add( 600.,  700., book(dummy,14,1,1));
      _h_deltaPhi_3J.add( 700.,  800., book(dummy,15,1,1));
      _h_deltaPhi_3J.add( 800., 1000., book(dummy,16,1,1));
      _h_deltaPhi_3J.add( 1000.,1200., book(dummy,17,1,1));
      _h_deltaPhi_3J.add( 1200.,4000., book(dummy,18,1,1));
    }

    /// Per-event analysis
    void analyze(const Event& event) {
      const Jets& jets = applyProjection<JetAlg>(event, "ANTIKT").jetsByPt(Cuts::absrap < 5. && Cuts::pT > 100*GeV);
      const Jets& lowjets = applyProjection<JetAlg>(event, "ANTIKT").jetsByPt(Cuts::absrap < 2.5 && Cuts::pT > 30*GeV);
      if (jets.size() < 2) vetoEvent;
      if (jets[0].absrap() > 2.5 || jets[1].absrap() > 2.5) vetoEvent;

      const double dphi = 180./M_PI*deltaPhi(jets[0].phi(), jets[1].phi());
      _h_deltaPhi_2J.fill(jets[0].pT(), dphi);    
      if (lowjets.size() > 2) _h_deltaPhi_3J.fill(jets[0].pT(), dphi);
    }

    /// Scale histograms
    void finalize() {
      int region_ptmax_2J = 0;
      double norm_finalize[9];
      for (Histo1DPtr histo_2J : _h_deltaPhi_2J.histos()) {
        norm_finalize[region_ptmax_2J] = histo_2J->integral();
        if (norm_finalize[region_ptmax_2J] != 0) scale(histo_2J, 1.0/norm_finalize[region_ptmax_2J]);
        region_ptmax_2J++;  
      }
      int region_ptmax_3J = 0;         
      for (Histo1DPtr histo_3J : _h_deltaPhi_3J.histos()) {       
        if (norm_finalize[region_ptmax_3J] != 0) scale(histo_3J, 1.0/norm_finalize[region_ptmax_3J]);
        region_ptmax_3J++;
      }
    }


  private:

    BinnedHistogram _h_deltaPhi_2J;
    BinnedHistogram _h_deltaPhi_3J;

  };
 

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2019_I1719955);

}
