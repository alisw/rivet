// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Inclusive jets at 8 TeV
  class ATLAS_2017_I1604271 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2017_I1604271);

    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs;
      declare(fs,"FinalState");
      FastJets fj04(fs, FastJets::ANTIKT, 0.4);
      FastJets fj06(fs, FastJets::ANTIKT, 0.6);
      fj04.useInvisibles();
      declare(fj04, "AntiKT04");
      fj06.useInvisibles();
      declare(fj06, "AntiKT06");

      // |y| and ystar bins
      const int nybins = 6;
      double ybins[nybins+1] = { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};

      // Book histograms
      // pT histograms
      for (size_t i = 0; i < nybins; ++i){ // loop over |y| bins
        {Histo1DPtr tmp; _pThistograms6.add(ybins[i], ybins[i+1], book(tmp, i+1,1,1));}
        {Histo1DPtr tmp; _pThistograms4.add(ybins[i], ybins[i+1], book(tmp, i+7,1,1));}
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const Jets& kt4Jets = applyProjection<FastJets>(event, "AntiKT04").jetsByPt(Cuts::pT>70*GeV && Cuts::absrap < 3.0);
      const Jets& kt6Jets = applyProjection<FastJets>(event, "AntiKT06").jetsByPt(Cuts::pT>70*GeV && Cuts::absrap < 3.0);

      int nJets4 = kt4Jets.size();
      int nJets6 = kt6Jets.size();

      // Inclusive jet selection
      for(int ijet=0;ijet<nJets4;++ijet){ // loop over jets
        FourMomentum jet = kt4Jets[ijet].momentum();
        // pT selection
        const double absy = jet.absrap();
        _pThistograms4.fill(absy,jet.pt()/GeV);
      }

      for(int ijet=0;ijet<nJets6;++ijet){ // loop over jets
        FourMomentum jet = kt6Jets[ijet].momentum();
        // pT selection
        const double absy = jet.absrap();
        _pThistograms6.fill(absy,jet.pt()/GeV);
      }


    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double xs_pb( crossSection() / picobarn );
      const double sumW( sumOfWeights() );
      const double xs_norm_factor( 0.5*xs_pb / sumW );
     
      MSG_DEBUG( "Cross-Section/pb     : " << xs_pb       );
      MSG_DEBUG( "ZH                   : " << crossSectionPerEvent()/ picobarn);
      MSG_DEBUG( "Sum of weights       : " << sumW        );
      MSG_DEBUG( "nEvents              : " << numEvents() );
      _pThistograms4.scale(xs_norm_factor, this);
      _pThistograms6.scale(xs_norm_factor, this);

    }

  private:

    // The inclusive pT spectrum for akt4 jets
    BinnedHistogram _pThistograms4;
    // The inclusive pT spectrum for akt6 jets
    BinnedHistogram _pThistograms6;

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1604271);

}
