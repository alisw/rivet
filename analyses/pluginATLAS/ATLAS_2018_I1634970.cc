// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class ATLAS_2018_I1634970 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2018_I1634970);

    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs;
      declare(fs,"FinalState");
      FastJets fj04(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::ALL, JetAlg::Invisibles::DECAY);
      declare(fj04, "AntiKT04");

      // |y| and ystar bins
      const int nybins          = 6;
      double ybins[nybins+1]     = { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};
      double ystarbins[nybins+1] = { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};

      // Book histograms
      // pT histograms
      for(size_t i=0;i<nybins;++i){// loop over |y| bins
        {Histo1DPtr tmp; _pThistograms.add(ybins[i], ybins[i+1], book(tmp, i+1,1,1));}
      }
      // mjj histograms
      for(size_t i=0;i<nybins;++i){// loop over ystar bins
        {Histo1DPtr tmp; _mjjhistograms.add(ystarbins[i], ystarbins[i+1], book(tmp, i+7,1,1));}
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const Jets& kt4Jets = apply<FastJets>(event, "AntiKT04").jetsByPt(Cuts::pT > 75*GeV && Cuts::absrap < 3.0);

      int nJets = kt4Jets.size();

      // Inclusive jet selection
      for(int ijet=0;ijet<nJets;++ijet){ // loop over jets
        FourMomentum jet = kt4Jets[ijet].momentum();
        // pT selection
        if(jet.pt()>100.0*GeV){
          // Fill distribution
          const double absy = jet.absrap();
          _pThistograms.fill(absy,jet.pt()/GeV);
        }
      }

      // Dijet selection
      if(nJets > 1){ // skip events with less than 2 jets passing pT>75GeV and |y|<3.0 cuts
        FourMomentum jet0  = kt4Jets[0].momentum(); 
        FourMomentum jet1  = kt4Jets[1].momentum();
        const double rap0  = jet0.rapidity();
        const double rap1  = jet1.rapidity();
        const double ystar = fabs(rap0-rap1)/2;
        const double mass  = (jet0 + jet1).mass(); 
        const double HT2   = jet0.pt()+jet1.pt();
        if(HT2>200*GeV && ystar<3.0){
          // Fill distribution
          _mjjhistograms.fill(ystar,mass/GeV);
        }
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
      _pThistograms.scale(xs_norm_factor, this);
      _mjjhistograms.scale(crossSectionPerEvent()/picobarn, this);

    }

  private:

    // The inclusive pT spectrum for akt4 jets
    BinnedHistogram _pThistograms;
    // The dijet mass spectrum for akt4 jets
    BinnedHistogram _mjjhistograms;

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2018_I1634970);

}
