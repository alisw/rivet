// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class ALICE_2015_I1357424 : public Analysis {
  public:

    ALICE_2015_I1357424()
      : Analysis("ALICE_2015_I1357424")
    {}


  public:

    void init() {
      const ChargedFinalState cfs(Cuts::absrap<0.5);
      declare(cfs, "CFS");
      //
      // plots from the paper
      book(_histPtPions          ,"d01-x01-y01");    // pions
      book(_histPtKaons          ,"d01-x01-y02");    // kaons
      book(_histPtProtons        ,"d01-x01-y03");    // protons
      book(_histPtKtoPi          ,"d02-x01-y01");  // K to pi ratio 
      book(_histPtPtoPi          ,"d03-x01-y01");  // p to pi ratio
      //
      // temp histos for ratios
      book(_histPtPionsR1        ,"TMP/pT_pi1", refData(2, 1, 1)); // pi histo compatible with more restricted kaon binning
      book(_histPtPionsR2        ,"TMP/pT_pi2", refData(3, 1, 1)); // pi histo compatible with more restricted proton binning
      book(_histPtKaonsR         ,"TMP/pT_K",   refData(2, 1, 1)); // K histo with more restricted binning
      book(_histPtProtonsR       ,"TMP/pT_p",   refData(3, 1, 1)); // p histo with more restricted binning
    }


    void analyze(const Event& event) {
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      for (const Particle& p : cfs.particles()) {
	// protections against mc generators decaying long-lived particles
	if ( !(p.hasAncestor(310)  || p.hasAncestor(-310)  ||     // K0s
	       p.hasAncestor(130)  || p.hasAncestor(-130)  ||     // K0l
	       p.hasAncestor(3322) || p.hasAncestor(-3322) ||     // Xi0
	       p.hasAncestor(3122) || p.hasAncestor(-3122) ||     // Lambda
	       p.hasAncestor(3222) || p.hasAncestor(-3222) ||     // Sigma+/-
	       p.hasAncestor(3312) || p.hasAncestor(-3312) ||     // Xi-/+ 
	   p.hasAncestor(3334) || p.hasAncestor(-3334) ))     // Omega-/+     
        {  
	  switch (abs(p.pid())) {
	  case 211: // pi+
	    _histPtPions->fill(p.pT()/GeV);
	    _histPtPionsR1->fill(p.pT()/GeV);
	    _histPtPionsR2->fill(p.pT()/GeV);
	    break;
	  case 2212: // proton
	    _histPtProtons->fill(p.pT()/GeV);
	    _histPtProtonsR->fill(p.pT()/GeV);
	    break;
	  case 321: // K+
	    _histPtKaons->fill(p.pT()/GeV);
	    _histPtKaonsR->fill(p.pT()/GeV);
	    break;
	  } // particle switch
	} // primary pi, K, p only
      } // particle loop
    }    

    void finalize() {
      divide(_histPtKaonsR,   _histPtPionsR1, _histPtKtoPi);
      divide(_histPtProtonsR, _histPtPionsR2, _histPtPtoPi);

      scale(_histPtPions,       1./sumOfWeights());
      scale(_histPtProtons,     1./sumOfWeights());
      scale(_histPtKaons,       1./sumOfWeights());
    }


  private:

    Histo1DPtr _histPtPions;
    Histo1DPtr _histPtProtons;
    Histo1DPtr _histPtKaons;

    Histo1DPtr _histPtPionsR1;
    Histo1DPtr _histPtPionsR2;
    Histo1DPtr _histPtProtonsR;
    Histo1DPtr _histPtKaonsR;

    Scatter2DPtr _histPtKtoPi;
    Scatter2DPtr _histPtPtoPi;
  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2015_I1357424);

}
