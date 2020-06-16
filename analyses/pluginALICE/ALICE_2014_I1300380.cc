// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
namespace Rivet {


  class ALICE_2014_I1300380 : public Analysis {
  public:

    ALICE_2014_I1300380()
      : Analysis("ALICE_2014_I1300380")
    {}


  public:

    void init() {
      const UnstableParticles cfs(Cuts::absrap<0.5);
      declare(cfs, "CFS");

      // Plots from the paper
      book(_histPtSigmaStarPlus        ,"d01-x01-y01");    // Sigma*+
      book(_histPtSigmaStarMinus       ,"d01-x01-y02");    // Sigma*- 
      book(_histPtSigmaStarPlusAnti    ,"d01-x01-y03");    // anti Sigma*-
      book(_histPtSigmaStarMinusAnti   ,"d01-x01-y04");    // anti Sigma*+
      book(_histPtXiStar               ,"d02-x01-y01");    // 0.5 * (xi star + anti xi star)
      book(_histAveragePt              ,"d03-x01-y01");  // <pT> profile
    }


    void analyze(const Event& event) {
      const UnstableParticles& cfs = apply<UnstableFinalState>(event, "CFS");
      for (const Particle& p : cfs.particles()) {
	// protections against mc generators decaying long-lived particles
	if ( !(p.hasAncestor(310)  || p.hasAncestor(-310)   || // K0s
	       p.hasAncestor(130)  || p.hasAncestor(-130)   ||     // K0l
	       p.hasAncestor(3322) || p.hasAncestor(-3322)  ||     // Xi0
	       p.hasAncestor(3122) || p.hasAncestor(-3122)  ||     // Lambda
	       p.hasAncestor(3222) || p.hasAncestor(-3222)  ||     // Sigma+/-
	       p.hasAncestor(3312) || p.hasAncestor(-3312)  ||     // Xi-/+
	       p.hasAncestor(3334) || p.hasAncestor(-3334)  ))     // Omega-/+     
	{   
	  int aid = p.abspid();
	  if (aid == 211  || // pi+ 
            aid == 321  || // K+
            aid == 313  || // K*(892)0
            aid == 2212 || // proton
            aid == 333  )  // phi(1020)
          {
	    _histAveragePt->fill(p.mass()/GeV, p.pT()/GeV);
	  }
	} // end if "rejection of long-lived particles"
      
      
        switch (p.pid()) {
	  case 3224:
	    _histPtSigmaStarPlus->fill(p.pT()/GeV);
	    _histAveragePt->fill(p.mass()/GeV, p.pT()/GeV);
	    break;
	  case -3224:
	    _histPtSigmaStarPlusAnti->fill(p.pT()/GeV);
	    _histAveragePt->fill(p.mass()/GeV, p.pT()/GeV);
	    break;
	  case 3114:
	    _histPtSigmaStarMinus->fill(p.pT()/GeV);
	    _histAveragePt->fill(p.mass()/GeV, p.pT()/GeV);
	    break;
	  case -3114:
	    _histPtSigmaStarMinusAnti->fill(p.pT()/GeV);
	    _histAveragePt->fill(p.mass()/GeV, p.pT()/GeV);
	    break;
	  case 3324:
	    _histPtXiStar->fill(p.pT()/GeV);
	    _histAveragePt->fill(p.mass()/GeV, p.pT()/GeV);
	    break;
	  case -3324:
	    _histPtXiStar->fill(p.pT()/GeV);
	    _histAveragePt->fill(p.mass()/GeV, p.pT()/GeV);
	    break;
	  case 3312:
	    _histAveragePt->fill(p.mass()/GeV, p.pT()/GeV);
	    break;
	  case -3312:
	    _histAveragePt->fill(p.mass()/GeV, p.pT()/GeV);
	    break;
	  case 3334:
	    _histAveragePt->fill(p.mass()/GeV, p.pT()/GeV);
	    break;
	  case -3334:
	    _histAveragePt->fill(p.mass()/GeV, p.pT()/GeV);
	    break;
        }
      }
    }


    void finalize() {
      scale(_histPtSigmaStarPlus,       1./sumOfWeights());
      scale(_histPtSigmaStarPlusAnti,   1./sumOfWeights());
      scale(_histPtSigmaStarMinus,      1./sumOfWeights());
      scale(_histPtSigmaStarMinusAnti,  1./sumOfWeights());
      scale(_histPtXiStar,              1./sumOfWeights()/ 2.);
    }


  private:
    // plots from the paper
    Histo1DPtr   _histPtSigmaStarPlus;
    Histo1DPtr   _histPtSigmaStarPlusAnti;
    Histo1DPtr   _histPtSigmaStarMinus;
    Histo1DPtr   _histPtSigmaStarMinusAnti;
    Histo1DPtr   _histPtXiStar;
    Profile1DPtr _histAveragePt;
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2014_I1300380);

}
