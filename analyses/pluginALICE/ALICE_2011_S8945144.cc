// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class ALICE_2011_S8945144 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2011_S8945144);


    void init() {
      const ChargedFinalState cfs((Cuts::etaIn(-15, 15)));
      declare(cfs, "CFS");

      book(_histPtPions        ,"d01-x01-y01");
      book(_histPtAntiPions    ,"d01-x01-y02");
      book(_histPtKaons        ,"d02-x01-y01");
      book(_histPtAntiKaons    ,"d02-x01-y02");
      book(_histPtProtons      ,"d03-x01-y01");
      book(_histPtAntiProtons  ,"d03-x01-y02");
      book(_histAveragePt      ,"d04-x01-y01");
    }


    void analyze(const Event& event) {
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      for (const Particle& p : cfs.particles()) {
        if(p.absrap()<0.5) {
          switch (p.pid()) {
            case 211:
              _histPtPions->fill(p.pT()/GeV);
              _histAveragePt->fill(p.mass()/GeV, p.pT()/GeV);
              break;
            case -211:
              _histPtAntiPions->fill(p.pT()/GeV);
              _histAveragePt->fill(p.mass()/GeV, p.pT()/GeV);
              break;
            case 2212:
              if ( !(p.hasAncestor(3322) ||                             // Xi0
                     p.hasAncestor(3122) || p.hasAncestor(-3122) ||     // Lambda
                     p.hasAncestor(3222) || p.hasAncestor(-3222) ||     // Sigma+/-
                     p.hasAncestor(3312) || p.hasAncestor(-3312) ) ) {  // Xi-/+
                _histPtProtons->fill(p.pT()/GeV);
                _histAveragePt->fill(p.mass()/GeV, p.pT()/GeV);
              }
              break;
            case -2212:
              if ( !(p.hasAncestor(3322) ||                             // Xi0
                     p.hasAncestor(3122) || p.hasAncestor(-3122) ||     // Lambda
                     p.hasAncestor(3222) || p.hasAncestor(-3222) ||     // Sigma+/-
                     p.hasAncestor(3312) || p.hasAncestor(-3312) ) ) {  // Xi-/+
                _histPtAntiProtons->fill(p.pT()/GeV);
                _histAveragePt->fill(p.mass()/GeV, p.pT()/GeV);
              }
              break;
            case 321:
              _histPtKaons->fill(p.pT()/GeV);
              _histAveragePt->fill(p.mass()/GeV, p.pT()/GeV);
              break;
            case -321:
              _histPtAntiKaons->fill(p.pT()/GeV);
              _histAveragePt->fill(p.mass()/GeV, p.pT()/GeV);
              break;
          }
        }
      }
    }


    void finalize() {
      scale(_histPtPions,       1./sumOfWeights());
      scale(_histPtProtons,     1./sumOfWeights());
      scale(_histPtKaons,       1./sumOfWeights());
      scale(_histPtAntiPions,   1./sumOfWeights());
      scale(_histPtAntiProtons, 1./sumOfWeights());
      scale(_histPtAntiKaons,   1./sumOfWeights());
    }


  private:

    Histo1DPtr _histPtPions;
    Histo1DPtr _histPtProtons;
    Histo1DPtr _histPtKaons;
    Histo1DPtr _histPtAntiPions;
    Histo1DPtr _histPtAntiProtons;
    Histo1DPtr _histPtAntiKaons;
    Profile1DPtr _histAveragePt;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(ALICE_2011_S8945144, ALICE_2011_I885104);

}
