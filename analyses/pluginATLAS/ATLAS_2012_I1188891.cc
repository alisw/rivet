// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Particle.hh"

namespace Rivet {


  class ATLAS_2012_I1188891 : public Analysis {
  public:


    ATLAS_2012_I1188891()
      : Analysis("ATLAS_2012_I1188891")
    {    }


    void init() {

      const FinalState fs;
      FastJets fj04(fs,  FastJets::ANTIKT, 0.4);
      declare(fj04, "AntiKT04");

      string histotitle[7]={"BBfraction","BCfraction","CCfraction","BUfraction","CUfraction","UUfraction","Total"};
      for (int i = 0 ; i < 7 ; i++){
        book(_h_temp[i] ,"TMP/"+histotitle[i],refData(1,1,1));
        if (i < 6) {
          book(_h_results[i], i+1, 1, 1);
        }
      }
    }


    void analyze(const Event& event) {

      double weight    = 1.0;
      double weight100 = 1.0 * 100.;  //to get results in %

      //keeps jets with pt>20 geV and ordered in decreasing pt
      Jets jetAr = apply<FastJets>(event, "AntiKT04").jetsByPt(20*GeV);

      int flav[2]={1,1};
      vector<FourMomentum> leadjets;

      //get b/c-hadrons
      vector<ConstGenParticlePtr> B_hadrons, C_hadrons;
      vector<ConstGenParticlePtr> allParticles = HepMCUtils::particles(event.genEvent());
      for (size_t i = 0; i < allParticles.size(); i++) {
        ConstGenParticlePtr p = allParticles.at(i);
        if(p->momentum().perp()*GeV < 5) continue;
        if ( (Rivet::PID::isHadron ( p->pdg_id() ) &&
              Rivet::PID::hasBottom( p->pdg_id() )    ) ) {
          B_hadrons.push_back(p);
        }
        if ( (Rivet::PID::isHadron( p->pdg_id() ) &&
              Rivet::PID::hasCharm( p->pdg_id() )    ) ) {
          C_hadrons.push_back(p);
        }
      }

      //select dijet
      for (const Jet& jet : jetAr) {

        const double pT   = jet.pT();
        const double absy = jet.absrap();

        bool isBjet = jet.bTagged();
        bool isCjet = jet.cTagged();

        int jetflav=1;
        if      (isBjet)jetflav=5;
        else if (isCjet)jetflav=4;

        if (absy <= 2.1 && leadjets.size() < 2) {

          if (pT > 500*GeV) continue;
          if ((leadjets.empty() && pT < 40*GeV) || pT < 20*GeV) continue;

          leadjets.push_back(jet.momentum());

          if (leadjets.size()==1) flav[0] = jetflav;
          if (leadjets.size()==2) flav[1] = jetflav;
        }
      }

      if (leadjets.size() < 2) vetoEvent;


      double pBinsLJ[7] = {40.,60.,80.,120.,160.,250.,500.};
      int    iPBinLJ = -1;

      for (int k = 0 ; k < 7 ; k++) {
        if (leadjets[0].pT() > pBinsLJ[k]*GeV) iPBinLJ=k;
        else break;
      }

      bool c_ljpt  = (iPBinLJ != -1);
      bool c_nljpt = leadjets[1].pT() > 20*GeV;
      bool c_dphi  = fabs( deltaPhi(leadjets[0],leadjets[1]) ) > 2.1;
      bool isDijet = c_ljpt & c_nljpt & c_dphi;
      if (!isDijet) vetoEvent;

      _h_temp[6]->fill(leadjets[0].pT(), weight);

      if (flav[0]==5 && flav[1]==5)                                  // BB dijet
        _h_temp[0]->fill(leadjets[0].pT(), weight100);

      if ((flav[0]==5 && flav[1]==4) || (flav[0]==4 && flav[1]==5))  // BC dijet
        _h_temp[1]->fill(leadjets[0].pT(), weight100);

      if (flav[0]==4 && flav[1]==4)                                  // CC dijet
        _h_temp[2]->fill(leadjets[0].pT(), weight100);

      if ((flav[0]==5 && flav[1]==1) || (flav[0]==1 && flav[1]==5))  // B-light dijet
        _h_temp[3]->fill(leadjets[0].pT(), weight100);

      if ((flav[0]==4 && flav[1]==1) || (flav[0]==1 && flav[1]==4))  // C-light dijet
        _h_temp[4]->fill(leadjets[0].pT(), weight100);

      if (flav[0]==1 && flav[1]==1)                                  // light-light dijet
        _h_temp[5]->fill(leadjets[0].pT(), weight100);
    }


    void finalize() {
      divide(_h_temp[0], _h_temp[6], _h_results[0]);
      divide(_h_temp[1], _h_temp[6], _h_results[1]);
      divide(_h_temp[2], _h_temp[6], _h_results[2]);
      divide(_h_temp[3], _h_temp[6], _h_results[3]);
      divide(_h_temp[4], _h_temp[6], _h_results[4]);
      divide(_h_temp[5], _h_temp[6], _h_results[5]);
    }

  private:

    Histo1DPtr   _h_temp[7];
    Scatter2DPtr _h_results[6];
  };


  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1188891);


}
