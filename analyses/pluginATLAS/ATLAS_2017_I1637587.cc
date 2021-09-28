#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

#include "fastjet/contrib/SoftDrop.hh"

namespace Rivet {

  /// @brief Soft drop mass at 13 TeV
  class ATLAS_2017_I1637587: public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2017_I1637587);

    /// Book cuts and projections
    void init() {
      // All final state particles
      const FinalState fs(Cuts::abseta < 5.0);

      FastJets jets(fs, FastJets::ANTIKT, 0.8, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jets, "jets");

      book(_h_Table1, 1,1,1);
      book(_h_Table2, 2,1,1);
      book(_h_Table3, 3,1,1);

      book(_h_Table4, 4,1,1);
      book(_h_Table5, 5,1,1);
      book(_h_Table6, 6,1,1);

      betas = { 0., 1., 2. };
      ptBins = { 600, 650, 700, 750, 800, 850, 900, 950, 1000, 2000 };
      rhoBins = { -4.5, -4.1, -3.7, -3.3, -2.9, -2.5, -2.1, -1.7, -1.3, -0.9, -0.5 };
    }

    void analyze(const Event& event) {

      const Jets& myJets = apply<FastJets>(event, "jets").jetsByPt(400*GeV);

      if (myJets.size() < 2)  vetoEvent;
      if (myJets[0].pT() > 1.5*myJets[1].pT())  vetoEvent;
      if (myJets[0].abseta() > 1.5 || myJets[1].abseta() > 1.5)  vetoEvent;

      for (size_t i = 0; i < 2; ++i) {
      	if (myJets[i].pT() < 600*GeV) continue;
        ClusterSequence cs_ca(myJets[i].constituents(), JetDefinition(fastjet::cambridge_algorithm, 0.8));
        PseudoJets myJet_ca = sorted_by_pt(cs_ca.inclusive_jets(400.0));
        if(myJet_ca.size()==0) continue;
        for (size_t ibeta = 0; ibeta < 3; ++ibeta) {
          fastjet::contrib::SoftDrop sd(betas[ibeta], 0.1); //beta, zcut
          PseudoJet sdJet = sd(myJet_ca[0]);
	        double rho2 = pow(sdJet.m()/myJets[i].pT(),2);
      	  double log10rho2 = log(rho2)/log(10.);

	        if (log10rho2 < -4.5) continue;
      	  if (ibeta==0)  _h_Table1->fill(log10rho2);
      	  if (ibeta==1)  _h_Table2->fill(log10rho2);
      	  if (ibeta==2)  _h_Table3->fill(log10rho2);

       	  if (ibeta==0)  _h_Table4->fill(return_bin(rho2, myJets[i].pT()));
          if (ibeta==1)  _h_Table5->fill(return_bin(rho2, myJets[i].pT()));
          if (ibeta==2)  _h_Table6->fill(return_bin(rho2, myJets[i].pT()));
        }
      }
    }

    void finalize() {
      //Normalization comes here.
      double norm0 = 0.;
      double norm1 = 0.;
      double norm2 = 0.;
      for (size_t i = 4; i <= 7; ++i) { //only normalize in the resummation region.
        norm0+=_h_Table1->bin(i).height();
      	norm1+=_h_Table2->bin(i).height();
      	norm2+=_h_Table3->bin(i).height();
      }

      _h_Table1->scaleW(1.0/norm0);
      _h_Table2->scaleW(1.0/norm1);
      _h_Table3->scaleW(1.0/norm2);


      ptNorm( _h_Table4 );
      ptNorm( _h_Table5 );
      ptNorm( _h_Table6 );
    }

    void ptNorm(Histo1DPtr ptBinnedHist) {
      for (size_t k = 0; k < 9; ++k){
        double normalization = 0;
        for (size_t j = 4; j <= 7; ++j) {
          normalization += ptBinnedHist->bin(k*10 + j).height();
        }
        if( normalization == 0 ) continue;

        for (size_t j = 0; j < 10; ++j) {
          ptBinnedHist->bin(k*10 + j).scaleW(1. / normalization);
        }
      }

      return;
    }

    size_t return_bin(double rho, double pT){
      if (pT < 600.)  return -1;
      if (rho < pow(10,-4.5))  return -1;

      size_t pTbin = 1;
      if (pT < 600) pTbin = 0; //should not happen
      else if (pT < 650) pTbin = 1;
      else if (pT < 700) pTbin = 2;
      else if (pT < 750) pTbin = 3;
      else if (pT < 800) pTbin = 4;
      else if (pT < 850) pTbin = 5;
      else if (pT < 900) pTbin = 6;
      else if (pT < 950) pTbin = 7;
      else if (pT < 1000) pTbin = 8;
      else pTbin = 9;

      size_t rhobin = 1;
      if (rho < pow(10,-4.5)) rhobin = 0; //this should not happen.
      else if (rho < pow(10,-4.1)) rhobin = 1;
      else if (rho < pow(10,-3.7)) rhobin = 2;
      else if (rho < pow(10,-3.3)) rhobin = 3;
      else if (rho < pow(10,-2.9)) rhobin = 4;
      else if (rho < pow(10,-2.5)) rhobin = 5;
      else if (rho < pow(10,-2.1)) rhobin = 6;
      else if (rho < pow(10,-1.7)) rhobin = 7;
      else if (rho < pow(10,-1.3)) rhobin = 8;
      else if (rho < pow(10,-0.9)) rhobin = 9;
      else if (rho < pow(10,-0.5)) rhobin = 10;
      else rhobin = 10;

      return (rhobin-1) + (pTbin-1)*10;
    }


  private:
    /// Histograms
    Histo1DPtr _h_Table1, _h_Table2, _h_Table3, _h_Table4, _h_Table5, _h_Table6;
    vector<double> betas, ptBins, rhoBins;

  };

 DECLARE_RIVET_PLUGIN(ATLAS_2017_I1637587);
}

