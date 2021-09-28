#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/Selector.hh"


namespace Rivet {


  /// @brief Soft-drop mass at 13 TeV
  class ATLAS_2019_I1772062: public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2019_I1772062);


    void getQuarkGluon(Rivet::Histo1DPtr hForward, Rivet::Histo1DPtr hCentral, Rivet::Histo1DPtr hQuark, Rivet::Histo1DPtr hGluon, int ptbin, string parName, size_t beta) {

      int nBins = rhoBins.size() - 1;
      if (parName == "rg" || parName == "trg") nBins = rgBins.size()-1;
      if ((parName == "zg" || parName == "tzg") && beta == 0) nBins = zgBinsBeta0.size()-1;
      if ((parName == "zg" || parName == "tzg") && beta != 0) nBins = zgBins.size()-1;

      double FGC = gluonFractionCentral[ptbin];
      double FGF = gluonFractionForward[ptbin];
      double FQC = 1.-FGC;
      double FQF = 1.-FGF;

      for (size_t i=0; i<hGluon->numBins(); i++) {
        double binCenter = hGluon->bin(i).midpoint();
        double gVal = 0., qVal = 0.;
        if ((FQF -  FQC) != 0.) {
          gVal = (FQF * hCentral->bin(ptbin*(nBins) + i).height() - FQC * hForward->bin(ptbin*(nBins) + i).height()) / (FQF - FQC);
          qVal = (FGF * hCentral->bin(ptbin*(nBins) + i).height() - FGC * hForward->bin(ptbin*(nBins) + i).height()) / (FQC - FQF);
        }
        hGluon->fill(binCenter, gVal);
        hQuark->fill(binCenter, qVal);
      }

      histNorm(hQuark, parName);
      histNorm(hGluon, parName);
    }

    void ptNorm(Rivet::Histo1DPtr ptBinnedHist, std::string var, size_t beta) {
      size_t varNormBin1 = 0;
      size_t varNormBin2 = 0;
      size_t nBins = 10;

      if (var=="m" || var=="tm") {
        varNormBin1 = normBin1;
        varNormBin2 = normBin2;
      }
      if (var=="zg" || var=="tzg") {
        if(beta==0){
          varNormBin2 = zgBinsBeta0.size()-1;
          nBins = zgBinsBeta0.size()-1;
        }
        else {
          varNormBin2 = zgBins.size()-1;
          nBins = zgBins.size()-1;
        }
      }
      if (var=="rg" || var=="trg") {
        varNormBin2 = rgBins.size()-1;
        nBins = rgBins.size()-1;
      }

      for (size_t k=0; k< ptBins.size()-1; ++k) {
        double normalization = 0;

        for (size_t j=varNormBin1; j<varNormBin2; ++j) {
          double binWidth = 1.;
          if (var=="m" || var=="tm") {
            binWidth = rhoBins[j+1] - rhoBins[j];
          }
          if (var=="zg" || var=="tzg") {
            if(beta==0){
              binWidth = zgBinsBeta0[j+1]- zgBinsBeta0[j];
            }
            else {
              binWidth = zgBins[j+1]- zgBins[j];
            }
          }

          if (var=="rg" || var=="trg") {
            binWidth = rgBins[j+1]- rgBins[j];
            if (j==nBins-1) ptBinnedHist->bin( k*nBins+j ).scaleW(2.);
            //normalization += ptBinnedHist->bin(k*nBins+j).height()*binWidth;
            normalization += ptBinnedHist->bin(k*(nBins) + j).height()*binWidth;
          }
          else{
            normalization += ptBinnedHist->bin(k*(nBins) + j).height()*binWidth;
          }
        }

        if (normalization == 0) continue;

        for (unsigned int j=0; j<nBins; j++) {
          if (var=="rg" || var=="trg") {
            ptBinnedHist->bin(k*(nBins) + j ).scaleW(1. / (normalization) );
          }
          else{
            ptBinnedHist->bin(k*(nBins) + j ).scaleW(1. / (normalization));
          }

        }
      }

      return;
    }

    void histNorm(Rivet::Histo1DPtr hist, std::string var) {
      if (var=="m" || var=="tm") {
        double norm = 0.;
        for (size_t i = normBin1; i < normBin2; i++) { //only normalize in the resummation region.
          norm+=hist->bin(i).area();
        }
        if (norm > 0.) {
          hist->scaleW(1.0/(norm));
        }
      }
      else if ( var=="zg" || var=="tzg") {
        normalize(hist);
      }
      else{
        normalize(hist);
      }
    }


    int return_bin(float pT, float rho, std::string whichvar, size_t beta) {
      // First thing's first
      if (pT < ptBins[0]) return -100;

      if (whichvar=="m" && rho < pow(10,-4.5)) return -100;
      if (whichvar=="tm" && rho < pow(10,-4.5)) return -100;

      if (whichvar=="zg" && rho <= 0) return -100;
      if (whichvar=="tzg" && rho <= 0) return -100;

      if (whichvar=="rg" && rho <= -1.2) return -100;
      if (whichvar=="trg" && rho <= -1.2) return -100;

      if (whichvar=="id" && rho <= 1) return -100;

      int pTbin = 1;
      if (pT < ptBins[0]) pTbin = 0; //should not happen
      else if (pT < ptBins[1]) pTbin = 1;
      else if (pT < ptBins[2]) pTbin = 2;
      else if (pT < ptBins[3]) pTbin = 3;
      else if (pT < ptBins[4]) pTbin = 4;
      else pTbin = 5;
      if (pTbin == 0) return -1;

      int rhobin = 1.;
      if ((whichvar=="m") || (whichvar=="tm"))
        {
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
          return rhobin*1. + (pTbin*1.-1.)*10.-1;
        }

      // zg
      else if (((whichvar=="zg")||(whichvar=="tzg")) && beta == 0)
        {
          if (rho < 0.10) return -10;
          else if (rho < 0.15) rhobin = 1;
          else if (rho < 0.20) rhobin = 2;
          else if (rho < 0.25) rhobin = 3;
          else if (rho < 0.30) rhobin = 4;
          else if (rho < 0.35) rhobin = 5;
          else if (rho < 0.40) rhobin = 6;
          else if (rho < 0.45) rhobin = 7;
          else if (rho < 0.50) rhobin = 8;
          else rhobin = 8;
          return rhobin*1. + (pTbin*1.-1.)*8.-1;
        }
      else if (((whichvar=="zg")||(whichvar=="tzg")) && beta != 0)
        {
          if (rho < 0.00) return -10;
          else if (rho < 0.05) rhobin = 1;
          else if (rho < 0.10) rhobin = 2;
          else if (rho < 0.15) rhobin = 3;
          else if (rho < 0.20) rhobin = 4;
          else if (rho < 0.25) rhobin = 5;
          else if (rho < 0.30) rhobin = 6;
          else if (rho < 0.35) rhobin = 7;
          else if (rho < 0.40) rhobin = 8;
          else if (rho < 0.45) rhobin = 9;
          else if (rho < 0.50) rhobin = 10;
          else rhobin = 10;
          return rhobin*1. + (pTbin*1.-1.)*10.-1;
        }

      //rg
      else if ((whichvar=="rg")||(whichvar=="trg"))
        {
          if (rho < -1.2) return -10;
          else if (rho < -1.0) rhobin = 1;
          else if (rho < -0.8) rhobin = 2;
          else if (rho < -0.6) rhobin = 3;
          else if (rho < -0.4) rhobin = 4;
          else if (rho < -0.2) rhobin = 5;
          else if (rho < -0.1) rhobin = 6;
          else rhobin = 6;
          return rhobin*1. + (pTbin*1.-1.)*6.-1;
        }

      return -100;
    }



    /// Book cuts and projections
    void init() {
      // All final state particles
      const FinalState fs(Cuts::abseta < 5.0);

      FastJets jets(fs, FastJets::ANTIKT, 0.8, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jets, "jets");

      ChargedFinalState tracks(Cuts::pT > 0.5*GeV && Cuts::abseta < 2.5);
      declare(tracks, "tracks");

      normBin1 = 2; normBin2 = 7;
      gluonFractionCentral = {0.75, 0.72, 0.66, 0.61, 0.54};
      gluonFractionForward = {0.70, 0.64, 0.57, 0.51, 0.43};
      rgBins = {-1.2, -1.0, -0.8, -0.6, -0.4, -0.2, -0.1};
      zgBins = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5};
      zgBinsBeta0 = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5};
      rhoBins = {-4.5, -4.1, -3.7, -3.3, -2.9, -2.5, -2.1, -1.7, -1.3, -0.9, -0.5};
      ptBins = {300, 400, 600, 800, 1000, 2000};



      book(_h["_h_Table1"], 1,1,1);   // Cluster, Rho, beta=0, inclusive pT, both jets
      book(_h["_h_Table2"], 2,1,1);   // Track, Rho, beta=0,   inclusive pT, both jets
      book(_h["_h_Table3"], 3,1,1);   // Cluster, Rho, beta=1, inclusive pT, both jets
      book(_h["_h_Table4"], 4,1,1);   // Track, Rho, beta=1,   inclusive pT, both jets
      book(_h["_h_Table5"], 5,1,1);   // Cluster, Rho, beta=2, inclusive pT, both jets
      book(_h["_h_Table6"], 6,1,1);   // Track, Rho, beta=2,   inclusive pT, both jets

      book(_h["_h_Table7"], 7,1,1);   // Cluster, zg, beta=0,  inclusive pT, both jets
      book(_h["_h_Table8"], 8,1,1);   // Track, zg, beta=0,    inclusive pT, both jets
      book(_h["_h_Table9"], 9,1,1);   // Cluster, zg, beta=1,  inclusive pT, both jets
      book(_h["_h_Table10"], 10,1,1);  // Track, zg, beta=1,    inclusive pT, both jets
      book(_h["_h_Table11"], 11,1,1);  // Cluster, zg, beta=2,  inclusive pT, both jets
      book(_h["_h_Table12"], 12,1,1);  // Track, zg, beta=2,    inclusive pT, both jets

      book(_h["_h_Table13"], 13,1,1);  // Cluster, rg, beta=0,  inclusive pT, both jets
      book(_h["_h_Table14"], 14,1,1);  // Track, rg, beta=0,    inclusive pT, both jets
      book(_h["_h_Table15"], 15,1,1);  // Cluster, rg, beta=1,  inclusive pT, both jets
      book(_h["_h_Table16"], 16,1,1);  // Track, rg, beta=1,    inclusive pT, both jets
      book(_h["_h_Table17"], 17,1,1);  // Cluster, rg, beta=2,  inclusive pT, both jets
      book(_h["_h_Table18"], 18,1,1);  // Track, rg, beta=2,    inclusive pT, both jets


      book(_h["_h_Table19"], 19,1,1);  // Cluster, Rho, beta=0, inclusive pT, central jet
      book(_h["_h_Table20"], 20,1,1);  // Track, Rho, beta=0,   inclusive pT, central jet
      book(_h["_h_Table21"], 21,1,1);  // Cluster, Rho, beta=1, inclusive pT, central jet
      book(_h["_h_Table22"], 22,1,1);  // Track, Rho, beta=1,   inclusive pT, central jet
      book(_h["_h_Table23"], 23,1,1);  // Cluster, Rho, beta=2, inclusive pT, central jet
      book(_h["_h_Table24"], 24,1,1);  // Track, Rho, beta=2,   inclusive pT, central jet

      book(_h["_h_Table25"], 25,1,1);  // Cluster, zg, beta=0,  inclusive pT, central jet
      book(_h["_h_Table26"], 26,1,1);  // Track, zg, beta=0,    inclusive pT, central jet
      book(_h["_h_Table27"], 27,1,1);  // Cluster, zg, beta=1,  inclusive pT, central jet
      book(_h["_h_Table28"], 28,1,1);  // Track, zg, beta=1,    inclusive pT, central jet
      book(_h["_h_Table29"], 29,1,1);  // Cluster, zg, beta=2,  inclusive pT, central jet
      book(_h["_h_Table30"], 30,1,1);  // Track, zg, beta=2,    inclusive pT, central jet

      book(_h["_h_Table31"], 31,1,1);  // Cluster, rg, beta=0,  inclusive pT, central jet
      book(_h["_h_Table32"], 32,1,1);  // Track, rg, beta=0,    inclusive pT, central jet
      book(_h["_h_Table33"], 33,1,1);  // Cluster, rg, beta=1,  inclusive pT, central jet
      book(_h["_h_Table34"], 34,1,1);  // Track, rg, beta=1,    inclusive pT, central jet
      book(_h["_h_Table35"], 35,1,1);  // Cluster, rg, beta=2,  inclusive pT, central jet
      book(_h["_h_Table36"], 36,1,1);  // Track, rg, beta=2,    inclusive pT, central jet


      book(_h["_h_Table37"], 37,1,1);  // Cluster, Rho, beta=0, inclusive pT, forward jet
      book(_h["_h_Table38"], 38,1,1);  // Track, Rho, beta=0,   inclusive pT, forward jet
      book(_h["_h_Table39"], 39,1,1);  // Cluster, Rho, beta=1, inclusive pT, forward jet
      book(_h["_h_Table40"], 40,1,1);  // Track, Rho, beta=1,   inclusive pT, forward jet
      book(_h["_h_Table41"], 41,1,1);  // Cluster, Rho, beta=2, inclusive pT, forward jet
      book(_h["_h_Table42"], 42,1,1);  // Track, Rho, beta=2,   inclusive pT, forward jet

      book(_h["_h_Table43"], 43,1,1);  // Cluster, zg, beta=0,  inclusive pT, forward jet
      book(_h["_h_Table44"], 44,1,1);  // Track, zg, beta=0,    inclusive pT, forward jet
      book(_h["_h_Table45"], 45,1,1);  // Cluster, zg, beta=1,  inclusive pT, forward jet
      book(_h["_h_Table46"], 46,1,1);  // Track, zg, beta=1,    inclusive pT, forward jet
      book(_h["_h_Table47"], 47,1,1);  // Cluster, zg, beta=2,  inclusive pT, forward jet
      book(_h["_h_Table48"], 48,1,1);  // Track, zg, beta=2,    inclusive pT, forward jet

      book(_h["_h_Table49"], 49,1,1);  // Cluster, rg, beta=0,  inclusive pT, forward jet
      book(_h["_h_Table50"], 50,1,1);  // Track, rg, beta=0,    inclusive pT, forward jet
      book(_h["_h_Table51"], 51,1,1);  // Cluster, rg, beta=1,  inclusive pT, forward jet
      book(_h["_h_Table52"], 52,1,1);  // Track, rg, beta=1,    inclusive pT, forward jet
      book(_h["_h_Table53"], 53,1,1);  // Cluster, rg, beta=2,  inclusive pT, forward jet
      book(_h["_h_Table54"], 54,1,1);  // Track, rg, beta=2,    inclusive pT, forward jet




      book(_h["_h_Table55"], 55,1,1);  // Cluster, Rho, beta=0, pT binned, both jets
      book(_h["_h_Table56"], 56,1,1);  // Track, Rho, beta=0,   pT binned, both jets
      book(_h["_h_Table57"], 57,1,1);  // Cluster, Rho, beta=1, pT binned, both jets
      book(_h["_h_Table58"], 58,1,1);  // Track, Rho, beta=1,   pT binned, both jets
      book(_h["_h_Table59"], 59,1,1);  // Cluster, Rho, beta=2, pT binned, both jets
      book(_h["_h_Table60"], 60,1,1);  // Track, Rho, beta=2,   pT binned, both jets

      book(_h["_h_Table61"], 61,1,1);  // Cluster, zg, beta=0,  pT binned, both jets
      book(_h["_h_Table62"], 62,1,1);  // Track, zg, beta=0,    pT binned, both jets
      book(_h["_h_Table63"], 63,1,1);  // Cluster, zg, beta=1,  pT binned, both jets
      book(_h["_h_Table64"], 64,1,1);  // Track, zg, beta=1,    pT binned, both jets
      book(_h["_h_Table65"], 65,1,1);  // Cluster, zg, beta=2,  pT binned, both jets
      book(_h["_h_Table66"], 66,1,1);  // Track, zg, beta=2,    pT binned, both jets

      book(_h["_h_Table67"], 67,1,1);  // Cluster, rg, beta=0,  pT binned, both jets
      book(_h["_h_Table68"], 68,1,1);  // Track, rg, beta=0,    pT binned, both jets
      book(_h["_h_Table69"], 69,1,1);  // Cluster, rg, beta=1,  pT binned, both jets
      book(_h["_h_Table70"], 70,1,1);  // Track, rg, beta=1,    pT binned, both jets
      book(_h["_h_Table71"], 71,1,1);  // Cluster, rg, beta=2,  pT binned, both jets
      book(_h["_h_Table72"], 72,1,1);  // Track, rg, beta=2,    pT binned, both jets


      book(_h["_h_Table73"], 73,1,1);  // Cluster, Rho, beta=0, pT binned, central jet
      book(_h["_h_Table74"], 74,1,1);  // Track, Rho, beta=0,   pT binned, central jet
      book(_h["_h_Table75"], 75,1,1);  // Cluster, Rho, beta=1, pT binned, central jet
      book(_h["_h_Table76"], 76,1,1);  // Track, Rho, beta=1,   pT binned, central jet
      book(_h["_h_Table77"], 77,1,1);  // Cluster, Rho, beta=2, pT binned, central jet
      book(_h["_h_Table78"], 78,1,1);  // Track, Rho, beta=2,   pT binned, central jet

      book(_h["_h_Table79"], 79,1,1);  // Cluster, zg, beta=0,  pT binned, central jet
      book(_h["_h_Table80"], 80,1,1);  // Track, zg, beta=0,    pT binned, central jet
      book(_h["_h_Table81"], 81,1,1);  // Cluster, zg, beta=1,  pT binned, central jet
      book(_h["_h_Table82"], 82,1,1); // Track, zg, beta=1,    pT binned, central jet
      book(_h["_h_Table83"], 83,1,1); // Cluster, zg, beta=2,  pT binned, central jet
      book(_h["_h_Table84"], 84,1,1); // Track, zg, beta=2,    pT binned, central jet

      book(_h["_h_Table85"], 85,1,1); // Cluster, rg, beta=0,  pT binned, central jet
      book(_h["_h_Table86"], 86,1,1); // Track, rg, beta=0,    pT binned, central jet
      book(_h["_h_Table87"], 87,1,1); // Cluster, rg, beta=1,  pT binned, central jet
      book(_h["_h_Table88"], 88,1,1); // Track, rg, beta=1,    pT binned, central jet
      book(_h["_h_Table89"], 89,1,1); // Cluster, rg, beta=2,  pT binned, central jet
      book(_h["_h_Table90"], 90,1,1); // Track, rg, beta=2,    pT binned, central jet


      book(_h["_h_Table91"], 91,1,1); // Cluster, Rho, beta=0, pT binned, forward jet
      book(_h["_h_Table92"], 92,1,1); // Track, Rho, beta=0,   pT binned, forward jet
      book(_h["_h_Table93"], 93,1,1); // Cluster, Rho, beta=1, pT binned, forward jet
      book(_h["_h_Table94"], 94,1,1); // Track, Rho, beta=1,   pT binned, forward jet
      book(_h["_h_Table95"], 95,1,1); // Cluster, Rho, beta=2, pT binned, forward jet
      book(_h["_h_Table96"], 96,1,1); // Track, Rho, beta=2,   pT binned, forward jet

      book(_h["_h_Table97"], 97,1,1); // Cluster, zg, beta=0,  pT binned, forward jet
      book(_h["_h_Table98"], 98,1,1); // Track, zg, beta=0,    pT binned, forward jet
      book(_h["_h_Table99"], 99,1,1); // Cluster, zg, beta=1,  pT binned, forward jet
      book(_h["_h_Table100"], 100,1,1); // Track, zg, beta=1,    pT binned, forward jet
      book(_h["_h_Table101"], 101,1,1); // Cluster, zg, beta=2,  pT binned, forward jet
      book(_h["_h_Table102"], 102,1,1); // Track, zg, beta=2,    pT binned, forward jet

      book(_h["_h_Table103"], 103,1,1); // Cluster, rg, beta=0,  pT binned, forward jet
      book(_h["_h_Table104"], 104,1,1); // Track, rg, beta=0,    pT binned, forward jet
      book(_h["_h_Table105"], 105,1,1); // Cluster, rg, beta=1,  pT binned, forward jet
      book(_h["_h_Table106"], 106,1,1); // Track, rg, beta=1,    pT binned, forward jet
      book(_h["_h_Table107"], 107,1,1); // Cluster, rg, beta=2,  pT binned, forward jet
      book(_h["_h_Table108"], 108,1,1); // Track, rg, beta=2,    pT binned, forward jet




      book(_h["_h_Table109"], 109,1,1);   // Cluster, Rho, beta=0, quark jet
      book(_h["_h_Table110"], 110,1,1);   // Track, Rho, beta=0,   quark jet
      book(_h["_h_Table111"], 111,1,1);   // Cluster, Rho, beta=1, quark jet
      book(_h["_h_Table112"], 112,1,1);   // Track, Rho, beta=1,   quark jet
      book(_h["_h_Table113"], 113,1,1);   // Cluster, Rho, beta=2, quark jet
      book(_h["_h_Table114"], 114,1,1);   // Track, Rho, beta=2,   quark jet

      book(_h["_h_Table115"], 115,1,1);   // Cluster, zg, beta=0,  quark jet
      book(_h["_h_Table116"], 116,1,1);   // Track, zg, beta=0,    quark jet
      book(_h["_h_Table117"], 117,1,1);   // Cluster, zg, beta=1,  quark jet
      book(_h["_h_Table118"], 118,1,1);  // Track, zg, beta=1,    quark jet
      book(_h["_h_Table119"], 119,1,1);  // Cluster, zg, beta=2,  quark jet
      book(_h["_h_Table120"], 120,1,1);  // Track, zg, beta=2,    quark jet

      book(_h["_h_Table121"], 121,1,1);  // Cluster, rg, beta=0,  quark jet
      book(_h["_h_Table122"], 122,1,1);  // Track, rg, beta=0,    quark jet
      book(_h["_h_Table123"], 123,1,1);  // Cluster, rg, beta=1,  quark jet
      book(_h["_h_Table124"], 124,1,1);  // Track, rg, beta=1,    quark jet
      book(_h["_h_Table125"], 125,1,1);  // Cluster, rg, beta=2,  quark jet
      book(_h["_h_Table126"], 126,1,1);  // Track, rg, beta=2,    quark jet


      book(_h["_h_Table127"], 127,1,1);   // Cluster, Rho, beta=0, gluon jet
      book(_h["_h_Table128"], 128,1,1);   // Track, Rho, beta=0,   gluon jet
      book(_h["_h_Table129"], 129,1,1);   // Cluster, Rho, beta=1, gluon jet
      book(_h["_h_Table130"], 130,1,1);   // Track, Rho, beta=1,   gluon jet
      book(_h["_h_Table131"], 131,1,1);   // Cluster, Rho, beta=2, gluon jet
      book(_h["_h_Table132"], 132,1,1);   // Track, Rho, beta=2,   gluon jet

      book(_h["_h_Table133"], 133,1,1);   // Cluster, zg, beta=0,  gluon jet
      book(_h["_h_Table134"], 134,1,1);   // Track, zg, beta=0,    gluon jet
      book(_h["_h_Table135"], 135,1,1);   // Cluster, zg, beta=1,  gluon jet
      book(_h["_h_Table136"], 136,1,1);  // Track, zg, beta=1,    gluon jet
      book(_h["_h_Table137"], 137,1,1);  // Cluster, zg, beta=2,  gluon jet
      book(_h["_h_Table138"], 138,1,1);  // Track, zg, beta=2,    gluon jet

      book(_h["_h_Table139"], 139,1,1);  // Cluster, rg, beta=0,  gluon jet
      book(_h["_h_Table140"], 140,1,1);  // Track, rg, beta=0,    gluon jet
      book(_h["_h_Table141"], 141,1,1);  // Cluster, rg, beta=1,  gluon jet
      book(_h["_h_Table142"], 142,1,1);  // Track, rg, beta=1,    gluon jet
      book(_h["_h_Table143"], 143,1,1);  // Cluster, rg, beta=2,  gluon jet
      book(_h["_h_Table144"], 144,1,1);  // Track, rg, beta=2,    gluon jet

    }


    void analyze(const Event& event) {

      const Jets& myJets = apply<FastJets>(event, "jets").jetsByPt(200*GeV);
      const Particles& tracks = apply<ChargedFinalState>(event, "tracks").particlesByPt();

      if (myJets.size() < 2)  vetoEvent;
      if (myJets[0].pT() > 1.5*myJets[1].pT())  vetoEvent;
      if (myJets[0].abseta() > 1.5 || myJets[1].abseta() > 1.5) vetoEvent;

      std::vector<bool> isCentral;
      isCentral.push_back(true);
      isCentral.push_back(false);
      if (myJets[0].abseta() > myJets[1].abseta()) {
        isCentral[0] = false;
        isCentral[1] = true;
      }

      for (size_t i = 0; i < 2; ++i) {
        if (myJets[i].pT() < 300*GeV) continue;

        vector<fastjet::PseudoJet> charged_constituents;
        for (const Particle& p : tracks) {
          const double dr = deltaR(myJets[i], p, PSEUDORAPIDITY);
          if (dr > 0.8) continue;
          if (abs(p.pid()) == 13) continue;
          charged_constituents.push_back(p);
        }

        fastjet::ClusterSequence cs_ca(myJets[i].constituents(), fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.8));
        vector<fastjet::PseudoJet> myJet_ca = fastjet::sorted_by_pt(cs_ca.inclusive_jets(10.0));

        fastjet::ClusterSequence cs_ca_charged(charged_constituents, fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.8));
        vector<fastjet::PseudoJet> myJet_ca_charged = fastjet::sorted_by_pt(cs_ca_charged.inclusive_jets(10.0));

        if (myJet_ca.size()==0) continue;
        if (myJet_ca_charged.size()==0) continue;

        // grooming parameters that are scanned.
        vector<size_t> betas = { 0, 1, 2 };
        for (size_t ibeta : betas) {

          fastjet::contrib::SoftDrop sd(double(ibeta), 0.1); //beta, zcut
          fastjet::PseudoJet sdJet = sd(myJet_ca[0]);
          fastjet::PseudoJet sdJet_charged = sd(myJet_ca_charged[0]);

          double rho2              = pow(sdJet.m()/myJets[i].pT(),2);
          double log10rho2         = log(rho2)/log(10.);
          double zg                = sdJet.structure_of<fastjet::contrib::SoftDrop>().symmetry();
          double rg                = sdJet.structure_of<fastjet::contrib::SoftDrop>().delta_R();
          rg = (rg > 0) ? log(rg) / log(10.) : -100;

          double rho2_charged      = pow(sdJet_charged.m()/myJet_ca_charged[0].pt(),2);
          double log10rho2_charged = log(rho2_charged)/log(10.);
          double zg_charged        = sdJet_charged.structure_of<fastjet::contrib::SoftDrop>().symmetry();
          double rg_charged             = sdJet.structure_of<fastjet::contrib::SoftDrop>().delta_R();
          rg_charged = (rg_charged > 0) ? log(rg_charged) / log(10.) : -100;

          double pt_log10rho2 = return_bin(myJets[i].pT()/GeV, rho2, "m", ibeta);
          double pt_zg = return_bin(myJets[i].pT()/GeV, zg, "zg", ibeta);
          double pt_rg = return_bin(myJets[i].pT()/GeV, rg, "rg", ibeta);

          double pt_log10rho2_charged = return_bin(myJets[i].pT()/GeV, rho2_charged, "tm", ibeta);
          double pt_zg_charged = return_bin(myJets[i].pT()/GeV,        zg_charged,   "tzg", ibeta);
          double pt_rg_charged = return_bin(myJets[i].pT()/GeV,        rg_charged,   "trg", ibeta);


          if (ibeta==0)  {
            _h["_h_Table1"]->fill(  log10rho2);
            _h["_h_Table2"]->fill(  log10rho2_charged);
            _h["_h_Table7"]->fill(  zg);
            _h["_h_Table8"]->fill(  zg_charged);
            if (rg > -1.2)
              _h["_h_Table13"]->fill(  rg);
            if (rg_charged > -1.2)
              _h["_h_Table14"]->fill(  rg_charged);

            _h["_h_Table55"]->fill(  pt_log10rho2);
            _h["_h_Table56"]->fill(  pt_log10rho2_charged);
            _h["_h_Table61"]->fill(  pt_zg);
            _h["_h_Table62"]->fill(  pt_zg_charged);
            _h["_h_Table67"]->fill(  pt_rg);
            _h["_h_Table68"]->fill(  pt_rg_charged);

            if (isCentral[i]) {
              _h["_h_Table19"]->fill(  log10rho2);
              _h["_h_Table20"]->fill(  log10rho2_charged);
              _h["_h_Table25"]->fill(  zg);
              _h["_h_Table26"]->fill(  zg_charged);
              if (rg > -1.2)
                _h["_h_Table31"]->fill(  rg);
              if (rg_charged > -1.2)
                _h["_h_Table32"]->fill(  rg_charged);

              _h["_h_Table73"]->fill(  pt_log10rho2);
              _h["_h_Table74"]->fill(  pt_log10rho2_charged);
              _h["_h_Table79"]->fill(  pt_zg);
              _h["_h_Table80"]->fill(  pt_zg_charged);
              _h["_h_Table85"]->fill(  pt_rg);
              _h["_h_Table86"]->fill(  pt_rg_charged);
            }

            if (!isCentral[i]) {
              _h["_h_Table37"]->fill(  log10rho2);
              _h["_h_Table38"]->fill(  log10rho2_charged);
              _h["_h_Table43"]->fill(  zg);
              _h["_h_Table44"]->fill(  zg_charged);
              if (rg > -1.2)
                _h["_h_Table49"]->fill(  rg);
              if (rg_charged > -1.2)
                _h["_h_Table50"]->fill(  rg_charged);

              _h["_h_Table91"]->fill(  pt_log10rho2);
              _h["_h_Table92"]->fill(  pt_log10rho2_charged);
              _h["_h_Table97"]->fill(  pt_zg);
              _h["_h_Table98"]->fill(  pt_zg_charged);
              _h["_h_Table103"]->fill(  pt_rg);
              _h["_h_Table104"]->fill(  pt_rg_charged);
            }
          }
          if (ibeta==1)  {
            _h["_h_Table3"]->fill(  log10rho2);
            _h["_h_Table4"]->fill(  log10rho2_charged);
            _h["_h_Table9"]->fill(  zg);
            _h["_h_Table10"]->fill(  zg_charged);
            if (rg > -1.2)
              _h["_h_Table15"]->fill(  rg);
            if (rg_charged > -1.2)
              _h["_h_Table16"]->fill(  rg_charged);

            _h["_h_Table57"]->fill(  pt_log10rho2);
            _h["_h_Table58"]->fill(  pt_log10rho2_charged);
            _h["_h_Table63"]->fill(  pt_zg);
            _h["_h_Table64"]->fill(  pt_zg_charged);
            _h["_h_Table69"]->fill(  pt_rg);
            _h["_h_Table70"]->fill(  pt_rg_charged);

            if (isCentral[i]) {
              _h["_h_Table21"]->fill(  log10rho2);
              _h["_h_Table22"]->fill(  log10rho2_charged);
              _h["_h_Table27"]->fill(  zg);
              _h["_h_Table28"]->fill(  zg_charged);
              if (rg > -1.2)
                _h["_h_Table33"]->fill(  rg);
              if (rg_charged > -1.2)
                _h["_h_Table34"]->fill(  rg_charged);

              _h["_h_Table75"]->fill(  pt_log10rho2);
              _h["_h_Table76"]->fill(  pt_log10rho2_charged);
              _h["_h_Table81"]->fill(  pt_zg);
              _h["_h_Table82"]->fill(  pt_zg_charged);
              _h["_h_Table87"]->fill(  pt_rg);
              _h["_h_Table88"]->fill(  pt_rg_charged);
            }

            if (!isCentral[i]) {
              _h["_h_Table39"]->fill(  log10rho2);
              _h["_h_Table40"]->fill(  log10rho2_charged);
              _h["_h_Table45"]->fill(  zg);
              _h["_h_Table46"]->fill(  zg_charged);
              if (rg > -1.2)
                _h["_h_Table51"]->fill(  rg);
              if (rg_charged > -1.2)
                _h["_h_Table52"]->fill(  rg_charged);

              _h["_h_Table93"]->fill(  pt_log10rho2);
              _h["_h_Table94"]->fill(  pt_log10rho2_charged);
              _h["_h_Table99"]->fill(  pt_zg);
              _h["_h_Table100"]->fill(  pt_zg_charged);
              _h["_h_Table105"]->fill(  pt_rg);
              _h["_h_Table106"]->fill(  pt_rg_charged);
            }
          }
          if (ibeta==2)  {
            _h["_h_Table5"]->fill(  log10rho2);
            _h["_h_Table6"]->fill(  log10rho2_charged);
            _h["_h_Table11"]->fill(  zg);
            _h["_h_Table12"]->fill(  zg_charged);
            if (rg > -1.2)
              _h["_h_Table17"]->fill(  rg);
            if (rg_charged > -1.2)
              _h["_h_Table18"]->fill(  rg_charged);

            _h["_h_Table59"]->fill(  pt_log10rho2);
            _h["_h_Table60"]->fill(  pt_log10rho2_charged);
            _h["_h_Table65"]->fill(  pt_zg);
            _h["_h_Table66"]->fill(  pt_zg_charged);
            _h["_h_Table71"]->fill(  pt_rg);
            _h["_h_Table72"]->fill(  pt_rg_charged);

            if (isCentral[i]) {
              _h["_h_Table23"]->fill(  log10rho2);
              _h["_h_Table24"]->fill(  log10rho2_charged);
              _h["_h_Table29"]->fill(  zg);
              _h["_h_Table30"]->fill(  zg_charged);
              if (rg > -1.2)
                _h["_h_Table35"]->fill(  rg);
              if (rg_charged > -1.2)
                _h["_h_Table36"]->fill(  rg_charged);

              _h["_h_Table77"]->fill(  pt_log10rho2);
              _h["_h_Table78"]->fill(  pt_log10rho2_charged);
              _h["_h_Table83"]->fill(  pt_zg);
              _h["_h_Table84"]->fill(  pt_zg_charged);
              _h["_h_Table89"]->fill(  pt_rg);
              _h["_h_Table90"]->fill(  pt_rg_charged);
            }

            if (!isCentral[i]) {
              _h["_h_Table41"]->fill(  log10rho2);
              _h["_h_Table42"]->fill(  log10rho2_charged);
              _h["_h_Table47"]->fill(  zg);
              _h["_h_Table48"]->fill(  zg_charged);
              if (rg > -1.2)
                _h["_h_Table53"]->fill(  rg);
              if (rg_charged > -1.2)
                _h["_h_Table54"]->fill(  rg_charged);

              _h["_h_Table95"]->fill(  pt_log10rho2);
              _h["_h_Table96"]->fill(  pt_log10rho2_charged);
              _h["_h_Table101"]->fill(  pt_zg);
              _h["_h_Table102"]->fill(  pt_zg_charged);
              _h["_h_Table107"]->fill(  pt_rg);
              _h["_h_Table108"]->fill(  pt_rg_charged);
            }
          }

        }
      }
    }

    void finalize() {
      //Normalization comes here.
      histNorm(_h["_h_Table1"], "m");
      histNorm(_h["_h_Table2"], "m");
      histNorm(_h["_h_Table3"], "m");
      histNorm(_h["_h_Table4"], "m");
      histNorm(_h["_h_Table5"], "m");
      histNorm(_h["_h_Table6"], "m");
      histNorm(_h["_h_Table7"], "zg");
      histNorm(_h["_h_Table8"], "zg");
      histNorm(_h["_h_Table9"], "zg");
      histNorm(_h["_h_Table10"], "zg");
      histNorm(_h["_h_Table11"], "zg");
      histNorm(_h["_h_Table12"], "zg");
      histNorm(_h["_h_Table13"], "rg");
      histNorm(_h["_h_Table14"], "rg");
      histNorm(_h["_h_Table15"], "rg");
      histNorm(_h["_h_Table16"], "rg");
      histNorm(_h["_h_Table17"], "rg");
      histNorm(_h["_h_Table18"], "rg");

      histNorm(_h["_h_Table19"], "m");
      histNorm(_h["_h_Table20"], "m");
      histNorm(_h["_h_Table21"], "m");
      histNorm(_h["_h_Table22"], "m");
      histNorm(_h["_h_Table23"], "m");
      histNorm(_h["_h_Table24"], "m");
      histNorm(_h["_h_Table25"], "zg");
      histNorm(_h["_h_Table26"], "zg");
      histNorm(_h["_h_Table27"], "zg");
      histNorm(_h["_h_Table28"], "zg");
      histNorm(_h["_h_Table29"], "zg");
      histNorm(_h["_h_Table30"], "zg");
      histNorm(_h["_h_Table31"], "rg");
      histNorm(_h["_h_Table32"], "rg");
      histNorm(_h["_h_Table33"], "rg");
      histNorm(_h["_h_Table34"], "rg");
      histNorm(_h["_h_Table35"], "rg");
      histNorm(_h["_h_Table36"], "rg");


      histNorm(_h["_h_Table37"], "m");
      histNorm(_h["_h_Table38"], "m");
      histNorm(_h["_h_Table39"], "m");
      histNorm(_h["_h_Table40"], "m");
      histNorm(_h["_h_Table41"], "m");
      histNorm(_h["_h_Table42"], "m");
      histNorm(_h["_h_Table43"], "zg");
      histNorm(_h["_h_Table44"], "zg");
      histNorm(_h["_h_Table45"], "zg");
      histNorm(_h["_h_Table46"], "zg");
      histNorm(_h["_h_Table47"], "zg");
      histNorm(_h["_h_Table48"], "zg");
      histNorm(_h["_h_Table49"], "rg");
      histNorm(_h["_h_Table50"], "rg");
      histNorm(_h["_h_Table51"], "rg");
      histNorm(_h["_h_Table52"], "rg");
      histNorm(_h["_h_Table53"], "rg");
      histNorm(_h["_h_Table54"], "rg");


      ptNorm(_h["_h_Table55"], "m", 0);
      ptNorm(_h["_h_Table56"], "m", 0);
      ptNorm(_h["_h_Table57"], "m", 1);
      ptNorm(_h["_h_Table58"], "m", 1);
      ptNorm(_h["_h_Table59"], "m", 2);
      ptNorm(_h["_h_Table60"], "m", 2);
      ptNorm(_h["_h_Table61"], "zg", 0);
      ptNorm(_h["_h_Table62"], "zg", 0);
      ptNorm(_h["_h_Table63"], "zg", 1);
      ptNorm(_h["_h_Table64"], "zg", 1);
      ptNorm(_h["_h_Table65"], "zg", 2);
      ptNorm(_h["_h_Table66"], "zg", 2);
      ptNorm(_h["_h_Table67"], "rg", 0);
      ptNorm(_h["_h_Table68"], "rg", 0);
      ptNorm(_h["_h_Table69"], "rg", 1);
      ptNorm(_h["_h_Table70"], "rg", 1);
      ptNorm(_h["_h_Table71"], "rg", 2);
      ptNorm(_h["_h_Table72"], "rg", 2);

      ptNorm(_h["_h_Table73"], "m", 0);
      ptNorm(_h["_h_Table74"], "m", 0);
      ptNorm(_h["_h_Table75"], "m", 1);
      ptNorm(_h["_h_Table76"], "m", 1);
      ptNorm(_h["_h_Table77"], "m", 2);
      ptNorm(_h["_h_Table78"], "m", 2);
      ptNorm(_h["_h_Table79"], "zg", 0);
      ptNorm(_h["_h_Table80"], "zg", 0);
      ptNorm(_h["_h_Table81"], "zg", 1);
      ptNorm(_h["_h_Table82"], "zg", 1);
      ptNorm(_h["_h_Table83"], "zg", 2);
      ptNorm(_h["_h_Table84"], "zg", 2);
      ptNorm(_h["_h_Table85"], "rg", 0);
      ptNorm(_h["_h_Table86"], "rg", 0);
      ptNorm(_h["_h_Table87"], "rg", 1);
      ptNorm(_h["_h_Table88"], "rg", 1);
      ptNorm(_h["_h_Table89"], "rg", 2);
      ptNorm(_h["_h_Table90"], "rg", 2);


      ptNorm(_h["_h_Table91"], "m", 0);
      ptNorm(_h["_h_Table92"], "m", 0);
      ptNorm(_h["_h_Table93"], "m", 1);
      ptNorm(_h["_h_Table94"], "m", 1);
      ptNorm(_h["_h_Table95"], "m", 2);
      ptNorm(_h["_h_Table96"], "m", 2);
      ptNorm(_h["_h_Table97"], "zg", 0);
      ptNorm(_h["_h_Table98"], "zg", 0);
      ptNorm(_h["_h_Table99"], "zg", 1);
      ptNorm(_h["_h_Table100"], "zg", 1);
      ptNorm(_h["_h_Table101"], "zg", 2);
      ptNorm(_h["_h_Table102"], "zg", 2);
      ptNorm(_h["_h_Table103"], "rg", 0);
      ptNorm(_h["_h_Table104"], "rg", 0);
      ptNorm(_h["_h_Table105"], "rg", 1);
      ptNorm(_h["_h_Table106"], "rg", 1);
      ptNorm(_h["_h_Table107"], "rg", 2);
      ptNorm(_h["_h_Table108"], "rg", 2);



      getQuarkGluon(_h["_h_Table91"], _h["_h_Table73"], _h["_h_Table109"], _h["_h_Table127"], 2, "m", 0);
      getQuarkGluon(_h["_h_Table92"], _h["_h_Table74"], _h["_h_Table110"], _h["_h_Table128"], 2, "m", 0);
      getQuarkGluon(_h["_h_Table93"], _h["_h_Table75"], _h["_h_Table111"], _h["_h_Table129"], 2, "m", 1);
      getQuarkGluon(_h["_h_Table94"], _h["_h_Table76"], _h["_h_Table112"], _h["_h_Table130"], 2, "m", 1);
      getQuarkGluon(_h["_h_Table95"], _h["_h_Table77"], _h["_h_Table113"], _h["_h_Table131"], 2, "m", 2);
      getQuarkGluon(_h["_h_Table96"], _h["_h_Table78"], _h["_h_Table114"], _h["_h_Table132"], 2, "m", 2);
      getQuarkGluon(_h["_h_Table97"], _h["_h_Table79"], _h["_h_Table115"], _h["_h_Table133"], 2, "zg", 0);
      getQuarkGluon(_h["_h_Table98"], _h["_h_Table80"], _h["_h_Table116"], _h["_h_Table134"], 2, "zg", 0);
      getQuarkGluon(_h["_h_Table99"], _h["_h_Table81"], _h["_h_Table117"], _h["_h_Table135"], 2, "zg", 1);
      getQuarkGluon(_h["_h_Table100"], _h["_h_Table82"], _h["_h_Table118"], _h["_h_Table136"], 2, "zg", 1);
      getQuarkGluon(_h["_h_Table101"], _h["_h_Table83"], _h["_h_Table119"], _h["_h_Table137"], 2, "zg", 2);
      getQuarkGluon(_h["_h_Table102"], _h["_h_Table84"], _h["_h_Table120"], _h["_h_Table138"], 2, "zg", 2);
      getQuarkGluon(_h["_h_Table103"], _h["_h_Table85"], _h["_h_Table121"], _h["_h_Table139"], 2, "rg", 0);
      getQuarkGluon(_h["_h_Table104"], _h["_h_Table86"], _h["_h_Table122"], _h["_h_Table140"], 2, "rg", 0);
      getQuarkGluon(_h["_h_Table105"], _h["_h_Table87"], _h["_h_Table123"], _h["_h_Table141"], 2, "rg", 1);
      getQuarkGluon(_h["_h_Table106"], _h["_h_Table88"], _h["_h_Table124"], _h["_h_Table142"], 2, "rg", 1);
      getQuarkGluon(_h["_h_Table107"], _h["_h_Table89"], _h["_h_Table125"], _h["_h_Table143"], 2, "rg", 2);
      getQuarkGluon(_h["_h_Table108"], _h["_h_Table90"], _h["_h_Table126"], _h["_h_Table144"], 2, "rg", 2);
    }

  protected:

    size_t normBin1, normBin2;
    vector<double> gluonFractionCentral, gluonFractionForward;
    vector<double> rgBins, zgBins, rhoBins, ptBins, zgBinsBeta0;

  private:
    map<string, Histo1DPtr> _h;

  };

  DECLARE_RIVET_PLUGIN(ATLAS_2019_I1772062);
}
