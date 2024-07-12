// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// CMS Azimuthal corellations at 13 TeV
  class CMS_2018_I1643640 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2018_I1643640);

    void init() {
      FinalState fs;
      FastJets akt(fs, FastJets::ANTIKT, 0.4);
      declare(akt, "antikT");

      Histo1DPtr tmp;

      _h_deltaPhi_2J_phi12.add( 200.,  300., book(tmp, 1, 1, 1));
      _h_deltaPhi_2J_phi12.add( 300.,  400., book(tmp, 2, 1, 1));
      _h_deltaPhi_2J_phi12.add( 400.,  500., book(tmp, 3, 1, 1));
      _h_deltaPhi_2J_phi12.add( 500.,  600., book(tmp, 4, 1, 1));
      _h_deltaPhi_2J_phi12.add( 600.,  700., book(tmp, 5, 1, 1));
      _h_deltaPhi_2J_phi12.add( 700.,  800., book(tmp, 6, 1, 1));
      _h_deltaPhi_2J_phi12.add( 800.,  1000., book(tmp, 7, 1, 1));
      _h_deltaPhi_2J_phi12.add( 1000., 1200., book(tmp, 8, 1, 1));
      _h_deltaPhi_2J_phi12.add( 1200., 7000., book(tmp, 9, 1, 1));

      _h_deltaPhi_3J_phi12.add( 200.,  300., book(tmp, 10, 1, 1));
      _h_deltaPhi_3J_phi12.add( 300.,  400., book(tmp, 11, 1, 1));
      _h_deltaPhi_3J_phi12.add( 400.,  500., book(tmp, 12, 1, 1));
      _h_deltaPhi_3J_phi12.add( 500.,  600., book(tmp, 13, 1, 1));
      _h_deltaPhi_3J_phi12.add( 600.,  700., book(tmp, 14, 1, 1));
      _h_deltaPhi_3J_phi12.add( 700.,  800., book(tmp, 15, 1, 1));
      _h_deltaPhi_3J_phi12.add( 800.,  1000., book(tmp, 16, 1, 1));
      _h_deltaPhi_3J_phi12.add( 1000., 7000., book(tmp, 17, 1, 1));

      _h_deltaPhi_4J_phi12.add( 200.,  300., book(tmp, 18, 1, 1));
      _h_deltaPhi_4J_phi12.add( 300.,  400., book(tmp, 19, 1, 1));
      _h_deltaPhi_4J_phi12.add( 400.,  500., book(tmp, 20, 1, 1));
      _h_deltaPhi_4J_phi12.add( 500.,  600., book(tmp, 21, 1, 1));
      _h_deltaPhi_4J_phi12.add( 600.,  700., book(tmp, 22, 1, 1));
      _h_deltaPhi_4J_phi12.add( 700.,  800., book(tmp, 23, 1, 1));
      _h_deltaPhi_4J_phi12.add( 800.,  1000., book(tmp, 24, 1, 1));
      _h_deltaPhi_4J_phi12.add( 1000., 7000., book(tmp, 25, 1, 1));

      _h_deltaPhi_3J_phimin2J.add( 200.,  300., book(tmp, 26, 1, 1));
      _h_deltaPhi_3J_phimin2J.add( 300.,  400., book(tmp, 27, 1, 1));
      _h_deltaPhi_3J_phimin2J.add( 400.,  500., book(tmp, 28, 1, 1));
      _h_deltaPhi_3J_phimin2J.add( 500.,  600., book(tmp, 29, 1, 1));
      _h_deltaPhi_3J_phimin2J.add( 600.,  700., book(tmp, 30, 1, 1));
      _h_deltaPhi_3J_phimin2J.add( 700.,  800., book(tmp, 31, 1, 1));
      _h_deltaPhi_3J_phimin2J.add( 800.,  1000., book(tmp, 32, 1, 1));
      _h_deltaPhi_3J_phimin2J.add( 1000., 7000., book(tmp, 33, 1, 1));

      _h_deltaPhi_4J_phimin2J.add( 200.,  300., book(tmp, 34, 1, 1));
      _h_deltaPhi_4J_phimin2J.add( 300.,  400., book(tmp, 35, 1, 1));
      _h_deltaPhi_4J_phimin2J.add( 400.,  500., book(tmp, 36, 1, 1));
      _h_deltaPhi_4J_phimin2J.add( 500.,  600., book(tmp, 37, 1, 1));
      _h_deltaPhi_4J_phimin2J.add( 600.,  700., book(tmp, 38, 1, 1));
      _h_deltaPhi_4J_phimin2J.add( 700.,  800., book(tmp, 39, 1, 1));
      _h_deltaPhi_4J_phimin2J.add( 800.,  1000., book(tmp, 40, 1, 1));
      _h_deltaPhi_4J_phimin2J.add( 1000., 7000., book(tmp, 41, 1, 1));

    }


    void analyze(const Event & event) {
      const Jets& jets = apply<JetAlg>(event, "antikT").jetsByPt();

      // 2 jet case and Delta_phi12
      if( jets.size() >= 2 ) {
        if ( (jets[0].pT() >= 200.*GeV)  &&  (jets[1].pT() >= 100.*GeV) ) {
          if ( (fabs(jets[0].rap()) <= 2.5)  &&  (fabs(jets[1].rap()) <= 2.5) ) {
            double dphi = deltaPhi(jets[0].phi(), jets[1].phi());
            _h_deltaPhi_2J_phi12.fill(jets[0].pT(), dphi, 1.0);
          }
        }
      }

      // 3 jet case and Delta_phi12
      if ( jets.size() >= 3 ) {
        if ( (jets[0].pT() >= 200.*GeV)  &&  (jets[1].pT() >= 100.*GeV)  && (jets[2].pT() >= 100.*GeV) ) {
          if ( (fabs(jets[0].rap()) <= 2.5)  &&  (fabs(jets[1].rap()) <= 2.5) &&  (fabs(jets[2].rap()) <= 2.5)) {
            double dphi = deltaPhi(jets[0].phi(), jets[1].phi());
            _h_deltaPhi_3J_phi12.fill(jets[0].pT(), dphi, 1.0);
          }
        }
      }

      // 4 jet case and Delta_phi12
      if ( jets.size() >= 4 ) {
        if ( (jets[0].pT() >= 200.*GeV)  &&  (jets[1].pT() >= 100.*GeV)  && (jets[2].pT() >= 100.*GeV)   && (jets[3].pT() >= 100.*GeV)) {
          if ( (fabs(jets[0].rap()) <= 2.5)  &&  (fabs(jets[1].rap()) <= 2.5) &&  (fabs(jets[2].rap()) <= 2.5) &&  (fabs(jets[3].rap()) <= 2.5)) {
            double dphi = deltaPhi(jets[0].phi(), jets[1].phi());
            _h_deltaPhi_4J_phi12.fill(jets[0].pT(), dphi, 1.0);
          }
        }
      }

      // 3 jet case and Delta_Phi_min2j
      if ( jets.size() >= 3 ) {
        if ( (jets[0].pT() >= 200.*GeV)  &&  (jets[1].pT() >= 100.*GeV)  && (jets[2].pT() >= 100.*GeV) ) {
          if ( (fabs(jets[0].rap()) <= 2.5)  &&  (fabs(jets[1].rap()) <= 2.5) &&  (fabs(jets[2].rap()) <= 2.5)) {
            double dphi01 = deltaPhi(jets[0].phi(), jets[1].phi());
            if (dphi01 >= PI/2. ){
              double dphi02 = deltaPhi(jets[0].phi(), jets[2].phi());
              double dphi12 = deltaPhi(jets[1].phi(), jets[2].phi());
              // evaluate DPhi2Jmin
              vector<double> Dphis2J{dphi01,dphi02,dphi12};
              double DPhi2Jmin = min(Dphis2J);
              // double Dphis2J[3] = {dphi01,dphi02,dphi12};
              // double DPhi2Jmin = Dphis2J[0];
              // for (int gg=1; gg<3; ++gg) { if (DPhi2Jmin>Dphis2J[gg]) DPhi2Jmin = Dphis2J[gg]; }
              _h_deltaPhi_3J_phimin2J.fill(jets[0].pT(), DPhi2Jmin, 1.0);
            }
          }
        }
      }

      // 4 jet case and Delta_Phi_min2j
      if ( jets.size() >= 4 ) {
        if ( (jets[0].pT() >= 200.*GeV)  &&  (jets[1].pT() >= 100.*GeV)  && (jets[2].pT() >= 100.*GeV)   && (jets[3].pT() >= 100.*GeV)) {
          if ( (fabs(jets[0].rap()) <= 2.5)  &&  (fabs(jets[1].rap()) <= 2.5) &&  (fabs(jets[2].rap()) <= 2.5) &&  (fabs(jets[3].rap()) <= 2.5)) {
            double dphi01 = deltaPhi(jets[0].phi(), jets[1].phi());
            if (dphi01 >= PI/2.) {
              double dphi02 = deltaPhi(jets[0].phi(), jets[2].phi());
              double dphi03 = deltaPhi(jets[0].phi(), jets[3].phi());
              double dphi12 = deltaPhi(jets[1].phi(), jets[2].phi());
              double dphi13 = deltaPhi(jets[1].phi(), jets[3].phi());
              double dphi23 = deltaPhi(jets[2].phi(), jets[3].phi());
              /// evaluate DPhi2Jmin
              // double Dphis2J[6]={dphi01,dphi02,dphi03,dphi12,dphi13,dphi23};
              // double DPhi2Jmin=Dphis2J[0];
              // for(int gg=1; gg<6; ++gg){ if(DPhi2Jmin>Dphis2J[gg]){DPhi2Jmin=Dphis2J[gg];} }
              vector<double> Dphis2J{dphi01,dphi02,dphi03,dphi12,dphi13,dphi23};
              double DPhi2Jmin = min(Dphis2J);
              _h_deltaPhi_4J_phimin2J.fill(jets[0].pT(), DPhi2Jmin, 1.0);
            }
          }
        }
      }
    }  // end analyze


    void finalize() {
      for (Histo1DPtr histo : _h_deltaPhi_2J_phi12.histos()) normalize(histo);
      for (Histo1DPtr histo : _h_deltaPhi_3J_phi12.histos()) normalize(histo);
      for (Histo1DPtr histo : _h_deltaPhi_4J_phi12.histos()) normalize(histo);
      for (Histo1DPtr histo : _h_deltaPhi_3J_phimin2J.histos()) normalize(histo);
      for (Histo1DPtr histo : _h_deltaPhi_4J_phimin2J.histos()) normalize(histo);
    }


  private:

    BinnedHistogram _h_deltaPhi_2J_phi12;
    BinnedHistogram _h_deltaPhi_3J_phi12;
    BinnedHistogram _h_deltaPhi_4J_phi12;
    BinnedHistogram _h_deltaPhi_3J_phimin2J;
    BinnedHistogram _h_deltaPhi_4J_phimin2J;

  };


  RIVET_DECLARE_PLUGIN(CMS_2018_I1643640);

}
