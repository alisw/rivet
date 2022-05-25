// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


   /// CMS ratio of 3-jet to 2-jet cross-sections
   class CMS_2011_S9088458 : public Analysis {
   public:

     RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2011_S9088458);


     void init() {
       FinalState fs;
       FastJets akt(fs, FastJets::ANTIKT, 0.5);
       declare(akt, "antikT");

       book(_h_tmp_dijet , "TMP/dijet", refData(1, 1, 1));
       book(_h_tmp_trijet, "TMP/trijet", refData(1, 1, 1));
       book(_h_r32, 1, 1, 1);
     }


     void analyze(const Event & event) {
       const double weight = 1.0;

       Jets highpT_jets;
       double HT = 0;
       for(const Jet & jet : apply<JetAlg>(event, "antikT").jetsByPt(50.0*GeV)) {
         if (jet.abseta() < 2.5) {
           highpT_jets.push_back(jet);
           HT += jet.pT();
         }
       }
       if (highpT_jets.size() < 2) vetoEvent;
       if (highpT_jets.size() >= 2) _h_tmp_dijet->fill(HT/TeV, weight);
       if (highpT_jets.size() >= 3) _h_tmp_trijet->fill(HT/TeV, weight);
     }


     void finalize() {
       divide(_h_tmp_trijet, _h_tmp_dijet, _h_r32);
     }


   private:

     /// @{
     Histo1DPtr _h_tmp_dijet, _h_tmp_trijet;
     Scatter2DPtr _h_r32;
     /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CMS_2011_S9088458, CMS_2011_I912560);

}
