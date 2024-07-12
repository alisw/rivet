#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PartonicTops.hh"

namespace Rivet {


  class CMS_2015_I1370682_PARTON : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2015_I1370682_PARTON);


    /// Book projections and histograms
    void init() {
      declare(PartonicTops(PartonicTops::DecayMode::E_MU, false), "LeptonicPartonTops");
      declare(PartonicTops(PartonicTops::DecayMode::HADRONIC),    "HadronicPartonTops");

      book(_hSL_topPt        , 15, 1, 1);
      book(_hSL_topPtTtbarSys, 16, 1, 1);
      book(_hSL_topY         , 17, 1, 1);
      book(_hSL_ttbarDelPhi  , 18, 1, 1);
      book(_hSL_topPtLead    , 19, 1, 1);
      book(_hSL_topPtSubLead , 20, 1, 1);
      book(_hSL_ttbarPt      , 21, 1, 1);
      book(_hSL_ttbarY       , 22, 1, 1);
      book(_hSL_ttbarMass    , 23, 1, 1);

      book(_hDL_topPt        , 24, 1, 1);
      book(_hDL_topPtTtbarSys, 25, 1, 1);
      book(_hDL_topY         , 26, 1, 1);
      book(_hDL_ttbarDelPhi  , 27, 1, 1);
      book(_hDL_topPtLead    , 28, 1, 1);
      book(_hDL_topPtSubLead , 29, 1, 1);
      book(_hDL_ttbarPt      , 30, 1, 1);
      book(_hDL_ttbarY       , 31, 1, 1);
      book(_hDL_ttbarMass    , 32, 1, 1);
      }


      void analyze(const Event& event) {

        // Do the analysis only for the ttbar full leptonic or semileptonic channel, without tau decay
        const Particles leptonicpartontops = apply<ParticleFinder>(event, "LeptonicPartonTops").particlesByPt();
        const Particles hadronicpartontops = apply<ParticleFinder>(event, "HadronicPartonTops").particlesByPt();
        const bool isSemilepton = (leptonicpartontops.size() == 1 && hadronicpartontops.size() == 1);
        const bool isDilepton = (leptonicpartontops.size() == 2 && hadronicpartontops.size() == 0);
        if (!isSemilepton && !isDilepton) vetoEvent;

        // Parton level at full phase space
        // Fill top quarks defined in the parton level, full phase space
        const FourMomentum t1P4 = leptonicpartontops[0];
        const FourMomentum t2P4 = isSemilepton ? hadronicpartontops[0] : leptonicpartontops[1];
        const double t1Pt = t1P4.pT(), t2Pt = t2P4.pT();
        const FourMomentum ttbarP4 = t1P4 + t2P4;
        const FourMomentum t1P4AtCM = LorentzTransform::mkFrameTransformFromBeta(ttbarP4.betaVec()).transform(t1P4);
        const double dPhi = deltaPhi(t1P4.phi(), t2P4.phi());

        const double weight = 1.0;
        if (isSemilepton) {
          _hSL_topPt->fill(t1Pt, weight);
          _hSL_topPt->fill(t2Pt, weight);
          _hSL_topPtTtbarSys->fill(t1P4AtCM.pT(), weight);
          _hSL_topY->fill(t1P4.rapidity(), weight);
          _hSL_topY->fill(t2P4.rapidity(), weight);
          _hSL_ttbarDelPhi->fill(dPhi, weight);
          _hSL_topPtLead->fill(std::max(t1Pt, t2Pt), weight);
          _hSL_topPtSubLead->fill(std::min(t1Pt, t2Pt), weight);
          _hSL_ttbarPt->fill(ttbarP4.pT(), weight);
          _hSL_ttbarY->fill(ttbarP4.rapidity(), weight);
          _hSL_ttbarMass->fill(ttbarP4.mass(), weight);
        } else { // if (isDilepton) {
          _hDL_topPt->fill(t1Pt, weight);
          _hDL_topPt->fill(t2Pt, weight);
          _hDL_topPtTtbarSys->fill(t1P4AtCM.pT(), weight);
          _hDL_topY->fill(t1P4.rapidity(), weight);
          _hDL_topY->fill(t2P4.rapidity(), weight);
          _hDL_ttbarDelPhi->fill(dPhi, weight);
          _hDL_topPtLead->fill(std::max(t1Pt, t2Pt), weight);
          _hDL_topPtSubLead->fill(std::min(t1Pt, t2Pt), weight);
          _hDL_ttbarPt->fill(ttbarP4.pT(), weight);
          _hDL_ttbarY->fill(ttbarP4.rapidity(), weight);
          _hDL_ttbarMass->fill(ttbarP4.mass(), weight);
        }
      }


      void finalize() {
        normalize(_hSL_topPt); normalize(_hSL_topPtTtbarSys); normalize(_hSL_topY); 
        normalize(_hSL_ttbarDelPhi); normalize(_hSL_topPtLead); normalize(_hSL_topPtSubLead); 
        normalize(_hSL_ttbarPt); normalize(_hSL_ttbarY); normalize(_hSL_ttbarMass);
        normalize(_hDL_topPt); normalize(_hDL_topPtTtbarSys); normalize(_hDL_topY); 
        normalize(_hDL_ttbarDelPhi); normalize(_hDL_topPtLead); normalize(_hDL_topPtSubLead); 
        normalize(_hDL_ttbarPt); normalize(_hDL_ttbarY); normalize(_hDL_ttbarMass);
      }


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _hSL_topPt, _hSL_topPtTtbarSys, _hSL_topY, _hSL_ttbarDelPhi, _hSL_topPtLead,
      _hSL_topPtSubLead, _hSL_ttbarPt, _hSL_ttbarY, _hSL_ttbarMass;
    Histo1DPtr _hDL_topPt, _hDL_topPtTtbarSys, _hDL_topY, _hDL_ttbarDelPhi, _hDL_topPtLead,
      _hDL_topPtSubLead, _hDL_ttbarPt, _hDL_ttbarY, _hDL_ttbarMass;
    //@}

  };


  RIVET_DECLARE_PLUGIN(CMS_2015_I1370682_PARTON);

}
