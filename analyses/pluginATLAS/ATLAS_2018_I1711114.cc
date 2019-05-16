#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// g -> bb at 13 TeV
  class ATLAS_2018_I1711114: public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2018_I1711114);


    /// Book cuts and projections
    void init() {
      // All final state particles
      FinalState fs(Cuts::abseta<5.0);

      ChargedFinalState cfs(Cuts::pT > 0.5*GeV && Cuts::abseta < 2.5);
      FastJets smallR_jets(cfs, FastJets::ANTIKT, 0.2, JetAlg::NO_MUONS, JetAlg::NO_INVISIBLES);
      declare(smallR_jets, "track_jets");

      FastJets largeR_jets(fs, FastJets::ANTIKT, 1.0, JetAlg::NO_MUONS, JetAlg::NO_INVISIBLES);
      declare(largeR_jets, "largeR_jets");

      _h_R   = bookHisto1D(1,1,1);
      _h_phi = bookHisto1D(2,1,1);
      _h_z   = bookHisto1D(3,1,1);
      _h_rho = bookHisto1D(4,1,1);
    }


    void analyze(const Event& event) {

      const double weight = event.weight();

      const PseudoJets& myJets = apply<FastJets>(event, "largeR_jets").pseudoJetsByPt(450*GeV);
      if (myJets.empty()) vetoEvent;

      fastjet::Filter trimmer(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.05));
      Jets myTrimmedJets;
      for (const Jet& jet : myJets)  myTrimmedJets += Jet(trimmer(jet));
      std::sort(myTrimmedJets.begin(), myTrimmedJets.end(), cmpMomByPt);
      if (myTrimmedJets[0].pT() < 450*GeV) vetoEvent;

      const Jets& myJets_charged = apply<FastJets>(event, "track_jets").jetsByPt(Cuts::pT > 10*GeV);
      PseudoJets pjs;
      for (auto tj : myTrimmedJets[0].pseudojet().constituents()) {
        fastjet::PseudoJet pj = tj;
        pj.set_user_index(-1); // dummmy
        pjs.push_back(pj);
      }
      for (size_t i = 0; i < myJets_charged.size(); ++i) {
        fastjet::PseudoJet pj = myJets_charged[i];
        pj *= 1e-20; // ghostify momentum
        pj.set_user_index(i);
        pjs.push_back(pj);
      }
      fastjet::ClusterSequence tagged_seq(pjs, fastjet::JetDefinition(fastjet::antikt_algorithm, 1.0));
      PseudoJets GAfatjet = fastjet::sorted_by_pt(tagged_seq.inclusive_jets(50.0));
      Jets associated_trackjets;
      for (auto pj : GAfatjet[0].constituents()) {
        if (pj.user_index() >= 0)  associated_trackjets += myJets_charged[pj.user_index()];
      }
      if (associated_trackjets.size() < 2) vetoEvent;
      std::sort(associated_trackjets.begin(), associated_trackjets.end(), cmpMomByPt);

      size_t nbtags = 0;
      for (unsigned int ij = 0; ij < 2; ++ij) {
        if (associated_trackjets[ij].bTagged(Cuts::pT > 5*GeV))  ++nbtags;
      }
      if (nbtags != 2) vetoEvent;

      // Now, time to compute the observables!
      const Vector3 b1v = associated_trackjets[0].p3();
      const Vector3 b2v = associated_trackjets[1].p3();
      const Vector3 plane1 = b1v.cross(b2v);
      const Vector3 gv = b1v + b2v;
      const Vector3 beam = Vector3(0,0,1);
      const Vector3 plane2 = beam.cross(gv);

      const double fTz = associated_trackjets[1].pT() / (associated_trackjets[0].pT()+associated_trackjets[1].pT());
      const double fTR = deltaR(associated_trackjets[0], associated_trackjets[1], RapScheme::YRAP); //N.B. this uses y and not eta
      const FourMomentum dijet = associated_trackjets[0].mom() + associated_trackjets[1].mom();
      const double fTrho = log(dijet.mass() / dijet.pT());
      const double fTphi = (plane1).angle(plane2)/M_PI;

      _h_R->fill(  fTR,   weight);
      _h_phi->fill(fTphi, weight);
      _h_z->fill(  fTz,   weight);
      _h_rho->fill(fTrho, weight);
    }


    /// Scale histos
    void finalize() {
      _h_R->normalize(  1.0, false);
      _h_phi->normalize(1.0, false);
      _h_z->normalize(  1.0, false);
      _h_rho->normalize(1.0, false);
    }


  private:

    /// Histograms
    Histo1DPtr _h_z, _h_R, _h_rho, _h_phi;

  };


  DECLARE_RIVET_PLUGIN(ATLAS_2018_I1711114);

}
