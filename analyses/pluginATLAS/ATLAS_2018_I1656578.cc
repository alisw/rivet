#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {


  /// ttbar l+jets cross sections at 13 TeV
  class ATLAS_2018_I1656578 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2018_I1656578);


    /// Book cuts and projections
    void init() {
      // Eta ranges
      Cut eta_full = (Cuts::abseta < 5.0);
      Cut lep_cuts = (Cuts::abseta < 2.5) && (Cuts::pT > 25*GeV);

      // All final state particles
      FinalState fs(eta_full);

      // Get photons to dress leptons
      IdentifiedFinalState all_photons(fs);
      all_photons.acceptIdPair(PID::PHOTON);

      PromptFinalState photons(Cuts::abspid == PID::PHOTON, true);
      declare(photons, "photons");

      // Projection to find the electrons
      PromptFinalState electrons(Cuts::abspid == PID::ELECTRON, true);

      DressedLeptons dressedelectrons(photons, electrons, 0.1, lep_cuts);
      declare(dressedelectrons, "elecs");

      DressedLeptons ewdressedelectrons(all_photons, electrons, 0.1, eta_full);

      // Projection to find the muons
      PromptFinalState muons(Cuts::abspid == PID::MUON, true);

      DressedLeptons dressedmuons(photons, muons, 0.1, lep_cuts);
      declare(dressedmuons, "muons");

      DressedLeptons ewdressedmuons(all_photons, muons, 0.1, eta_full);

      // Projection to find MET
      declare(MissingMomentum(fs), "MET");

      // Jet clustering.
      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(ewdressedelectrons);
      vfs.addVetoOnThisFinalState(ewdressedmuons);
      FastJets jets(vfs, FastJets::ANTIKT, 0.4, JetAlg::Muons::ALL, JetAlg::Invisibles::DECAY);
      declare(jets, "jets");


      book(_h["absPout_inc"],                 114, 1, 1);
      book(_h["absPout_inc_norm"],            115, 1, 1);
      book(_h["ptpseudotophadron_r1"],        98, 1, 1);
      book(_h["ptpseudotophadron_r1_norm"],   99, 1, 1);
      book(_h["ptttbar_r1"],                  100, 1, 1);
      book(_h["ptttbar_r1_norm"],             101, 1, 1);
      book(_h["absPout_r1"],                  96, 1, 1);
      book(_h["absPout_r1_norm"],             97, 1, 1);
      book(_h["ptpseudotophadron_r2"],        110, 1, 1);
      book(_h["ptpseudotophadron_r2_norm"],   111, 1, 1);
      book(_h["ptttbar_r2"],                  112, 1, 1);
      book(_h["ptttbar_r2_norm"],             113, 1, 1);
      book(_h["absPout_r2"],                  108, 1, 1);
      book(_h["absPout_r2_norm"],             109, 1, 1);
      book(_h["ptpseudotophadron_r3"],        104, 1, 1);
      book(_h["ptpseudotophadron_r3_norm"],   105, 1, 1);
      book(_h["ptttbar_r3"],                  106, 1, 1);
      book(_h["ptttbar_r3_norm"],             107, 1, 1);
      book(_h["absPout_r3"],                  102, 1, 1);
      book(_h["absPout_r3_norm"],             103, 1, 1);

    }


    void analyze(const Event& event) {

      // Get the selected objects, using the projections.
      vector<DressedLepton> electrons = apply<DressedLeptons>(event, "elecs").dressedLeptons();
      vector<DressedLepton> muons     = apply<DressedLeptons>(event, "muons").dressedLeptons();
      const Jets& jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);

      const Vector3 met = apply<MissingMomentum>(event, "MET").vectorMPT();

      Jets bjets, lightjets;
      for (Jet jet : jets) {
        bool b_tagged = jet.bTagged(Cuts::pT > 5*GeV);
        if ( b_tagged && bjets.size() < 2) bjets += jet;
        else lightjets += jet;
      }

      bool single_electron = (electrons.size() == 1) && (muons.empty());
      bool single_muon = (muons.size() == 1) && (electrons.empty());


      DressedLepton *lepton = NULL;
      if (single_electron)   lepton = &electrons[0];
      else if (single_muon)  lepton = &muons[0];

      if (!single_electron && !single_muon)  vetoEvent;

      bool num_b_tagged_jets = (bjets.size() == 2);
      if (!num_b_tagged_jets)  vetoEvent;

      if (jets.size() < 4)  vetoEvent;

      bool reg_4jex2bin = (jets.size() == 4);
      bool reg_5jex2bin = (jets.size() == 5);
      bool reg_6jin2bin = (jets.size() >= 6);


      FourMomentum pbjet1; //Momentum of bjet1
      FourMomentum pbjet2; //Momentum of bjet

      if ( deltaR(bjets[0], *lepton) <= deltaR(bjets[1], *lepton) ) {
        pbjet1 = bjets[0].momentum();
        pbjet2 = bjets[1].momentum();
      } else {
        pbjet1 = bjets[1].momentum();
        pbjet2 = bjets[0].momentum();
      }

      double bestWmass = 1000.0*TeV;
      double mWPDG = 80.399*GeV;
      int Wj1index = -1, Wj2index = -1;
      for (unsigned int i = 0; i < (lightjets.size() - 1); ++i) {
        for (unsigned int j = i + 1; j < lightjets.size(); ++j) {
          double wmass = (lightjets[i].momentum() + lightjets[j].momentum()).mass();
          if (fabs(wmass - mWPDG) < fabs(bestWmass - mWPDG)) {
            bestWmass = wmass;
            Wj1index = i;
            Wj2index = j;
          }
        }
      }

      FourMomentum pjet1 = lightjets[Wj1index].momentum();
      FourMomentum pjet2 = lightjets[Wj2index].momentum();

      // compute hadronic W boson
      FourMomentum pWhadron = pjet1 + pjet2;
      double pz = computeneutrinoz(lepton->momentum(), met);
      FourMomentum ppseudoneutrino( sqrt(sqr(met.x()) + sqr(met.y()) + sqr(pz)), met.x(), met.y(), pz);

      //compute leptonic, hadronic, combined pseudo-top
      FourMomentum ppseudotoplepton = lepton->momentum() + ppseudoneutrino + pbjet1;
      FourMomentum ppseudotophadron = pbjet2 + pWhadron;
      FourMomentum pttbar = ppseudotoplepton + ppseudotophadron;

      Vector3 z_versor(0,0,1);
      Vector3 vpseudotophadron = ppseudotophadron.vector3();
      Vector3 vpseudotoplepton = ppseudotoplepton.vector3();
      // Variables
      double absPout = fabs(vpseudotophadron.dot((vpseudotoplepton.cross(z_versor))/(vpseudotoplepton.cross(z_versor).mod())));

      //pseudotop hadrons and leptons fill histogram
      if (reg_4jex2bin) {
        _h["ptpseudotophadron_r1"]->fill(ppseudotophadron.pt()); //pT of pseudo top hadron
        _h["ptpseudotophadron_r1_norm"]->fill(ppseudotophadron.pt()); //pT of pseudo top hadron
        _h["ptttbar_r1"]->fill(pttbar.pt()); //fill pT of ttbar in combined channel
        _h["ptttbar_r1_norm"]->fill(pttbar.pt()); //fill pT of ttbar in combined channel
        _h["absPout_r1"]->fill(absPout);
        _h["absPout_r1_norm"]->fill(absPout);
      }
      if (reg_5jex2bin) {
        _h["ptpseudotophadron_r2"]->fill(ppseudotophadron.pt()); //pT of pseudo top hadron
        _h["ptpseudotophadron_r2_norm"]->fill(ppseudotophadron.pt()); //pT of pseudo top hadron
        _h["ptttbar_r2"]->fill(pttbar.pt()); //fill pT of ttbar in combined channel
        _h["ptttbar_r2_norm"]->fill(pttbar.pt()); //fill pT of ttbar in combined channel
        _h["absPout_r2"]->fill(absPout);
        _h["absPout_r2_norm"]->fill(absPout);
      }
      if (reg_6jin2bin) {
        _h["ptpseudotophadron_r3"]->fill(ppseudotophadron.pt()); //pT of pseudo top hadron
        _h["ptpseudotophadron_r3_norm"]->fill(ppseudotophadron.pt()); //pT of pseudo top hadron
        _h["ptttbar_r3"]->fill(pttbar.pt()); //fill pT of ttbar in combined channel
        _h["ptttbar_r3_norm"]->fill(pttbar.pt()); //fill pT of ttbar in combined channel
        _h["absPout_r3"]->fill(absPout);
        _h["absPout_r3_norm"]->fill(absPout);
      }
      _h["absPout_inc"]->fill(absPout);
      _h["absPout_inc_norm"]->fill(absPout);
    }


    void finalize() {
      // Normalize to cross-section
      const double sf = (crossSection() / sumOfWeights());
      for (auto hist : _h) {
        scale(hist.second, sf);
        // Normalized distributions
        if (hist.first.find("_norm") != string::npos)  normalize(hist.second);
      }
    }


    double computeneutrinoz(const FourMomentum& lepton, const Vector3 &met) const {
      // computing z component of neutrino momentum given lepton and met
      double pzneutrino;
      double m_W = 80.399; // in GeV, given in the paper
      double k = (( sqr( m_W ) - sqr( lepton.mass() ) ) / 2 ) + (lepton.px() * met.x() + lepton.py() * met.y());
      double a = sqr ( lepton.E() )- sqr ( lepton.pz() );
      double b = -2*k*lepton.pz();
      double c = sqr( lepton.E() ) * sqr( met.mod() ) - sqr( k );
      double discriminant = sqr(b) - 4 * a * c;
      double quad[2] = { (- b - sqrt(discriminant)) / (2 * a), (- b + sqrt(discriminant)) / (2 * a) }; //two possible quadratic solns
      if (discriminant < 0)  pzneutrino = - b / (2 * a); //if the discriminant is negative
      else { //if the discriminant is greater than or equal to zero, take the soln with smallest absolute value
        double absquad[2];
        for (int n=0; n<2; ++n)  absquad[n] = fabs(quad[n]);
        if (absquad[0] < absquad[1])  pzneutrino = quad[0];
        else                          pzneutrino = quad[1];
      }
      return pzneutrino;
    }


  private:

    /// @name Objects that are used by the event selection decisions
    map<string, Histo1DPtr> _h;

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2018_I1656578);


}
