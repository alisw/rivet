// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {

  /// @brief Measurements of angular distance and momentum ratio distributions in three-jet and Z + two-jet final states in pp collisions
  class CMS_2021_I1847230 : public Analysis {
  public:

    CMS_2021_I1847230 ()
      : Analysis("CMS_2021_I1847230")
    {}
    void init() {

      _mode = 0;
      if ( getOption("MODE") == "QCD8TeV" ) _mode = 1;
      else if ( getOption("MODE") == "QCD13TeV" ) _mode = 2;
      else if ( getOption("MODE") == "ZJet" ) _mode = 3;
      if (_mode == 1) {
         _jr = 0.5;
         book(_h1, "d01-x01-y01");
         book(_h2, "d02-x01-y01");
         book(_h3, "d03-x01-y01");
         book(_h4, "d04-x01-y01");
       }
       if (_mode == 2) {
         _jr = 0.4;
         book(_h1, "d05-x01-y01");
         book(_h2, "d06-x01-y01");
         book(_h3, "d07-x01-y01");
         book(_h4, "d08-x01-y01");
       }
       if (_mode == 1 or _mode == 2) {
         const FastJets jets(FinalState(), FastJets::ANTIKT, _jr);
         declare(jets, "jets");
       }
       if (_mode == 3) {
         FinalState fs(Cuts::abseta < 2.4 and Cuts::pT > 100*MeV);
         declare(fs, "FS");
  
         ZFinder zfinder(fs, Cuts::abseta < 5. and Cuts::pT > 30*GeV, PID::MUON, 70*GeV, 110*GeV,
              0.2, ZFinder::ChargedLeptons::PROMPT, ZFinder::ClusterPhotons::NODECAY, 
              ZFinder::AddPhotons::NO, 91.2*GeV);
  
         declare(zfinder, "ZFinder");
         declare(FastJets(zfinder.remainingFinalState(), FastJets::ANTIKT, 0.5), "JetsAK5_zj");
  
         book(_h1, "d09-x01-y01");  
         book(_h2, "d10-x01-y01");  
         book(_h3, "d11-x01-y01");  
         book(_h4, "d12-x01-y01"); 
 
         book(_ZJw_gen, "TMP/ZJw_gen");
       }
    }

    void analyze(const Event& event) {

      if (_mode == 1 or _mode ==2) {
        const Jets& jets = apply<JetAlg>(event, "jets").jetsByPt(Cuts::pT > 30.0*GeV);
        if (jets.size() < 3) vetoEvent;
        const FourMomentum jet1 = jets[0].momentum();
        const FourMomentum jet2 = jets[1].momentum();
        const FourMomentum jet3 = jets[2].momentum();
        
        if (jet1.pT() < 510.0*GeV) vetoEvent;
        if (jet1.absrapidity() > 2.5 or jet2.absrapidity() > 2.5) vetoEvent;
        const double del_phi12 = mapAngle0ToPi(jet2.phi() - jet1.phi());
        if (abs(del_phi12 - M_PI) > 1.0) vetoEvent;
        const double jet3_pt_jet2_pt = jet3.pT()/jet2.pT();
        if (!inRange(jet3_pt_jet2_pt, 0.1, 0.9)) vetoEvent;
        const double del_r23 = deltaR(jet3.rapidity(), jet3.phi(), jet2.rapidity(), jet2.phi());
        if (!inRange(del_r23, _jr+0.1, 1.5)) vetoEvent;
  
        if (del_r23 < 1.0) _h1->fill(jet3_pt_jet2_pt);
        if (del_r23 > 1.0) _h2->fill(jet3_pt_jet2_pt);
        if (jet3_pt_jet2_pt < 0.3) _h3->fill(del_r23);
        if (jet3_pt_jet2_pt > 0.6) _h4->fill(del_r23); 
      }

      if (_mode == 3) { 
        const ZFinder& zfinder = apply<ZFinder>(event, "ZFinder");
        if (zfinder.bosons().size() != 1) vetoEvent;
        const Particle& z = zfinder.bosons()[0];
        const Particles leptons = sortBy(zfinder.constituents(), cmpMomByPt);
        if (leptons[0].pT() < 25.0*GeV or leptons[1].pT() < 10.0*GeV or z.pT() < 80.0*GeV) vetoEvent;
        if (leptons[0].absrapidity() > 2.1 or leptons[1].absrapidity() > 2.4) vetoEvent;
  
        const PseudoJets& psjetsAK5_zj = apply<FastJets>(event, "JetsAK5_zj").pseudoJetsByPt(20.0*GeV);
  
        if (psjetsAK5_zj.empty()) vetoEvent;
        
        const fastjet::PseudoJet& j0 = psjetsAK5_zj[0];
        const FourMomentum jmom0(j0.e(), j0.px(), j0.py(), j0.pz());
        
        if (jmom0.absrapidity() > 1.0 or jmom0.pT() < 80.0*GeV) vetoEvent;
        if (!(deltaPhi(z, jmom0) > 2.0 and deltaR(leptons[0], jmom0) > 0.5 and deltaR(leptons[1], jmom0) > 0.5)) vetoEvent;
        
        _ZJw_gen ->fill(); 

        if(psjetsAK5_zj.size() < 2) vetoEvent;
        
        const fastjet::PseudoJet& j1 = psjetsAK5_zj[1];
        const FourMomentum jmom1(j1.e(), j1.px(), j1.py(), j1.pz());
        if(deltaR(leptons[0], jmom1) < 0.5 or deltaR(leptons[1], jmom1) < 0.5 or jmom1.absrapidity() > 2.4) vetoEvent;
        
        const double dR_gen_Jj = deltaR(jmom0, jmom1);
        if (!inRange(dR_gen_Jj, 0.5, 1.5)) vetoEvent;
        const double rPt_gen_Jj = jmom1.pT()/jmom0.pT(); 
         
        if (dR_gen_Jj < 1.0) _h1->fill(rPt_gen_Jj);
        if (dR_gen_Jj > 1.0) _h2->fill(rPt_gen_Jj);
        if (rPt_gen_Jj < 0.3) _h3->fill(dR_gen_Jj);
        if (rPt_gen_Jj > 0.6) _h4->fill(dR_gen_Jj);
      }
    }

    void finalize() {
      if (_mode == 1 or _mode == 2) {
        normalize(_h1);
        normalize(_h2);
        normalize(_h3);
        normalize(_h4);
        //scale(_h1, 1.0/_h1->effNumEntries());
        //scale(_h2, 1.0/_h2->effNumEntries());
        //scale(_h3, 1.0/_h3->effNumEntries());
        //scale(_h4, 1.0/_h4->effNumEntries());
      }
      if (_mode == 3) {
        for (size_t i = 0; i < _h1->numBins(); i++) {
          _h1->bin(i).scaleW(_h1->bin(i).width());
        }
        for (size_t i = 0; i < _h2->numBins(); i++) {
          _h2->bin(i).scaleW(_h2->bin(i).width());
        }
        for (size_t i = 0; i < _h3->numBins(); i++) {
          _h3->bin(i).scaleW(_h3->bin(i).width());
        }
        for (size_t i = 0; i < _h4->numBins(); i++) {
          _h4->bin(i).scaleW(_h4->bin(i).width());
        }
        scale(_h1, 1.0/ *_ZJw_gen);
        scale(_h2, 1.0/ *_ZJw_gen);
        scale(_h3, 1.0/ *_ZJw_gen);
        scale(_h4, 1.0/ *_ZJw_gen);
      }
    }

  private:

    Histo1DPtr  _h1;
    Histo1DPtr  _h2;
    Histo1DPtr  _h3;
    Histo1DPtr  _h4;
    
    CounterPtr _ZJw_gen;

    double _jr; 

  protected:

    size_t _mode;

  };

  DECLARE_RIVET_PLUGIN(CMS_2021_I1847230);
}
