// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/WFinder.hh"

namespace Rivet {


  class ATLAS_2017_I1517194 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    ///@brief: Electroweak Wjj production at 8 TeV
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2017_I1517194);
    //@}

    /// @name Analysis methods
    //@{
    /// Book histograms and initialise projections before the run
    void init() {

      // Get options from the new option system
      _mode = 0;
      if ( getOption("LMODE") == "EL" ) _mode = 0;
      if ( getOption("LMODE") == "MU" ) _mode = 1;

      const FinalState fs;
      // W Selection
      WFinder wfinder(fs, Cuts::rap < 2.5 && Cuts::pT >= 25*GeV, _mode? PID::MUON : PID::ELECTRON, 0*GeV, 13*TeV, 0*GeV, 0.1);
      declare(wfinder, "WFinder");

      FastJets jets( wfinder.remainingFinalState(), FastJets::ANTIKT, 0.4, JetAlg::Muons::DECAY, JetAlg::Invisibles::ALL);
      declare(jets, "Jets_w");

      MissingMomentum missmom(FinalState(Cuts::eta < 5.0));
      declare(missmom, "mm");

      const vector<string> phase_spaces = { "highmass15", "antiLC", "signal10",
                                            "highmass10", "inclusive", "highmass20",
                                            "antiLCantiJC", "antiJC", "signal",  };
      const vector<string> variables = { "dijetmass", "dijetpt", "dphi12", "dy12", "j1pt", "JC", "LC", "ngapjets" };
      size_t hepdataid = 10;
      for (size_t group = 0; group < 4; ++group) {
        for ( size_t ps=0; ps < phase_spaces.size(); ++ps ) {
          for ( size_t var=0; var < variables.size(); ++var ) {
            if (group < 2) {
              if ((ps == 0 || ps == 2 || ps == 3 || ps == 5) && var == 0)  continue;
              if ((ps == 1 || ps == 2 || ps  > 5) && var  > 4)  continue;
              if (group == 1 && ps == 7 && var == 3)  continue;
            }
            else {
              if (ps == 1 || ps == 4 || ps > 5)  continue;
              if ((ps == 0 || ps == 5) && var  < 2)  continue;
              if (ps == 2 && var  > 4)  continue;
              if ((ps == 0 || ps == 3 || ps == 5) && var == 5)  continue;
              if (group == 2) {
                if ((ps == 0 || ps == 5) && var == 3)  continue;
                if (ps == 3 && var == 1)  continue;
              }
              else {
                if ((ps == 0 || ps == 3 || ps == 5) && var == 6)  continue;
                if ((ps == 2 || ps == 3) && var  < 2)  continue;
              }
            }

            ++hepdataid;

            string label = variables[var]+"_"+phase_spaces[ps];
            //string pre = group > 1? "ew_" : "";
            //std::cout << "rivet -- " << pre << label << suff << " " << hepdataid << std::endl;
            string suff = group % 2? "" : "_norm";
            if (group > 1)  book(_hists["ew_" + label + suff], hepdataid, 1, 1);
            else            book(_hists[label + suff], hepdataid, 1, 1);
          }
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      FourMomentum boson, lepton, neutrino;
      const WFinder& wfinder = apply<WFinder>(event, "WFinder");
      const FastJets& jetpro = apply<FastJets>(event, "Jets_w");
      const MissingMomentum& missmom = apply<MissingMomentum>(event, "mm");

      if ( wfinder.bosons().size() != 1 ) { vetoEvent; }

      boson    = wfinder.bosons().front().momentum();
      lepton   = wfinder.constituentLeptons().front().momentum();
      neutrino = wfinder.constituentNeutrinos().front().momentum();

      vector<FourMomentum> jets;
      for (const Jet& jet : jetpro.jetsByPt(30*GeV)) {
        if ( fabs(jet.momentum().rapidity()) > 4.4 ) continue;
        if ( fabs(deltaR(jet, lepton)) < 0.3 ) continue;
        jets.push_back(jet.momentum());
      }

      if (jets.size() < 2)  vetoEvent;

      double mT = sqrt( 2.0*lepton.pT()*neutrino.Et()*( 1.0-cos( lepton.phi()-neutrino.phi() ) ) );
      double DeltaRap = fabs( jets[0].rapidity()-jets[1].rapidity() );
      double dijet_mass = FourMomentum(jets[0]+jets[1]).mass();
      size_t nojets = jets.size();

      if (missmom.vectorEt().mod() < 25*GeV)  vetoEvent;
      if (jets[0].pT() < 80*GeV)  vetoEvent;
      if (jets[1].pT() < 60*GeV)  vetoEvent;
      if (dijet_mass < 500*GeV)  vetoEvent;
      if (DeltaRap < 2.0)  vetoEvent;
      if (mT < 40*GeV)  vetoEvent;

      // By now, the event has passed all VBF cuts
      double DijetPt = FourMomentum(jets[0]+jets[1]).pT();
      double dphi = fabs(jets[0].phi() - jets[1].phi());
      double DeltaPhi = ( dphi<=pi ) ? dphi/pi : (2.*pi-dphi)/pi;
      double jet_1_rap = jets[0].rapidity();
      double jet_2_rap = jets[1].rapidity();

      // Njets in Gap Control Regions info
      int njetsingap = 0;
      bool firstIsForward = jet_1_rap > jet_2_rap;
      int jF = (firstIsForward) ? 0 : 1; // sets most positive jet to be forward index
      int jB = (firstIsForward) ? 1 : 0; // sets most negative jet to be backward index
      for (size_t j = 2; j < nojets; ++j) {
        if ( (jets[j].rapidity()<jets[jF].rapidity()) && (jets[j].rapidity()>jets[jB].rapidity()) ) njetsingap++;
      }

      // Third+ Jet Centrality Cut (Vetos any event that has a jet between raps of the two leading jets)
      bool passJC = false;
      std::vector<double> JCvals;
      JCvals.clear();
      if ( nojets < 3) passJC = true;
      else if ( nojets > 2 ) {
        passJC = true;
        for (size_t j = 2; j < nojets; ++j) {
          double jet_3_rap = jets[j].rapidity();
          double jet3gap = fabs(( jet_3_rap - ((jet_1_rap + jet_2_rap)/2.0))/(jet_1_rap - jet_2_rap));
          JCvals.push_back(jet3gap);
          if ( jet3gap < 0.4 ) { passJC = false; }
        }
      }

      double lepton_rap = lepton.rapidity();
      double lep_cent = fabs((lepton_rap - ((jet_1_rap + jet_2_rap)/2.0) )/(jet_1_rap - jet_2_rap));
      bool passLC = (lep_cent < 0.4);

      map<string, bool> phaseSpaces;
      phaseSpaces["inclusive"] = true;
      phaseSpaces["highmass10"] = (dijet_mass>1000.0*GeV);
      phaseSpaces["highmass15"] = (dijet_mass>1500.0*GeV);
      phaseSpaces["highmass20"] = (dijet_mass>2000.0*GeV);
      phaseSpaces["antiLC"] = ( !passLC && passJC );
      phaseSpaces["antiJC"] = ( passLC && !passJC );
      phaseSpaces["antiLCantiJC"] = ( !passLC && !passJC );
      phaseSpaces["signal"] = ( passLC && passJC );
      phaseSpaces["signal10"] = ( (dijet_mass>1000.0*GeV) && passLC && passJC );

      for (const auto& ps : phaseSpaces) {
        if (!ps.second)  continue;
        fillHisto("dijetmass_"+ps.first, dijet_mass);
        fillHisto("dijetpt_"+ps.first, DijetPt);
        fillHisto("dy12_"+ps.first, fabs(jet_1_rap-jet_2_rap));
        fillHisto("dphi12_"+ps.first, DeltaPhi);
        fillHisto("j1pt_"+ps.first, jets[0].pT());
        if (ps.first == "inclusive" || ps.first.find("highmass") != string::npos) {
          fillHisto("LC_"+ps.first, lep_cent);
          fillHisto("ngapjets_"+ps.first, njetsingap);
          for (auto& jc : JCvals) {
            fillHisto("JC_"+ps.first, jc);
          }
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double factor = crossSection()/sumOfWeights()/femtobarn;  // refData is in fb
      for (const auto& key_hist : _hists) {
        scale(key_hist.second, factor);
        if (key_hist.first.find("_norm") != string::npos) normalize(key_hist.second);
      }
    }

    //@}


    void fillHisto(const string& label, const double value) {
      if (_hists.find(label) != _hists.end()) { // QCD+EW, absolute
        _hists[label]->fill(value);
      }
      if (_hists.find(label + "_norm") != _hists.end()) { // QCD+EW, normalised
        _hists[label + "_norm"]->fill(value);
      }
      if (_hists.find("ew_" + label) != _hists.end()) { // EW-only, absolute
        _hists["ew_" + label]->fill(value);
      }
      if (_hists.find("ew_" + label + "_norm") != _hists.end()) { // EW-only, normalised
        _hists["ew_" + label + "_norm"]->fill(value);
      }
    }


  protected:

    size_t _mode;

  private:

    map<string, Histo1DPtr> _hists;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1517194);


}
