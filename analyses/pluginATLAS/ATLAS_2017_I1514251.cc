// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {
  
  /// Z + jets in pp at 13 TeV 
  /// @note This base class contains a "mode" variable for combined, e, and mu channel derived classes
  class ATLAS_2017_I1514251 : public Analysis {
  public:
    
    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2017_I1514251);
    
    /// Book histograms and initialise projections before the run
    void init() {

      // Get options from the new option system
      // default to combined.
      _mode = 2;
      if ( getOption("LMODE") == "EL" ) _mode = 0;
      if ( getOption("LMODE") == "MU" ) _mode = 1;
      if ( getOption("LMODE") == "EMU" ) _mode = 2;

      const FinalState fs;
      
      Cut cuts = (Cuts::pT > 25*GeV) && (Cuts::abseta < 2.5);
      
      ZFinder zeefinder(fs, cuts, PID::ELECTRON, 71*GeV, 111*GeV);
      ZFinder zmumufinder(fs, cuts, PID::MUON, 71*GeV, 111*GeV);
      declare(zeefinder, "zeefinder");
      declare(zmumufinder, "zmumufinder");
      
      // Define veto FS in order to prevent Z-decay products entering the jet algorithm
      VetoedFinalState had_fs;
      had_fs.addVetoOnThisFinalState(zeefinder);
      had_fs.addVetoOnThisFinalState(zmumufinder);
      FastJets jets(had_fs, FastJets::ANTIKT, 0.4);
      jets.useInvisibles(true);
      declare(jets, "jets");
      
      // individual channels
      _h_Njets_excl  = bookHisto1D(  _mode + 1, 1, 1);
      _h_Njets       = bookHisto1D(  _mode + 4, 1, 1);
      _h_Njets_Ratio = bookScatter2D(_mode + 7, 1, 1, true);

      _h_leading_jet_pT_eq1jet = bookHisto1D(_mode + 10, 1, 1);
      _h_leading_jet_pT        = bookHisto1D(_mode + 13, 1, 1);
      _h_leading_jet_pT_2jet   = bookHisto1D(_mode + 16, 1, 1);
      _h_leading_jet_pT_3jet   = bookHisto1D(_mode + 19, 1, 1);
      _h_leading_jet_pT_4jet   = bookHisto1D(_mode + 22, 1, 1);
      _h_leading_jet_rap       = bookHisto1D(_mode + 25, 1, 1);
      _h_HT                    = bookHisto1D(_mode + 28, 1, 1);
      _h_jet_dphi              = bookHisto1D(_mode + 31, 1, 1);
      _h_jet_mass              = bookHisto1D(_mode + 34, 1, 1);

    }
    
    
    /// Perform the per-event analysis
    void analyze(const Event& event) {
      
      const double weight = event.weight(); 
      
      const ZFinder& zeefinder = apply<ZFinder>(event, "zeefinder");
      const ZFinder& zmumufinder = apply<ZFinder>(event, "zmumufinder");
      
      const Particles& zees = zeefinder.bosons();
      const Particles& zmumus = zmumufinder.bosons();

      //Veto Z->mumu in electron mode, and vice versa:      
      if (_mode==0 && (zees.size()!=1 || zmumus.size() ) )  vetoEvent;

      if (_mode==1 && (zees.size() || zmumus.size()!=1 ) )  vetoEvent;

      if (zees.size() + zmumus.size() != 1) {
        // Running in combined mode, we did not find exactly one Z. Not good.
        MSG_DEBUG("Did not find exactly one good Z candidate");
        vetoEvent;
      }

      // Find the (dressed!) leptons
      const Particles& leptons = zees.size() ? zeefinder.constituents() : zmumufinder.constituents();
      if (leptons.size() != 2) vetoEvent;

      Jets jets =  apply<JetAlg>(event, "jets").jetsByPt(Cuts::pT > 30*GeV && Cuts::absrap < 2.5);
      
      bool veto = false;
      foreach(const Jet& j, jets)  {
        foreach(const Particle& l, leptons) { veto |= deltaR(j, l) < 0.4; }
      }
      if (veto) vetoEvent;
      
      double HT=0;
      foreach(const Particle& l, leptons) { HT += l.pT(); }

      const size_t Njets = jets.size();
      _h_Njets_excl->fill(Njets, weight);	
      for(size_t i = 0; i <= Njets; ++i) { _h_Njets->fill(i, weight);	}

      if (Njets < 1)  vetoEvent;

      
      for(size_t i = 0; i < Njets; ++i) { HT += jets[i].pT(); }
      const double pT = jets[0].pT();
      const double rap = jets[0].rapidity();
      
      _h_HT->fill(HT, weight);      
      _h_leading_jet_rap->fill(fabs(rap), weight);
      _h_leading_jet_pT->fill(pT, weight);
      if (Njets == 1)  _h_leading_jet_pT_eq1jet->fill(pT, weight);
      if (Njets > 1) {
        _h_leading_jet_pT_2jet->fill(pT, weight);
        _h_jet_dphi->fill( deltaPhi(jets[0], jets[1]), weight );
        _h_jet_mass->fill( (jets[0].momentum()+jets[1].momentum()).mass() , weight );
      }
      
      if (Njets > 2)  _h_leading_jet_pT_3jet->fill(pT, weight);
      if (Njets > 3)  _h_leading_jet_pT_4jet->fill(pT, weight);

    }

    void finalize() {
      for (size_t i = 0; i < _h_Njets->numBins()-2; ++i) {
        double  n = _h_Njets->bin(i + 1).sumW();
        double dN = _h_Njets->bin(i + 1).sumW2();
        double  d = _h_Njets->bin(i).sumW();
        double dD = _h_Njets->bin(i).sumW2();
        double r = safediv(n, d);
        double e = sqrt( safediv(r * (1 - r), d) );
        if ( _h_Njets->effNumEntries() != _h_Njets->numEntries() ) {
          // use F. James's approximation for weighted events:
          e = sqrt( safediv((1 - 2 * r) * dN + r * r * dD, d * d) );
        }
        _h_Njets_Ratio->point(i).setY(r, e);
      }

      // when running in combined mode, need to average to get lepton xsec
      double normfac = crossSectionPerEvent();
      if (_mode == 2) normfac = 0.5*normfac;

      scale(_h_Njets,                  normfac );
      scale(_h_Njets_excl,             normfac );
      scale(_h_HT,                     normfac );
      scale(_h_leading_jet_rap,        normfac );
      scale(_h_leading_jet_pT,         normfac );
      scale(_h_leading_jet_pT_eq1jet,  normfac );
      scale(_h_leading_jet_pT_2jet,    normfac );
      scale(_h_leading_jet_pT_3jet,    normfac );
      scale(_h_leading_jet_pT_4jet,    normfac );
      scale(_h_jet_dphi,               normfac );
      scale(_h_jet_mass,               normfac );

    }

    //@}


  protected:

    size_t _mode;


  private:

    Scatter2DPtr _h_Njets_Ratio;
    Histo1DPtr   _h_Njets;
    Scatter2DPtr _h_Njets_excl_Ratio;
    Histo1DPtr   _h_Njets_excl;
    Histo1DPtr   _h_HT;
    Histo1DPtr   _h_leading_jet_rap;
    Histo1DPtr   _h_leading_jet_pT;
    Histo1DPtr   _h_leading_jet_pT_eq1jet;
    Histo1DPtr   _h_leading_jet_pT_2jet;
    Histo1DPtr   _h_leading_jet_pT_3jet;
    Histo1DPtr   _h_leading_jet_pT_4jet;
    Histo1DPtr   _h_jet_dphi;
    Histo1DPtr   _h_jet_mass;
  
  };

  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1514251);

}
