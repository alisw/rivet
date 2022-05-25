// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  /// Z + jets in pp at 7 TeV (combined channel / base class)
  /// @note This base class contains a "mode" variable for combined, e, and mu channel derived classes
  class ATLAS_2013_I1230812 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2013_I1230812);
    //@}


    /// Book histograms and initialise projections before the run
    void init() {

      // Get options from the new option system
      _mode = 0;
      if ( getOption("LMODE") == "EL" ) _mode = 1;
      if ( getOption("LMODE") == "MU" ) _mode = 2;

      // Determine the e/mu decay channels used (NB Prompt leptons only).
      /// @todo Note that Zs are accepted with any rapidity: the cuts are on the e/mu: is this correct?
      Cut pt20 = Cuts::pT >= 20*GeV;
      Cut eta_e = _mode? Cuts::abseta < 1.37 || Cuts::absetaIn(1.52, 2.47) : Cuts::abseta < 2.5;
      Cut eta_m = _mode? Cuts::abseta < 2.4 : Cuts::abseta < 2.5;
      ZFinder zfinder_el(FinalState(eta_e), pt20, PID::ELECTRON, 66*GeV, 116*GeV);
      ZFinder zfinder_mu(FinalState(eta_m), pt20, PID::MUON, 66*GeV, 116*GeV);
      declare(zfinder_el, "zfinder_el");
      declare(zfinder_mu, "zfinder_mu");

      // Define veto FS in order to prevent Z-decay products entering the jet algorithm
      VetoedFinalState had_fs;
      had_fs.addVetoOnThisFinalState(getProjection<ZFinder>("zfinder_el"));
      had_fs.addVetoOnThisFinalState(getProjection<ZFinder>("zfinder_mu"));
      FastJets jets(had_fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::ALL, JetAlg::Invisibles::ALL);
      declare(jets, "jets");

      book(_h_njet_incl              ,  1, 1, _mode+1);
      book(_h_njet_incl_ratio        ,  2, 1, _mode+1, true);
      book(_h_njet_excl              ,  3, 1, _mode+1);
      book(_h_njet_excl_ratio        ,  4, 1, _mode+1, true);
      book(_h_njet_excl_pt150        ,  5, 1, _mode+1);
      book(_h_njet_excl_pt150_ratio  ,  6, 1, _mode+1, true);
      book(_h_njet_excl_vbf          ,  7, 1, _mode+1);
      book(_h_njet_excl_vbf_ratio    ,  8, 1, _mode+1, true);
      book(_h_ptlead                 ,  9, 1, _mode+1);
      book(_h_ptseclead              , 10, 1, _mode+1);
      book(_h_ptthirdlead            , 11, 1, _mode+1);
      book(_h_ptfourthlead           , 12, 1, _mode+1);
      book(_h_ptlead_excl            , 13, 1, _mode+1);
      book(_h_pt_ratio               , 14, 1, _mode+1);
      book(_h_pt_z                   , 15, 1, _mode+1);
      book(_h_pt_z_excl              , 16, 1, _mode+1);
      book(_h_ylead                  , 17, 1, _mode+1);
      book(_h_yseclead               , 18, 1, _mode+1);
      book(_h_ythirdlead             , 19, 1, _mode+1);
      book(_h_yfourthlead            , 20, 1, _mode+1);
      book(_h_deltay                 , 21, 1, _mode+1);
      book(_h_mass                   , 22, 1, _mode+1);
      book(_h_deltaphi               , 23, 1, _mode+1);
      book(_h_deltaR                 , 24, 1, _mode+1);
      book(_h_ptthirdlead_vbf        , 25, 1, _mode+1);
      book(_h_ythirdlead_vbf         , 26, 1, _mode+1);
      book(_h_ht                     , 27, 1, _mode+1);
      book(_h_st                     , 28, 1, _mode+1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      FourMomentum z, lp, lm;
      const ZFinder& zfinder_el = apply<ZFinder>(event, "zfinder_el");
      const ZFinder& zfinder_mu = apply<ZFinder>(event, "zfinder_mu");

      bool e_ok = zfinder_el.constituents().size() == 2 && zfinder_mu.constituents().size() ==0;
      bool m_ok = zfinder_el.constituents().size() == 0 && zfinder_mu.constituents().size() ==2;

      if (_mode == 0 &&  !e_ok && !m_ok ) vetoEvent;
      if (_mode == 1 && !e_ok) vetoEvent;
      if (_mode == 2 && !m_ok) vetoEvent;

      if (zfinder_el.constituents().size() == 2) {
        z = zfinder_el.boson().momentum();
        lp = zfinder_el.constituents()[0].momentum();
        lm = zfinder_el.constituents()[1].momentum();
      }
      else if (zfinder_mu.constituents().size() == 2) {
        z = zfinder_mu.boson().momentum();
        lp = zfinder_mu.constituents()[0].momentum();
        lm = zfinder_mu.constituents()[1].momentum();
      }
      else  vetoEvent;

      if (deltaR(lp, lm) < 0.2) vetoEvent;

      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 30*GeV && Cuts::absrap < 4.4);
      ifilter_discard(jets, deltaRLess(lp, 0.5));
      ifilter_discard(jets, deltaRLess(lm, 0.5));

      // Fill jet multiplicities
      for (size_t ijet = 0; ijet <= jets.size(); ++ijet) {
        _h_njet_incl->fill(ijet);
      }
      _h_njet_excl->fill(jets.size());

      // Require at least one jet
      if (jets.size() >= 1) {
        // Leading jet histos
        const double ptlead   = jets[0].pT()/GeV;
        const double yabslead = fabs(jets[0].rapidity());
        const double ptz   = z.pT()/GeV;
        _h_ptlead->fill(ptlead);
        _h_ylead ->fill(yabslead);
        _h_pt_z  ->fill(ptz);
        // Fill jet multiplicities
        if (ptlead > 150)  _h_njet_excl_pt150->fill(jets.size());

        // Loop over selected jets, fill inclusive distributions
        double st = 0;
        double ht = lp.pT()/GeV + lm.pT()/GeV;
        for (size_t ijet = 0; ijet < jets.size(); ++ijet) {
          ht += jets[ijet].pT()/GeV;
          st += jets[ijet].pT()/GeV;
        }
        _h_ht->fill(ht);
        _h_st->fill(st);

        // Require exactly one jet
        if (jets.size() == 1) {
          _h_ptlead_excl->fill(ptlead);
          _h_pt_z_excl  ->fill(ptz);
        }
      }


      // Require at least two jets
      if (jets.size() >= 2) {
        // Second jet histos
        const double ptlead      = jets[0].pT()/GeV;
        const double pt2ndlead   = jets[1].pT()/GeV;
        const double ptratio     = pt2ndlead/ptlead;
        const double yabs2ndlead = fabs(jets[1].rapidity());
        _h_ptseclead->fill(pt2ndlead);
        _h_yseclead->fill( yabs2ndlead);
        _h_pt_ratio->fill( ptratio);

        // Dijet histos
        const double deltaphi = fabs(deltaPhi(jets[1], jets[0]));
        const double deltarap = fabs(jets[0].rapidity() - jets[1].rapidity()) ;
        const double deltar   = fabs(deltaR(jets[0], jets[1], RAPIDITY));
        const double mass     = (jets[0].momentum() + jets[1].momentum()).mass()/GeV;
        _h_mass->fill(    mass);
        _h_deltay->fill(  deltarap);
        _h_deltaphi->fill(deltaphi);
        _h_deltaR->fill(  deltar);

        if (mass > 350 && deltarap > 3)  _h_njet_excl_vbf->fill(jets.size());
      }

      // Require at least three jets
      if (jets.size() >= 3) {
        // Third jet histos
        const double pt3rdlead   = jets[2].pT()/GeV;
        const double yabs3rdlead = fabs(jets[2].rapidity());
        _h_ptthirdlead->fill(pt3rdlead);
        _h_ythirdlead->fill( yabs3rdlead);

        //Histos after VBF preselection
        const double deltarap = fabs(jets[0].rapidity() - jets[1].rapidity()) ;
        const double mass     = (jets[0].momentum() + jets[1].momentum()).mass();
        if (mass > 350 && deltarap > 3) {
          _h_ptthirdlead_vbf->fill(pt3rdlead);
          _h_ythirdlead_vbf->fill( yabs3rdlead);
        }
      }

      // Require at least four jets
      if (jets.size() >= 4) {
        // Fourth jet histos
        const double pt4thlead   = jets[3].pT()/GeV;
        const double yabs4thlead = fabs(jets[3].rapidity());
        _h_ptfourthlead->fill(pt4thlead);
        _h_yfourthlead->fill( yabs4thlead);
      }
    }

    /// @name Ratio calculator util functions
    //@{

    /// Calculate the efficiency error, being careful about div-by-zero
    double err_incl(const HistoBin1D &M, const HistoBin1D &N, bool hasWeights) {
      double r = safediv(M.sumW(), N.sumW());
      if (hasWeights) { // use F. James's approximation for weighted events
        return sqrt( safediv((1 - 2 * r) * M.sumW2() + r * r * N.sumW2(), N.sumW() * N.sumW()) );
      }
      return sqrt( safediv(r * (1 - r), N.sumW()) );
    }

    /// Calculate the ratio error, being careful about div-by-zero
    double err_excl(const HistoBin1D &A, const HistoBin1D &B) {
      double r = safediv(A.sumW(), B.sumW());
      double dAsquared = safediv(A.sumW2(), A.sumW() * A.sumW()); // squared relative error of A
      double dBsquared = safediv(B.sumW2(), B.sumW() * B.sumW()); // squared relative error of B
      return r * sqrt(dAsquared + dBsquared);
    }

    //@}


    void finalize() {
      bool hasWeights = _h_njet_incl->effNumEntries() != _h_njet_incl->numEntries();
      for (size_t i = 0; i < 6; ++i) {
        _h_njet_incl_ratio->point(i).setY(safediv(_h_njet_incl->bin(i + 1).sumW(), _h_njet_incl->bin(i).sumW()),
                                          err_incl(_h_njet_incl->bin(i + 1), _h_njet_incl->bin(i), hasWeights));
        _h_njet_excl_ratio->point(i).setY(safediv(_h_njet_excl->bin(i + 1).sumW(), _h_njet_excl->bin(i).sumW()),
                                          err_excl(_h_njet_excl->bin(i + 1), _h_njet_excl->bin(i)));
        if (i >= 1) {
          _h_njet_excl_pt150_ratio->point(i - 1).setY(safediv(_h_njet_excl_pt150->bin(i).sumW(), _h_njet_excl_pt150->bin(i - 1).sumW()),
                                                      err_excl(_h_njet_excl_pt150->bin(i), _h_njet_excl_pt150->bin(i - 1)));
          if (i >= 2) {
            _h_njet_excl_vbf_ratio->point(i - 2).setY(safediv(_h_njet_excl_vbf->bin(i).sumW(), _h_njet_excl_vbf->bin(i - 1).sumW()),
                                                      err_excl(_h_njet_excl_vbf->bin(i), _h_njet_excl_vbf->bin(i - 1)));
          }
        }
      }

      double sf = _mode? 1.0 : 0.5;
      const double xs = sf * crossSectionPerEvent()/picobarn;

      scale(_h_njet_incl, xs); scale(_h_njet_excl, xs); scale(_h_njet_excl_pt150, xs); 
      scale(_h_njet_excl_vbf, xs); scale(_h_ptlead, xs); scale(_h_ptseclead, xs); 
      scale(_h_ptthirdlead, xs); scale(_h_ptfourthlead, xs); scale(_h_ptlead_excl, xs);
      scale(_h_pt_ratio, xs); scale(_h_pt_z, xs); scale(_h_pt_z_excl, xs);
      scale(_h_ylead, xs); scale(_h_yseclead, xs); scale(_h_ythirdlead, xs); 
      scale(_h_yfourthlead, xs); scale(_h_deltay, xs); scale(_h_mass, xs); 
      scale(_h_deltaphi, xs); scale(_h_deltaR, xs); scale(_h_ptthirdlead_vbf, xs); 
      scale(_h_ythirdlead_vbf, xs); scale(_h_ht, xs); scale(_h_st, xs);
    }

    //@}


  protected:

    size_t _mode;


  private:

    Scatter2DPtr _h_njet_incl_ratio;
    Scatter2DPtr _h_njet_excl_ratio;
    Scatter2DPtr _h_njet_excl_pt150_ratio;
    Scatter2DPtr _h_njet_excl_vbf_ratio;
    Histo1DPtr _h_njet_incl;
    Histo1DPtr _h_njet_excl;
    Histo1DPtr _h_njet_excl_pt150;
    Histo1DPtr _h_njet_excl_vbf;
    Histo1DPtr _h_ptlead;
    Histo1DPtr _h_ptseclead;
    Histo1DPtr _h_ptthirdlead;
    Histo1DPtr _h_ptfourthlead;
    Histo1DPtr _h_ptlead_excl;
    Histo1DPtr _h_pt_ratio;
    Histo1DPtr _h_pt_z;
    Histo1DPtr _h_pt_z_excl;
    Histo1DPtr _h_ylead;
    Histo1DPtr _h_yseclead;
    Histo1DPtr _h_ythirdlead;
    Histo1DPtr _h_yfourthlead;
    Histo1DPtr _h_deltay;
    Histo1DPtr _h_mass;
    Histo1DPtr _h_deltaphi;
    Histo1DPtr _h_deltaR;
    Histo1DPtr _h_ptthirdlead_vbf;
    Histo1DPtr _h_ythirdlead_vbf;
    Histo1DPtr _h_ht;
    Histo1DPtr _h_st;
  };


  RIVET_DECLARE_PLUGIN(ATLAS_2013_I1230812);

}
