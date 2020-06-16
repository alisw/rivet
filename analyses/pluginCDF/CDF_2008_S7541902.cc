// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include <algorithm>

namespace Rivet {


  /// @brief CDF jet pT and multiplicity distributions in W + jets events
  ///
  /// This CDF analysis provides jet pT distributions for 4 jet multiplicity bins
  /// as well as the jet multiplicity distribution in W + jets events.
  /// e-Print: arXiv:0711.4044 [hep-ex]
  class CDF_2008_S7541902 : public Analysis {
  public:

    /// Constructor
    CDF_2008_S7541902()
      : Analysis("CDF_2008_S7541902"),
        _electronETCut(20.0*GeV), _electronETACut(1.1),
        _eTmissCut(30.0*GeV), _mTCut(20.0*GeV),
        _jetEtCutA(20.0*GeV),  _jetEtCutB(25.0*GeV), _jetETA(2.0)
    {    }


    /// @name Analysis methods
    //@{

    void init() {
      // Set up projections
      // Basic FS
      FinalState fs((Cuts::etaIn(-3.6, 3.6)));
      declare(fs, "FS");

      // Create a final state with any e-nu pair with invariant mass 65 -> 95 GeV and ET > 20 (W decay products)
      vector<pair<PdgId,PdgId> > vids;
      vids += make_pair(PID::ELECTRON, PID::NU_EBAR);
      vids += make_pair(PID::POSITRON, PID::NU_E);
      FinalState fs2((Cuts::etaIn(-3.6, 3.6) && Cuts::pT >=  20*GeV));
      InvMassFinalState invfs(fs2, vids, 65*GeV, 95*GeV);
      declare(invfs, "INVFS");

      // Make a final state without the W decay products for jet clustering
      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(invfs);
      declare(vfs, "VFS");
      declare(FastJets(vfs, FastJets::CDFJETCLU, 0.4), "Jets");

      // Book histograms
      for (int i = 0 ; i < 4 ; ++i) {
        book(_histJetEt[i] ,1+i, 1, 1);
        book(_histJetMultRatio[i], 5, 1, i+1, true);
        /// @todo These would be better off as YODA::Counter until finalize()
        book(_histJetMult[i] ,6+i, 1, 1); // _sumW is essentially the 0th "histo" counter
      }

      book(_sumW,"sumW");
    }


    /// Do the analysis
    void analyze(const Event& event) {
      // Get the W decay products (electron and neutrino)
      const InvMassFinalState& invMassFinalState = apply<InvMassFinalState>(event, "INVFS");
      const Particles&  wDecayProducts = invMassFinalState.particles();

      FourMomentum electronP, neutrinoP;
      bool gotElectron(false), gotNeutrino(false);
      for (const Particle& p : wDecayProducts) {
        FourMomentum p4 = p.momentum();
        if (p4.Et() > _electronETCut && fabs(p4.eta()) < _electronETACut && p.abspid() == PID::ELECTRON) {
          electronP = p4;
          gotElectron = true;
        }
        else if (p4.Et() > _eTmissCut && p.abspid() == PID::NU_E) {
          neutrinoP = p4;
          gotNeutrino = true;
        }
      }

      // Veto event if the electron or MET cuts fail
      if (!gotElectron || !gotNeutrino) vetoEvent;

      // Veto event if the MTR cut fails
      double mT2 = 2.0 * ( electronP.pT()*neutrinoP.pT() - electronP.px()*neutrinoP.px() - electronP.py()*neutrinoP.py() );
      if (sqrt(mT2) < _mTCut ) vetoEvent;

      // Get the jets
      const JetAlg& jetProj = apply<FastJets>(event, "Jets");
      Jets theJets = jetProj.jets(cmpMomByEt, Cuts::Et > _jetEtCutA);
      size_t njetsA(0), njetsB(0);
      for (const Jet& j : theJets) {
        const FourMomentum pj = j.momentum();
        if (fabs(pj.rapidity()) < _jetETA) {
          // Fill differential histograms for top 4 jets with Et > 20
          if (njetsA < 4 && pj.Et() > _jetEtCutA) {
            ++njetsA;
            _histJetEt[njetsA-1]->fill(pj.Et());
          }
          // Count number of jets with Et > 25 (for multiplicity histograms)
          if (pj.Et() > _jetEtCutB) ++njetsB;
        }
      }

      // Increment event counter
      _sumW->fill();

      // Jet multiplicity
      for (size_t i = 1; i <= njetsB; ++i) {
        /// @todo This isn't really a histogram: replace with a YODA::Counter when we have one!
        _histJetMult[i-1]->fill(1960.);
        if (i == 4) break;
      }
    }



    /// Finalize
    void finalize() {

      // Fill the 0th ratio histogram specially
      /// @todo This special case for 1-to-0 will disappear if we use Counters for all mults including 0.
      if (_sumW->val() > 0) {
        const YODA::Histo1D::Bin& b0 = _histJetMult[0]->bin(0);
        double ratio = b0.area()/dbl(*_sumW);
        double frac_err = 1/dbl(*_sumW); ///< This 1/sqrt{N} error treatment isn't right for weighted events: use YODA::Counter
        if (b0.area() > 0) frac_err = sqrt( sqr(frac_err) + sqr(b0.areaErr()/b0.area()) );
        _histJetMultRatio[0]->point(0).setY(ratio, ratio*frac_err);
      }

      // Loop over the non-zero multiplicities
      for (size_t i = 0; i < 3; ++i) {
        const YODA::Histo1D::Bin& b1 = _histJetMult[i]->bin(0);
        const YODA::Histo1D::Bin& b2 = _histJetMult[i+1]->bin(0);
        if (b1.area() == 0.0) continue;
        double ratio = b2.area()/b1.area();
        double frac_err = b1.areaErr()/b1.area();
        if (b2.area() > 0) frac_err = sqrt( sqr(frac_err) + sqr(b2.areaErr()/b2.area()) );
        _histJetMultRatio[i+1]->point(0).setY(ratio, ratio*frac_err);
      }

      // Normalize the non-ratio histograms
      for (size_t i = 0; i < 4; ++i) {
        scale(_histJetEt[i], crossSection()/picobarn/sumOfWeights());
        scale(_histJetMult[i], crossSection()/picobarn/sumOfWeights());
      }

    }

    //@}


  private:

    /// @name Cuts
    //@{

    /// Cut on the electron ET:
    double _electronETCut;
    /// Cut on the electron ETA:
    double _electronETACut;
    /// Cut on the missing ET
    double _eTmissCut;
    /// Cut on the transverse mass squared
    double _mTCut;
    /// Cut on the jet ET for differential cross sections
    double _jetEtCutA;
    /// Cut on the jet ET for jet multiplicity
    double _jetEtCutB;
    /// Cut on the jet ETA
    double _jetETA;

    //@}

    /// @name Histograms
    //@{
    Histo1DPtr _histJetEt[4];
    Histo1DPtr _histJetMultNorm;
    Scatter2DPtr _histJetMultRatio[4];
    Histo1DPtr _histJetMult[4];
    CounterPtr _sumW;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CDF_2008_S7541902);

}
