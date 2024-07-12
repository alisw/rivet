// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/TauFinder.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {


  /// @brief MC validation analysis for tau decays involving photons
  class MC_TAUS_PHOTONS : public Analysis {
    public:

      RIVET_DEFAULT_ANALYSIS_CTOR(MC_TAUS_PHOTONS);

      /// Book projections and histograms
      void init() {
        TauFinder taus(TauFinder::DecayMode::ANY, Cuts::pT > 500*MeV);
        declare(taus, "Taus");

        IdentifiedFinalState photons(Cuts::pT > 500*MeV, PID::PHOTON);
        declare(photons, "Photons");

        book(_nPhotonsEl, "nPhotonsEl", 20, 0, 20);
        book(_nPhotonsMu, "nPhotonsMu", 20, 0, 20);
        book(_nPhotonsHad, "nPhotonsHad", 20, 0, 20);

        book(_tauMassEl, "tauMassEl", 50, 0, 5);
        book(_tauMassMu, "tauMassMu", 50, 0, 5);
        book(_tauMassHad, "tauMassHad", 50, 0, 5);

        book(_pFracPhotonsEl, "pFracPhotonsEl", 50, 0, 2);
        book(_pFracPhotonsMu, "pFracPhotonsMu", 50, 0, 2);
        book(_pFracPhotonsHad, "pFracPhotonsHad", 50, 0, 2);

        book(_logPFracPhotonsEl, "logPFracPhotonsEl", 50, -5, 0);
        book(_logPFracPhotonsMu, "logPFracPhotonsMu", 50, -5, 0);
        book(_logPFracPhotonsHad, "logPFracPhotonsHad", 50, -5, 0);

        book(_logPFracNotPhotonsEl, "logPFracNotPhotonsEl", 50, -5, 0);
        book(_logPFracNotPhotonsMu, "logPFracNotPhotonsMu", 50, -5, 0);
        book(_logPFracNotPhotonsHad, "logPFracNotPhotonsHad", 50, -5, 0);

        book(_restFramePhotonsEnergyEl, "RestFramePhotonsEnergyEl", 50, 0, 3);
        book(_restFramePhotonsEnergyMu, "RestFramePhotonsEnergyMu", 50, 0, 3);
        book(_restFramePhotonsEnergyHad, "RestFramePhotonsEnergyHad", 50, 0, 3);

        return;
      }


      /// Per-event analysis
      void analyze(const Event& event) {
        const Particles taus = apply<TauFinder>(event, "Taus").particlesByPt();
        const Particles photons = apply<IdentifiedFinalState>(event, "Photons").particlesByPt();

        for (const auto& tau : taus) {
          Particles descendants = tau.stableDescendants();
          bool hasHad = false, hasEl = false, hasMu = false;

          for (const auto& descendant : descendants) {
            if (isHadron(descendant)) {
              hasHad = true;
            } else if (abs(descendant.pid()) == PID::ELECTRON) {
              hasEl = true;
            } else if (abs(descendant.pid()) == PID::MUON) {
              hasMu = true;
            }

            continue;
          }

          // decaymode: hadronic = 0, electron = 1, muon = 2
          int decaymode = -1;
          if (hasHad) decaymode = 0;
          else if (hasEl) decaymode = 1;
          else if (hasMu) decaymode = 2;
          assert(decaymode >= 0);

          FourMomentum taumom;
          FourMomentum photonmom;
          int nphotons = 0;
          for (const auto& descendant : descendants) {
            if (descendant.pid() == PID::PHOTON) {
              photonmom += descendant.mom();
              nphotons++;
            }

            taumom += descendant.mom();

            continue;
          }

          LorentzTransform boost =
            LorentzTransform::mkFrameTransformFromBeta(taumom.betaVec());

          float taumass = taumom.mass();

          float photonfrac = photonmom.p() / tau.p();
          if (photonfrac >= 1.0)
            photonfrac = 1.0 - 1e-6;

          if (decaymode == 0) {
            _nPhotonsHad->fill(nphotons);
            _tauMassHad->fill(taumass / GeV);

            _pFracPhotonsHad->fill(photonfrac);
            _logPFracPhotonsHad->fill(log(photonfrac));
            _logPFracNotPhotonsHad->fill(log(1.0 - photonfrac));
            _restFramePhotonsEnergyHad->fill(boost(photonmom).E() / GeV);

          } else if (decaymode == 1) {
            _nPhotonsEl->fill(nphotons);
            _tauMassEl->fill(taumass / GeV);

            _pFracPhotonsEl->fill(photonfrac);
            _logPFracPhotonsEl->fill(log(photonfrac));
            _logPFracNotPhotonsEl->fill(log(1.0 - photonfrac));
            _restFramePhotonsEnergyEl->fill(boost(photonmom).E() / GeV);

          } else if (decaymode == 2) {
            _nPhotonsMu->fill(nphotons);
            _tauMassMu->fill(taumass / GeV);

            _pFracPhotonsMu->fill(photonfrac);
            _logPFracPhotonsMu->fill(log(photonfrac));
            _logPFracNotPhotonsMu->fill(log(1.0 - photonfrac));
            _restFramePhotonsEnergyMu->fill(boost(photonmom).E() / GeV);

          }

          continue;
        }

        return;
      }


      /// Normalisations etc.
      void finalize() {
        _nPhotonsEl->normalize();
        _nPhotonsMu->normalize();
        _nPhotonsHad->normalize();

        _tauMassEl->normalize();
        _tauMassMu->normalize();
        _tauMassHad->normalize();

        _pFracPhotonsEl->normalize();
        _pFracPhotonsMu->normalize();
        _pFracPhotonsHad->normalize();

        _logPFracPhotonsEl->normalize();
        _logPFracPhotonsMu->normalize();
        _logPFracPhotonsHad->normalize();

        _logPFracNotPhotonsEl->normalize();
        _logPFracNotPhotonsMu->normalize();
        _logPFracNotPhotonsHad->normalize();

        _restFramePhotonsEnergyEl->normalize();
        _restFramePhotonsEnergyMu->normalize();
        _restFramePhotonsEnergyHad->normalize();

        return;
      }

    private:

      Histo1DPtr
          _nPhotonsEl
        , _nPhotonsMu
        , _nPhotonsHad
        , _tauMassEl
        , _tauMassMu
        , _tauMassHad
        , _pFracPhotonsEl
        , _pFracPhotonsMu
        , _pFracPhotonsHad
        , _logPFracPhotonsEl
        , _logPFracPhotonsMu
        , _logPFracPhotonsHad
        , _logPFracNotPhotonsEl
        , _logPFracNotPhotonsMu
        , _logPFracNotPhotonsHad
        , _restFramePhotonsEnergyEl
        , _restFramePhotonsEnergyMu
        , _restFramePhotonsEnergyHad
        ;

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MC_TAUS_PHOTONS);

} 
