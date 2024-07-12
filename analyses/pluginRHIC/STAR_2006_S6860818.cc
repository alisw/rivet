// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// STAR strange particle spectra in pp at 200 GeV
  class STAR_2006_S6860818 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(STAR_2006_S6860818);

    /// Book projections and histograms
    void init() {
      ChargedFinalState bbc1(Cuts::etaIn(-5.0, -3.5)); // beam-beam-counter trigger
      ChargedFinalState bbc2(Cuts::etaIn( 3.5,  5.0)); // beam-beam-counter trigger
      declare(bbc1, "BBC1");
      declare(bbc2, "BBC2");

      UnstableParticles ufs(Cuts::abseta < 2.5);
      declare(ufs, "UFS");

      book(_h_pT_k0s        ,1, 1, 1);
      book(_h_pT_kminus     ,1, 2, 1);
      book(_h_pT_kplus      ,1, 3, 1);
      book(_h_pT_lambda     ,1, 4, 1);
      book(_h_pT_lambdabar  ,1, 5, 1);
      book(_h_pT_ximinus    ,1, 6, 1);
      book(_h_pT_xiplus     ,1, 7, 1);
      //book(_h_pT_omega      ,1, 8, 1);
      book(_h_antibaryon_baryon_ratio, 2, 1, 1);
      book(_h_lambar_lam, 2, 2, 1);
      book(_h_xiplus_ximinus, 2, 3, 1);
      book(_h_pT_vs_mass    ,3, 1, 1);

      for (size_t i = 0; i < 4; i++) {
        book(_nWeightedBaryon[i], "TMP/nWeightedBaryon"+to_str(i));
        book(_nWeightedAntiBaryon[i], "TMP/nWeightedAntiBaryon"+to_str(i));
      }
      book(_sumWeightSelected, "sumWselected");
    }


    /// Do the analysis
    void analyze(const Event& event) {
      const ChargedFinalState& bbc1 = apply<ChargedFinalState>(event, "BBC1");
      const ChargedFinalState& bbc2 = apply<ChargedFinalState>(event, "BBC2");
      if (bbc1.size() < 1 || bbc2.size() < 1) {
        MSG_DEBUG("Failed beam-beam-counter trigger");
        vetoEvent;
      }
      _sumWeightSelected->fill();

      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
        if (p.absrap() < 0.5) {
          const PdgId pid = p.pid();
          const double pT = p.pT() / GeV;
          switch (abs(pid)) {
          case PID::PIPLUS:
            if (pid < 0) _h_pT_vs_mass->fill(0.1396, pT);
            break;
          case PID::PROTON:
            if (pid < 0) _h_pT_vs_mass->fill(0.9383, pT);
            if (pT > 0.4) {
              pid > 0 ? _nWeightedBaryon[0]->fill() : _nWeightedAntiBaryon[0]->fill();
            }
            break;
          case PID::K0S:
            if (pT > 0.2) {
              _h_pT_k0s->fill(pT, 1.0/pT);
            }
            _h_pT_vs_mass->fill(0.5056, pT);
            break;
          case PID::K0L:
            _h_pT_vs_mass->fill(0.5056, pT);
            break;
          case 113: // rho0(770)
            _h_pT_vs_mass->fill(0.7755, pT);
            break;
          case 313: // K0*(892)
            _h_pT_vs_mass->fill(0.8960, pT);
            break;
          case 333: // phi(1020)
            _h_pT_vs_mass->fill(1.0190, pT);
            break;
          case 3214: // Sigma(1385)
            _h_pT_vs_mass->fill(1.3840, pT);
            break;
          case 3124: // Lambda(1520)
            _h_pT_vs_mass->fill(1.5200, pT);
            break;
          case PID::KPLUS:
            if (pid < 0) _h_pT_vs_mass->fill(0.4856, pT);
            if (pT > 0.2) {
              pid > 0 ? _h_pT_kplus->fill(pT, 1.0/pT) : _h_pT_kminus->fill(pT, 1.0/pT);
            }
            break;
          case PID::LAMBDA:
            pid > 0 ? _h_pT_vs_mass->fill(1.1050, pT) : _h_pT_vs_mass->fill(1.1250, pT);
            if (pT > 0.3) {
              pid > 0 ? _h_pT_lambda->fill(pT, 1.0/pT) : _h_pT_lambdabar->fill(pT, 1.0/pT);
              pid > 0 ? _nWeightedBaryon[1]->fill() : _nWeightedAntiBaryon[1]->fill();
            }
            break;
          case PID::XIMINUS:
            pid > 0 ? _h_pT_vs_mass->fill(1.3120, pT) : _h_pT_vs_mass->fill(1.3320, pT);
            if (pT > 0.5) {
              pid > 0 ? _h_pT_ximinus->fill(pT, 1.0/pT) : _h_pT_xiplus->fill(pT, 1.0/pT);
              pid > 0 ? _nWeightedBaryon[2]->fill() : _nWeightedAntiBaryon[2]->fill();
            }
            break;
          case PID::OMEGAMINUS:
            _h_pT_vs_mass->fill(1.6720, pT);
            if (pT > 0.5) {
              //_h_pT_omega->fill(pT, 1.0/pT);
              pid > 0 ? _nWeightedBaryon[3]->fill() : _nWeightedAntiBaryon[3]->fill();
            }
            break;
          }

        }
      }
    }


    /// Finalize
    void finalize() {
      std::vector<Point2D> points;
      for (size_t i=0 ; i<4 ; i++) {
        if (_nWeightedBaryon[i]->val()==0 || _nWeightedAntiBaryon[i]->val()==0) {
          points.push_back(Point2D(i,0,0.5,0));
        } else {
          double y  = safediv(_nWeightedAntiBaryon[i]->val(), _nWeightedBaryon[i]->val(), 0.);
          double dy = sqrt( safediv(1., _nWeightedAntiBaryon[i]->numEntries(), 0.) + safediv(1., _nWeightedBaryon[i]->numEntries(), 0.) );
          points.push_back(Point2D(i,y,0.5,y*dy));
        }
      }
      _h_antibaryon_baryon_ratio->addPoints( points );

      divide(_h_pT_lambdabar,_h_pT_lambda, _h_lambar_lam);
      divide(_h_pT_xiplus,_h_pT_ximinus, _h_xiplus_ximinus);

      const YODA::Scatter1D factor = (1./(2.0 * M_PI)) / *_sumWeightSelected;
      scale(_h_pT_k0s,       factor);
      scale(_h_pT_kminus,    factor);
      scale(_h_pT_kplus,     factor);
      scale(_h_pT_lambda,    factor);
      scale(_h_pT_lambdabar, factor);
      scale(_h_pT_ximinus,   factor);
      scale(_h_pT_xiplus,    factor);
      //scale(_h_pT_omega,     1./(2*M_PI*_sumWeightSelected));
      MSG_DEBUG("sumOfWeights()     = " << sumOfWeights());
      MSG_DEBUG("_sumWeightSelected = " << _sumWeightSelected->val());
    }

  private:

    CounterPtr _sumWeightSelected;
    array<CounterPtr, 4> _nWeightedBaryon;
    array<CounterPtr, 4> _nWeightedAntiBaryon;

    Histo1DPtr _h_pT_k0s, _h_pT_kminus, _h_pT_kplus, _h_pT_lambda, _h_pT_lambdabar, _h_pT_ximinus, _h_pT_xiplus;
    //Histo1DPtr _h_pT_omega;
    Scatter2DPtr _h_antibaryon_baryon_ratio;
    Profile1DPtr  _h_pT_vs_mass;
    Scatter2DPtr _h_lambar_lam;
    Scatter2DPtr _h_xiplus_ximinus;

  };


  RIVET_DECLARE_ALIASED_PLUGIN(STAR_2006_S6860818, STAR_2006_I722757);

}
