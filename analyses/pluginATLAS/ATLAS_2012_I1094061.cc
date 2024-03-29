// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// Two-particle correlation functions in pp collisions at 900 GeV and 7 TeV
  class ATLAS_2012_I1094061 : public Analysis {


    /// Container for a pair of foreground and background histos, divided at the end of the analysis
    struct HistoPair{
      enum HistoType { FOREGROUND, BACKGROUND };

      HistoPair() : _ana(nullptr)
      {  }

      void init(int ds, int xaxis, int yaxis, ATLAS_2012_I1094061* ana) {
        _ana = ana;
        const string hcode = ana->makeAxisCode(ds, xaxis, yaxis);
        _h_foreground = ana->bookHisto1D("tmpForeground_" + hcode, ana->refData(ds, xaxis, yaxis));
        _h_background = ana->bookHisto1D("tmpBackground_" + hcode, ana->refData(ds, xaxis, yaxis));
        _s_final = ana->bookScatter2D(ds, xaxis, yaxis, true);
      }

      void fillForeground(double value, double weight) {
        _h_foreground->fill(value, weight);
        _h_foreground->fill(-value, weight);
      }

      void fillBackground(double value, double weight) {
        _h_background->fill(value, weight);
        _h_background->fill(-value, weight);
      }

      void fill(double value, double weight, HistoType type) {
        if (type == FOREGROUND) {
          fillForeground(value, weight);
        } else { // type == BACKGROUND
          fillBackground(value, weight);
        }
      }

      void finalize(double wgtSum, double bgWeight, double avNTracks) {
        _h_foreground->scaleW(1/wgtSum);
        _h_background->scaleW(1/bgWeight);
        _ana->divide(_h_foreground, _h_background, _s_final);
        for (Point2D& p : _s_final->points()) {
          p.setY(p.y() - (avNTracks-1));
        }
      }

    private:

      ATLAS_2012_I1094061 *_ana;

      Histo1DPtr _h_foreground;
      Histo1DPtr _h_background;
      Scatter2DPtr _s_final;

    };


  public:

    /// Constructor
    ATLAS_2012_I1094061()
      : Analysis("ATLAS_2012_I1094061"),
        _minpT(100.*MeV), _nVersions(5), _version(0),
        _etaCut(2.), _phiCut(0.5*M_PI),
        _historyInclusive(_nVersions, ParticleVector()), _historyN20(_nVersions, ParticleVector()),
        _historyInclusiveWgts(_nVersions, 0.), _historyN20Wgts(_nVersions, 0.),
        _particleCountInclusive(0.), _particleCountN20(0.),
        _weightInclusive(0.), _weightN20(0.),
        _bgWeightInclusive(0.), _bgWeightN20(0.)
    {    }


    /// @name Analysis methods
    //@{

    void init() {

      const ChargedFinalState cfs(-2.5, 2.5, _minpT);
      declare(cfs, "ChargedParticles");

      // Only do the multiplicity > 20 plots for 7 TeV collisions
      _doN20 = (fabs(sqrtS() - 7000.*GeV) < 0.1*GeV);

      int yaxis = (_doN20) ? 2: 1;
      _hp_DEta_0_pi.init(1, 1, yaxis, this);
      _hp_DEta_0_pi2.init(2, 1, yaxis, this);
      _hp_DEta_pi2_pi.init(3, 1, yaxis, this);
      _hp_DPhi_0_2.init(4, 1, yaxis, this);
      _hp_DPhi_2_5.init(5, 1, yaxis, this);

      if (_doN20) {
        yaxis = 3;
        _hp_N20_DEta_0_pi.init(1, 1, yaxis, this);
        _hp_N20_DEta_0_pi2.init(2, 1, yaxis, this);
        _hp_N20_DEta_pi2_pi.init(3, 1, yaxis, this);
        _hp_N20_DPhi_0_2.init(4, 1, yaxis, this);
        _hp_N20_DPhi_2_5.init(5, 1, yaxis, this);
      }
    }


    void analyze(const Event& evt) {

      const ChargedFinalState& cfsProj = apply<ChargedFinalState>(evt, "ChargedParticles");

      Particles chargedParticles = cfsProj.particles();
      if (chargedParticles.size() < 2) vetoEvent;
      const bool hasN20 = (_doN20 && chargedParticles.size() >= 20);

      const double dMultiplicity = (double) chargedParticles.size();
      const double multiplicityWeightIncr = dMultiplicity * evt.weight();

      _weightInclusive += evt.weight();
      _particleCountInclusive += multiplicityWeightIncr;
      if (hasN20) {
        _weightN20 += evt.weight();
        _particleCountN20 += multiplicityWeightIncr;
      }

      double fgWeight = 2.*evt.weight() / dMultiplicity;

      for (Particles::const_iterator p1 = chargedParticles.begin(); p1 != chargedParticles.end(); ++p1) {
        Particles::const_iterator p2 = p1;
        ++p2;

        // Fill the foreground distributions
        while (p2 != chargedParticles.end()) {
          fillHistosInclusive(*p1, *p2, fgWeight, HistoPair::FOREGROUND);
          if (hasN20)
            fillHistosN20(*p1, *p2, fgWeight, HistoPair::FOREGROUND);
          ++p2;
        }

        // Loop over the history of particles from previous events and fill the background
        // by correlating those particles with the current event
        for (size_t version = 0; version != _nVersions; ++version) {

          const Particles& bgParticles = _historyInclusive[version];
          double bgWeight = evt.weight() * _historyInclusiveWgts[version];
          for (Particles::const_iterator p2 = bgParticles.begin(); p2 != bgParticles.end(); ++p2) {
            fillHistosInclusive(*p1, *p2, bgWeight, HistoPair::BACKGROUND);
            _bgWeightInclusive += bgWeight;
          }

          if (!hasN20) continue;

          const Particles& bgParticlesN20 = _historyN20[version];
          bgWeight = evt.weight() * _historyN20Wgts[version];

          for (Particles::const_iterator p2 = bgParticlesN20.begin(); p2 != bgParticlesN20.end(); ++p2) {
            fillHistosN20(*p1, *p2, bgWeight, HistoPair::BACKGROUND);
            _bgWeightN20 += bgWeight;
          }

        }

      }

      // Overwrite the history for the version count number
      _historyInclusive[_version] = chargedParticles;
      _historyInclusiveWgts[_version] = evt.weight();
      if (hasN20) {
        _historyN20[_version] = chargedParticles;
        _historyN20Wgts[_version] = evt.weight();
      }
      ++_version;
      if (_version == _nVersions) _version = 0;

    }


    void finalize() {

      const double avMultiplicity = _particleCountInclusive / _weightInclusive;
      _hp_DEta_0_pi.finalize(_weightInclusive,  _bgWeightInclusive, avMultiplicity);
      _hp_DEta_0_pi2.finalize(_weightInclusive, _bgWeightInclusive, avMultiplicity);
      _hp_DEta_pi2_pi.finalize(_weightInclusive,_bgWeightInclusive, avMultiplicity);
      _hp_DPhi_0_2.finalize(_weightInclusive, _bgWeightInclusive, avMultiplicity);
      _hp_DPhi_2_5.finalize(_weightInclusive, _bgWeightInclusive, avMultiplicity);

      if (_doN20) {
        const double avMultiplicityN20 = _particleCountN20 / _weightN20;
        _hp_N20_DEta_0_pi.finalize(_weightN20,   _bgWeightN20, avMultiplicityN20);
        _hp_N20_DEta_0_pi2.finalize(_weightN20,  _bgWeightN20, avMultiplicityN20);
        _hp_N20_DEta_pi2_pi.finalize(_weightN20, _bgWeightN20, avMultiplicityN20);
        _hp_N20_DPhi_0_2.finalize(_weightN20, _bgWeightN20, avMultiplicityN20);
        _hp_N20_DPhi_2_5.finalize(_weightN20, _bgWeightN20, avMultiplicityN20);
      }

    }

    //@}


    void fillHistos(const Particle &p1, const Particle &p2, double weight, HistoPair::HistoType type, bool inclusive) {

      const double dEta = fabs(p1.eta() - p2.eta());
      const double dPhi = mapAngle0ToPi(p1.phi() - p2.phi());
      const double dPhiShift = TWOPI - dPhi;

      HistoPair& dEta_0_pi   = (inclusive)? _hp_DEta_0_pi   :_hp_N20_DEta_0_pi;
      HistoPair& dPhi_0_2    = (inclusive)? _hp_DPhi_0_2    :_hp_N20_DPhi_0_2;
      HistoPair& dPhi_2_5    = (inclusive)? _hp_DPhi_2_5    :_hp_N20_DPhi_2_5;
      HistoPair& dEta_0_pi2  = (inclusive)? _hp_DEta_0_pi2  :_hp_N20_DEta_0_pi2;
      HistoPair& dEta_pi2_pi = (inclusive)? _hp_DEta_pi2_pi :_hp_N20_DEta_pi2_pi;

      dEta_0_pi.fill(dEta, weight, type);

      if (dEta < _etaCut) {
        dPhi_0_2.fill(dPhi, weight, type);
        dPhi_0_2.fill(dPhiShift, weight, type);
      } else {
        dPhi_2_5.fill(dPhi, weight, type);
        dPhi_2_5.fill(dPhiShift, weight, type);
      }

      if (dPhi < _phiCut) {
        dEta_0_pi2.fill(dEta, weight, type);
      } else {
        dEta_pi2_pi.fill(dEta, weight, type);
      }

    }


    void fillHistosInclusive(const Particle &p1, const Particle &p2, double weight, HistoPair::HistoType type) {
      fillHistos(p1, p2, weight, type, true);
    }

    void fillHistosN20(const Particle &p1, const Particle &p2, double weight, HistoPair::HistoType type) {
      fillHistos(p1, p2, weight, type, false);
    }


  private:

    /// Cut values
    double _minpT;

    /// History versions
    size_t _nVersions, _version;

    /// Cut values
    double _etaCut, _phiCut;

    /// Vectors of particles from _nVersions previous events, to construct the background correlation.
    vector<Particles> _historyInclusive, _historyN20;

    /// History-event weights
    vector<double> _historyInclusiveWgts, _historyN20Wgts;

    double _particleCountInclusive, _particleCountN20;
    double _weightInclusive, _weightN20;
    double _bgWeightInclusive, _bgWeightN20;

    bool _doN20;

    HistoPair _hp_DEta_0_pi, _hp_DEta_0_pi2, _hp_DEta_pi2_pi;
    HistoPair _hp_DPhi_0_2, _hp_DPhi_2_5;
    HistoPair _hp_N20_DEta_0_pi, _hp_N20_DEta_0_pi2, _hp_N20_DEta_pi2_pi;
    HistoPair _hp_N20_DPhi_0_2, _hp_N20_DPhi_2_5;

  };


  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1094061);

}
