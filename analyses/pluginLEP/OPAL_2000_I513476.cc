// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ParisiTensor.hh"
#include "Rivet/Projections/Hemispheres.hh"

namespace Rivet {


  /// @brief  event shapes at 172, 183, 189
  class OPAL_2000_I513476 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(OPAL_2000_I513476);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      // Projections
      declare(Beam(), "Beams");
      const FinalState fs;
      declare(fs, "FS");
      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      declare(FastJets(fs, FastJets::DURHAM, 0.7), "DurhamJets");
      declare(Sphericity(fs), "Sphericity");
      declare(ParisiTensor(fs), "Parisi");
      const Thrust thrust(fs);
      declare(thrust, "Thrust");
      declare(Hemispheres(thrust), "Hemispheres");

      // Book histograms
      int ih=-1;
      if (fuzzyEquals(sqrtS()/GeV, 172)) {
        ih = 1;
      } else if (fuzzyEquals(sqrtS()/GeV, 183)) {
        ih = 2;
      } else if (fuzzyEquals(sqrtS()/GeV, 189)) {
        ih = 3;
      }
      else {
        throw Error("Invalid CMS energy for OPAL_2000_I513476");
      }
      book(_h_thrust    ,  1,1,ih);
      book(_h_major     ,  2,1,ih);
      book(_h_minor     ,  3,1,ih);
      book(_h_aplanarity,  4,1,ih);
      book(_h_oblateness,  5,1,ih);
      book(_h_C         ,  6,1,ih);
      book(_h_rhoH      ,  7,1,ih);
      book(_h_sphericity,  8,1,ih);
      book(_h_totalB    ,  9,1,ih);
      book(_h_wideB     , 10,1,ih);
      book(_h_y23       , 11,1,ih);
      book(_h_mult      , 13,1,ih);
      book(_h_pTin      , 15,1,ih);
      book(_h_pTout     , 16,1,ih);
      book(_h_y         , 17,1,ih);
      book(_h_x         , 18,1,ih);
      book(_h_xi        , 19,1,ih);
      book(_sumW,"/TMP/sumW");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      if (cfs.size() < 2) vetoEvent;

      _sumW->fill();
      
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      
      // Thrust related
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      _h_thrust    ->fill(thrust.thrust()     );
      _h_major     ->fill(thrust.thrustMajor());
      _h_minor     ->fill(thrust.thrustMinor());
      _h_oblateness->fill(thrust.oblateness() );

      // Sphericity related
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");
      _h_sphericity->fill(sphericity.sphericity());
      _h_aplanarity->fill(sphericity.aplanarity());
      
      // C parameter
      const ParisiTensor& parisi = apply<ParisiTensor>(event, "Parisi");
      _h_C->fill(parisi.C());

      // Hemispheres
      const Hemispheres& hemi = apply<Hemispheres>(event, "Hemispheres");

      _h_rhoH  ->fill(hemi.scaledMhigh());
      _h_wideB ->fill(hemi.Bmax());
      _h_totalB->fill(hemi.Bsum());
      
      // Jets
      const FastJets& durjet = apply<FastJets>(event, "DurhamJets");
      const double y23 = durjet.clusterSeq()->exclusive_ymerge_max(2);
      _h_y23->fill(y23);

      // charged particles
      _h_mult->fill(cfs.particles().size());
      for (const Particle& p : cfs.particles()) {
        const Vector3 mom3  = p.p3();
        const double energy = p.E();
        const double pTinT  = dot(mom3, thrust.thrustMajorAxis());
        const double pToutT = dot(mom3, thrust.thrustMinorAxis());
      	_h_pTin ->fill(fabs(pTinT/GeV) );
      	_h_pTout->fill(fabs(pToutT/GeV));
        const double momT = dot(thrust.thrustAxis(), mom3);
        const double rapidityT = 0.5 * std::log((energy + momT) / (energy - momT));
      	_h_y->fill(fabs(rapidityT));
        const double mom = mom3.mod();
        const double scaledMom = mom/meanBeamMom;
        const double logInvScaledMom = -std::log(scaledMom);
        _h_xi->fill(logInvScaledMom);
        _h_x ->fill(scaledMom      );
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_thrust    ,1./ *_sumW);
      scale(_h_major     ,1./ *_sumW);
      scale(_h_minor     ,1./ *_sumW);
      scale(_h_aplanarity,1./ *_sumW);
      scale(_h_oblateness,1./ *_sumW);
      scale(_h_C         ,1./ *_sumW);
      scale(_h_rhoH      ,1./ *_sumW);
      scale(_h_sphericity,1./ *_sumW);
      scale(_h_totalB    ,1./ *_sumW);
      scale(_h_wideB     ,1./ *_sumW);
      scale(_h_y23       ,1./ *_sumW);
      scale(_h_mult      ,200./ *_sumW);
      scale(_h_pTin      ,1./ *_sumW);
      scale(_h_pTout     ,1./ *_sumW);
      scale(_h_y         ,1./ *_sumW);
      scale(_h_x         ,1./ *_sumW);
      scale(_h_xi        ,1./ *_sumW);
      // mean multiplicity
      double nch     = _h_mult->xMean();
      double nch_err = _h_mult->xStdErr();
      Scatter2DPtr m_ch;
      book(m_ch,14,1,1);
      m_ch->addPoint(sqrtS(),nch,0.,nch_err);
      // mean ptIn
      double pTin     = _h_pTin->xMean();
      double pTin_err = _h_pTin->xStdErr();
      Scatter2DPtr m_pTin;
      book(m_pTin,20,1,1);
      m_pTin->addPoint(sqrtS(),pTin,0.,pTin_err);
      // mean ptOut
      double pTout     = _h_pTout->xMean();
      double pTout_err = _h_pTout->xStdErr();
      Scatter2DPtr m_pTout;
      book(m_pTout,20,1,2);
      m_pTout->addPoint(sqrtS(),pTout,0.,pTout_err);
      // mean y
      double y     = _h_y->xMean();
      double y_err = _h_y->xStdErr();
      Scatter2DPtr m_y;
      book(m_y,20,1,3);
      m_y->addPoint(sqrtS(),y,0.,y_err);
      // mean x
      double x     = _h_x->xMean();
      double x_err = _h_x->xStdErr();
      Scatter2DPtr m_x;
      book(m_x,20,1,4);
      m_x->addPoint(sqrtS(),x,0.,x_err);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_thrust,_h_major,_h_minor,_h_aplanarity,_h_oblateness,_h_C,_h_rhoH,_h_sphericity;
    Histo1DPtr _h_totalB,_h_wideB,_h_y23,_h_mult,_h_pTin,_h_pTout,_h_y,_h_x,_h_xi;
    CounterPtr _sumW;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(OPAL_2000_I513476);


}
