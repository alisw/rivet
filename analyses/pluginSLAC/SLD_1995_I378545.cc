// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include "Rivet/Projections/ParisiTensor.hh"
//#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Event shapes at MZ
  class SLD_1995_I378545 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SLD_1995_I378545);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      const FinalState fs;
      declare(fs, "FS");
      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      Sphericity sphere(fs);
      declare(sphere, "Sphericity");
      declare(Thrust    (fs), "Thrust"    );
      declare(Hemispheres(sphere), "Hemispheres");
      declare(ParisiTensor(fs), "Parisi");
      // declare(FastJets(fs, FastJets::DURHAM, 0.7), "DurhamJets");
      // declare(FastJets(fs, FastJets::JADE  , 0.7), "JadeJets"  );
      // histograms
      book(_histThrust     , 2, 1, 1);
      book(_histMJetHigh   , 3, 1, 1);
      book(_histBT         , 4, 1, 1);
      book(_histBW         , 5, 1, 1);
      book(_histOblateness , 6, 1, 1);
      book(_histC          , 7, 1, 1);
      // book(_histy23JADE    , 8, 1, 1);
      // book(_histy23Durham  ,12, 1, 1);
      book(_histEEC        ,14, 1, 1);
      book(_histAEEC       ,15, 1, 1);
      book(_histJCEF       ,16, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      if ( fs.particles().size() < 2) vetoEvent;
      // thrust
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      _histThrust    ->fill(1.-thrust.thrust());
      _histOblateness->fill(thrust.oblateness() );
      // visible energy
      double Evis = 0.0;
      for (const Particle& p : fs.particles()) {
        Evis += p.E();
      }
      double Evis2 = sqr(Evis);
      // (A)EEC
      // Need iterators since second loop starts at current outer loop iterator, i.e. no "foreach" here!
      for (Particles::const_iterator p_i = fs.particles().begin(); p_i != fs.particles().end(); ++p_i) {
        for (Particles::const_iterator p_j = p_i; p_j != fs.particles().end(); ++p_j) {
          const Vector3 mom3_i = p_i->momentum().p3();
          const Vector3 mom3_j = p_j->momentum().p3();
          const double energy_i = p_i->momentum().E();
          const double energy_j = p_j->momentum().E();
          const double thetaij = 180.*mom3_i.unit().angle(mom3_j.unit())/M_PI;
          double eec = (energy_i*energy_j) / Evis2;
	  if(p_i != p_j)  eec *= 2.;
          _histEEC->fill(           thetaij, eec);
          if (thetaij <90.){
	    _histAEEC->fill(thetaij, -eec);
	  }
          else {
	    _histAEEC->fill(180.-thetaij, eec);
	  }
        }
      }
      // hemisphere related
      const Hemispheres& hemi = apply<Hemispheres>(event, "Hemispheres");
      _histMJetHigh->fill(hemi.scaledM2high());
      _histBW->fill(hemi.Bmax() );
      _histBT->fill(hemi.Bsum() );
      // C-parameter
      const ParisiTensor& parisi = apply<ParisiTensor>(event, "Parisi");
      _histC->fill(parisi.C());
      // jets
      // const FastJets&  durjet = apply<FastJets>(event, "DurhamJets");
      // const FastJets& jadejet = apply<FastJets>(event, "JadeJets");
      // if (durjet .clusterSeq()) _histy23Durham->fill( durjet.clusterSeq()->exclusive_ymerge_max(2));
      // if (jadejet.clusterSeq()) _histy23JADE  ->fill(jadejet.clusterSeq()->exclusive_ymerge_max(2));
      // jet cone
      Vector3 jetAxis=thrust.thrustAxis();
      if(hemi.highMassDirection()) jetAxis *=-1.;
      for(const Particle & p : fs.particles()) {
	const double thetaij = 180.*jetAxis.angle(p.p3().unit())/M_PI;
	double jcef = p.E()/ Evis;
	_histJCEF->fill(          thetaij,jcef);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_histThrust     , 1./sumOfWeights());
      scale(_histMJetHigh   , 1./sumOfWeights());
      scale(_histBT         , 1./sumOfWeights());
      scale(_histBW         , 1./sumOfWeights());
      scale(_histOblateness , 1./sumOfWeights());
      scale(_histC          , 1./sumOfWeights());
      // scale(_histy23JADE    , 1./sumOfWeights());
      // scale(_histy23Durham  , 1./sumOfWeights());
      scale(_histEEC        , 180./M_PI/sumOfWeights());
      scale(_histAEEC       , 180./M_PI/sumOfWeights());
      scale(_histJCEF       , 180./M_PI/sumOfWeights());
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr  _histThrust, _histMJetHigh, _histBT, _histBW;
    Histo1DPtr _histOblateness, _histC, _histy23JADE, _histy23Durham;
    Histo1DPtr _histEEC, _histAEEC, _histJCEF;

    ///@}


  };


  RIVET_DECLARE_PLUGIN(SLD_1995_I378545);

}
