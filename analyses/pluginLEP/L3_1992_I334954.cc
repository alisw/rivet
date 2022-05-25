// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include "Rivet/Projections/ParisiTensor.hh"

namespace Rivet {


  /// @brief Event shapes at MZ
  class L3_1992_I334954 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(L3_1992_I334954);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      const FinalState fs;
      declare(fs, "FS");
      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      declare(Sphericity(fs), "Sphericity");
      Thrust thrust(fs);
      declare(thrust, "Thrust"    );
      declare(Hemispheres(thrust), "Hemispheres");
      declare(ParisiTensor(fs), "Parisi");
      declare(FastJets(cfs, FastJets::DURHAM, 0.7), "DurhamJets");
      declare(FastJets(cfs, FastJets::JADE  , 0.7), "JadeJets"  );
      // histograms
      book(_histThrust    ,  1, 1, 1);
      book(_histMajor     ,  2, 1, 1);
      book(_histMinor     ,  3, 1, 1);
      book(_histOblateness,  4, 1, 1);
      book(_histJade      ,  6, 1, 1);
      book(_histDurham    ,  7, 1, 1);
      book(_histSphericity, 10, 1, 1);
      book(_histAplanarity, 11, 1, 1);
      book(_histC         , 12, 1, 1);
      book(_histD         , 13, 1, 1);
      book(_histMJetHeavy , 14, 1, 1);
      book(_histMJetLight , 15, 1, 1);
      book(_histMult      , 16, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // 5 charged particles
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      if(cfs.particles().size()<5) vetoEvent;
      // charged particle mult
      _histMult->fill(cfs.particles().size());
      // Sphericity related
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");
      _histSphericity->fill(sphericity.sphericity());
      _histAplanarity->fill(sphericity.aplanarity());
      // thrust related
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      _histThrust    ->fill(thrust.thrust());
      _histMajor     ->fill(thrust.thrustMajor());
      _histMinor     ->fill(thrust.thrustMinor());
      _histOblateness->fill(thrust.oblateness());
      // hemisphere related
      const Hemispheres& hemi = apply<Hemispheres>(event, "Hemispheres");
      _histMJetHeavy->fill(hemi.scaledM2high());
      _histMJetLight->fill(hemi.scaledM2low());
      // C and D
      const ParisiTensor& parisi = apply<ParisiTensor>(event, "Parisi");
      _histC->fill(parisi.C());
      _histD->fill(parisi.D());
      // jet rates
      const FastJets& durjet = apply<FastJets>(event, "DurhamJets");
      double y23 = durjet.clusterSeq()->exclusive_ymerge_max(2);
      _histDurham->fill(y23);
      const FastJets& jadejet = apply<FastJets>(event, "JadeJets");
      y23 = jadejet.clusterSeq()->exclusive_ymerge_max(2);
      _histJade->fill(y23);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_histThrust    , 1./sumOfWeights());
      scale(_histMajor     , 1./sumOfWeights());
      scale(_histMinor     , 1./sumOfWeights());
      scale(_histOblateness, 1./sumOfWeights());
      scale(_histJade      , 1./sumOfWeights());
      scale(_histDurham    , 1./sumOfWeights());
      scale(_histSphericity, 1./sumOfWeights());
      scale(_histAplanarity, 1./sumOfWeights());
      scale(_histC         , 1./sumOfWeights());
      scale(_histD         , 1./sumOfWeights());
      scale(_histMJetHeavy , 1./sumOfWeights());
      scale(_histMJetLight , 1./sumOfWeights());
      // percentage and bin width
      scale(_histMult, 100.*2./sumOfWeights());
      // scale(, 1./sumOfWeights());
      // scale(, 1./sumOfWeights());
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _histThrust, _histMajor, _histMinor, _histOblateness;
    Histo1DPtr _histJade,_histDurham;
    Histo1DPtr _histSphericity, _histAplanarity;
    Histo1DPtr _histC, _histD;
    Histo1DPtr _histMJetHeavy, _histMJetLight, _histMult;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(L3_1992_I334954);

}
