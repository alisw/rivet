// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief psi(2S) -> J/psi gamma gamma
  class BESIII_2012_I1097066 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2012_I1097066);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==100443);
      declare(ufs, "UFS");
      DecayedParticles PSI(ufs);
      PSI.addStable(PID::JPSI);
      declare(PSI, "PSI");
      declare(Beam(), "Beams");
      // book histograms
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // get the axis, direction of incoming electron
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      Vector3 axis;
      if(beams.first.pid()>0) axis = beams.first .momentum().p3().unit();
      else                    axis = beams.second.momentum().p3().unit();
      // find the J/psi decays
      static const map<PdgId,unsigned int> & mode = { { 443,1},{ 22,2}};
      DecayedParticles PSI = apply<DecayedParticles>(event, "PSI");
      if( PSI.decaying().size()!=1) vetoEvent;
      if(!PSI.modeMatches(0,3,mode)) vetoEvent;
      if( PSI.decaying()[0].children().size()!=3) vetoEvent;
      // particles
      const Particle  & JPsi = PSI.decayProducts()[0].at(443)[0];
      const Particles & gam  = PSI.decayProducts()[0].at( 22);
      Vector3 norm = gam[0].p3().cross(gam[1].p3()).unit();
      _h[0]->fill(abs(axis.dot(norm)));
      if(JPsi.children().size()!=2) vetoEvent;
      if(JPsi.children()[0].pid()!=-JPsi.children()[1].pid()) vetoEvent;
      if(JPsi.children()[0].abspid()!=PID::EMINUS &&
	 JPsi.children()[0].abspid()!=PID::MUON) vetoEvent;
      Particle lm = JPsi.children()[0];
      Particle lp = JPsi.children()[1];
      if(lm.pid()<0) swap(lm,lp);
      Vector3 axis2=JPsi.momentum().p3().unit();
      LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(JPsi.momentum().betaVec());
      Vector3 axis3 = boost.transform(lm.momentum()).p3().unit();
      _h[1]->fill(abs(axis2.dot(axis3)));
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2012_I1097066);

}
