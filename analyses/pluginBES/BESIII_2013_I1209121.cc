// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi -> gamma eta eta
  class BESIII_2013_I1209121 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2013_I1209121);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==443);
      declare(ufs, "UFS");
      DecayedParticles PSI(ufs);
      PSI.addStable(PID::ETA);
      declare(PSI, "PSI");
      declare(Beam(), "Beams");
      // book histograms
      for(unsigned int ix=0;ix<4;++ix)
	book(_h[ix],1,1,1+ix);
    }

    // angle cuts due regions of BES calorimeter
    bool vetoPhoton(const double & cTheta) {
      return cTheta>0.92 || (cTheta>0.8 && cTheta<0.86);
    }
    
    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // get the axis, direction of incoming electron
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      Vector3 axis;
      if(beams.first.pid()>0)
	axis = beams.first .momentum().p3().unit();
      else
	axis = beams.second.momentum().p3().unit();
      // find the J/psi decays
      static const map<PdgId,unsigned int> & mode = { { 221,2},{ 22,1}};
      DecayedParticles PSI = apply<DecayedParticles>(event, "PSI");
      if( PSI.decaying().size()!=1) vetoEvent;
      if(!PSI.modeMatches(0,3,mode)) vetoEvent;
      // particles
      const Particles & eta = PSI.decayProducts()[0].at(221);
      const Particle  & gam = PSI.decayProducts()[0].at( 22)[0];
      _h[0]->fill((eta[0].momentum()+eta[1].momentum()).mass());
      double cTheta = axis.dot(gam.p3().unit());
      if(vetoPhoton(abs(cTheta))) vetoEvent;
      _h[1]->fill(cTheta);
      // remaining angles
      LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(PSI.decaying()[0].momentum().betaVec());
      FourMomentum pGamma = boost1.transform(gam.momentum());
      FourMomentum pEtaEta= boost1.transform(eta[0].momentum()+eta[1].momentum());
      Vector3 e1z = pGamma.p3().unit();
      Vector3 e1y = e1z.cross(axis).unit();
      Vector3 e1x = e1y.cross(e1z).unit();
      LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pEtaEta.betaVec());
      Vector3 axis2 = boost2.transform(boost1.transform(eta[0].momentum())).p3().unit();
      _h[2]->fill(e1z.dot(axis2));
      _h[3]->fill(atan2(e1y.dot(axis2),e1x.dot(axis2)));
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<4;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2013_I1209121);

}
