// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi -> gamma KS0 KS0
  class BESIII_2018_I1689296 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1689296);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==443);
      declare(ufs, "UFS");
      DecayedParticles PSI(ufs);
      PSI.addStable(PID::K0S);
      declare(PSI, "PSI");
      declare(Beam(), "Beams");
      // hisotgrams
      for(unsigned int ix=0;ix<3;++ix) {
	if(ix<2) book(_h_mass[ix],1,1,1+ix);
	book(_h_angle[ix],2,1,1+ix);
      }
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
      static const map<PdgId,unsigned int> & mode = { { 310,2},{ 22,1}};
      DecayedParticles PSI = apply<DecayedParticles>(event, "PSI");
      if( PSI.decaying().size()!=1) vetoEvent;
      if(!PSI.modeMatches(0,3,mode)) vetoEvent;
      // particles
      const Particles & K0  = PSI.decayProducts()[0].at(310);
      const Particle  & gam = PSI.decayProducts()[0].at( 22)[0];
      double mKK = (K0[0].momentum()+K0[1].momentum()).mass();
      _h_mass[0]->fill(mKK);
      for(unsigned int ix=0;ix<2;++ix)
	_h_mass[1]->fill((gam.momentum()+K0[ix].momentum()).mass());
      double cTheta = axis.dot(gam.p3().unit());
      if(vetoPhoton(abs(cTheta))) vetoEvent;
      _h_angle[0]->fill(cTheta);
      // remaining angles
      LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(PSI.decaying()[0].momentum().betaVec());
      FourMomentum pGamma = boost1.transform(gam.momentum());
      FourMomentum pKK    = boost1.transform(K0[0].momentum()+K0[1].momentum());
      Vector3 e1z = pGamma.p3().unit();
      Vector3 e1y = e1z.cross(axis).unit();
      Vector3 e1x = e1y.cross(e1z).unit();
      LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pKK.betaVec());
      Vector3 axis2 = boost2.transform(boost1.transform(K0[0].momentum())).p3().unit();
      _h_angle[1]->fill(e1z.dot(axis2));
      double phi = atan2(e1y.dot(axis2),e1x.dot(axis2));
      _h_angle[2]->fill(phi);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix) {
	if(ix<2) normalize(_h_mass[ix],1.,false);
	normalize(_h_angle[ix],1.,false);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[2], _h_angle[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2018_I1689296);

}
