// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi -> gamma p pbar
  class BESIII_2012_I1079921 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2012_I1079921);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==443 or
						Cuts::abspid==100443);
      declare(ufs, "UFS");
      DecayedParticles PSI(ufs);
      declare(PSI, "PSI");
      declare(Beam(), "Beams");
      // book histos
      if (isCompatibleWithSqrtS(3.1,0.001)) {
	for(unsigned int ix=0;ix<4;++ix)
	  book(_h[ix],1,1,1+ix);
      }
      else if(isCompatibleWithSqrtS(3.686,0.001)) {
	book(_h[0],2,1,1);
      }
      else {
	cerr << "testing problem " << sqrtS() << "\n";
	throw Error("Unexpected sqrtS ! Only 3.1 and 3.686 GeV are supported");
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
      static const map<PdgId,unsigned int> & mode = { { 2212,1}, { -2212,1},{ 22,1}};
      DecayedParticles PSI = apply<DecayedParticles>(event, "PSI");
      if( PSI.decaying().size()!=1) vetoEvent;
      if(sqrtS()>3.2 && PSI.decaying()[0].pid()==443) vetoEvent;
      if(!PSI.modeMatches(0,3,mode)) vetoEvent;
      const Particle & pp   = PSI.decayProducts()[0].at( 2212)[0];
      const Particle & pbar = PSI.decayProducts()[0].at(-2212)[0];
      const Particle & gam  = PSI.decayProducts()[0].at(   22)[0];
      double mass = (pp.momentum()+pbar.momentum()).mass()-pp.mass()-pbar.mass();
      _h[0]->fill(mass);
      if(!_h[1] || mass>0.05) return;
      double cTheta = axis.dot(gam.p3().unit());
      // photon angle
      if(vetoPhoton(abs(cTheta))) vetoEvent;
      _h[1]->fill(cTheta);
      // remaining angles
      LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(PSI.decaying()[0].momentum().betaVec());
      FourMomentum pGamma    = boost1.transform(gam.momentum());
      FourMomentum pppbar = boost1.transform(pp.momentum()+pbar.momentum());
      Vector3 e1z = pGamma.p3().unit();
      Vector3 e1y = e1z.cross(axis).unit();
      Vector3 e1x = e1y.cross(e1z).unit();
      LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pppbar.betaVec());
      Vector3 axis2 = boost2.transform(boost1.transform(pp.momentum())).p3().unit();
      _h[2]->fill(e1z.dot(axis2));
      double phi = atan2(axis2.dot(e1y),axis2.dot(e1x))/M_PI*180.;
      _h[3]->fill(phi);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<4;++ix)
	if(_h[ix]) normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2012_I1079921);

}
