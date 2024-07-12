// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B_c - > jpsi pi+ and pi+pi+pi-
  class LHCB_2012_I1097092 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2012_I1097092);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==541);
      declare(ufs, "UFS");
      DecayedParticles BC(ufs);
      BC.addStable( PID::PI0);
      BC.addStable( PID::K0S);
      BC.addStable( PID::JPSI);
      declare(BC, "BC");
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_mass  [ix],1+ix,1,1);
	book(_h_ctheta[ix],3,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 443,1}, { 211,1} };
      static const map<PdgId,unsigned int> & mode1CC = { { 443,1}, {-211,1} };
      static const map<PdgId,unsigned int> & mode2   = { { 443,1}, { 211,2}, {-211,1} };
      static const map<PdgId,unsigned int> & mode2CC = { { 443,1}, {-211,2}, { 211,1} };
      DecayedParticles BC = apply<DecayedParticles>(event, "BC");
      // loop over particles
      for(unsigned int ix=0;ix<BC.decaying().size();++ix) {
	int sign = BC.decaying()[ix].pid()/BC.decaying()[ix].abspid();
	int imode=-1;
	if ((sign== 1 && BC.modeMatches(ix,2,mode1  )) ||
	    (sign==-1 && BC.modeMatches(ix,2,mode1CC))) {
	  imode=0;
	}
	else if ((sign== 1 && BC.modeMatches(ix,4,mode2  )) ||
		 (sign==-1 && BC.modeMatches(ix,4,mode2CC))) {
	  const Particle  & pim = BC.decayProducts()[ix].at(-sign*211)[0];
	  const Particles & pip = BC.decayProducts()[ix].at( sign*211);
	  _h_mass[0]->fill((pim.momentum()+pip[0].momentum()+pip[1].momentum()).mass());
	  _h_mass[1]->fill((pim.momentum()+pip[0].momentum()).mass());
	  _h_mass[1]->fill((pim.momentum()+pip[1].momentum()).mass());
	  imode=1;
	}
	else
	  continue;
	// check the J/psi decay mode
	const Particle & jpsi = BC.decayProducts()[ix].at(443)[0];
	if(jpsi.children().size()!=2 || jpsi.children()[0].pid()!=-jpsi.children()[1].pid() ||
	   jpsi.children()[0].abspid()!=13) continue;
	Particle muon = jpsi.children()[0].pid()==13 ? jpsi.children()[0] : jpsi.children()[1];
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(BC.decaying()[ix].momentum().betaVec());
	FourMomentum pJPsi = boost1.transform(jpsi.momentum());
	FourMomentum pMu   = boost1.transform(muon.momentum());
	// to j/psi rest frame
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pJPsi.betaVec());
	Vector3 axis = pJPsi.p3().unit();
	FourMomentum pp = boost2.transform(pMu);
	// calculate angle
	double cTheta = pp.p3().unit().dot(axis);
	_h_ctheta[imode]->fill(cTheta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_h_mass  [ix]);
	normalize(_h_ctheta[ix]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[2];
    Histo1DPtr _h_ctheta[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2012_I1097092);

}
