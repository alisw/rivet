// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B0 -> psi(2S) K+ pi-
  class BELLE_2013_I1239347 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2013_I1239347);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      DecayedParticles B0(ufs);
      B0.addStable(100443);
      declare(B0, "B0");
      // histograms
      for(unsigned int ix=0;ix<8;++ix)
	book(_h_mass[ix],1,1,1+ix);
      book(_h_mass[8],2,1,1);
      for(unsigned int ix=0;ix<2;++ix)
	book(_h_angle[ix],3,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 321,1},{-211,1}, { 100443,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-321,1},{ 211,1}, { 100443,1}};
      DecayedParticles B0 = apply<DecayedParticles>(event, "B0");
      // loop over particles
      for(unsigned int ix=0;ix<B0.decaying().size();++ix) {
      	int sign = 1;
      	if (B0.decaying()[ix].pid()>0 && B0.modeMatches(ix,3,mode)) {
      	  sign=1;
      	}
      	else if  (B0.decaying()[ix].pid()<0 && B0.modeMatches(ix,3,modeCC)) {
      	  sign=-1;
      	}
      	else
      	  continue;
      	const Particle & Kp  = B0.decayProducts()[ix].at( 321*sign)[0];
      	const Particle & pim = B0.decayProducts()[ix].at(-211*sign)[0];
	const Particle & psi = B0.decayProducts()[ix].at( 100443  )[0];
	double mKpi   = (Kp .momentum()+pim.momentum()).mass();
	double m2Psipi= (psi.momentum()+pim.momentum()).mass2();
	if(m2Psipi<19.)                       _h_mass[0]->fill(sqr(mKpi));
	else if(m2Psipi>=19. && m2Psipi<20.5) _h_mass[1]->fill(sqr(mKpi));
	else if(m2Psipi>=20.5)                _h_mass[2]->fill(sqr(mKpi));

	if(mKpi<0.796) {
	  _h_mass[3]->fill(m2Psipi);
	  _h_mass[8]->fill(m2Psipi);
	}
	else if(mKpi>=0.796&&mKpi<0.996)
	  _h_mass[4]->fill(m2Psipi);
	else if(mKpi>=0.996&&mKpi<1.332) {
	  _h_mass[5]->fill(m2Psipi);
	  _h_mass[8]->fill(m2Psipi);
	}
	else if(mKpi>=1.332&&mKpi<1.532)
	  _h_mass[6]->fill(m2Psipi);
	else if(mKpi>=1.532) {
	  _h_mass[7]->fill(m2Psipi);
	  _h_mass[8]->fill(m2Psipi);
	}
	// need leptonic psi' decay for angular dists
	if(psi.children().size()!=2|| psi.children()[0].pid()!=-psi.children()[1].pid() ||
	   (psi.children()[0].abspid()!=11 &&
	    psi.children()[0].abspid()!=11)) vetoEvent;
	Particle lm = psi.children()[0];
	Particle lp = psi.children()[1];
	if(lm.pid()<0) swap(lm,lp);
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(B0.decaying()[ix].momentum().betaVec());
	FourMomentum pKstar = boost1.transform(Kp.momentum()+pim.momentum());
	FourMomentum pPsi   = boost1.transform(psi.momentum());
	// trans vector in K* frame
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pKstar.betaVec());
	FourMomentum pKp = boost2.transform(boost1.transform(Kp.momentum()));
	Vector3 axis1 = pKstar.p3().unit();
	double cTheta1 = axis1.dot(pKp.p3().unit());
	Vector3 trans1 = pKp.p3() - cTheta1*pKp.p3().mod()*axis1;
	// leptons in psi' frame
	LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pPsi.betaVec());
	FourMomentum plm = boost3.transform(boost1.transform(lm.momentum()));
	double cTheta2 = axis1.dot(plm.p3().unit());
	Vector3 trans2 = plm.p3() - cTheta2*plm.p3().mod()*axis1;
	_h_angle[0]->fill(cTheta2);
	// angle between planes
	double chi = atan2(trans1.cross(trans2).dot(axis1),trans1.dot(trans2));
	_h_angle[1]->fill(chi);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<9;++ix)
	normalize(_h_mass[ix],1.,false);
      normalize(_h_angle[0],1.,false);
      normalize(_h_angle[1],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[9],_h_angle[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2013_I1239347);

}
