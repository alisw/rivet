// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> D* omega pi
  class BELLE_2015_I1369998 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2015_I1369998);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      DecayedParticles B0(ufs);
      B0.addStable( 413);
      B0.addStable(-413);
      B0.addStable( 223);
      declare(B0, "B0");
      for(unsigned int ix=0;ix<4;++ix)
	for(unsigned int iy=0;iy<6;++iy)
	  book(_h[ix][iy],1+ix,1,1+iy);
    }

    void findChildren(const Particle & p, Particles & pim, Particles & pip,
		      Particles & pi0, unsigned int &ncount) {
      for( const Particle &child : p.children()) {
	if(child.pid()==PID::PIPLUS) {
	  pip.push_back(child);
	  ncount+=1;
	}
	else if(child.pid()==PID::PIMINUS) {
	  pim.push_back(child);
	  ncount+=1;
	}
	else if(child.pid()==PID::PI0) {
	  pi0.push_back(child);
	  ncount+=1;
	}
	else if(child.children().empty()) {
	  ncount+=1;
	}
    	else
    	  findChildren(child,pim,pip,pi0,ncount);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 413,1},{ 223,1}, {-211,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-413,1},{ 223,1}, { 211,1}};
      DecayedParticles B0 = apply<DecayedParticles>(event, "B0");
      // loop over particles
      for(unsigned int ix=0;ix<B0.decaying().size();++ix) {
	int sign = 1;
      	if(B0.decaying()[ix].pid()<0 && B0.modeMatches(ix,3,mode))
	  sign =  1;
	else if(B0.decaying()[ix].pid()>0 && B0.modeMatches(ix,3,modeCC))
	  sign = -1; 
	else 
      	  continue;
	const Particle & Dstar = B0.decayProducts()[ix].at( sign*413)[0];
	const Particle & omega = B0.decayProducts()[ix].at(      223)[0];
	const Particle & pim1  = B0.decayProducts()[ix].at(-sign*211)[0];
	// mass hists, no cuts
	double mOmegaPi2 = (omega.momentum()+pim1.momentum()).mass2();
	_h[0][0]->fill(mOmegaPi2);
	double mDstarpi2 = (Dstar.momentum()+pim1.momentum()).mass2();
	_h[1][0]->fill(mDstarpi2);
	// check the no of decay products
	if(Dstar.children().size()!=2 || omega.children().size()!=3)
	  continue;
	// find the children of the D* meson
	Particle D0,pip1;
	if(Dstar.children()[0].pid()==sign*211 &&
	   Dstar.children()[1].pid()==sign*421) {
	  pip1 = Dstar.children()[0];
	  D0   = Dstar.children()[1];
	}
	else if(Dstar.children()[1].pid()==sign*211 &&
		Dstar.children()[0].pid()==sign*421) {
	  pip1 = Dstar.children()[1];
	  D0   = Dstar.children()[0];
	}
	else
	  continue;
	// children of the omega
	unsigned int ncount=0;
	Particles pip,pim,pi0;
	findChildren(omega,pim,pip,pi0,ncount);
	if( ncount!=3 || !(pim.size()==1 && pip.size()==1 && pi0.size()==1)) continue;
	// first bottom to the B frame
	LorentzTransform boostB = LorentzTransform::mkFrameTransformFromBeta(B0.decaying()[ix].momentum().betaVec());
	FourMomentum pOmega = boostB.transform(omega.momentum());
	FourMomentum pDstar = boostB.transform(Dstar.momentum());
	FourMomentum pD     = boostB.transform(D0   .momentum());
	FourMomentum ppim1  = boostB.transform(pim1 .momentum());
	FourMomentum ppim2  = boostB.transform(pim[0].momentum());
	FourMomentum ppip1  = boostB.transform(pip1 .momentum());
	FourMomentum ppip2  = boostB.transform(pip[0].momentum());
	// ---------------------- First set of angles --------------------------------------
	// first the angles for D* (pi omega)
	LorentzTransform boostD = LorentzTransform::mkFrameTransformFromBeta(pDstar.betaVec());
	Vector3 axisD    = boostD.transform(pD   ).p3().unit();
	Vector3 axispip1 = boostD.transform(ppip1).p3().unit();
	Vector3 axisDstar = (pOmega+ppim1).p3().unit();
	double cBeta1 = axisDstar.dot(axisD);
	_h[0][3]->fill(cBeta1);
	LorentzTransform boostWpi = LorentzTransform::mkFrameTransformFromBeta((pOmega+ppim1).betaVec());
	FourMomentum pOmega2 = boostWpi.transform(pOmega);
	Vector3 axisW   = pOmega2.p3().unit();
	Vector3 axisWpi = (pOmega+ppim1).p3().unit();
	double cXi1 = axisWpi.dot(axisW);
	_h[0][1]->fill(cXi1);
	// now angle between the two planes
	Vector3 transW = axisW-cXi1*axisWpi;
	Vector3 transD = axisD-cBeta1*axisDstar;
	double psi1 = atan2(transW.cross(transD).dot(axisDstar), transW.dot(transD));
	_h[0][5]->fill(psi1);
	// normal to omega decay plane
	LorentzTransform boostW = LorentzTransform::mkFrameTransformFromBeta(pOmega2.betaVec());
	FourMomentum ppim3  = boostW.transform(boostWpi.transform(ppim2));
	FourMomentum ppip3  = boostW.transform(boostWpi.transform(ppip2));
	Vector3 nW = ppim3.p3().cross(ppip3.p3()).unit();
	// boost B decay products to omega rest frame
	FourMomentum pOmegaPi = boostW.transform(boostWpi.transform(pOmega+ppim1));
	FourMomentum pDstar2  = boostW.transform(boostWpi.transform(pDstar));
	Vector3 axisWpi2 = pOmegaPi.p3().unit();
	double cTheta1 = axisWpi2.dot(nW);
	transW = nW-cTheta1*axisWpi2;
	transD = pDstar2.p3().unit()-pDstar2.p3().unit().dot(axisWpi2)*axisWpi2;
	double phi1 = atan2(transW.cross(transD).dot(axisWpi2), transW.dot(transD));
	_h[0][2]->fill(cTheta1);
	_h[0][4]->fill(phi1);
	// ---------------------- Second set of angles --------------------------------------
	// boost to D* pi frame
	LorentzTransform boostDpi = LorentzTransform::mkFrameTransformFromBeta((pDstar+ppim1).betaVec());
	pDstar2 = boostDpi.transform(pDstar);
	pOmega2 = boostDpi.transform(pOmega);
	axisW     = pOmega2.p3().unit();
	axisDstar = pDstar2.p3().unit();
	double cXi2 = axisW.dot(axisDstar);
	_h[1][1]->fill(cXi2);
	// boost to D* rest frame
	LorentzTransform boostDstar = LorentzTransform::mkFrameTransformFromBeta(pDstar2.betaVec());
	axisW = boostDstar.transform(pOmega2).p3().unit();
	Vector3 axisDSpi= boostDstar.transform(boostDpi.transform(pDstar+ppim1)).p3().unit();
        axisD= boostDstar.transform(boostDpi.transform(pD)).p3().unit();
	double cBeta2 = axisD.dot(axisDSpi);
	_h[1][3]->fill(cBeta2);
	transW = axisW-axisW.dot(axisDSpi)*axisDSpi;
	transD = axisD-cBeta2*axisDSpi;
	double psi2 = atan2(transW.cross(transD).dot(axisDSpi), transW.dot(transD));
	_h[1][5]->fill(psi2);
	// boost to omega frame
	boostW = LorentzTransform::mkFrameTransformFromBeta(pOmega.betaVec());
	ppim3  = boostW.transform(ppim2);
	ppip3  = boostW.transform(ppip2);
	nW = ppim3.p3().cross(ppip3.p3()).unit();
	axisDSpi  = boostW.transform(pDstar+ppim1).p3().unit();
	axisDstar = boostW.transform(pDstar).p3().unit();
	double cTheta2 = axisDSpi.dot(nW);
	_h[1][2]->fill(cTheta2);
	transW = nW-cTheta2*axisDSpi;
	transD = axisDstar.unit()-axisDstar.dot(axisDSpi)*axisDSpi;
	double phi2 = atan2(transW.cross(transD).dot(axisDSpi), transW.dot(transD));
	_h[1][4]->fill(psi2);
	// restricted plots
	if(abs(cTheta1)>.5) {
	  _h[2][0]->fill(mOmegaPi2);
	}
	else {
	  _h[2][1]->fill(mOmegaPi2);
	  _h[2][3]->fill(cBeta1);
	  _h[2][5]->fill(psi1);
	  _h[3][1]->fill(mDstarpi2);
	  _h[3][3]->fill(cTheta2);
	  _h[3][5]->fill(phi2);
	}
	if(cXi2>-.4) {
	  _h[2][2]->fill(cBeta1);
	  _h[2][4]->fill(psi1);
	  _h[3][0]->fill(mDstarpi2);
	  _h[3][2]->fill(cTheta2);
	  _h[3][4]->fill(phi2);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<4;++ix)
	for(unsigned int iy=0;iy<6;++iy)
	  normalize(_h[ix][iy],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4][6];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2015_I1369998);

}
