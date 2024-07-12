// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Bbar0 -> J/[psi K- pi+
  class BELLE_2014_I1312626 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2014_I1312626);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      DecayedParticles B0(ufs);
      B0.addStable( 443);
      declare(B0, "B0");
      // histograms
      for(unsigned int ix=0;ix<7;++ix)
	book(_h_mass[ix],1,1,1+ix);
      for(unsigned int ix=0;ix<2;++ix)
	book(_h_angle[ix],2,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 211,1}, {-321,1}, { 443,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-211,1}, { 321,1}, { 443,1}};
      DecayedParticles B0 = apply<DecayedParticles>(event, "B0");
      for(unsigned int ix=0;ix<B0.decaying().size();++ix) {
	int sign=1;
	if (B0.decaying()[ix].pid()<0 && B0.modeMatches(ix,3,mode  )) sign = 1;
	if (B0.decaying()[ix].pid()>0 && B0.modeMatches(ix,3,modeCC)) sign =-1;
	else continue;
      	const Particle & pip  = B0.decayProducts()[ix].at( sign*211)[0];
       	const Particle & Km   = B0.decayProducts()[ix].at(-sign*321)[0];
       	const Particle & JPsi = B0.decayProducts()[ix].at( 443)[0];
	double mpiJ2 = (pip.momentum()+JPsi.momentum()).mass2();
	double mKpi2 = (pip.momentum()+Km  .momentum()).mass2();
	if(mpiJ2<16.)       _h_mass[0]->fill(mKpi2);
	else if(mpiJ2<19.)  _h_mass[1]->fill(mKpi2);
	else                _h_mass[2]->fill(mKpi2);
	if(mKpi2<1.2)       _h_mass[3]->fill(mpiJ2);
	else if(mKpi2<2.05) _h_mass[4]->fill(mpiJ2);
	else if(mKpi2<3.2 ) _h_mass[5]->fill(mpiJ2);
	else                _h_mass[6]->fill(mpiJ2);
	// mass cuts for angular variables
	if(mKpi2<1.2 || mpiJ2<16. || mpiJ2>19.) continue;
	// helicity angle find J.psi leptonic children
	if(JPsi.children().size()!=2) continue;
	if(JPsi.children()[0].pid()!=-JPsi.children()[1].pid()) continue;
	if(JPsi.children()[0].abspid()!=PID::EMINUS &&
	   JPsi.children()[0].abspid()!=PID::MUON) continue;
	Particle lm = JPsi.children()[0];
	Particle lp = JPsi.children()[1];
	if(lm.pid()<0) swap(lm,lp);
       	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(B0.decaying()[ix].momentum().betaVec());
       	FourMomentum ppsi = boost1.transform(JPsi.momentum());
      	Vector3 axis1 = -ppsi.p3().unit();
      	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(ppsi.betaVec());
      	FourMomentum plp = boost2.transform(boost1.transform(lp.momentum()));
      	double cL = plp.p3().unit().dot(axis1);
	_h_angle[0]->fill(cL);
	Vector3 LTrans = plp.p3() - cL*plp.p3().mod()*axis1;
	FourMomentum pKpi = pip.momentum()+Km.momentum();
      	LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pKpi.betaVec());
      	FourMomentum ppi = boost3.transform(boost1.transform(pip.momentum()));
	double cPi = ppi.p3().unit().dot(axis1);
	Vector3 PTrans = ppi.p3() - cPi*ppi.p3().mod()*axis1;
        double phi = atan2(LTrans.cross(PTrans).dot(axis1), LTrans.dot(PTrans));
	_h_angle[1]->fill(phi);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<7;++ix)
	normalize(_h_mass[ix],1.,false);
      for(unsigned int ix=0;ix<2;++ix)
	normalize(_h_angle[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[7],_h_angle[2];
    /// @}
  };


  RIVET_DECLARE_PLUGIN(BELLE_2014_I1312626);

}
