// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Bs0 -> pi+ pi- J/psi
  class BELLE_2011_I889524 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2011_I889524);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==531);
      declare(ufs, "UFS");
      DecayedParticles BS0(ufs);
      BS0.addStable( 443);
      declare(BS0, "BS0");
      // histograms
      book(_h_mass,1,1,1);
      for(unsigned int ix=0;ix<2;++ix)
	book(_h_cTheta[ix],2,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 211,1}, {-211,1}, { 443,1}};
      DecayedParticles BS0 = apply<DecayedParticles>(event, "BS0");
      for(unsigned int ix=0;ix<BS0.decaying().size();++ix) {
       	if (!BS0.modeMatches(ix,3,mode)) continue;
      	const Particle & pip  = BS0.decayProducts()[ix].at( 211)[0];
      	const Particle & pim  = BS0.decayProducts()[ix].at(-211)[0];
      	const Particle & JPsi = BS0.decayProducts()[ix].at( 443)[0];
	double mpipi = (pip.momentum()+pim.momentum()).mass();
      	_h_mass->fill(mpipi);
	// helicity angle find J.psi leptonic children
	if(JPsi.children().size()!=2) continue;
	if(JPsi.children()[0].pid()!=-JPsi.children()[1].pid()) continue;
	if(JPsi.children()[0].abspid()!=PID::EMINUS &&
	   JPsi.children()[0].abspid()!=PID::MUON) continue;
	Particle lm = JPsi.children()[0];
	Particle lp = JPsi.children()[1];
	if(lm.pid()<0) swap(lm,lp);
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(BS0.decaying()[ix].momentum().betaVec());
	FourMomentum ppsi = boost1.transform(JPsi.momentum());
	Vector3 axis1 = ppsi.p3().unit();
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(ppsi.betaVec());
	FourMomentum plp = boost2.transform(boost1.transform(lp.momentum()));
	double cTheta = plp.p3().unit().dot(axis1);
	if     (mpipi>0.8 && mpipi<1.16) _h_cTheta[0]->fill(cTheta);
	else if(mpipi>1.3 && mpipi<1.5)  _h_cTheta[1]->fill(cTheta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_mass,1.,false);
      for(unsigned int ix=0;ix<2;++ix)
	normalize(_h_cTheta[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass,_h_cTheta[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2011_I889524);

}
