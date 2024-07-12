// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi -> gamma eta eta'
  class BESIII_2022_I2135117 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2135117);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==443);
      declare(ufs, "UFS");
      DecayedParticles PSI(ufs);
      PSI.addStable(PID::ETA);
      PSI.addStable(PID::ETAPRIME);
      declare(PSI, "PSI");
      // histos
      for(unsigned int ix=0;ix<3;++ix) {
	book(_h_angle[ix],2,1,1+ix);
	book(_h_mass [ix],1,1,1+ix);
      }
      book(_h_angle[3],1,1,4);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // find the J/psi decays
      static const map<PdgId,unsigned int> & mode = { { 22,1},{ 221,1},{ 331,1}};
      DecayedParticles PSI = apply<DecayedParticles>(event, "PSI");
      for(unsigned int ix=0;ix<PSI.decaying().size();++ix) {
	if(!PSI.modeMatches(ix,3,mode)) continue;
	const Particle  & eta  = PSI.decayProducts()[ix].at(221)[0];
	const Particle  & etap = PSI.decayProducts()[ix].at(331)[0];
	const Particle  & gam  = PSI.decayProducts()[ix].at( 22)[0];
	double mEE = (eta.momentum()+etap.momentum()).mass();
	double mEG = (gam.momentum()+eta .momentum()).mass();
	if(abs(mEG-1.019461)<0.04) continue;
	_h_mass[0]->fill(mEE);
	_h_mass[1]->fill(mEG);
	_h_mass[2]->fill((gam.momentum()+etap.momentum()).mass());
	// angles
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(PSI.decaying()[0].momentum().betaVec());
	FourMomentum pGamma = boost1.transform(gam.momentum());
	FourMomentum pEE    = boost1.transform(eta.momentum()+etap.momentum());
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pEE.betaVec());
	Vector3 axis2 = boost2.transform(boost1.transform(eta.momentum())).p3().unit();
	double cTheta = pGamma.p3().unit().dot(axis2);
	_h_angle[3]->fill(cTheta);
	if( mEE>1.5 && mEE<1.7)
	  _h_angle[0]->fill(cTheta);
	else if(mEE>1.7 && mEE<2.)
	  _h_angle[1]->fill(cTheta);
	else if(mEE>2. && mEE<3.2)
	  _h_angle[2]->fill(cTheta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<4;++ix) {
	normalize(_h_angle[ix],1.,false);
	if(ix<3) normalize(_h_mass[ix],1.,false);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[3],_h_angle[4];
    /// @}
  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2135117);

}
