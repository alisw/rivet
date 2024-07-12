// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Upsilon(1S) -> gamma pi+pi- K+K-
  class BABAR_2018_I1667191 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2018_I1667191);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==200553 or
						Cuts::abspid==100553);
      declare(ufs, "UFS");
      DecayedParticles UPS(ufs);
      UPS.addStable( 553);
      declare(UPS, "UPS");
      // histograms
      for(unsigned int ix=0;ix<2;++ix)
	book(_h_mass[ix],1+ix,1,1);
      book(_h_pi,3,1,1);
      for(unsigned int ix=0;ix<3;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  book(_h_angle[ix][iy],4+ix,1,1+iy);
    }


   /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& gamma,
			   Particles & pip, Particles & pim,
			   Particles & Kp , Particles & Km,unsigned int & nstable) {
      for(const Particle & p: mother.children()) {
	if     (p.pid()== 211) pip.push_back(p);
	else if(p.pid()==-211) pim.push_back(p);
	else if(p.pid()== 321) Kp .push_back(p);
	else if(p.pid()==-321) Km .push_back(p);
	else if(p.pid()==22) gamma.push_back(p);
	else if(p.children().empty())
	  nstable+=1;
	else {
	  findDecayProducts(p, gamma,pip,pim,Kp,Km,nstable);
	  
	}
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 553,1},{ 211,1}, {-211,1}};
      DecayedParticles UPS = apply<DecayedParticles>(event, "UPS");
      // loop over particles
      for(unsigned int ix=0;ix<UPS.decaying().size();++ix) {
	// check pi+pi- upslion(1S) decay mode
      	if (!UPS.modeMatches(ix,3,mode)) continue;
	const Particle  & ups1 = UPS.decayProducts()[ix].at( 553)[0];
	const Particle  & pips = UPS.decayProducts()[ix].at( 211)[0];
	const Particle  & pims = UPS.decayProducts()[ix].at(-211)[0];
	// boost to rest frame
	LorentzTransform boost;
	if (UPS.decaying()[ix].p3().mod() > 1*MeV)
	  boost = LorentzTransform::mkFrameTransformFromBeta(UPS.decaying()[ix].momentum().betaVec());
	FourMomentum ppipi = boost.transform(pips.momentum()+pims.momentum());
	Vector3 axis1 = ppipi.p3().unit();
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(ppipi);
	FourMomentum ppi = boost1.transform(boost.transform(pips.momentum()));
	_h_pi->fill(abs(ppi.p3().unit().dot(axis1)));
	unsigned int nstable=0;
	Particles gamma,pip,pim,Kp,Km;
	findDecayProducts(ups1, gamma,pip,pim,Kp,Km,nstable);
	if(gamma.size()!=1 || nstable!=0) continue;
	// gamma pi+pi-
	if(Kp.empty()&&Km.empty()&&pip.size()==1&&pim.size()==1) {
	  ppipi = pip[0].momentum()+pim[0].momentum();
	  double mpipi = ppipi.mass();
	  _h_mass[0]->fill(mpipi);
	  axis1 *=-1;
	  LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(boost.transform(ups1.momentum()).betaVec());
	  FourMomentum pgamma = boost2.transform(boost.transform(gamma[0].momentum()));
	  Vector3 axis2 = pgamma.p3().unit();
	  double cGamma = axis2.dot(axis1);
	  ppipi = boost2.transform(boost1.transform(ppipi));
	  LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(ppipi.betaVec());
	  Vector3 axis3 = boost3.transform(boost2.transform(boost1.transform(pip[0].momentum()))).p3().unit();
	  double cH = axis3.dot(axis2);
	  int iloc=-1;
	  if(mpipi>0.6&&mpipi<1.)          iloc=0;
	  else if(mpipi>1.092&&mpipi<1.46) iloc=1;
	  if(iloc>=0) {
	    _h_angle[iloc][0]->fill(cGamma);
	    _h_angle[iloc][1]->fill(cH);
	  }
	}
	// gamma K+K-
	else if (pip.empty()&&pim.empty()&&Kp.size()==1&&Km.size()==1) {
	  FourMomentum pKK = Kp[0].momentum()+Km[0].momentum(); 
	  double mKK = pKK.mass();
	  _h_mass[1]->fill(mKK);
	  axis1 *=-1;
	  LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(boost.transform(ups1.momentum()).betaVec());
	  FourMomentum pgamma = boost2.transform(boost.transform(gamma[0].momentum()));
	  Vector3 axis2 = pgamma.p3().unit();
	  double cGamma = axis2.dot(axis1);
	  pKK = boost2.transform(boost1.transform(pKK));
	  LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pKK.betaVec());
	  Vector3 axis3 = boost3.transform(boost2.transform(boost1.transform(Kp[0].momentum()))).p3().unit();
	  double cH = axis3.dot(axis2);
	  if(mKK>1.424 && mKK<1.62) {
	    _h_angle[2][0]->fill(cGamma);
	    _h_angle[2][1]->fill(cH);
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	normalize(_h_mass[ix],1.,false);
      normalize(_h_pi,1.,false);
      for(unsigned int ix=0;ix<3;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  normalize(_h_angle[ix][iy],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[2];
    Histo1DPtr _h_pi;
    Histo1DPtr _h_angle[3][2];
    /// @}

  };


  RIVET_DECLARE_PLUGIN(BABAR_2018_I1667191);

}
