// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief  D0 -> KS) K+/- pi-/+
  class LHCB_2016_I1394391 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2016_I1394391);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==421);
      declare(ufs, "UFS");
      DecayedParticles D0(ufs);
      D0.addStable(PID::PI0);
      D0.addStable(PID::K0S);
      D0.addStable(PID::ETA);
      D0.addStable(PID::ETAPRIME);
      declare(D0, "D0");
      // histograms
      book(_h_Kmpip,1,1,1);
      book(_h_K0pip,1,1,2);
      book(_h_K0Km ,1,1,3);
      book(_h_Kppim,2,1,1);
      book(_h_K0pim,2,1,2);
      book(_h_K0Kp ,2,1,3);
      book(_dalitz [0],"dalitz_1",50,0.3,2.0,50,0.3,2.);
      book(_dalitz [1],"dalitz_2",50,0.3,2.0,50,0.3,2.);
    }

    double efficiency(const double & x, const double &y) {
      double X=x-2., Y=y-1.;
      static const double E0 = 5.8096, Ex = -3.645, Ey = -3.174, Ex2=  0.831,
	Exy = 2.131, Ey2 = 4.43, Ex3 = -0.427, Ex2y = 2.65, Exy2 = 1.50, Ey3 = -3.92;
      return E0 + Ex*X + Ey*Y + Ex2*sqr(X) + Ey2*sqr(Y) + Exy*X*Y +
	Ex3*pow(X,3) + Ex2y*sqr(X)*Y + Exy2*X*sqr(Y) + Ey3*pow(Y,3);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 321,1},{-211,1}, { 310,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-321,1},{ 211,1}, { 310,1}};
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for(unsigned int ix=0;ix<D0.decaying().size();++ix) {
	if( !D0.modeMatches(ix,3,mode  ) &&
	    !D0.modeMatches(ix,3,modeCC) ) continue;
	const Particles & K0 = D0.decayProducts()[ix].at(310);
	int sign = D0.decaying()[ix].pid()/421;
	const Particles & pip= D0.decayProducts()[ix].find( sign*211) == D0.decayProducts()[ix].end() ?
	  Particles() : D0.decayProducts()[ix].at( sign*211);
	const Particles & pim= D0.decayProducts()[ix].find(-sign*211) == D0.decayProducts()[ix].end() ?
	  Particles() : D0.decayProducts()[ix].at(-sign*211);
	const Particles & Kp = D0.decayProducts()[ix].find( sign*321) == D0.decayProducts()[ix].end() ?
	  Particles() : D0.decayProducts()[ix].at( sign*321);
	const Particles & Km = D0.decayProducts()[ix].find(-sign*321) == D0.decayProducts()[ix].end() ?
	  Particles() : D0.decayProducts()[ix].at(-sign*321);
	// K0S K- pi+
	if( Km.size()==1 && pip.size()==1) {
	  double mK0pip = (K0[0].momentum()+pip[0].momentum() ).mass2();
	  double mKmpip = (Km[0].momentum()+pip[0].momentum() ).mass2();
	  double mKK    = (K0[0].momentum()+Km [0].momentum() ).mass2();
	  double eff = efficiency(mKK,mK0pip);
	  _h_K0Km ->fill(mKK   ,eff);
	  _h_K0pip->fill(mK0pip,eff);
	  _h_Kmpip->fill(mKmpip,eff);
	  _dalitz[0]->fill(mKmpip,mK0pip); 
	}
	// K0S K+ pi-
	else if( Kp.size()==1 && pim.size()==1) {
	  double mK0pim = (K0[0].momentum()+pim[0].momentum() ).mass2();
	  double mKppim = (Kp[0].momentum()+pim[0].momentum() ).mass2();
	  double mKK    = (K0[0].momentum()+Kp [0].momentum() ).mass2();
	  double eff = efficiency(mKK,mK0pim);
	  _h_K0Kp ->fill(mKK   ,eff);
	  _h_K0pim->fill(mK0pim,eff);
	  _h_Kppim->fill(mKppim,eff);
	  _dalitz[1]->fill(mKppim,mK0pim); 
	}
      }
    }


    /// Normalise histograms etc., after the runbook
    void finalize() {
      normalize(_h_Kmpip);
      normalize(_h_K0pip);
      normalize(_h_K0Km );
      normalize(_h_Kppim);
      normalize(_h_K0pim);
      normalize(_h_K0Kp );
      normalize(_dalitz [0]);
      normalize(_dalitz [1]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kmpip, _h_K0pip, _h_K0Km;
    Histo1DPtr _h_Kppim, _h_K0pim, _h_K0Kp;
    Histo2DPtr _dalitz[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2016_I1394391);

}
