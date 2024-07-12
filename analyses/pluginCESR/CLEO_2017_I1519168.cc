// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D0 -> K+ K- pi+ pi- and 2pi+2pi-
  class CLEO_2017_I1519168 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_2017_I1519168);


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
      for(unsigned int ix=0;ix<8;++ix)
	book(_h[ix   ],1,1,1+ix);
      for(unsigned int ix=0;ix<6;++ix)
	book(_h[ix+ 8],2,1,1+ix);
      for(unsigned int ix=0;ix<4;++ix)
	book(_h[ix+14],3,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // define the decay mode
      static const map<PdgId,unsigned int> & mode1   = { { 211,2}, { -211,2}};
      static const map<PdgId,unsigned int> & mode2   = { { 321,1}, { -321,1}, { 211,1}, { -211,1}};
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for(unsigned int ix=0;ix<D0.decaying().size();++ix) {
	int sign = D0.decaying()[ix].pid()/421;
	if ( D0.modeMatches(ix,4,mode1)) {
	  const Particles & pip= D0.decayProducts()[ix].at( sign*211);
	  const Particles & pim= D0.decayProducts()[ix].at(-sign*211);
	  bool KSveto=false;
	  set<double> mpm; 
	  for(unsigned int ix=0;ix<2;++ix) {
	    for(unsigned int iy=0;iy<2;++iy) {
	      double m2 = (pip[ix].momentum()+pim[iy].momentum()).mass2();
	      double m = sqrt(m2);
	      mpm.insert(m2);
	      if(abs(m-0.497611)<0.0165) KSveto=true;
	    }
	  }
	  if(KSveto) continue;
	  _h[0]->fill(*mpm.begin());
	  _h[1]->fill(*mpm.rbegin());
	  for(const double & m2 : mpm) _h[2]->fill(m2);
	  FourMomentum ppp = pip[0].momentum()+pip[1].momentum();
	  _h[3]->fill(ppp.mass2());
	  FourMomentum pmm = pim[0].momentum()+pim[1].momentum();
	  double m2ppm[2] = {(ppp+pim[0].momentum()).mass2(),(ppp+pim[1].momentum()).mass2()};
	  if(m2ppm[0]>m2ppm[1]) swap(m2ppm[0],m2ppm[1]);
	  _h[4]->fill(m2ppm[0]);
	  _h[5]->fill(m2ppm[1]);
	  double m2mmp[2] = {(pmm+pip[0].momentum()).mass2(),(pmm+pip[1].momentum()).mass2()};
	  if(m2mmp[0]>m2mmp[1]) swap(m2mmp[0],m2mmp[1]);
	  _h[6]->fill(m2mmp[0]);
	  _h[7]->fill(m2ppm[1]);
	}
	else if ( D0.modeMatches(ix,4,mode2)) {
	  const Particles & Kp = D0.decayProducts()[ix].at( sign*321);
	  const Particles & Km = D0.decayProducts()[ix].at(-sign*321);
	  const Particles & pip= D0.decayProducts()[ix].at( sign*211);
	  const Particles & pim= D0.decayProducts()[ix].at(-sign*211);
	  double mpipi = (pip[0].momentum()+pim[0].momentum()).mass();
	  if(abs(mpipi-0.497611)<0.0165) continue;
	  _h[ 8]->fill((Kp [0].momentum()+Km [0].momentum()).mass2());
	  _h[ 9]->fill((Kp [0].momentum()+pip[0].momentum()).mass2());
	  _h[10]->fill((Kp [0].momentum()+pim[0].momentum()).mass2());
	  _h[11]->fill((Km [0].momentum()+pip[0].momentum()).mass2());
	  _h[12]->fill((Km [0].momentum()+pim[0].momentum()).mass2());
	  _h[13]->fill(sqr(mpipi));
	  _h[14]->fill((Kp [0].momentum()+Km [0].momentum()+pip[0].momentum()).mass2());
	  _h[15]->fill((Kp [0].momentum()+Km [0].momentum()+pim[0].momentum()).mass2());
	  _h[16]->fill((Kp [0].momentum()+pip[0].momentum()+pim[0].momentum()).mass2());
	  _h[17]->fill((Km [0].momentum()+pip[0].momentum()+pim[0].momentum()).mass2());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<18;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[18];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_2017_I1519168);

}
