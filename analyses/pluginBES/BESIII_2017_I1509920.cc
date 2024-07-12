// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief psi(2S) -> e+e- chi_c and chi_c  -> e+e- Jpsi$
  class BESIII_2017_I1509920 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2017_I1509920);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(Beam(), "Beams");
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==PID::PSI2S);
      DecayedParticles psi(ufs);
      psi.addStable(20443);
      psi.addStable(445);
      declare(psi, "PSI");
      ufs = UnstableParticles(Cuts::abspid==445 or Cuts::abspid==20443);
      DecayedParticles chi(ufs);
      chi.addStable(PID::JPSI);
      declare(chi, "CHI");
      // histograms
      for(unsigned int ix=0;ix<2;++ix)
	for(unsigned int iy=0;iy<4;++iy)
	  book(_h[ix][iy],1+ix,1,1+iy);
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
      // decaying particles
      static const map<PdgId,unsigned int> & mode1   = { {20443,1},{ 11,1}, { -11,1}};
      static const map<PdgId,unsigned int> & mode2   = { {  445,1},{ 11,1}, { -11,1}};
      static const map<PdgId,unsigned int> & mode3   = { {  443,1},{ 11,1}, { -11,1}};
      unsigned int iproj=0;
      for(const DecayedParticles & in : {apply<DecayedParticles>(event, "PSI"),apply<DecayedParticles>(event, "CHI")} ) {
	for(unsigned int ix=0;ix<in.decaying().size();++ix) {
	  int imode=-1,ichi=0;
	  if(iproj==0) {
	    if(in.modeMatches(ix,3,mode1)) {
	      imode=0;
	      ichi=20443;
	    }
	    else if(in.modeMatches(ix,3,mode2)) {
	      imode=1;
	      ichi=445;
	    }
	    else
	      continue;
	  }
	  else {
	    if(in.modeMatches(ix,3,mode3)) {
	      imode = in.decaying()[ix].pid()==20443 ? 2 : 3;
	      ichi=443;
	    }
	    else
	      continue;
	  }
	  // extract particles
	  const Particle & em  = in.decayProducts()[ix].at(  11)[0];
	  const Particle & ep  = in.decayProducts()[ix].at( -11)[0];
	  const Particle & out = in.decayProducts()[ix].at(ichi)[0];
	  FourMomentum qq = ep.momentum()+em.momentum();
	  double q = qq.mass();
	  _h[0][imode]->fill(q);
	  if(iproj==0) {
	    _h[1][imode]->fill(axis.dot((em.momentum()+ep.momentum()).p3().unit()));
	  }
	  else if(iproj==1) {
	    // boost everything to decaying particle frame
	    LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(in.decaying()[ix].momentum().betaVec());
	    FourMomentum pout = boost1.transform(out.momentum());
	    FourMomentum pe   = boost1.transform(em.momentum());
	    qq = boost1.transform(qq);
	    Vector3 axis1 = pout.p3().unit();
	    LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(qq.betaVec());
	    pe = boost2.transform(pe);
	    _h[1][imode]->fill(pe.p3().unit().dot(axis1));
	  }
	}
	++iproj;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	for(unsigned int iy=0;iy<4;++iy)
	  normalize(_h[ix][iy],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2][4];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2017_I1509920);

}
