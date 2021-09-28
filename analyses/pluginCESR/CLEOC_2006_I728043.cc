// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief eta, eta' and phi in D0, D+, Ds decays
  class CLEOC_2006_I728043 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CLEOC_2006_I728043);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(UnstableParticles(),"UFS");
      // histograms
      unsigned int imin(0),imax(3);
      if(fuzzyEquals(sqrtS(),3.77,1e-2))
	imax=2;
      else if(fuzzyEquals(sqrtS(),4.17,1e-2))
	imin=2;
      else
	MSG_ERROR("Invalid CMS energy in CLEOC_2006_I728043");
      for(unsigned int ix=imin;ix<imax;++ix) {
	std::ostringstream title;
	title << "TMP/n_D_" << ix;
	book(_n_D[ix],title.str());
	book(_br_eta     [ix],1,1,ix+1);
	book(_br_etaPrime[ix],2,1,ix+1);
	book(_br_phi     [ix],3,1,ix+1);
	book(_s_eta      [ix],4,1,ix+1);
	book(_s_phi      [ix],5,1,ix+1);
      }
    }

    void fillHistos(const Particle & Dmeson, const LorentzTransform & boost) {
      Particles ssbar;
      unsigned int iMeson=0;
      if(Dmeson.abspid()==421)
	iMeson = 1;
      else if(Dmeson.abspid()==431)
	iMeson = 2;
      _n_D[iMeson]->fill();
      findDecayProducts(Dmeson,ssbar);
      for(const Particle & dec : ssbar) {
	FourMomentum p = boost.transform(dec.momentum());
	double mom=p.p3().mod();
	if(dec.pid()==221) {
	  _br_eta[iMeson]->fill(0.5);
	  _s_eta[iMeson]->fill(mom);
	}
	else if(dec.pid()==331) {
	  _br_etaPrime[iMeson]->fill(0.5);
	}
	else {
	  _br_phi[iMeson]->fill(0.5);
	  _s_phi[iMeson] ->fill(mom);
	}
      }
    }
    
    void findDecayProducts(const Particle & mother, Particles & ssbar) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
      	if (id==221 || id==331 || id==333) {
	  ssbar.push_back(p);
	  findDecayProducts(p,ssbar);
	}
	else if ( !p.children().empty() ) {
	  findDecayProducts(p,ssbar);
	}
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // find psi(3770)
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particles psi = ufs.particles(Cuts::pid==30443);
      // D_s
      if(psi.empty()) {
	LorentzTransform boost;
	for(const Particle & Dmeson : apply<UnstableParticles>("UFS",event).particles(Cuts::abspid==431)) {
	  fillHistos(Dmeson,boost);
	}
      }
      // D0 D+
      else {
	for(const Particle& p : psi) {
	  // boost to rest frame
	  LorentzTransform boost;
	  if (p.p3().mod() > 1*MeV)
	    boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	  // loop over D0 and D+ children
	  for(const Particle & Dmeson : p.children()) {
	    if(Dmeson.abspid()!=411 && Dmeson.abspid()!=421) continue;
	    fillHistos(Dmeson,boost);
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      unsigned int imin(0),imax(3);
      if(fuzzyEquals(sqrtS(),3.77,1e-2))
	imax=2;
      else if(fuzzyEquals(sqrtS(),4.17,1e-2))
	imin=2;
      else
	MSG_ERROR("Invalid CMS energy in CLEOC_2006_I728043");
      for(unsigned int ix=imin;ix<imax;++ix) {
	if(_n_D[ix]->effNumEntries()<=0.) continue;
      	scale(_br_eta     [ix], 100./ *_n_D[ix]);
      	scale(_br_etaPrime[ix], 100./ *_n_D[ix]);
      	scale(_br_phi     [ix], 100./ *_n_D[ix]);
      	scale(_s_eta      [ix], 100./ *_n_D[ix]);
      	scale(_s_phi      [ix], 100./ *_n_D[ix]);
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    CounterPtr _n_D[3];
    Histo1DPtr _br_eta[3],_br_etaPrime[3],_br_phi[3];
    Histo1DPtr _s_eta[3], _s_phi[3];
    ///@}


  };


  DECLARE_RIVET_PLUGIN(CLEOC_2006_I728043);

}
