// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

namespace Rivet {


  /// @brief MC decay M -> M l+ l-
  class MC_Meson_Meson_Leptons_Decay : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_Meson_Meson_Leptons_Decay);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      
      // Initialise and register projections
      declare(UnstableParticles(),"UFS");

      // book the Histos
      // pi0 dalitz
      bookHistos(111,22,11,0.140);
      // eta dalitz
      bookHistos(221,22,11,0.55);
      bookHistos(221,22,13,0.55);
      // eta' dalitz
      bookHistos(331,22,11,0.96);
      bookHistos(331,22,13,0.96);
      // omega -> pi0
      bookHistos(223,111,11,0.8);
      bookHistos(223,111,13,0.8);
      // phi -> pi0
      bookHistos(333,111,11,1.1);
      bookHistos(333,111,13,1.1);
      // phi -> eta
      bookHistos(333,221,11,1.1);
      bookHistos(333,221,13,1.1);
      // J/psi dalitz
      bookHistos(443,22,11,3.1);
      bookHistos(443,22,13,3.1);
      // B -> s gamma
      bookHistos(511,313,11,5.3);
      bookHistos(511,313,13,5.3);
    }

    void bookHistos(int id1, int id2, int il, double dM) {
      if(abs(id2)%10==3 || id2==22) {
	_incoming_P.push_back(id1);
	_outgoingV.push_back(id2);
	_outgoingf_P.push_back(il);
	std::ostringstream title;
	title << "h2_" << abs(id1);
	if(id1>0) title << "p";
	else      title << "m";
	title << "_" << abs(id2);
	if(id2>0) title << "p";
	else                    title << "m";
	title << "_" << il << "_";
	_mff_P   .push_back(Histo1DPtr());
	book(_mff_P.back(), title.str()+"mff"   , 100, 0., dM);
	_mVf   .push_back(Histo1DPtr());
	book(_mVf.back(), title.str()+"mVf"   , 100, 0., dM);
	_mVfbar.push_back(Histo1DPtr());
	book(_mVfbar.back(), title.str()+"mVfbar", 100, 0., dM);
      }
      else {
	_incomingV.push_back(id1);
	_outgoingP.push_back(id2);
	_outgoingf_V.push_back(il);
	std::ostringstream title;
	title << "h_" << abs(id1);
	if(id1>0) title << "p";
	else      title << "m";
	title << "_" << abs(id2);
	if(id2>0) title << "p";
	else                    title << "m";
	title << "_" << il << "_";
	_mff_V   .push_back(Histo1DPtr());
	book(_mff_V.back(), title.str()+"mff"   , 100, 0., dM);
	_mPf   .push_back(Histo1DPtr());
	book(_mPf.back(), title.str()+"mPf"   , 100, 0., dM);
	_mPfbar.push_back(Histo1DPtr());
	book(_mPfbar.back(), title.str()+"mPfbar", 100, 0., dM);
      }
    }

    void findDecayProducts(const Particle & mother,
			   unsigned int & nstable,
			   Particles& lp, Particles& lm,
			   Particles& scalar,
			   Particles& vector) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
	if ( id == PID::EMINUS || id == PID::MUON ) {
       	  lm.push_back(p);
       	  ++nstable;
       	}
	else if (id == PID::EPLUS || id == PID::ANTIMUON) {
	  lp.push_back(p);
	  ++nstable;
	}
	else if (abs(id)%10==1 && PID::isMeson(id)) {
	  scalar.push_back(p);
	  ++nstable;
	}
	else if ((abs(id)%10==3 && PID::isMeson(id)) ||
		 id==PID::PHOTON ) {
	  vector.push_back(p);
	  ++nstable;
	}
	else if ( !p.children().empty() ) {
	  findDecayProducts(p,nstable,lp,lm,scalar,vector);
	}
	else
	  ++nstable;
      }
    }
    
    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over unstable particles
      for(const Particle& iMeson : apply<UnstableParticles>(event, "UFS").particles()) {
	// only consider scalar/vector mesons
	long pid = iMeson.pid();
	if(!PID::isMeson(pid)) continue;
	if(abs(pid)%10!=3 and abs(pid)%10!=1 ) continue;
	Particles lp,lm,scalar,vector;
	unsigned int nstable(0);
	findDecayProducts(iMeson,nstable,lp,lm,scalar,vector);
	if(nstable!=3 || lp.size()!=1 || lm.size()!=1 || lp[0].pid()!=-lm[0].pid()) continue;
	if(scalar.size()==1) {
	  // check if we already have this decay
	  unsigned int ix=0; bool found(false); 
	  while(!found&&ix<_incomingV.size()) {
	    if(_incomingV[ix]==pid && _outgoingP[ix]==scalar[0].pid() &&
	       _outgoingf_V[ix]==lm[0].pid()) {
	      found=true;
	    }
	    else {
	      ++ix;
	    }
	  }
	  // create a new graph if needed
	  if(!found) {
	    MSG_WARNING("MC_Meson_Meson_Leptons_Decay S" << iMeson.pid() << " " << scalar[0].pid() << " "
			<< iMeson.mass() << "\n");
	    continue;
	  }
	  // add the results to the histogram
	  _mff_V   [ix]->fill((lm    [0].momentum()+lp[0].momentum()).mass());
	  _mPf   [ix]->fill((scalar[0].momentum()+lm[0].momentum()).mass());
	  _mPfbar[ix]->fill((scalar[0].momentum()+lp[0].momentum()).mass());
	}
	else if(vector.size()==1) {
	  // check if we already have this decay
	  unsigned int ix=0; bool found(false); 
	  while(!found&&ix<_incoming_P.size()) {
	    if(_incoming_P[ix]==pid && _outgoingV[ix]==vector[0].pid() &&
	       _outgoingf_P[ix]==lm[0].pid()) {
	      found=true;
	    }
	    else {
	      ++ix;
	    }
	  }
	  // create a new graph if needed
	  if(!found) {
	    MSG_WARNING("MC_Meson_Meson_Leptons_Decay V" << iMeson.pid() << " " << vector[0].pid() << " "
			<< iMeson.mass() << "\n");
	    continue;
	  }
	  // add the results to the histogram
	  _mff_P   [ix]->fill((lm    [0].momentum()+lp[0].momentum()).mass());
	  _mVf   [ix]->fill((vector[0].momentum()+lm[0].momentum()).mass());
	  _mVfbar[ix]->fill((vector[0].momentum()+lp[0].momentum()).mass());
	}
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      
      // normalize to unity V->P
      for(unsigned int ix=0;ix<_mff_V.size();++ix) {
        normalize(_mff_V);
        normalize(_mPf);
        normalize(_mPfbar);
      }
      // normalize to unity P->V
      for(unsigned int ix=0;ix<_mff_P.size();++ix) {
        normalize(_mff_P);
        normalize(_mVf);
        normalize(_mVfbar);
      }

    }

    //@}



    /// @name Histograms for V -> P
    //@{
    /**
     *  PDG codes of the incoming particles
     */
    vector<long> _incomingV;
    
    /**
     *  PDG codes of the outgoing pseudoscalar mesons
     */
    vector<long> _outgoingP;
    
    /**
     *  PDG codes of the outgoing fermion
     */
    vector<long> _outgoingf_V;
    
    /**
     *  Histograms for the mass of the fermion-antifermion pair
     */
    vector<Histo1DPtr> _mff_V;
    
    /**
     *  Histograms for the masses of the pseudoscalar and the fermion
     */
    vector<Histo1DPtr> _mPf;
    
    /**
     *  Histograms for the masses of the pseudoscalar and the antifermion
     */
    vector<Histo1DPtr> _mPfbar;
    //@}
    
    /// @name Histograms P->V
    //@{
    /**
     *  PDG codes of the incoming_P particles
     */
    vector<long> _incoming_P;
    
    /**
     *  PDG codes of the outgoing vector mesons
     */
    vector<long> _outgoingV;
    
    /**
     *  PDG codes of the outgoing fermion
     */
    vector<long> _outgoingf_P;
    
    /**
     *  Histograms for the mass of the fermion-antifermion pair
     */
    vector<Histo1DPtr> _mff_P;
    
    /**
     *  Histograms for the masses of the vector and the fermion
     */
    vector<Histo1DPtr> _mVf;
    
    /**
     *  Histograms for the masses of the vector and the antifermion
     */
    vector<Histo1DPtr> _mVfbar;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_Meson_Meson_Leptons_Decay);


}
