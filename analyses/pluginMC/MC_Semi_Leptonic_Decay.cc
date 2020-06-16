// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief MC Semi-leptonic decay analysis
  class MC_Semi_Leptonic_Decay : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_Semi_Leptonic_Decay);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // B decays
      bookHistos( 511,  -413,-11,5.3);
      bookHistos(-511,   413, 11,5.3);
      bookHistos( 511,  -413,-13,5.3);
      bookHistos(-511,   413, 13,5.3);
      bookHistos( 521,  -423,-11,5.3);
      bookHistos(-521,   423, 11,5.3);
      bookHistos( 521,  -423,-13,5.3);
      bookHistos(-521,   423, 13,5.3);
      bookHistos( 511,  -411,-11,5.3);
      bookHistos(-511,   411, 11,5.3);
      bookHistos( 511,  -411,-13,5.3);
      bookHistos(-511,   411, 13,5.3);
      bookHistos( 521,  -421,-11,5.3);
      bookHistos(-521,   421, 11,5.3);
      bookHistos( 521,  -421,-13,5.3);
      bookHistos(-521,   421, 13,5.3);
      bookHistos( 511,  -415,-11,5.3);
      bookHistos(-511,   415, 11,5.3);
      bookHistos( 511,  -415,-13,5.3);
      bookHistos(-511,   415, 13,5.3);
      bookHistos( 521,  -425,-11,5.3);
      bookHistos(-521,   425, 11,5.3);
      bookHistos( 521,  -425,-13,5.3);
      bookHistos(-521,   425, 13,5.3);
      bookHistos( 511,-10411,-11,5.3);
      bookHistos(-511, 10411, 11,5.3);
      bookHistos( 511,-10411,-13,5.3);
      bookHistos(-511, 10411, 13,5.3);
      bookHistos( 521,-10421,-11,5.3);
      bookHistos(-521, 10421, 11,5.3);
      bookHistos( 521,-10421,-13,5.3);
      bookHistos(-521, 10421, 13,5.3);
      bookHistos( 511,-10413,-11,5.3);
      bookHistos(-511, 10413, 11,5.3);
      bookHistos( 511,-10413,-13,5.3);
      bookHistos(-511, 10413, 13,5.3);
      bookHistos( 521,-10423,-11,5.3);
      bookHistos(-521, 10423, 11,5.3);
      bookHistos( 521,-10423,-13,5.3);
      bookHistos(-521, 10423, 13,5.3);
      bookHistos( 511,-20413,-11,5.3);
      bookHistos(-511, 20413, 11,5.3);
      bookHistos( 511,-20413,-13,5.3);
      bookHistos(-511, 20413, 13,5.3);
      bookHistos( 521,-20423,-11,5.3);
      bookHistos(-521, 20423, 11,5.3);
      bookHistos( 521,-20423,-13,5.3);
      bookHistos(-521, 20423, 13,5.3);

      bookHistos( 511,  -213,-11,5.3);
      bookHistos(-511,   213, 11,5.3);
      bookHistos( 511,  -213,-13,5.3);
      bookHistos(-511,   213, 13,5.3);
      bookHistos( 521,   113,-11,5.3);
      bookHistos(-521,   113, 11,5.3);
      bookHistos( 521,   113,-13,5.3);
      bookHistos(-521,   113, 13,5.3);
      bookHistos( 521,   223,-11,5.3);
      bookHistos(-521,   223, 11,5.3);
      bookHistos( 521,   223,-13,5.3);
      bookHistos(-521,   223, 13,5.3);
      bookHistos( 511,  -211,-11,5.3);
      bookHistos(-511,   211, 11,5.3);
      bookHistos( 511,  -211,-13,5.3);
      bookHistos(-511,   211, 13,5.3);
      bookHistos( 521,   111,-11,5.3);
      bookHistos(-521,   111, 11,5.3);
      bookHistos( 521,   111,-13,5.3);
      bookHistos(-521,   111, 13,5.3);
      bookHistos( 521,   221,-11,5.3);
      bookHistos(-521,   221, 11,5.3);
      bookHistos( 521,   221,-13,5.3);
      bookHistos(-521,   221, 13,5.3);
      bookHistos( 521,   331,-11,5.3);
      bookHistos(-521,   331, 11,5.3);
      bookHistos( 521,   331,-13,5.3);
      bookHistos(-521,   331, 13,5.3);
      
      // D decays
      bookHistos( 411,-311,-11,1.9);
      bookHistos(-411, 311, 11,1.9);
      bookHistos( 411,-311,-13,1.9);
      bookHistos(-411, 311, 13,1.9);
      bookHistos( 421,-321,-11,1.9);
      bookHistos(-421, 321, 11,1.9);
      bookHistos( 421,-321,-13,1.9);
      bookHistos(-421, 321, 13,1.9);      
      bookHistos( 411,-313,-11,1.9);
      bookHistos(-411, 313, 11,1.9);
      bookHistos( 411,-313,-13,1.9);
      bookHistos(-411, 313, 13,1.9);
      bookHistos( 421,-323,-11,1.9);
      bookHistos(-421, 323, 11,1.9);
      bookHistos( 421,-323,-13,1.9);
      bookHistos(-421, 323, 13,1.9);
    
      bookHistos( 411,-10313,-11,1.9);
      bookHistos(-411, 10313, 11,1.9);
      bookHistos( 411,-10313,-13,1.9);
      bookHistos(-411, 10313, 13,1.9);
      bookHistos( 421,-10323,-11,1.9);
      bookHistos(-421, 10323, 11,1.9);
      bookHistos( 421,-10323,-13,1.9);
      bookHistos(-421, 10323, 13,1.9);
      
      bookHistos( 411,-315,-11,1.9);
      bookHistos(-411, 315, 11,1.9);
      bookHistos( 411,-315,-13,1.9);
      bookHistos(-411, 315, 13,1.9);
      bookHistos( 421,-325,-11,1.9);
      bookHistos(-421, 325, 11,1.9);
      bookHistos( 421,-325,-13,1.9);
      bookHistos(-421, 325, 13,1.9);

      bookHistos( 411, 111,-11,1.9);
      bookHistos(-411, 111, 11,1.9);
      bookHistos( 411, 111,-13,1.9);
      bookHistos(-411, 111, 13,1.9);
      bookHistos( 411, 221,-11,1.9);
      bookHistos(-411, 221, 11,1.9);
      bookHistos( 411, 221,-13,1.9);
      bookHistos(-411, 221, 13,1.9);
      bookHistos( 411, 331,-11,1.9);
      bookHistos(-411, 331, 11,1.9);
      bookHistos( 411, 331,-13,1.9);
      bookHistos(-411, 331, 13,1.9);
      bookHistos( 421,-211,-11,1.9);
      bookHistos(-421, 211, 11,1.9);
      bookHistos( 421,-211,-13,1.9);
      bookHistos(-421, 211, 13,1.9);

      bookHistos( 411, 113,-11,1.9);
      bookHistos(-411, 113, 11,1.9);
      bookHistos( 411, 113,-13,1.9);
      bookHistos(-411, 113, 13,1.9);
      bookHistos( 411, 223,-11,1.9);
      bookHistos(-411, 223, 11,1.9);
      bookHistos( 411, 223,-13,1.9);
      bookHistos(-411, 223, 13,1.9);
      bookHistos( 421,-213,-11,1.9);
      bookHistos(-421, 213, 11,1.9);
      bookHistos( 421,-213,-13,1.9);
      bookHistos(-421, 213, 13,1.9);
      
      // D_s decays
      bookHistos( 431, 221,-11,1.9);
      bookHistos(-431, 221, 11,1.9);
      bookHistos( 431, 221,-13,1.9);
      bookHistos(-431, 221, 13,1.9);
      bookHistos( 431, 331,-11,1.9);
      bookHistos(-431, 331, 11,1.9);
      bookHistos( 431, 331,-13,1.9);
      bookHistos(-431, 331, 13,1.9);
      bookHistos( 431, 333,-11,1.9);
      bookHistos(-431, 333, 11,1.9);
      bookHistos( 431, 333,-13,1.9);
      bookHistos(-431, 333, 13,1.9);
      
      bookHistos( 431, 311,-11,1.9);
      bookHistos(-431,-311, 11,1.9);
      bookHistos( 431, 311,-13,1.9);
      bookHistos(-431,-311, 13,1.9);
      bookHistos( 431, 313,-11,1.9);
      bookHistos(-431,-313, 11,1.9);
      bookHistos( 431, 313,-13,1.9);
      bookHistos(-431,-313, 13,1.9);
    }

    void bookHistos(int id1, int id2, int ilep, double mass) {
      _incoming.push_back(id1);
      _outgoing.push_back(id2);
      _outgoingL.push_back(ilep);
      std::ostringstream title;
      title << "h_" << abs(id1);
      if(id1>0) title << "p";
      else     title << "m";
      title << "_" << abs(id2);
      if(id2>0) title << "p";
      else                 title << "m";
      title << "_" << abs(ilep);
      if(ilep>0) title << "p";
      else       title << "m";
      title << "_";
      _energy.push_back(Histo1DPtr());
      book(_energy.back(), title.str()+"energy",200,0.0,0.5*mass/MeV);
      _scale .push_back(Histo1DPtr());
      book(_scale.back(), title.str()+"scale"  ,200,0.0,mass/MeV);
    }

    void findDecayProducts(const Particle & mother,
			   unsigned int & nstable,
			   Particles& lp, Particles& lm,
			   Particles& nu, Particles& nub,
			   Particles& out) {
      for(const Particle & p : mother.children()) {
	int id = p.pid();
     	if ( id == PID::EMINUS || id == PID::MUON ) {
      	  lm .push_back(p);
      	  ++nstable;
      	}
      	else if (id == PID::EPLUS || id == PID::ANTIMUON) {
      	  lp .push_back(p);
      	  ++nstable;
      	}
	else if ( id == PID::NU_E || id == PID::NU_EBAR ) {
      	  nu .push_back(p);
      	  ++nstable;
      	}
      	else if (id == PID::NU_MU || id == PID::NU_MUBAR ) {
      	  nub.push_back(p);
      	  ++nstable;
      	}
	else if (PID::isMeson(id)) {
	  out.push_back(p);
	  ++nstable;
	}
      	else if ( !p.children().empty() ) {
      	  findDecayProducts(p,nstable,lp,lm,nu,nub,out);
      	}
      	else
      	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over unstable particles
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles()) {
	int id = meson.pid();
	// spin 0 mesons
	if(!PID::isMeson(id)) continue;
	if(abs(id)%10!=1) continue;
	unsigned int nstable(0);
	Particles lp, lm, nu, nub, out;
	findDecayProducts(meson,nstable,lp,lm,nu,nub,out);
	if(nstable!=3 || out.size()!=1) continue;
	int ilep=0;
	FourMomentum plep,pmnu=out[0].momentum();
	double me2(0.);
	if( lp.size()==1 && nu.size()==1 && out.size()==1 ) {
	  if(nu[0].pid()  != -lp[0].pid()+1) continue;
	  ilep =  lp[0].pid();
	  plep = nu[0].momentum()+lp[0].momentum();
	  pmnu += nu[0].momentum();
	  me2 = lp[0].mass2();
	}
	else if( lm.size()==1 && nub.size()==1 && out.size()==1 ) {
	  if(nub[0].pid() != -lm[0].pid()-1) continue;
	  ilep =  lm[0].pid();
	  plep = nub[0].momentum()+lm[0].momentum();
	  pmnu += nub[0].momentum();
	  me2 = lm[0].mass2();
	}
	else
	  continue;
	// check if histos already exist
	unsigned int iloc=0;
	bool found(false);
	while(!found&&iloc<_incoming.size()) {
	  if(_incoming[iloc] == id  &&
	     _outgoing[iloc] == out[0].pid() &&
	     ilep==_outgoingL[iloc]) found=true; 
	  else ++iloc;
	}
	if(!found) {
	  MSG_WARNING("MC_Semi_Leptonic_Decay" << id << " " << out[0].pid() << " " << " " << ilep << " "
		      << meson.mass() << "\n");
	  continue;
	}
	// add the results to the histogram
	_scale[iloc]->fill(plep.mass()/MeV);
	double ee = 0.5/meson.mass()*(meson.mass2()-pmnu.mass2()+me2);
	_energy[iloc]->fill(ee/MeV);
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<_energy.size();++ix) {
	normalize(_energy);
	normalize(_scale );
      }
    }

    //@}


    /// @name Histograms
    //@{
    /**
     *  PDG codes of the decaying mesons
     */ 
    vector<long> _incoming;
    
    /**
     *  PDG codes of the decay products
     */
    vector<long> _outgoing;  
    
    /**
     *  Identidies of the leptons
     */
    vector<long> _outgoingL;
    
    /**
     *  Histograms
     */
    //@{
    /**
     *  The lepton energy
     */
    vector<Histo1DPtr> _energy;
    
    /**
     *  The \f$q\f$ value
     */
    vector<Histo1DPtr> _scale;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_Semi_Leptonic_Decay);


}
