// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class MC_Onium_PiPi_Decay : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MC_Onium_PiPi_Decay);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(),"UFS");
      // psi 2S
      bookHistos(100443,   443,0.6);
      // psi(3770)
      bookHistos( 30443,   443,0.7);
      // Upsilon (4S)
      bookHistos(300553,   553,1.2);
      bookHistos(300553,100553,0.6);
      // Upsilon (3S)
      bookHistos(200553,   553,0.9);
      bookHistos(200553,100553,0.4);
      // Upsilon (2S)
      bookHistos(100553,   553,0.6);
    }

    void bookHistos(int id1, int id2, double deltaM) {
      double twompi = 0.378;
      _incoming.push_back(id1);
      _outgoing.push_back(id2);
      std::ostringstream title;
      title << "h_" << id1 << "_" << id2 << "_";
      _mpipi.push_back(make_pair(Histo1DPtr(), Histo1DPtr()));
      book(_mpipi.back().first, title.str()+"mpippim",100,twompi/GeV,deltaM/GeV);
      book(_mpipi.back().second, title.str()+"mpi0pi0",100,twompi/GeV,deltaM/GeV);
      _hel.push_back(make_pair(Histo1DPtr(), Histo1DPtr()));
      book(_hel.back().first, title.str()+"hpippim",100,-1.,1.);
      book(_hel.back().second, title.str()+"hpi0pi0",100, 0.,1.);
    }
    
    void findDecayProducts(const Particle & mother,
			   unsigned int & nstable,
			   Particles& pip, Particles& pim,
			   Particles& pi0, Particles & onium) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
      	if ( id == PID::PIMINUS) {
	  pim.push_back(p);
	  ++nstable;
	}
       	else if (id == PID::PIPLUS) {
       	  pip.push_back(p);
       	  ++nstable;
       	}
       	else if (id == PID::PI0) {
       	  pi0.push_back(p);
       	  ++nstable;
       	}
	else if (abs(id)%1000==443 || abs(id)%1000==553) {
	  onium.push_back(p);
	  ++nstable;
	}
	else if ( !p.children().empty() ) {
	  findDecayProducts(p,nstable,pip,pim,pi0,onium);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over unstable particles
      for(const Particle& vMeson : apply<UnstableParticles>(event, "UFS").particles()) {
	int id = vMeson.pid();
	if(id%1000!=443 && id%1000!=553) continue;
	unsigned int nstable(0);
	Particles pip, pim, pi0, onium;
	findDecayProducts(vMeson,nstable,pip,pim,pi0,onium);
	// check for onium
	if(onium.size() !=1 || nstable !=3) continue;
	// check for pipi
	if( ! ((pip.size()==1 && pim.size() ==1) || pi0.size()==2)) continue;
	// check if histos already made
	unsigned int iloc=0; bool found(false);
	while(!found&&iloc<_incoming.size()) {
	  if(_incoming[iloc]==vMeson.pid()&&_outgoing[iloc]==onium[0].pid()) found=true;
	  else ++iloc;
	}
	// if histos not made, make them
	if(!found) {
	  MSG_WARNING("MC_Onium_PiPi_Decay" << vMeson.pid() << " " << onium[0].pid() << " "
		      << vMeson.mass()-onium[0].mass() << "\n");
	  continue;
	}
	// boost to rest frame of the pion pair
	FourMomentum q = vMeson.momentum()-onium[0].momentum();
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(q.betaVec());
	FourMomentum qp = onium[0].momentum();
	FourMomentum ppi= pip.size()==1 ? pip[0].momentum() : pi0[0].momentum();
	qp  = boost.transform(qp);
	ppi = boost.transform(ppi);
	double cX=-ppi.p3().unit().dot(qp.p3().unit());
	if(pi0.size()==2) {
	  _mpipi[iloc].second->fill(q.mass());
	  _hel  [iloc].second->fill(abs(cX));
	}
	else {
	  _mpipi[iloc].first->fill(q.mass());
	  _hel  [iloc].first->fill(cX);
	}
      }
    }
    
    
    /// Normalise histograms etc., after the run
    void finalize() {

      // normalize to unity
      for(unsigned int ix=0;ix<_mpipi.size();++ix) {
	normalize(_mpipi[ix].first );
	normalize(_mpipi[ix].second);
	normalize(_hel[ix].first );
	normalize(_hel[ix].second);
      }
    }

    //@}

    /**
     *  Incoming onium states
     */
    vector<long> _incoming;
    
    /**
     *  Outgoing onium states
     */
    vector<long> _outgoing;
    
    /**
     *  Histograms for the \f$\pi^+\pi^-\f$ masses
     */
    vector<pair<Histo1DPtr,Histo1DPtr> > _mpipi;
    
    /**
     *  Histmgrams for the helicity angles
     */
    vector<pair<Histo1DPtr,Histo1DPtr> > _hel;

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MC_Onium_PiPi_Decay);


}
