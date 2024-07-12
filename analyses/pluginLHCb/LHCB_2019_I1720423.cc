// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D+ -> K+K+K-
  class LHCB_2019_I1720423 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2019_I1720423);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_h_KpKm[0],1,1,1);
      book(_h_KpKp   ,1,1,2);
      book(_h_KpKm[1],1,1,3);
      book(_h_KpKm[2],1,1,4);
      book(_dalitz, "dalitz",50,0.9,1.8,50,1.1,1.9);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & Kp , Particles & Km) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::PIPLUS || id == PID::PIMINUS ||
	     id == PID::K0S   || id == PID::K0L ||
	     id == PID::PI0 ) {
	  ++nstable;
	}
	else if (id == PID::KPLUS) {
	  Kp.push_back(p);
	  ++nstable;
	}
	else if (id == PID::KMINUS) {
	  Km.push_back(p);
	  ++nstable;
	}
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, Kp, Km);
	}
	else
	  ++nstable;
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid== 411 )) {
	unsigned int nstable(0);
	Particles Kp, Km;
	findDecayProducts(meson, nstable, Kp, Km);
	if(nstable !=3) continue;
	if(meson.pid()<0) {
	  swap(Km,Kp);
	}
	if (Km.size()==1&&Kp.size()==2) {
	  double m1 = (Km[0].momentum()+Kp[0].momentum()).mass2();
	  double m2 = (Km[0].momentum()+Kp[1].momentum()).mass2();
	  double m3 = (Kp[0].momentum()+Kp[1].momentum()).mass2();
	  if(m1>m2) swap(m1,m2);
	  _dalitz->fill(m1,m2);
	  _h_KpKm[0]->fill(m1);
	  _h_KpKm[0]->fill(m2);
	  _h_KpKm[2]->fill(m1);
	  _h_KpKm[1]->fill(m2);
	  _h_KpKp->fill(m3);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix)
	normalize(_h_KpKm[ix]);
      normalize(_h_KpKp);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_KpKm[3],_h_KpKp;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2019_I1720423);

}
