// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Excited Lambda_c spectra
  class CLEOII_1994_I381696 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOII_1994_I381696);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      // book histos
      book(_h_2595,5,1,1);
      book(_h_2625,6,1,1);
      book(_r_2595,4,1,1);
      book(_r_2625,4,1,2);
      book(_c_lam,"TMP/c_lam");
    }

   

    void findDecayProducts(const Particle& mother, unsigned int & npip,
			   unsigned int & npim, unsigned int & nlam) {
      for(const Particle & p : mother.children()) {
	if(p.abspid()  == 4122) ++nlam;
	else if(p.pid()==  211) ++npip;
	else if(p.pid()== -211) ++npim;
	else if(!p.children().empty())
	  findDecayProducts(p,npip,npim,nlam);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const int id2595 = 14122;
      static const int id2625 = 4124;
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double Emax = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==id2595 or
					     Cuts::abspid==id2625)) {
	unsigned int nlam(0),npip(0),npim(0);
	findDecayProducts(p,npip,npim,nlam);
	bool isDecay = npip==1 && npim==1 && nlam==1;
	// spectrum
	if(p.abspid()==id2595) {
	  double Pmax = sqrt(sqr(Emax)-sqr(2.595));
	  double xp = p.momentum().p3().mod()/Pmax;
	  _h_2595->fill(xp);
	  if(isDecay) _r_2595->fill(10.55);
	}
	else {
	  double Pmax = sqrt(sqr(Emax)-sqr(2.625));
	  double xp = p.momentum().p3().mod()/Pmax;
	  _h_2625->fill(xp);
	  if(isDecay) _r_2625->fill(10.55);
	}
      }
      unsigned int nlam = ufs.particles(Cuts::abspid==4122).size();
      _c_lam->fill(nlam);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_2595);
      normalize(_h_2625);
      scale(_r_2595, 0.06/ *_c_lam);
      scale(_r_2625, 0.06/ *_c_lam);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_2595,_h_2625;
    Histo1DPtr _r_2595,_r_2625;
    CounterPtr _c_lam;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(CLEOII_1994_I381696);

}
