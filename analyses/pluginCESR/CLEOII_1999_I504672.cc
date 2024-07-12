// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Charged particle multplicity in Upsilon 4s decays
  class CLEOII_1999_I504672 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOII_1999_I504672);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_n_tot  , 1, 1, 1);
      book(_n_lep  , 1, 1, 2);
      book(_n_nolep, 1, 1, 3);
      book(_h_n    , 2, 1, 1);
      book(_n_Ups_All  ,"/TMP/n_Ups_All");
      book(_n_Ups_Lep  ,"/TMP/n_Ups_Lep");
      book(_n_Ups_NoLep,"/TMP/n_Ups_NoLep");
    }

    void findChildren(const Particle & p,int & nCharged, int & nLep) {
      bool isBottom = PID::isBottomHadron(p.pid());
      for( const Particle &child : p.children()) {
	if(child.children().empty()) {
	  if(PID::isCharged(child.pid())           ) ++nCharged;
	  if(isBottom && (child.abspid()==11||child.abspid()==13)) ++nLep;
	}
	else
	  findChildren(child,nCharged,nLep);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for (const Particle& p :  apply<FinalState>(event, "UFS").particles(Cuts::pid==300553)) {
	_n_Ups_All->fill();
	int nCharged(0),nLep(0);
	findChildren(p,nCharged,nLep);
	_h_n->fill(nCharged);
	_n_tot->fill(10.1,nCharged);
	if(nLep==2) {
	  _n_lep  ->fill(10.1,nCharged);
	  _n_Ups_Lep->fill();
	}
	else {
	  _n_nolep->fill(10.1,nCharged);
	  _n_Ups_NoLep->fill();
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_n    , 1./ *_n_Ups_All  );
      scale(_n_tot  , 1./ *_n_Ups_All  );
      scale(_n_lep  , 1./ *_n_Ups_Lep  );
      scale(_n_nolep, 1./ *_n_Ups_NoLep);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_n;
    Histo1DPtr _n_tot,_n_lep,_n_nolep;
    CounterPtr _n_Ups_All,_n_Ups_Lep,_n_Ups_NoLep;
    
    //@}


  };


  RIVET_DECLARE_PLUGIN(CLEOII_1999_I504672);

}
