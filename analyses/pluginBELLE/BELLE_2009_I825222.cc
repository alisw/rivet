// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> Xs gamma
  class BELLE_2009_I825222 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2009_I825222);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::abspid==521 || Cuts::abspid==511), "UFS");
      // Book histograms
      book(_h_br, 1, 1, 1);
      book(_p_E,  1, 1, 2);
      book(_p_E2,"TMP/E2",refData(1,1,3));
      book(_nBottom, "TMP/BottomCounter");
    }

    void findDecayProducts(const Particle& mother,
                           unsigned int& nK0, unsigned int& nKp,
			   unsigned int& nKm) {
      for (const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::KPLUS )      ++nKp;
        else if (id == PID::KMINUS ) ++nKm;
        else if (id == PID::K0S)     ++nK0;
        else if (id == PID::PI0 || id == PID::PIPLUS || id == PID::PIMINUS) {
	  continue;
        }
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nK0, nKp, nKm);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Loop over bottoms
      for (const Particle& bottom : apply<UnstableParticles>(event, "UFS").particles()) {
	// remove mixing entries etc
	if(bottom.children()[0].abspid()==bottom.abspid()) continue;
        _nBottom->fill();
	FourMomentum pgamma(0.,0.,0.,0.);
	unsigned int ngamma = 0;
        for (const Particle & child : bottom.children()) {
	  if (child.pid() == PID::PHOTON) {
            ngamma += 1;
            pgamma += child.momentum();
          }
	}
	if (ngamma != 1) continue;
        unsigned int nK0(0),nKp(0),nKm(0);
        FourMomentum p_tot(0,0,0,0);
        findDecayProducts(bottom, nK0, nKp, nKm);
        unsigned int nk = nKp-nKm+nK0;
        if (nk % 2 == 1) {
	  const LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(bottom.momentum().betaVec());
	  double eGamma = boost.transform(pgamma).E();
	  for(auto bin : _h_br->bins()) {
	    if(eGamma>bin.xMin()) {
	      _h_br->fill(bin.xMid());
	      _p_E ->fill(bin.xMid(),eGamma);
	      _p_E2->fill(bin.xMid(),sqr(eGamma));
	    }
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 1e4 for br ormalization and 0.1 for bin width
      scale(_h_br, 1e3/_nBottom->sumW());
      // dispersion
      Scatter2DPtr dispersion;
      book(dispersion,1,1,3);
      for(unsigned int ix=0;ix<_p_E2->numBins();++ix) {
	double val = _p_E2->bins()[ix].mean()-sqr(_p_E->bins()[ix].mean());
	double err = val*sqrt(sqr(_p_E2->bins()[ix].stdErr()/_p_E2->bins()[ix].mean())+
			      4.*sqr(_p_E->bins()[ix].stdErr()/_p_E->bins()[ix].mean()));
	double dx = 0.5*_p_E2->bins()[ix].xWidth();
	dispersion->addPoint(_p_E2->bins()[ix].xMid(),val,make_pair(dx,dx),make_pair(err,err));
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_br;
    Profile1DPtr _p_E,_p_E2;
    CounterPtr _nBottom;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2009_I825222);

}
