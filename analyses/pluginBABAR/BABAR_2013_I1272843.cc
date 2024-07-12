// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> Xs l+l-
  class BABAR_2013_I1272843 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2013_I1272843);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::abspid==521 || Cuts::abspid==511), "UFS");
      // Book histograms
      for(unsigned int ix=0;ix<2;++ix) book(_p_A[ix],1,1+ix,4);
      for(unsigned int iy=0;iy<3;++iy) {
	book(_h_mX[iy],2,1,1+iy);
	for(unsigned int ix=0;ix<2;++ix)
	  book(_h_q2[ix][iy],1,1+ix,1+iy);
      }
      book(_nBottom, "TMP/BottomCounter");
    }

    void findDecayProducts(bool & charm, const Particle& mother,
                           unsigned int& nK0, unsigned int& nKp,
			   unsigned int& nKm, Particles & lp, Particles & lm) {
      for (const Particle & p : mother.children()) {
        int id = p.pid();
	if(PID::isCharmHadron(p.pid())) charm = true;
	else if ( id == PID::POSITRON || id == PID::ANTIMUON) lp.push_back(p);
	else if ( id == PID::ELECTRON || id == PID::MUON    ) lm.push_back(p);
        else if ( id == PID::KPLUS )    ++nKp;
        else if (id == PID::KMINUS )    ++nKm;
        else if (id == PID::K0S)        ++nK0;
        else if (id == PID::PI0 || id == PID::PIPLUS || id == PID::PIMINUS) {
	  continue;
        }
        else if ( !p.children().empty() ) {
          findDecayProducts(charm, p, nK0, nKp, nKm,lp,lm);
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
	bool charm = false;
	Particles lp,lm;
	unsigned int nK0(0),nKp(0),nKm(0);
	findDecayProducts(charm,bottom, nK0, nKp, nKm,lp,lm);
	if(charm) continue;
        unsigned int nk = nKp-nKm+nK0;
	if( nk % 2 == 0) continue;
	if (lp.size()!=1 || lm.size()!=1 || lp[0].pid()!=-lm[0].pid()) continue;
	FourMomentum pl = (lp[0].momentum()+lm[0].momentum());
	FourMomentum pX=bottom.momentum()-pl;
	double q2 = pl.mass2();
	double mX = pX.mass();
	if(lm[0].pid()==11) {
	  _h_q2[0][0]->fill(q2);
	  _h_q2[1][0]->fill(q2);
	  _h_mX[0]->fill(mX);
	}
	else {
	  _h_q2[0][1]->fill(q2);
	  _h_q2[1][1]->fill(q2);
	  _h_mX[1]->fill(mX);
	}
	_h_q2[0][2]->fill(q2);
	_h_q2[1][2]->fill(q2);
	_h_mX[2]->fill(mX);
	// A_CP doesn't include charmonium regions
	if(q2<6.8 || (q2>10.1 && q2<12.9) || q2>14.2) {
	  double wgt = bottom.pid()>0 ? -1 : 1;
	  _p_A[0]->fill(q2,wgt);
	  _p_A[1]->fill(q2,wgt);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int iy=0;iy<3;++iy) {
	scale(_h_mX[iy],1e6/ *_nBottom);
	for(unsigned int ix=0;ix<2;++ix)
	  scale(_h_q2[ix][iy],1e6/ *_nBottom);
      }
    }
    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_q2[2][3],_h_mX[3];
    Profile1DPtr _p_A[2];
    CounterPtr _nBottom;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2013_I1272843);

}
