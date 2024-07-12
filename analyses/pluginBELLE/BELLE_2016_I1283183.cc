// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> Xs l+l-
  class BELLE_2016_I1283183 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2016_I1283183);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::abspid==521 || Cuts::abspid==511), "UFS");
      // book the profile hist
      book(_p,2,1,1);
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
	bool charm = false;
	Particles lp,lm;
	unsigned int nK0(0),nKp(0),nKm(0);
	findDecayProducts(charm,bottom, nK0, nKp, nKm,lp,lm);
	if(charm) continue;
        unsigned int nk = nKp-nKm+nK0;
	if( nk % 2 == 0) continue;
	if (lp.size()!=1 || lm.size()!=1 || lp[0].pid()!=-lm[0].pid()) continue;
	if(bottom.pid()>0) swap(lp,lm);
	double q2 = (lp[0].momentum()+lm[0].momentum()).mass2();
	// veto region valid for muons but not electrons
	if(lm[0].pid()==PID::ELECTRON && ( (q2>7.3 && q2<8.1) || (q2>11.8 && q2<12.5) )) continue;
	// first boost to bottom frame
	const LorentzTransform boost  = LorentzTransform::mkFrameTransformFromBeta(bottom.momentum().betaVec());
	FourMomentum plp = boost.transform(lp[0] .momentum());
	FourMomentum plm = boost.transform(lm[0] .momentum());
	FourMomentum pB  = boost.transform(bottom.momentum());
	const LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta((plp+plm).betaVec());
	plp = boost2.transform(plp);
	pB  = boost .transform(pB );
	double cTheta = plp.p3().unit().dot(pB.p3().unit());
	_p->fill(q2, cTheta>0 ? 1 : -1);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
    }

    /// @}


    /// @name Histograms
    /// @{
    Profile1DPtr _p;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2016_I1283183);

}
