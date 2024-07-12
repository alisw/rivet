// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

#define I_KNOW_THE_INITIAL_QUARKS_PROJECTION_IS_DODGY_BUT_NEED_TO_USE_IT
#include "Rivet/Projections/InitialQuarks.hh"

namespace Rivet {

  /// @brief B* production at LEP1
  class DELPHI_1995_I395026 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DELPHI_1995_I395026);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // // Initialise and register projections
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(InitialQuarks(), "IQF");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_ctheta1, 5,1,1);
      book(_h_ctheta2, "/TMP/ctheta",20,-1.,1.);
      book(_h_z      , 4,1,1);
      book(_c_hadron , "/TMP/chadron");
      book(_c_bottom , "/TMP/cbottom");
      book(_c_bStar  , "/TMP/cbStar ");
      book(_c_B      , "/TMP/cB     ");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");

      int flavour = 0;
      const InitialQuarks& iqf = apply<InitialQuarks>(event, "IQF");

      // If we only have two quarks (qqbar), just take the flavour.
      // If we have more than two quarks, look for the highest energetic q-qbar pair.
      if (iqf.particles().size() == 2) {
        flavour = iqf.particles().front().abspid();
      }
      else {
        map<int, double> quarkmap;
        for (const Particle& p : iqf.particles()) {
          if (quarkmap[p.pid()] < p.E()) {
            quarkmap[p.pid()] = p.E();
          }
        }
        double maxenergy = 0.;
        for (int i = 1; i <= 5; ++i) {
          if (quarkmap[i]+quarkmap[-i] > maxenergy) {
            flavour = i;
          }
        }
      }
      if     (flavour==5) _c_bottom->fill();
      _c_hadron->fill();
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      
      // loop over the particles
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==513 or Cuts::abspid==523 or
									       Cuts::abspid==511 or Cuts::abspid==521)) {
	int sign = p.pid()/p.abspid();
	// count number of Bs not from mixing or B*
	if(p.abspid()==511 || p.abspid()==521) {
	  if(p.parents()[0].abspid()==p.abspid()) continue;
	  if(p.parents()[0].abspid()==513 || p.parents()[0].abspid()==523) continue;
	  _c_B->fill(); 
	}
	// B*
	else {
	  _c_bStar->fill(); 
	  double xE = p.momentum().t()/meanBeamMom;
	  _h_z->fill(xE);
	  Particle decay;
	  if(p.children().size()!=2) continue;
	  int mid = p.abspid()-2;
	  if(p.children()[0].pid()==sign*mid && 
	     p.children()[1].pid()==22) {
	    decay = p.children()[1];
	  }
	  else if(p.children()[1].pid()==sign*mid && 
		  p.children()[0].pid()==22) {
	    decay = p.children()[0];
	  }
	  else
	    continue;
	  LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	  Vector3 e1z = p.p3().unit();	
	  FourMomentum pp = boost.transform(decay.momentum());
	  Vector3 axis1 = boost.transform(decay.momentum()).p3().unit();
	  double ctheta = e1z.dot(axis1);
	  _h_ctheta1->fill(ctheta);
	  _h_ctheta2->fill(ctheta);
	}
      }
    }
    
    pair<double,double> calcRho(Histo1DPtr hist) {
      if(hist->numEntries()==0.) return make_pair(0.,0.);
      double sum1(0.),sum2(0.);
      for (auto bin : hist->bins() ) {
	double Oi = bin.area();
	if(Oi==0.) continue;
	double ai = 0.125*( -bin.xMin()*(3.+sqr(bin.xMin())) + bin.xMax()*(3.+sqr(bin.xMax())));
	double bi = 0.375*( -bin.xMin()*(1.-sqr(bin.xMin())) + bin.xMax()*(1.-sqr(bin.xMax())));
	double Ei = bin.areaErr();
	sum1 += sqr(bi/Ei);
	sum2 += bi/sqr(Ei)*(Oi-ai);
      }
      return make_pair(sum2/sum1,sqrt(1./sum1));
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      // spectrum
      scale(_h_z      ,1./_c_hadron->val());
      // polarization
      scale(_h_ctheta1,1./_c_hadron->val());
      normalize(_h_ctheta2);
      pair<double,double> rho = calcRho(_h_ctheta2);
      Scatter2DPtr h_rho;
      book(h_rho, 3,1,1);
      h_rho->addPoint(0.5, rho.first, make_pair(0.5,0.5),
		      make_pair(rho.second,rho.second) );
      // no of B* per hadronic Z
      double val = _c_bStar->val()/_c_hadron->val();
      double err = val*sqrt(sqr(_c_bStar->err()/_c_bStar->val())+sqr(_c_hadron->err()/_c_hadron->val()));
      Scatter2DPtr h_nBS;
      book(h_nBS,2,1,1);
      h_nBS->addPoint(1.,val,make_pair(0.5,0.5),make_pair(err,err) );
      // no of B* per b bbar
      val = _c_bStar->val()/_c_bottom->val();
      err = val*sqrt(sqr(_c_bStar->err()/_c_bStar->val())+sqr(_c_bottom->err()/_c_bottom->val()));
      Scatter2DPtr h1;
      book(h1,1,1,1);
      h1->addPoint(1.,val,make_pair(0.5,0.5),make_pair(err,err) );
      Counter ctemp = *_c_bStar+*_c_B;
      // no of B*/B+B*
      val = _c_bStar->val()/ctemp.val();
      err = val*sqrt(sqr(_c_bStar->err()/_c_bStar->val())+sqr(ctemp.err()/ctemp.val()));
      h1->addPoint(2.,val,make_pair(0.5,0.5),make_pair(err,err) );
      // average x_E
      val = _h_z->xMean();
      err = _h_z->xStdErr();
      h1->addPoint(3.,val,make_pair(0.5,0.5),make_pair(err,err) );
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_ctheta1, _h_ctheta2, _h_z;
    CounterPtr _c_hadron,_c_bottom,_c_bStar,_c_B;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(DELPHI_1995_I395026);


}
