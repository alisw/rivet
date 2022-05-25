// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief chi_b(2S) decays
  class CLEOII_1992_I32611 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOII_1992_I32611);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(),"UFS");
      // Book histograms
      // averages
      _h_N_aver  = {Profile1DPtr(),Profile1DPtr(),Profile1DPtr()};
      _h_R2_aver = {Profile1DPtr(),Profile1DPtr(),Profile1DPtr()};
      book(_h_N_aver [0],1,3,1);
      book(_h_N_aver [1],1,2,1);
      book(_h_N_aver [2],1,1,1);
      book(_h_R2_aver[0],1,3,2);
      book(_h_R2_aver[1],1,2,2);
      book(_h_R2_aver[2],1,1,2);
      // dists
      _h_N  = {Histo1DPtr(),Histo1DPtr(),Histo1DPtr()};
      _h_R2 = {Histo1DPtr(),Histo1DPtr(),Histo1DPtr()};
      book(_h_N [0],2,1,1);
      book(_h_N [1],2,1,2);
      book(_h_N [2],2,1,3);
      book(_h_R2[0],3,1,1);
      book(_h_R2[1],3,1,2);
      book(_h_R2[2],3,1,3);
    }

    void findDecayProducts(Particle parent, Particles & children, unsigned int & nCharged) {
      for(const Particle & p: parent.children()) {
	if(p.children().empty()) {
	  if(isCharged(p)) ++nCharged;
	  children.push_back(p);
	}
	else
	  findDecayProducts(p,children,nCharged);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles chib = apply<UnstableParticles>(event,"UFS").particles(Cuts::pid==110551 or
								       Cuts::pid==120553 or
								       Cuts::pid==100555);
      for(const Particle & p : chib) {
	unsigned int iHist = (p.pid()%10)/2;
	unsigned int nCharged(0);
	Particles children;
	findDecayProducts(p,children,nCharged);
	// ncharged
	_h_N     [iHist]->fill(nCharged);
	_h_N_aver[iHist]->fill(0.5,nCharged);
	// R_2
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	vector<FourMomentum> mom;
	mom.reserve(children.size());
	for(const Particle & p2: children) {
	  mom.push_back(boost.transform(p2.momentum()));
	}
	// compute R2
	double H0(0.),H2(0.);
	for(const FourMomentum & p1:mom) {
	  double mod1 = p1.p3().mod();
	  Vector3 axis = p1.p3().unit();
	  for(const FourMomentum & p2:mom) {
	    double mod2 = p2.p3().mod();
	    double cTheta = axis.dot(p2.p3().unit());
	    H0 += mod1*mod2;
	    H2 += mod1*mod2*0.5*(3.*sqr(cTheta)-1.);
	  }
	}
	double R2=H2/H0;
	_h_R2     [iHist]->fill(R2);
	_h_R2_aver[iHist]->fill(0.5,R2);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix) {
	normalize( _h_N [ix]);
	normalize( _h_R2[ix]);
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    vector<Histo1DPtr> _h_N,_h_R2;
    vector<Profile1DPtr> _h_N_aver,_h_R2_aver;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(CLEOII_1992_I32611);

}
