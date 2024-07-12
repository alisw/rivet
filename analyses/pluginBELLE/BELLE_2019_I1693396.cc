// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B0 -> D*- semileptonic
  class BELLE_2019_I1693396 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2019_I1693396);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(Cuts::pid==511), "UFS");

	// Book histograms
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h[0][ix], "TMP/h_w_"     +toString(ix+1), refData(1, 1, ix+1));
	book(_h[1][ix], "TMP/h_costhl_"+toString(ix+1), refData(3, 1, ix+1));
	book(_h[2][ix], "TMP/h_costhv_"+toString(ix+1), refData(2, 1, ix+1));
	book(_h[3][ix], "TMP/h_chi_"   +toString(ix+1), refData(4, 1, ix+1));
      }
    }

    /// Perform the per-event analysis
    bool analyzeDecay(Particle mother, vector<int> ids) {
      // There is no point in looking for decays with less particles than to be analysed
      if (mother.children().size() == ids.size()) {
        bool decayfound = true;
        for (int id : ids) {
          if (!contains(mother, id)) decayfound = false;
        }
        return decayfound;
      }
      return false;
    }

    bool contains(Particle& mother, int id) {
      return any(mother.children(), HasPID(id));
    }

    double recoilW(const Particle& mother) {
      FourMomentum lepton, neutrino, meson, q;
      for(const Particle& c : mother.children()) {
        if (c.isNeutrino()) neutrino=c.mom();
        if (c.isLepton() &! c.isNeutrino()) lepton =c.mom();
        if (c.isHadron()) meson=c.mom();
      }
      q = lepton + neutrino; //no hadron before
      double mb2= mother.mom()*mother.mom();
      double mD2 = meson*meson;
      return (mb2 + mD2 - q*q )/ (2. * sqrt(mb2) * sqrt(mD2) );
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      FourMomentum pl, pnu, pB, pD, pDs, ppi;
      // Iterate of B0bar mesons
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
        pB = p.momentum();
        // Find semileptonic decays
	int iloc=-1;
        if (analyzeDecay(p, {PID::DSTARMINUS,12,-11}))      iloc = 0;
	else if(analyzeDecay(p, {PID::DSTARMINUS,14,-13}) ) iloc = 1;
	else continue;
	_h[0][iloc]->fill(recoilW(p));
        // Get the necessary momenta for the angles
        bool foundDdecay=false;
        for (const Particle & c : p.children()) {
          if ( (c.pid() == PID::DSTARMINUS)  && (analyzeDecay(c, {PID::PIMINUS, PID::D0BAR}) || analyzeDecay(c, {PID::PI0, PID::DMINUS})) ) {
	    foundDdecay=true;
	    pDs = c.momentum();
	    for (const Particle & dc : c.children()) {
	      if (dc.hasCharm()) pD = dc.momentum(); 
	      else ppi = dc.momentum(); 
	    }
	  }
          if (c.pid() == -11 || c.pid() == -13) pl  = c.momentum();
          if (c.pid() ==  12 || c.pid() ==  14) pnu = c.momentum();
	}
	// This is the angle analysis
	if (!foundDdecay) continue;
	// First boost all relevant momenta into the B-rest frame
	const LorentzTransform B_boost = LorentzTransform::mkFrameTransformFromBeta(pB.betaVec());
        // Momenta in B rest frame:
        FourMomentum lv_brest_Dstar = B_boost.transform(pDs);
        FourMomentum lv_brest_w     = B_boost.transform(pB - pDs);
        FourMomentum lv_brest_D     = B_boost.transform(pD);
        FourMomentum lv_brest_lep   = B_boost.transform(pl);
            
        const LorentzTransform Ds_boost = LorentzTransform::mkFrameTransformFromBeta(lv_brest_Dstar.betaVec());
        FourMomentum lv_Dstarrest_D     = Ds_boost.transform(lv_brest_D);
        const LorentzTransform W_boost  = LorentzTransform::mkFrameTransformFromBeta(lv_brest_w.betaVec());
        FourMomentum lv_wrest_lep       = W_boost.transform(lv_brest_lep);

        double cos_thetaV = cos(lv_brest_Dstar.p3().angle(lv_Dstarrest_D.p3()));
        _h[2][iloc]->fill(cos_thetaV);
            
        double cos_thetaL = cos(lv_brest_w.p3().angle(lv_wrest_lep.p3()));
        _h[1][iloc]->fill(cos_thetaL);

        Vector3 LTrans = lv_wrest_lep.p3()   - cos_thetaL*lv_wrest_lep.p3().perp()*lv_brest_w.p3().unit();
        Vector3 VTrans = lv_Dstarrest_D.p3() - cos_thetaV*lv_Dstarrest_D.p3().perp()*lv_brest_Dstar.p3().unit();
        float chi = atan2(LTrans.cross(VTrans).dot(lv_brest_w.p3().unit()), LTrans.dot(VTrans));
	_h[3][iloc]->fill(chi);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // efficiencies
      vector<double> eff[4][2]={{{2.72,5.72,7.7,9.1,10.03,10.61,10.74,10.67,10.23,9.1},
				 {2.68,5.66,7.66,9.05,9.91,10.43,10.6,10.52,10.04,9.14}},
				{{3.12,3.97,5.73,7.96,9.31,9.85,10.23,10.59,11.06,11.21},
				 {3.16,3.52,5.19,7.59,9.1,9.78,10.27,10.43,11.0,11.36}},
				{{11.72,11.52,11.35,10.88,10.2,9.34,8.29,7.16,6.05,4.82},
				 {11.54,11.43,11.14,10.74,10.09,9.29,8.25,7.1,5.97,4.72}},
				{{8.6,8.74,8.96,9.3,9.81,9.82,9.33,9.0,8.77,8.59},
				 {8.51,8.67,8.82,9.15,9.7,9.73,9.2,8.83,8.62,8.54}}};
      vector<double> efe[4][2]={{{0.02,0.02,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03},
				 {0.02,0.02,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03}},
				{{0.03,0.02,0.02,0.03,0.03,0.03,0.03,0.03,0.03,0.03},
				 {0.03,0.02,0.02,0.03,0.03,0.03,0.03,0.03,0.03,0.03}},
				{{0.03,0.03,0.03,0.04,0.04,0.04,0.03,0.03,0.02,0.02},
				 {0.03,0.03,0.03,0.04,0.04,0.04,0.03,0.03,0.02,0.02}},
				{{0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03},
				 {0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03}}};
      // response matricesdouble
      double response[4][2][10][10] = {{{{0.803,0.053,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
					 {0.197,0.778,0.098,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
					 {0.0,0.168,0.717,0.126,0.002,0.0,0.0,0.0,0.0,0.0},
					 {0.0,0.001,0.182,0.667,0.149,0.006,0.0,0.0,0.0,0.0},
					 {0.0,0.0,0.004,0.199,0.626,0.167,0.011,0.0,0.0,0.0},
					 {0.0,0.0,0.0,0.009,0.207,0.592,0.177,0.015,0.0,0.0},
					 {0.0,0.0,0.0,0.0,0.016,0.215,0.575,0.183,0.018,0.0},
					 {0.0,0.0,0.0,0.0,0.0,0.021,0.213,0.567,0.186,0.017},
					 {0.0,0.0,0.0,0.0,0.0,0.0,0.024,0.214,0.598,0.186},
					 {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.022,0.198,0.797}},
					{{0.961,0.024,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
					 {0.038,0.952,0.027,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
					 {0.0,0.021,0.948,0.041,0.001,0.0,0.0,0.0,0.0,0.0},
					 {0.0,0.001,0.023,0.918,0.067,0.003,0.001,0.001,0.001,0.0},
					 {0.0,0.001,0.001,0.04,0.871,0.097,0.005,0.001,0.001,0.0},
					 {0.0,0.0,0.0,0.001,0.06,0.817,0.129,0.006,0.001,0.0},
					 {0.0,0.0,0.0,0.0,0.001,0.082,0.758,0.164,0.007,0.001},
					 {0.0,0.0,0.0,0.0,0.0,0.001,0.106,0.698,0.196,0.008},
					 {0.0,0.0,0.0,0.0,0.0,0.0,0.001,0.128,0.657,0.212},
					 {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.002,0.137,0.777}}},
				       {{{0.918,0.077,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
					 {0.082,0.806,0.095,0.001,0.0,0.0,0.0,0.0,0.0,0.0},
					 {0.0,0.115,0.761,0.101,0.002,0.0,0.0,0.0,0.0,0.0},
					 {0.0,0.001,0.141,0.735,0.105,0.002,0.0,0.0,0.0,0.0},
					 {0.0,0.0,0.002,0.16,0.719,0.1,0.001,0.0,0.0,0.0},
					 {0.0,0.0,0.0,0.003,0.17,0.722,0.093,0.001,0.0,0.0},
					 {0.0,0.0,0.0,0.0,0.003,0.173,0.738,0.08,0.001,0.0},
					 {0.0,0.0,0.0,0.0,0.0,0.002,0.166,0.771,0.072,0.0},
					 {0.0,0.0,0.0,0.0,0.0,0.0,0.001,0.147,0.819,0.064},
					 {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.001,0.108,0.936}},
					{{0.659,0.129,0.011,0.003,0.002,0.002,0.002,0.004,0.013,0.144},
					 {0.151,0.691,0.132,0.012,0.004,0.002,0.002,0.002,0.004,0.016},
					 {0.015,0.141,0.697,0.147,0.016,0.005,0.002,0.002,0.002,0.005},
					 {0.005,0.012,0.134,0.671,0.162,0.018,0.005,0.002,0.002,0.002},
					 {0.002,0.004,0.013,0.14,0.634,0.155,0.016,0.004,0.002,0.002},
					 {0.002,0.002,0.004,0.015,0.155,0.633,0.141,0.013,0.004,0.003},
					 {0.002,0.002,0.002,0.004,0.018,0.163,0.67,0.136,0.012,0.004},
					 {0.005,0.002,0.002,0.002,0.005,0.015,0.147,0.695,0.14,0.015},
					 {0.016,0.004,0.002,0.002,0.002,0.004,0.013,0.132,0.691,0.15},
					 {0.142,0.013,0.003,0.002,0.002,0.002,0.003,0.012,0.13,0.659}}},
				       {{{0.812,0.051,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
					 {0.188,0.784,0.096,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
					 {0.0,0.164,0.728,0.126,0.002,0.0,0.0,0.0,0.0,0.0},
					 {0.0,0.001,0.172,0.676,0.149,0.006,0.0,0.0,0.0,0.0},
					 {0.0,0.0,0.004,0.19,0.631,0.165,0.01,0.0,0.0,0.0},
					 {0.0,0.0,0.0,0.008,0.203,0.6,0.181,0.016,0.0,0.0},
					 {0.0,0.0,0.0,0.0,0.014,0.209,0.578,0.187,0.019,0.0},
					 {0.0,0.0,0.0,0.0,0.0,0.02,0.209,0.573,0.195,0.017},
					 {0.0,0.0,0.0,0.0,0.0,0.0,0.022,0.205,0.6,0.195},
					 {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.019,0.186,0.788}},
					{{0.959,0.022,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
					 {0.039,0.955,0.012,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
					 {0.0,0.021,0.96,0.022,0.001,0.0,0.0,0.0,0.0,0.0},
					 {0.001,0.001,0.026,0.931,0.043,0.001,0.0,0.0,0.0,0.0},
					 {0.0,0.0,0.001,0.047,0.889,0.07,0.002,0.0,0.0,0.0},
					 {0.0,0.001,0.0,0.0,0.067,0.837,0.103,0.002,0.001,0.0},
					 {0.0,0.0,0.0,0.0,0.0,0.091,0.778,0.138,0.003,0.0},
					 {0.0,0.0,0.0,0.0,0.0,0.0,0.117,0.715,0.174,0.004},
					 {0.0,0.0,0.0,0.0,0.0,0.0,0.001,0.142,0.672,0.193},
					 {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.002,0.151,0.803}}},
				       {{{0.918,0.077,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
					 {0.082,0.805,0.091,0.001,0.0,0.0,0.0,0.0,0.0,0.0},
					 {0.0,0.117,0.763,0.101,0.002,0.0,0.0,0.0,0.0,0.0},
					 {0.0,0.001,0.142,0.735,0.103,0.002,0.0,0.0,0.0,0.0},
					 {0.0,0.0,0.003,0.159,0.723,0.098,0.001,0.0,0.0,0.0},
					 {0.0,0.0,0.0,0.004,0.169,0.726,0.091,0.001,0.0,0.0},
					 {0.0,0.0,0.0,0.0,0.004,0.172,0.745,0.082,0.001,0.0},
					 {0.0,0.0,0.0,0.0,0.0,0.002,0.161,0.771,0.074,0.0},
					 {0.0,0.0,0.0,0.0,0.0,0.0,0.001,0.145,0.817,0.066},
					 {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.107,0.934}},
					{{0.653,0.129,0.012,0.004,0.003,0.002,0.002,0.004,0.014,0.144},
					 {0.152,0.686,0.13,0.013,0.004,0.003,0.002,0.002,0.005,0.017},
					 {0.016,0.143,0.693,0.147,0.016,0.006,0.003,0.002,0.003,0.005},
					 {0.005,0.013,0.138,0.667,0.16,0.018,0.005,0.002,0.002,0.003},
					 {0.003,0.004,0.013,0.142,0.63,0.156,0.015,0.004,0.002,0.002},
					 {0.002,0.002,0.004,0.015,0.158,0.629,0.142,0.013,0.004,0.003},
					 {0.003,0.002,0.002,0.005,0.018,0.164,0.667,0.138,0.013,0.005},
					 {0.005,0.003,0.002,0.003,0.006,0.016,0.148,0.692,0.141,0.016},
					 {0.017,0.004,0.002,0.002,0.003,0.005,0.013,0.131,0.686,0.152},
					 {0.144,0.014,0.004,0.002,0.002,0.002,0.004,0.012,0.129,0.654}}}};
      // correct the values
      for(unsigned int ix=0;ix<4;++ix) {
	for(unsigned int iy=0;iy<2;++iy) {
	  Scatter2DPtr corrected;
	  book(corrected,ix+1,1,iy+1);
	  // first extract values and errors applying efficiency
	  Vector<10> val,err;
	  for(unsigned int ibin=0;ibin<_h[ix][iy]->bins().size();++ibin) {
	    val[ibin] = eff[ix][iy][ibin]/100. * _h[ix][iy]->bins()[ibin].area();
	    err[ibin] = sqr(eff[ix][iy][ibin]/100. * _h[ix][iy]->bins()[ibin].areaErr());
	                sqr(efe[ix][iy][ibin]/100. * _h[ix][iy]->bins()[ibin].area());
	  }
	  // put response into a matrix
	  Matrix<10> R,R2;
	  for(unsigned int i1=0;i1<10;++i1) {
	    for(unsigned int i2=0;i2<10;++i2) {
	      R .set(i1,i2,    response[ix][iy][i1][i2]);
	      R2.set(i1,i2,sqr(response[ix][iy][i1][i2]));
	    }
	  }
	  // multiply to get value and error^2
	  val = multiply(R ,val);
	  err = multiply(R2,err);
	  // total for normalization
	  double total=0.;
	  for(unsigned int i1=0;i1<10;++i1) total+=val[i1];
	  // finally the output scatter
	  for(unsigned int ibin=0;ibin<_h[ix][iy]->bins().size();++ibin) {
	    double dx = 0.5*_h[ix][iy]->bins()[ibin].xWidth();
	    double dy = sqrt(err[ibin])/total/2./dx;
	    corrected->addPoint(_h[ix][iy]->bins()[ibin].xMid(),
				val[ibin]/total/2./dx,
				make_pair(dx,dx),
				make_pair(dy,dy));
	  }
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4][2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2019_I1693396);

}
