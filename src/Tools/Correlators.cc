// -*- C++ -*-
#include "Rivet/Tools/Correlators.hh"

namespace Rivet {


  // Constructor
  Correlators::Correlators(const ParticleFinder& fsp, int nMaxIn,
    int pMaxIn, vector<double> pTbinEdgesIn) :
      nMax(nMaxIn + 1), pMax(pMaxIn + 1), pTbinEdges(pTbinEdgesIn) {
    setName("Correlators");
    declareProjection(fsp, "FS");
    isPtDiff   = !pTbinEdges.empty();
    if (isPtDiff) {
      vector<double>::iterator underflow = pTbinEdges.begin();
      pTbinEdges.insert(underflow,pTbinEdges[0]-1);
    }
    setToZero();
  }

  // Alternative constructor.
  Correlators::Correlators(const ParticleFinder& fsp, int nMaxIn,
    int pMaxIn, const YODA::Scatter2D hIn) : nMax(nMaxIn + 1), pMax(pMaxIn + 1) {
    for (auto b : hIn.points()) pTbinEdges.push_back(b.xMin());
    pTbinEdges.push_back(hIn.points().back().xMax());
    setName("Correlators");
    declareProjection(fsp, "FS");
    isPtDiff   = !pTbinEdges.empty();
    if (isPtDiff) {
      vector<double>::iterator underflow = pTbinEdges.begin();
      pTbinEdges.insert(underflow,pTbinEdges[0]-1);
    }
    setToZero();
  }


  // Set all elements in vectors to zero
  void Correlators::setToZero(){
    vector< complex<double> > pTmp(pMax, _ZERO);
    Vec2D qTmp(nMax, pTmp);
    qVec = qTmp;
    if (isPtDiff) {
      pVec.erase(pVec.begin(), pVec.end());
      for (double pT : pTbinEdges)
         pVec.insert(pair<double, Vec2D>(pT, qVec));
    }
  }

  // Functions for output:
  const pair<double,double> Correlators::intCorrelator(vector<int> n) const {
    // Create vector of zeros for normalisation and vector of initial
    // powers
    int m = n.size();
    vector<int> powers(m, 1);
    vector<int> zeros(m, 0);
    complex<double> num = recCorr(m, n, powers, false);
    complex<double> den = recCorr(m, zeros, powers, false);
    pair<double, double> ret;
    ret.second = (den.real() < _TINY) ? 0. : den.real();
    ret.first = num.real();
    return ret;
  }

  const vector<pair<double,double>> Correlators::pTBinnedCorrelators(vector<int> n,
    bool overflow) const {
    // Create vector of zeros for normalisation and vector of initial
    // powers
    if (!isPtDiff)
      cout << "You must book the correlator with a binning if you want to"
	      " extract binned correlators! Failing." << endl;
    int m = n.size();
    vector<int> powers(m, 1);
    vector<int> zeros(m, 0);
    vector<pair<double,double>> ret;
    for (double pT : pTbinEdges){
      complex<double> num = recCorr(m, n, powers, true, pT);
      complex<double> den = recCorr(m, zeros, powers, true, pT);
      pair<double, double> tmp;
      tmp.second = (den.real() < _TINY) ? 0. : den.real();
      tmp.first = num.real();
      ret.push_back(tmp);
    }
    if (!overflow)
      return vector<pair<double, double> > (ret.begin() + 1, ret.end() - 1);
    return ret;
  }

  // M-particle correlation with eta-gap
  const pair<double,double> Correlators::intCorrelatorGap(const Correlators& other,
    vector<int> n1, vector<int> n2) const {
    // Create vectors of zeros for normalisation and vectors of initial
    // powers
    int m1 = n1.size();
    int m2 = n2.size();
    // Return if too few particles in event.
    vector<int> zero1(m1, 0);
    vector<int> zero2(m2, 0);
    vector<int> p1(m1, 1);
    vector<int> p2(m2, 1);
    complex<double> num1 = recCorr(m1, n1, p1, false);
    complex<double> den1 = recCorr(m1, zero1, p1, false);
    complex<double> num2 = other.recCorr(m2, n2, p2, false);
    complex<double> den2 = other.recCorr(m2, zero2, p2, false);
    complex<double> num  = num1 * num2;
    complex<double> den  = den1 * den2;
    pair<double, double> ret;
    ret.second = (den1.real() < _TINY || den2.real() < _TINY)  ? 0. :
	   den.real();
    ret.first = num.real();
    return ret;
  }

  // M-particle correlation with eta-gap
  const vector<pair<double,double>> Correlators::pTBinnedCorrelatorsGap(
    const Correlators& other, vector<int> n1, vector<int> n2,
    bool overflow) const {
    if (!isPtDiff)
      cout << "You must book the correlator with a binning if you want to"
	      " extract binned correlators! Failing." << endl;
    // Create vectors of zeros for normalisation and vectors of initial
    // powers
    int m1 = n1.size();
    int m2 = n2.size();
    vector<int> zero1(m1, 0);
    vector<int> zero2(m2, 0);
    vector<int> p1(m1, 1);
    vector<int> p2(m2, 1);
    vector<pair<double,double>> ret;
    for (double pT : pTbinEdges) {
      complex<double> num1 = recCorr(m1, n1, p1, true, pT);
      complex<double> den1 = recCorr(m1, zero1, p1, true, pT);
      complex<double> num2 = other.recCorr(m2, n2, p2, false);
      complex<double> den2 = other.recCorr(m2, zero2, p2, false);
      complex<double> num  = num1 * num2;
      complex<double> den  = den1 * den2;
      pair<double, double> tmp;
      tmp.second = (den1.real() < _TINY || den2.real() < _TINY)
        ? 0. : den.real();
      tmp.first = num.real();
      ret.push_back(tmp);
    }
    //If we don't want to include underflow/overflow, remove them here
    if (!overflow)
      return vector<pair<double, double> > (ret.begin() + 1, ret.end() - 1);
    return ret;
  }


  // Project function. Loops over array and calculates Q vectors
  void Correlators::project(const Event& e) {
    setToZero();
    // @TODO: Weight could be implemented to account for non-uniform
    // acceptance if detector simulation is needed. If not, the weight
    // can be unity. Note that this weight is not the MC event weight, which
    // should be used when filling histograms (as usual).
    const double w = 1.0;;
    const Particles& parts = applyProjection<ParticleFinder>(e, "FS").particles();
    // Check that we have at least two particles in the event
    if (parts.size() > 2)  {
      for(const Particle& p : parts)
        fillCorrelators(p, w);
    }
  }

  // Calculate correlators from one particle
  void Correlators::fillCorrelators(const Particle& p, const double& weight = 1.) {
    for (int iN = 0; iN < nMax; ++iN)
    for (int iP = 0; iP < pMax; ++iP) {
      double real = cos(iN * p.phi());
      double imag = sin(iN * p.phi());
      complex<double> expi(real, imag);
      complex<double> tmp = pow(weight, iP) * expi;
      qVec[iN][iP] += tmp;
      if (isPtDiff) {
        map<double, Vec2D>::iterator pTitr = pVec.lower_bound(p.pT());
        // Move to the correct bin.
        if (pTitr != pVec.begin()) pTitr--;
        pTitr->second[iN][iP] += tmp;
      }
    }
  }

  // Two-particle correlator Eq. (19) p. 6 in Generic Fr. paper.
  const complex<double> Correlators::twoPartCorr(int n1, int n2, int p1,
    int p2, double pT, bool useP) const {
    complex<double> tmp1 = (!useP) ? getQ(n1, p1) :
                           getP(n1, p1, pT);
    complex<double> tmp2 = getQ(n2, p2);
    complex<double> tmp3 = (!useP) ? getQ(n1+n2, p1+p2) :
                           getP(n1+n2, p1+p2, pT);
    complex<double> sum  = tmp1 * tmp2 - tmp3;
    return sum;
  }

  // Find correlators by recursion. Order = M (# of particles),
  // n's are harmonics, p's are the powers of the weights
  const complex<double> Correlators::recCorr(int order, vector<int> n,
    vector<int> p, bool useP, double pT) const {
    // Sanity checks
    int nUsed = 0;
    for (int i = 0, N = n.size(); i < N; ++i) nUsed += n[i];
    if (nMax < nUsed)
      cout <<"Requested n = " << nUsed << ", nMax = " << nMax << endl;
    if (int(p.size()) > pMax)
      cout << "Requested p = " << p.size() << ", pMax = " << pMax << endl;
    // If order is 1, then return Q/p vector (important when dealing
    // with gaps and one side has only one particle
    if (order < 2)
      return (!useP) ? getQ(n[0], p[0]) : getP(n[0], p[0], pT);
    // Return 2-p correlator explicitly.
    if ( order < 3 )
      return twoPartCorr(n[0], n[1], p[0], p[1], pT, useP);

    // Else find nth order harmonics by recursion
    // at order M - 1
    int orderNow           = order - 1;
    int nNow               = n[orderNow];
    int pNow               = p[orderNow];
    complex<double> recNow = getQ(n[orderNow], p[orderNow]);
    recNow                *= recCorr(orderNow, n, p, useP, pT);
    for (int i = 0; i < orderNow; ++i){
      vector<int> tmpN, tmpP;
      for (int j = 0; j < orderNow; ++j){
        tmpN.push_back(n[j]);
        tmpP.push_back(p[j]);
      }
      tmpN[i] += nNow;
      tmpP[i] += pNow;
      recNow -= recCorr(orderNow, tmpN, tmpP, useP, pT);
    }
    return recNow;
  }

}
