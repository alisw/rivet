// -*- C++ -*-
#include "Rivet/Projections/DISRapidityGap.hh"

namespace Rivet {


  int DISRapidityGap::compare(const Projection& p) const {
    return mkNamedPCmp(p, "DISKIN") || mkNamedPCmp(p, "DISFS");
  }


  void DISRapidityGap::project(const Event& e) {
    const DISKinematics& dk =
      apply<DISKinematics>(e, "DISKIN");
    const Particles& p =
      apply<DISFinalState>(e, "DISFS").particles(cmpMomByEta);
    findgap(p, dk);
  }

  void DISRapidityGap::clearAll() {
    _M2X = _M2Y = _t = _gap = 0.;
    _gapUpp = _gapLow = -8.;
    _ePpzX_HCM = _eMpzX_HCM =_ePpzX_LAB =
      _eMpzX_LAB = _ePpzX_XCM = _eMpzX_XCM = 0.;
    _momX_HCM.setPE(0., 0., 0., 0.);
    _momY_HCM.setPE(0., 0., 0., 0.);
    _momX_XCM.setPE(0., 0., 0., 0.);
    _momY_XCM.setPE(0., 0., 0., 0.);
    _momX_LAB.setPE(0., 0., 0., 0.);
    _momY_LAB.setPE(0., 0., 0., 0.);
    _pX_HCM.clear();
    _pY_HCM.clear();
    _pX_XCM.clear();
    _pY_XCM.clear();
    _pX_LAB.clear();
    _pY_LAB.clear();
  }

  void DISRapidityGap::findgap(const Particles& particles,
                               const DISKinematics& diskin) {

    clearAll();

    // Begin by finding largest gap and gapedges between all final
    // state particles in HCM frame.
    int nP  = particles.size();
    int dir = diskin.orientation();
    for (int i = 0; i < nP-1; ++i){
      double tmpGap = abs(particles[i+1].eta() - particles[i].eta());
      if (tmpGap > _gap) {
        _gap    = tmpGap;
        _gapLow = (dir > 0) ? particles[i].eta() : dir * particles[i+1].eta();
        _gapUpp = (dir > 0) ? particles[i+1].eta() : dir * particles[i].eta();
      }
    }

    // Define the two systems X and Y.
    Particles tmp_pX, tmp_pY;
    foreach (const Particle& ip, particles) {
      if (dir * ip.eta() > _gapLow) tmp_pX.push_back(ip);
      else tmp_pY.push_back(ip);
    }

    Particles pX, pY;
    pX = (dir < 0) ? tmp_pY : tmp_pX;
    pY = (dir < 0) ? tmp_pX : tmp_pY;

    // Find variables related to HCM frame.
    // Note that HCM has photon along +z, as opposed to
    // H1 where proton is along +z. This results in a sign change
    // as compared to H1 papers!
      
    // X - side
    FourMomentum momX;
    for (const Particle& jp : pX) {
      momX  += jp.momentum();
      _ePpzX_HCM += jp.E() - jp.pz(); // Sign + => -
      _eMpzX_HCM += jp.E() + jp.pz(); // Sign - => +
    }
    _momX_HCM = momX;
    _pX_HCM   = pX;
    _M2X      = _momX_HCM.mass2();
  
    // Y - side
    FourMomentum momY;
    for (const Particle& kp : pY) momY += kp.momentum();
    _momY_HCM = momY;
    _pY_HCM   = pY;
    _M2Y      = _momY_HCM.mass2();

    // Find variables related to LAB frame
    const LorentzTransform hcmboost   = diskin.boostHCM();
    const LorentzTransform hcminverse = hcmboost.inverse();
    _momX_LAB = hcminverse.transform(_momX_HCM);
    _momY_LAB = hcminverse.transform(_momY_HCM);

    // Find momenta in XCM frame. Note that it is HCM frame that is
    // boosted, resulting in a sign change later!
    const bool doXCM = (momX.betaVec().mod2() < 1.);
    if (doXCM) {
      const LorentzTransform xcmboost =
        LorentzTransform::mkFrameTransformFromBeta(momX.betaVec());
      _momX_XCM = xcmboost.transform(momX);
      _momY_XCM = xcmboost.transform(momY);
    }

    for (const Particle& jp : pX) {
      // Boost from HCM to LAB. 
      FourMomentum lab = hcminverse.transform(jp.momentum());
      _ePpzX_LAB += lab.E() + dir * lab.pz();
      _eMpzX_LAB += lab.E() - dir * lab.pz();
      Particle plab = jp;
      plab.setMomentum(lab);
      _pX_LAB.push_back(plab);
      // Set XCM. Note that since HCM frame is boosted to XCM frame,
      // we have a sign change
      if (doXCM) {
        const LorentzTransform xcmboost =
          LorentzTransform::mkFrameTransformFromBeta(_momX_HCM.betaVec());
        FourMomentum xcm = xcmboost.transform(jp.momentum());
        _ePpzX_XCM += xcm.E() - xcm.pz(); // Sign + => -
        _eMpzX_XCM += xcm.E() + xcm.pz(); // Sign - => +
        Particle pxcm = jp;
        pxcm.setMomentum(xcm);
        _pX_XCM.push_back(pxcm);
      }
    }

    for ( const Particle& jp : pY ) {
      // Boost from HCM to LAB
      FourMomentum lab = hcminverse.transform(jp.momentum());
      Particle plab = jp;
      plab.setMomentum(lab);
      _pY_LAB.push_back(plab);
      // Boost from HCM to XCM
      if (doXCM) {
        const LorentzTransform xcmboost =
          LorentzTransform::mkFrameTransformFromBeta(_momX_HCM.betaVec());
        FourMomentum xcm = xcmboost.transform(jp.momentum());
        Particle pxcm = jp;
        pxcm.setMomentum(xcm);
        _pY_XCM.push_back(pxcm);
      }
    }

    // Find t: Currently can only handle gap on proton side.
    // @TODO: Expand to also handle gap on photon side
    // Boost p from LAB to HCM frame to find t.
    FourMomentum proton = hcmboost.transform(diskin.beamHadron().momentum());
    FourMomentum pPom = proton - _momY_HCM;
    _t = pPom * pPom;

  }
}
