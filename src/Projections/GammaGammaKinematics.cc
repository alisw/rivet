// -*- C++ -*-
#include "Rivet/Projections/GammaGammaKinematics.hh"
#include "Rivet/Math/Constants.hh"

namespace Rivet {


  void GammaGammaKinematics::project(const Event& e) {
    // Find appropriate GammaGamma leptons
    const GammaGammaLeptons& gglep = applyProjection<GammaGammaLeptons>(e, "Lepton");
    if ( gglep.failed() ) {
      fail();
      return;
    }
    _inLepton  = gglep. in();
    _outLepton = gglep.out();

    // // Get the GammaGamma lepton and store some of its properties
    // const FourMomentum pHad = _inHadron.momentum();
    const pair<FourMomentum,FourMomentum> pLepIn = make_pair(_inLepton.first.momentum(),
							     _inLepton.second.momentum());
    const pair<FourMomentum,FourMomentum> pLepOut = make_pair(_outLepton.first.momentum(),
							      _outLepton.second.momentum());
    const pair<FourMomentum,FourMomentum> pGamma  = make_pair(pLepIn.first   - pLepOut.first,
							      pLepIn.second  - pLepOut.second);
    const FourMomentum tothad = pGamma.first+pGamma.second;
    _theQ2 = make_pair(-pGamma.first.mass2(),-pGamma.second.mass2());
    _theW2 = tothad.mass2();
    // _theX = Q2()/(2.0 * pGamma * pHad);
    // _theY = (pGamma * pHad) / (pLepIn * pHad);
    // _theS = invariant(pLepIn + pHad);

    // // Calculate boost vector for boost into HCM-system
    // LorentzTransform tmp;
    // tmp.setBetaVec(-tothad.boostVector());

    // // Rotate so the photon is in x-z plane in HCM rest frame
    // FourMomentum pGammaHCM = tmp.transform(pGamma);
    // tmp.preMult(Matrix3(Vector3::mkZ(), -pGammaHCM.azimuthalAngle()));
    // pGammaHCM = tmp.transform(pGamma);
    // assert(isZero(dot(pGammaHCM.vector3(), Vector3::mkY())));

    // // Rotate so the photon is along the positive z-axis
    // const double rot_angle = pGammaHCM.polarAngle() * (pGammaHCM.px() >= 0 ? -1 : 1);
    // tmp.preMult(Matrix3(Vector3::mkY(), rot_angle));
    // // Check that final HCM photon lies along +ve z as expected
    // pGammaHCM = tmp.transform(pGamma);
    // assert(isZero(dot(pGammaHCM.vector3(), Vector3::mkX()), 1e-3));
    // assert(isZero(dot(pGammaHCM.vector3(), Vector3::mkY()), 1e-3));
    // assert(isZero(angle(pGammaHCM.vector3(), Vector3::mkZ()), 1e-3));

    // // Finally rotate so outgoing lepton at phi = 0
    // FourMomentum pLepOutHCM = tmp.transform(pLepOut);
    // tmp.preMult(Matrix3(Vector3::mkZ(), -pLepOutHCM.azimuthalAngle()));
    // assert(isZero(tmp.transform(pLepOut).azimuthalAngle()));
    // _hcm = tmp;

    // // Boost to Breit frame (use opposite convention for photon --- along *minus* z)
    // tmp.preMult(Matrix3(Vector3::mkX(), PI));
    // const double bz = 1 - 2*x();
    // _breit = LorentzTransform::mkObjTransformFromBeta(Vector3::mkZ() * bz).combine(tmp);
    // assert(isZero(angle(_breit.transform(pGamma).vector3(), -Vector3::mkZ()), 1e-3));
    // assert(isZero(_breit.transform(pLepOut).azimuthalAngle(), 1e-3));
  }


  CmpState GammaGammaKinematics::compare(const Projection & p) const {
    const GammaGammaKinematics& other = pcast<GammaGammaKinematics>(p);
    return mkNamedPCmp(other, "Lepton");
  }


}
