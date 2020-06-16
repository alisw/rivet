#include "Rivet/Tools/RivetMT2.hh"
#include "Rivet/Math/Vector4.hh"
#include "lester_mt2_bisect.hh"

namespace Rivet {

    double mT2Sq(const FourMomentum& a, const FourMomentum& b, const Vector3& ptmiss,
                 double invisiblesMass, double invisiblesMass2) {
      return asymm_mt2_lester_bisect::get_mT2_Sq(a.mass(), a.px(), a.py(),
                                                 b.mass(), b.px(), b.py(),
                                                 ptmiss.x(), ptmiss.y(),
                                                 invisiblesMass,
                                                 (invisiblesMass2 < 0 ? invisiblesMass : invisiblesMass2),
                                                 0.0, true);
    }
    
}
