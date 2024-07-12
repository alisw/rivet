// -*- C++ -*-
#include "Rivet/Math/Vector4.hh"
#include "Rivet/Tools/RivetFastJet.hh"

namespace Rivet {


  FourVector::operator fastjet::PseudoJet () const {
    return fastjet::PseudoJet(x(), y(), z(), t());
  }


}
