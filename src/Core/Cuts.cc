#include <Rivet/Cuts.hh>
#include <boost/make_shared.hpp>

// headers for converters
#include <Rivet/Particle.hh>
#include <Rivet/Jet.hh>
#include <Rivet/Math/Vectors.hh>
#include <fastjet/PseudoJet.hh>
#include <HepMC/SimpleVector.h>

// todo Sort out what can go into anonymous namespace{}

namespace Rivet {


  // Base for all wrapper classes that translate ClassToCheck to Cuttable
  class CuttableBase {
  public:
    virtual double getValue(Cuts::Quantity) const = 0;
    virtual ~CuttableBase() {}
  };


  // Cuttables can be directly passed to @ref _accept
  template <>
  bool CutBase::accept<CuttableBase>(const CuttableBase & t) const {
    return _accept(t);
  }

  // bool operator==(const Cut & a, const Cut & b) {
  //   return *a == b;
  // }


  // Open Cut singleton
  class Open_Cut : public CutBase {
  public:
    bool operator==(const Cut & c) const {
      shared_ptr<Open_Cut> cc = dynamic_pointer_cast<Open_Cut>(c);
      return bool(cc);
    }
  protected:
    // open cut accepts everything
    bool _accept(const CuttableBase &) const { return true; }
  };

  const Cut & Cuts::open() {
    // only ever need one static open cut object
    static const Cut open = boost::make_shared<Open_Cut>();
    return open;
  }





  // Cut constructor for >=
  class Cut_GtrEq : public CutBase {
  public:
    Cut_GtrEq(const Cuts::Quantity qty, const double low) : qty_(qty), low_(low) {}
    bool operator==(const Cut & c) const {
      shared_ptr<Cut_GtrEq> cc = dynamic_pointer_cast<Cut_GtrEq>(c);
      return cc && qty_ == cc->qty_  &&  low_ == cc->low_;
    }
  protected:
    bool _accept(const CuttableBase & o) const { return o.getValue(qty_) >= low_; }
  private:
    Cuts::Quantity qty_;
    double low_;
  };

  // Cut constructor for <
  class Cut_Less : public CutBase {
  public:
    Cut_Less(const Cuts::Quantity qty, const double high) : qty_(qty), high_(high) {}
    bool operator==(const Cut & c) const {
      shared_ptr<Cut_Less> cc = dynamic_pointer_cast<Cut_Less>(c);
      return cc && qty_ == cc->qty_  &&  high_ == cc->high_;
    }
  protected:
    bool _accept(const CuttableBase & o) const { return o.getValue(qty_) < high_; }
  private:
    Cuts::Quantity qty_;
    double high_;
  };

  // Cut constructor for >=
  class Cut_Gtr : public CutBase {
  public:
    Cut_Gtr(const Cuts::Quantity qty, const double low) : qty_(qty), low_(low) {}
    bool operator==(const Cut & c) const {
      shared_ptr<Cut_Gtr> cc = dynamic_pointer_cast<Cut_Gtr>(c);
      return cc && qty_ == cc->qty_  &&  low_ == cc->low_;
    }
  protected:
    bool _accept(const CuttableBase & o) const { return o.getValue(qty_) > low_; }
  private:
    Cuts::Quantity qty_;
    double low_;
  };

  // Cut constructor for <
  class Cut_LessEq : public CutBase {
  public:
    Cut_LessEq(const Cuts::Quantity qty, const double high) : qty_(qty), high_(high) {}
    bool operator==(const Cut & c) const {
      shared_ptr<Cut_LessEq> cc = dynamic_pointer_cast<Cut_LessEq>(c);
      return cc && qty_ == cc->qty_  &&  high_ == cc->high_;
    }
  protected:
    bool _accept(const CuttableBase & o) const { return o.getValue(qty_) <= high_; }
  private:
    Cuts::Quantity qty_;
    double high_;
  };


  template <typename T>
  inline Cut make_cut(T t) {
    return boost::make_shared<T>(t);
  }

  Cut operator < (Cuts::Quantity qty, double n) {
    return make_cut(Cut_Less(qty, n));
  }

  Cut operator >= (Cuts::Quantity qty, double n) {
    return make_cut(Cut_GtrEq(qty, n));
  }

  Cut operator <= (Cuts::Quantity qty, double n) {
    return make_cut(Cut_LessEq(qty, n));
  }

  Cut operator > (Cuts::Quantity qty, double n) {
    return make_cut(Cut_Gtr(qty, n));
  }

  Cut Cuts::range(Cuts::Quantity qty, double m, double n) {
    if (m > n) std::swap(m,n);
    return (qty >= m) & (qty < n);
  }


  //////////////
  /// Combiners

  /// AND, OR, NOT, and XOR objects for combining cuts

  class CutsOr : public CutBase {
  public:
    CutsOr(const Cut & c1, const Cut & c2) : cut1(c1), cut2(c2) {}
    bool operator==(const Cut & c) const {
      shared_ptr<CutsOr> cc = dynamic_pointer_cast<CutsOr>(c);
      return cc && (   ( cut1 == cc->cut1  &&  cut2 == cc->cut2 )
                       || ( cut1 == cc->cut2  &&  cut2 == cc->cut1 ));
    }
  protected:
    bool _accept(const CuttableBase & o) const {
      return cut1->accept(o) || cut2->accept(o);
    }
  private:
    const Cut cut1;
    const Cut cut2;
  };

  class CutsAnd : public CutBase {
  public:
    CutsAnd(const Cut & c1, const Cut & c2) : cut1(c1), cut2(c2) {}
    bool operator==(const Cut & c) const {
      shared_ptr<CutsAnd> cc = dynamic_pointer_cast<CutsAnd>(c);
      return cc && (   ( cut1 == cc->cut1  &&  cut2 == cc->cut2 )
                       || ( cut1 == cc->cut2  &&  cut2 == cc->cut1 ));
    }
  protected:
    bool _accept(const CuttableBase & o) const {
      return cut1->accept(o) && cut2->accept(o);
    }
  private:
    const Cut cut1;
    const Cut cut2;
  };

  class CutInvert : public CutBase {
  public:
    CutInvert(const Cut & c1) : cut(c1) {}
    bool operator==(const Cut & c) const {
      shared_ptr<CutInvert> cc = dynamic_pointer_cast<CutInvert>(c);
      return cc && cut == cc->cut;
    }
  protected:
    bool _accept(const CuttableBase & o) const {
      return !cut->accept(o);
    }
  private:
    const Cut cut;
  };

  class CutsXor : public CutBase {
  public:
    CutsXor(const Cut & c1, const Cut & c2) : cut1(c1), cut2(c2) {}
    bool operator==(const Cut & c) const {
      shared_ptr<CutsXor> cc = dynamic_pointer_cast<CutsXor>(c);
      return cc && (   ( cut1 == cc->cut1  &&  cut2 == cc->cut2 )
                       || ( cut1 == cc->cut2  &&  cut2 == cc->cut1 ));
    }
  protected:
    bool _accept(const CuttableBase & o) const {
      bool A_and_B = cut1->accept(o) && cut2->accept(o);
      bool A_or_B  = cut1->accept(o) || cut2->accept(o);
      return A_or_B && (! A_and_B);
    }
  private:
    const Cut cut1;
    const Cut cut2;
  };

  ////////////
  ///Operators

  Cut operator && (const Cut & aptr, const Cut & bptr) {
    return make_cut(CutsAnd(aptr, bptr));
  }

  Cut operator || (const Cut & aptr, const Cut & bptr) {
    return make_cut(CutsOr(aptr, bptr));
  }

  Cut operator ! (const Cut & cptr) {
    return make_cut(CutInvert(cptr));
  }


  Cut operator & (const Cut & aptr, const Cut & bptr) {
    return make_cut(CutsAnd(aptr, bptr));
  }

  Cut operator | (const Cut & aptr, const Cut & bptr) {
    return make_cut(CutsOr(aptr, bptr));
  }

  Cut operator ~ (const Cut & cptr) {
    return make_cut(CutInvert(cptr));
  }

  Cut operator ^ (const Cut & aptr, const Cut & bptr) {
    return make_cut(CutsXor(aptr, bptr));
  }

  ///////////////////////
  /// Cuts

  // Non-functional Cuttable template class. Must be specialized below.
  template <typename T> class Cuttable;

  // Non-cuttables need to be wrapped into a Cuttable first
  #define SPECIALISE_ACCEPT(TYPENAME)                           \
    template <>                                                 \
    bool CutBase::accept<TYPENAME>(const TYPENAME & t) const {  \
      return _accept(Cuttable<TYPENAME>(t));                    \
    }                                                           \


  void qty_not_found() {
    throw Exception("Missing implementation for a Cuts::Quantity.");
  }


  template<>
  class Cuttable<Particle> : public CuttableBase {
  public:
    Cuttable(const Particle& p) : p_(p) {}
    double getValue(Cuts::Quantity qty) const {
      switch ( qty ) {
      case Cuts::pT:     return p_.momentum().pT();
      case Cuts::Et:     return p_.momentum().Et();
      case Cuts::mass:   return p_.momentum().mass();
      case Cuts::rap:    return p_.momentum().rapidity();
      case Cuts::absrap: return std::abs(p_.momentum().rapidity());
      case Cuts::eta:    return p_.momentum().pseudorapidity();
      case Cuts::abseta: return std::abs(p_.momentum().pseudorapidity());
      case Cuts::phi:    return p_.momentum().phi();
      default: qty_not_found();
      }
      return -999.;
    }
  private:
    const Particle & p_;
  };
  SPECIALISE_ACCEPT(Particle)


  template<>
  class Cuttable<FourMomentum> : public CuttableBase {
  public:
    Cuttable(const FourMomentum& fm) : fm_(fm) {}
    double getValue(Cuts::Quantity qty) const {
      switch ( qty ) {
      case Cuts::pT:     return fm_.pT();
      case Cuts::Et:     return fm_.Et();
      case Cuts::mass:   return fm_.mass();
      case Cuts::rap:    return fm_.rapidity();
      case Cuts::absrap: return std::abs(fm_.rapidity());
      case Cuts::eta:    return fm_.pseudorapidity();
      case Cuts::abseta: return std::abs(fm_.pseudorapidity());
      case Cuts::phi:    return fm_.phi();
      default: qty_not_found();
      }
      return -999.;
    }
  private:
    const FourMomentum & fm_;
  };
  SPECIALISE_ACCEPT(FourMomentum)


  template<>
  class Cuttable<Jet> : public CuttableBase {
  public:
    Cuttable(const Jet& jet) : jet_(jet) {}
    double getValue(Cuts::Quantity qty) const {
      switch ( qty ) {
      case Cuts::pT:     return jet_.momentum().pT();
      case Cuts::Et:     return jet_.momentum().Et();
      case Cuts::mass:   return jet_.momentum().mass();
      case Cuts::rap:    return jet_.momentum().rapidity();
      case Cuts::absrap: return std::abs(jet_.momentum().rapidity());
      case Cuts::eta:    return jet_.momentum().pseudorapidity();
      case Cuts::abseta: return std::abs(jet_.momentum().pseudorapidity());
      case Cuts::phi:    return jet_.momentum().phi();
      default: qty_not_found();
      }
      return -999.;
    }
  private:
    const Jet & jet_;
  };
  SPECIALISE_ACCEPT(Jet)


  template<>
  class Cuttable<fastjet::PseudoJet> : public CuttableBase {
  public:
    Cuttable(const fastjet::PseudoJet& pjet) : pjet_(pjet) {}
    double getValue(Cuts::Quantity qty) const {
      switch ( qty ) {
      case Cuts::pT:     return pjet_.perp();
      case Cuts::Et:     return pjet_.Et();
      case Cuts::mass:   return pjet_.m();
      case Cuts::rap:    return pjet_.rap();
      case Cuts::absrap: return std::abs(pjet_.rap());
      case Cuts::eta:    return pjet_.eta();
      case Cuts::abseta: return std::abs(pjet_.eta());
      case Cuts::phi:    return pjet_.phi();
      default: qty_not_found();
      }
      return -999.;
    }
  private:
    const fastjet::PseudoJet & pjet_;
  };
  SPECIALISE_ACCEPT(fastjet::PseudoJet)


  template<>
  class Cuttable<HepMC::FourVector> : public CuttableBase {
  public:
    Cuttable(const HepMC::FourVector & vec) : vec_(vec) {}
    double getValue(Cuts::Quantity qty) const {
      switch ( qty ) {
      case Cuts::pT:     return vec_.perp();
        /// @todo case Cuts::Et:     return vec_.perp();
      case Cuts::mass:   return vec_.m();
      case Cuts::rap:    return 0.5*std::log((vec_.t()+vec_.z())/(vec_.t()-vec_.z()));
      case Cuts::absrap: return std::abs(getValue(Cuts::rap));
      case Cuts::eta:    return vec_.pseudoRapidity();
      case Cuts::abseta: return std::abs(vec_.pseudoRapidity());
      case Cuts::phi:    return vec_.phi();
      default: qty_not_found();
      }
      return -999.;
    }
  private:
    const HepMC::FourVector & vec_;
  };
  SPECIALISE_ACCEPT(HepMC::FourVector)


}
