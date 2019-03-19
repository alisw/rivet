#ifndef RIVET_MATH_VECTOR4
#define RIVET_MATH_VECTOR4

#include "Rivet/Math/MathHeader.hh"
#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Math/VectorN.hh"
#include "Rivet/Math/Vector3.hh"

namespace Rivet {


  class FourVector;
  class FourMomentum;
  class LorentzTransform;
  typedef FourVector Vector4;
  FourVector transform(const LorentzTransform& lt, const FourVector& v4);


  /// @brief Specialisation of VectorN to a general (non-momentum) Lorentz 4-vector.
  ///
  /// @todo Add composite set/mk methods from different coord systems
  class FourVector : public Vector<4> {
    friend FourVector multiply(const double a, const FourVector& v);
    friend FourVector multiply(const FourVector& v, const double a);
    friend FourVector add(const FourVector& a, const FourVector& b);
    friend FourVector transform(const LorentzTransform& lt, const FourVector& v4);

  public:

    FourVector() : Vector<4>() { }

    template<typename V4>
    FourVector(const V4& other) {
      this->setT(other.t());
      this->setX(other.x());
      this->setY(other.y());
      this->setZ(other.z());
    }

    FourVector(const Vector<4>& other)
      : Vector<4>(other) { }

    FourVector(const double t, const double x, const double y, const double z) {
      this->setT(t);
      this->setX(x);
      this->setY(y);
      this->setZ(z);
    }

    virtual ~FourVector() { }

  public:

    double t() const { return get(0); }
    double t2() const { return sqr(t()); }
    FourVector& setT(const double t) { set(0, t); return *this; }

    double x() const { return get(1); }
    double x2() const { return sqr(x()); }
    FourVector& setX(const double x) { set(1, x); return *this; }

    double y() const { return get(2); }
    double y2() const { return sqr(y()); }
    FourVector& setY(const double y) { set(2, y); return *this; }

    double z() const { return get(3); }
    double z2() const { return sqr(z()); }
    FourVector& setZ(const double z) { set(3, z); return *this; }

    double invariant() const {
      // Done this way for numerical precision
      return (t() + z())*(t() - z()) - x()*x() - y()*y();
    }

    bool isNull() const {
      return Rivet::isZero(invariant());
    }

    /// Angle between this vector and another
    double angle(const FourVector& v) const {
      return vector3().angle( v.vector3() );
    }
    /// Angle between this vector and another (3-vector)
    double angle(const Vector3& v3) const {
      return vector3().angle(v3);
    }

    /// @brief Mod-square of the projection of the 3-vector on to the \f$ x-y \f$ plane
    /// This is a more efficient function than @c polarRadius, as it avoids the square root.
    /// Use it if you only need the squared value, or e.g. an ordering by magnitude.
    double polarRadius2() const {
      return vector3().polarRadius2();
    }
    /// Synonym for polarRadius2
    double perp2() const {
      return vector3().perp2();
    }
    /// Synonym for polarRadius2
    double rho2() const {
      return vector3().rho2();
    }

    /// Magnitude of projection of 3-vector on to the \f$ x-y \f$ plane
    double polarRadius() const {
      return vector3().polarRadius();
    }
    /// Synonym for polarRadius
    double perp() const {
      return vector3().perp();
    }
    /// Synonym for polarRadius
    double rho() const {
      return vector3().rho();
    }

    /// Projection of 3-vector on to the \f$ x-y \f$ plane
    Vector3 polarVec() const {
      return vector3().polarVec();
    }
    /// Synonym for polarVec
    Vector3 perpVec() const {
      return vector3().perpVec();
    }
    /// Synonym for polarVec
    Vector3 rhoVec() const {
      return vector3().rhoVec();
    }

    /// Angle subtended by the 3-vector's projection in x-y and the x-axis.
    double azimuthalAngle(const PhiMapping mapping=ZERO_2PI) const {
      return vector3().azimuthalAngle(mapping);
    }
    /// Synonym for azimuthalAngle.
    double phi(const PhiMapping mapping=ZERO_2PI) const {
      return vector3().phi(mapping);
    }

    /// Angle subtended by the 3-vector and the z-axis.
    double polarAngle() const {
      return vector3().polarAngle();
    }
    /// Synonym for polarAngle.
    double theta() const {
      return vector3().theta();
    }

    /// Pseudorapidity (defined purely by the 3-vector components)
    double pseudorapidity() const {
      return vector3().pseudorapidity();
    }
    /// Synonym for pseudorapidity.
    double eta() const {
      return vector3().eta();
    }

    /// Get the \f$ |\eta| \f$ directly.
    double abspseudorapidity() const { return fabs(eta()); }
    /// Get the \f$ |\eta| \f$ directly (alias).
    double abseta() const { return fabs(eta()); }

    /// Get the spatial part of the 4-vector as a 3-vector.
    Vector3 vector3() const {
      return Vector3(get(1), get(2), get(3));
    }

    /// Implicit cast to a 3-vector
    operator Vector3 () const { return vector3(); }


  public:

    /// Contract two 4-vectors, with metric signature (+ - - -).
    double contract(const FourVector& v) const {
      const double result = t()*v.t() - x()*v.x() - y()*v.y() - z()*v.z();
      return result;
    }

    /// Contract two 4-vectors, with metric signature (+ - - -).
    double dot(const FourVector& v) const {
      return contract(v);
    }

    /// Contract two 4-vectors, with metric signature (+ - - -).
    double operator*(const FourVector& v) const {
      return contract(v);
    }

    /// Multiply by a scalar.
    FourVector& operator*=(double a) {
      _vec = multiply(a, *this)._vec;
      return *this;
    }

    /// Divide by a scalar.
    FourVector& operator/=(double a) {
      _vec = multiply(1.0/a, *this)._vec;
      return *this;
    }

    /// Add to this 4-vector.
    FourVector& operator+=(const FourVector& v) {
      _vec = add(*this, v)._vec;
      return *this;
    }

    /// Subtract from this 4-vector. NB time as well as space components are subtracted.
    FourVector& operator-=(const FourVector& v) {
      _vec = add(*this, -v)._vec;
      return *this;
    }

    /// Multiply all components (space and time) by -1.
    FourVector operator-() const {
      FourVector result;
      result._vec = -_vec;
      return result;
    }

    /// Multiply space components only by -1.
    FourVector reverse() const {
      FourVector result = -*this;
      result.setT(-result.t());
      return result;
    }

  };


  /// Contract two 4-vectors, with metric signature (+ - - -).
  inline double contract(const FourVector& a, const FourVector& b) {
    return a.contract(b);
  }

  /// Contract two 4-vectors, with metric signature (+ - - -).
  inline double dot(const FourVector& a, const FourVector& b) {
    return contract(a, b);
  }

  inline FourVector multiply(const double a, const FourVector& v) {
    FourVector result;
    result._vec = a * v._vec;
    return result;
  }

  inline FourVector multiply(const FourVector& v, const double a) {
    return multiply(a, v);
  }

  inline FourVector operator*(const double a, const FourVector& v) {
    return multiply(a, v);
  }

  inline FourVector operator*(const FourVector& v, const double a) {
    return multiply(a, v);
  }

  inline FourVector operator/(const FourVector& v, const double a) {
    return multiply(1.0/a, v);
  }

  inline FourVector add(const FourVector& a, const FourVector& b) {
    FourVector result;
    result._vec = a._vec + b._vec;
    return result;
  }

  inline FourVector operator+(const FourVector& a, const FourVector& b) {
    return add(a, b);
  }

  inline FourVector operator-(const FourVector& a, const FourVector& b) {
    return add(a, -b);
  }

  /// Calculate the Lorentz self-invariant of a 4-vector.
  /// \f$ v_\mu v^\mu = g_{\mu\nu} x^\mu x^\nu \f$.
  inline double invariant(const FourVector& lv) {
    return lv.invariant();
  }

  /// Angle (in radians) between spatial parts of two Lorentz vectors.
  inline double angle(const FourVector& a, const FourVector& b) {
    return a.angle(b);
  }

  /// Angle (in radians) between spatial parts of two Lorentz vectors.
  inline double angle(const Vector3& a, const FourVector& b) {
    return angle( a, b.vector3() );
  }

  /// Angle (in radians) between spatial parts of two Lorentz vectors.
  inline double angle(const FourVector& a, const Vector3& b) {
    return a.angle(b);
  }


  ////////////////////////////////////////////////


  /// Specialized version of the FourVector with momentum/energy functionality.
  class FourMomentum : public FourVector {
    friend FourMomentum multiply(const double a, const FourMomentum& v);
    friend FourMomentum multiply(const FourMomentum& v, const double a);
    friend FourMomentum add(const FourMomentum& a, const FourMomentum& b);
    friend FourMomentum transform(const LorentzTransform& lt, const FourMomentum& v4);

  public:
    FourMomentum() { }

    template<typename V4>
    FourMomentum(const V4& other) {
      this->setE(other.t());
      this->setPx(other.x());
      this->setPy(other.y());
      this->setPz(other.z());
    }

    FourMomentum(const Vector<4>& other)
      : FourVector(other) { }

    FourMomentum(const double E, const double px, const double py, const double pz) {
      this->setE(E);
      this->setPx(px);
      this->setPy(py);
      this->setPz(pz);
    }

    ~FourMomentum() {}

  public:


    /// @name Coordinate setters
    //@{

    /// Set energy \f$ E \f$ (time component of momentum).
    FourMomentum& setE(double E) {
      setT(E);
      return *this;
    }

    /// Set x-component of momentum \f$ p_x \f$.
    FourMomentum& setPx(double px) {
      setX(px);
      return *this;
    }

    /// Set y-component of momentum \f$ p_y \f$.
    FourMomentum& setPy(double py) {
      setY(py);
      return *this;
    }

    /// Set z-component of momentum \f$ p_z \f$.
    FourMomentum& setPz(double pz) {
      setZ(pz);
      return *this;
    }


    /// Set the p coordinates and energy simultaneously
    FourMomentum& setPE(double px, double py, double pz, double E) {
      if (E < 0)
        throw std::invalid_argument("Negative energy given as argument: " + to_str(E));
      setPx(px); setPy(py); setPz(pz); setE(E);
      return *this;
    }
    /// Alias for setPE
    FourMomentum& setXYZE(double px, double py, double pz, double E) {
      return setPE(px, py, pz, E);
    }
    // /// Near-alias with switched arg order
    // FourMomentum& setEP(double E, double px, double py, double pz) {
    //   return setPE(px, py, pz, E);
    // }
    // /// Alias for setEP
    // FourMomentum& setEXYZ(double E, double px, double py, double pz) {
    //   return setEP(E, px, py, pz);
    // }


    /// Set the p coordinates and mass simultaneously
    FourMomentum& setPM(double px, double py, double pz, double mass) {
      if (mass < 0)
        throw std::invalid_argument("Negative mass given as argument: " + to_str(mass));
      const double E = sqrt( sqr(mass) + sqr(px) + sqr(py) + sqr(pz) );
      // setPx(px); setPy(py); setPz(pz); setE(E);
      return setPE(px, py, pz, E);
    }
    /// Alias for setPM
    FourMomentum& setXYZM(double px, double py, double pz, double mass) {
      return setPM(px, py, pz, mass);
    }


    /// Set the vector state from (eta,phi,energy) coordinates and the mass
    ///
    /// eta = -ln(tan(theta/2))
    /// -> theta = 2 atan(exp(-eta))
    FourMomentum& setEtaPhiME(double eta, double phi, double mass, double E) {
      if (mass < 0)
        throw std::invalid_argument("Negative mass given as argument");
      if (E < 0)
        throw std::invalid_argument("Negative energy given as argument");
      const double theta = 2 * atan(exp(-eta));
      if (theta < 0 || theta > M_PI)
        throw std::domain_error("Polar angle outside 0..pi in calculation");
      setThetaPhiME(theta, phi, mass, E);
      return *this;
    }

    /// Set the vector state from (eta,phi,pT) coordinates and the mass
    ///
    /// eta = -ln(tan(theta/2))
    /// -> theta = 2 atan(exp(-eta))
    FourMomentum& setEtaPhiMPt(double eta, double phi, double mass, double pt) {
      if (mass < 0)
        throw std::invalid_argument("Negative mass given as argument");
      if (pt < 0)
        throw std::invalid_argument("Negative transverse momentum given as argument");
      const double theta = 2 * atan(exp(-eta));
      if (theta < 0 || theta > M_PI)
        throw std::domain_error("Polar angle outside 0..pi in calculation");
      const double p = pt / sin(theta);
      const double E = sqrt( sqr(p) + sqr(mass) );
      setThetaPhiME(theta, phi, mass, E);
      return *this;
    }

    /// Set the vector state from (y,phi,energy) coordinates and the mass
    ///
    /// y = 0.5 * ln((E+pz)/(E-pz))
    /// -> (E^2 - pz^2) exp(2y) = (E+pz)^2
    ///  & (E^2 - pz^2) exp(-2y) = (E-pz)^2
    /// -> E = sqrt(pt^2 + m^2) cosh(y)
    /// -> pz = sqrt(pt^2 + m^2) sinh(y)
    /// -> sqrt(pt^2 + m^2) = E / cosh(y)
    FourMomentum& setRapPhiME(double y, double phi, double mass, double E) {
      if (mass < 0)
        throw std::invalid_argument("Negative mass given as argument");
      if (E < 0)
        throw std::invalid_argument("Negative energy given as argument");
      const double sqrt_pt2_m2 = E / cosh(y);
      const double pt = sqrt( sqr(sqrt_pt2_m2) - sqr(mass) );
      if (pt < 0)
        throw std::domain_error("Negative transverse momentum in calculation");
      const double pz = sqrt_pt2_m2 * sinh(y);
      const double px = pt * cos(phi);
      const double py = pt * sin(phi);
      setPE(px, py, pz, E);
      return *this;
    }

    /// Set the vector state from (y,phi,pT) coordinates and the mass
    ///
    /// y = 0.5 * ln((E+pz)/(E-pz))
    /// -> E = sqrt(pt^2 + m^2) cosh(y)  [see above]
    FourMomentum& setRapPhiMPt(double y, double phi, double mass, double pt) {
      if (mass < 0)
        throw std::invalid_argument("Negative mass given as argument");
      if (pt < 0)
        throw std::invalid_argument("Negative transverse mass given as argument");
      const double E = sqrt( sqr(pt) + sqr(mass) ) * cosh(y);
      if (E < 0)
        throw std::domain_error("Negative energy in calculation");
      setRapPhiME(y, phi, mass, E);
      return *this;
    }

    /// Set the vector state from (theta,phi,energy) coordinates and the mass
    ///
    /// p = sqrt(E^2 - mass^2)
    /// pz = p cos(theta)
    /// pt = p sin(theta)
    FourMomentum& setThetaPhiME(double theta, double phi, double mass, double E) {
      if (theta < 0 || theta > M_PI)
        throw std::invalid_argument("Polar angle outside 0..pi given as argument");
      if (mass < 0)
        throw std::invalid_argument("Negative mass given as argument");
      if (E < 0)
        throw std::invalid_argument("Negative energy given as argument");
      const double p = sqrt( sqr(E) - sqr(mass) );
      const double pz = p * cos(theta);
      const double pt = p * sin(theta);
      if (pt < 0)
        throw std::invalid_argument("Negative transverse momentum in calculation");
      const double px = pt * cos(phi);
      const double py = pt * sin(phi);
      setPE(px, py, pz, E);
      return *this;
    }

    /// Set the vector state from (theta,phi,pT) coordinates and the mass
    ///
    /// p = pt / sin(theta)
    /// pz = p cos(theta)
    /// E = sqrt(p^2 + mass^2)
    FourMomentum& setThetaPhiMPt(double theta, double phi, double mass, double pt) {
      if (theta < 0 || theta > M_PI)
        throw std::invalid_argument("Polar angle outside 0..pi given as argument");
      if (mass < 0)
        throw std::invalid_argument("Negative mass given as argument");
      if (pt < 0)
        throw std::invalid_argument("Negative transverse momentum given as argument");
      const double p = pt / sin(theta);
      const double px = pt * cos(phi);
      const double py = pt * sin(phi);
      const double pz = p * cos(theta);
      const double E = sqrt( sqr(p) + sqr(mass) );
      setPE(px, py, pz, E);
      return *this;
    }

    /// Set the vector state from (pT,phi,energy) coordinates and the mass
    ///
    /// pz = sqrt(E^2 - mass^2 - pt^2)
    FourMomentum& setPtPhiME(double pt, double phi, double mass, double E) {
      if (pt < 0)
        throw std::invalid_argument("Negative transverse momentum given as argument");
      if (mass < 0)
        throw std::invalid_argument("Negative mass given as argument");
      if (E < 0)
        throw std::invalid_argument("Negative energy given as argument");
      const double px = pt * cos(phi);
      const double py = pt * sin(phi);
      const double pz = sqrt(sqr(E) - sqr(mass) - sqr(pt));
      setPE(px, py, pz, E);
      return *this;
    }

    //@}


    /// @name Accessors
    //@{

    /// Get energy \f$ E \f$ (time component of momentum).
    double E() const { return t(); }
    /// Get energy-squared \f$ E^2 \f$.
    double E2() const { return t2(); }

    /// Get x-component of momentum \f$ p_x \f$.
    double px() const { return x(); }
    /// Get x-squared \f$ p_x^2 \f$.
    double px2() const { return x2(); }

    /// Get y-component of momentum \f$ p_y \f$.
    double py() const { return y(); }
    /// Get y-squared \f$ p_y^2 \f$.
    double py2() const { return y2(); }

    /// Get z-component of momentum \f$ p_z \f$.
    double pz() const { return z(); }
    /// Get z-squared \f$ p_z^2 \f$.
    double pz2() const { return z2(); }


    /// @brief Get the mass \f$ m = \sqrt{E^2 - p^2} \f$ (the Lorentz self-invariant).
    ///
    /// For spacelike momenta, the mass will be -sqrt(|mass2|).
    double mass() const {
      // assert(Rivet::isZero(mass2()) || mass2() > 0);
      // if (Rivet::isZero(mass2())) {
      //   return 0.0;
      // } else {
      //   return sqrt(mass2());
      // }
      return sign(mass2()) * sqrt(fabs(mass2()));
    }

    /// Get the squared mass \f$ m^2 = E^2 - p^2 \f$ (the Lorentz self-invariant).
    double mass2() const {
      return invariant();
    }


    /// Get 3-momentum part, \f$ p \f$.
    Vector3 p3() const { return vector3(); }

    /// Get the modulus of the 3-momentum
    double p() const {
      return p3().mod();
    }

    /// Get the modulus-squared of the 3-momentum
    double p2() const {
      return p3().mod2();
    }


    /// Calculate the rapidity.
    double rapidity() const {
      return 0.5 * std::log( (E() + pz()) / (E() - pz()) );
    }
    /// Alias for rapidity.
    double rap() const {
      return rapidity();
    }

    /// Absolute rapidity.
    double absrapidity() const {
      return fabs(rapidity());
    }
    /// Absolute rapidity.
    double absrap() const {
      return fabs(rap());
    }

    /// Calculate the transverse momentum vector \f$ \vec{p}_T \f$.
    Vector3 pTvec() const {
      return p3().polarVec();
    }
    /// Synonym for pTvec
    Vector3 ptvec() const {
      return pTvec();
    }

    /// Calculate the squared transverse momentum \f$ p_T^2 \f$.
    double pT2() const {
      return vector3().polarRadius2();
    }
    /// Calculate the squared transverse momentum \f$ p_T^2 \f$.
    double pt2() const {
      return vector3().polarRadius2();
    }

    /// Calculate the transverse momentum \f$ p_T \f$.
    double pT() const {
      return sqrt(pT2());
    }
    /// Calculate the transverse momentum \f$ p_T \f$.
    double pt() const {
      return sqrt(pT2());
    }

    /// Calculate the transverse energy \f$ E_T^2 = E^2 \sin^2{\theta} \f$.
    double Et2() const {
      return Et() * Et();
    }
    /// Calculate the transverse energy \f$ E_T = E \sin{\theta} \f$.
    double Et() const {
      return E() * sin(polarAngle());
    }

    //@}


    /// @name Lorentz boost factors and vectors
    //@{

    /// Calculate the boost factor \f$ \gamma \f$.
    /// @note \f$ \gamma = E/mc^2 \f$ so we rely on the c=1 convention
    double gamma() const {
      return sqrt(E2()/mass2());
    }

    /// Calculate the boost vector \f$ \vec{\gamma} \f$.
    /// @note \f$ \gamma = E/mc^2 \f$ so we rely on the c=1 convention
    Vector3 gammaVec() const {
      return gamma() * p3().unit();
    }

    /// Calculate the boost factor \f$ \beta \f$.
    /// @note \f$ \beta = pc/E \f$ so we rely on the c=1 convention
    double beta() const {
      return p()/E();
    }

    /// Calculate the boost vector \f$ \vec{\beta} \f$.
    /// @note \f$ \beta = pc/E \f$ so we rely on the c=1 convention
    Vector3 betaVec() const {
      // return Vector3(px()/E(), py()/E(), pz()/E());
      return p3()/E();
    }

    /// @brief Deprecated alias for betaVec
    /// @deprecated This will be removed; use betaVec() instead
    Vector3 boostVector() const { return betaVec(); }

    //@}


    ////////////////////////////////////////


    /// @name Sorting helpers
    //@{

    /// Struct for sorting by increasing energy
    struct byEAscending {
      bool operator()(const FourMomentum& left, const FourMomentum& right) const{
        const double pt2left = left.E();
        const double pt2right = right.E();
        return pt2left < pt2right;
      }

      bool operator()(const FourMomentum* left, const FourMomentum* right) const{
        return (*this)(*left, *right);
      }
    };


    /// Struct for sorting by decreasing energy
    struct byEDescending {
      bool operator()(const FourMomentum& left, const FourMomentum& right) const{
        return byEAscending()(right, left);
      }

      bool operator()(const FourMomentum* left, const FourVector* right) const{
        return (*this)(*left, *right);
      }
    };

    //@}


    ////////////////////////////////////////


    /// @name Arithmetic operators
    //@{

    /// Multiply by a scalar
    FourMomentum& operator*=(double a) {
      _vec = multiply(a, *this)._vec;
      return *this;
    }

    /// Divide by a scalar
    FourMomentum& operator/=(double a) {
      _vec = multiply(1.0/a, *this)._vec;
      return *this;
    }

    /// Add to this 4-vector. NB time as well as space components are added.
    FourMomentum& operator+=(const FourMomentum& v) {
      _vec = add(*this, v)._vec;
      return *this;
    }

    /// Subtract from this 4-vector. NB time as well as space components are subtracted.
    FourMomentum& operator-=(const FourMomentum& v) {
      _vec = add(*this, -v)._vec;
      return *this;
    }

    /// Multiply all components (time and space) by -1.
    FourMomentum operator-() const {
      FourMomentum result;
      result._vec = -_vec;
      return result;
    }

    /// Multiply space components only by -1.
    FourMomentum reverse() const {
      FourMomentum result = -*this;
      result.setE(-result.E());
      return result;
    }

    //@}


    ////////////////////////////////////////


    /// @name Factory functions
    //@{

    /// Make a vector from (px,py,pz,E) coordinates
    static FourMomentum mkXYZE(double px, double py, double pz, double E) {
      return FourMomentum().setPE(px, py, pz, E);
    }

    /// Make a vector from (px,py,pz) coordinates and the mass
    static FourMomentum mkXYZM(double px, double py, double pz, double mass) {
      return FourMomentum().setPM(px, py, pz, mass);
    }

    /// Make a vector from (eta,phi,energy) coordinates and the mass
    static FourMomentum mkEtaPhiME(double eta, double phi, double mass, double E) {
      return FourMomentum().setEtaPhiME(eta, phi, mass, E);
    }

    /// Make a vector from (eta,phi,pT) coordinates and the mass
    static FourMomentum mkEtaPhiMPt(double eta, double phi, double mass, double pt) {
      return FourMomentum().setEtaPhiMPt(eta, phi, mass, pt);
    }

    /// Make a vector from (y,phi,energy) coordinates and the mass
    static FourMomentum mkRapPhiME(double y, double phi, double mass, double E) {
      return FourMomentum().setRapPhiME(y, phi, mass, E);
    }

    /// Make a vector from (y,phi,pT) coordinates and the mass
    static FourMomentum mkRapPhiMPt(double y, double phi, double mass, double pt) {
      return FourMomentum().setRapPhiMPt(y, phi, mass, pt);
    }

    /// Make a vector from (theta,phi,energy) coordinates and the mass
    static FourMomentum mkThetaPhiME(double theta, double phi, double mass, double E) {
      return FourMomentum().setThetaPhiME(theta, phi, mass, E);
    }

    /// Make a vector from (theta,phi,pT) coordinates and the mass
    static FourMomentum mkThetaPhiMPt(double theta, double phi, double mass, double pt) {
      return FourMomentum().setThetaPhiMPt(theta, phi, mass, pt);
    }

    /// Make a vector from (pT,phi,energy) coordinates and the mass
    static FourMomentum mkPtPhiME(double pt, double phi, double mass, double E) {
      return FourMomentum().setPtPhiME(pt, phi, mass, E);
    }

    //@}


  };



  inline FourMomentum multiply(const double a, const FourMomentum& v) {
    FourMomentum result;
    result._vec = a * v._vec;
    return result;
  }

  inline FourMomentum multiply(const FourMomentum& v, const double a) {
    return multiply(a, v);
  }

  inline FourMomentum operator*(const double a, const FourMomentum& v) {
    return multiply(a, v);
  }

  inline FourMomentum operator*(const FourMomentum& v, const double a) {
    return multiply(a, v);
  }

  inline FourMomentum operator/(const FourMomentum& v, const double a) {
    return multiply(1.0/a, v);
  }

  inline FourMomentum add(const FourMomentum& a, const FourMomentum& b) {
    FourMomentum result;
    result._vec = a._vec + b._vec;
    return result;
  }

  inline FourMomentum operator+(const FourMomentum& a, const FourMomentum& b) {
    return add(a, b);
  }

  inline FourMomentum operator-(const FourMomentum& a, const FourMomentum& b) {
    return add(a, -b);
  }


  //////////////////////////////////////////////////////


  /// @name \f$ \Delta R \f$ calculations from 4-vectors
  //@{

  /// @brief Calculate the squared 2D rapidity-azimuthal ("eta-phi") distance between two four-vectors.
  ///
  /// There is a scheme ambiguity for momentum-type four vectors as to whether
  /// the pseudorapidity (a purely geometric concept) or the rapidity (a
  /// relativistic energy-momentum quantity) is to be used: this can be chosen
  /// via the optional scheme parameter. Use of this scheme option is
  /// discouraged in this case since @c RAPIDITY is only a valid option for
  /// vectors whose type is really the FourMomentum derived class.
  inline double deltaR2(const FourVector& a, const FourVector& b,
                       RapScheme scheme=PSEUDORAPIDITY) {
    switch (scheme) {
    case PSEUDORAPIDITY :
      return deltaR2(a.vector3(), b.vector3());
    case RAPIDITY:
      {
        const FourMomentum* ma = dynamic_cast<const FourMomentum*>(&a);
        const FourMomentum* mb = dynamic_cast<const FourMomentum*>(&b);
        if (!ma || !mb) {
          string err = "deltaR with scheme RAPIDITY can only be called with FourMomentum objects, not FourVectors";
          throw std::runtime_error(err);
        }
        return deltaR2(*ma, *mb, scheme);
      }
    default:
      throw std::runtime_error("The specified deltaR scheme is not yet implemented");
    }
  }

  /// @brief Calculate the 2D rapidity-azimuthal ("eta-phi") distance between two four-vectors.
  ///
  /// There is a scheme ambiguity for momentum-type four vectors as to whether
  /// the pseudorapidity (a purely geometric concept) or the rapidity (a
  /// relativistic energy-momentum quantity) is to be used: this can be chosen
  /// via the optional scheme parameter. Use of this scheme option is
  /// discouraged in this case since @c RAPIDITY is only a valid option for
  /// vectors whose type is really the FourMomentum derived class.
  inline double deltaR(const FourVector& a, const FourVector& b,
                       RapScheme scheme=PSEUDORAPIDITY) {
    return sqrt(deltaR2(a, b, scheme));
  }



  /// @brief Calculate the squared 2D rapidity-azimuthal ("eta-phi") distance between two four-vectors.
  ///
  /// There is a scheme ambiguity for momentum-type four vectors
  /// as to whether the pseudorapidity (a purely geometric concept) or the
  /// rapidity (a relativistic energy-momentum quantity) is to be used: this can
  /// be chosen via the optional scheme parameter.
  inline double deltaR2(const FourVector& v,
                       double eta2, double phi2,
                       RapScheme scheme=PSEUDORAPIDITY) {
    switch (scheme) {
    case PSEUDORAPIDITY :
      return deltaR2(v.vector3(), eta2, phi2);
    case RAPIDITY:
      {
        const FourMomentum* mv = dynamic_cast<const FourMomentum*>(&v);
        if (!mv) {
          string err = "deltaR with scheme RAPIDITY can only be called with FourMomentum objects, not FourVectors";
          throw std::runtime_error(err);
        }
        return deltaR2(*mv, eta2, phi2, scheme);
      }
    default:
      throw std::runtime_error("The specified deltaR scheme is not yet implemented");
    }
  }

  /// @brief Calculate the 2D rapidity-azimuthal ("eta-phi") distance between two four-vectors.
  ///
  /// There is a scheme ambiguity for momentum-type four vectors
  /// as to whether the pseudorapidity (a purely geometric concept) or the
  /// rapidity (a relativistic energy-momentum quantity) is to be used: this can
  /// be chosen via the optional scheme parameter.
  inline double deltaR(const FourVector& v,
                       double eta2, double phi2,
                       RapScheme scheme=PSEUDORAPIDITY) {
    return sqrt(deltaR2(v, eta2, phi2, scheme));
  }


  /// @brief Calculate the squared 2D rapidity-azimuthal ("eta-phi") distance between two four-vectors.
  ///
  /// There is a scheme ambiguity for momentum-type four vectors
  /// as to whether the pseudorapidity (a purely geometric concept) or the
  /// rapidity (a relativistic energy-momentum quantity) is to be used: this can
  /// be chosen via the optional scheme parameter.
  inline double deltaR2(double eta1, double phi1,
                        const FourVector& v,
                        RapScheme scheme=PSEUDORAPIDITY) {
    switch (scheme) {
    case PSEUDORAPIDITY :
      return deltaR2(eta1, phi1, v.vector3());
    case RAPIDITY:
      {
        const FourMomentum* mv = dynamic_cast<const FourMomentum*>(&v);
        if (!mv) {
          string err = "deltaR with scheme RAPIDITY can only be called with FourMomentum objects, not FourVectors";
          throw std::runtime_error(err);
        }
        return deltaR2(eta1, phi1, *mv, scheme);
      }
    default:
      throw std::runtime_error("The specified deltaR scheme is not yet implemented");
    }
  }

  /// @brief Calculate the 2D rapidity-azimuthal ("eta-phi") distance between two four-vectors.
  ///
  /// There is a scheme ambiguity for momentum-type four vectors
  /// as to whether the pseudorapidity (a purely geometric concept) or the
  /// rapidity (a relativistic energy-momentum quantity) is to be used: this can
  /// be chosen via the optional scheme parameter.
  inline double deltaR(double eta1, double phi1,
                       const FourVector& v,
                       RapScheme scheme=PSEUDORAPIDITY) {
    return sqrt(deltaR2(eta1, phi1, v, scheme));
  }


  /// @brief Calculate the squared 2D rapidity-azimuthal ("eta-phi") distance between two four-vectors.
  ///
  /// There is a scheme ambiguity for momentum-type four vectors
  /// as to whether the pseudorapidity (a purely geometric concept) or the
  /// rapidity (a relativistic energy-momentum quantity) is to be used: this can
  /// be chosen via the optional scheme parameter.
  inline double deltaR2(const FourMomentum& a, const FourMomentum& b,
                       RapScheme scheme=PSEUDORAPIDITY) {
    switch (scheme) {
    case PSEUDORAPIDITY:
      return deltaR2(a.vector3(), b.vector3());
    case RAPIDITY:
      return deltaR2(a.rapidity(), a.azimuthalAngle(), b.rapidity(), b.azimuthalAngle());
    default:
      throw std::runtime_error("The specified deltaR scheme is not yet implemented");
    }
  }

  /// @brief Calculate the 2D rapidity-azimuthal ("eta-phi") distance between two four-vectors.
  ///
  /// There is a scheme ambiguity for momentum-type four vectors
  /// as to whether the pseudorapidity (a purely geometric concept) or the
  /// rapidity (a relativistic energy-momentum quantity) is to be used: this can
  /// be chosen via the optional scheme parameter.
  inline double deltaR(const FourMomentum& a, const FourMomentum& b,
                       RapScheme scheme=PSEUDORAPIDITY) {
    return sqrt(deltaR2(a, b, scheme));
  }


  /// @brief Calculate the squared 2D rapidity-azimuthal ("eta-phi") distance between two four-vectors.
  /// There is a scheme ambiguity for momentum-type four vectors
  /// as to whether the pseudorapidity (a purely geometric concept) or the
  /// rapidity (a relativistic energy-momentum quantity) is to be used: this can
  /// be chosen via the optional scheme parameter.
  inline double deltaR2(const FourMomentum& v,
                        double eta2, double phi2,
                        RapScheme scheme=PSEUDORAPIDITY) {
    switch (scheme) {
    case PSEUDORAPIDITY:
      return deltaR2(v.vector3(), eta2, phi2);
    case RAPIDITY:
      return deltaR2(v.rapidity(), v.azimuthalAngle(), eta2, phi2);
    default:
      throw std::runtime_error("The specified deltaR scheme is not yet implemented");
    }
  }

  /// @brief Calculate the 2D rapidity-azimuthal ("eta-phi") distance between two four-vectors.
  /// There is a scheme ambiguity for momentum-type four vectors
  /// as to whether the pseudorapidity (a purely geometric concept) or the
  /// rapidity (a relativistic energy-momentum quantity) is to be used: this can
  /// be chosen via the optional scheme parameter.
  inline double deltaR(const FourMomentum& v,
                       double eta2, double phi2,
                       RapScheme scheme=PSEUDORAPIDITY) {
    return sqrt(deltaR2(v, eta2, phi2, scheme));
  }


  /// @brief Calculate the squared 2D rapidity-azimuthal ("eta-phi") distance between two four-vectors.
  /// There is a scheme ambiguity for momentum-type four vectors
  /// as to whether the pseudorapidity (a purely geometric concept) or the
  /// rapidity (a relativistic energy-momentum quantity) is to be used: this can
  /// be chosen via the optional scheme parameter.
  inline double deltaR2(double eta1, double phi1,
                        const FourMomentum& v,
                        RapScheme scheme=PSEUDORAPIDITY) {
    switch (scheme) {
    case PSEUDORAPIDITY:
      return deltaR2(eta1, phi1, v.vector3());
    case RAPIDITY:
      return deltaR2(eta1, phi1, v.rapidity(), v.azimuthalAngle());
    default:
      throw std::runtime_error("The specified deltaR scheme is not yet implemented");
    }
  }

  /// @brief Calculate the 2D rapidity-azimuthal ("eta-phi") distance between two four-vectors.
  /// There is a scheme ambiguity for momentum-type four vectors
  /// as to whether the pseudorapidity (a purely geometric concept) or the
  /// rapidity (a relativistic energy-momentum quantity) is to be used: this can
  /// be chosen via the optional scheme parameter.
  inline double deltaR(double eta1, double phi1,
                        const FourMomentum& v,
                        RapScheme scheme=PSEUDORAPIDITY) {
    return sqrt(deltaR2(eta1, phi1, v, scheme));
  }


  /// @brief Calculate the squared 2D rapidity-azimuthal ("eta-phi") distance between two four-vectors.
  /// There is a scheme ambiguity for momentum-type four vectors
  /// as to whether the pseudorapidity (a purely geometric concept) or the
  /// rapidity (a relativistic energy-momentum quantity) is to be used: this can
  /// be chosen via the optional scheme parameter.
  inline double deltaR2(const FourMomentum& a, const FourVector& b,
                        RapScheme scheme=PSEUDORAPIDITY) {
    switch (scheme) {
    case PSEUDORAPIDITY:
      return deltaR2(a.vector3(), b.vector3());
    case RAPIDITY:
      return deltaR2(a.rapidity(), a.azimuthalAngle(), FourMomentum(b).rapidity(), b.azimuthalAngle());
    default:
      throw std::runtime_error("The specified deltaR scheme is not yet implemented");
    }
  }

  /// @brief Calculate the 2D rapidity-azimuthal ("eta-phi") distance between two four-vectors.
  /// There is a scheme ambiguity for momentum-type four vectors
  /// as to whether the pseudorapidity (a purely geometric concept) or the
  /// rapidity (a relativistic energy-momentum quantity) is to be used: this can
  /// be chosen via the optional scheme parameter.
  inline double deltaR(const FourMomentum& a, const FourVector& b,
                       RapScheme scheme=PSEUDORAPIDITY) {
    return sqrt(deltaR2(a, b, scheme));
  }


  /// @brief Calculate the 2D rapidity-azimuthal ("eta-phi") distance between two four-vectors.
  /// There is a scheme ambiguity for momentum-type four vectors
  /// as to whether the pseudorapidity (a purely geometric concept) or the
  /// rapidity (a relativistic energy-momentum quantity) is to be used: this can
  /// be chosen via the optional scheme parameter.
  inline double deltaR2(const FourVector& a, const FourMomentum& b,
                        RapScheme scheme=PSEUDORAPIDITY) {
    return deltaR2(b, a, scheme); //< note reversed args
  }

  /// @brief Calculate the 2D rapidity-azimuthal ("eta-phi") distance between two four-vectors.
  /// There is a scheme ambiguity for momentum-type four vectors
  /// as to whether the pseudorapidity (a purely geometric concept) or the
  /// rapidity (a relativistic energy-momentum quantity) is to be used: this can
  /// be chosen via the optional scheme parameter.
  inline double deltaR(const FourVector& a, const FourMomentum& b,
                       RapScheme scheme=PSEUDORAPIDITY) {
    return deltaR(b, a, scheme); //< note reversed args
  }


  /// @brief Calculate the 2D rapidity-azimuthal ("eta-phi") distance between a
  /// three-vector and a four-vector.
  inline double deltaR2(const FourMomentum& a, const Vector3& b) {
    return deltaR2(a.vector3(), b);
  }

  /// @brief Calculate the 2D rapidity-azimuthal ("eta-phi") distance between a
  /// three-vector and a four-vector.
  inline double deltaR(const FourMomentum& a, const Vector3& b) {
    return deltaR(a.vector3(), b);
  }

  /// @brief Calculate the 2D rapidity-azimuthal ("eta-phi") distance between a
  /// three-vector and a four-vector.
  inline double deltaR2(const Vector3& a, const FourMomentum& b) {
    return deltaR2(a, b.vector3());
  }

  /// @brief Calculate the 2D rapidity-azimuthal ("eta-phi") distance between a
  /// three-vector and a four-vector.
  inline double deltaR(const Vector3& a, const FourMomentum& b) {
    return deltaR(a, b.vector3());
  }

  /// @brief Calculate the 2D rapidity-azimuthal ("eta-phi") distance between a
  /// three-vector and a four-vector.
  inline double deltaR2(const FourVector& a, const Vector3& b) {
    return deltaR2(a.vector3(), b);
  }

  /// @brief Calculate the 2D rapidity-azimuthal ("eta-phi") distance between a
  /// three-vector and a four-vector.
  inline double deltaR(const FourVector& a, const Vector3& b) {
    return deltaR(a.vector3(), b);
  }

  /// @brief Calculate the 2D rapidity-azimuthal ("eta-phi") distance between a
  /// three-vector and a four-vector.
  inline double deltaR2(const Vector3& a, const FourVector& b) {
    return deltaR2(a, b.vector3());
  }

  /// @brief Calculate the 2D rapidity-azimuthal ("eta-phi") distance between a
  /// three-vector and a four-vector.
  inline double deltaR(const Vector3& a, const FourVector& b) {
    return deltaR(a, b.vector3());
  }

  //@}


  //////////////////////////////////////////////////////


  /// @name \f$ \Delta phi \f$ calculations from 4-vectors
  //@{

  /// Calculate the difference in azimuthal angle between two vectors.
  inline double deltaPhi(const FourMomentum& a, const FourMomentum& b, bool sign=false) {
    return deltaPhi(a.vector3(), b.vector3(), sign);
  }

  /// Calculate the difference in azimuthal angle between two vectors.
  inline double deltaPhi(const FourMomentum& v, double phi2, bool sign=false) {
    return deltaPhi(v.vector3(), phi2, sign);
  }

  /// Calculate the difference in azimuthal angle between two vectors.
  inline double deltaPhi(double phi1, const FourMomentum& v, bool sign=false) {
    return deltaPhi(phi1, v.vector3(), sign);
  }

  /// Calculate the difference in azimuthal angle between two vectors.
  inline double deltaPhi(const FourVector& a, const FourVector& b, bool sign=false) {
    return deltaPhi(a.vector3(), b.vector3(), sign);
  }

  /// Calculate the difference in azimuthal angle between two vectors.
  inline double deltaPhi(const FourVector& v, double phi2, bool sign=false) {
    return deltaPhi(v.vector3(), phi2, sign);
  }

  /// Calculate the difference in azimuthal angle between two vectors.
  inline double deltaPhi(double phi1, const FourVector& v, bool sign=false) {
    return deltaPhi(phi1, v.vector3(), sign);
  }

  /// Calculate the difference in azimuthal angle between two vectors.
  inline double deltaPhi(const FourVector& a, const FourMomentum& b, bool sign=false) {
    return deltaPhi(a.vector3(), b.vector3(), sign);
  }

  /// Calculate the difference in azimuthal angle between two vectors.
  inline double deltaPhi(const FourMomentum& a, const FourVector& b, bool sign=false) {
    return deltaPhi(a.vector3(), b.vector3(), sign);
  }

  /// Calculate the difference in azimuthal angle between two vectors.
  inline double deltaPhi(const FourVector& a, const Vector3& b, bool sign=false) {
    return deltaPhi(a.vector3(), b, sign);
  }

  /// Calculate the difference in azimuthal angle between two vectors.
  inline double deltaPhi(const Vector3& a, const FourVector& b, bool sign=false) {
    return deltaPhi(a, b.vector3(), sign);
  }

  /// Calculate the difference in azimuthal angle between two vectors.
  inline double deltaPhi(const FourMomentum& a, const Vector3& b, bool sign=false) {
    return deltaPhi(a.vector3(), b, sign);
  }

  /// Calculate the difference in azimuthal angle between two vectors.
  inline double deltaPhi(const Vector3& a, const FourMomentum& b, bool sign=false) {
    return deltaPhi(a, b.vector3(), sign);
  }

  //@}


  //////////////////////////////////////////////////////


  /// @name \f$ |\Delta eta| \f$ calculations from 4-vectors
  //@{

  /// Calculate the difference in pseudorapidity between two vectors.
  inline double deltaEta(const FourMomentum& a, const FourMomentum& b) {
    return deltaEta(a.vector3(), b.vector3());
  }

  /// Calculate the difference in pseudorapidity between two vectors.
  inline double deltaEta(const FourMomentum& v, double eta2) {
    return deltaEta(v.vector3(), eta2);
  }

  /// Calculate the difference in pseudorapidity between two vectors.
  inline double deltaEta(double eta1, const FourMomentum& v) {
    return deltaEta(eta1, v.vector3());
  }

  /// Calculate the difference in pseudorapidity between two vectors.
  inline double deltaEta(const FourVector& a, const FourVector& b) {
    return deltaEta(a.vector3(), b.vector3());
  }

  /// Calculate the difference in pseudorapidity between two vectors.
  inline double deltaEta(const FourVector& v, double eta2) {
    return deltaEta(v.vector3(), eta2);
  }

  /// Calculate the difference in pseudorapidity between two vectors.
  inline double deltaEta(double eta1, const FourVector& v) {
    return deltaEta(eta1, v.vector3());
  }

  /// Calculate the difference in pseudorapidity between two vectors.
  inline double deltaEta(const FourVector& a, const FourMomentum& b) {
    return deltaEta(a.vector3(), b.vector3());
  }

  /// Calculate the difference in pseudorapidity between two vectors.
  inline double deltaEta(const FourMomentum& a, const FourVector& b) {
    return deltaEta(a.vector3(), b.vector3());
  }

  /// Calculate the difference in pseudorapidity between two vectors.
  inline double deltaEta(const FourVector& a, const Vector3& b) {
    return deltaEta(a.vector3(), b);
  }

  /// Calculate the difference in pseudorapidity between two vectors.
  inline double deltaEta(const Vector3& a, const FourVector& b) {
    return deltaEta(a, b.vector3());
  }

  /// Calculate the difference in pseudorapidity between two vectors.
  inline double deltaEta(const FourMomentum& a, const Vector3& b) {
    return deltaEta(a.vector3(), b);
  }

  /// Calculate the difference in pseudorapidity between two vectors.
  inline double deltaEta(const Vector3& a, const FourMomentum& b) {
    return deltaEta(a, b.vector3());
  }

  //@}


  /// @name \f$ |\Delta y| \f$ calculations from 4-momentum vectors
  //@{

  /// Calculate the difference in rapidity between two 4-momentum vectors.
  inline double deltaRap(const FourMomentum& a, const FourMomentum& b) {
    return deltaRap(a.rapidity(), b.rapidity());
  }

  /// Calculate the difference in rapidity between two 4-momentum vectors.
  inline double deltaRap(const FourMomentum& v, double y2) {
    return deltaRap(v.rapidity(), y2);
  }

  /// Calculate the difference in rapidity between two 4-momentum vectors.
  inline double deltaRap(double y1, const FourMomentum& v) {
    return deltaRap(y1, v.rapidity());
  }

  //@}


  //////////////////////////////////////////////////////


  /// @name 4-vector comparison functions (for sorting)
  //@{

  /// Comparison to give a sorting by decreasing pT
  inline bool cmpMomByPt(const FourMomentum& a, const FourMomentum& b) {
    return a.pt() > b.pt();
  }
  /// Comparison to give a sorting by increasing pT
  inline bool cmpMomByAscPt(const FourMomentum& a, const FourMomentum& b) {
    return a.pt() < b.pt();
  }

  /// Comparison to give a sorting by decreasing 3-momentum magnitude |p|
  inline bool cmpMomByP(const FourMomentum& a, const FourMomentum& b) {
    return a.vector3().mod() > b.vector3().mod();
  }
  /// Comparison to give a sorting by increasing 3-momentum magnitude |p|
  inline bool cmpMomByAscP(const FourMomentum& a, const FourMomentum& b) {
    return a.vector3().mod() < b.vector3().mod();
  }

  /// Comparison to give a sorting by decreasing transverse energy
  inline bool cmpMomByEt(const FourMomentum& a, const FourMomentum& b) {
    return a.Et() > b.Et();
  }
  /// Comparison to give a sorting by increasing transverse energy
  inline bool cmpMomByAscEt(const FourMomentum& a, const FourMomentum& b) {
    return a.Et() < b.Et();
  }

  /// Comparison to give a sorting by decreasing energy
  inline bool cmpMomByE(const FourMomentum& a, const FourMomentum& b) {
    return a.E() > b.E();
  }
  /// Comparison to give a sorting by increasing energy
  inline bool cmpMomByAscE(const FourMomentum& a, const FourMomentum& b) {
    return a.E() < b.E();
  }

  /// Comparison to give a sorting by decreasing mass
  inline bool cmpMomByMass(const FourMomentum& a, const FourMomentum& b) {
    return a.mass() > b.mass();
  }
  /// Comparison to give a sorting by increasing mass
  inline bool cmpMomByAscMass(const FourMomentum& a, const FourMomentum& b) {
    return a.mass() < b.mass();
  }

  /// Comparison to give a sorting by increasing eta (pseudorapidity)
  inline bool cmpMomByEta(const FourMomentum& a, const FourMomentum& b) {
    return a.eta() < b.eta();
  }

  /// Comparison to give a sorting by decreasing eta (pseudorapidity)
  inline bool cmpMomByDescEta(const FourMomentum& a, const FourMomentum& b) {
    return a.pseudorapidity() > b.pseudorapidity();
  }

  /// Comparison to give a sorting by increasing absolute eta (pseudorapidity)
  inline bool cmpMomByAbsEta(const FourMomentum& a, const FourMomentum& b) {
    return fabs(a.eta()) < fabs(b.eta());
  }

  /// Comparison to give a sorting by increasing absolute eta (pseudorapidity)
  inline bool cmpMomByDescAbsEta(const FourMomentum& a, const FourMomentum& b) {
    return fabs(a.eta()) > fabs(b.eta());
  }

  /// Comparison to give a sorting by increasing rapidity
  inline bool cmpMomByRap(const FourMomentum& a, const FourMomentum& b) {
    return a.rapidity() < b.rapidity();
  }

  /// Comparison to give a sorting by decreasing rapidity
  inline bool cmpMomByDescRap(const FourMomentum& a, const FourMomentum& b) {
    return a.rapidity() > b.rapidity();
  }

  /// Comparison to give a sorting by increasing absolute rapidity
  inline bool cmpMomByAbsRap(const FourMomentum& a, const FourMomentum& b) {
    return fabs(a.rapidity()) < fabs(b.rapidity());
  }

  /// Comparison to give a sorting by decreasing absolute rapidity
  inline bool cmpMomByDescAbsRap(const FourMomentum& a, const FourMomentum& b) {
    return fabs(a.rapidity()) > fabs(b.rapidity());
  }

  /// @todo Add sorting by phi [0..2PI]


  /// Sort a container of momenta by cmp and return by reference for non-const inputs
  template<typename MOMS, typename CMP>
  inline MOMS& sortBy(MOMS& pbs, const CMP& cmp) {
    std::sort(pbs.begin(), pbs.end(), cmp);
    return pbs;
  }
  /// Sort a container of momenta by cmp and return by value for const inputs
  template<typename MOMS, typename CMP>
  inline MOMS sortBy(const MOMS& pbs, const CMP& cmp) {
    MOMS rtn = pbs;
    std::sort(rtn.begin(), rtn.end(), cmp);
    return rtn;
  }

  /// Sort a container of momenta by pT (decreasing) and return by reference for non-const inputs
  template<typename MOMS>
  inline MOMS& sortByPt(MOMS& pbs) {
    return sortBy(pbs, cmpMomByPt);
  }
  /// Sort a container of momenta by pT (decreasing) and return by value for const inputs
  template<typename MOMS>
  inline MOMS sortByPt(const MOMS& pbs) {
    return sortBy(pbs, cmpMomByPt);
  }

  /// Sort a container of momenta by E (decreasing) and return by reference for non-const inputs
  template<typename MOMS>
  inline MOMS& sortByE(MOMS& pbs) {
    return sortBy(pbs, cmpMomByE);
  }
  /// Sort a container of momenta by E (decreasing) and return by value for const inputs
  template<typename MOMS>
  inline MOMS sortByE(const MOMS& pbs) {
    return sortBy(pbs, cmpMomByE);
  }

  /// Sort a container of momenta by Et (decreasing) and return by reference for non-const inputs
  template<typename MOMS>
  inline MOMS& sortByEt(MOMS& pbs) {
    return sortBy(pbs, cmpMomByEt);
  }
  /// Sort a container of momenta by Et (decreasing) and return by value for const inputs
  template<typename MOMS>
  inline MOMS sortByEt(const MOMS& pbs) {
    return sortBy(pbs, cmpMomByEt);
  }

  //@}


  /// @name MT calculation
  //@{

  /// Calculate transverse mass of a visible and an invisible 3-vector
  inline double mT(const Vector3& vis, const Vector3& invis) {
    // return sqrt(2*vis.perp()*invis.perp() * (1 - cos(deltaPhi(vis, invis))) );
    return mT(vis.perp(), invis.perp(), deltaPhi(vis, invis));
  }

  /// Calculate transverse mass of a visible and an invisible 4-vector
  inline double mT(const FourMomentum& vis, const FourMomentum& invis) {
    return mT(vis.p3(), invis.p3());
  }

  /// Calculate transverse mass of a visible 4-vector and an invisible 3-vector
  inline double mT(const FourMomentum& vis, const Vector3& invis) {
    return mT(vis.p3(), invis);
  }

  /// Calculate transverse mass of a visible 4-vector and an invisible 3-vector
  inline double mT(const Vector3& vis, const FourMomentum& invis) {
    return mT(vis, invis.p3());
  }

  //@}


  //////////////////////////////////////////////////////


  /// @name 4-vector string representations
  //@{

  /// Render a 4-vector as a string.
  inline std::string toString(const FourVector& lv) {
    ostringstream out;
    out << "("  << (fabs(lv.t()) < 1E-30 ? 0.0 : lv.t())
        << "; " << (fabs(lv.x()) < 1E-30 ? 0.0 : lv.x())
        << ", " << (fabs(lv.y()) < 1E-30 ? 0.0 : lv.y())
        << ", " << (fabs(lv.z()) < 1E-30 ? 0.0 : lv.z())
        << ")";
    return out.str();
  }

  /// Write a 4-vector to an ostream.
  inline std::ostream& operator<<(std::ostream& out, const FourVector& lv) {
    out << toString(lv);
    return out;
  }

  //@}


  /// @name Typedefs of vector types to short names
  /// @todo Switch canonical and alias names
  //@{
  //typedef FourVector V4; //< generic
  typedef FourVector X4; //< spatial
  typedef FourMomentum P4; //< momentum
  //@}

  /// @name Typedefs for lists of vector types
  //@{
  typedef std::vector<FourVector> FourVectors;
  typedef std::vector<FourMomentum> FourMomenta;
  typedef std::vector<X4> X4s;
  typedef std::vector<P4> P4a;
  //@}


}

#endif
