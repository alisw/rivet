#ifndef RIVET_MATH_VECTOR2
#define RIVET_MATH_VECTOR2

#include "Rivet/Math/MathConstants.hh"
#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Math/VectorN.hh"

namespace Rivet {


  class Vector2;
  typedef Vector2 TwoVector;
  //class Matrix2;

  Vector2 multiply(const double, const Vector2&);
  Vector2 multiply(const Vector2&, const double);
  Vector2 add(const Vector2&, const Vector2&);
  Vector2 operator*(const double, const Vector2&);
  Vector2 operator*(const Vector2&, const double);
  Vector2 operator/(const Vector2&, const double);
  Vector2 operator+(const Vector2&, const Vector2&);
  Vector2 operator-(const Vector2&, const Vector2&);


  /// @brief Two-dimensional specialisation of Vector.
  class Vector2 : public Vector<2> {

    // friend class Matrix2;
    friend Vector2 multiply(const double, const Vector2&);
    friend Vector2 multiply(const Vector2&, const double);
    friend Vector2 add(const Vector2&, const Vector2&);
    friend Vector2 subtract(const Vector2&, const Vector2&);

  public:
    Vector2() : Vector<2>() { }

    template<typename V2>
    Vector2(const V2& other) {
      this->setX(other.x());
      this->setY(other.y());
    }

    Vector2(const Vector<2>& other) {
      this->setX(other.get(0));
      this->setY(other.get(1));
    }

    Vector2(double x, double y) {
      this->setX(x);
      this->setY(y);
    }

    ~Vector2() { }


  public:

    static Vector2 mkX() { return Vector2(1,0); }
    static Vector2 mkY() { return Vector2(0,1); }


  public:

    double x() const { return get(0); }
    double y() const { return get(1); }
    Vector2& setX(double x) { set(0, x); return *this; }
    Vector2& setY(double y) { set(1, y); return *this; }


    /// Dot-product with another vector
    double dot(const Vector2& v) const {
      return _vec.dot(v._vec);
    }

    /// Angle in radians to another vector
    double angle(const Vector2& v) const {
      const double localDotOther = unit().dot(v.unit());
      if (localDotOther > 1.0) return 0.0;
      if (localDotOther < -1.0) return M_PI;
      return acos(localDotOther);
    }


    /// Unit-normalized version of this vector.
    Vector2 unitVec() const {
      double md = mod();
      if ( md <= 0.0 ) return Vector2();
      else return *this * 1.0/md;
    }

    /// Synonym for unitVec
    Vector2 unit() const {
      return unitVec();
    }


  public:
    Vector2& operator*=(const double a) {
      _vec = multiply(a, *this)._vec;
      return *this;
    }

    Vector2& operator/=(const double a) {
      _vec = multiply(1.0/a, *this)._vec;
      return *this;
    }

    Vector2& operator+=(const Vector2& v) {
      _vec = add(*this, v)._vec;
      return *this;
    }

    Vector2& operator-=(const Vector2& v) {
      _vec = subtract(*this, v)._vec;
      return *this;
    }

    Vector2 operator-() const {
      Vector2 rtn;
      rtn._vec = -_vec;
      return rtn;
    }

  };



  inline double dot(const Vector2& a, const Vector2& b) {
    return a.dot(b);
  }

  inline Vector2 multiply(const double a, const Vector2& v) {
    Vector2 result;
    result._vec = a * v._vec;
    return result;
  }

  inline Vector2 multiply(const Vector2& v, const double a) {
    return multiply(a, v);
  }

  inline Vector2 operator*(const double a, const Vector2& v) {
    return multiply(a, v);
  }

  inline Vector2 operator*(const Vector2& v, const double a) {
    return multiply(a, v);
  }

  inline Vector2 operator/(const Vector2& v, const double a) {
    return multiply(1.0/a, v);
  }

  inline Vector2 add(const Vector2& a, const Vector2& b) {
    Vector2 result;
    result._vec = a._vec + b._vec;
    return result;
  }

  inline Vector2 subtract(const Vector2& a, const Vector2& b) {
    Vector2 result;
    result._vec = a._vec - b._vec;
    return result;
  }

  inline Vector2 operator+(const Vector2& a, const Vector2& b) {
    return add(a, b);
  }

  inline Vector2 operator-(const Vector2& a, const Vector2& b) {
    return subtract(a, b);
  }

  /// Angle (in radians) between two 2-vectors.
  inline double angle(const Vector2& a, const Vector2& b) {
    return a.angle(b);
  }


}

#endif
