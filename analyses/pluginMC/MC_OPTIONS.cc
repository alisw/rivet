// -*- C++ -*-
#include "Rivet/Analysis.hh"

namespace Rivet {
  using std::istream;
  using std::ostream;

  // Example analysis to show how to use options in an analysis
  // Example of a custom class to be read in from an option.
  class A {
  public:
    A() : a(-1.0) {}
  private:
    double a;
    // Custom class must be streamable.
    friend istream& operator>> (istream& is, A& a);
    friend ostream& operator<< (ostream& os, const A& a);
  };
  
  // Custom class must be streamable.
  istream& operator>> (istream& is, A& a) {
    is >> a.a;
    return is;
  }
  ostream& operator<< (ostream& os, const A& a) {
    os << a.a;
    return os;
  }

  class MC_OPTIONS : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MC_OPTIONS);

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {


      // Parameters read in.
      // A double.
      double f = getOption<double>("foo", 1.0);
      // A string.
      string s = getOption<string>("bar", "");
      // A custom object. 
      A a = getOption<A>("baz", A());

      cout << "foo = " << f << endl;
      cout << "bar = " << s << endl;
      cout << "baz = " << a << endl;
      value = f;
      book(h, "hist",10,0,10);
    }

    // Perform the per-event analysis
    void analyze(const Event& event) {
      h->fill(value);

    }


    /// Finalize
    void finalize() {
    
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr h;
    //@}
    double value;
    
  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MC_OPTIONS);

}
