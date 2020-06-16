// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  class CMS_2018_I1653948 : public Analysis {
  public:

    CMS_2018_I1653948()
      : Analysis("CMS_2018_I1653948"), _xi_hf_cut(1E-6), _xi_castor_cut(1E-7)
    {    }


    /// Book projections and histograms
    void init() {
      declare(FinalState(),"FS");
      book(_h_xsec, 1, 1, 1);
    }


    /// Analyze each event
    void analyze(const Event& event) {

      const FinalState& fs = applyProjection<FinalState>(event, "FS");
      if (fs.size() < 3) vetoEvent; // veto on elastic events
      const Particles particlesByRapidity = fs.particles(cmpMomByRap);
      const size_t num_particles = particlesByRapidity.size();

      vector<double> gaps;
      vector<double> midpoints;

      for (size_t ip = 1; ip < num_particles; ++ip) {
        const Particle& p1 = particlesByRapidity[ip-1];
        const Particle& p2 = particlesByRapidity[ip];
        const double gap = p2.momentum().rapidity() - p1.momentum().rapidity();
        const double mid = (p2.momentum().rapidity() + p1.momentum().rapidity()) / 2.;
        gaps.push_back(gap);
        midpoints.push_back(mid);
      }

      int imid = std::distance(gaps.begin(), max_element(gaps.begin(), gaps.end()));
      double gapcenter = midpoints[imid];

      FourMomentum MxFourVector(0.,0.,0.,0.);
      FourMomentum MyFourVector(0.,0.,0.,0.);

      for (const Particle& p : fs.particles(cmpMomByEta)) {
        if (p.momentum().rapidity() < gapcenter) {
          MxFourVector += p.momentum();
        } else {
          MyFourVector += p.momentum();
        }
      }

      double Mx = MxFourVector.mass();
      double My = MyFourVector.mass();

      double xix = (Mx * Mx) / (sqrtS()/GeV * sqrtS()/GeV);
      double xiy = (My * My) / (sqrtS()/GeV * sqrtS()/GeV);
      double xi  = max(xix, xiy);

      if (xi > _xi_hf_cut) _h_xsec->fill(0.5);
      if (xix > _xi_castor_cut || xiy > _xi_hf_cut) _h_xsec->fill(1.5);
    }


    /// Normalizations, etc.
    void finalize() {
      scale(_h_xsec, crossSection()/millibarn/sumOfWeights());
    }


  private:

    Histo1DPtr _h_xsec;
    double _xi_hf_cut;
    double _xi_castor_cut;

  };


  DECLARE_RIVET_PLUGIN(CMS_2018_I1653948);

}
