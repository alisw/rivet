// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// D0 direct photon pair production
  class D0_2010_S8570965 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(D0_2010_S8570965);


    void init() {
      FinalState fs;
      declare(fs, "FS");

      IdentifiedFinalState ifs(Cuts::abseta < 0.9 && Cuts::pT > 20*GeV);
      ifs.acceptId(PID::PHOTON);
      declare(ifs, "IFS");

      book(_h_M ,1, 1, 1);
      book(_h_pT ,2, 1, 1);
      book(_h_dPhi ,3, 1, 1);
      book(_h_costheta ,4, 1, 1);

      std::pair<double, double> M_ranges[] = { std::make_pair(30.0, 50.0),
                                               std::make_pair(50.0, 80.0),
                                               std::make_pair(80.0, 350.0) };

      for (size_t i = 0; i < 3; ++i) {
        Histo1DPtr a,b,c;
        _h_pT_M.add(M_ranges[i].first, M_ranges[i].second, book(a, 5+3*i, 1, 1));
        _h_dPhi_M.add(M_ranges[i].first, M_ranges[i].second, book(b, 6+3*i, 1, 1));
        _h_costheta_M.add(M_ranges[i].first, M_ranges[i].second, book(c, 7+3*i, 1, 1));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = 1.0;

      Particles photons = apply<IdentifiedFinalState>(event, "IFS").particlesByPt();
      if (photons.size() < 2 ||
          (photons[0].pT() < 21.0*GeV)) {
        vetoEvent;
      }

      // Isolate photons with ET_sum in cone
      Particles isolated_photons;
      Particles fs = apply<FinalState>(event, "FS").particles();
      for (const Particle& photon : photons) {
        double eta_P = photon.eta();
        double phi_P = photon.phi();
        double Etsum=0.0;
        for (const Particle& p : fs) {
          if (HepMCUtils::uniqueId(p.genParticle()) != HepMCUtils::uniqueId(photon.genParticle()) &&
              deltaR(eta_P, phi_P, p.eta(), p.phi()) < 0.4) {
            Etsum += p.Et();
          }
        }
        if (Etsum < 2.5*GeV) {
          isolated_photons.push_back(photon);
        }
      }

      if (isolated_photons.size() != 2) {
        vetoEvent;
      }
      std::sort(isolated_photons.begin(), isolated_photons.end(), cmpMomByPt);

      FourMomentum y1=isolated_photons[0].momentum();
      FourMomentum y2=isolated_photons[1].momentum();
      if (deltaR(y1, y2)<0.4) {
        vetoEvent;
      }

      FourMomentum yy=y1+y2;
      double Myy = yy.mass()/GeV;
      if (Myy<30.0 || Myy>350.0) {
        vetoEvent;
      }

      double pTyy = yy.pT()/GeV;
      if (Myy<pTyy) {
        vetoEvent;
      }

      double dPhiyy = mapAngle0ToPi(y1.phi()-y2.phi());
      if (dPhiyy<0.5*M_PI) {
        vetoEvent;
      }

      double costhetayy = fabs(tanh((y1.eta()-y2.eta())/2.0));

      _h_M->fill(Myy, weight);
      _h_pT->fill(pTyy, weight);
      _h_dPhi->fill(dPhiyy, weight);
      _h_costheta->fill(costhetayy, weight);

      _h_pT_M.fill(Myy, pTyy, weight);
      _h_dPhi_M.fill(Myy, dPhiyy, weight);
      _h_costheta_M.fill(Myy, costhetayy, weight);
    }


    void finalize() {

      scale(_h_M, crossSection()/sumOfWeights());
      scale(_h_pT, crossSection()/sumOfWeights());
      scale(_h_dPhi, crossSection()/sumOfWeights());
      scale(_h_costheta, crossSection()/sumOfWeights());

      _h_pT_M.scale(crossSection()/sumOfWeights(), this);
      _h_dPhi_M.scale(crossSection()/sumOfWeights(), this);
      _h_costheta_M.scale(crossSection()/sumOfWeights(), this);

    }


  private:

    Histo1DPtr _h_M;
    Histo1DPtr _h_pT;
    Histo1DPtr _h_dPhi;
    Histo1DPtr _h_costheta;
    BinnedHistogram _h_pT_M;
    BinnedHistogram _h_dPhi_M;
    BinnedHistogram _h_costheta_M;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(D0_2010_S8570965, D0_2010_I846997);

}
