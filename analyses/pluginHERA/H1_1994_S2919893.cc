// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Math/Constants.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"

namespace Rivet {


  /// @brief H1 energy flow and charged particle spectra
  ///
  /// @author Peter Richardson
  /// Based on the equivalent HZTool analysis
  class H1_1994_S2919893 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(H1_1994_S2919893);


    /// @name Analysis methods
    //@{

    /// Initialise projections and histograms
    void init() {
      // Projections
      declare(DISLepton(), "Lepton");
      declare(DISKinematics(), "Kinematics");
      declare(FinalState(), "FS");

      // Histos
      book(_histEnergyFlowLowX ,1, 1, 1);
      book(_histEnergyFlowHighX ,1, 1, 2);

      book(_histEECLowX ,2, 1, 1);
      book(_histEECHighX ,2, 1, 2);

      book(_histSpectraW77 ,3, 1, 1);
      book(_histSpectraW122 ,3, 1, 2);
      book(_histSpectraW169 ,3, 1, 3);
      book(_histSpectraW117 ,3, 1, 4);

      book(_histPT2 ,4, 1, 1);

      book(_w77 .first, "TMP/w77_1");
      book(_w122.first, "TMP/w122_1");
      book(_w169.first, "TMP/w169_1");
      book(_w117.first, "TMP/w117_1");
      book(_wEnergy.first, "TMP/wEnergy_1");

      book(_w77 .second, "TMP/w77_2");
      book(_w122.second, "TMP/w122_2");
      book(_w169.second, "TMP/w169_2");
      book(_w117.second, "TMP/w117_2");
      book(_wEnergy.second, "TMP/wEnergy_2");
    }


    /// Analyse each event
    void analyze(const Event& event) {

      // Get the DIS kinematics
      const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");
      if ( dk.failed() ) vetoEvent;
      const double x  = dk.x();
      const double w2 = dk.W2();
      const double w = sqrt(w2);

      // Momentum of the scattered lepton
      const DISLepton& dl = apply<DISLepton>(event,"Lepton");
      if ( dl.failed() ) return;
      const FourMomentum leptonMom = dl.out();
      const double ptel = leptonMom.pT();
      const double enel = leptonMom.E();
      const double thel = leptonMom.angle(dk.beamHadron().mom())/degree;

      // Extract the particles other than the lepton
      const FinalState& fs = apply<FinalState>(event, "FS");
      Particles particles;
      particles.reserve(fs.particles().size());
      ConstGenParticlePtr dislepGP = dl.out().genParticle();
      for(const Particle& p: fs.particles()) {
        ConstGenParticlePtr loopGP = p.genParticle();
        if (loopGP == dislepGP) continue;
        particles.push_back(p);
      }

      // Cut on the forward energy
      double efwd = 0.0;
      for (const Particle& p : particles) {
        const double th = p.angle(dk.beamHadron())/degree;
        if (inRange(th, 4.4, 15)) efwd += p.E();
      }

      // Apply the cuts
      // Lepton energy and angle, w2 and forward energy
      MSG_DEBUG("enel/GeV = " << enel/GeV << ", thel = " << thel
                << ", w2 = " << w2 << ", efwd/GeV = " << efwd/GeV);
      bool cut = enel/GeV > 14. && thel > 157. && thel < 172.5 && w2 >= 3000. && efwd/GeV > 0.5;
      if (!cut) vetoEvent;

      // Weight of the event
      (x < 1e-3 ? _wEnergy.first : _wEnergy.second)->fill();

      // Boost to hadronic CM
      const LorentzTransform hcmboost = dk.boostHCM();
      // Loop over the particles
      long ncharged(0);
      for (size_t ip1 = 0; ip1 < particles.size(); ++ip1) {
        const Particle& p = particles[ip1];

        const double th = p.angle(dk.beamHadron().momentum()) / degree;
        // Boost momentum to lab
        const FourMomentum hcmMom = hcmboost.transform(p.momentum());
        // Angular cut
        if (th <= 4.4) continue;

        // Energy flow histogram
        const double et = fabs(hcmMom.Et());
        const double eta = hcmMom.eta();
        (x < 1e-3 ? _histEnergyFlowLowX : _histEnergyFlowHighX)->fill(eta, et);
        if (PID::charge3(p.pid()) != 0) {
          /// @todo Use units in w comparisons... what are the units?
          if (w > 50. && w <= 200.) {
            double xf= 2 * hcmMom.z() / w;
            double pt2 = hcmMom.pT2();
            if (w > 50. && w <= 100.) {
              _histSpectraW77 ->fill(xf);
            } else if (w > 100. && w <= 150.) {
              _histSpectraW122->fill(xf);
            } else if (w > 150. && w <= 200.) {
              _histSpectraW169->fill(xf);
            }
            _histSpectraW117->fill(xf);
            /// @todo Is this profile meant to be filled with 2 weight factors?
            _histPT2->fill(xf, pt2/GeV2);
            ++ncharged;
          }
        }


        // Energy-energy correlation
        if (th <= 8.) continue;
        double phi1 = p.phi(ZERO_2PI);
        double eta1 = p.eta();
        double et1 = fabs(p.momentum().Et());
        for (size_t ip2 = ip1+1; ip2 < particles.size(); ++ip2) {
          const Particle& p2 = particles[ip2];

          //double th2 = beamAngle(p2.momentum(), order);
          double th2 = p2.angle(dk.beamHadron().momentum()) / degree;
          if (th2 <= 8.) continue;
          double phi2 = p2.phi(ZERO_2PI);

          /// @todo Use angle function
          double deltaphi = phi1 - phi2;
          if (fabs(deltaphi) > PI) deltaphi = fabs(fabs(deltaphi) - TWOPI);
          double eta2 = p2.eta();
          double omega = sqrt(sqr(eta1-eta2) + sqr(deltaphi));
          double et2 = fabs(p2.momentum().Et());
          double wt = et1*et2 / sqr(ptel);
          (x < 1e-3 ? _histEECLowX : _histEECHighX)->fill(omega, wt);
        }
      }

      // Factors for normalization
      if (w > 50. && w <= 200.) {
        if (w <= 100.) {
          _w77.first ->fill(ncharged);
          _w77.second->fill();
        } else if (w <= 150.) {
          _w122.first ->fill(ncharged);
          _w122.second->fill();
        } else {
          _w169.first ->fill(ncharged);
          _w169.second->fill();
        }
        _w117.first ->fill(ncharged);
        _w117.second->fill();
      }
    }


    // Normalize inclusive single particle distributions to the average number of charged particles per event.
    void finalize() {
      normalize(_histSpectraW77,  *_w77.first/ *_w77.second);
      normalize(_histSpectraW122,  *_w122.first/ *_w122.second);
      normalize(_histSpectraW169,  *_w169.first/ *_w169.second);
      normalize(_histSpectraW117,  *_w117.first/ *_w117.second);

      scale(_histEnergyFlowLowX , 1./ *_wEnergy.first );
      scale(_histEnergyFlowHighX, 1./ *_wEnergy.second);

      scale(_histEECLowX , 1./ *_wEnergy.first );
      scale(_histEECHighX, 1./ *_wEnergy.second);
    }

    //@}


  private:

    /// Polar angle with right direction of the beam
    inline double beamAngle(const FourVector& v, bool order) {
      double thel = v.polarAngle()/degree;
      if (thel < 0) thel += 180.;
      if (!order) thel = 180 - thel;
      return thel;
    }

    /// @name Histograms
    /// @{
    Histo1DPtr _histEnergyFlowLowX, _histEnergyFlowHighX;
    Histo1DPtr _histEECLowX, _histEECHighX;
    Histo1DPtr _histSpectraW77, _histSpectraW122, _histSpectraW169, _histSpectraW117;
    Profile1DPtr _histPT2;
    /// @}

    /// @name Storage of weights to calculate averages for normalisation
    /// @{
    pair<CounterPtr,CounterPtr> _w77, _w122, _w169, _w117, _wEnergy;
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(H1_1994_S2919893, H1_1994_I372350);

}
