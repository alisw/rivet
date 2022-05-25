// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/ParisiTensor.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// ALEPH jet rates and event shapes at LEP 1 and 2
  class ALEPH_2004_S5765862 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(ALEPH_2004_S5765862);


    void init() {
      _initialisedJets    = true;
      _initialisedSpectra = true;
      // TODO: According to the paper they seem to discard neutral particles
      //       between 1 and 2 GeV. That correction is included in the systematic
      //       uncertainties and overly complicated to program, so we ignore it.
      const FinalState fs;
      declare(fs, "FS");
      FastJets durhamjets(fs, FastJets::DURHAM, 0.7, JetAlg::Muons::ALL, JetAlg::Invisibles::ALL);
      declare(durhamjets, "DurhamJets");

      const Thrust thrust(fs);
      declare(thrust, "Thrust");
      declare(Sphericity(fs), "Sphericity");
      declare(ParisiTensor(fs), "Parisi");
      declare(Hemispheres(thrust), "Hemispheres");

      const ChargedFinalState cfs;
      declare(Beam(), "Beams");
      declare(cfs, "CFS");

      // Histos
      // offset for the event shapes and jets
      int offset = 0;
      switch (int(sqrtS()/GeV + 0.5)) {
        case 91: offset = 0; break;
        case 133: offset = 1; break;
        case 161: offset = 2; break;
        case 172: offset = 3; break;
        case 183: offset = 4; break;
        case 189: offset = 5; break;
        case 200: offset = 6; break;
        case 206: offset = 7; break;
        default:
          _initialisedJets = false;
      }
      // event shapes
      if(_initialisedJets) {
        book(_h_thrust ,offset+54, 1, 1);
        book(_h_heavyjetmass ,offset+62, 1, 1);
        book(_h_totaljetbroadening ,offset+70, 1, 1);
        book(_h_widejetbroadening ,offset+78, 1, 1);
        book(_h_cparameter ,offset+86, 1, 1);
        book(_h_thrustmajor ,offset+94, 1, 1);
        book(_h_thrustminor ,offset+102, 1, 1);
        book(_h_jetmassdifference ,offset+110, 1, 1);
        book(_h_aplanarity ,offset+118, 1, 1);
        if ( offset != 0 )
          book(_h_planarity, offset+125, 1, 1);
        book(_h_oblateness ,offset+133, 1, 1);
        book(_h_sphericity ,offset+141, 1, 1);

        // Durham n->m jet resolutions
        book(_h_y_Durham[0] ,offset+149, 1, 1);   // y12 d149 ... d156
        book(_h_y_Durham[1] ,offset+157, 1, 1);   // y23 d157 ... d164
        if (offset < 6) { // there is no y34, y45 and y56 for 200 gev
          book(_h_y_Durham[2] ,offset+165, 1, 1); // y34 d165 ... d172, but not 171
          book(_h_y_Durham[3] ,offset+173, 1, 1); // y45 d173 ... d179
          book(_h_y_Durham[4] ,offset+180, 1, 1); // y56 d180 ... d186
        }
        else if (offset == 6) {
          _h_y_Durham[2] = Histo1DPtr();
          _h_y_Durham[3] = Histo1DPtr();
          _h_y_Durham[4] = Histo1DPtr();
        }
        else if (offset == 7) {
          book(_h_y_Durham[2] ,172, 1, 1);
          book(_h_y_Durham[3] ,179, 1, 1);
          book(_h_y_Durham[4] ,186, 1, 1);
        }

        // Durham n-jet fractions
        book(_h_R_Durham[0] ,offset+187, 1, 1); // R1 d187 ... d194
        book(_h_R_Durham[1] ,offset+195, 1, 1); // R2 d195 ... d202
        book(_h_R_Durham[2] ,offset+203, 1, 1); // R3 d203 ... d210
        book(_h_R_Durham[3] ,offset+211, 1, 1); // R4 d211 ... d218
        book(_h_R_Durham[4] ,offset+219, 1, 1); // R5 d219 ... d226
        book(_h_R_Durham[5] ,offset+227, 1, 1); // R>=6 d227 ... d234
      }
      // offset for the charged particle distributions
      offset = 0;
      switch (int(sqrtS() + 0.5)) {
        case 133: offset = 0; break;
        case 161: offset = 1; break;
        case 172: offset = 2; break;
        case 183: offset = 3; break;
        case 189: offset = 4; break;
        case 196: offset = 5; break;
        case 200: offset = 6; break;
        case 206: offset = 7; break;
        default:
          _initialisedSpectra = false;
      }
      if (_initialisedSpectra) {
        book(_h_xp , 2+offset, 1, 1);
        book(_h_xi ,11+offset, 1, 1);
        book(_h_xe ,19+offset, 1, 1);
        book(_h_pTin  ,27+offset, 1, 1);
        if (offset == 7)
          book(_h_pTout, 35, 1, 1);
        book(_h_rapidityT ,36+offset, 1, 1);
        book(_h_rapidityS ,44+offset, 1, 1);
      }
      book(_weightedTotalChargedPartNum, "_weightedTotalChargedPartNum");
      if (!_initialisedSpectra && !_initialisedJets) {
        MSG_WARNING("CoM energy of events sqrt(s) = " << sqrtS()/GeV
                    << " doesn't match any available analysis energy .");
      }

      book(mult, 1, 1, 1);
    }


    void analyze(const Event& e) {

      const Thrust& thrust = apply<Thrust>(e, "Thrust");
      const Sphericity& sphericity = apply<Sphericity>(e, "Sphericity");

      if(_initialisedJets) {
        bool LEP1 = isCompatibleWithSqrtS(91.2*GeV,0.01);
        // event shapes
        double thr = LEP1 ? thrust.thrust() : 1.0 - thrust.thrust();
        _h_thrust->fill(thr);
        _h_thrustmajor->fill(thrust.thrustMajor());
        if(LEP1)
          _h_thrustminor->fill(log(thrust.thrustMinor()));
        else
          _h_thrustminor->fill(thrust.thrustMinor());
        _h_oblateness->fill(thrust.oblateness());

        const Hemispheres& hemi = apply<Hemispheres>(e, "Hemispheres");
        _h_heavyjetmass->fill(hemi.scaledM2high());
        _h_jetmassdifference->fill(hemi.scaledM2diff());
        _h_totaljetbroadening->fill(hemi.Bsum());
        _h_widejetbroadening->fill(hemi.Bmax());

        const ParisiTensor& parisi = apply<ParisiTensor>(e, "Parisi");
        _h_cparameter->fill(parisi.C());

        _h_aplanarity->fill(sphericity.aplanarity());
        if(_h_planarity)
          _h_planarity->fill(sphericity.planarity());
        _h_sphericity->fill(sphericity.sphericity());

        // Jet rates
        const FastJets& durjet = apply<FastJets>(e, "DurhamJets");
        double log10e = log10(exp(1.));
        if (durjet.clusterSeq()) {
          double logynm1=0.;
          double logyn;
          for (size_t i=0; i<5; ++i) {
            double yn = durjet.clusterSeq()->exclusive_ymerge_max(i+1);
            if (yn<=0.0) continue;
            logyn = -log(yn);
            if (_h_y_Durham[i]) {
              _h_y_Durham[i]->fill(logyn);
            }
            if(!LEP1) logyn *= log10e;
            for (size_t j = 0; j < _h_R_Durham[i]->numBins(); ++j) {
              double val   = _h_R_Durham[i]->bin(j).xMin();
              double width = _h_R_Durham[i]->bin(j).xWidth();
              if(-val<=logynm1) break;
              if(-val<logyn) {
                _h_R_Durham[i]->fill(val+0.5*width, width);
              }
            }
            logynm1 = logyn;
          }
          for (size_t j = 0; j < _h_R_Durham[5]->numBins(); ++j) {
            double val   = _h_R_Durham[5]->bin(j).xMin();
            double width = _h_R_Durham[5]->bin(j).xWidth();
            if(-val<=logynm1) break;
            _h_R_Durham[5]->fill(val+0.5*width, width);
          }
        }
        if( !_initialisedSpectra) {
          const ChargedFinalState& cfs = apply<ChargedFinalState>(e, "CFS");
          const size_t numParticles = cfs.particles().size();
          _weightedTotalChargedPartNum->fill(numParticles);
        }
      }

      // charged particle distributions
      if(_initialisedSpectra) {
        const ChargedFinalState& cfs = apply<ChargedFinalState>(e, "CFS");
        const size_t numParticles = cfs.particles().size();
        _weightedTotalChargedPartNum->fill(numParticles);
        const ParticlePair& beams = apply<Beam>(e, "Beams").beams();
        const double meanBeamMom = ( beams.first.p3().mod() +
                                     beams.second.p3().mod() ) / 2.0;
        for (const Particle& p : cfs.particles()) {
          const double xp = p.p3().mod()/meanBeamMom;
          _h_xp->fill(xp   );
          const double logxp = -std::log(xp);
          _h_xi->fill(logxp);
          const double xe = p.E()/meanBeamMom;
          _h_xe->fill(xe   );
          const double pTinT  = dot(p.p3(), thrust.thrustMajorAxis());
          const double pToutT = dot(p.p3(), thrust.thrustMinorAxis());
          _h_pTin->fill(fabs(pTinT/GeV));
          if(_h_pTout) _h_pTout->fill(fabs(pToutT/GeV));
          const double momT = dot(thrust.thrustAxis()        ,p.p3());
          const double rapidityT = 0.5 * std::log((p.E() + momT) /
                                                  (p.E() - momT));
          _h_rapidityT->fill(fabs(rapidityT));
          const double momS = dot(sphericity.sphericityAxis(),p.p3());
          const double rapidityS = 0.5 * std::log((p.E() + momS) /
                                                  (p.E() - momS));
          _h_rapidityS->fill(fabs(rapidityS));
        }
      }
    }

    void finalize() {
      if(!_initialisedJets && !_initialisedSpectra) return;

      if (_initialisedJets) {
        normalize(_h_thrust);
        normalize(_h_heavyjetmass);
        normalize(_h_totaljetbroadening);
        normalize(_h_widejetbroadening);
        normalize(_h_cparameter);
        normalize(_h_thrustmajor);
        normalize(_h_thrustminor);
        normalize(_h_jetmassdifference);
        normalize(_h_aplanarity);
        if(_h_planarity) normalize(_h_planarity);
        normalize(_h_oblateness);
        normalize(_h_sphericity);

        for (size_t n=0; n<6; ++n) {
          scale(_h_R_Durham[n], 1./sumOfWeights());
        }

        for (size_t n = 0; n < 5; ++n) {
          if (_h_y_Durham[n]) {
            scale(_h_y_Durham[n], 1.0/sumOfWeights());
          }
        }
      }

      Histo1D temphisto(refData(1, 1, 1));
      const double avgNumParts = dbl(*_weightedTotalChargedPartNum) / sumOfWeights();


      for (size_t b = 0; b < temphisto.numBins(); b++) {
        const double x  = temphisto.bin(b).xMid();
        const double ex = temphisto.bin(b).xWidth()/2.;
        if (inRange(sqrtS()/GeV, x-ex, x+ex)) {
          mult->addPoint(x, avgNumParts, ex, 0.);
        }
      }

      if (_initialisedSpectra) {
        normalize(_h_xp, avgNumParts);
        normalize(_h_xi, avgNumParts);
        normalize(_h_xe, avgNumParts);
        normalize(_h_pTin , avgNumParts);
        if (_h_pTout) normalize(_h_pTout, avgNumParts);
        normalize(_h_rapidityT, avgNumParts);
        normalize(_h_rapidityS, avgNumParts);
      }
    }


  private:

    bool _initialisedJets = false;
    bool _initialisedSpectra = false;

    Scatter2DPtr mult;
    Histo1DPtr _h_xp;
    Histo1DPtr _h_xi;
    Histo1DPtr _h_xe;
    Histo1DPtr _h_pTin;
    Histo1DPtr _h_pTout;
    Histo1DPtr _h_rapidityT;
    Histo1DPtr _h_rapidityS;
    Histo1DPtr _h_thrust;
    Histo1DPtr _h_heavyjetmass;
    Histo1DPtr _h_totaljetbroadening;
    Histo1DPtr _h_widejetbroadening;
    Histo1DPtr _h_cparameter;
    Histo1DPtr _h_thrustmajor;
    Histo1DPtr _h_thrustminor;
    Histo1DPtr _h_jetmassdifference;
    Histo1DPtr _h_aplanarity;
    Histo1DPtr _h_planarity;
    Histo1DPtr _h_oblateness;
    Histo1DPtr _h_sphericity;

    Histo1DPtr _h_R_Durham[6];
    Histo1DPtr _h_y_Durham[5];

    CounterPtr _weightedTotalChargedPartNum;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(ALEPH_2004_S5765862, ALEPH_2004_I636645);

}
