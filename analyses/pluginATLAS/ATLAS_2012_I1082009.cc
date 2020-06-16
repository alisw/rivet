// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  class ATLAS_2012_I1082009 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    ATLAS_2012_I1082009()
      : Analysis("ATLAS_2012_I1082009")
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Input for the jets: No neutrinos, no muons
      VetoedFinalState veto;
      veto.addVetoPairId(PID::MUON);
      veto.vetoNeutrinos();
      FastJets jets(veto, FastJets::ANTIKT, 0.6);
      declare(jets, "jets");
      // unstable final-state for D*
      declare(UnstableParticles(), "UFS");

      book(_weight25_30, "_weight_25_30");
      book(_weight30_40, "_weight_30_40");
      book(_weight40_50, "_weight_40_50");
      book(_weight50_60, "_weight_50_60");
      book(_weight60_70, "_weight_60_70");
      book(_weight25_70, "_weight_25_70");

      book(_h_pt25_30 , 8,1,1);
      book(_h_pt30_40 , 9,1,1);
      book(_h_pt40_50 ,10,1,1);
      book(_h_pt50_60 ,11,1,1);
      book(_h_pt60_70 ,12,1,1);
      book(_h_pt25_70 ,13,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // get the jets
      Jets jets;
      for (const Jet& jet : apply<FastJets>(event, "jets").jetsByPt(25.0*GeV)) {
        if ( jet.abseta() < 2.5 ) jets.push_back(jet);
      }
      // get the D* mesons
      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");
      Particles Dstar;
      for (const Particle& p : ufs.particles()) {
        const int id = p.abspid();
        if(id==413) Dstar.push_back(p);
      }

      // loop over the jobs
      for (const Jet& jet : jets ) {
        double perp = jet.perp();
        bool found = false;
        double z(0.);
        if(perp<25.||perp>70.) continue;
        for(const Particle & p : Dstar) {
          if(p.perp()<7.5) continue;
          if(deltaR(p, jet.momentum())<0.6) {
            Vector3 axis = jet.p3().unit();
            z = axis.dot(p.p3())/jet.E();
            if(z<0.3) continue;
            found = true;
            break;
          }
        }
        _weight25_70->fill();
        if(found) _h_pt25_70->fill(z);
        if(perp>=25.&&perp<30.) {
          _weight25_30->fill();
          if(found) _h_pt25_30->fill(z);
        }
        else if(perp>=30.&&perp<40.) {
          _weight30_40->fill();
          if(found) _h_pt30_40->fill(z);
        }
        else if(perp>=40.&&perp<50.) {
          _weight40_50->fill();
          if(found) _h_pt40_50->fill(z);
        }
        else if(perp>=50.&&perp<60.) {
          _weight50_60->fill();
          if(found) _h_pt50_60->fill(z);
        }
        else if(perp>=60.&&perp<70.) {
          _weight60_70->fill();
          if(found) _h_pt60_70->fill(z);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_pt25_30,1./ *_weight25_30);
      scale(_h_pt30_40,1./ *_weight30_40);
      scale(_h_pt40_50,1./ *_weight40_50);
      scale(_h_pt50_60,1./ *_weight50_60);
      scale(_h_pt60_70,1./ *_weight60_70);
      scale(_h_pt25_70,1./ *_weight25_70);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    CounterPtr _weight25_30,_weight30_40,_weight40_50;
    CounterPtr _weight50_60,_weight60_70,_weight25_70;

    Histo1DPtr _h_pt25_30;
    Histo1DPtr _h_pt30_40;
    Histo1DPtr _h_pt40_50;
    Histo1DPtr _h_pt50_60;
    Histo1DPtr _h_pt60_70;
    Histo1DPtr _h_pt25_70;
    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1082009);

}
