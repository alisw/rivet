// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/RivetHTT.hh"

namespace Rivet {


  /// Example of how to execute the HEPTopTagger interface
  class EXAMPLE_HTT : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(EXAMPLE_HTT);


    /// Book histograms and initialise projections before the run
    void init() {

      // Declare the fat jet
      FastJets fatjets(FinalState(Cuts::abseta < 2.5), FastJets::CA, 1.5,
                       JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(fatjets, "FatJets");

      // Declare histograms
      book(_h_deltaMtop, "DeltaMtop", 20, 0, 40.);
      book(_h_Mpruned, "Mpruned", 40, 0.0, 400.0);
      book(_h_Munfiltered, "Munfiltered", 40, 0.0, 400.0);
      book(_h_Mcand, "Mcand", 40, 0.0, 400.0);
      book(_h_DMpruned, "DMpruned", 20, 0.0, 200.0);
      book(_h_DMunfiltered, "DMunfiltered", 20, 0.0, 200.0);

      // Set parameters for HTT
      HTT::InputParameters parameters;
      parameters.fw = getOption<double>("fw", 0.15);
      parameters.mass_drop = getOption<double>("mass_drop", 0.8);
      parameters.max_subjet_mass = getOption<double>("max_subjet_mass", 30)*GeV;
      parameters.filt_N = getOption<unsigned int>("filt_N", 5);

      parameters.Mtop_min = getOption<double>("Mtop_min", 150.)*GeV;
      parameters.Mtop_max = getOption<double>("Mtop_max", 200.)*GeV;

      parameters.mass_ratio_range_min = (1 - parameters.fw) * 80.379 / 172.9;
      parameters.mass_ratio_range_max = (1 + parameters.fw) * 80.379 / 172.9;

      parameters.m23cut = getOption<double>("m23cut", 0.35);
      parameters.m13cutmin = getOption<double>("m13cutmin", 0.2);
      parameters.m13cutmax = getOption<double>("m13cutmax", 1.3);

      parameters.prune_zcut = getOption<double>("prune_zcut", 0.1);
      parameters.prune_rcut = getOption<double>("prune_rcut", .5);

      _tagger.setParams(parameters);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      Jets fatjets = apply<JetAlg>(event, "FatJets").jetsByPt();
      if (fatjets.size() > 0 && fatjets[0].pT() > 200*GeV) {
        _tagger.calc(fatjets[0]);

        // Fill the histograms if the Jet is tagged
        if (_tagger.isTopTagged()) {
            _h_deltaMtop->fill(_tagger.deltaTopMass());
            _h_Mpruned->fill(_tagger.prunedMass());
            _h_Munfiltered->fill(_tagger.unfilteredMass());
            _h_Mcand->fill(_tagger.topJet().mass());
            _h_DMpruned->fill(fabs(_tagger.prunedMass() - _tagger.topJet().mass()));
            _h_DMunfiltered->fill(fabs(_tagger.unfilteredMass() - _tagger.topJet().mass()));
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale({_h_deltaMtop, _h_Mpruned, _h_Munfiltered, _h_Mcand, _h_DMpruned, _h_DMunfiltered}, crossSection()/picobarn/sumW());
    }


  private:

    // Tagger class
    HTT _tagger;

    // HTT distributions
    Histo1DPtr _h_Mpruned, _h_Munfiltered, _h_Mcand, _h_DMpruned, _h_DMunfiltered, _h_deltaMtop;

  };


  RIVET_DECLARE_PLUGIN(EXAMPLE_HTT);

}
