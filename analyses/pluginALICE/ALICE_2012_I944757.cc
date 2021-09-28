#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Charm production at central rapidity in pp at 7 TeV
  class ALICE_2012_I944757 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2012_I944757);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(Cuts::absrap < 0.5), "UFS");

      // Book histograms
      book(_h_D0,     1, 1, 1);
      book(_h_Dplus,  2, 1, 1);
      book(_h_Dstarp, 3, 1, 1);
      book(_h_integ,  4, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

        /*PDG code IDs used inside the for loop: 421 = D0, 411 = D+, 413 = D*+ */

        for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
          if (p.fromBottom())  continue;
          if (p.abspid() == 421) {
            _h_D0->fill(p.pT()/GeV);
            _h_integ->fill(1);
          }
          else if (p.abspid() == 411) {
            _h_Dplus->fill(p.pT()/GeV);
            _h_integ->fill(2);
          }
          else if (p.abspid()== 413) {
            _h_Dstarp->fill(p.pT()/GeV);
            _h_integ->fill(3);
          }
        }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_D0, crossSection()/(microbarn*2*sumOfWeights())); // norm to cross section
      scale(_h_Dplus, crossSection()/(microbarn*2*sumOfWeights())); // norm to cross section
      scale(_h_Dstarp, crossSection()/(microbarn*2*sumOfWeights())); // norm to cross section
      scale(_h_integ, crossSection()/(microbarn*2*sumOfWeights())); // norm to cross section
      /* Obtained cross sections data at this point consider both particles and antiparticles
      hence the added factor 2 in the normalization solves the issue (as done in the paper) */
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_D0, _h_Dplus, _h_Dstarp, _h_integ;
    //@}


  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2012_I944757);
}
