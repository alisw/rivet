// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Math/Constants.hh"
#include "Rivet/Math/Units.hh"

namespace Rivet {


  class LHCB_2010_S8758301 : public Analysis {
  public:

  	const double MIN_PT = 0.0001; // [GeV/c]


    /// Constructor
    LHCB_2010_S8758301()
      : Analysis("LHCB_2010_S8758301"), _mode(0),
      sumKs0_all(0),
      sumKs0_outup(0), sumKs0_outdwn(0),
      sum_low_pt_loss(0), sum_high_pt_loss(0)
    {  }


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
    	// only interested in Ks0 particles
      declare(UnstableParticles(Cuts::pid == 310), "UFS");

      if (getOption("DropPlots") == "SIG") _mode = 1;
    	if (getOption("DropPlots") == "D2SIG") _mode = 2;

    	if ((_mode & 1) == 0) {	//SIG
    		book(_h_K0s_pt_30    ,1,1,1);
    		book(_h_K0s_pt_35    ,1,1,2);
    		book(_h_K0s_pt_40    ,1,1,3);
    	};

    	if ((_mode & 2) == 0) { //D2SIG
    		book(_h_K0s_pt_y_30  ,2,1,1);
    		book(_h_K0s_pt_y_35  ,2,1,2);
    		book(_h_K0s_pt_y_40  ,2,1,3);
    	};

      book(_h_K0s_pt_y_all ,3,1,1);

      book(sumKs0_30, "TMP/sumKs0_30");
      book(sumKs0_35, "TMP/sumKs0_35");
      book(sumKs0_40, "TMP/sumKs0_40");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      // safe to do this as container has only Ks0
      sumKs0_all += ufs.particles().size();

      for (const Particle& p : ufs.particles()) {

        double y = p.absrapidity(); // symmetric LHCb (factor 1/2 in finalize!)
        double pT = p.pT();

        if (pT < MIN_PT) {
          sum_low_pt_loss ++; // just for debug since no inferior limit on pT
          MSG_DEBUG("Small pT K^0_S: " << pT << " GeV/c.");
        }
        if (pT > 1.6) {
          sum_high_pt_loss ++;
          continue; // no need to flow overflow bin
        }
        if ((y > 2.5) && (y < 4.0)) {
          _h_K0s_pt_y_all->fill(pT);
          if ((y > 2.5) && (y < 3.0)) {
            if ((_mode & 2) == 0) _h_K0s_pt_y_30->fill(pT);
            if ((_mode & 1) == 0) _h_K0s_pt_30->fill(pT);
            sumKs0_30->fill();
          } else if ((y > 3.0) && (y < 3.5)) {
            if ((_mode & 2) == 0) _h_K0s_pt_y_35->fill(pT);
            if ((_mode & 1) == 0) _h_K0s_pt_35->fill(pT);
            sumKs0_35->fill();
          } else if ((y > 3.5) && (y < 4.0)) {
            if ((_mode & 2) == 0) _h_K0s_pt_y_40->fill(pT);
            if ((_mode & 1) == 0) _h_K0s_pt_40->fill(pT);
            sumKs0_40->fill();
          }
        } else if (y < 2.5) {
          sumKs0_outdwn ++;
        } else if (y > 4.0) {
          sumKs0_outup ++;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      MSG_DEBUG("Total number Ks0: " << sumKs0_all << endl
                << "Sum of weights: " << sumOfWeights() << endl
                << "Weight Ks0 (2.5 < y < 3.0): " <<  sumKs0_30 ->sumW()<< endl
                << "Weight Ks0 (3.0 < y < 3.5): " << sumKs0_35->sumW() << endl
                << "Weight Ks0 (3.5 < y < 4.0): " << sumKs0_40->sumW() << endl
                << "Nb. Ks0 (y > 4.0): " << sumKs0_outup << endl
                << "Nb. Ks0 (y < 2.5): " << sumKs0_outdwn << endl
                << "Nb. Ks0 (pT < " << (MIN_PT/MeV) << " MeV/c): " << sum_low_pt_loss << endl
                << "Nb. Ks0 (pT > 1.6 GeV/c): " << sum_high_pt_loss << endl
                << "Cross-section [mb]: " << crossSection()/millibarn << endl
                << "Nb. events: " << numEvents());
      // Compute cross-section; multiply by bin width for correct scaling
      // cross-section given by Rivet in pb (symmetric LHCb!)
      double xsection_factor = 0.5 * crossSection()/sumOfWeights();
      // Divide by bin area for consistent scaling, xsection in mub
      scale(_h_K0s_pt_30, 0.1*xsection_factor/microbarn);
      scale(_h_K0s_pt_35, 0.1*xsection_factor/microbarn);
      scale(_h_K0s_pt_40, 0.1*xsection_factor/microbarn);
      // Already divided by dy (rapidity window width), xsection in mb
      scale(_h_K0s_pt_y_30, xsection_factor/millibarn);
      scale(_h_K0s_pt_y_35, xsection_factor/millibarn);
      scale(_h_K0s_pt_y_40, xsection_factor/millibarn);
      // factorize to phase-space volume (area)
      scale(_h_K0s_pt_y_all, xsection_factor/1.5/1.6/millibarn);
    }

    /// @}


    /// Mode switch
    size_t _mode;

    /// @name Histograms
    /// @{
    Histo1DPtr _h_K0s_pt_y_30;         // histogram for 2.5 < y < 3.0 (d2sigma)
    Histo1DPtr _h_K0s_pt_y_35;         // histogram for 3.0 < y < 3.5 (d2sigma)
    Histo1DPtr _h_K0s_pt_y_40;         // histogram for 3.5 < y < 4.0 (d2sigma)
    Histo1DPtr _h_K0s_pt_30;           // histogram for 2.5 < y < 3.0 (sigma)
    Histo1DPtr _h_K0s_pt_35;           // histogram for 3.0 < y < 3.5 (sigma)
    Histo1DPtr _h_K0s_pt_40;           // histogram for 3.5 < y < 4.0 (sigma)
    Histo1DPtr _h_K0s_pt_y_all;        // histogram for 2.5 < y < 4.0 (d2sigma)
    CounterPtr sumKs0_30;                           // Sum of weights 2.5 < y < 3.0
    CounterPtr sumKs0_35;                           // Sum of weights 3.0 < y < 3.5
    CounterPtr sumKs0_40;                           // Sum of weights 3.5 < y < 4.0
    // Various counters mainly for debugging and comparisons between different generators
    size_t sumKs0_all;                          // Nb of all Ks0 generated
    size_t sumKs0_outup;                        // Nb of mesons with y > 4.0
    size_t sumKs0_outdwn;                       // Nb of mesons with y < 2.5
    size_t sum_low_pt_loss;                     // Nb of mesons with very low pT (indicates when units are mixed-up)
    size_t sum_high_pt_loss;                    // Nb of mesons with pT > 1.6 GeV/c
    /// @}
  };


  RIVET_DECLARE_ALIASED_PLUGIN(LHCB_2010_S8758301, LHCB_2010_I865584);

}
