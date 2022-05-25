// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DISKinematics.hh"

namespace Rivet {


  /// @brief   F2 structure function in DIS e+ p (ZEUS)
  
  class ZEUS_1996_I420332 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ZEUS_1996_I420332);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      // Initialise and register projections
      declare(DISLepton(), "Lepton");
      declare(DISKinematics(), "Kinematics");
		
	
       Histo1DPtr dummy;
	 _h_f2.add( 3.2,    4.0,book(dummy,1,1,1)); 
	 _h_f2.add( 4.0,    5.0,book(dummy,2,1,1));
	 _h_f2.add( 5.0,    7.0,book(dummy,3,1,1));
	 _h_f2.add( 7.0,    9.0,book(dummy,5,1,1));
	 _h_f2.add( 9.0,    11.0,book(dummy,6,1,1));
	 _h_f2.add( 11.0,   13.0,book(dummy,8,1,1));
	 _h_f2.add( 13.0,   16.0,book(dummy,9,1,1));
	 _h_f2.add( 16.0,   20.0,book(dummy,11,1,1));
	 _h_f2.add( 20.,    32,  book(dummy,13,1,1));
	 _h_f2.add( 32.,    40., book(dummy,15,1,1));
	 _h_f2.add( 40.,    50., book(dummy,17,1,1));
	 _h_f2.add( 50.,    65., book(dummy,18,1,1));
	 _h_f2.add( 65.,    85., book(dummy,19,1,1));
	 _h_f2.add( 85.,    110.,book(dummy,20,1,1));
	 _h_f2.add( 110.,   140.,book(dummy,21,1,1));
	 _h_f2.add( 140.,   185.,book(dummy,22,1,1));
	 _h_f2.add( 185.,   240.,book(dummy,23,1,1));
	 _h_f2.add( 240.,   310.,book(dummy,24,1,1));
	 _h_f2.add( 310.,   410.,book(dummy,25,1,1));
	 _h_f2.add( 410.,   530.,book(dummy,26,1,1));
	 _h_f2.add( 530.,   710.,book(dummy,27,1,1));
	 _h_f2.add( 710.,   900.,book(dummy,28,1,1));
	 _h_f2.add( 900.,  1300.,book(dummy,29,1,1));
	 _h_f2.add( 1300., 1800.,book(dummy,30,1,1));
	 _h_f2.add( 1800., 2500.,book(dummy,31,1,1));
	 _h_f2.add( 2500., 3500.,book(dummy,32,1,1));
	 _h_f2.add( 3500., 15000.,book(dummy,33,1,1));

/*
        book(_hist_Q2_10, "Q2_10",100,1., 11.0);
        book(_hist_Q2_100, "Q2_100",100,10., 100.0);
        book(_hist_Q2_1000, "Q2_1000",100,100., 1000.0);
        book(_hist_Q2_2000, "Q2_2000",100,800., 5000.0);
        book(_hist_Q2_3000, "Q2_3000",100,3000., 10000.0);
*/

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const DISKinematics& dk = applyProjection<DISKinematics>(event, "Kinematics");
      //const DISLepton& dl = applyProjection<DISLepton>(event,"Lepton");

      // Get the DIS kinematics
      double x  = dk.x();
      double y = dk.y();
      double Q2 = dk.Q2()/GeV;
	
	// Flux factor
	const double alpha = 7.29927e-3;
	double F = x*pow(Q2,2.)/(2.0*M_PI*pow(alpha,2.)*(1.0+pow((1.-y),2.)));
/*
	_hist_Q2_10-> fill(Q2) ;
	_hist_Q2_100-> fill(Q2) ;
	_hist_Q2_1000-> fill(Q2) ;
	_hist_Q2_2000-> fill(Q2) ;
	_hist_Q2_3000-> fill(Q2) ;
*/
	_h_f2.fill(Q2,x,F); // wypelniamy histogram x w skali Q2


    }


    /// Normalise histograms etc., after the run
    void finalize() {
	double gev2nb =0.389e6;
      double scalefactor=crossSection()/nanobarn/sumOfWeights()/gev2nb ;
      // with _h_f2.scale also q2 bin width is scaled
      _h_f2.scale(scalefactor, this);

    }

    ///@}


    /// @name Histograms
    ///@{
	BinnedHistogram _h_f2;
      Histo1DPtr _hist_Q2_10,_hist_Q2_100,_hist_Q2_1000,_hist_Q2_2000,_hist_Q2_3000;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(ZEUS_1996_I420332);

}
