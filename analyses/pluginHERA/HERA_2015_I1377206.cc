// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief Measurement of sigma_red (F2) of H1 and ZEUS at different beam energies
  class HERA_2015_I1377206 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(HERA_2015_I1377206);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(Cuts::abseta < 5 && Cuts::pT > 100*MeV), "FS");
      declare(DISLepton(), "Lepton");
      declare(DISKinematics(), "Kinematics");
		
	
      Histo1DPtr dummy;
      string beamOpt = getOption<string>("BEAM","NONE");
      
      if (beamOpt == "NONE") {
        const ParticlePair& beam = beams();
        // cout << " beam id "<< beam.first.pid() << " " << beam.second.pid() << " sqrts " << sqrtS() << endl;
        if( beam.first.pid() == PID::POSITRON || beam.second.pid() == PID::POSITRON ) { 
          positron = true ;
        }
        else {
          positron = false ;
        }
      }
      else {
        if( beamOpt == "EMINUS" )
          positron = false;
        else if( beamOpt == "EPLUS" )
          positron = true;
        else {
          MSG_ERROR("Beam option error. You have specified an unsupported beam.");
          return;
        }
      }
        
      double eps = 0.01 ;
      // NC e+ p at sqrts=318
      if (isCompatibleWithSqrtS(318., eps) && positron  ) {
        // cout << " NC e+ p sqrts = " << sqrtS() << endl ;
        _h_sigred.add( 0.1,     0.15, book(dummy,1,1,1)); // Q2=0.15
        _h_sigred.add( 0.15,     0.2, book(dummy,1,1,2)); // Q2=0.2
        _h_sigred.add( 0.2,      0.3, book(dummy,1,1,3)); // Q2=0.25
        _h_sigred.add( 0.3,      0.4, book(dummy,1,1,4)); // Q2=0.35
        _h_sigred.add( 0.4,      0.45, book(dummy,1,1,5)); // Q2=0.4
        _h_sigred.add( 0.45,     0.6, book(dummy,1,1,6)); // Q2=0.5
        _h_sigred.add( 0.6,      0.7, book(dummy,1,1,7)); // Q2=0.65
        _h_sigred.add( 0.7,      1.0, book(dummy,1,1,8)); // Q2=0.85
        _h_sigred.add( 1.1,      1.3, book(dummy,1,1,9)); // Q2=1.2
        _h_sigred.add( 1.3,      1.7, book(dummy,1,1,10)); // Q2=1.5
        _h_sigred.add( 1.7,      2.3, book(dummy,1,1,11)); // Q2=2
        _h_sigred.add( 2.3,      3.1, book(dummy,1,1,12)); // Q2=2.7
        _h_sigred.add( 3.1,      3.8, book(dummy,1,1,13)); // Q2=3.5
        _h_sigred.add( 3.8,      5.3, book(dummy,1,1,14)); // Q2=4.5
        _h_sigred.add( 5.3,      8.0, book(dummy,1,1,15)); // Q2=6.5
        _h_sigred.add( 8.0,      9.1, book(dummy,1,1,16)); // Q2=8.5
        _h_sigred.add( 9.1,      11., book(dummy,1,1,17)); // Q2=10
        _h_sigred.add( 11.,      13., book(dummy,1,1,18)); // Q2=12
        _h_sigred.add( 13.,     17.4, book(dummy,1,1,19)); // Q2=15
        _h_sigred.add( 17.4,    19.1, book(dummy,1,1,20)); // Q2=18
        _h_sigred.add( 19.1,    25.8, book(dummy,1,1,21)); // Q2=22
        _h_sigred.add( 25.8,     28., book(dummy,1,1,22)); // Q2=27
        _h_sigred.add( 30.,      42., book(dummy,1,1,23)); // Q2=35
        _h_sigred.add( 42.,      49., book(dummy,1,1,24)); // Q2=45
        _h_sigred.add( 54.,      65., book(dummy,1,1,25)); // Q2=60
        _h_sigred.add( 65.,      75., book(dummy,1,1,26)); // Q2=70
        _h_sigred.add( 75.,     108., book(dummy,1,1,27)); // Q2=90
        _h_sigred.add( 108.,    134., book(dummy,1,1,28)); // Q2=120
        _h_sigred.add( 134.,    180., book(dummy,1,1,29)); // Q2=150
        _h_sigred.add( 180.,    225., book(dummy,1,1,30)); // Q2=200
        _h_sigred.add( 225.,    280., book(dummy,1,1,31)); // Q2=250
        _h_sigred.add( 280.,    325., book(dummy,1,1,32)); // Q2=300
        _h_sigred.add( 355.,    455., book(dummy,1,1,33)); // Q2=400
        _h_sigred.add( 460.,    545., book(dummy,1,1,34)); // Q2=500
        _h_sigred.add( 560.,    765., book(dummy,1,1,35)); // Q2=650
        _h_sigred.add( 770.,    835., book(dummy,1,1,36)); // Q2=800
        _h_sigred.add( 900.,   1120., book(dummy,1,1,37)); // Q2=1000
        _h_sigred.add( 1120.,  1295., book(dummy,1,1,38)); // Q2=1200
        _h_sigred.add( 1300.,  1755., book(dummy,1,1,39)); // Q2=1500
        _h_sigred.add( 1800.,  2270., book(dummy,1,1,40)); // Q2=2000
        _h_sigred.add( 2500.,  3685., book(dummy,1,1,41)); // Q2=3000
        _h_sigred.add( 4000.,  6520., book(dummy,1,1,42)); // Q2=5000
        _h_sigred.add( 7000.,  9275., book(dummy,1,1,43)); // Q2=8000
        _h_sigred.add( 10000.,15000., book(dummy,1,1,44)); // Q2=12000
        _h_sigred.add( 17000.,24770., book(dummy,1,1,45)); // Q2=20000
        _h_sigred.add( 25000.,42000., book(dummy,1,1,46)); // Q2=30000 
        // CC e+ p at sqrts=318
        // cout << " CC e+ p sqrts = " << sqrtS() << endl ;
        _h_sigred_cc.add( 280.,    325., book(dummy,6,1,1)); // Q2=300
        _h_sigred_cc.add( 460.,    545., book(dummy,6,1,2)); // Q2=500
        _h_sigred_cc.add( 900.,   1120., book(dummy,6,1,3)); // Q2=1000
        _h_sigred_cc.add( 1300.,  1755., book(dummy,6,1,4)); // Q2=1500
        _h_sigred_cc.add( 1800.,  2270., book(dummy,6,1,5)); // Q2=2000
        _h_sigred_cc.add( 2500.,  3685., book(dummy,6,1,6)); // Q2=3000
        _h_sigred_cc.add( 4000.,  6520., book(dummy,6,1,7)); // Q2=5000
        _h_sigred_cc.add( 7000.,  9275., book(dummy,6,1,8)); // Q2=8000
        _h_sigred_cc.add( 10000.,20000., book(dummy,6,1,9)); // Q2=15000
        _h_sigred_cc.add( 20000.,42000., book(dummy,6,1,10)); // Q2=30000
       }
      // NC e+ p at sqrts=300
      else if (isCompatibleWithSqrtS(300., eps) && positron  ) {
        // cout << " NC e+ p sqrts = " << sqrtS() << endl ;
         _h_sigred.add( 0.01,     0.05, book(dummy,2,1,1)); // Q2=0.045
        _h_sigred.add( 0.05,     0.07, book(dummy,2,1,2)); // Q2=0.065
        _h_sigred.add( 0.07,     0.09, book(dummy,2,1,3)); // Q2=0.085
        _h_sigred.add( 0.09,     0.12, book(dummy,2,1,4)); // Q2=0.11
        _h_sigred.add( 0.12,     0.18, book(dummy,2,1,5)); // Q2=0.15
        _h_sigred.add( 0.18,     0.22, book(dummy,2,1,6)); // Q2=0.2
        _h_sigred.add( 0.22,     0.32, book(dummy,2,1,7)); // Q2=0.25
        _h_sigred.add( 0.32,     0.4,  book(dummy,2,1,8)); // Q2=0.35
        _h_sigred.add( 0.4,      0.45, book(dummy,2,1,9)); // Q2=0.4
        _h_sigred.add( 0.45,     0.6, book(dummy,2,1,10)); // Q2=0.5
        _h_sigred.add( 0.6,      0.7, book(dummy,2,1,11)); // Q2=0.65
        _h_sigred.add( 0.7,      1.0, book(dummy,2,1,12)); // Q2=0.85
        _h_sigred.add( 1.1,      1.3, book(dummy,2,1,13)); // Q2=1.2
        _h_sigred.add( 1.3,      1.7, book(dummy,2,1,14)); // Q2=1.5
        _h_sigred.add( 1.7,      2.3, book(dummy,2,1,15)); // Q2=2
        _h_sigred.add( 2.3,      3.1, book(dummy,2,1,16)); // Q2=2.7
        _h_sigred.add( 3.1,      3.8, book(dummy,2,1,17)); // Q2=3.5
        _h_sigred.add( 3.8,      5.3, book(dummy,2,1,18)); // Q2=4.5
        _h_sigred.add( 5.3,      8.0, book(dummy,2,1,19)); // Q2=6.5
        _h_sigred.add( 8.0,      9.1, book(dummy,2,1,20)); // Q2=8.5
        _h_sigred.add( 9.1,      11., book(dummy,2,1,21)); // Q2=10
        _h_sigred.add( 11.,      13., book(dummy,2,1,22)); // Q2=12
        _h_sigred.add( 13.,     17.4, book(dummy,2,1,23)); // Q2=15
        _h_sigred.add( 17.4,    19.1, book(dummy,2,1,24)); // Q2=18
        _h_sigred.add( 19.1,    25.8, book(dummy,2,1,25)); // Q2=22
        _h_sigred.add( 25.8,     28., book(dummy,2,1,26)); // Q2=27
        _h_sigred.add( 30.,      42., book(dummy,2,1,27)); // Q2=35
        _h_sigred.add( 42.,      49., book(dummy,2,1,28)); // Q2=45
        _h_sigred.add( 54.,      65., book(dummy,2,1,29)); // Q2=60
        _h_sigred.add( 65.,      75., book(dummy,2,1,30)); // Q2=70
        _h_sigred.add( 75.,     108., book(dummy,2,1,31)); // Q2=90
        _h_sigred.add( 108.,    134., book(dummy,2,1,32)); // Q2=120
        _h_sigred.add( 134.,    180., book(dummy,2,1,33)); // Q2=150
        _h_sigred.add( 180.,    225., book(dummy,2,1,34)); // Q2=200
        _h_sigred.add( 225.,    280., book(dummy,2,1,35)); // Q2=250
        _h_sigred.add( 280.,    325., book(dummy,2,1,36)); // Q2=300
        _h_sigred.add( 355.,    455., book(dummy,2,1,37)); // Q2=400
        _h_sigred.add( 460.,    545., book(dummy,2,1,38)); // Q2=500
        _h_sigred.add( 560.,    765., book(dummy,2,1,39)); // Q2=650
        _h_sigred.add( 770.,    835., book(dummy,2,1,40)); // Q2=800
        _h_sigred.add( 900.,   1120., book(dummy,2,1,41)); // Q2=1000
        _h_sigred.add( 1120.,  1295., book(dummy,2,1,42)); // Q2=1200
        _h_sigred.add( 1300.,  1755., book(dummy,2,1,43)); // Q2=1500
        _h_sigred.add( 1800.,  2270., book(dummy,2,1,44)); // Q2=2000
        _h_sigred.add( 2500.,  3685., book(dummy,2,1,45)); // Q2=3000
        _h_sigred.add( 4000.,  6520., book(dummy,2,1,46)); // Q2=5000
        _h_sigred.add( 7000.,  9275., book(dummy,2,1,47)); // Q2=8000
        _h_sigred.add( 10000.,15000., book(dummy,2,1,48)); // Q2=12000
        _h_sigred.add( 17000.,24770., book(dummy,2,1,49)); // Q2=20000
        _h_sigred.add( 25000.,42000., book(dummy,2,1,50)); // Q2=30000
      }
      else if (isCompatibleWithSqrtS(251., eps) && positron  ) {
        // NC e+ p at sqrts=251
        // cout << " NC e+ p sqrts = " << sqrtS() << endl ;
        _h_sigred.add( 1.0,      1.7, book(dummy,3,1,1)); // Q2=1.5
        _h_sigred.add( 1.7,      2.3, book(dummy,3,1,2)); // Q2=2
        _h_sigred.add( 2.3,      3.1, book(dummy,3,1,3)); // Q2=2.5
        _h_sigred.add( 3.1,      3.8, book(dummy,3,1,4)); // Q2=3.5
        _h_sigred.add( 3.8,      5.3, book(dummy,3,1,5)); // Q2=5
        _h_sigred.add( 5.3,      8.0, book(dummy,3,1,6)); // Q2=6.5
        _h_sigred.add( 8.0,      9.1, book(dummy,3,1,7)); // Q2=8.5
        _h_sigred.add( 11.,      13., book(dummy,3,1,8)); // Q2=12
        _h_sigred.add( 13.,     17.4, book(dummy,3,1,9)); // Q2=15
        _h_sigred.add( 17.4,    22.1, book(dummy,3,1,10)); // Q2=20
        _h_sigred.add( 22.1,    28. , book(dummy,3,1,11)); // Q2=25
        _h_sigred.add( 30.,      42., book(dummy,3,1,12)); // Q2=35
        _h_sigred.add( 42.,      49., book(dummy,3,1,13)); // Q2=45
        _h_sigred.add( 54.,      65., book(dummy,3,1,14)); // Q2=60
        _h_sigred.add( 75.,     108., book(dummy,3,1,15)); // Q2=90
        _h_sigred.add( 108.,    134., book(dummy,3,1,16)); // Q2=120
        _h_sigred.add( 134.,    180., book(dummy,3,1,17)); // Q2=150
        _h_sigred.add( 180.,    225., book(dummy,3,1,18)); // Q2=200
        _h_sigred.add( 225.,    280., book(dummy,3,1,19)); // Q2=250
        _h_sigred.add( 280.,    325., book(dummy,3,1,20)); // Q2=300
        _h_sigred.add( 355.,    455., book(dummy,3,1,21)); // Q2=400
        _h_sigred.add( 460.,    545., book(dummy,3,1,22)); // Q2=500
        _h_sigred.add( 560.,    765., book(dummy,3,1,23)); // Q2=650
        _h_sigred.add( 770.,    835., book(dummy,3,1,24)); // Q2=800
      }
      else if (isCompatibleWithSqrtS(225., eps) && positron  ) {
        // NC e+ p at sqrts=225
        // cout << " NC e+ p sqrts = " << sqrtS() << endl ;
        _h_sigred.add( 1.0,      1.7, book(dummy,4,1,1)); // Q2=1.5
        _h_sigred.add( 1.7,      2.3, book(dummy,4,1,2)); // Q2=2
        _h_sigred.add( 2.3,      3.1, book(dummy,4,1,3)); // Q2=2.5
        _h_sigred.add( 3.1,      3.8, book(dummy,4,1,4)); // Q2=3.5
        _h_sigred.add( 3.8,      5.3, book(dummy,4,1,5)); // Q2=5
        _h_sigred.add( 5.3,      8.0, book(dummy,4,1,6)); // Q2=6.5
        _h_sigred.add( 8.0,      9.1, book(dummy,4,1,7)); // Q2=8.5
        _h_sigred.add( 11.,      13., book(dummy,4,1,8)); // Q2=12
        _h_sigred.add( 13.,     17.4, book(dummy,4,1,9)); // Q2=15
        _h_sigred.add( 17.4,    22.1, book(dummy,4,1,10)); // Q2=20
        _h_sigred.add( 22.1,    28. , book(dummy,4,1,11)); // Q2=25
        _h_sigred.add( 30.,      42., book(dummy,4,1,12)); // Q2=35
        _h_sigred.add( 42.,      49., book(dummy,4,1,13)); // Q2=45
        _h_sigred.add( 54.,      65., book(dummy,4,1,14)); // Q2=60
        _h_sigred.add( 75.,     108., book(dummy,4,1,15)); // Q2=90
        _h_sigred.add( 108.,    134., book(dummy,4,1,16)); // Q2=120
        _h_sigred.add( 134.,    180., book(dummy,4,1,17)); // Q2=150
        _h_sigred.add( 180.,    225., book(dummy,4,1,18)); // Q2=200
        _h_sigred.add( 225.,    280., book(dummy,4,1,19)); // Q2=250
        _h_sigred.add( 280.,    325., book(dummy,4,1,20)); // Q2=300
        _h_sigred.add( 355.,    455., book(dummy,4,1,21)); // Q2=400
        _h_sigred.add( 460.,    545., book(dummy,4,1,22)); // Q2=500
        _h_sigred.add( 560.,    765., book(dummy,4,1,23)); // Q2=650
        _h_sigred.add( 770.,    835., book(dummy,4,1,24)); // Q2=800
      }
      else if (isCompatibleWithSqrtS(225., eps) && positron  ) {
        // NC e- p at sqrts=318
        // cout << " NC e- p sqrts = " << sqrtS() << endl ;
        _h_sigred.add( 54.,      65., book(dummy,5,1,1)); // Q2=60
        _h_sigred.add( 75.,     108., book(dummy,5,1,2)); // Q2=90
        _h_sigred.add( 108.,    134., book(dummy,5,1,3)); // Q2=120
        _h_sigred.add( 134.,    180., book(dummy,5,1,4)); // Q2=150
        _h_sigred.add( 180.,    225., book(dummy,5,1,5)); // Q2=200
        _h_sigred.add( 225.,    280., book(dummy,5,1,6)); // Q2=250
        _h_sigred.add( 280.,    325., book(dummy,5,1,7)); // Q2=300
        _h_sigred.add( 355.,    455., book(dummy,5,1,8)); // Q2=400
        _h_sigred.add( 460.,    545., book(dummy,5,1,9)); // Q2=500
        _h_sigred.add( 560.,    765., book(dummy,5,1,10)); // Q2=650
        _h_sigred.add( 770.,    835., book(dummy,5,1,11)); // Q2=800
        _h_sigred.add( 900.,   1120., book(dummy,5,1,12)); // Q2=1000
        _h_sigred.add( 1120.,  1295., book(dummy,5,1,13)); // Q2=1200
        _h_sigred.add( 1300.,  1755., book(dummy,5,1,14)); // Q2=1500
        _h_sigred.add( 1800.,  2270., book(dummy,5,1,15)); // Q2=2000
        _h_sigred.add( 2500.,  3685., book(dummy,5,1,16)); // Q2=3000
        _h_sigred.add( 4000.,  6520., book(dummy,5,1,17)); // Q2=5000
        _h_sigred.add( 7000.,  9275., book(dummy,5,1,18)); // Q2=8000
        _h_sigred.add( 10000.,15000., book(dummy,5,1,19)); // Q2=12000
        _h_sigred.add( 17000.,24770., book(dummy,5,1,20)); // Q2=20000
        _h_sigred.add( 25000.,42000., book(dummy,5,1,21)); // Q2=30000
        _h_sigred.add( 42000.,70000., book(dummy,5,1,22)); // Q2=50000
        // CC e- p at sqrts=318
        // cout << " CC e- p sqrts = " << sqrtS() << endl ;
        _h_sigred_cc.add( 280.,    325., book(dummy,7,1,1)); // Q2=300
        _h_sigred_cc.add( 460.,    545., book(dummy,7,1,2)); // Q2=500
        _h_sigred_cc.add( 900.,   1120., book(dummy,7,1,3)); // Q2=1000
        _h_sigred_cc.add( 1300.,  1755., book(dummy,7,1,4)); // Q2=1500
        _h_sigred_cc.add( 1800.,  2270., book(dummy,7,1,5)); // Q2=2000
        _h_sigred_cc.add( 2500.,  3685., book(dummy,7,1,6)); // Q2=3000
        _h_sigred_cc.add( 4000.,  6520., book(dummy,7,1,7)); // Q2=5000
        _h_sigred_cc.add( 7000.,  9275., book(dummy,7,1,8)); // Q2=8000
        _h_sigred_cc.add( 10000.,20000., book(dummy,7,1,9)); // Q2=15000
        _h_sigred_cc.add( 20000.,42000., book(dummy,7,1,10)); // Q2=30000

      }

/*
        book(_hist_Q2_10, "Q2_10",100,1., 11.0);
        book(_hist_Q2_100, "Q2_100",100,10., 100.0);
        book(_hist_Q2_1000, "Q2_1000",100,100., 1000.0);
        book(_hist_Q2_2000, "Q2_2000",100,800., 5000.0);
        book(_hist_Q2_3000, "Q2_3000",100,3000., 10000.0);
*/


    }


    void analyze(const Event& event) {

      /// @todo Do the event by event analysis here
      const DISKinematics& dk = applyProjection<DISKinematics>(event, "Kinematics");
      const DISLepton& dl = applyProjection<DISLepton>(event,"Lepton");

      // Get the DIS kinematics
      double x  = dk.x();
      double y = dk.y();
      double Q2 = dk.Q2()/GeV;
	
      // Flux factor
      const double alpha = 7.29927e-3;
      // GF = 1.16638e-5 Fermi constant
      const double GF2 = 1.16638e-5*1.16638e-5;
      // MW = 80.385 W-boson mass
      const double MW2 = 80.385 * 80.385;
/*
	_hist_Q2_10-> fill(Q2) ;
	_hist_Q2_100-> fill(Q2) ;
	_hist_Q2_1000-> fill(Q2) ;
	_hist_Q2_2000-> fill(Q2) ;
	_hist_Q2_3000-> fill(Q2) ;
*/

      if (PID::isNeutrino(dl.out().abspid()) ) {
        // fill histo for CC      
        double F = 2.0*M_PI*x/GF2 * pow((MW2 + Q2)/MW2,2);
        _h_sigred_cc.fill(Q2,x,F); // fill histogram x,Q2 
      }
      else { 
        // fill histo for NC      
        double F = x*pow(Q2,2.)/(2.0*M_PI*pow(alpha,2.)*(1.0+pow((1.-y),2.)));
        _h_sigred.fill(Q2,x,F); // fill histogram x,Q2 
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double gev2nb =0.389e6;
      const double scalefactor=crossSection()/nanobarn/sumOfWeights()/gev2nb ;
      // with _h_sigred.scale also q2 bin width is scaled
      _h_sigred.scale(scalefactor, this);
      _h_sigred_cc.scale(scalefactor, this);
    }

    //@}


    /// @name Histograms
    //@{
      BinnedHistogram _h_sigred, _h_sigred_cc;
      Histo1DPtr _hist_Q2_10,_hist_Q2_100,_hist_Q2_1000,_hist_Q2_2000,_hist_Q2_3000;
      bool positron ;
    //@}


  };

  RIVET_DECLARE_PLUGIN(HERA_2015_I1377206);
}
