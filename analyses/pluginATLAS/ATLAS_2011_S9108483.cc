// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Math/Constants.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/NonHadronicFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  class ATLAS_2011_S9108483 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2011_S9108483);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // get the non-hadronic final-state particles
      double etaMax = 2.5;
      const NonHadronicFinalState nhfs((Cuts::etaIn(-etaMax,etaMax) && Cuts::pT >= 13.*GeV));
      declare(nhfs,"NHFS");
      // select the charged ones
      const ChargedFinalState cfs(nhfs);
      declare(cfs,"CFS");
      // and then veto electrons, and taus to be safe
      VetoedFinalState vfs(cfs);
      vfs.addVetoPairId(PID::ELECTRON);

      declare(vfs,"VFS");

      /// Book histograms
      book(_count_trigger   ,"count_trigger"  , 1, 0., 1.);
      book(_count_event     ,"count_selection", 1, 0., 1.);
      book(_count_quality   ,"count_quality"  , 1, 0., 1.);
      book(_count_beta      ,"count_beta"     , 1, 0., 1.);
      book(_count_90  ,"count_90" , 1, 0., 1.);
      book(_count_110 ,"count_110", 1, 0., 1.);
      book(_count_120 ,"count_120", 1, 0., 1.);
      book(_count_130 ,"count_130", 1, 0., 1.);

      book(_hist_beta ,"beta",1000, 0.,   2.);
      book(_hist_time ,"time",1000, -50,  50.);
      book(_hist_mass ,"mass",  60, 5., 305.);
    }


    double rndGauss(double sigma, double mean) {
      double r = sqrt(-2.0*log(rand()/static_cast<double>(RAND_MAX)));
      double phi = rand()/static_cast<double>(RAND_MAX)*2.0*pi;
      return mean+sigma*r*sin(phi);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // smearing parameters
      // time measurement (in ns)
//       const double tsmear=5.*0.7;
      const double tsmear=0.7;
      // sagita error
      const double csag  =1.1E-4;
      // multiple scattering
      const double cms   =2.0E-2;
      // muon chamber radius (in metres)
      const double radius = 10.e3;
      // convert to ns
      const double tr = radius/c_light;
      // get the charged final-state particles
      Particles charged =
        apply<VetoedFinalState>(event,"VFS").particles();
      // need at least two candidates
      if(charged.size()<2) vetoEvent;
      // number passing trigger
      _count_trigger->fill(0.5);
      // Z mass veto
      for ( const Particle & mu1 : charged ) {
        for ( const Particle & mu2 : charged ) {
          double mass = (mu1.momentum()+mu2.momentum()).mass();
          double diff = abs(mass-91.18);
          if(diff<10.) vetoEvent;
        }
      }
      // number passing first event selection
      _count_event->fill(0.5);
      // now find the candidates
      // loop over the particles and find muons and heavy charged particles
      map<double,Particle> muonCandidates;
      for (const Particle& mu : charged) {
        // calculate the smeared momentum
        double pT     = mu.pT2();
        double pmag   = sqrt(pT+sqr(mu.pz()));
        double deltap = sqrt( sqr(csag*sqr(pmag)) + sqr(cms*mu.E()/GeV));
        double psmear = rndGauss(deltap,pmag);
        // keep particles with pT>40
        if(psmear/pmag*mu.perp()<40.*GeV||
           psmear/pmag*mu.perp()>1000.*GeV) continue;
        muonCandidates.insert(make_pair(psmear,mu));
      }
      // require two candidates
      if(muonCandidates.size()<2) vetoEvent;
      // number passing "quality" cut
      _count_quality->fill(0.5);
      // now do the time of flight
      bool filled = false;
      for(map<double,Particle>::const_iterator it=muonCandidates.begin();
          it!=muonCandidates.end();++it) {
        // true magnitude and pT of momentum
        double pT     = it->second.pT2();
        double pmag   = sqrt(pT+sqr(it->second.pz()));
        pT = sqrt(pT);
        // true time difference in ns
        double deltaT  =tr *(it->second.E()-pmag)/pT;
        // smear it
        deltaT = rndGauss(tsmear,deltaT);
        // beta
        double beta = 1./(1.+deltaT/tr*pT/pmag);
        _hist_beta->fill(beta);
        _hist_time->fill(deltaT);
        // beta cut
        if(beta<0.95) continue;
        // mass
        double mass = 2.*pT*it->first*deltaT/tr*(1.+0.5*deltaT/tr*pT/pmag);
        if(mass<0.) continue;
        mass = sqrt(mass);
        filled = true;
        _hist_mass->fill(mass);
        if(mass>90. ) {
          _count_90 ->fill(0.5);
          if(mass>110.) {
            _count_110->fill(0.5);
            if(mass>120.) {
              _count_120->fill(0.5);
              if(mass>130.) {
                _count_130->fill(0.5);
              }
            }
          }
        }
      }
      if(!filled) vetoEvent;
      // number passing beta cut
      _count_beta->fill(0.5);
    }

    //@}

    void finalize() {
      double fact = crossSection()/sumOfWeights()*37;
      MSG_WARNING("testing " << crossSection() << " " << sumOfWeights() << " " << fact);
      scale(_hist_beta,fact);
      scale(_hist_time,fact);
      scale(_hist_mass,fact);
      scale( _count_trigger, fact);
      scale( _count_event, fact);
      scale( _count_quality, fact);
      scale( _count_beta, fact);
      scale( _count_90, fact);
      scale( _count_110, fact);
      scale( _count_120, fact);
      scale( _count_130, fact);
    }

  private:

    /// @name Histograms
    //@{
    Histo1DPtr _hist_beta;
    Histo1DPtr _hist_time;
    Histo1DPtr _hist_mass;
    Histo1DPtr _count_trigger;
    Histo1DPtr _count_event;
    Histo1DPtr _count_quality;
    Histo1DPtr _count_beta;
    Histo1DPtr _count_90;
    Histo1DPtr _count_110;
    Histo1DPtr _count_120;
    Histo1DPtr _count_130;
    //@}

  };


  RIVET_DECLARE_ALIASED_PLUGIN(ATLAS_2011_S9108483, ATLAS_2011_I914491);

}
