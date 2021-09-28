// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Sphericity.hh"


namespace Rivet {

  /// @brief Multijet event shapes at 13 TeV
  class ATLAS_2020_I1808726 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2020_I1808726);

    double xs1 = 0.0;
    double xs2 = 0.0;
    double xs3 = 0.0;

    /// Initialization, called once before running
    void init() {

      //Jet collection (excluding muons and neutrinos)
      const FinalState fs(Cuts::abseta < 4.5);

      FastJets jets(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jets, "Jets");

      // Book histograms
      //Jet multiplicity
      book(_h["njet_h1"], 73, 1, 1);
      book(_h["njet_h2"], 74, 1, 1);
      book(_h["njet_h3"], 75, 1, 1);

      //Transverse thrust
      book(_h["transThrust_j3_h1"], 1, 1, 1);
      book(_h["transThrust_j3_h2"], 5, 1, 1);
      book(_h["transThrust_j3_h3"], 9, 1, 1);
      book(_h["transThrust_j4_h1"], 2, 1, 1);
      book(_h["transThrust_j4_h2"], 6, 1, 1);
      book(_h["transThrust_j4_h3"], 10, 1, 1);
      book(_h["transThrust_j5_h1"], 3, 1, 1);
      book(_h["transThrust_j5_h2"], 7, 1, 1);
      book(_h["transThrust_j5_h3"], 11, 1, 1);
      book(_h["transThrust_j6_h1"], 4, 1, 1);
      book(_h["transThrust_j6_h2"], 8, 1, 1);
      book(_h["transThrust_j6_h3"], 12, 1, 1);

      //Thrust minor
      book(_h["transMinor_j3_h1"], 13, 1, 1);
      book(_h["transMinor_j3_h2"], 17, 1, 1);
      book(_h["transMinor_j3_h3"], 21, 1, 1);
      book(_h["transMinor_j4_h1"], 14, 1, 1);
      book(_h["transMinor_j4_h2"], 18, 1, 1);
      book(_h["transMinor_j4_h3"], 22, 1, 1);
      book(_h["transMinor_j5_h1"], 15, 1, 1);
      book(_h["transMinor_j5_h2"], 19, 1, 1);
      book(_h["transMinor_j5_h3"], 23, 1, 1);
      book(_h["transMinor_j6_h1"], 16, 1, 1);
      book(_h["transMinor_j6_h2"], 20, 1, 1);
      book(_h["transMinor_j6_h3"], 24, 1, 1);

      //Transverse sphericity
      book(_h["transSphericity_j3_h1"], 25, 1, 1);
      book(_h["transSphericity_j3_h2"], 29, 1, 1);
      book(_h["transSphericity_j3_h3"], 33, 1, 1);
      book(_h["transSphericity_j4_h1"], 26, 1, 1);
      book(_h["transSphericity_j4_h2"], 30, 1, 1);
      book(_h["transSphericity_j4_h3"], 34, 1, 1);
      book(_h["transSphericity_j5_h1"], 27, 1, 1);
      book(_h["transSphericity_j5_h2"], 31, 1, 1);
      book(_h["transSphericity_j5_h3"], 35, 1, 1);
      book(_h["transSphericity_j6_h1"], 28, 1, 1);
      book(_h["transSphericity_j6_h2"], 32, 1, 1);
      book(_h["transSphericity_j6_h3"], 36, 1, 1);

      //Aplanarity
      book(_h["aplanarity_j3_h1"], 37, 1, 1);
      book(_h["aplanarity_j3_h2"], 41, 1, 1);
      book(_h["aplanarity_j3_h3"], 45, 1, 1);
      book(_h["aplanarity_j4_h1"], 38, 1, 1);
      book(_h["aplanarity_j4_h2"], 42, 1, 1);
      book(_h["aplanarity_j4_h3"], 46, 1, 1);
      book(_h["aplanarity_j5_h1"], 39, 1, 1);
      book(_h["aplanarity_j5_h2"], 43, 1, 1);
      book(_h["aplanarity_j5_h3"], 47, 1, 1);
      book(_h["aplanarity_j6_h1"], 40, 1, 1);
      book(_h["aplanarity_j6_h2"], 44, 1, 1);
      book(_h["aplanarity_j6_h3"], 48, 1, 1);

      //C
      book(_h["C_j3_h1"], 49, 1, 1);
      book(_h["C_j3_h2"], 53, 1, 1);
      book(_h["C_j3_h3"], 57, 1, 1);
      book(_h["C_j4_h1"], 50, 1, 1);
      book(_h["C_j4_h2"], 54, 1, 1);
      book(_h["C_j4_h3"], 58, 1, 1);
      book(_h["C_j5_h1"], 51, 1, 1);
      book(_h["C_j5_h2"], 55, 1, 1);
      book(_h["C_j5_h3"], 59, 1, 1);
      book(_h["C_j6_h1"], 52, 1, 1);
      book(_h["C_j6_h2"], 56, 1, 1);
      book(_h["C_j6_h3"], 60, 1, 1);

      //D
      book(_h["D_j3_h1"], 61, 1, 1);
      book(_h["D_j3_h2"], 65, 1, 1);
      book(_h["D_j3_h3"], 69, 1, 1);
      book(_h["D_j4_h1"], 62, 1, 1);
      book(_h["D_j4_h2"], 66, 1, 1);
      book(_h["D_j4_h3"], 70, 1, 1);
      book(_h["D_j5_h1"], 63, 1, 1);
      book(_h["D_j5_h2"], 67, 1, 1);
      book(_h["D_j5_h3"], 71, 1, 1);
      book(_h["D_j6_h1"], 64, 1, 1);
      book(_h["D_j6_h2"], 68, 1, 1);
      book(_h["D_j6_h3"], 72, 1, 1);
    }
    
    void analyze(const Event& event) {
      
      const Jets& jets = apply<FastJets>(event, "Jets").jetsByPt(7.0*GeV);
      
      //Select jets passing kinematic cuts
      std::vector<const Jet*> goodJets; goodJets.clear();
      std::vector<Vector3> momenta2; momenta2.clear();
      std::vector<Vector3> momenta3; momenta3.clear();

      //foreach (const Jet& j, jets) {
      for (const Jet& j : jets) {
        if (j.abseta() < 2.4 && j.pt() > 100.0*GeV){
	         goodJets.push_back(&j);

             Vector3 jet2 = j.p3();
             jet2.setZ(0.0);
              momenta2.push_back(jet2);

	          Vector3 jet3 = j.p3();
	          momenta3.push_back(jet3);
	       }
        }
      
      //Dijet event selection
      if (goodJets.size() < 2) vetoEvent;
      double ht2 = goodJets[0]->pt()+goodJets[1]->pt();
      if (ht2 <= 1000.0*GeV) vetoEvent;

      //Jet multiplicity
      if (ht2 > 1000.0*GeV && ht2 < 1500.0*GeV) _h["njet_h1"]->fill(goodJets.size());
      if (ht2 > 1500.0*GeV && ht2 < 2000.0*GeV) _h["njet_h2"]->fill(goodJets.size());
      if (ht2 > 2000.0*GeV) _h["njet_h3"]->fill(goodJets.size());
      
      //Thrust calculation
      Thrust thrust;
      thrust.calc(momenta2);
      const double transThrust  = 1.0 - thrust.thrust();
      const double transMinor = thrust.thrustMajor();
      
      //Linearized sphericity calculation (2D)
      double a11 = 0.0; double a22 = 0.0;
      double a12 = 0.0;
      double modSum2 = 0.0;

      for (size_t k = 0; k < momenta2.size(); ++k) {
        modSum2 += momenta2[k].mod();
	      a11 += momenta2[k].x()*momenta2[k].x()/momenta2[k].mod();
        a22 += momenta2[k].y()*momenta2[k].y()/momenta2[k].mod();
        a12 += momenta2[k].x()*momenta2[k].y()/momenta2[k].mod();
      }

      double trc2 = (a11+a22)/modSum2;
      double det2 = (a11*a22-a12*a12)/pow(modSum2,2);

      double eigen21 = (trc2+sqrt(pow(trc2,2)-4*det2))/2;
      double eigen22 = (trc2-sqrt(pow(trc2,2)-4*det2))/2;
      double transSphericity = 2*eigen22/(eigen21+eigen22);

      //Linearized sphericity calculation (3D)
      double b11 = 0.0; double b12 = 0.0; double b13 = 0.0;
      double b22 = 0.0; double b23 = 0.0;
      double b33 = 0.0;
      double modSum3 = 0.0;

      for (size_t k = 0; k < momenta3.size(); ++k){
        modSum3 += momenta3[k].mod();
        b11 += momenta3[k].x()*momenta3[k].x()/momenta3[k].mod();
        b22 += momenta3[k].y()*momenta3[k].y()/momenta3[k].mod();
        b33 += momenta3[k].z()*momenta3[k].z()/momenta3[k].mod();
        b12 += momenta3[k].x()*momenta3[k].y()/momenta3[k].mod();
        b13 += momenta3[k].x()*momenta3[k].z()/momenta3[k].mod();
        b23 += momenta3[k].y()*momenta3[k].z()/momenta3[k].mod();
      }

      Matrix3 sph3;
      sph3.set(0,0, b11/modSum3); sph3.set(0,1, b12/modSum3); sph3.set(0,2, b13/modSum3);
      sph3.set(1,0, b12/modSum3); sph3.set(1,1, b22/modSum3); sph3.set(1,2, b23/modSum3);
      sph3.set(2,0, b13/modSum3); sph3.set(2,1, b23/modSum3); sph3.set(2,2, b33/modSum3);

      double q = sph3.trace()/3.;
      double p1 = sph3.get(0,1)*sph3.get(0,1) + sph3.get(0,2)*sph3.get(0,2) + sph3.get(1,2)*sph3.get(1,2);
      double p2 = (sph3.get(0,0)-q)*(sph3.get(0,0)-q) + (sph3.get(1,1)-q)*(sph3.get(1,1)-q) + (sph3.get(2,2)-q)*(sph3.get(2,2)-q) + 2*p1;
      double p = sqrt(p2/6.);
      
      Matrix3 I3 = Matrix3::mkIdentity();
      double r = ( 1./p * (sph3 - q*I3)).det()/2.;
      
      double phi(0);
      if (r <= -1) phi = M_PI / 3.;
      else if (r >= 1) phi = 0;
      else phi = acos(r) / 3.;

      double eigen31 = q + 2 * p * cos(phi);
      double eigen33 = q + 2 * p * cos(phi + (2*M_PI/3.));
      double eigen32 = 3 * q - eigen31 - eigen33;

      double aplanarity = (3./2)*eigen33;
      double C = 3*(eigen31*eigen32 + eigen31*eigen33 + eigen32*eigen33);
      double D = 27*eigen31*eigen32*eigen33;

      //Fill event-shape histograms
      if (ht2 > 1000.0*GeV && ht2 < 1500.0*GeV){

      	if (goodJets.size() == 3){ 
	       _h["transThrust_j3_h1"]->fill(transThrust); _h["transMinor_j3_h1"]->fill(transMinor);
	       _h["transSphericity_j3_h1"]->fill(transSphericity); _h["aplanarity_j3_h1"]->fill(aplanarity);
	       _h["C_j3_h1"]->fill(C); _h["D_j3_h1"]->fill(D);
	      }

	   if (goodJets.size() == 4){
          _h["transThrust_j4_h1"]->fill(transThrust); _h["transMinor_j4_h1"]->fill(transMinor);
          _h["transSphericity_j4_h1"]->fill(transSphericity); _h["aplanarity_j4_h1"]->fill(aplanarity);
          _h["C_j4_h1"]->fill(C); _h["D_j4_h1"]->fill(D);
	     }

	   if (goodJets.size() == 5){
          _h["transThrust_j5_h1"]->fill(transThrust); _h["transMinor_j5_h1"]->fill(transMinor);
          _h["transSphericity_j5_h1"]->fill(transSphericity); _h["aplanarity_j5_h1"]->fill(aplanarity);
          _h["C_j5_h1"]->fill(C); _h["D_j5_h1"]->fill(D);
     	}

	   if (goodJets.size() >= 6){
          _h["transThrust_j6_h1"]->fill(transThrust); _h["transMinor_j6_h1"]->fill(transMinor);
          _h["transSphericity_j6_h1"]->fill(transSphericity); _h["aplanarity_j6_h1"]->fill(aplanarity);
          _h["C_j6_h1"]->fill(C); _h["D_j6_h1"]->fill(D);
	    }
      }


      if (ht2 > 1500.0*GeV && ht2 < 2000.0*GeV){

        if (goodJets.size() == 3){
          _h["transThrust_j3_h2"]->fill(transThrust); _h["transMinor_j3_h2"]->fill(transMinor);
          _h["transSphericity_j3_h2"]->fill(transSphericity); _h["aplanarity_j3_h2"]->fill(aplanarity);
          _h["C_j3_h2"]->fill(C); _h["D_j3_h2"]->fill(D);
        }

        if (goodJets.size() == 4){
          _h["transThrust_j4_h2"]->fill(transThrust); _h["transMinor_j4_h2"]->fill(transMinor);
          _h["transSphericity_j4_h2"]->fill(transSphericity); _h["aplanarity_j4_h2"]->fill(aplanarity);
          _h["C_j4_h2"]->fill(C); _h["D_j4_h2"]->fill(D);
        }

        if (goodJets.size() == 5){
          _h["transThrust_j5_h2"]->fill(transThrust); _h["transMinor_j5_h2"]->fill(transMinor);
          _h["transSphericity_j5_h2"]->fill(transSphericity); _h["aplanarity_j5_h2"]->fill(aplanarity);
          _h["C_j5_h2"]->fill(C); _h["D_j5_h2"]->fill(D);
        }

        if (goodJets.size() >= 6){
          _h["transThrust_j6_h2"]->fill(transThrust); _h["transMinor_j6_h2"]->fill(transMinor);
          _h["transSphericity_j6_h2"]->fill(transSphericity); _h["aplanarity_j6_h2"]->fill(aplanarity);
          _h["C_j6_h2"]->fill(C); _h["D_j6_h2"]->fill(D);
        }
      }

      if (ht2 > 2000.0*GeV){

        if (goodJets.size() == 3){
          _h["transThrust_j3_h3"]->fill(transThrust); _h["transMinor_j3_h3"]->fill(transMinor);
          _h["transSphericity_j3_h3"]->fill(transSphericity); _h["aplanarity_j3_h3"]->fill(aplanarity);
          _h["C_j3_h3"]->fill(C); _h["D_j3_h3"]->fill(D);
        }

        if (goodJets.size() == 4){
          _h["transThrust_j4_h3"]->fill(transThrust); _h["transMinor_j4_h3"]->fill(transMinor);
          _h["transSphericity_j4_h3"]->fill(transSphericity); _h["aplanarity_j4_h3"]->fill(aplanarity);
          _h["C_j4_h3"]->fill(C); _h["D_j4_h3"]->fill(D);
        }

        if (goodJets.size() == 5){
          _h["transThrust_j5_h3"]->fill(transThrust); _h["transMinor_j5_h3"]->fill(transMinor);
          _h["transSphericity_j5_h3"]->fill(transSphericity); _h["aplanarity_j5_h3"]->fill(aplanarity);
          _h["C_j5_h3"]->fill(C); _h["D_j5_h3"]->fill(D);
        }

        if (goodJets.size() >= 6){
          _h["transThrust_j6_h3"]->fill(transThrust); _h["transMinor_j6_h3"]->fill(transMinor);
          _h["transSphericity_j6_h3"]->fill(transSphericity); _h["aplanarity_j6_h3"]->fill(aplanarity);
          _h["C_j6_h3"]->fill(C); _h["D_j6_h3"]->fill(D);
        }
      }
    }


    void finalize() {
    
      const double xs1 = _h["njet_h1"]->sumW();
      const double xs2 = _h["njet_h2"]->sumW();
      const double xs3 = _h["njet_h3"]->sumW();
      
      for (auto& hist : _h) {
        if (hist.first.find("njet_") != string::npos)  scale(hist.second, crossSectionPerEvent()/picobarn);
        else if (hist.first.find("_h1") != string::npos) scale(hist.second, 1.0/xs1);
        else if (hist.first.find("_h2") != string::npos) scale(hist.second, 1.0/xs2);
        else if (hist.first.find("_h3") != string::npos) scale(hist.second, 1.0/xs3);
      }
    }

  private:
  
    //Jet multiplicity
    map<string,Histo1DPtr> _h;

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2020_I1808726);
}
