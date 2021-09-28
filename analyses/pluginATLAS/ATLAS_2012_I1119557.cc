// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  namespace { // unnamed namespace
    // Forward declarations of calculator functions: implementations at bottom of file
    double getWidth(const Jet& jet);
    /// @todo Re-enable eccentricity calculation
    // double getEcc(const Jet& jet);
    double getPFlow(const Jet& jet);
    double getAngularity(const Jet& jet);
  }


  class ATLAS_2012_I1119557 : public Analysis {
  public:

    ATLAS_2012_I1119557()
      : Analysis("ATLAS_2012_I1119557")
    {    }


    void init() {
      const FinalState fs;
      declare(fs, "FinalState");

      FastJets fj06(fs, FastJets::ANTIKT, 0.6);
      declare(fj06, "AntiKT06");
      FastJets fj10(fs, FastJets::ANTIKT, 1.0);
      declare(fj10, "AntiKT10");

      for (size_t alg = 0; alg < 2; ++alg) {
        book(_hs_mass[alg]  ,1, alg+1, 1);
        book(_hs_width[alg] ,2, alg+1, 1);
        /// @todo Commented eccentricity out for now: reinstate
        //book( _hs_eccentricity[alg] ,3, alg+1, 1);
      }
      book(_h_planarFlow ,4, 2, 1);
      book(_h_angularity ,5, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Jets jetAr[2];
      jetAr[0] = apply<FastJets>(event, "AntiKT06").jetsByPt(Cuts::pT > 300*GeV && Cuts::abseta < 2.0);
      jetAr[1] = apply<FastJets>(event, "AntiKT10").jetsByPt(Cuts::pT > 300*GeV && Cuts::abseta < 2.0);

      for (size_t alg = 0; alg < 2; ++alg) {
        // Require at least one jet
        if (jetAr[alg].size() < 1) continue;

        // The leading jet
        const Jet& jet = jetAr[alg][0];
        const double m   = jet.mass();
        const double eta = jet.eta();

        _hs_mass[alg]->fill(m/GeV);
        _hs_width[alg]->fill(getWidth(jet));
        /// @todo Commented eccentricity out for now: reinstate
        // if (fabs(eta) < 0.7 && m > 100*GeV) _hs_eccentricity[alg]->fill(getEcc(jet));

        if (fabs(eta) < 0.7) {
          if (alg == 0 && inRange(m/GeV, 100., 130.)) _h_angularity->fill(getAngularity(jet));
          if (alg == 1 && inRange(m/GeV, 130., 210.)) _h_planarFlow->fill(getPFlow(jet));
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (size_t alg = 0; alg < 2; ++alg) {
        normalize(_hs_mass[alg]);
        normalize(_hs_width[alg]);
        /// @todo Commented eccentricity out for now: reinstate
        // normalize(_hs_eccentricity[alg]);
      }
      normalize(_h_planarFlow);
      normalize(_h_angularity);
    }


  private:

    Histo1DPtr _hs_mass[2];
    Histo1DPtr _hs_width[2];
    /// @todo Commented eccentricity out for now: reinstate
    // Histo1DPtr _hs_eccentricity[2];
    Histo1DPtr _h_planarFlow;
    Histo1DPtr _h_angularity;
  };


  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1119557);



  namespace { // unnamed namespace

    // Adapted code from Lily
    /// @todo Convert to use the Rivet rotation matrix code (should be simpler)
    FourMomentum RotateAxes(const Rivet::FourMomentum& p, double M[3][3]){
      double px_rot = M[0][0]*(p.px()) + M[0][1]*(p.py()) + M[0][2]*(p.pz());
      double py_rot = M[1][0]*(p.px()) + M[1][1]*(p.py()) + M[1][2]*(p.pz());
      double pz_rot = M[2][0]*(p.px()) + M[2][1]*(p.py()) + M[2][2]*(p.pz());
      return FourMomentum(p.E(), px_rot, py_rot, pz_rot);
    }

    // Untouched code from Lily
    /// @todo Convert to use the Rivet rotation matrix code (should be simpler)
    void CalcRotationMatrix(double nvec[3],double rot_mat[3][3]){
      // clear momentum tensor
      for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
          rot_mat[i][j] = 0.;
        }
      }
      double mag3 = sqrt(nvec[0]*nvec[0] + nvec[1]*nvec[1]+ nvec[2]*nvec[2]);
      double mag2 = sqrt(nvec[0]*nvec[0] + nvec[1]*nvec[1]);
      /// @todo cout is not a valid response to a numerical error! Is the error condition reported?!? Assert added by AB for Rivet 1.8.2
      assert(mag3 > 0);
      if (mag3 <= 0) {
        cout << "rotation axis is null" << '\n';
        return;
      }

      double ctheta0 = nvec[2]/mag3;
      double stheta0 = mag2/mag3;
      double cphi0 = (mag2 > 0.) ? nvec[0]/mag2 : 0;
      double sphi0 = (mag2 > 0.) ? nvec[1]/mag2 : 0;

      rot_mat[0][0] = -ctheta0*cphi0;
      rot_mat[0][1] = -ctheta0*sphi0;
      rot_mat[0][2] = stheta0;
      rot_mat[1][0] = sphi0;
      rot_mat[1][1] = -1.*cphi0;
      rot_mat[1][2] = 0.;
      rot_mat[2][0] = stheta0*cphi0;
      rot_mat[2][1] = stheta0*sphi0;
      rot_mat[2][2] = ctheta0;
    }


    /// Jet width calculation
    double getWidth(const Jet& jet) {
      const double phi_jet = jet.phi();
      const double eta_jet = jet.eta();
      double width(0), pTsum(0);
      for (const Particle& p : jet.particles()) {
        double pT = p.pT();
        double eta = p.eta();
        double phi = p.phi();
        width += sqrt(pow(phi_jet - phi,2) + pow(eta_jet - eta, 2)) * pT;
        pTsum += pT;
      }
      const double rtn = (pTsum != 0.0) ? width/pTsum : -1;
      return rtn;
    }


    /// Eccentricity calculation, copied and adapted from Lily's code
    /// @todo Re-enable eccentricity calculation
    // double getEcc(const Jet& jet) {
    //   vector<double> phis;
    //   vector<double> etas;
    //   vector<double> energies;

    //   double etaSum(0), phiSum(0), eTot(0);
    //   for (const Particle& p : jet.particles()) {
    //     const double E = p.E();
    //     const double eta = p.eta();

    //     energies.push_back(E);
    //     etas.push_back(jet.eta() - eta);

    //     eTot   += E;
    //     etaSum += eta * E;

    //     /// @todo Replace with the Rivet deltaPhi function (or the mapAngleTo* functions)
    //     double dPhi = jet.phi() - p.phi();
    //     //if DPhi does not lie within 0 < DPhi < PI take 2*PI off DPhi
    //     //this covers cases where DPhi is greater than PI
    //     if( fabs( dPhi - TWOPI ) < fabs(dPhi) ) dPhi -= TWOPI;
    //     //if DPhi does not lie within -PI < DPhi < 0 add 2*PI to DPhi
    //     //this covers cases where DPhi is less than -PI
    //     else if( fabs(dPhi + TWOPI) < fabs(dPhi) ) dPhi += TWOPI;
    //     phis.push_back(dPhi);

    //     phiSum += dPhi * E;
    //   }

    //   //these are the "pull" away from the jet axis
    //   etaSum = etaSum/eTot;
    //   phiSum = phiSum/eTot;

    //   // now for every cluster we alter its position by moving it:
    //   // away from the new axis if it is in the direction of -ve pull
    //   // closer to the new axis if it is in the direction of +ve pull
    //   // the effect of this will be that the new energy weighted center will be on the old jet axis.
    //   double little_x(0), little_y(0);
    //   for (size_t k = 0; k < jet.particles().size(); ++k) {
    //     little_x+= etas[k]-etaSum;
    //     little_y+= phis[k]-phiSum;
    //     etas[k] = etas[k]-etaSum;
    //     phis[k] = phis[k]-phiSum;
    //   }

    //   double x1(0), x2(0);
    //   for (size_t i = 0; i < jet.particles().size(); ++i) {
    //     x1 += 2. * energies[i]* etas[i] * phis[i]; // this is 2*X*Y
    //     x2 += energies[i]*(phis[i] * phis[i] - etas[i] * etas[i] ); // this is X^2 - Y^2
    //   }

    //   // Variance calculations
    //   double theta = .5*atan2(x1, x2);
    //   double sinTheta =sin(theta);
    //   double cosTheta = cos(theta);
    //   double theta2 = theta + 0.5*PI;
    //   double sinThetaPrime = sin(theta2);
    //   double cosThetaPrime = cos(theta2);

    //   double varX(0), varY(0);
    //   for (size_t i = 0; i < jet.particles().size(); i++) {
    //     const double x = sinTheta*etas[i] + cosTheta*phis[i];
    //     const double y = sinThetaPrime*etas[i] + cosThetaPrime*phis[i];
    //     varX += energies[i]* sqr(x);
    //     varY += energies[i]* sqr(y);
    //   }
    //   const double varianceMax = max(varX, varY);
    //   const double varianceMin = min(varX, varY);
    //   const double ecc = (varianceMax != 0.0) ? 1 - varianceMin/varianceMax : -1;
    //   return ecc;
    // }


    /// Planar flow calculation, copied and adapted from Lily's code
    double getPFlow(const Jet& jet) {
      const double phi0 = jet.phi();
      const double eta0 = jet.eta();

      double nref[3]; ///< @todo 3-vector to rotate x to? Use Rivet vector classes
      nref[0] = cos(phi0)/cosh(eta0);
      nref[1] = sin(phi0)/cosh(eta0);
      nref[2] = tanh(eta0);

      // Rotation matrix to align with nref
      double rotationMatrix[3][3];
      CalcRotationMatrix(nref, rotationMatrix);

      double iw00(0.), iw01(0.), iw11(0.), iw10(0.);
      for (const Particle& p : jet.particles()) {
        double a = 1./(p.E()*jet.mass());
        FourMomentum rotclus = RotateAxes(p.momentum(), rotationMatrix);
        iw00 += a*pow(rotclus.px(), 2);
        iw01 += a*rotclus.px()*rotclus.py();
        iw10 += a*rotclus.py()*rotclus.px();
        iw11 += a*pow(rotclus.py(), 2);
      }

      const double det = iw00*iw11 - iw01*iw10;
      const double trace = iw00 + iw11;

      const double pf = (trace != 0.0) ? (4.0*det)/sqr(trace) : -1;
      return pf;
    }


    /// Angularity calculation, copied and adapted from Lily's code
    double getAngularity(const Jet& jet) {
      double sum_a = 0.;
      // a can take any value < 2 (e.g. 1,0,-0.5 etc) for infrared safety
      const double a = -2.;
      for (const Particle& p : jet.particles()) {
        double e_i       = p.E();
        double theta_i   = jet.momentum().angle(p.momentum());
        double e_theta_i = e_i * pow(sin(theta_i), a) * pow(1-cos(theta_i), 1-a);
        sum_a += e_theta_i;
      }
      const double rtn = (jet.mass() != 0.0 && !std::isnan(sum_a)) ? sum_a/jet.mass() : -1;
      return rtn;
    }

  }


}
