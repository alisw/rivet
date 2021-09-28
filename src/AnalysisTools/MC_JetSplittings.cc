// -*- C++ -*-
#include "Rivet/Analyses/MC_JetSplittings.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {




  MC_JetSplittings::MC_JetSplittings(const string& name,
                                     size_t njet,
                                     const string& jetpro_name)
    : Analysis(name), m_njet(njet), m_jetpro_name(jetpro_name),
      _h_log10_d(njet), _h_log10_R(njet+1)
  {

  }


  // Book histograms
  void MC_JetSplittings::init() {
    const double sqrts = sqrtS() ? sqrtS() : 14000.*GeV;

    for (size_t i = 0; i < m_njet; ++i) {
      string dname = "log10_d_" + to_str(i) + to_str(i+1);
      book(_h_log10_d[i] ,dname, 100, 0.2, log10(0.5*sqrts/GeV));
      string Rname = "log10_R_" + to_str(i);
      book(_h_log10_R[i], Rname, 50, 0.2, log10(0.5*sqrts/GeV));
    }
    string Rname = "log10_R_" + to_str(m_njet);
    book(_h_log10_R[m_njet], Rname, 50, 0.2, log10(0.5*sqrts/GeV));
  }



  // Do the analysis
  void MC_JetSplittings::analyze(const Event & e) {
    const double weight = 1.0;

    const FastJets& jetpro = apply<FastJets>(e, m_jetpro_name);
    const auto seq = jetpro.clusterSeq();
    if (!seq) vetoEvent; //< the cseq is the whole point in this sort of analysis!!

    // Jet resolutions and integrated jet rates
    double previous_dij = 10.0;
    for (size_t i = 0; i < min(m_njet,(size_t)seq->n_particles()); ++i) {
      const double d_ij2 = seq->exclusive_dmerge_max(i);
      if (d_ij2 <= 0) continue; ///< @todo Is < 0 possible? Feels like no; I should check ;-)
      // Jet resolution i -> j
      const double d_ij = log10(sqrt(d_ij2));

      // Fill differential jet resolution
      _h_log10_d[i]->fill(d_ij, weight);

      // Fill integrated jet resolution
      for (size_t ibin = 0; ibin < _h_log10_R[i]->numPoints(); ++ibin) {
        Point2D& dp = _h_log10_R[i]->point(ibin);
        if (dp.x() > d_ij && dp.x() < previous_dij) {
          dp.setY(dp.y() + weight);
        }
      }
      previous_dij = d_ij;
    }
    // One remaining integrated jet resolution
    for (size_t ibin = 0; ibin<_h_log10_R[m_njet]->numPoints(); ++ibin) {
      Point2D & dp = _h_log10_R[m_njet]->point(ibin);
      if (dp.x() < previous_dij) {
        dp.setY(dp.y() + weight);
      }
    }

  }


  // Finalize
  void MC_JetSplittings::finalize() {
    const double xsec_unitw = crossSection()/picobarn/sumOfWeights();
    for (size_t i = 0; i < m_njet; ++i) {
      scale(_h_log10_d[i], xsec_unitw);
      for (size_t ibin = 0; ibin<_h_log10_R[i]->numPoints(); ++ibin) {
        Point2D& dp = _h_log10_R[i]->point(ibin);
        dp.setY(dp.y()*xsec_unitw);
      }
    }

    for (size_t ibin = 0; ibin < _h_log10_R[m_njet]->numPoints(); ++ibin) {
      Point2D& dp =_h_log10_R[m_njet]->point(ibin);
      dp.setY(dp.y()*xsec_unitw);
    }
  }


}
