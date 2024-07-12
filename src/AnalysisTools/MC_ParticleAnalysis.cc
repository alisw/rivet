// -*- C++ -*-
#include "Rivet/Analyses/MC_ParticleAnalysis.hh"

namespace Rivet {




  MC_ParticleAnalysis::MC_ParticleAnalysis(const string& name,
                                           size_t nparticles,
                                           const string& particle_name)
    : Analysis(name),
      _nparts(nparticles), _pname(particle_name),
      _h_pt(nparticles),
      _h_eta(nparticles), _h_eta_plus(nparticles), _h_eta_minus(nparticles),
      _h_rap(nparticles), _h_rap_plus(nparticles), _h_rap_minus(nparticles),
      tmpeta(nparticles), tmprap(nparticles)
  {

  }



  // Book histograms
  void MC_ParticleAnalysis::init() {

    for (size_t i = 0; i < _nparts; ++i) {
      book(tmpeta[i], _pname + "_eta_pmratio_" + to_str(i+1));
      book(tmprap[i], _pname + "_y_pmratio_" + to_str(i+1));

      const string ptname = _pname + "_pt_" + to_str(i+1);
      const double ptmax = 1.0/(double(i)+2.0) * (sqrtS()>0.?sqrtS():14000.)/GeV/2.0;
      const int nbins_pt = 100/(i+1);
      book(_h_pt[i] ,ptname, logspace(nbins_pt, 1.0, ptmax));

      const string etaname = _pname + "_eta_" + to_str(i+1);
      book(_h_eta[i] ,etaname, i > 1 ? 25 : 50, -5.0, 5.0);
      book(_h_eta_plus[i], "_" + etaname + "_plus", i > 1 ? 15 : 25, 0, 5);
      book(_h_eta_minus[i], "_" + etaname + "_minus", i > 1 ? 15 : 25, 0, 5);

      const string rapname = _pname + "_y_" + to_str(i+1);
      book(_h_rap[i] ,rapname, i > 1 ? 25 : 50, -5.0, 5.0);
      book(_h_rap_plus[i], "_" + rapname + "_plus", i > 1 ? 15 : 25, 0, 5);
      book(_h_rap_minus[i], "_" + rapname + "_minus", i > 1 ? 15 : 25, 0, 5);

      for (size_t j = i+1; j < min(size_t(3), _nparts); ++j) {
        const pair<size_t, size_t> ij = std::make_pair(i, j);

        string detaname = _pname + "s_deta_" + to_str(i+1) + to_str(j+1);
        Histo1DPtr tmpeta;
        book(tmpeta, detaname, 25, -5.0, 5.0);
        _h_deta.insert(make_pair(ij, tmpeta));

        string dphiname = _pname + "s_dphi_" + to_str(i+1) + to_str(j+1);
        Histo1DPtr tmpphi;
        book(tmpphi, dphiname, 25, 0.0, M_PI);
        _h_dphi.insert(make_pair(ij, tmpphi));

        string dRname = _pname + "s_dR_" + to_str(i+1) + to_str(j+1);
        Histo1DPtr tmpR;
        book(tmpR, dRname, 25, 0.0, 5.0);
        _h_dR.insert(make_pair(ij, tmpR));
      }
    }

    book(_h_multi_exclusive ,_pname + "_multi_exclusive", _nparts+3, -0.5, _nparts+3-0.5);
    book(_h_multi_inclusive ,_pname + "_multi_inclusive", _nparts+3, -0.5, _nparts+3-0.5);
    book(_h_multi_ratio, _pname + "_multi_ratio");

    book(_h_multi_exclusive_prompt ,_pname + "_multi_exclusive_prompt", _nparts+3, -0.5, _nparts+3-0.5);
    book(_h_multi_inclusive_prompt ,_pname + "_multi_inclusive_prompt", _nparts+3, -0.5, _nparts+3-0.5);
    book(_h_multi_ratio_prompt, _pname + "_multi_ratio_prompt");
  }


  // Do the analysis
  void MC_ParticleAnalysis::_analyze(const Event& event, const Particles& particles) {
    Particles promptparticles;
    for (const Particle& p : particles)
      if (p.isPrompt()) promptparticles += p;

    for (size_t i = 0; i < _nparts; ++i) {
      if (particles.size() < i+1) continue;
      _h_pt[i]->fill(particles[i].pt()/GeV);

      // Eta
      const double eta_i = particles[i].eta();
      _h_eta[i]->fill(eta_i);
      (eta_i > 0.0 ? _h_eta_plus : _h_eta_minus)[i]->fill(fabs(eta_i));

      // Rapidity
      const double rap_i = particles[i].rapidity();
      _h_rap[i]->fill(rap_i);
      (rap_i > 0.0 ? _h_rap_plus : _h_rap_minus)[i]->fill(fabs(rap_i));

      // Inter-particle properties
      for (size_t j = i+1; j < min(size_t(3),_nparts); ++j) {
        if (particles.size() < j+1) continue;
        std::pair<size_t, size_t> ij = std::make_pair(i, j);
        double deta = particles[i].eta() - particles[j].eta();
        double dphi = deltaPhi(particles[i].momentum(), particles[j].momentum());
        double dR = deltaR(particles[i].momentum(), particles[j].momentum());
        _h_deta[ij]->fill(deta);
        _h_dphi[ij]->fill(dphi);
        _h_dR[ij]->fill(dR);
      }
    }

    // Multiplicities
    _h_multi_exclusive->fill(particles.size());
    _h_multi_exclusive_prompt->fill(promptparticles.size());
    for (size_t i = 0; i < _nparts+2; ++i) {
      if (particles.size() >= i) _h_multi_inclusive->fill(i);
      if (promptparticles.size() >= i) _h_multi_inclusive_prompt->fill(i);
    }

  }


  // Finalize
  void MC_ParticleAnalysis::finalize() {
    const double scaling = crossSection()/sumOfWeights();
    for (size_t i = 0; i < _nparts; ++i) {
      scale(_h_pt[i], scaling);
      scale(_h_eta[i], scaling);
      scale(_h_rap[i], scaling);

      // Create eta/rapidity ratio plots
      divide(_h_eta_plus[i], _h_eta_minus[i], tmpeta[i]);
      divide(_h_rap_plus[i], _h_rap_minus[i], tmprap[i]);
    }

    // Scale the d{eta,phi,R} histograms
    typedef map<pair<size_t, size_t>, Histo1DPtr> HistMap;
    for (HistMap::value_type& it : _h_deta) scale(it.second, scaling);
    for (HistMap::value_type& it : _h_dphi) scale(it.second, scaling);
    for (HistMap::value_type& it : _h_dR) scale(it.second, scaling);

    // Fill inclusive multi ratios
    for (size_t i = 0; i < _h_multi_inclusive->numBins()-1; ++i) {
      _h_multi_ratio->addPoint(i+1, 0, 0.5, 0);
      if (_h_multi_inclusive->bin(i).sumW() > 0.0) {
        const double ratio = _h_multi_inclusive->bin(i+1).sumW() / _h_multi_inclusive->bin(i).sumW();
        const double relerr_i = _h_multi_inclusive->bin(i).relErr();
        const double relerr_j = _h_multi_inclusive->bin(i+1).relErr();
        const double err = ratio * (relerr_i + relerr_j);
        _h_multi_ratio->point(i).setY(ratio, err);
      }
    }
    for (size_t i = 0; i < _h_multi_inclusive_prompt->numBins()-1; ++i) {
      _h_multi_ratio_prompt->addPoint(i+1, 0, 0.5, 0);
      if (_h_multi_inclusive_prompt->bin(i).sumW() > 0.0) {
        const double ratio = _h_multi_inclusive_prompt->bin(i+1).sumW() / _h_multi_inclusive_prompt->bin(i).sumW();
        const double relerr_i = _h_multi_inclusive_prompt->bin(i).relErr();
        const double relerr_j = _h_multi_inclusive_prompt->bin(i+1).relErr();
        const double err = ratio * (relerr_i + relerr_j);
        _h_multi_ratio_prompt->point(i).setY(ratio, err);
      }
    }

    scale(_h_multi_exclusive, scaling);
    scale(_h_multi_exclusive_prompt, scaling);
    scale(_h_multi_inclusive, scaling);
    scale(_h_multi_inclusive_prompt, scaling);
  }


}
