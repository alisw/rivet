// -*- C++ -*-

#include "../Core/zstr/zstr.hpp"
#include "HepMC3/GenCrossSection.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/ReaderAsciiHepMC2.h"
#include "HepMC3/ReaderFactory.h"
#include "Rivet/Tools/RivetHepMC.hh"
#include <cassert>

namespace Rivet {
  namespace HepMCUtils {


    ConstGenParticlePtr getParticlePtr(const RivetHepMC::GenParticle& gp) {
      return gp.shared_from_this();
    }

    std::vector<ConstGenParticlePtr> particles(ConstGenEventPtr ge) {
      assert(ge != nullptr);
      return ge->particles();
    }

    std::vector<ConstGenParticlePtr> particles(const GenEvent* ge) {
      assert(ge != nullptr);
      return ge->particles();
    }

    std::vector<ConstGenVertexPtr> vertices(ConstGenEventPtr ge) {
      assert(ge != nullptr);
      return ge->vertices();
    }

    std::vector<ConstGenVertexPtr> vertices(const GenEvent* ge) {
      assert(ge != nullptr);
      return ge->vertices();
    }

    std::vector<ConstGenParticlePtr> particles(ConstGenVertexPtr gv, const Relatives& relo) {
      return relo(gv);
    }

    std::vector<ConstGenParticlePtr> particles(ConstGenParticlePtr gp, const Relatives& relo) {
      return relo(gp);
    }

    int particles_size(ConstGenEventPtr ge) {
      return particles(ge).size();
    }

    int particles_size(const GenEvent* ge) {
      return particles(ge).size();
    }

    int uniqueId(ConstGenParticlePtr gp) {
      return (gp != nullptr) ? gp->id() : 0;
    }


    std::pair<ConstGenParticlePtr, ConstGenParticlePtr> beams(const GenEvent* ge) {
      assert(ge != nullptr);
      std::vector<ConstGenParticlePtr> beamlist = ge->beams();
      if (beamlist.size() < 2) {
        std::cerr << "CANNOT FIND ANY BEAMS!" << std::endl;
        return std::pair<ConstGenParticlePtr, ConstGenParticlePtr>();
      }
      return std::make_pair(beamlist[0], beamlist[1]);
    }


    bool readEvent(std::shared_ptr<HepMC_IO_type> io, std::shared_ptr<GenEvent> evt) {
      return io->read_event(*evt) && !io->failed();
      /// @todo Any problem due to these?! Factored failure return-lines are nicer if we can have them
      // if (!io->read_event(*evt)) return false;
      // if (io->failed()) return false;
      /// @todo Check that this is working when reading from a MEV-unit file... or should the reader auto-convert evt is GEV and io is MEV?
      evt->set_units(HepMC3::Units::GEV, HepMC3::Units::MM);
    }


    shared_ptr<HepMC_IO_type> makeReader(std::string filename, std::shared_ptr<std::istream>& istrp, std::string* errm) {
      shared_ptr<HepMC_IO_type> ret;

      #ifdef HAVE_LIBZ
      if (filename == "-")
        istrp = make_shared<zstr::istream>(std::cin);
      else
        istrp = make_shared<zstr::ifstream>(filename.c_str());
      std::istream& istr = *istrp;
      #else
      if (filename != "-") istrp = make_shared<std::ifstream>(filename.c_str());
      std::istream& istr = filename == "-" ? std::cin : *istrp;
      #endif

      ret = RivetHepMC::deduce_reader(istr);

      // Check that everything was ok.
      if (ret) {
        if (ret->failed()) {
          if (errm) *errm = "Problems reading from HepMC file. ";
          ret = shared_ptr<HepMC_IO_type>();
        }
        return ret;
      }
      if (!ret && filename == "-") {
        if (errm) *errm += "Problems reading HepMC from stdin. No header found. ";
        return shared_ptr<HepMC_IO_type>();
      }

      return ret;
    }


    void strip(GenEvent& ge, const set<long>& stripid) {
      //      std::cerr << "Stripping event " << ge.event_number() << std::endl;
      vector<HepMC3::GenParticlePtr> allparticles = ge.particles();
      for (auto& p : allparticles) {
        if (!p->production_vertex() || !p->end_vertex() || stripid.count(p->pid()) == 0 || p->production_vertex()->id() == 0)
          continue;
        // std::cout << "Removing particle " << p->id() << " (" << p->pid() << ")" << std::endl;
        HepMC3::GenVertexPtr vp = p->production_vertex();
        HepMC3::GenVertexPtr ve = p->end_vertex();
        if (!vp || !ve) continue;
        if (vp == ve) continue;
        // Check if the vertices would leave particles with the sam
        // production as decay vertex - we don't want that.
        if ((vp->particles_out().size() == 1 && vp->particles_out()[0] == p) ||
            (ve->particles_in().size() == 1 && ve->particles_in()[0] == p)) {
          bool loop = false;
          for (auto pi : vp->particles_in())
            for (auto po : ve->particles_out())
              if (pi == po) loop = true;
          if (loop) continue;
        }
        if (vp->particles_in().size() == 1 && (vp->particles_in()[0]->pid() > 21 && vp->particles_in()[0]->pid() < 30)) continue;

        vp->remove_particle_out(p);
        ve->remove_particle_in(p);

        if (ve->particles_in().empty()) {
          auto prem = ve->particles_out();
          for (auto po : prem) vp->add_particle_out(po);
          ge.remove_vertex(ve);
        } else if (vp->particles_out().empty()) {
          auto prem = vp->particles_in();
          for (auto pi : prem) ve->add_particle_in(pi);
          ge.remove_vertex(vp);
        }
        ge.remove_particle(p);
      }
    }


    pair<double, double> crossSection(const GenEvent& ge, size_t index) {
      if (!ge.cross_section()) {
        printf("Cross-section not set for GenEvent! Will return dummy value.\n");
        return make_pair(0.0, 0.0);
      }
      // Work-around since access functions are not const.
      HepMC3::GenCrossSection xs = *ge.cross_section();
      try {
        return make_pair(xs.xsec(index), xs.xsec_err(index));
      } catch (...) {
        // Probably a HepMC2 input file which only
        // has a single (nominal) cross-section
      }
      return make_pair(xs.xsec(0), xs.xsec_err(0));
    }


    vector<string> weightNames(const GenEvent& ge) {
      vector<string> ret;
      try {
        ret = ge.weight_names();
      } catch (...) { return vector<string>(); }
      return ret;
    }


    std::valarray<double> weights(const GenEvent& ge) {
      // std::valarray<double> rtn(ge.weights().size());
      // for (size_t i = 0; i < ge.weights().size(); ++i) rtn[i] = ge.weights()[i];
      // return rtn;
      return std::valarray<double>(ge.weights().data(), ge.weights().size());
    }


  }
}
