// -*- C++ -*-

#include "../Core/zstr/zstr.hpp"
#include "HepMC3/GenCrossSection.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/ReaderAsciiHepMC2.h"
#include "HepMC3/ReaderFactory.h"
#include "Rivet/Tools/ReaderCompressedAscii.hh"
#include "Rivet/Tools/RivetHepMC.hh"
#include <cassert>

/// @todo Delete this and require > 3.2 at some point
#if HEPMC3_VERSION_CODE < 3002000

namespace HepMC3 {


  /// Deduce the type of input stream based on its content and return an appropriate Reader
  std::shared_ptr<Reader> deduce_reader(std::istream& stream) {
    std::vector<std::string> head;
    head.push_back("");
    size_t back = 0;
    size_t backnonempty = 0;
    while ((back < 200 && backnonempty < 100) && stream) {
      char c = stream.get();
      back++;
      if (c == '\n') {
        if (head.back().length() != 0) head.push_back("");
      } else {
        head.back() += c;
        backnonempty++;
      }
    }
    if (!stream) {
      printf("Info in deduce_reader: input stream is too short or invalid.\n");
      return shared_ptr<Reader>(nullptr);
    }

    for (size_t i = 0; i < back; i++) stream.unget();

    if (strncmp(head.at(0).c_str(), "HepMC::Version", 14) == 0 && strncmp(head.at(1).c_str(), "HepMC::Asciiv3", 14) == 0) {
      printf("Info in deduce_reader: Attempt ReaderAscii\n");
      return std::shared_ptr<Reader>((Reader*)(new ReaderAscii(stream)));
    }

    if (strncmp(head.at(0).c_str(), "HepMC::Version", 14) == 0 && strncmp(head.at(1).c_str(), "HepMC::IO_GenEvent", 18) == 0) {
      printf("Info in deduce_reader: Attempt ReaderAsciiHepMC2\n");
      return std::shared_ptr<Reader>((Reader*)(new ReaderAsciiHepMC2(stream)));
    }

    #if HEPMC3_VERSION_CODE >= 3001002
    if (strncmp(head.at(0).c_str(), "<LesHouchesEvents", 17) == 0) {
      printf("Info in deduce_reader: Attempt ReaderLHEF\n");
      return std::shared_ptr<Reader>((Reader*)(new ReaderLHEF(stream)));
    }
    printf("Info in deduce_reader: Attempt ReaderHEPEVT\n");
    std::stringstream st_e(head.at(0).c_str());
    char attr = ' ';
    bool HEPEVT = true;
    int m_i, m_p;
    while (true) {
      if (!(st_e >> attr)) {
        HEPEVT = false;
        break;
      }
      if (attr == ' ') continue;
      if (attr != 'E') {
        HEPEVT = false;
        break;
      }
      HEPEVT = static_cast<bool>(st_e >> m_i >> m_p);
      break;
    }
    if (HEPEVT) return std::shared_ptr<Reader>((Reader*)(new ReaderHEPEVT(stream)));
    printf("Info in deduce_reader: All attempts failed\n");
    #endif
    return shared_ptr<Reader>(nullptr);
  }


}

#endif



namespace Rivet {
  namespace HepMCUtils {


    ConstGenParticlePtr getParticlePtr(const RivetHepMC::GenParticle& gp) {
      return gp.shared_from_this();
    }

    std::vector<ConstGenParticlePtr> particles(ConstGenEventPtr ge) {
      return ge->particles();
    }

    std::vector<ConstGenParticlePtr> particles(const GenEvent* ge) {
      assert(ge);
      return ge->particles();
    }

    std::vector<ConstGenVertexPtr> vertices(ConstGenEventPtr ge) {
      return ge->vertices();
    }

    std::vector<ConstGenVertexPtr> vertices(const GenEvent* ge) {
      assert(ge);
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
      return gp->id();
    }


    std::pair<ConstGenParticlePtr, ConstGenParticlePtr> beams(const GenEvent* ge) {
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

      // First scan forward and check if there is some hint as to what
      // kind of file we are looking att.
      int ntry = 10;
      std::string header;
      int filetype = -1;
      while (ntry) {
        std::getline(istr, header);
        if (header.empty()) continue;
        if (header.substr(0, 34) == "HepMC::Asciiv3-START_EVENT_LISTING") {
          filetype = 3;
          break;
        }
        if (header.substr(0, 44) == "HepMC::CompressedAsciiv3-START_EVENT_LISTING") {
          filetype = 4;
          break;
        }
        if (header.substr(0, 38) == "HepMC::IO_GenEvent-START_EVENT_LISTING") {
          filetype = 2;
          break;
        }
        ntry -= 1;
      }

      if (filetype == 3) {
        ret = make_shared<RivetHepMC::ReaderAscii>(istr);
      } else if (filetype == 4) {
        ret = make_shared<Rivet::ReaderCompressedAscii>(istr);
      } else if (filetype == 2) {
        ret = make_shared<RivetHepMC::ReaderAsciiHepMC2>(istr);
      }

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

      // Now we try to reopen the file and see if we can read something.
      if (errm) *errm += "Could not deduce file format. Will ask HepMC3 to try. ";
      ret = RivetHepMC::deduce_reader(filename);

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


    pair<double, double> crossSection(const GenEvent& ge) {
      // Work-around since access functions are not const.
      HepMC3::GenCrossSection xs = *ge.cross_section();
      return make_pair(xs.xsec(), xs.xsec_err());
    }


    vector<string> weightNames(const GenEvent& ge) {
      vector<string> ret;
      try {
        #if HEPMC3_VERSION_CODE >= 3002000
        ret = ge.weight_names();
        #else
        ret = ge.weight_names("");
        #endif
      } catch (...) { return vector<string>(); }
      return ret;
    }


    std::valarray<double> weights(const GenEvent& ge) {
      return std::valarray<double>(&ge.weights()[0], ge.weights().size());
    }


  }
}
