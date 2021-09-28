// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
///
/// @file WriterCompressedAscii.cc
/// @brief Implementation of \b class WriterCompressedAscii
///
#include "Rivet/Tools/WriterCompressedAscii.hh"

#include "HepMC3/Version.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Units.h"
#include <cstring>
#include <algorithm>//min max for VS2017

#if HEPMC3_VERSION_CODE >= 3002000
#define ERROR(x) HEPMC3_ERROR(x)
#define WARNING(x) HEPMC3_WARNING(x)
#define DEBUG(x,y) HEPMC3_DEBUG(x,y)
#endif

namespace Rivet {

using namespace HepMC3;


WriterCompressedAscii::WriterCompressedAscii(const std::string &filename, shared_ptr<GenRunInfo> run)
  : m_use_integers(false),
    m_file(filename),
    m_stream(&m_file),
    m_precision_phi(0.0001),
    m_precision_eta(0.0001),
    m_precision_e(0.001),
    m_precision_m(0.000001),
    m_precision(5),
    m_current(0) {
  set_run_info(run);
  if ( !m_file.is_open() ) {
      ERROR( "WriterCompressedAscii: could not open output file: "<<filename )
  } else {
    m_file << "HepMC::Version " << version() << std::endl;
    m_file << "HepMC::Asciiv3-START_EVENT_LISTING" << std::endl;
    if ( run_info() ) write_run_info();
  }
}


WriterCompressedAscii::WriterCompressedAscii(std::ostream &stream, shared_ptr<GenRunInfo> run)
  : m_use_integers(false),
    m_file(),
    m_stream(&stream),
    m_precision_phi(0.0001),
    m_precision_eta(0.0001),
    m_precision_e(0.001),
    m_precision_m(0.000001),
    m_precision(5),
    m_current(0) {
  set_run_info(run);
  (*m_stream) << "HepMC::Version " << version() << std::endl;
  (*m_stream) << "HepMC::CompressedAsciiv3-START_EVENT_LISTING" << std::endl;
  if ( run_info() ) write_run_info();
}


WriterCompressedAscii::~WriterCompressedAscii() {
    close();
}


void WriterCompressedAscii::write_event(const GenEvent &evt) {
  
  if ( !m_stripid.empty() ) {
    GenEvent e = evt;
    // cout << "#beams " << e.beams().size() << endl;
    // for ( auto bp : evt.beams() ) {
    //   GenParticlePtr nb = e.particles()[bp->id() - 1];
    //   cout << "Beam: " << bp->id() << " " << bp->pid() << " "
    //        << (bp->production_vertex()? bp->production_vertex()->id(): 999) << endl;
    //   cout << "NewBeam: " << nb->id() << " " << nb->pid() << " "
    //        << (nb->production_vertex()? nb->production_vertex()->id(): 999) << endl;
    // }

    strip(e);
    set<long> saveid;
    swap(m_stripid, saveid);
    write_event(e);
    swap(saveid, m_stripid);
    return;
  }

  m_current = &evt;
  m_masses.clear();
  os.clear();
  os.str(std::string());
  if ( !run_info() ) {
    set_run_info(evt.run_info());
    write_run_info();
  } else {
    if ( evt.run_info() && run_info() != evt.run_info() ) {
      WARNING( "WriterCompressedAscii::write_event: GenEvents contain "
               "different GenRunInfo objects from - only the "
               "first such object will be serialized." )
        }
  }

  // Write event info
  os << "E " << evt.event_number()
     << " " << evt.vertices().size()
     << " " << evt.particles().size();
  // Write event position if not zero
  const FourVector &pos = evt.event_pos();
  if ( !pos.is_zero() ) write_position(pos);
  os << endl;

  // Write conversions made for double -> ints
  if ( m_use_integers )
    os << "C " << m_precision_phi << " " << m_precision_eta
       << " " << m_precision_e << " " << m_precision_m << endl;
  else
    os << "C 0 0 0 0" << endl;

  // Write weight values if present
  if ( evt.weights().size() ) {
    os << "W";
    for (auto  w: evt.weights())
      os << " " << w;
    os << endl;
  }

  // Write attributes
  for ( auto vt1: evt.attributes() ) {
    for ( auto vt2: vt1.second ) {

      string st;
      bool status = vt2.second->to_string(st);

      if( !status ) {
        WARNING( "WriterCompressedAscii::write_event: problem "
                 "serializing attribute: "<<vt1.first )
          }
      else {
        os << "A " << vt2.first
           << " " << vt1.first
           << " " << escape(st) << endl;
      }
    }
  }

  set<ConstGenVertexPtr> done;

  // Print particles
  for(ConstGenParticlePtr p: evt.particles() ) {

    // Check to see if we need to write a vertex first
    ConstGenVertexPtr v = p->production_vertex();

    if ( !v ) cout << "WARMING particle " << p->id()
                   << " has no productions vertex!" << endl;
    else if ( v->id() < 0 && done.insert(v).second )
      write_vertex(v);
    
    write_particle(p);
  }

  bool orphans = false;
  for(ConstGenVertexPtr v: evt.vertices() ) {
    if ( done.insert(v).second ) {
      orphans = true;
      write_vertex(v);
    }
  }
  if ( orphans ) cout << "WARMING found unconnected vertices" << endl;
  
  (*m_stream) << os.str();

}

string WriterCompressedAscii::escape(const string& s)  const {
  string ret;
  ret.reserve( s.length()*2 );
  for ( string::const_iterator it = s.begin(); it != s.end(); ++it ) {
    switch ( *it ) {
    case '\\': ret += "\\\\"; break;
    case '\n': ret += "\\|"; break;
    default: ret += *it;
    }
  }
  return ret;
}

void WriterCompressedAscii::write_vertex(ConstGenVertexPtr v) {

  os << "V " << v->id() << " " << v->status() << " ";

  bool first = true;
  for ( auto p : v->particles_in() ) {
    os << (first? '[': ',') << p->id();
    first = false;
  }
  os << "]";
  
  const FourVector &pos = v->position();
  if ( !pos.is_zero() ) write_position(pos);

  os << endl;
  
}


void WriterCompressedAscii::write_run_info() {

  // If no run info object set, create a dummy one.
  if ( !run_info() ) set_run_info(make_shared<GenRunInfo>());

  vector<string> names = run_info()->weight_names();

  if ( !names.empty() ) {
    string out = names[0];
    for ( int i = 1, N = names.size(); i < N; ++i )
      out += "\n" + names[i];
    os << "W " << escape(out) << endl;
  }
  
  for ( int i = 0, N = run_info()->tools().size(); i < N; ++i  ) {
    string out = "T " + run_info()->tools()[i].name + "\n"
      + run_info()->tools()[i].version + "\n"
      + run_info()->tools()[i].description;
    os << escape(out) << endl;
  }


  for ( auto att: run_info()->attributes() ) {
    string st;
    if ( ! att.second->to_string(st) ) {
      WARNING ("WriterCompressedAscii::write_run_info: problem serializing attribute: "<< att.first )
        }
    else {
      os << "A " << att.first << " " << escape(st) << endl;

    }
  }
}

void WriterCompressedAscii::
write_particle(ConstGenParticlePtr p) {

  ConstGenVertexPtr vp = p->production_vertex();

  os << "P " << p->id()
     << " "  << (vp? vp->id(): 0)
     << " "  << p->pid();
  write_momentum(p->momentum());
  write_mass(p);
  os << " " << p->status() << endl;

}

void WriterCompressedAscii::close() {
  std::ofstream* ofs = dynamic_cast<std::ofstream*>(m_stream);
  if (ofs && !ofs->is_open()) return;
  (*m_stream) << "HepMC::CompressedAsciiv3-END_EVENT_LISTING" << endl << endl;
  if (ofs) ofs->close();
}

double WriterCompressedAscii::psrap(const FourVector & p) const {
  static const double MAXETA = 100.0;
  static const double MAXLOG = exp(-MAXETA);
  double nom = p.p3mod() + abs(p.pz());
  if ( nom <= 0.0 ) return 0.0;
  double den = max(p.perp(), nom*MAXLOG);
  return p.pz() > 0? log(nom/den): -log(nom/den);
}

void WriterCompressedAscii::write_momentum(FourVector p) {

  Units::convert(p, m_current->momentum_unit(), Units::GEV);
  
  if ( m_use_integers ) {
    long ie = long(round(p.e()/precision_e()));
    // Avoid zero momentum particles
    if ( ie == 0 && p.e() != 0.0 )
      os << " " << p.e()/precision_e();
    else
      os << " " << ie;
    os << " " << long(round(psrap(p)/precision_eta()))
       << " " << long(round(p.phi()/(M_PI*precision_phi())));
    return;
  }

  std::ostringstream osse;
  osse << std::scientific << setprecision(precision())
       << " " << p.px()
       << " " << p.py()
       << " " << p.pz()
       << " " << p.e();
  os << osse.str();

}

void WriterCompressedAscii::write_mass(ConstGenParticlePtr p) {

  double m = p->generated_mass();
  if ( m_current->momentum_unit() != Units::GEV ) m /= 1000.0;
  if ( m_use_integers ) {
    long im = long(round(m/precision_m()));
    auto pm = m_masses.find(p->pid());
    if ( pm == m_masses.end() || pm->second != im ) {
      os << " " << im;
      m_masses[p->pid()] = im;
    } else {
      os << " *";
    }
    return;
  }
  std::ostringstream osse;
  osse << std::scientific << setprecision(precision())
       << " " << m;
  os << osse.str();

}

void WriterCompressedAscii::write_position(FourVector pos) {

  Units::convert(pos, m_current->length_unit(), Units::MM);
  std::ostringstream osse;
  osse << std::scientific << setprecision(precision());

  if ( m_use_integers ) {
    osse << " @ " << long(psrap(pos)/precision_eta())
         << " " << long(pos.phi()/(M_PI*precision_phi()))
         << " " << pos.p3mod()
         << " " << pos.t();
  } else {
    osse << " @ " << pos.x()
         << " " << pos.y()
         << " " << pos.z()
         << " " << pos.t();
  }
  os << osse.str();
}

void WriterCompressedAscii::strip(GenEvent & e) {
  HepMCUtils::strip(e, m_stripid);
}

}
