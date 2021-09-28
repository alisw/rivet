// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
///
/// @file ReaderCompressedAscii.cc
/// @brief Implementation of \b class ReaderCompressedAscii
///
#include "Rivet/Tools/ReaderCompressedAscii.hh"

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Units.h"
#include <cstring>
#include <sstream>

#if HEPMC3_VERSION_CODE >= 3002000
#define ERROR(x) HEPMC3_ERROR(x)
#define WARNING(x) HEPMC3_WARNING(x)
#define DEBUG(x,y) HEPMC3_DEBUG(x,y)
#endif

namespace Rivet {

using namespace HepMC3;


ReaderCompressedAscii::ReaderCompressedAscii(const string &filename)
  : m_file(filename), m_stream(0), m_evt(0), m_precision_phi(0.001),
    m_precision_eta(0.001), m_precision_e(0.001), m_precision_m(0.000001),
    m_using_integers(false) {
  if( !m_file.is_open() ) {
    ERROR( "ReaderCompressedAscii: could not open input file: "<<filename )
      }
  m_stream = &m_file;
  set_run_info(make_shared<GenRunInfo>());
}


// Ctor for reading from stdin
ReaderCompressedAscii::ReaderCompressedAscii(std::istream & stream)
  : m_stream(&stream), m_evt(0), m_precision_phi(0.001),
    m_precision_eta(0.001), m_precision_e(0.001), m_precision_m(0.000001),
    m_using_integers(false) {
  if( !m_stream ) {
    ERROR( "ReaderCompressedAscii: could not open input stream " )
      }
  set_run_info(make_shared<GenRunInfo>());
}



ReaderCompressedAscii::~ReaderCompressedAscii() { }


bool ReaderCompressedAscii::read_event(GenEvent &evt) {

  if ( failed() ) return false;

  bool               parsed_event_header    = false;
  bool               is_parsing_successful  = true;
  pair<int,int> vertices_and_particles(0,0);
  std::string line;

  m_evt = &evt;

  evt.clear();
  evt.set_run_info(run_info());
  m_masses.clear();
  m_vertices.clear();
  m_particles.clear();
  m_ppvx.clear();
  m_vpin.clear();
  //
  // Parse event, vertex and particle information
  //
  while(!failed()) {

    std::getline(*m_stream, line);
    is.clear();
    is.str(line);
is.get(); // Remove the first character from the stream.
    if ( line.empty() ) continue;

    if ( line.substr(0, 5) == "HepMC" ) {
      if ( line.substr(0, 14) != "HepMC::Version" &&
           line.substr(0, 24) != "HepMC::CompressedAsciiv3" ) {
            WARNING( "ReaderCompressedAscii: found unsupported expression "
                     "in header. Will close the input." )
              std::cout << line << std::endl;
      }
      continue;
    }

    switch( line[0] ) {
    case 'E':
      vertices_and_particles = parse_event_information();
      if (vertices_and_particles.second < 0) {
        is_parsing_successful = false;
      } else {
        is_parsing_successful = true;
        parsed_event_header   = true;
      }
      break;
    case 'V':
      is_parsing_successful = parse_vertex_information();
      break;
    case 'P':
      is_parsing_successful = parse_particle_information();
      break;
    case 'W':
      if ( parsed_event_header )
        is_parsing_successful = parse_weight_values();
      else
        is_parsing_successful = parse_weight_names();
      break;
    case 'C':
      is_parsing_successful = parse_precision();
      break;
    case 'T':
      is_parsing_successful = parse_tool();
      break;
    case 'A':
      if ( parsed_event_header )
        is_parsing_successful = parse_attribute();
      else
        is_parsing_successful = parse_run_attribute();
      break;
    default:
      WARNING( "ReaderCompressedAscii: skipping unrecognised prefix: "
               << line[0] )
        is_parsing_successful = true;
      break;
    }

    if( !is_parsing_successful ) break;

    // Check for next event
    if ( parsed_event_header &&
         ( m_stream->peek() == 'E' || m_stream->peek() == 'H') )break;
  }

    // Set the production vertex for all particles.
    for ( int ip = 0, Np = m_particles.size(); ip < Np; ++ip )
      if ( m_ppvx[ip] && m_vertices[m_ppvx[ip]] )
        m_vertices[m_ppvx[ip]]->add_particle_out(m_particles[ip]);

    // Add the incoming particles to all vertices
    for ( auto iv : m_vertices )
      for ( auto ip : m_vpin[iv.first] ) iv.second->add_particle_in(m_particles[ip - 1]);

    // When all particles and vertices are connected we add all of them to the event.
    for ( auto p : m_particles ) evt.add_particle(p);
    for ( auto v : m_vertices ) evt.add_vertex(v.second);

  // Check if all particles and vertices were parsed
  if ((int)m_evt->particles().size() > vertices_and_particles.second ) {
    ERROR( "ReaderCompressedAscii: too many particles were parsed" )
      printf("%zu  vs  %i expected\n",m_evt->particles().size(),vertices_and_particles.second );
    is_parsing_successful = false;
    }
  if ((int)m_evt->particles().size() < vertices_and_particles.second ) {
        ERROR( "ReaderCompressedAscii: too few  particles were parsed" )
          printf("%zu  vs  %i expected\n",m_evt->particles().size(),vertices_and_particles.second );
        is_parsing_successful = false;
  }

  if ((int)m_evt->vertices().size()  > vertices_and_particles.first) {
    ERROR( "ReaderCompressedAscii: too many vertices were parsed" )
      printf("%zu  vs  %i expected\n",m_evt->vertices().size(),vertices_and_particles.first );
    is_parsing_successful =  false;
    }

  if ((int)m_evt->vertices().size()  < vertices_and_particles.first) {
    ERROR( "ReaderCompressedAscii: too few vertices were parsed" )
      printf("%zu  vs  %i expected\n",m_evt->vertices().size(),vertices_and_particles.first );
    is_parsing_successful =  false;
  }    
  // Check if there were errors during parsing
  if( !is_parsing_successful ) {
    ERROR( "ReaderCompressedAscii: event parsing failed. Returning empty event" )
    DEBUG( 1, "Parsing failed at line:" << endl << line )
      
    m_evt->clear();

    return false;
  }
  return true;
}


pair<int,int> ReaderCompressedAscii::parse_event_information() {
  static const pair<int,int>  err(-1,-1);
  pair<int,int>               ret(-1,-1);
  int                         event_no = 0;

  // event number
  if ( !(is >> event_no) ) return err;
  m_evt->set_event_number(event_no);

  // num_vertices
  if ( !(is>> ret.first) ) return err;

  // num_particles
  if ( !(is>> ret.second) ) return err;

  if ( !read_position() ) return err;

  DEBUG( 10, "ReaderCompressedAscii: E: "<<event_no<<" ("<<ret.first<<"V, "<<ret.second<<"P)" )

    return ret;
}


bool ReaderCompressedAscii::parse_weight_values() {

  vector<double> wts;
  double w;
  while ( is >> w ) wts.push_back(w);
  if ( run_info() && run_info()->weight_names().size() &&
       run_info()->weight_names().size() != wts.size() )
    throw std::logic_error(
      "ReaderCompressedAscii::parse_weight_values: "
      "The number of weights ("+
      std::to_string((long long int)(wts.size()))+
      ") does not match the  number weight names("
      + std::to_string((long long int)(run_info()->weight_names().size()))+
      ") in the GenRunInfo object");

  m_evt->weights() = wts;

  return true;
}


bool ReaderCompressedAscii::parse_precision() {
  if ( !(is >> m_precision_phi >>  m_precision_eta
         >>  m_precision_e >> m_precision_m) ) return false;
  m_using_integers = ( m_precision_phi > 0.0 );
  return true;
}

bool ReaderCompressedAscii::parse_vertex_information() {
  GenVertexPtr  data = make_shared<GenVertex>();

  int id = 0;
  if ( !(is >> id) ) return false;

  int status = 0;
  if  ( !(is >> status) ) return false;
  data->set_status( status );

  std::string incoming;
  if ( !(is >> incoming) ) return false;
  std::string::size_type i = std::string::npos;
  while ( ( i = incoming.find_first_of("[,]") ) != std::string::npos )
    incoming[i] = ' ';
  std::istringstream isin(incoming);
  int pin = 0;
  vector<int> vpin;
  while ( isin >> pin ) vpin.push_back(pin);

  if ( !read_position(data) ) return false;

  m_vertices[-id] = data;
  m_vpin[-id] = vpin;

  return true;
}


bool ReaderCompressedAscii::parse_particle_information() {
  GenParticlePtr  data = make_shared<GenParticle>();

  int id = 0;
  if ( !(is >> id) ) return false;

  int ivp = 0;
  if ( !(is >> ivp) ) return false;

  int pdgid = 0;
  if ( !(is >> pdgid) ) return false;
  data->set_pid(pdgid);

  if ( !read_momentum(data) ) return false;

  int status = 0;
  if ( !(is >> status) ) return false;
  data->set_status(status);

  m_particles.push_back(data);
  m_ppvx.push_back(-ivp);

  return true;
}


bool ReaderCompressedAscii::parse_attribute() {
  int id = 0;
  if ( !(is >> id ) ) return false;

  string name;
  if ( !(is >> name) ) return false;
  is.get();
  string contents;
  if ( !std::getline(is, contents) ) return false;
  shared_ptr<Attribute> att =
    make_shared<StringAttribute>(StringAttribute(unescape(contents)));

  m_evt->add_attribute(name, att, id);

  return true;
}

bool ReaderCompressedAscii::parse_run_attribute() {

  string name;
  if ( !(is >> name) ) return false;
  is.get();
  string contents;
  if ( !std::getline(is, contents) ) return false;
  shared_ptr<Attribute> att =
    make_shared<StringAttribute>(StringAttribute(unescape(contents)));

  run_info()->add_attribute(name, att);

  return true;

}


bool ReaderCompressedAscii::parse_weight_names() {

  vector<string> names;
  string name;
  while ( is >> name ) names.push_back(name);

  run_info()->set_weight_names(names);

  return true;

}

bool ReaderCompressedAscii::parse_tool() {

  std::string line;
  if ( !(is >> line) ) return false;
  line = unescape(line);

  GenRunInfo::ToolInfo tool;

  std::string::size_type pos = line.find("\n");
  tool.name = line.substr(0, pos);

  line = line.substr(pos + 1);
  pos = line.find("\n");
  tool.version = line.substr(0, pos);

  tool.description = line.substr(pos + 1);
  run_info()->tools().push_back(tool);

  return true;

}

bool ReaderCompressedAscii::read_position(GenVertexPtr v) {
  string at;
  if ( !(is >> at) ) return true;
  if ( at != "@" ) return false;

  if ( !m_using_integers ) {
    double x = 0.0, y = 0.0, z = 0.0, t = 0.0;
    if ( !(is >> x >> y >> z >> t) ) return false;
    FourVector pos(x, y, z, t);
    Units::convert(pos, Units::MM, m_evt->length_unit());
    v->set_position(pos);
    return true;
  }

  long ieta = 0;
  long iphi = 0;
  double p3mod = 0.0;
  double t = 0.0;
  if ( !(is >> ieta >> iphi >> p3mod >> t) ) return false;

  double eta = double(ieta)*m_precision_eta;
  double phi = double(iphi)*m_precision_phi*M_PI;
  double pt = p3mod/cosh(eta);
  FourVector pos(pt*cos(phi), pt*sin(phi), p3mod*tanh(eta), t);

  Units::convert(pos, Units::MM, m_evt->length_unit());
  v->set_position(pos);

  return true;

}

bool ReaderCompressedAscii::read_momentum(GenParticlePtr p) {
  if ( !m_using_integers ) {
    double px = 0.0, py = 0.0, pz = 0.0, e = 0.0, m = 0.0;
    if ( !(is >> px >> py >> pz >> e >> m) ) return false;
    FourVector pp(px, py, pz, e);
  
    if ( m_evt->momentum_unit() != Units::GEV ) {
      m *= 1000.0;
      Units::convert(pp, Units::GEV, m_evt->momentum_unit());
    }
    p->set_momentum(pp);
    p->set_generated_mass(m);
    return true;
  }
    
  long iphi = 0;
  long ieta = 0;
  double ie = 0;
  std::string sm;
  if ( !(is >> ie >> ieta >> iphi >> sm) ) return false;

  double m = 0.0;
  if ( sm == "*" ) {
    m = m_masses[p->pid()]*m_precision_m;
  } else {
    m = (m_masses[p->pid()] = stol(sm))*m_precision_m;
  }

  double e = double(ie)*m_precision_e;
  double m2 = ( m >= 0.0? m*m: -m*m );
  double p3mod = sqrt(max(e*e - m2, 0.0));
  double eta = double(ieta)*m_precision_eta;
  double phi = double(iphi)*m_precision_phi*M_PI;
  double pt = abs(eta) < 100.0? p3mod/cosh(eta): 0.0;
  FourVector pp(pt*cos(phi), pt*sin(phi), p3mod*tanh(eta), e);
  
  if ( m_evt->momentum_unit() != Units::GEV ) {
    m *= 1000.0;
    Units::convert(pp, Units::GEV, m_evt->momentum_unit());
  }

  p->set_momentum(pp);
  p->set_generated_mass(m);

  return true;
}

string ReaderCompressedAscii::unescape(const string& s) {
  string ret;
  ret.reserve(s.length());
  for ( string::const_iterator it = s.begin(); it != s.end(); ++it ) {
    if ( *it == '\\' ) {
      ++it;
      if ( *it == '|' )
        ret += '\n';
      else
        ret += *it;
    } else
      ret += *it;
  }

  return ret;
}


void ReaderCompressedAscii::close() {
    if( !m_file.is_open()) return;
    m_file.close();
}


} // namespace HepMC3
