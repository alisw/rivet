// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_READERCOMPRESSEDASCII_H
#define HEPMC3_READERCOMPRESSEDASCII_H
///
/// @file  ReaderCompressedAscii.h
/// @brief Definition of class \b ReaderCompressedAscii
///
/// @class HepMC3::ReaderCompressedAscii
/// @brief GenEvent I/O parsing for structured text files
///
/// @ingroup IO
///
#include "Rivet/Tools/RivetHepMC.hh"
#include "HepMC3/Reader.h"
#include "HepMC3/GenEvent.h"
#include <string>
#include <fstream>
#include <istream>

namespace Rivet {


class ReaderCompressedAscii : public HepMC3::Reader {

public:

  typedef HepMC3::GenParticlePtr GenParticlePtr;
  typedef HepMC3::GenVertexPtr GenVertexPtr;

public:

  /// @brief Constructor
  /// @warning If file already exists, it will be cleared before writing
  ReaderCompressedAscii(const std::string& filename);

  /// The ctor to read from stdin
  ReaderCompressedAscii(std::istream &);

  /// @brief Destructor
  ~ReaderCompressedAscii();

  /// @brief Load event from file
  ///
  /// @param[out] evt Event to be filled
  bool read_event(GenEvent& evt);

  /// @brief Return status of the stream
  bool failed() {
    return !(*m_stream);
  }

  /// @brief Close file stream
  void close();

private:

  /// @brief Unsecape '\' and '\n' characters in string
  std::string unescape(const std::string& s);

  /// @name Read helpers
  //@{

  /// @brief Parse event
  ///
  /// Helper routine for parsing event information
  /// @return vertices count and particles count for verification
  std::pair<int,int> parse_event_information();

  /// @brief Parse weight value lines
  ///
  /// Helper routine for parsing weight value information
  bool parse_weight_values();

  /// @brief Parse precision
  ///
  /// Helper routine for parsing precision information
  bool parse_precision();

  /// @brief Parse vertex
  ///
  /// Helper routine for parsing single event information
  bool parse_vertex_information();

  /// @brief Parse particle
  ///
  /// Helper routine for parsing single particle information
  bool parse_particle_information();

  /// @brief Parse attribute
  ///
  /// Helper routine for parsing single attribute information
  bool parse_attribute();

  /// @brief Parse run-level attribute.
  ///
  /// Helper routine for parsing single attribute information
  bool parse_run_attribute();

  /// @brief Parse run-level weight names.
  ///
  /// Helper routine for parsing a line with information about
  /// weight names.
  bool parse_weight_names();

  /// @brief Parse run-level tool information.
  ///
  /// Helper routine for parsing a line with information about
  /// tools being used.
  bool parse_tool();

  /// @brief Read position information
  ///
  /// Reads position information from the current line and sets the
  /// information in the given vertex. If no vertex is given the root
  /// vertiex is assumed.
  bool read_position(GenVertexPtr v = GenVertexPtr());

  /// @brief Read momentum information
  ///
  /// Reads momentum and mass information from the current line and
  /// sets the information in the given particle.
  bool read_momentum(GenParticlePtr p);

  //@}


private:

  std::ifstream m_file;       //!< Input file
  std::istream* m_stream;     //!< The stream being read from 

  std::istringstream is;      //!< A stream to read from the current line.

  GenEvent * m_evt;           //!< The event being read in.
  
  double m_precision_phi;     //!< Input precision in phi
  double m_precision_eta;     //!< Input precision in eta
  double m_precision_e;       //!< Input precision energy
  double m_precision_m;       //!< Input precision mass
  bool m_using_integers;      //!< Reading integers

  map<long,long> m_masses;    //!< Keep track of masses being read.

  /// Keep track of read particles
  vector<GenParticlePtr>  m_particles;
  /// Keep track of read particles
  vector<int>  m_ppvx;
  /// Keep track of read vertices
  map<int,GenVertexPtr> m_vertices;
  /// Keep track of incoming particles to vertices
  map<int, std::vector<int> > m_vpin;

  /** @brief Store attributes global to the run being written/read. */
  std::map< std::string, shared_ptr<HepMC3::Attribute> > m_global_attributes;

};


} // namespace HepMC3

#endif
