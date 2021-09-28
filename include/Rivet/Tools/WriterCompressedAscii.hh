// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_WRITERCOMPRESSEDASCII_H
#define HEPMC3_WRITERCOMPRESSEDASCII_H
///
/// @file  WriterCompressedAscii.h
/// @brief Definition of class \b WriterCompressedAscii
///
/// @class HepMC3::WriterCompressedAscii
/// @brief GenEvent I/O serialization for structured text files
///
/// @ingroup IO
///
#include "Rivet/Tools/RivetHepMC.hh"
#include "HepMC3/Writer.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenRunInfo.h"
#include <string>
#include <fstream>

namespace Rivet {

class WriterCompressedAscii : public HepMC3::Writer {

public:

  typedef HepMC3::GenRunInfo GenRunInfo;
  typedef HepMC3::FourVector FourVector;
  typedef HepMC3::ConstGenVertexPtr ConstGenVertexPtr;
  typedef HepMC3::ConstGenParticlePtr ConstGenParticlePtr;

  /// @brief Constructor
  /// @warning If file already exists, it will be cleared before writing
  WriterCompressedAscii(const std::string& filename,
                        shared_ptr<GenRunInfo> run = shared_ptr<GenRunInfo>());

  /// @brief Constructor from ostream
  WriterCompressedAscii(std::ostream& stream,
                        shared_ptr<GenRunInfo> run = shared_ptr<GenRunInfo>());
  
  /// @brief Destructor
  ~WriterCompressedAscii();

  /// @brief Write event to file
  ///
  /// @param[in] evt Event to be serialized
  void write_event(const GenEvent& evt);

  /// @brief Write the GenRunInfo object to file.
  void write_run_info();

  /// @brief Return status of the stream
  bool failed() { return (bool)m_file.rdstate(); }

  /// @brief Close file stream
  void close();

  /// @brief Use cartesian coordinates
  ///
  /// Momenta and positions will be written out as doubles using
  /// standard cartesian coordinates.
  void use_doubles() {
    m_use_integers = false;
  }

  /// @brief Use cylindical coordinates
  ///
  /// Momenta and positions will be written out as integers using
  /// eta-phi coordinates.
  void use_integers() {
    m_use_integers = true;
  }

  /// @brief Add a particle id to be stripped
  ///
  /// Specify the PDG id of a (unobservable) particle which will be
  /// attempted to be removed from the event befor writing.
  void add_stripid(long pdgid) {
    m_stripid.insert(pdgid);
  }

    /// @brief Remove onobservable particles
    ///
    /// Go through the event and try to take away all intermediate
    /// particles specified in add_stripid()
  void strip(GenEvent & e);

    /// @brief Set output double precision
    ///
    /// General output precision for double
    void set_precision(int prec) {
        m_precision = prec;
    }

    /// @brief Set output precision in phi
    ///
    /// Azimuth angles will be written out as integers corresponding to
    /// this precision
    void set_precision_phi(double prec) {
        m_precision_phi = prec;
    }

    /// @brief Set output precision in eta
    ///
    /// Pseudo repisities will be written out as integers corresponding to
    /// this precision
    void set_precision_eta(double prec) {
        m_precision_eta = prec;
    }

    /// @brief Set output precision in energy
    ///
    /// Energies will be written out as integers corresponding to
    /// this precision (in GeV)
    void set_precision_e(double prec) {
        m_precision_e = prec;
    }

    /// @brief Set output precision in mass
    ///
    /// Masses will be written out as integers corresponding to
    /// this precision (in GeV)
    void set_precision_m(double prec) {
        m_precision_m = prec;
    }

    /// @brief Return output precision for doubles.
    int precision() const {
      return m_precision;
    }

    /// @brief Return output precision for azimuth angles.
    double precision_phi() const {
      return m_precision_phi;
    }

    /// @brief Return output precision for pseudo rapisities.
    double precision_eta() const {
      return m_precision_eta;
    }

    /// @brief Return output precision for energies.
    double precision_e() const {
      return m_precision_e;
    }

    /// @brief Return output precision for masses.
    double precision_m() const {
      return m_precision_m;
    }

    /// @brief Internal function to calculate the pseudo rapidity.
    double psrap(const FourVector & p) const;

private:

    /// @brief Escape '\' and '\n' characters in string
    std::string escape(const std::string& s)  const;

    //@}


    /// @name Write helpers
    //@{

    /// @brief Inline function for writing positions
    void write_position(FourVector pos);

    /// @brief Inline function for writing momenta
    void write_momentum(FourVector p);

    /// @brief Inline function for writing mass of a particle
    void write_mass(ConstGenParticlePtr p);

    /// @brief Write vertex
    ///
    /// Helper routine for writing single vertex to file
    void write_vertex(ConstGenVertexPtr v);

    /// @brief Write particle
    ///
    /// Helper routine for writing single particle to file
  void write_particle(ConstGenParticlePtr p);

    /// @brief Helper function to access the root vertex
    ConstGenVertexPtr rootvertex() {
      vector<ConstGenParticlePtr> beams = m_current->beams();
      if ( beams.empty() ) return ConstGenVertexPtr();
      return beams[0]->production_vertex();
    }
    //@}

private:

  bool m_use_integers;        //!< Compress by using intergers and
                              //!  cylindrical coordinates

  std::ofstream m_file;       //!< Output file
  std::ostream* m_stream;     //!< Output stream

  double m_precision_phi;     //!< Output precision in phi
  double m_precision_eta;     //!< Output precision in eta
  double m_precision_e;       //!< Output precision energy
  double m_precision_m;       //!< Output precision mass
  int m_precision;            //!< General double output precision
  set<long> m_stripid;        //!< Strip matching intermediate particles
  map<long,long> m_masses;    //!< Keep track of masses being written.

  const GenEvent * m_current; //!< The event being written
  std::ostringstream os;      //!< Internal stream to write an entire event.

};


} // namespace HepMC3

#endif
