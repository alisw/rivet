#include "Rivet/Tools/RivetHepMC.hh"
#include "Rivet/Tools/WriterCompressedAscii.hh"
#include "../src/Core/zstr/zstr.hpp"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/WriterAscii.h"
using namespace std;

int main(int argc, char** argv) {

  string ifile;
  string ofile;
  bool etaphi = false;
  bool strip = false;
  bool userivet = false;
  bool help = false;
  double pphi = 0.0;
  double peta = 0.0;
  double pe = 0.0;
  double pm = 0.0;
  double po = 0.0;

  for ( int iarg = 1; iarg < argc; ++iarg ) {
    string arg = argv[iarg];
    if ( arg == "-s" )
      userivet = strip = etaphi = true;
    else if ( arg == "-S" )
      userivet = strip = true;
    else if ( arg == "-c" )
      userivet = etaphi = true;
    else if ( arg == "-r" )
      userivet = true;
    else if ( arg == "-p" ) {
      if ( iarg + 1 >= argc ) {
        help = true;
        break;
      }
      arg = argv[++iarg];
      if ( arg.substr(0,4) == "phi=" ) pphi = stod(arg.substr(4));
      if ( arg.substr(0,4) == "eta=" ) peta = stod(arg.substr(4));
      if ( arg.substr(0,2) == "e=" ) pe = stod(arg.substr(2));
      if ( arg.substr(0,2) == "m=" ) pm = stod(arg.substr(2));
      if ( arg.substr(0,2) == "p=" ) po = stod(arg.substr(2));
      userivet = etaphi = true;
    }
    else if ( arg == "-h" )
      help = true;
    else {
      if ( ifile.empty() ) ifile = arg;
      else if ( ofile.empty() ) ofile = arg;
      else {
        cout << "Unknown argument '" << arg << "'" << endl;
        help = true;
      }
    }
  }
  if ( ofile.empty() ) help = true;
  if ( help )  {
    cout << "Usage: " << argv[0]
         << " [options] <input-hepmcfile> <output-hepmcfile>" << endl;
    cout << "  where options are one or more of" << endl
         << "   -r: rivet internal hepmc output" << endl
         << "   -c: compressed hepmc output (implies -r)" << endl
         << "   -s: strips hepmc from unobservable (some) particles (implies -c)" << endl
         << "   -p type=prec: precision in compressed hepmc output (implies -c)" << endl
         << "      allowed types are phi, eta, e and m. " << endl
         << "   -h: write this message and exit." << endl;

    return 1;
  }


  
  std::shared_ptr<std::istream> istr;
  shared_ptr<ostream> output;
  if ( ofile.substr(ofile.length() - 3) == ".gz" ||
       ofile.substr(ofile.length() - 6) == ".hepmz" )
    output = make_shared<Rivet::zstr::ofstream>(ofile);
  else
    output = make_shared<ofstream>(ofile);

  shared_ptr<Rivet::HepMC_IO_type>
    reader = Rivet::HepMCUtils::makeReader(ifile, istr);
  
  shared_ptr<Rivet::RivetHepMC::GenEvent>
    evt = make_shared<Rivet::RivetHepMC::GenEvent>();

  shared_ptr<HepMC3::Writer> writer;
  if ( userivet ) {
    auto compressed = make_shared<Rivet::WriterCompressedAscii>(*output);
    if ( etaphi ) {
      compressed->use_integers();
      if ( pphi > 0.0 ) compressed->set_precision_phi(pphi);
      if ( peta > 0.0 ) compressed->set_precision_eta(peta);
      if ( pe > 0.0 ) compressed->set_precision_e(pe);
      if ( pm > 0.0 ) compressed->set_precision_m(pm);
      if ( po > 0.0 ) compressed->set_precision(po);
    }
    if ( strip ) {
      compressed->add_stripid(21);
      compressed->add_stripid(-1);
      compressed->add_stripid(1);
      compressed->add_stripid(-2);
      compressed->add_stripid(2);
      compressed->add_stripid(-3);
      compressed->add_stripid(3);
    }
    writer = compressed;
  } else {
    writer = make_shared<HepMC3::WriterAscii>(*output);
  }

  while(reader && Rivet::HepMCUtils::readEvent(reader, evt) ) {
    writer->write_event(*evt);
  }
  
  return 0;
}


