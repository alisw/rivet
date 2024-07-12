#include "Rivet/Tools/RivetHepMC.hh"
#include "Rivet/Tools/WriterCompressedAscii.hh"
#include "../src/Core/zstr/zstr.hpp"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/WriterAscii.h"
using namespace std;

int main(int argc, char** argv) {

  if (argc < 3) {
    cerr << "Usage: " << argv[0]
         << " <input-hepmcfile> <output-hepmcfile>" << endl;
    return 1;
  }

  int outputmode = 0;
  if ( argc >= 4 ) outputmode = atoi(argv[3]);
  Rivet::zstr::ifstream input(argv[1]);
  Rivet::zstr::ofstream output(argv[2]);

  std::shared_ptr<Rivet::HepMC_IO_type>
    reader = Rivet::HepMCUtils::makeReader(input);
  
  std::shared_ptr<Rivet::RivetHepMC::GenEvent>
    evt = make_shared<Rivet::RivetHepMC::GenEvent>();

  shared_ptr<HepMC3::Writer> writer;
  if ( outputmode == 0 )
    writer = make_shared<HepMC3::WriterAscii>(output);
  else {
    auto compressed = make_shared<Rivet::WriterCompressedAscii>(output);
    if ( outputmode >= 2 ) compressed->use_integers();
    if ( abs(outputmode) == 3 ) {
      compressed->add_stripid(21);
      compressed->add_stripid(-1);
      compressed->add_stripid(1);
      compressed->add_stripid(-2);
      compressed->add_stripid(2);
      compressed->add_stripid(-3);
      compressed->add_stripid(3);
    }
    if ( abs(outputmode) == 4 ) {
      compressed->add_stripid(-22);
    }
    writer = compressed;
  }

  while(reader && Rivet::HepMCUtils::readEvent(reader, evt) ) {
    //    cout << evt->event_number() << endl;
    writer->write_event(*evt);
  }
  
  return 0;
}


