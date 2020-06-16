#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Tools/RivetHepMC.hh"
#include <fstream>

using namespace std;

int main(int argc, char* argv[]) {
  assert(argc > 1);

  Rivet::AnalysisHandler ah;
  Rivet::Log::setLevel("Rivet", Rivet::Log::DEBUG);

  // Specify the analyses to be used
  ah.addAnalysis("EXAMPLE");
  ah.addAnalyses({{ "MC_JETS", "EXAMPLE_CUTS", "EXAMPLE_SMEAR" }});

  shared_ptr<std::istream> file;
  shared_ptr<Rivet::HepMC_IO_type> reader = Rivet::HepMCUtils::makeReader("testApi.hepmc", file);
  std::shared_ptr<Rivet::GenEvent> evt = make_shared<Rivet::GenEvent>();
  double sum_of_weights = 0.0;
  while ( Rivet::HepMCUtils::readEvent(reader, evt) ) {
    // Analyse current event
    ah.analyze(*evt);
    sum_of_weights += evt->weights()[0];
  }

  ah.setCrossSection(make_pair(1.0, 0.1));

  ah.finalize();
  ah.writeData("out.yoda");

  return 0;
}
