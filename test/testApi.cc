#include "Rivet/AnalysisHandler.hh"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"

using namespace std;

int main() {
  Rivet::AnalysisHandler ah;
  Rivet::Log::setLevel("Rivet", Rivet::Log::TRACE);

  // Specify the analyses to be used
  ah.addAnalysis("EXAMPLE");
  vector<string> moreanalyses(1, "MC_JETS");
  ah.addAnalyses(moreanalyses);
  ah.addAnalysis("EXAMPLE_CUTS");

  std::istream* file = new std::fstream("testApi.hepmc", std::ios::in);
  HepMC::IO_GenEvent hepmcio(*file);
  HepMC::GenEvent* evt = hepmcio.read_next_event();
  double sum_of_weights = 0.0;
  while (evt) {
    // Analyse current event
    ah.analyze(*evt);
    sum_of_weights += evt->weights()[0];

    // Clean up and get next event
    delete evt; evt = 0;
    hepmcio >> evt;
  }
  delete file; file = 0;

  ah.setCrossSection(1.0);
  ah.setSumOfWeights(sum_of_weights); ///< Not necessary, but allowed
  ah.finalize();
  ah.writeData("out.yoda");

  return 0;
}
