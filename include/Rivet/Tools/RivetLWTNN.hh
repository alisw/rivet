// -*- C++ -*-
#ifndef RIVET_RivetLWTNN_HH
#define RIVET_RivetLWTNN_HH

#include "Rivet/Tools/RivetPaths.hh"
#include "lwtnn/LightweightNeuralNetwork.hh"
#include "lwtnn/LightweightGraph.hh"
#include "lwtnn/Exceptions.hh"
#include "lwtnn/parse_json.hh"
#include <fstream>

namespace Rivet {
  using namespace std;


  /// Read a LWT DNN config from the JSON path
  lwt::JSONConfig readLWTNNConfig(const string& jsonpath) {
    ifstream input;
    try {
      // Note: a failed read here may fail quietly, and cause the filestream to
      // go bad, making it look like the hepmc event-read has failed.
      input = std::ifstream(jsonpath);
      return lwt::parse_json(input);
    } catch (lwt::LightweightNNException &e) {
      input.close();
      throw IOError("Error loading LWTNN JSON config");
    }
  }


  /// @brief Read a LWT Graph config from the JSON path
  ///
  /// Note graph here means "not linear" rather than a GNN
  lwt::GraphConfig readLWTNNGraphConfig(const string& jsonpath) {
    ifstream input;
    try {
      // Note: a failed read here may fail quietly, and cause the filestream to
      // go bad, making it look like the hepmc event-read has failed.
      input = std::ifstream(jsonpath);
      return lwt::parse_json_graph(input);
    } catch (lwt::LightweightNNException &e) {
      input.close();
      throw IOError("Error loading LWTNN JSON config");
    }
  }

  /// Make a LWT DNN from the JSON config object
  std::unique_ptr<lwt::LightweightNeuralNetwork> mkLWTNN(const lwt::JSONConfig& jsonconfig) {
    try {
      return std::make_unique<lwt::LightweightNeuralNetwork>(jsonconfig.inputs, jsonconfig.layers, jsonconfig.outputs);
    } catch (lwt::LightweightNNException &e) {
      throw IOError("Error initialising from LWTNN JSON config");
    }
  }

  /// @brief Make a LWT Graph from the JSON config object
  ///
  /// Note graph here means "not linear" rather than a GNN
  std::unique_ptr<lwt::LightweightGraph> mkGraphLWTNN(const lwt::GraphConfig& graphconfig) {
    try {
      return std::make_unique<lwt::LightweightGraph>(graphconfig);
    } catch (lwt::LightweightNNException &e) {
      throw IOError("Error initialising from LWTNN JSON config");
    }
  }


  /// Make a LWT DNN direct from the JSON config path
  std::unique_ptr<lwt::LightweightNeuralNetwork> mkLWTNN(const string& jsonpath) {
    lwt::JSONConfig config = readLWTNNConfig(jsonpath);
    return mkLWTNN(config);
  }

  /// @brief Make a LWT graph direct from the JSON config path
  ///
  /// Note graph here means "not linear" rather than a GNN
  std::unique_ptr<lwt::LightweightGraph> mkGraphLWTNN(const string& jsonpath) {
    lwt::GraphConfig config = readLWTNNGraphConfig(jsonpath);
    return mkGraphLWTNN(config);
  }

}

#endif
