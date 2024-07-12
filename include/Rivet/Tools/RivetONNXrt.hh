// -*- C++ -*-
#ifndef RIVET_RivetONNXrt_HH
#define RIVET_RivetONNXrt_HH

#include <iostream>

#include "Rivet/Tools/RivetPaths.hh"
#include "onnxruntime_cxx_api.h"



namespace Rivet {
  using namespace std;
    /// @brief Simple object to automatically take care of basic ONNX networks
    ///
    /// Assumes one input/output node (note a node is not a neuron - a node is a single 
    /// tensor of arbitrary dimension size)
    /// See examples/EXAMPLE_ONNX.cc for how to use this.
    class RivetONNXrt{

      public:
      /// constructor
      RivetONNXrt(const string& filename, const string& runname = "RivetONNXrt"){
        //Set some ORT variables that need to be kept in memory
        _env = std::make_unique<Ort::Env>(ORT_LOGGING_LEVEL_WARNING,runname.c_str() );
        Ort::SessionOptions sessionopts; //todo - check this is allowed to go out of scope.
        _session = std::make_unique<Ort::Session> (*_env, filename.c_str(), sessionopts);

        //Get network hyper-params and store them (input, output shape, etc.) in the class.
        getNetworkInfo();

        MSG_DEBUG(*this);
      }

      //Default constructor with no args causes problems - delete it.
      RivetONNXrt() = delete;

      /// given an input vector, populate an output vector
      void compute(std::vector<float> &inputs, std::vector<float>& outputs){
        //Create ONNXrt inputs
        // create input tensor object from data values
        auto memory_info = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
        auto input_tensor = Ort::Value::CreateTensor<float>(memory_info, inputs.data(),
                                                              inputs.size(), _inputNodeDims.data(), 2);//TODO: Understand this magic number

        //Work around for stupid pointer stuff
        const char* temp_inputNodeName = _inputNodeName.c_str();
        const char* temp_outputNodeName = _outputNodeName.c_str();
        
        auto output_tensors =
          _session->Run(Ort::RunOptions{nullptr}, &temp_inputNodeName, &input_tensor, 
                              1, &temp_outputNodeName, 1);//"magic" 1's reflect the number of input/output nodes.

        float* floatarr = output_tensors.front().GetTensorMutableData<float>();

        outputs.clear();
        // TODO (longer-term): Generalise for different shape output arrays.
        outputs.assign(floatarr, floatarr+_outputNodeDims[0]);                                                     
      }

      void getNetworkInfo(){
        Ort::AllocatorWithDefaultOptions allocator;
        //Magic 0's are the fact we only have 1 input/output node
        auto input_name = _session->GetInputNameAllocated(0, allocator);
        _inputNodeName = input_name.get();
        auto in_type_info = _session->GetInputTypeInfo(0);
        auto in_tensor_info = in_type_info.GetTensorTypeAndShapeInfo();
        _inType = in_tensor_info.GetElementType();//TODO: Use this for SFINAE
        _inputNodeDims = in_tensor_info.GetShape();
        // Check for -1's: This is an artifact of batch size issues.
        // TODO: It's interesting that this is problematic in C++ and not in python.
        // I'd like to know why.
        for (auto& i : _inputNodeDims){
          if (i < 0)
            i = abs(i);
        }

        auto output_name = _session->GetOutputNameAllocated(0, allocator);
        _outputNodeName = output_name.get();
        auto out_type_info = _session->GetInputTypeInfo(0);
        auto out_tensor_info = in_type_info.GetTensorTypeAndShapeInfo();
        _outType = out_tensor_info.GetElementType();//TODO: Use this for SFINAE
        _outputNodeDims = out_tensor_info.GetShape();
        // Check for -1's: This is an artifact of batch size issues.
        // TODO: It's interesting that this is problematic in C++ and not in python.
        // I'd like to know why.
        for (auto& i : _outputNodeDims){
          if (i < 0)
            i = abs(i);
        }

        //Do some basic sanity checks:
        if (_session->GetInputCount() != 1 || _session->GetInputCount() != 1){
          throw("RivetONNXrt class cannot deal with multiple input/output nodes");
        }
      }

      /// Printing function for debugging.
      friend ostream& operator <<(std::ostream& os, const RivetONNXrt& rort){
        os << "RivetONNXrt Network Summary: \n";
        os << "Input name: " << rort._inputNodeName << "; Output name: " << rort._outputNodeName;
        os << "\nInput dimensions: (";
        for (size_t i = 0; i < rort._inputNodeDims.size()-1; ++i){
          os << rort._inputNodeDims[i] << ", ";
        }
        os << rort._inputNodeDims[rort._inputNodeDims.size() - 1] << ")\n";
        os << "Input Type (in ONNX enum form): " << rort._inType << "\n";
        os << "\nOutput dimensions: (";
        for (size_t i = 0; i < rort._outputNodeDims.size()-1; ++i){
          os << rort._outputNodeDims[i] << ", ";
        }
        os << rort._outputNodeDims[rort._outputNodeDims.size() - 1] << ")\n";
        os << "Output Type (in ONNX enum form): " << rort._outType << "\n";
        return os;
      }

      /// Logger
      Log& getLog() const {
        string logname = "Rivet.RivetONNXrt";
        return Log::getLog(logname);
      }

      private:
      /// ORT things that need to be preserved
      std::unique_ptr<Ort::Env> _env;
      std::unique_ptr<Ort::Session> _session;

      /// Useful Info about the object, may as well be cached in member vars, will be needed
      /// multiple times.
      std::vector<int64_t> _inputNodeDims;//I don't like this int64_t stuff but ORT insisted.
      std::vector<int64_t> _outputNodeDims;
      ONNXTensorElementDataType _inType;
      ONNXTensorElementDataType _outType;

      string _inputNodeName;
      string _outputNodeName;  
  };
}


#endif