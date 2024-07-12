// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2016 Benoit Steiner <benoit.steiner.goog@gmail.com>
// Copyright (C) 2019 Rasmus Munk Larsen <rmlarsen@google.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_GENERIC_TYPE_CASTING_H
#define EIGEN_GENERIC_TYPE_CASTING_H

namespace RivetEigen {

namespace internal {

template<>
struct scalar_cast_op<float, RivetEigen::half> {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_cast_op)
  typedef RivetEigen::half result_type;
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE RivetEigen::half operator() (const float& a) const {
    #if (defined(EIGEN_HAS_CUDA_FP16) && defined(EIGEN_CUDA_ARCH) && EIGEN_CUDA_ARCH >= 300) || \
      (defined(EIGEN_HAS_HIP_FP16) && defined(EIGEN_HIP_DEVICE_COMPILE))
      return __float2half(a);
    #else
      return RivetEigen::half(a);
    #endif
  }
};

template<>
struct functor_traits<scalar_cast_op<float, RivetEigen::half> >
{ enum { Cost = NumTraits<float>::AddCost, PacketAccess = false }; };


template<>
struct scalar_cast_op<int, RivetEigen::half> {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_cast_op)
  typedef RivetEigen::half result_type;
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE RivetEigen::half operator() (const int& a) const {
    #if (defined(EIGEN_HAS_CUDA_FP16) && defined(EIGEN_CUDA_ARCH) && EIGEN_CUDA_ARCH >= 300) || \
      (defined(EIGEN_HAS_HIP_FP16) && defined(EIGEN_HIP_DEVICE_COMPILE))
      return __float2half(static_cast<float>(a));
    #else
      return RivetEigen::half(static_cast<float>(a));
    #endif
  }
};

template<>
struct functor_traits<scalar_cast_op<int, RivetEigen::half> >
{ enum { Cost = NumTraits<float>::AddCost, PacketAccess = false }; };


template<>
struct scalar_cast_op<RivetEigen::half, float> {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_cast_op)
  typedef float result_type;
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE float operator() (const RivetEigen::half& a) const {
    #if (defined(EIGEN_HAS_CUDA_FP16) && defined(EIGEN_CUDA_ARCH) && EIGEN_CUDA_ARCH >= 300) || \
      (defined(EIGEN_HAS_HIP_FP16) && defined(EIGEN_HIP_DEVICE_COMPILE))
      return __half2float(a);
    #else
      return static_cast<float>(a);
    #endif
  }
};

template<>
struct functor_traits<scalar_cast_op<RivetEigen::half, float> >
{ enum { Cost = NumTraits<float>::AddCost, PacketAccess = false }; };


template<>
struct scalar_cast_op<float, RivetEigen::bfloat16> {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_cast_op)
  typedef RivetEigen::bfloat16 result_type;
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE RivetEigen::bfloat16 operator() (const float& a) const {
    return RivetEigen::bfloat16(a);
  }
};

template<>
struct functor_traits<scalar_cast_op<float, RivetEigen::bfloat16> >
{ enum { Cost = NumTraits<float>::AddCost, PacketAccess = false }; };


template<>
struct scalar_cast_op<int, RivetEigen::bfloat16> {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_cast_op)
  typedef RivetEigen::bfloat16 result_type;
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE RivetEigen::bfloat16 operator() (const int& a) const {
    return RivetEigen::bfloat16(static_cast<float>(a));
  }
};

template<>
struct functor_traits<scalar_cast_op<int, RivetEigen::bfloat16> >
{ enum { Cost = NumTraits<float>::AddCost, PacketAccess = false }; };


template<>
struct scalar_cast_op<RivetEigen::bfloat16, float> {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_cast_op)
  typedef float result_type;
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE float operator() (const RivetEigen::bfloat16& a) const {
    return static_cast<float>(a);
  }
};

template<>
struct functor_traits<scalar_cast_op<RivetEigen::bfloat16, float> >
{ enum { Cost = NumTraits<float>::AddCost, PacketAccess = false }; };


}
}

#endif  // EIGEN_GENERIC_TYPE_CASTING_H
