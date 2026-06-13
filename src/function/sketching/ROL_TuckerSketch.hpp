// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TUCKERSKETCH_HPP
#define ROL_TUCKERSKETCH_HPP

#include "ROL_Sketch.hpp"
#include "ROL_Vector.hpp"
#include "ROL_StdVector.hpp"
#include "/projects/TuckerMPI/src/serial/Tucker.hpp"

namespace ROL {

/** @ingroup func_group
    \class TuckerSketch
    \brief Provides a sketching interface using Tucker decomposition.
*/
template <class Real>
class TuckerSketch : public Sketch<Real> {
private:
  // TuckerMPI storage
  Tucker::Tensor<Real>* historyTensor_;
  int stateDim_;
  int timeDim_;

public:
  TuckerSketch(const Vector<Real> &x, int ncol, int rank,
               Real orthTol = 1e-8, int orthIt = 2, bool truncate = false,
               unsigned dom_seed = 0, unsigned rng_seed = 0)
    : Sketch<Real>(x, ncol, rank, orthTol, orthIt, truncate, dom_seed, rng_seed),
      historyTensor_(nullPtr) {
    stateDim_ = x.dimension();
    timeDim_ = ncol;

    // Initialize the history tensor.
    // Mode 0: Time, Mode 1: State
    historyTensor_ = new Tucker::Tensor<Real>(timeDim_, stateDim_);
    this->reset(true);
  }

  ~TuckerSketch() {
    if (historyTensor_ != nullPtr) {
      delete historyTensor_;
      historyTensor_ = nullPtr;
    }
  }

  void reset(bool randomize = true) override {
    // In the pseudo-online approach, we just zero out the tensor.
    // 'randomize' is ignored as Tucker decomposition is deterministic
    // based on the data provided.
    for (int i = 0; i < timeDim_; ++i) {
      for (int j = 0; j < stateDim_; ++j) {
        (*historyTensor_)(i, j) = static_cast<Real>(0);
      }
    }
  }

  int advance(Real nu, const Vector<Real> &h, int col, Real eta = 1.0) override {
    if (col >= timeDim_ || col < 0) return 1;

    const Real scale = nu * eta;

    // To be truly general for any ROL::Vector (e.g. PDE Opt vectors wrapping Tpetra),
    // we extract values using the basis vectors.
    // While O(N*dim), this is the only way to avoid assuming a concrete subclass.
    for (int j = 0; j < stateDim_; ++j) {
      Ptr<Vector<Real>> bj = h.basis(j);
      if (bj == nullPtr) {
        // If basis is not implemented, this Vector type is not supported.
        return 1;
      }
      (*historyTensor_)(col, j) = h.dot(*bj) * scale;
    }
    return 0;
  }

  int reconstruct(Vector<Real> &a, const int col) override {
    if (col >= timeDim_ || col < 0) return 2;

    try {
      Tucker::SizeArray reducedI(2);
      reducedI[0] = this->rank_;
      reducedI[1] = this->rank_;

      const Tucker::TuckerTensor<Real>* tt = Tucker::STHOSVD(historyTensor_, &reducedI);
      if (tt == nullPtr) return 5;

      // Reconstruct the vector at slice 'col'.
      // a = G(col, :) * U_1
      const Tucker::Matrix<Real>* U1 = tt->U[1];
      const Tucker::Tensor<Real>* G = tt->G;

      // Use a temporary StdVector to build the result.
      Ptr<StdVector<Real>> tmp = makePtr<StdVector<Real>>(stateDim_);
      for (int j = 0; j < stateDim_; ++j) {
        Real sum = 0;
        for (int k = 0; k < this->rank_; ++k) {
          sum += (*G)(col, k) * (*U1)(j, k);
        }
        (*tmp)[j] = sum;
      }

      // Copy the reconstructed result back into the general ROL::Vector.
      a.set(*tmp);

      delete tt;
    } catch (const std::exception& e) {
      return 5;
    }

    return 0;
  }

  void setRank(int rank) override {
    Sketch<Real>::setRank(rank);
  }
};

} // namespace ROL

#endif // ROL_TUCKERSKETCH_HPP
