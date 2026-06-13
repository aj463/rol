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

#include "Tucker.hpp"

#if defined(__has_include)
#  if __has_include("ROL_TpetraMultiVector.hpp") && __has_include("Tpetra_MultiVector.hpp") && __has_include("Teuchos_CommHelpers.hpp")
#    include "ROL_TpetraMultiVector.hpp"
#    include "Teuchos_CommHelpers.hpp"
#    define ROL_TUCKERSKETCH_HAS_TPETRA 1
#  endif
#endif
#ifndef ROL_TUCKERSKETCH_HAS_TPETRA
#  define ROL_TUCKERSKETCH_HAS_TPETRA 0
#endif

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace ROL {

/** @ingroup func_group
    \class TuckerSketch
    \brief Provides a sketching interface using Tucker decomposition.
*/
template <class Real>
class TuckerSketch : public Sketch<Real> {
private:
  // TuckerMPI storage
  std::vector<Real> history_;
  Tucker::Tensor<Real>* historyTensor_;
  const Tucker::TuckerTensor<Real>* factorization_;
  std::vector<Ptr<Vector<Real>>> basis_;
  int stateDim_;
  int timeDim_;
  int maxRank_;
  int stateRank_;
  int timeRank_;
  bool tpetraFastPath_;

  // Preallocated workspaces to prevent heap allocations in hot loops
  mutable std::vector<Real> coeffWorkspace_;
  mutable std::vector<Real> tempWorkspace_;

#if ROL_TUCKERSKETCH_HAS_TPETRA
  std::vector<Real> localColumnWorkspace_;
  std::vector<Real> globalColumnWorkspace_;
#endif

  int index(const int state, const int time) const {
    return state + time * stateDim_;
  }

  Real matrixEntry(const Tucker::Matrix<Real>* mat,
                   const int row, const int col) const {
    return mat->data()[row + col * mat->nrows()];
  }

  Real coreEntry(const Tucker::Tensor<Real>* core,
                 const int row, const int col) const {
    return core->data()[row + col * core->size(0)];
  }

  void setTuckerRank(const int rank) {
    const int requestedRank = std::max(rank, 1);
    const int rankBound = std::max(maxRank_, 1);
    const int effectiveRank = std::min(requestedRank, rankBound);
    stateRank_ = std::min(effectiveRank, stateDim_);
    timeRank_  = std::min(effectiveRank, timeDim_);
  }

  void clearFactorization() {
    if (factorization_ != nullptr) {
      Tucker::MemoryManager::safe_delete(factorization_);
      factorization_ = nullptr;
    }
  }

  void allocateTensor() {
    if (historyTensor_ != nullptr) return;

    Tucker::SizeArray tensorSize(2);
    tensorSize[0] = stateDim_;
    tensorSize[1] = timeDim_;
    historyTensor_ = Tucker::MemoryManager::safe_new<Tucker::Tensor<Real>>(tensorSize);
  }

  void copyHistoryToTensor() {
    allocateTensor();
    std::copy(history_.begin(), history_.end(), historyTensor_->data());
  }

  int computeFactorization() {
    if (factorization_ != nullptr) return 0;
    if (stateRank_ <= 0 || timeRank_ <= 0) return 5;

    copyHistoryToTensor();
    Tucker::SizeArray reducedI(2);
    reducedI[0] = stateRank_;
    reducedI[1] = timeRank_;
    factorization_ = Tucker::STHOSVD(historyTensor_, &reducedI);
    return (factorization_ == nullptr ? 5 : 0);
  }

  void buildBasis(const Vector<Real> &x) {
    basis_.resize(stateDim_);
    for (int i = 0; i < stateDim_; ++i) {
      basis_[i] = x.basis(i);
      ROL_TEST_FOR_EXCEPTION(basis_[i] == nullPtr, std::invalid_argument,
        "ROL::TuckerSketch requires Vector::basis to extract and rebuild coordinates.");
    }
  }

  bool hasTpetraFastPath(const Vector<Real> &x) const {
#if ROL_TUCKERSKETCH_HAS_TPETRA
    return (dynamic_cast<const TpetraMultiVector<Real>*>(&x) != nullptr);
#else
    ROL_UNUSED(x);
    return false;
#endif
  }

#if ROL_TUCKERSKETCH_HAS_TPETRA
  bool copyFromTpetraMultiVector(Real nu, const Vector<Real> &h, const int col) {
    const TpetraMultiVector<Real>* hv = dynamic_cast<const TpetraMultiVector<Real>*>(&h);
    if (hv == nullptr) return false;

    const Ptr<const Tpetra::MultiVector<Real>> mv = hv->getVector();
    const int globalLength = static_cast<int>(mv->getGlobalLength());
    const int numVectors   = static_cast<int>(mv->getNumVectors());
    if (stateDim_ != globalLength * numVectors) return false;

    std::fill(localColumnWorkspace_.begin(), localColumnWorkspace_.end(), static_cast<Real>(0));
    std::fill(globalColumnWorkspace_.begin(), globalColumnWorkspace_.end(), static_cast<Real>(0));

    const auto view = mv->getLocalViewHost(Tpetra::Access::ReadOnly);
    const auto map  = mv->getMap();
    const int localLength = static_cast<int>(mv->getLocalLength());

    for (int j = 0; j < numVectors; ++j) {
      for (int i = 0; i < localLength; ++i) {
        const int gid = static_cast<int>(map->getGlobalElement(i) - map->getIndexBase());
        if (gid < 0 || gid >= globalLength) return false;
        localColumnWorkspace_[gid + j * globalLength] = nu * view(i,j);
      }
    }

    // WARNING: This assumes the global state fits easily within local memory!
    Teuchos::reduceAll<int,Real>(*map->getComm(), Teuchos::REDUCE_SUM, stateDim_,
                                 &localColumnWorkspace_[0], &globalColumnWorkspace_[0]);

    for (int i = 0; i < stateDim_; ++i) {
      history_[index(i, col)] += globalColumnWorkspace_[i];
    }
    return true;
  }

  bool copyToTpetraMultiVector(Vector<Real> &a,
                               const std::vector<Real> &coeff) const {
    TpetraMultiVector<Real>* av = dynamic_cast<TpetraMultiVector<Real>*>(&a);
    if (av == nullptr) return false;

    const Ptr<Tpetra::MultiVector<Real>> mv = av->getVector();
    const int globalLength = static_cast<int>(mv->getGlobalLength());
    const int numVectors   = static_cast<int>(mv->getNumVectors());
    if (stateDim_ != globalLength * numVectors) return false;

    auto view = mv->getLocalViewHost(Tpetra::Access::ReadWrite);
    const auto map  = mv->getMap();
    const int localLength = static_cast<int>(mv->getLocalLength());

    for (int j = 0; j < numVectors; ++j) {
      for (int i = 0; i < localLength; ++i) {
        const int gid = static_cast<int>(map->getGlobalElement(i) - map->getIndexBase());
        if (gid < 0 || gid >= globalLength) return false;
        view(i,j) = coeff[gid + j * globalLength];
      }
    }
    return true;
  }
#endif

  void computeReconstructionColumn(std::vector<Real> &coeff,
                                   const int col) const {
    const Tucker::Matrix<Real>* U0 = factorization_->U[0];
    const Tucker::Matrix<Real>* U1 = factorization_->U[1];
    const Tucker::Tensor<Real>* G  = factorization_->G;
    ROL_TEST_FOR_EXCEPTION(U0 == nullptr || U1 == nullptr || G == nullptr,
      std::logic_error, "ROL::TuckerSketch has an invalid Tucker factorization.");

    std::fill(coeff.begin(), coeff.end(), static_cast<Real>(0));
    std::fill(tempWorkspace_.begin(), tempWorkspace_.begin() + stateRank_, static_cast<Real>(0));

    // Cache-friendly tensor contraction (G is column-major: row varies fastest)
    for (int q = 0; q < timeRank_; ++q) {
      const Real u1_val = matrixEntry(U1, col, q);
      for (int p = 0; p < stateRank_; ++p) { // Inner loop traverses G's rows sequentially
        tempWorkspace_[p] += coreEntry(G, p, q) * u1_val;
      }
    }

    // Cache-friendly outer contraction (U0 is column-major: row varies fastest)
    for (int p = 0; p < stateRank_; ++p) {
      const Real temp_val = tempWorkspace_[p];
      for (int i = 0; i < stateDim_; ++i) { // Inner loop traverses U0's rows sequentially
        coeff[i] += matrixEntry(U0, i, p) * temp_val;
      }
    }
  }

  void reconstructWithBasis(Vector<Real> &a,
                            const std::vector<Real> &coeff) const {
    ROL_TEST_FOR_EXCEPTION(basis_.size() != static_cast<size_t>(stateDim_),
      std::logic_error,
      "ROL::TuckerSketch does not have a basis representation for reconstruction.");

    a.zero();
    for (int i = 0; i < stateDim_; ++i) {
      if (coeff[i] != static_cast<Real>(0)) {
        a.axpy(coeff[i], *basis_[i]);
      }
    }
  }

public:
  // Disallow copy/move to prevent double free of Tucker-owned pointers.
  TuckerSketch(const TuckerSketch&) = delete;
  TuckerSketch& operator=(const TuckerSketch&) = delete;
  TuckerSketch(TuckerSketch&&) = delete;
  TuckerSketch& operator=(TuckerSketch&&) = delete;

  TuckerSketch(const Vector<Real> &x, int ncol, int rank,
               Real orthTol = 1e-8, int orthIt = 2, bool truncate = false,
               unsigned dom_seed = 0, unsigned rng_seed = 0)
    : Sketch<Real>(x, ncol, rank, orthTol, orthIt, truncate, dom_seed, rng_seed),
      history_(),
      historyTensor_(nullptr),
      factorization_(nullptr),
      stateDim_(x.dimension()),
      timeDim_(ncol),
      maxRank_(std::min(stateDim_, timeDim_)),
      stateRank_(0),
      timeRank_(0),
      tpetraFastPath_(false) {
    ROL_TEST_FOR_EXCEPTION(stateDim_ <= 0, std::invalid_argument,
      "ROL::TuckerSketch requires a vector with positive dimension.");
    ROL_TEST_FOR_EXCEPTION(timeDim_ <= 0, std::invalid_argument,
      "ROL::TuckerSketch requires a positive number of columns.");

    history_.assign(static_cast<size_t>(stateDim_) * static_cast<size_t>(timeDim_),
                    static_cast<Real>(0));

    // Allocate generic workspaces once up front
    coeffWorkspace_.assign(stateDim_, static_cast<Real>(0));
    tempWorkspace_.assign(maxRank_, static_cast<Real>(0));

    tpetraFastPath_ = hasTpetraFastPath(x);
#if ROL_TUCKERSKETCH_HAS_TPETRA
    if (tpetraFastPath_) {
      localColumnWorkspace_.assign(stateDim_, static_cast<Real>(0));
      globalColumnWorkspace_.assign(stateDim_, static_cast<Real>(0));
    }
#endif

    setTuckerRank(rank);
    if (!tpetraFastPath_) {
      buildBasis(x);
    }
    this->reset(true);
  }

  ~TuckerSketch() override {
    clearFactorization();
    if (historyTensor_ != nullptr) {
      Tucker::MemoryManager::safe_delete(historyTensor_);
    }
  }

  void reset(bool randomize = true) override {
    ROL_UNUSED(randomize);
    clearFactorization();
    std::fill(history_.begin(), history_.end(), static_cast<Real>(0));
  }

  int advance(Real nu, const Vector<Real> &h, int col, Real eta = 1.0) override {
    if (col >= timeDim_ || col < 0) return 1;
    if (h.dimension() != stateDim_) return 1;

    clearFactorization();
    if (eta != static_cast<Real>(1)) {
      for (auto& val : history_) {
        val *= eta;
      }
    }

    if (tpetraFastPath_) {
#if ROL_TUCKERSKETCH_HAS_TPETRA
      return (copyFromTpetraMultiVector(nu, h, col) ? 0 : 1);
#else
      return 1;
#endif
    }

    for (int j = 0; j < stateDim_; ++j) {
      history_[index(j, col)] += nu * h.dot(*basis_[j]);
    }
    return 0;
  }

  int reconstruct(Vector<Real> &a, const int col) override {
    if (col >= timeDim_ || col < 0) return 2;
    if (a.dimension() != stateDim_) return 5;

    try {
      const int info = computeFactorization();
      if (info != 0) return info;

      computeReconstructionColumn(coeffWorkspace_, col);

      if (tpetraFastPath_) {
#if ROL_TUCKERSKETCH_HAS_TPETRA
        return (copyToTpetraMultiVector(a, coeffWorkspace_) ? 0 : 5);
#else
        return 5;
#endif
      }
      reconstructWithBasis(a, coeffWorkspace_);
    } catch (const std::exception&) {
      return 5;
    }

    return 0;
  }

  void setRank(int rank) override {
    setTuckerRank(rank);
    reset(true);
  }

  void update(void) override {
    reset(true);
  }
};

} // namespace ROL

#endif // ROL_TUCKERSKETCH_HPP
