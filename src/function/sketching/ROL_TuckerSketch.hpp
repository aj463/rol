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
#include "Tucker_StreamingTuckerTensor.hpp"

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
#include <cmath>
#include <cstddef>
#include <memory>
#include <limits>
#include <stdexcept>
#include <vector>

namespace ROL {

/**
    \class TuckerSketch
    \brief Provides a sketching interface using Tucker decomposition.
*/
template <class Real>
class TuckerSketch : public Sketch<Real> {
private:
  enum StatusCode {
    STATUS_SUCCESS              = 0,
    STATUS_ADVANCE_INPUT_ERROR  = 1,
    STATUS_RECONSTRUCT_ERROR    = 2,
    STATUS_FACTORIZATION_ERROR  = 5
  };

  // Generic deleter for all Tucker objects allocated via Tucker::MemoryManager
  struct TuckerDeleter {
    template <typename T>
    void operator()(T* ptr) const {
      if (ptr != nullptr) {
        Tucker::MemoryManager::safe_delete(ptr);
      }
    }
  };

  using TensorPtr          = std::unique_ptr<Tucker::Tensor<Real>, TuckerDeleter>;
  using StreamingTensorPtr = std::unique_ptr<Tucker::StreamingTuckerTensor<Real>, TuckerDeleter>;
  using TuckerVectorPtr    = std::unique_ptr<Tucker::Vector<Real>, TuckerDeleter>;

  // Variables ordered precisely for safe initialization in the constructor
  int stateDim_;
  int timeDim_;
  int stateRank_;
  int timeRank_;
  int streamedColumns_;
  int streamDirection_;
  Real epsilon_;
  bool tpetraFastPath_;

  // TuckerMPI storage utilizing RAII
  StreamingTensorPtr streamingFactorization_;
  TuckerVectorPtr epsilonList_;

  // Generic coordinate basis used when no direct vector-storage path is available.
  std::vector<Ptr<Vector<Real>>> basis_;

  // Reused workspaces for reconstruction.
  std::vector<Real> coeffWorkspace_;
  std::vector<Real> tempWorkspace_;

#if ROL_TUCKERSKETCH_HAS_TPETRA
  // Chunk workspace used while reducing distributed Tpetra data into the dense Tucker slice.
  std::vector<Real> tpetraReduceWorkspace_;
#endif

  static Real matrixEntry(const Tucker::Matrix<Real>* mat, const int row, const int col) {
    return mat->data()[row + col * mat->nrows()];
  }

  static Real coreEntry(const Tucker::Tensor<Real>* core, const int row, const int col) {
    return core->data()[row + col * core->size(0)];
  }

  static bool hasTpetraFastPath(const Vector<Real>& x) {
#if ROL_TUCKERSKETCH_HAS_TPETRA
    return (dynamic_cast<const TpetraMultiVector<Real>*>(&x) != nullptr);
#else
    ROL_UNUSED(x);
    return false;
#endif
  }

  bool validColumn(const int col) const {
    return (0 <= col && col < timeDim_);
  }

  int expectedInputColumn() const {
    return (streamDirection_ < 0 ? timeDim_ - 1 - streamedColumns_
                                 : streamedColumns_);
  }

  int streamedColumn(const int col) const {
    return (streamDirection_ < 0 ? timeDim_ - 1 - col : col);
  }

  bool updateStreamDirection(const int col) {
    if (streamedColumns_ == 0) {
      if (col == 0) {
        streamDirection_ = 1;
      }
      else if (col == timeDim_ - 1) {
        streamDirection_ = -1;
      }
      else {
        return false;
      }
    }
    return (col == expectedInputColumn());
  }

  size_t tpetraReduceChunkSize() const {
    // Keep the extra communication buffer bounded. The dense Tucker slice is
    // still required by the TuckerMPI streaming update, but the gather no longer
    // allocates another full global receive vector on each rank.
    static constexpr size_t maxEntries = static_cast<size_t>(1) << 20;
    return std::min(static_cast<size_t>(stateDim_), maxEntries);
  }

  void clearFactorization() {
    streamingFactorization_.reset();
    streamedColumns_ = 0;
    streamDirection_ = 0;
    stateRank_ = 0;
    timeRank_  = 0;
  }

  void updateEpsilonList() {
    if (!epsilonList_) {
      epsilonList_.reset(Tucker::MemoryManager::safe_new<Tucker::Vector<Real>>(2));
    }
    const Real epsilonPerMode = epsilon_ / std::sqrt(static_cast<Real>(2));
    (*epsilonList_)[0] = epsilonPerMode;
    (*epsilonList_)[1] = epsilonPerMode;
  }

  TensorPtr createSliceTensor() const {
    Tucker::SizeArray sliceSize(2);
    sliceSize[0] = stateDim_;
    sliceSize[1] = 1;
    return TensorPtr(Tucker::MemoryManager::safe_new<Tucker::Tensor<Real>>(sliceSize));
  }

  const Tucker::TuckerTensor<Real>* factorization() const {
    return (streamingFactorization_ == nullptr ? nullptr : streamingFactorization_->factorization);
  }

  int updateRanksFromFactorization() {
    const Tucker::TuckerTensor<Real>* F = factorization();
    if (F == nullptr || F->U[0] == nullptr || F->U[1] == nullptr) {
      return STATUS_FACTORIZATION_ERROR;
    }

    stateRank_ = F->U[0]->ncols();
    timeRank_  = F->U[1]->ncols();
    if (stateRank_ <= 0 || timeRank_ <= 0) {
      return STATUS_FACTORIZATION_ERROR;
    }

    if (tempWorkspace_.size() < static_cast<size_t>(stateRank_)) {
      tempWorkspace_.resize(stateRank_);
    }
    return STATUS_SUCCESS;
  }

  int initializeStreamingFactorization(const Tucker::TuckerTensor<Real>* initial,
                                       const Tucker::Tensor<Real>& slice) {
    if (initial == nullptr) {
      return STATUS_FACTORIZATION_ERROR;
    }

    streamingFactorization_.reset(
      Tucker::MemoryManager::safe_new<Tucker::StreamingTuckerTensor<Real>>(initial));

    for (int i = 0; i < streamingFactorization_->N - 1; ++i) {
      streamingFactorization_->Gram[i] = nullptr;
    }

    streamingFactorization_->isvd->initializeFactors(streamingFactorization_->factorization);
    streamingFactorization_->Xnorm2 = slice.norm2();

    const int ndims = streamingFactorization_->N;
    const Tucker::TuckerTensor<Real>* F = streamingFactorization_->factorization;

    for (int n = 0; n < ndims - 1; ++n) {
      streamingFactorization_->squared_errors[n] = static_cast<Real>(0);

      for (int i = F->G->size(n); i < F->U[n]->nrows(); ++i) {
        if (F->eigenvalues != nullptr && F->eigenvalues[n] != nullptr) {
          streamingFactorization_->squared_errors[n] += std::abs(F->eigenvalues[n][i]);
        }
        else if (F->singularValues != nullptr && F->singularValues[n] != nullptr) {
          streamingFactorization_->squared_errors[n] += F->singularValues[n][i] * F->singularValues[n][i];
        }
      }
    }

    const Real isvdError = streamingFactorization_->isvd->getErrorNorm();
    streamingFactorization_->squared_errors[ndims - 1] = isvdError * isvdError;

    return updateRanksFromFactorization();
  }

  void buildBasis(const Vector<Real>& x) {
    basis_.resize(stateDim_);
    for (int i = 0; i < stateDim_; ++i) {
      basis_[i] = x.basis(i);
      ROL_TEST_FOR_EXCEPTION(basis_[i] == nullPtr, std::invalid_argument,
        "ROL::TuckerSketch requires Vector::basis to extract and rebuild coordinates.");
    }
  }

  void computeTemporalCoefficients(const int col) {
    const Tucker::TuckerTensor<Real>* F = factorization();
    ROL_TEST_FOR_EXCEPTION(F == nullptr,
      std::logic_error, "ROL::TuckerSketch has no Tucker factorization.");

    const Tucker::Matrix<Real>* U1 = F->U[1];
    const Tucker::Tensor<Real>* G  = F->G;
    ROL_TEST_FOR_EXCEPTION(U1 == nullptr || G == nullptr,
      std::logic_error, "ROL::TuckerSketch has an invalid Tucker factorization.");

    std::fill(tempWorkspace_.begin(), tempWorkspace_.begin() + stateRank_, static_cast<Real>(0));

    for (int q = 0; q < timeRank_; ++q) {
      const Real u1Value = matrixEntry(U1, col, q);
      for (int p = 0; p < stateRank_; ++p) {
        tempWorkspace_[p] += coreEntry(G, p, q) * u1Value;
      }
    }
  }

#if ROL_TUCKERSKETCH_HAS_TPETRA
  bool tpetraLayoutIsCompatible(const Tpetra::MultiVector<Real>& mv,
                                int& globalLength,
                                int& numVectors) const {
    const auto rawGlobalLength = mv.getGlobalLength();
    const auto rawNumVectors   = mv.getNumVectors();

    if (rawGlobalLength > static_cast<decltype(rawGlobalLength)>(std::numeric_limits<int>::max()) ||
        rawNumVectors   > static_cast<decltype(rawNumVectors)>(std::numeric_limits<int>::max())) {
      return false;
    }

    globalLength = static_cast<int>(rawGlobalLength);
    numVectors   = static_cast<int>(rawNumVectors);

    const size_t expectedDim =
      static_cast<size_t>(globalLength) * static_cast<size_t>(numVectors);
    return (expectedDim == static_cast<size_t>(stateDim_));
  }

  bool scatterLocalTpetraEntriesIntoSlice(Tucker::Tensor<Real>& slice,
                                          const Real nu,
                                          const Tpetra::MultiVector<Real>& mv,
                                          const int globalLength,
                                          const int numVectors) const {
    std::fill(slice.data(), slice.data() + stateDim_, static_cast<Real>(0));

    const auto view = mv.getLocalViewHost(Tpetra::Access::ReadOnly);
    const auto map  = mv.getMap();
    const int localLength = static_cast<int>(mv.getLocalLength());

    for (int j = 0; j < numVectors; ++j) {
      for (int i = 0; i < localLength; ++i) {
        const int gid = static_cast<int>(map->getGlobalElement(i) - map->getIndexBase());
        if (gid < 0 || gid >= globalLength) {
          return false;
        }
        slice.data()[gid + j * globalLength] = nu * view(i, j);
      }
    }
    return true;
  }

  bool reduceTpetraSliceInChunks(Tucker::Tensor<Real>& slice,
                                 const Teuchos::Comm<int>& comm) {
    if (tpetraReduceWorkspace_.empty()) {
      tpetraReduceWorkspace_.assign(tpetraReduceChunkSize(), static_cast<Real>(0));
    }

    const size_t chunkSize = tpetraReduceWorkspace_.size();
    for (size_t offset = 0; offset < static_cast<size_t>(stateDim_); offset += chunkSize) {
      const size_t count = std::min(chunkSize, static_cast<size_t>(stateDim_) - offset);

      Teuchos::reduceAll<int, Real>(comm, Teuchos::REDUCE_SUM,
                                    static_cast<int>(count),
                                    slice.data() + offset,
                                    tpetraReduceWorkspace_.data());

      std::copy(tpetraReduceWorkspace_.begin(),
                tpetraReduceWorkspace_.begin() + count,
                slice.data() + offset);
    }
    return true;
  }

  bool copyFromTpetraMultiVector(Tucker::Tensor<Real>& slice,
                                 const Real nu,
                                 const Vector<Real>& h) {
    const TpetraMultiVector<Real>* hv = dynamic_cast<const TpetraMultiVector<Real>*>(&h);
    if (hv == nullptr) {
      return false;
    }

    const Ptr<const Tpetra::MultiVector<Real>> mv = hv->getVector();
    int globalLength = 0;
    int numVectors = 0;
    if (!tpetraLayoutIsCompatible(*mv, globalLength, numVectors)) {
      return false;
    }

    if (!scatterLocalTpetraEntriesIntoSlice(slice, nu, *mv, globalLength, numVectors)) {
      return false;
    }

    return reduceTpetraSliceInChunks(slice, *mv->getMap()->getComm());
  }

  bool reconstructTpetraMultiVector(Vector<Real>& a, const int col) {
    TpetraMultiVector<Real>* av = dynamic_cast<TpetraMultiVector<Real>*>(&a);
    if (av == nullptr) {
      return false;
    }

    const Ptr<Tpetra::MultiVector<Real>> mv = av->getVector();
    int globalLength = 0;
    int numVectors = 0;
    if (!tpetraLayoutIsCompatible(*mv, globalLength, numVectors)) {
      return false;
    }

    const Tucker::TuckerTensor<Real>* F = factorization();
    ROL_TEST_FOR_EXCEPTION(F == nullptr || F->U[0] == nullptr,
      std::logic_error, "ROL::TuckerSketch has an invalid Tucker factorization.");

    computeTemporalCoefficients(col);

    auto view = mv->getLocalViewHost(Tpetra::Access::ReadWrite);
    const auto map = mv->getMap();
    const int localLength = static_cast<int>(mv->getLocalLength());

    for (int j = 0; j < numVectors; ++j) {
      for (int i = 0; i < localLength; ++i) {
        const int gid = static_cast<int>(map->getGlobalElement(i) - map->getIndexBase());
        if (gid < 0 || gid >= globalLength) {
          return false;
        }

        const int row = gid + j * globalLength;
        Real value = static_cast<Real>(0);
        for (int p = 0; p < stateRank_; ++p) {
          value += matrixEntry(F->U[0], row, p) * tempWorkspace_[p];
        }
        view(i, j) = value;
      }
    }
    return true;
  }
#endif

  bool fillSliceTensor(Tucker::Tensor<Real>& slice,
                       const Real nu,
                       const Vector<Real>& h) const {
    for (int j = 0; j < stateDim_; ++j) {
      slice.data()[j] = nu * h.dot(*basis_[j]);
    }
    return true;
  }

  void computeReconstructionColumn(std::vector<Real>& coeff, const int col) {
    const Tucker::TuckerTensor<Real>* F = factorization();
    ROL_TEST_FOR_EXCEPTION(F == nullptr,
      std::logic_error, "ROL::TuckerSketch has no Tucker factorization.");

    const Tucker::Matrix<Real>* U0 = F->U[0];
    ROL_TEST_FOR_EXCEPTION(U0 == nullptr,
      std::logic_error, "ROL::TuckerSketch has an invalid Tucker factorization.");

    std::fill(coeff.begin(), coeff.end(), static_cast<Real>(0));
    computeTemporalCoefficients(col);

    // Compute coeff = U0 * temp.
    for (int p = 0; p < stateRank_; ++p) {
      const Real tempValue = tempWorkspace_[p];
      for (int i = 0; i < stateDim_; ++i) {
        coeff[i] += matrixEntry(U0, i, p) * tempValue;
      }
    }
  }

  void reconstructWithBasis(Vector<Real>& a, const std::vector<Real>& coeff) const {
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
  TuckerSketch(const TuckerSketch&) = delete;
  TuckerSketch& operator=(const TuckerSketch&) = delete;
  TuckerSketch(TuckerSketch&&) = delete;
  TuckerSketch& operator=(TuckerSketch&&) = delete;

  TuckerSketch(const Vector<Real>& x,
               int ncol,
               int rank,
               Real orthTol = 1e-8,
               int orthIt = 2,
               bool truncate = false,
               unsigned dom_seed = 0,
               unsigned rng_seed = 0)
    : TuckerSketch(x, ncol, rank, orthTol, orthIt, truncate,
                   dom_seed, rng_seed, static_cast<Real>(1e-6)) {}

  TuckerSketch(const Vector<Real>& x,
               int ncol,
               int rank,
               Real orthTol,
               int orthIt,
               bool truncate,
               unsigned dom_seed,
               unsigned rng_seed,
               Real epsilon)
    : Sketch<Real>(x, ncol, rank, orthTol, orthIt, truncate, dom_seed, rng_seed),
      stateDim_(x.dimension()),
      timeDim_(ncol),
      stateRank_(0),
      timeRank_(0),
      streamedColumns_(0),
      streamDirection_(0),
      epsilon_(epsilon),
      tpetraFastPath_(hasTpetraFastPath(x)),
      streamingFactorization_(nullptr),
      epsilonList_(nullptr) {

    ROL_TEST_FOR_EXCEPTION(stateDim_ <= 0, std::invalid_argument,
      "ROL::TuckerSketch requires a vector with positive dimension.");
    ROL_TEST_FOR_EXCEPTION(timeDim_ <= 0, std::invalid_argument,
      "ROL::TuckerSketch requires a positive number of columns.");
    ROL_TEST_FOR_EXCEPTION(epsilon_ < static_cast<Real>(0), std::invalid_argument,
      "ROL::TuckerSketch requires a nonnegative STHOSVD tolerance.");

#if ROL_TUCKERSKETCH_HAS_TPETRA
    if (tpetraFastPath_) {
      tpetraReduceWorkspace_.assign(tpetraReduceChunkSize(), static_cast<Real>(0));
    }
#endif

    ROL_UNUSED(rank);
    updateEpsilonList();

    if (!tpetraFastPath_) {
      coeffWorkspace_.assign(stateDim_, static_cast<Real>(0));
      buildBasis(x);
    }
  }

  ~TuckerSketch() override = default;

  void reset(bool randomize = true) override {
    ROL_UNUSED(randomize);
    clearFactorization();
  }

  int advance(Real nu, const Vector<Real>& h, int col, Real eta = static_cast<Real>(1)) override {
    if (!validColumn(col) || h.dimension() != stateDim_) {
      return STATUS_ADVANCE_INPUT_ERROR;
    }
    if (eta != static_cast<Real>(1) || !updateStreamDirection(col)) {
      return STATUS_ADVANCE_INPUT_ERROR;
    }
    if (epsilon_ < static_cast<Real>(0)) {
      return STATUS_FACTORIZATION_ERROR;
    }

    try {
      updateEpsilonList();
      TensorPtr slice = createSliceTensor();

      bool filled = false;
      if (tpetraFastPath_) {
#if ROL_TUCKERSKETCH_HAS_TPETRA
        filled = copyFromTpetraMultiVector(*slice, nu, h);
#endif
      }
      else {
        filled = fillSliceTensor(*slice, nu, h);
      }

      if (!filled) {
        return STATUS_ADVANCE_INPUT_ERROR;
      }

      if (streamedColumns_ == 0) {
        // NOTE: This follows the existing TuckerMPI ownership convention used by
        // StreamingTuckerTensor: the constructor assumes responsibility for the
        // factorization returned by STHOSVD.
        const Tucker::TuckerTensor<Real>* initial = Tucker::STHOSVD(slice.get(), epsilon_);
        const int info = initializeStreamingFactorization(initial, *slice);
        if (info != STATUS_SUCCESS) {
          clearFactorization();
          return info;
        }
      }
      else {
        if (!streamingFactorization_) {
          return STATUS_FACTORIZATION_ERROR;
        }

        Tucker::Tensor<Real>* streamingSlice = slice.release();
        Tucker::StreamingSTHOSVDUpdate(streamingFactorization_.get(), streamingSlice, epsilonList_.get());
        // StreamingSTHOSVDUpdate consumes streamingSlice as workspace.

        if (updateRanksFromFactorization() != STATUS_SUCCESS) {
          clearFactorization();
          return STATUS_FACTORIZATION_ERROR;
        }
      }

      ++streamedColumns_;
    }
    catch (const std::exception&) {
      clearFactorization();
      return STATUS_FACTORIZATION_ERROR;
    }

    return STATUS_SUCCESS;
  }

  int reconstruct(Vector<Real>& a, const int col) override {
    const int scol = streamedColumn(col);
    if (!validColumn(col) || scol >= streamedColumns_) {
      return STATUS_RECONSTRUCT_ERROR;
    }
    if (a.dimension() != stateDim_) {
      return STATUS_FACTORIZATION_ERROR;
    }

    try {
      if (updateRanksFromFactorization() != STATUS_SUCCESS) {
        return STATUS_FACTORIZATION_ERROR;
      }

      if (tpetraFastPath_) {
#if ROL_TUCKERSKETCH_HAS_TPETRA
        return reconstructTpetraMultiVector(a, scol) ? STATUS_SUCCESS : STATUS_FACTORIZATION_ERROR;
#else
        return STATUS_FACTORIZATION_ERROR;
#endif
      }

      computeReconstructionColumn(coeffWorkspace_, scol);
      reconstructWithBasis(a, coeffWorkspace_);
    }
    catch (const std::exception&) {
      return STATUS_FACTORIZATION_ERROR;
    }

    return STATUS_SUCCESS;
  }

  void setRank(int rank) override {
    ROL_UNUSED(rank);
  }

  void setTolerance(Real epsilon) override {
    ROL_TEST_FOR_EXCEPTION(epsilon < static_cast<Real>(0), std::invalid_argument,
      "ROL::TuckerSketch requires a nonnegative STHOSVD tolerance.");
    epsilon_ = epsilon;
    updateEpsilonList();
    clearFactorization();
  }

  void scaleTolerance(Real factor) override {
    ROL_TEST_FOR_EXCEPTION(factor <= static_cast<Real>(0), std::invalid_argument,
      "ROL::TuckerSketch requires a positive tolerance scale factor.");
    epsilon_ = factor * epsilon_;
    updateEpsilonList();
    clearFactorization();
  }

  Real getTolerance(void) const override {
    return epsilon_;
  }

  void update(void) override {
    reset(true);
  }
};

} // namespace ROL

#endif // ROL_TUCKERSKETCH_HPP
