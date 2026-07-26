#include "spinsim/ensemble/EnsembleRunner.hpp"

#include <atomic>
#include <exception>
#include <format>
#include <mutex>
#include <stdexcept>
#include <thread>

namespace spinsim {

void VectorRecorder::record(double time, std::span<const double> values) {
  times_.push_back(time);
  values_.insert(values_.end(), values.begin(), values.end());
}

void VectorRecorder::clear() {
  times_.clear();
  values_.clear();
}

LocalRunner::LocalRunner(std::size_t threads)
    : threads_(threads > 0 ? threads : std::thread::hardware_concurrency()) {
  if (threads_ == 0) threads_ = 1;
}

EnsembleResult LocalRunner::run(std::size_t realizations,
                                std::span<const std::string> columns,
                                const RealizationJob& job) {
  if (realizations == 0) {
    throw std::invalid_argument("an ensemble needs at least one realisation");
  }

  // Every realisation's output is kept, then folded in index order below. The
  // storage is realisations * points * columns doubles -- a few hundred
  // kilobytes for the runs this code was written for, and the price of an
  // answer that does not depend on the thread count.
  std::vector<VectorRecorder> perRealization(realizations);

  std::atomic<std::size_t> next{0};
  std::mutex failureMutex;
  std::exception_ptr failure;

  const std::size_t workers = std::min(threads_, realizations);
  {
    std::vector<std::jthread> pool;
    pool.reserve(workers);
    for (std::size_t w = 0; w < workers; ++w) {
      pool.emplace_back([&] {
        while (true) {
          const std::size_t index = next.fetch_add(1);
          if (index >= realizations) return;
          try {
            job(index, perRealization[index]);
          } catch (...) {
            const std::lock_guard lock(failureMutex);
            if (!failure) failure = std::current_exception();
            // Stop handing out work; the run is already lost.
            next.store(realizations);
            return;
          }
        }
      });
    }
  }  // jthreads join here

  if (failure) std::rethrow_exception(failure);

  EnsembleResult result;
  result.times = perRealization.front().times();
  result.columns.assign(columns.begin(), columns.end());

  const std::size_t points = result.times.size();
  const std::size_t width = points * result.columns.size();
  result.statistics = Accumulator(width);

  for (std::size_t i = 0; i < realizations; ++i) {
    const VectorRecorder& recorded = perRealization[i];
    if (recorded.rows() != points || recorded.values().size() != width) {
      throw std::runtime_error(std::format(
          "realisation {} recorded {} points and {} values; expected {} and {}",
          i, recorded.rows(), recorded.values().size(), points, width));
    }
    result.statistics.add(recorded.values());
  }

  return result;
}

}  // namespace spinsim
