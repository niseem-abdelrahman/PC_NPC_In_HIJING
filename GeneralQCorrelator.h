#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace qcorr {

using Mask = std::uint64_t;

class GeneralQCorrelator {
public:
  enum class Subevent : std::uint8_t { Any = 0, A = 1, B = 2 };

  struct Track {
    double phi{0.0};
    double eta{0.0};
    double pt{0.0};
    double acc_weight{1.0};
    double gmc_weight{1.0};
    Mask categories{0};
  };

  struct Leg {
    int harmonic{0};
    Mask required_categories{0};
    Subevent subevent{Subevent::Any};
    bool use_gmc_weight{false};
  };

  struct Result {
    std::complex<double> numerator;
    double denominator;

    Result()
      : numerator(0.0, 0.0),
        denominator(0.0)
    {}

    Result(const std::complex<double>& n, double d)
      : numerator(n),
        denominator(d)
    {}

    [[nodiscard]] bool valid(double eps = 0.0) const {
      return std::abs(denominator) > eps;
    }

    [[nodiscard]] std::complex<double> value() const {
      if (!valid()) {
        return std::complex<double>(0.0, 0.0);
      }
      return numerator / denominator;
    }

    [[nodiscard]] double real() const { return value().real(); }
    [[nodiscard]] double imag() const { return value().imag(); }
  };

  static constexpr Mask kRef = 1ULL << 0;
  static constexpr Mask kPoi1 = 1ULL << 1;
  static constexpr Mask kPoi2 = 1ULL << 2;

  explicit GeneralQCorrelator(double eta_gap = 0.0)
    : eta_gap_(eta_gap)
  {}

  void clear_event() {
    tracks_.clear();
    clear_block_cache();
    clear_phase_cache();
  }

  void reserve(std::size_t n) {
    tracks_.reserve(n);
  }

  void reserve_block_cache(std::size_t n) const {
    event_block_cache_.reserve(n);
  }

  void set_eta_gap(double eta_gap) {
    eta_gap_ = eta_gap;
    clear_block_cache();
  }

  [[nodiscard]] double eta_gap() const { return eta_gap_; }

  void add_track(const Track& track) {
    tracks_.push_back(track);
    clear_block_cache();
    clear_phase_cache();
  }

  void clear_block_cache() const {
    event_block_cache_.clear();
  }

  [[nodiscard]] Result correlation(const std::vector<Leg>& legs) const {
    if (legs.size() < 2 || legs.size() > 6) {
      throw std::invalid_argument("GeneralQCorrelator::correlation expects 2 <= k <= 6 legs.");
    }

    const auto& partitions = partitions_for(static_cast<int>(legs.size()));
    std::array<BlockCacheEntry, 1U << 6> local_cache{};

    std::complex<double> numerator(0.0, 0.0);
    double denominator = 0.0;

    for (const auto& partition : partitions) {
      std::complex<double> part_num(partition.coefficient, 0.0);
      double part_den = partition.coefficient;

      bool zero_num = false;
      bool zero_den = false;

      for (const unsigned block_mask : partition.blocks) {
        const auto sums = block_sums(legs, block_mask, local_cache);
        part_num *= sums.first;
        part_den *= sums.second;

        if (sums.first.real() == 0.0 && sums.first.imag() == 0.0) {
          zero_num = true;
        }
        if (sums.second == 0.0) {
          zero_den = true;
        }
        if (zero_num && zero_den) {
          break;
        }
      }

      numerator += part_num;
      denominator += part_den;
    }

    return Result(numerator, denominator);
  }

  [[nodiscard]] Result same_event(std::vector<Leg> legs) const {
    for (auto& leg : legs) {
      leg.subevent = Subevent::Any;
    }
    return correlation(legs);
  }

  [[nodiscard]] Result two_subevent(const std::vector<Leg>& legs) const {
    return correlation(legs);
  }

  [[nodiscard]] std::size_t size() const { return tracks_.size(); }

private:
  struct Partition {
    std::vector<unsigned> blocks;
    double coefficient;

    Partition()
      : coefficient(0.0)
    {}

    Partition(const std::vector<unsigned>& b, double c)
      : blocks(b),
        coefficient(c)
    {}
  };

  struct BlockCacheEntry {
    bool ready;
    std::complex<double> numerator;
    double denominator;

    BlockCacheEntry()
      : ready(false),
        numerator(0.0, 0.0),
        denominator(0.0)
    {}
  };

  struct EventBlockKey {
    Mask required_categories{0};
    int harmonic_sum{0};
    int block_size{0};
    int gmc_power{0};
    Subevent subevent{Subevent::Any};

    bool operator==(const EventBlockKey& other) const noexcept {
      return required_categories == other.required_categories &&
             harmonic_sum == other.harmonic_sum &&
             block_size == other.block_size &&
             gmc_power == other.gmc_power &&
             subevent == other.subevent;
    }
  };

  struct EventBlockKeyHash {
    std::size_t operator()(const EventBlockKey& key) const noexcept {
      std::size_t h = static_cast<std::size_t>(key.required_categories);
      mix(h, static_cast<std::size_t>(static_cast<std::int64_t>(key.harmonic_sum) + 1024));
      mix(h, static_cast<std::size_t>(key.block_size));
      mix(h, static_cast<std::size_t>(key.gmc_power));
      mix(h, static_cast<std::size_t>(key.subevent));
      return h;
    }

    static void mix(std::size_t& seed, std::size_t value) noexcept {
      seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    }
  };

  struct BlockSum {
    std::complex<double> numerator;
    double denominator;

    BlockSum()
      : numerator(0.0, 0.0),
        denominator(0.0)
    {}

    BlockSum(const std::complex<double>& n, double d)
      : numerator(n),
        denominator(d)
    {}
  };

  double eta_gap_{0.0};
  std::vector<Track> tracks_{};
  mutable std::unordered_map<EventBlockKey, BlockSum, EventBlockKeyHash> event_block_cache_{};
  mutable std::unordered_map<int, std::vector<std::complex<double>>> phase_cache_{};

  void clear_phase_cache() const {
    phase_cache_.clear();
  }

  static int popcount(unsigned x) {
    int count = 0;
    while (x != 0U) {
      count += static_cast<int>(x & 1U);
      x >>= 1U;
    }
    return count;
  }

  static double factorial_int(int n) {
    double out = 1.0;
    for (int i = 2; i <= n; ++i) {
      out *= static_cast<double>(i);
    }
    return out;
  }

  static double block_coefficient(int block_size) {
    const double mag = factorial_int(block_size - 1);
    return ((block_size - 1) % 2 == 0) ? mag : -mag;
  }

  static void build_partitions_recursive(int k,
                                         int idx,
                                         int nblocks,
                                         std::array<int, 6>& assignment,
                                         std::vector<Partition>& out) {
    if (idx == k) {
      std::vector<unsigned> blocks(static_cast<std::size_t>(nblocks), 0U);
      for (int i = 0; i < k; ++i) {
        blocks[static_cast<std::size_t>(assignment[i])] |= (1U << i);
      }
      std::sort(blocks.begin(), blocks.end());

      double coeff = 1.0;
      for (const unsigned block : blocks) {
        coeff *= block_coefficient(popcount(block));
      }
      out.push_back(Partition(blocks, coeff));
      return;
    }

    for (int b = 0; b <= nblocks; ++b) {
      assignment[static_cast<std::size_t>(idx)] = b;
      build_partitions_recursive(k, idx + 1, nblocks + (b == nblocks ? 1 : 0), assignment, out);
    }
  }

  static std::vector<Partition> make_partitions_for(int k) {
    std::vector<Partition> out;
    std::array<int, 6> assignment{};
    assignment.fill(0);
    assignment[0] = 0;
    build_partitions_recursive(k, 1, 1, assignment, out);
    return out;
  }

  static const std::vector<Partition>& partitions_for(int k) {
    static const std::array<std::vector<Partition>, 7> all = []() {
      std::array<std::vector<Partition>, 7> tmp{};
      for (int order = 2; order <= 6; ++order) {
        tmp[static_cast<std::size_t>(order)] = make_partitions_for(order);
      }
      return tmp;
    }();
    return all[static_cast<std::size_t>(k)];
  }

  [[nodiscard]] bool track_in_subevent(const Track& track, Subevent subevent) const {
    switch (subevent) {
    case Subevent::Any:
      return true;
    case Subevent::A:
      return track.eta > (+0.5 * eta_gap_);
    case Subevent::B:
      return track.eta < (-0.5 * eta_gap_);
    }
    return false;
  }

  static Subevent merged_subevent(const std::vector<Leg>& legs,
                                  unsigned block_mask,
                                  bool& compatible) {
    compatible = true;
    Subevent out = Subevent::Any;

    for (std::size_t i = 0; i < legs.size(); ++i) {
      if ((block_mask & (1U << i)) == 0U) {
        continue;
      }
      const Subevent current = legs[i].subevent;
      if (current == Subevent::Any) {
        continue;
      }
      if (out == Subevent::Any) {
        out = current;
        continue;
      }
      if (out != current) {
        compatible = false;
        return Subevent::Any;
      }
    }

    return out;
  }

  const std::vector<std::complex<double>>& phase_factors(int harmonic_sum) const {
    auto it = phase_cache_.find(harmonic_sum);
    if (it != phase_cache_.end()) {
      return it->second;
    }

    std::vector<std::complex<double>> values;
    values.reserve(tracks_.size());
    for (const auto& track : tracks_) {
      const double phase = static_cast<double>(harmonic_sum) * track.phi;
      values.emplace_back(std::cos(phase), std::sin(phase));
    }

    auto inserted = phase_cache_.emplace(harmonic_sum, std::move(values));
    return inserted.first->second;
  }

  static double integer_power(double x, int p) {
    double out = 1.0;
    for (int i = 0; i < p; ++i) {
      out *= x;
    }
    return out;
  }

  std::pair<std::complex<double>, double> block_sums(
      const std::vector<Leg>& legs,
      unsigned block_mask,
      std::array<BlockCacheEntry, 1U << 6>& local_cache) const {
    auto& local_entry = local_cache[block_mask];
    if (local_entry.ready) {
      return std::make_pair(local_entry.numerator, local_entry.denominator);
    }

    bool compatible = true;
    const Subevent block_subevent = merged_subevent(legs, block_mask, compatible);
    if (!compatible) {
      local_entry.ready = true;
      local_entry.numerator = std::complex<double>(0.0, 0.0);
      local_entry.denominator = 0.0;
      return std::make_pair(local_entry.numerator, local_entry.denominator);
    }

    Mask required_categories = 0;
    int harmonic_sum = 0;
    int block_size = 0;
    int gmc_power = 0;

    for (std::size_t i = 0; i < legs.size(); ++i) {
      if ((block_mask & (1U << i)) == 0U) {
        continue;
      }

      required_categories |= legs[i].required_categories;
      harmonic_sum += legs[i].harmonic;
      ++block_size;

      if (legs[i].use_gmc_weight) {
        ++gmc_power;
      }
    }

    const EventBlockKey key{required_categories, harmonic_sum, block_size, gmc_power, block_subevent};
    const auto found = event_block_cache_.find(key);
    if (found != event_block_cache_.end()) {
      local_entry.ready = true;
      local_entry.numerator = found->second.numerator;
      local_entry.denominator = found->second.denominator;
      return std::make_pair(local_entry.numerator, local_entry.denominator);
    }

    std::complex<double> numerator(0.0, 0.0);
    double denominator = 0.0;
    const auto& phases = phase_factors(harmonic_sum);

    for (std::size_t it = 0; it < tracks_.size(); ++it) {
      const auto& track = tracks_[it];
      if ((track.categories & required_categories) != required_categories) {
        continue;
      }
      if (!track_in_subevent(track, block_subevent)) {
        continue;
      }

      const double acc_prod = integer_power(track.acc_weight, block_size);
      const double num_weight = acc_prod * integer_power(track.gmc_weight, gmc_power);

      numerator += num_weight * phases[it];
      denominator += acc_prod;
    }

    local_entry.ready = true;
    local_entry.numerator = numerator;
    local_entry.denominator = denominator;

    event_block_cache_.emplace(key, BlockSum(numerator, denominator));

    return std::make_pair(local_entry.numerator, local_entry.denominator);
  }
};

}  // namespace qcorr
